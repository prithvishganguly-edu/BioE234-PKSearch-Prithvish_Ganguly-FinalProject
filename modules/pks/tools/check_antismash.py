import re
import time
import requests

# Maps antiSMASH internal substrate codes to human-readable names.
_AT_SUBSTRATE_NAMES = {
    "mal":        "Malonyl-CoA",
    "mmal":       "Methylmalonyl-CoA",
    "emal":       "Ethylmalonyl-CoA",
    "hmal":       "Hydroxymalonyl-CoA",
    "mxmal":      "Methoxymalonyl-CoA",
    "butmal":     "Butylmalonyl-CoA",
    "hxmal":      "Hexylmalonyl-CoA",
    "isobut":     "Isobutyryl-CoA",
    "prop":       "Propionyl-CoA",
    "Acetyl-CoA": "Acetyl-CoA",
    "DHCH":       "DHCH (cyclohexane carboxyl)",
    "trans-1,2-CPDA": "trans-1,2-CPDA",
    "cemal":      "Chloroethylmalonyl-CoA",
    "2metbut":    "2-Methylbutyryl-CoA",
    "CHC-CoA":    "CHC-CoA",
}

# Maps antiSMASH aSDomain type to short name for domain order strings.
_DOMAIN_SHORT = {
    "PKS_KS": "KS", "PKS_AT": "AT", "PKS_KR": "KR",
    "PKS_DH": "DH", "PKS_ER": "ER", "PKS_PP": "ACP",
    "PKS_TE": "TE", "Trans-AT_docking": "docking",
    "AMP-binding": "A",  "PCP": "T", "Epimerization": "E",
    "Condensation_LCL": "C", "Condensation_DCL": "C",
}


def _parse_location_start(loc_str: str) -> int:
    """Extract start position from antiSMASH location string like '[3:1272](+)'."""
    m = re.search(r"\[(\d+):", loc_str)
    return int(m.group(1)) if m else 0


class CheckAntiSmash:
    """
    Description:
        Checks the status of an antiSMASH job and parses the results for PKS domains if complete.
    Input:
        job_id (str): The Job ID returned by submit_antismash.
        wait (bool): If True, poll until the job completes (up to timeout_seconds). Default False.
        timeout_seconds (int): Max seconds to wait when wait=True. Default 300.
        expected_domains (list): Optional list of expected domain short names per gene,
            e.g. [["KS","AT","KR","ACP"]] for a single-gene construct. When provided,
            returns a validation dict comparing expected vs detected.
    Output:
        dict: Structured summary including domain_order, domain_predictions,
              predicted_polymer, mibig_protein_hits, pks_clusters, and
              optionally validation.
    """

    def initiate(self) -> None:
        self.status_base = "https://antismash.secondarymetabolites.org/api/v1.0/status"
        self.results_base = "https://antismash.secondarymetabolites.org/upload"

    def run(self, job_id: str, wait: bool = False, timeout_seconds: int = 300,
            expected_domains: list = None) -> dict:
        """Check status and parse the resulting JSON for PKS domains."""
        if not job_id:
            raise ValueError("Job ID cannot be empty.")

        # 1. Check Job Status (with optional polling loop)
        status_url = f"{self.status_base}/{job_id}"
        deadline = time.time() + timeout_seconds
        poll_interval = 15

        while True:
            try:
                status_res = requests.get(status_url)
                status_res.raise_for_status()
            except requests.exceptions.RequestException as e:
                raise ValueError(f"Failed to check status: {e}")

            status_data = status_res.json()
            current_status = status_data.get("status", "unknown")

            if current_status in ("completed", "done"):
                break
            if current_status.startswith("failed"):
                return {"status": current_status, "message": "The antiSMASH job failed. Check the error message and resubmit."}
            if not wait or time.time() >= deadline:
                return {
                    "status": current_status,
                    "message": f"Job is still processing ({current_status}). Try again in 1-2 minutes, or call with wait=True to poll automatically.",
                }
            time.sleep(poll_interval)

        # 2. Fetch Results JSON
        data = None
        for candidate in [
            f"{job_id}.json",
            status_data.get("filename", "").replace(".fasta", ".json").replace(".gb", ".json"),
            "user_pks.json",
        ]:
            if not candidate or candidate == ".json":
                continue
            try:
                r = requests.get(f"{self.results_base}/{job_id}/{candidate}")
                if r.status_code == 200:
                    data = r.json()
                    break
            except requests.exceptions.RequestException:
                continue
        if data is None:
            return {"error": "Job completed, but failed to fetch results JSON."}

        # 3. Parse PKS Data
        parsed_results = {
            "status": "completed",
            "visualization_url": f"{self.results_base}/{job_id}/index.html",
            "genes": {},           # per-gene domain order + details
            "domain_predictions": {},
            "predicted_polymer": None,
            "pks_clusters": [],
        }

        records = data.get("records", [])
        if not records:
            return parsed_results

        rec = records[0]
        nrps_pks = rec.get("modules", {}).get("antismash.modules.nrps_pks", {})

        # 3a. Build per-gene domain order from aSDomain features (position-sorted).
        #     This covers KS, AT, KR, DH, ER, ACP, TE for all CDS in the record.
        asd_features = [f for f in rec.get("features", []) if f.get("type") == "aSDomain"]
        asd_features.sort(key=lambda f: _parse_location_start(f.get("location", "")))

        domain_pred_map = nrps_pks.get("domain_predictions", {})

        for feat in asd_features:
            q = feat.get("qualifiers", {})
            domain_type = q.get("aSDomain", ["?"])[0]
            domain_id   = q.get("domain_id", [""])[0]
            locus_tag   = q.get("locus_tag", ["unknown"])[0]
            short       = _DOMAIN_SHORT.get(domain_type, domain_type)
            evalue      = q.get("evalue", [None])[0]
            score       = q.get("score", [None])[0]
            subtype     = q.get("domain_subtypes", [""])[0]

            if locus_tag not in parsed_results["genes"]:
                parsed_results["genes"][locus_tag] = {
                    "domain_order": [],
                    "domain_details": [],
                }

            detail = {
                "domain": short,
                "domain_type": domain_type,
                "domain_id": domain_id,
                "location": feat.get("location"),
                "evalue": evalue,
                "score": score,
            }
            if subtype:
                detail["subtype"] = subtype

            # Attach AT/KR predictions if available
            pred = domain_pred_map.get(domain_id, {})
            if "signature" in pred:
                top = sorted(
                    pred["signature"].get("predictions", {}).items(),
                    key=lambda x: x[1][2] if isinstance(x[1], list) else 0,
                    reverse=True,
                )
                if top:
                    raw = top[0][0]
                    detail["AT_substrate"] = _AT_SUBSTRATE_NAMES.get(raw, raw)
                    detail["AT_substrate_code"] = raw
                    detail["AT_confidence"] = round(top[0][1][2], 1) if isinstance(top[0][1], list) else None
            if "kr_stereochem" in pred:
                detail["KR_stereochemistry"] = pred["kr_stereochem"].get("prediction")
            if "kr_activity" in pred:
                detail["KR_activity"] = pred["kr_activity"].get("prediction")

            parsed_results["genes"][locus_tag]["domain_order"].append(short)
            parsed_results["genes"][locus_tag]["domain_details"].append(detail)

            # Also populate flat domain_predictions for backward compatibility
            if any(k in detail for k in ("AT_substrate", "KR_stereochemistry")):
                entry = {k: v for k, v in detail.items()
                         if k in ("AT_substrate", "AT_substrate_code", "AT_confidence",
                                  "KR_stereochemistry", "KR_activity")}
                parsed_results["domain_predictions"][domain_id] = entry

        # Build domain_order string per gene
        for gene_data in parsed_results["genes"].values():
            gene_data["domain_order_string"] = "-".join(gene_data["domain_order"])

        # 3b. Expected vs actual validation
        if expected_domains is not None:
            gene_list = list(parsed_results["genes"].values())
            validation = []
            for i, expected in enumerate(expected_domains):
                if i < len(gene_list):
                    detected = gene_list[i]["domain_order"]
                else:
                    detected = []
                exp_set = set(expected)
                det_set = set(detected)
                validation.append({
                    "expected": expected,
                    "detected": detected,
                    "domain_order_string": "-".join(detected),
                    "missing": sorted(exp_set - det_set),
                    "unexpected": sorted(det_set - exp_set),
                    "match": expected == detected,
                })
            parsed_results["validation"] = validation

        # 3c. Predicted polymer SMILES
        region_preds = nrps_pks.get("region_predictions", {})
        if region_preds:
            first = next(iter(region_preds.values()))
            if first:
                parsed_results["predicted_polymer"] = {
                    "polymer": first[0].get("polymer"),
                    "smiles": first[0].get("smiles"),
                }

        # 3d. KnownClusterBlast — protein-level MIBiG matches (always present)
        clusterblast = rec.get("modules", {}).get("antismash.modules.clusterblast", {})
        knowncluster = clusterblast.get("knowncluster", {})
        mibig_protein_hits = []
        for _region, cds_dict in knowncluster.get("mibig_entries", {}).items():
            for gene_id, matches in cds_dict.items():
                for match in matches[:3]:
                    if len(match) >= 6:
                        mibig_protein_hits.append({
                            "gene": gene_id,
                            "protein_accession": match[0],
                            "protein_name": match[1],
                            "bgc_accession": match[2],
                            "product_type": match[4],
                            "similarity_pct": match[5],
                        })
        if mibig_protein_hits:
            parsed_results["mibig_protein_hits"] = sorted(
                mibig_protein_hits, key=lambda x: x["similarity_pct"], reverse=True
            )[:5]

        # 3e. BGC region clusters (≥10 kb constructs only)
        kcb_results = {r["region_number"]: r for r in knowncluster.get("results", [])}
        for region in rec.get("regions", []):
            products = region.get("products", [])
            if any("pks" in p.lower() or "nrps" in p.lower() for p in products):
                region_num = region.get("region_number", region.get("idx", 1))
                cluster_info = {
                    "type": products,
                    "location": f"Base pairs {region.get('start')} to {region.get('end')}",
                    "pathway_modules": [],
                    "known_cluster_hits": [],
                }
                for mod in nrps_pks.get("modules", []):
                    if mod.get("start", 0) >= region.get("start", 0) and mod.get("end", 0) <= region.get("end", 0):
                        cluster_info["pathway_modules"].append({
                            "module_number": mod.get("module_number", "Unknown"),
                            "domains": [d.get("name") for d in mod.get("domains", [])],
                            "predicted_substrate": mod.get("predictions", {}).get("specificity", "Unknown"),
                        })
                for entry in kcb_results.get(region_num, {}).get("ranking", [])[:3]:
                    if len(entry) >= 2:
                        cdict, sdict = entry[0], entry[1]
                        cluster_info["known_cluster_hits"].append({
                            "bgc_accession": cdict.get("accession"),
                            "description": cdict.get("description"),
                            "similarity_pct": sdict.get("similarity"),
                            "hits": sdict.get("hits"),
                        })
                parsed_results["pks_clusters"].append(cluster_info)

        return parsed_results


_instance = CheckAntiSmash()
_instance.initiate()
check_antismash = _instance.run
