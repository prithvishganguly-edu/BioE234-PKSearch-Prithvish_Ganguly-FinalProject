import requests

class CheckAntiSmash:
    """
    Description:
        Checks the status of an antiSMASH job and parses the results for PKS domains if complete.
    Input:
        job_id (str): The Job ID returned by submit_antismash.
    Output:
        dict: A structured summary of the pathway modules and known cluster hits.
    """

    def initiate(self) -> None:
        self.status_base = "https://antismash.secondarymetabolites.org/api/v1.0/status"
        self.results_base = "https://antismash.secondarymetabolites.org/upload"

    def run(self, job_id: str) -> dict:
        """Check status and parse the resulting JSON for PKS clusters."""
        if not job_id:
            raise ValueError("Job ID cannot be empty.")

        # 1. Check Job Status
        status_url = f"{self.status_base}/{job_id}"
        try:
            status_res = requests.get(status_url)
            status_res.raise_for_status()
        except requests.exceptions.RequestException as e:
            raise ValueError(f"Failed to check status: {e}")
            
        status_data = status_res.json()
        if status_data.get("status") not in ("completed", "done"):
             return {
                 "status": status_data.get("status", "unknown"), 
                 "message": "The job is still processing. Please try again in 1-2 minutes."
             }

        # 2. Fetch Results JSON — filename depends on what was uploaded,
        #    so try <job_id>.json first, then scan the index for any .json file.
        data = None
        for candidate in [f"{job_id}.json", status_data.get("filename", "").replace(".fasta", ".json").replace(".gb", ".json"), "user_pks.json"]:
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
            "pks_clusters": [],
            "domain_predictions": {},
            "predicted_polymer": None,
        }

        records = data.get("records", [])
        if not records:
            return parsed_results

        rec = records[0]
        nrps_pks = rec.get("modules", {}).get("antismash.modules.nrps_pks", {})

        # Always surface domain-level predictions regardless of BGC region calls.
        # This is critical for short constructs (<10 kb) that antiSMASH won't assign
        # to a region but still annotates at the domain level.
        domain_preds = nrps_pks.get("domain_predictions", {})
        for domain_id, pred in domain_preds.items():
            summary = {}
            if "signature" in pred:
                top = sorted(pred["signature"].get("predictions", {}).items(),
                             key=lambda x: x[1][2] if isinstance(x[1], list) else 0,
                             reverse=True)
                if top:
                    summary["AT_substrate"] = top[0][0]
                    summary["AT_confidence"] = round(top[0][1][2], 1) if isinstance(top[0][1], list) else None
            if "kr_stereochem" in pred:
                summary["KR_stereochemistry"] = pred["kr_stereochem"].get("prediction")
            if "kr_activity" in pred:
                summary["KR_activity"] = pred["kr_activity"].get("prediction")
            if summary:
                parsed_results["domain_predictions"][domain_id] = summary

        # Predicted polymer SMILES from region_predictions
        region_preds = nrps_pks.get("region_predictions", {})
        if region_preds:
            first = next(iter(region_preds.values()))
            if first:
                parsed_results["predicted_polymer"] = {
                    "polymer": first[0].get("polymer"),
                    "smiles": first[0].get("smiles"),
                }

        # BGC region clusters (populated for full-length constructs)
        for region in rec.get("regions", []):
            products = region.get("products", [])
            if any("pks" in p.lower() or "nrps" in p.lower() for p in products):
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
                known_hits = region.get("knownclusterblast", {}).get("hits", [])
                for hit in known_hits[:3]:
                    cluster_info["known_cluster_hits"].append({
                        "name": hit.get("name"),
                        "similarity": f"{hit.get('similarity_score', 0)}%",
                    })
                parsed_results["pks_clusters"].append(cluster_info)

        return parsed_results

_instance = CheckAntiSmash()
_instance.initiate()
check_antismash = _instance.run