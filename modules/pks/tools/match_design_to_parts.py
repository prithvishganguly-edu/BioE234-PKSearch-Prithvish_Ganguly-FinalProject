from __future__ import annotations


class MatchDesignToParts:
    """
    Description:
        Takes a PKS design from RetroTide or TridentSynth and automatically
        finds matching natural biological parts in ClusterCAD, returning
        amino acid sequences for each matched domain.

    Input:
        design (dict): A single design object from RetroTide or TridentSynth.
        source (str): "retrotide" or "tridentsynth".
        max_matches_per_module (int): Max ClusterCAD matches per module (default 3).

    Output:
        dict: Module-by-module matches with AA sequences from ClusterCAD.

    Tests:
        - Case:
            Input: design={}, source="retrotide"
            Expected Output: ValueError
            Description: Empty design raises error
        - Case:
            Input: design={"modules": [...]}, source="invalid"
            Expected Output: ValueError
            Description: Invalid source raises error
    """

    _VALID_SOURCES = ("retrotide", "tridentsynth")

    def initiate(self) -> None:
        from modules.pks.tools.clustercad_search_domains import ClusterCADSearchDomains
        from modules.pks.tools.clustercad_domain_lookup import ClusterCADDomainLookup

        self._searcher = ClusterCADSearchDomains()
        self._searcher.initiate()

        self._domain_lookup = ClusterCADDomainLookup()
        self._domain_lookup.initiate()

    # ------------------------------------------------------------------
    # Module normalization
    # ------------------------------------------------------------------

    @staticmethod
    def _normalize_retrotide_modules(design: dict) -> list[dict]:
        raw_modules = design.get("modules")
        if not raw_modules:
            raise ValueError("RetroTide design has no 'modules' key or it is empty.")

        normalized = []
        for mod in raw_modules:
            domains = mod.get("domains", {})
            domain_types = [k for k in domains.keys() if k not in ("ACP", "KS")]

            at_info = domains.get("AT", {})
            at_substrate = at_info.get("substrate") if isinstance(at_info, dict) else None

            normalized.append({
                "loading": mod.get("loading", False),
                "domain_types": domain_types,
                "at_substrate": at_substrate,
                "raw_domains": domains,
            })
        return normalized

    @staticmethod
    def _normalize_tridentsynth_modules(design: dict) -> list[dict]:
        if "pks_modules" in design:
            raw_modules = design["pks_modules"]
        elif "result" in design and "pks_modules" in design.get("result", {}):
            raw_modules = design["result"]["pks_modules"]
        elif isinstance(design, list):
            raw_modules = design
        else:
            raise ValueError(
                "TridentSynth design has no 'pks_modules' key. "
                "Pass the full TridentSynth result dict or the pks_modules list."
            )

        if not raw_modules:
            raise ValueError("TridentSynth design has no modules.")

        normalized = []
        for mod in raw_modules:
            domains_list = mod.get("domains", [])
            domain_types = [
                d["domain"] for d in domains_list
                if d.get("domain") not in ("ACP", "KS", "KSq")
            ]

            at_substrate = None
            for d in domains_list:
                if d.get("domain") == "AT" and d.get("substrate"):
                    at_substrate = d["substrate"]
                    break

            module_type = mod.get("module_type", "")
            is_loading = "starter" in module_type.lower() or "load" in module_type.lower()

            raw_domains = {}
            for d in domains_list:
                raw_domains[d["domain"]] = {"substrate": d.get("substrate")} if d.get("substrate") else {}

            normalized.append({
                "loading": is_loading,
                "domain_types": domain_types,
                "at_substrate": at_substrate,
                "raw_domains": raw_domains,
            })
        return normalized

    # ------------------------------------------------------------------
    # ClusterCAD search + AA fetch
    # ------------------------------------------------------------------

    def _find_parts_for_module(
        self, module: dict, max_matches: int, warnings: list[str], module_index: int,
    ) -> list[dict]:
        search_kwargs = {
            "max_results": max_matches,
            "loading_module_only": module["loading"],
        }

        if module["at_substrate"]:
            search_kwargs["domain_type"] = "AT"
            search_kwargs["annotation_contains"] = module["at_substrate"]

        reductive_domains = [d for d in module["domain_types"] if d in ("KR", "DH", "ER")]
        if reductive_domains:
            search_kwargs["domain_types"] = reductive_domains

        try:
            matches = self._searcher.run(**search_kwargs)
        except Exception as exc:
            warnings.append(f"Module {module_index}: ClusterCAD search failed: {exc}")
            return []

        if not matches:
            warnings.append(f"Module {module_index}: no ClusterCAD matches found.")
            return []

        results = []
        for match in matches:
            domains_with_seq = []
            for dom in match.get("all_domains", []):
                domain_entry = {
                    "domain_type": dom["domain_type"],
                    "domain_id": dom["domain_id"],
                    "annotation": dom.get("annotation", ""),
                    "aa_sequence": None,
                }
                try:
                    lookup = self._domain_lookup.run(dom["domain_id"])
                    domain_entry["aa_sequence"] = lookup.get("AAsequence")
                except Exception as exc:
                    warnings.append(
                        f"Module {module_index}: failed to fetch AA for "
                        f"domain {dom['domain_id']}: {exc}"
                    )
                domains_with_seq.append(domain_entry)

            results.append({
                "accession": match.get("accession", ""),
                "cluster_name": match.get("description", ""),
                "subunit_name": match.get("subunit_name", ""),
                "subunit_id": match.get("subunit_id"),
                "module_label": match.get("module", ""),
                "domains": domains_with_seq,
            })

        return results

    # ------------------------------------------------------------------
    # Main entry point
    # ------------------------------------------------------------------

    def run(
        self,
        design: dict,
        source: str,
        max_matches_per_module: int = 3,
    ) -> dict:
        if not design:
            raise ValueError("design must be a non-empty dict.")

        source = source.strip().lower()
        if source not in self._VALID_SOURCES:
            raise ValueError(
                f"source must be one of {self._VALID_SOURCES}, got {source!r}."
            )

        if not isinstance(max_matches_per_module, int) or max_matches_per_module < 1:
            raise ValueError("max_matches_per_module must be a positive integer.")

        if source == "retrotide":
            normalized = self._normalize_retrotide_modules(design)
        else:
            normalized = self._normalize_tridentsynth_modules(design)

        warnings: list[str] = []
        module_matches = []

        for idx, module in enumerate(normalized):
            natural = self._find_parts_for_module(
                module, max_matches_per_module, warnings, idx,
            )
            module_matches.append({
                "module_index": idx,
                "loading": module["loading"],
                "design_domains": module["raw_domains"],
                "natural_matches": natural,
            })

        modules_with_matches = sum(
            1 for m in module_matches if m["natural_matches"]
        )

        return {
            "source": source,
            "total_modules": len(normalized),
            "modules_with_matches": modules_with_matches,
            "module_matches": module_matches,
            "warnings": warnings,
        }


_instance = MatchDesignToParts()
_instance.initiate()
match_design_to_parts = _instance.run
