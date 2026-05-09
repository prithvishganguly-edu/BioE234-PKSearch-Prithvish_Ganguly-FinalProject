from __future__ import annotations


_VALID_REDUCTIVE = {"KR", "DH", "ER"}


class FindPKSModuleParts:
    """
    Description:
        Finds matching natural PKS modules in the ClusterCAD database for
        a single module with specified domain properties.  Call this once
        per module in a design.  Accepts simple scalar inputs (loading
        status, AT substrate, reduction domains present) and returns
        cluster matches with amino acid sequences for all domains.

    Input:
        loading            (bool): Whether this is a loading/starter module.
        at_substrate       (str):  AT substrate name, e.g. "Malonyl-CoA"
                                   or "" if unspecified.
        reductive_domains  (str):  Comma-separated active reduction domains
                                   present, e.g. "KR,DH,ER" or "KR" or "".
        max_matches        (int):  Max ClusterCAD matches to return (default 3).

    Output:
        dict with keys:
            - "loading"            (bool)
            - "at_substrate"       (str | None)
            - "reductive_domains"  (list[str])
            - "total_matches"      (int)
            - "matches"            (list[dict])
            - "warnings"           (list[str])

    Tests:
        - Case:
            Input: loading=True, at_substrate="Malonyl-CoA"
            Expected Output: dict with total_matches >= 0 and matches list
            Description: Search for loading modules with Malonyl-CoA AT.
        - Case:
            Input: loading=False, reductive_domains="KR,DH,ER"
            Expected Output: dict with matches containing domains with aa_sequence
            Description: Search for extension modules with full reductive loop.
        - Case:
            Input: loading=True, max_matches=0
            Expected Output: ValueError raised
            Description: max_matches must be >= 1.
        - Case:
            Input: loading=True, reductive_domains="ZZ"
            Expected Output: ValueError raised
            Description: Unrecognized domain type.
    """

    def initiate(self) -> None:
        from modules.pks.tools.clustercad_search_domains import ClusterCADSearchDomains
        from modules.pks.tools.clustercad_domain_lookup import ClusterCADDomainLookup

        self._searcher = ClusterCADSearchDomains()
        self._searcher.initiate()

        self._domain_lookup = ClusterCADDomainLookup()
        self._domain_lookup.initiate()

    @staticmethod
    def _parse_reductive_domains(raw: str) -> list[str]:
        if not raw or not raw.strip():
            return []
        parts = [d.strip().upper() for d in raw.split(",") if d.strip()]
        for d in parts:
            if d not in _VALID_REDUCTIVE:
                raise ValueError(
                    f"Unrecognized reductive domain {d!r}. "
                    f"Valid values: {sorted(_VALID_REDUCTIVE)}"
                )
        return parts

    def run(
        self,
        loading: bool,
        at_substrate: str = "",
        reductive_domains: str = "",
        max_matches: int = 3,
    ) -> dict:
        if not isinstance(max_matches, int) or max_matches < 1:
            raise ValueError("max_matches must be a positive integer >= 1.")

        parsed_domains = self._parse_reductive_domains(reductive_domains)
        at_sub = at_substrate.strip() if at_substrate else None

        search_kwargs = {
            "max_results": max_matches,
            "loading_module_only": bool(loading),
        }

        if at_sub:
            search_kwargs["domain_type"] = "AT"
            search_kwargs["annotation_contains"] = at_sub

        if parsed_domains:
            search_kwargs["domain_types"] = parsed_domains

        warnings: list[str] = []

        try:
            raw_matches = self._searcher.run(**search_kwargs)
        except Exception as exc:
            warnings.append(f"ClusterCAD search failed: {exc}")
            raw_matches = []

        matches = []
        for match in raw_matches:
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
                        f"Failed to fetch AA for domain {dom['domain_id']}: {exc}"
                    )
                domains_with_seq.append(domain_entry)

            matches.append({
                "accession": match.get("accession", ""),
                "cluster_name": match.get("description", ""),
                "subunit_name": match.get("subunit_name", ""),
                "subunit_id": match.get("subunit_id"),
                "module_label": match.get("module", ""),
                "domains": domains_with_seq,
            })

        return {
            "loading": bool(loading),
            "at_substrate": at_sub,
            "reductive_domains": parsed_domains,
            "total_matches": len(matches),
            "matches": matches,
            "warnings": warnings,
        }


_instance = FindPKSModuleParts()
_instance.initiate()
find_pks_module_parts = _instance.run
