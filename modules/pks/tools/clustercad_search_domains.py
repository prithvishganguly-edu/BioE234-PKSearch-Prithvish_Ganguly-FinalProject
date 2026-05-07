import requests
import json
from bs4 import BeautifulSoup
from pathlib import Path


class ClusterCADSearchDomains:
    """
    Description:
        Searches all PKS clusters in the ClusterCAD database for modules
        containing domains matching specified criteria. Builds a local cache
        of all cluster architectures on first run for fast subsequent searches.

    Input:
        domain_type (str): Domain type to search for (e.g. 'AT', 'KS', 'KR',
                           'DH', 'ER', 'ACP'). Case insensitive.
        annotation_contains (str): Optional text to search for in domain
                                   annotations. Case insensitive.
        domain_types (list): Optional list of domain types that must ALL be
                             present in the module (e.g. ['KR', 'DH', 'ER']).
        exclude_annotation (str): Optional text to exclude domains with this
                                  annotation (e.g. 'inactive').
        active_only (bool): If True, only return active domains.
        loading_module_only (bool): If True, only return loading modules.
        cluster_description_contains (str): Filter by cluster name.
        min_modules (int): Only return clusters with at least this many modules.
        max_modules (int): Only return clusters with at most this many modules.
        reviewed_only (bool): If True, only search reviewed clusters.
        max_results (int): Maximum number of results to return. Defaults to 10.
        force_refresh (bool): If True, rebuilds the cache.

    Output:
        list: List of matching results with cluster, subunit, module, and
              domain details.

    Tests:
        - Case:
            Input: domain_type="AT", annotation_contains="loading"
            Expected Output: list of dicts with loading AT domains
            Description: Returns modules with loading AT domains.
        - Case:
            Input: domain_types=["KR", "DH"], active_only=True
            Expected Output: list of dicts with active KR and DH domains
            Description: Returns modules with both active KR and DH domains.
        - Case:
            Input: domain_type="", domain_types=[]
            Expected Output: ValueError
            Description: No search criteria raises ValueError.
    """

    BASE_URL     = "https://clustercad.jbei.org"
    CACHE_FILE   = Path(__file__).parent / "clustercad_cache.json"
    REVIEWED_URL = "https://clustercad.jbei.org/pks/"

    def initiate(self) -> None:
        self.session = requests.Session()
        self.session.headers.update({
            "User-Agent": "Mozilla/5.0"
        })

        if self.CACHE_FILE.exists():
            print(
                f"[clustercad_search_domains] Loading cache from {self.CACHE_FILE}",
                flush=True
            )
            with open(self.CACHE_FILE, "r") as f:
                self.cache = json.load(f)
            print(
                f"[clustercad_search_domains] Cache loaded: "
                f"{len(self.cache)} clusters",
                flush=True
            )
        else:
            print(
                "[clustercad_search_domains] No cache found — building cache now "
                "(this may take ~70 seconds)...",
                flush=True
            )
            self.cache = self._build_cache()

        self.reviewed_accessions = self._get_reviewed_accessions()

    def run(
        self,
        domain_type                  : str  = "",
        annotation_contains          : str  = "",
        domain_types                 : list = None,
        exclude_annotation           : str  = "",
        active_only                  : bool = False,
        loading_module_only          : bool = False,
        cluster_description_contains : str  = "",
        min_modules                  : int  = 0,
        max_modules                  : int  = 9999,
        reviewed_only                : bool = False,
        max_results                  : int  = 10,
        force_refresh                : bool = False,
    ) -> list:
        """Search all PKS clusters for modules matching domain criteria."""

        # handle string inputs from Gemini
        if isinstance(active_only, str):
            active_only = active_only.strip().lower() == "true"
        if isinstance(loading_module_only, str):
            loading_module_only = loading_module_only.strip().lower() == "true"
        if isinstance(reviewed_only, str):
            reviewed_only = reviewed_only.strip().lower() == "true"
        if isinstance(force_refresh, str):
            force_refresh = force_refresh.strip().lower() == "true"
        if isinstance(max_results, str):
            max_results = int(max_results)
        if isinstance(min_modules, str):
            min_modules = int(min_modules)
        if isinstance(max_modules, str):
            max_modules = int(max_modules)
        if domain_types is None:
            domain_types = []

        # require at least one search criterion
        if not domain_type.strip() and not domain_types:
            raise ValueError(
                "Provide at least one of: domain_type (str) or "
                "domain_types (list of str)."
            )

        if not 1 <= max_results <= 531:
            raise ValueError("max_results must be between 1 and 531.")

        if force_refresh:
            print(
                "[clustercad_search_domains] Force refresh — rebuilding cache...",
                flush=True
            )
            self.cache = self._build_cache()
            self.reviewed_accessions = self._get_reviewed_accessions()

        # synonym dictionary
        synonyms = {
            "malonyl"              : "mal",
            "malonyl-coa"          : "mal",
            "malonate"             : "mal",
            "methylmalonyl"        : "mmal",
            "methylmalonyl-coa"    : "mmal",
            "methyl malonyl"       : "mmal",
            "ethylmalonyl"         : "emal",
            "ethylmalonyl-coa"     : "emal",
            "ethyl malonyl"        : "emal",
            "butyryl"              : "butmal",
            "butyryl-coa"          : "butmal",
            "butylmalonyl"         : "butmal",
            "butylmalonyl-coa"     : "butmal",
            "butyl malonyl"        : "butmal",
            "hexylmalonyl"         : "hxmal",
            "hexylmalonyl-coa"     : "hxmal",
            "hexyl malonyl"        : "hxmal",
            "hexanoyl"             : "hexmal",
            "hexanoyl-coa"         : "hexmal",
            "hydroxymalonyl"       : "hmal",
            "hydroxymalonyl-coa"   : "hmal",
            "hydroxy malonyl"      : "hmal",
            "methoxymalonyl"       : "mxmal",
            "methoxymalonyl-coa"   : "mxmal",
            "methoxy malonyl"      : "mxmal",
            "isobutyryl"           : "isobut",
            "isobutyryl-coa"       : "isobut",
            "propionyl"            : "prop",
            "propionyl-coa"        : "prop",
            "propanoyl"            : "prop",
            "acetyl"               : "Acetyl-CoA",
            "acetyl-coa"           : "Acetyl-CoA",
            "pyruvate"             : "pyr",
            "pyruvyl"              : "pyr",
            "isobutylmalonyl"      : "isobutmal",
            "isobutylmalonyl-coa"  : "isobutmal",
            "cyclohexane"          : "CHC-CoA",
            "chc"                  : "CHC-CoA",
            "cyclohexanoyl"        : "CHC-CoA",
            "3-methylhexylmalonyl" : "3-me-hexmal",
            "methyl hexyl malonyl" : "3-me-hexmal",
            "2-oxobutylmalonyl"    : "2-oxobutmal",
            "oxobutyl"             : "2-oxobutmal",
            "2-methylbutyryl"      : "2metbut",
            "methylbutyryl"        : "2metbut",
            "cpda"                 : "trans-1,2-CPDA",
            "cyclopropane"         : "trans-1,2-CPDA",
            "dcp"                  : "DCP",
            "dichloropyrrole"      : "DCP",
        }

        # translate annotation synonyms
        annotation_lower = annotation_contains.strip().lower()
        if annotation_lower in synonyms:
            translated = synonyms[annotation_lower]
            print(
                f"[clustercad_search_domains] Translating "
                f"'{annotation_contains}' → '{translated}'",
                flush=True
            )
            annotation_contains = translated
            annotation_lower    = translated.lower()

        # translate exclude_annotation synonyms
        exclude_lower = exclude_annotation.strip().lower()
        if exclude_lower in synonyms:
            exclude_annotation = synonyms[exclude_lower]
            exclude_lower      = exclude_annotation.lower()

        domain_type_upper  = domain_type.strip().upper()
        domain_types_upper = [d.strip().upper() for d in domain_types]
        desc_lower         = cluster_description_contains.strip().lower()

        results = []

        for cluster in self.cache:
            if len(results) >= max_results:
                break

            accession   = cluster["accession"]
            description = cluster["description"]

            # filter by reviewed only
            if reviewed_only and accession not in self.reviewed_accessions:
                continue

            # filter by cluster description
            if desc_lower and desc_lower not in description.lower():
                continue

            # count total modules in this cluster
            total_modules = sum(
                len(s.get("modules", []))
                for s in cluster.get("subunits", [])
            )

            # filter by module count
            if total_modules < min_modules or total_modules > max_modules:
                continue

            for subunit in cluster.get("subunits", []):
                if len(results) >= max_results:
                    break

                subunit_name = subunit["subunit_name"]
                subunit_id   = subunit["subunit_id"]

                for module in subunit.get("modules", []):
                    if len(results) >= max_results:
                        break

                    module_num  = module["module"]
                    all_domains = module.get("domains", [])

                    # filter loading modules only
                    if loading_module_only and module_num != "module 0":
                        continue

                    # filter by domain_types — ALL must be present
                    if domain_types_upper:
                        present = {d["domain_type"].upper() for d in all_domains}
                        if not all(dt in present for dt in domain_types_upper):
                            continue

                    # find matching domains for primary domain_type search
                    if domain_type_upper:
                        matching_domains = []
                        for domain in all_domains:
                            if domain["domain_type"].upper() != domain_type_upper:
                                continue
                            if annotation_lower and \
                               annotation_lower not in domain["annotation"].lower():
                                continue
                            if exclude_lower and \
                               exclude_lower in domain["annotation"].lower():
                                continue
                            if active_only and \
                               "inactive" in domain["annotation"].lower():
                                continue
                            matching_domains.append(domain)

                        if not matching_domains:
                            continue
                    else:
                        # no primary domain_type — return all domains
                        matching_domains = all_domains

                    results.append({
                        "accession"        : accession,
                        "description"      : description,
                        "subunit_name"     : subunit_name,
                        "subunit_id"       : subunit_id,
                        "module"           : module_num,
                        "total_modules"    : total_modules,
                        "product_smiles"   : module.get("product_smiles", ""),
                        "all_domains"      : all_domains,
                        "matching_domains" : matching_domains,
                    })

        if not results:
            raise ValueError(
                f"No modules found matching your criteria. "
                f"Valid annotation terms: mal, mmal, emal, butmal, hmal, hxmal, "
                f"mxmal, isobut, prop, pyr, Acetyl-CoA, CHC-CoA, DCP, isobutmal, "
                f"3-me-hexmal, 2-oxobutmal, 2metbut, trans-1,2-CPDA. "
                f"Try broadening your search criteria."
            )

        return results

    def _get_reviewed_accessions(self) -> set:
        """Fetch the set of reviewed cluster accessions."""
        try:
            r     = self.session.get(self.REVIEWED_URL, timeout=10)
            soup  = BeautifulSoup(r.text, "html.parser")
            table = soup.find("table", {"id": "clusterTable"})
            if not table:
                return set()
            return {
                row.find_all("td")[0].get_text(strip=True)
                for row in table.select("tbody tr")
                if len(row.find_all("td")) >= 1
            }
        except Exception:
            return set()

    def _build_cache(self) -> list:
        """Download all cluster architectures and save to cache file."""

        r = self.session.get(f"{self.BASE_URL}/pks/all/", timeout=10)
        if r.status_code != 200:
            raise ValueError(
                f"Could not fetch cluster list: status {r.status_code}"
            )

        soup  = BeautifulSoup(r.text, "html.parser")
        table = soup.find("table", {"id": "clusterTable"})
        if not table:
            raise ValueError("Could not find cluster table.")

        rows = table.select("tbody tr")
        print(
            f"[clustercad_search_domains] Downloading {len(rows)} clusters...",
            flush=True
        )

        cache = []
        for i, row in enumerate(rows):
            cols = row.find_all("td")
            if len(cols) < 2:
                continue

            accession   = cols[0].get_text(strip=True)
            description = cols[1].get_text(strip=True)

            try:
                subunits = self._get_subunits(accession)
            except Exception as e:
                print(
                    f"[clustercad_search_domains] WARNING: Could not fetch "
                    f"{accession}: {e}",
                    flush=True
                )
                subunits = []

            cache.append({
                "accession"  : accession,
                "description": description,
                "subunits"   : subunits,
            })

            if (i + 1) % 50 == 0:
                print(
                    f"[clustercad_search_domains] Progress: "
                    f"{i+1}/{len(rows)} clusters",
                    flush=True
                )

        with open(self.CACHE_FILE, "w") as f:
            json.dump(cache, f, indent=2)

        print(
            f"[clustercad_search_domains] Cache built and saved to "
            f"{self.CACHE_FILE}",
            flush=True
        )

        return cache

    def _get_subunits(self, mibig_accession: str) -> list:
            """Fetch subunit/domain architecture for a single cluster."""

            url      = f"{self.BASE_URL}/pks/{mibig_accession}"
            response = self.session.get(url, timeout=10)

            if response.status_code != 200:
                return []

            soup     = BeautifulSoup(response.text, "html.parser")
            subunits = []

            for subunit_tag in soup.find_all("a", attrs={"data-subunitid": True}):
                subunit_id   = subunit_tag["data-subunitid"]
                subunit_name = subunit_tag.get_text(strip=True)
                parent_ul    = subunit_tag.find_parent("ul")

                if not parent_ul:
                    continue

                modules = []
                for li in parent_ul.find_all("li", class_="list-group-item"):
                    for col_div in li.find_all(
                        "div", class_=lambda c: c and "col-md" in c
                    ):
                        module_num = None
                        for text in col_div.stripped_strings:
                            if text.startswith("module"):
                                module_num = text
                                break

                        if not module_num:
                            continue

                        domains = []
                        for btn in col_div.find_all(
                            "button", attrs={"data-domainid": True}
                        ):
                            domains.append({
                                "domain_type": btn.get("data-domain", ""),
                                "domain_id"  : btn.get("data-domainid", ""),
                                "annotation" : btn.get("title", ""),
                            })

                        img            = col_div.find("img", attrs={"data-smiles": True})
                        product_smiles = img["data-smiles"] if img else ""

                        if domains:
                            modules.append({
                                "module"        : module_num,
                                "domains"       : domains,
                                "product_smiles": product_smiles,
                            })

                subunits.append({
                    "subunit_name": subunit_name,
                    "subunit_id"  : subunit_id,
                    "modules"     : modules,
                })

            return subunits