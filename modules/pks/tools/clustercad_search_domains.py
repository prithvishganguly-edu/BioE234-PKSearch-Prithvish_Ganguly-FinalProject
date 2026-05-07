import requests
import json
import os
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
                                   annotations (e.g. 'butyryl', 'malonyl',
                                   'loading'). Case insensitive.
        loading_module_only (bool): If True, only return loading modules
                                    (module 0). Defaults to False.
        max_results (int): Maximum number of results to return. Defaults to 10.
        force_refresh (bool): If True, re-downloads all cluster data and
                              rebuilds the cache. Defaults to False.

    Output:
        list: List of matching results, each with cluster accession,
              description, subunit name, module number, and matching domains.

    Tests:
        - Case:
            Input: domain_type="AT", annotation_contains="loading"
            Expected Output: list of dicts with loading AT domains
            Description: Returns modules with loading AT domains.
        - Case:
            Input: domain_type="KR", annotation_contains="B1"
            Expected Output: list of dicts with KR type B1 domains
            Description: Returns modules with KR type B1 domains.
        - Case:
            Input: domain_type=""
            Expected Output: ValueError
            Description: Empty domain type raises ValueError.
    """

    BASE_URL    = "https://clustercad.jbei.org"
    CACHE_FILE  = Path(__file__).parent / "clustercad_cache.json"

    def initiate(self) -> None:
        self.session = requests.Session()
        self.session.headers.update({
            "User-Agent": "Mozilla/5.0"
        })

        # load cache if it exists, otherwise build it
        if self.CACHE_FILE.exists():
            print(
                f"[clustercad_search_domains] Loading cache from {self.CACHE_FILE}",
                flush=True
            )
            with open(self.CACHE_FILE, "r") as f:
                self.cache = json.load(f)
            print(
                f"[clustercad_search_domains] Cache loaded: {len(self.cache)} clusters",
                flush=True
            )
        else:
            print(
                "[clustercad_search_domains] No cache found — building cache now "
                "(this may take ~70 seconds)...",
                flush=True
            )
            self.cache = self._build_cache()

    def run(
        self,
        domain_type           : str,
        annotation_contains   : str  = "",
        loading_module_only   : bool = False,
        max_results           : int  = 10,
        force_refresh         : bool = False,
    ) -> list:
        """Search all PKS clusters for modules matching domain criteria."""

        # handle string inputs from Gemini
        if isinstance(loading_module_only, str):
            loading_module_only = loading_module_only.strip().lower() == "true"

        if isinstance(force_refresh, str):
            force_refresh = force_refresh.strip().lower() == "true"

        if isinstance(max_results, str):
            max_results = int(max_results)

        if not domain_type or not domain_type.strip():
            raise ValueError(
                "domain_type must be a non-empty string (e.g. 'AT', 'KS', 'KR')."
            )

        if not 1 <= max_results <= 531:
            raise ValueError("max_results must be between 1 and 531.")

        # rebuild cache if requested
        if force_refresh:
            print(
                "[clustercad_search_domains] Force refresh requested — "
                "rebuilding cache...",
                flush=True
            )
            self.cache = self._build_cache()

        # synonym dictionary — translates common names to ClusterCAD annotation terms
        synonyms = {
            # malonyl-CoA variants
            "malonyl"           : "mal",
            "malonyl-coa"       : "mal",
            "malonate"          : "mal",

            # methylmalonyl-CoA variants
            "methylmalonyl"     : "mmal",
            "methylmalonyl-coa" : "mmal",
            "methyl malonyl"    : "mmal",

            # ethylmalonyl-CoA variants
            "ethylmalonyl"      : "emal",
            "ethylmalonyl-coa"  : "emal",
            "ethyl malonyl"     : "emal",

            # butylmalonyl-CoA variants
            "butyryl"           : "butmal",
            "butyryl-coa"       : "butmal",
            "butylmalonyl"      : "butmal",
            "butylmalonyl-coa"  : "butmal",
            "butyl malonyl"     : "butmal",

            # hexylmalonyl-CoA variants
            "hexylmalonyl"      : "hxmal",
            "hexylmalonyl-coa"  : "hxmal",
            "hexyl malonyl"     : "hxmal",
            "hexanoyl"          : "hexmal",
            "hexanoyl-coa"      : "hexmal",

            # hydroxymalonyl-CoA variants
            "hydroxymalonyl"    : "hmal",
            "hydroxymalonyl-coa": "hmal",
            "hydroxy malonyl"   : "hmal",

            # methoxymalonyl-CoA variants
            "methoxymalonyl"    : "mxmal",
            "methoxymalonyl-coa": "mxmal",
            "methoxy malonyl"   : "mxmal",

            # isobutyryl-CoA variants
            "isobutyryl"        : "isobut",
            "isobutyryl-coa"    : "isobut",

            # propionyl-CoA variants
            "propionyl"         : "prop",
            "propionyl-coa"     : "prop",
            "propanoyl"         : "prop",

            # acetyl-CoA variants
            "acetyl"            : "Acetyl-CoA",
            "acetyl-coa"        : "Acetyl-CoA",

            # pyruvate variants
            "pyruvate"          : "pyr",
            "pyruvyl"           : "pyr",

            # isobutylmalonyl variants
            "isobutylmalonyl"   : "isobutmal",
            "isobutylmalonyl-coa": "isobutmal",

            # cyclohexane carboxyl variants
            "cyclohexane"       : "CHC-CoA",
            "chc"               : "CHC-CoA",
            "cyclohexanoyl"     : "CHC-CoA",

            # 3-methylhexylmalonyl variants
            "3-methylhexylmalonyl" : "3-me-hexmal",
            "methyl hexyl malonyl" : "3-me-hexmal",

            # 2-oxobutylmalonyl variants
            "2-oxobutylmalonyl" : "2-oxobutmal",
            "oxobutyl"          : "2-oxobutmal",

            # 2-methylbutyryl variants
            "2-methylbutyryl"   : "2metbut",
            "methylbutyryl"     : "2metbut",

            # trans-1,2-CPDA variants
            "cpda"              : "trans-1,2-CPDA",
            "cyclopropane"      : "trans-1,2-CPDA",

            # DCP variants
            "dcp"               : "DCP",
            "dichloropyrrole"   : "DCP",
        }

        # translate annotation if a synonym is found
        annotation_lower = annotation_contains.strip().lower()
        if annotation_lower in synonyms:
            translated = synonyms[annotation_lower]
            print(
                f"[clustercad_search_domains] Translating '{annotation_contains}' "
                f"→ '{translated}'",
                flush=True
            )
            annotation_contains = translated
            annotation_lower    = translated.lower()

        domain_type_upper = domain_type.strip().upper()

        results = []

        for cluster in self.cache:
            if len(results) >= max_results:
                break

            accession   = cluster["accession"]
            description = cluster["description"]

            for subunit in cluster.get("subunits", []):
                if len(results) >= max_results:
                    break

                subunit_name = subunit["subunit_name"]
                subunit_id   = subunit["subunit_id"]

                for module in subunit.get("modules", []):
                    if len(results) >= max_results:
                        break

                    module_num = module["module"]

                    # filter loading modules only
                    if loading_module_only and module_num != "module 0":
                        continue

                    # find matching domains in this module
                    matching_domains = []
                    for domain in module.get("domains", []):
                        if domain["domain_type"].upper() != domain_type_upper:
                            continue

                        # filter by annotation if specified
                        if annotation_lower and \
                           annotation_lower not in domain["annotation"].lower():
                            continue

                        matching_domains.append(domain)

                    if matching_domains:
                        results.append({
                            "accession"       : accession,
                            "description"     : description,
                            "subunit_name"    : subunit_name,
                            "subunit_id"      : subunit_id,
                            "module"          : module_num,
                            "product_smiles"  : module.get("product_smiles", ""),
                            "matching_domains" : matching_domains,
                        })

        if not results:
            raise ValueError(
                f"No modules found with domain type '{domain_type}'"
                + (f" and annotation containing '{annotation_contains}'"
                   if annotation_contains else "")
                + f". Valid annotation terms include: mal, mmal, emal, butmal, "
                  f"hmal, hxmal, mxmal, isobut, prop, pyr, Acetyl-CoA, CHC-CoA, "
                  f"DCP, isobutmal, 3-me-hexmal, 2-oxobutmal, 2metbut, "
                  f"trans-1,2-CPDA. Try broadening your search criteria."
            )

        return results
        
    def _build_cache(self) -> list:
        """Download all cluster architectures and save to cache file."""

        # get full cluster list
        r = self.session.get(
            f"{self.BASE_URL}/pks/all/", timeout=10
        )
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

            # get subunit/domain architecture for this cluster
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

            # print progress every 50 clusters
            if (i + 1) % 50 == 0:
                print(
                    f"[clustercad_search_domains] Progress: {i+1}/{len(rows)} clusters",
                    flush=True
                )

        # save to cache file
        with open(self.CACHE_FILE, "w") as f:
            json.dump(cache, f, indent=2)

        print(
            f"[clustercad_search_domains] Cache built and saved to {self.CACHE_FILE}",
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

                    img = col_div.find("img", attrs={"data-smiles": True})
                    product_smiles = img["data-smiles"] if img else ""

                    if domains:
                        modules.append({
                            "module"         : module_num,
                            "domains"        : domains,
                            "product_smiles" : product_smiles,
                        })

            subunits.append({
                "subunit_name": subunit_name,
                "subunit_id"  : subunit_id,
                "modules"     : modules,
            })

        return subunits