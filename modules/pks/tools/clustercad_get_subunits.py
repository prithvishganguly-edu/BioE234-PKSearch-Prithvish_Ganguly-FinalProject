import requests
from bs4 import BeautifulSoup


class ClusterCADGetSubunits:
    """
    Description:
        Retrieves all subunits and their domain architectures for a specific
        PKS cluster from the ClusterCAD database using its MIBiG accession number.

    Input:
        mibig_accession (str): MIBiG accession number of the cluster
                               (e.g. 'BGC0001492.1' for Abyssomicin).

    Output:
        list: A list of dicts, each representing a subunit with its name,
              subunit ID, and list of modules with domain architectures.

    Tests:
        - Case:
            Input: mibig_accession="BGC0001492.1"
            Expected Output: list of dicts with 'subunit_name', 'subunit_id', 'modules'
            Description: Returns subunits and domain architecture for Abyssomicin.
        - Case:
            Input: mibig_accession=""
            Expected Output: ValueError
            Description: Empty accession raises ValueError.
        - Case:
            Input: mibig_accession="INVALID"
            Expected Output: ValueError
            Description: Invalid accession raises ValueError.
    """

    BASE_URL = "https://clustercad.jbei.org"

    def initiate(self) -> None:
        self.session = requests.Session()
        self.session.headers.update({
            "User-Agent": "Mozilla/5.0"
        })

    def run(self, mibig_accession: str) -> list:
            """Return all subunits and domain architectures for a PKS cluster."""

            if not isinstance(mibig_accession, str) or not mibig_accession.strip():
                raise ValueError("mibig_accession must be a non-empty string.")

            mibig_accession = mibig_accession.strip()

            if not mibig_accession.upper().startswith("BGC"):
                raise ValueError(
                    f"Invalid MIBiG accession '{mibig_accession}'. "
                    "Accessions should start with 'BGC', e.g. 'BGC0001492.1'."
                )

            url = f"{self.BASE_URL}/pks/{mibig_accession}"
            response = self.session.get(url, timeout=10)

            if response.status_code == 404:
                raise ValueError(
                    f"Cluster '{mibig_accession}' not found in ClusterCAD. "
                    "Use clustercad_list_clusters to find valid accession numbers."
                )

            if response.status_code != 200:
                raise ValueError(
                    f"ClusterCAD request failed with status code {response.status_code}."
                )

            soup = BeautifulSoup(response.text, "html.parser")

            subunits = []

            # each subunit has an <a> tag with data-subunitid
            for subunit_tag in soup.find_all("a", attrs={"data-subunitid": True}):
                subunit_id   = subunit_tag["data-subunitid"]
                subunit_name = subunit_tag.get_text(strip=True)

                # get the parent <ul> which contains all modules for this subunit
                parent_ul = subunit_tag.find_parent("ul")
                if not parent_ul:
                    continue

                modules = []

                # each module is a <div> with col classes containing module text
                for li in parent_ul.find_all("li", class_="list-group-item"):
                    # each col div inside the li is a separate module
                    for col_div in li.find_all("div", class_=lambda c: c and "col-md" in c):
                        # get module number
                        module_num = None
                        for text in col_div.stripped_strings:
                            if text.startswith("module"):
                                module_num = text
                                break

                        if not module_num:
                            continue

                        # get all domains in this module from buttons
                        domains = []
                        for btn in col_div.find_all("button", attrs={"data-domainid": True}):
                            domain = {
                                "domain_type": btn.get("data-domain", ""),
                                "domain_id"  : btn.get("data-domainid", ""),
                                "annotation" : btn.get("title", ""),
                            }
                            domains.append(domain)

                        # get predicted product SMILES for this module if available
                        img = col_div.find("img", attrs={"data-smiles": True})
                        product_smiles = img["data-smiles"] if img else ""

                        if domains:
                            modules.append({
                                "module"         : module_num,
                                "domains"        : domains,
                                "product_smiles" : product_smiles,
                            })

                subunits.append({
                    "subunit_name" : subunit_name,
                    "subunit_id"   : subunit_id,
                    "modules"      : modules,
                })

            if not subunits:
                raise ValueError(
                    f"No subunits found for cluster '{mibig_accession}'. "
                    "The ClusterCAD page structure may have changed."
                )

            return subunits