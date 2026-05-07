import requests
from bs4 import BeautifulSoup


class ClusterCADClusterDetails:
    """
    Description:
        Retrieves details about a specific PKS cluster from the ClusterCAD
        database using its MIBiG accession number.

    Input:
        mibig_accession (str): MIBiG accession number of the cluster
                               (e.g. 'BGC0001491.1' or 'BGC0001491').

    Output:
        dict: Cluster details including description, MIBiG accession,
              subunit count, module count, and URL.

    Tests:
        - Case:
            Input: mibig_accession="BGC0001492.1"
            Expected Output: dict with 'description' containing 'Abyssomicin'
            Description: Returns details for the Abyssomicin PKS cluster.
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

    def run(self, mibig_accession: str) -> dict:
        """Return details for a specific PKS cluster from ClusterCAD."""

        if not isinstance(mibig_accession, str) or not mibig_accession.strip():
            raise ValueError("mibig_accession must be a non-empty string.")

        mibig_accession = mibig_accession.strip()

        if not mibig_accession.upper().startswith("BGC"):
            raise ValueError(
                f"Invalid MIBiG accession '{mibig_accession}'. "
                "Accessions should start with 'BGC', e.g. 'BGC0001491.1'."
            )

        # search the cluster list page to find matching row
        list_url = f"{self.BASE_URL}/pks/all/"
        response = self.session.get(list_url, timeout=10)

        if response.status_code != 200:
            raise ValueError(
                f"ClusterCAD request failed with status code {response.status_code}."
            )

        soup = BeautifulSoup(response.text, "html.parser")
        table = soup.find("table", {"id": "clusterTable"})
        if not table:
            raise ValueError(
                "Could not find cluster table. The ClusterCAD page structure may have changed."
            )

        # search for the matching row by accession
        for row in table.select("tbody tr"):
            cols = row.find_all("td")
            if len(cols) < 4:
                continue

            accession = cols[0].get_text(strip=True)

            # match with or without version suffix (e.g. BGC0000055 or BGC0000055.1)
            if accession.upper() == mibig_accession.upper() or \
               accession.upper().split(".")[0] == mibig_accession.upper().split(".")[0]:

                description  = cols[1].get_text(strip=True)
                subunit_count = int(cols[2].get_text(strip=True))
                module_count  = int(cols[3].get_text(strip=True))
                cluster_url   = f"{self.BASE_URL}{row.get('data-href', '')}"

                return {
                    "accession": accession,
                    "description": description,
                    "subunit_count": subunit_count,
                    "module_count": module_count,
                    "url": cluster_url,
                }

        raise ValueError(
            f"Cluster '{mibig_accession}' not found in ClusterCAD. "
            "Use clustercad_list_clusters to find valid accession numbers."
        )