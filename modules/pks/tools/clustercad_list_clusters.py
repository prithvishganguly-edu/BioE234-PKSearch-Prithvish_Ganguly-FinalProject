import requests
from bs4 import BeautifulSoup

#checking
class ClusterCADListClusters:
    """
    Description:
        Retrieves a list of PKS clusters from the ClusterCAD database.

    Input:
        reviewed_only (bool): If True, return only manually reviewed clusters.
                              If False, return all clusters including machine-annotated ones.
        max_results (int): Maximum number of clusters to return. Defaults to 20.

    Output:
        list: A list of dicts, each with 'accession' and 'description' keys.

    Tests:
        - Case:
            Input: reviewed_only=True, max_results=5
            Expected Output: list of up to 5 dicts with 'accession' and 'description'
            Description: Returns up to 5 reviewed PKS clusters.
        - Case:
            Input: reviewed_only=False, max_results=10
            Expected Output: list of up to 10 dicts with 'accession' and 'description'
            Description: Returns up to 10 clusters including unreviewed.
    """

    BASE_URL = "https://clustercad.jbei.org"

    def initiate(self) -> None:
        self.session = requests.Session()
        self.session.headers.update({
            "User-Agent": "Mozilla/5.0"
        })

    def run(self, reviewed_only: bool = True, max_results: int = 20) -> list:
        """Return a list of PKS clusters from ClusterCAD."""

        # handle string inputs from Gemini
        if isinstance(reviewed_only, str):
            if reviewed_only.strip().lower() == "true":
                reviewed_only = True
            elif reviewed_only.strip().lower() == "false":
                reviewed_only = False
            else:
                raise ValueError("reviewed_only must be true or false.")

        if not isinstance(reviewed_only, bool):
            raise ValueError("reviewed_only must be a boolean (True or False).")

        if isinstance(max_results, str):
            max_results = int(max_results)

        if not 1 <= max_results <= 500:
            raise ValueError("max_results must be between 1 and 500.")

        # choose endpoint based on reviewed_only flag
        url = f"{self.BASE_URL}/pks/" if reviewed_only else f"{self.BASE_URL}/pks/all/"

        response = self.session.get(url, timeout=10)

        if response.status_code != 200:
            raise ValueError(
                f"ClusterCAD request failed with status code {response.status_code}."
            )

        soup = BeautifulSoup(response.text, "html.parser")

        # target the specific table by its id
        table = soup.find("table", {"id": "clusterTable"})
        if not table:
            raise ValueError(
                "Could not find cluster table. The ClusterCAD page structure may have changed."
            )

        clusters = []

        for row in table.select("tbody tr"):
            if len(clusters) >= max_results:
                break

            cols = row.find_all("td")
            if len(cols) < 2:
                continue

            # col 0 = MIBiG accession, col 1 = description
            accession = cols[0].get_text(strip=True)
            description = cols[1].get_text(strip=True)

            if accession:
                clusters.append({
                    "accession": accession,
                    "description": description,
                })

        if not clusters:
            raise ValueError(
                "No clusters found. The ClusterCAD page structure may have changed."
            )

        return clusters