import requests


class ClusterCADDomainLookup:
    """
    Description:
        Retrieves the amino acid sequence and annotations for a specific PKS
        domain from the ClusterCAD database using its domain ID. Domain IDs
        can be found using the clustercad_get_subunits tool.

    Input:
        domain_id (int): The numeric domain ID from ClusterCAD
                         (e.g. 27717 for an AT domain in Abyssomicin).

    Output:
        dict: Domain details including name, start/stop positions,
              annotations, and amino acid sequence.

    Tests:
        - Case:
            Input: domain_id=27717
            Expected Output: dict with 'name', 'start', 'stop', 'annotations', 'AAsequence'
            Description: Returns details for AT domain 27717 in Abyssomicin.
        - Case:
            Input: domain_id=-1
            Expected Output: ValueError
            Description: Negative domain ID raises ValueError.
        - Case:
            Input: domain_id=0
            Expected Output: ValueError
            Description: Zero domain ID raises ValueError.
    """

    BASE_URL = "https://clustercad.jbei.org"

    def initiate(self) -> None:
        self.session = requests.Session()
        self.session.headers.update({
            "User-Agent": "Mozilla/5.0",
            "X-Requested-With": "XMLHttpRequest",  # required for AJAX endpoint
        })

    def run(self, domain_id: int) -> dict:
        """Return amino acid sequence and details for a specific PKS domain."""

        # handle string inputs from Gemini
        if isinstance(domain_id, str):
            try:
                domain_id = int(domain_id.strip())
            except ValueError:
                raise ValueError(
                    f"domain_id must be an integer, got '{domain_id}'."
                )

        if not isinstance(domain_id, int):
            raise ValueError("domain_id must be an integer.")

        if domain_id <= 0:
            raise ValueError("domain_id must be a positive integer.")

        url = f"{self.BASE_URL}/pks/domainLookup/"
        response = self.session.get(url, params={"domainid": domain_id}, timeout=10)

        if response.status_code == 404:
            raise ValueError(
                f"Domain ID {domain_id} not found in ClusterCAD. "
                "Use clustercad_get_subunits to find valid domain IDs."
            )

        if response.status_code != 200:
            raise ValueError(
                f"ClusterCAD request failed with status code {response.status_code}."
            )

        try:
            data = response.json()
        except Exception:
            raise ValueError(
                f"Could not parse response for domain ID {domain_id}. "
                "The ClusterCAD page structure may have changed."
            )

        return {
            "domain_id"  : domain_id,
            "name"       : data.get("name", ""),
            "start"      : data.get("start", ""),
            "stop"       : data.get("stop", ""),
            "annotations": data.get("annotations", ""),
            "AAsequence" : data.get("AAsequence", ""),
        }