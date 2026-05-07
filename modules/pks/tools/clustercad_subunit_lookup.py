import requests


class ClusterCADSubunitLookup:
    """
    Description:
        Retrieves the amino acid sequence, nucleotide sequence, and details
        for a specific PKS subunit from the ClusterCAD database using its
        subunit ID. Subunit IDs can be found using the clustercad_get_subunits tool.

    Input:
        subunit_id (int): The numeric subunit ID from ClusterCAD
                          (e.g. 24119 for AbsB1 in Abyssomicin).

    Output:
        dict: Subunit details including name, start/stop positions,
              GenBank accession, amino acid sequence, and nucleotide sequence.

    Tests:
        - Case:
            Input: subunit_id=24119
            Expected Output: dict with 'name' containing 'AbsB1'
            Description: Returns details for AbsB1 subunit in Abyssomicin.
        - Case:
            Input: subunit_id=-1
            Expected Output: ValueError
            Description: Negative subunit ID raises ValueError.
        - Case:
            Input: subunit_id=0
            Expected Output: ValueError
            Description: Zero subunit ID raises ValueError.
    """

    BASE_URL = "https://clustercad.jbei.org"

    def initiate(self) -> None:
        self.session = requests.Session()
        self.session.headers.update({
            "User-Agent": "Mozilla/5.0",
            "X-Requested-With": "XMLHttpRequest",  # required for AJAX endpoint
        })

    def run(self, subunit_id: int) -> dict:
        """Return amino acid and nucleotide sequences for a specific PKS subunit."""

        # handle string inputs from Gemini
        if isinstance(subunit_id, str):
            try:
                subunit_id = int(subunit_id.strip())
            except ValueError:
                raise ValueError(
                    f"subunit_id must be an integer, got '{subunit_id}'."
                )

        if not isinstance(subunit_id, int):
            raise ValueError("subunit_id must be an integer.")

        if subunit_id <= 0:
            raise ValueError("subunit_id must be a positive integer.")

        url = f"{self.BASE_URL}/pks/subunitLookup/"
        response = self.session.get(
            url, params={"subunitid": subunit_id}, timeout=10
        )

        if response.status_code == 404:
            raise ValueError(
                f"Subunit ID {subunit_id} not found in ClusterCAD. "
                "Use clustercad_get_subunits to find valid subunit IDs."
            )

        if response.status_code != 200:
            raise ValueError(
                f"ClusterCAD request failed with status code {response.status_code}."
            )

        try:
            data = response.json()
        except Exception:
            raise ValueError(
                f"Could not parse response for subunit ID {subunit_id}. "
                "The ClusterCAD page structure may have changed."
            )

        return {
            "subunit_id"          : subunit_id,
            "name"                : data.get("name", ""),
            "start"               : data.get("start", ""),
            "stop"                : data.get("stop", ""),
            "genbank_accession"   : data.get("genbankAccession", ""),
            "AAsequence"          : data.get("AAsequence", ""),
            "DNAsequence"         : data.get("DNAsequence", ""),
        }