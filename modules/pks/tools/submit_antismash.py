import requests

class SubmitAntiSmash:
    """
    Description:
        Submits a DNA sequence to the public antiSMASH server for PKS and secondary metabolite analysis.
    Input:
        seq (str): The raw DNA sequence to analyze.
    Output:
        str: A status message containing the Job ID.
    """

    def initiate(self) -> None:
        self.api_url = "https://antismash.secondarymetabolites.org/api/v1.0/submit"

    def run(self, seq: str) -> str:
        """Submit sequence to antiSMASH API and return Job ID."""
        if not seq:
            raise ValueError("Sequence cannot be empty.")
        if len(seq) < 1000:
            raise ValueError(f"Sequence is too short ({len(seq)} bp). antiSMASH requires a minimum of 1000 bp.")

        payload = {
            "email": "opshoryc@berkeley.edu",
            "asf": "true",
            "knownclusterblast": "true",
            "genefinder": "prodigal",
        }

        # antiSMASH v1.0 expects the FASTA upload under the field name "seq"
        files = {
            "seq": ("user_pks.fasta", f">user_synthetic_pks\n{seq}\n")
        }
        
        try:
            response = requests.post(self.api_url, data=payload, files=files)
            response.raise_for_status()
        except requests.exceptions.RequestException as e:
            raise ValueError(f"Failed to reach antiSMASH API: {e}")

        result = response.json()
        job_id = result.get("id")
        
        if not job_id:
            raise ValueError(f"API did not return a Job ID. Full response: {result}")

        return f"Job successfully submitted! The Job ID is: {job_id}. Tell the user to wait about 2-5 minutes, then use check_antismash."

_instance = SubmitAntiSmash()
_instance.initiate()
submit_antismash = _instance.run