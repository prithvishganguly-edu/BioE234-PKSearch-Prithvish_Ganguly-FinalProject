import os
import requests

class SubmitAntiSmash:
    """
    Description:
        Submits a DNA sequence or GenBank file to the public antiSMASH server for PKS analysis.
    Input:
        seq (str): Raw DNA sequence string to analyze. Must be at least 1000 bp.
        filepath (str): Path to a GenBank (.gb) file to upload directly. When provided,
                        seq is ignored and the CDS annotations are used as-is (no Prodigal).
    Output:
        str: A status message containing the Job ID.
    """

    def initiate(self) -> None:
        self.api_url = "https://antismash.secondarymetabolites.org/api/v1.0/submit"

    def run(self, seq: str = "", filepath: str = "") -> str:
        """Submit a sequence or GenBank file to antiSMASH and return the Job ID."""

        payload = {
            "email": "opshoryc@berkeley.edu",
            "asf": "true",
            "knownclusterblast": "true",
        }

        if filepath:
            # Upload a GenBank file directly — antiSMASH uses existing CDS annotations,
            # no gene prediction needed.
            if not os.path.isfile(filepath):
                raise ValueError(f"File not found: {filepath}")
            ext = os.path.splitext(filepath)[1].lower()
            if ext not in (".gb", ".gbk", ".genbank"):
                raise ValueError(f"filepath must be a GenBank file (.gb/.gbk), got: {ext}")
            with open(filepath, "rb") as f:
                fname = os.path.basename(filepath)
                files = {"seq": (fname, f.read())}
        else:
            # Submit raw DNA as FASTA — use Prodigal for gene prediction.
            if not seq:
                raise ValueError("Provide either seq (DNA string) or filepath (GenBank path).")
            if len(seq) < 1000:
                raise ValueError(
                    f"Sequence is too short ({len(seq)} bp). antiSMASH requires a minimum of 1000 bp."
                )
            payload["genefinder"] = "prodigal"
            files = {"seq": ("user_pks.fasta", f">user_synthetic_pks\n{seq}\n")}

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
