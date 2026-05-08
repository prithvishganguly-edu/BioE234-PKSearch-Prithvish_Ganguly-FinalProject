import os
import requests

class SubmitAntiSmash:
    """
    Description:
        Submits a DNA sequence, GenBank file, or NCBI accession to the public
        antiSMASH server for PKS domain analysis. Provide exactly one of the
        three input parameters. Active Site Finder and KnownClusterBlast are
        always enabled.
    Input:
        seq (str): Raw DNA sequence string (≥1000 bp). Prodigal used for gene
                   prediction. Mutually exclusive with filepath and ncbi.
        filepath (str): Path to a GenBank (.gb) file. Uses existing CDS
                        annotations directly — no Prodigal. Preferred for
                        in-silico designs from reverse_translate.
        ncbi (str): NCBI nucleotide accession (e.g. 'AM420293', 'NC_003888').
                    antiSMASH fetches the record server-side. Best for
                    sequenced clones deposited on NCBI.
    Output:
        str: A status message containing the Job ID for use with check_antismash.
    """

    def initiate(self) -> None:
        self.api_url = "https://antismash.secondarymetabolites.org/api/v1.0/submit"

    def run(self, seq: str = "", filepath: str = "", ncbi: str = "") -> str:
        """Submit to antiSMASH and return the Job ID."""

        # Validate — exactly one input must be provided
        provided = sum(bool(x) for x in (seq, filepath, ncbi))
        if provided == 0:
            raise ValueError("Provide one of: seq (DNA string), filepath (GenBank path), or ncbi (accession).")
        if provided > 1:
            raise ValueError("Provide only one of seq, filepath, or ncbi — not multiple.")

        payload = {
            "email": "opshoryc@berkeley.edu",
            "asf": "true",
            "knownclusterblast": "true",
        }

        if ncbi:
            # Let antiSMASH fetch the record directly from NCBI.
            # NCBI records already have CDS annotations — no gene finder needed.
            payload["ncbi"] = ncbi.strip()
            try:
                response = requests.post(self.api_url, data=payload)
                response.raise_for_status()
            except requests.exceptions.RequestException as e:
                raise ValueError(f"Failed to reach antiSMASH API: {e}")

        elif filepath:
            # Upload a GenBank file — uses existing CDS annotations, no Prodigal.
            ext = os.path.splitext(filepath)[1].lower()
            if ext not in (".gb", ".gbk", ".genbank"):
                raise ValueError(f"filepath must be a GenBank file (.gb/.gbk), got: {ext}")
            if not os.path.isfile(filepath):
                raise ValueError(f"File not found: {filepath}")
            with open(filepath, "rb") as f:
                fname = os.path.basename(filepath)
                files = {"seq": (fname, f.read())}
            try:
                response = requests.post(self.api_url, data=payload, files=files)
                response.raise_for_status()
            except requests.exceptions.RequestException as e:
                raise ValueError(f"Failed to reach antiSMASH API: {e}")

        else:
            # Submit raw DNA as FASTA — use Prodigal for gene prediction.
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
