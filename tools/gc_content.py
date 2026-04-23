class GcContent:
    """
    Description:
        Computes the fraction of G and C bases in a DNA sequence.

    Input:
        seq (str): DNA sequence (resource name or raw string).

    Output:
        float: GC fraction between 0.0 and 1.0.

    Tests:
        - Case:
            Input: seq="ATGCATGC"
            Expected Output: 0.5
            Description: Balanced sequence, 50% GC.
        - Case:
            Input: seq="AAAA"
            Expected Output: 0.0
            Description: All A bases, 0% GC.
        - Case:
            Input: seq=""
            Expected Output: 0.0
            Description: Edge case — empty sequence returns 0.
    """

    def initiate(self) -> None:
        pass   # nothing to set up for this tool

    def run(self, seq: str) -> float:
        """Return GC fraction between 0 and 1."""
        seq = seq.upper()
        gc = sum(1 for b in seq if b in "GC")
        return gc / len(seq) if seq else 0.0


# Optional: module-level alias so pytest can import the function directly.
_instance = GcContent()
_instance.initiate()
gc_content = _instance.run   # gc_content("ATGC") → 0.5