from __future__ import annotations

"""
MCP-callable wrapper for the RetroTide PKS retrobiosynthesis tool.

Place this file at:
    modules/<your_module>/tools/retrotide.py
"""


class Retrotide:
    """
    Description:
        Given a target chemical structure as a SMILES string, uses the
        RetroTide algorithm to propose one or more chimeric type-I modular
        PKS designs that can synthesise the compound's carbon backbone
        without any post-PKS enzymatic decoration.  Each design is ranked by
        chemical similarity (Tanimoto on Morgan fingerprints) between the
        predicted PKS product and the target.

        Use this tool when the user asks for a product that can be made
        *purely* by a PKS — i.e. no downstream tailoring enzymes are needed.

    Input:
        smiles       (str):  SMILES string for the target polyketide molecule.
        max_designs  (int):  Maximum number of PKS designs to return
                             (default 5, max 20).
        similarity   (str):  Fingerprint metric for ranking.  Choices:
                             "tanimoto" (default), "dice".

    Output:
        list[dict]: Ordered list (best first) of PKS designs.  Each dict has:
            - "rank"            (int)   : 1-based rank (1 = most similar)
            - "similarity"      (float) : Tanimoto / Dice score vs. target
            - "product_smiles"  (str)   : SMILES of the PKS product
            - "modules"         (list)  : Ordered list of module descriptors.
              Each module dict contains:
                  "starter"     (bool)         : True only for the loading module
                  "domains"     (list[str])     : Active domain labels, e.g.
                                                 ["KS","AT","KR","ACP"]
                  "extender"    (str | None)   : Extender unit CoA name, e.g.
                                                 "malonyl-CoA"
                  "at_specificity" (str | None): AT domain substrate
            - "exact_match"     (bool)  : True when product_smiles equals target

    Tests:
        - Case:
            Input: smiles="CCCC(=O)O", max_designs=3, similarity="tanimoto"
            Expected Output: list with 1-3 dicts each having "rank", "similarity",
                             "product_smiles", "modules", "exact_match"
            Description: Short-chain acid; RetroTide returns at least one design.

        - Case:
            Input: smiles="", max_designs=5, similarity="tanimoto"
            Expected Output: ValueError raised
            Description: Empty SMILES should raise ValueError.

        - Case:
            Input: smiles="CCCC(=O)O", max_designs=0, similarity="tanimoto"
            Expected Output: ValueError raised
            Description: max_designs must be >= 1.
    """

    # ------------------------------------------------------------------ #
    #  Allowed similarity metrics accepted by RetroTide's designPKS()      #
    # ------------------------------------------------------------------ #
    _VALID_METRICS = {"tanimoto", "dice"}
    _MAX_DESIGNS_HARD_LIMIT = 20

    def initiate(self) -> None:
        self._import_error = None
        self._designPKS = None
        self._MolFromSmiles = None
        try:
            from retrotide.retrotide import designPKS
            from rdkit.Chem import MolFromSmiles
            self._designPKS = designPKS
            self._MolFromSmiles = MolFromSmiles
        except ImportError as exc:
            self._import_error = exc  # store, do NOT raise`

    # ------------------------------------------------------------------ #
    def run(
        self,
        smiles: str,
        max_designs: int = 5,
        similarity: str = "tanimoto",
    ) -> list[dict]:
        """
        Call RetroTide and return serialisable PKS design objects.

        Parameters
        ----------
        smiles : str
            SMILES of the target polyketide (must be non-empty and parseable
            by RDKit).
        max_designs : int
            How many top-ranked designs to return (1 – 20).
        similarity : str
            Ranking metric: "tanimoto" or "dice".

        Returns
        -------
        list[dict]
            JSON-serialisable list of PKS design dicts, best first.
        """
        # ---- Input validation ---------------------------------------- #
        smiles = smiles.strip()
        if not smiles:
            raise ValueError("smiles must be a non-empty SMILES string.")

        mol = self._MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(
                f"RDKit could not parse the provided SMILES: {smiles!r}"
            )

        if not isinstance(max_designs, int) or max_designs < 1:
            raise ValueError("max_designs must be a positive integer >= 1.")
        max_designs = min(max_designs, self._MAX_DESIGNS_HARD_LIMIT)

        similarity = similarity.lower().strip()
        if similarity not in self._VALID_METRICS:
            raise ValueError(
                f"similarity must be one of {sorted(self._VALID_METRICS)}, "
                f"got {similarity!r}."
            )

        # ---- Call RetroTide ------------------------------------------ #
        # designPKS signature (from JBEI/RetroTide):
        #   designPKS(
        #       targetSMILES: str,
        #       similarityMetric: str = "tanimoto",   # "tanimoto" | "dice"
        #       numDesigns:      int  = 10,
        #   ) -> list[PKSDesign]
        #
        # Each PKSDesign has attributes:
        #   .similarity    float          – score vs. target
        #   .productSMILES str            – SMILES of PKS product
        #   .modules       list[Module]   – ordered list of PKS modules
        #
        # Each Module has attributes:
        #   .isLoading     bool           – True for loading / starter module
        #   .domains       list[str]      – e.g. ["KS","AT","KR","ACP"]
        #   .extender      str | None     – e.g. "malonyl-CoA"
        #   .atSpecificity str | None     – AT domain substrate label
        raw_designs = self._designPKS(
            smiles,
            similarityMetric=similarity,
            numDesigns=max_designs,
        )

        # ---- Serialise ----------------------------------------------- #
        results: list[dict] = []
        for rank, design in enumerate(raw_designs[:max_designs], start=1):
            modules_out = []
            for mod in design.modules:
                modules_out.append(
                    {
                        "starter":        getattr(mod, "isLoading", False),
                        "domains":        list(getattr(mod, "domains", [])),
                        "extender":       getattr(mod, "extender", None),
                        "at_specificity": getattr(mod, "atSpecificity", None),
                    }
                )

            product_smiles: str = getattr(design, "productSMILES", "")
            results.append(
                {
                    "rank":           rank,
                    "similarity":     round(float(design.similarity), 6),
                    "product_smiles": product_smiles,
                    "modules":        modules_out,
                    "exact_match":    (product_smiles == smiles),
                }
            )

        return results


# --------------------------------------------------------------------------- #
#  Module-level alias so pytest can import the function directly.              #
# --------------------------------------------------------------------------- #
_instance = Retrotide()
try:
    _instance.initiate()
except Exception as e:
    _instance._import_error = e
retrotide_designer = _instance.run