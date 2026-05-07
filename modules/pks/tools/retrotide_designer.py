from __future__ import annotations

"""
MCP-callable wrapper for the RetroTide PKS retrobiosynthesis tool.
"""


class Retrotide:
    """
    Description:
        Given a target chemical structure as a SMILES string, uses the
        RetroTide algorithm to propose chimeric type-I modular PKS designs
        that can synthesise the compound's carbon backbone without any
        post-PKS enzymatic decoration.  Each design is ranked by chemical
        similarity between the predicted PKS product and the target.

        Use this tool when the user asks for a product that can be made
        *purely* by a PKS — i.e. no downstream tailoring enzymes are needed.

    Input:
        smiles       (str):  SMILES string for the target polyketide molecule.
        max_designs  (int):  Maximum number of PKS designs to return per round
                             (default 5, max 25).
        similarity   (str):  Similarity metric for ranking.  Choices:
                             "atompairs" (default), "atomatompath".

    Output:
        list[dict]: Ordered list (best first) of PKS designs.  Each dict has:
            - "rank"            (int)   : 1-based rank (1 = most similar)
            - "similarity"      (float) : similarity score vs. target
            - "product_smiles"  (str)   : SMILES of the PKS product
            - "modules"         (list)  : Ordered list of module descriptors.
              Each module dict contains:
                  "loading"     (bool)         : True only for the loading module
                  "domains"     (dict)         : Domain labels and their params,
                                                 e.g. {"AT": {"substrate": "malonyl-CoA"}}
            - "exact_match"     (bool)  : True when product_smiles equals target

    Tests:
        - Case:
            Input: smiles="CCCC(=O)O", max_designs=3, similarity="atompairs"
            Expected Output: list with dicts each having "rank", "similarity",
                             "product_smiles", "modules", "exact_match"
            Description: Short-chain acid; RetroTide returns at least one design.

        - Case:
            Input: smiles="", max_designs=5, similarity="atompairs"
            Expected Output: ValueError raised
            Description: Empty SMILES should raise ValueError.

        - Case:
            Input: smiles="CCCC(=O)O", max_designs=0, similarity="atompairs"
            Expected Output: ValueError raised
            Description: max_designs must be >= 1.
    """

    _VALID_METRICS = {"atompairs", "atomatompath"}
    _MAX_DESIGNS_HARD_LIMIT = 25

    def initiate(self) -> None:
        self._import_error = None
        self._designPKS = None
        self._MolFromSmiles = None
        self._MolToSmiles = None
        try:
            from retrotide.retrotide import designPKS
            from rdkit.Chem import MolFromSmiles, MolToSmiles
            self._designPKS = designPKS
            self._MolFromSmiles = MolFromSmiles
            self._MolToSmiles = MolToSmiles
            self._patch_bcs_domain()
        except ImportError as exc:
            self._import_error = exc

    @staticmethod
    def _patch_bcs_domain() -> None:
        """Patch bcs to work around two upstream bugs:

        1. Some Domain subclass instances (especially DH) are constructed
           without 'active' being set, crashing __repr__ / __hash__.
        2. designPKS generates module combos missing from structureDB,
           causing KeyError inside Cluster.computeProduct.
        """
        from bcs.bcs import Domain, Cluster
        import retrotide.retrotide as _rt
        if hasattr(Domain, "_patched"):
            return

        def _domain_getattr(self, name):
            if name == "active":
                self.active = True
                return True
            raise AttributeError(
                f"'{type(self).__name__}' object has no attribute {name!r}"
            )

        Domain.__getattr__ = _domain_getattr

        _orig_compute = Cluster.computeProduct

        def _safe_compute(self, structureDB, chain=None):
            try:
                return _orig_compute(self, structureDB, chain=chain)
            except (KeyError, IndexError):
                return None

        Cluster.computeProduct = _safe_compute

        _orig_compare = _rt.compareToTarget

        def _safe_compare(structure, target, similarity="atompairs",
                          targetpathintegers=None):
            if structure is None:
                return 0.0
            return _orig_compare(structure, target, similarity,
                                 targetpathintegers)

        _rt.compareToTarget = _safe_compare

        Domain._patched = True

    def run(
        self,
        smiles: str,
        max_designs: int = 5,
        similarity: str = "atompairs",
    ) -> list[dict]:
        if self._import_error is not None:
            raise RuntimeError(
                f"RetroTide dependencies not available: {self._import_error}"
            )

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

        raw_rounds = self._designPKS(
            mol,
            maxDesignsPerRound=max_designs,
            similarity=similarity,
        )

        all_designs = []
        for rnd in raw_rounds:
            for cluster, score, product_mol in rnd:
                product_smi = self._MolToSmiles(product_mol) if product_mol else ""
                modules_out = []
                for mod in cluster.modules:
                    domains_dict = {}
                    raw_domains = getattr(mod, "domains", {})
                    for domain_cls, domain_obj in raw_domains.items():
                        domain_name = domain_cls.__name__ if hasattr(domain_cls, "__name__") else str(domain_cls)
                        if hasattr(domain_obj, "design"):
                            params = domain_obj.design()
                            params.pop("active", None)
                        elif isinstance(domain_obj, dict):
                            params = dict(domain_obj)
                        else:
                            params = str(domain_obj)
                        domains_dict[domain_name] = params
                    modules_out.append({
                        "loading": getattr(mod, "loading", False),
                        "domains": domains_dict,
                    })

                all_designs.append({
                    "similarity": round(float(score), 6),
                    "product_smiles": product_smi,
                    "modules": modules_out,
                    "exact_match": (product_smi == smiles),
                })

        all_designs.sort(key=lambda d: d["similarity"], reverse=True)
        all_designs = all_designs[:max_designs]

        results = []
        for rank, design in enumerate(all_designs, start=1):
            design["rank"] = rank
            results.append(design)

        return results


_instance = Retrotide()
try:
    _instance.initiate()
except Exception as e:
    _instance._import_error = e
retrotide_designer = _instance.run
