from __future__ import annotations

import requests
from rdkit.Chem import MolFromSmiles, MolToSmiles


_PUBCHEM_URL = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    "/compound/name/{name}/property/IsomericSMILES,MolecularFormula,IUPACName/JSON"
)


class ResolveSmiles:
    """
    Description:
        Resolve a chemical name (common name, trade name, or IUPAC name) to a
        canonical SMILES string by querying PubChem.  If the input is already a
        valid SMILES, it is returned immediately without a network call.

        Use this tool when a user provides a chemical name instead of a SMILES
        string, or when RDKit fails to parse a SMILES and you want to look it
        up by name.

    Input:
        query  (str):  Chemical name or SMILES string to resolve.

    Output:
        dict with keys:
            - "smiles"            (str) : Canonical SMILES string
            - "molecular_formula" (str) : Molecular formula (e.g. "C37H67NO13")
            - "iupac_name"        (str) : IUPAC systematic name
            - "source"            (str) : "rdkit" if already valid SMILES,
                                          "pubchem" if resolved via PubChem
            - "cid"               (int | None) : PubChem Compound ID, or None
                                                 if the input was already SMILES

    Tests:
        - Case:
            Input: query="erythromycin"
            Expected Output: dict with "smiles" containing a non-empty SMILES,
                             "source" == "pubchem", "cid" is an integer
            Description: Common antibiotic name resolved via PubChem.

        - Case:
            Input: query="CCCC(=O)O"
            Expected Output: dict with "smiles" == "CCCC(=O)O",
                             "source" == "rdkit", "cid" is None
            Description: Already-valid SMILES short-circuits without network call.

        - Case:
            Input: query=""
            Expected Output: ValueError raised
            Description: Empty query should raise ValueError.

        - Case:
            Input: query="xyzzy_not_a_chemical"
            Expected Output: ValueError raised
            Description: Unknown name should raise ValueError.
    """

    def initiate(self) -> None:
        pass

    def run(self, query: str) -> dict:
        query = query.strip()
        if not query:
            raise ValueError("query must be a non-empty string.")

        mol = MolFromSmiles(query)
        if mol is not None:
            return {
                "smiles": MolToSmiles(mol),
                "molecular_formula": "",
                "iupac_name": "",
                "source": "rdkit",
                "cid": None,
            }

        try:
            resp = requests.get(
                _PUBCHEM_URL.format(name=requests.utils.quote(query, safe="")),
                timeout=15,
            )
        except requests.RequestException as exc:
            raise RuntimeError(f"PubChem request failed: {exc}") from exc

        if resp.status_code == 404:
            raise ValueError(
                f"Could not resolve {query!r} — not found in PubChem."
            )
        resp.raise_for_status()

        data = resp.json()
        props = data["PropertyTable"]["Properties"][0]

        return {
            "smiles": props.get("IsomericSMILES") or props.get("SMILES", ""),
            "molecular_formula": props.get("MolecularFormula", ""),
            "iupac_name": props.get("IUPACName", ""),
            "source": "pubchem",
            "cid": props.get("CID"),
        }


_instance = ResolveSmiles()
_instance.initiate()
resolve_smiles = _instance.run
