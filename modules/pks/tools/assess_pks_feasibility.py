from __future__ import annotations


class AssessPksFeasibility:
    """
    Description:
        Pre-screens a target molecule to determine whether it is amenable
        to type-I modular PKS biosynthesis.  Returns a structured feasibility
        score (0-1) with diagnostic checks so the caller can decide whether
        to proceed with pks_design_retrotide or fall back to tridentsynth.

    Input:
        smiles (str): SMILES string of the target molecule.

    Output:
        dict: Feasibility assessment with score, checks, molecular
              properties, and a recommendation string.

    Tests:
        - Case:
            Input: "CC(O)CC(=O)CC(O)CC(=O)O"
            Expected Output: feasible=True, score >= 0.6
            Description: Linear oxygenated chain — classic polyketide backbone
        - Case:
            Input: "c1ccccc1"
            Expected Output: feasible=False
            Description: Benzene — aromatic, no oxygens, not PKS-accessible
        - Case:
            Input: ""
            Expected Output: ValueError
            Description: Empty SMILES string
    """

    _PASS = "pass"
    _WARN = "warn"
    _FAIL = "fail"

    _STATUS_MULTIPLIER = {_PASS: 1.0, _WARN: 0.5, _FAIL: 0.0}

    _FEASIBILITY_THRESHOLD = 0.6

    def initiate(self) -> None:
        self._import_error = None
        try:
            from rdkit.Chem import MolFromSmiles, MolToSmiles, MolFromSmarts
            from rdkit.Chem import Descriptors, rdMolDescriptors

            self._MolFromSmiles = MolFromSmiles
            self._MolToSmiles = MolToSmiles
            self._MolFromSmarts = MolFromSmarts
            self._MolWt = Descriptors.MolWt
            self._CalcNumAtomStereoCenters = rdMolDescriptors.CalcNumAtomStereoCenters
            self._CalcNumRings = rdMolDescriptors.CalcNumRings
            self._CalcNumAromaticRings = rdMolDescriptors.CalcNumAromaticRings
        except ImportError as exc:
            self._import_error = exc

    # ------------------------------------------------------------------
    # Individual checks
    # ------------------------------------------------------------------

    def _check_molecular_weight(self, mw: float) -> tuple[str, str]:
        if 100 <= mw <= 2000:
            return self._PASS, f"MW {mw:.1f} Da is within typical polyketide range (100-2000)."
        if 80 <= mw < 100 or 2000 < mw <= 3000:
            return self._WARN, f"MW {mw:.1f} Da is borderline for type-I PKS products."
        return self._FAIL, f"MW {mw:.1f} Da is outside feasible PKS range."

    def _check_carbon_chain(self, num_c: int) -> tuple[str, str]:
        if num_c >= 6:
            return self._PASS, f"{num_c} carbons — sufficient chain length for PKS extension."
        if 4 <= num_c < 6:
            return self._WARN, f"{num_c} carbons — short chain; limited PKS module space."
        return self._FAIL, f"{num_c} carbons — too few for meaningful PKS assembly."

    def _check_heteroatoms(self, num_n: int, num_s: int, num_hal: int) -> tuple[str, str]:
        if num_n == 0 and num_s == 0 and num_hal == 0:
            return self._PASS, "Contains only C, H, O — ideal for PKS."
        if num_hal > 0:
            return self._FAIL, f"Contains {num_hal} halogen(s) — requires specialized halogenase tailoring."
        parts = []
        if num_n > 0:
            parts.append(f"{num_n} nitrogen(s)")
        if num_s > 0:
            parts.append(f"{num_s} sulfur(s)")
        return self._WARN, f"Contains {', '.join(parts)} — post-PKS tailoring enzymes needed."

    def _check_aromatic_rings(self, num_ar: int) -> tuple[str, str]:
        if num_ar == 0:
            return self._PASS, "No aromatic rings — compatible with type-I modular PKS."
        if num_ar <= 2:
            return self._WARN, f"{num_ar} aromatic ring(s) — may be accessible via aromatic starter units."
        return self._FAIL, f"{num_ar} aromatic rings — type-I modular PKS cannot form polyaromatic scaffolds."

    def _check_oxygen_ratio(self, num_o: int, num_c: int) -> tuple[str, str]:
        if num_c == 0:
            return self._FAIL, "No carbons — cannot assess oxygen ratio."
        ratio = num_o / num_c
        if 0.15 <= ratio <= 0.50:
            return self._PASS, f"O/C ratio {ratio:.2f} — typical polyketide oxygenation pattern."
        if 0.10 <= ratio < 0.15 or 0.50 < ratio <= 0.75:
            return self._WARN, f"O/C ratio {ratio:.2f} — atypical but possible oxygenation level."
        if num_o == 0:
            return self._FAIL, "No oxygens — polyketides require oxygenated intermediates."
        return self._FAIL, f"O/C ratio {ratio:.2f} — far outside typical polyketide range."

    def _check_functional_groups(self, mol) -> tuple[str, str]:
        complex_patterns = {
            "epoxide": "[#6]1[#6][O]1",
            "aziridine": "[#6]1[#6][N]1",
            "cyclopropane": "[C]1[C][C]1",
            "nitro": "[N+](=O)[O-]",
            "azide": "[N-]=[N+]=[N-]",
            "phosphate": "[P](=O)([O])([O])[O]",
        }
        found_complex = []
        for name, smarts in complex_patterns.items():
            pat = self._MolFromSmarts(smarts)
            if pat is not None and mol.HasSubstructMatch(pat):
                found_complex.append(name)

        if found_complex:
            return self._FAIL, f"Contains complex groups ({', '.join(found_complex)}) not directly accessible by PKS."

        moderate_patterns = {
            "amide": "[NX3][CX3](=[OX1])",
            "amine": "[NX3;!$(NC=O)]",
            "thioether": "[#6][SX2][#6]",
            "thiol": "[SX2H]",
        }
        found_moderate = []
        for name, smarts in moderate_patterns.items():
            pat = self._MolFromSmarts(smarts)
            if pat is not None and mol.HasSubstructMatch(pat):
                found_moderate.append(name)

        if found_moderate:
            return self._WARN, f"Contains groups ({', '.join(found_moderate)}) that require post-PKS tailoring."

        return self._PASS, "Functional groups are PKS-compatible (ketones, alcohols, ethers, esters)."

    def _check_ring_complexity(self, num_rings: int) -> tuple[str, str]:
        if num_rings <= 1:
            return self._PASS, f"{num_rings} ring(s) — macrolactone or linear chain, typical PKS product."
        if num_rings <= 3:
            return self._WARN, f"{num_rings} rings — moderate ring complexity; may need tailoring cyclases."
        return self._FAIL, f"{num_rings} rings — high ring complexity beyond typical PKS capability."

    # ------------------------------------------------------------------
    # Molecular property extraction
    # ------------------------------------------------------------------

    @staticmethod
    def _count_elements(mol) -> dict:
        counts = {"C": 0, "O": 0, "N": 0, "S": 0, "halogens": 0}
        halogens = {"F", "Cl", "Br", "I"}
        for atom in mol.GetAtoms():
            sym = atom.GetSymbol()
            if sym in counts:
                counts[sym] += 1
            elif sym in halogens:
                counts["halogens"] += 1
        return counts

    # ------------------------------------------------------------------
    # Main entry point
    # ------------------------------------------------------------------

    def run(self, smiles: str) -> dict:
        if self._import_error is not None:
            raise RuntimeError(f"RDKit not available: {self._import_error}")

        if isinstance(smiles, str):
            smiles = smiles.strip()
        if not smiles:
            raise ValueError("smiles must be a non-empty string.")

        mol = self._MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"RDKit could not parse SMILES: {smiles!r}")

        mw = self._MolWt(mol)
        elements = self._count_elements(mol)
        num_stereo = self._CalcNumAtomStereoCenters(mol)
        num_rings = self._CalcNumRings(mol)
        num_aromatic = self._CalcNumAromaticRings(mol)

        mol_props = {
            "molecular_weight": round(mw, 2),
            "num_carbons": elements["C"],
            "num_oxygens": elements["O"],
            "num_nitrogens": elements["N"],
            "num_sulfurs": elements["S"],
            "num_halogens": elements["halogens"],
            "num_stereocenters": num_stereo,
            "num_rings": num_rings,
            "num_aromatic_rings": num_aromatic,
        }

        checks_spec = [
            ("molecular_weight", 0.12, self._check_molecular_weight(mw)),
            ("carbon_chain_length", 0.15, self._check_carbon_chain(elements["C"])),
            ("heteroatom_content", 0.14, self._check_heteroatoms(elements["N"], elements["S"], elements["halogens"])),
            ("aromatic_rings", 0.12, self._check_aromatic_rings(num_aromatic)),
            ("oxygen_ratio", 0.16, self._check_oxygen_ratio(elements["O"], elements["C"])),
            ("functional_groups", 0.21, self._check_functional_groups(mol)),
            ("ring_complexity", 0.10, self._check_ring_complexity(num_rings)),
        ]

        checks = []
        score = 0.0
        for name, weight, (status, detail) in checks_spec:
            contribution = weight * self._STATUS_MULTIPLIER[status]
            score += contribution
            checks.append({
                "name": name,
                "status": status,
                "detail": detail,
                "weight": weight,
            })

        score = round(score, 4)
        feasible = score >= self._FEASIBILITY_THRESHOLD

        if score >= 0.8:
            rec = "Strong PKS target — proceed with pks_design_retrotide."
        elif score >= self._FEASIBILITY_THRESHOLD:
            rec = ("Moderate PKS target — some features may need post-PKS "
                   "tailoring. Consider also running tridentsynth.")
        else:
            rec = ("Poor PKS target — this compound does not match type-I "
                   "modular PKS characteristics. Use tridentsynth for hybrid "
                   "PKS + tailoring pathways.")

        return {
            "feasible": feasible,
            "score": score,
            "molecular_properties": mol_props,
            "checks": checks,
            "recommendation": rec,
        }


_instance = AssessPksFeasibility()
_instance.initiate()
assess_pks_feasibility = _instance.run
