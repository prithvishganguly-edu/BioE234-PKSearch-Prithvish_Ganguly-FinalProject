"""MCP tool wrapper for TridentSynth query preparation.

This tool does not directly execute a remote synthesis job. Instead, it validates,
normalizes, and auto-fills a TridentSynth-style query so the MCP client can hand
Gemini a clean, structured representation of what should be submitted.
"""

from __future__ import annotations

from typing import Any, Optional


class TridentSynth:
    """
    Description:
        Validate and prepare a TridentSynth synthesis query for a target molecule.
        The tool mirrors the public TridentSynth query form: target SMILES,
        PKS/Bio/Chem strategy selection, optional biological and chemical step
        counts, PKS release mechanism, PKS starter/extender substrates, and
        optional DORAnet C/N/O atom filters.

        If some optional parameters are left blank, the tool can auto-fill a
        reasonable exploratory configuration so the assistant can still proceed.

    Input:
        target_smiles (str): One-molecule SMILES string for the desired target.
        use_pks (bool, optional): Include PKS synthesis search. Default: True.
        use_bio (bool, optional): Include biological tailoring / retrosynthesis.
            Default: False.
        use_chem (bool, optional): Include chemistry tailoring / retrosynthesis.
            Default: False.
        max_bio_steps (int, optional): Maximum biological synthesis steps.
            Allowed: 1-3.
        max_chem_steps (int, optional): Maximum chemical synthesis steps.
            Allowed: 1-3.
        pks_release_mechanism (str, optional): "thiolysis" or "cyclization".
        pks_starters (list[str], optional): PKS starter substrates.
        pks_extenders (list[str], optional): PKS extender substrates.
        max_carbon (int, optional): DORAnet max carbon atom filter.
        max_nitrogen (int, optional): DORAnet max nitrogen atom filter.
        max_oxygen (int, optional): DORAnet max oxygen atom filter.
        auto_optimize_unspecified (bool, optional): If True, fill in omitted
            optional fields with a compact exploratory setup. Default: True.

    Output:
        dict: JSON-serializable summary containing a normalized TridentSynth
        query, the inferred strategy label, any auto-filled defaults, warnings,
        and the public submission URL.

    Tests:
        - Case:
            Input: target_smiles="CCCCCCC", use_pks=True, use_bio=True
            Expected Output: result["ok"] == True
            Description: Valid PKS+Bio query with auto-filled biology steps.
        - Case:
            Input: target_smiles="CC=C(C)C(=O)O", use_pks=True, use_bio=False, use_chem=False
            Expected Output: result["normalized_query"]["strategies"] == ["PKS"]
            Description: PKS-only query is allowed.
        - Case:
            Input: target_smiles="CC.CC"
            Expected Exception: ValueError
            Description: TridentSynth accepts one molecule at a time, so "." is rejected.
        - Case:
            Input: target_smiles="CCCC", use_pks=False, use_bio=False, use_chem=False
            Expected Exception: ValueError
            Description: At least one synthesis strategy must be selected.
        - Case:
            Input: target_smiles="CCCC", max_bio_steps=4
            Expected Exception: ValueError
            Description: Bio/Chem step limits must be in the 1-3 range.
    """

    allowed_release_mechanisms: set[str]
    allowed_starters: dict[str, str]
    allowed_extenders: dict[str, str]

    def initiate(self) -> None:
        self.allowed_release_mechanisms = {"thiolysis", "cyclization"}
        self.allowed_starters = {
            "malonyl-coa": "Malonyl-CoA",
            "methylmalonyl-coa": "Methylmalonyl-CoA",
            "allylmalonyl-coa": "Allylmalonyl-CoA",
            "methoxymalonyl-coa": "Methoxymalonyl-CoA",
            "hydroxymalonyl-coa": "Hydroxymalonyl-CoA",
        }
        self.allowed_extenders = {
            "malonyl-coa": "Malonyl-CoA",
            "methylmalonyl-coa": "Methylmalonyl-CoA",
            "ethylmalonyl-coa": "Ethylmalonyl-CoA",
            "methoxymalonyl-coa": "Methoxymalonyl-CoA",
            "allylmalonyl-coa": "Allylmalonyl-CoA",
        }

    def _canonicalize_choices(
        self,
        values: Optional[list[str]],
        allowed_map: dict[str, str],
        field_name: str,
    ) -> list[str]:
        if not values:
            return []

        cleaned: list[str] = []
        seen: set[str] = set()
        for value in values:
            key = str(value).strip().lower()
            if key not in allowed_map:
                valid = ", ".join(sorted(allowed_map.values()))
                raise ValueError(f"Invalid {field_name}: {value!r}. Allowed values: {valid}")
            canonical = allowed_map[key]
            if canonical not in seen:
                seen.add(canonical)
                cleaned.append(canonical)
        return cleaned

    def _validate_step_count(self, value: Optional[int], field_name: str) -> Optional[int]:
        if value is None:
            return None
        if value not in (1, 2, 3):
            raise ValueError(f"{field_name} must be 1, 2, or 3 — got {value}")
        return value

    def _validate_atom_limit(self, value: Optional[int], field_name: str) -> Optional[int]:
        if value is None:
            return None
        if value < 0:
            raise ValueError(f"{field_name} must be >= 0")
        return value

    def _infer_release_mechanism(self, target_smiles: str) -> str:
        # Lightweight heuristic for omitted release mechanism:
        # if the target looks ring-like and contains O or N, cyclization may be more sensible.
        has_ring_digit = any(ch.isdigit() for ch in target_smiles)
        has_o_or_n = any(ch in target_smiles for ch in ("O", "N", "o", "n"))
        if has_ring_digit and has_o_or_n:
            return "cyclization"
        return "thiolysis"

    def run(
        self,
        target_smiles: str,
        use_pks: bool = True,
        use_bio: bool = False,
        use_chem: bool = False,
        max_bio_steps: Optional[int] = None,
        max_chem_steps: Optional[int] = None,
        pks_release_mechanism: Optional[str] = None,
        pks_starters: Optional[list[str]] = None,
        pks_extenders: Optional[list[str]] = None,
        max_carbon: Optional[int] = None,
        max_nitrogen: Optional[int] = None,
        max_oxygen: Optional[int] = None,
        auto_optimize_unspecified: bool = True,
    ) -> dict[str, Any]:
        target_smiles = target_smiles.strip()
        if not target_smiles:
            raise ValueError("target_smiles must not be empty")
        if "." in target_smiles:
            raise ValueError(
                "TridentSynth accepts one molecule at a time. Remove '.' and submit a single target molecule."
            )
        if not any((use_pks, use_bio, use_chem)):
            raise ValueError("At least one synthesis strategy must be selected: use_pks, use_bio, or use_chem")

        max_bio_steps = self._validate_step_count(max_bio_steps, "max_bio_steps")
        max_chem_steps = self._validate_step_count(max_chem_steps, "max_chem_steps")
        max_carbon = self._validate_atom_limit(max_carbon, "max_carbon")
        max_nitrogen = self._validate_atom_limit(max_nitrogen, "max_nitrogen")
        max_oxygen = self._validate_atom_limit(max_oxygen, "max_oxygen")

        auto_filled: dict[str, Any] = {}
        warnings: list[str] = []

        if auto_optimize_unspecified:
            if use_bio and max_bio_steps is None:
                max_bio_steps = 1
                auto_filled["max_bio_steps"] = 1
            if use_chem and max_chem_steps is None:
                max_chem_steps = 1
                auto_filled["max_chem_steps"] = 1
            if use_pks and pks_release_mechanism is None:
                pks_release_mechanism = self._infer_release_mechanism(target_smiles)
                auto_filled["pks_release_mechanism"] = pks_release_mechanism
            if use_pks and not pks_starters:
                pks_starters = ["Malonyl-CoA", "Methylmalonyl-CoA"]
                auto_filled["pks_starters"] = pks_starters
            if use_pks and not pks_extenders:
                pks_extenders = ["Malonyl-CoA", "Methylmalonyl-CoA"]
                auto_filled["pks_extenders"] = pks_extenders

        if pks_release_mechanism is not None:
            pks_release_mechanism = pks_release_mechanism.strip().lower()
            if pks_release_mechanism not in self.allowed_release_mechanisms:
                allowed = ", ".join(sorted(self.allowed_release_mechanisms))
                raise ValueError(
                    f"Invalid pks_release_mechanism: {pks_release_mechanism!r}. Allowed values: {allowed}"
                )

        starters = self._canonicalize_choices(pks_starters, self.allowed_starters, "pks_starters")
        extenders = self._canonicalize_choices(pks_extenders, self.allowed_extenders, "pks_extenders")

        if not use_pks and (pks_release_mechanism or starters or extenders):
            warnings.append(
                "PKS-specific parameters were provided while use_pks=False; they are retained in the query summary but would be ignored by a non-PKS run."
            )

        strategies: list[str] = []
        if use_pks:
            strategies.append("PKS")
        if use_bio:
            strategies.append("Bio")
        if use_chem:
            strategies.append("Chem")

        strategy_label = " + ".join(strategies)

        normalized_query = {
            "target_smiles": target_smiles,
            "strategies": strategies,
            "max_bio_steps": max_bio_steps,
            "max_chem_steps": max_chem_steps,
            "pks_release_mechanism": pks_release_mechanism,
            "pks_starters": starters,
            "pks_extenders": extenders,
            "doranet_max_atoms": {
                "C": max_carbon,
                "N": max_nitrogen,
                "O": max_oxygen,
            },
        }

        notes: list[str] = [
            "TridentSynth neutralizes charged structures and removes stereochemistry before running.",
            "The public TridentSynth site accepts one target molecule per submission.",
        ]
        if auto_filled:
            notes.append("Some optional inputs were auto-filled to create a compact exploratory query.")

        return {
            "ok": True,
            "tool": "TridentSynth",
            "submission_url": "https://tridentsynth.lbl.gov/run_TridentSynth/",
            "strategy_label": strategy_label,
            "normalized_query": normalized_query,
            "auto_filled": auto_filled,
            "warnings": warnings,
            "notes": notes,
            "next_step": "Submit the normalized query to the public TridentSynth web form or extend this wrapper with a site-specific HTTP client if a documented API becomes available.",
        }


_instance = TridentSynth()
_instance.initiate()
tridentsynth = _instance.run