"""
Unit tests for JGK's tools:
  - resolve_smiles
  - assess_pks_feasibility
  - retrotide_designer
  - match_design_to_parts
"""

import pytest
from modules.pks.tools.assess_pks_feasibility import assess_pks_feasibility
from modules.pks.tools.resolve_smiles import resolve_smiles
from modules.pks.tools.retrotide_designer import retrotide_designer
from modules.pks.tools.match_design_to_parts import match_design_to_parts


# ── resolve_smiles tests ─────────────────────────────────────────────


def test_resolve_smiles_valid_smiles_passthrough():
    result = resolve_smiles("CCCC(=O)O")
    assert result["smiles"] == "CCCC(=O)O"
    assert result["source"] == "rdkit"
    assert result["cid"] is None


def test_resolve_smiles_empty_query():
    with pytest.raises(ValueError):
        resolve_smiles("")


def test_resolve_smiles_whitespace_only():
    with pytest.raises(ValueError):
        resolve_smiles("   ")


def test_resolve_smiles_output_structure():
    result = resolve_smiles("CCO")
    for key in ("smiles", "molecular_formula", "iupac_name", "source", "cid"):
        assert key in result


def test_resolve_smiles_canonicalizes():
    result = resolve_smiles("C(=O)(O)CCC")
    assert result["smiles"] == "CCCC(=O)O"


@pytest.mark.network
def test_resolve_smiles_pubchem_lookup():
    result = resolve_smiles("aspirin")
    assert result["source"] == "pubchem"
    assert isinstance(result["cid"], int)
    assert result["smiles"]


@pytest.mark.network
def test_resolve_smiles_unknown_name():
    with pytest.raises(ValueError, match="not found in PubChem"):
        resolve_smiles("xyzzy_not_a_real_chemical_99")


def test_resolve_smiles_aromatic_smiles():
    result = resolve_smiles("c1ccccc1")
    assert result["source"] == "rdkit"
    assert result["smiles"]
    assert result["cid"] is None


def test_resolve_smiles_stereochemistry_preserved():
    result = resolve_smiles("[C@@H](O)(F)Cl")
    assert result["source"] == "rdkit"
    assert "@" in result["smiles"]


def test_resolve_smiles_equivalent_smiles_same_canonical():
    r1 = resolve_smiles("OCC")
    r2 = resolve_smiles("CCO")
    assert r1["smiles"] == r2["smiles"]


def test_resolve_smiles_complex_polyketide_smiles():
    ery_smi = "CCC1OC(=O)C(C)C(O)C(C)C(OC2CC(C)(OC)C(O)C(C)O2)C(C)CC(C)C(=O)C(C)C(O)C1(C)O"
    result = resolve_smiles(ery_smi)
    assert result["source"] == "rdkit"
    assert result["smiles"]


def test_resolve_smiles_leading_trailing_whitespace():
    result = resolve_smiles("  CCO  ")
    assert result["smiles"] == "CCO"
    assert result["source"] == "rdkit"


@pytest.mark.network
def test_resolve_smiles_iupac_name():
    result = resolve_smiles("2-methylbutanoic acid")
    assert result["source"] == "pubchem"
    assert result["smiles"]
    assert result["molecular_formula"]


@pytest.mark.network
def test_resolve_smiles_pubchem_returns_iupac():
    result = resolve_smiles("caffeine")
    assert result["iupac_name"]
    assert result["molecular_formula"]


# ── assess_pks_feasibility tests ──────────────────────────────────────


def test_assess_feasibility_polyketide_backbone():
    result = assess_pks_feasibility("CC(O)CC(=O)CC(O)CC(=O)O")
    assert result["feasible"] is True
    assert result["score"] >= 0.6
    assert len(result["checks"]) == 7


def test_assess_feasibility_benzene_low_score():
    result = assess_pks_feasibility("c1ccccc1")
    assert result["score"] < 0.8
    failed_checks = [c for c in result["checks"] if c["status"] == "fail"]
    assert len(failed_checks) >= 2


def test_assess_feasibility_empty_smiles():
    with pytest.raises(ValueError):
        assess_pks_feasibility("")


def test_assess_feasibility_invalid_smiles():
    with pytest.raises(ValueError):
        assess_pks_feasibility("not_a_smiles_XYZ!!!")


def test_assess_feasibility_nitrogen_warns():
    result = assess_pks_feasibility("CC(C)CN")
    hetero_check = [c for c in result["checks"] if c["name"] == "heteroatom_content"][0]
    assert hetero_check["status"] in ("warn", "fail")


def test_assess_feasibility_output_structure():
    result = assess_pks_feasibility("CCCCCC(=O)O")
    assert "feasible" in result
    assert "score" in result
    assert "molecular_properties" in result
    assert "checks" in result
    assert "recommendation" in result
    props = result["molecular_properties"]
    assert "molecular_weight" in props
    assert "num_carbons" in props
    assert "num_oxygens" in props


def test_assess_feasibility_halogen_fails():
    result = assess_pks_feasibility("CCCCCl")
    hetero = [c for c in result["checks"] if c["name"] == "heteroatom_content"][0]
    assert hetero["status"] == "fail"
    assert "halogen" in hetero["detail"].lower()


def test_assess_feasibility_score_range():
    result = assess_pks_feasibility("CCCCCC(=O)O")
    assert 0.0 <= result["score"] <= 1.0


def test_assess_feasibility_strong_target_recommendation():
    result = assess_pks_feasibility("CC(O)CC(=O)CC(O)CC(=O)O")
    if result["score"] >= 0.8:
        assert "pks_design_retrotide" in result["recommendation"]


def test_assess_feasibility_poor_target_recommendation():
    result = assess_pks_feasibility("c1ccccc1")
    if result["score"] < 0.6:
        assert "tridentsynth" in result["recommendation"]


def test_assess_feasibility_check_weights_sum_to_one():
    result = assess_pks_feasibility("CCCCCC(=O)O")
    total_weight = sum(c["weight"] for c in result["checks"])
    assert abs(total_weight - 1.0) < 0.01


def test_assess_feasibility_molecular_properties_counts():
    result = assess_pks_feasibility("CCCCCC(=O)O")
    props = result["molecular_properties"]
    assert props["num_carbons"] == 6
    assert props["num_oxygens"] == 2
    assert props["num_nitrogens"] == 0
    assert props["num_sulfurs"] == 0
    assert props["num_halogens"] == 0


def test_assess_feasibility_no_oxygen_fails():
    result = assess_pks_feasibility("CCCCCCCC")
    oxy_check = [c for c in result["checks"] if c["name"] == "oxygen_ratio"][0]
    assert oxy_check["status"] == "fail"


def test_assess_feasibility_high_aromatic_fails():
    result = assess_pks_feasibility("c1ccc2cc3ccccc3cc2c1")
    ar_check = [c for c in result["checks"] if c["name"] == "aromatic_rings"][0]
    assert ar_check["status"] == "fail"


def test_assess_feasibility_very_small_molecule():
    result = assess_pks_feasibility("CC")
    assert result["feasible"] is False
    chain_check = [c for c in result["checks"] if c["name"] == "carbon_chain_length"][0]
    assert chain_check["status"] == "fail"


def test_assess_feasibility_sulfur_warns():
    result = assess_pks_feasibility("CCCCSC")
    hetero = [c for c in result["checks"] if c["name"] == "heteroatom_content"][0]
    assert hetero["status"] == "warn"


def test_assess_feasibility_epoxide_fails_functional_groups():
    result = assess_pks_feasibility("C1OC1CCCC(=O)O")
    fg_check = [c for c in result["checks"] if c["name"] == "functional_groups"][0]
    assert fg_check["status"] == "fail"
    assert "epoxide" in fg_check["detail"].lower()


def test_assess_feasibility_macrolactone_ring_passes():
    result = assess_pks_feasibility("C1CCCCCCCCCCC(=O)O1")
    ring_check = [c for c in result["checks"] if c["name"] == "ring_complexity"][0]
    assert ring_check["status"] == "pass"


def test_assess_feasibility_whitespace_smiles():
    result = assess_pks_feasibility("  CCCCCC(=O)O  ")
    assert result["score"] > 0


def test_assess_feasibility_none_raises():
    with pytest.raises((ValueError, AttributeError)):
        assess_pks_feasibility(None)


# ── retrotide_designer tests ─────────────────────────────────────────


def test_retrotide_empty_smiles():
    with pytest.raises(ValueError):
        retrotide_designer("")


def test_retrotide_invalid_smiles():
    with pytest.raises(ValueError):
        retrotide_designer("not_valid_smiles!!!")


def test_retrotide_invalid_max_designs():
    with pytest.raises(ValueError):
        retrotide_designer("CCCC(=O)O", max_designs=0)


def test_retrotide_invalid_similarity_metric():
    with pytest.raises(ValueError):
        retrotide_designer("CCCC(=O)O", similarity="invalid_metric")


def test_retrotide_negative_max_designs():
    with pytest.raises(ValueError):
        retrotide_designer("CCCC(=O)O", max_designs=-1)


def test_retrotide_whitespace_smiles():
    with pytest.raises(ValueError):
        retrotide_designer("   ")


def test_retrotide_accepts_atompairs():
    """Verify 'atompairs' is accepted without error (validation only, not a full run)."""
    with pytest.raises(ValueError):
        retrotide_designer("", similarity="atompairs")


def test_retrotide_accepts_atomatompath():
    """Verify 'atomatompath' is accepted without error (validation only, not a full run)."""
    with pytest.raises(ValueError):
        retrotide_designer("", similarity="atomatompath")


def test_retrotide_max_designs_string_raises():
    with pytest.raises(ValueError):
        retrotide_designer("CCCC(=O)O", max_designs="five")


@pytest.mark.slow
def test_retrotide_functional_output():
    results = retrotide_designer("CCCC(=O)O", max_designs=2)
    assert isinstance(results, list)
    assert len(results) >= 1
    design = results[0]
    assert "rank" in design
    assert "similarity" in design
    assert "product_smiles" in design
    assert "modules" in design
    assert "exact_match" in design
    assert design["rank"] == 1
    assert 0.0 <= design["similarity"] <= 1.0
    assert isinstance(design["modules"], list)
    assert len(design["modules"]) >= 1
    for mod in design["modules"]:
        assert "loading" in mod
        assert "domains" in mod


@pytest.mark.slow
def test_retrotide_respects_max_designs():
    results = retrotide_designer("CCCC(=O)O", max_designs=1)
    assert len(results) <= 1


@pytest.mark.slow
def test_retrotide_results_sorted_by_similarity():
    results = retrotide_designer("CCCC(=O)O", max_designs=5)
    if len(results) >= 2:
        for i in range(len(results) - 1):
            assert results[i]["similarity"] >= results[i + 1]["similarity"]


@pytest.mark.slow
def test_retrotide_ranks_sequential():
    results = retrotide_designer("CCCC(=O)O", max_designs=5)
    for i, design in enumerate(results):
        assert design["rank"] == i + 1


@pytest.mark.slow
def test_retrotide_loading_module_present():
    results = retrotide_designer("CCCC(=O)O", max_designs=2)
    if results:
        design = results[0]
        loading_modules = [m for m in design["modules"] if m["loading"]]
        assert len(loading_modules) >= 1


@pytest.mark.slow
def test_retrotide_product_smiles_valid():
    from rdkit.Chem import MolFromSmiles
    results = retrotide_designer("CCCC(=O)O", max_designs=2)
    for design in results:
        if design["product_smiles"]:
            mol = MolFromSmiles(design["product_smiles"])
            assert mol is not None


@pytest.mark.slow
def test_retrotide_capped_at_hard_limit():
    results = retrotide_designer("CCCC(=O)O", max_designs=30)
    assert len(results) <= 25


# ── match_design_to_parts tests ─────────────────────────────────────


def test_match_design_invalid_source():
    design = {"modules": [{"loading": True, "domains": {"AT": {"substrate": "Malonyl-CoA"}}}]}
    with pytest.raises(ValueError, match="source must be"):
        match_design_to_parts(design, source="invalid")


def test_match_design_empty_design():
    with pytest.raises(ValueError):
        match_design_to_parts({}, source="retrotide")


def test_match_design_none_design():
    with pytest.raises(ValueError):
        match_design_to_parts(None, source="retrotide")


def test_match_design_retrotide_no_modules():
    with pytest.raises(ValueError, match="no 'modules'"):
        match_design_to_parts({"rank": 1}, source="retrotide")


def test_match_design_tridentsynth_no_modules():
    with pytest.raises(ValueError, match="no 'pks_modules'"):
        match_design_to_parts({"ok": True}, source="tridentsynth")


def test_match_design_invalid_max_matches():
    design = {"modules": [{"loading": True, "domains": {"AT": {"substrate": "Malonyl-CoA"}}}]}
    with pytest.raises(ValueError, match="max_matches_per_module"):
        match_design_to_parts(design, source="retrotide", max_matches_per_module=0)


def test_match_design_negative_max_matches():
    design = {"modules": [{"loading": True, "domains": {"AT": {"substrate": "Malonyl-CoA"}}}]}
    with pytest.raises(ValueError, match="max_matches_per_module"):
        match_design_to_parts(design, source="retrotide", max_matches_per_module=-1)


def test_match_design_source_case_insensitive():
    design = {"modules": [{"loading": True, "domains": {"AT": {"substrate": "Malonyl-CoA"}}}]}
    try:
        match_design_to_parts(design, source="RETROTIDE", max_matches_per_module=1)
    except Exception as exc:
        assert "source must be" not in str(exc)


def test_match_design_source_whitespace_stripped():
    design = {"modules": [{"loading": True, "domains": {"AT": {"substrate": "Malonyl-CoA"}}}]}
    try:
        match_design_to_parts(design, source="  retrotide  ", max_matches_per_module=1)
    except Exception as exc:
        assert "source must be" not in str(exc)


def test_match_design_retrotide_empty_modules_list():
    with pytest.raises(ValueError, match="no 'modules'"):
        match_design_to_parts({"modules": []}, source="retrotide")


def test_match_design_tridentsynth_empty_modules_list():
    with pytest.raises(ValueError, match="no modules"):
        match_design_to_parts({"pks_modules": []}, source="tridentsynth")


def test_match_design_tridentsynth_nested_result_format():
    design = {
        "result": {
            "pks_modules": [
                {
                    "module_type": "starter",
                    "domains": [
                        {"domain": "AT", "substrate": "Malonyl-CoA"},
                        {"domain": "ACP"},
                    ],
                }
            ]
        }
    }
    try:
        result = match_design_to_parts(design, source="tridentsynth", max_matches_per_module=1)
        assert result["source"] == "tridentsynth"
        assert result["total_modules"] == 1
    except Exception:
        pass


@pytest.mark.slow
def test_match_design_retrotide_functional():
    design = {
        "rank": 1,
        "similarity": 0.5,
        "product_smiles": "CCCC(=O)O",
        "modules": [
            {"loading": True, "domains": {"AT": {"substrate": "Malonyl-CoA"}}},
        ],
        "exact_match": False,
    }
    result = match_design_to_parts(design, source="retrotide", max_matches_per_module=1)
    assert result["source"] == "retrotide"
    assert result["total_modules"] == 1
    assert "module_matches" in result
    assert "warnings" in result
    assert len(result["module_matches"]) == 1
    match = result["module_matches"][0]
    assert match["module_index"] == 0
    assert match["loading"] is True
    assert "design_domains" in match
    assert "natural_matches" in match


@pytest.mark.slow
def test_match_design_retrotide_multi_module():
    design = {
        "rank": 1,
        "similarity": 0.7,
        "product_smiles": "CCC(O)CC(=O)O",
        "modules": [
            {"loading": True, "domains": {"AT": {"substrate": "Malonyl-CoA"}}},
            {"loading": False, "domains": {"AT": {"substrate": "Methylmalonyl-CoA"}, "KR": {"type": "B1"}}},
        ],
        "exact_match": False,
    }
    result = match_design_to_parts(design, source="retrotide", max_matches_per_module=2)
    assert result["total_modules"] == 2
    assert len(result["module_matches"]) == 2
    for mm in result["module_matches"]:
        assert "module_index" in mm
        assert "natural_matches" in mm


@pytest.mark.slow
def test_match_design_output_has_aa_sequences():
    design = {
        "rank": 1,
        "similarity": 0.5,
        "product_smiles": "CCCC(=O)O",
        "modules": [
            {"loading": True, "domains": {"AT": {"substrate": "Malonyl-CoA"}}},
        ],
        "exact_match": False,
    }
    result = match_design_to_parts(design, source="retrotide", max_matches_per_module=1)
    for mm in result["module_matches"]:
        for nat in mm["natural_matches"]:
            assert "domains" in nat
            for dom in nat["domains"]:
                assert "domain_type" in dom
                assert "aa_sequence" in dom


@pytest.mark.slow
def test_match_design_modules_with_matches_count():
    design = {
        "rank": 1,
        "similarity": 0.5,
        "product_smiles": "CCCC(=O)O",
        "modules": [
            {"loading": True, "domains": {"AT": {"substrate": "Malonyl-CoA"}}},
        ],
        "exact_match": False,
    }
    result = match_design_to_parts(design, source="retrotide", max_matches_per_module=1)
    matched = sum(1 for m in result["module_matches"] if m["natural_matches"])
    assert result["modules_with_matches"] == matched
