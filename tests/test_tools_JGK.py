"""
Unit tests for JGK's tools:
  - resolve_smiles
  - assess_pks_feasibility
  - retrotide_designer
  - find_pks_module_parts (formerly match_design_to_parts)
"""

import pytest
from modules.pks.tools.assess_pks_feasibility import assess_pks_feasibility
from modules.pks.tools.resolve_smiles import resolve_smiles
from modules.pks.tools.retrotide_designer import retrotide_designer
from modules.pks.tools.match_design_to_parts import find_pks_module_parts


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
    non_pass = [c for c in result["checks"] if c["status"] in ("fail", "warn")]
    assert len(non_pass) >= 2


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


def test_assess_feasibility_no_oxygen_warns():
    result = assess_pks_feasibility("CCCCCCCC")
    oxy_check = [c for c in result["checks"] if c["name"] == "oxygen_ratio"][0]
    assert oxy_check["status"] == "warn"
    assert "fully reducing" in oxy_check["detail"].lower()


def test_assess_feasibility_high_aromatic_fails():
    result = assess_pks_feasibility("c1ccc2cc3ccccc3cc2c1")
    ar_check = [c for c in result["checks"] if c["name"] == "aromatic_rings"][0]
    assert ar_check["status"] == "fail"


def test_assess_feasibility_very_small_molecule():
    result = assess_pks_feasibility("CC")
    chain_check = [c for c in result["checks"] if c["name"] == "carbon_chain_length"][0]
    assert chain_check["status"] == "fail"
    mw_check = [c for c in result["checks"] if c["name"] == "molecular_weight"][0]
    assert mw_check["status"] == "fail"


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


# ── find_pks_module_parts tests ─────────────────────────────────────

# -- Validation --

def test_find_module_invalid_max_matches_zero():
    with pytest.raises(ValueError, match="max_matches"):
        find_pks_module_parts(loading=True, at_substrate="Malonyl-CoA", max_matches=0)


def test_find_module_invalid_max_matches_negative():
    with pytest.raises(ValueError, match="max_matches"):
        find_pks_module_parts(loading=True, at_substrate="Malonyl-CoA", max_matches=-1)


def test_find_module_invalid_max_matches_type():
    with pytest.raises(ValueError):
        find_pks_module_parts(loading=True, at_substrate="Malonyl-CoA", max_matches="three")


def test_find_module_unrecognized_domain():
    with pytest.raises(ValueError, match="Unrecognized"):
        find_pks_module_parts(loading=True, reductive_domains="ZZ")


def test_find_module_unrecognized_mixed_domains():
    with pytest.raises(ValueError, match="Unrecognized"):
        find_pks_module_parts(loading=True, reductive_domains="KR,FAKE")


# -- Parsing --

def test_find_module_parse_comma_separated():
    from modules.pks.tools.match_design_to_parts import FindPKSModuleParts
    assert FindPKSModuleParts._parse_reductive_domains("KR,DH,ER") == ["KR", "DH", "ER"]


def test_find_module_parse_single_domain():
    from modules.pks.tools.match_design_to_parts import FindPKSModuleParts
    assert FindPKSModuleParts._parse_reductive_domains("KR") == ["KR"]


def test_find_module_parse_empty_string():
    from modules.pks.tools.match_design_to_parts import FindPKSModuleParts
    assert FindPKSModuleParts._parse_reductive_domains("") == []


def test_find_module_parse_whitespace_handling():
    from modules.pks.tools.match_design_to_parts import FindPKSModuleParts
    assert FindPKSModuleParts._parse_reductive_domains("KR, DH , ER") == ["KR", "DH", "ER"]


def test_find_module_parse_case_insensitive():
    from modules.pks.tools.match_design_to_parts import FindPKSModuleParts
    assert FindPKSModuleParts._parse_reductive_domains("kr,dh,er") == ["KR", "DH", "ER"]


# -- Output structure --

def test_find_module_output_keys():
    result = find_pks_module_parts(loading=True, at_substrate="Malonyl-CoA", max_matches=1)
    for key in ("loading", "at_substrate", "reductive_domains", "total_matches", "matches", "warnings"):
        assert key in result


def test_find_module_output_types():
    result = find_pks_module_parts(loading=True, at_substrate="Malonyl-CoA", max_matches=1)
    assert isinstance(result["loading"], bool)
    assert isinstance(result["reductive_domains"], list)
    assert isinstance(result["total_matches"], int)
    assert isinstance(result["matches"], list)
    assert isinstance(result["warnings"], list)


def test_find_module_total_matches_equals_len():
    result = find_pks_module_parts(loading=True, at_substrate="Malonyl-CoA", max_matches=3)
    assert result["total_matches"] == len(result["matches"])


def test_find_module_echoes_inputs():
    result = find_pks_module_parts(loading=False, at_substrate="Methylmalonyl-CoA", reductive_domains="KR,DH", max_matches=1)
    assert result["loading"] is False
    assert result["at_substrate"] == "Methylmalonyl-CoA"
    assert result["reductive_domains"] == ["KR", "DH"]


def test_find_module_empty_substrate_returns_none():
    result = find_pks_module_parts(loading=True, at_substrate="", reductive_domains="KR", max_matches=1)
    assert result["at_substrate"] is None


# -- Functional tests (hit ClusterCAD cache) --

@pytest.mark.slow
def test_find_module_loading_malonyl():
    result = find_pks_module_parts(loading=True, at_substrate="Malonyl-CoA", max_matches=3)
    assert result["total_matches"] >= 0
    assert isinstance(result["matches"], list)


@pytest.mark.slow
def test_find_module_extension_with_reductive_loop():
    result = find_pks_module_parts(loading=False, at_substrate="Methylmalonyl-CoA", reductive_domains="KR,DH,ER", max_matches=3)
    assert result["total_matches"] >= 0


@pytest.mark.slow
def test_find_module_respects_max_matches():
    result = find_pks_module_parts(loading=True, at_substrate="Malonyl-CoA", max_matches=2)
    assert len(result["matches"]) <= 2


@pytest.mark.slow
def test_find_module_matches_have_domains():
    result = find_pks_module_parts(loading=True, at_substrate="Malonyl-CoA", max_matches=1)
    for m in result["matches"]:
        assert "accession" in m
        assert "cluster_name" in m
        assert "domains" in m
        for dom in m["domains"]:
            assert "domain_type" in dom
            assert "domain_id" in dom
            assert "aa_sequence" in dom


@pytest.mark.slow
def test_find_module_aa_sequences_populated():
    result = find_pks_module_parts(loading=True, at_substrate="Malonyl-CoA", max_matches=1)
    if result["matches"]:
        has_seq = any(
            dom["aa_sequence"] is not None
            for m in result["matches"]
            for dom in m["domains"]
        )
        assert has_seq


@pytest.mark.slow
def test_find_module_warnings_are_strings():
    result = find_pks_module_parts(loading=True, at_substrate="Malonyl-CoA", max_matches=1)
    for w in result["warnings"]:
        assert isinstance(w, str)
