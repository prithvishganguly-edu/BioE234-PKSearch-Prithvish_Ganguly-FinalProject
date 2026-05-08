"""
Unit tests for the seq_basics example tools and pks tools.

:: Each tool is now a class following the Python Function Object Pattern
(initiate / run).  Tests cover both the canonical class interface AND the
module-level alias  (for example: `reverse_complement = _instance.run`)  so that
direct imports continue to work for students who prefer that style.
"""

import os
import pytest
from modules.pks.tools.reverse_translate import reverse_translate
from modules.pks.tools.submit_antismash import submit_antismash
from modules.pks.tools.check_antismash import check_antismash
from modules.pks.tools.assess_pks_feasibility import assess_pks_feasibility
from modules.pks.tools.resolve_smiles import resolve_smiles
from modules.pks.tools.retrotide_designer import retrotide_designer
from modules.pks.tools.match_design_to_parts import match_design_to_parts


# ── reverse_translate tests ──────────────────────────────────────────


def test_reverse_translate():
    result = reverse_translate("MKV", host="e_coli", filename="test_output.gb")
    assert result["status"] == "success"
    assert "dna_sequence" in result
    assert result["dna_sequence"].startswith("ATG")


def test_reverse_translate_includes_length():
    result = reverse_translate("MKV", host="e_coli", filename="test_output.gb")
    assert "dna_length_bp" in result
    assert result["dna_length_bp"] == len(result["dna_sequence"])


def test_reverse_translate_short_sequence_warns():
    result = reverse_translate("MK", host="e_coli", filename="test_output.gb")
    assert "warning" in result
    assert "1000" in result["warning"]


def test_reverse_translate_long_sequence_no_warning():
    # 334 aa × 3 = 1002 bp — just over the threshold
    long_aa = "MKVLSAFG" * 42
    result = reverse_translate(long_aa, host="e_coli", filename="test_output.gb")
    assert "warning" not in result


def test_reverse_translate_genbank_has_cds():
    from Bio import SeqIO
    result = reverse_translate("MKVL", host="e_coli", filename="test_output.gb")
    rec = SeqIO.read(result["file_saved_at"], "genbank")
    cds_features = [f for f in rec.features if f.type == "CDS"]
    assert len(cds_features) == 1
    assert "translation" in cds_features[0].qualifiers


def test_reverse_translate_empty_raises():
    with pytest.raises(ValueError):
        reverse_translate("", host="e_coli")


def test_reverse_translate_invalid_aa_raises():
    with pytest.raises(ValueError):
        reverse_translate("MKV!XZ", host="e_coli")


# ── submit_antismash tests (no network) ──────────────────────────────


def test_submit_antismash_empty_seq_raises():
    with pytest.raises(ValueError, match="Provide either"):
        submit_antismash()


def test_submit_antismash_short_seq_raises():
    with pytest.raises(ValueError, match="1000 bp"):
        submit_antismash(seq="ATGAAATAA")


def test_submit_antismash_filepath_missing_raises():
    with pytest.raises(ValueError, match="File not found"):
        submit_antismash(filepath="nonexistent_file.gb")


def test_submit_antismash_filepath_wrong_extension_raises():
    with pytest.raises(ValueError, match="GenBank"):
        submit_antismash(filepath="modules/pks/data/test_output.fasta")


# ── check_antismash tests (no network) ──────────────────────────────


def test_check_antismash_empty_job_id_raises():
    with pytest.raises(ValueError, match="cannot be empty"):
        check_antismash("")


@pytest.mark.network
def test_check_antismash_unknown_job_id_raises():
    with pytest.raises(ValueError, match="Failed to check status"):
        check_antismash("bacteria-00000000-0000-0000-0000-000000000000")


"""
from modules.seq_basics.tools.translate import translate
from modules.seq_basics.tools.reverse_complement import reverse_complement


def test_reverse_complement_basic():
    assert reverse_complement("ATGC") == "GCAT"


def test_reverse_complement_ambiguity_codes():
    assert reverse_complement("ATRYSWKMN")


def test_translate_basic():
    assert translate("ATGGCT") == "MA"


def test_translate_frame_validation():
    with pytest.raises(ValueError):
        translate("ATGGCT", frame=0)
    with pytest.raises(ValueError):
        translate("ATGGCT", frame=4)


def test_translate_with_coordinates_and_frame():
    # sequence: A ATG GCT AAA
    # start=1 → ATGGCTAAA
    # frame=1 → ATG GCT AAA → M A K
    assert translate("AATGGCTAAA", start=1, end=None, frame=1) == "MAK"
"""


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
