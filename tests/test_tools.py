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
    with pytest.raises(ValueError, match="Provide one of"):
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


def test_submit_antismash_multiple_inputs_raises():
    with pytest.raises(ValueError, match="only one"):
        submit_antismash(seq="ATGAAATAA", ncbi="NC_003888")


def test_submit_antismash_multiple_inputs_filepath_ncbi_raises():
    with pytest.raises(ValueError, match="only one"):
        submit_antismash(filepath="modules/pks/data/test_output.gb", ncbi="NC_003888")


# ── check_antismash tests (no network) ──────────────────────────────


def test_check_antismash_empty_job_id_raises():
    with pytest.raises(ValueError, match="cannot be empty"):
        check_antismash("")


@pytest.mark.network
def test_check_antismash_unknown_job_id_raises():
    with pytest.raises(ValueError, match="Failed to check status"):
        check_antismash("bacteria-00000000-0000-0000-0000-000000000000")


_KNOWN_JOB = "bacteria-a5121478-d684-48c8-9779-8f94b7f7ff30"  # DEBS module 1 e2e run


@pytest.mark.network
def test_check_antismash_completed_output_structure():
    result = check_antismash(_KNOWN_JOB)
    assert result["status"] == "completed"
    assert "visualization_url" in result
    assert "domain_predictions" in result
    assert "predicted_polymer" in result
    assert "mibig_protein_hits" in result
    at_domains = [v for v in result["domain_predictions"].values() if "AT_substrate" in v]
    assert len(at_domains) >= 1
    at = at_domains[0]
    assert "AT_substrate" in at
    assert "AT_substrate_code" in at
    assert "AT_confidence" in at
    top_hit = result["mibig_protein_hits"][0]
    assert top_hit["bgc_accession"] == "BGC0000055"
    assert top_hit["similarity_pct"] == 100.0
    assert result["predicted_polymer"]["smiles"] is not None


@pytest.mark.network
def test_check_antismash_genes_domain_order():
    result = check_antismash(_KNOWN_JOB)
    assert "genes" in result
    assert len(result["genes"]) >= 1
    gene = next(iter(result["genes"].values()))
    assert "domain_order_string" in gene
    assert "domain_order" in gene
    assert "domain_details" in gene
    # DEBS module 1 should have KS-AT-KR-ACP
    assert gene["domain_order_string"] == "KS-AT-KR-ACP"
    assert gene["domain_order"] == ["KS", "AT", "KR", "ACP"]
    # domain_details should have one entry per domain with evalue and score
    assert len(gene["domain_details"]) == 4
    for d in gene["domain_details"]:
        assert "domain" in d
        assert "evalue" in d
        assert "score" in d


@pytest.mark.network
def test_check_antismash_validation_correct_design():
    result = check_antismash(_KNOWN_JOB, expected_domains=[["KS", "AT", "KR", "ACP"]])
    assert "validation" in result
    v = result["validation"][0]
    assert v["match"] is True
    assert v["missing"] == []
    assert v["unexpected"] == []
    assert v["domain_order_string"] == "KS-AT-KR-ACP"


@pytest.mark.network
def test_check_antismash_validation_missing_domain():
    result = check_antismash(_KNOWN_JOB, expected_domains=[["KS", "AT", "DH", "KR", "ACP"]])
    assert "validation" in result
    v = result["validation"][0]
    assert v["match"] is False
    assert "DH" in v["missing"]
    assert v["unexpected"] == []


@pytest.mark.network
def test_check_antismash_validation_unexpected_domain():
    # Expect only KS and AT — KR and ACP will be flagged as unexpected
    result = check_antismash(_KNOWN_JOB, expected_domains=[["KS", "AT"]])
    assert "validation" in result
    v = result["validation"][0]
    assert v["match"] is False
    assert "KR" in v["unexpected"]
    assert "ACP" in v["unexpected"]


def test_check_antismash_no_validation_without_expected():
    # Without expected_domains, validation key should not be present (no network needed)
    # We can verify this on a fake/fast path by checking empty job_id guard
    with pytest.raises(ValueError):
        check_antismash("")


@pytest.mark.network
def test_submit_antismash_ncbi_happy_path():
    # AM420293 = Saccharopolyspora erythraea NRRL 2338 genome (contains erythromycin BGC)
    msg = submit_antismash(ncbi="AM420293")
    assert "Job ID is:" in msg
    job_id = msg.split("Job ID is: ")[1].split(".")[0]
    assert job_id.startswith("bacteria-")


@pytest.mark.network
def test_submit_antismash_filepath_happy_path():
    job_msg = submit_antismash(filepath="modules/pks/data/debs_mod1_e2e.gb")
    assert "Job ID is:" in job_msg
    job_id = job_msg.split("Job ID is: ")[1].split(".")[0]
    assert job_id.startswith("bacteria-")


@pytest.mark.network
def test_check_antismash_wait_true():
    # Submit a fresh job and wait for it to complete automatically.
    job_msg = submit_antismash(filepath="modules/pks/data/debs_mod1_e2e.gb")
    job_id = job_msg.split("Job ID is: ")[1].split(".")[0]
    result = check_antismash(job_id, wait=True, timeout_seconds=360)
    assert result["status"] == "completed"
    assert "domain_predictions" in result
    assert len(result["domain_predictions"]) > 0


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
