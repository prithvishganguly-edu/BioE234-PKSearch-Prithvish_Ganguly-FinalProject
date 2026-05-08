"""
Pytest suite for the check_antismash tool.
Fast tests cover input validation. Network tests verify parsing against
a known-completed antiSMASH job (DEBS module 1 e2e run).
"""

import pytest
from modules.pks.tools.check_antismash import check_antismash
from modules.pks.tools.submit_antismash import submit_antismash

# Known-completed job: DEBS module 1, S. coelicolor codon-optimized, GenBank submission.
# Expected architecture: KS-AT-KR-ACP | AT=Methylmalonyl-CoA | KR=B2 active
_KNOWN_JOB = "bacteria-a5121478-d684-48c8-9779-8f94b7f7ff30"


# ── Input validation (no network) ────────────────────────────────────

def test_empty_job_id_raises():
    with pytest.raises(ValueError, match="cannot be empty"):
        check_antismash("")


def test_no_validation_key_without_expected_domains():
    with pytest.raises(ValueError):
        check_antismash("")


# ── Output structure (network) ────────────────────────────────────────

@pytest.mark.network
def test_unknown_job_id_raises():
    with pytest.raises(ValueError, match="Failed to check status"):
        check_antismash("bacteria-00000000-0000-0000-0000-000000000000")


@pytest.mark.network
def test_completed_job_output_structure():
    result = check_antismash(_KNOWN_JOB)
    assert result["status"] == "completed"
    assert "visualization_url" in result
    assert "genes" in result
    assert "domain_predictions" in result
    assert "predicted_polymer" in result
    assert "mibig_protein_hits" in result
    assert "pks_clusters" in result


@pytest.mark.network
def test_domain_predictions_at_keys():
    result = check_antismash(_KNOWN_JOB)
    at_domains = [v for v in result["domain_predictions"].values() if "AT_substrate" in v]
    assert len(at_domains) >= 1
    at = at_domains[0]
    assert at["AT_substrate"] == "Methylmalonyl-CoA"
    assert at["AT_substrate_code"] == "mmal"
    assert at["AT_confidence"] == 100.0


@pytest.mark.network
def test_genes_domain_order_string():
    result = check_antismash(_KNOWN_JOB)
    assert len(result["genes"]) >= 1
    gene = next(iter(result["genes"].values()))
    assert gene["domain_order_string"] == "KS-AT-KR-ACP"
    assert gene["domain_order"] == ["KS", "AT", "KR", "ACP"]


@pytest.mark.network
def test_domain_details_structure():
    result = check_antismash(_KNOWN_JOB)
    gene = next(iter(result["genes"].values()))
    assert len(gene["domain_details"]) == 4
    for d in gene["domain_details"]:
        assert "domain" in d
        assert "evalue" in d
        assert "score" in d


@pytest.mark.network
def test_mibig_protein_hits_top_match():
    result = check_antismash(_KNOWN_JOB)
    top = result["mibig_protein_hits"][0]
    assert top["bgc_accession"] == "BGC0000055"
    assert top["similarity_pct"] == 100.0


@pytest.mark.network
def test_predicted_polymer_present():
    result = check_antismash(_KNOWN_JOB)
    assert result["predicted_polymer"]["polymer"] == "(Me-ohmal)"
    assert result["predicted_polymer"]["smiles"] is not None


# ── Validation / expected_domains (network) ───────────────────────────

@pytest.mark.network
def test_validation_correct_design():
    result = check_antismash(_KNOWN_JOB, expected_domains=[["KS", "AT", "KR", "ACP"]])
    v = result["validation"][0]
    assert v["match"] is True
    assert v["missing"] == []
    assert v["unexpected"] == []


@pytest.mark.network
def test_validation_missing_domain():
    result = check_antismash(_KNOWN_JOB, expected_domains=[["KS", "AT", "DH", "KR", "ACP"]])
    v = result["validation"][0]
    assert v["match"] is False
    assert "DH" in v["missing"]
    assert v["unexpected"] == []


@pytest.mark.network
def test_validation_unexpected_domain():
    result = check_antismash(_KNOWN_JOB, expected_domains=[["KS", "AT"]])
    v = result["validation"][0]
    assert v["match"] is False
    assert "KR" in v["unexpected"]
    assert "ACP" in v["unexpected"]


# ── wait=True polling (network) ───────────────────────────────────────

@pytest.mark.network
def test_wait_true_returns_completed():
    job_msg = submit_antismash(filepath="modules/pks/data/debs_mod1_e2e.gb")
    job_id = job_msg.split("Job ID is: ")[1].split(".")[0]
    result = check_antismash(job_id, wait=True, timeout_seconds=360)
    assert result["status"] == "completed"
    assert len(result["domain_predictions"]) > 0
