"""
Pytest suite for the submit_antismash tool.
Fast tests cover input validation. Network tests verify actual submission
against the live antiSMASH API.
"""

import pytest
from modules.pks.tools.submit_antismash import submit_antismash


# ── Input validation (no network) ────────────────────────────────────

def test_no_inputs_raises():
    with pytest.raises(ValueError, match="Provide one of"):
        submit_antismash()


def test_multiple_inputs_seq_ncbi_raises():
    with pytest.raises(ValueError, match="only one"):
        submit_antismash(seq="ATGAAATAA", ncbi="NC_003888")


def test_multiple_inputs_filepath_ncbi_raises():
    with pytest.raises(ValueError, match="only one"):
        submit_antismash(filepath="modules/pks/data/test_output.gb", ncbi="NC_003888")


def test_seq_too_short_raises():
    with pytest.raises(ValueError, match="1000 bp"):
        submit_antismash(seq="ATGAAATAA")


def test_filepath_missing_raises():
    with pytest.raises(ValueError, match="File not found"):
        submit_antismash(filepath="nonexistent_file.gb")


def test_filepath_wrong_extension_raises():
    with pytest.raises(ValueError, match="GenBank"):
        submit_antismash(filepath="modules/pks/data/test_output.fasta")


# ── Live API submissions (network) ───────────────────────────────────

@pytest.mark.network
def test_filepath_submission():
    msg = submit_antismash(filepath="modules/pks/data/debs_mod1_e2e.gb")
    assert "Job ID is:" in msg
    assert msg.split("Job ID is: ")[1].split(".")[0].startswith("bacteria-")


@pytest.mark.network
def test_ncbi_submission():
    # AM420293 = Saccharopolyspora erythraea NRRL 2338 (contains erythromycin BGC)
    msg = submit_antismash(ncbi="AM420293")
    assert "Job ID is:" in msg
    assert msg.split("Job ID is: ")[1].split(".")[0].startswith("bacteria-")
