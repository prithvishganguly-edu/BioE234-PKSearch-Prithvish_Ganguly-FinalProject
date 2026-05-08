"""
Pytest suite for the reverse_translate tool.
Tests cover: output structure, dna_length_bp, short-sequence warning,
CDS annotation, host parsing, and error handling.
"""

import pytest
from Bio import SeqIO
from modules.pks.tools.reverse_translate import reverse_translate


def test_reverse_translate_basic():
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
    long_aa = "MKVLSAFG" * 42   # 336 aa = 1008 bp — just over the 1000 bp threshold
    result = reverse_translate(long_aa, host="e_coli", filename="test_output.gb")
    assert "warning" not in result


def test_reverse_translate_genbank_has_cds():
    result = reverse_translate("MKVL", host="e_coli", filename="test_output.gb")
    rec = SeqIO.read(result["file_saved_at"], "genbank")
    cds_features = [f for f in rec.features if f.type == "CDS"]
    assert len(cds_features) == 1
    assert "translation" in cds_features[0].qualifiers


def test_reverse_translate_scoelicolor_host():
    result = reverse_translate("MKVL", host="s_coelicolor", filename="test_output.gb")
    assert result["status"] == "success"
    assert "s_coelicolor" in result["message"]


def test_reverse_translate_salbus_host():
    result = reverse_translate("MKVL", host="s_albus", filename="test_output.gb")
    assert result["status"] == "success"


def test_reverse_translate_file_saved():
    import os
    result = reverse_translate("MKVL", host="e_coli", filename="test_output.gb")
    assert os.path.isfile(result["file_saved_at"])


def test_reverse_translate_empty_raises():
    with pytest.raises(ValueError):
        reverse_translate("", host="e_coli")


def test_reverse_translate_invalid_aa_raises():
    with pytest.raises(ValueError):
        reverse_translate("MKV!XZ", host="e_coli")
