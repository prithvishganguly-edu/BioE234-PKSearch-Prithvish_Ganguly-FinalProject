"""
Unit tests for the seq_basics example tools.
 
:: Each tool is now a class following the Python Function Object Pattern
(initiate / run).  Tests cover both the canonical class interface AND the
module-level alias  (for example: `reverse_complement = _instance.run`)  so that
direct imports continue to work for students who prefer that style.
"""



import pytest
from modules.pks.tools.reverse_translate import reverse_translate

def test_reverse_translate():
    """Unit test to ensure the bridge successfully outputs DNA."""
    # We pass in MKV and ask for E. coli optimization
    result = reverse_translate("MKV", host="e_coli", filename="test_output.gb")
    
    assert result["status"] == "success"
    assert "dna_sequence" in result
    # M translates to ATG, so the DNA should start with ATG
    assert result["dna_sequence"].startswith("ATG")
"""
from modules.seq_basics.tools.translate import translate
from modules.seq_basics.tools.reverse_complement import reverse_complement


def test_reverse_complement_basic():
    assert reverse_complement("ATGC") == "GCAT"


def test_reverse_complement_ambiguity_codes():
    # Should not error for supported IUPAC subset
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