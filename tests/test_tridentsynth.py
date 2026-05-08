from __future__ import annotations

import importlib.util
import os
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
TRIDENTSYNTH_PATH = ROOT / "modules" / "pks" / "tools" / "tridentsynth.py"

spec = importlib.util.spec_from_file_location("tridentsynth_module", TRIDENTSYNTH_PATH)
assert spec is not None
assert spec.loader is not None

tridentsynth_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tridentsynth_module)

TridentSynth = tridentsynth_module.TridentSynth


@pytest.fixture
def tool():
    t = TridentSynth()
    t.initiate()
    return t


def test_build_payload_hexane_pks_chem(tool):
    payload, submitted_query = tool._build_payload(
        target_smiles="CCCCCC",
        use_pks=True,
        use_bio=False,
        use_chem=True,
        max_bio_steps=None,
        max_chem_steps=1,
        pks_release_mechanism=None,
        pks_starters=None,
        pks_extenders=None,
        max_carbon=None,
        max_nitrogen=None,
        max_oxygen=None,
        auto_optimize_unspecified=True,
    )

    payload_preview = tool._payload_preview(payload)

    assert submitted_query["target_smiles"] == "CCCCCC"
    assert submitted_query["use_pks"] is True
    assert submitted_query["use_bio"] is False
    assert submitted_query["use_chem"] is True

    assert payload_preview["smiles"] == "CCCCCC"
    assert payload_preview["synthesisStrategy_pks"] == "on"
    assert "synthesisStrategy_bio" not in payload_preview
    assert payload_preview["synthesisStrategy_chem"] == "on"

    assert payload_preview["rangebio"] == "1"
    assert payload_preview["rangechem"] == "1"

    assert payload_preview["releaseMechanism"] == "thiolysis"
    assert payload_preview["pksStarters[]"] == ["mal", "mmal"]
    assert payload_preview["pksExtenders[]"] == ["mal", "mmal"]

    assert payload_preview["maxAtomsC"] == ""
    assert payload_preview["maxAtomsN"] == ""
    assert payload_preview["maxAtomsO"] == ""


def test_build_payload_respects_explicit_pks_options(tool):
    payload, submitted_query = tool._build_payload(
        target_smiles="CCCCCC",
        use_pks=True,
        use_bio=False,
        use_chem=True,
        max_bio_steps=None,
        max_chem_steps=1,
        pks_release_mechanism="thiolysis",
        pks_starters=["Malonyl-CoA"],
        pks_extenders=["Malonyl-CoA", "Methylmalonyl-CoA"],
        max_carbon=8,
        max_nitrogen=0,
        max_oxygen=2,
        auto_optimize_unspecified=True,
    )

    payload_preview = tool._payload_preview(payload)

    assert submitted_query["pks_release_mechanism"] == "thiolysis"
    assert submitted_query["pks_starters"] == ["mal"]
    assert submitted_query["pks_extenders"] == ["mal", "mmal"]

    assert payload_preview["pksStarters[]"] == "mal"
    assert payload_preview["pksExtenders[]"] == ["mal", "mmal"]

    assert payload_preview["maxAtomsC"] == "8"
    assert payload_preview["maxAtomsN"] == "0"
    assert payload_preview["maxAtomsO"] == "2"


def test_build_payload_bio_only(tool):
    payload, submitted_query = tool._build_payload(
        target_smiles="CC(C(O)=O)=O",
        use_pks=False,
        use_bio=True,
        use_chem=False,
        max_bio_steps=1,
        max_chem_steps=None,
        pks_release_mechanism=None,
        pks_starters=None,
        pks_extenders=None,
        max_carbon=None,
        max_nitrogen=None,
        max_oxygen=None,
        auto_optimize_unspecified=True,
    )

    payload_preview = tool._payload_preview(payload)

    assert submitted_query["target_smiles"] == "CC(C(O)=O)=O"
    assert submitted_query["use_pks"] is False
    assert submitted_query["use_bio"] is True
    assert submitted_query["use_chem"] is False

    assert payload_preview["smiles"] == "CC(C(O)=O)=O"
    assert "synthesisStrategy_pks" not in payload_preview
    assert payload_preview["synthesisStrategy_bio"] == "on"
    assert "synthesisStrategy_chem" not in payload_preview
    assert "releaseMechanism" not in payload_preview
    assert "pksStarters[]" not in payload_preview
    assert "pksExtenders[]" not in payload_preview


def test_build_payload_rejects_multi_molecule_smiles(tool):
    with pytest.raises(ValueError, match="one molecule at a time"):
        tool._build_payload(
            target_smiles="C.CC",
            use_pks=True,
            use_bio=False,
            use_chem=True,
            max_bio_steps=None,
            max_chem_steps=1,
            pks_release_mechanism=None,
            pks_starters=None,
            pks_extenders=None,
            max_carbon=None,
            max_nitrogen=None,
            max_oxygen=None,
            auto_optimize_unspecified=True,
        )


def test_build_payload_requires_at_least_one_strategy(tool):
    with pytest.raises(ValueError, match="At least one synthesis strategy"):
        tool._build_payload(
            target_smiles="CCCCCC",
            use_pks=False,
            use_bio=False,
            use_chem=False,
            max_bio_steps=None,
            max_chem_steps=None,
            pks_release_mechanism=None,
            pks_starters=None,
            pks_extenders=None,
            max_carbon=None,
            max_nitrogen=None,
            max_oxygen=None,
            auto_optimize_unspecified=True,
        )


def test_build_payload_rejects_invalid_step_counts(tool):
    with pytest.raises(ValueError, match="max_chem_steps must be 1, 2, or 3"):
        tool._build_payload(
            target_smiles="CCCCCC",
            use_pks=True,
            use_bio=False,
            use_chem=True,
            max_bio_steps=None,
            max_chem_steps=5,
            pks_release_mechanism=None,
            pks_starters=None,
            pks_extenders=None,
            max_carbon=None,
            max_nitrogen=None,
            max_oxygen=None,
            auto_optimize_unspecified=True,
        )


def test_reaction_to_structures(tool):
    reaction = "CCCCCCC(=O)O>>CCCCCC.O=C=O"

    parsed = tool._reaction_to_structures(reaction)

    assert parsed == {
        "reactants": ["CCCCCCC(=O)O"],
        "products": ["CCCCCC", "O=C=O"],
    }


def test_unique_ordered_filters_bad_smiles(tool):
    values = [
        "i",
        "CCCCCC",
        "CCCCCC",
        "O=C=O",
        None,
        "",
        "info",
    ]

    assert tool._unique_ordered(values) == ["CCCCCC", "O=C=O"]


def test_extract_pks_modules_only_first_candidate_set(tool):
    fake_result_text = """
    MODULE 1 (Loading module)
    KSq
    AT (substrate: Methylmalonyl-CoA)
    ACP

    MODULE 2 (Extension module)
    KS
    AT (substrate: Malonyl-CoA)
    KR
    DH
    ER
    ACP

    MODULE 3 (Extension module)
    KS
    AT (substrate: Malonyl-CoA)
    KR
    DH
    ER
    ACP

    Domain legend

    MODULE 1 (Loading module)
    KSq
    AT (substrate: Malonyl-CoA)
    ACP
    """

    modules = tool._extract_pks_modules(fake_result_text)

    assert len(modules) == 3

    assert modules[0]["module_number"] == 1
    assert modules[0]["module_type"] == "Loading module"
    assert modules[0]["domains"] == [
        {"domain": "KSq", "substrate": None},
        {"domain": "AT", "substrate": "Methylmalonyl-CoA"},
        {"domain": "ACP", "substrate": None},
    ]

    assert modules[1]["module_number"] == 2
    assert "AT substrate Malonyl-CoA" in modules[1]["text_summary"]

    assert modules[2]["module_number"] == 3
    assert "KR" in modules[2]["text_summary"]
    assert "DH" in modules[2]["text_summary"]
    assert "ER" in modules[2]["text_summary"]


def test_parse_result_page_extracts_expected_fields(tool):
    fake_html = """
    <html>
      <body>
        <h2>Synthesis parameters</h2>
        <p>Target SMILES</p>
        <p>CCCCCC</p>
        <p>Target Name</p>
        <p>hexane</p>
        <p>Pathway Sequence</p>
        <p>pks, chem</p>
        <p>PKS Termination Step</p>
        <p>thiolysis</p>
        <p>PKS Extenders</p>
        <p>Malonyl-CoA, Methylmalonyl-CoA</p>
        <p>PKS Starters</p>
        <p>Malonyl-CoA, Methylmalonyl-CoA</p>
        <p># Bio Steps</p>
        <p>1</p>
        <p># Chem Steps</p>
        <p>1</p>
        <p>Job Id</p>
        <p>test-task-id</p>

        <h2>MODULE 1 (Loading module)</h2>
        <p>KSq</p>
        <p>AT (substrate: Methylmalonyl-CoA)</p>
        <p>ACP</p>

        <h2>MODULE 2 (Extension module)</h2>
        <p>KS</p>
        <p>AT (substrate: Malonyl-CoA)</p>
        <p>KR</p>
        <p>DH</p>
        <p>ER</p>
        <p>ACP</p>

        <h2>Full pathway design #1</h2>
        <p>PKS product</p>
        <p>CCCCCCC(=O)O</p>
        <p>similarity to target</p>
        <p>0.67</p>

        <p>Post-PKS product</p>
        <p>CCCCCC</p>
        <p>similarity to target</p>
        <p>1.0</p>

        <p>Reactions (SMILES)</p>
        <p>CCCCCCC(=O)O>>CCCCCC.O=C=O</p>

        <p>Reaction rules</p>
        <p>Carboxylic Acids Decarboxylation</p>

        <p>Reaction enthalpies (kcal/mol)</p>
        <p>-6.19 kcal/mol</p>

        <p>Step feasibilities</p>
        <p>0.99</p>

        <p>Net feasibility</p>
        <p>0.99</p>

        <h2>Full pathway design #2</h2>
      </body>
    </html>
    """

    parsed = tool._parse_result_page(fake_html, task_id="test-task-id")

    assert parsed["task_id"] == "test-task-id"
    assert parsed["synthesis_parameters"]["target_smiles"] == "CCCCCC"
    assert parsed["synthesis_parameters"]["target_name"] == "hexane"
    assert parsed["synthesis_parameters"]["pathway_sequence"] == "pks, chem"
    assert parsed["synthesis_parameters"]["pks_termination_step"] == "thiolysis"

    assert len(parsed["pks_modules"]) == 2
    assert parsed["pks_modules"][0]["module_type"] == "Loading module"

    best = parsed["best_pathway"]

    assert best["pks_product_smiles"] == "CCCCCCC(=O)O"
    assert best["pks_similarity_to_target"] == 0.67
    assert best["post_pks_product_smiles"] == "CCCCCC"
    assert best["post_pks_similarity_to_target"] == 1.0

    assert best["reaction_smiles"] == ["CCCCCCC(=O)O>>CCCCCC.O=C=O"]
    assert best["reaction_structures"] == [
        {
            "reactants": ["CCCCCCC(=O)O"],
            "products": ["CCCCCC", "O=C=O"],
        }
    ]

    assert best["reaction_rule_names"] == ["Carboxylic Acids Decarboxylation"]
    assert best["reaction_enthalpies"] == ["-6.19 kcal/mol"]
    assert best["pathway_structures_smiles"] == [
        "CCCCCCC(=O)O",
        "CCCCCC",
        "O=C=O",
    ]


def test_add_selected_steps_fixes_target_and_selected_steps(tool):
    parsed = {
        "synthesis_parameters": {
            "target_smiles": "CCCCCCC",
            "pathway_sequence": "pks, chem",
            "bio_steps": "1",
            "chem_steps": "1",
            "pks_termination_step": "thiolysis",
            "pks_starters": "Malonyl-CoA, Methylmalonyl-CoA",
            "pks_extenders": "Malonyl-CoA, Methylmalonyl-CoA",
        },
        "pks_modules": [
            {
                "text_summary": (
                    "Module 1 (Loading module): KSq, "
                    "AT substrate Methylmalonyl-CoA, ACP"
                )
            }
        ],
        "best_pathway": {
            "pks_product_smiles": "CCCCCCC(=O)O",
            "pks_similarity_to_target": 0.67,
            "post_pks_product_smiles": "CCCCCCC",
            "post_pks_similarity_to_target": 1.0,
            "reaction_smiles": ["CCCCCCC(=O)O>>CCCCCC.O=C=O"],
            "reaction_structures": [
                {
                    "reactants": ["CCCCCCC(=O)O"],
                    "products": ["CCCCCC", "O=C=O"],
                }
            ],
            "reaction_rule_names": ["Carboxylic Acids Decarboxylation"],
            "reaction_enthalpies": ["-6.19 kcal/mol"],
            "pathway_structures_smiles": ["CCCCCCC", "CCCCCCC(=O)O"],
        },
    }

    submitted_query = {
        "target_smiles": "CCCCCC",
        "use_pks": True,
        "use_bio": False,
        "use_chem": True,
    }

    cleaned = tool._add_selected_steps(parsed, submitted_query)

    assert cleaned["synthesis_parameters"]["target_smiles"] == "CCCCCC"
    assert cleaned["synthesis_parameters"]["pathway_sequence"] == "PKS, Chemical synthesis"

    assert cleaned["selected_steps"] == {
        "bio_steps": None,
        "chem_steps": "1",
    }

    best = cleaned["best_pathway"]

    assert best["post_pks_product_smiles"] == "CCCCCC"
    assert best["pathway_structures_smiles"] == [
        "CCCCCCC(=O)O",
        "CCCCCC",
        "O=C=O",
    ]

    if "text_summary" in cleaned:
        assert "Target SMILES: CCCCCC" in cleaned["text_summary"]
        assert "Bio steps: 1" not in cleaned["text_summary"]


def test_live_submission_only_optional(tool):
    if os.getenv("RUN_LIVE_TRIDENTSYNTH") != "1":
        pytest.skip("Set RUN_LIVE_TRIDENTSYNTH=1 to run live TridentSynth submission test.")

    result = tool.run(
        target_smiles="CCCCCC",
        use_pks=True,
        use_bio=False,
        use_chem=True,
        max_chem_steps=1,
        wait_for_completion=False,
    )

    assert result["ok"] is True
    assert result["status"] == "submitted"
    assert result["task_id"]
    assert result["submitted_query"]["target_smiles"] == "CCCCCC"