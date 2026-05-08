"""
Tests for SearchPKS — combined SBSPKS intermediate + MIBiG polyketide search.

Test strategy
-------------
Unit tests (no network, no disk)
    The ``tool`` fixture patches ``_load_or_build_combined_index`` to return
    a small hand-crafted TEST_INDEX (5 entries) covering all cases needed:
    SBSPKS final product, SBSPKS intermediate, SBSPKS starter unit, and a
    MIBiG entry with a BGC accession.  RDKit runs for real; no HTTP goes out.

Integration test (live network)
    Decorated with @pytest.mark.skipif so it is skipped when the SBSPKS
    server or GitHub is unreachable, keeping CI green on flaky networks.
"""

from __future__ import annotations

import importlib.util
import os

import pytest
from unittest.mock import MagicMock, patch

# ---------------------------------------------------------------------------
# SMILES constants — all confirmed from live servers
# ---------------------------------------------------------------------------

PIKROMYCIN_SMILES = (
    "CCC1OC(=O)C(C)C(=O)C(C)C(OC2OC(C)CC(C2O)N(C)C)"
    "C(C)CC(C)C(=O)C=CC1(C)O"
)
ERYTHROMYCIN_SMILES = (
    "CCC1C(C(C(C(=O)C(CC(C(C(C(C(C(=O)O1)C)"
    "OC2CC(C(C(O2)C)O)(C)O)C)OC3C(C(CC(O3)C)N(C)C)O)(C)O)C)C)O)C"
)
# Confirmed from display_smiles.cgi?smile=erythromycin.9.smiles (Step 2 test)
ERYTHROMYCIN_MODULE2_SMILES = "CCC(O)C(C)C(S)=O"
# Confirmed from display_smiles.cgi?smile=erythromycin.11.smiles (Step 2 test)
ERYTHROMYCIN_STARTER_SMILES = "CCC(O)=O"
GLUCOSE_SMILES = "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"

# ---------------------------------------------------------------------------
# Small test index injected into every unit test via the fixture.
# Covers all index entry types needed by the full test suite.
# ---------------------------------------------------------------------------

TEST_INDEX = [
    # 1. SBSPKS final product — pikromycin (self-match anchor)
    {
        "smiles":          PIKROMYCIN_SMILES,
        "source":          "sbspks",
        "is_intermediate": False,
        "module_number":   None,
        "pathway_name":    "Pikromycin",
        "organism":        "Streptomyces venezuelae",
        "compound_name":   "Pikromycin",
        "bgc_accession":   None,
        "path_key":        "pikromycin",
    },
    # 2. SBSPKS final product — erythromycin
    {
        "smiles":          ERYTHROMYCIN_SMILES,
        "source":          "sbspks",
        "is_intermediate": False,
        "module_number":   None,
        "pathway_name":    "Erythromycin",
        "organism":        "Saccharopolyspora erythraea",
        "compound_name":   "Erythromycin",
        "bgc_accession":   None,
        "path_key":        "erythromycin",
    },
    # 3. SBSPKS intermediate — erythromycin module 2 (key test case)
    {
        "smiles":          ERYTHROMYCIN_MODULE2_SMILES,
        "source":          "sbspks",
        "is_intermediate": True,
        "module_number":   2,
        "pathway_name":    "Erythromycin",
        "organism":        "Saccharopolyspora erythraea",
        "compound_name":   "Erythromycin module 2 intermediate",
        "bgc_accession":   None,
        "path_key":        "erythromycin",
    },
    # 4. SBSPKS starter unit — module_number=0
    {
        "smiles":          ERYTHROMYCIN_STARTER_SMILES,
        "source":          "sbspks",
        "is_intermediate": True,
        "module_number":   0,
        "pathway_name":    "Erythromycin",
        "organism":        "Saccharopolyspora erythraea",
        "compound_name":   "Erythromycin (starter unit)",
        "bgc_accession":   None,
        "path_key":        "erythromycin",
    },
    # 5. MIBiG entry — erythromycin A with BGC accession (antiSMASH hook)
    {
        "smiles":          ERYTHROMYCIN_SMILES,
        "source":          "mibig",
        "is_intermediate": False,
        "module_number":   None,
        "pathway_name":    "erythromycin A",
        "organism":        "Saccharopolyspora erythraea NRRL 2338",
        "compound_name":   "erythromycin A",
        "bgc_accession":   "BGC0000055",
    },
]


# ---------------------------------------------------------------------------
# Import helpers
# ---------------------------------------------------------------------------

def _import_module():
    """Load the search_pks module directly, bypassing modules/__init__."""
    spec = importlib.util.spec_from_file_location(
        "search_pks",
        os.path.join(
            os.path.dirname(__file__),
            "..", "modules", "pks", "tools", "search_pks.py",
        ),
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _import_search_pks():
    """Return the SearchPKS class from the module."""
    return _import_module().SearchPKS


# ---------------------------------------------------------------------------
# Mock response helpers
# ---------------------------------------------------------------------------

def _make_mock_pathway_response(path_key: str) -> MagicMock:
    """Return a mock make_reaction.cgi page with one tailoring and one PKS edge."""
    mock_resp = MagicMock()
    mock_resp.raise_for_status.return_value = None
    mock_resp.text = (
        f"{{ data: {{ id: '{path_key}_2-{path_key}_1', source: '{path_key}_2',"
        f" target: '{path_key}_1',"
        f" label: '+methylmalonate\\n{path_key}PKS\\n(KS condensation + KR reduction)',"
        f" href1: '&&moduleno=2&&' }}, classes: 'wrapped' }},"
        f"{{ data: {{ id: '{path_key}_1-{path_key}', source: '{path_key}_1',"
        f" target: '{path_key}', label: 'PikC\\n(Hydroxylation)',"
        f" href1: '/tailoring.php' }}, classes: 'wrapped' }},"
    )
    return mock_resp


# ---------------------------------------------------------------------------
# Fixture
# ---------------------------------------------------------------------------

@pytest.fixture()
def tool():
    """
    Return a fully initialised SearchPKS that searches TEST_INDEX only.

    ``_load_or_build_combined_index`` is patched to return TEST_INDEX so
    no network or disk I/O occurs during setup.  RDKit runs for real so
    fingerprint computation and Tanimoto scores are genuine.  A mock HTTP
    session is wired up so ``pathway_search`` mode can return plausible
    ``pathway_steps`` without hitting the real server.
    """
    SearchPKS = _import_search_pks()
    instance = SearchPKS()

    with patch.object(
        instance, "_load_or_build_combined_index", return_value=TEST_INDEX
    ):
        instance.initiate()

    assert len(instance._combined_index) == len(TEST_INDEX)
    assert len(instance._index_fps) == len(TEST_INDEX)

    mock_session = MagicMock()
    mock_session.get.side_effect = lambda url, **kw: (
        _make_mock_pathway_response(url.split("path=")[-1])
        if "make_reaction.cgi" in url
        else MagicMock(raise_for_status=MagicMock(), text="")
    )
    instance._session = mock_session

    return instance


# ---------------------------------------------------------------------------
# TestCombinedIndex — index structure
# ---------------------------------------------------------------------------

class TestCombinedIndex:
    def test_index_contains_both_sources(self, tool):
        """Combined index must have entries from both sbspks and mibig."""
        sources = {e["source"] for e in tool._combined_index}
        assert "sbspks" in sources
        assert "mibig" in sources

    def test_index_contains_intermediates(self, tool):
        """At least one entry must be marked is_intermediate=True."""
        assert any(e["is_intermediate"] for e in tool._combined_index)

    def test_index_fingerprints_length_matches(self, tool):
        """One fingerprint slot per index entry (some may be None for bad SMILES)."""
        assert len(tool._index_fps) == len(tool._combined_index)

    def test_index_fingerprints_are_valid(self, tool):
        """All TEST_INDEX entries have parseable SMILES — no None fingerprints."""
        assert all(fp is not None for fp in tool._index_fps)


# ---------------------------------------------------------------------------
# TestKnownPolyketide — basic search behaviour on well-known SMILES
# ---------------------------------------------------------------------------

class TestKnownPolyketide:
    def test_pikromycin_returns_self_as_top_hit(self, tool):
        """
        Querying pikromycin's own SMILES must return it as the top hit
        with score exactly 1.0.
        """
        result = tool.run(
            query_smiles=PIKROMYCIN_SMILES,
            search_type="reaction_search",
            max_results=3,
            similarity_threshold=0.3,
        )
        assert result["query_smiles"] == PIKROMYCIN_SMILES
        assert result["search_type_used"] == "reaction_search"
        assert len(result["results"]) >= 1

        top = result["results"][0]
        assert top["compound_name"] == "Pikromycin"
        assert abs(top["similarity_score"] - 1.0) < 1e-6
        assert top["source"] == "sbspks"
        assert top["pathway_name"] == "Pikromycin"
        assert top["organism"] == "Streptomyces venezuelae"

    def test_result_has_new_schema_fields(self, tool):
        """Every result dict must contain all fields defined in the new schema."""
        result = tool.run(
            query_smiles=PIKROMYCIN_SMILES,
            search_type="reaction_search",
            max_results=1,
            similarity_threshold=0.0,
        )
        required = {
            "compound_name", "smiles", "similarity_score", "source",
            "is_intermediate", "module_number", "pathway_name",
            "organism", "bgc_url", "engineering_hint",
            "engineering_recommendation", "pathway_steps",
        }
        for entry in result["results"]:
            assert required.issubset(entry.keys()), (
                f"Missing keys: {required - entry.keys()}"
            )
            assert isinstance(entry["similarity_score"], float)
            assert 0.0 <= entry["similarity_score"] <= 1.0

    def test_total_hits_present_and_valid(self, tool):
        """Output must include total_hits as a non-negative integer."""
        result = tool.run(
            query_smiles=PIKROMYCIN_SMILES,
            max_results=1,
            similarity_threshold=0.0,
        )
        assert "total_hits" in result
        assert isinstance(result["total_hits"], int)
        assert result["total_hits"] >= 1

    def test_total_hits_geq_results_length(self, tool):
        """total_hits counts all threshold-passing entries, not just max_results."""
        result = tool.run(
            query_smiles=ERYTHROMYCIN_SMILES,
            max_results=1,
            similarity_threshold=0.0,
        )
        assert result["total_hits"] >= len(result["results"])

    def test_results_sorted_descending(self, tool):
        """Results must be ordered by similarity_score, highest first."""
        result = tool.run(
            query_smiles=ERYTHROMYCIN_SMILES,
            max_results=10,
            similarity_threshold=0.0,
        )
        scores = [r["similarity_score"] for r in result["results"]]
        assert scores == sorted(scores, reverse=True)

    def test_max_results_respected(self, tool):
        """Never return more entries than max_results."""
        result = tool.run(
            query_smiles=ERYTHROMYCIN_SMILES,
            max_results=2,
            similarity_threshold=0.0,
        )
        assert len(result["results"]) <= 2


# ---------------------------------------------------------------------------
# TestIntermediateHits — the key non-redundant feature vs ClusterCAD
# ---------------------------------------------------------------------------

class TestIntermediateHits:
    def test_returns_at_least_one_intermediate_hit(self, tool):
        """
        Querying the erythromycin module-2 intermediate SMILES against the
        test index must return at least one is_intermediate=True hit.
        """
        result = tool.run(
            query_smiles=ERYTHROMYCIN_MODULE2_SMILES,
            max_results=5,
            similarity_threshold=0.0,
        )
        intermediate_hits = [r for r in result["results"] if r["is_intermediate"]]
        assert len(intermediate_hits) >= 1, (
            "Expected at least one intermediate hit for a query that exactly "
            "matches module 2 of the erythromycin pathway."
        )

    def test_intermediate_hit_has_non_null_module_number(self, tool):
        """
        Every is_intermediate=True hit with module_number > 0 must have a
        non-None, positive integer module_number.
        """
        result = tool.run(
            query_smiles=ERYTHROMYCIN_MODULE2_SMILES,
            max_results=5,
            similarity_threshold=0.0,
        )
        pks_intermediates = [
            r for r in result["results"]
            if r["is_intermediate"] and r["module_number"] != 0
        ]
        assert len(pks_intermediates) >= 1
        for hit in pks_intermediates:
            assert hit["module_number"] is not None
            assert isinstance(hit["module_number"], int)
            assert hit["module_number"] > 0

    def test_intermediate_hit_has_engineering_hint(self, tool):
        """
        Every is_intermediate=True hit must have a non-None engineering_hint
        containing at least the pathway name.
        """
        result = tool.run(
            query_smiles=ERYTHROMYCIN_MODULE2_SMILES,
            max_results=5,
            similarity_threshold=0.0,
        )
        intermediate_hits = [r for r in result["results"] if r["is_intermediate"]]
        assert len(intermediate_hits) >= 1
        for hit in intermediate_hits:
            assert hit["engineering_hint"] is not None, (
                f"engineering_hint is None for intermediate hit: {hit['compound_name']}"
            )
            assert len(hit["engineering_hint"]) > 20

    def test_final_product_has_no_engineering_hint(self, tool):
        """
        is_intermediate=False hits must always have engineering_hint=None.
        """
        result = tool.run(
            query_smiles=PIKROMYCIN_SMILES,
            max_results=1,
            similarity_threshold=0.9,
        )
        final_hits = [r for r in result["results"] if not r["is_intermediate"]]
        assert len(final_hits) >= 1
        for hit in final_hits:
            assert hit["engineering_hint"] is None

    def test_starter_unit_hit_has_module_number_zero(self, tool):
        """
        Starter-unit entries (module_number=0) must appear as intermediates
        with engineering_hint referencing the starter substrate.
        """
        result = tool.run(
            query_smiles=ERYTHROMYCIN_STARTER_SMILES,
            max_results=5,
            similarity_threshold=0.0,
        )
        starter_hits = [
            r for r in result["results"]
            if r["is_intermediate"] and r["module_number"] == 0
        ]
        assert len(starter_hits) >= 1
        for hit in starter_hits:
            assert hit["engineering_hint"] is not None
            assert "starter" in hit["engineering_hint"].lower()


# ---------------------------------------------------------------------------
# TestMibigHits — BGC accession for antiSMASH validator integration
# ---------------------------------------------------------------------------

class TestMibigHits:
    def test_polyketide_returns_mibig_hit(self, tool):
        """
        Querying erythromycin must return at least one source='mibig' result.
        """
        result = tool.run(
            query_smiles=ERYTHROMYCIN_SMILES,
            max_results=5,
            similarity_threshold=0.9,
        )
        mibig_hits = [r for r in result["results"] if r["source"] == "mibig"]
        assert len(mibig_hits) >= 1

    def test_mibig_hit_has_bgc_url(self, tool):
        """
        Every mibig-source result must carry a non-None bgc_url that
        points to the MIBiG reference page.
        """
        result = tool.run(
            query_smiles=ERYTHROMYCIN_SMILES,
            max_results=5,
            similarity_threshold=0.9,
        )
        mibig_hits = [r for r in result["results"] if r["source"] == "mibig"]
        assert len(mibig_hits) >= 1
        for hit in mibig_hits:
            assert hit["bgc_url"] is not None
            assert "mibig.secondarymetabolites.org/go/BGC" in hit["bgc_url"], (
                f"Expected MIBiG URL, got: {hit['bgc_url']}"
            )

    def test_sbspks_hit_has_null_bgc_url(self, tool):
        """
        Every sbspks-source result must have bgc_url=None.
        """
        result = tool.run(
            query_smiles=PIKROMYCIN_SMILES,
            max_results=3,
            similarity_threshold=0.9,
        )
        sbspks_hits = [r for r in result["results"] if r["source"] == "sbspks"]
        assert len(sbspks_hits) >= 1
        for hit in sbspks_hits:
            assert hit["bgc_url"] is None


# ---------------------------------------------------------------------------
# TestNonPolyketide — glucose should not score above threshold
# ---------------------------------------------------------------------------

class TestNonPolyketide:
    def test_glucose_returns_no_hits_above_0_6(self, tool):
        """
        Glucose has no structural similarity to macrolide polyketides.
        No results should appear above a 0.6 Tanimoto threshold.
        """
        result = tool.run(
            query_smiles=GLUCOSE_SMILES,
            max_results=5,
            similarity_threshold=0.6,
        )
        assert len(result["results"]) == 0
        assert len(result["warnings"]) >= 1

    def test_glucose_low_threshold_has_low_scores(self, tool):
        """
        Even at threshold=0.0, all glucose scores must be well below 0.5.
        """
        result = tool.run(
            query_smiles=GLUCOSE_SMILES,
            max_results=5,
            similarity_threshold=0.0,
        )
        for entry in result["results"]:
            assert entry["similarity_score"] < 0.5, (
                f"Unexpectedly high score {entry['similarity_score']} "
                f"for glucose vs {entry['compound_name']}"
            )


# ---------------------------------------------------------------------------
# TestInvalidSmiles — must raise before any computation
# ---------------------------------------------------------------------------

class TestInvalidSmiles:
    @pytest.mark.parametrize("bad_smiles", [
        "not-a-smiles!!",
        "INVALID",
        "C(C)(C)(C)(C)(C)(C)",  # valence violation — RDKit rejects
        "",
    ])
    def test_invalid_smiles_raises_value_error(self, tool, bad_smiles):
        """Any unparseable or empty SMILES must raise ValueError immediately."""
        with pytest.raises(ValueError):
            tool.run(
                query_smiles=bad_smiles,
                max_results=5,
                similarity_threshold=0.5,
            )


# ---------------------------------------------------------------------------
# TestInputValidation — parameter range checks
# ---------------------------------------------------------------------------

class TestInputValidation:
    def test_zero_max_results_raises(self, tool):
        with pytest.raises(ValueError, match="max_results"):
            tool.run(PIKROMYCIN_SMILES, max_results=0)

    def test_threshold_above_one_raises(self, tool):
        with pytest.raises(ValueError, match="similarity_threshold"):
            tool.run(PIKROMYCIN_SMILES, similarity_threshold=1.5)

    def test_threshold_below_zero_raises(self, tool):
        with pytest.raises(ValueError, match="similarity_threshold"):
            tool.run(PIKROMYCIN_SMILES, similarity_threshold=-0.1)

    def test_invalid_search_type_raises(self, tool):
        with pytest.raises(ValueError, match="search_type"):
            tool.run(PIKROMYCIN_SMILES, search_type="magic_search")


# ---------------------------------------------------------------------------
# TestSearchTypes — reaction_search vs pathway_search behaviour
# ---------------------------------------------------------------------------

class TestSearchTypes:
    def test_pathway_search_adds_steps_for_sbspks_hits(self, tool):
        """
        pathway_search must attach pathway_steps to every sbspks-source hit.
        """
        result = tool.run(
            query_smiles=ERYTHROMYCIN_SMILES,
            search_type="pathway_search",
            max_results=5,
            similarity_threshold=0.9,
        )
        assert result["search_type_used"] == "pathway_search"
        sbspks_hits = [r for r in result["results"] if r["source"] == "sbspks"]
        assert len(sbspks_hits) >= 1
        for hit in sbspks_hits:
            assert "pathway_steps" in hit, (
                f"pathway_steps missing for sbspks hit: {hit['compound_name']}"
            )
            assert isinstance(hit["pathway_steps"], list)

    def test_mibig_hits_have_empty_pathway_steps(self, tool):
        """
        MIBiG hits must have pathway_steps=[] — MIBiG has no SBSPKS pathway page.
        """
        result = tool.run(
            query_smiles=ERYTHROMYCIN_SMILES,
            search_type="reaction_search",
            max_results=5,
            similarity_threshold=0.9,
        )
        mibig_hits = [r for r in result["results"] if r["source"] == "mibig"]
        assert len(mibig_hits) >= 1
        for hit in mibig_hits:
            assert "pathway_steps" in hit
            assert hit["pathway_steps"] == []

    def test_reaction_search_always_includes_pathway_steps_field(self, tool):
        """
        pathway_steps must always be present on every result regardless of
        search_type — SBSPKS hits get steps fetched, MIBiG hits get empty list.
        """
        result = tool.run(
            query_smiles=ERYTHROMYCIN_SMILES,
            search_type="reaction_search",
            max_results=5,
            similarity_threshold=0.0,
        )
        for hit in result["results"]:
            assert "pathway_steps" in hit
            assert isinstance(hit["pathway_steps"], list)


# ---------------------------------------------------------------------------
# TestEngineeringHintFunction — unit test for _generate_engineering_hint
# ---------------------------------------------------------------------------

class TestEngineeringRecommendation:
    def test_every_result_has_recommendation(self, tool):
        """engineering_recommendation must be present and non-empty on every hit."""
        result = tool.run(
            query_smiles=ERYTHROMYCIN_SMILES,
            max_results=5,
            similarity_threshold=0.0,
        )
        for hit in result["results"]:
            assert "engineering_recommendation" in hit
            assert isinstance(hit["engineering_recommendation"], str)
            assert len(hit["engineering_recommendation"]) > 20

    def test_exact_match_recommendation_says_no_engineering(self, tool):
        """Exact match (sim=1.0) recommendation must say no engineering needed."""
        result = tool.run(
            query_smiles=PIKROMYCIN_SMILES,
            max_results=1,
            similarity_threshold=0.9,
        )
        top = result["results"][0]
        assert top["similarity_score"] == 1.0
        assert "no engineering" in top["engineering_recommendation"].lower()

    def test_intermediate_recommendation_mentions_te_domain(self, tool):
        """Intermediate hit recommendation must mention TE domain relocation."""
        result = tool.run(
            query_smiles=ERYTHROMYCIN_MODULE2_SMILES,
            max_results=5,
            similarity_threshold=0.0,
        )
        intermediate_hits = [r for r in result["results"] if r["is_intermediate"] and r["module_number"] and r["module_number"] > 0]
        assert len(intermediate_hits) >= 1
        for hit in intermediate_hits:
            assert "te" in hit["engineering_recommendation"].lower() or "thioesterase" in hit["engineering_recommendation"].lower()

    def test_mibig_hits_have_empty_pathway_steps(self, tool):
        """MIBiG hits must always have pathway_steps=[] since they have no SBSPKS page."""
        result = tool.run(
            query_smiles=ERYTHROMYCIN_SMILES,
            max_results=5,
            similarity_threshold=0.9,
        )
        mibig_hits = [r for r in result["results"] if r["source"] == "mibig"]
        assert len(mibig_hits) >= 1
        for hit in mibig_hits:
            assert hit["pathway_steps"] == []


class TestEngineeringHintFunction:
    @pytest.fixture(autouse=True)
    def load_hint_fn(self):
        self._hint = _import_module()._generate_engineering_hint

    def test_regular_intermediate_mentions_module_and_pathway(self):
        entry = {
            "is_intermediate": True, "module_number": 4,
            "pathway_name": "Radicicol",
        }
        hint = self._hint(entry)
        assert "module 4" in hint
        assert "Radicicol" in hint
        assert "TE domain" in hint or "thioesterase" in hint.lower()

    def test_starter_unit_mentions_starter(self):
        entry = {
            "is_intermediate": True, "module_number": 0,
            "pathway_name": "Erythromycin",
        }
        hint = self._hint(entry)
        assert hint is not None
        assert "starter" in hint.lower()
        assert "Erythromycin" in hint

    def test_final_product_returns_none(self):
        entry = {
            "is_intermediate": False, "module_number": None,
            "pathway_name": "Pikromycin",
        }
        assert self._hint(entry) is None

    def test_mibig_final_product_returns_none(self):
        entry = {
            "is_intermediate": False, "module_number": None,
            "pathway_name": "rapamycin", "source": "mibig",
        }
        assert self._hint(entry) is None


# ---------------------------------------------------------------------------
# TestTailoringReactionParsing — unit tests for _fetch_pathway_data
# (internal helper still used by pathway_search mode)
# ---------------------------------------------------------------------------

class TestTailoringReactionParsing:
    def test_pks_elongation_steps_excluded_from_tailoring(self, tool):
        """
        _fetch_pathway_data must not include KS condensation steps in
        tailoring_reactions — only genuine post-PKS modifications.
        """
        html = (
            "{ data: { id: 'ery_2-ery_1', source: 'ery_2', target: 'ery_1',"
            " label: '+methylmalonate\\neryAI\\n(KS condensation + KR reduction)',"
            " href1: '' }, classes: 'wrapped' },"
            "{ data: { id: 'ery_1-ery', source: 'ery_1', target: 'ery',"
            " label: 'eryF\\n(Hydroxylation)', href1: '' }, classes: 'wrapped' },"
        )
        mock_resp = MagicMock()
        mock_resp.raise_for_status.return_value = None
        mock_resp.text = html
        tool._session.get.side_effect = None
        tool._session.get.return_value = mock_resp

        data = tool._fetch_pathway_data("erythromycin", [])
        assert data["tailoring_reactions"] == ["eryF (Hydroxylation)"]
        assert len(data["all_steps"]) == 2

    def test_multiple_tailoring_reactions_all_captured(self, tool):
        """All post-PKS tailoring edges must appear in tailoring_reactions."""
        html = (
            "{ data: { id: 'e_3-e_2', source: 'e_3', target: 'e_2',"
            " label: 'eryF\\n(Hydroxylation)', href1: '' }, classes: 'wrapped' },"
            "{ data: { id: 'e_2-e_1', source: 'e_2', target: 'e_1',"
            " label: 'eryB\\n(O-glycosylation)', href1: '' }, classes: 'wrapped' },"
            "{ data: { id: 'e_1-e', source: 'e_1', target: 'e',"
            " label: 'eryC\\n(O-glycosylation)', href1: '' }, classes: 'wrapped' },"
        )
        mock_resp = MagicMock()
        mock_resp.raise_for_status.return_value = None
        mock_resp.text = html
        tool._session.get.side_effect = None
        tool._session.get.return_value = mock_resp

        data = tool._fetch_pathway_data("erythromycin", [])
        assert len(data["tailoring_reactions"]) == 3


# ---------------------------------------------------------------------------
# Integration test — skipped when live servers are unreachable
# ---------------------------------------------------------------------------

def _servers_reachable() -> bool:
    try:
        import requests
        r = requests.head(
            "http://202.54.226.228/~pksdb/retro/make_reaction.cgi",
            timeout=5,
        )
        return r.status_code < 500
    except Exception:
        return False


def _cache_exists() -> bool:
    """Return True only if the combined index cache is built and has sufficient entries."""
    import json as _json
    cache = os.path.join(
        os.path.dirname(__file__),
        "..", "modules", "pks", "data", "combined_index.json",
    )
    if not os.path.exists(cache) or os.path.getsize(cache) < 100_000:
        return False
    try:
        with open(cache) as f:
            data = _json.load(f)
        # Need at least 1000 entries and at least one sbspks intermediate
        has_intermediates = any(
            e.get("source") == "sbspks" and e.get("is_intermediate")
            for e in data[:500]
        )
        return len(data) >= 1000 and has_intermediates
    except Exception:
        return False


@pytest.mark.skipif(
    not _cache_exists(),
    reason="Combined index cache not built — run initiate() first to build it",
)
class TestLiveIntegration:
    def test_pikromycin_live_search_returns_correct_schema(self):
        """
        Live test: load the full combined index from disk cache (built in
        previous steps) and verify pikromycin returns a correct result with
        all new schema fields present, is_intermediate correct, and at least
        one MIBiG hit with a BGC accession.
        """
        SearchPKS = _import_search_pks()
        instance = SearchPKS()
        instance.initiate()

        result = instance.run(
            query_smiles=PIKROMYCIN_SMILES,
            search_type="reaction_search",
            max_results=10,
            similarity_threshold=0.5,
        )

        compound_names = [r["compound_name"] for r in result["results"]]
        assert "Pikromycin" in compound_names

        pik = next(r for r in result["results"] if r["compound_name"] == "Pikromycin")
        assert abs(pik["similarity_score"] - 1.0) < 1e-4
        assert pik["source"] == "sbspks"
        assert pik["is_intermediate"] is False
        assert pik["module_number"] is None
        assert pik["bgc_url"] is None
        assert pik["engineering_hint"] is None
        assert pik["organism"] == "Streptomyces venezuelae"
        assert "pathway_steps" in pik
        assert "engineering_recommendation" in pik

        mibig_hits = [r for r in result["results"] if r["source"] == "mibig"]
        assert len(mibig_hits) >= 1
        for hit in mibig_hits:
            assert hit["bgc_url"] is not None
            assert "mibig.secondarymetabolites.org/go/BGC" in hit["bgc_url"]

    def test_live_search_returns_intermediate_hits(self):
        """
        Querying a known intermediate SMILES against the full index must
        surface at least one is_intermediate=True hit with engineering_hint.
        """
        SearchPKS = _import_search_pks()
        instance = SearchPKS()
        instance.initiate()

        result = instance.run(
            query_smiles=ERYTHROMYCIN_MODULE2_SMILES,
            search_type="reaction_search",
            max_results=10,
            similarity_threshold=0.5,
        )
        intermediate_hits = [r for r in result["results"] if r["is_intermediate"]]
        assert len(intermediate_hits) >= 1
        for hit in intermediate_hits:
            assert hit["engineering_hint"] is not None
            assert hit["module_number"] is not None
