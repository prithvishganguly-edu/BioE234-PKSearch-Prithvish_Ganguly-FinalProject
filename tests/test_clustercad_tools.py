"""
Tests for ClusterCAD MCP tools.
Run with: pytest tests/test_clustercad_tools.py -v
"""

import pytest
from modules.pks.tools.clustercad_list_clusters import ClusterCADListClusters
from modules.pks.tools.clustercad_cluster_details import ClusterCADClusterDetails
from modules.pks.tools.clustercad_get_subunits import ClusterCADGetSubunits
from modules.pks.tools.clustercad_domain_lookup import ClusterCADDomainLookup
from modules.pks.tools.clustercad_subunit_lookup import ClusterCADSubunitLookup
from modules.pks.tools.clustercad_search_domains import ClusterCADSearchDomains


# ---------------------------------------------------------------------------
# Fixtures — instantiate and initiate each tool once per test session
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def list_clusters():
    c = ClusterCADListClusters()
    c.initiate()
    return c

@pytest.fixture(scope="session")
def cluster_details():
    c = ClusterCADClusterDetails()
    c.initiate()
    return c

@pytest.fixture(scope="session")
def get_subunits():
    c = ClusterCADGetSubunits()
    c.initiate()
    return c

@pytest.fixture(scope="session")
def domain_lookup():
    c = ClusterCADDomainLookup()
    c.initiate()
    return c

@pytest.fixture(scope="session")
def subunit_lookup():
    c = ClusterCADSubunitLookup()
    c.initiate()
    return c

@pytest.fixture(scope="session")
def search_domains():
    c = ClusterCADSearchDomains()
    c.initiate()
    return c


# ---------------------------------------------------------------------------
# Tool 1: clustercad_list_clusters
# ---------------------------------------------------------------------------

class TestListClusters:

    def test_returns_list(self, list_clusters):
        results = list_clusters.run()
        assert isinstance(results, list)

    def test_default_returns_results(self, list_clusters):
        results = list_clusters.run()
        assert len(results) > 0

    def test_each_result_has_required_keys(self, list_clusters):
        results = list_clusters.run()
        for r in results:
            assert "accession" in r
            assert "description" in r

    def test_accession_format(self, list_clusters):
        results = list_clusters.run()
        for r in results:
            assert r["accession"].startswith("BGC")

    def test_max_results_respected(self, list_clusters):
        results = list_clusters.run(max_results=5)
        assert len(results) <= 5

    def test_reviewed_only_true(self, list_clusters):
        results = list_clusters.run(reviewed_only=True, max_results=20)
        assert len(results) > 0

    def test_reviewed_only_false_returns_more(self, list_clusters):
        reviewed = list_clusters.run(reviewed_only=True, max_results=50)
        all_clusters = list_clusters.run(reviewed_only=False, max_results=50)
        assert len(all_clusters) >= len(reviewed)

    def test_invalid_max_results_raises(self, list_clusters):
        with pytest.raises(ValueError):
            list_clusters.run(max_results=0)

    def test_abyssomicin_in_results(self, list_clusters):
        results = list_clusters.run(reviewed_only=False, max_results=100)
        descriptions = [r["description"].lower() for r in results]
        assert any("abyssomicin" in d for d in descriptions)


# ---------------------------------------------------------------------------
# Tool 2: clustercad_cluster_details
# ---------------------------------------------------------------------------

class TestClusterDetails:

    ACCESSION = "BGC0001492.1"  # Abyssomicin

    def test_returns_dict(self, cluster_details):
        result = cluster_details.run(self.ACCESSION)
        assert isinstance(result, dict)

    def test_required_keys(self, cluster_details):
        result = cluster_details.run(self.ACCESSION)
        assert "accession" in result
        assert "description" in result
        assert "subunit_count" in result
        assert "module_count" in result
        assert "url" in result

    def test_correct_description(self, cluster_details):
        result = cluster_details.run(self.ACCESSION)
        assert "abyssomicin" in result["description"].lower()

    def test_correct_subunit_count(self, cluster_details):
        result = cluster_details.run(self.ACCESSION)
        assert result["subunit_count"] == 3

    def test_correct_module_count(self, cluster_details):
        result = cluster_details.run(self.ACCESSION)
        assert result["module_count"] == 7

    def test_url_format(self, cluster_details):
        result = cluster_details.run(self.ACCESSION)
        assert result["url"].startswith("https://clustercad.jbei.org")

    def test_empty_accession_raises(self, cluster_details):
        with pytest.raises(ValueError):
            cluster_details.run("")

    def test_invalid_accession_raises(self, cluster_details):
        with pytest.raises(ValueError):
            cluster_details.run("INVALID123")

    def test_nonexistent_accession_raises(self, cluster_details):
        with pytest.raises(ValueError):
            cluster_details.run("BGC9999999.1")

    def test_accession_without_version(self, cluster_details):
        result = cluster_details.run("BGC0001492")
        assert result["accession"] is not None


# ---------------------------------------------------------------------------
# Tool 3: clustercad_get_subunits
# ---------------------------------------------------------------------------

class TestGetSubunits:

    ACCESSION = "BGC0001492.1"  # Abyssomicin — 3 subunits, 7 modules

    def test_returns_list(self, get_subunits):
        result = get_subunits.run(self.ACCESSION)
        assert isinstance(result, list)

    def test_correct_subunit_count(self, get_subunits):
        result = get_subunits.run(self.ACCESSION)
        assert len(result) == 3

    def test_subunit_has_required_keys(self, get_subunits):
        result = get_subunits.run(self.ACCESSION)
        for subunit in result:
            assert "subunit_name" in subunit
            assert "subunit_id" in subunit
            assert "modules" in subunit

    def test_module_has_required_keys(self, get_subunits):
        result = get_subunits.run(self.ACCESSION)
        for subunit in result:
            for module in subunit["modules"]:
                assert "module" in module
                assert "domains" in module
                assert "product_smiles" in module

    def test_domain_has_required_keys(self, get_subunits):
        result = get_subunits.run(self.ACCESSION)
        for subunit in result:
            for module in subunit["modules"]:
                for domain in module["domains"]:
                    assert "domain_type" in domain
                    assert "domain_id" in domain
                    assert "annotation" in domain

    def test_correct_total_modules(self, get_subunits):
        result = get_subunits.run(self.ACCESSION)
        total = sum(len(s["modules"]) for s in result)
        assert total == 7

    def test_first_subunit_name(self, get_subunits):
        result = get_subunits.run(self.ACCESSION)
        assert result[0]["subunit_name"] == "AbsB1"

    def test_loading_module_is_module_0(self, get_subunits):
        result = get_subunits.run(self.ACCESSION)
        first_module = result[0]["modules"][0]["module"]
        assert first_module == "module 0"

    def test_product_smiles_not_empty(self, get_subunits):
        result = get_subunits.run(self.ACCESSION)
        for subunit in result:
            for module in subunit["modules"]:
                assert module["product_smiles"] != ""

    def test_empty_accession_raises(self, get_subunits):
        with pytest.raises(ValueError):
            get_subunits.run("")

    def test_invalid_accession_raises(self, get_subunits):
        with pytest.raises(ValueError):
            get_subunits.run("INVALID")


# ---------------------------------------------------------------------------
# Tool 4: clustercad_domain_lookup
# ---------------------------------------------------------------------------

class TestDomainLookup:

    DOMAIN_ID = 27717  # AT domain in Abyssomicin AbsB1 module 0

    def test_returns_dict(self, domain_lookup):
        result = domain_lookup.run(self.DOMAIN_ID)
        assert isinstance(result, dict)

    def test_required_keys(self, domain_lookup):
        result = domain_lookup.run(self.DOMAIN_ID)
        assert "domain_id" in result
        assert "name" in result
        assert "start" in result
        assert "stop" in result
        assert "annotations" in result
        assert "AAsequence" in result

    def test_correct_domain_id(self, domain_lookup):
        result = domain_lookup.run(self.DOMAIN_ID)
        assert result["domain_id"] == self.DOMAIN_ID

    def test_name_contains_AT(self, domain_lookup):
        result = domain_lookup.run(self.DOMAIN_ID)
        assert "AT" in result["name"]

    def test_aa_sequence_not_empty(self, domain_lookup):
        result = domain_lookup.run(self.DOMAIN_ID)
        assert len(result["AAsequence"]) > 0

    def test_aa_sequence_is_protein(self, domain_lookup):
        result = domain_lookup.run(self.DOMAIN_ID)
        valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
        assert all(c in valid_aa for c in result["AAsequence"])

    def test_start_stop_are_numeric(self, domain_lookup):
        result = domain_lookup.run(self.DOMAIN_ID)
        assert result["start"].isdigit()
        assert result["stop"].isdigit()

    def test_zero_domain_id_raises(self, domain_lookup):
        with pytest.raises(ValueError):
            domain_lookup.run(0)

    def test_negative_domain_id_raises(self, domain_lookup):
        with pytest.raises(ValueError):
            domain_lookup.run(-1)

    def test_string_domain_id_converted(self, domain_lookup):
        result = domain_lookup.run("27717")
        assert result["domain_id"] == self.DOMAIN_ID


# ---------------------------------------------------------------------------
# Tool 5: clustercad_subunit_lookup
# ---------------------------------------------------------------------------

class TestSubunitLookup:

    SUBUNIT_ID = 24119  # AbsB1 in Abyssomicin

    def test_returns_dict(self, subunit_lookup):
        result = subunit_lookup.run(self.SUBUNIT_ID)
        assert isinstance(result, dict)

    def test_required_keys(self, subunit_lookup):
        result = subunit_lookup.run(self.SUBUNIT_ID)
        assert "subunit_id" in result
        assert "name" in result
        assert "start" in result
        assert "stop" in result
        assert "genbank_accession" in result
        assert "AAsequence" in result
        assert "DNAsequence" in result

    def test_correct_name(self, subunit_lookup):
        result = subunit_lookup.run(self.SUBUNIT_ID)
        assert "AbsB1" in result["name"]

    def test_aa_sequence_not_empty(self, subunit_lookup):
        result = subunit_lookup.run(self.SUBUNIT_ID)
        assert len(result["AAsequence"]) > 0

    def test_dna_sequence_not_empty(self, subunit_lookup):
        result = subunit_lookup.run(self.SUBUNIT_ID)
        assert len(result["DNAsequence"]) > 0

    def test_dna_sequence_is_nucleotide(self, subunit_lookup):
        result = subunit_lookup.run(self.SUBUNIT_ID)
        valid_bases = set("ATGCatgc")
        assert all(c in valid_bases for c in result["DNAsequence"])

    def test_dna_length_is_3x_aa(self, subunit_lookup):
        result = subunit_lookup.run(self.SUBUNIT_ID)
        aa_len  = len(result["AAsequence"])
        dna_len = len(result["DNAsequence"])
        assert abs(dna_len - aa_len * 3) <= 3  # allow for stop codon

    def test_genbank_accession_format(self, subunit_lookup):
        result = subunit_lookup.run(self.SUBUNIT_ID)
        assert "." in result["genbank_accession"]  # e.g. ARE67853.1

    def test_zero_subunit_id_raises(self, subunit_lookup):
        with pytest.raises(ValueError):
            subunit_lookup.run(0)

    def test_negative_subunit_id_raises(self, subunit_lookup):
        with pytest.raises(ValueError):
            subunit_lookup.run(-1)

    def test_string_subunit_id_converted(self, subunit_lookup):
        result = subunit_lookup.run("24119")
        assert result["subunit_id"] == self.SUBUNIT_ID


# ---------------------------------------------------------------------------
# Tool 6: clustercad_search_domains
# ---------------------------------------------------------------------------

class TestSearchDomains:

    def test_returns_list(self, search_domains):
        results = search_domains.run(domain_type="AT")
        assert isinstance(results, list)

    def test_search_by_domain_type(self, search_domains):
        results = search_domains.run(domain_type="AT", max_results=5)
        assert len(results) > 0
        for r in results:
            domain_types = [d["domain_type"] for d in r["matching_domains"]]
            assert "AT" in domain_types

    def test_synonym_translation_butyryl(self, search_domains):
        results = search_domains.run(
            domain_type="AT",
            annotation_contains="butyryl",
            max_results=3
        )
        assert len(results) > 0
        for r in results:
            annotations = [d["annotation"] for d in r["matching_domains"]]
            assert any("butmal" in a for a in annotations)

    def test_synonym_translation_malonyl(self, search_domains):
        results = search_domains.run(
            domain_type="AT",
            annotation_contains="malonyl",
            max_results=3
        )
        assert len(results) > 0

    def test_loading_module_only(self, search_domains):
        results = search_domains.run(
            domain_type="AT",
            loading_module_only=True,
            max_results=5
        )
        assert len(results) > 0
        for r in results:
            assert r["module"] == "module 0"

    def test_domain_types_all_present(self, search_domains):
        results = search_domains.run(
            domain_types=["KR", "DH", "ER"],
            max_results=5
        )
        assert len(results) > 0
        for r in results:
            present = {d["domain_type"] for d in r["all_domains"]}
            assert "KR" in present
            assert "DH" in present
            assert "ER" in present

    def test_active_only_no_inactive(self, search_domains):
        results = search_domains.run(
            domain_type="KR",
            active_only=True,
            max_results=5
        )
        assert len(results) > 0
        for r in results:
            for d in r["matching_domains"]:
                assert "inactive" not in d["annotation"].lower()

    def test_exclude_annotation(self, search_domains):
        results = search_domains.run(
            domain_type="DH",
            exclude_annotation="inactive",
            max_results=5
        )
        assert len(results) > 0
        for r in results:
            for d in r["matching_domains"]:
                assert "inactive" not in d["annotation"].lower()

    def test_cluster_description_filter(self, search_domains):
        results = search_domains.run(
            domain_type="AT",
            cluster_description_contains="Abyssomicin",
            max_results=10
        )
        assert len(results) > 0
        for r in results:
            assert "abyssomicin" in r["description"].lower()

    def test_min_modules_filter(self, search_domains):
        results = search_domains.run(
            domain_type="AT",
            min_modules=10,
            max_results=5
        )
        assert len(results) > 0
        for r in results:
            assert r["total_modules"] >= 10

    def test_max_modules_filter(self, search_domains):
        results = search_domains.run(
            domain_type="AT",
            max_modules=5,
            max_results=5
        )
        assert len(results) > 0
        for r in results:
            assert r["total_modules"] <= 5

    def test_reviewed_only(self, search_domains):
        results = search_domains.run(
            domain_type="AT",
            reviewed_only=True,
            max_results=5
        )
        assert len(results) > 0

    def test_max_results_respected(self, search_domains):
        results = search_domains.run(
            domain_type="AT",
            max_results=3
        )
        assert len(results) <= 3

    def test_result_has_required_keys(self, search_domains):
        results = search_domains.run(domain_type="AT", max_results=3)
        for r in results:
            assert "accession" in r
            assert "description" in r
            assert "subunit_name" in r
            assert "subunit_id" in r
            assert "module" in r
            assert "total_modules" in r
            assert "product_smiles" in r
            assert "all_domains" in r
            assert "matching_domains" in r

    def test_no_criteria_raises(self, search_domains):
        with pytest.raises(ValueError):
            search_domains.run(domain_type="", domain_types=[])

    def test_invalid_max_results_raises(self, search_domains):
        with pytest.raises(ValueError):
            search_domains.run(domain_type="AT", max_results=0)

    def test_nonexistent_annotation_raises(self, search_domains):
        with pytest.raises(ValueError):
            search_domains.run(
                domain_type="AT",
                annotation_contains="zzznomatch999"
            )