# PKS Module — search_pks

**Author:** Prithvish Ganguly  
**Course:** BioE 234, Spring 2026  

---

## What this tool does

`search_pks` answers the question: *"I want to make compound X using a PKS — what existing biosynthetic pathway in nature is closest to what I need?"*

Given a target molecule as a SMILES string, it searches a combined database of 4,088 polyketide structures and returns the most structurally similar entries, ranked by Tanimoto similarity score. Crucially, it searches not just final natural products but every biosynthetic intermediate along each PKS assembly line. If your target matches a mid-pathway intermediate, the tool tells you exactly which module to truncate to release that structure early.

---

## Files

| File | Description |
|------|-------------|
| `tools/search_pks.py` | Python implementation — `SearchPKS` class with `initiate()` and `run()` |
| `tools/search_pks.json` | C9 JSON wrapper — tool schema, parameters, examples |
| `tools/prompts.json` | 8 test prompts for Gemini tool-calling evaluation |
| `tools/test_search_pks.py` | 43 pytest tests covering all functionality |
| `SKILL.md` | Gemini guidance — when and how to call this tool |
| `THEORY.md` | Algorithm explanation and biological background |
| `data/` | Cache files (auto-generated, not committed) |

---

## How to use it

The tool requires a SMILES string. If you have a compound name, call `resolve_smiles` first:

```
User: I want to make erythromycin with a PKS. What's closest in nature?

Gemini calls: resolve_smiles("erythromycin") → SMILES
Gemini calls: search_pks(query_smiles=SMILES, similarity_threshold=0.4)
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `query_smiles` | string | required | SMILES of the target molecule |
| `search_type` | string | `reaction_search` | `reaction_search` or `pathway_search` |
| `similarity_threshold` | float | `0.6` | Minimum Tanimoto score (0.0–1.0) |
| `max_results` | int | `5` | Maximum hits to return |

### Output fields (per result)

| Field | Description |
|-------|-------------|
| `compound_name` | Name of the matching compound |
| `organism` | Producing organism (all strains joined for MIBiG duplicates) |
| `similarity_score` | Tanimoto coefficient (0.0–1.0) |
| `source` | `sbspks` or `mibig` |
| `is_intermediate` | `true` if match is a mid-pathway intermediate |
| `module_number` | Which PKS module this intermediate is from |
| `engineering_hint` | TE domain relocation hint (intermediates only) |
| `engineering_recommendation` | Plain-English engineering guidance for every hit |
| `bgc_url` | Full MIBiG URL (e.g. `https://mibig.secondarymetabolites.org/go/BGC0000055`) |
| `all_bgc_urls` | All BGC URLs when multiple strains produce the same compound (MIBiG only) |
| `producing_strains` | All producing strains when MIBiG has multiple entries for same compound |
| `pathway_steps` | Biosynthetic assembly line steps — fetched in parallel for all SBSPKS hits |

---

## Database

The combined index has two sources:

**SBSPKS v2** (~2,440 entries)  
Every biosynthetic intermediate and final product from ~225 PKS/NRPS pathways. Intermediates are the key feature — they let you find a pathway that builds your target molecule partway through its assembly line, even if no natural product exactly matches your target.

**MIBiG** (~1,648 entries)  
Curated biosynthetic gene cluster database. Entries include `bgc_accession` fields usable directly with antiSMASH for genome mining.

The index is cached to `modules/pks/data/combined_index.json` and rebuilt automatically with a 30-day TTL.

---

## Running the tests

```bash
pytest tests/test_search_pks.py -v
```

To run only fast tests (no network/database):
```bash
pytest tests/test_search_pks.py -v -m "not slow"
```
