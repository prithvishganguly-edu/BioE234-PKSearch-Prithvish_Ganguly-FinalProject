# PKS Module — reverse_translate, submit_antismash, check_antismash

**Author:** Opshory Choudhury  
**Course:** BioE 234, Spring 2026  

---

## What these tools do

This set of three tools forms the **assembly and validation** stage of the PKS design pipeline. Once a module architecture has been designed (via RetroTide or ClusterCAD) and an amino acid sequence is in hand, these tools convert it to codon-optimized DNA and verify that the resulting construct actually encodes the expected PKS domain architecture.

```
amino acid sequence
        │
        ▼
  reverse_translate   →  codon-optimized GenBank file
        │
        ▼
  submit_antismash    →  antiSMASH job ID
        │
        ▼
  check_antismash     →  domain annotation + structure prediction + MIBiG hits
```

---

## Files

| File | Description |
|------|-------------|
| `tools/reverse_translate.py` | Codon optimization via DnaChisel, outputs GenBank |
| `tools/reverse_translate.json` | C9 JSON wrapper |
| `tools/submit_antismash.py` | POSTs a DNA sequence to the antiSMASH v1.0 API |
| `tools/submit_antismash.json` | C9 JSON wrapper |
| `tools/check_antismash.py` | Polls antiSMASH for results and parses PKS domain data |
| `tools/check_antismash.json` | C9 JSON wrapper |
| `tools/prompts.json` | Test prompts for Gemini tool-calling evaluation |
| `data/` | Output GenBank files (auto-generated) |

---

## Tool 1 — `reverse_translate`

Converts an amino acid sequence into a codon-optimized DNA sequence for a target host organism and saves it as a GenBank file.

Uses **DnaChisel** with `EnforceTranslation` (locks the amino acid sequence) and `CodonOptimize` (maximizes codon usage for the target host).

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `aa_sequence` | string | required | Amino acid sequence (single-letter codes) |
| `host` | string | `e_coli` | Target organism for codon optimization |
| `filename` | string | `synthetic_pks.gb` | Output GenBank filename, saved to `data/` |

### Supported hosts

| Input alias | Organism |
|-------------|----------|
| `e_coli`, `ecoli` | *Escherichia coli* |
| `s_coelicolor`, `streptomyces` | *Streptomyces coelicolor* |
| `s_albus` | *Streptomyces albus* |
| `s_venezuelae` | *Streptomyces venezuelae* |
| `p_putida`, `pseudomonas` | *Pseudomonas putida* |
| `s_cerevisiae`, `yeast` | *Saccharomyces cerevisiae* |

If an unrecognized host is provided, falls back to *E. coli* and notifies the user.

### Output

```json
{
  "status": "success",
  "dna_sequence": "ATGGTCGCC...",
  "file_saved_at": "modules/pks/data/synthetic_pks.gb",
  "message": "Sequence optimized for s_coelicolor and saved to synthetic_pks.gb."
}
```

### Example usage

```
User: Reverse translate this amino acid sequence MKVL for S. coelicolor.

Gemini calls: reverse_translate(aa_sequence="MKVL", host="s_coelicolor", filename="mkvl_scoel.gb")
```

---

## Tool 2 — `submit_antismash`

Submits a DNA sequence to the public **antiSMASH 7** server for biosynthetic gene cluster analysis. Returns a job ID for status polling.

The sequence is wrapped in FASTA format and uploaded to the antiSMASH v1.0 REST API. Two analyses are always enabled:
- **Active Site Finder (ASF)** — annotates catalytic residues within each PKS domain
- **KnownClusterBlast** — compares your cluster against all MIBiG known clusters and reports percent similarity

### Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `seq` | string | Raw DNA sequence string to analyze (minimum 1000 bp) |

### Output

```
"Job successfully submitted! The Job ID is: bacteria-xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx."
```

### Example usage

```
User: Submit this DNA to antiSMASH: ATGAAACGT...

Gemini calls: submit_antismash(seq="ATGAAACGT...")
→ Returns job ID, then immediately calls check_antismash
```

> **Note:** The tool accepts raw DNA strings only. NCBI accession numbers and MIBiG BGC accessions are not currently supported as direct inputs.

---

## Tool 3 — `check_antismash`

Polls the antiSMASH server for job status and, when complete, parses the results JSON to extract PKS domain architecture, substrate predictions, polymer SMILES, and KnownClusterBlast hits.

### Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `job_id` | string | The job ID returned by `submit_antismash` |

### Output (job still running)

```json
{
  "status": "queued",
  "message": "The job is still processing. Please try again in 1-2 minutes."
}
```

### Output (job complete)

```json
{
  "status": "completed",
  "visualization_url": "https://antismash.secondarymetabolites.org/upload/<job_id>/index.html",
  "domain_predictions": {
    "nrpspksdomains_ctg1_1_PKS_AT.1": {
      "AT_substrate": "mmal",
      "AT_confidence": 100.0
    },
    "nrpspksdomains_ctg1_1_PKS_KR.1": {
      "KR_stereochemistry": "B2",
      "KR_activity": "active"
    }
  },
  "predicted_polymer": {
    "polymer": "(Me-ohmal)",
    "smiles": "C(C)C(O)C(=O)O"
  },
  "pks_clusters": []
}
```

### Output fields explained

| Field | Description |
|-------|-------------|
| `visualization_url` | Direct link to the antiSMASH results page |
| `domain_predictions` | Per-domain annotations keyed by antiSMASH domain ID |
| `domain_predictions[*].AT_substrate` | Predicted extender unit (e.g. `mmal` = methylmalonyl-CoA) |
| `domain_predictions[*].AT_confidence` | Confidence score 0–100 |
| `domain_predictions[*].KR_stereochemistry` | KR stereo type: A1, A2, B1, B2, C1, or C2 |
| `domain_predictions[*].KR_activity` | `"active"` or `"inactive"` |
| `predicted_polymer.polymer` | Shorthand name of the predicted chain extension product |
| `predicted_polymer.smiles` | SMILES of that product |
| `pks_clusters` | BGC region hits — only populated for constructs ≥10 kb |

> **Note:** `domain_predictions` and `predicted_polymer` are always returned for any construct where PKS domains are detected, even short single-module fragments that don't trigger BGC region calling. `pks_clusters` will be empty for constructs under ~10 kb.

### Example usage

```
User: Check antiSMASH job bacteria-1abc9db1-...

Gemini calls: check_antismash(job_id="bacteria-1abc9db1-4f2e-4ab3-a5dd-21210cb2f4b8")
```

### Validation workflow (AI behavior)

When results come back, Gemini should:
1. Recall the original module architecture from RetroTide or TridentSynth
2. Compare `domain_predictions` against what was requested:
   - Does `AT_substrate` match the extender unit RetroTide specified?
   - Does `KR_stereochemistry` match the KR type (A/B)?
   - Are any expected domains (DH, ER) absent?
3. Interpret `predicted_polymer` — does the SMILES match the expected chain extension for that module?
4. Flag any mismatches explicitly
5. Provide the `visualization_url` so the user can inspect the annotated map

---

## Running the tests

```bash
pytest tests/test_tools.py -v -k "reverse_translate"
```

---

## Known API quirks fixed during development

| Bug | Fix |
|-----|-----|
| antiSMASH endpoint was on dead `/api/v2.0/` path | Updated to `/api/v1.0/` |
| Submission field was `seqfile`; API expects `seq` | Corrected field name |
| Status check compared against `"completed"`; API returns `"done"` | Added `"done"` to accepted statuses |
| `genefinder` defaulted to `none` → "all records skipped" for FASTA input | Added `genefinder: prodigal` to payload |
| Result JSON named after uploaded file, not job ID | `check_antismash` now tries multiple candidate filenames |
| `reverse_translate` emitted `misc_feature`; antiSMASH ignores it | Changed to `CDS` feature with embedded translation |
| `check_antismash` only parsed BGC regions → empty output for short constructs | Now always parses `nrps_pks` domain predictions and polymer SMILES directly |
| Sequences under 1000 bp submitted silently and failed server-side | Added client-side length check with clear error message |
