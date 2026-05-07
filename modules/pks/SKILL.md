# pks — Skill Guidance for Gemini

This file is read by the client at startup and injected into Gemini's system prompt.
It gives Gemini domain knowledge to use the PKS tools correctly and present their
results fully.

---

## What this module does

The `pks` module provides tools for polyketide synthase (PKS) design and chemical
name resolution, supporting retrobiosynthesis workflows.

---

## Tools and when to use them

### `pks_design_retrotide`
Proposes chimeric type-I modular PKS designs for a target molecule using the
RetroTide retrobiosynthesis algorithm.

Use when the user asks to:
- "design a PKS for [molecule]"
- "what PKS can make [molecule]"
- "find a biosynthetic pathway using only a PKS"

**Input:** A SMILES string. If the user provides a chemical name instead, call
`resolve_smiles` first to convert it.

**Presenting results — IMPORTANT:**
Always display the **full module architecture** for each design. Each design contains:
- `rank` and `similarity` score
- `product_smiles` — the SMILES of what the PKS produces
- `modules` — the ordered list of PKS modules. **Always show every module** with:
  - `loading` — whether this is the starter/loading module
  - `domains` — the catalytic domains (AT, KR, DH, ER) and their parameters

Present the modules as a table or structured list so the user can see the full
domain architecture of each proposed PKS. Do NOT omit or summarize the modules.

**Domain key:**
- **AT** (Acyltransferase) — selects the extender unit; `substrate` indicates which
  (e.g. Malonyl-CoA, Methylmalonyl-CoA, butmal)
- **KR** (Ketoreductase) — reduces the beta-keto group; `type` indicates
  stereospecificity (A1, A2, B1, B2, C1, C2)
- **DH** (Dehydratase) — removes water to form a double bond; `type` indicates
  stereochemistry (E or Z)
- **ER** (Enoylreductase) — reduces the double bond; `type` indicates
  stereospecificity (D or L)
- A **loading module** (`loading: true`) is always first and uses a starter unit

### `pks_search_sbspks`
Searches a combined database of 4,088 polyketide structures (2,440 SBSPKS biosynthetic
intermediates + 1,648 MIBiG polyketide compounds) for entries structurally similar to a
query SMILES. Returns ranked hits with Tanimoto similarity scores and engineering hints.

Use when the user asks to:
- "find natural products similar to [molecule]"
- "what PKS pathway builds something like [molecule]"
- "I want to make [molecule] with a PKS, what's closest in nature"
- "search for related polyketides"
- "does this molecule appear as a PKS intermediate"

**Workflow — always follow this order:**
1. If the user gives a compound name (not a SMILES), call `resolve_smiles` first.
2. Then call `pks_search_sbspks` with that SMILES.
3. Never pass a compound name directly to `pks_search_sbspks`.

**Parameters:**
- `query_smiles` — SMILES of the target (use resolve_smiles first if given a name)
- `search_type` — "reaction_search" (default) or "pathway_search" (includes full biosynthetic steps)
- `similarity_threshold` — minimum Tanimoto score; default 0.6, lower to 0.3–0.4 for distant scaffolds
- `max_results` — max hits to return, default 5

**Interpreting results:**
- `similarity_score` 1.0 = exact match already in database
- `similarity_score` >0.7 = very similar scaffold, same compound family
- `similarity_score` 0.4–0.7 = related polyketide family, worth engineering
- `similarity_score` <0.4 = distantly related; mention this caveat to the user
- `is_intermediate=true` means the match is a mid-pathway intermediate, not a final product;
  the `engineering_hint` field explains how to engineer early chain release from that pathway
- `bgc_accession` is present only for MIBiG hits and can be used with antiSMASH
- If no hits appear, suggest lowering `similarity_threshold` to 0.3

### `resolve_smiles`
Converts a chemical name (common, trade, or IUPAC) to a canonical SMILES string
via PubChem. If the input is already valid SMILES, returns immediately.

Use when:
- The user provides a chemical name instead of SMILES (e.g. "erythromycin")
- Another tool requires SMILES input but the user gave a name
- You need to verify or canonicalize a SMILES string

Returns: `smiles`, `molecular_formula`, `iupac_name`, `source`, and `cid`.

### `tridentsynth`
Submits a target SMILES to the TridentSynth synthesis planner and returns
biosynthetic pathway results. Use for broader synthesis planning beyond pure PKS.

### `dna_gc_content`
Computes the GC content (fraction of G and C bases) of a DNA sequence.

---

## Workflow guidance

For a typical PKS design request:
1. If the user gives a name, call `resolve_smiles` to get the SMILES
2. Call `pks_design_retrotide` with the SMILES
3. Present the full results including all module domains
4. Optionally call `tridentsynth` if the user wants alternative pathways
