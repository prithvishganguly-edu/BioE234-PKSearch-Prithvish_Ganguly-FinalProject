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
