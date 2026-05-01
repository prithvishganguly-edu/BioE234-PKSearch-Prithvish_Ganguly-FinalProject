---
name: pks
description: The `pks` module provides tools for polyketide synthase (PKS) design and chemical name resolution, supporting retrobiosynthesis workflows.
---

# pks — Skill Guidance for Gemini

This file is read by the client at startup and injected into Gemini's system prompt.
It gives Gemini domain knowledge to use the PKS tools correctly and present their
results fully.

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
TridentSynth Live Pathway Tool

Description:
This MCP tool connects Gemini to the live TridentSynth web server to generate PKS-based synthesis pathway predictions for a target molecule. The user provides a target SMILES string and selects which synthesis strategies to use, including PKS assembly, biological tailoring, and/or chemical tailoring. The tool submits the query to TridentSynth, waits for the completed results page, parses the best pathway, and returns the key results as text.

Inputs:
- target_smiles: Target molecule as a single SMILES string.
- use_pks: Whether to include PKS assembly.
- use_bio: Whether to include biological tailoring steps.
- use_chem: Whether to include chemical tailoring steps.
- max_bio_steps: Maximum number of biological steps, usually 1–3.
- max_chem_steps: Maximum number of chemical steps, usually 1–3.
- pks_release_mechanism: PKS termination method, such as thiolysis or cyclization.
- pks_starters: Optional PKS starter substrates.
- pks_extenders: Optional PKS extender substrates.
- max_carbon, max_nitrogen, max_oxygen: Optional atom-count filters.
- wait_for_completion: Whether to wait for the result page and parse the completed pathway.

Outputs:
- Job status and task ID.
- Synthesis parameters used by TridentSynth.
- Best pathway found.
- PKS modules used, including module type, domains, and AT substrate choices.
- PKS product SMILES and similarity to the target.
- Post-PKS product SMILES and similarity to the target.
- Reaction SMILES for the best post-PKS transformation.
- Reaction rule name, such as carboxylic acid decarboxylation.
- Reaction enthalpy and feasibility values when available.
- Pathway structures returned as SMILES strings.

Example Use:
User: Use the TridentSynth tool to find a PKS + chemical synthesis pathway for heptane.

Example Output:
The tool submits heptane as CCCCCCC, runs a PKS + chemical search, and returns the best pathway. For example, it may identify a PKS-generated carboxylic acid intermediate that undergoes decarboxylation to form heptane. The output includes the PKS modules, reaction SMILES, pathway structures, reaction rule, enthalpy, and target similarity.

Why This Tool Is Useful:
This tool lets Gemini directly access a specialized PKS retrosynthesis platform instead of only reasoning from general knowledge. It turns TridentSynth from a standalone web tool into an interactive MCP skill that can automatically run pathway searches, extract the most relevant pathway, and summarize the result in plain text for downstream analysis or design decisions.

Limitations:
- The tool depends on the live TridentSynth website being available.
- Runtime depends on how long the TridentSynth server takes to complete the job.
- The parser extracts text from the results page, so major website layout changes could require updates.
- The tool predicts computationally plausible pathways; results should still be biologically and experimentally validated.

### `dna_gc_content`
Computes the GC content (fraction of G and C bases) of a DNA sequence.

---

## Workflow guidance

For a typical PKS design request:
1. If the user gives a name, call `resolve_smiles` to get the SMILES
2. Call `pks_design_retrotide` with the SMILES
3. Present the full results including all module domains
4. Optionally call `tridentsynth` if the user wants alternative pathways
