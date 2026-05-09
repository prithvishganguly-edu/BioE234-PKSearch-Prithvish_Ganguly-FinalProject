# Theory — TridentSynth MCP Tool

**Author:** Shalen Ardeshna  
**Course:** BioE 134, Spring 2026  

---

## Overview

The `tridentsynth` MCP tool connects Gemini to the live TridentSynth web server so that a user can submit a target molecule and receive a PKS-based synthesis pathway prediction.

The tool does not implement a full retrosynthesis algorithm locally. Instead, it acts as an MCP integration layer around TridentSynth. It handles input validation, live job submission, result polling, explicit result parsing, output cleanup, and deterministic summarization.

The main purpose is to make a specialized PKS pathway design platform usable from natural-language prompts.

---

## Biological background

Polyketide synthases, or PKSs, are modular biosynthetic enzymes that assemble complex carbon skeletons from small acyl-CoA building blocks. A modular type I PKS works like an assembly line: each module adds one extender unit to a growing chain and may chemically modify the chain before passing it to the next module.

Common PKS building blocks include:

| Building block | Role |
|---------------|------|
| Malonyl-CoA | Common extender unit that adds a two-carbon unit |
| Methylmalonyl-CoA | Extender unit that introduces methyl branching |
| Substituted malonyl-CoAs | Less common extender units that introduce additional side-chain diversity |

Common PKS domains include:

| Domain | Function |
|--------|----------|
| KS | Ketosynthase; catalyzes chain extension |
| KSq | Decarboxylating KS-like loading domain |
| AT | Acyltransferase; selects the starter or extender substrate |
| ACP | Acyl carrier protein; carries the growing chain |
| KR | Ketoreductase; reduces beta-keto groups to alcohols |
| DH | Dehydratase; forms alkenes by dehydration |
| ER | Enoylreductase; reduces alkenes to saturated bonds |

The AT substrate choice and reductive domain pattern strongly influence the structure of the PKS product.

---

## Why TridentSynth is useful

Designing a PKS pathway is difficult because many different module combinations can produce related carbon skeletons. A target molecule may require decisions about:

- starter unit selection
- extender unit selection
- number of modules
- reductive domain pattern
- chain release mechanism
- biological or chemical tailoring after PKS assembly

TridentSynth is useful because it can search across multiple synthesis strategies:

1. PKS assembly
2. Biological tailoring
3. Chemical tailoring

This matters because the best pathway may not be PKS-only. A PKS may produce a close intermediate that requires a biological or chemical step to reach the final target.

---

## Tool workflow

The MCP tool follows this workflow:

```text
target SMILES
      │
      ▼
validate and normalize input
      │
      ▼
construct TridentSynth form payload
      │
      ▼
submit live job to TridentSynth
      │
      ▼
extract task ID
      │
      ▼
poll result page
      │
      ▼
parse best pathway
      │
      ▼
return structured JSON and text_summary
```

The tool submits to the live TridentSynth form using field names used by the website:

```text
smiles
synthesisStrategy_pks
synthesisStrategy_bio
synthesisStrategy_chem
rangebio
rangechem
releaseMechanism
pksStarters[]
pksExtenders[]
maxAtomsC
maxAtomsN
maxAtomsO
```

After submission, the tool extracts the task ID and polls:

```text
/run_TridentSynth/results/<task_id>/
```

until the result is complete or the timeout is reached.

---

## Input validation

The tool performs several checks before submitting a job.

### Single-molecule target

TridentSynth accepts one target molecule at a time. Therefore, the tool rejects SMILES strings containing a dot, such as:

```text
C.CC
```

A dot usually indicates multiple disconnected molecular species, which would make the target ambiguous.

### Strategy selection

At least one synthesis strategy must be selected:

- `use_pks`
- `use_bio`
- `use_chem`

If all three are false, the tool raises an error because there is no pathway search to run.

### Step limits

Biological and chemical step counts must be:

```text
1, 2, or 3
```

This matches the expected TridentSynth search controls and keeps searches compact.

### Atom filters

Optional atom filters such as `max_carbon`, `max_nitrogen`, and `max_oxygen` must be nonnegative when provided.

---

## Auto-filled search settings

When `auto_optimize_unspecified=True`, the tool fills blank optional parameters with a compact exploratory setup.

| Setting | Auto-filled value |
|---------|------------------|
| `max_bio_steps` | `1` |
| `max_chem_steps` | `1` |
| `pks_release_mechanism` | inferred from target |
| `pks_starters` | `mal`, `mmal` |
| `pks_extenders` | `mal`, `mmal` |

These defaults allow natural-language prompts to work without requiring the user to manually specify every TridentSynth option.

---

## PKS release mechanism

The PKS release mechanism controls how the completed chain leaves the assembly line.

| Release mechanism | Meaning |
|------------------|---------|
| `thiolysis` | Releases a mostly linear product, often as a carboxylic acid or related structure |
| `cyclization` | Releases the chain through ring formation |

The tool uses a simple heuristic to infer the release mechanism. If the target appears ring-like and contains heteroatoms, it may choose `cyclization`; otherwise, it defaults to `thiolysis`. The user can override this with `pks_release_mechanism`.

---

## Parsing strategy

The TridentSynth result page contains multiple related structures, so explicit parsing is important. Gemini may otherwise confuse:

- target SMILES
- PKS product SMILES
- final post-PKS product SMILES
- intermediate pathway structures

To avoid this, the tool assigns each field from a specific source:

| Output field | Source |
|-------------|--------|
| `target_smiles` | submitted user input |
| `pks_product_smiles` | PKS product label in the best pathway block |
| `post_pks_product_smiles` | final product label or target-matching reaction product |
| `reaction_smiles` | reaction SMILES section of the best pathway block |
| `pks_modules` | first/best PKS module set |
| `pathway_structures_smiles` | parsed PKS product plus reaction reactants/products |

The first pathway block is treated as the best pathway. This prevents later candidate pathways from being mixed into the primary output.

---

## Output interpretation

The tool returns structured JSON plus a deterministic `text_summary`.

Important top-level result fields include:

| Field | Meaning |
|------|---------|
| `synthesis_parameters` | TridentSynth settings parsed from the result page |
| `selected_steps` | Bio/Chem step counts that were actually selected |
| `pks_modules` | Parsed PKS module architecture |
| `best_pathway` | Parsed best pathway result |
| `text_summary` | Clean non-redundant summary for Gemini to display |

Important `best_pathway` fields include:

| Field | Meaning |
|------|---------|
| `pks_product_smiles` | Product generated by PKS before tailoring |
| `pks_similarity_to_target` | Similarity between PKS product and target |
| `post_pks_product_smiles` | Final product after Bio/Chem tailoring |
| `post_pks_similarity_to_target` | Similarity between final product and target |
| `reaction_smiles` | Reaction SMILES for the post-PKS transformation |
| `reaction_rule_names` | Human-readable transformation rule names |
| `reaction_enthalpies` | Reaction enthalpy values when available |
| `pathway_structures_smiles` | Unique SMILES structures in the selected pathway |

---

## Similarity scores

TridentSynth reports similarity scores comparing products to the target molecule. The tool preserves these as:

```text
pks_similarity_to_target
post_pks_similarity_to_target
```

These values distinguish how close the initial PKS product is from how close the final product is after tailoring.

For example:

```text
PKS similarity to target: 0.67
Final product similarity to target: 1.0
```

This means the PKS product is only partially similar to the target, but a post-PKS transformation reaches the final target.

---

## Reaction SMILES

Reaction SMILES represent a transformation in the format:

```text
reactants >> products
```

For example:

```text
CCCCCCC(=O)O>>CCCCCC.O=C=O
```

This represents decarboxylation of a carboxylic acid intermediate to form hexane and carbon dioxide.

The tool parses this into:

```json
{
  "reactants": ["CCCCCCC(=O)O"],
  "products": ["CCCCCC", "O=C=O"]
}
```

This makes the pathway easier for Gemini and downstream tools to interpret.

---

## Deterministic summarization

The `text_summary` field is designed to be displayed directly by Gemini. It reduces:

- repeated sections
- incorrect target SMILES
- confusion between PKS product and final product
- showing Bio steps when Bio was not selected
- unnecessary raw JSON interpretation

This is important because large tool outputs can contain many related structures. A deterministic summary makes the user-facing response more reliable.

---

## Example pathway interpretation

For a hexane target:

```text
Target SMILES: CCCCCC
```

TridentSynth may return a pathway where a PKS creates a carboxylic acid intermediate:

```text
CCCCCCC(=O)O
```

Then a chemical step performs decarboxylation:

```text
CCCCCCC(=O)O>>CCCCCC.O=C=O
```

This means:

1. The PKS produces a carboxylic acid intermediate.
2. A chemical reaction removes carbon dioxide.
3. The final product is hexane.
4. The final product similarity to target is 1.0.

This shows why combined PKS plus chemical synthesis can succeed even when PKS assembly alone does not directly produce the final target.

---

## GUI design rationale

The Streamlit GUI is not part of the TridentSynth chemistry algorithm, but it improves usability.

The GUI:

1. loads Gemini
2. starts the MCP client
3. lists available tools
4. sends user prompts to Gemini
5. lets Gemini call MCP tools
6. displays the final answer

The GUI uses `st.cache_resource` to keep a persistent MCP connection alive, avoiding server restarts after every user message. This makes the demo smoother and reduces overhead when running multiple prompts.

---

## Testing strategy

The pytest file tests the tool without depending on the live TridentSynth website by default. Most tests use fake HTML or internal helper functions, which keeps tests fast and stable.

The tests cover:

| Area | Reason |
|------|--------|
| Payload construction | Ensures correct TridentSynth form fields |
| Strategy selection | Ensures PKS/Bio/Chem flags are correctly submitted |
| Input validation | Prevents invalid SMILES and invalid step counts |
| Reaction parsing | Ensures reaction SMILES are split correctly |
| PKS module parsing | Ensures modules and domains are extracted |
| Result-page parsing | Ensures best pathway fields are populated |
| Target correction | Prevents target/intermediate SMILES mix-ups |
| Selected steps | Prevents unselected Bio/Chem steps from appearing |
| Optional live test | Allows real website submission when explicitly enabled |

The optional live test is skipped unless `RUN_LIVE_TRIDENTSYNTH=1` is set. This avoids flaky grading caused by network or server availability issues.

---

## Assumptions

The tool assumes:

1. TridentSynth keeps the same core form field names.
2. The result page continues to contain labels such as `PKS product`, `Post-PKS product`, `Reactions (SMILES)`, and `Reaction rules`.
3. The first pathway block is the best pathway.
4. The submitted target SMILES is the authoritative target identity.
5. If the final product has similarity 1.0 or appears as a reaction product, it should be treated as the target product.

---

## Known issues fixed during development

| Issue | Fix |
|------|-----|
| Tool originally returned only a results link | Added result-page polling and parsing |
| TridentSynth returned only a task ID after submission | Added task ID extraction and result URL construction |
| Tool was slower than manual website use | Switched to direct polling of `/results/<task_id>/` |
| Gemini confused target SMILES with intermediate SMILES | Forced `target_smiles` to use submitted input |
| Gemini showed Bio steps when Bio was not selected | Added `selected_steps` based on selected strategies |
| Gemini repeated the same pathway information | Added a clean `text_summary` |
| Pathway structures included bad strings | Added SMILES-like filtering and de-duplication |
| Parser returned multiple candidate pathways | Limited parsing to the first/best pathway block |
| Tests could have depended on the live website | Made live test optional |

---

## Limitations

- The tool depends on the live TridentSynth website being available.
- Runtime depends on TridentSynth server speed.
- Complex targets may take longer or time out.
- Website layout changes may require parser updates.
- Similarity and feasibility values are computational predictions.
- A predicted pathway is not proof that the pathway will work experimentally.
- PKS module designs should be reviewed before construct design.
- Post-PKS chemical transformations may not be biologically implementable without additional engineering.

---

## Summary

The `tridentsynth` MCP tool is a live wrapper around a specialized PKS pathway design platform. It allows Gemini to submit target molecules to TridentSynth, retrieve the best pathway, and return interpretable PKS design information.

The main technical contribution is not a new retrosynthesis algorithm, but a robust MCP integration layer. This layer handles validation, live submission, result polling, explicit parsing, output cleanup, deterministic summarization, and GUI-based interaction.

This makes TridentSynth usable as part of an AI-assisted PKS design workflow and connects natural-language prompts to real pathway design predictions.