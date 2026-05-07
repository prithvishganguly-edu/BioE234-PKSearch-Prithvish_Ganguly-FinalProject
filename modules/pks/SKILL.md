---
name: pks
description: The `pks` module provides tools for polyketide synthase (PKS) design and chemical name resolution, supporting retrobiosynthesis workflows.
---

# pks — Skill Guidance for Gemini

This file is read by the client at startup and injected into Gemini's system prompt.
It gives Gemini domain knowledge to use the PKS tools correctly and present their
results fully.

## What is a PKS?
A Type I modular PKS is a biological assembly line that builds complex natural
products step by step. Each module adds one building block (extender unit) to
a growing chain, and each domain within a module performs a specific chemical
transformation:
- **KS** (Ketosynthase): Condenses the growing chain with the extender unit
- **AT** (Acyltransferase): Selects and loads the extender unit (e.g. malonyl-CoA, methylmalonyl-CoA)
- **KR** (Ketoreductase): Reduces the ketone group (active/inactive, type A/B)
- **DH** (Dehydratase): Dehydrates the hydroxyl group (active/inactive)
- **ER** (Enoylreductase): Fully reduces the double bond (active/inactive)
- **ACP** (Acyl Carrier Protein): Carries the growing chain between domains
- **TE** (Thioesterase): Releases the final product

---

## Tools and when to use them

### `resolve_smiles`
Converts a chemical name (common, trade, or IUPAC) to a canonical SMILES string
via PubChem. If the input is already valid SMILES, returns immediately.

Use when:
- The user provides a chemical name instead of SMILES (e.g. "erythromycin")
- Another tool requires SMILES input but the user gave a name
- You need to verify or canonicalize a SMILES string

Returns: `smiles`, `molecular_formula`, `iupac_name`, `source`, and `cid`.

**Always call this before any tool that requires SMILES input.**

---

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
TridentSynth Live Pathway Tool

Connects Gemini to the live TridentSynth web server to generate PKS-based synthesis
pathway predictions for a target molecule.

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
- Best pathway found including PKS modules, domains, and AT substrate choices.
- PKS product SMILES and similarity to the target.
- Post-PKS product SMILES and reaction SMILES for best post-PKS transformation.
- Reaction rule name, enthalpy, and feasibility values when available.

---

### `dna_gc_content`
Computes the GC content (fraction of G and C bases) of a DNA sequence.

---

### ClusterCAD Tools
You have 6 tools that work together to explore the ClusterCAD database of 531
real PKS clusters. These tools find natural biological parts that can be used
to implement PKS designs proposed by RetroTide or TridentSynth.

#### `clustercad_list_clusters`
Lists PKS clusters from ClusterCAD with their MIBiG accessions and descriptions.

Use when:
- The user asks to browse available PKS clusters
- You need a MIBiG accession number for a named cluster

Parameters: `reviewed_only` (default True), `max_results` (default 20)

#### `clustercad_cluster_details`
Returns summary information (subunit count, module count, URL) for a specific cluster.

Use when:
- You need a quick overview of a cluster before diving into its architecture

Parameters: `mibig_accession` (e.g. "BGC0001492.1")

#### `clustercad_get_subunits`
Returns the full domain architecture for every subunit and module in a cluster,
including domain IDs, annotations, and intermediate SMILES for each module.

Use when:
- You need to see the full domain layout of a cluster
- You need domain IDs for `clustercad_domain_lookup`
- You need subunit IDs for `clustercad_subunit_lookup`
- You want to trace the intermediate chain SMILES step by step

Parameters: `mibig_accession`

#### `clustercad_domain_lookup`
Returns the amino acid sequence and positional details for a specific domain.

Use when:
- The user asks for the amino acid sequence of a specific domain
- You need the sequence of a particular AT, KR, DH, or other domain

Parameters: `domain_id` (integer, found via `clustercad_get_subunits`)

#### `clustercad_subunit_lookup`
Returns both the amino acid AND nucleotide (DNA) sequence for a full subunit,
plus its GenBank accession.

Use when:
- The user asks for the DNA or protein sequence of a full subunit
- You need a GenBank accession for a subunit

Parameters: `subunit_id` (integer, found via `clustercad_get_subunits`)

#### `clustercad_search_domains`
Searches all 531 PKS clusters instantly (using a local cache) for modules
matching specified domain type, substrate annotation, and other criteria.

Use when:
- The user wants to find natural PKS modules with specific properties
- RetroTide proposes a domain architecture and you want real natural examples
- The user asks "which PKS clusters use X substrate?"

Parameters:
- `domain_type`: e.g. "AT", "KR", "DH", "ER"
- `annotation_contains`: substrate or activity (auto-translated — see below)
- `domain_types`: list of domain types ALL required in module e.g. ["KR","DH","ER"]
- `exclude_annotation`: e.g. "inactive"
- `active_only`: True/False
- `loading_module_only`: True/False
- `cluster_description_contains`: filter by cluster name
- `min_modules` / `max_modules`: filter by cluster complexity
- `reviewed_only`: True/False
- `max_results`: default 10
- `force_refresh`: rebuilds cache if needed

**Substrate name auto-translation (no need to know ClusterCAD terms):**
| Common Name | ClusterCAD Term |
|------------|----------------|
| malonyl-CoA | mal |
| methylmalonyl-CoA | mmal |
| ethylmalonyl-CoA | emal |
| butyryl / butylmalonyl-CoA | butmal |
| hexylmalonyl-CoA | hxmal |
| hydroxymalonyl-CoA | hmal |
| methoxymalonyl-CoA | mxmal |
| isobutyryl-CoA | isobut |
| propionyl-CoA | prop |
| acetyl-CoA | Acetyl-CoA |
| pyruvate | pyr |

---

## Key Rules — ALWAYS Follow These

1. **NEVER ask the user for MIBiG accession numbers, domain IDs, or subunit IDs** — look them up yourself using the tools
2. **ALWAYS call `resolve_smiles` first** if the user provides a chemical name instead of SMILES
3. **ALWAYS call `clustercad_list_clusters` first** if you need an accession number for a named cluster
4. **ALWAYS call `clustercad_get_subunits` first** to find domain IDs before calling `clustercad_domain_lookup`
5. **ALWAYS call `clustercad_get_subunits` first** to find subunit IDs before calling `clustercad_subunit_lookup`
6. **Chain tools proactively** — never stop halfway and ask the user for information you can retrieve yourself

---

## Workflow Guidance

### "Design a PKS for [molecule]"
1. `resolve_smiles(molecule name)` → get SMILES
2. `pks_design_retrotide(smiles)` → get proposed domain architecture
3. `clustercad_search_domains(domain_type, annotation)` → find natural examples of each proposed module
4. Present full results including all module domains

### "Tell me about the Erythromycin PKS"
1. `clustercad_list_clusters(reviewed_only=True)` → find accession
2. `clustercad_cluster_details(accession)` → get subunit/module count
3. `clustercad_get_subunits(accession)` → get full domain architecture

### "Find PKS modules that load butyryl-CoA"
1. `clustercad_search_domains(domain_type="AT", annotation_contains="butyryl", loading_module_only=True)`
   → auto-translates butyryl → butmal, returns matching modules instantly

### "Give me the DNA sequence of a subunit"
1. `clustercad_list_clusters()` → find accession
2. `clustercad_get_subunits(accession)` → find subunit_id
3. `clustercad_subunit_lookup(subunit_id)` → returns full DNA + AA sequence

### "Find modules with fully reducing loops"
1. `clustercad_search_domains(domain_types=["KR", "DH", "ER"], active_only=True)`
   → returns modules with all three reductive domains active

### "Find alternative pathways"
1. `resolve_smiles(molecule name)` → get SMILES
2. `tridentsynth(target_smiles)` → get PKS + tailoring pathway

---

## How ClusterCAD Complements RetroTide and TridentSynth
- **RetroTide / TridentSynth** propose a PKS domain architecture for a target molecule
- **ClusterCAD tools** find natural examples of each proposed module from 531 real clusters
- **Together** they provide verified biological parts for PKS engineering instead of relying on purely computational designs