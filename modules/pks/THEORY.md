# Theory — PKS Pipeline Algorithms and Biological Background

---

## 1. What is a Polyketide Synthase (PKS)?

A Type I modular PKS is a biological assembly line encoded by a cluster of genes. Each **module** in the assembly line adds one two-carbon unit (extender unit) to a growing chain, and each **domain** within a module performs a specific chemical transformation:

- **KS** (Ketosynthase) — condenses the growing chain with the next extender unit
- **AT** (Acyltransferase) — selects which extender unit to add (malonyl-CoA, methylmalonyl-CoA, etc.)
- **KR** (Ketoreductase) — optionally reduces the beta-keto group to a hydroxyl
- **DH** (Dehydratase) — optionally dehydrates the hydroxyl to a double bond
- **ER** (Enoylreductase) — optionally reduces the double bond fully
- **ACP** (Acyl Carrier Protein) — carries the growing chain between domains
- **TE** (Thioesterase) — releases the final product, usually by macrolactonization

The pattern of active/inactive domains across all modules determines the final product's carbon skeleton and stereochemistry. Well-known PKS natural products include erythromycin, rapamycin, epothilone, and pikromycin.

---

## 2. Why Search Biosynthetic Intermediates?

Most existing tools (e.g. ClusterCAD) only search **final products** — the compound released after the last module and TE domain. This misses a key engineering opportunity.

If your target molecule matches a **mid-pathway intermediate** — the chain state after module 3 of erythromycin, for example — then you can engineer that PKS to release the chain early by:
1. Relocating the TE domain to after the matching module
2. Deleting all downstream modules

`search_pks` searches all ~2,440 SBSPKS intermediates (every node in every pathway graph) so these truncation opportunities are surfaced automatically.

---

## 3. Chemical Similarity: Morgan Fingerprints and Tanimoto Coefficient

### Morgan Fingerprints (Extended Connectivity Fingerprints, ECFP)

A Morgan fingerprint encodes the local chemical environment around each atom in a molecule as a fixed-length bit vector. The algorithm:

1. Assign each atom an initial identifier based on its atomic number, charge, and degree
2. Iteratively update each atom's identifier by hashing it together with its neighbors' identifiers
3. After `radius` iterations (radius=2 → ECFP4), collect all generated identifiers
4. Map them into a bit vector of length `nBits` (2048 in this tool) using modular hashing

The result is a 2048-bit binary fingerprint where each bit represents the presence or absence of a particular local chemical substructure. This implementation uses RDKit's `MorganGenerator` with radius=2 and fpSize=2048.

### Tanimoto Coefficient

Given two fingerprints A and B, the Tanimoto (Jaccard) similarity is:

```
Tanimoto(A, B) = |A ∩ B| / |A ∪ B|
               = (bits set in both) / (bits set in either)
```

- Range: 0.0 (nothing in common) to 1.0 (identical)
- Scores >0.7 generally indicate the same compound family
- Scores 0.4–0.7 indicate related scaffolds worth engineering
- Scores <0.4 may reflect shared chemical class rather than true structural similarity

This is the industry-standard method for virtual screening and compound similarity search in drug discovery.

---

## 4. Database Construction

### SBSPKS v2
The SBSPKS server exposes pathway data via `make_reaction.cgi`, which returns the full biosynthetic graph for each pathway — every intermediate node with its SMILES. During `initiate()`, the tool:
1. Fetches all ~225 pathway pages
2. Extracts intermediate nodes and fetches their SMILES via `display_smiles.cgi`
3. Computes Morgan fingerprints for all ~2,440 entries
4. Caches to `modules/pks/data/sbspks_intermediate_index.json`

Note: The SBSPKS direct structure-search CGI endpoints (`search_similar1.cgi`, `search_functional.cgi`) have a persistent server-side permission bug and cannot be used. All similarity search is done locally.

### MIBiG
The MIBiG database is downloaded as a ZIP from the ClusterCAD GitHub mirror, parsed for all polyketide entries with SMILES, and cached to `modules/pks/data/mibig_cache.json`. MIBiG entries include `bgc_accession` fields (e.g. BGC0000055) usable directly with antiSMASH.

### Combined Index
Both sources are merged into `modules/pks/data/combined_index.json` (~4,088 entries total). The cache has a 30-day TTL and is rebuilt automatically on expiry.

---

## 5. Engineering Hints

When a query matches a biosynthetic intermediate (not a final product), `search_pks` generates an engineering hint explaining how to truncate that pathway to release the target structure:

> *"Target matches the module 7 intermediate of the Megalomycin pathway. Consider engineering this PKS to release the chain early by relocating the thioesterase (TE) domain after module 7, or by deleting all modules downstream of module 7."*

This is based on the well-established PKS engineering principle that TE domain relocation is sufficient to trigger early chain release at the new truncation point (Weissman & Leadlay, Nature Reviews Microbiology, 2005).

---

## 6. Limitations

- Similarity search uses 2D fingerprints — stereochemistry is not fully captured
- SBSPKS intermediates are predicted by the assembly-line logic, not experimentally verified
- Novel scaffolds with unusual extender units (e.g. methoxy-malonyl) may score low against all database entries
- The MIBiG entries include only compounds where SMILES are available; some BGC products lack structural data
- SBSPKS structure-search CGI endpoints (`search_similar1.cgi`, `search_functional.cgi`) have a persistent server-side permission bug and cannot be used for searching; however, `make_reaction.cgi` (pathway data) and `display_smiles.cgi` (SMILES retrieval) work correctly and are used for index building and pathway step fetching
- MIBiG hits do not include assembly line steps since MIBiG only stores the final compound and BGC metadata

---

## 7. Chemical Name Resolution (`resolve_smiles`) - Johanna Kann

Before any computational analysis can begin, the target molecule must be represented as a SMILES (Simplified Molecular Input Line Entry System) string — a compact, text-based notation that encodes molecular structure as a sequence of atoms, bonds, branches, and ring closures. For example, acetic acid is `CC(=O)O` and benzene is `c1ccccc1`.

Users often know a compound by its common name (e.g. "erythromycin"), trade name, or IUPAC systematic name rather than its SMILES. `resolve_smiles` bridges this gap using a two-step strategy:

1. **Local SMILES validation** — The input string is first parsed by RDKit's `MolFromSmiles`. If it produces a valid molecule object, the input is already SMILES and is canonicalized (rewritten into RDKit's unique canonical form) and returned immediately with no network call. Canonicalization ensures that equivalent SMILES strings (e.g. `OCC` and `CCO`) map to the same representation.

2. **PubChem REST lookup** — If the input is not valid SMILES, it is treated as a chemical name and submitted to the PubChem PUG REST API, which resolves the name against PubChem's database of >110 million compounds. The API returns the isomeric SMILES (preserving stereochemistry), molecular formula, and IUPAC name.

This tool is the mandatory first step in the pipeline: all downstream tools (`assess_pks_feasibility`, `pks_design_retrotide`, `pks_search_sbspks`) require SMILES input.

---

## 8. PKS Feasibility Pre-screening (`assess_pks_feasibility`) -- Johanna Kann

Not every molecule is a viable target for type-I modular PKS biosynthesis. `assess_pks_feasibility` applies seven weighted heuristic checks derived from the known biosynthetic capabilities and constraints of type-I modular PKS systems:

| Check | Weight | Rationale |
|-------|--------|-----------|
| **Molecular weight** | 0.12 | Known type-I PKS products range from ~100–2000 Da. Products outside this range require an impractical number of modules or are too small for meaningful PKS assembly. |
| **Carbon chain length** | 0.15 | Each PKS extension module adds a two-carbon unit; at least 6 carbons are needed for a minimal 3-module assembly line. |
| **Heteroatom content** | 0.14 | PKS assembles C/H/O backbones natively. Nitrogen, sulfur, and halogens require post-PKS tailoring enzymes (aminotransferases, halogenases) that fall outside the core PKS machinery. |
| **Aromatic rings** | 0.12 | Type-I modular PKS produces linear or macrocyclic scaffolds. Aromatic rings are the domain of type-II (iterative) PKS and cannot be formed by the modular assembly line, though aromatic starter units are possible. |
| **Oxygen/carbon ratio** | 0.16 | Polyketide intermediates typically have O/C ratios of 0.15–0.50, reflecting the ketone, hydroxyl, and ester functionalities introduced by KR/DH domains and TE-mediated lactonization. |
| **Functional groups** | 0.21 | PKS-compatible functional groups include ketones, hydroxyl groups, alkenes, ethers, and esters. Groups like epoxides, azides, nitro groups, and phosphates are not accessible by PKS domains. |
| **Ring complexity** | 0.10 | Single macrolactone rings are typical PKS products; 2–3 rings may be achievable with tailoring cyclases, but higher ring counts exceed typical PKS capability. |

Each check returns a status of **pass** (multiplier 1.0), **warn** (0.5), or **fail** (0.0). The final score is the sum of each check's weight multiplied by its status multiplier. Molecules scoring ≥0.6 are considered feasible PKS targets; scores ≥0.8 are strong candidates. The score guides which downstream design tools to prioritize.

---

## 9. Retrobiosynthetic PKS Design (`retrotide_designer`) -- Johanna Kann

### The RetroTide Algorithm

RetroTide performs **retrobiosynthesis**: given a target polyketide structure, it works backwards to propose chimeric type-I modular PKS assembly lines whose predicted product best matches the target. The algorithm, developed by Hagen et al. (2016), uses a combinatorial approach:

1. **Target decomposition** — The target SMILES is parsed and its carbon backbone is analyzed to determine the number of extension modules needed and the pattern of oxidation states at each position (ketone, hydroxyl, alkene, or fully reduced methylene).

2. **Module enumeration** — For each module position, RetroTide enumerates possible domain configurations from a library of characterized PKS domains. The key choices per module are:
   - **AT substrate** — which extender unit (malonyl-CoA, methylmalonyl-CoA, ethylmalonyl-CoA, etc.) introduces the correct branching pattern
   - **KR type and activity** — whether the beta-keto group is reduced and with which stereochemistry (A-type vs. B-type, yielding L- or D-configured alcohols)
   - **DH activity** — whether dehydration occurs, producing an enoyl intermediate
   - **ER activity** — whether the double bond is fully reduced

3. **Product prediction** — Each candidate assembly line is fed through a structure prediction engine (BCS — Biosynthetic Cluster Simulator) that computes the product SMILES by simulating the sequential condensation, reduction, and release reactions.

4. **Similarity ranking** — The predicted product of each design is compared to the target using chemical similarity metrics. Two scoring options are available:
   - **Atom-pair similarity** (default) — counts the frequency of all pairs of atoms at given topological distances in the molecular graph
   - **Atom-atom-path similarity** — considers the full atom-to-atom paths rather than just pairwise distances, providing a more path-sensitive comparison

Designs are ranked by similarity score (1.0 = exact match) and the top candidates are returned. Each result includes the full module-by-module domain architecture, enabling direct implementation.

### Chimeric PKS Engineering -- Johanna Kann

The designs RetroTide proposes are "chimeric" — they combine domains sourced from different natural PKS clusters. This is biologically grounded: domain swapping experiments have demonstrated that individual domains (particularly AT domains) can be exchanged between clusters while retaining function (Ruan et al., 1997; Reeves et al., 2001). RetroTide selects domains based on their characterized substrate specificity and reduction activity, drawing from a curated library of experimentally validated PKS domains.

---

## 10. Mapping Designs to Natural Parts (`find_pks_module_parts`) -- Johanna Kann

A PKS design from RetroTide (or TridentSynth) specifies an abstract domain architecture — e.g., "module 2 needs an AT that loads methylmalonyl-CoA, an active B-type KR, and an active DH." To build this in the lab, each abstract domain must be matched to a real amino acid sequence from a characterized natural PKS.

`find_pks_module_parts` automates this mapping one module at a time using the ClusterCAD database, which contains 531 experimentally characterized PKS clusters with full domain annotations. The tool accepts simple scalar inputs — whether the module is a loading module, its AT substrate, and which reductive domains (KR, DH, ER) are active — so it can be called independently for each module in any design:

1. **Per-module search** — Given a module's properties (loading status, AT substrate specificity, active reductive domains), the tool queries ClusterCAD's domain index for natural modules with matching characteristics. This approach works with designs from any upstream tool without requiring format-specific normalization.

2. **Amino acid sequence retrieval** — For each matching natural module, the tool fetches the amino acid sequence of every domain from ClusterCAD's API. These sequences are the starting point for gene synthesis or domain swapping constructs.

3. **Result assembly** — The output provides a ranked list of natural module matches with their source cluster, subunit, and per-domain amino acid sequences. This gives the engineer multiple options for each position in the chimeric PKS, enabling selection based on factors like phylogenetic compatibility or expression host preferences.

This step connects computational design to physical DNA construction, closing the gap between in silico pathway prediction and wet-lab implementation.

---

---

## 11. Codon Optimization and DNA Synthesis (`reverse_translate`) — Opshory Choudhury

### Why Codon Usage Matters

All organisms use the same 20 amino acids but the 64-codon genetic code is degenerate — most amino acids are encoded by 2–6 synonymous codons. However, cells do not use all codons equally. Each organism maintains a characteristic **codon usage bias** driven by the relative abundance of cognate tRNAs. Rare codons slow ribosome elongation, cause premature translation termination, increase the probability of misincorporation errors, and reduce overall protein yield.

For heterologous expression of PKS genes — which already encode very large, multi-domain proteins — codon usage mismatches are a major source of failure. A DEBS module 1 gene optimized for *E. coli* performs poorly in *Streptomyces coelicolor* because the two organisms have radically different GC contents (~51% vs ~72%) and codon preferences. Proper codon optimization for the target host is therefore essential before DNA synthesis.

### The DnaChisel Algorithm

`reverse_translate` uses **DnaChisel** (Zulkower & Rosser, 2020) to convert an amino acid sequence into a codon-optimized DNA sequence for a chosen host. The optimization proceeds in two steps:

1. **Constraint enforcement — `EnforceTranslation`**
   The amino acid sequence is locked in place. DnaChisel builds an initial DNA string using a naive per-residue codon lookup, then treats the correct protein sequence as a hard constraint during all subsequent optimization. This guarantees that no optimization step alters the encoded protein.

2. **Objective maximization — `CodonOptimize`**
   Given the host organism's codon usage table (sourced from the Codon Usage Database via NCBI taxon ID), DnaChisel maximizes the Codon Adaptation Index (CAI) of the DNA sequence. CAI is defined as the geometric mean of the relative synonymous codon usage (RSCU) values across all codons in the sequence:

   ```
   CAI = exp( (1/L) × Σ ln(w_i) )
   ```

   where L is the number of codons and w_i is the RSCU of codon i relative to the most used synonymous codon for that amino acid. A CAI of 1.0 means every codon is the most preferred codon for that amino acid in that organism.

   DnaChisel optimizes CAI using a local search over synonymous codon substitutions, respecting the `EnforceTranslation` constraint throughout.

### Supported Hosts

| Alias | Organism | NCBI Taxon ID | Notable use case |
|-------|----------|---------------|-----------------|
| `e_coli` | *Escherichia coli* | built-in | Rapid prototyping, solubility testing |
| `s_coelicolor` | *Streptomyces coelicolor* | 100226 | Native PKS host, high GC content |
| `s_albus` | *Streptomyces albus* | 318586 | Clean genetic background for PKS expression |
| `s_venezuelae` | *Streptomyces venezuelae* | 54571 | Fast-growing Streptomyces, pikromycin producer |
| `p_putida` | *Pseudomonas putida* | 160488 | Robust heterologous expression |
| `s_cerevisiae` | *Saccharomyces cerevisiae* | built-in | Eukaryotic expression, post-translational modifications |

### Why CDS Annotation Matters

The GenBank output includes a `CDS` feature covering the full coding sequence, with an embedded `translation` qualifier. This is critical for downstream antiSMASH analysis: antiSMASH's domain annotation pipeline requires annotated CDS features to identify where genes begin and end. Without a CDS feature, antiSMASH falls back to Prodigal for ab initio gene prediction — which performs poorly on short single-ORF constructs and on high-GC Streptomyces-optimized sequences that may lack canonical Shine-Dalgarno sequences.

---

## 12. Submitting Constructs to antiSMASH (`submit_antismash`) — Opshory Choudhury

### What antiSMASH Does

**antiSMASH** (antibiotics and Secondary Metabolite Analysis SHell) is the field-standard platform for identifying and annotating biosynthetic gene clusters (BGCs) in genomic sequences (Blin et al., 2023). For PKS engineering, it provides two key services:

1. **Domain annotation** — Each ORF is scanned against a curated library of profile Hidden Markov Models (pHMMs) built from experimentally characterized PKS domains. Every detected domain (KS, AT, KR, DH, ER, ACP, TE, etc.) is reported with its position, HMM score, and E-value. Low E-values (e.g. 1.9×10⁻¹⁷³ for a KS domain) indicate confident detection; partial or truncated domains have elevated E-values and reduced bit scores, serving as an immediate flag for assembly errors.

2. **Product structure prediction** — After domain detection, antiSMASH's NRPS/PKS module predicts the core polyketide scaffold by walking the modules in order and simulating the chemical transformations at each step (see Section 13).

### The Three Submission Modes

`submit_antismash` supports three input modes that differ in how gene boundaries are determined:

| Mode | Input | Gene finding | Best for |
|------|-------|-------------|---------|
| `seq` | Raw DNA string | Prodigal ab initio | Raw sequences without annotations |
| `filepath` | GenBank `.gb` file | Uses existing CDS annotations | In-silico designs from `reverse_translate` |
| `ncbi` | NCBI accession | Uses existing NCBI annotations | Sequenced clones deposited on NCBI |

The `filepath` and `ncbi` modes are preferred because they provide precise gene boundaries, ensuring domain annotation begins at the correct reading frame. Prodigal (`seq` mode) requires a minimum of ~1000 bp and performs gene prediction de novo, which can miss ORFs in short constructs or misidentify the start codon.

### Always-On Analyses

Every submission enables two additional analyses:

- **Active Site Finder (ASF)** — Scans each detected domain for the catalytic residue motifs that determine activity. For KR domains, it distinguishes active (NADPH-binding Lys and Ser motifs present) from inactive (degenerate motifs) variants. For AT domains, it identifies the catalytic Ser that forms the acyl-enzyme intermediate. Missing active site residues indicate a domain that will be non-functional even if it is structurally intact.

- **KnownClusterBlast** — Compares every protein in the submitted sequence against the MIBiG database of experimentally characterized BGCs using DIAMOND protein search. Results show which known clusters share homologous proteins with your construct and at what percent identity. For a codon-optimized DEBS construct, this should return a high-identity hit to BGC0000055 (erythromycin), confirming the amino acid sequence is correct regardless of synonymous codon changes.

---

## 13. Parsing and Validating antiSMASH Results (`check_antismash`) — Opshory Choudhury

### AT Substrate Prediction

The AT domain determines which extender unit is loaded onto the ACP — the fundamental choice that sets the branching pattern and functional group at each position in the polyketide chain. antiSMASH predicts AT substrate specificity using two independent methods whose results are reported side by side:

1. **Signature-based prediction (Rausch et al., 2005)** — A set of 24 residue positions in the AT active site, collectively called the "AT signature," are extracted from the aligned domain sequence. These positions are compared against a reference library of AT domains with known substrate specificity using a scoring matrix. The method returns a ranked list of substrate predictions with percent identity scores. Substrates predicted at ≥85% identity are highly reliable; those at 60–85% should be treated as tentative.

2. **Minowa profile HMMs (Minowa et al., 2007)** — Each AT domain is scored against a library of substrate-specific pHMMs. This method is complementary to signature-based prediction and particularly useful for unusual extender units (e.g. ethylmalonyl-CoA, methoxymalonyl-CoA) that are underrepresented in the signature library.

`check_antismash` reports the top signature-based prediction as `AT_substrate` (translated to a human-readable name like "Methylmalonyl-CoA") along with the raw antiSMASH code (`AT_substrate_code`, e.g. "mmal") for direct comparison against RetroTide's substrate nomenclature.

### KR Stereochemistry Prediction

The KR domain determines the stereochemistry of the hydroxyl group introduced at the beta-carbon. antiSMASH classifies KR domains into types based on the pattern of conserved residues in two diagnostic sequence motifs:

- **Type A KR** — produces L-configured (S) hydroxyl groups; associated with the LDD motif in the "lid" region
- **Type B KR** — produces D-configured (R) hydroxyl groups; associated with a different lid sequence
- **Subtypes A1, A2, B1, B2** — further distinguish whether the KR is C2-methylated (A1/B1) or unmethylated (A2/B2), refining the stereochemical prediction
- **Type C** — non-stereospecific or inactive; associated with degenerate motifs

The prediction accuracy is ~70–80% (Keatinge-Clay, 2007). Mismatches between expected and detected KR type are flagged in the `validation` output and indicate that the incorporated natural KR domain may produce the wrong stereoisomer.

### PKS Structure Prediction

After domain annotation, antiSMASH predicts the core polyketide scaffold by simulating the assembly line chemistry module by module:

1. The loading module (or starter unit) is identified and its starter CoA thioester is recorded
2. Each extension module's AT substrate adds its extender unit to the growing chain
3. The reductive loop (KR/DH/ER) determines the oxidation state at each new beta-carbon:
   - No reductive domains → beta-keto (C=O)
   - KR only → beta-hydroxy (C-OH, stereochemistry from KR type)
   - KR + DH → alpha,beta-unsaturated (C=C, E or Z from DH type)
   - KR + DH + ER → fully reduced methylene (C-C)
4. The terminal TE domain determines cyclization (macrolactonization for most PKS products) or hydrolysis

The result is reported as a polymer shorthand (e.g. `(Me-ohmal)` for a methylated, hydroxylated malonate extension) and as a SMILES string. This predicted SMILES is the first approximation of what the engineered PKS will produce — useful for comparing against the RetroTide design intent.

### Domain Order String and Design Validation

`check_antismash` parses all `aSDomain` features from the antiSMASH output, sorts them by their genomic position, and builds an ordered domain string per gene (e.g. `KS-AT-KR-ACP`). This string is the most compact summary of what the annotated construct actually encodes.

When `expected_domains` is provided, the tool computes a structured diff:
- **Missing** domains — present in the design but not detected by antiSMASH; indicates a frame-shift, synthesis error, truncated ORF, or domain below the detection threshold
- **Unexpected** domains — detected by antiSMASH but not in the design; indicates an unintended ORF, a retained upstream/downstream sequence, or a mis-assembled construct

This automated comparison replaces the manual inspection of antiSMASH HTML output and provides a machine-readable validation signal directly in the pipeline.


---

## 14. ClusterCAD Database Tools (`clustercad_*`) - Dennis Wu

### What is ClusterCAD?

ClusterCAD is a curated database of 531 experimentally characterized Type I modular PKS clusters, each with full domain annotations, substrate specificities, subunit sequences, and predicted biosynthetic intermediates. It is maintained by the Joint BioEnergy Institute (JBEI) and provides a ground-truth reference for PKS domain architectures derived from real organisms.

### Role in the Pipeline

The ClusterCAD tools serve as a **parts database** for PKS engineering. Once RetroTide or TridentSynth proposes an abstract domain architecture for a target molecule, the ClusterCAD tools find real, experimentally validated biological parts that implement each proposed module. This closes the gap between computational design and physical DNA construction.


RetroTide / TridentSynth → proposes abstract domain architecture
↓
clustercad_search_domains → finds natural modules matching each step
↓
clustercad_get_subunits → retrieves full domain architecture and IDs
↓
clustercad_domain_lookup / clustercad_subunit_lookup → retrieves AA and DNA sequences
↓
reverse_translate + submit_antismash → constructs and validates the chimeric PKS

### The Six ClusterCAD Tools

#### `clustercad_list_clusters`
Returns a browsable list of PKS clusters from ClusterCAD with their MIBiG accession numbers and descriptions. Used as the entry point when the user references a cluster by name rather than accession number. Supports filtering by reviewed status and result count.

#### `clustercad_cluster_details`
Returns a summary of a specific PKS cluster given its MIBiG accession number, including subunit count, module count, and a link to the ClusterCAD web page. Used to get a quick overview before retrieving the full domain architecture.

#### `clustercad_get_subunits`
Returns the complete domain architecture for every subunit and module in a PKS cluster. For each module, it provides:
- All domain types present (KS, AT, KR, DH, ER, ACP, TE)
- Domain IDs for downstream sequence retrieval
- Substrate annotations (e.g. "substrate mal, loading", "type B1, active")
- Predicted intermediate SMILES — the growing chain state after each module

This is the key tool for tracing how a natural PKS builds its product step by step, and for identifying which module produces an intermediate matching the engineering target.

#### `clustercad_domain_lookup`
Returns the amino acid sequence and positional coordinates of a single domain given its domain ID. Domain IDs are obtained from `clustercad_get_subunits`. This tool provides the sequence-level information needed to design domain-swapping constructs for chimeric PKS engineering.

#### `clustercad_subunit_lookup`
Returns both the amino acid sequence and the full nucleotide (DNA) sequence for an entire PKS subunit, along with its GenBank accession number. Subunit IDs are obtained from `clustercad_get_subunits`. This tool is the final step before gene synthesis — providing the exact DNA sequence to order or codon-optimize for a given expression host.

#### `clustercad_search_domains`
Searches all 531 PKS clusters simultaneously for modules matching specified domain criteria. Uses a locally cached JSON index (built once at first run, ~70 seconds) for instant subsequent queries. Supports filtering by:
- Domain type (AT, KS, KR, DH, ER, ACP, TE)
- Substrate annotation (with automatic synonym translation — e.g. "butyryl" → "butmal")
- Domain combination (e.g. modules containing KR + DH + ER all active)
- Active/inactive status
- Loading module vs. extension module
- Cluster name, module count, and reviewed status

### Substrate Annotation Vocabulary

ClusterCAD uses abbreviated substrate names in its domain annotations. The search tool automatically translates common chemical names to ClusterCAD terms:

| Common Name | ClusterCAD Term |
|------------|----------------|
| Malonyl-CoA | mal |
| Methylmalonyl-CoA | mmal |
| Ethylmalonyl-CoA | emal |
| Butyryl / Butylmalonyl-CoA | butmal |
| Hexylmalonyl-CoA | hxmal |
| Hydroxymalonyl-CoA | hmal |
| Methoxymalonyl-CoA | mxmal |
| Isobutyryl-CoA | isobut |
| Propionyl-CoA | prop |
| Acetyl-CoA | Acetyl-CoA |
| Pyruvate | pyr |

### Intermediate SMILES and Engineering Opportunities

A key feature of the ClusterCAD tools is access to the predicted intermediate SMILES at each module — the chain state after that module has completed its condensation and reductive cycle. This enables the same truncation engineering strategy described in Section 2 and Section 5 for SBSPKS: if a target molecule matches a mid-pathway intermediate in a natural PKS cluster, the engineer can:

1. Relocate the TE domain to after the matching module
2. Delete all downstream modules
3. Retrieve the subunit DNA sequence using `clustercad_subunit_lookup` for direct gene synthesis

This approach is grounded in the same principle as SBSPKS intermediate search (Section 2) but uses the curated, experimentally validated ClusterCAD dataset rather than the computationally predicted SBSPKS intermediates.

---

## 15. Live TridentSynth Pathway Search (`tridentsynth`) — Shalen Ardeshna

### What TridentSynth Adds to the Pipeline

`tridentsynth` connects the MCP pipeline to the live TridentSynth web server. While tools like `retrotide_designer` generate abstract PKS module designs and `search_pks` searches known natural products/intermediates, TridentSynth searches for pathways that combine:

1. **PKS assembly**
2. **Biological tailoring**
3. **Chemical tailoring**

This is useful because many target molecules are not direct PKS products. A PKS may produce a close intermediate that must then be converted to the target through one biological or chemical step.

### Live Job Submission

The tool accepts a target SMILES string and builds a TridentSynth-compatible form payload using the same field names as the website:

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

The main user-controlled parameters are whether to use PKS, Bio, and/or Chem synthesis; the maximum number of Bio/Chem steps; the PKS release mechanism; and the PKS starter/extender substrates. If optional values are not provided, the tool auto-fills a compact exploratory setup using one Bio step, one Chem step, and common malonyl/methylmalonyl starter and extender units.

### Result Polling and Parsing

After submission, TridentSynth returns a task ID. The tool polls:

```text
/run_TridentSynth/results/<task_id>/
```

until the result page is complete or the timeout is reached.

The parser extracts the first/best pathway block and returns:

- synthesis parameters
- PKS modules
- PKS product SMILES
- post-PKS product SMILES
- reaction SMILES
- reaction rule names
- reaction enthalpy
- similarity scores
- pathway structure SMILES

Only the first/best pathway block is used so that alternate candidate pathways do not get mixed into the main result.

### Preventing SMILES Misinterpretation

A TridentSynth result page can contain several related molecules: the submitted target, the PKS product, final post-PKS product, byproducts, and alternate pathway products. These can be easy for Gemini to confuse if it reads the raw JSON directly.

To avoid this, the tool explicitly assigns each output field from a specific result-page section. It also creates a deterministic `text_summary` that Gemini can display directly. This prevents common errors such as reporting an intermediate as the target, showing Bio steps when Bio was not selected, or repeating the same pathway information multiple times.

### Example: Hexane Pathway

For hexane:

```text
CCCCCC
```

TridentSynth may predict a PKS-generated carboxylic acid intermediate followed by chemical decarboxylation:

```text
CCCCCCC(=O)O>>CCCCCC.O=C=O
```

This means the PKS creates a carboxylic acid precursor, the chemical step removes carbon dioxide, and the final product is hexane. The final product can therefore match the target even when the direct PKS product only partially matches it.

### Limitations

- The tool depends on the live TridentSynth website being available.
- Runtime depends on TridentSynth server speed.
- Complex targets may take longer or time out.
- If the TridentSynth page layout changes, the parser may need to be updated.
- TridentSynth pathways are computational predictions and still require manual review and experimental validation.
- Post-PKS chemical transformations may not be directly biologically implementable without additional engineering.

---

## References

- Keatinge-Clay, A.T. (2012). The structures of type I polyketide synthases. *Natural Product Reports*, 29, 1050.
- Weissman, K.J. & Leadlay, P.F. (2005). Combinatorial biosynthesis of reduced polyketides. *Nature Reviews Microbiology*, 3, 925–936.
- Rogers, D. & Hahn, M. (2010). Extended-connectivity fingerprints. *Journal of Chemical Information and Modeling*, 50(5), 742–754.
- Tanimoto, T.T. (1958). An elementary mathematical theory of classification and prediction. IBM Technical Report.
- MIBiG Consortium (2023). MIBiG 3.0: a community-driven effort to annotate experimentally validated biosynthetic gene clusters. *Nucleic Acids Research*, 51, D603–D610.
- Weininger, D. (1988). SMILES, a chemical language and information system. *Journal of Chemical Information and Computer Sciences*, 28(1), 31–36.
- Kim, S. et al. (2021). PubChem in 2021: new data content and improved web interfaces. *Nucleic Acids Research*, 49, D1388–D1395.
- Hagen, A. et al. (2016). Engineering a polyketide synthase for in vitro production of adipic acid. *ACS Synthetic Biology*, 5(1), 21–27.
- Reeves, C.D. et al. (2001). Alteration of the substrate specificity of a modular polyketide synthase acyltransferase domain through site-directed mutagenesis. *Biochemistry*, 40(51), 15464–15470.
- Eng, C.H. et al. (2018). ClusterCAD: a computational platform for type I modular polyketide synthase design. *Nucleic Acids Research*, 46(D1), D509–D515.
- Zulkower, V. & Rosser, S. (2020). DNA Chisel, a versatile sequence optimizer. *Bioinformatics*, 36(16), 4508–4509.
- Blin, K. et al. (2023). antiSMASH 7.0: new and improved predictions for detection, regulation, and visualisation. *Nucleic Acids Research*, 51(W1), W46–W50.
- Rausch, C. et al. (2005). Specificity prediction of adenylation domains in nonribosomal peptide synthetases (NRPS) using transductive support vector machines (TSVMs). *Nucleic Acids Research*, 33(18), 5799–5808.
- Minowa, Y. et al. (2007). Comprehensive analysis of distinctive polyketide and nonribosomal peptide structural motifs encoded in microbial genomes. *Journal of Molecular Biology*, 368(5), 1500–1517.
- Keatinge-Clay, A.T. (2007). A tylosin ketoreductase reveals how chirality is determined in polyketides. *Chemistry & Biology*, 14(8), 898–908.
