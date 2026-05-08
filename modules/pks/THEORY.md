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
- SBSPKS live endpoints are unavailable due to a server-side bug; the index reflects the last successful cache build

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

## 10. Mapping Designs to Natural Parts (`match_design_to_parts`) -- Johanna Kann

A PKS design from RetroTide (or TridentSynth) specifies an abstract domain architecture — e.g., "module 2 needs an AT that loads methylmalonyl-CoA, an active B-type KR, and an active DH." To build this in the lab, each abstract domain must be matched to a real amino acid sequence from a characterized natural PKS.

`match_design_to_parts` automates this mapping using the ClusterCAD database, which contains 531 experimentally characterized PKS clusters with full domain annotations:

1. **Design normalization** — The tool accepts designs from either RetroTide or TridentSynth, which use different output formats. It normalizes both into a common representation: an ordered list of modules, each with its domain types and AT substrate annotation.

2. **Module-by-module search** — For each module in the design, the tool queries ClusterCAD's domain index for natural modules with matching characteristics:
   - AT substrate specificity (e.g., methylmalonyl-CoA)
   - Presence of required reductive domains (KR, DH, ER)
   - Loading module status (starter modules use different domain architectures)

3. **Amino acid sequence retrieval** — For each matching natural module, the tool fetches the amino acid sequence of every domain from ClusterCAD's API. These sequences are the starting point for gene synthesis or domain swapping constructs.

4. **Result assembly** — The output provides, for each design module, a ranked list of natural module matches with their source cluster, subunit, and per-domain amino acid sequences. This gives the engineer multiple options for each position in the chimeric PKS, enabling selection based on factors like phylogenetic compatibility or expression host preferences.

This step connects computational design to physical DNA construction, closing the gap between in silico pathway prediction and wet-lab implementation.

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
