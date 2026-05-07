# Theory — search_pks Algorithm and Biological Background

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

## References

- Keatinge-Clay, A.T. (2012). The structures of type I polyketide synthases. *Natural Product Reports*, 29, 1050.
- Weissman, K.J. & Leadlay, P.F. (2005). Combinatorial biosynthesis of reduced polyketides. *Nature Reviews Microbiology*, 3, 925–936.
- Rogers, D. & Hahn, M. (2010). Extended-connectivity fingerprints. *Journal of Chemical Information and Modeling*, 50(5), 742–754.
- Tanimoto, T.T. (1958). An elementary mathematical theory of classification and prediction. IBM Technical Report.
- MIBiG Consortium (2023). MIBiG 3.0: a community-driven effort to annotate experimentally validated biosynthetic gene clusters. *Nucleic Acids Research*, 51, D603–D610.
