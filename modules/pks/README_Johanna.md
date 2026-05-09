# PKS De-Novo Design — resolve_smiles, assess_pks_feasibility, retrotide_designer, and find_pks_module_parts

**Author:** Johanna Kann 
**Course:** BioE 234, Spring 2026  

---

## What these tools do

This set of three tools forms the **preliminary design** stage of the PKS design pipeline. 

First, users input a natural language query for a molecule -- ranging from heptane to erythromycin to nitrogen. The first step calls **resolve_smiles** to convert this to a SMILES format, which queries PubChem to convert from a natural language name to a SMILES, which is a format required by downstream tools. The idea of adding this tool framework is to make sure the LLM doesn't hallucinate a SMILES for the molecule.

Next, the platform uses **assess_pks_feasibility** to determine if the inputted molecule could actuallly be made by a PKS. It scores the molecule based on a variety of qualities, which are described in more detail below, and rates its PKS-ability from 0 to 1. This helps direct whether the query goes to retrotide (pure PKS designer), Search PKS (natural PKS searcher), or tridentsynth (PKS + biology + chemistry). If the molecule is a great candidate for a PKS, the LLM prioritizes retrotide; if not, it prioritizes tridentsynth.

Next, if it is directed towards retrotide, the SMILES is inputed by the LLM into the **retrotide_designer** tool.This tool utilizes retrotide, a previously developed platform for retro-biosynthesis of PKSs to design. It is downloaded from github using the requirements.txt pip install required when starting up your virtual environment. The output of this is a series of 5 ranked PKS designs, with similarity scores of the products compared to the input molecule.

Finally, match_design_to_parts was added to help bridge betweek the PKS designs outputted by retrotide and tridentsynth, and the ClusterCAD query platforms. This is the natural next step for the PKS scientist -- not that we have a theoretical, de-novo design of a PKS, we can use CluserCAD to look at natural PKS designs and see which ones we can combine together to make a novel PKS to test in the lab.

```
natural language molecule
        │
        ▼
  resolve_smiles   →  SMILES format of molecule
        │
        ▼
  assess_pks_feasibility   →  rating from 0 to 1 on how PKS-able the inputted molecule is
        │
        ▼
  retrotide_designer     →  de-novo design of a theoretical PKS
         │
        ▼
  find_pks_module_parts   → searches ClusterCAD for real biological parts matching each module
```
---
## Novelty Aspect

While Retrotide itself was a previously developed PKS retrobiosynthesis platform, the software itself has been clunky for scientists to use, and often only as been able to be used for idea generation. The novelty in this section of the project primarily lies in the other tools, and how they make retrotide more useful for scientists in real life. Assess PKS feasibiltiy is entirely novel, and it helps scientists who may be new to the field evaluate if their desired molecule could realistically be produced by a PKS before they go down the rabithole of trying to design one. Additionally, the output of retrotide is formatted to be more useful to the user than the original software, and more useable. Finally, the match/find PKS modules is also novel in design, and helps easily bridge the users design process by allowing them to directly design a PKS they can test in lab by connecting retrotide outputs to tangiable, translatable laboratory designs, such as amino acid sequences through ClusterCad.

---

## Tool Summaries

### `resolve_smiles`

Molecules can be referred to in many ways — by a common name like "erythromycin," an IUPAC name, or a SMILES string (a text-based encoding of chemical structure). All of the downstream tools in this pipeline require SMILES as input. `resolve_smiles` handles this conversion automatically: if the user types a chemical name, it queries PubChem (a public database of over 110 million compounds) to look up the correct SMILES. If the input is already valid SMILES, it simply verifies and standardizes it. This prevents the AI from guessing or hallucinating a molecular structure.

**Input:** A chemical name or SMILES string (e.g., `"erythromycin"` or `"CCCC(=O)O"`)

**Output:** The canonical SMILES, molecular formula, and IUPAC name.

---

### `assess_pks_feasibility`

Not every molecule can be made by a type-I modular polyketide synthase (PKS). PKS assembly lines build molecules from two-carbon units and are best suited for linear or macrocyclic structures rich in carbon and oxygen. Molecules that are too small, too aromatic, or contain atoms like nitrogen or halogens are poor candidates. `assess_pks_feasibility` runs seven quick diagnostic checks — molecular weight, carbon chain length, heteroatom content, aromatic rings, oxygen-to-carbon ratio, functional group complexity, and ring count — and produces a score from 0 to 1. This score tells the AI which design tool to prioritize: high scores favor the pure-PKS designer (RetroTide), while low scores steer toward hybrid approaches (TridentSynth) that combine PKS with other enzymatic or chemical steps.

**Input:** A SMILES string.

**Output:** A feasibility score (0–1), pass/warn/fail status for each of the seven checks, and a recommendation on which design tool to try next.

---

### `retrotide_designer`

This tool performs retrobiosynthesis: given a target molecule, it works backwards to propose PKS assembly lines whose predicted product best matches the target. It uses the RetroTide algorithm, which enumerates possible combinations of PKS domains — choosing which extender unit each module loads (via the AT domain), whether each beta-keto group is reduced (KR), dehydrated (DH), or fully reduced (ER) — and then simulates each candidate assembly line to predict its product. The predicted products are ranked by chemical similarity to the target, and the top designs are returned with their full module-by-module domain architecture.

**Input:** A SMILES string, plus optional parameters for how many designs to return (default 5) and which similarity metric to use.

**Output:** A ranked list of PKS designs, each with a similarity score, predicted product SMILES, and the ordered list of modules showing every domain and its configuration.

---

### `find_pks_module_parts`

A RetroTide design is a blueprint — it tells you what each module should do (e.g., "load methylmalonyl-CoA, reduce with a B-type KR, dehydrate with DH") but not which real protein sequences to use. To actually build the PKS in the lab, each module in the design needs to be matched to a real amino acid sequence from a characterized natural PKS. `find_pks_module_parts` searches ClusterCAD, a database of 531 experimentally characterized PKS gene clusters, for natural modules that match a given domain profile. It accepts simple inputs — whether the module is a starter module, what substrate the AT domain loads, and which reductive domains are active — and returns matching natural modules with the amino acid sequence of every domain. These sequences are the starting point for gene synthesis.

**Input:** Whether the module is a loading module (true/false), the AT substrate name (e.g., "Malonyl-CoA"), and which reductive domains are active (e.g., "KR,DH,ER"). Called once per module in the design.

**Output:** A list of matching natural PKS modules from ClusterCAD, each with cluster accession, subunit name, and amino acid sequences for all domains.
