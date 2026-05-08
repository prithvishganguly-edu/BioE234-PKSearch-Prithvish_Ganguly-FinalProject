# PKS De-Novo Design — resolve_smiles, assess_pks_feasibility, retrotide_designer, and match_design_to_parts

**Author:** Johanna Kann 
**Course:** BioE 234, Spring 2026  

---

## What these tools do

This set of three tools forms the **preliminary design** stage of the PKS design pipeline. First, users input a natural language query for a molecule -- ranging from heptane to erythromycin to nitrogen. The first step is to convert this to a SMILES format, which queries PubChem to 

 utilizes retrotide, a previously developed platform for retro-biosynthesis of PKSs to design a 
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
  match_design_to_parts   → provides framework to plug output of retrotide/tridentsynth into clustercad tools
```

---
