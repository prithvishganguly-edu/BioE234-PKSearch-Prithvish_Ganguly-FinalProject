# PKS Module — tridentsynth and Streamlit GUI

**Author:** Shalen Ardeshna  
**Course:** BioE 134, Spring 2026  

---

## What this contribution does

This contribution adds two connected pieces to the PKS design MCP project:

1. **`tridentsynth`** — a live MCP tool that submits target molecules to the TridentSynth web server and parses the best PKS-based synthesis pathway.
2. **`app_streamlit.py`** — a browser-based GUI for interacting with Gemini and the MCP tools without using the terminal-only client.

The goal is to let Gemini move from a natural-language request to a real TridentSynth pathway search, then return useful pathway information in a clean text summary.

---

## Files

| File | Description |
|------|-------------|
| `modules/pks/tools/tridentsynth.py` | Python implementation for live TridentSynth submission, polling, parsing, and result cleanup |
| `modules/pks/tools/tridentsynth.json` | C9 JSON wrapper for MCP tool registration |
| `modules/pks/tools/prompts.json` | Prompt examples for Gemini tool-calling evaluation |
| `tests/test_tridentsynth.py` | Pytest coverage for payload construction, parsing, validation, and cleanup |
| `app_streamlit.py` | Streamlit GUI for browser-based Gemini + MCP interaction |
| `modules/pks/tools/THEORY.md` | Background and design explanation |

---

## Tool 1 — `tridentsynth`

`tridentsynth` answers the question:

> Given a target molecule, what PKS-based pathway can TridentSynth propose to make it?

The tool accepts a target SMILES string and optional synthesis settings, submits the job to TridentSynth, waits for the result page, and parses the best pathway.

### Main parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `target_smiles` | string | required | Single target-molecule SMILES string |
| `use_pks` | bool | `true` | Include PKS assembly |
| `use_bio` | bool | `false` | Include biological tailoring |
| `use_chem` | bool | `false` | Include chemical tailoring |
| `max_bio_steps` | int | `1` when auto-filled | Maximum biological tailoring steps |
| `max_chem_steps` | int | `1` when auto-filled | Maximum chemical tailoring steps |
| `pks_release_mechanism` | string | inferred | `thiolysis` or `cyclization` |
| `pks_starters` | array | auto-filled | Optional PKS starter substrates |
| `pks_extenders` | array | auto-filled | Optional PKS extender substrates |
| `wait_for_completion` | bool | `true` | Polls until results are complete |
| `timeout_seconds` | int | `600` | Maximum wait time |
| `poll_seconds` | int | `5` | Time between polling attempts |

### Main outputs

| Field | Description |
|-------|-------------|
| `status` | `submitted`, `completed`, or `timed_out` |
| `task_id` | TridentSynth job/task ID |
| `submitted_query` | Normalized query submitted to TridentSynth |
| `result.text_summary` | Clean text summary for Gemini to display directly |
| `result.synthesis_parameters` | Target, strategy, PKS release, starters/extenders, and step settings |
| `result.selected_steps` | Bio/Chem step counts actually selected |
| `result.pks_modules` | Parsed PKS module architecture |
| `result.best_pathway` | PKS product, final product, reaction SMILES, reaction rules, enthalpy, feasibility, and pathway structures |

---

## Example usage

```text
User: Use TridentSynth to find a PKS plus chemical synthesis pathway for hexane. The target SMILES is CCCCCC.
```

Example tool call:

```text
tridentsynth(
    target_smiles="CCCCCC",
    use_pks=true,
    use_bio=false,
    use_chem=true,
    max_chem_steps=1,
    wait_for_completion=true,
    timeout_seconds=300,
    poll_seconds=5
)
```

Example output summary:

```text
TridentSynth result
- Target SMILES: CCCCCC
- Strategy: PKS, Chemical synthesis
- PKS termination: thiolysis
- Selected steps: Chemical steps: 1

PKS modules:
- Module 1 (Loading module): KSq, AT substrate Methylmalonyl-CoA, ACP
- Module 2 (Extension module): KS, AT substrate Malonyl-CoA, KR, DH, ER, ACP
- Module 3 (Extension module): KS, AT substrate Malonyl-CoA, KR, DH, ER, ACP

Best pathway:
- PKS product SMILES: CCCCCCC(=O)O
- Final product SMILES: CCCCCC
- Final product similarity to target: 1.0
- Reaction SMILES: CCCCCCC(=O)O>>CCCCCC.O=C=O
- Reaction details: rule = Carboxylic Acids Decarboxylation, enthalpy = -6.19 kcal/mol
- Pathway structures SMILES: CCCCCCC(=O)O, CCCCCC, O=C=O
```

---

## Tool 2 — Streamlit GUI

The Streamlit GUI gives the project a browser chat interface instead of requiring users to run the terminal-only Gemini client.

The GUI:

1. Connects to Gemini.
2. Starts and maintains an MCP connection.
3. Lists registered MCP tools.
4. Lets the user type prompts in a browser.
5. Lets Gemini call MCP tools automatically.
6. Displays the final response in the browser.

The GUI uses `st.cache_resource` to keep the MCP connection alive across messages. Use the **Restart MCP connection** button after changing tool files or if the MCP connection gets stuck.

---

## Running the GUI

From the repo root:

```bash
source .venv/bin/activate
streamlit run app_streamlit.py
```

The app usually opens at:

```text
http://localhost:8501
```

---

## Running tests

Run the TridentSynth tests with:

```bash
python -m pytest tests/test_tridentsynth.py
```

Expected normal result:

```text
11 passed, 1 skipped
```

The skipped test is the optional live TridentSynth submission test. To run it:

```bash
RUN_LIVE_TRIDENTSYNTH=1 python -m pytest tests/test_tridentsynth.py
```

---

## Test coverage

The pytest file covers:

| Area | Description |
|------|-------------|
| Payload construction | Confirms correct TridentSynth form fields |
| Strategy selection | Checks PKS, Bio, and Chem options |
| Input validation | Rejects invalid SMILES and invalid step counts |
| PKS options | Tests release mechanism, starters, and extenders |
| Reaction parsing | Splits reaction SMILES into reactants and products |
| PKS module parsing | Extracts module/domain architecture |
| Result parsing | Extracts best pathway fields from result HTML |
| Target correction | Prevents target SMILES from being confused with intermediates |
| Selected steps | Prevents unselected Bio/Chem steps from appearing |
| Optional live test | Allows real website submission only when explicitly enabled |

---

## Notes for users

- `tridentsynth` depends on the live TridentSynth website being available.
- Large pathway searches may take longer than simple examples.
- The returned pathway is computationally predicted and should be experimentally validated before use.
- The tool is intended to return TridentSynth pathway predictions, not to replace manual biological review.

---

## Summary

This contribution connects the MCP project to a real external PKS pathway design platform and adds a browser-based interface for using it. The `tridentsynth` tool runs live pathway searches and returns cleaned pathway results, while the Streamlit GUI makes the system easier to demo and interact with.
