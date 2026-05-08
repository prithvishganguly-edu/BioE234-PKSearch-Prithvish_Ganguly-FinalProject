# TridentSynth Tool and GUI Contribution

## Overview

This contribution adds two main features to the PKS design MCP project:

1. A live `tridentsynth` MCP tool that connects Gemini to the TridentSynth web server.
2. A Streamlit GUI that lets users interact with the MCP + Gemini system through a browser instead of the terminal.

Together, these additions make the project more usable for PKS pathway design because users can submit target molecules, run pathway searches, and view summarized results in a graphical chat interface.

---

## TridentSynth MCP Tool

### Purpose

The `tridentsynth` tool allows Gemini to run live TridentSynth pathway searches from a target SMILES string. Instead of only reasoning from general chemistry knowledge, Gemini can submit a real query to TridentSynth and return the best predicted PKS-based synthesis pathway.

### What the Tool Does

The tool:

1. Takes a target molecule as a SMILES string.
2. Builds a TridentSynth-compatible form submission.
3. Selects synthesis strategies such as PKS, biological synthesis, and/or chemical synthesis.
4. Submits the job to the live TridentSynth web server.
5. Waits for the completed results page.
6. Parses the best pathway.
7. Returns the result as structured data and a clean text summary.

### Inputs

Important inputs include:

- `target_smiles`: Target molecule as a SMILES string.
- `use_pks`: Whether to use PKS assembly.
- `use_bio`: Whether to use biological tailoring steps.
- `use_chem`: Whether to use chemical tailoring steps.
- `max_bio_steps`: Maximum number of biological steps.
- `max_chem_steps`: Maximum number of chemical steps.
- `pks_release_mechanism`: PKS release method, such as `thiolysis` or `cyclization`.
- `pks_starters`: Optional PKS starter substrates.
- `pks_extenders`: Optional PKS extender substrates.
- `wait_for_completion`: Whether to wait for the full completed result.

### Outputs

The tool returns:

- job status
- TridentSynth task ID
- submitted query
- synthesis parameters
- selected steps
- PKS modules
- PKS product SMILES
- final post-PKS product SMILES
- similarity scores
- reaction SMILES
- reaction rule names
- reaction enthalpy
- pathway structures as SMILES
- a clean `text_summary` for Gemini to display directly

### Example Prompt

```text
Use the tridentsynth tool with target_smiles CCCCCC, use_pks true, use_bio false, use_chem true, max_chem_steps 1, wait_for_completion true, timeout_seconds 300, poll_seconds 5.
```

### Example Output

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

### Why This Matters

This tool makes the project more powerful because it connects the MCP system to an external pathway design platform. Gemini can now run a specialized TridentSynth search, extract the best pathway, and summarize relevant PKS design information in plain text.

---

## Streamlit GUI

### Purpose

The Streamlit GUI provides a browser-based interface for the MCP + Gemini system. Before this, users had to interact with Gemini through the terminal. The GUI makes the project easier to demo and easier for users to interact with.

### What the GUI Does

The GUI:

1. Starts a local Streamlit web app.
2. Connects to Gemini.
3. Connects to the MCP server.
4. Lists registered MCP tools.
5. Lets the user type prompts into a chat interface.
6. Allows Gemini to call tools automatically.
7. Displays the final response in the browser.

### Main File

The GUI is implemented in:

```text
app_streamlit.py
```

It should be placed in the repo root, at the same level as:

```text
client_gemini.py
server.py
requirements.txt
modules/
```

### Running the GUI

From the repo root:

```bash
cd ~/Documents/GitHub/2026-bioe234-final-project-pks
source .venv/bin/activate
streamlit run app_streamlit.py
```

This opens a local browser app, usually at:

```text
http://localhost:8501
```

### Persistent MCP Connection

The GUI uses `st.cache_resource` to keep a persistent MCP connection alive. This prevents the app from restarting the MCP server for every message, making the GUI faster and smoother.

The sidebar includes a restart button that can be used after changing tool files:

```text
Restart MCP connection
```

### Why This Matters

The GUI makes the project more accessible and presentation-ready. Instead of typing prompts in a terminal, users can interact with the system through a clean chat interface. This is especially useful for demonstrating TridentSynth because the user can run a pathway search and view the parsed results directly in the browser.

---

## Requirements Added

This contribution requires the following packages:

```text
streamlit
requests
beautifulsoup4
pytest
```

These should be included in `requirements.txt`.

Install all requirements with:

```bash
python -m pip install -r requirements.txt
```

---

## Testing

The TridentSynth tool includes a pytest file:

```text
tests/test_tridentsynth.py
```

The tests check:

- payload construction
- PKS/Bio/Chem strategy selection
- invalid SMILES validation
- reaction SMILES parsing
- PKS module parsing
- result-page parsing
- target SMILES correction
- selected step handling

Run the tests with:

```bash
python -m pytest tests/test_tridentsynth.py
```

Expected result:

```text
11 passed, 1 skipped
```

The skipped test is the optional live TridentSynth test. It only runs if:

```bash
RUN_LIVE_TRIDENTSYNTH=1
```

is set.

---

## Summary of Contribution

This contribution adds a live TridentSynth pathway design tool and a graphical interface for the MCP project.

The TridentSynth tool expands the system from static/local PKS utilities to a live pathway search workflow. The GUI improves usability by allowing users to interact with Gemini and the MCP tools through a browser. Together, these additions make the project easier to use, easier to demo, and more useful for PKS pathway design.
