"""
HTTP/SSE transport server for connecting to Claude.ai → Settings → Integrations.

Usage:
    python server_http.py           # starts on http://0.0.0.0:8000
    python server_http.py --port 9000

Then expose it publicly (required by Claude.ai):
    ngrok http 8000
    # Use the https://.... URL ngrok prints as your integration URL in Claude.ai.

The MCP endpoint path is:  /mcp
Full URL example:  https://abc123.ngrok-free.app/mcp
"""

from __future__ import annotations

import sys
import argparse

from fastmcp import FastMCP
from modules import register_all

mcp = FastMCP(
    "BioE234 MCP Starter",
    instructions="""
You are a polyketide synthase (PKS) and natural products biosynthesis assistant.
You have access to a suite of tools for searching, designing, and analyzing PKS pathways.

## Tool chaining rules

1. **If the user gives a chemical name** (e.g. "erythromycin", "rapamycin") instead of a SMILES string:
   - ALWAYS call `resolve_smiles` first to convert the name to a canonical SMILES.
   - Then use that SMILES in any subsequent tool calls.
   - Never ask the user to provide a SMILES manually if they gave a name.

2. **If the user wants to find similar polyketides or pathways** (e.g. "find PKS for X", "what is similar to X", "search for X"):
   - Call `resolve_smiles` if input is a name, then call `search_pks` with the resulting SMILES.
   - Use search_type="reaction_search" by default; use "pathway_search" if the user asks for pathway steps.

3. **If the user wants to design a PKS** (e.g. "design a PKS to make X", "what modules do I need for X"):
   - Call `resolve_smiles` if needed, then call `pks_design_retrotide`.

4. **If the user asks about ClusterCAD clusters or domains**:
   - Use `clustercad_list_clusters`, `clustercad_cluster_details`, `clustercad_get_subunits`, or `clustercad_domain_lookup` as appropriate.

5. **If the user provides a DNA sequence and asks about GC content**:
   - Call `dna_gc_content`.

6. **If the user asks for a reverse complement or translation**:
   - Call `dna_reverse_complement` or `dna_translate`.

## General behavior
- Always call tools rather than answering from general knowledge when a tool is available.
- Show the user the key results in plain English after each tool call — don't just dump raw JSON.
- If a search returns engineering hints for intermediate matches, highlight them clearly.
- If resolve_smiles fails (unknown compound), tell the user and ask them to provide a SMILES directly.
"""
)

print("[server_http] Starting BioE234 MCP server (HTTP transport)...", file=sys.stderr)

register_all(mcp)

print("[server_http] All modules registered. Server ready.", file=sys.stderr)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--port", type=int, default=8000)
    parser.add_argument("--host", default="0.0.0.0")
    args = parser.parse_args()

    print(f"[server_http] Listening on http://{args.host}:{args.port}/mcp", file=sys.stderr)
    mcp.run(transport="http", host=args.host, port=args.port)
