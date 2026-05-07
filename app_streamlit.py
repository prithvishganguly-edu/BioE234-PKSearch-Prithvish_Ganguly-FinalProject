from __future__ import annotations

import asyncio
import json
import re
import time
from pathlib import Path
from typing import Any, Dict

import streamlit as st
from dotenv import load_dotenv
from fastmcp import Client
from google import genai
from google.genai import errors, types


# ----------------------------
# Helper functions
# ----------------------------

def _load_skill_context(modules_dir: Path) -> str:
    skill_texts = []

    if not modules_dir.exists():
        return ""

    for module_dir in sorted(modules_dir.iterdir()):
        if not module_dir.is_dir() or module_dir.name.startswith("_"):
            continue

        skill_file = module_dir / "SKILL.md"
        if skill_file.exists():
            skill_texts.append(skill_file.read_text())

    return "\n\n---\n\n".join(skill_texts)


def _strip_ctx_from_schema(schema: dict) -> dict:
    schema = dict(schema or {})
    props = dict(schema.get("properties", {}))
    props.pop("ctx", None)
    schema["properties"] = props

    if "required" in schema:
        schema["required"] = [r for r in schema["required"] if r != "ctx"]

    return schema


def _build_system_content(mcp_tools, mcp_resources, skill_context: str = "") -> types.Content:
    tools_json = [
        {
            "name": t.name,
            "description": getattr(t, "description", ""),
            "input_schema": getattr(t, "inputSchema", None),
        }
        for t in mcp_tools
    ]

    resources_json = []
    for r in mcp_resources:
        uri_obj = getattr(r, "uri", None) or getattr(r, "name", None)
        uri_str = str(uri_obj) if uri_obj is not None else None
        resource_name = uri_str.split("/")[-1] if uri_str else None

        resources_json.append(
            {
                "uri": uri_str,
                "name": resource_name,
                "description": getattr(r, "description", ""),
            }
        )

    payload = {
        "mcp_tools": tools_json,
        "mcp_resources": resources_json,
        "instruction": (
            "You are a GUI assistant connected to MCP tools. "
            "Use MCP tools when they are relevant. "
            "After each tool call, summarize the result clearly in plain English. "
            "For TridentSynth results, if result.text_summary is present, display that text_summary directly. "
            "Do not reinterpret target SMILES, post-PKS product SMILES, pathway structures, selected steps, "
            "or reaction details from other fields. Do not add duplicate sections or repeat the same information."
            "PKS modules, best reaction SMILES, pathway structure SMILES, reaction rule, "
            "reaction enthalpy, similarity scores, and selected steps."
        ),
    }

    system_text = "SYSTEM CONTEXT:\n" + json.dumps(payload, indent=2)

    if skill_context:
        system_text += "\n\n--- SKILL GUIDANCE ---\n\n" + skill_context

    return types.Content(
        role="model",
        parts=[types.Part.from_text(text=system_text)],
    )


def _mcp_tool_to_fn_declaration(tool: Any) -> types.FunctionDeclaration:
    params: Dict[str, Any] = getattr(tool, "inputSchema", None) or {
        "type": "object",
        "properties": {},
    }

    params = _strip_ctx_from_schema(params)

    desc = (getattr(tool, "description", None) or "").strip()
    if not desc:
        desc = f"MCP tool: {tool.name}"

    return types.FunctionDeclaration(
        name=tool.name,
        description=desc,
        parameters_json_schema=params,
    )


def _messages_to_contents(messages: list[dict[str, str]]) -> list[types.Content]:
    contents: list[types.Content] = []

    for message in messages:
        role = "user" if message["role"] == "user" else "model"
        contents.append(
            types.Content(
                role=role,
                parts=[types.Part.from_text(text=message["content"])],
            )
        )

    return contents


def _safe_generate(
    gemini: genai.Client,
    *,
    model: str,
    contents: list[types.Content],
    config: types.GenerateContentConfig,
    retries: int = 3,
    backoff_seconds: int = 2,
):
    for attempt in range(retries):
        try:
            return gemini.models.generate_content(
                model=model,
                contents=contents,
                config=config,
            )

        except errors.ServerError as e:
            msg = str(e)

            if ("503" in msg or "UNAVAILABLE" in msg) and attempt < retries - 1:
                wait = backoff_seconds * (2**attempt)
                time.sleep(wait)
                continue

            raise

        except errors.ClientError as e:
            msg = str(e)

            if ("429" in msg or "RESOURCE_EXHAUSTED" in msg) and attempt < retries - 1:
                match = re.search(r"retry[^\d]*(\d+(?:\.\d+)?)\s*s", msg, re.IGNORECASE)
                wait = float(match.group(1)) + 1 if match else backoff_seconds * (2**attempt)
                wait = max(wait, 2)
                time.sleep(wait)
                continue

            raise


async def _run_tool_loop(
    *,
    mcp: Client,
    initial_resp,
    contents: list[types.Content],
    gemini: genai.Client,
    model: str,
    config: types.GenerateContentConfig,
) -> str:
    resp = initial_resp
    contents = list(contents)

    while True:
        function_calls = resp.function_calls or []

        if not function_calls:
            return resp.text or "[No text response]"

        fc_content = resp.candidates[0].content
        fr_parts: list[types.Part] = []

        for fc in function_calls:
            tool_name = fc.name
            tool_args = dict(fc.args or {})

            try:
                tool_result = await mcp.call_tool(tool_name, tool_args)

                if isinstance(tool_result, list):
                    result_data = "\n".join(
                        getattr(item, "text", str(item)) for item in tool_result
                    )
                elif hasattr(tool_result, "content"):
                    result_data = "\n".join(
                        getattr(item, "text", str(item)) for item in tool_result.content
                    )
                else:
                    result_data = str(tool_result)

                fn_response = {"result": result_data}

            except Exception as e:
                fn_response = {"error": str(e)}

            fr_parts.append(
                types.Part.from_function_response(
                    name=tool_name,
                    response=fn_response,
                )
            )

        fr_content = types.Content(role="user", parts=fr_parts)
        contents.extend([fc_content, fr_content])

        resp = _safe_generate(
            gemini,
            model=model,
            contents=contents,
            config=config,
        )


# ----------------------------
# Persistent MCP + Gemini object
# ----------------------------

class PersistentGeminiMCP:
    def __init__(self):
        load_dotenv()

        self.project_root = Path(__file__).parent
        self.model = "gemini-2.5-flash"
        self.gemini = genai.Client()

        self.loop = asyncio.new_event_loop()
        asyncio.set_event_loop(self.loop)

        server_path = str(self.project_root / "server.py")
        self.mcp = Client(server_path)

        self.loop.run_until_complete(self._initialize())

    async def _initialize(self):
        await self.mcp.__aenter__()

        self.mcp_tools = await self.mcp.list_tools()
        self.mcp_resources = await self.mcp.list_resources()

        fn_decls = [_mcp_tool_to_fn_declaration(t) for t in self.mcp_tools]

        if fn_decls:
            tool_obj = types.Tool(function_declarations=fn_decls)
            self.config = types.GenerateContentConfig(tools=[tool_obj])
        else:
            self.config = types.GenerateContentConfig()

        modules_dir = self.project_root / "modules"
        skill_context = _load_skill_context(modules_dir)

        self.system_content = _build_system_content(
            mcp_tools=self.mcp_tools,
            mcp_resources=self.mcp_resources,
            skill_context=skill_context,
        )

    def ask(self, user_text: str, chat_history: list[dict[str, str]]) -> str:
        return self.loop.run_until_complete(
            self._ask_async(
                user_text=user_text,
                chat_history=chat_history,
            )
        )

    async def _ask_async(self, user_text: str, chat_history: list[dict[str, str]]) -> str:
        history_contents = _messages_to_contents(chat_history)

        user_content = types.Content(
            role="user",
            parts=[types.Part.from_text(text=user_text)],
        )

        contents = [self.system_content, *history_contents, user_content]

        initial_resp = _safe_generate(
            self.gemini,
            model=self.model,
            contents=contents,
            config=self.config,
        )

        final_text = await _run_tool_loop(
            mcp=self.mcp,
            initial_resp=initial_resp,
            contents=contents,
            gemini=self.gemini,
            model=self.model,
            config=self.config,
        )

        return final_text

    def close(self):
        try:
            self.loop.run_until_complete(self.mcp.__aexit__(None, None, None))
        except Exception:
            pass


@st.cache_resource
def get_assistant() -> PersistentGeminiMCP:
    return PersistentGeminiMCP()


# ----------------------------
# Streamlit UI
# ----------------------------

st.set_page_config(
    page_title="PKS Gemini MCP Assistant",
    page_icon="🧬",
    layout="wide",
)

st.title("🧬 PKS Gemini MCP Assistant")

with st.sidebar:
    st.header("Controls")

    if st.button("Restart MCP connection"):
        try:
            assistant = get_assistant()
            assistant.close()
        except Exception:
            pass

        st.cache_resource.clear()
        st.session_state.messages = []
        st.rerun()

    if st.button("Clear chat"):
        st.session_state.messages = []
        st.rerun()

    st.divider()

    try:
        assistant = get_assistant()
        tool_names = [t.name for t in assistant.mcp_tools]

        st.header("Registered tools")
        for name in tool_names:
            st.write(f"- {name}")

    except Exception as e:
        st.error(f"Could not load MCP tools: {e}")


if "messages" not in st.session_state:
    st.session_state.messages = []


for message in st.session_state.messages:
    with st.chat_message(message["role"]):
        st.markdown(message["content"])


prompt = st.chat_input("Ask Gemini to run your MCP tools...")

if prompt:
    st.session_state.messages.append(
        {
            "role": "user",
            "content": prompt,
        }
    )

    with st.chat_message("user"):
        st.markdown(prompt)

    with st.chat_message("assistant"):
        with st.spinner("Gemini is thinking and may call MCP tools..."):
            try:
                assistant = get_assistant()

                response = assistant.ask(
                    user_text=prompt,
                    chat_history=st.session_state.messages[:-1],
                )

            except Exception as e:
                response = f"Error: {e}"

        st.markdown(response)

    st.session_state.messages.append(
        {
            "role": "assistant",
            "content": response,
        }
    )