from __future__ import annotations

import re
import time
from dataclasses import dataclass
from typing import Any, Optional
from urllib.parse import urljoin

import requests
from bs4 import BeautifulSoup, Tag


@dataclass
class _FieldCandidate:
    tag: Tag
    name: str
    field_type: str
    context: str
    value: str | None = None


class TridentSynth:
    """
    Description:
        Submit a live TridentSynth job through the public TridentSynth web form,
        then optionally poll the results page until completion.

        This tool discovers the current live HTML form fields dynamically instead
        of hard-coding field names, so it is more resilient to small frontend changes.

    Input:
        target_smiles (str): Single target-molecule SMILES string.
        use_pks (bool, optional): Include PKS synthesis strategy. Default: True.
        use_bio (bool, optional): Include Bio synthesis strategy. Default: False.
        use_chem (bool, optional): Include Chem synthesis strategy. Default: False.
        max_bio_steps (int, optional): Max biological steps. Allowed: 1-3.
        max_chem_steps (int, optional): Max chemical steps. Allowed: 1-3.
        pks_release_mechanism (str, optional): "thiolysis" or "cyclization".
        pks_starters (list[str], optional): PKS starter substrates.
        pks_extenders (list[str], optional): PKS extender substrates.
        max_carbon (int, optional): Optional DORAnet max carbon filter.
        max_nitrogen (int, optional): Optional DORAnet max nitrogen filter.
        max_oxygen (int, optional): Optional DORAnet max oxygen filter.
        auto_optimize_unspecified (bool, optional): Auto-fill omitted optional
            settings with a compact exploratory configuration. Default: True.
        wait_for_completion (bool, optional): Poll the results page until the
            job finishes. Default: True.
        poll_seconds (int, optional): Seconds between polling attempts. Default: 10.
        timeout_seconds (int, optional): Maximum total polling time. Default: 300.
        acknowledge_controlled_substance_warning (bool, optional): If the site
            shows its controlled-substance acknowledgement screen, continue only
            when this is True. Default: False.

    Output:
        dict: JSON-serializable summary containing the submitted query, parsed job
        info, results URL, and parsed results if available.

    Tests:
        - Case:
            Input: target_smiles="CCCCCCC", use_pks=True, use_bio=True
            Expected Output: result["ok"] == True
            Description: Valid PKS+Bio live submission setup.
        - Case:
            Input: target_smiles="CC=C(C)C(=O)O", use_pks=True, use_bio=False, use_chem=False
            Expected Output: result["normalized_query"]["use_pks"] == True
            Description: PKS-only live submission setup.
        - Case:
            Input: target_smiles="CC.CC"
            Expected Exception: ValueError
            Description: TridentSynth accepts one molecule at a time.
        - Case:
            Input: target_smiles="CCCC", use_pks=False, use_bio=False, use_chem=False
            Expected Exception: ValueError
            Description: At least one synthesis strategy must be selected.
        - Case:
            Input: target_smiles="CCCC", max_bio_steps=4
            Expected Exception: ValueError
            Description: Bio/Chem step limits must be in the 1-3 range.
    """

    def initiate(self) -> None:
        self.base_url = "https://tridentsynth.lbl.gov/run_TridentSynth/"
        self.default_headers = {
            "User-Agent": "BioE234-TridentSynth-MCP/1.0 (+research use)",
        }
        self.allowed_release_mechanisms = {"thiolysis", "cyclization"}
        self.allowed_starters = {
            "malonyl-coa": "Malonyl-CoA",
            "methylmalonyl-coa": "Methylmalonyl-CoA",
            "allylmalonyl-coa": "Allylmalonyl-CoA",
            "methoxymalonyl-coa": "Methoxymalonyl-CoA",
            "hydroxymalonyl-coa": "Hydroxymalonyl-CoA",
        }
        self.allowed_extenders = {
            "malonyl-coa": "Malonyl-CoA",
            "methylmalonyl-coa": "Methylmalonyl-CoA",
            "ethylmalonyl-coa": "Ethylmalonyl-CoA",
            "methoxymalonyl-coa": "Methoxymalonyl-CoA",
            "allylmalonyl-coa": "Allylmalonyl-CoA",
        }

    def _session(self) -> requests.Session:
        session = requests.Session()
        session.headers.update(self.default_headers)
        return session

    def _clean(self, text: str) -> str:
        return re.sub(r"\s+", " ", text).strip().lower()

    def _normalize_choice_list(
        self,
        values: Optional[list[str]],
        allowed_map: dict[str, str],
        field_name: str,
    ) -> list[str]:
        if not values:
            return []

        out: list[str] = []
        seen: set[str] = set()
        for value in values:
            key = self._clean(str(value))
            if key not in allowed_map:
                valid = ", ".join(sorted(allowed_map.values()))
                raise ValueError(f"Invalid {field_name}: {value!r}. Allowed values: {valid}")
            canonical = allowed_map[key]
            if canonical not in seen:
                seen.add(canonical)
                out.append(canonical)
        return out

    def _validate_step_count(self, value: Optional[int], field_name: str) -> Optional[int]:
        if value is None:
            return None
        if value not in (1, 2, 3):
            raise ValueError(f"{field_name} must be 1, 2, or 3 — got {value}")
        return value

    def _validate_atom_limit(self, value: Optional[int], field_name: str) -> Optional[int]:
        if value is None:
            return None
        if value < 0:
            raise ValueError(f"{field_name} must be >= 0")
        return value

    def _infer_release_mechanism(self, target_smiles: str) -> str:
        has_ring_digit = any(ch.isdigit() for ch in target_smiles)
        has_o_or_n = any(ch in target_smiles for ch in ("O", "N", "o", "n"))
        return "cyclization" if has_ring_digit and has_o_or_n else "thiolysis"

    def _get_soup(self, session: requests.Session, url: str) -> BeautifulSoup:
        response = session.get(url, timeout=60)
        response.raise_for_status()
        return BeautifulSoup(response.text, "html.parser")

    def _choose_submission_form(self, soup: BeautifulSoup) -> Tag:
        forms = soup.find_all("form")
        if not forms:
            raise RuntimeError("No HTML form was found on the TridentSynth page.")

        def score(form: Tag) -> int:
            text = self._clean(form.get_text(" ", strip=True))
            s = 0
            if "run" in text:
                s += 5
            if "target molecule" in text or "smiles" in text:
                s += 8
            if "pks" in text:
                s += 2
            if form.find(["input", "select", "textarea"]):
                s += 2
            return s

        return max(forms, key=score)

    def _context_text(self, tag: Tag) -> str:
        parts: list[str] = []
        attrs = [
            tag.get("name", ""),
            tag.get("id", ""),
            tag.get("placeholder", ""),
            tag.get("aria-label", ""),
            tag.get("title", ""),
            tag.get("value", "") if tag.name == "input" else "",
        ]
        parts.extend(str(x) for x in attrs if x)

        tag_id = tag.get("id")
        if tag_id:
            parent = tag.find_parent()
            if parent is not None:
                label = parent.find("label", attrs={"for": tag_id})
                if label:
                    parts.append(label.get_text(" ", strip=True))

        parent = tag.find_parent(["label", "div", "td", "th", "p", "section"])
        if parent:
            parts.append(parent.get_text(" ", strip=True))

        previous = []
        for sib in list(tag.previous_strings)[-4:]:
            txt = str(sib).strip()
            if txt:
                previous.append(txt)
        parts.extend(previous)

        return self._clean(" ".join(parts))

    def _gather_fields(self, form: Tag) -> list[_FieldCandidate]:
        fields: list[_FieldCandidate] = []
        for tag in form.find_all(["input", "select", "textarea"]):
            name = tag.get("name")
            if not name:
                continue
            field_type = tag.get("type", tag.name).lower()
            fields.append(
                _FieldCandidate(
                    tag=tag,
                    name=name,
                    field_type=field_type,
                    context=self._context_text(tag),
                    value=tag.get("value"),
                )
            )
        return fields

    def _find_best(
        self,
        fields: list[_FieldCandidate],
        must_have: list[str],
        optional: Optional[list[str]] = None,
        field_types: Optional[list[str]] = None,
    ) -> Optional[_FieldCandidate]:
        optional = optional or []
        best: Optional[_FieldCandidate] = None
        best_score = -10**9

        for field in fields:
            if field_types and field.field_type not in field_types:
                continue
            text = field.context
            if not all(token in text for token in must_have):
                continue
            score = 0
            for token in must_have:
                if token in text:
                    score += 8
            for token in optional:
                if token in text:
                    score += 2
            if field.field_type == "hidden":
                score -= 10
            if score > best_score:
                best = field
                best_score = score
        return best

    def _find_checkbox_by_value_or_context(
        self,
        fields: list[_FieldCandidate],
        value_keywords: list[str],
        group_keywords: Optional[list[str]] = None,
    ) -> Optional[_FieldCandidate]:
        group_keywords = group_keywords or []
        best: Optional[_FieldCandidate] = None
        best_score = -10**9

        for field in fields:
            if field.field_type not in {"checkbox", "radio"}:
                continue
            hay = self._clean(f"{field.context} {field.value or ''}")
            if not all(tok in hay for tok in value_keywords):
                continue
            score = 0
            for tok in value_keywords:
                if tok in hay:
                    score += 8
            for tok in group_keywords:
                if tok in hay:
                    score += 2
            if score > best_score:
                best = field
                best_score = score
        return best

    def _default_payload_items(self, form: Tag) -> list[tuple[str, str]]:
        items: list[tuple[str, str]] = []
        for tag in form.find_all(["input", "select", "textarea"]):
            name = tag.get("name")
            if not name:
                continue

            if tag.name == "input":
                field_type = tag.get("type", "text").lower()
                if field_type == "hidden":
                    items.append((name, tag.get("value", "")))
                elif field_type in {"checkbox", "radio"} and tag.has_attr("checked"):
                    items.append((name, tag.get("value", "on")))
                elif field_type not in {"submit", "button", "file", "image"} and tag.get("value"):
                    items.append((name, tag.get("value", "")))

            elif tag.name == "textarea":
                if tag.text:
                    items.append((name, tag.text))

            elif tag.name == "select":
                selected = tag.find("option", selected=True)
                if selected is not None and selected.get("value") is not None:
                    items.append((name, selected.get("value", "")))

        return items

    def _replace_single(self, items: list[tuple[str, str]], name: str, value: str) -> None:
        items[:] = [(k, v) for (k, v) in items if k != name]
        items.append((name, value))

    def _append_multi(self, items: list[tuple[str, str]], name: str, value: str) -> None:
        items.append((name, value))

    def _set_strategy_checkboxes(
        self,
        items: list[tuple[str, str]],
        fields: list[_FieldCandidate],
        use_pks: bool,
        use_bio: bool,
        use_chem: bool,
    ) -> dict[str, str]:
        chosen: dict[str, str] = {}
        for enabled, key in [(use_pks, "pks"), (use_bio, "bio"), (use_chem, "chem")]:
            if not enabled:
                continue
            field = self._find_checkbox_by_value_or_context(fields, [key], ["strategy", "synthesis"])
            if field is None:
                raise RuntimeError(f"Could not find the {key.upper()} strategy field on the live TridentSynth form.")
            self._append_multi(items, field.name, field.value or "on")
            chosen[key] = field.value or "on"
        return chosen

    def _set_text_or_select(
        self,
        items: list[tuple[str, str]],
        fields: list[_FieldCandidate],
        value: Optional[Any],
        must_have: list[str],
        optional: Optional[list[str]] = None,
        field_types: Optional[list[str]] = None,
        required: bool = False,
    ) -> Optional[str]:
        if value is None:
            return None
        field = self._find_best(fields, must_have, optional=optional, field_types=field_types)
        if field is None:
            if required:
                raise RuntimeError(f"Could not find form field for {' '.join(must_have)}")
            return None
        self._replace_single(items, field.name, str(value))
        return field.name

    def _set_checkbox_group(
        self,
        items: list[tuple[str, str]],
        fields: list[_FieldCandidate],
        selections: list[str],
        group_keyword: str,
    ) -> dict[str, list[str]]:
        chosen: dict[str, list[str]] = {}
        if not selections:
            return chosen

        matched_fields: list[_FieldCandidate] = []
        for selection in selections:
            normalized = self._clean(selection)
            keywords = [normalized]
            field = self._find_checkbox_by_value_or_context(fields, keywords, [group_keyword])
            if field is None:
                split_keywords = [tok for tok in re.split(r"[^a-z0-9]+", normalized) if tok]
                field = self._find_checkbox_by_value_or_context(fields, split_keywords, [group_keyword])
            if field is None:
                raise RuntimeError(f"Could not find checkbox for {group_keyword} choice {selection!r}")
            matched_fields.append(field)

        names = {field.name for field in matched_fields}
        items[:] = [(k, v) for (k, v) in items if k not in names]

        for field in matched_fields:
            self._append_multi(items, field.name, field.value or "on")
            chosen.setdefault(field.name, []).append(field.value or "on")

        return chosen

    def _submit_form(
        self,
        session: requests.Session,
        base_soup: BeautifulSoup,
        items: list[tuple[str, str]],
    ) -> requests.Response:
        form = self._choose_submission_form(base_soup)
        action = form.get("action") or self.base_url
        method = (form.get("method") or "post").lower()
        url = urljoin(self.base_url, action)

        if method == "get":
            response = session.get(url, params=items, timeout=120)
        else:
            response = session.post(url, data=items, timeout=120)

        response.raise_for_status()
        return response

    def _handle_controlled_substance_warning(
        self,
        session: requests.Session,
        response: requests.Response,
        acknowledge: bool,
    ) -> Optional[requests.Response]:
        soup = BeautifulSoup(response.text, "html.parser")
        text = soup.get_text(" ", strip=True)

        if "Controlled Substance Warning" not in text:
            return response

        if not acknowledge:
            return None

        forms = soup.find_all("form")
        for form in forms:
            form_text = self._clean(form.get_text(" ", strip=True))
            if "continue anyway" not in form_text:
                continue

            items = self._default_payload_items(form)

            buttons = form.find_all(["button", "input"])
            for button in buttons:
                button_text = self._clean(button.get_text(" ", strip=True) or button.get("value", ""))
                if "continue" in button_text and button.get("name"):
                    items.append((button["name"], button.get("value", "Continue Anyway")))
                    break

            action = urljoin(response.url, form.get("action") or response.url)
            method = (form.get("method") or "post").lower()

            if method == "get":
                next_response = session.get(action, params=items, timeout=120)
            else:
                next_response = session.post(action, data=items, timeout=120)

            next_response.raise_for_status()
            return next_response

        return response

    def _extract_job_info(self, soup: BeautifulSoup, current_url: str) -> dict[str, Any]:
        text = soup.get_text("\n", strip=True)
        lines = [line.strip() for line in text.splitlines() if line.strip()]

        job_id = None
        results_link = None
        status_hint = None

        for idx, line in enumerate(lines):
            low = line.lower()
            if low == "job id" and idx + 1 < len(lines):
                job_id = lines[idx + 1]
            elif low == "job status" and idx + 1 < len(lines):
                status_hint = lines[idx + 1]

        for a in soup.find_all("a", href=True):
            href = a["href"]
            full = urljoin(current_url, href)
            if "/run_TridentSynth/" not in full:
                continue
            if any(bad in full for bad in ["tutorial", "about", "publications", "common_precursors"]):
                continue
            if full.rstrip("/") != self.base_url.rstrip("/"):
                results_link = full
                break

        return {
            "job_id": job_id,
            "results_url": results_link,
            "submission_status_hint": status_hint,
            "raw_page_excerpt": text[:1500],
        }

    def _extract_summary_value(self, lines: list[str], label: str) -> Optional[str]:
        for idx, line in enumerate(lines):
            if line.strip().lower() == label.lower() and idx + 1 < len(lines):
                return lines[idx + 1]
        return None

    def _parse_results(self, soup: BeautifulSoup, url: str) -> dict[str, Any]:
        text = soup.get_text("\n", strip=True)
        lines = [line.strip() for line in text.splitlines() if line.strip()]

        status_message = None
        for phrase in [
            "Pathways to target found!",
            "The target can be synthesized by PKS assembly alone",
            "Closest reachable product",
            "closest reachable product",
            "Job submitted. Please wait.",
        ]:
            if phrase in text:
                status_message = phrase
                break

        summary = {
            "task_id": self._extract_summary_value(lines, "Task ID"),
            "target_smiles": self._extract_summary_value(lines, "Target SMILES"),
            "target_name": self._extract_summary_value(lines, "Target Name"),
            "pathway_sequence": self._extract_summary_value(lines, "Pathway Sequence"),
            "pks_termination_step": self._extract_summary_value(lines, "PKS Termination Step"),
            "pks_extenders": self._extract_summary_value(lines, "PKS Extenders"),
            "pks_starters": self._extract_summary_value(lines, "PKS Starters"),
            "bio_steps": self._extract_summary_value(lines, "# Bio Steps"),
            "chem_steps": self._extract_summary_value(lines, "# Chem Steps"),
            "job_id": self._extract_summary_value(lines, "Job Id"),
        }

        top_design: dict[str, Any] = {
            "pks_product": self._extract_summary_value(lines, "PKS product"),
            "post_pks_product": self._extract_summary_value(lines, "Post-PKS product"),
        }

        similarity_values = []
        for idx, line in enumerate(lines):
            if line.lower().startswith("similarity to target") and idx + 1 < len(lines):
                similarity_values.append(lines[idx + 1])

        if similarity_values:
            top_design["pks_similarity_to_target"] = similarity_values[0]
        if len(similarity_values) > 1:
            top_design["post_pks_similarity_to_target"] = similarity_values[1]

        net_feasibility = self._extract_summary_value(lines, "Net feasibility")
        if net_feasibility is not None:
            top_design["top_pathway_net_feasibility"] = net_feasibility

        rxn_smiles: list[str] = []
        collect = False
        for line in lines:
            if line == "Reactions (SMILES)":
                collect = True
                continue
            if collect:
                if line.startswith("Reaction rules") or line.startswith("Step feasibilities"):
                    break
                if ">>" in line:
                    rxn_smiles.append(line.strip("•* ").strip())

        if rxn_smiles:
            top_design["top_pathway_reaction_smiles"] = rxn_smiles

        return {
            "results_url": url,
            "status_message": status_message,
            "job_summary": summary,
            "top_design": top_design,
            "raw_text_excerpt": text[:4000],
        }

    def _is_complete(self, parsed_results: dict[str, Any]) -> bool:
        msg = (parsed_results.get("status_message") or "").lower()
        if not msg:
            return False
        if "please wait" in msg:
            return False
        return True

    def run(
        self,
        target_smiles: str,
        use_pks: bool = True,
        use_bio: bool = False,
        use_chem: bool = False,
        max_bio_steps: Optional[int] = None,
        max_chem_steps: Optional[int] = None,
        pks_release_mechanism: Optional[str] = None,
        pks_starters: Optional[list[str]] = None,
        pks_extenders: Optional[list[str]] = None,
        max_carbon: Optional[int] = None,
        max_nitrogen: Optional[int] = None,
        max_oxygen: Optional[int] = None,
        auto_optimize_unspecified: bool = True,
        wait_for_completion: bool = True,
        poll_seconds: int = 10,
        timeout_seconds: int = 300,
        acknowledge_controlled_substance_warning: bool = False,
    ) -> dict[str, Any]:
        target_smiles = target_smiles.strip()
        if not target_smiles:
            raise ValueError("target_smiles must not be empty")
        if "." in target_smiles:
            raise ValueError(
                "TridentSynth accepts one molecule at a time. Remove '.' and submit a single target molecule."
            )
        if not any((use_pks, use_bio, use_chem)):
            raise ValueError("At least one synthesis strategy must be selected: use_pks, use_bio, or use_chem")

        max_bio_steps = self._validate_step_count(max_bio_steps, "max_bio_steps")
        max_chem_steps = self._validate_step_count(max_chem_steps, "max_chem_steps")
        max_carbon = self._validate_atom_limit(max_carbon, "max_carbon")
        max_nitrogen = self._validate_atom_limit(max_nitrogen, "max_nitrogen")
        max_oxygen = self._validate_atom_limit(max_oxygen, "max_oxygen")

        starters = self._normalize_choice_list(pks_starters, self.allowed_starters, "pks_starters")
        extenders = self._normalize_choice_list(pks_extenders, self.allowed_extenders, "pks_extenders")

        auto_filled: dict[str, Any] = {}
        if auto_optimize_unspecified:
            if use_bio and max_bio_steps is None:
                max_bio_steps = 1
                auto_filled["max_bio_steps"] = 1
            if use_chem and max_chem_steps is None:
                max_chem_steps = 1
                auto_filled["max_chem_steps"] = 1
            if use_pks and pks_release_mechanism is None:
                pks_release_mechanism = self._infer_release_mechanism(target_smiles)
                auto_filled["pks_release_mechanism"] = pks_release_mechanism
            if use_pks and not starters:
                starters = ["Malonyl-CoA", "Methylmalonyl-CoA"]
                auto_filled["pks_starters"] = starters
            if use_pks and not extenders:
                extenders = ["Malonyl-CoA", "Methylmalonyl-CoA"]
                auto_filled["pks_extenders"] = extenders

        if pks_release_mechanism is not None:
            pks_release_mechanism = self._clean(pks_release_mechanism)
            if pks_release_mechanism not in self.allowed_release_mechanisms:
                allowed = ", ".join(sorted(self.allowed_release_mechanisms))
                raise ValueError(
                    f"Invalid pks_release_mechanism: {pks_release_mechanism!r}. Allowed values: {allowed}"
                )

        session = self._session()
        base_soup = self._get_soup(session, self.base_url)
        form = self._choose_submission_form(base_soup)
        fields = self._gather_fields(form)
        items = self._default_payload_items(form)

        target_field = (
            self._find_best(fields, ["target", "molecule"], ["smiles"], ["text", "search", "textarea"])
            or self._find_best(fields, ["smiles"], ["target"], ["text", "search", "textarea"])
            or self._find_best(fields, ["molecule"], ["target"], ["text", "search", "textarea"])
        )
        if target_field is None:
            raise RuntimeError("Could not find the target SMILES field on the live TridentSynth form.")
        self._replace_single(items, target_field.name, target_smiles)

        strategy_fields = self._set_strategy_checkboxes(items, fields, use_pks, use_bio, use_chem)

        self._set_text_or_select(
            items, fields, max_bio_steps, ["biological", "steps"], ["bio"], ["select", "number", "text"]
        )
        self._set_text_or_select(
            items, fields, max_chem_steps, ["chemical", "steps"], ["chem"], ["select", "number", "text"]
        )
        self._set_text_or_select(
            items, fields, pks_release_mechanism, ["release", "mechanism"], ["termination", "pks"], ["select", "radio", "text"]
        )
        starter_fields = self._set_checkbox_group(items, fields, starters, "starter")
        extender_fields = self._set_checkbox_group(items, fields, extenders, "extender")
        self._set_text_or_select(items, fields, max_carbon, ["carbon"], ["max", "atom"], ["select", "number", "text"])
        self._set_text_or_select(items, fields, max_nitrogen, ["nitrogen"], ["max", "atom"], ["select", "number", "text"])
        self._set_text_or_select(items, fields, max_oxygen, ["oxygen"], ["max", "atom"], ["select", "number", "text"])

        submit_response = self._submit_form(session, base_soup, items)
        submit_response = self._handle_controlled_substance_warning(
            session,
            submit_response,
            acknowledge=acknowledge_controlled_substance_warning,
        )

        if submit_response is None:
            return {
                "ok": False,
                "tool": "TridentSynth",
                "submission_url": self.base_url,
                "warning": (
                    "TridentSynth returned a controlled-substance acknowledgement screen. "
                    "Re-run with acknowledge_controlled_substance_warning=True only if you "
                    "have a legitimate research purpose."
                ),
                "normalized_query": {
                    "target_smiles": target_smiles,
                    "use_pks": use_pks,
                    "use_bio": use_bio,
                    "use_chem": use_chem,
                    "max_bio_steps": max_bio_steps,
                    "max_chem_steps": max_chem_steps,
                    "pks_release_mechanism": pks_release_mechanism,
                    "pks_starters": starters,
                    "pks_extenders": extenders,
                    "max_carbon": max_carbon,
                    "max_nitrogen": max_nitrogen,
                    "max_oxygen": max_oxygen,
                },
            }

        submit_soup = BeautifulSoup(submit_response.text, "html.parser")
        job_info = self._extract_job_info(submit_soup, submit_response.url)

        out: dict[str, Any] = {
            "ok": True,
            "tool": "TridentSynth",
            "submission_url": self.base_url,
            "normalized_query": {
                "target_smiles": target_smiles,
                "use_pks": use_pks,
                "use_bio": use_bio,
                "use_chem": use_chem,
                "max_bio_steps": max_bio_steps,
                "max_chem_steps": max_chem_steps,
                "pks_release_mechanism": pks_release_mechanism,
                "pks_starters": starters,
                "pks_extenders": extenders,
                "max_carbon": max_carbon,
                "max_nitrogen": max_nitrogen,
                "max_oxygen": max_oxygen,
            },
            "auto_filled": auto_filled,
            "live_form_fields_used": {
                "target_smiles": target_field.name,
                "strategies": strategy_fields,
                "pks_starters": starter_fields,
                "pks_extenders": extender_fields,
            },
            "job_id": job_info.get("job_id"),
            "results_url": job_info.get("results_url"),
            "submission_status_hint": job_info.get("submission_status_hint"),
            "submission_page_excerpt": job_info.get("raw_page_excerpt"),
        }

        results_url = job_info.get("results_url")
        if not wait_for_completion or not results_url:
            out["status"] = "submitted"
            if not results_url:
                out["warning"] = (
                    "The job appears to have been submitted, but no results link could be parsed "
                    "from the response page."
                )
            return out

        deadline = time.time() + max(timeout_seconds, 1)
        last_parsed: dict[str, Any] | None = None

        while time.time() <= deadline:
            result_response = session.get(results_url, timeout=120)
            result_response.raise_for_status()
            result_soup = BeautifulSoup(result_response.text, "html.parser")
            last_parsed = self._parse_results(result_soup, results_url)

            if self._is_complete(last_parsed):
                out["status"] = "completed"
                out["results"] = last_parsed
                return out

            time.sleep(max(poll_seconds, 1))

        out["status"] = "submitted_but_still_running"
        if last_parsed is not None:
            out["results"] = last_parsed
        out["warning"] = (
            f"Timed out after {timeout_seconds} seconds while polling the results page. "
            "Use the returned results_url to check later."
        )
        return out


_instance = TridentSynth()
_instance.initiate()
tridentsynth = _instance.run