from __future__ import annotations

import re
import time
from typing import Any, Optional
from urllib.parse import urljoin

import requests
from bs4 import BeautifulSoup


class TridentSynth:
    """
    Live TridentSynth runner.

    Submits a job to the public TridentSynth site, waits for completion, parses the
    best pathway, and returns pathway structures, PKS modules, reaction rules, and
    other result information as text/JSON.
    """

    def initiate(self) -> None:
        self.base_url = "https://tridentsynth.lbl.gov/run_TridentSynth/"
        self.headers = {
            "User-Agent": "BioE234-TridentSynth-MCP/1.0"
        }

        self.allowed_release_mechanisms = {"thiolysis", "cyclization"}

        self.starter_map = {
            "malonyl-coa": "mal",
            "malonyl coa": "mal",
            "mal": "mal",
            "methylmalonyl-coa": "mmal",
            "methylmalonyl coa": "mmal",
            "mmal": "mmal",
            "allylmalonyl-coa": "allylmal",
            "allylmalonyl coa": "allylmal",
            "allylmal": "allylmal",
            "methoxymalonyl-coa": "mxmal",
            "methoxymalonyl coa": "mxmal",
            "mxmal": "mxmal",
            "hydroxymalonyl-coa": "hmal",
            "hydroxymalonyl coa": "hmal",
            "hmal": "hmal",
        }

        self.extender_map = {
            "malonyl-coa": "mal",
            "malonyl coa": "mal",
            "mal": "mal",
            "methylmalonyl-coa": "mmal",
            "methylmalonyl coa": "mmal",
            "mmal": "mmal",
            "ethylmalonyl-coa": "emal",
            "ethylmalonyl coa": "emal",
            "emal": "emal",
            "methoxymalonyl-coa": "mxmal",
            "methoxymalonyl coa": "mxmal",
            "mxmal": "mxmal",
            "allylmalonyl-coa": "allylmal",
            "allylmalonyl coa": "allylmal",
            "allylmal": "allylmal",
        }

    def _session(self) -> requests.Session:
        session = requests.Session()
        session.headers.update(self.headers)
        return session

    def _clean(self, value: Any) -> str:
        return re.sub(r"\s+", " ", str(value)).strip().lower()

    def _validate_step_count(self, value: Optional[int], field_name: str) -> Optional[int]:
        if value is None:
            return None
        if value not in (1, 2, 3):
            raise ValueError(f"{field_name} must be 1, 2, or 3.")
        return value

    def _validate_atom_limit(self, value: Optional[int], field_name: str) -> Optional[int]:
        if value is None:
            return None
        if value < 0:
            raise ValueError(f"{field_name} must be >= 0.")
        return value

    def _normalize_choices(
        self,
        values: Optional[list[str]],
        mapping: dict[str, str],
        field_name: str,
    ) -> list[str]:
        if not values:
            return []

        normalized: list[str] = []
        for raw in values:
            key = self._clean(raw)
            if key not in mapping:
                allowed = ", ".join(sorted(mapping.keys()))
                raise ValueError(f"Invalid {field_name}: {raw!r}. Allowed values include: {allowed}")

            short_value = mapping[key]
            if short_value not in normalized:
                normalized.append(short_value)

        return normalized

    def _infer_release_mechanism(self, target_smiles: str) -> str:
        has_ring = any(ch.isdigit() for ch in target_smiles)
        has_heteroatom = any(ch in target_smiles for ch in ["O", "N", "o", "n"])
        return "cyclization" if has_ring and has_heteroatom else "thiolysis"

    def _payload_preview(self, payload: list[tuple[str, str]]) -> dict[str, Any]:
        preview: dict[str, Any] = {}

        for key, value in payload:
            if key in preview:
                if isinstance(preview[key], list):
                    preview[key].append(value)
                else:
                    preview[key] = [preview[key], value]
            else:
                preview[key] = value

        return preview

    def _build_payload(
        self,
        target_smiles: str,
        use_pks: bool,
        use_bio: bool,
        use_chem: bool,
        max_bio_steps: Optional[int],
        max_chem_steps: Optional[int],
        pks_release_mechanism: Optional[str],
        pks_starters: Optional[list[str]],
        pks_extenders: Optional[list[str]],
        max_carbon: Optional[int],
        max_nitrogen: Optional[int],
        max_oxygen: Optional[int],
        auto_optimize_unspecified: bool,
    ) -> tuple[list[tuple[str, str]], dict[str, Any]]:
        target_smiles = target_smiles.strip()

        if not target_smiles:
            raise ValueError("target_smiles must not be empty.")
        if "." in target_smiles:
            raise ValueError("TridentSynth accepts one molecule at a time. Remove '.' from the SMILES.")
        if not any([use_pks, use_bio, use_chem]):
            raise ValueError("At least one synthesis strategy must be selected.")

        max_bio_steps = self._validate_step_count(max_bio_steps, "max_bio_steps")
        max_chem_steps = self._validate_step_count(max_chem_steps, "max_chem_steps")
        max_carbon = self._validate_atom_limit(max_carbon, "max_carbon")
        max_nitrogen = self._validate_atom_limit(max_nitrogen, "max_nitrogen")
        max_oxygen = self._validate_atom_limit(max_oxygen, "max_oxygen")

        starters = self._normalize_choices(pks_starters, self.starter_map, "pks_starters")
        extenders = self._normalize_choices(pks_extenders, self.extender_map, "pks_extenders")

        auto_filled: dict[str, Any] = {}

        if auto_optimize_unspecified:
            if max_bio_steps is None:
                max_bio_steps = 1
                auto_filled["max_bio_steps"] = 1

            if max_chem_steps is None:
                max_chem_steps = 1
                auto_filled["max_chem_steps"] = 1

            if use_pks and pks_release_mechanism is None:
                pks_release_mechanism = self._infer_release_mechanism(target_smiles)
                auto_filled["pks_release_mechanism"] = pks_release_mechanism

            if use_pks and not starters:
                starters = ["mal", "mmal"]
                auto_filled["pks_starters"] = starters

            if use_pks and not extenders:
                extenders = ["mal", "mmal"]
                auto_filled["pks_extenders"] = extenders

        if pks_release_mechanism is not None:
            pks_release_mechanism = self._clean(pks_release_mechanism)
            if pks_release_mechanism not in self.allowed_release_mechanisms:
                raise ValueError("pks_release_mechanism must be 'thiolysis' or 'cyclization'.")

        payload: list[tuple[str, str]] = []

        payload.append(("smiles", target_smiles))

        if use_pks:
            payload.append(("synthesisStrategy_pks", "on"))
        if use_bio:
            payload.append(("synthesisStrategy_bio", "on"))
        if use_chem:
            payload.append(("synthesisStrategy_chem", "on"))

        # Exact live TridentSynth field names from browser Form Data.
        payload.append(("rangebio", str(max_bio_steps if max_bio_steps is not None else 1)))
        payload.append(("rangechem", str(max_chem_steps if max_chem_steps is not None else 1)))

        if use_pks:
            payload.append(("releaseMechanism", pks_release_mechanism or "thiolysis"))

            for starter in starters:
                payload.append(("pksStarters[]", starter))

            for extender in extenders:
                payload.append(("pksExtenders[]", extender))

        # Live form includes these even when blank.
        payload.append(("maxAtomsC", "" if max_carbon is None else str(max_carbon)))
        payload.append(("maxAtomsN", "" if max_nitrogen is None else str(max_nitrogen)))
        payload.append(("maxAtomsO", "" if max_oxygen is None else str(max_oxygen)))

        normalized_query = {
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
            "auto_filled": auto_filled,
            "payload_preview": self._payload_preview(payload),
        }

        return payload, normalized_query

    def _submit_payload(
        self,
        session: requests.Session,
        payload: list[tuple[str, str]],
        acknowledge_controlled_substance_warning: bool,
    ) -> requests.Response:
        session.get(self.base_url, timeout=60)
        response = session.post(self.base_url, data=payload, timeout=120)

        if response.ok and "Controlled Substance Warning" in response.text:
            if not acknowledge_controlled_substance_warning:
                raise RuntimeError(
                    "TridentSynth returned its controlled-substance warning. "
                    "The tool stopped because acknowledge_controlled_substance_warning=False."
                )

            warning_soup = BeautifulSoup(response.text, "html.parser")
            warning_form = warning_soup.find("form")
            if warning_form is None:
                raise RuntimeError("Controlled-substance warning appeared, but no continue form was found.")

            warning_payload: list[tuple[str, str]] = []
            for tag in warning_form.find_all(["input", "button"]):
                name = tag.get("name")
                value = tag.get("value", "")
                if name:
                    warning_payload.append((name, value))

            action = warning_form.get("action") or response.url
            method = (warning_form.get("method") or "post").lower()
            continue_url = urljoin(response.url, action)

            if method == "get":
                response = session.get(continue_url, params=warning_payload, timeout=120)
            else:
                response = session.post(continue_url, data=warning_payload, timeout=120)

        if not response.ok:
            raise RuntimeError(
                f"TridentSynth submission failed with HTTP {response.status_code}.\n"
                f"Payload preview: {self._payload_preview(payload)}\n"
                f"Response excerpt: {response.text[:1500]}"
            )

        return response

    def _extract_task_id(self, response: requests.Response) -> Optional[str]:
        try:
            payload = response.json()
            if isinstance(payload, dict):
                return payload.get("task_id") or payload.get("job_id")
        except Exception:
            pass

        text = response.text
        match = re.search(r'"(?:task_id|job_id)"\s*:\s*"([^"]+)"', text)
        if match:
            return match.group(1)

        soup = BeautifulSoup(text, "html.parser")
        page_text = soup.get_text("\n", strip=True)
        lines = [line.strip() for line in page_text.splitlines() if line.strip()]

        for i, line in enumerate(lines):
            if line.lower().rstrip(":") in {"job id", "task id"} and i + 1 < len(lines):
                return lines[i + 1]

        return None

    def _candidate_result_urls(self, task_id: str) -> list[str]:
        return [
            urljoin(self.base_url, f"results/{task_id}/"),
        ]

    def _fetch_result_page(
        self,
        session: requests.Session,
        task_id: str,
        payload: list[tuple[str, str]],
    ) -> tuple[Optional[str], Optional[str]]:
        url = urljoin(self.base_url, f"results/{task_id}/")

        try:
            response = session.get(
                url,
                timeout=30,
                headers={
                    "Referer": self.base_url,
                    "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
                },
            )
        except Exception:
            return None, None

        if response.status_code != 200:
            return None, None

        text = response.text

        if (
            task_id in text
            or "TridentSynth job summary" in text
            or "Pathways to target found" in text
            or "Job submitted. Please wait" in text
            or "PKS product" in text
            or "Post-PKS product" in text
            or "Synthesis parameters" in text
            or "Full pathway design" in text
        ):
            return url, text

        return None, None

    def _is_result_complete(self, html: str) -> bool:
        text = BeautifulSoup(html, "html.parser").get_text(" ", strip=True).lower()

        if "job submitted. please wait" in text:
            return False
        if "pending" in text or "running" in text or "queued" in text:
            return False

        completion_phrases = [
            "pathways to target found",
            "closest reachable product",
            "the target can be synthesized by pks assembly alone",
            "full pathway design",
            "pks product",
            "post-pks product",
            "synthesis parameters",
        ]

        return any(phrase in text for phrase in completion_phrases)

    def _text_lines(self, soup: BeautifulSoup) -> list[str]:
        text = soup.get_text("\n", strip=True)
        return [line.strip() for line in text.splitlines() if line.strip()]

    def _value_after_label(self, lines: list[str], label: str) -> Optional[str]:
        label_clean = self._clean(label).rstrip(":")

        for i, line in enumerate(lines):
            if self._clean(line).rstrip(":") == label_clean:
                if i + 1 < len(lines):
                    return lines[i + 1].strip("` ")

        return None

    def _looks_like_smiles(self, value: str) -> bool:
        if not value:
            return False

        value = value.strip("` ").strip()

        if len(value) < 2:
            return False

        if value.lower() in {"i", "info", "none", "null", "nan"}:
            return False

        return bool(re.search(r"[CONFPSIBrcnos\[\]\(\)=#@+\-\\/0-9]", value))

    def _unique_ordered(self, values: list[Optional[str]]) -> list[str]:
        out: list[str] = []

        for value in values:
            if not value:
                continue

            value = value.strip("` ").strip()

            if not self._looks_like_smiles(value):
                continue

            if value not in out:
                out.append(value)

        return out

    def _extract_best_pathway_block(self, text: str) -> str:
        patterns = [
            r"Full pathway design #1(.*?)(Full pathway design #2|$)",
            r"Pathway 1(.*?)(Pathway 2|$)",
            r"Best pathway(.*?)(Full pathway design #2|Pathway 2|$)",
            r"Top pathway(.*?)(Full pathway design #2|Pathway 2|$)",
        ]

        for pattern in patterns:
            match = re.search(pattern, text, flags=re.IGNORECASE | re.DOTALL)
            if match:
                return match.group(1)

        return text

    def _extract_product_smiles(self, text: str, label: str) -> Optional[str]:
        pattern = rf"{re.escape(label)}.*?`([^`]+)`"
        match = re.search(pattern, text, flags=re.IGNORECASE | re.DOTALL)
        if match:
            value = match.group(1).strip()
            return value if self._looks_like_smiles(value) else None

        lines = [line.strip() for line in text.splitlines() if line.strip()]
        label_clean = self._clean(label)

        for i, line in enumerate(lines):
            if self._clean(line).rstrip(":") == label_clean.rstrip(":"):
                if i + 1 < len(lines):
                    value = lines[i + 1].strip("` ")
                    return value if self._looks_like_smiles(value) else None

        return None

    def _extract_float_after_label(self, text: str, label: str) -> Optional[float]:
        pattern = rf"{re.escape(label)}.*?([0-9]+(?:\.[0-9]+)?)"
        match = re.search(pattern, text, flags=re.IGNORECASE | re.DOTALL)
        if not match:
            return None

        try:
            return float(match.group(1))
        except ValueError:
            return None

    def _extract_reaction_smiles(self, text: str) -> list[str]:
        cleaned: list[str] = []

        section_match = re.search(
            r"Reactions \(SMILES\)(.*?)(Reaction rules|Step feasibilities|Net feasibility|Full pathway design #2|Pathway 2|$)",
            text,
            flags=re.IGNORECASE | re.DOTALL,
        )

        search_text = section_match.group(1) if section_match else text

        candidates = re.findall(r"`([^`]*>>[^`]*)`", search_text)

        if not candidates:
            candidates = re.findall(
                r"([A-Za-z0-9@\+\-\[\]\(\)=#$\\/%\.:]+>>[A-Za-z0-9@\+\-\[\]\(\)=#$\\/%\.:]+)",
                search_text,
            )

        for item in candidates:
            item = item.strip("` ").strip()
            if ">>" not in item:
                continue
            if item not in cleaned:
                cleaned.append(item)

        return cleaned

    def _reaction_to_structures(self, reaction_smiles: str) -> dict[str, list[str]]:
        if ">>" not in reaction_smiles:
            return {"reactants": [], "products": []}

        left, right = reaction_smiles.split(">>", 1)

        return {
            "reactants": [s for s in left.split(".") if self._looks_like_smiles(s)],
            "products": [s for s in right.split(".") if self._looks_like_smiles(s)],
        }

    def _extract_pks_modules(self, text: str) -> list[dict[str, Any]]:
        modules: list[dict[str, Any]] = []

        # Only parse before downstream post-PKS pathway sections.
        pre_pathway_text = re.split(
            r"Post-PKS pathways|Pathway 1|Reaction rules|Reaction enthalpies",
            text,
            maxsplit=1,
            flags=re.IGNORECASE,
        )[0]

        pattern = re.compile(
            r"MODULE\s+(\d+)\s+\((.*?)\)(.*?)(?=MODULE\s+\d+\s+\(|Domain legend|Post-PKS pathways|$)",
            flags=re.IGNORECASE | re.DOTALL,
        )

        seen_module_numbers: set[int] = set()

        for match in pattern.finditer(pre_pathway_text):
            module_number = int(match.group(1))

            # If Module 1, 2, etc. appears again, the page has moved to another candidate design.
            # Stop so we only return the first/best PKS module set.
            if module_number in seen_module_numbers:
                break

            seen_module_numbers.add(module_number)

            module_type = re.sub(r"\s+", " ", match.group(2)).strip()
            module_body = match.group(3)

            domains: list[dict[str, Optional[str]]] = []

            domain_candidates = re.findall(
                r"\b(KSq|KS|AT|KR|DH|ER|ACP)\b(?:\s*\(substrate:\s*([^)]+)\))?",
                module_body,
                flags=re.IGNORECASE,
            )

            for domain, substrate in domain_candidates:
                canonical_domain = "KSq" if domain.lower() == "ksq" else domain.upper()
                clean_substrate = substrate.strip() if substrate else None

                domains.append(
                    {
                        "domain": canonical_domain,
                        "substrate": clean_substrate,
                    }
                )

            domain_text_parts = []
            for item in domains:
                if item["substrate"]:
                    domain_text_parts.append(f"{item['domain']} substrate {item['substrate']}")
                else:
                    domain_text_parts.append(str(item["domain"]))

            modules.append(
                {
                    "module_number": module_number,
                    "module_type": module_type,
                    "domains": domains,
                    "text_summary": f"Module {module_number} ({module_type}): " + ", ".join(domain_text_parts),
                }
            )

        return modules

    def _extract_domain_legend(self, text: str) -> Optional[str]:
        match = re.search(
            r"Domain legend\s*(.*?)(Post-PKS pathways|Pathway 1|$)",
            text,
            flags=re.IGNORECASE | re.DOTALL,
        )
        if not match:
            return None

        legend = re.sub(r"\s+", " ", match.group(1)).strip()
        return legend or None

    def _extract_reaction_rule_names(self, text: str) -> list[str]:
        match = re.search(
            r"Reaction rules\s*(.*?)(Reaction enthalpies|Step feasibilities|Net feasibility|Full pathway design #2|Pathway 2|$)",
            text,
            flags=re.IGNORECASE | re.DOTALL,
        )
        if not match:
            return []

        section = match.group(1)
        lines = [line.strip() for line in section.splitlines() if line.strip()]

        cleaned: list[str] = []
        for line in lines:
            if line.lower() in {"reaction rules", "none"}:
                continue
            if line not in cleaned:
                cleaned.append(line)

        return cleaned

    def _extract_reaction_enthalpies(self, text: str) -> list[str]:
        match = re.search(
            r"Reaction enthalpies\s*\(kcal/mol\)\s*(.*?)(Step feasibilities|Net feasibility|Full pathway design #2|Pathway 2|$)",
            text,
            flags=re.IGNORECASE | re.DOTALL,
        )
        if not match:
            return []

        section = match.group(1)
        enthalpies = re.findall(r"-?\d+(?:\.\d+)?\s*kcal/mol", section)

        return list(dict.fromkeys(enthalpies))

    def _parse_result_page(self, html: str, task_id: Optional[str]) -> dict[str, Any]:
        soup = BeautifulSoup(html, "html.parser")
        lines = self._text_lines(soup)
        full_text = soup.get_text("\n", strip=True)
        best_block = self._extract_best_pathway_block(full_text)

        target_smiles = self._value_after_label(lines, "Target SMILES")
        target_name = self._value_after_label(lines, "Target Name")
        pathway_sequence = self._value_after_label(lines, "Pathway Sequence")
        pks_termination_step = self._value_after_label(lines, "PKS Termination Step")
        pks_extenders = self._value_after_label(lines, "PKS Extenders")
        pks_starters = self._value_after_label(lines, "PKS Starters")
        bio_steps = self._value_after_label(lines, "# Bio Steps")
        chem_steps = self._value_after_label(lines, "# Chem Steps")
        job_id = self._value_after_label(lines, "Job Id") or task_id

        pks_product = self._extract_product_smiles(best_block, "PKS product")
        post_pks_product = self._extract_product_smiles(best_block, "Post-PKS product")

        similarity_matches = re.findall(
            r"similarity to target.*?([0-9]+(?:\.[0-9]+)?)",
            best_block,
            flags=re.IGNORECASE | re.DOTALL,
        )

        pks_similarity = float(similarity_matches[0]) if len(similarity_matches) >= 1 else None
        post_pks_similarity = float(similarity_matches[1]) if len(similarity_matches) >= 2 else None

        net_feasibility = self._extract_float_after_label(best_block, "Net feasibility")

        reaction_smiles = self._extract_reaction_smiles(best_block)

        # Fallback safety: for pages where pathway sections are not clearly separated,
        # only keep the first/top reaction so alternative reactions do not flood output.
        if len(reaction_smiles) > 1:
            reaction_smiles = reaction_smiles[:1]

        reaction_structures = [self._reaction_to_structures(rxn) for rxn in reaction_smiles]

        reaction_rule_ids = self._unique_ordered(re.findall(r"rule\d+_\d+", best_block))
        reaction_rule_names = self._extract_reaction_rule_names(best_block)
        reaction_enthalpies = self._extract_reaction_enthalpies(best_block)

        pks_modules = self._extract_pks_modules(full_text)
        domain_legend = self._extract_domain_legend(full_text)

        step_feasibilities: list[float] = []
        step_section = re.search(
            r"Step feasibilities(.*?)(Pathway \d+|Full pathway design #\d+|$)",
            best_block,
            flags=re.IGNORECASE | re.DOTALL,
        )
        if step_section:
            for raw in re.findall(r"(?<!rule)(?<!\d)(0\.\d+|1\.0+|1)(?!\d)", step_section.group(1)):
                try:
                    step_feasibilities.append(float(raw))
                except ValueError:
                    pass

        structure_values: list[Optional[str]] = [pks_product]

        for rxn_struct in reaction_structures:
            structure_values.extend(rxn_struct["reactants"])
            structure_values.extend(rxn_struct["products"])

        structure_values.append(post_pks_product)
        structure_values.append(target_smiles)

        pathway_structures_smiles = self._unique_ordered(structure_values)

        if "Pathways to target found!" in full_text:
            status_message = "Pathways to target found."
        elif "Closest reachable product" in full_text or "closest reachable product" in full_text:
            status_message = "Exact target not reached; closest reachable product shown."
        elif "The target can be synthesized by PKS assembly alone" in full_text:
            status_message = "Target can be synthesized by PKS assembly alone."
        elif "Job submitted. Please wait" in full_text:
            status_message = "Job still running."
        else:
            status_message = "Result page parsed."

        return {
            "task_id": job_id,
            "status_message": status_message,
            "synthesis_parameters": {
                "target_smiles": target_smiles,
                "target_name": target_name,
                "pathway_sequence": pathway_sequence,
                "pks_termination_step": pks_termination_step,
                "pks_extenders": pks_extenders,
                "pks_starters": pks_starters,
                "bio_steps": bio_steps,
                "chem_steps": chem_steps,
            },
            "pks_modules": pks_modules,
            "domain_legend": domain_legend,
            "best_pathway": {
                "pks_product_smiles": pks_product,
                "pks_similarity_to_target": pks_similarity,
                "post_pks_product_smiles": post_pks_product,
                "post_pks_similarity_to_target": post_pks_similarity,
                "net_feasibility": net_feasibility,
                "reaction_smiles": reaction_smiles,
                "reaction_structures": reaction_structures,
                "reaction_rule_ids": reaction_rule_ids,
                "reaction_rule_names": reaction_rule_names,
                "reaction_enthalpies": reaction_enthalpies,
                "step_feasibilities": step_feasibilities,
                "pathway_structures_smiles": pathway_structures_smiles,
            },
        }

    def _add_selected_steps(self, parsed: dict[str, Any], submitted_query: dict[str, Any]) -> dict[str, Any]:
        params = parsed.get("synthesis_parameters", {})

        parsed["selected_steps"] = {
            "bio_steps": params.get("bio_steps") if submitted_query.get("use_bio") else None,
            "chem_steps": params.get("chem_steps") if submitted_query.get("use_chem") else None,
        }

        return parsed

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
        timeout_seconds: int = 600,
        poll_seconds: int = 5,
        acknowledge_controlled_substance_warning: bool = False,
    ) -> dict[str, Any]:
        session = self._session()

        payload, normalized_query = self._build_payload(
            target_smiles=target_smiles,
            use_pks=use_pks,
            use_bio=use_bio,
            use_chem=use_chem,
            max_bio_steps=max_bio_steps,
            max_chem_steps=max_chem_steps,
            pks_release_mechanism=pks_release_mechanism,
            pks_starters=pks_starters,
            pks_extenders=pks_extenders,
            max_carbon=max_carbon,
            max_nitrogen=max_nitrogen,
            max_oxygen=max_oxygen,
            auto_optimize_unspecified=auto_optimize_unspecified,
        )

        submit_response = self._submit_payload(
            session=session,
            payload=payload,
            acknowledge_controlled_substance_warning=acknowledge_controlled_substance_warning,
        )

        task_id = self._extract_task_id(submit_response)

        if not task_id:
            if self._is_result_complete(submit_response.text):
                parsed = self._add_selected_steps(
                    self._parse_result_page(submit_response.text, task_id=None),
                    normalized_query,
                )
                return {
                    "ok": True,
                    "status": "completed",
                    "submitted_query": normalized_query,
                    "result": parsed,
                }

            raise RuntimeError(
                "TridentSynth submission succeeded, but no task_id or result page could be parsed.\n"
                f"Payload preview: {self._payload_preview(payload)}\n"
                f"Response excerpt: {submit_response.text[:1500]}"
            )

        if not wait_for_completion:
            return {
                "ok": True,
                "status": "submitted",
                "task_id": task_id,
                "submitted_query": normalized_query,
                "message": "Job submitted successfully. Re-run with wait_for_completion=True to parse the best pathway.",
            }

        deadline = time.time() + max(timeout_seconds, 1)
        last_html: Optional[str] = None

        while time.time() <= deadline:
            result_url, result_html = self._fetch_result_page(session, task_id, payload)

            if result_html:
                last_html = result_html

                if self._is_result_complete(result_html):
                    parsed = self._add_selected_steps(
                        self._parse_result_page(result_html, task_id=task_id),
                        normalized_query,
                    )
                    return {
                        "ok": True,
                        "status": "completed",
                        "task_id": task_id,
                        "submitted_query": normalized_query,
                        "result": parsed,
                    }

            time.sleep(max(poll_seconds, 1))

        timed_out_result = None
        if last_html:
            timed_out_result = self._add_selected_steps(
                self._parse_result_page(last_html, task_id=task_id),
                normalized_query,
            )

        return {
            "ok": False,
            "status": "timed_out",
            "task_id": task_id,
            "submitted_query": normalized_query,
            "partial_result": timed_out_result,
            "message": f"Timed out after {timeout_seconds} seconds while waiting for TridentSynth results.",
        }


_instance = TridentSynth()
_instance.initiate()
tridentsynth = _instance.run