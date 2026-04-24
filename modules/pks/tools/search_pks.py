"""
SearchPKS — SBSPKS v2 structural-similarity and pathway search.

Implementation notes
--------------------
The SBSPKS v2 server (http://202.54.226.228/~pksdb/retro) exposes three CGI
endpoints intended for structure search:

  search_similar1.cgi  — "Search similar chemical structure"
  search_similar.cgi   — linked from display_smiles pages
  search_functional.cgi — "Search tailoring reactions"

All three fail with the server-side error::

    cannot make dir TMP_1 at search_similar1.cgi line 148.

This is a persistent permission bug on the remote server that prevents any
temp-directory creation for the Perl CGI scripts.  The following endpoints DO
work correctly and are used here instead:

  display_smiles.cgi?smile=<path_key>.smiles  — returns compound SMILES
  make_reaction.cgi?path=<path_key>           — returns full pathway data
  keyword.cgi  (POST, field "text")           — keyword-based pathway search

Strategy
--------
1.  A static catalog of all ~225 SBSPKS compound names / path keys is embedded
    in this module.
2.  SMILES are fetched from display_smiles.cgi and cached to disk; two known
    SMILES (erythromycin, pikromycin) are hard-wired as a fall-back seed.
3.  RDKit Morgan fingerprints + Tanimoto coefficient provide local similarity.
4.  Top hits are enriched with live pathway data from make_reaction.cgi.
"""

from __future__ import annotations

import io
import json
import os
import re
import time
import urllib.request
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Any, Dict, List, Literal, Optional

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_BASE_URL = "http://202.54.226.228/~pksdb/retro"
_REQUEST_TIMEOUT = 12          # seconds per individual HTTP request
_MAX_WORKERS = 10              # parallel SMILES-fetch threads
_CACHE_TTL = 30 * 24 * 3600   # disk-cache TTL: 30 days

# Path to the on-disk SMILES cache, co-located with this file
_CACHE_FILE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..", "data", "sbspks_smiles_cache.json",
)

# Path to the on-disk MIBiG compound cache
_MIBIG_CACHE_FILE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..", "data", "mibig_cache.json",
)

_MIBIG_ZIP_URL = (
    "https://github.com/mibig-secmet/mibig-json"
    "/archive/refs/heads/master.zip"
)

# ---------------------------------------------------------------------------
# Compound catalog: (path_key, display_name, pathway_type)
# Extracted verbatim from the SBSPKS v2 sidebar navigation.
# ---------------------------------------------------------------------------

_CATALOG: list[tuple[str, str, str]] = [
    # --- Modular PKS ---
    ("9methylstreptimidone", "9-Methyl-streptimidone",  "Modular PKS"),
    ("ambruticin",           "Ambruticin",              "Modular PKS"),
    ("amphotericin",         "Amphotericin",            "Modular PKS"),
    ("ansamitocin",          "Ansamitocin",             "Modular PKS"),
    ("apoptolidin",          "Apoptolidin",             "Modular PKS"),
    ("ascomycin",            "Ascomycin",               "Modular PKS"),
    ("aureothin",            "Aureothin",               "Modular PKS"),
    ("avermectin",           "Avermectin",              "Modular PKS"),
    ("bafilomycin",          "Bafilomycin",             "Modular PKS"),
    ("BE14106",              "BE-14106",                "Modular PKS"),
    ("borrelidin",           "Borrelidin",              "Modular PKS"),
    ("chalcomycin",          "Chalcomycin",             "Modular PKS"),
    ("chlorothricin",        "Chlorothricin",           "Modular PKS"),
    ("chondramide",          "Chondramide",             "Modular PKS"),
    ("concanamycin",         "Concanamycin",            "Modular PKS"),
    ("coronafacicacid",      "Coronafacic acid",        "Modular PKS"),
    ("cryptophycin",         "Cryptophycin",            "Modular PKS"),
    ("curacin",              "Curacin",                 "Modular PKS"),
    ("cycloheximide",        "Cycloheximide",           "Modular PKS"),
    ("cylindrospermopsin",   "Cylindrospermopsin",      "Modular PKS"),
    ("cystothiazole",        "CystothiazoleA",          "Modular PKS"),
    ("erythromycin",         "Erythromycin",            "Modular PKS"),
    ("FD891",                "FD-891",                  "Modular PKS"),
    ("fk506",                "FK506",                   "Modular PKS"),
    ("fostriecin",           "Fostriecin",              "Modular PKS"),
    ("fr008",                "Fr008",                   "Modular PKS"),
    ("geldanamycin",         "Geldanamycin",            "Modular PKS"),
    ("gephyronic",           "Gephyronic acid",         "Modular PKS"),
    ("halstoctacosanolide",  "Halstoctacosanolide",     "Modular PKS"),
    ("herbimycin",           "Herbimycin",              "Modular PKS"),
    ("herboxidiene",         "Herboxidiene",            "Modular PKS"),
    ("jerangolid",           "Jerangolid",              "Modular PKS"),
    ("kijanimicin",          "Kijanimicin",             "Modular PKS"),
    ("lankamycin",           "Lankamycin",              "Modular PKS"),
    ("lasalocid",            "Lasalocid",               "Modular PKS"),
    ("lipomycin",            "Lipomycin",               "Modular PKS"),
    ("megalomycin",          "Megalomycin",             "Modular PKS"),
    ("meilingmycin",         "Meilingmycin D",          "Modular PKS"),
    ("meridamycin",          "Meridamycin",             "Modular PKS"),
    ("methymycin",           "Methymycin",              "Modular PKS"),
    ("ML449",                "ML-449",                  "Modular PKS"),
    ("monensin",             "Monensin",                "Modular PKS"),
    ("mycinamicin",          "Mycinamycin",             "Modular PKS"),
    ("mycolactone",          "Mycolactone",             "Modular PKS"),
    ("myxalamide",           "Myxalamide",              "Modular PKS"),
    ("nanchangmycin",        "Nanchangmycin",           "Modular PKS"),
    ("narbomycin",           "Narbomycin",              "Modular PKS"),
    ("neoaureothin",         "Neoaureothin",            "Modular PKS"),
    ("niddamycin",           "Niddamycin",              "Modular PKS"),
    ("nystatin",             "Nystatin",                "Modular PKS"),
    ("oleandomycin",         "Oleandomycin",            "Modular PKS"),
    ("phoslactomycin",       "Phoslactomycin",          "Modular PKS"),
    ("piericidin",           "Piericidin",              "Modular PKS"),
    ("pikromycin",           "Pikromycin",              "Modular PKS"),
    ("pyoluteorin",          "Pyoluteorin",             "Modular PKS"),
    ("pyrrolomycin",         "Pyrrolomycin",            "Modular PKS"),
    ("rifamycin",            "Rifamycin",               "Modular PKS"),
    ("RK682",                "RK-682",                  "Modular PKS"),
    ("salinomycin",          "Salinomycin",             "Modular PKS"),
    ("soraphen",             "Soraphen",                "Modular PKS"),
    ("spinosad",             "Spinosad",                "Modular PKS"),
    ("spirangein",           "Spirangein",              "Modular PKS"),
    ("stigmatellin",         "Stigmatellin",            "Modular PKS"),
    ("tautomycetin",         "Tautomycetin",            "Modular PKS"),
    ("tautomycin",           "Tautomycin",              "Modular PKS"),
    ("tetrocarcin",          "Tetrocarcin A",           "Modular PKS"),
    ("tetronomycin",         "Tetronomycin",            "Modular PKS"),
    ("vicenistatin",         "Vicenistatin",            "Modular PKS"),
    # --- Trans-AT PKS ---
    ("bacilleane",           "Bacilleane",              "Trans-AT PKS"),
    ("bryostatin",           "Bryostatin",              "Trans-AT PKS"),
    ("chivosazol",           "Chivosazol",              "Trans-AT PKS"),
    ("difficidin",           "Difficidin",              "Trans-AT PKS"),
    ("disorazol",            "Disorazol",               "Trans-AT PKS"),
    ("isomigrastatin",       "ISO-Migrastatin",         "Trans-AT PKS"),
    ("kirromycin",           "Kirromycin",              "Trans-AT PKS"),
    ("lactimidomycin",       "Lactimidomycin",          "Trans-AT PKS"),
    ("macrolactin",          "Macrolactin",             "Trans-AT PKS"),
    ("mupirocin",            "Mupirocin",               "Trans-AT PKS"),
    ("myxovirescin",         "MyxovirescinA",           "Trans-AT PKS"),
    ("rhizoxin",             "Rhizoxin",                "Trans-AT PKS"),
    ("virginiamycin",        "Virginiamycin",           "Trans-AT PKS"),
    # --- Iterative PKS ---
    ("aflatoxin",            "Aflatoxin",               "Iterative PKS"),
    ("avilamycin",           "Avilamycin",              "Iterative PKS"),
    ("bikaverin",            "Bikaverin",               "Iterative PKS"),
    ("compactin",            "Compactin",               "Iterative PKS"),
    ("fumonisin",            "Fumonisin",               "Iterative PKS"),
    ("griseofulvin",         "Griseofulvin",            "Iterative PKS"),
    ("hypothemycin",         "Hypothemycin",            "Iterative PKS"),
    ("lovastatin",           "Lovastatin",              "Iterative PKS"),
    ("msas1",                "MSAS (G.loz)",            "Iterative PKS"),
    ("msas2",                "MSAS (A.par)",            "Iterative PKS"),
    ("msas3",                "MSAS (P.gri)",            "Iterative PKS"),
    ("msas4",                "MSAS (A.ter)",            "Iterative PKS"),
    ("msas5",                "MSAS (P.pat)",            "Iterative PKS"),
    ("naphthopyrone",        "Naphthopyrone",           "Iterative PKS"),
    ("radicicol",            "Radicicol",               "Iterative PKS"),
    ("sterigmatocystin",     "Sterigmatocystin",        "Iterative PKS"),
    ("THN1",                 "THN (A.fum)",             "Iterative PKS"),
    ("THN2",                 "THN (C.lag)",             "Iterative PKS"),
    ("THN3",                 "THN (E.der)",             "Iterative PKS"),
    ("viridicatumtoxin",     "Viridicatumtoxin",        "Iterative PKS"),
    # --- NRPS ---
    ("A500359",              "A-500359",                "NRPS"),
    ("A40926",               "A40926",                  "NRPS"),
    ("A47934",               "A47934",                  "NRPS"),
    ("A54145",               "A54145",                  "NRPS"),
    ("acinetobactin",        "Acinetobactin",           "NRPS"),
    ("actinomycin",          "Actinomycin",             "NRPS"),
    ("acv",                  "ACV",                     "NRPS"),
    ("aeruginoside",         "Aeruginoside",            "NRPS"),
    ("aeruginosin",          "Aeruginosin",             "NRPS"),
    ("anabaenopeptilide",    "Anabaenopeptilide",       "NRPS"),
    ("anabaenopeptin",       "Anabaenopeptin",          "NRPS"),
    ("anguibactin",          "Anguibactin",             "NRPS"),
    ("anthramycin",          "Anthramycin",             "NRPS"),
    ("apicidin",             "Apicidin",                "NRPS"),
    ("arthrofactin",         "Arthrofactin",            "NRPS"),
    ("bacillibactin",        "Bacillibactin",           "NRPS"),
    ("bacitracin",           "Bacitracin",              "NRPS"),
    ("balhimycin",           "Balhimycin",              "NRPS"),
    ("capreomycin",          "Capreomycin",             "NRPS"),
    ("cda",                  "CDA",                     "NRPS"),
    ("cephalosporin",        "Cephalosporin",           "NRPS"),
    ("cephamycin",           "Cephamycin-C",            "NRPS"),
    ("cereulide",            "Cereulide",               "NRPS"),
    ("chloramphenicol",      "Chloramphenicol",         "NRPS"),
    ("chloroeremomycin",     "Chloroeremomycin",        "NRPS"),
    ("complestatin",         "Complestatin",            "NRPS"),
    ("cyanopeptolin",        "Cyanopeptolin 984",       "NRPS"),
    ("cyclosporin",          "Cyclosporin",             "NRPS"),
    ("daptomycin",           "Daptomycin",              "NRPS"),
    ("echinomycin",          "Echinomycin",             "NRPS"),
    ("enduracidin",          "Enduracidin",             "NRPS"),
    ("enniatin",             "Enniatin",                "NRPS"),
    ("enterobactin",         "Enterobactin",            "NRPS"),
    ("ergovaline",           "Ergovaline",              "NRPS"),
    ("exochelin",            "Exochelin",               "NRPS"),
    ("fengycin",             "Fengycin",                "NRPS"),
    ("friulimicin",          "Friulimicin",             "NRPS"),
    ("fumitremorgin",        "Fumitremorgin",           "NRPS"),
    ("fusaricidin",          "Fusaricidin",             "NRPS"),
    ("gliotoxin",            "Gliotoxin",               "NRPS"),
    ("glycopeptilide",       "Glycopeptilide",          "NRPS"),
    ("gramicidin",           "Gramicidin",              "NRPS"),
    ("HC_toxin",             "HC-Toxin",                "NRPS"),
    ("hormaomycin",          "Hormaomycin",             "NRPS"),
    ("indigoidine",          "Indigoidine",             "NRPS"),
    ("kutzneride",           "Kutzneride",              "NRPS"),
    ("lichenycin",           "Lichenycin",              "NRPS"),
    ("lyngbyatoxin",         "Lyngbyatoxin",            "NRPS"),
    ("mannopeptimycins",     "Mannopeptimycins",        "NRPS"),
    ("massetolide",          "Massetolide",             "NRPS"),
    ("methylpendomycin",     "Methylpendomycin",        "NRPS"),
    ("nostocyclopeptide",    "Nostocyclopeptide",       "NRPS"),
    ("pacidamycin",          "Pacidamycin 3",           "NRPS"),
    ("paebacillibactin",     "Paebacillibactin",        "NRPS"),
    ("paenibactin",          "Paenibactin",             "NRPS"),
    ("phosphinothricin",     "Phosphinothricin",        "NRPS"),
    ("polymyxin",            "Polymyxin",               "NRPS"),
    ("pristinamycin",        "Pristinamycin",           "NRPS"),
    ("putisolvin",           "Putisolvin II",           "NRPS"),
    ("pyochelin",            "Pyochelin",               "NRPS"),
    ("pyoverdine",           "Pyoverdine",              "NRPS"),
    ("safracin",             "Safracin B",              "NRPS"),
    ("saframycin",           "Saframycin A",            "NRPS"),
    ("sibiromycin",          "Sibiromycin",             "NRPS"),
    ("spumigins",            "Spumigin A",              "NRPS"),
    ("surfactin",            "Surfactin",               "NRPS"),
    ("syringopeptin",        "Syringopeptin",           "NRPS"),
    ("taromycin",            "Taromycin A",             "NRPS"),
    ("teicoplanin",          "Teicoplanin",             "NRPS"),
    ("thaxtomin",            "Thaxtomin",               "NRPS"),
    ("thiocoraline",         "Thiocoraline",            "NRPS"),
    ("tyrocidine",           "Tyrocidine A",            "NRPS"),
    ("valinomycin",          "Valinomycin",             "NRPS"),
    ("vanchrobactin",        "Vanchrobactin",           "NRPS"),
    ("vibriobactin",         "Vibriobactin",            "NRPS"),
    ("vicibactin",           "Vicibactin",              "NRPS"),
    ("viomycin",             "Viomycin",                "NRPS"),
    ("xantholysin",          "Xantholysin A",           "NRPS"),
    # --- PKS-NRPS Hybrid ---
    ("ajudazol",             "Ajudazol",                "PKS-NRPS Hybrid"),
    ("anatoxin",             "Anatoxin",                "PKS-NRPS Hybrid"),
    ("barbamide",            "Barbamide",               "PKS-NRPS Hybrid"),
    ("caerulomycin",         "Caerulomycin",            "PKS-NRPS Hybrid"),
    ("chondrochlorens",      "Chondrochlorens",         "PKS-NRPS Hybrid"),
    ("colibactin",           "Colibactin",              "PKS-NRPS Hybrid"),
    ("collismycina",         "Collismycin A",           "PKS-NRPS Hybrid"),
    ("coronatine",           "Coronatine",              "PKS-NRPS Hybrid"),
    ("crocain",              "Crocain",                 "PKS-NRPS Hybrid"),
    ("epothilone",           "Epothilone",              "PKS-NRPS Hybrid"),
    ("fk520",                "FK520",                   "PKS-NRPS Hybrid"),
    ("guadinomine",          "Guadinomine",             "PKS-NRPS Hybrid"),
    ("hectochlorin",         "Hectochlorin",            "PKS-NRPS Hybrid"),
    ("indanomycin",          "Indanomycin",             "PKS-NRPS Hybrid"),
    ("iturin",               "Iturin",                  "PKS-NRPS Hybrid"),
    ("jamaicamide",          "Jamaicamide-A",           "PKS-NRPS Hybrid"),
    ("lankacidin",           "Lankacidin",              "PKS-NRPS Hybrid"),
    ("leinamycin",           "Leinamycin",              "PKS-NRPS Hybrid"),
    ("lymphostin",           "Lymphostin",              "PKS-NRPS Hybrid"),
    ("maduropeptin",         "Maduropeptin",            "PKS-NRPS Hybrid"),
    ("melithiazol",          "Melithiazol",             "PKS-NRPS Hybrid"),
    ("microcystin",          "Microcystin",             "PKS-NRPS Hybrid"),
    ("mycobactin",           "Mycobactin",              "PKS-NRPS Hybrid"),
    ("mycosubtilin",         "Mycosubtilin",            "PKS-NRPS Hybrid"),
    ("myxalamids",           "Myxalamid S",             "PKS-NRPS Hybrid"),
    ("myxochromide",         "Myxochromide",            "PKS-NRPS Hybrid"),
    ("myxothiazol",          "Myxothiazol",             "PKS-NRPS Hybrid"),
    ("naphthomycin",         "Naphthomycin A",          "PKS-NRPS Hybrid"),
    ("nodularin",            "Nodularin",               "PKS-NRPS Hybrid"),
    ("nostopeptilide",       "Nostopeptilide",          "PKS-NRPS Hybrid"),
    ("oxazolomycin",         "Oxazolomycin",            "PKS-NRPS Hybrid"),
    ("pederin",              "Pederin",                 "PKS-NRPS Hybrid"),
    ("pellasoren",           "Pellasoren",              "PKS-NRPS Hybrid"),
    ("prodigiosin",          "Prodigiosin",             "PKS-NRPS Hybrid"),
    ("psymberin",            "Psymberin",               "PKS-NRPS Hybrid"),
    ("rapamycin",            "Rapamycin",               "PKS-NRPS Hybrid"),
    ("simocyclinone",        "Simocyclinone",           "PKS-NRPS Hybrid"),
    ("streptolydigin",       "Streptolydigin",          "PKS-NRPS Hybrid"),
    ("tenellin",             "Tenellin",                "PKS-NRPS Hybrid"),
    ("thailandamide",        "Thailandamide",           "PKS-NRPS Hybrid"),
    ("thuggacin",            "Thuggacin A",             "PKS-NRPS Hybrid"),
    ("tirandamycin",         "Tirandamycin 3",          "PKS-NRPS Hybrid"),
    ("tubulysin",            "Tubulysin",               "PKS-NRPS Hybrid"),
    ("yersiniabactin",       "Yersiniabactin",          "PKS-NRPS Hybrid"),
    ("zearalenone",          "Zearalenone",             "PKS-NRPS Hybrid"),
    ("zorbamycin",           "Zorbamycin",              "PKS-NRPS Hybrid"),
    ("zwittermicin",         "Zwittermicin",            "PKS-NRPS Hybrid"),
]

# SMILES confirmed directly from the live server (display_smiles.cgi)
_SEED_SMILES: dict[str, str] = {
    # erythromycin.smiles — confirmed from SBSPKS display_smiles.cgi 2026-04-23
    "erythromycin": (
        "CCC1C(C(C(C(=O)C(CC(C(C(C(C(C(=O)O1)C)"
        "OC2CC(C(C(O2)C)O)(C)O)C)OC3C(C(CC(O3)C)N(C)C)O)(C)O)C)C)O)C"
    ),
    # pikromycin — confirmed from the SBSPKS example form
    "pikromycin": (
        "CCC1OC(=O)C(C)C(=O)C(C)C(OC2OC(C)CC(C2O)N(C)C)"
        "C(C)CC(C)C(=O)C=CC1(C)O"
    ),
}

# Curated organism mapping for common SBSPKS compounds
_ORGANISM_MAP: dict[str, str] = {
    "erythromycin":     "Saccharopolyspora erythraea",
    "pikromycin":       "Streptomyces venezuelae",
    "methymycin":       "Streptomyces venezuelae",
    "narbomycin":       "Streptomyces venezuelae",
    "niddamycin":       "Streptomyces caelestis",
    "oleandomycin":     "Streptomyces antibioticus",
    "megalomycin":      "Micromonospora megalomicea",
    "epothilone":       "Sorangium cellulosum",
    "rapamycin":        "Streptomyces hygroscopicus",
    "fk506":            "Streptomyces tsukubaensis",
    "fk520":            "Streptomyces hygroscopicus",
    "ascomycin":        "Streptomyces hygroscopicus",
    "rifamycin":        "Amycolatopsis mediterranei",
    "lovastatin":       "Aspergillus terreus",
    "compactin":        "Penicillium citrinum",
    "griseofulvin":     "Penicillium griseofulvum",
    "avermectin":       "Streptomyces avermitilis",
    "spinosad":         "Saccharopolyspora spinosa",
    "mupirocin":        "Pseudomonas fluorescens",
    "bryostatin":       "Candidatus Endobugula sertula",
    "difficidin":       "Bacillus subtilis",
    "aflatoxin":        "Aspergillus flavus",
    "bikaverin":        "Fusarium moniliforme",
    "fumonisin":        "Fusarium verticillioides",
    "soraphen":         "Sorangium cellulosum",
    "amphotericin":     "Streptomyces nodosus",
    "nystatin":         "Streptomyces noursei",
    "cyclosporin":      "Tolypocladium inflatum",
    "daptomycin":       "Streptomyces roseosporus",
    "surfactin":        "Bacillus subtilis",
    "geldanamycin":     "Streptomyces hygroscopicus",
    "herbimycin":       "Streptomyces hygroscopicus",
    "hypothemycin":     "Hypomyces subiculosus",
    "radicicol":        "Pochonia chlamydosporia",
    "salinomycin":      "Streptomyces albus",
    "monensin":         "Streptomyces cinnamonensis",
    "lasalocid":        "Streptomyces lasaliensis",
    "nanchangmycin":    "Streptomyces nanchangensis",
    "cystothiazole":    "Cystobacter fuscus",
    "epothilone":       "Sorangium cellulosum",
    "myxothiazol":      "Myxococcus fulvus",
    "rhizoxin":         "Burkholderia rhizoxinica",
    "mupirocin":        "Pseudomonas fluorescens",
}

# PKS module reaction keywords — used to exclude non-tailoring steps
_PKS_MODULE_KEYWORDS = frozenset({
    "KS condensation", "Loading Module", "KR reduction",
    "DH ", " DH", "ER ", " ER", "ACP", "AT domain",
    "propionate", "malonate", "methylmalonate", "methoxymalonyl",
    "+methylmalonate", "+malonate", "+propionate",
})


class SearchPKS:
    """
    Description:
        Queries the SBSPKS v2 Chemical Space database to find polyketides,
        NRPSs, and PKS-NRPS hybrids that are structurally similar to a query
        molecule, and retrieves their associated PKS/NRPS biosynthetic pathway
        context.

        Because the SBSPKS live search CGI endpoints have a persistent
        server-side directory-creation bug, similarity is computed locally
        using RDKit Morgan fingerprints (radius=2, 2048 bits) and the Tanimoto
        coefficient.  SMILES for the ~225 SBSPKS compounds are fetched from the
        working display_smiles.cgi endpoint and cached to disk.

    Input:
        query_smiles (str): SMILES string of the target molecule.
        search_type  (str): "reaction_search" — compare structure to all
                            SBSPKS final products; "pathway_search" — same
                            structural search but returns full pathway graphs
                            for every hit.
        max_results  (int): Maximum number of similar compounds to return.
        similarity_threshold (float): Minimum Tanimoto cutoff (0.0–1.0).

    Output:
        Dict with keys "query_smiles", "results", "search_type_used",
        "warnings".  See run() docstring for the full schema.

    Tests:
        - Case:
            Input: query_smiles="CCC1OC(=O)C(C)C(=O)C(C)C(OC2OC(C)CC(C2O)N(C)C)C(C)CC(C)C(=O)C=CC1(C)O",
                   search_type="reaction_search", max_results=3,
                   similarity_threshold=0.3
            Expected Output: results list containing at least one entry with
                             similarity_score >= 0.3 (pikromycin or related
                             macrolide should be the top hit)
            Description: Known PKS macrolide query; SBSPKS should contain
                         related macrolides.

        - Case:
            Input: query_smiles="OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
                   search_type="reaction_search", max_results=5,
                   similarity_threshold=0.6
            Expected Output: results list is empty OR all scores < 0.6
            Description: Glucose is not a polyketide; no close SBSPKS hits
                         expected above 0.6.

        - Case:
            Input: query_smiles="not-a-smiles!!",
                   search_type="reaction_search", max_results=5,
                   similarity_threshold=0.5
            Expected Output: ValueError raised
            Description: Invalid SMILES must raise ValueError.
    """

    # ------------------------------------------------------------------ #
    #  initiate                                                            #
    # ------------------------------------------------------------------ #

    def initiate(self) -> None:
        """
        One-time setup: import heavy dependencies, create HTTP session, and
        load (or build) the in-memory SMILES catalog used for similarity
        comparisons.

        If a fresh disk cache exists it is loaded immediately.  Otherwise
        SMILES are fetched from SBSPKS in parallel and the result is written
        to disk for future calls.  The seed SMILES for erythromycin and
        pikromycin are always available even when the server is unreachable.
        """
        self._import_error: Exception | None = None

        # Lazy-import heavy libraries so module load never crashes
        try:
            import requests
            from rdkit.Chem import MolFromSmiles, RDKFingerprint
            from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect
            from rdkit import DataStructs
            self._requests = requests
            self._MolFromSmiles = MolFromSmiles
            self._GetMorganFP = GetMorganFingerprintAsBitVect
            self._DataStructs = DataStructs
        except ImportError as exc:
            self._import_error = exc

        # HTTP session with a descriptive User-Agent
        if self._import_error is None:
            sess = self._requests.Session()
            sess.headers.update({"User-Agent": "PKSearch-BioE234/1.0 (educational)"})
            self._session = sess
        else:
            self._session = None

        # Build path_key → {display_name, pathway_type} lookup
        self._catalog_meta: dict[str, dict] = {
            pk: {"name": name, "pathway_type": ptype}
            for pk, name, ptype in _CATALOG
        }

        # Load SMILES catalog: disk cache → fall through to live fetch
        self._smiles_cache: dict[str, str] = dict(_SEED_SMILES)
        if self._import_error is None:
            self._smiles_cache.update(self._load_smiles_cache())

    # ------------------------------------------------------------------ #
    #  run                                                                 #
    # ------------------------------------------------------------------ #

    def run(
        self,
        query_smiles: str,
        search_type: Literal["reaction_search", "pathway_search"] = "reaction_search",
        max_results: int = 5,
        similarity_threshold: Optional[float] = 0.6,
    ) -> Dict[str, Any]:
        """
        Queries SBSPKS v2 for structurally similar polyketides and biosynthetic
        context using local RDKit Tanimoto similarity.

        Args:
            query_smiles (str): SMILES string of the target molecule.
            search_type (str): "reaction_search" compares structure to stored
                final-product SMILES for all ~225 SBSPKS compounds and returns
                ranked hits with basic metadata.  "pathway_search" performs the
                same structural search and additionally returns the full
                biosynthetic pathway graph (all reaction steps and intermediates)
                for each hit.
            max_results (int): Maximum number of similar compounds to return.
            similarity_threshold (float): Minimum Tanimoto score (0.0–1.0);
                compounds below this value are excluded.

        Returns:
            Dict[str, Any]:
                - "query_smiles" (str): The input SMILES.
                - "results" (List[Dict]): Ranked list of similar compounds, each
                    containing:
                    - "name" (str): Compound display name from SBSPKS.
                    - "smiles" (str): SMILES of the matched compound.
                    - "similarity_score" (float): Tanimoto coefficient (Morgan
                        fingerprint, radius=2, 2048 bits) against query.
                    - "pathway_name" (str): Pathway identifier (same as name).
                    - "pathway_type" (str): Modular PKS / Trans-AT PKS /
                        Iterative PKS / NRPS / PKS-NRPS Hybrid.
                    - "organism" (str): Source organism (curated map or
                        "Unknown").
                    - "tailoring_reactions" (List[str]): Post-PKS tailoring
                        steps parsed from the SBSPKS pathway graph
                        (e.g. "eryF (Hydroxylation)").
                    - "pathway_steps" (List[str]): Full reaction step labels
                        (present only when search_type="pathway_search").
                    - "sbspks_url" (str): Direct URL to the SBSPKS pathway page.
                - "search_type_used" (str): Which search mode was executed.
                - "warnings" (List[str]): Non-fatal issues encountered (e.g.
                    parse errors, network timeouts).

        Raises:
            ValueError: If query_smiles is empty or RDKit cannot parse it.
            RuntimeError: If required libraries (rdkit, requests) are not
                installed.
        """
        warnings: List[str] = []

        # ---- Dependency check ---------------------------------------- #
        if self._import_error is not None:
            raise RuntimeError(
                f"Required library not available: {self._import_error}. "
                "Ensure rdkit and requests are installed."
            )

        # ---- Input validation ---------------------------------------- #
        query_smiles = query_smiles.strip()
        if not query_smiles:
            raise ValueError("query_smiles must be a non-empty SMILES string.")

        query_mol = self._MolFromSmiles(query_smiles)
        if query_mol is None:
            raise ValueError(
                f"RDKit could not parse the provided SMILES: {query_smiles!r}. "
                "Ensure the string is a valid SMILES (not InChI, name, or CAS)."
            )

        if not isinstance(max_results, int) or max_results < 1:
            raise ValueError("max_results must be a positive integer.")

        if similarity_threshold is not None:
            if not (0.0 <= similarity_threshold <= 1.0):
                raise ValueError("similarity_threshold must be between 0.0 and 1.0.")

        valid_types = {"reaction_search", "pathway_search"}
        if search_type not in valid_types:
            raise ValueError(
                f"search_type must be one of {sorted(valid_types)}, "
                f"got {search_type!r}."
            )

        # ---- Compute query fingerprint -------------------------------- #
        query_fp = self._GetMorganFP(query_mol, radius=2, nBits=2048)

        # ---- Fetch any missing SMILES from the live server ------------ #
        # Only attempt for compounds not already in the cache
        missing = [
            pk for pk, _, _ in _CATALOG
            if pk not in self._smiles_cache
        ]
        if missing:
            newly_fetched = self._parallel_fetch_smiles(missing)
            self._smiles_cache.update(newly_fetched)
            if newly_fetched:
                self._save_smiles_cache(self._smiles_cache)
            if len(newly_fetched) < len(missing):
                warnings.append(
                    f"Could not fetch SMILES for "
                    f"{len(missing) - len(newly_fetched)} compounds "
                    f"(server timeout or parse error); those were skipped."
                )

        # ---- Compute Tanimoto for every compound with a known SMILES -- #
        scored: list[tuple[float, str]] = []  # (score, path_key)
        for path_key, smiles_str in self._smiles_cache.items():
            if path_key not in self._catalog_meta:
                continue
            mol = self._MolFromSmiles(smiles_str)
            if mol is None:
                continue
            fp = self._GetMorganFP(mol, radius=2, nBits=2048)
            score = self._DataStructs.TanimotoSimilarity(query_fp, fp)
            if similarity_threshold is None or score >= similarity_threshold:
                scored.append((score, path_key))

        # Sort descending by similarity, then take top N
        scored.sort(key=lambda x: x[0], reverse=True)
        top_hits = scored[:max_results]

        if not top_hits:
            warnings.append(
                "No compounds in the SBSPKS catalog met the similarity "
                f"threshold of {similarity_threshold}."
            )

        # ---- Enrich top hits with pathway metadata -------------------- #
        results: List[Dict[str, Any]] = []
        for score, path_key in top_hits:
            meta = self._catalog_meta[path_key]
            pathway_data = self._fetch_pathway_data(path_key, warnings)

            entry: Dict[str, Any] = {
                "name":               meta["name"],
                "smiles":             self._smiles_cache.get(path_key, ""),
                "similarity_score":   round(score, 4),
                "pathway_name":       meta["name"],
                "pathway_type":       meta["pathway_type"],
                "organism":           _ORGANISM_MAP.get(path_key, "Unknown"),
                "tailoring_reactions": pathway_data.get("tailoring_reactions", []),
                "sbspks_url": (
                    f"{_BASE_URL}/make_reaction.cgi?path={path_key}"
                ),
            }
            if search_type == "pathway_search":
                entry["pathway_steps"] = pathway_data.get("all_steps", [])

            results.append(entry)

        return {
            "query_smiles":     query_smiles,
            "results":          results,
            "search_type_used": search_type,
            "warnings":         warnings,
        }

    # ------------------------------------------------------------------ #
    #  Private helpers                                                     #
    # ------------------------------------------------------------------ #

    def _fetch_smiles_for_compound(self, path_key: str) -> str | None:
        """
        Fetch the final-product SMILES for one SBSPKS compound.

        Hits the display_smiles.cgi endpoint which embeds the SMILES in a
        hidden form field::

            <input type="hidden" name="hid_smile" value="<SMILES>" />

        Returns the SMILES string on success, or None on any error.
        """
        url = f"{_BASE_URL}/display_smiles.cgi?smile={path_key}.smiles"
        try:
            resp = self._session.get(url, timeout=_REQUEST_TIMEOUT)
            resp.raise_for_status()
        except Exception:
            return None

        # The SMILES is in the hidden field hid_smile
        m = re.search(r'name="hid_smile"\s+value="([^"]+)"', resp.text)
        if m:
            return m.group(1).strip()
        # Fallback: look for the plain-text SMILES paragraph
        m2 = re.search(
            r"<p[^>]*>\s*([A-Za-z0-9@\[\]()=#\-\+%\.\\\/]{15,})\s*&nbsp",
            resp.text,
        )
        if m2:
            return m2.group(1).strip()
        return None

    def _parallel_fetch_smiles(self, path_keys: list[str]) -> dict[str, str]:
        """
        Fetch SMILES for a list of path_keys concurrently and return a
        {path_key: smiles} dict for every successful fetch.
        """
        result: dict[str, str] = {}
        with ThreadPoolExecutor(max_workers=_MAX_WORKERS) as pool:
            future_map = {
                pool.submit(self._fetch_smiles_for_compound, pk): pk
                for pk in path_keys
            }
            for future in as_completed(future_map, timeout=60):
                pk = future_map[future]
                try:
                    smiles = future.result(timeout=_REQUEST_TIMEOUT + 2)
                    if smiles:
                        result[pk] = smiles
                except Exception:
                    pass
        return result

    def _fetch_pathway_data(
        self, path_key: str, warnings: List[str]
    ) -> dict[str, Any]:
        """
        Fetch the biosynthetic pathway page for a compound and parse it into
        structured data.

        The make_reaction.cgi endpoint embeds Cytoscape.js graph data as
        inline JavaScript.  Edge objects have the structure::

            { data: { id: '...', source: '...', target: '...',
                      label: '<enzyme>\\n(<reaction_type>)', href1: '...' },
              classes: 'wrapped' }

        Tailoring reactions are identified as edges whose label does NOT
        contain keywords characteristic of PKS elongation modules (e.g.
        "KS condensation", "Loading Module").

        Returns a dict with keys "tailoring_reactions" (List[str]) and
        "all_steps" (List[str]).  Returns empty lists on any network error.
        """
        url = f"{_BASE_URL}/make_reaction.cgi?path={path_key}"
        try:
            resp = self._session.get(url, timeout=_REQUEST_TIMEOUT)
            resp.raise_for_status()
        except Exception as exc:
            warnings.append(
                f"Could not fetch pathway data for {path_key}: {exc}"
            )
            return {"tailoring_reactions": [], "all_steps": []}

        # Extract all edge label strings from the embedded Cytoscape JS
        raw_labels = re.findall(r"label:\s*'([^']+)'", resp.text)

        all_steps: List[str] = []
        tailoring: List[str] = []

        for raw in raw_labels:
            # Convert JS \n escapes to a readable form
            clean = raw.replace("\\n", " ").strip()
            if not clean:
                continue
            all_steps.append(clean)

            # Classify: tailoring vs. PKS-module reaction
            lower = clean.lower()
            is_pks_module = any(
                kw.lower() in lower for kw in _PKS_MODULE_KEYWORDS
            )
            if not is_pks_module:
                tailoring.append(clean)

        return {"tailoring_reactions": tailoring, "all_steps": all_steps}

    # ------------------------------------------------------------------ #
    #  Disk cache helpers                                                  #
    # ------------------------------------------------------------------ #

    def _load_smiles_cache(self) -> dict[str, str]:
        """
        Load the on-disk SMILES cache if it exists and is not stale.

        Returns an empty dict if the file is missing, unreadable, or older
        than _CACHE_TTL seconds.
        """
        cache_path = _CACHE_FILE
        if not os.path.exists(cache_path):
            return {}
        age = time.time() - os.path.getmtime(cache_path)
        if age > _CACHE_TTL:
            return {}
        try:
            with open(cache_path, "r", encoding="utf-8") as fh:
                data = json.load(fh)
            if isinstance(data, dict):
                return {k: v for k, v in data.items() if isinstance(v, str)}
        except Exception:
            pass
        return {}

    def _save_smiles_cache(self, cache: dict[str, str]) -> None:
        """
        Persist the SMILES cache to disk.  Silently swallows I/O errors so
        a read-only filesystem never causes run() to fail.
        """
        cache_path = _CACHE_FILE
        try:
            os.makedirs(os.path.dirname(cache_path), exist_ok=True)
            with open(cache_path, "w", encoding="utf-8") as fh:
                json.dump(cache, fh, indent=2)
        except Exception:
            pass


# ---------------------------------------------------------------------------
# MIBiG index builder (Step 3)
# ---------------------------------------------------------------------------

def build_mibig_index(
    cache_path: str = _MIBIG_CACHE_FILE,
) -> list[dict]:
    """
    Download, parse, and cache the MIBiG polyketide compound index.

    MIBiG (Minimum Information about a Biosynthetic Gene Cluster) is a
    curated database of experimentally validated biosynthetic gene clusters.
    The full dataset is available as a ZIP archive of JSON files from GitHub
    (mibig-secmet/mibig-json, master branch).

    Each BGC JSON file has this structure::

        {
          "cluster": {
            "mibig_accession": "BGC0000055",
            "biosyn_class": ["Polyketide", "Saccharide"],
            "organism_name": "Saccharopolyspora erythraea NRRL 2338",
            "loci": { "accession": "AM420293.1" },
            "compounds": [
              {
                "compound": "erythromycin A",
                "chem_struct": "CC[C@@H]1...",   // SMILES
                ...
              },
              ...
            ]
          }
        }

    This function filters for entries whose ``biosyn_class`` list contains
    ``"Polyketide"``, then yields one dict per compound within each matching
    BGC that has a ``chem_struct`` (SMILES) value.  A single BGC can produce
    multiple compounds (e.g. erythromycins A, B, C, D all come from
    BGC0000055), so the number of returned dicts exceeds the number of BGCs.

    Results are cached to ``cache_path`` as a JSON file.  If a fresh cache
    (< 30 days old) exists, it is loaded directly without re-downloading.

    Args:
        cache_path (str): Absolute path for the on-disk JSON cache file.
            Defaults to ``modules/pks/data/mibig_cache.json``.

    Returns:
        list[dict]: One dict per polyketide compound with a known SMILES.
            Each dict contains:

            - ``compound_name`` (str): Compound name from MIBiG, e.g.
              ``"erythromycin A"``.
            - ``smiles`` (str): SMILES string from the ``chem_struct`` field.
            - ``bgc_accession`` (str): MIBiG BGC accession, e.g.
              ``"BGC0000055"``.  This is the ID the antiSMASH validator
              (Opshory's component) consumes directly.
            - ``organism`` (str): Source organism name, e.g.
              ``"Saccharopolyspora erythraea NRRL 2338"``.
            - ``genbank_accession`` (str | None): GenBank locus accession
              for the BGC, e.g. ``"AM420293.1"``, or None if not present.

    Raises:
        RuntimeError: If the MIBiG ZIP cannot be downloaded and no cache
            exists on disk.
    """
    # --- Load from cache if fresh ----------------------------------------
    if os.path.exists(cache_path):
        age = time.time() - os.path.getmtime(cache_path)
        if age < _CACHE_TTL:
            try:
                with open(cache_path, "r", encoding="utf-8") as fh:
                    cached = json.load(fh)
                if isinstance(cached, list) and cached:
                    return cached
            except Exception:
                pass  # corrupt cache — fall through to re-download

    # --- Download ZIP archive -------------------------------------------
    try:
        with urllib.request.urlopen(_MIBIG_ZIP_URL, timeout=60) as resp:
            zip_bytes = resp.read()
    except Exception as exc:
        raise RuntimeError(
            f"Failed to download MIBiG dataset from {_MIBIG_ZIP_URL}: {exc}. "
            "Check network connectivity."
        ) from exc

    # --- Parse JSON files from the ZIP -----------------------------------
    entries: list[dict] = []

    with zipfile.ZipFile(io.BytesIO(zip_bytes)) as zf:
        # JSON BGC files live under <root>/data/BGCxxxxxxx.json
        bgc_files = [
            n for n in zf.namelist()
            if n.endswith(".json") and "/data/" in n
        ]
        for filename in bgc_files:
            try:
                cluster_data = json.loads(zf.read(filename))
            except Exception:
                continue  # skip malformed JSON

            cluster = cluster_data.get("cluster", {})

            # Filter: must include Polyketide in biosynthetic class
            if "Polyketide" not in cluster.get("biosyn_class", []):
                continue

            bgc_accession: str = cluster.get("mibig_accession", "")
            organism: str = cluster.get("organism_name", "Unknown")
            genbank_accession: str | None = (
                cluster.get("loci", {}).get("accession") or None
            )

            # One entry per compound that has a SMILES string
            for compound in cluster.get("compounds", []):
                smiles = compound.get("chem_struct", "").strip()
                if not smiles:
                    continue
                compound_name = compound.get("compound", "Unknown compound")
                entries.append({
                    "compound_name":      compound_name,
                    "smiles":             smiles,
                    "bgc_accession":      bgc_accession,
                    "organism":           organism,
                    "genbank_accession":  genbank_accession,
                })

    if not entries:
        raise RuntimeError(
            "MIBiG ZIP was downloaded but no polyketide entries with SMILES "
            "were parsed.  The ZIP structure may have changed."
        )

    # --- Save to disk cache ----------------------------------------------
    try:
        os.makedirs(os.path.dirname(cache_path), exist_ok=True)
        with open(cache_path, "w", encoding="utf-8") as fh:
            json.dump(entries, fh, indent=2)
    except Exception:
        pass  # non-fatal: cache write failure never breaks the caller

    return entries


# ---------------------------------------------------------------------------
# Module-level convenience alias (matches the framework pattern)
# ---------------------------------------------------------------------------

_instance = SearchPKS()
_instance.initiate()
search_pks = _instance.run


# ---------------------------------------------------------------------------
# Standalone intermediate-extraction function (Step 2)
# ---------------------------------------------------------------------------

def extract_intermediates_from_pathway(html: str) -> list[dict]:
    """
    Parse the Cytoscape.js graph embedded in a make_reaction.cgi pathway page
    and return metadata for every biosynthetic node (intermediates + final
    product).

    The SBSPKS pathway pages embed graph data as inline JavaScript.  Each node
    looks like::

        { data: { id: 'erythromycin_11',
                  href:'/~pksdb/retro/display_smiles.cgi?smile=erythromycin.11.smiles'}},

    Each edge looks like::

        { data: { id: 'erythromycin_11-erythromycin_10',
                  source: 'erythromycin_11', target: 'erythromycin_10',
                  label: '+propionate\\neryAI\\n(Loading Module)',
                  href1: '...&&moduleno=1&&...' },
          classes: 'wrapped' },

    SMILES strings are NOT embedded in the HTML.  Each node's ``href`` field
    points to ``display_smiles.cgi?smile=<smiles_file>`` where
    ``<smiles_file>`` is e.g. ``erythromycin.11.smiles``.  The caller must
    fetch those separately.

    The node numbering is REVERSE order through the pathway: the node with the
    highest numeric suffix is the earliest intermediate (starter unit), and the
    unnumbered node (id == path_key) is the final product.  The
    ``module_number`` for each node is read from the ``moduleno`` parameter in
    the incoming edge's ``href1`` field — this is the authoritative step number
    as annotated in the SBSPKS database.

    Args:
        html (str): Raw HTML text from a make_reaction.cgi?path=<path_key>
            request.

    Returns:
        list[dict]: One dict per node, sorted by node_suffix descending
            (earliest intermediate first, final product last).  Each dict
            contains:

            - ``node_id`` (str): Full Cytoscape node ID, e.g.
              ``"erythromycin_11"``.
            - ``path_key`` (str): Parent pathway key, e.g.
              ``"erythromycin"``.
            - ``smiles_file`` (str): The ``smile=`` parameter value from the
              node href, e.g. ``"erythromycin.11.smiles"``.  Pass as
              ``display_smiles.cgi?smile=<smiles_file>`` to fetch the SMILES.
            - ``is_final_product`` (bool): True only for the unnumbered node
              whose id equals the path_key.
            - ``node_suffix`` (int | None): The numeric suffix from the node
              id (11, 10, … 1), or None for the final product node.
            - ``module_number`` (int | None): The PKS module step number from
              the incoming edge's ``moduleno=`` URL parameter (1-based), or
              None if no incoming edge was found (i.e. the first/source node).
            - ``step_label`` (str): Human-readable label of the edge that
              produces this node (e.g. ``"+methylmalonate eryAI
              (KS condensation + KR reduction)"``), or ``""`` for the first
              node which has no incoming edge.
    """
    # --- 1. Parse all node lines -----------------------------------------
    # Nodes have id + href but NO source/target/label fields.
    # Pattern: { data: { id: '<id>', href:'...?smile=<smiles_file>'}}
    node_pattern = re.compile(
        r"\{ data: \{ id: '([^']+)', href:'[^?]+\?smile=([^']+)'\}\}"
    )

    # --- 2. Parse all edge lines -----------------------------------------
    # Edges have source, target, label, href1 and end with classes: 'wrapped'
    # Extract: target node id, step label, moduleno from href1
    edge_pattern = re.compile(
        r"\{ data: \{ id: '[^']+', source: '([^']+)', target: '([^']+)'"
        r" ,label: '([^']+)', href1: '([^']*)'\},classes:"
    )
    moduleno_pattern = re.compile(r"moduleno=(\d+)")

    # Build: target_node_id → (step_label, module_number)
    incoming: dict[str, tuple[str, int | None]] = {}
    for m in edge_pattern.finditer(html):
        _source, target, raw_label, href1 = m.groups()
        clean_label = raw_label.replace("\\n", " ").strip()
        mod_m = moduleno_pattern.search(href1)
        module_number = int(mod_m.group(1)) if mod_m else None
        incoming[target] = (clean_label, module_number)

    # --- 3. Build result list --------------------------------------------
    suffix_pattern = re.compile(r"^(.+)_(\d+)$")

    results: list[dict] = []
    seen_ids: set[str] = set()

    for m in node_pattern.finditer(html):
        node_id, smiles_file = m.groups()

        # de-duplicate (Cytoscape sometimes repeats nodes)
        if node_id in seen_ids:
            continue
        seen_ids.add(node_id)

        # Determine path_key and numeric suffix
        suffix_m = suffix_pattern.match(node_id)
        if suffix_m:
            path_key = suffix_m.group(1)
            node_suffix: int | None = int(suffix_m.group(2))
            is_final_product = False
        else:
            # No suffix — this IS the path_key, i.e. the final product
            path_key = node_id
            node_suffix = None
            is_final_product = True

        step_label, module_number = incoming.get(node_id, ("", None))

        results.append({
            "node_id":          node_id,
            "path_key":         path_key,
            "smiles_file":      smiles_file,
            "is_final_product": is_final_product,
            "node_suffix":      node_suffix,
            "module_number":    module_number,
            "step_label":       step_label,
        })

    # Sort: highest suffix first (earliest in pathway), final product last
    results.sort(
        key=lambda d: (d["node_suffix"] is None, -(d["node_suffix"] or 0))
    )
    return results
