"""Documentation server (Model Context Protocol).

``python -m keta.mcp`` (or the ``keta-mcp`` command) starts a stdio MCP server
through which an MCP client searches and reads the keta documentation: the Python API
user manual, the calculation principles, the executed examples, and the
signature and docstring of every class, method and function of the installed
module. The server only reads documentation; it runs no calculation and no
code supplied by the client.

Requires the optional dependency ``mcp`` (``pip install "keta[mcp]"``).

Client configuration, for example::

    {"mcpServers": {"keta": {"command": "keta-mcp"}}}

The documentation corpus (``docs_corpus.json.gz``, built with the manuals) is
part of the package; the functions below can also be used without MCP.
"""

from __future__ import annotations

import gzip
import inspect
import json
import math
import re
import unicodedata
from collections import Counter
from functools import lru_cache
from importlib import resources

__all__ = ["search_documentation", "read_section", "list_contents", "list_examples", "read_example",
           "describe_api", "main"]

LANGUAGES = ("en", "fr")


@lru_cache(maxsize=1)
def _corpus() -> dict:
    with resources.files("keta").joinpath("docs_corpus.json.gz").open("rb") as f:
        return json.loads(gzip.decompress(f.read()).decode("utf-8"))


def _words(text: str) -> list[str]:
    t = unicodedata.normalize("NFKD", text.lower())
    t = "".join(c for c in t if not unicodedata.combining(c))
    return re.findall(r"[a-z0-9_]+", t)


@lru_cache(maxsize=2)
def _index(lang: str):
    """BM25 index of the sections of one language: (sections, term frequencies, document frequencies, mean length)."""
    secs = [s for s in _corpus()["sections"] if s["lang"] == lang]
    tfs = []
    df: Counter = Counter()
    for s in secs:
        # the title and the chapter title weigh three times the body
        words = _words(s["text"]) + 3 * (_words(s["title"]) + _words(s["chapter_title"]))
        tf = Counter(words)
        tfs.append((tf, len(words)))
        df.update(tf.keys())
    mean = sum(n for _, n in tfs) / max(len(tfs), 1)
    return secs, tfs, df, mean


def _check_lang(language: str) -> str:
    language = (language or "en").lower()
    if language not in LANGUAGES:
        raise ValueError(f"language must be one of {LANGUAGES}")
    return language


def _excerpt(text: str, terms: set[str], width: int = 300) -> str:
    """A window of the text around the first occurrence of a query term."""
    flat = re.sub(r"\s+", " ", text)
    pos = next((m.start() for m in re.finditer(r"\w+", flat) if (_words(m.group(0)) or [""])[0] in terms), 0)
    start = max(0, pos - width // 3)
    out = flat[start:start + width].strip()
    return ("… " if start else "") + out + (" …" if start + width < len(flat) else "")


def search_documentation(query: str, language: str = "en", limit: int = 8) -> list[dict]:
    """Search the user manual and the calculation principles.

    Returns the best matching sections (id, document, chapter, title, url, excerpt), most relevant
    first. Read a whole section with ``read_section(id)``.
    """
    lang = _check_lang(language)
    secs, tfs, df, mean = _index(lang)
    terms = set(_words(query))
    if not terms:
        return []
    n = len(secs)
    k1, b = 1.4, 0.75
    scores = []
    for s, (tf, length) in zip(secs, tfs):
        score = 0.0
        for t in terms:
            if t in tf:
                idf = math.log(1 + (n - df[t] + 0.5) / (df[t] + 0.5))
                score += idf * tf[t] * (k1 + 1) / (tf[t] + k1 * (1 - b + b * length / mean))
        if score > 0:
            scores.append((score, s))
    scores.sort(key=lambda x: -x[0])
    return [dict(id=s["id"], document=s["document"], chapter=f'{s["chapter"]}. {s["chapter_title"]}',
                 title=s["title"], url=s["url"], excerpt=_excerpt(s["text"], terms))
            for _, s in scores[:max(1, int(limit))]]


def read_section(section_id: str) -> dict:
    """Full text (Markdown, formulas in LaTeX) of a section, e.g. ``en/manual/14.2``.

    A chapter number alone (``en/manual/14``) returns the chapter introduction followed by the list of
    its sections.
    """
    secs = {s["id"]: s for s in _corpus()["sections"]}
    s = secs.get(section_id.strip())
    if s is None:
        raise KeyError(f"unknown section {section_id!r}; use search_documentation or list_contents")
    out = dict(id=s["id"], document=s["document"], chapter=f'{s["chapter"]}. {s["chapter_title"]}',
               title=s["title"], url=s["url"], text=s["text"])
    if s["number"] == s["chapter"]:
        prefix = s["id"].rsplit("/", 1)[0] + "/" + s["chapter"] + "."
        out["sections"] = [dict(id=x["id"], title=x["title"]) for x in _corpus()["sections"]
                           if x["id"].startswith(prefix)]
    return out


def list_contents(language: str = "en", document: str = "") -> list[dict]:
    """Table of contents: chapters with their numbered sections (``document``: 'manual', 'principles' or both)."""
    lang = _check_lang(language)
    out: list[dict] = []
    for s in _corpus()["sections"]:
        if s["lang"] != lang or (document and f"/{document}/" not in s["id"]):
            continue
        if s["number"] == s["chapter"]:
            out.append(dict(id=s["id"], document=s["document"], chapter=s["chapter"], title=s["chapter_title"],
                            sections=[]))
        elif out:
            out[-1]["sections"].append(dict(id=s["id"], title=f'{s["number"]} {s["title"]}'))
    return out


def list_examples(language: str = "en", query: str = "") -> list[dict]:
    """The runnable examples of the user manual (name and first comment line), optionally filtered by words."""
    lang = _check_lang(language)
    terms = set(_words(query))
    out = []
    for e in _corpus()["examples"]:
        if e["lang"] != lang:
            continue
        if terms and not terms & set(_words(e["code"]) + _words(e["name"])):
            continue
        first = next((l.lstrip("# ").strip() for l in e["code"].splitlines() if l.startswith("#") and len(l) > 2),
                     "")
        out.append(dict(name=e["name"], summary=first))
    return out


def read_example(name: str, language: str = "en") -> dict:
    """Code of an example of the user manual and the output it produces."""
    lang = _check_lang(language)
    stem = name.strip().removesuffix(".py")
    for e in _corpus()["examples"]:
        if e["lang"] == lang and e["name"] == stem + ".py":
            return dict(name=e["name"], code=e["code"], output=e["output"])
    raise KeyError(f"unknown example {name!r}; use list_examples")


def describe_api(name: str) -> dict:
    """Signature and docstring of a class, method or function of the installed module.

    ``name`` may be ``Model``, ``Model.add_beam``, ``keta.modal_combination`` or a bare method name
    (``add_beam``), which is looked up in every class. For a class, the list of its public members is
    returned too.
    """
    import keta

    path = name.strip().removeprefix("keta.")
    found = []
    obj = keta
    try:
        for part in path.split("."):
            obj = getattr(obj, part)
        found.append((path, obj))
    except AttributeError:
        if "." not in path:
            for cname, cls in inspect.getmembers(keta, inspect.isclass):
                if hasattr(cls, path):
                    found.append((f"{cname}.{path}", getattr(cls, path)))
    if not found:
        raise KeyError(f"no attribute {name!r} in keta")
    out = []
    for qual, o in found[:10]:
        item = dict(name=f"keta.{qual}", doc=inspect.getdoc(o) or "")
        if inspect.isclass(o):
            item["members"] = sorted(n for n in dir(o) if not n.startswith("_"))
        out.append(item)
    return dict(version=keta.__version__, matches=out)


INSTRUCTIONS = (
    "keta is a finite-element kernel for bar structures and bridges, written in C++20 and driven from Python "
    "(import keta). SI units (m, N, Pa, kg, rad), global Z axis upwards. Use search_documentation to find the "
    "relevant sections of the user manual (every command, its parameters and an executed example) and of the "
    "calculation principles, read_section to read them, list_examples / read_example for complete scripts with "
    "their output, and describe_api for the exact signature of a method in the installed version. Documentation "
    "exists in English ('en') and French ('fr')."
)


def build_server():
    """The MCP server (mcp >= 2: MCPServer; mcp 1.x: FastMCP)."""
    try:
        from mcp.server.mcpserver import MCPServer as Server
    except ImportError:
        try:
            from mcp.server.fastmcp import FastMCP as Server
        except ImportError as e:
            raise SystemExit('the MCP server needs the "mcp" package: pip install "keta[mcp]"') from e
    server = Server("keta", instructions=INSTRUCTIONS)
    for f in (search_documentation, read_section, list_contents, list_examples, read_example, describe_api):
        server.tool()(f)
    return server


def main() -> None:
    build_server().run("stdio")


if __name__ == "__main__":
    main()
