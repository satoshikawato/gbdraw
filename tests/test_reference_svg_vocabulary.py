"""SVG 1.1 vocabulary of the reference outputs.

Every SVG in ``tests/reference_outputs/`` is well-formed XML with an SVG ``svg``
root. Each element is an SVG 1.1 element in the SVG namespace that its parent's
content model allows, and each attribute, namespace declarations included, is
one the SVG 1.1 DTD declares for that element or matches an entry of
``SVG11_EXTENSIONS``. The DTD is the W3C file, unchanged and with its notice:
``tests/fixtures/svg11-flat-20110816.dtd``.

``SVG11_EXTENSIONS`` is a registered design-rule allowlist
(``tools/web-design-rule-guards.json``, SVG11-VOCABULARY, kind ``py-set``):
adding an entry is an authority-only change that needs the Owner's approval,
and removing one is free.
"""

from __future__ import annotations

import fnmatch
import hashlib
import re
import xml.etree.ElementTree as ET
from pathlib import Path

TESTS_DIR = Path(__file__).resolve().parent
REFERENCE_DIR = TESTS_DIR / "reference_outputs"
SVG11_DTD = TESTS_DIR / "fixtures" / "svg11-flat-20110816.dtd"
# https://www.w3.org/Graphics/SVG/1.1/DTD/svg11-flat-20110816.dtd
SVG11_DTD_SHA256 = "7c22266d65655a7d359e9e79a107d156e64b29a87491313be2092bd66cab3b9b"
SVG_NAMESPACE = "http://www.w3.org/2000/svg"
# The prefix the DTD gives the attributes of each namespace.
ATTRIBUTE_PREFIXES = {
    "http://www.w3.org/1999/xlink": "xlink:",
    "http://www.w3.org/XML/1998/namespace": "xml:",
}

# Names (fnmatch patterns) outside SVG 1.1 that gbdraw writes on purpose.
SVG11_EXTENSIONS: frozenset[str] = frozenset(
    {
        # Custom data attributes (SVG 2): the interactive metadata that the Web
        # editor reads to find a feature, record, label or match by its identity.
        "data-[a-z]*",
        # WAI-ARIA (SVG 2): label halos, leader lines and backings are hidden
        # from assistive technology (gbdraw/render/drawers/label_backing.py).
        "aria-hidden",
        # svgwrite declares the XML Events namespace on the root; nothing uses it.
        "xmlns:ev",
    }
)

_DTD_TOKEN = re.compile(
    r"<!--.*?-->"
    r"|<!\[\s*%(?P<condition>[\w.:-]+);\s*\["
    r"|(?P<end>\]\]>)"
    r"|<!ENTITY\s+%\s+(?P<entity>[\w.:-]+)\s+(?P<value>\"[^\"]*\"|'[^']*')\s*>"
    r"|<!(?P<kind>ELEMENT|ATTLIST)\s+(?P<declaration>.*?)>",
    re.S,
)
_ATTRIBUTE_DEFINITION = re.compile(
    r"([\w:.-]+)\s+(?:\([^)]*\)|\S+)\s+(?:#REQUIRED|#IMPLIED|(?:#FIXED\s+)?(?:\"[^\"]*\"|'[^']*'))"
)


def _read_svg11_dtd() -> tuple[dict[str, set[str]], dict[str, set[str]]]:
    """Each SVG 1.1 element's allowed child elements and attributes.

    As XML 1.0 reads a DTD: the first declaration of a parameter entity wins,
    entities expand where they are used, and a conditional section counts only
    when its keyword expands to INCLUDE.
    """
    entities: dict[str, str] = {}
    included: list[bool] = []
    children: dict[str, set[str]] = {}
    attributes: dict[str, set[str]] = {}

    def expand(text: str) -> str:
        return re.sub(r"%([\w.:-]+);", lambda match: expand(entities[match[1]]), text)

    for token in _DTD_TOKEN.finditer(SVG11_DTD.read_text(encoding="ascii")):
        if token["condition"]:
            included.append(
                all(included) and expand(f"%{token['condition']};").strip() == "INCLUDE"
            )
        elif token["end"]:
            included.pop()
        elif token["entity"] and all(included):
            entities.setdefault(token["entity"], token["value"][1:-1])
        elif token["kind"] and all(included):
            name, body = expand(token["declaration"]).split(None, 1)
            if token["kind"] == "ELEMENT":
                children[name] = set(re.findall(r"[\w:.-]+", body))
            else:
                attributes.setdefault(name, set()).update(
                    _ATTRIBUTE_DEFINITION.findall(body)
                )
    return children, attributes


def _vocabulary_outside_svg11(
    path: Path, children: dict[str, set[str]], attributes: dict[str, set[str]]
) -> tuple[set[str], set[str]]:
    """The ``parent > element`` and ``element @attribute`` uses that neither SVG 1.1
    nor ``SVG11_EXTENSIONS`` allows, and the extension entries that the file uses."""
    found: set[str] = set()
    used: set[str] = set()
    parents: list[str] = []
    declarations: list[str] = []
    for event, item in ET.iterparse(path, events=("start-ns", "start", "end")):
        if event == "start-ns":
            prefix, _uri = item
            declarations.append(f"xmlns:{prefix}" if prefix else "xmlns")
        elif event == "end":
            parents.pop()
        else:
            tag = item.tag.removeprefix(f"{{{SVG_NAMESPACE}}}")
            if tag not in (children.get(parents[-1], set()) if parents else {"svg"}):
                found.add(f"{parents[-1] if parents else '(root)'} > {tag}")
            for name in [*declarations, *item.attrib]:
                uri, brace, local = name[1:].partition("}")
                qualified = (
                    ATTRIBUTE_PREFIXES.get(uri, f"{{{uri}}}") + local if brace else name
                )
                if qualified in attributes.get(tag, set()):
                    continue
                extensions = {
                    pattern
                    for pattern in SVG11_EXTENSIONS
                    if fnmatch.fnmatchcase(qualified, pattern)
                }
                used |= extensions
                if not extensions:
                    found.add(f"{tag} @{qualified}")
            declarations = []
            parents.append(tag)
    return found, used


def test_the_dtd_is_the_w3c_file() -> None:
    digest = hashlib.sha256(SVG11_DTD.read_bytes()).hexdigest()
    assert digest == SVG11_DTD_SHA256, (
        "the SVG 1.1 DTD differs from the W3C file; list additions in SVG11_EXTENSIONS"
    )


def test_reference_outputs_use_svg11_vocabulary_and_listed_extensions() -> None:
    children, attributes = _read_svg11_dtd()
    paths = sorted(REFERENCE_DIR.glob("*.svg"))
    assert paths, f"no SVG in {REFERENCE_DIR}"
    outside: dict[str, str] = {}
    used: set[str] = set()
    for path in paths:
        found, extensions = _vocabulary_outside_svg11(path, children, attributes)
        used |= extensions
        for use in found:
            outside.setdefault(use, path.name)
    assert not outside, (
        "reference SVGs use vocabulary that is neither SVG 1.1 nor in SVG11_EXTENSIONS "
        "(first file in brackets):\n"
        + "\n".join(f"{use} [{name}]" for use, name in sorted(outside.items()))
    )
    assert used == SVG11_EXTENSIONS, (
        f"no reference SVG uses {sorted(SVG11_EXTENSIONS - used)}; remove it from SVG11_EXTENSIONS"
    )
