"""G-G(1): Python validation producers keep a Web-visible failure meaning.

Every ``raise ValidationError(...)`` in the package either carries a producer
``diagnostic=`` or has a message the Web adapter classifies. The remaining
sites are counted by ``UNPROVEN_BASELINE``, and the adapter's message tables by
``NATIVE_TABLE_BASELINE``. Both are registered design-rule allowlists
(``tools/web-design-rule-guards.json``, R6): an observed count may not exceed
its baseline, lowering a baseline is optional, and raising one is an
authority-only change.
"""

from __future__ import annotations

import ast
import itertools
import re
import subprocess
import warnings
from collections import Counter
from functools import lru_cache
from pathlib import Path

from gbdraw.exceptions import ValidationError
from gbdraw.web_support import error_adapter
from gbdraw.web_support.error_adapter import serialize_web_error

REPO_ROOT = Path(__file__).resolve().parents[1]
PYTHON_HELPERS = REPO_ROOT / "gbdraw" / "web" / "js" / "app" / "python-helpers.js"
WEB_WORDING = REPO_ROOT / "gbdraw" / "web" / "js" / "utils" / "error-normalization.js"
UNCLASSIFIED = {"UNKNOWN", "VALIDATION_UNCLASSIFIED"}
# Placeholder values for interpolations: empty, numbers, private text, a #index.
FILLS = ("", "0", "1", "PRIVATE", "#1")

# Sizes of the adapter's message-classification tables (R6, kind
# `py-count-map`). A producer migrated to diagnostic= removes its row.
NATIVE_TABLE_BASELINE = {"_EXACT": 24, "_TEMPLATES": 20, "_CONSTRAINTS": 22}

# The ValidationError sites of the whole package that neither carry
# diagnostic= nor classify (R6, kind `py-count`): typed API and request
# contract checks, Session decoding, comparison invariants, depth and
# annotation input checks, CLI adapter and table reader messages, and internal
# renderer state. One total, so a site that moves to another module does not
# change it.
UNPROVEN_BASELINE = 1180


def _package_sources() -> list[Path]:
    listed = subprocess.run(
        ["git", "ls-files", "gbdraw/*.py"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=True,
    ).stdout.split()
    return [REPO_ROOT / path for path in listed]


_HOLE = object()


def _parts(node: ast.AST) -> list | None:
    """Return a message as literal parts and interpolation holes."""

    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return [node.value]
    if isinstance(node, ast.JoinedStr):
        return [value.value if isinstance(value, ast.Constant) else _HOLE for value in node.values]
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Add):
        left, right = _parts(node.left), _parts(node.right)
        return None if left is None or right is None else left + right
    if (
        isinstance(node, ast.Call)
        and isinstance(node.func, ast.Attribute)
        and node.func.attr == "format"
    ):
        template = _parts(node.func.value)
        if template is None or _HOLE in template:
            return None
        pieces = re.split(r"\{[^{}]*\}", "".join(template))
        return [item for piece in pieces for item in (piece, _HOLE)][:-1]
    return None


def _renderings(node: ast.AST) -> list[str]:
    """Fill each hole with each placeholder (all combinations for <= 3 holes)."""

    parts = _parts(node)
    if parts is None:
        return []
    holes = parts.count(_HOLE)
    choices = itertools.product(FILLS, repeat=holes) if holes <= 3 else ((fill,) * holes for fill in FILLS)
    messages = []
    for choice in choices:
        values = iter(choice)
        messages.append("".join(next(values) if part is _HOLE else part for part in parts))
    return messages


def _call_name(call: ast.Call) -> str | None:
    if isinstance(call.func, ast.Name):
        return call.func.id
    if isinstance(call.func, ast.Attribute):
        return call.func.attr
    return None


@lru_cache(maxsize=None)
def _trees() -> tuple[tuple[str, ast.Module], ...]:
    trees = [
        (str(path.relative_to(REPO_ROOT)), ast.parse(path.read_text(encoding="utf-8")))
        for path in _package_sources()
    ]
    helpers = PYTHON_HELPERS.read_text(encoding="utf-8")
    source = helpers[helpers.index("`") + 1 : helpers.rindex("`;")]
    trees.append((str(PYTHON_HELPERS.relative_to(REPO_ROOT)), ast.parse(source)))
    return tuple(trees)


def _validation_raise_sites():
    # The Worker's embedded Python helpers are producers too.
    for path, tree in _trees():
        for node in ast.walk(tree):
            if (
                isinstance(node, ast.Raise)
                and isinstance(node.exc, ast.Call)
                and _call_name(node.exc) == "ValidationError"
            ):
                yield path, node.lineno, node.exc


def _classifies(message: str) -> bool:
    payload = serialize_web_error(
        ValidationError(message), operation="generate", stage="request-validation"
    )
    return payload["code"] not in UNCLASSIFIED


def _site_status(call: ast.Call) -> str:
    # ``diagnostic=None`` states no meaning, so it cannot satisfy the ratchet.
    if any(
        keyword.arg == "diagnostic"
        and not (isinstance(keyword.value, ast.Constant) and keyword.value.value is None)
        for keyword in call.keywords
    ):
        return "diagnostic"
    if call.args:
        if any(_classifies(message) for message in _renderings(call.args[0])):
            return "classified"
    return "unproven"


def _unproven_by_file() -> Counter:
    counts: Counter = Counter()
    for path, _line, call in _validation_raise_sites():
        if _site_status(call) == "unproven":
            counts[path] += 1
    return counts


def test_unproven_validation_producers_only_shrink():
    actual = _unproven_by_file()
    total = sum(actual.values())
    assert total <= UNPROVEN_BASELINE, (
        f"{total} ValidationError sites classify as UNKNOWN or VALIDATION_UNCLASSIFIED, "
        f"above UNPROVEN_BASELINE {UNPROVEN_BASELINE}. Raise the new ones with diagnostic= "
        f"instead. Sites per file (compare with dev): {dict(sorted(actual.items()))}"
    )
    if total < UNPROVEN_BASELINE:
        warnings.warn(f"UNPROVEN_BASELINE may be lowered: {UNPROVEN_BASELINE} -> {total}")


def _constant_items(node: ast.Dict) -> dict[str, object]:
    return {
        key.value: value.value
        for key, value in zip(node.keys, node.values)
        if isinstance(key, ast.Constant) and isinstance(value, ast.Constant)
    }


def _literal_diagnostics():
    for path, _line, call in _validation_raise_sites():
        for keyword in call.keywords:
            if keyword.arg == "diagnostic" and isinstance(keyword.value, ast.Dict):
                yield path, _constant_items(keyword.value)
    # Producer helpers that build diagnostics return dict literals with a code.
    for path, tree in _trees():
        if path.endswith("error_adapter.py") or path.endswith(".js"):
            continue
        for node in ast.walk(tree):
            if not isinstance(node, (ast.FunctionDef, ast.Assign)):
                continue
            name = node.name if isinstance(node, ast.FunctionDef) else ast.unparse(node.targets[0])
            if "diagnostic" not in name and "invalid" not in name.lower() and name != "_UNREADABLE":
                continue
            for child in ast.walk(node):
                if isinstance(child, ast.Dict) and any(
                    isinstance(key, ast.Constant) and key.value == "code" for key in child.keys
                ):
                    yield path, _constant_items(child)


def test_literal_diagnostics_use_the_published_vocabulary():
    seen = list(_literal_diagnostics())
    assert len(seen) >= 10
    for path, diagnostic in seen:
        assert diagnostic.get("code") in error_adapter.DIAGNOSTIC_CODES, (path, diagnostic)
        if "reason" in diagnostic:
            assert diagnostic["reason"] in error_adapter.DIAGNOSTIC_REASONS, (path, diagnostic)
        if "field" in diagnostic:
            assert diagnostic["field"] in error_adapter.FIELDS, (path, diagnostic)
        payload = serialize_web_error(
            ValidationError("PRIVATE", diagnostic=diagnostic),
            operation="generate",
            stage="request-validation",
        )
        assert payload["code"] == diagnostic["code"]
        assert "PRIVATE" not in str(payload)


def _web_block(start: str) -> str:
    text = WEB_WORDING.read_text(encoding="utf-8")
    begin = text.index(start)
    return text[begin : text.index("});", begin)]


def test_python_diagnostic_vocabulary_matches_the_web_wording_owner():
    web_codes = set(re.findall(r"\b([A-Z][A-Z0-9_]+): \[", _web_block("const DEFINITIONS")))
    web_reasons = set(re.findall(r"\b([A-Z][A-Z0-9_]+): '", _web_block("const REASONS")))
    text = WEB_WORDING.read_text(encoding="utf-8")
    context_block = text[text.index("const contextFor") : text.index("export const normalizeUserFacingError")]
    web_context_keys = (
        set(re.findall(r"'([A-Za-z]+)'", context_block))
        | set(re.findall(r"value\.([A-Za-z]+)", context_block))
    )
    web_fields = set(re.findall(r"[a-zA-Z_]+", text[text.index("const FIELDS") : text.index("const REASONS")]))
    assert error_adapter.DIAGNOSTIC_CODES <= web_codes
    assert error_adapter.DIAGNOSTIC_REASONS <= web_reasons
    assert error_adapter._DIAGNOSTIC_INTEGER_KEYS | {"configPath"} <= web_context_keys
    assert error_adapter.FIELDS <= web_fields
    # The engine-stage fallback (B9) and its bounded exception-class context.
    assert "RENDER_FAILED" in web_codes
    assert "exceptionType" in web_context_keys
    web_exception_types = set(re.findall(r"[A-Za-z]+", text[text.index("const EXCEPTION_TYPES") : text.index("const ORDINAL_LABELS")]))
    assert error_adapter.EXCEPTION_TYPE_NAMES <= web_exception_types


def test_native_message_tables_only_shrink():
    actual = {name: len(getattr(error_adapter, name)) for name in ("_EXACT", "_TEMPLATES", "_CONSTRAINTS")}
    grown = {name: size for name, size in actual.items() if size > NATIVE_TABLE_BASELINE.get(name, 0)}
    assert not grown, (
        "The adapter's English classification tables may only shrink; raise new "
        f"failures with diagnostic= (above NATIVE_TABLE_BASELINE: {grown})"
    )
    lowerable = {name: size for name, size in actual.items() if size < NATIVE_TABLE_BASELINE.get(name, 0)}
    if lowerable:
        warnings.warn(f"NATIVE_TABLE_BASELINE may be lowered: {lowerable}")


def _producer_messages() -> set[str]:
    messages: set[str] = set()
    for _path, tree in _trees():
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call) or not node.args:
                continue
            name = _call_name(node)
            if not name or not (name.endswith("Error") or name == "Exception"):
                continue
            messages.update(_renderings(node.args[0]))
    return messages


def test_native_classification_rows_all_have_a_producer():
    messages = _producer_messages()
    dead_exact = [message for message in error_adapter._EXACT if message not in messages]
    dead_templates = [
        template
        for template, _code, _reason in error_adapter._TEMPLATES
        if not any(re.fullmatch(template, message) for message in messages)
    ]
    dead_constraints = [
        clause
        for clause in error_adapter._CONSTRAINTS
        if not any(
            re.fullmatch(r"([\s\S]+?) " + re.escape(clause) + r"(?:\.| or null)?", message)
            for message in messages
        )
    ]
    assert not dead_exact, dead_exact
    assert not dead_templates, dead_templates
    assert not dead_constraints, dead_constraints
