"""BLAST comparison tables: one reader for every entry point and record-ID binding.

CO-05 and N-03 (PD-OI-077): every reader keeps the first 12 typed columns and a
missing or malformed table fails instead of being skipped. CO-06 (PD-OI-074):
table IDs that contradict the endpoint records fail; unknown IDs warn.
"""

from __future__ import annotations

import ast
import json
import logging
import math
import random
import re
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import gbdraw
from gbdraw.analysis.conservation import load_conservation_sources
from gbdraw.api import LinearComparison, LinearDiagramOptions
from gbdraw.api.config import apply_config_overrides
from gbdraw.api.diagram import assemble_linear_diagram_from_records
from gbdraw.api.record_planning import resolve_linear_options
from gbdraw.exceptions import ComparisonIdentityError, ValidationError
from gbdraw.diagrams.linear.assemble import _read_comparison_files
from gbdraw.io.comparisons import COMPARISON_COLUMNS
from gbdraw.mode_profiles import ComparisonThresholds
from gbdraw.session_request_codec import (
    CANONICAL_REQUEST_SCHEMA,
    _decode_comparisons,
    decode_canonical_request,
)
from gbdraw.web_support.error_adapter import serialize_web_error
from gbdraw.api.assembly_inputs import LinearAssemblyInputs

_VECTORS = json.loads(
    (Path(__file__).parent / "fixtures" / "comparison_outfmt6_table_cases.json").read_text(
        encoding="utf-8"
    )
)
_PERMISSIVE = SimpleNamespace(evalue=10.0, bitscore=0.0, identity=0.0, alignment_length=0)
_THRESHOLDS = ComparisonThresholds(evalue=10.0, bitscore=0.0, identity=0.0, alignment_length=0)


def _table_text(lines: list[object]) -> str:
    return "".join(
        (line if isinstance(line, str) else "\t".join(line)) + "\n" for line in lines
    )


def _write_table(path: Path, lines: list[object]) -> Path:
    path.write_text(_table_text(lines), encoding="utf-8")
    return path


def _two_records() -> list[SeqRecord]:
    records = [SeqRecord(Seq("A" * 3000), id=record_id) for record_id in ("SA", "SB")]
    for record in records:
        record.annotations["molecule_type"] = "DNA"
    return records


def _read_with_shared_reader(path: Path) -> pd.DataFrame:
    from gbdraw.io.comparisons import read_comparison_table

    return read_comparison_table(path)


def _read_with_blast_option(path: Path) -> pd.DataFrame:
    return _read_comparison_files([str(path)], _THRESHOLDS)[0]


def _read_with_canonical_codec(path: Path) -> pd.DataFrame:
    decoded = _decode_comparisons(
        [
            {
                "kind": "nucleotideBlast",
                "resourceId": "comparison-nucleotide-1",
                "queryRecordIndex": 0,
                "subjectRecordIndex": 1,
            }
        ],
        mode="linear",
        schema=CANONICAL_REQUEST_SCHEMA,
        resource_paths={"comparison-nucleotide-1": path},
        resource_names={},
    )
    return decoded["linear_comparisons"][0].matches


def _read_with_comparisons_table(path: Path) -> pd.DataFrame:
    manifest = path.parent / "comparisons.tsv"
    manifest.write_text(f"blast\tquery\tsubject\n{path.name}\t#1\t#2\n", encoding="utf-8")
    options = resolve_linear_options(
        LinearDiagramOptions(comparison_table_file=str(manifest)),
        records=_two_records(),
        layout=None,
    )
    return options.linear_comparisons[0].matches


def _read_with_conservation(path: Path) -> pd.DataFrame:
    result = load_conservation_sources(blast_config=_PERMISSIVE, conservation_files=[str(path)])
    if result.skipped_sources:
        # Similarity rings keep their source index and report the reader error.
        raise ValidationError(str(result.skipped_sources[0].skip_reason))
    return result.sources[0].dataframe


_ENTRY_POINTS = {
    "shared-reader": _read_with_shared_reader,
    "-b (_read_comparison_files)": _read_with_blast_option,
    "canonical codec (Web)": _read_with_canonical_codec,
    "--comparisons_table": _read_with_comparisons_table,
    "conservation ring": _read_with_conservation,
}
# A similarity ring skips an unreadable source and reports its text, so only the
# entry points that raise the reader error carry its Web diagnostic.
_DIAGNOSTIC_ENTRY_POINTS = frozenset(_ENTRY_POINTS) - {"conservation ring"}


def _web(error: BaseException) -> dict:
    return serialize_web_error(error, operation="generate", stage="render")


@pytest.mark.parametrize("entry_point", list(_ENTRY_POINTS))
@pytest.mark.parametrize(
    "case", _VECTORS["cases"], ids=[case["name"] for case in _VECTORS["cases"]]
)
def test_every_comparison_reader_follows_the_shared_outfmt_vectors(
    tmp_path: Path, case: dict, entry_point: str
) -> None:
    path = _write_table(tmp_path / "pair.tsv", case["lines"])
    read = _ENTRY_POINTS[entry_point]
    expect = case["expect"]
    if "error" in expect:
        with pytest.raises(ValidationError, match=expect["error"]) as caught:
            read(path)
        if entry_point in _DIAGNOSTIC_ENTRY_POINTS:
            payload = _web(caught.value)
            diagnostic = expect["diagnostic"]
            if entry_point == "canonical codec (Web)":
                # The Web request names the record pair of the table (CI-04).
                diagnostic = {**diagnostic, "context": {
                    **diagnostic["context"], "queryRecordIndex": 0, "subjectRecordIndex": 1,
                }}
            assert {key: payload[key] for key in ("code", "context")} == diagnostic
        return

    frame = read(path)

    assert list(frame.columns[: len(COMPARISON_COLUMNS)]) == list(COMPARISON_COLUMNS)
    assert len(frame) == expect["rowCount"]
    assert isinstance(frame.index, pd.RangeIndex)
    if not expect["rowCount"]:
        return
    first = frame.loc[0, list(COMPARISON_COLUMNS)].tolist()
    assert first[:2] == expect["firstRow"][:2]
    for actual, expected in zip(first[2:], expect["firstRow"][2:], strict=True):
        assert math.isclose(float(actual), float(expected), rel_tol=1e-12)
    for column in ("alignment_length", "qstart", "qend", "sstart", "send"):
        assert pd.api.types.is_integer_dtype(frame[column])


def test_extra_columns_are_reported_in_an_info_log(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    case = next(case for case in _VECTORS["cases"] if case["name"] == "14 columns (std qlen slen)")
    path = _write_table(tmp_path / "pair.tsv", case["lines"])

    with caplog.at_level(logging.INFO, logger="gbdraw.io.comparisons"):
        _read_with_blast_option(path)

    assert "ignoring 2 extra column(s)" in caplog.text


def test_normalized_dataframe_uses_first_twelve_positional_columns() -> None:
    from gbdraw.io.comparisons import normalize_comparison_dataframe

    frame = pd.DataFrame(
        [["SA", "SB", "95.0", "200", "0", "0", "101", "300", "1101", "1300", "1e-50", "300", "3000"]]
    )

    normalized = normalize_comparison_dataframe(frame)

    assert normalized.loc[0, ["query", "identity", "qstart", "evalue"]].tolist() == [
        "SA",
        95.0,
        101,
        pytest.approx(1e-50),
    ]
    with pytest.raises(ValidationError, match=r"row 1, column 7 \(qstart\): 'x' is not an integer"):
        normalize_comparison_dataframe(frame.replace({"101": "x"}))
    boolean = frame.copy()
    boolean[2] = True
    with pytest.raises(ValidationError, match=r"row 1, column 3 \(identity\): 'True' is not a finite number"):
        normalize_comparison_dataframe(boolean)


@pytest.mark.parametrize(
    ("value", "context"),
    [
        (
            pd.DataFrame([["SA", "SB", "95.0", "200", "0", "0", "x", "300", "1101", "1300", "1e-50", "300"]]),
            {"reason": "INTEGER", "row": 1, "column": 7},
        ),
        (pd.DataFrame([["SA", "SB", "95.0"]]), {"reason": "FIELDS", "columnCount": 12}),
        ("PRIVATE_TABLE", {"field": "comparison"}),
    ],
    ids=["invalid value", "too few columns", "not a DataFrame"],
)
def test_normalized_dataframe_failures_carry_a_comparison_diagnostic(value, context) -> None:
    from gbdraw.io.comparisons import normalize_comparison_dataframe

    with pytest.raises(ValidationError) as caught:
        normalize_comparison_dataframe(value)

    payload = _web(caught.value)
    assert payload["code"] == "COMPARISON_INPUT"
    assert payload["context"] == context
    assert "PRIVATE" not in json.dumps(payload)


def test_the_b_option_rejects_a_missing_file_instead_of_shifting_later_tables(
    tmp_path: Path,
) -> None:
    valid = _write_table(tmp_path / "SB_SC.tsv", _VECTORS["cases"][0]["lines"])
    missing = tmp_path / "SA_SB.tsv"

    with pytest.raises(ValidationError, match=r"SA_SB\.tsv: the comparison file does not exist") as caught:
        _read_comparison_files([str(missing), str(valid)], _THRESHOLDS)
    assert _web(caught.value)["code"] == "INPUT_UNREADABLE"
    assert _web(caught.value)["context"] == {"field": "comparison"}


def test_the_b_option_rejects_a_malformed_file_instead_of_shifting_later_tables(
    tmp_path: Path,
) -> None:
    valid = _write_table(tmp_path / "SB_SC.tsv", _VECTORS["cases"][0]["lines"])
    malformed = _write_table(tmp_path / "SA_SB.tsv", [["SA", "SB", "95.0"]])

    with pytest.raises(ValidationError, match=r"SA_SB\.tsv: line 1: expected at least 12"):
        _read_comparison_files([str(malformed), str(valid)], _THRESHOLDS)


@pytest.fixture(scope="module")
def web_helpers() -> dict[str, object]:
    helpers_js = Path(gbdraw.__file__).parent / "web" / "js" / "app" / "python-helpers.js"
    source = helpers_js.read_text(encoding="utf-8")
    namespace: dict[str, object] = {}
    exec(source.split("`", 1)[1].rsplit("`", 1)[0], namespace)
    return namespace


def test_web_generate_reports_a_malformed_blast_table_line_as_a_comparison_diagnostic(
    tmp_path: Path, web_helpers: dict[str, object]
) -> None:
    """The Worker's Generate entry point returns the reader's diagnostic, not UNKNOWN."""

    from gbdraw.api import InMemoryRecordSource, LinearDiagramRequest, RecordInput
    from gbdraw.session_request_codec import encode_canonical_request

    request = LinearDiagramRequest(
        records=tuple(
            RecordInput(source=InMemoryRecordSource(record)) for record in _pair_records("SA", "SB")
        ),
        options=LinearDiagramOptions(linear_comparisons=(LinearComparison(0, 1, _hits("SA", "SB")),)),
    )
    encoded = encode_canonical_request(request)
    workspace = tmp_path / "gbdraw-web-render-1"
    resources = workspace / "resources"
    resources.mkdir(parents=True)
    (workspace / ".gbdraw-worker-render-workspace").touch()
    resource_paths = {}
    for index, resource in enumerate(encoded.resources, start=1):
        target = resources / f"{index:04d}.bin"
        target.write_bytes(resource.content)
        resource_paths[resource.resource_id] = str(target)
    # The browser uploads a BLAST table as a nucleotideBlast resource.
    comparison = encoded.payload["comparisons"][0]
    assert comparison["kind"] == "nucleotideBlast"  # a BLAST table is stored as one (OV-3512)
    case = next(
        case for case in _VECTORS["cases"] if case["name"] == "non-integer coordinate after comments"
    )
    Path(resource_paths[comparison["resourceId"]]).write_text(
        _table_text(case["lines"]), encoding="utf-8"
    )

    result = web_helpers["run_canonical_request_wrapper"](
        json.dumps(encoded.payload), json.dumps(resource_paths), str(workspace)
    )

    error = result["error"]
    assert (error["code"], error["operation"]) == ("COMPARISON_INPUT", "generate")
    assert error["context"] == {
        "reason": "INTEGER", "row": 3, "column": 7, "queryRecordIndex": 0, "subjectRecordIndex": 1,
    }
    assert "101.5" not in json.dumps(error)
    assert not workspace.exists()


def test_an_uploaded_gzip_blast_table_decodes_as_the_plain_one(tmp_path: Path) -> None:
    """The page sends a chosen BLAST table as chosen, so a gzip one reaches the
    Worker's Generate (and a Session's entry) as gzip bytes: the request decoder
    expands it as the CLI reads the file."""

    import gzip as gzip_module

    from gbdraw.api import InMemoryRecordSource, LinearDiagramRequest, RecordInput
    from gbdraw.session_request_codec import decode_canonical_request, encode_canonical_request

    request = LinearDiagramRequest(
        records=tuple(
            RecordInput(source=InMemoryRecordSource(record)) for record in _pair_records("SA", "SB")
        ),
        options=LinearDiagramOptions(linear_comparisons=(LinearComparison(0, 1, _hits("SA", "SB")),)),
    )
    encoded = encode_canonical_request(request)
    table = _table_text([[str(cell) for cell in _hits("SA", "SB").iloc[0]]]).encode("utf-8")
    decoded = []
    for name, content in (("plain", table), ("gzip", gzip_module.compress(table))):
        resources = tmp_path / name
        resources.mkdir()
        resource_paths = {}
        for index, resource in enumerate(encoded.resources, start=1):
            target = resources / f"{index:04d}.bin"
            target.write_bytes(resource.content)
            resource_paths[resource.resource_id] = str(target)
        payload = json.loads(json.dumps(encoded.payload))
        comparison = payload["comparisons"][0]
        assert comparison["kind"] == "nucleotideBlast"  # a BLAST table is stored as one (OV-3512)
        Path(resource_paths[comparison["resourceId"]]).write_bytes(content)
        decoded.append(decode_canonical_request(
            payload, resource_paths=resource_paths, output_directory=str(tmp_path / f"{name}-out")
        ).options.linear_comparisons[0])
    assert decoded[1].matches.equals(decoded[0].matches)
    assert decoded[1].search_frame_text == decoded[0].search_frame_text == table.decode("utf-8")


# A mapping of a table's columns onto COMPARISON_COLUMNS, as each language spells it. A reader takes positions or
# looks up a foreign header; it is seen in both languages as a loop, comprehension or callback over the columns (or a
# slice of them) that
# - indexes by each entry's position or keeps the position (`fields[index]`, `{column: index}`, `new Map(... [c, i])`),
# - or looks each entry up in a header (`header.index(column)`, `columns.indexOf(column)`, `index.get(column)`,
#   `get_loc`, `findIndex`);
# and as the columns subscripted by a position (`COMPARISON_COLUMNS[index]`), a field found by a column's position
# (`fields[COMPARISON_COLUMNS.index("evalue")]`), and, in Python, pandas column names (`names=`, `.columns =`,
# `set_axis(`, `rename(columns=...)`) and `dict(enumerate(COMPARISON_COLUMNS))`. An imported alias or a namespace
# (`schema.COMPARISON_COLUMNS`) counts as the constant. A consumer is not a mapping: a value subscripted by a
# canonical name (`row[column]`, `LABELS[column]`, a TSV writer), joining, membership, length, a fixed entry or slice,
# `zip(COMPARISON_COLUMNS, row)` and a frame built with `columns=COMPARISON_COLUMNS` (consumers of values already in
# that order; eleven on dev outside gbdraw/io/comparisons.py).
_PY_LOOKUP_METHODS = frozenset({"index", "get", "get_loc", "get_indexer"})
_JS_LOOKUP_METHODS = r"(?:indexOf|lastIndexOf|findIndex|findLastIndex|get)"
_JS_IMPORT_ALIAS = re.compile(r"\bCOMPARISON_COLUMNS\s+as\s+([\w$]+)")
_JS_IDENTIFIER = re.compile(r"[A-Za-z_$][\w$]*")
_JS_STATEMENT_START = re.compile(r"^[A-Za-z_$]", re.MULTILINE)
_JS_FUNCTION_DECLARATION = re.compile(
    r"(?:export\s+)?(?:default\s+)?(?:async\s+)?function\s*\*?\s*([\w$]+)"
    r"|(?:export\s+)?(?:const|let|var)\s+([\w$]+)\s*=\s*(?:async\s*)?(?:function\b|\([^)]*\)\s*=>|[\w$]+\s*=>)"
    r"|(?:export\s+)?(?:default\s+)?class\s+([\w$]+)"
)


def _python_column_mappings(tree: ast.Module) -> list[int]:
    """Line numbers of each Python mapping onto COMPARISON_COLUMNS (or an imported alias of it)."""

    names = {"COMPARISON_COLUMNS"} | {
        alias.asname
        for node in ast.walk(tree)
        if isinstance(node, ast.ImportFrom)
        for alias in node.names
        if alias.name == "COMPARISON_COLUMNS" and alias.asname
    }

    def columns(node: ast.AST | None) -> bool:
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id in {"list", "tuple"}:
            return bool(node.args) and columns(node.args[0])
        if isinstance(node, ast.Subscript) and isinstance(node.slice, ast.Slice):
            return columns(node.value)
        if isinstance(node, (ast.List, ast.Tuple)):
            return any(isinstance(item, ast.Starred) and columns(item.value) for item in node.elts)
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Add):
            return columns(node.left) or columns(node.right)
        return (isinstance(node, ast.Name) and node.id in names) or (
            isinstance(node, ast.Attribute) and node.attr == "COMPARISON_COLUMNS"
        )

    def bound(node: ast.AST) -> set[str]:
        return {name.id for name in ast.walk(node) if isinstance(name, ast.Name)}

    def iteration(iterable: ast.AST, target: ast.AST) -> tuple[set[str], set[str]] | None:
        """(positions, entries) bound by a loop over the columns, or None for another loop."""

        call = iterable if isinstance(iterable, ast.Call) and isinstance(iterable.func, ast.Name) else None
        if call is not None and call.func.id == "enumerate" and call.args and columns(call.args[0]):
            if isinstance(target, ast.Tuple) and len(target.elts) == 2:
                return bound(target.elts[0]), bound(target.elts[1])
            return bound(target), set()
        if call is not None and call.func.id == "reversed" and call.args and columns(call.args[0]):
            return set(), bound(target)
        return (set(), bound(target)) if columns(iterable) else None

    def reads(nodes: list[ast.AST], positions: set[str], entries: set[str], kept: list[ast.AST]) -> bool:
        """Whether the code indexes by a position, looks an entry up, or keeps a position in `kept` values."""

        for node in (inner for outer in nodes for inner in ast.walk(outer)):
            if isinstance(node, ast.Subscript) and bound(node.slice) & positions:
                return True
            if (
                isinstance(node, ast.Call)
                and isinstance(node.func, ast.Attribute)
                and node.func.attr in _PY_LOOKUP_METHODS
                and any(bound(arg) & entries for arg in node.args)
            ):
                return True
        return any(bound(value) & positions for value in kept)

    def stored(body: list[ast.stmt]) -> list[ast.AST]:
        """The values a loop body stores under a key (`index_of[column] = position`)."""

        return [
            node.value
            for statement in body
            for node in ast.walk(statement)
            if isinstance(node, ast.Assign) and any(isinstance(target, ast.Subscript) for target in node.targets)
        ]

    lines = []
    for node in ast.walk(tree):
        found = False
        if isinstance(node, (ast.For, ast.AsyncFor)):
            loop = iteration(node.iter, node.target)
            found = loop is not None and reads(list(node.body), *loop, stored(node.body))
        elif isinstance(node, (ast.ListComp, ast.SetComp, ast.GeneratorExp, ast.DictComp)):
            results = [node.key, node.value] if isinstance(node, ast.DictComp) else [node.elt]
            for generator in node.generators:
                loop = iteration(generator.iter, generator.target)
                found = found or (loop is not None and reads(results, *loop, results))
        elif isinstance(node, ast.Subscript):
            index = node.slice
            fixed = isinstance(index, (ast.Constant, ast.Slice)) or (
                isinstance(index, ast.UnaryOp) and isinstance(index.operand, ast.Constant)
            )
            found = (columns(node.value) and not fixed) or any(
                isinstance(call, ast.Call)
                and isinstance(call.func, ast.Attribute)
                and call.func.attr == "index"
                and columns(call.func.value)
                for call in ast.walk(index)
            )
        elif isinstance(node, ast.Call):
            function = node.func.id if isinstance(node.func, ast.Name) else getattr(node.func, "attr", "")
            keywords = {keyword.arg: keyword.value for keyword in node.keywords}
            first = node.args[0] if node.args else None
            found = (
                (
                    function == "dict"
                    and isinstance(first, ast.Call)
                    and isinstance(first.func, ast.Name)
                    and first.func.id == "enumerate"
                    and bool(first.args)
                    and columns(first.args[0])
                )
                or (function == "set_axis" and (columns(first) or columns(keywords.get("labels"))))
                or (
                    function == "rename"
                    and "columns" in keywords
                    and any(columns(inner) for inner in ast.walk(keywords["columns"]))
                )
                or columns(keywords.get("names"))
            )
        elif isinstance(node, ast.Assign):
            found = columns(node.value) and any(
                isinstance(target, ast.Attribute) and target.attr == "columns" for target in node.targets
            )
        if found:
            lines.append(node.lineno)
    return sorted(lines)


def _python_top_level_name(tree: ast.Module, line: int) -> str:
    """The top-level function (or `Class.method`) holding a line, or `<module>`."""

    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            if node.lineno <= line <= (node.end_lineno or node.lineno):
                if isinstance(node, ast.ClassDef):
                    for member in node.body:
                        if isinstance(member, (ast.FunctionDef, ast.AsyncFunctionDef)) and (
                            member.lineno <= line <= (member.end_lineno or member.lineno)
                        ):
                            return f"{node.name}.{member.name}"
                return node.name
    return "<module>"


def _js_literal_end(text: str, start: int) -> int:
    """The offset after a string or regular expression literal that opens at `start`."""

    quote, index, in_class = text[start], start + 1, False
    while index < len(text):
        char = text[index]
        if char == "\\":
            index += 2
            continue
        if char == "\n" and quote != "`":
            return index  # an unterminated literal ends at its line
        if quote == "/" and char in "[]":
            in_class = char == "["
        elif char == quote and not in_class:
            return index + 1
        index += 1
    return len(text)


def _js_without_comments(text: str) -> tuple[str, str]:
    """The text with comments blanked, and with comments and literals blanked; offsets and lines are kept.

    A string or regular expression literal is skipped, so a `/*` or `//` inside one (`'text/*'`) is not a comment.
    """

    code, skeleton = list(text), list(text)
    index, length, previous = 0, len(text), ""
    while index < length:
        char = text[index]
        # A slash after an operator, a bracket, `return` or `typeof`, or at the start opens a regular expression.
        regex = char == "/" and text[index + 1 : index + 2] not in {"/", "*"} and (
            previous in "(,=:[!&|?{};+-*%<>~^"
            or re.search(r"\b(?:return|typeof)\s*$", text[max(0, index - 12) : index]) is not None
        )
        if char in "'\"`" or regex:
            end = _js_literal_end(text, index)
            for blank in range(index + 1, end - 1):
                if text[blank] != "\n":
                    skeleton[blank] = " "
            index, previous = end, "a"
        elif text.startswith("//", index) or text.startswith("/*", index):
            end = text.find("\n" if text[index + 1] == "/" else "*/", index + 2)
            end = length if end < 0 else end + (0 if text[index + 1] == "/" else 2)
            for blank in range(index, end):
                if text[blank] != "\n":
                    code[blank] = skeleton[blank] = " "
            index = end
        else:
            previous = previous if char.isspace() else char
            index += 1
    return "".join(code), "".join(skeleton)


def _balanced(text: str, start: int, opening: str, closing: str) -> str:
    """The text from `start` (an opening bracket) to its closing bracket, or to the end."""

    depth = 0
    for index in range(start, len(text)):
        if text[index] == opening:
            depth += 1
        elif text[index] == closing:
            depth -= 1
            if depth == 0:
                return text[start : index + 1]
    return text[start:]


def _js_reads(body: str, positions: list[str], entries: list[str], *, keeps: bool) -> bool:
    """Whether code indexes by a position, looks an entry up, or (`keeps`) uses a position at all."""

    for name in positions:
        word = rf"(?<![\w$.]){re.escape(name)}(?![\w$])"
        if re.search(rf"\[\s*{word}\s*\]", body) or (keeps and re.search(word, body)):
            return True
    for match in re.finditer(rf"\.\s*{_JS_LOOKUP_METHODS}\s*\(", body):
        arguments = _balanced(body, match.end() - 1, "(", ")")
        if any(re.search(rf"(?<![\w$.]){re.escape(name)}(?![\w$])", arguments) for name in entries):
            return True
    return False


def _js_column_mappings(code: str) -> list[int]:
    """Offsets of each JavaScript mapping onto COMPARISON_COLUMNS (or an alias of it) in comment-free code."""

    names = ["COMPARISON_COLUMNS", *_JS_IMPORT_ALIAS.findall(code)]
    constant = "(?:" + "|".join(rf"(?<![\w$]){re.escape(name)}(?![\w$])" for name in names) + ")"
    columns = rf"{constant}(?:\s*\.\s*slice\s*\([^)]*\))?"
    offsets = []
    callbacks = (
        rf"(?:{columns}\s*\.\s*(map|flatMap|forEach|reduce|reduceRight)\s*"
        rf"|Array\s*\.\s*from\s*\(\s*{columns}\s*,)"
    )
    for match in re.finditer(callbacks, code):
        if match.group(1) and code[match.end() : match.end() + 1] != "(":
            continue
        # The callback from an opening parenthesis: the call's own, or one put before Array.from's second argument.
        if match.group(1):
            callback = _balanced(code, match.end(), "(", ")")
        else:
            callback = _balanced("(" + code[match.end() :], 0, "(", ")")
        head = re.match(r"\(\s*(?:async\s+)?(?:function\s*[\w$]*\s*)?(\([^)]*\)|[\w$]+)", callback)
        parameters = _JS_IDENTIFIER.findall(head.group(1)) if head else []
        if match.group(1) in {"reduce", "reduceRight"}:
            parameters = parameters[1:]
        body = callback[head.end() :] if head else callback
        if _js_reads(body, parameters[1:2], parameters[:1], keeps=True):
            offsets.append(match.start())
    loop = rf"for\s*\(\s*(?:const|let|var)\s+([^;]+?)\s+of\s+{columns}(\s*\.\s*entries\s*\(\s*\))?\s*\)"
    for match in re.finditer(loop, code):
        rest = code[match.end() :]
        brace = re.match(r"\s*\{", rest)
        body = _balanced(code, match.end() + brace.end() - 1, "{", "}") if brace else rest.split(";", 1)[0]
        bound = _JS_IDENTIFIER.findall(match.group(1))
        positions, entries = (bound[:1], bound[1:]) if match.group(2) else ([], bound)
        stores = any(
            re.search(rf"(?:\]\s*=(?!=)|\.\s*set\s*\()[^;]*(?<![\w$.]){re.escape(name)}(?![\w$])", body)
            for name in positions
        )
        if stores or _js_reads(body, positions, entries, keeps=False):
            offsets.append(match.start())
    offsets.extend(match.start() for match in re.finditer(rf"{constant}\s*\[\s*(?!\d)[^\]\s]", code))
    offsets.extend(match.start() for match in re.finditer(rf"\[\s*{constant}\s*\.\s*indexOf\s*\(", code))
    return sorted(set(offsets))


def _js_top_level_name(skeleton: str, offset: int) -> str:
    """The top-level function or class whose statement holds an offset, or `<module>`."""

    starts = [match.start() for match in _JS_STATEMENT_START.finditer(skeleton) if match.start() <= offset]
    if not starts:
        return "<module>"
    declaration = _JS_FUNCTION_DECLARATION.match(skeleton, starts[-1])
    return next((name for name in declaration.groups() if name), "<module>") if declaration else "<module>"


def _column_mappings(path: str, text: str) -> list[str]:
    """`path:function` of each mapping onto COMPARISON_COLUMNS in one source file (`<module>` outside a function)."""

    if path.endswith(".py"):
        tree = ast.parse(text)
        return [f"{path}:{_python_top_level_name(tree, line)}" for line in _python_column_mappings(tree)]
    code, skeleton = _js_without_comments(text)
    return [f"{path}:{_js_top_level_name(skeleton, offset)}" for offset in _js_column_mappings(code)]


def _package_column_mappings() -> dict[str, int]:
    package_root = Path(gbdraw.__file__).parent
    sources = [*package_root.rglob("*.py"), *(package_root / "web" / "js").rglob("*.js")]
    readers: dict[str, int] = {}
    for path in sorted(sources):
        text = path.read_text(encoding="utf-8")
        if "COMPARISON_COLUMNS" not in text:  # an alias or a namespace access names it too
            continue
        relative = path.relative_to(package_root.parent).as_posix()
        for site in _column_mappings(relative, text):
            readers[site] = readers.get(site, 0) + 1
    return readers


_OWNER_MODULE = "gbdraw/io/comparisons.py"


def _unregistered_column_mappings(readers: dict[str, int], registry: dict) -> dict[str, int]:
    """Mappings outside the owner module and outside every owner function of a twin of the owner module.

    A twin of the owner module is a row that is not pending and names a gbdraw/io/comparisons.py function as its
    Python owner; an owner of any other row is not exempt.
    """

    owners = {
        owner
        for row in registry["twins"]
        if row.get("pending") is not True
        and any(owner.startswith(f"{_OWNER_MODULE}:") for owner in row["pythonOwner"])
        for owner in (*row["pythonOwner"], *row["jsOwner"])
    }
    return {
        site: count
        for site, count in readers.items()
        if not site.startswith(f"{_OWNER_MODULE}:") and site not in owners
    }


def test_outfmt_table_columns_have_one_reader() -> None:
    """R9: gbdraw.io.comparisons maps table columns onto COMPARISON_COLUMNS; a second mapping is a registered twin.

    A mapping elsewhere lies in a function that a row of tools/twins.json pairing with gbdraw/io/comparisons.py (not
    pending) names as an owner.
    """

    readers = _package_column_mappings()
    registry = json.loads((Path(__file__).parents[1] / "tools" / "twins.json").read_text(encoding="utf-8"))
    assert _unregistered_column_mappings(readers, registry) == {}, (
        "a mapping onto COMPARISON_COLUMNS outside gbdraw/io/comparisons.py (file:function): call the owner, "
        "or register the twin in tools/twins.json"
    )
    # The file reader (its column names, and the position of an invalid cell's column) and the DataFrame normalizer.
    assert {site: count for site, count in readers.items() if site.startswith(f"{_OWNER_MODULE}:")} == {
        f"{_OWNER_MODULE}:read_comparison_table": 2,
        f"{_OWNER_MODULE}:normalize_comparison_dataframe": 1,
    }


def test_losatp_outfmt6_keeps_its_strict_grammar_errors() -> None:
    """The LOSATP parser leaves the outfmt 6 field grammar to gbdraw.io.comparisons, with the same errors."""

    from gbdraw.analysis.protein_colinearity import parse_losatp_outfmt6
    from gbdraw.exceptions import ParseError

    good = "q1\ts1\t91.5\t42\t1\t0\t1\t42\t2\t43\t1e-20\t150"
    cases = {
        good.replace("\t150", ""): "LOSATP blastp output does not match outfmt 6 columns.",
        good.replace("\t42\t1\t", "\t42.0\t1\t", 1): "LOSATP blastp output contains non-numeric outfmt 6 fields.",
        good.replace("1e-20", "1e999"): "LOSATP blastp output contains non-numeric outfmt 6 fields.",
    }
    for line, message in cases.items():
        with pytest.raises(ParseError) as raised:
            parse_losatp_outfmt6(f"{good}\n{line}\n")
        assert str(raised.value) == message
    assert len(parse_losatp_outfmt6(f"{good}\n{good}\n")) == 2


_SYNTHETIC_JS = """\
import { COMPARISON_COLUMNS, OTHER } from '../contracts/schema.generated.js';
import { COMPARISON_COLUMNS as BLAST_COLUMNS } from '../contracts/schema.generated.js';
import * as schema from '../contracts/schema.generated.js';

/** Not a mapping: COMPARISON_COLUMNS.map((column) => columns.indexOf(column)) in a comment. */
const ACCEPT = '.txt,text/*,.gz';
const QUOTE = /["'/*]/g;
const HEADER = COMPARISON_COLUMNS.join('\\t');
const labels = COMPARISON_COLUMNS.map((column) => LABELS[column]);
const known = (column) => COMPARISON_COLUMNS.includes(column) && COMPARISON_COLUMNS.length > 0;
const writeRow = (row) => COMPARISON_COLUMNS.map((column) => String(row?.[column])).join('\\t');

export const registeredReorder = (columns) => {
  const positions = COMPARISON_COLUMNS.map((column) => columns.indexOf(column));
  return positions;
};

/** `/*` in the string above does not hide this one. */
function headerMap(index) {
  return COMPARISON_COLUMNS.map((column) => index.get(column));
}

const namedFields = (fields) => Object.fromEntries(fields.map((value, index) => [COMPARISON_COLUMNS[index], value]));

const looped = (fields) => {
  const named = {};
  for (const [index, column] of COMPARISON_COLUMNS.entries()) {
    named[column] = fields[index];
  }
  return named;
};

const aliased = (columns) => BLAST_COLUMNS.reduce(
  (found, column) => found.concat(columns.findIndex((c) => c === column)), []);

const positionMap = () => new Map(schema.COMPARISON_COLUMNS.map((column, index) => [column, index]));

const evalueOf = (fields) => fields[COMPARISON_COLUMNS.indexOf('evalue')];

const fromSlice = (fields) => Array.from(COMPARISON_COLUMNS.slice(2), (column, index) => Number(fields[index + 2]));

const INDEX = {};
COMPARISON_COLUMNS.forEach((column, index) => { INDEX[column] = index; });
"""

_SYNTHETIC_PY = """\
import pandas as pd
from gbdraw.io.comparisons import COMPARISON_COLUMNS


class Reader:
    def read(self, path):
        return pd.read_csv(path, sep="\\\\t", names=list(COMPARISON_COLUMNS))


def relabel(frame):
    frame.columns = [*COMPARISON_COLUMNS]
    return frame


ORDER = list(COMPARISON_COLUMNS)
HEADER = "\\\\t".join(COMPARISON_COLUMNS)
IDS = COMPARISON_COLUMNS[:2]
FIRST = COMPARISON_COLUMNS[0]
POSITIONS = {column: index for index, column in enumerate(COMPARISON_COLUMNS)}


def not_mappings(column, frame, row):
    empty = {name: [] for name in COMPARISON_COLUMNS}
    known = column in COMPARISON_COLUMNS, COMPARISON_COLUMNS.index(column) + 1, len(COMPARISON_COLUMNS)
    rows = [dict(zip(COMPARISON_COLUMNS, values)) for values in frame.itertuples(index=False)]
    written = "\\\\t".join(str(row[name]) for name in COMPARISON_COLUMNS)
    labels = [LABELS[name] for name in COMPARISON_COLUMNS]
    series = [frame[name] for name in COMPARISON_COLUMNS]
    for order, name in enumerate(COMPARISON_COLUMNS):
        check(frame[name], order)
    built = pd.DataFrame(columns=COMPARISON_COLUMNS), pd.DataFrame(rows, columns=COMPARISON_COLUMNS)
    return empty, known, rows, written, labels, series, built


def positions_of(header):
    return [header.index(column) for column in COMPARISON_COLUMNS]


def by_header(index):
    return tuple(index.get(column) for column in COMPARISON_COLUMNS)


def numbered():
    return dict(enumerate(COMPARISON_COLUMNS))


def by_position(fields):
    named = {}
    for index in range(len(fields)):
        named[COMPARISON_COLUMNS[index]] = fields[index]
    return named


def outer(rows):
    def inner(row):
        named = {}
        for index, column in enumerate(COMPARISON_COLUMNS):
            named[column] = row[index]
        return named

    return [inner(row) for row in rows]


def strict(fields):
    return all(valid(column, fields[index]) for index, column in enumerate(COMPARISON_COLUMNS[2:], start=2))


def evalue_of(fields):
    return fields[COMPARISON_COLUMNS.index("evalue")]


def renamed(frame):
    return frame.rename(columns=dict(zip(frame.columns, COMPARISON_COLUMNS)))


def aliased(header):
    from gbdraw.io.comparisons import COMPARISON_COLUMNS as COLUMNS

    index_of = {}
    for position, column in enumerate(COLUMNS):
        index_of[column] = position
    return index_of
"""


def test_column_mapping_detector_names_each_mapping_by_function() -> None:
    """The one-reader guard sees each spelling of a mapping, by file and top-level function, and reads the twin rows."""

    readers = {}
    for path, text in (("gbdraw/web/js/x.js", _SYNTHETIC_JS), ("gbdraw/y.py", _SYNTHETIC_PY)):
        for site in _column_mappings(path, text):
            readers[site] = readers.get(site, 0) + 1
    mappings = {
        "gbdraw/web/js/x.js:registeredReorder": 1,
        "gbdraw/web/js/x.js:headerMap": 1,
        "gbdraw/web/js/x.js:namedFields": 1,
        "gbdraw/web/js/x.js:looped": 1,
        "gbdraw/web/js/x.js:aliased": 1,
        "gbdraw/web/js/x.js:positionMap": 1,
        "gbdraw/web/js/x.js:evalueOf": 1,
        "gbdraw/web/js/x.js:fromSlice": 1,
        "gbdraw/web/js/x.js:<module>": 1,
        "gbdraw/y.py:Reader.read": 1,
        "gbdraw/y.py:relabel": 1,
        "gbdraw/y.py:<module>": 1,
        "gbdraw/y.py:positions_of": 1,
        "gbdraw/y.py:by_header": 1,
        "gbdraw/y.py:numbered": 1,
        "gbdraw/y.py:by_position": 1,
        "gbdraw/y.py:outer": 1,
        "gbdraw/y.py:strict": 1,
        "gbdraw/y.py:evalue_of": 1,
        "gbdraw/y.py:renamed": 1,
        "gbdraw/y.py:aliased": 1,
    }
    assert readers == mappings

    registry = {"twins": [
        # A twin of the owner module: its owners are exempt.
        {"pythonOwner": [f"{_OWNER_MODULE}:read_comparison_table"],
         "jsOwner": ["gbdraw/web/js/x.js:registeredReorder"]},
        # Pending, or a twin of another module: its owners are not.
        {"pythonOwner": [f"{_OWNER_MODULE}:read_comparison_table"],
         "jsOwner": ["gbdraw/web/js/x.js:namedFields"], "pending": True},
        {"pythonOwner": ["gbdraw/y.py:Reader.read"], "jsOwner": ["gbdraw/web/js/x.js:headerMap"]},
    ]}
    unregistered = {site: 1 for site in mappings if site != "gbdraw/web/js/x.js:registeredReorder"}
    assert _unregistered_column_mappings(
        {**readers, f"{_OWNER_MODULE}:read_comparison_table": 1}, registry
    ) == unregistered


# --- CO-06: record-ID binding --------------------------------------------------


def _random_sequence(seed: int, length: int = 3000) -> str:
    generator = random.Random(seed)
    return "".join(generator.choice("ACGT") for _ in range(length))


def _pair_records(*ids: str) -> list[SeqRecord]:
    records = []
    for index, record_id in enumerate(ids):
        record = SeqRecord(Seq(_random_sequence(index)), id=record_id, name=record_id.split(".")[0])
        record.annotations["molecule_type"] = "DNA"
        records.append(record)
    return records


def _hits(query: str, subject: str) -> pd.DataFrame:
    row = [query, subject, 99.0, 1000, 0, 0, 2001, 3000, 1, 1000, 1e-100, 1800.0]
    return pd.DataFrame([row], columns=COMPARISON_COLUMNS)


def _render_svg(records: list[SeqRecord], comparisons: list[LinearComparison]) -> str:
    return assemble_linear_diagram_from_records(
        records,
        LinearAssemblyInputs(
            cfg=apply_config_overrides(
                None,
                {
                    "labels.linear.scope": "none",
                    "canvas.show_gc": False,
                    "canvas.show_skew": False,
                },
            ),
            linear_comparisons=comparisons,
            legend="none",
        ),
    ).tostring()


def test_table_matching_the_endpoint_records_draws_with_endpoint_metadata() -> None:
    svg = _render_svg(_pair_records("R2", "R3"), [LinearComparison(0, 1, _hits("R2", "R3"))])

    assert 'data-query-record-id="R2"' in svg
    assert 'data-subject-record-id="R3"' in svg
    assert 'data-qstart="2001"' in svg


def test_version_suffix_differences_are_tolerated() -> None:
    svg = _render_svg(
        _pair_records("R2.1", "R3.2"), [LinearComparison(0, 1, _hits("R2.3", "R3"))]
    )

    assert 'data-query-record-id="R2.1"' in svg
    assert 'data-subject-record-id="R3.2"' in svg


def test_swapped_query_and_subject_table_is_rejected() -> None:
    with pytest.raises(
        ComparisonIdentityError, match=r"the query column names 'R3', which is the subject record"
    ) as error:
        _render_svg(_pair_records("R2", "R3"), [LinearComparison(0, 1, _hits("R3", "R2"))])

    assert error.value.reason == "RECORD_ID"
    # CI-04: the Web names the correction and the column instead of the generic
    # comparison-identity text.
    payload = _web(error.value)
    assert payload["code"] == "COMPARISON_INPUT"
    assert payload["context"] == {
        "reason": "RECORD_ID", "column": 1, "queryRecordIndex": 0, "subjectRecordIndex": 1,
    }


def test_table_naming_another_displayed_record_is_rejected() -> None:
    records = _pair_records("R1", "R2", "R3")
    with pytest.raises(ComparisonIdentityError, match=r"the subject column names 'R3', which is record #3 'R3'"):
        _render_svg(records, [LinearComparison(0, 1, _hits("R1", "R3"))])


def test_distinct_versions_of_one_accession_stay_distinct() -> None:
    records = _pair_records("X.1", "X.2", "Y")
    with pytest.raises(ComparisonIdentityError, match=r"the query column names 'X.1', which is record #1"):
        _render_svg(records, [LinearComparison(1, 2, _hits("X.1", "Y"))])


def test_unknown_table_ids_keep_positional_placement_with_a_warning(
    caplog: pytest.LogCaptureFixture,
) -> None:
    with caplog.at_level(logging.WARNING, logger="gbdraw.linear_comparison"):
        svg = _render_svg(
            _pair_records("R2", "R3"), [LinearComparison(0, 1, _hits("contig_A", "contig_B"))]
        )

    assert 'data-query-record-id="R2"' in svg
    assert 'data-subject-record-id="R3"' in svg
    assert 'data-query-record-index="0"' in svg
    assert "1 row(s) use sequence IDs that match no displayed record ('contig_A', 'contig_B')" in caplog.text


def test_unknown_table_ids_are_returned_as_structured_warnings() -> None:
    """PD-OI-074: the Web carries the warning the CLI logs (logging is off in the browser)."""

    def build(*tables: pd.DataFrame):
        return assemble_linear_diagram_from_records(
            _pair_records("R2", "R3"),
            LinearAssemblyInputs(
                cfg=apply_config_overrides(
                    None,
                    {
                        "labels.linear.scope": "none",
                        "canvas.show_gc": False,
                        "canvas.show_skew": False,
                    },
                ),
                linear_comparisons=[LinearComparison(0, 1, table) for table in tables],
                legend="none",
            ),
            _return_build_result=True,
        )

    warnings = build(_hits("contig_A", "contig_B")).metadata.comparison_warnings
    assert len(warnings) == 1
    warning = warnings[0]
    assert warning.code == "comparison_record_id_unmatched"
    assert (warning.query_record_index, warning.subject_record_index) == (0, 1)
    assert (warning.query_record_id, warning.subject_record_id) == ("R2", "R3")
    assert warning.row_count == 1
    assert warning.example_ids == ("contig_A", "contig_B")
    assert warning.message.startswith("Comparison between query record #1 'R2' and subject record #2 'R3'")

    # Matching IDs, including version-suffix differences, warn about nothing.
    assert build(_hits("R2", "R3")).metadata.comparison_warnings == ()


def test_uploaded_rows_merged_with_source_bound_rows_are_still_checked() -> None:
    records = _pair_records("R2", "R3")
    # Source-bound rows are checked by feature identity, so their table IDs are exempt.
    bound = _hits("R3", "R2").assign(
        query_feature_index=0,
        subject_feature_index=0,
        query_feature_svg_id="f_query",
        subject_feature_svg_id="f_subject",
    )
    svg = _render_svg(records, [LinearComparison(0, 1, bound)])
    assert 'data-query-record-id="R2"' in svg

    # Sources for one pair are merged; the unbound uploaded rows keep the ID check.
    with pytest.raises(ComparisonIdentityError, match=r"the query column names 'R3'"):
        _render_svg(
            records,
            [LinearComparison(0, 1, bound), LinearComparison(0, 1, _hits("R3", "R2"))],
        )


def _write_genbank(path: Path, record: SeqRecord) -> Path:
    SeqIO.write([record], path, "genbank")
    return path


def _write_hits(path: Path, frame: pd.DataFrame) -> Path:
    frame.to_csv(path, sep="\t", header=False, index=False)
    return path


def test_web_request_and_cli_reject_the_same_swapped_table(tmp_path: Path, gbdraw_runner) -> None:
    from gbdraw.api.request_render import render_request
    from gbdraw.api import InMemoryRecordSource, LinearDiagramRequest, RecordInput
    from gbdraw.session_request_codec import encode_canonical_request

    records = _pair_records("R2", "R3")
    swapped = _write_hits(tmp_path / "R3_R2.tsv", _hits("R3", "R2"))

    request = LinearDiagramRequest(
        records=tuple(RecordInput(source=InMemoryRecordSource(record)) for record in records),
        options=LinearDiagramOptions(linear_comparisons=(LinearComparison(0, 1, _hits("R2", "R3")),)),
    )
    encoded = encode_canonical_request(request)
    resource_paths = {}
    for resource in encoded.resources:
        target = tmp_path / resource.name
        target.write_bytes(resource.content)
        resource_paths[resource.resource_id] = target
    comparison = encoded.payload["comparisons"][0]
    assert comparison["kind"] == "nucleotideBlast"  # a BLAST table is stored as one (OV-3512)
    resource_paths[comparison["resourceId"]] = swapped
    decoded = decode_canonical_request(
        encoded.payload, resource_paths=resource_paths, output_directory=tmp_path / "out"
    )
    with pytest.raises(ComparisonIdentityError, match=r"the query column names 'R3', which is the subject record") as web_error:
        render_request(decoded)

    gbk_files = [
        _write_genbank(tmp_path / f"{record.id}.gb", record) for record in records
    ]
    returncode, output, _svg = gbdraw_runner.run(
        "linear", gbk_files, "swapped", tmp_path, blast_files=[swapped]
    )

    assert returncode != 0
    assert str(web_error.value) in output


def test_cli_logs_the_comparison_thresholds_and_files_it_reads(tmp_path: Path, gbdraw_runner) -> None:
    records = _pair_records("SA", "SB")
    gbk_files = [_write_genbank(tmp_path / f"{record.id}.gb", record) for record in records]
    hits = _write_hits(tmp_path / "SA_SB.tsv", _hits("SA", "SB"))

    returncode, output, _svg = gbdraw_runner.run(
        "linear", gbk_files, "logged", tmp_path, blast_files=[hits], extra_args=["--evalue", "1e-5"]
    )

    assert returncode == 0, output
    assert output.count("INFO: BLAST output visualization settings: e-value threshold: 1e-05;") == 1
    assert output.count("INFO: Loading comparison file(s)...") == 1
    assert output.count(f"INFO: Loading {hits}") == 1


def test_cli_exits_non_zero_for_a_missing_comparison_file(tmp_path: Path, gbdraw_runner) -> None:
    records = _pair_records("SA", "SB", "SC")
    gbk_files = [_write_genbank(tmp_path / f"{record.id}.gb", record) for record in records]
    valid = _write_hits(tmp_path / "SB_SC.tsv", _hits("SB", "SC"))

    returncode, output, svg_path = gbdraw_runner.run(
        "linear", gbk_files, "missing", tmp_path, blast_files=[tmp_path / "SA_SB.tsv", valid]
    )

    assert returncode != 0
    assert re.search(r"SA_SB\.tsv: the comparison file does not exist", output)
    assert not svg_path.exists()
