"""tools/fix_metrics.py on a small repository built here: no network, no use of this clone's history."""

from __future__ import annotations

import datetime as dt
import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

TOOL = Path(__file__).resolve().parents[1] / "tools" / "fix_metrics.py"
_SPEC = importlib.util.spec_from_file_location("fix_metrics", TOOL)
assert _SPEC and _SPEC.loader
fix_metrics = importlib.util.module_from_spec(_SPEC)
sys.modules["fix_metrics"] = fix_metrics
_SPEC.loader.exec_module(fix_metrics)

JST = dt.timezone(dt.timedelta(hours=9))
MONDAY = dt.datetime(2026, 11, 2, 10, 0, tzinfo=JST)  # day 0
MODULE = "gbdraw/mod.py"


class Repo:
    def __init__(self, path: Path) -> None:
        self.path = path
        self.number = 0
        self.git("init", "-q", "-b", "dev")

    def git(self, *args: str, day: float = 0.0) -> str:
        stamp = (MONDAY + dt.timedelta(days=day)).isoformat()
        env = {**os.environ, "GIT_AUTHOR_DATE": stamp, "GIT_COMMITTER_DATE": stamp, "GIT_AUTHOR_NAME": "t",
               "GIT_AUTHOR_EMAIL": "t@example.com", "GIT_COMMITTER_NAME": "t", "GIT_COMMITTER_EMAIL": "t@example.com",
               "GIT_CONFIG_GLOBAL": os.devnull, "GIT_CONFIG_NOSYSTEM": "1"}
        return subprocess.run(["git", "-C", str(self.path), *args], check=True, capture_output=True, text=True,
                              env=env).stdout

    def write(self, files: dict[str, str]) -> None:
        for name, text in files.items():
            target = self.path / name
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(text, encoding="utf-8")
            self.git("add", name)

    def branch(self, name: str, files: dict[str, str], day: float) -> None:
        self.git("switch", "-q", "-c", name, "dev")
        self.write(files)
        self.git("commit", "-q", "-m", f"work on {name}", day=day - 0.1)
        self.git("switch", "-q", "dev")

    def pull_request(self, name: str, title: str, files: dict[str, str], day: float) -> str:
        """Branch, commit and merge as GitHub does; returns '#<number>'."""
        self.number += 1
        self.branch(name, files, day)
        self.git("merge", "-q", "--no-ff", name, "-m", f"Merge pull request #{self.number} from owner/{name}",
                 "-m", title, day=day)
        return f"#{self.number}"


def module(lines: dict[int, str]) -> str:
    """Ten statement lines; ``lines`` replaces some of them."""
    return "".join(lines.get(n, f"value_{n} = {n}") + "\n" for n in range(1, 11))


def twins(rows: int, without_vector: int) -> str:
    return json.dumps({"twins": [{"id": f"concern:c{n}", "vector": None if n < without_vector else "tests/fixtures/v.json"}
                                 for n in range(rows)]})


@pytest.fixture(scope="module")
def repo(tmp_path_factory: pytest.TempPathFactory) -> Repo:
    """
    Day 0  #1 feat: the module (10 lines), two twin rows (one without a vector), both baselines.
    Day 1  #2 fix: rewrites lines 1-4 (#1's, 1 day old).
    Day 2  #3 fix: rewrites lines 1-3 (#2's).
    Day 3  #4 fix: rewrites lines 1-3 (#3's); #3 rewrote #2 there: a third fix (#2 .. #4: 2 days).
           Its comment and blank lines are not counted.
    Day 4  #5 fix: a test only (no production change).
    Day 12 #6 fix: rewrites lines 7-9 (#1's, 12 days old) and drops #4's comment line; a twin row added;
           the Python baseline falls.
    """
    repo = Repo(tmp_path_factory.mktemp("history"))
    repo.write({"README.md": "x\n"})
    repo.git("commit", "-q", "-m", "root", day=-1)
    repo.pull_request("feat/module", "Add the module", {
        MODULE: module({}), "tools/twins.json": twins(2, 1),
        "tests/test_web_error_producer_coverage.py": "UNPROVEN_BASELINE = 10\n",
        "tests/web/error-producer-coverage.test.mjs": "const UNCLASSIFIED_THROW_BASELINE = 7;\n",
    }, day=0)
    repo.pull_request("fix/first", "Keep the label clear of the ring", {
        MODULE: module({n: f"value_{n} = {n} + 1" for n in range(1, 5)})}, day=1)
    repo.pull_request("fix/second", "Restore the colors after Undo", {
        MODULE: module({**{n: f"value_{n} = {n} + 2" for n in range(1, 4)}, 4: "value_4 = 4 + 1"})}, day=2)
    third = {**{n: f"value_{n} = {n} + 3" for n in range(1, 4)}, 4: "value_4 = 4 + 1"}
    repo.pull_request("fix/third", "Name the slot", {MODULE: module(third) + "# note\n\n"}, day=3)
    repo.pull_request("fix/tests-only", "Pin the slot name", {"tests/test_slot.py": "def test_x():\n    pass\n"},
                      day=4)
    later = {**third, **{n: f"value_{n} = {n} * 2" for n in range(7, 10)}}
    repo.pull_request("fix/later", "Draw the axis", {
        MODULE: module(later) + "\n", "tools/twins.json": twins(3, 1),
        "tests/test_web_error_producer_coverage.py": "UNPROVEN_BASELINE = 9\n",
    }, day=12)
    return repo


def test_week_counts_rework_lineage_and_size(repo: Repo) -> None:
    report = fix_metrics.period_report(fix_metrics.Git(repo.path), "dev", fix_metrics.day_start("2026-11-02"),
                                       fix_metrics.day_start("2026-11-09"), jobs=2)
    row = report.total
    assert (row.fix_prs, row.with_rewritten_lines, row.rework) == (4, 3, 3)
    assert (row.rewrite_recent_fix, row.third_fixes) == (2, 1)
    assert row.no_production_change == 1
    assert (row.state_path_titles, row.layout_titles) == (1, 1)
    entries = {fix.pr: fix for fix in report.fixes}
    third = entries["#4"].third_fix
    assert [(t.path, t.rewrites, t.which_rewrote, t.rewrites_lines) for t in third] == [(MODULE, "#3", "#2", [1, 2, 3])]
    assert third[0].span_days == 2.0
    assert [(r.of, r.lines) for r in entries["#3"].rewrites_recent_fix] == [("#2", [1, 2, 3])]
    assert entries["#2"].rewrites_recent_fix == []  # #1 is not a fix
    assert entries["#4"].net_production_lines == 2  # a comment line and a blank line added
    assert [w.label for w in report.weeks] == ["2026-11-02"]


def test_lineage_stays_in_one_file(tmp_path: Path) -> None:
    """#4 rewrites #3 in a.py, but #3 rewrote #2 only in b.py: no third fix."""
    repo = Repo(tmp_path)
    repo.write({"README.md": "x\n"})
    repo.git("commit", "-q", "-m", "root", day=-1)
    a, b = "gbdraw/a.py", "gbdraw/b.py"
    first_three = {n: f"value_{n} = 'a'" for n in range(1, 4)}
    repo.pull_request("feat/two", "Add two modules", {a: module({}), b: module({})}, day=0)
    repo.pull_request("fix/b", "Fix b", {b: module(first_three)}, day=1)
    repo.pull_request("fix/a-and-b", "Fix a and b", {a: module(first_three),
                                                     b: module({n: f"value_{n} = 'b'" for n in range(1, 4)})}, day=2)
    # #4 also rewrites one line of #3 in b.py: fewer than three lines is no rewrite of #3 there.
    repo.pull_request("fix/a", "Fix a", {a: module({n: f"value_{n} = 'c'" for n in range(1, 4)}),
                                         b: module({1: "value_1 = 'd'", 2: "value_2 = 'b'", 3: "value_3 = 'b'"})},
                      day=3)
    report = fix_metrics.period_report(fix_metrics.Git(repo.path), "dev", fix_metrics.day_start("2026-11-02"),
                                       fix_metrics.day_start("2026-11-09"), jobs=1)
    last = report.fixes[-1]
    assert [(r.path, r.of) for r in last.rewrites_recent_fix] == [(a, "#3")]
    assert last.third_fix == [] and report.total.third_fixes == 0


def test_lineage_spans_seven_days_and_half_the_lines_is_rework(tmp_path: Path) -> None:
    """#4 rewrites #3 (4 days), #3 rewrote #2 (5.5 days), but #2 merged 9.5 days before #4: no third fix.
    Half of #4's lines are recent: rework."""
    repo = Repo(tmp_path)
    repo.write({"README.md": "x\n"})
    repo.git("commit", "-q", "-m", "root", day=-1)
    a = "gbdraw/a.py"
    repo.pull_request("feat/a", "Add a", {a: module({})}, day=0)
    repo.pull_request("fix/one", "One", {a: module({n: f"value_{n} = 'one'" for n in range(1, 4)})}, day=0.5)
    repo.pull_request("fix/two", "Two", {a: module({n: f"value_{n} = 'two'" for n in range(1, 4)})}, day=6)
    repo.pull_request("fix/three", "Three", {a: module({n: f"value_{n} = 'three'" for n in [1, 2, 3, 7, 8, 9]})},
                      day=10)
    report = fix_metrics.period_report(fix_metrics.Git(repo.path), "dev", fix_metrics.day_start("2026-11-02"),
                                       fix_metrics.day_start("2026-11-16"), jobs=1)
    last = report.fixes[-1]
    assert [(r.of, r.age_days) for r in last.rewrites_recent_fix] == [("#3", 4.0)]
    assert report.fixes[1].third_fix == [] and last.third_fix == []
    assert (last.rewritten_lines, last.recent_lines, report.weeks[-1].rework) == (6, 3, 1)


def test_range_reads_registry_and_baselines_at_both_ends(repo: Repo) -> None:
    report = fix_metrics.period_report(fix_metrics.Git(repo.path), "dev", fix_metrics.day_start("2026-11-03"),
                                       fix_metrics.day_start("2026-11-15"), jobs=1)
    assert [w.label for w in report.weeks] == ["2026-11-02", "2026-11-09"]
    assert report.twins_start == fix_metrics.Registry(rows=2, without_vector=1)
    assert [w.twins for w in report.weeks] == [fix_metrics.Registry(rows=2, without_vector=1),
                                               fix_metrics.Registry(rows=3, without_vector=1)]
    assert report.twin_changes == [fix_metrics.TwinChange("#6", "2026-11-14 10:00", "2026-11-09", ["concern:c2"], [])]
    assert [(w.twin_rows_added, w.twin_rows_removed) for w in report.weeks] == [(0, 0), (1, 0)]
    assert report.error_baselines_start == {"python_unclassified_validation_errors": 10,
                                            "javascript_unclassified_throws": 7}
    assert report.total.error_baselines["python_unclassified_validation_errors"] == 9
    later = report.fixes[-1]
    assert (later.pr, later.rewritten_lines, later.recent_lines, later.rewrites_recent_fix) == ("#6", 3, 0, [])
    assert report.weeks[1].rework == 0


def test_branch_names_the_lines_that_make_it_a_third_fix(repo: Repo) -> None:
    repo.branch("fix/fourth", {MODULE: module({**{n: f"value_{n} = {n} + 4" for n in range(1, 4)},
                                                4: "value_4 = 4 + 1", **{n: f"value_{n} = {n} * 2" for n in range(7, 10)}})
                               + "# note\n\n"}, day=13)
    git = fix_metrics.Git(repo.path)
    at = int((MONDAY + dt.timedelta(days=4)).timestamp())
    report = fix_metrics.branch_report(git, "dev", "fix/fourth", at)
    assert report.merged_cleanly and report.kind == "fix"
    assert [(r.of, r.lines) for r in report.rewrites_recent_fix] == [("#4", [1, 2, 3])]
    assert [(t.rewrites, t.which_rewrote) for t in report.third_fix] == [("#4", "#3")]
    text = fix_metrics.branch_text(report)
    assert "Third fix in a row on the same lines: yes" in text
    assert "lines 1-3 of #4 (Name the slot)" in text

    later = fix_metrics.branch_report(git, "dev", "fix/fourth", int((MONDAY + dt.timedelta(days=11)).timestamp()))
    assert later.rewrites_recent_fix == [] and later.third_fix == []
    assert "Rewrites 3+ lines of a fix merged in the last 7 days: no" in fix_metrics.branch_text(later)


def test_branch_that_conflicts_is_judged_from_its_merge_base(repo: Repo) -> None:
    repo.git("switch", "-q", "-c", "fix/conflict", "dev~2")  # before #6 rewrote lines 7-9
    repo.write({MODULE: module({**{n: f"value_{n} = {n} + 3" for n in range(1, 4)}, 4: "value_4 = 4 + 1",
                                7: "value_7 = -7"}) + "# note\n\n"})
    repo.git("commit", "-q", "-m", "conflict", day=12.5)
    repo.git("switch", "-q", "dev")
    report = fix_metrics.branch_report(fix_metrics.Git(repo.path), "dev", "fix/conflict",
                                       int((MONDAY + dt.timedelta(days=13)).timestamp()))
    assert not report.merged_cleanly
    assert report.rewritten_lines == 1 and report.rewrites_recent_fix == []


def test_command_line_writes_markdown_and_json(repo: Repo, tmp_path: Path) -> None:
    out = tmp_path / "week.json"
    done = subprocess.run([sys.executable, str(TOOL), "--repo", str(repo.path), "--ref", "dev", "--week", "2026-11-02",
                           "--json", str(out)], check=True, capture_output=True, text=True)
    assert "| 2026-11-02 | 4 | 3 | 3 (100%) | 2 (50%) | 1 (25%) |" in done.stdout
    assert "- #4 2026-11-05 10:00 Name the slot: `gbdraw/mod.py` rewrites 3 lines of #3" in done.stdout
    data = json.loads(out.read_text(encoding="utf-8"))
    assert data["total"]["third_fixes"] == 1 and data["fixes"][2]["third_fix"][0]["which_rewrote"] == "#2"
    refused = subprocess.run([sys.executable, str(TOOL), "--repo", str(repo.path), "--week", "2026-11-03"],
                             capture_output=True, text=True)
    assert refused.returncode == 2 and "is not a Monday" in refused.stderr


@pytest.mark.parametrize(("body", "expected"), [
    ("## Boundary analysis\n- Cause class: state path (P7); guard: tests/web/x.test.mjs\n", "state path"),
    ("- **Cause class:** Duplicate-implementation, guard: none because ...", "duplicate implementation"),
    ("Cause class: lookup by name or position", "lookup by name or position"),
    ("Cause class: a typo", "unrecognised"),
    ("No line here", "none"),
])
def test_cause_class_line(body: str, expected: str) -> None:
    assert fix_metrics.cause_class(body) == expected


def test_production_paths_and_comment_lines() -> None:
    assert fix_metrics.is_production("gbdraw/web/js/app/app-setup.js")
    assert not fix_metrics.is_production("gbdraw/web/js/contracts/session.js")
    assert not fix_metrics.is_production("tests/test_x.py")
    assert not fix_metrics.is_production("gbdraw/data/fonts/a.ttf")
    assert fix_metrics.is_comment_or_blank("   // a comment") and fix_metrics.is_comment_or_blank("  ")
    assert not fix_metrics.is_comment_or_blank("x = 1  # trailing")
