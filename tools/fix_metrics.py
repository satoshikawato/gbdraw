#!/usr/bin/env python3
"""Count how fix pull requests rework recent changes, from the first-parent history of a branch.

The history is read with git only (no network), unless ``--github`` asks for the two
values that only GitHub holds. Nothing is written to the repository.

    python tools/fix_metrics.py                                  # the last full week (Monday JST)
    python tools/fix_metrics.py --week 2026-11-09
    python tools/fix_metrics.py --since 2026-10-12 --until 2026-11-15 --json report.json
    python tools/fix_metrics.py --branch fix/my-change            # before opening a pull request

A pull request is a first-parent merge "Merge pull request #N from <owner>/<kind>/<name>";
a fix pull request is one whose kind is ``fix``. A rewritten line is a line of production
code (``PRODUCTION_SUFFIXES`` under ``gbdraw/``, minus ``NOT_PRODUCTION``) that the merge
deletes or changes and that is neither blank nor a comment. ``git blame`` at the merge's
first parent names the pull request that wrote it; its age is the time between the two
merges. ``git blame`` names the last change to a line, which is not always the change that
caused a defect; counted the same way every week, the trend still shows.

- Rework: half or more of a fix's rewritten lines are at most ``WINDOW_DAYS`` old.
- Rewrite of a recent fix: in one file, a fix rewrites at least ``LINEAGE_LINES`` lines
  that another fix merged at most ``WINDOW_DAYS`` earlier wrote.
- Third fix: fix F rewrites a recent fix G that way, G rewrote a fix H that way in the
  same file, and H merged at most ``WINDOW_DAYS`` before F.

``--branch`` judges a work branch the same way, as the merge of it into the reference
branch would be (``git merge-tree``), at the current time or ``--at``.
"""

from __future__ import annotations

import argparse
import dataclasses
import datetime as dt
import json
import os
import re
import statistics
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Iterable

JST = dt.timezone(dt.timedelta(hours=9))
DAY = 86400
WINDOW_DAYS = 7
WINDOW = WINDOW_DAYS * DAY
LINEAGE_LINES = 3
REWORK_SHARE = 0.5
PRODUCTION_SUFFIXES = (".py", ".js", ".mjs", ".html", ".css", ".toml")
NOT_PRODUCTION = (
    "gbdraw/web/vendor/",
    "gbdraw/web/gallery/",
    "gbdraw/web/js/contracts/",
    "gbdraw/web/js/utils/i18n-catalogs/",
    "gbdraw/web/css/tailwind.css",
)
# The line starts that open or continue a comment, by file suffix: `*args` in Python and a `#private`
# field in JavaScript are code.
COMMENT_PREFIXES = {".py": ("#",), ".toml": ("#",), ".js": ("//", "/*", "*"), ".mjs": ("//", "/*", "*"),
                    ".css": ("/*", "*"), ".html": ("<!--",)}
# A file a merge deletes or rewrites wholesale says nothing about rework; its lines are not blamed.
MAX_REMOVED_LINES = 600
DIRECT = "direct"

TWIN_REGISTRY = "tools/twins.json"
# The two unclassified-error ratchets: one total, or (before 2026-10-12) a table of counts per file.
ERROR_BASELINES = (
    ("python_unclassified_validation_errors", "tests/test_web_error_producer_coverage.py",
     re.compile(r"^UNPROVEN_BASELINE\b[^=]*=\s*(\d+|\{.*?^\})", re.M | re.S)),
    ("javascript_unclassified_throws", "tests/web/error-producer-coverage.test.mjs",
     re.compile(r"^const UNCLASSIFIED_THROW_BASELINE\s*=\s*(\d+|\{.*?^\})", re.M | re.S)),
)
BASELINE_ENTRY = re.compile(r":\s*\(?\s*(\d+)")
# Fix titles that name a state path or a layout rule (words of the title, any case).
STATE_PATH_TITLE = re.compile(
    r"\b(?:sessions?|(?:re)?load(?:s|ed|ing)?|sav(?:e|es|ed|ing)|undo(?:es|ne)?|redo|history|"
    r"reset(?:s|ting)?|switch(?:es|ed|ing)?|generate[sd]?|replay(?:s|ed|ing)?|restor(?:e|es|ed|ing))\b",
    re.I,
)
LAYOUT_TITLE = re.compile(
    r"\b(?:overlap\w*|clear of|room|margins?|padding|spacing|space\w*\b.*\bapart|gap between|overflow\w*|"
    r"cross(?:es|ing)|collid\w*|collision\w*|clip(?:s|ped|ping)?|touch(?:es|ing)|make way|"
    r"keep\b.*\b(?:off|outside|inside))\b",
    re.I,
)
# The cause class line of a fix pull request's Boundary analysis: "Cause class: <class>; ...".
CAUSE_CLASS_LINE = re.compile(r"^[ \t]*(?:[-*][ \t]+)?\**cause class\**[ \t]*:[ \t]*(.+)$", re.I | re.M)
CAUSE_CLASSES = (
    "duplicate implementation",
    "state path",
    "undefined edge",
    "combination",
    "single-site logic error",
    "lookup by name or position",
    "performance or memory",
    "environment",
    "not a product defect",
    "unknown",
)
MERGE_SUBJECT = re.compile(r"^Merge pull request #(\d+) from ([^/\s]+)/(\S+)")
HUNK = re.compile(r"^@@ -(\d+)(?:,\d+)? \+\d+(?:,\d+)? @@")
BLAME_HEADER = re.compile(r"^([0-9a-f]{40}) \d+ (\d+)")


class Git:
    def __init__(self, repo: Path, env: dict[str, str] | None = None) -> None:
        self.repo = repo
        self.env = {**os.environ, **env} if env else None

    def attempt(self, *args: str) -> tuple[int, str, str]:
        done = subprocess.run(["git", "-C", str(self.repo), "-c", "core.quotepath=false", *args],
                              capture_output=True, check=False, env=self.env)
        return done.returncode, done.stdout.decode("utf-8", errors="replace"), done.stderr.decode(errors="replace")

    def run(self, *args: str) -> str:
        code, out, err = self.attempt(*args)
        if code:
            raise RuntimeError(f"git {' '.join(args)}: {err.strip()}")
        return out

    def show(self, revision: str, path: str) -> str | None:
        code, text, _ = self.attempt("show", f"{revision}:{path}")
        return text if code == 0 else None


@dataclasses.dataclass(frozen=True)
class Merge:
    """One first-parent merge of the reference branch; ``key`` is ``#<number>`` for a pull request."""

    key: str
    branch: str
    kind: str
    sha: str
    base: str
    time: int
    title: str


@dataclasses.dataclass(frozen=True)
class Change:
    """What a merge did to production code: the rewritten lines by file and writer, and its size."""

    origins: dict[str, dict[str, list[int]]]
    touched: tuple[str, ...]
    added: int
    removed: int
    skipped: tuple[str, ...]


@dataclasses.dataclass(frozen=True)
class Commit:
    sha: str
    time: int


def is_production(path: str) -> bool:
    return (path.startswith("gbdraw/") and path.endswith(PRODUCTION_SUFFIXES)
            and not path.startswith(NOT_PRODUCTION))


def is_comment_or_blank(path: str, text: str) -> bool:
    stripped = text.strip()
    return not stripped or stripped.startswith(COMMENT_PREFIXES.get(Path(path).suffix, ()))


def kind_of(branch_path: str) -> str:
    """``fix/name`` -> ``fix``; a branch with no slash is ``other``."""
    return branch_path.split("/", 1)[0] if "/" in branch_path else "other"


class History:
    """The first-parent history of ``ref``: its pull requests and which one brought each commit."""

    def __init__(self, git: Git, ref: str) -> None:
        self.git = git
        self.tip = git.run("rev-parse", "--verify", f"{ref}^{{commit}}").strip()
        self.chain: list[Commit] = []
        self.merges: list[Merge] = []
        log = git.run("log", "--first-parent", "--format=%H%x1f%P%x1f%ct%x1f%s%x1f%b%x1e", self.tip)
        for record in reversed(log.split("\x1e")):
            record = record.strip("\n")
            if not record:
                continue
            sha, parents, time, subject, body = record.split("\x1f", 4)
            self.chain.append(Commit(sha, int(time)))
            parent_list = parents.split()
            if len(parent_list) < 2:
                continue
            title = next((line.strip() for line in body.splitlines() if line.strip()), "")
            match = MERGE_SUBJECT.match(subject)
            if match:
                key, branch, kind = f"#{match.group(1)}", f"{match.group(2)}/{match.group(3)}", kind_of(match.group(3))
            else:
                key, branch, kind = f"m{sha[:9]}", subject, "other"
            self.merges.append(Merge(key, branch, kind, sha, parent_list[0], int(time), title))
        self.by_key = {merge.key: merge for merge in self.merges}
        self.writer = self._writers()

    def _writers(self) -> dict[str, str]:
        """Commit -> key of the merge that brought it to the first-parent line (``DIRECT`` if none)."""
        parents: dict[str, list[str]] = {}
        for line in self.git.run("rev-list", "--parents", self.tip).splitlines():
            sha, *rest = line.split()
            parents[sha] = rest
        merge_key = {merge.sha: merge.key for merge in self.merges}
        writer: dict[str, str] = {}
        for commit in self.chain:
            key = merge_key.get(commit.sha, DIRECT)
            writer[commit.sha] = key
            stack = list(parents[commit.sha][1:])
            while stack:
                sha = stack.pop()
                if sha in writer:
                    continue
                writer[sha] = key
                stack.extend(parents.get(sha, []))
        return writer

    def commit_before(self, time: int) -> Commit | None:
        """The last first-parent commit made before ``time``."""
        found = None
        for commit in self.chain:
            if commit.time >= time:
                break
            found = commit
        return found


def parse_diff(text: str) -> tuple[dict[str, list[tuple[int, str]]], dict[str, list[int]]]:
    """``git diff -U0`` -> production files' deleted lines (old number, text) and [added, deleted]."""
    removed: dict[str, list[tuple[int, str]]] = {}
    sizes: dict[str, list[int]] = {}
    current: str | None = None
    old = 0
    in_hunk = False
    for line in text.split("\n"):
        if line.startswith("diff --git "):
            match = re.match(r"diff --git a/(.*) b/(.*)$", line)
            path = match.group(2) if match else None
            current = path if path and is_production(path) else None
            if current:
                removed.setdefault(current, [])
                sizes.setdefault(current, [0, 0])
            in_hunk = False
        elif current is None:
            continue
        elif line.startswith("@@"):
            hunk = HUNK.match(line)
            old = int(hunk.group(1)) if hunk else 0
            in_hunk = True
        elif not in_hunk or line.startswith("\\"):
            continue
        elif line.startswith("-"):
            removed[current].append((old, line[1:]))
            sizes[current][1] += 1
            old += 1
        elif line.startswith("+"):
            sizes[current][0] += 1
    return removed, sizes


def line_ranges(numbers: Iterable[int]) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    for number in sorted(numbers):
        if out and number == out[-1][1] + 1:
            out[-1] = (out[-1][0], number)
        else:
            out.append((number, number))
    return out


@dataclasses.dataclass(frozen=True)
class Rewrite:
    """In ``path``, the fix ``of`` (merged ``age_days`` before) wrote ``lines`` that this change rewrites."""

    path: str
    of: str
    lines: list[int]
    age_days: float


@dataclasses.dataclass(frozen=True)
class ThirdFix:
    path: str
    rewrites: str
    rewrites_lines: list[int]
    which_rewrote: str
    which_rewrote_lines: list[int]
    span_days: float


class Analyser:
    """Diffs and blames each merge once (``change`` may run in several threads) and judges the changes."""

    def __init__(self, history: History) -> None:
        self.history = history
        self.git = history.git
        self._changes: dict[str, Change] = {}

    def change(self, merge: Merge) -> Change:
        cached = self._changes.get(merge.sha)
        if cached is None:
            cached = self._changes[merge.sha] = self.analyse(merge.base, merge.sha)
        return cached

    def analyse(self, base: str, head: str) -> Change:
        """The change from ``base`` to ``head`` (commits or trees), blamed at ``base``."""
        removed, sizes = parse_diff(self.git.run("diff", "-U0", "--no-renames", base, head, "--", "gbdraw/"))
        origins: dict[str, dict[str, list[int]]] = {}
        skipped = []
        for path, lines in removed.items():
            if len(lines) > MAX_REMOVED_LINES:
                skipped.append(path)
                continue
            numbers = [number for number, text in lines if not is_comment_or_blank(path, text)]
            if not numbers:
                continue
            by_writer: dict[str, list[int]] = defaultdict(list)
            for number, sha in self._blame(base, path, numbers).items():
                by_writer[self.history.writer.get(sha, DIRECT)].append(number)
            origins[path] = {key: sorted(found) for key, found in by_writer.items()}
        return Change(
            origins=origins,
            touched=tuple(sorted(sizes)),
            added=sum(size[0] for size in sizes.values()),
            removed=sum(size[1] for size in sizes.values()),
            skipped=tuple(sorted(skipped)),
        )

    def _blame(self, base: str, path: str, numbers: list[int]) -> dict[int, str]:
        args = ["blame", "--porcelain"]
        for first, last in line_ranges(numbers):
            args += ["-L", f"{first},{last}"]
        found: dict[int, str] = {}
        for line in self.git.run(*args, base, "--", path).splitlines():
            match = BLAME_HEADER.match(line)
            if match:
                found[int(match.group(2))] = match.group(1)
        return found


    def recent_fix_rewrites(self, change: Change, time: int) -> list[Rewrite]:
        found = []
        for path, by_writer in sorted(change.origins.items()):
            for key, lines in sorted(by_writer.items()):
                origin = self.history.by_key.get(key)
                if origin is None or origin.kind != "fix" or len(lines) < LINEAGE_LINES:
                    continue
                age = time - origin.time
                if 0 <= age <= WINDOW:
                    found.append(Rewrite(path, key, lines, round(age / DAY, 2)))
        return found

    def third_fixes(self, change: Change, time: int) -> list[ThirdFix]:
        found = []
        for first in self.recent_fix_rewrites(change, time):
            middle = self.history.by_key[first.of]
            for second in self.recent_fix_rewrites(self.change(middle), middle.time):
                oldest = self.history.by_key[second.of]
                if second.path == first.path and time - oldest.time <= WINDOW:
                    found.append(ThirdFix(first.path, first.of, first.lines, second.of, second.lines,
                                          round((time - oldest.time) / DAY, 2)))
        return found

    def rework(self, change: Change, time: int) -> tuple[int, int]:
        """(rewritten lines, those at most WINDOW_DAYS old)."""
        total = recent = 0
        for by_writer in change.origins.values():
            for key, lines in by_writer.items():
                total += len(lines)
                origin = self.history.by_key.get(key)
                if origin is not None and 0 <= time - origin.time <= WINDOW:
                    recent += len(lines)
        return total, recent


def day_start(text: str) -> int:
    return int(dt.datetime.strptime(text, "%Y-%m-%d").replace(tzinfo=JST).timestamp())


def week_start(time: int) -> dt.date:
    day = dt.datetime.fromtimestamp(time, JST).date()
    return day - dt.timedelta(days=day.weekday())


def jst(time: int) -> str:
    return dt.datetime.fromtimestamp(time, JST).strftime("%Y-%m-%d %H:%M")


def weeks_of(since: int, until: int) -> list[str]:
    """The Mondays of the weeks that ``since`` .. ``until`` touches."""
    monday, found = week_start(since), []
    while day_start(monday.isoformat()) < until:
        found.append(monday.isoformat())
        monday += dt.timedelta(days=7)
    return found


# ------------------------------------------------------------------ period report

@dataclasses.dataclass(frozen=True)
class FixEntry:
    """One fix pull request of the period and what it rewrote."""

    pr: str
    merged: str
    week: str
    branch: str
    title: str
    rewritten_lines: int
    recent_lines: int
    net_production_lines: int
    production_files: int
    skipped_files: list[str]
    state_path_title: bool
    layout_title: bool
    rewrites_recent_fix: list[Rewrite]
    third_fix: list[ThirdFix]
    cause_class: str | None = None


@dataclasses.dataclass(frozen=True)
class Registry:
    rows: int
    without_vector: int


@dataclasses.dataclass(frozen=True)
class Row:
    """The counts of one week, or of the whole period (``label`` "total"); the registry and the
    baselines as the last first-parent commit before its end (``end_commit``) has them."""

    label: str
    fix_prs: int
    with_rewritten_lines: int
    rework: int
    rewrite_recent_fix: int
    third_fixes: int
    net_increase: int
    median_net_lines: float | None
    no_production_change: int
    state_path_titles: int
    layout_titles: int
    twin_rows_added: int
    twin_rows_removed: int
    end_commit: str | None
    twins: Registry | None
    error_baselines: dict[str, int | None]
    cause_classes: dict[str, int] | None = None
    tests_runs: int | None = None
    tests_median_minutes: float | None = None


@dataclasses.dataclass(frozen=True)
class TwinChange:
    """The rows of ``TWIN_REGISTRY`` that one first-parent commit added and removed."""

    pr: str
    merged: str
    week: str
    added: list[str]
    removed: list[str]


@dataclasses.dataclass(frozen=True)
class PeriodReport:
    ref: str
    tip: str
    since: str
    until: str
    start_commit: str | None
    weeks: list[Row]
    total: Row
    twins_start: Registry | None
    error_baselines_start: dict[str, int | None]
    twin_changes: list[TwinChange]
    fixes: list[FixEntry]


def summarize(label: str, fixes: list[FixEntry], twins: list[TwinChange], git: Git, end: Commit | None) -> Row:
    rewrote = [fix for fix in fixes if fix.rewritten_lines]
    nets = [fix.net_production_lines for fix in fixes]
    return Row(
        label=label,
        fix_prs=len(fixes),
        with_rewritten_lines=len(rewrote),
        rework=sum(1 for fix in rewrote if fix.recent_lines / fix.rewritten_lines >= REWORK_SHARE),
        rewrite_recent_fix=sum(1 for fix in fixes if fix.rewrites_recent_fix),
        third_fixes=sum(1 for fix in fixes if fix.third_fix),
        net_increase=sum(1 for net in nets if net > 0),
        median_net_lines=statistics.median(nets) if nets else None,
        no_production_change=sum(1 for fix in fixes if not fix.production_files),
        state_path_titles=sum(1 for fix in fixes if fix.state_path_title),
        layout_titles=sum(1 for fix in fixes if fix.layout_title),
        twin_rows_added=sum(len(change.added) for change in twins),
        twin_rows_removed=sum(len(change.removed) for change in twins),
        end_commit=end.sha if end else None,
        twins=registry_at(git, end),
        error_baselines=baselines_at(git, end),
    )


def twin_ids(git: Git, revision: str) -> set[str]:
    text = git.show(revision, TWIN_REGISTRY)
    return set() if text is None else {row["id"] for row in json.loads(text)["twins"]}


def twin_changes(history: History, since: int, until: int) -> list[TwinChange]:
    keys = {merge.sha: merge.key for merge in history.merges}
    found = []
    log = history.git.run("log", "--first-parent", "--format=%H %ct", history.tip, "--", TWIN_REGISTRY)
    for line in reversed(log.splitlines()):
        sha, stamp = line.split()
        time = int(stamp)
        if since <= time < until:
            before, after = twin_ids(history.git, f"{sha}^1"), twin_ids(history.git, sha)
            if before != after:
                found.append(TwinChange(keys.get(sha, DIRECT), jst(time), week_start(time).isoformat(),
                                        sorted(after - before), sorted(before - after)))
    return found


def registry_at(git: Git, commit: Commit | None) -> Registry | None:
    text = git.show(commit.sha, TWIN_REGISTRY) if commit else None
    if text is None:
        return None
    rows = json.loads(text)["twins"]
    return Registry(rows=len(rows), without_vector=sum(1 for row in rows if row.get("vector") is None))


def baselines_at(git: Git, commit: Commit | None) -> dict[str, int | None]:
    values: dict[str, int | None] = {}
    for name, path, pattern in ERROR_BASELINES:
        text = git.show(commit.sha, path) if commit else None
        match = pattern.search(text) if text else None
        value = match.group(1) if match else None
        values[name] = None if value is None else int(value) if value.isdigit() else sum(
            int(entry) for entry in BASELINE_ENTRY.findall(value))
    return values


def period_report(git: Git, ref: str, since: int, until: int, jobs: int) -> PeriodReport:
    history = History(git, ref)
    analyser = Analyser(history)
    fixes = [merge for merge in history.merges if merge.kind == "fix" and since <= merge.time < until]
    # A third fix's middle fix merged at most WINDOW before it: analyse those too, in parallel.
    lineage = [merge for merge in history.merges if merge.kind == "fix" and since - WINDOW <= merge.time < until]
    with ThreadPoolExecutor(max(1, jobs)) as pool:
        list(pool.map(analyser.change, lineage))
    entries = []
    for merge in fixes:
        change = analyser.change(merge)
        rewritten, recent = analyser.rework(change, merge.time)
        entries.append(FixEntry(
            pr=merge.key, merged=jst(merge.time), week=week_start(merge.time).isoformat(), branch=merge.branch,
            title=merge.title, rewritten_lines=rewritten, recent_lines=recent,
            net_production_lines=change.added - change.removed, production_files=len(change.touched),
            skipped_files=list(change.skipped), state_path_title=bool(STATE_PATH_TITLE.search(merge.title)),
            layout_title=bool(LAYOUT_TITLE.search(merge.title)),
            rewrites_recent_fix=analyser.recent_fix_rewrites(change, merge.time),
            third_fix=analyser.third_fixes(change, merge.time),
        ))
    start = history.commit_before(since)
    twins = twin_changes(history, since, until)
    weeks = weeks_of(since, until)
    return PeriodReport(
        ref=ref, tip=history.tip,
        since=dt.datetime.fromtimestamp(since, JST).date().isoformat(),
        until=dt.datetime.fromtimestamp(until - DAY, JST).date().isoformat(),
        start_commit=start.sha if start else None,
        weeks=[summarize(week, [e for e in entries if e.week == week], [c for c in twins if c.week == week], git,
                         history.commit_before(min(day_start(week) + 7 * DAY, until))) for week in weeks],
        total=summarize("total", entries, twins, git, history.commit_before(until)),
        twins_start=registry_at(git, start), error_baselines_start=baselines_at(git, start),
        twin_changes=twins, fixes=entries,
    )


# ------------------------------------------------------------------ GitHub (only with --github)

def github_repository(git: Git) -> str:
    url = git.run("remote", "get-url", "origin").strip()
    match = re.search(r"github\.com[:/]([^/]+/[^/]+?)(?:\.git)?$", url)
    if not match:
        raise RuntimeError(f"origin is not a GitHub repository: {url}")
    return match.group(1)


def gh(*args: str) -> str:
    return subprocess.run(["gh", *args], capture_output=True, check=True, text=True).stdout


def cause_class(body: str) -> str:
    """The class that a pull request body's "Cause class:" line names, ``none``, or ``unrecognised``."""
    match = CAUSE_CLASS_LINE.search(body)
    if not match:
        return "none"
    value = re.split(r"[;,(]", match.group(1).strip().strip("*`").lower(), maxsplit=1)[0]
    value = re.sub(r"[\s_-]+", " ", value).strip()
    return next((name for name in CAUSE_CLASSES if value == name.replace("-", " ")), "unrecognised")


def pull_request_bodies(repository: str, numbers: list[int]) -> dict[int, str]:
    owner, name = repository.split("/")
    bodies: dict[int, str] = {}
    for start in range(0, len(numbers), 50):
        chunk = numbers[start:start + 50]
        fields = " ".join(f"p{number}: pullRequest(number: {number}) {{ body }}" for number in chunk)
        query = f'query {{ repository(owner: "{owner}", name: "{name}") {{ {fields} }} }}'
        data = json.loads(gh("api", "graphql", "-f", f"query={query}"))["data"]["repository"]
        bodies.update({number: (data.get(f"p{number}") or {}).get("body") or "" for number in chunk})
    return bodies


def tests_run_minutes(repository: str, since: int, until: int) -> list[tuple[str, float]]:
    """(week, minutes from start to last update) of each successful pull_request run of Tests."""
    found = []
    for day in range(since, until, DAY):  # a day per query: GitHub lists at most 1,000 runs per filtered query
        span = "..".join(dt.datetime.fromtimestamp(t, dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
                         for t in (day, min(day + DAY, until) - 1))
        text = gh("api", "-X", "GET", f"repos/{repository}/actions/workflows/test.yml/runs", "--paginate",
                  "-f", "event=pull_request", "-f", "status=success", "-f", f"created={span}", "-f", "per_page=100",
                  "--jq", ".workflow_runs[] | [.run_started_at, .updated_at] | @tsv")
        for line in text.splitlines():
            started, updated = (dt.datetime.fromisoformat(value.replace("Z", "+00:00")) for value in line.split("\t"))
            found.append((week_start(int(started.timestamp())).isoformat(), (updated - started).total_seconds() / 60))
    return found


def with_github(report: PeriodReport, git: Git, since: int, until: int) -> PeriodReport:
    repository = github_repository(git)
    # A fix's key is always "#<number>": only a "Merge pull request" subject gives a kind.
    bodies = pull_request_bodies(repository, [int(fix.pr[1:]) for fix in report.fixes])
    fixes = [dataclasses.replace(fix, cause_class=cause_class(bodies[int(fix.pr[1:])])) for fix in report.fixes]
    runs = tests_run_minutes(repository, since, until)

    def add(row: Row, label: str | None) -> Row:
        minutes = [m for week, m in runs if label is None or week == label]
        classes = Counter(str(fix.cause_class) for fix in fixes if label is None or fix.week == label)
        return dataclasses.replace(row, cause_classes=dict(sorted(classes.items())), tests_runs=len(minutes),
                                   tests_median_minutes=round(statistics.median(minutes), 1) if minutes else None)

    return dataclasses.replace(report, fixes=fixes, weeks=[add(row, row.label) for row in report.weeks],
                               total=add(report.total, None))


# ------------------------------------------------------------------ --branch

@dataclasses.dataclass(frozen=True)
class BranchReport:
    branch: str
    kind: str
    ref: str
    tip: str
    base: str
    merged_cleanly: bool
    at: str
    rewritten_lines: int
    recent_lines: int
    rewrites_recent_fix: list[Rewrite]
    third_fix: list[ThirdFix]
    titles: dict[str, str]


def branch_report(git: Git, ref: str, branch: str, at: int) -> BranchReport:
    # `git merge-tree --write-tree` writes the merged tree: keep its objects in a scratch folder that reads
    # the repository's objects as an alternate, so nothing is written to the repository.
    objects = git.run("rev-parse", "--path-format=absolute", "--git-path", "objects").strip()
    with tempfile.TemporaryDirectory(prefix="fix-metrics-") as scratch:
        return _branch_report(Git(git.repo, {"GIT_OBJECT_DIRECTORY": scratch,
                                             "GIT_ALTERNATE_OBJECT_DIRECTORIES": objects}), ref, branch, at)


def _branch_report(git: Git, ref: str, branch: str, at: int) -> BranchReport:
    history = History(git, ref)
    analyser = Analyser(history)
    head = git.run("rev-parse", "--verify", f"{branch}^{{commit}}").strip()
    code, out, err = git.attempt("merge-tree", "--write-tree", "--no-messages", history.tip, head)
    if code not in (0, 1):  # 1 is a conflict
        raise RuntimeError(f"git merge-tree {history.tip} {head}: {err.strip()}")
    if code == 0:
        base, merged = history.tip, out.split()[0]
    else:  # a conflict: judge the branch against the commit it forked from
        base, merged = git.run("merge-base", history.tip, head).strip(), head
    change = analyser.analyse(base, merged)
    rewritten, recent = analyser.rework(change, at)
    rewrites, thirds = analyser.recent_fix_rewrites(change, at), analyser.third_fixes(change, at)
    named = {r.of for r in rewrites} | {t.which_rewrote for t in thirds}
    return BranchReport(
        branch=branch, kind=kind_of(branch.removeprefix("refs/heads/").removeprefix("origin/")), ref=ref,
        tip=history.tip, base=base, merged_cleanly=code == 0, at=jst(at), rewritten_lines=rewritten,
        recent_lines=recent, rewrites_recent_fix=rewrites, third_fix=thirds,
        titles={key: history.by_key[key].title for key in sorted(named)},
    )


def ranges_text(lines: list[int]) -> str:
    return ", ".join(str(a) if a == b else f"{a}-{b}" for a, b in line_ranges(lines))


def branch_text(report: BranchReport) -> str:
    def pr(key: str) -> str:
        return f"{key} ({report.titles[key]})" if report.titles.get(key) else key

    how = "merged with" if report.merged_cleanly else "conflicts with, so judged from its merge base with"
    out = [
        f"{report.branch} ({how} {report.ref} {report.tip[:9]}, at {report.at} JST): {report.rewritten_lines} "
        f"rewritten production lines, {report.recent_lines} of them at most {WINDOW_DAYS} days old.",
        f"Rewrites {LINEAGE_LINES}+ lines of a fix merged in the last {WINDOW_DAYS} days: "
        f"{'yes' if report.rewrites_recent_fix else 'no'}",
    ]
    out += [f"  - {r.path}: lines {ranges_text(r.lines)} of {pr(r.of)}, merged {r.age_days} days ago"
            for r in report.rewrites_recent_fix]
    out.append(f"Third fix in a row in the same file: {'yes' if report.third_fix else 'no'}")
    out += [f"  - {t.path}: rewrites lines {ranges_text(t.rewrites_lines)} of {t.rewrites}, which rewrote lines "
            f"{ranges_text(t.which_rewrote_lines)} of {pr(t.which_rewrote)}; {t.span_days} days from "
            f"{t.which_rewrote} to now" for t in report.third_fix]
    return "\n".join(out) + "\n"


# ------------------------------------------------------------------ Markdown

def _or_dash(value: int | None) -> str:
    return "-" if value is None else str(value)


def count(part: int, whole: int) -> str:
    return f"{part} ({part / whole:.0%})" if whole else str(part)


def period_text(report: PeriodReport) -> str:
    github = report.total.tests_runs is not None
    head = ["week", "fix PRs", "rewrote lines", "rework", "rewrite a recent fix", "third fix",
            "net increase (median lines)", "no production change", "state-path title", "layout title"]
    head += ["Tests median minutes (runs)"] if github else []
    lines = [
        f"# Fix metrics {report.since} .. {report.until} (JST)",
        "",
        f"First-parent merges of `{report.ref}` at `{report.tip[:9]}`; fix PRs are `<owner>/fix/...` branches. "
        "Definitions: `tools/fix_metrics.py`.",
        "",
        "| " + " | ".join(head) + " |",
        "|" + "---|" * len(head),
    ]
    for row in [*report.weeks, report.total]:
        cells = [
            row.label, str(row.fix_prs), str(row.with_rewritten_lines),
            count(row.rework, row.with_rewritten_lines), count(row.rewrite_recent_fix, row.fix_prs),
            count(row.third_fixes, row.fix_prs),
            f"{count(row.net_increase, row.fix_prs)} ({'-' if row.median_net_lines is None else f'{row.median_net_lines:g}'})",
            count(row.no_production_change, row.fix_prs), count(row.state_path_titles, row.fix_prs),
            str(row.layout_titles),
        ]
        cells += [f"{row.tests_median_minutes} ({row.tests_runs})"] if github else []
        lines.append("| " + " | ".join(cells) + " |")
    if report.total.cause_classes is not None:
        lines += ["", "Cause class lines of the fix PRs: "
                  + ", ".join(f"{name} {n}" for name, n in report.total.cause_classes.items()) + "."]
    head = ["at the end of", "commit", f"rows of `{TWIN_REGISTRY}`", "without a shared vector", "rows added / removed",
            *(f"`{path}` baseline" for _, path, _ in ERROR_BASELINES)]
    lines += ["", "| " + " | ".join(head) + " |", "|" + "---|" * len(head)]

    def registry_cells(label: str, commit: str | None, twins: Registry | None, added: str,
                       baselines: dict[str, int | None]) -> str:
        cells = [label, f"`{(commit or '-')[:9]}`", _or_dash(twins.rows if twins else None),
                 _or_dash(twins.without_vector if twins else None), added,
                 *(_or_dash(baselines[name]) for name, _, _ in ERROR_BASELINES)]
        return "| " + " | ".join(cells) + " |"

    lines.append(registry_cells("start", report.start_commit, report.twins_start, "-", report.error_baselines_start))
    lines += [registry_cells(row.label, row.end_commit, row.twins, f"+{row.twin_rows_added} / -{row.twin_rows_removed}",
                             row.error_baselines) for row in [*report.weeks, report.total]]
    lines += ["", "## Twin rows added and removed", ""]
    lines += [f"- {change.pr} {change.merged}: " + "; ".join(
        part for part in (f"+ {', '.join(change.added)}" if change.added else "",
                          f"- {', '.join(change.removed)}" if change.removed else "") if part)
        for change in report.twin_changes]
    lines += ["", "## Third fixes", ""]
    lines += [f"- {fix.pr} {fix.merged} {fix.title}: `{t.path}` rewrites {len(t.rewrites_lines)} lines of "
              f"{t.rewrites}, which rewrote {len(t.which_rewrote_lines)} lines of {t.which_rewrote} "
              f"({t.span_days} days)" for fix in report.fixes for t in fix.third_fix]
    lines += ["", f"## Rewrites of a fix merged in the previous {WINDOW_DAYS} days", ""]
    lines += [f"- {fix.pr} {fix.merged} {fix.title}: " + "; ".join(
        f"`{r.path}` {len(r.lines)} lines of {r.of} ({r.age_days} days)" for r in fix.rewrites_recent_fix)
        for fix in report.fixes if fix.rewrites_recent_fix]
    return "\n".join(lines) + "\n"


def parse_args(argv: list[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=(__doc__ or "").split("\n\n")[0])
    parser.add_argument("--repo", default=".", help="the git checkout to read (default: the current directory)")
    parser.add_argument("--ref", default="origin/dev", help="the branch whose first-parent history is read")
    parser.add_argument("--week", help="a Monday (YYYY-MM-DD, JST); default: the last full week")
    parser.add_argument("--since", help="first day of a range (YYYY-MM-DD, JST)")
    parser.add_argument("--until", help="last day of a range, included (YYYY-MM-DD, JST)")
    parser.add_argument("--branch", help="judge this work branch instead of reporting on a period")
    parser.add_argument("--at", help="judge --branch at this time (YYYY-MM-DDTHH:MM, JST); default: now")
    parser.add_argument("--github", action="store_true",
                        help="also read the fix PRs' cause class lines and the Tests run times with gh")
    parser.add_argument("--json", type=Path, help="also write the report as JSON to this file")
    parser.add_argument("--jobs", type=int, default=4, help="git processes at once (default 4)")
    args = parser.parse_args(argv)
    if args.branch and (args.week or args.since or args.until or args.github):
        parser.error("--branch takes no period and no --github")
    if bool(args.since) != bool(args.until) or (args.since and args.week):
        parser.error("give --since and --until together, or --week")
    if args.week and dt.date.fromisoformat(args.week).weekday() != 0:
        parser.error(f"--week {args.week} is not a Monday")
    return args


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    git = Git(Path(args.repo))
    now = int(dt.datetime.now(JST).timestamp())
    report: PeriodReport | BranchReport
    if args.branch:
        at = int(dt.datetime.strptime(args.at, "%Y-%m-%dT%H:%M").replace(tzinfo=JST).timestamp()) if args.at else now
        report = branch_report(git, args.ref, args.branch, at)
        text = branch_text(report)
    else:
        if args.since:
            since, until = day_start(args.since), day_start(args.until) + DAY
        else:
            monday = dt.date.fromisoformat(args.week) if args.week else week_start(now) - dt.timedelta(days=7)
            since = day_start(monday.isoformat())
            until = since + 7 * DAY
        report = period_report(git, args.ref, since, until, args.jobs)
        if args.github:
            report = with_github(report, git, since, until)
        text = period_text(report)
    if args.json:
        args.json.write_text(json.dumps(dataclasses.asdict(report), indent=1, ensure_ascii=False) + "\n",
                             encoding="utf-8")
    sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
