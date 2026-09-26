"""Emit pinned, read-only S00 source/authority/ref evidence from this clone."""
import hashlib
import importlib.metadata
import json
from pathlib import Path
import re
import subprocess
import sys

PLAN = Path("docs/internal/issue-601-bug15-bug19-implementation-20260926")
START = "be4000ede35ea5692d7a88d46c0130ffec992f5f"
BASE = "af5d942af60353dda199aa487da9152a3576b3fe"
OLD = "2edc00aebc74e01003da643dfc957b513d5dcfe5"
BRANCHES = {
    "fix/issue-601-bug15-bug19": START,
    "fix/issue-601-export-output-20260926": "43836eb77798924bc826d0064d3cc5e219379d91",
    "fix/issue-598-alignment-direction-reset-20260926": "d5edc4bb00b0d793d94361202c31a270f78a4aa0",
    "fix/issue-602-linear-live-edit-20260926": "94faf8a98eddf823f5ebfed859d324e258daa689",
}
AUTHORITY = ["docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md", "tools/web-product-impact-map.json",
             "tools/web-product-decisions.json", "docs/internal/PRODUCT_IMPACT_RATCHET.md",
             "docs/internal/WEB_CHANGE_POLICY.md", "docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md"]
SOURCES = ["gbdraw/exceptions.py", "gbdraw/render/groups/linear/pairwise_match.py",
           "gbdraw/web/js/services/error-normalization.js", "gbdraw/web/js/app/python-helpers.js",
           "gbdraw/web/js/workers/diagram-generation-worker.js", "gbdraw/web/js/services/diagram-generation.js",
           "gbdraw/web/js/app/run-analysis.js", "gbdraw/web/js/app/similarity-alignment.js",
           "gbdraw/web/js/app/app-setup.js", "gbdraw/web/index.html", "gbdraw/web/js/app/rule-matching.js",
           "gbdraw/web_support/rule_matching.py", "gbdraw/web/js/app/feature-editor/rule-actions.js",
           "gbdraw/web/js/app/feature-editor/label-actions.js", "gbdraw/web/js/app/feature-editor/label-override-table.js",
           "gbdraw/web/js/app/file-imports.js", "gbdraw/web/js/app/watchers.js", "gbdraw/web/js/app/feature-visibility.js",
           "gbdraw/web/js/app/feature-editor/visibility-actions.js", "gbdraw/web/js/services/export.js",
           "gbdraw/web/js/services/pdf-fonts.js", "gbdraw/web/js/app/feature-search/search-core.js",
           "gbdraw/web/js/services/standalone-interactivity-assets.js", "gbdraw/features/colors.py",
           "gbdraw/labels/filtering.py", "gbdraw/features/visibility.py", "gbdraw/web_support/request_render.py"]


def git(*args):
    return subprocess.check_output(["git", *args], text=True).strip()


def digest(ref, path):
    return hashlib.sha256(subprocess.check_output(["git", "show", f"{ref}:{path}"])).hexdigest()


fields = dict(zip(["Concern", "Scenario revision", "Choice", "Rationale", "Must preserve", "May retire", "Accepted residual risk", "Owner", "Decision date"],
                  ["concern", "scenarioRevision", "choice", "rationale", "mustPreserve", "mayRetire", "acceptedResidualRisk", "owner", "decisionDate"]))
receipts = []
for filename in ["DECISION_01_ERROR_DISCLOSURE.md", "DECISION_02_REGEX_EDIT_RECOVERY.md"]:
    text = (PLAN / "decisions" / filename).read_text()
    receipt = re.search(r"~~~text\n(PRODUCT_DECISION\n.*?)\n~~~", text, re.S).group(1)
    expected = {fields[k]: int(v) if k == "Scenario revision" else v for k, v in (line.split(": ", 1) for line in receipt.splitlines()[1:])}
    actual = json.loads(re.search(r"~~~json\n(.*?)\n~~~", text, re.S).group(1))
    assert expected == actual and len(actual) == 9
    sha = hashlib.sha256(receipt.encode()).hexdigest()
    assert sha in text
    receipts.append({"file": filename, "receiptSha256": sha, "fieldsEqual": 9, "concern": actual["concern"], "choice": actual["choice"]})

contract = git("show", f"{BASE}:{AUTHORITY[0]}")
pd_revisions = []
for section in re.split(r"(?=^### PD-OI-)", contract, flags=re.M):
    match = re.match(r"### (PD-OI-\d+):", section)
    if match:
        revision = re.search(r"- Scenario revision: `(\d+)`", section)
        pd_revisions.append({"id": match[1], "scenarioRevision": int(revision[1]) if revision else None})
decisions = json.loads(git("show", f"{BASE}:{AUTHORITY[2]}"))
assert "satoshikawato" in decisions["maintainerLogins"]
product_map = json.loads(git("show", f"{BASE}:{AUTHORITY[1]}"))
inventory = {
    "startSha": START, "inspectedDevSha": BASE, "priorInvestigationSha": OLD,
    "python": sys.version.split()[0], "node": subprocess.check_output(["node", "--version"], text=True).strip(),
    "packages": {name: importlib.metadata.version(name) for name in ["biopython", "pandas", "svgwrite", "pytest", "playwright"]},
    "unchangedTrees": {path: {ref: git("rev-parse", f"{ref}:{path}") for ref in [OLD, BASE, START]} for path in ["gbdraw", "tests", "tools", ".github"]},
    "sourceSha256": {path: digest(START, path) for path in SOURCES},
    "authoritySha256": {path: digest(BASE, path) for path in AUTHORITY},
    "contractRevision": int(re.search(r"Contract revision: `(\d+)`", contract)[1]),
    "pdRevisions": pd_revisions, "activeBdDecisions": decisions["decisions"],
    "mapConcernKeys": re.findall(r'"key":\s*"([^"]+)"', json.dumps(product_map)),
    "receiptChecks": receipts,
    "disclosureConcernsPresentOnDev": {key: key in contract for key in ["web.errors.diagnostic-disclosure", "web.errors.user-facing-diagnostic-disclosure", "web.rules.rejected-pattern-edit-recovery"]},
    "relatedBranches": {},
    "runAnalysisErrorReturns": [{"line": i, "source": line.strip()} for i, line in enumerate(git("show", f"{START}:{SOURCES[6]}").splitlines(), 1) if re.search(r"return.*status: 'error'", line)],
    "originalAuditHistory": git("log", "--all", "--format=%H", "--", "docs/internal/GUI_AUDIT_DEV_20260926.md").splitlines(),
}
for branch, sha in BRANCHES.items():
    files = git("ls-tree", "-r", "--name-only", sha, "docs/internal").splitlines()
    results = [p for p in files if re.search(r"issue-(598|601|602).*(?:RESULT\.md|/(?:results|evidence)/S\d\d\.md)$", p)]
    inventory["relatedBranches"][branch] = {
        "sha": sha, "mergeBaseWithDev": git("merge-base", BASE, sha),
        "unmergedRuntimePaths": git("diff", "--name-only", f"{BASE}...{sha}", "--", "gbdraw", "tests", "tools", ".github").splitlines(),
        "isAncestorOfDev": subprocess.run(["git", "merge-base", "--is-ancestor", sha, BASE]).returncode == 0,
        "results": {path: {"sha256": digest(sha, path), "lastCommit": git("log", "-1", "--format=%H", sha, "--", path)} for path in results},
    }
assert all(len(set(trees.values())) == 1 for trees in inventory["unchangedTrees"].values())
assert len(inventory["runAnalysisErrorReturns"]) == 9
print(json.dumps(inventory, ensure_ascii=False, indent=2))
