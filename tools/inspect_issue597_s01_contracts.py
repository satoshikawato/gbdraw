"""Verify merged issue 597 Product receipts without modifying authority."""

import hashlib
import json
import re
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
plan = root / "docs/internal/issue-597-input-session-implementation-20260926"
contract = (root / "docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md").read_text()
fields = {
    "Concern": "concern",
    "Scenario revision": "scenarioRevision",
    "Choice": "choice",
    "Rationale": "rationale",
    "Must preserve": "mustPreserve",
    "May retire": "mayRetire",
    "Accepted residual risk": "acceptedResidualRisk",
    "Owner": "owner",
    "Decision date": "decisionDate",
}
result = {
    "head": subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=root, text=True
    ).strip(),
    "dev": subprocess.check_output(
        ["git", "rev-parse", "origin/dev"], cwd=root, text=True
    ).strip(),
    "mergeSha": "af5d942af60353dda199aa487da9152a3576b3fe",
    "revision": int(re.search(r"Contract revision: `(\d+)`", contract)[1]),
    "receipts": [],
}
assert result["revision"] == 21
for ref in ["HEAD", "origin/dev"]:
    subprocess.run(
        ["git", "merge-base", "--is-ancestor", result["mergeSha"], ref],
        cwd=root,
        check=True,
    )
assert (
    subprocess.check_output(
        [
            "git",
            "show",
            "origin/dev:docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md",
        ],
        cwd=root,
    )
    == (root / "docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md").read_bytes()
)
for id, name in [
    ("PD-OI-044", "02_RECORD_DISCOVERY.md"),
    ("PD-OI-045", "03_SESSION_OPERATIONS.md"),
]:
    text = (plan / "decisions" / name).read_text()
    human = re.search(r"```text\nPRODUCT_DECISION\n(.*?)\n```", text, re.S)[1]
    receipt = {
        fields[k]: int(v) if k == "Scenario revision" else v
        for k, v in (line.split(": ", 1) for line in human.splitlines())
    }
    serialized = json.loads(re.search(r"```json\n(.*?)\n```", text, re.S)[1])
    section = re.search(
        r"^### " + id + r":.*?(?=^### PD-OI-|^## Acceptance contract catalog|\Z)",
        contract,
        re.S | re.M,
    )[0]
    authority = json.loads(re.search(r"```json\n(.*?)\n```", section, re.S)[1])
    assert len(receipt) == 9 and receipt == serialized == authority
    result["receipts"].append(
        {
            "id": id,
            "fieldsEqual": list(receipt),
            "sha256": hashlib.sha256(text.encode()).hexdigest(),
            "approvedOutcomeDigest": re.search(
                r"承認対象資料の SHA-256: `([a-f0-9]{64})`", text
            )[1],
        }
    )
print(json.dumps(result, indent=2))
