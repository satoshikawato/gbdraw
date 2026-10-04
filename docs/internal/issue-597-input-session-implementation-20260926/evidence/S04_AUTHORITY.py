"""Reproduce S04 receipt evidence; do not modify or replace a policy checker."""

import hashlib
import json
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
PLAN = ROOT / "docs/internal/issue-597-input-session-implementation-20260926"
MERGE = "af5d942af60353dda199aa487da9152a3576b3fe"
BASE = "c922fc38aac78da9be83342c09ac0164ecef6ff6"
CONTRACT = "docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md"
FIELDS = {
    "Concern": "concern", "Scenario revision": "scenarioRevision", "Choice": "choice",
    "Rationale": "rationale", "Must preserve": "mustPreserve", "May retire": "mayRetire",
    "Accepted residual risk": "acceptedResidualRisk", "Owner": "owner",
    "Decision date": "decisionDate",
}


def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT)


for target in (BASE, "HEAD"):
    subprocess.run(["git", "merge-base", "--is-ancestor", MERGE, target], cwd=ROOT, check=True)
current = (ROOT / CONTRACT).read_bytes()
assert current == git("show", f"{BASE}:{CONTRACT}")
historical = git("show", f"{MERGE}:{CONTRACT}")
assert int(re.search(rb"Contract revision: `(\d+)`", historical)[1]) == 21
assert int(re.search(rb"Contract revision: `(\d+)`", current)[1]) == 22
receipts = []
for decision_id, filename in (
    ("PD-OI-044", "02_RECORD_DISCOVERY.md"),
    ("PD-OI-045", "03_SESSION_OPERATIONS.md"),
):
    text = (PLAN / "decisions" / filename).read_text()
    human = re.search(r"```text\nPRODUCT_DECISION\n(.*?)\n```", text, re.S)[1]
    receipt = {
        FIELDS[key]: int(value) if key == "Scenario revision" else value
        for key, value in (line.split(": ", 1) for line in human.splitlines())
    }
    serialized = json.loads(re.search(r"```json\n(.*?)\n```", text, re.S)[1])
    assert len(receipt) == 9 and receipt == serialized
    for contract in (historical.decode(), current.decode()):
        section = re.search(
            rf"^### {decision_id}:.*?(?=^### PD-OI-|^## Acceptance contract catalog|\Z)",
            contract, re.S | re.M,
        )[0]
        authority = json.loads(re.search(r"```json\n(.*?)\n```", section, re.S)[1])
        assert receipt == authority
    receipts.append({
        "id": decision_id, "fieldsEqual": list(receipt),
        "receiptSha256": hashlib.sha256(text.encode()).hexdigest(),
    })
print(json.dumps({
    "head": git("rev-parse", "HEAD").decode().strip(), "base": BASE,
    "authorityMerge": MERGE, "historicalRevision": 21, "currentRevision": 22,
    "baseContractBytesEqual": True, "receipts": receipts,
    "contractSha256": hashlib.sha256(current).hexdigest(),
}, indent=2))
