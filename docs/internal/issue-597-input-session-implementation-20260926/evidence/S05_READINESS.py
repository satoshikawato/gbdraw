"""Read-only S05 authority evidence, not a policy checker or runtime permission."""

import argparse
import hashlib
import json
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
PLAN = "docs/internal/issue-597-input-session-implementation-20260926"
CONTRACT = "docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md"
AUTHORITY_MERGE = "af5d942af60353dda199aa487da9152a3576b3fe"
PREREQUISITES = {
    "S01": "15bcbca89392cbd63a88fe80dc44b71ad4868061",
    "S03": "e3684d496df7a3b2a7492443cad50168f6959b36",
    "S04Integration": "eed70f05d04f638a5a4ab6c6567982d45e774596",
    "S04": "b092499d27705092d84037d5acd51032177c619a",
}
FIELDS = {
    "Concern": "concern", "Scenario revision": "scenarioRevision", "Choice": "choice",
    "Rationale": "rationale", "Must preserve": "mustPreserve", "May retire": "mayRetire",
    "Accepted residual risk": "acceptedResidualRisk", "Owner": "owner",
    "Decision date": "decisionDate",
}


def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT)


def blob(ref, path):
    return git("show", f"{ref}:{path}")


def digest(value):
    return hashlib.sha256(value).hexdigest()


def ancestor(commit, ref):
    result = subprocess.run(["git", "merge-base", "--is-ancestor", commit, ref], cwd=ROOT)
    assert result.returncode in (0, 1)
    return result.returncode == 0


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--head", default="HEAD")
parser.add_argument("--dev", default="origin/dev")
args = parser.parse_args()
head = git("rev-parse", args.head).decode().strip()
dev = git("rev-parse", args.dev).decode().strip()
ancestry = {name: {"sha": sha, "headAncestor": ancestor(sha, head)}
            for name, sha in PREREQUISITES.items()}
assert all(row["headAncestor"] for row in ancestry.values())
assert ancestor(AUTHORITY_MERGE, head) and ancestor(AUTHORITY_MERGE, dev)
contracts = {
    "historical": blob(AUTHORITY_MERGE, CONTRACT),
    "head": blob(head, CONTRACT),
    "dev": blob(dev, CONTRACT),
}
assert contracts["head"] == contracts["dev"]
receipts = []
for decision_id, name in (("PD-OI-044", "02_RECORD_DISCOVERY.md"),
                          ("PD-OI-045", "03_SESSION_OPERATIONS.md")):
    receipt_bytes = blob(head, f"{PLAN}/decisions/{name}")
    text = receipt_bytes.decode()
    human = re.search(r"```text\nPRODUCT_DECISION\n(.*?)\n```", text, re.S)[1]
    receipt = {
        FIELDS[key]: int(value) if key == "Scenario revision" else value
        for key, value in (line.split(": ", 1) for line in human.splitlines())
    }
    serialized = json.loads(re.search(r"```json\n(.*?)\n```", text, re.S)[1])
    assert len(receipt) == 9 and receipt == serialized
    for contract in contracts.values():
        section = re.search(
            rf"^### {decision_id}:.*?(?=^### PD-OI-|^## Acceptance contract catalog|\Z)",
            contract.decode(), re.S | re.M,
        )[0]
        authority = json.loads(re.search(r"```json\n(.*?)\n```", section, re.S)[1])
        assert receipt == authority
    receipts.append({"id": decision_id, "fieldsEqual": list(receipt),
                     "receiptSha256": digest(receipt_bytes)})

policy_path = "tools/web-change-policy.json"
policy_bytes = blob(dev, policy_path)
policy = json.loads(policy_bytes)
permissions = [
    {"kind": "constructor", "capability": "Diagram Worker",
     "subject": "services/session-import-client.js",
     "present": "services/session-import-client.js"
     in policy["allowedPrivilegedOwners"]["Diagram Worker"]},
    {"kind": "importer", "target": "services/session-file.js",
     "subject": "workers/session-import-worker.js",
     "present": "workers/session-import-worker.js"
     in policy["allowedPrivilegedImporters"]["services/session-file.js"]},
]
map_path = "tools/web-product-impact-map.json"
mapped = json.loads(blob(dev, map_path))
mapped_evidence = []
for concern in mapped["concerns"]:
    for contract in concern["contracts"]:
        path = contract["ref"].split("::", 1)[0]
        base_bytes = blob(dev, path)
        mapped_evidence.append({
            "concern": concern["key"], "ref": contract["ref"],
            "execution": contract["execution"],
            "devSha256": digest(base_bytes),
            "headBytesEqualDev": blob(head, path) == base_bytes,
        })
missing = [row for row in permissions if not row["present"]]
print(json.dumps({
    "scope": "S05 prerequisite evidence only; no runtime admission",
    "head": head, "dev": dev, "prerequisites": ancestry,
    "productAuthorityMerge": AUTHORITY_MERGE,
    "productAuthorityAncestry": {"head": True, "dev": True},
    "contractRevisions": {
        name: int(re.search(rb"Contract revision: `(\d+)`", value)[1])
        for name, value in contracts.items()
    },
    "contractBytesEqualLatestDev": True,
    "contractSha256": digest(contracts["dev"]), "receipts": receipts,
    "permissions": permissions, "missingPermissions": missing,
    "policyBytesEqualLatestDev": blob(head, policy_path) == policy_bytes,
    "policySha256": digest(policy_bytes),
    "productMapBytesEqualLatestDev": blob(head, map_path) == blob(dev, map_path),
    "mappedEvidence": mapped_evidence,
    "readiness": ("BLOCKED_MISSING_PRIVILEGED_PERMISSION" if missing
                  else "PERMISSION_SUBJECTS_PRESENT_ONLY; verify mapped evidence and integrate dev"),
    "runtimeStarted": False,
}, indent=2))
raise SystemExit(2 if missing else 0)
