#!/usr/bin/env bash
set -uo pipefail

# Claude Code Stop hook.
# It blocks one completion attempt when relevant Web files changed and the
# repository guard command fails. On a repeated stop-hook continuation it reports
# the failure without creating an infinite loop; CI and pre-push remain the final
# enforcement boundaries.

INPUT="$(cat)"
PROJECT_DIR="${CLAUDE_PROJECT_DIR:-}"

if [[ -z "$PROJECT_DIR" ]]; then
  PROJECT_DIR="$(python - "$INPUT" <<'PY'
import json
import sys
try:
    payload = json.loads(sys.argv[1])
except Exception:
    payload = {}
print(payload.get("cwd", ""))
PY
)"
fi

if [[ -z "$PROJECT_DIR" || ! -d "$PROJECT_DIR/.git" ]]; then
  exit 0
fi

STOP_HOOK_ACTIVE="$(python - "$INPUT" <<'PY'
import json
import sys
try:
    payload = json.loads(sys.argv[1])
except Exception:
    payload = {}
print("true" if payload.get("stop_hook_active") else "false")
PY
)"

cd "$PROJECT_DIR" || exit 0

CHANGED="$(
  {
    git diff --name-only --cached
    git diff --name-only
    git ls-files --others --exclude-standard
  } 2>/dev/null | sort -u
)"

if ! printf '%s\n' "$CHANGED" | grep -Eq \
  '^(gbdraw/web/js/|tests/web/|AGENTS\.md$|CLAUDE\.md$|gbdraw/web/CLAUDE\.md$|\.agents/skills/refactor-gbdraw-web-safely/|\.claude/hooks/|\.claude/settings\.json$|package\.json$|package-lock\.json$|\.github/workflows/web-refactor-guardrails\.yml$)'; then
  exit 0
fi

OUTPUT_FILE="$(mktemp)"
trap 'rm -f "$OUTPUT_FILE"' EXIT

if npm run test:web:refactor-guards >"$OUTPUT_FILE" 2>&1; then
  exit 0
fi

python - "$OUTPUT_FILE" "$STOP_HOOK_ACTIVE" <<'PY'
import json
import pathlib
import sys

path = pathlib.Path(sys.argv[1])
stop_hook_active = sys.argv[2] == "true"
output = path.read_text(encoding="utf-8", errors="replace")
tail = output[-7000:]
message = (
    "Web refactor guardrails failed. Fix the failures before declaring the task "
    "complete.\n\n" + tail
)

if stop_hook_active:
    # Avoid an unbounded Stop-hook loop. The failure is still shown to Claude and
    # the user; pre-push and CI remain blocking enforcement.
    print(json.dumps({
        "hookSpecificOutput": {
            "hookEventName": "Stop",
            "additionalContext": message
        }
    }))
else:
    print(json.dumps({
        "decision": "block",
        "reason": message
    }))
PY
