#!/usr/bin/env bash
set -euo pipefail

ROOT="$(git rev-parse --show-toplevel)"
cd "$ROOT"

chmod +x .githooks/pre-push
chmod +x .claude/hooks/run-web-refactor-guards.sh
git config core.hooksPath .githooks

echo "Installed repository Git hooks from .githooks/"
echo "Claude Code project hook is configured in .claude/settings.json"
