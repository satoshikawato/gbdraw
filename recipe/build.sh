#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# --- 1. Build and sync the browser wheel artifact used by the web UI ---
echo "Building outer wheel for browser-wheel sync..."
WHEEL_BUILD_DIR=$(mktemp -d "${TMPDIR:-/tmp}/gbdraw-build-XXXXXX")
trap 'rm -rf "$WHEEL_BUILD_DIR"' EXIT

$PYTHON -m build --wheel --no-isolation --outdir "$WHEEL_BUILD_DIR"

OUTER_WHEEL=$(find "$WHEEL_BUILD_DIR" -maxdepth 1 -name 'gbdraw-*.whl' | head -n 1)
if [ -z "$OUTER_WHEEL" ]; then
    echo "Could not find built outer wheel in $WHEEL_BUILD_DIR" >&2
    exit 1
fi

echo "Synchronizing browser wheel from $OUTER_WHEEL..."
$PYTHON tools/sync_browser_wheel.py "$OUTER_WHEEL"

# --- 2. Standard installation (for CLI usage) ---
$PYTHON -m pip install . --no-deps --ignore-installed -vv

# --- 3. Copy font files ---
mkdir -p $PREFIX/fonts
cp gbdraw/data/*.ttf $PREFIX/fonts/

# --- 4. Force copy the synced web directory to site-packages ---
# Ensure the web directory (including index.html and the wheel) is copied to the installation path
echo "Copying web directory to site-packages..."
mkdir -p $SP_DIR/gbdraw/web
cp -r gbdraw/web/* $SP_DIR/gbdraw/web/
