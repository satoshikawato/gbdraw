#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# --- 1. Prepare the browser wheel bundled with the offline GUI ---
echo "Preparing browser wheel..."
$PYTHON tools/prepare_browser_wheel.py

# --- 2. Standard installation (for CLI usage) ---
$PYTHON -m pip install . --no-deps --ignore-installed -vv

# --- 3. Copy font files ---
mkdir -p $PREFIX/fonts
cp gbdraw/data/*.ttf $PREFIX/fonts/

# --- 4. Force copy local web runtime assets to site-packages ---
echo "Copying local web runtime assets to site-packages..."
mkdir -p $SP_DIR/gbdraw/web
cp gbdraw/web/index.html $SP_DIR/gbdraw/web/
cp gbdraw/web/open-source-notices.html $SP_DIR/gbdraw/web/
cp gbdraw/web/gbdraw-*.whl $SP_DIR/gbdraw/web/
for web_asset_dir in assets js presets vendor wasm; do
    cp -r gbdraw/web/$web_asset_dir $SP_DIR/gbdraw/web/
done
