#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# --- 1. Build the wheel file for Pyodide (WebAssembly) ---
echo "Building wheel for Pyodide..."
# Build the wheel from the current directory and output it to gbdraw/web/
$PYTHON -m pip wheel . --no-deps --no-build-isolation --wheel-dir gbdraw/web

# --- 2. Dynamically retrieve the filename and update index.html ---
# Find the generated wheel file (e.g., gbdraw-0.8.0-py3-none-any.whl)
WHEEL_FILE=$(ls gbdraw/web/gbdraw-*.whl | head -n 1)
WHEEL_NAME=$(basename $WHEEL_FILE)
echo "Generated wheel name: $WHEEL_NAME"

# Update index.html to use the correct wheel filename
# This replaces 'const GBDRAW_WHEEL_NAME = "...";' with the actual filename
sed -i "s/const GBDRAW_WHEEL_NAME = \".*\";/const GBDRAW_WHEEL_NAME = \"$WHEEL_NAME\";/" gbdraw/web/index.html

# --- 3. Standard installation (for CLI usage) ---
$PYTHON -m pip install . --no-deps --ignore-installed -vv

# --- 4. Copy font files ---
mkdir -p $PREFIX/fonts
cp gbdraw/data/*.ttf $PREFIX/fonts/

# --- 5. Force copy the web directory to site-packages ---
# Ensure the web directory (including index.html and the wheel) is copied to the installation path
echo "Copying web directory to site-packages..."
mkdir -p $SP_DIR/gbdraw/web
cp -r gbdraw/web/* $SP_DIR/gbdraw/web/
