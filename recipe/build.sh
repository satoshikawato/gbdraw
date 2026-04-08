#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# --- 1. Synchronize the browser wheel bundled with the offline GUI ---
echo "Synchronizing browser wheel..."
$PYTHON setup.py build_py

# --- 2. Standard installation (for CLI usage) ---
$PYTHON -m pip install . --no-deps --ignore-installed -vv

# --- 3. Copy font files ---
mkdir -p $PREFIX/fonts
cp gbdraw/data/*.ttf $PREFIX/fonts/

# --- 4. Force copy the web directory to site-packages ---
# Ensure the web directory (including index.html and the wheel) is copied to the installation path
echo "Copying web directory to site-packages..."
mkdir -p $SP_DIR/gbdraw/web
cp -r gbdraw/web/* $SP_DIR/gbdraw/web/
