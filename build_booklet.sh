#!/bin/bash
# Compile the CASSIA Instruction Booklet and copy the PDF to the repo root.
# Requires the Instruction_Booklet/ directory (Overleaf-linked git, not tracked here).
set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BOOKLET_DIR="$SCRIPT_DIR/Instruction_Booklet"

if [ ! -d "$BOOKLET_DIR" ]; then
    echo "ERROR: Instruction_Booklet/ not found. Clone it from Overleaf first." >&2
    exit 1
fi

cd "$BOOKLET_DIR"
latexmk -pdf -f -interaction=nonstopmode main.tex
cp main.pdf "$SCRIPT_DIR/CASSIA_Instruction_Booklet.pdf"
echo "PDF written to: $SCRIPT_DIR/CASSIA_Instruction_Booklet.pdf"
