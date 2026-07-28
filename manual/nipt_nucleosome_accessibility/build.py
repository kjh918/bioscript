#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from pathlib import Path

def create_structure():
    """Creates the directory structure and initial files for the pipeline."""
    base_dir = Path("nipt_pipeline")
    
    structure = {
        "src": {
            "__init__.py": "",
            "core": {
                "__init__.py": "",
                "fragment.py": "# Fragment class definition\n",
                "utils.py": "# BAM and BED parsing utilities\n"
            },
            "accessibility": {
                "__init__.py": "",
                "calculator.py": "# Accessibility calculation logic\n"
            },
            "cnv": {
                "__init__.py": "",
                "calculator.py": "# CNV Binning, RD/FS calculation, and Z-score logic\n",
                "normalization.py": "# GC-bias and baseline normalization logic\n"
            }
        },
        "main.py": "# Main CLI entry point\n",
        "requirements.txt": "pysam>=0.21.0\nnumpy\nscipy\n",
        "README.md": "# NIPT Pipeline\n\nCombined Accessibility and CNV analysis.\n"
    }

    def build_tree(current_path, tree):
        for name, content in tree.items():
            path = current_path / name
            if isinstance(content, dict):
                path.mkdir(parents=True, exist_ok=True)
                build_tree(path, content)
            else:
                if not path.exists():
                    with open(path, "w") as f:
                        f.write(content)
                    print(f"Created: {path}")

    base_dir.mkdir(parents=True, exist_ok=True)
    build_tree(base_dir, structure)
    print(f"\nProject structure successfully created in './{base_dir}'")

if __name__ == "__main__":
    create_structure()