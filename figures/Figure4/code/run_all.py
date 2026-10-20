"""Rebuild every Figure 4 panel from ../data/ into ../panels/.   Usage:  python code/run_all.py
(prep_data.py, which regenerates ../data/ from the raw inputs, is NOT run here; run it separately if needed.)"""
import os, subprocess, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
for s in sorted(HERE.glob('panel_*.py')):
    print(f'== {s.name}', flush=True)
    subprocess.run([os.path.abspath(sys.executable), str(s)], cwd=HERE, check=True)
