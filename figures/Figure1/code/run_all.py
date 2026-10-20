"""Rebuild every Figure 1 panel from ../data/ :  python code/run_all.py   (from the Figure1/ folder or anywhere)
Requires only python 3 + pandas, numpy, matplotlib (scipy not needed). Data tables are produced by prep_data.py."""
import os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = [f"panel_{p}.py" for p in "ABCDEFGHIJKL"] + ["panel_JKL_track_labels.py"]
for s in SCRIPTS:
    r = subprocess.run([sys.executable, os.path.join(HERE, s)], cwd=HERE)
    if r.returncode: sys.exit(f"{s} failed")
print(f"{len(SCRIPTS)} panel files written to", os.path.normpath(os.path.join(HERE, "..", "panels")))
