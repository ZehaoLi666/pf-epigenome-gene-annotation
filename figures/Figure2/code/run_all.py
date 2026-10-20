"""Rebuild every Figure 2 panel from ../data (plain python; numpy/pandas/matplotlib)."""
import os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))
for s in ["panel_A.py", "panel_B.py", "panel_C.py"]:
    print("running", s, flush=True)
    subprocess.run([os.path.realpath(sys.executable), os.path.join(HERE, s)], cwd=HERE, check=True)
print("done ->", os.path.normpath(os.path.join(HERE, "..", "panels")))
