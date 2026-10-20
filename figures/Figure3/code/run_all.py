"""Render every Figure 3 panel from ../data/ into ../panels/.

    python code/run_all.py        (from the Figure3 folder, or from anywhere)

Each panel script is run in its own Python process, exactly as if run by hand.
prep_data.py is NOT run here (it needs upstream result tables and local bigWig tracks);
run it only to regenerate ../data/.
"""
import os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = ["panel_A.py", "panel_B.py", "panel_C.py", "legend_key.py", "panel_D.py", "panel_E.py", "panel_F.py", "panel_G.py"]

if __name__ == "__main__":
    for s in SCRIPTS:
        print(f"running {s}", flush=True)
        subprocess.run([sys.executable, os.path.join(HERE, s)], cwd=HERE, check=True, env=dict(os.environ, PYTHONPATH=HERE))
    print("done:", os.path.normpath(os.path.join(HERE, "..", "panels")))
