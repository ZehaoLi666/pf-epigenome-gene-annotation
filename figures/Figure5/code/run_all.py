"""Re-draw every Figure 5 panel from ../data (run prep_data.py first only if inputs change)."""
import os, runpy
HERE = os.path.dirname(os.path.abspath(__file__))
for s in ["panel_A.py", "panel_B.py", "panel_C.py", "panel_D.py", "legend_AD.py", "panel_E.py", "panel_F.py",
          "panel_G.py", "panel_H.py", "panel_I.py", "panel_J.py"]:
    runpy.run_path(os.path.join(HERE, s), run_name="__main__"); print("done", s)
