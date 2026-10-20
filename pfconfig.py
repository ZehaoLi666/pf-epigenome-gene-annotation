"""Central path configuration for all scripts in this repository.

Every script that needs an input outside the repository imports its paths from here, e.g.

    from pfconfig import LAB_DIR, result
    gff = f"{LAB_DIR}/PlasmoDB/PlasmoDB-68/PlasmoDB-68_Pfalciparum3D7.gff"
    genes = result("gene_boundaries_v3_all_genes.tsv")

Values are read from a small file of `KEY: value` lines (YAML subset; PyYAML not needed):
  1. the file named by the environment variable PF_CONFIG, else
  2. <repository>/config/paths.yaml  (copy config/paths.example.yaml and edit it).
Any key can also be overridden by an environment variable of the same name prefixed with PF_
(e.g. PF_RESULTS_DIR=/data/pf_results).

Keys
  LAB_DIR         local data root: PlasmoDB releases, bigWig tracks ("lncRNA project tracks/"),
                  IGV_longread_transcript_models/, ...  (layout: see README, "Input data layout")
  RESULTS_DIR     folder holding the intermediate result tables deposited with the paper
                  (result(name) -> RESULTS_DIR/name)
  HPCC_BIGDATA    cluster work area used by the HPCC job scripts (large intermediate files)
  HPCC_LAB        cluster folder with raw data (SRA downloads, PlasmoDB genome files)
  HPCC_HOME       cluster home directory (conda environments)
  SUBMISSION_DIR  manuscript working folder (only referenced in a few descriptive strings)
  FIGURES_DIR     folder containing Figure1..Figure5/data (default: <repo>/figures)
  TABLES_OUT_DIR  output folder of the supplementary-table builders (default: <repo>/supplementary/tables/out)
  WORK_DIR        scratch folder for analysis outputs (default: current directory)
"""
import os

REPO = os.path.dirname(os.path.abspath(__file__))
_DEFAULTS = {
    "LAB_DIR": os.path.join(REPO, "data", "lab"),
    "RESULTS_DIR": os.path.join(REPO, "data", "results"),
    "HPCC_BIGDATA": os.path.join(REPO, "data", "hpcc_bigdata"),
    "HPCC_LAB": os.path.join(REPO, "data", "hpcc_lab"),
    "HPCC_HOME": os.path.join(REPO, "data", "hpcc_home"),
    "SUBMISSION_DIR": os.path.join(REPO, "data", "submission"),
    "FIGURES_DIR": os.path.join(REPO, "figures"),
    "TABLES_OUT_DIR": os.path.join(REPO, "supplementary", "tables", "out"),
    "WORK_DIR": os.getcwd(),
}


def _read(path):
    out = {}
    if not os.path.exists(path):
        return out
    for line in open(path, encoding="utf-8"):
        line = line.split("#", 1)[0].strip()
        if not line or ":" not in line:
            continue
        k, v = line.split(":", 1)
        v = v.strip().strip('"').strip("'")
        if v:
            out[k.strip()] = os.path.expanduser(os.path.expandvars(v))
    return out


_CFG = dict(_DEFAULTS)
_CFG.update(_read(os.environ.get("PF_CONFIG", os.path.join(REPO, "config", "paths.yaml"))))
for _k in list(_CFG):
    if os.environ.get("PF_" + _k):
        _CFG[_k] = os.environ["PF_" + _k]

LAB_DIR = _CFG["LAB_DIR"]
RESULTS_DIR = _CFG["RESULTS_DIR"]
HPCC_BIGDATA = _CFG["HPCC_BIGDATA"]
HPCC_LAB = _CFG["HPCC_LAB"]
HPCC_HOME = _CFG["HPCC_HOME"]
SUBMISSION_DIR = _CFG["SUBMISSION_DIR"]
FIGURES_DIR = _CFG["FIGURES_DIR"]
TABLES_OUT_DIR = _CFG["TABLES_OUT_DIR"]
WORK_DIR = _CFG["WORK_DIR"]


def result(name):
    """Path of a deposited intermediate result table (RESULTS_DIR/name)."""
    return os.path.join(RESULTS_DIR, name)


if __name__ == "__main__":
    for k in _DEFAULTS:
        print(f"{k:15s} {_CFG[k]}")
