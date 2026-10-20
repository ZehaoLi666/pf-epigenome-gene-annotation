# Supplementary tables

`build_S1.py` … `build_S7_S9.py` assemble Tables S1–S9 (one Excel workbook per table with a README sheet, plus one TSV
per sheet) from the intermediate result tables (`RESULTS_DIR`) and the figure source tables (`FIGURES_DIR`).
`common.py` holds the writer and the coordinate-convention text and is executed by each builder.
Output goes to `TABLES_OUT_DIR`. Requires pandas and openpyxl.

Supplementary figures are assembled by hand from the source tables listed in their legends; no plotting code
for them is included.
