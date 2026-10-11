"""Assemble REPORT.md from code/report_template.md, results/tables.md sections and code/report_text.md.

    python3 cs581/basemeth/code/fill_report.py
"""
import os
import re

B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
tables = open(os.path.join(B, "results", "tables.md")).read()
text = open(os.path.join(B, "code", "report_text.md")).read()
tpl = open(os.path.join(B, "code", "report_template.md")).read()


def section(title_start):
    """Body of the tables.md section whose heading starts with title_start (heading line dropped)."""
    m = re.search(r"^## " + re.escape(title_start) + r".*?$(.*?)(?=^## |\Z)", tables, re.S | re.M)
    return m.group(1).strip() if m else "(not available)"


parts = dict(re.findall(r"^@@(\w+)\n(.*?)(?=^@@|\Z)", text, re.S | re.M))
fill = {
    "MAINTABLE": section("Final MAGUS"),
    "PAIRED": section("Paired differences"),
    "PRIMARY": section("Pre-registered primary"),
    "SUBSET": "Subset-level SP error (%), the 25 subset alignments of each method scored against the reference "
              "restricted to the subset, pooled:\n\n" + section("Subset-level"),
    "COST": "Total CPU seconds of the 25 subset alignments per dataset (mean over the set); ratio to L-INS-i.\n\n"
            + section("Subset-alignment cost") + "\n\nControl (merge of MAGUS's own subsets vs re-aligned L-INS-i "
            "subsets, SP error %):\n\n" + section("Control"),
    "TREES": section("Trees"),
    "BASE": parts.get("BASE", "").strip() + "\n\n" + section("Whole-dataset"),
    "METHODS": open(os.path.join(B, "REPORT_methods.md")).read().strip(),
}
for k in ("VERDICT", "DNA", "PLAN", "RISKS"):
    fill[k] = parts.get(k, "").strip()
for k, v in fill.items():
    tpl = tpl.replace("__{}__".format(k), v)
open(os.path.join(B, "REPORT.md"), "w").write(tpl)
