#!/usr/bin/env python3
"""
Split Bloomington stock table into allele / insertion categories.

• Reads:  Bloomington_single_insertion_StockHasPaperRef_AlleleOrInsertion_NoBalancers_stocks.csv
• Writes: Allele/ Insertion Categories/<Classification>.csv

Author: Aadish (helped by ChatGPT), 2025-04-25
"""

# ───────── imports ─────────
from pathlib import Path
import re
import pandas as pd

# ───────── config ─────────
INPUT_CSV  = Path("Bloomington_single_insertion_StockHasPaperRef_AlleleOrInsertion_NoBalancers_stocks.csv")
OUT_DIR    = Path("Allele/ Insertion Categories")   # will be created if absent
# which genotype column to look at
GENO_COLS  = ["Genotype"]

# ordered (pattern, class-label) list  ➜ first match wins
RULES: list[tuple[re.Pattern, str]] = [
    (re.compile(r"Mi\{PT-GFSTF", re.I),                               "MiMIC protein trap"),
    (re.compile(r"(Mi\{Trojan-[^}]+\}|-TG4\.\d)", re.I),              "Trojan MiMIC line"),
    (re.compile(r"Mi\{(?:MIC)?", re.I),                               "MiMIC insertion"),

    (re.compile(r"PBac\{WHr?\}", re.I),                               "PBac-WH insertion"),
    (re.compile(r"PBac\{PB\}",  re.I),                                "PBac-PB insertion"),
    (re.compile(r"PBac\{RB\}",  re.I),                                "PBac-RB insertion"),

    (re.compile(r"PBac\{IT\.GAL4\}", re.I),                           "GAL4 driver (PBac-IT.GAL4)"),

    (re.compile(r"TI\{(?:2A-GAL4|TI|GAL4)\}", re.I),                  "Targeted GAL4 insertion"),
    (re.compile(r"TI\{(?:CRIMIC\.TG4\.1|CRIMIC\.TG4\.0|T-GEM\.1)\}", re.I),
                                                                     "CRIMIC-GAL4 trap"),
    (re.compile(r"TI\{(?:2A-lexA|lexA::UnkAD)\}", re.I),              "Targeted LexA driver"),

    (re.compile(r"PBac\{.*GFP.*\}", re.I),                            "Protein trap (GFP-tagged)"),
    (re.compile(r"M\{.*GFP.*\}",   re.I),                             "BAC-GFP transgene"),
    (re.compile(r"M\{.*mCherry.*\}", re.I),                           "BAC-mCherry reporter"),

    (re.compile(r"PBac\{3HPy\[\+\]\}", re.I),                         "PBac-Other insertion"),

    (re.compile(r"P\{.*UAS.*", re.I),                                 "UAS construct"),     # needs both checks below
    (re.compile(r"P\{",        re.I),                                 "Other P-element construct"),

    (re.compile(r"^[^{]*$",    re.I),                                 "Classical allele"),  # no “{” at all
]

# ───────── helpers ─────────
def pick_genotype(row) -> str:
    """Return the first non-NA genotype string."""
    for col in GENO_COLS:
        if col in row and pd.notna(row[col]):
            return str(row[col])
    return ""

def classify(genotype: str) -> str | None:
    """Return the first matching class label, or None."""
    # special case: distinguish UAS vs other P{…}
    if "P{" in genotype:
        if "UAS" in genotype:
            return "UAS construct"
        # fall through to other P-element (handled later)
    for pat, label in RULES:
        if pat.search(genotype):
            return label
    return None   # should be rare

# ───────── main ─────────
def main() -> None:
    df = pd.read_csv(INPUT_CSV)
    df = df.loc[:, ~df.columns.str.contains(r"^Unnamed")]
    df["__genotype"] = df.apply(pick_genotype, axis=1)
    df["Category"]   = df["__genotype"].map(classify)

    # ensure output directory exists
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # write one file per category
    for cat, sub in df.groupby("Category", sort=False):
        if cat is None:
            cat = "Unclassified"
        outfile = OUT_DIR / f"{cat}.csv"
        sub.drop(columns="__genotype").to_csv(outfile, index=False)
        print(f"wrote {len(sub):>4} rows ➜ {outfile}")

    if "FBst" in df.columns:
        # count distinct StockIDs in each category
        summary_df = (
            df.groupby("Category")["FBst"]
            .nunique()
            .reset_index(name="Count")
        )
    else:
        # fallback to row‐counts
        summary_df = (
            df.groupby("Category")
            .size()
            .reset_index(name="Count")
        )

        # ───── write summary to Excel ─────
    excel_path = "Category_counts.xlsx"
    with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
        summary_df.to_excel(writer, index=False, sheet_name="Counts")
    print(f"\nWrote summary → {excel_path}")

if __name__ == "__main__":
    main()
