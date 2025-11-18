#!/usr/bin/env python3
"""
Make a reference-centric Excel from a CSV with a 'Supporting_PMCIDs' column.

- Sheet 1 ("source"): exact contents of the input CSV
- Sheet 2 ("refs"): each PMCID exploded to its own row, retaining all
  original row metadata, plus author, journal, year, and title fetched
  from Europe PMC.

Usage:
  python make_refs_excel.py /path/to/input.csv

Requires:
  pip install pandas requests openpyxl
"""
import argparse
import pathlib
import time
import re
from typing import Dict, Optional, Tuple

import pandas as pd
import requests


EUROPE_PMC_SEARCH = (
    "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
)  # ?query=EXT_ID:PMCxxxx&format=json&pageSize=1

import os
import time
import requests
from typing import Optional, Tuple

NCBI_IDCONV = "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/"
NCBI_ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
NCBI_API_KEY = os.getenv("NCBI_API_KEY")  # optional but recommended for higher rate limits

def fetch_pmc_metadata(pmcid: str, retries: int = 3, backoff: float = 0.8
                      ) -> Optional[Tuple[str, str, str, str, str]]:
    """
    NCBI-based resolver:
      PMCID --(idconv)--> PMID --(esummary)--> (PMID, first_author_lastname, journal, year, title)
    """
    pmcid = (pmcid or "").strip()
    if not pmcid:
        return None

    # --- 1) PMCID -> PMID via ID Converter ---
    idconv_params = {"ids": pmcid, "format": "json"}
    if NCBI_API_KEY:
        idconv_params["api_key"] = NCBI_API_KEY

    pmid = None
    for attempt in range(retries):
        try:
            r = requests.get(NCBI_IDCONV, params=idconv_params, timeout=15)
            r.raise_for_status()
            data = r.json()
            recs = data.get("records") or []
            if not recs:
                return None
            pmid = recs[0].get("pmid")
            break
        except requests.RequestException:
            if attempt == retries - 1:
                return None
            time.sleep(backoff * (2 ** attempt))

    if not pmid:
        return None

    # --- 2) PMID -> metadata via ESummary ---
    esum_params = {"db": "pubmed", "id": pmid, "retmode": "json"}
    if NCBI_API_KEY:
        esum_params["api_key"] = NCBI_API_KEY

    for attempt in range(retries):
        try:
            r = requests.get(NCBI_ESUMMARY, params=esum_params, timeout=15)
            r.raise_for_status()
            js = r.json()
            uid = (js.get("result") or {}).get("uids", [None])[0]
            if not uid:
                return None
            rec = js["result"][uid]

            title = rec.get("title")
            journal = rec.get("fulljournalname") or rec.get("source")
            pubdate = rec.get("pubdate") or ""
            year = None
            for tok in str(pubdate).split():
                if tok[:4].isdigit():
                    year = tok[:4]
                    break

            authors = rec.get("authors") or []
            first_author_last = None
            if authors and isinstance(authors, list):
                nm = authors[0].get("name")  # "Smith J" or "Smith, John"
                if nm:
                    # Always force last name only
                    if "," in nm:
                        first_author_last = nm.split(",")[0].strip()
                    else:
                        first_author_last = nm.split()[0].strip()

            return pmid, first_author_last, journal, year, title
        except requests.RequestException:
            if attempt == retries - 1:
                return None
            time.sleep(backoff * (2 ** attempt))

def build_refs_sheet(df: pd.DataFrame, pmcid_col: str = "Supporting_PMCIDs") -> pd.DataFrame:
    if pmcid_col not in df.columns:
        raise ValueError(
            f"Expected column '{pmcid_col}' not found. Columns present: {list(df.columns)}"
        )

    # Split by ';' and explode to one PMCID per row
    work = df.copy()
    work[pmcid_col] = work[pmcid_col].fillna("").astype(str)
    def _normalize_pmcid(value: str) -> str:
        s = str(value).strip()
        if not s:
            return ""
        # Remove common prefixes and whitespace/punctuation
        s = s.replace("PMCID:", "").replace("pmcid:", "")
        s = s.replace("PMC ", "PMC").replace(" ", "")
        s = re.sub(r"[^A-Za-z0-9]", "", s)
        # Ensure PMC prefix if only digits are provided
        if s and not s.upper().startswith("PMC"):
            if s.isdigit():
                s = "PMC" + s
        # Normalize case to 'PMC' prefix
        if s.upper().startswith("PMC"):
            s = "PMC" + re.sub(r"^PMC", "", s.upper())
        return s

    work["PMCID"] = (
        work[pmcid_col]
        .str.split(r"\s*;\s*")
        .apply(lambda xs: [_normalize_pmcid(x) for x in xs if _normalize_pmcid(x)])
    )
    work = work.explode("PMCID", ignore_index=True)
    # Drop rows that had no PMCID
    work = work[work["PMCID"].notna() & (work["PMCID"] != "")]
    if work.empty:
        # Still provide structure if nothing to explode
        work = df.copy()
        work["PMCID"] = pd.NA

    # Fetch metadata for unique PMCIDs (minimize network calls)
    unique_pmcs = sorted(set(x for x in work["PMCID"].dropna()))
    meta_map: Dict[str, Tuple[str, str, str, str]] = {}
    for pmc in unique_pmcs:
        meta = fetch_pmc_metadata(pmc)
        if meta:
            meta_map[pmc] = meta

    # Map into columns
    work["Ref_PMID"] = work["PMCID"].map(lambda x: (meta_map.get(x) or [None]*5)[0])
    work["Ref_FirstAuthorLast"] = work["PMCID"].map(lambda x: (meta_map.get(x) or [None]*5)[1])
    work["Ref_Journal"] = work["PMCID"].map(lambda x: (meta_map.get(x) or [None]*5)[2])
    work["Ref_Year"] = work["PMCID"].map(lambda x: (meta_map.get(x) or [None]*5)[3])
    work["Ref_Title"] = work["PMCID"].map(lambda x: (meta_map.get(x) or [None]*5)[4])

    # New column: endnote_ref
    work["endnote_ref"] = work.apply(
        lambda r: f"{{{r['Ref_PMID']}, {r['Ref_FirstAuthorLast']}, {r['Ref_Year']}}}"
        if pd.notna(r["Ref_PMID"]) and pd.notna(r["Ref_FirstAuthorLast"]) and pd.notna(r["Ref_Year"])
        else None,
        axis=1
    )


    # Column order: PMCID + ref fields first, and drop the original pmcid_col on the refs sheet
    ref_cols = ["PMCID", "Ref_PMID", "Ref_FirstAuthorLast", "Ref_Journal", "Ref_Year", "endnote_ref", "Ref_Title"]
    base_cols = [c for c in df.columns if c not in ref_cols and c != pmcid_col]
    final_cols = [c for c in ref_cols if c in work.columns] + [
        c for c in base_cols if c in work.columns
    ]
    # Include any new columns created during processing (excluding pmcid_col)
    leftover = [c for c in work.columns if c not in final_cols and c != pmcid_col]
    final_cols = final_cols + leftover

    return work[final_cols]


def main():
    ap = argparse.ArgumentParser(
        description="Create a reference-centric Excel from a CSV with Supporting_PMCIDs."
    )
    ap.add_argument("csv_path", type=pathlib.Path, help="Path to input CSV")
    ap.add_argument(
        "--pmcid-col",
        default="Supporting_PMCIDs",
        help="Name of the column containing semicolon-separated PMCIDs (default: Supporting_PMCIDs)",
    )
    ap.add_argument(
        "--outfile",
        type=pathlib.Path,
        default=None,
        help="Optional output .xlsx path. Defaults to '<input>_with_refs.xlsx' next to the CSV.",
    )
    args = ap.parse_args()

    if not args.csv_path.exists():
        raise SystemExit(f"Input not found: {args.csv_path}")

    # Read CSV (let pandas infer; if you need a specific encoding, add encoding='utf-8')
    df = pd.read_csv(args.csv_path)

    # Build sheets
    refs_df = build_refs_sheet(df, pmcid_col=args.pmcid_col)

    # Output path
    if args.outfile is not None:
        out_xlsx = args.outfile
    else:
        out_xlsx = args.csv_path.with_name(args.csv_path.stem + "_with_refs.xlsx")

    # Write Excel
    with pd.ExcelWriter(out_xlsx, engine="openpyxl") as writer:
        df.to_excel(writer, sheet_name="source", index=False)
        refs_df.to_excel(writer, sheet_name="refs", index=False)

    print(f"✓ Wrote {out_xlsx}")


if __name__ == "__main__":
    main()

