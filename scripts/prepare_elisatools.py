"""Convert the ELISAtools example SoftMax Pro exports into pipeline inputs.

Source: CRAN package ELISAtools 0.1.8 (Feng Feng, MIT license), ``inst/extdata``,
kept unmodified in ``data/raw/elisatools_feng2019/``.

Five real sandwich-ELISA plates in two kit-lot batches. The plate map is
rebuilt from SoftMax's own ``Group:`` blocks, which are authoritative. The
package's ``AnExp_2plate.txt`` annotation does not match these plates, and
cannot represent the swapped 3000/1500 standard rows on Assay_11_and_12
plate 1.

* standards: ``Group: Standards`` (nominal level on the first row of each triplicate).
  The 0 pg/mL standard is typed ``blank``;
* unknowns: ``Group: Unk_Dilution``, dilution from its ``Dilution`` column (1:10).
  The P-labels are template slots, not patient IDs; the same label on two
  plates is not the same sample. They are therefore prefixed with the plate ID;
* ``Group: Control``: the kit QC control, the same material on every plate.
  It is the bridge for inter-plate normalization (``QC_CONTROL``);
* operator masks (``Masked`` in any group) are kept in an ``operator_masked``
  column. They are *not* applied, so the pipeline's own outlier QC can be
  compared with the analyst's judgement.

Usage::

    python scripts/prepare_elisatools.py
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from immunoassay.io import COLUMNS, ROWS, normalize_well, read_softmax_export

REPO = Path(__file__).resolve().parents[1]
RAW = REPO / "data/raw/elisatools_feng2019"
OUT = REPO / "data/processed/elisatools_feng2019"
BATCH = {"Assay_2": "Batch1", "Assay_3_and_4": "Batch1", "Assay_11_and_12": "Batch2"}


def _group_rows(df: pd.DataFrame) -> list[tuple[str, str, bool, pd.Series]]:
    """(sample, well, masked, row) for every well row of a SoftMax group table."""
    out = []
    well_col = "Wells"
    value_col = next(c for c in ("Value", "Values") if c in df.columns)
    for _, row in df.iterrows():
        well = normalize_well(str(row[well_col]))
        if well:
            out.append((str(row["Sample"]), well, row[value_col] == "Masked", row))
    return out


def main() -> None:
    (OUT / "plates").mkdir(parents=True, exist_ok=True)
    (OUT / "layouts").mkdir(parents=True, exist_ok=True)
    manifest = []
    for stem, batch in BATCH.items():
        for sp in read_softmax_export(RAW / f"{stem}.txt"):
            pid = sp.plate["plate_id"].iloc[0]
            layout = []
            std = sp.groups["Standards"]
            level = None
            for sample, well, masked, row in _group_rows(std):
                conc = row["Concentration pg/mL"]
                level = float(conc) if conc not in ("", None) else level
                typ = "blank" if level == 0 else "standard"
                layout.append(
                    {
                        "well": well,
                        "type": typ,
                        "sample_id": "" if typ == "blank" else f"STD{sample}",
                        "concentration": "" if typ == "blank" else level,
                        "dilution": 1,
                        "operator_masked": masked,
                    }
                )
            dil_by_sample = {}
            for sample, well, masked, row in _group_rows(sp.groups["Unk_Dilution"]):
                if row["Dilution"] not in ("", None):
                    dil_by_sample[sample] = float(row["Dilution"])
                layout.append(
                    {
                        "well": well,
                        "type": "sample",
                        "sample_id": f"{pid}:{sample}",
                        "concentration": "",
                        "dilution": dil_by_sample[sample],
                        "operator_masked": masked,
                    }
                )
            for _, well, masked, _ in _group_rows(sp.groups["Control"]):
                layout.append(
                    {
                        "well": well,
                        "type": "control",
                        "sample_id": "QC_CONTROL",
                        "concentration": "",
                        "dilution": 1,
                        "operator_masked": masked,
                    }
                )
            lay = pd.DataFrame(layout)
            if lay["well"].duplicated().any():
                raise RuntimeError(f"{pid}: duplicated wells in layout")

            grid = sp.plate.pivot(index="row", columns="column", values="od")
            lines = [
                f"Source,ELISAtools 0.1.8 extdata {stem}.txt plate {sp.index + 1} (MIT)",
                "Value,OD450 - OD620",
                "",
                "," + ",".join(map(str, COLUMNS)),
            ]
            lines += [f"{r}," + ",".join(f"{grid.loc[r, c]:.4f}" for c in COLUMNS) for r in ROWS]
            (OUT / "plates" / f"{pid}.csv").write_text("\n".join(lines) + "\n")
            lay.to_csv(OUT / "layouts" / f"{pid}_layout.csv", index=False)
            manifest.append(
                {
                    "plate_id": pid,
                    "analyte": "analyte (undisclosed)",
                    "plate_file": f"plates/{pid}.csv",
                    "layout_file": f"layouts/{pid}_layout.csv",
                    "units": "pg/mL",
                    "batch": batch,
                }
            )
            print(f"{pid}: {len(lay)} wells, {int(lay.operator_masked.sum())} masked by operator")
    pd.DataFrame(manifest).to_csv(OUT / "manifest.csv", index=False)


if __name__ == "__main__":
    main()
