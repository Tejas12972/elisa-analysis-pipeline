# ELISAtools example plates (Feng et al.): 5 real SoftMax Pro sandwich-ELISA plates

## Source

- **Package:** ELISAtools, "ELISA Data Analysis with Batch Correction", version 0.1.8
  (CRAN, published 2025-04-06). Files were taken from `inst/extdata/`.
- **Author / maintainer:** Feng Feng (Boston University, ffeng@BU.edu)
- **Retrieved from:** https://github.com/cran/ELISAtools, the read-only CRAN mirror, at
  commit `d28128c9c29c52a9e52c7896e500010536cdcd0b`. The Bioconductor git (`git.bioconductor.org/packages/ELISAtools`) and
  `github.com/bioc/ELISAtools` no longer serve this package, so the CRAN copy is the live one.
  Upstream development happens at https://github.com/BULQI/ELISAtools.
- **Associated paper:** Feng F, Thompson MP, Thomas BE, Duffy ER, Kim J, Kurosawa S, Tashjian JY,
  Wei Y, Andry C, Stearns-Kurosawa DJ. *A Computational Solution To Improve Biomarker
  Reproducibility During Long-term Projects.* bioRxiv 2018, doi:10.1101/483800
  (v1 CC-BY, v2 CC-BY-NC). The data files are covered by the package's MIT license, not
  by the preprint's license.
- **Date retrieved:** 2026-09-29

## License

The DESCRIPTION file says `License: MIT + file LICENSE`. The LICENSE file (copied here unmodified) reads:

```
YEAR: 2019
COPYRIGHT HOLDER: Feng Feng
```

That is the CRAN convention for the standard MIT license:

> Copyright (c) 2019 Feng Feng
>
> Permission is hereby granted, free of charge, to any person obtaining a copy of this software
> and associated documentation files (the "Software"), to deal in the Software without
> restriction, including without limitation the rights to use, copy, modify, merge, publish,
> distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the
> Software is furnished to do so, subject to the following conditions: The above copyright notice
> and this permission notice shall be included in all copies or substantial portions of the
> Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND ...

You may redistribute these files if you keep this notice. `LICENSE` and `DESCRIPTION` are kept
next to the data for that reason.

## Files (all unmodified)

| File | What it is |
|---|---|
| `Assay_2.txt` | SoftMax Pro text export. One plate ("assay 2"). Saved 10/4/2017. |
| `Assay_3_and_4.txt` | SoftMax Pro text export. Two plates ("assay 3 and 4"). Saved 10/4/2017. |
| `Assay_11_and_12.txt` | SoftMax Pro text export. Two plates ("assay 11 and 12 12202017"). Saved 12/21/2017. |
| `design.txt` | ELISAtools run/batch design table (TSV, CRLF line endings). |
| `AnExp_2plate.txt` | ELISAtools plate-map "annotation" grid. It does not match the real layouts; see "Oddities". |
| `stdConc.txt` | ELISAtools standard-concentration table: s1..s8 = 3000, 1500, 750, 375, 187.5, 93.8, 46.9, 0. |
| `DESCRIPTION`, `LICENSE` | Package metadata and license stub, kept for attribution. |

Not copied: the package's other files. `plate2_27_dataSdf.txt` has no matching layout, and
its values do not fit the annotation it is paired with in `design_plate2_27.txt`. Also left out:
`annote*.txt` (a different, unrelated layout), `tutorial.R`, and a `.docx` tutorial.

## Assay and analyte

- **Assay type:** sandwich ELISA with a TMB/HRP readout. It is read at **450 nm with a 620 nm
  reference**: the header line says `Absorbance ... 2 wavelengths "450 620"`.
- **Analyte: not documented anywhere in the package.** The standard units are **pg/mL** (from the
  SoftMax group header `Concentration pg/mL`). The associated paper measured plasma biomarkers
  from cancer patients over about 10 months. Its headline analyte is called only "PF", and there
  were four other biomarkers. The example files are plausibly from that project, but nothing
  says which analyte they are. **Treat the analyte as unknown.**
- **Samples:** unknowns `P1`..`P15`, each in triplicate, all at a **1:10 dilution**. The dilution
  factor sits in the SoftMax `Unk_Dilution` group, column `Dilution` = 10. The paper's context
  suggests human plasma, but the files do not say so. There is one kit/QC **Control** in
  triplicate on every plate.
- **No biological groups.** The files contain no treatment, case or control labels, only P-numbers.
- **Batches:** `design.txt` puts `Assay_2` plus `Assay_3_and_4` in **Batch1** (3 plates) and
  `Assay_11_and_12` in **Batch2** (2 plates). In the paper's sense, a batch is one ELISA
  kit lot. The date column (9/18/2009) is the same on every row and looks like a placeholder.
  The export filenames date the runs to 2017.

## File format (SoftMax Pro "Plate + Groups" text export)

- Tab-delimited with **CRLF** line endings, encoded **Latin-1**. `Temperature(¡C)` contains byte 0xA1,
  which is not valid UTF-8. Every line is padded with trailing tabs to 98 fields.
- Line 1 is `##BLOCKS= N`. N is the number of blocks, counting both plate blocks and group blocks.
- **Plate block**, one per plate:
  - `Plate:\tPlate1\t1.3\tTimeFormat\tEndpoint\tAbsorbance\tRaw\tFALSE\t1\t...\t2\t450 620 \t1\t12\t96\t1\t8`
  - Header line: `\tTemperature(¡C)\tA1\tA2 ... H12`, which is **one row in long/row-major order,
    96 wells**, not an 8x12 grid.
  - Data line 1: `\t<temp>\t<96 values>`. These are the **450 nm** ODs.
  - A blank line.
  - Data line 2: `\t<temp>\t<96 values>`. These are the **620 nm** reference ODs, all about 0.035.
  - A blank line, then `~End`.
  - The software's reported value for each well is `OD450 - OD620`. For example, Assay_2 A1:
    2.2416 - 0.0372 = 2.204, which is the value shown in the Standards group. ELISAtools'
    `read.plate()` does the same subtraction: it treats the second line as the "blank".
- **Group blocks** follow each plate. They hold SoftMax's own analysis and **carry the real plate
  map**:
  - `Group: Standards`. Columns: `Sample`, `Concentration pg/mL`, `BackCalcConc`, `Wells`,
    `Value`, `MeanValue`, `SD`, `CV`. Each standard's nominal concentration appears on the first
    row of its triplicate. `Masked` marks a well the operator excluded. `Range?` means the
    back-calculated value fell outside the curve.
  - `Group: Unk_Dilution`. Columns: `Sample`, `Wells`, `Value`, `R`, `Result`, `MeanResult`,
    `SD`, `CV`, `Dilution`, `AdjResult`. `R` flags wells outside the standard range.
  - `Group: Control`: the QC control triplicate.
  - `Group: Unknowns`, only in `Assay_11_and_12`: `Blank Wells Diluent Only`, which are the
    diluent-only blank wells.
  - Each group is followed by `Group Column` formula definitions and `Group Summaries`, then `~End`.
- The last line holds the original filename and save date.

## Plate layouts

These come from the SoftMax group blocks, which are authoritative. They are not taken from
`AnExp_2plate.txt`.

**Standards.** Columns 1-3, rows A-H, in triplicate: 3000, 1500, 750, 375, 187.5, 93.75, 46.875 and
0 pg/mL. The 0 pg/mL standard (row H) is the zero standard.

| Plate | Standards | Unknowns (P1-P15, triplicate, 1:10) | Control | Blanks |
|---|---|---|---|---|
| Assay_2 p1 | A-H rows, cols 1-3. Masked: C1, E1, F1, G2 | P1=B4-6, P2=E4-6, P3=C4-6, P4=F4-6, P5=D4-6, P6=G4-6, P7=D7-9, P8=B7-9, P9=E7-9, P10=C7-9, P11=F7-9, P12=C11+D10+D11, P13=G7-9, P14=E10+E11+F10, P15=B10+B11+C10 | F11, G10, G11 | Only the H1-3 zero standard. Row A cols 4-12, row H cols 4-12 and col 12 (OD about 0.06) are **unassigned**, presumably empty or diluent. |
| Assay_3_and_4 p1 | Same layout. Nothing masked. | Same layout as Assay_2 | F11, G10, G11 | Same as above |
| Assay_3_and_4 p2 | Same layout | Same layout | F11, G10, G11 | Same as above |
| Assay_11_and_12 p1 | **Row A = 1500 and row B = 3000 (swapped)**. Masked: A3, H1 | P1-P11 as above. P12=E10-G10, P13=G7-9, P14=B11-D11, P15=B10-D10 | E11-G11 | **Explicit**: `Blank Wells Diluent Only` = A4-A12, B12-G12, H4-H12 (24 wells) |
| Assay_11_and_12 p2 | Normal order (A=3000). Masked: H1 | Only **P1-P14**, and in a different arrangement: P5=G4-6, P6=D7-9, ... P14=B10-D10 | E11-G11 | A4-A12, B12, C12, D4-D6, D12-G12, H4-H12 (27 wells) |

Other points:

- Unknown ODs are low: about 0.04-0.2, against a top standard of about 2.2-2.5. Most samples
  read near the bottom of the curve. In Assay_11_and_12 plate 1, several are below the curve
  (`Range?`). The Control sits mid-curve at about 0.7-0.9 OD.
- The P-number labels are probably not the same patient from plate to plate. SoftMax reports very
  different values for the same label on different plates, for example P3 at 824 and 1165 pg/mL
  on the two Assay_3_and_4 plates. The labels look like layout slots in a reused template.
  **Nothing links sample identity across plates**, so treat the unknowns as anonymous.

## Oddities and gotchas

1. **`AnExp_2plate.txt` is wrong for these plates.** It maps p1 to B4-6, p2 to C4-6, P3 to D4-6,
   and so on. SoftMax says P1 = B4-6, P2 = E4-6, P3 = C4-6. The annotation also puts Control in
   E11-G11 on every plate, while Assay_2 and Assay_3_and_4 have it at F11, G10 and G11. And
   the annotation cannot show the swapped standards in Assay_11_and_12 plate 1. For a correct
   parse, build the layout from the SoftMax `Group:` blocks.
2. Inconsistent case in the annotation: `p1`, `p2` versus `P3`...`P15`.
3. Latin-1 encoding, CRLF line endings, trailing tab padding.
4. ODs are stored row-major on one line rather than as an 8x12 grid. Each plate has two
   wavelength lines.
5. Masked wells show as `Masked` in the group blocks. The raw OD lines still hold their values.
   In Assay_2, masked C1 has a raw OD450 of 0.694 against 0.940 and 0.779 for its partners.
6. Blank wells are labelled on only 2 of the 5 plates.

## What this dataset is good for

It is the best-documented real multi-plate data here for **standard-curve fitting, masked-well
and outlier QC, back-calculation with a 1:10 dilution, and inter-plate/inter-lot (batch)
comparison**. It has 5 plates in 2 kit-lot batches, with a repeated QC Control that works as
a Levey-Jennings / inter-assay CV anchor. It is **not** usable for group statistics or biological
interpretation, because it has no groups and the analyte is unknown.
