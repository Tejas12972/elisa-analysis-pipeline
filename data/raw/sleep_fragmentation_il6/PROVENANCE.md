# Mouse serum IL-6 sandwich ELISA after acute sleep fragmentation (4 plates, 11 groups)

## Source

- **Title (Zenodo/Dryad record):** "Gene expression and ELISA data". The README calls it
  *"Data for the article 'Inflammation from Sleep Fragmentation Starts in the Periphery Rather
  than Brain in Male Mice'"* and lists the journal as "Scientific Report", submitted.
  Final publication status and title were **not verified**.
- **Authors:** Van-Thuan Nguyen, Cameron Fields, Noah T. Ashley (Dept. of Biology, Western
  Kentucky University). Funding: NIH R15GM117534.
- **DOI:** 10.5061/dryad.tdz08kq3p (Dryad). Downloaded from the Zenodo mirror:
  https://zenodo.org/records/7627300, published 2023-02-09. Dryad's own download endpoints
  returned 401/403 to scripted requests.
- **Date retrieved:** 2026-09-29

## License

**CC0 1.0 Universal (public domain dedication).** Zenodo metadata: `license: cc-zero`. From the
README:

> License: Use of these data is covered by the following license: Title: CC0 1.0 Universal
> (CC0 1.0) Specification: https://creativecommons.org/publicdomain/zero/1.0/; the authors
> respectfully request to be contacted by researchers interested in the re-use of these data so
> that the possibility of collaboration can be discussed.

CC0 needs no attribution, but cite the DOI anyway.

## Files (unmodified)

| File | What it is |
|---|---|
| `Male_NSF_ASF_IL-6_ELISA_result.xlsx` | Three sheets: raw plate ODs with the authors' calculations, analyzed concentrations, and intra/inter-assay CVs. |
| `README.md` | The authors' README for the whole deposit. It also covers corticosterone and qPCR files, which were not downloaded. |

Not copied from the same deposit: `CORT-result.xlsx` (corticosterone, a competitive ELISA, holding
concentrations plus duplicate ODs), `Gene_Expression_results.xlsx` and
`Raw_data_Gene_Expression_ELISA_results.xlsx` (qPCR plus summary tables). They are available at
the DOI if a CORT-vs-IL-6 analysis is wanted later.

## Experiment

- **Animals:** male mice. **Groups:** 2 sleep treatments x time since treatment start:
  - `NSF` (non-sleep-fragmentation control) at 0, 1, 2, 6, 12, 24 h
  - `ASF` (acute sleep fragmentation) at 1, 2, 6, 12, 24 h. There is no ASF 0 h; NSF 0 h is
    the shared baseline.
  - About 10 mice per group, 110 sample IDs in total, with IDs like `M_ASF6h-3`. The README gives
    n = 9-10 per group after outlier removal.
- **Matrix:** trunk-blood serum.
- **Kit:** BioLegend **ELISA MAX Deluxe Set Mouse IL-6**, cat. 431304, a sandwich ELISA.
  The README gives intra- and inter-assay CV as 8.24% and 7.43%. The reader wavelength is not
  stated. The kit reads at 450 nm, and the values here are single-wavelength ODs.
- **Authors' statistics:** in GraphPad Prism, two-way ANOVA (treatment x time) and one-way
  ANOVA against 0 h, with Tukey and Bonferroni post-hoc tests. Outliers were removed at |z| > 2.
- **Expected biology,** from the authors' analyzed sheet: serum IL-6 rises sharply in ASF mice
  at **6 h**. The ASF 6 h group averages about 26 pg/mL (range 11-42) against about 1-6 pg/mL for NSF and other
  time points. Smaller rises appear at ASF 12 h and ASF 2 h.

## Sheet `IL-6 Raw and SD`: the raw data

This sheet is **four 96-well plates stacked vertically**, each written as a plate-reader grid. Excel
column C = plate column 1 ... column N = plate column 12, and Excel column B holds the row letter.

| Plate | Header row (1..12) | Grid rows (A..H) | Samples on plate |
|---|---|---|---|
| 1 | 1 | 2-9 | NSF0h-1..10, NSF1h-1..10, ASF1h-1..10 |
| 2 | 27 | 28-35 | NSF2h, ASF2h, NSF6h (x10 each) |
| 3 | 58 | 59-66 | ASF6h, NSF12h, ASF12h (x10 each) |
| 4 | 85 | 86-93 | NSF24h, ASF24h (x10 each) |

The label "Plate N" sits in column A partway down each grid, at rows 5, 31, 61 and 89.

**Standards.** Plate columns 1-2, in duplicate. Rows A-G = 500, 250, 125, 62.5, 31.3, 15.6 and
7.8 pg/mL. **Row H cols 1-2 is the 0 pg/mL standard.** The `IL-6 Intra and inter assay` sheet
lists H1/H2 as the eighth "Standard". Nominal concentrations are written next to each grid in
columns Q (label, e.g. `S500`) and R (value). Column S holds the mean OD by formula, and
column V holds mean OD minus blank.

**Blank the authors subtracted (`S0` cell).** It is **not** H1/H2:
- Plates 1-3: `=AVERAGE(K9:L9)`, meaning plate wells **H9, H10**. Those wells hold no sample.
- Plate 4: `=AVERAGE(K90:L90)`, meaning wells **E9, E10** (0.057 and 0.101). On plate 4 all of
  cols 9-10 and E-H 7-8 are unassigned, with ODs about 0.05. Plate 4's blank wells differ from the
  other plates, and the E9/E10 pair has poor precision.

**Samples.** Duplicate wells in adjacent column pairs (3-4, 5-6, 7-8, 9-10), running down the
rows. The formula cells in the tables under each grid map every sample to its wells. For
example, `M_NSF0h-1: =AVERAGE(E2:F2)` means wells A3/A4. The layout is column-major within
column pairs:

- **Plate 1:** NSF0h-1..8 = rows A-H of cols 3-4. NSF0h-9, -10 = A, B of cols 5-6. NSF1h-1..6 = C-H
  of cols 5-6. NSF1h-7..10 = A-D of cols 7-8. ASF1h-1..4 = E-H of cols 7-8. ASF1h-5..10 = A-F of
  cols 9-10. **Unused:** G9, G10, H9, H10 (the blank), and all of cols 11-12.
- **Plates 2 and 3:** the same pattern for three groups of 10: group 1 = A-H 3-4 plus A-B 5-6,
  group 2 = C-H 5-6 plus A-D 7-8, group 3 = E-H 7-8 plus A-F 9-10.
- **Plate 4:** two groups of 10. Group 1 = A-H 3-4 plus A-B 5-6. Group 2 = C-H 5-6 plus A-D 7-8.

Columns 11-12 on plates 1-3 read about 0.040. That is below the zero standard (about 0.08-0.19),
so these are likely empty or uncoated wells. Do not use them as blanks. Plate 3 A12 = 0.454 is a
stray value in an unused well. On plate 4, cols 11-12 read about 0.06-0.14 and are also unassigned.

**Dilution factors.** These are only implicit, in the concentration formulas. Most samples use
`=(OD - blank - a)/b`, which is neat. Samples flagged with `*` in the adjacent column use
`=2*(...)` or `=3*(...)`, which means a **2x or 3x dilution**. The README never states this, so
treat it as an inference from the formulas. On plate 1, NSF0h samples 2-8 and 10 apply the `3*` to
the **OD** (`D13 = 3*(C13-blank)`) rather than to the concentration. That is equivalent only
because the authors used a linear fit, and it is still an unusual thing to do.

| Plate | Samples at 2x (formula `2*`) | 3x |
|---|---|---|
| 1 | none | NSF0h-2, -3, -4, -5, -6, -7, -8, -10 (the `3*` is on OD) |
| 2 | ASF2h-1, -2, -7; NSF2h-4, -5, -9; NSF6h-1, -4, -7 | none |
| 3 | ASF6h-1, -3, -6; NSF12h-1, -3, -7; ASF12h-3, -6, -8, -9, -10 | none |
| 4 | NSF24h-4, -5, -7, -9; ASF24h-1, -3, -4 | ASF24h-2, -6 |

**The authors' curve fitting.** They used a **linear** fit per plate, of the form
`conc = (OD_net - a) / b`, with the coefficients typed in by hand:
- plate 1: a = 0.0032, b = 0.003. The NSF0h column uses a = **0.032** instead, which looks like a
  typo.
- plate 2: a = 0.0097, b = 0.0144.
- plate 3: a = 0.1167 (0.117 for the NSF12h column), b = 0.0104. Here the blank is **not**
  subtracted first.
- plate 4: a = 0.1038, b = 0.0103, also without blank subtraction.

Some sample concentrations are hard-coded to `0` rather than computed: plate 1 NSF0h-3, -6, -7
and -8. A 4PL/5PL refit from the raw ODs is a real methodological improvement here.

## Other sheets

- `IL-6 Analyzed and graph`: column A = sample ID and column B = concentration in pg/mL, which
  is **long format with group encoded in the ID**, plus pivoted NSF/ASF x time tables before and
  after outlier removal. This is the fallback concentration table for the statistics stage.
- `IL-6 Intra and inter assay`: long table with columns plate, sample ID, duplicate OD1, OD2, mean,
  SD and %CV. It includes the 8 standards per plate, H1/H2 being the zero. It is a convenient
  cross-check that the well map above is right: NSF0h-1 = 0.104 and 0.120 = A3/A4. It also
  contains a separate re-run inter-assay block (ASF6h-1..10, NSF12h-1..) in columns N-S.

## Oddities and gotchas

1. The file is an analysis workbook built from the plate export, not a reader export itself. Most
   cells are formulas. Use `openpyxl` with `data_only=False` to recover the well mapping. The raw
   OD cells are literal numbers.
2. The blank wells differ between plates 1-3 (H9/H10) and plate 4 (E9/E10).
3. Standard ODs vary a lot between plates. The 500 pg/mL top standard reads 1.55 on plate 1 and
   2.1-2.5 on plates 2-4. That is a good real-world inter-plate normalization case.
4. The design is confounded: each plate holds different groups, so plate is confounded with
   time point. Plate 1 has 0 h and 1 h, plate 2 has 2 h and NSF 6 h, plate 3 has ASF 6 h and
   12 h, plate 4 has 24 h. ASF 6 h and NSF 6 h are on **different plates**, so the headline
   comparison crosses plates.
5. The dilution factors come only from formula multipliers, and the 3x factor on plate 1 is
   applied to OD.
6. Many NSF samples fall below the 7.8 pg/mL standard (OD_net near 0). An LLOQ/censoring policy
   matters.
