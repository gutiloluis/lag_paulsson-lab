# CRISPRi library design (lLAG08 / lLAG10)

Code that designed the two mismatch-CRISPRi sgRNA libraries for *Bacillus subtilis*
PY79 used in this project (ordered May 2024).

## Contents

| Path | Origin |
|---|---|
| `generate_sgrnas.py` | **Third party** — van Gestel et al. (see Attribution) |
| `2024-03-15_Bacillus_CRISPRi_Library/generate_spacers.py` | Original |
| `2024-03-15_Bacillus_CRISPRi_Library/2024-04-19_generate_sgrnas_essentials_nonessentials.ipynb` | Original — the notebook that produced the ordered libraries |
| `2024-08-13_Library-Sampling-Simulation/Library-Sampling-Simulation.ipynb` | Original — sampling/representation simulation |
| `environment.yml` | conda env `library_design` (`bioframe` additionally installed with pip) |

Input data (genome, published supplementary tables, model weights) and all outputs
live outside this repository, under the working directory referenced as `headpath`
in the notebook.

## Attribution

### `generate_sgrnas.py` — third-party code

| | |
|---|---|
| Source | <https://github.com/jordivangestel/sgRNAs-for-mismatch-CRISPRi> |
| Authors | Jordi van Gestel, John S. Hawkins, Horia Todor, Carol A. Gross |
| Retrieved | 2024-03-18 (added here in commit `9d1b783`) |
| Licence | None stated upstream (no LICENSE file) |

Published as part of:

> van Gestel J, Hawkins JS, Todor H, Gross CA. Computational pipeline for designing
> guide RNAs for mismatch-CRISPRi. *STAR Protocols* 2(2):100521 (2021).
> doi:10.1016/j.xpro.2021.100521

**Local modification.** One change from upstream, applied at two places in
`construct_perfect_sgrnas()`: the feature selector `feature.type == "CDS"` was
replaced by `'locus_tag' in feature.qualifiers`, so that tRNA and rRNA features are
scanned as well as coding sequences. The original lines are retained as comments.
Apart from that change and the provenance header, the file is byte-for-byte upstream.

This copy is retained for reproducibility of
the library design and is not redistributed under any claim of ownership.

### Acknowledgements

We thank Jason M. Peters (University of Wisconsin–Madison) for advice on the design
of this library — in particular the recommendation to select mismatched sgRNAs
empirically from the measured data of Hawkins et al. (2020) — and Ryan D. Ward for
sharing mismatch-CRISPRi guide-design code.

### How to cite

- van Gestel J, Hawkins JS, Todor H, Gross CA. *STAR Protoc* 2(2):100521 (2021).
  doi:10.1016/j.xpro.2021.100521 — design pipeline
- Hawkins JS, Silvis MR, Koo BM, Peters JM, Osadnik H, Jost M, Hearne CC,
  Weissman JS, Todor H, Gross CA. *Cell Syst* 11(5):523–535.e9 (2020).
  doi:10.1016/j.cels.2020.09.009 — measured mismatch sgRNA fitness, controls,
  and the mismatch-efficacy model
- Koo BM, Kritikos G, Farelli JD, Todor H, Tong K, Kimsey H, Wapinski I,
  Galardini M, Cabal A, Peters JM, Hachmann AB, Rudner DZ, Allen KN, Typas A,
  Gross CA. *Cell Syst* 4(3):291–305.e7 (2017). doi:10.1016/j.cels.2016.12.013 —
  auxotrophic gene definitions

## Input data

| File | Source |
|---|---|
| `PY79_BSGAtlas_RATT_manually_curated_modified.gb` | PY79 genome, BSGAtlas annotation transferred with RATT, manually curated; commas in `locus_tag`/`label` replaced by semicolons (`change_commas_to_semicolon_in_locus_tag`) |
| `Table_S3_1-s2.0-S2405471220303665-mmc4.xlsx` | Hawkins et al. 2020, Table S3 — sheets `relative fitness (bsu data)` and `relative fitness (bsu controls)` |
| `KO_S4_mmc5.xlsx` | Koo et al. 2017, Table S4, sheet D — auxotrophs |
| `model_param_bsub.csv` | Hawkins et al. 2020, Table S2 ("Linear model parameters and weights from trained models", related to Figure 2), sheet `trained model params (bsu)`, converted to the CSV layout `generate_sgrnas.py` expects (`y-intercept` → `intercept`, `mismatch at nucleotide N` → `N`, `mismatch X -> Y` → `XY`, `parent GC` → `GC_content`) |

**Which trained model.** Table S2 contains three fits of the same 33-parameter
linear model, trained on *B. subtilis* GFP data, on *E. coli* GFP data, or on both.
This library used the *B. subtilis* fit, matching the organism. The `model_param.csv`
distributed with `generate_sgrnas.py` is the species-averaged ("both") fit, which is
also the fit behind the `relative fitness (predicted)` column of Table S3; predicted
efficacies here therefore sit about 0.04 above that column on average.

## What was designed how

Guides were obtained by two different routes:

- **From measured data** — essential genes and controls were taken directly from
  Hawkins et al. 2020 Table S3 (measured relative fitness), not designed de novo.
- **From `generate_sgrnas.py`** — auxotrophs, non-essential genes, tRNAs and rRNAs
  were found in the genome by the van Gestel script (`--step find,filter`, with
  `--off_target_upper 10` and the offset filter effectively disabled at
  `--offset_upper 1000000`). Mismatched variants and their predicted efficacies were
  generated with `--step mismatch` for auxotrophs, tRNAs and rRNAs; non-essential
  genes received perfect-match guides only.

Ordered spacers:

| Library | Set | Spacers | Genes | Route |
|---|---|---|---|---|
| **opLAG_1** (essential) | essential | 5,983 | 300 | Hawkins et al. 2020 |
| | auxotroph | 1,700 | 85 | `generate_sgrnas.py` |
| | tRNA | 200 | 20 | `generate_sgrnas.py` |
| | rRNA | 100 | 5 | `generate_sgrnas.py` |
| | control | 500 | — | Hawkins et al. 2020 |
| | **total** | **8,483** | | |
| **opLAG_2** (non-essential) | non-essential | 10,871 | 3,680 | `generate_sgrnas.py` |
| | control | 500 | — | Hawkins et al. 2020 |
| | **total** | **11,371** | | |

The same 500 control spacers were included in both order files.

## Selection, as implemented in `generate_spacers.py`

Everything below `generate_sgrnas.py` is original to this repository:

1. **Off-target re-analysis.** Each spacer's 9-nt PAM-proximal seed is located
   genome-wide on both strands (`find_all_hits_seeds_parallel`), hits are mapped to
   genes with `bioframe.overlap`, and each hit is characterised by the number of
   mismatches outside the seed and by whether it binds the non-template strand.
   Spacers with more than 10 off-targets are discarded.
2. **Biofilm exclusion.** Spacers with any off-target in `sinR`, `abrB` or `ywcC`
   are removed, to avoid biofilm-inducing knockdowns.
3. **Cloning filters.** Spacers containing BbsI (`GAAGAC`), BsaI (`GGTCTC`) or PvuI
   (`CGATCG`) sites, on either strand, are removed; the adapter-flanked order
   sequence is re-checked for BsaI and PvuI. Duplicated spacers are dropped entirely.
4. **Efficacy binning.** Predicted (and, for essential genes, measured) efficacy is
   binned into ten bins of width 0.1 (`BIN_DICT`). Guides are selected per gene and
   per bin, ranked by off-target count, off-target strand, number of mismatches, and
   then sequence-quality tie-breakers (`gc_downstream_8_9`, `agg_pam`,
   `n_t_in_start`, `offset_gt_1000`). Empty bins are filled from the nearest
   populated bin (`generate_bin_sample_single_gene`). Guides per gene: 20 for
   essential genes (two passes, measured and predicted) and auxotrophs, 20 for
   rRNAs, 10 for tRNAs, 3 for non-essential genes.
5. **Order sequences.** Selected spacers are flanked with cloning adapters
   (`add_adapter_seqs`) and a pool-specific suffix, and written to
   `df_ess_to_order_all.csv` / `df_non_to_order_all.csv` and the simplified
   two-column order files.
