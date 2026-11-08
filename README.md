# protocol-ladder

Code and results for *How far validation protocols transfer between datasets: a
data-geometry audit of blast furnace slag viscosity models*.

Twenty-two model configurations are evaluated under four validation regimes on
three datasets. Everything the paper reports can be regenerated from the files
here, and `python verify.py` fails the build if any of 28 quoted numbers moves.

## The ladder

| Rung | Rule | Where it is defined |
|---|---|---|
| T1 | resubstitution, train = test | plant record only |
| T2 | row-wise random K-fold | all three datasets |
| T3 | source-grouped, every record of one sample held out | all three datasets |
| T4 | strictly forward in time | plant record only |

T1 and T4 need a record that carries a temperature curve and a time stamp, so
they are defined only on the plant record. The two literature-derived sets
contain at most one record per composition, so only T2 and T3 apply there, and
their T3 holds out a block of the composition domain (the compilation) or a
whole literature source (both). **The four rungs are therefore analogous across
datasets, not identical**, and the paper says so in Section 2.1.

## What reproduces what

| Paper item | Artefact in `outputs/` | Regenerate with |
|---|---|---|
| Plant record, 22 configurations × 4 tiers (88 rows) | `g4_plant_ladder.csv` | `scripts/run_ladder_plant.py`, needs the plant record, about 10 min |
| Wide-range database (88 rows) | `g4_gpvisc_ladder.csv` | `scripts/run_full_grid_gpvisc.py`, public data only, about 10 min |
| Literature compilation (88 rows) | `g4_compilation_ladder.csv` | shipped; assembly script not included, see Data below |
| One configuration with per-fold predictions | `metrics_rebuilt.csv` | `scripts/run_ladder_gpvisc.py`, about 2 min |
| Sampling-density manipulation (Table S4) | `spacing_experiment.csv`, `spacing_part1_draws.csv` | `scripts/spacing_experiment.py`, public data only |

## Key numbers pinned by `verify.py`

- Median paired loss T2 → T3: **0.7644** (plant), **0.0086** (wide-range), **0.7724** (compilation)
- Plant ladder, median R² by rung: +0.9773 / +0.7274 / +0.0534 / −1.5302
- Plant ladder, Extra-Trees on F3: RMSE 0.0000 / 0.2163 / 0.4804 / 0.5575
- Wide-range ladder, Extra-Trees on F3: RMSE 0.4312 / 0.5381 / 0.6262 / 0.9255

The 89-fold contrast the paper quotes is the ratio of the first two entries in
the first line, 0.7644 / 0.0086. The compilation is reported separately because
its grouped rung measures extrapolation, not leakage.

## Sampling-density manipulation (Table S4)

`outputs/spacing_experiment.csv` and `outputs/spacing_part1_draws.csv` back the
controlled check quoted in Section 1.3 and tabulated as Table S4 of the
Supplementary Information. Regenerate both with
`python scripts/spacing_experiment.py` (public data only, a few minutes).

| Condition | Placement of the five points | Spacing (K) | Repeats | Paired loss | sd |
|---|---|---|---|---|---|
| tight | consecutive | 32.0 | 30 | 0.2100 | 0.0849 |
| wide | quantile-spaced | 79.5 | 30 | 0.1185 | 0.0435 |

The difference is 0.0915 with SE 0.0174 and t = 5.25; a companion 100-draw
observational sweep gives Spearman rho = -0.304 (p = 0.0021). Draws are
independent between conditions, so the test is two-sample rather than paired:
the manipulation fixes the placement rule, not the identity of the compositions.
Extrapolating the observed slope to the 9 K plant spacing predicts 0.273, only
36% of the plant loss 0.7644, so Table S4 supports the direction of the effect
and not a quantitative transfer function.

## Quick start

```bash
pip install -r requirements.txt
python scripts/download_public.py --check   # verifies the public database
python verify.py                            # 28 checks against the shipped CSVs
```

## Data

Three datasets, not equally available.

- **Wide-range database** (10,739 records, 1,367 groups). Public. Fetched by
  `scripts/download_public.py` from the gpvisc release (DOI
  10.5281/zenodo.13843250, version 0.3.5) and checked against a recorded
  SHA-256. Nothing is repackaged.
- **Literature compilation** (242 records, 241 compositions). Redistributed as a
  machine-readable table in which every row carries its source label. The
  per-configuration results are shipped here; the assembly script is not,
  because the underlying workbooks belong to their original authors.
- **Plant record** (967 records, 110 slag samples, 1450–1530 °C). Proprietary,
  held under a data-use agreement, and **not** redistributed. What is shipped is
  its schema and the per-configuration results, so the plant owner can
  regenerate every table from their own copy without editing any code.

Set `DSTAR_DATA`, `DSTAR_OUT`, `DSTAR_PLANT_CSV` or `DSTAR_GPVISC_CSV` to
relocate any of these.

## Note on double-blind review

The manuscript is reviewed double-blind. This repository carries author names
and must not be cited in the blinded version. For the review copy, use an
anonymised mirror and cite that instead; the named address and a DOI belong in
the accepted version only.

## Layout

```
paths.py                central path resolution, overridable by environment variable
verify.py               21 checks against the shipped results
scripts/
  download_public.py    fetch + checksum the public database
  run_ladder_plant.py   22 configurations × 4 tiers on the plant record
  run_full_grid_gpvisc.py  the same grid on the wide-range database
  run_ladder_gpvisc.py  one configuration, with per-fold predictions
  res_common.py         shared loading and metrics
  make_figures_resolution.py
  spacing_experiment.py  sampling-density manipulation (Table S4)
outputs/                per-configuration results
                        spacing_experiment.csv, spacing_part1_draws.csv
data/schema/            expected schema of the plant record
```

## Licence

MIT, see `LICENSE`. The plant record is not covered by this licence and is not
included.

## Contact

Jinyu Wang, Hongwei Guo (corresponding), Helan Liang (corresponding),
Soochow University.
