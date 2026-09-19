# ahr v1

_Generated from `manifest.toml` and `metrics.json`. Do not edit._

Released 2026-09-06 · signature 1

**Why this version.** first release: binary XGBoost on Morgan, MACCS and RDKit descriptors over the Tox21 aryl-hydrocarbon-receptor screen

## Outputs

| column | dtype | range | meaning |
|---|---|---|---|
| `ahr_agonist` | float32 | 0.0–1.0 | P(activates the aryl hydrocarbon receptor signalling pathway in the Tox21 qHTS reporter). The readout is transcriptional output. |
| `ahr_cytotox` | float32 | 0.0–1.0 | P(reduces viability in the counter-screen over the same library). A compound scoring high on both readouts activated the reporter in a cell that was also dying |

Missing values: NaN when RDKit cannot parse the input SMILES

## Performance — `scaffold-shuffle-5seed@1`

Protocol `scaffold-shuffle-5seed@1` — Bemis-Murcko scaffold split with scaffold groups permuted by seed, so distinct seeds give distinct test sets. Five seeds; report mean and standard deviation over the held-out test folds.

Evaluated 2026-09-06 on n_train=5208, n_val=12, n_test=1724.

### `ahr_agonist`

| metric | mean | std | per seed |
|---|---|---|---|
| auc_roc | 0.691 | 0.113 | 0.577, 0.814, 0.583, 0.641, 0.839 |
| auprc | 0.192 | 0.219 | 0.006, 0.436, 0.010, 0.024, 0.482 |
| mcc | 0.153 | 0.177 | -0.010, 0.350, 0.040, 0.000, 0.386 |
| brier | 0.174 | 0.061 | 0.230, 0.100, 0.231, 0.211, 0.100 |

### `ahr_cytotox`

Measured on n_train=4753, n_val=12, n_test=1690.

| metric | mean | std | per seed |
|---|---|---|---|
| auc_roc | 0.794 | 0.067 | 0.722, 0.880, 0.747, 0.751, 0.872 |
| auprc | 0.217 | 0.123 | 0.090, 0.348, 0.203, 0.077, 0.366 |
| mcc | 0.246 | 0.075 | 0.190, 0.352, 0.225, 0.151, 0.312 |
| brier | 0.163 | 0.068 | 0.231, 0.066, 0.185, 0.230, 0.101 |

## Performance — `scaffold-balanced-5seed@1`

Protocol `scaffold-balanced-5seed@1` — Bemis-Murcko scaffold split with scaffold groups permuted by seed and each group placed in the fold it overfills least, so a group larger than a fold's capacity settles in train instead of starving that fold. Same fold fractions, seeds and metrics as scaffold-shuffle-5seed@1; only the packing differs. Five seeds; report mean and standard deviation over the held-out test folds.

Evaluated 2026-09-19 on n_train=5679, n_val=506, n_test=759.

### `ahr_agonist`

| metric | mean | std | per seed |
|---|---|---|---|
| auc_roc | 0.859 | 0.021 | 0.843, 0.902, 0.851, 0.850, 0.850 |
| auprc | 0.505 | 0.067 | 0.429, 0.611, 0.482, 0.551, 0.453 |
| mcc | 0.415 | 0.037 | 0.406, 0.481, 0.372, 0.428, 0.391 |
| brier | 0.123 | 0.039 | 0.154, 0.075, 0.090, 0.115, 0.181 |

### `ahr_cytotox`

Measured on n_train=5338, n_val=431, n_test=686.

| metric | mean | std | per seed |
|---|---|---|---|
| auc_roc | 0.831 | 0.021 | 0.832, 0.824, 0.847, 0.796, 0.858 |
| auprc | 0.358 | 0.036 | 0.410, 0.366, 0.304, 0.337, 0.375 |
| mcc | 0.319 | 0.059 | 0.373, 0.373, 0.240, 0.256, 0.355 |
| brier | 0.079 | 0.018 | 0.061, 0.064, 0.098, 0.103, 0.068 |

> Comparable only with metrics carrying the same protocol id.

## Data

Tox21 aryl hydrocarbon receptor qHTS — PubChem BioAssay AID 743085 (qHTS assay to identify small molecule that activate the aryl hydrocarbon receptor (AhR) signaling pathway) and AID 743086 (qHTS assay to identify small molecule that activate the aryl hydrocarbon receptor (AhR) signaling pathway - cell viability counter screen), rows called Active or Inactive, one row per compound labelled by majority call across its assay records. Retrieved 2026-09-06, licensed public-domain, redistributed here.

`6944` compounds, positive rate `0.084`, table SHA-256 `a51cfd84d52a7c2d…`

Regenerate and check for upstream drift with `python -m vp_ahr.data fetch --verify`.

## Model

xgboost-binary on `combo3` features. Shipped weights: one model per output, each on every compound its endpoint labels minus a 10% scaffold carve used for early stopping

`weights.joblib` SHA-256 `a0b0da88f1580ebe…`

## Provenance

Environment: python 3.11.11, rdkit 2026.03.5, xgboost 3.2.0.

Reproducibility is to this dataset hash and this environment, not bit-exact: the sources are live endpoints and RDKit descriptor values move between releases.
