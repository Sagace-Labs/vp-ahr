# vp-ahr

Predicts activation of the aryl hydrocarbon receptor from a SMILES string — a
molecular initiating event in which a planar aromatic ligand drives the
AhR/ARNT complex to the xenobiotic response element and induces CYP1A1 and
CYP1A2.

## Install

    pip install vp-ahr

## Use

    from vp_ahr import predict
    predict(["CC(=O)Oc1ccccc1C(=O)O"])   # -> DataFrame[ahr_agonist, ahr_cytotox]

Returns one row per input and one column per declared output. `ahr_agonist` is
the probability of activating the receptor's transcriptional output;
`ahr_cytotox` is the probability of reducing viability in the counter-screen
over the same library. A compound scoring high on both activated the reporter
in a cell that was also dying. The readout is transcription.
Unparseable SMILES come back as NaN. Pin a version with `predict(smiles,
version="v1")`; list what is available with `versions()`.

## Current version

**v1**, signature 1, measured under protocols `scaffold-shuffle-5seed@1`
and `scaffold-balanced-5seed@1`. The full record — metrics per output,
protocol and seed, dataset hash, environment — is in
[`src/vp_ahr/versions/v1/CARD.md`](src/vp_ahr/versions/v1/CARD.md).

## Data

PubChem BioAssay AID 743085, the Tox21 qHTS screen for activators of the aryl
hydrocarbon receptor signalling pathway, and AID 743086, its cell-viability
counter-screen over the same library. Both are reduced to one row per compound
labelled by the majority call across its assay records, retrieved 2026-09-06 and
redistributed here as a United States government work in the public domain.
Rebuild and check for upstream drift with
`python -m vp_ahr.data fetch --verify`; see [`data/README.md`](data/README.md)
for the expected layout.

## Retrain

    python -m vp_ahr.train --version v2 --reason "why this version exists"
    python -m vp_ahr.evaluate --version v2

`train` fits one deployment model per output on the whole dataset and writes a
new version directory; `evaluate` refits per seed under the protocol and
records what those held-out models scored. Reproducibility is to the recorded
dataset hash and environment, which can change.

## Licence

Code is Apache-2.0 ([`LICENSE`](LICENSE)). The bundled dataset is in the public
domain ([`LICENSE-DATA`](LICENSE-DATA)).

## Cite

See [`CITATION.cff`](CITATION.cff).
