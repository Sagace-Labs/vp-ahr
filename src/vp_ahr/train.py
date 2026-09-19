"""Build a new released version.

    python -m vp_ahr.train --version v2 --reason "why this version exists"

Writes ``versions/<version>/manifest.toml`` and ``weights.joblib``.

The shipped weights are fit on the whole dataset, with a scaffold carve held
out only for early stopping.
"""

from __future__ import annotations

import argparse
import platform
import shutil
import sys
from datetime import date
from pathlib import Path

from vp_ahr import contract as ahr_contract
from vp_ahr import data as ahr_data
from vp_ahr import model as ahr_model
from vp_ahr.target import CYTOTOX, TARGET
from vp_core import fingerprints

__all__ = ["build_version", "main"]

VERSIONS_DIR = Path(__file__).resolve().parent / "versions"

DEFAULT_PROTOCOL = "scaffold-shuffle-5seed@1"

#: Which label column supplies each declared output.
OUTPUT_LABELS: dict[str, str] = {
    "ahr_agonist": "label",
    "ahr_cytotox": "cytotox",
}


def _provenance() -> dict:
    """The environment that produced the version."""
    import rdkit
    import xgboost

    return {
        "python": platform.python_version(),
        "rdkit": rdkit.__version__,
        "xgboost": xgboost.__version__,
    }


def build_version(
    version: str,
    *,
    reason: str,
    protocol: str = DEFAULT_PROTOCOL,
    supersedes: str | None = None,
    seed: int = 0,
) -> Path:
    """Fit the deployment models and write the version directory."""
    from vp_core import protocols
    from vp_core.splits import scaffold_train_val

    protocols.get(protocol)  # fail early on an unknown protocol

    directory = VERSIONS_DIR / version
    if directory.exists():
        raise FileExistsError(
            f"{directory} already exists. Released versions are immutable — "
            "publish a new version instead of editing this one."
        )

    table = ahr_data.load()
    smiles = table["smiles"].tolist()
    X = fingerprints.featurize(smiles, ahr_model.FEATURES)

    fitted = {}
    for output in ahr_contract.column_names():
        column = OUTPUT_LABELS[output]
        rows = ahr_data.labelled(table, column)
        y = table[column].to_numpy()[rows].astype(int)
        subset = [smiles[i] for i in rows]

        # Deployment fit: everything this endpoint labels, with a small scaffold
        # carve that stops boosting.
        train_idx, val_idx = scaffold_train_val(subset, val_frac=0.10, seed=seed)
        fitted[output] = ahr_model.fit(
            X[rows][train_idx],
            y[train_idx],
            X[rows][val_idx],
            y[val_idx],
            seed=seed,
        )
        print(
            f"  {output}: {len(rows)} labelled compounds, "
            f"positive rate {y.mean():.3f}, "
            f"best iteration {getattr(fitted[output], 'best_iteration', None)}",
            file=sys.stderr,
        )

    directory.mkdir(parents=True)
    try:
        return _write_version(
            directory,
            fitted,
            table,
            version=version,
            reason=reason,
            protocol=protocol,
            supersedes=supersedes,
        )
    except Exception:
        shutil.rmtree(directory, ignore_errors=True)
        raise


def _write_version(
    directory: Path,
    fitted: dict,
    table,
    *,
    version: str,
    reason: str,
    protocol: str,
    supersedes: str | None,
) -> Path:
    import joblib

    import vp_core
    from vp_core import dataset, hashing, manifest

    weights_path = directory / "weights.joblib"
    joblib.dump(fitted, weights_path)

    labels = list(ahr_data.LABELS)
    record = {
        "schema": manifest.SCHEMA_VERSION,
        "pathway": "ahr",
        "version": version,
        "released": date.today().isoformat(),
        "supersedes": supersedes,
        "reason": reason,
        "signature": ahr_contract.as_manifest_table(),
        "dataset": {
            "name": "Tox21 aryl hydrocarbon receptor qHTS",
            "source": (
                f"PubChem BioAssay AID {TARGET.pubchem_aid} ({TARGET.assay_name}) "
                f"and AID {CYTOTOX.pubchem_aid} ({CYTOTOX.assay_name}), rows called "
                f"Active or Inactive, one row per compound labelled by majority "
                f"call across its assay records"
            ),
            "url": (
                "https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/"
                f"{TARGET.pubchem_aid}/concise/CSV"
            ),
            "retrieved": "2026-09-06",
            "licence": "public-domain",
            "redistributable": True,
            "path": "data/ahr_tox21.parquet",
            "labels": labels,
            "sha256": dataset.dataset_hash(table, labels=labels),
            "n_rows": len(table),
            "base_rate": round(float(table["label"].mean()), 6),
            "fetch": "python -m vp_ahr.data fetch --verify",
        },
        "model": {
            "family": "xgboost-binary",
            "features": ahr_model.FEATURES,
            "fit": (
                "one model per output, each on every compound its endpoint labels "
                "minus a 10% scaffold carve used for early stopping"
            ),
            "weights": "weights.joblib",
            "sha256": hashing.sha256_file(weights_path),
        },
        "protocol": {
            "id": protocol,
            "provider": "vp-core",
            "core_version": vp_core.__version__,
        },
        "provenance": _provenance(),
    }

    problems = manifest.validate(record, version_dir=directory)
    if problems:
        raise ValueError(f"refusing to write an invalid manifest: {'; '.join(problems)}")
    manifest.write(directory / "manifest.toml", record)

    print(
        f"wrote {directory}\n"
        f"  {len(table)} compounds, positive rate {table['label'].mean():.3f}\n"
        f"  next: python -m vp_ahr.evaluate --version {version}",
        file=sys.stderr,
    )
    return directory


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m vp_ahr.train")
    parser.add_argument("--version", required=True, help="new version name, e.g. v2")
    parser.add_argument("--reason", required=True, help="why this version exists")
    parser.add_argument("--protocol", default=DEFAULT_PROTOCOL)
    parser.add_argument("--supersedes", default=None)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args(argv)

    build_version(
        args.version,
        reason=args.reason,
        protocol=args.protocol,
        supersedes=args.supersedes,
        seed=args.seed,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
