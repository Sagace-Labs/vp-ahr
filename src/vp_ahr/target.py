"""The endpoints this pathway predicts.

Frozen so the data sources are auditable.

Both assays were verified live against PubChem on 2026-09-06, and screen the
same library:

    AID 743085   qHTS assay to identify small molecule that activate the aryl
                 hydrocarbon receptor (AhR) signaling pathway
                                                          10486 substances
    AID 743086   the same, cell viability counter screen   10486 substances
"""

from __future__ import annotations

from dataclasses import dataclass

__all__ = ["CYTOTOX", "TARGET", "TARGETS", "Endpoint", "all_names", "get"]


@dataclass(frozen=True)
class Endpoint:
    """One molecular initiating event, identified by the assay that reads it out."""

    name: str
    pathway: str
    mie: str
    pubchem_aid: int
    assay_name: str


TARGET = Endpoint(
    name="AHR",
    pathway="aryl hydrocarbon receptor / xenobiotic induction",
    mie=(
        "binding of the aryl hydrocarbon receptor by a planar aromatic ligand, which "
        "drives the AhR/ARNT complex to the xenobiotic response element and induces "
        "CYP1A1 and CYP1A2"
    ),
    pubchem_aid=743085,
    assay_name=(
        "qHTS assay to identify small molecule that activate the aryl hydrocarbon "
        "receptor (AhR) signaling pathway"
    ),
)

CYTOTOX = Endpoint(
    name="VIABILITY",
    pathway="aryl hydrocarbon receptor / xenobiotic induction",
    mie="loss of cell viability, which registers on the reporter readout as a consequence",
    pubchem_aid=743086,
    assay_name=(
        "qHTS assay to identify small molecule that activate the aryl hydrocarbon "
        "receptor (AhR) signaling pathway - cell viability counter screen"
    ),
)

TARGETS: dict[str, Endpoint] = {"AHR": TARGET, "VIABILITY": CYTOTOX}


def get(name: str) -> Endpoint:
    try:
        return TARGETS[name.upper()]
    except KeyError:
        raise KeyError(f"unknown target {name!r}; known: {sorted(TARGETS)}") from None


def all_names() -> list[str]:
    return sorted(TARGETS)
