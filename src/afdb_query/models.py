"""Result objects and helpers for afdb-query."""

from __future__ import annotations

from dataclasses import dataclass, field


def confidence_url(model_url: str) -> str:
    """Derive the per-residue confidence-JSON URL from a model (CIF) URL.

    AFDB names files ``...-model_vN.cif`` and ``...-confidence_vN.json`` in the
    same directory, so the per-residue pLDDT URL is a pure string transform of
    the model URL (verified against the live API).
    """
    url = model_url.replace("-model_", "-confidence_")
    if url.endswith(".bcif"):
        return url[: -len(".bcif")] + ".json"
    if url.endswith(".cif"):
        return url[: -len(".cif")] + ".json"
    return url


@dataclass(frozen=True)
class Plddt:
    """Per-residue pLDDT for one structure.

    ``scores`` and ``residue_numbers`` are parallel lists. ``raw`` is the full
    confidence-JSON document (escape hatch).
    """

    scores: list[float]
    residue_numbers: list[int]
    raw: dict = field(repr=False)

    @classmethod
    def from_dict(cls, data: dict) -> "Plddt":
        return cls(
            scores=data["confidenceScore"],
            residue_numbers=data["residueNumber"],
            raw=data,
        )

    def first(self, n: int) -> list[float]:
        """First ``n`` per-residue pLDDT values, or all of them if fewer than ``n``.

        Never pads and never raises on short structures: returns ``scores[:n]``.
        """
        return self.scores[:n]
