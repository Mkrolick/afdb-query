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


@dataclass(frozen=True)
class Structure:
    """One AFDB structure match for a queried sequence.

    Thin typed wrapper over the endpoint's ``summary`` dict. ``raw`` is the full
    summary (escape hatch). ``plddt()`` lazily fetches per-residue pLDDT.
    """

    raw: dict
    _client: "AlphaFold" = field(repr=False, compare=False)  # noqa: F821
    _cache: dict = field(default_factory=dict, repr=False, compare=False)

    @property
    def model_identifier(self) -> str | None:
        return self.raw.get("model_identifier")

    @property
    def model_url(self) -> str | None:
        return self.raw.get("model_url")

    @property
    def global_plddt(self) -> float | None:
        return self.raw.get("confidence_avg_local_score")

    @property
    def sequence_identity(self) -> float | None:
        return self.raw.get("sequence_identity")

    @property
    def coverage(self) -> float | None:
        return self.raw.get("coverage")

    @property
    def uniprot_accession(self) -> str | None:
        for entity in self.raw.get("entities") or []:
            if entity.get("identifier_category") == "UNIPROT":
                return entity.get("identifier")
        return None

    @property
    def description(self) -> str | None:
        for entity in self.raw.get("entities") or []:
            if entity.get("description"):
                return entity["description"]
        return None

    def plddt(self) -> Plddt:
        """Tier 2: per-residue pLDDT for this structure (fetched once, then cached)."""
        if "plddt" not in self._cache:
            data = self._client._fetch_confidence(self.model_url)
            self._cache["plddt"] = Plddt.from_dict(data)
        return self._cache["plddt"]
