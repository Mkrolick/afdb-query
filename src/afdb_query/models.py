"""Result objects and helpers for afdb-query."""

from __future__ import annotations

from dataclasses import dataclass, field

from .plddt import is_contiguous as _is_contiguous
from .plddt import mean_plddt as _mean_plddt
from .selection import is_monomer as _is_monomer


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
    def from_dict(cls, data: dict) -> Plddt:
        return cls(
            scores=data["confidenceScore"],
            residue_numbers=data["residueNumber"],
            raw=data,
        )

    def mean(self, start: int | None = None, stop: int | None = None) -> float | None:
        """Mean pLDDT over ``scores[start:stop]``; ``None`` when that slice is empty.

        Delegates to :func:`afdb_query.plddt.mean_plddt` rather than re-implementing the
        average, so this and every other mean in a codebase are the same computation.
        Note it may differ in the last digit from :attr:`Structure.global_plddt`, which
        AFDB rounds server-side.
        """
        return _mean_plddt(self.scores, start, stop)

    @property
    def is_contiguous(self) -> bool:
        """Whether ``residue_numbers`` is exactly ``1..N``, so positions are residues."""
        return _is_contiguous(self.residue_numbers)


@dataclass(frozen=True)
class Structure:
    """One AFDB structure match for a queried sequence.

    Thin typed wrapper over the endpoint's ``summary`` dict. ``raw`` is the full
    summary (escape hatch). ``plddt()`` lazily fetches per-residue pLDDT.
    """

    raw: dict
    _client: AlphaFold = field(repr=False, compare=False)  # noqa: F821
    _cache: dict = field(default_factory=dict, repr=False, compare=False)

    @property
    def model_identifier(self) -> str | None:
        return self.raw.get("model_identifier")

    @property
    def model_url(self) -> str | None:
        return self.raw.get("model_url")

    @property
    def global_plddt(self) -> float | None:
        """AFDB's ``confidence_avg_local_score``: mean pLDDT over the WHOLE deposited model.

        For a complex that average spans every chain, not only the one matching the
        query, so it is not comparable with a monomer's. Check :attr:`oligomeric_state`
        (or ``afdb_query.is_monomer``) before comparing this across structures.
        """
        return self.raw.get("confidence_avg_local_score")

    @property
    def oligomeric_state(self) -> str | None:
        """``"MONOMER"``, ``"HOMODIMER"``, ``"HETERODIMER"``, ... or None if unset."""
        return self.raw.get("oligomeric_state")

    @property
    def is_monomer(self) -> bool:
        """Whether this is a single-chain model of one entity. See ``selection.is_monomer``."""
        return _is_monomer(self.raw)

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
            data = self._client.fetch_confidence(self.model_url)
            self._cache["plddt"] = Plddt.from_dict(data)
        return self._cache["plddt"]
