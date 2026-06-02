"""Result objects and helpers for afdb-query."""

from __future__ import annotations


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
