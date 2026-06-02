def test_public_exports():
    import afdb_query as m

    for name in [
        "AlphaFold",
        "Structure",
        "Plddt",
        "filter_reason",
        "confidence_url",
        "AFDBError",
        "InvalidSequenceError",
    ]:
        assert hasattr(m, name), f"missing export: {name}"
