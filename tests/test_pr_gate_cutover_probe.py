def test_pr_gate_blocks_a_failed_fast_job() -> None:
    """Intentional failure for the isolated hosted branch-protection probe."""
    raise AssertionError("intentional PR / gate cutover probe; do not merge")
