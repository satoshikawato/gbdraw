"""The error rule of the warm-equals-cold tests (``error_rule_failures``)."""

from __future__ import annotations

from tests.utils.render_stage_equivalence import (
    RenderOutcome,
    StepResult,
    equivalence_failures,
    error_rule_failures,
)

_OK = RenderOutcome(reply={"svg": "x"}, error=None, metrics={})
_INVALID = RenderOutcome(reply=None, error=("ValidationError", "bad input"), metrics={}, user_input_error=True)
_DEFECT = RenderOutcome(reply=None, error=("KeyError", "'x'"), metrics={})


def _step(label: str, outcome: RenderOutcome) -> StepResult:
    return StepResult(label, (), outcome, outcome)


def test_an_identical_user_input_error_passes() -> None:
    assert not error_rule_failures(("base", _OK), [("a", _OK), ("b", _INVALID)])


def test_an_identical_other_error_fails() -> None:
    failures = equivalence_failures([_step("base", _OK), _step("a", _OK), _step("b", _DEFECT)])
    assert failures == ["b: render failed with ('KeyError', \"'x'\")"]


def test_the_base_must_render() -> None:
    assert error_rule_failures(("base", _INVALID), [("a", _OK)]) == [
        "base: the base state must render, got ('ValidationError', 'bad input')"
    ]
    assert error_rule_failures(("base", None), [("a", _OK)]) == ["base: the base state must render, got None"]


def test_every_step_failing_fails() -> None:
    assert error_rule_failures(("base", _OK), [("a", _INVALID), ("b", _INVALID)]) == [
        "all 2 steps failed with a user-input validation error"
    ]
