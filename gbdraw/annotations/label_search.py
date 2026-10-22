"""Where an annotation label goes along its region when its default place conflicts (OV-845, both modes).

A label keeps its default place when nothing conflicts there. Otherwise it
takes the candidate place with the fewest conflicts, nearest the default,
the smaller one on a tie. Each mode supplies its candidates (places that keep
the label's centre over one of its own marks) and its conflict count.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable


def nearest_fewest_conflicts(default: float, candidates: Iterable[float], conflicts: Callable[[float], int]) -> float:
    """``default`` when it has no conflict; else the candidate with the fewest conflicts, nearest ``default``,
    the smaller on a tie. Candidates are tried nearest first, so the search stops at the first clear one."""

    best_count, best = -1, default
    for value in sorted({default, *candidates}, key=lambda value: (abs(value - default), value)):
        count = conflicts(value)
        if count == 0:
            return value
        if best_count < 0 or count < best_count:
            best_count, best = count, value
    return best


__all__ = ["nearest_fewest_conflicts"]
