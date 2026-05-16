from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

ScalarUnit = Literal["px", "factor"]


@dataclass(frozen=True)
class ScalarSpec:
    """A scalar value with an explicit unit."""

    value: float
    unit: ScalarUnit = "factor"

    @classmethod
    def parse(cls, raw: str) -> "ScalarSpec":
        s = str(raw).strip()
        if not s:
            raise ValueError("empty scalar")
        if s.endswith("px"):
            return cls(value=float(s[:-2]), unit="px")
        if s.endswith("%"):
            return cls(value=float(s[:-1]) / 100.0, unit="factor")
        return cls(value=float(s), unit="factor")

    def resolve(self, reference: float) -> float:
        return self.value if self.unit == "px" else (self.value * reference)


__all__ = ["ScalarSpec", "ScalarUnit"]
