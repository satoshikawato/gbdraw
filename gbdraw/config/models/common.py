from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping


@dataclass(frozen=True)
class ShortLongFloatConfig:
    short: float
    long: float

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "ShortLongFloatConfig":
        return cls(short=float(d["short"]), long=float(d["long"]))

    def for_length_param(self, length_param: str) -> float:
        return getattr(self, length_param)


__all__ = ["ShortLongFloatConfig"]


