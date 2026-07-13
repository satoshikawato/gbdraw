from __future__ import annotations

from dataclasses import dataclass


@dataclass
class CircularTrackSlotParseError(ValueError):
    message: str
    raw: str

    def __str__(self) -> str:  # pragma: no cover
        return f"{self.message}: {self.raw!r}"


def parse_bool(raw: object) -> bool:
    s = str(raw).strip().lower()
    if s in {"1", "true", "yes", "on"}:
        return True
    if s in {"0", "false", "no", "off"}:
        return False
    raise ValueError(f"invalid bool: {raw!r}")


def split_kv_list(raw: str) -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = []
    for part in raw.split(","):
        part = part.strip()
        if not part:
            continue
        if "=" not in part:
            raise ValueError(f"expected key=value, got: {part!r}")
        key, value = part.split("=", 1)
        out.append((key.strip(), value.strip()))
    return out


def normalize_dinucleotide_skew_color_params(params: dict) -> dict:
    if "positive_color" not in params and "high_color" in params:
        params["positive_color"] = params["high_color"]
    if "negative_color" not in params and "low_color" in params:
        params["negative_color"] = params["low_color"]
    params.pop("high_color", None)
    params.pop("low_color", None)
    return params


def strip_inline_comment(raw: str) -> str:
    for marker in (" #", "\t#"):
        index = raw.find(marker)
        if index >= 0:
            return raw[:index].strip()
    return raw


__all__ = [
    "CircularTrackSlotParseError",
    "parse_bool",
    "split_kv_list",
]
