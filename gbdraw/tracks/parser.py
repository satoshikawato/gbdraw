from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable

from .spec import (
    CircularTrackPlacement,
    LayoutMode,
    LinearTrackPlacement,
    ScalarSpec,
    TrackKind,
    TrackSpec,
)


@dataclass(frozen=True)
class TrackSpecParseError(ValueError):
    message: str
    raw: str

    def __str__(self) -> str:  # pragma: no cover
        return f"{self.message}: {self.raw!r}"


def _parse_bool(raw: str) -> bool:
    s = str(raw).strip().lower()
    if s in {"1", "true", "yes", "on"}:
        return True
    if s in {"0", "false", "no", "off"}:
        return False
    raise ValueError(f"invalid bool: {raw!r}")


def _split_kv_list(raw: str) -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = []
    for part in raw.split(","):
        part = part.strip()
        if not part:
            continue
        if "=" not in part:
            raise ValueError(f"expected key=value, got: {part!r}")
        k, v = part.split("=", 1)
        out.append((k.strip(), v.strip()))
    return out


def parse_track_spec(raw: str, *, mode: LayoutMode) -> TrackSpec:
    """
    Parse a single track rule string into a `TrackSpec`.

    Minimal grammar (subject to change):
      <kind>[:<id>] [@ <k>=<v>(,<k>=<v>)*]

    Examples:
      - "features:genes@z=10"
      - "gc_content@r=0.78,w=0.08"
      - "legend@show=false"
      - "labels@show=true"
    """
    original = raw
    s = str(raw).strip()
    if not s or s.startswith("#"):
        raise TrackSpecParseError("empty/comment line", original)

    # Strip inline comments.
    if "#" in s:
        s = s.split("#", 1)[0].strip()

    if "@" in s:
        head, opts = s.split("@", 1)
        opts = opts.strip()
    else:
        head, opts = s, ""

    head = head.strip()
    if not head:
        raise TrackSpecParseError("missing track head", original)

    if ":" in head:
        kind_raw, id_raw = head.split(":", 1)
        kind = kind_raw.strip()
        track_id = id_raw.strip()
    else:
        kind = head.strip()
        track_id = kind

    if not kind:
        raise TrackSpecParseError("missing kind", original)
    if not track_id:
        track_id = kind

    # NOTE: we keep this permissive for now; unknown kinds can be represented as "custom".
    kind_norm = kind.lower()
    kind_typed: TrackKind = kind_norm if kind_norm in TrackKind.__args__ else "custom"  # type: ignore[attr-defined]

    placement: CircularTrackPlacement | LinearTrackPlacement | None = None
    show = True
    params: dict[str, Any] = {}

    if opts:
        try:
            for k, v in _split_kv_list(opts):
                key = k.lower()
                if key in {"id"}:
                    track_id = v
                elif key in {"show", "visible"}:
                    show = _parse_bool(v)
                elif key in {"z", "z_index", "zindex"}:
                    params["z"] = int(v)
                elif mode == "circular":
                    if placement is None or not isinstance(placement, CircularTrackPlacement):
                        placement = CircularTrackPlacement()
                    if key in {"r", "radius"}:
                        placement = CircularTrackPlacement(
                            radius=ScalarSpec.parse(v),
                            inner_radius=placement.inner_radius,
                            outer_radius=placement.outer_radius,
                            width=placement.width,
                            z=placement.z,
                        )
                    elif key in {"ri", "inner", "inner_radius"}:
                        placement = CircularTrackPlacement(
                            radius=placement.radius,
                            inner_radius=ScalarSpec.parse(v),
                            outer_radius=placement.outer_radius,
                            width=placement.width,
                            z=placement.z,
                        )
                    elif key in {"ro", "outer", "outer_radius"}:
                        placement = CircularTrackPlacement(
                            radius=placement.radius,
                            inner_radius=placement.inner_radius,
                            outer_radius=ScalarSpec.parse(v),
                            width=placement.width,
                            z=placement.z,
                        )
                    elif key in {"w", "width"}:
                        placement = CircularTrackPlacement(
                            radius=placement.radius,
                            inner_radius=placement.inner_radius,
                            outer_radius=placement.outer_radius,
                            width=ScalarSpec.parse(v),
                            z=placement.z,
                        )
                    else:
                        params[key] = v
                else:  # linear
                    if placement is None or not isinstance(placement, LinearTrackPlacement):
                        placement = LinearTrackPlacement()
                    if key in {"y"}:
                        placement = LinearTrackPlacement(
                            y=ScalarSpec.parse(v),
                            height=placement.height,
                            z=placement.z,
                        )
                    elif key in {"h", "height"}:
                        placement = LinearTrackPlacement(
                            y=placement.y,
                            height=ScalarSpec.parse(v),
                            z=placement.z,
                        )
                    else:
                        params[key] = v

            # Apply z if present.
            if "z" in params:
                z = int(params.pop("z"))
                if mode == "circular":
                    placement = (placement if isinstance(placement, CircularTrackPlacement) else CircularTrackPlacement())
                    placement = CircularTrackPlacement(
                        radius=placement.radius,
                        inner_radius=placement.inner_radius,
                        outer_radius=placement.outer_radius,
                        width=placement.width,
                        z=z,
                    )
                else:
                    placement = (placement if isinstance(placement, LinearTrackPlacement) else LinearTrackPlacement())
                    placement = LinearTrackPlacement(y=placement.y, height=placement.height, z=z)
        except Exception as exc:
            raise TrackSpecParseError(str(exc), original) from exc

    return TrackSpec(
        id=track_id,
        kind=kind_typed,
        mode=mode,
        placement=placement,
        show=show,
        params=params or None,
    )


def parse_track_specs(specs: Iterable[str], *, mode: LayoutMode) -> list[TrackSpec]:
    """Parse multiple rule strings into a list of `TrackSpec`."""
    out: list[TrackSpec] = []
    seen: set[str] = set()
    for raw in specs:
        s = str(raw).strip()
        if not s or s.startswith("#"):
            continue
        ts = parse_track_spec(s, mode=mode)
        if ts.id in seen:
            raise TrackSpecParseError("duplicate track id", ts.id)
        seen.add(ts.id)
        out.append(ts)
    return out


__all__ = [
    "TrackSpecParseError",
    "parse_track_spec",
    "parse_track_specs",
]


