#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from pathlib import Path

SPACING = {
    "0": "0",
    "0.5": "0.125rem",
    "1": "0.25rem",
    "1.5": "0.375rem",
    "2": "0.5rem",
    "2.5": "0.625rem",
    "3": "0.75rem",
    "4": "1rem",
    "5": "1.25rem",
    "6": "1.5rem",
    "7": "1.75rem",
    "8": "2rem",
    "10": "2.5rem",
    "12": "3rem",
    "14": "3.5rem",
    "16": "4rem",
    "20": "5rem",
    "24": "6rem",
    "32": "8rem",
    "40": "10rem",
    "48": "12rem",
    "72": "18rem",
    "80": "20rem",
}

FRACTIONS = {
    "1/2": "50%",
    "1/3": "33.333333%",
    "1/4": "25%",
    "full": "100%",
}

COLORS = {
    "transparent": "transparent",
    "white": "#ffffff",
    "black": "#000000",
    "slate-50": "#f8fafc",
    "slate-100": "#f1f5f9",
    "slate-200": "#e2e8f0",
    "slate-300": "#cbd5e1",
    "slate-400": "#94a3b8",
    "slate-500": "#64748b",
    "slate-600": "#475569",
    "slate-700": "#334155",
    "slate-800": "#1e293b",
    "blue-50": "#eff6ff",
    "blue-100": "#dbeafe",
    "blue-200": "#bfdbfe",
    "blue-300": "#93c5fd",
    "blue-400": "#60a5fa",
    "blue-500": "#3b82f6",
    "blue-600": "#2563eb",
    "blue-700": "#1d4ed8",
    "blue-800": "#1e40af",
    "green-50": "#f0fdf4",
    "green-100": "#dcfce7",
    "green-200": "#bbf7d0",
    "green-300": "#86efac",
    "green-500": "#22c55e",
    "green-600": "#16a34a",
    "green-700": "#15803d",
    "purple-50": "#faf5ff",
    "purple-100": "#f3e8ff",
    "purple-200": "#e9d5ff",
    "purple-300": "#d8b4fe",
    "purple-600": "#9333ea",
    "purple-700": "#7e22ce",
    "purple-800": "#6b21a8",
    "red-50": "#fef2f2",
    "red-100": "#fee2e2",
    "red-200": "#fecaca",
    "red-400": "#f87171",
    "red-500": "#ef4444",
    "red-600": "#dc2626",
    "red-700": "#b91c1c",
    "red-800": "#991b1b",
    "amber-700": "#b45309",
    "teal-100": "#ccfbf1",
    "teal-200": "#99f6e4",
    "teal-800": "#115e59",
    "gray-100": "#f3f4f6",
    "gray-300": "#d1d5db",
    "indigo-200": "#c7d2fe",
    "indigo-600": "#4f46e5",
    "indigo-700": "#4338ca",
    "indigo-800": "#3730a3",
}

TEXT_SIZES = {
    "xs": "font-size:0.75rem;line-height:1rem",
    "sm": "font-size:0.875rem;line-height:1.25rem",
    "base": "font-size:1rem;line-height:1.5rem",
    "lg": "font-size:1.125rem;line-height:1.75rem",
    "xl": "font-size:1.25rem;line-height:1.75rem",
    "2xl": "font-size:1.5rem;line-height:2rem",
    "5xl": "font-size:3rem;line-height:1",
    "6xl": "font-size:3.75rem;line-height:1",
}

FONT_WEIGHTS = {
    "normal": "400",
    "medium": "500",
    "semibold": "600",
    "bold": "700",
}

LEADING = {
    "none": "1",
    "tight": "1.25",
    "relaxed": "1.625",
}

TRACKING = {
    "tight": "-0.025em",
    "wide": "0.025em",
    "wider": "0.05em",
}

TRANSITIONS = {
    "transition-all": "transition-property:all;transition-timing-function:cubic-bezier(.4,0,.2,1);transition-duration:0.15s",
    "transition-colors": "transition-property:color,background-color,border-color,text-decoration-color,fill,stroke;transition-timing-function:cubic-bezier(.4,0,.2,1);transition-duration:0.15s",
    "transition-opacity": "transition-property:opacity;transition-timing-function:cubic-bezier(.4,0,.2,1);transition-duration:0.15s",
    "transition-transform": "transition-property:transform;transition-timing-function:cubic-bezier(.4,0,.2,1);transition-duration:0.15s",
}

SHADOWS = {
    "shadow": "box-shadow:0 1px 3px 0 rgb(0 0 0 / 0.1),0 1px 2px -1px rgb(0 0 0 / 0.1)",
    "shadow-sm": "box-shadow:0 1px 2px 0 rgb(0 0 0 / 0.05)",
    "shadow-sm/50": "box-shadow:0 1px 2px 0 rgb(0 0 0 / 0.08)",
    "shadow-md": "box-shadow:0 4px 6px -1px rgb(0 0 0 / 0.1),0 2px 4px -2px rgb(0 0 0 / 0.1)",
    "shadow-lg": "box-shadow:0 10px 15px -3px rgb(0 0 0 / 0.1),0 4px 6px -4px rgb(0 0 0 / 0.1)",
    "shadow-xl": "box-shadow:0 20px 25px -5px rgb(0 0 0 / 0.1),0 8px 10px -6px rgb(0 0 0 / 0.1)",
    "shadow-2xl": "box-shadow:0 25px 50px -12px rgb(0 0 0 / 0.25)",
    "shadow-none": "box-shadow:0 0 #0000",
}

ROUNDED = {
    "rounded": "border-radius:0.25rem",
    "rounded-md": "border-radius:0.375rem",
    "rounded-lg": "border-radius:0.5rem",
    "rounded-xl": "border-radius:0.75rem",
    "rounded-3xl": "border-radius:1.5rem",
    "rounded-full": "border-radius:9999px",
    "rounded-l-lg": "border-top-left-radius:0.5rem;border-bottom-left-radius:0.5rem",
}

DISPLAY = {
    "block": "display:block",
    "inline-block": "display:inline-block",
    "inline-flex": "display:inline-flex",
    "flex": "display:flex",
    "grid": "display:grid",
    "hidden": "display:none",
}

POSITION = {
    "relative": "position:relative",
    "absolute": "position:absolute",
    "fixed": "position:fixed",
    "sticky": "position:sticky",
}

SIMPLE = {
    "flex-col": "flex-direction:column",
    "flex-wrap": "flex-wrap:wrap",
    "flex-1": "flex:1 1 0%",
    "flex-grow": "flex-grow:1",
    "shrink-0": "flex-shrink:0",
    "items-center": "align-items:center",
    "items-start": "align-items:flex-start",
    "items-end": "align-items:flex-end",
    "justify-center": "justify-content:center",
    "justify-between": "justify-content:space-between",
    "justify-start": "justify-content:flex-start",
    "overflow-hidden": "overflow:hidden",
    "overflow-auto": "overflow:auto",
    "overflow-y-auto": "overflow-y:auto",
    "text-left": "text-align:left",
    "text-center": "text-align:center",
    "font-mono": "font-family:ui-monospace,SFMono-Regular,Menlo,Monaco,Consolas,Liberation Mono,Courier New,monospace",
    "uppercase": "text-transform:uppercase",
    "italic": "font-style:italic",
    "truncate": "overflow:hidden;text-overflow:ellipsis;white-space:nowrap",
    "whitespace-nowrap": "white-space:nowrap",
    "whitespace-pre-wrap": "white-space:pre-wrap",
    "break-all": "word-break:break-all",
    "select-none": "user-select:none",
    "select-all": "user-select:all",
    "cursor-pointer": "cursor:pointer",
    "cursor-not-allowed": "cursor:not-allowed",
    "cursor-col-resize": "cursor:col-resize",
    "cursor-help": "cursor:help",
    "pointer-events-none": "pointer-events:none",
    "origin-top": "transform-origin:top",
    "outline-none": "outline:2px solid transparent;outline-offset:2px",
    "underline": "text-decoration-line:underline",
    "antialiased": "-webkit-font-smoothing:antialiased;-moz-osx-font-smoothing:grayscale",
    "leading-none": "line-height:1",
    "leading-tight": "line-height:1.25",
    "leading-relaxed": "line-height:1.625",
    "align-text-bottom": "vertical-align:text-bottom",
    "border": "border-width:1px;border-style:solid",
    "border-0": "border-width:0",
    "border-2": "border-width:2px",
    "border-4": "border-width:4px",
    "border-b": "border-bottom-width:1px",
    "border-b-4": "border-bottom-width:4px",
    "border-l": "border-left-width:1px",
    "border-l-4": "border-left-width:4px",
    "border-t": "border-top-width:1px",
    "border-none": "border-style:none",
    "border-dashed": "border-style:dashed",
    "border-t-transparent": "border-top-color:transparent",
    "list-none": "list-style-type:none",
    "animate-spin": "animation:tw-spin 1s linear infinite",
    "animate-pulse": "animation:tw-pulse 2s cubic-bezier(0.4,0,0.6,1) infinite",
    "backdrop-blur-sm": "backdrop-filter:blur(4px)",
    "backdrop-blur-md": "backdrop-filter:blur(12px)",
    "bg-gradient-to-r": "background-image:linear-gradient(to right,var(--tw-gradient-stops))",
    "from-blue-600": "--tw-gradient-from:#2563eb;--tw-gradient-to:rgb(37 99 235 / 0);--tw-gradient-stops:var(--tw-gradient-from),var(--tw-gradient-to)",
    "to-indigo-600": "--tw-gradient-to:#4f46e5",
    "from-blue-700": "--tw-gradient-from:#1d4ed8;--tw-gradient-to:rgb(29 78 216 / 0);--tw-gradient-stops:var(--tw-gradient-from),var(--tw-gradient-to)",
    "to-indigo-700": "--tw-gradient-to:#4338ca",
    "bg-contain": "background-size:contain",
    "bg-no-repeat": "background-repeat:no-repeat",
    "ring-0": "box-shadow:0 0 #0000",
}

VARIANTS = {"hover", "focus", "active", "disabled", "group-hover"}

CLASS_TOKEN_RE = re.compile(r"^-?[A-Za-z0-9][A-Za-z0-9_:/\-\[\]\.()%#]+$")


class BuildError(RuntimeError):
    pass


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("-i", dest="input_css", required=True)
    parser.add_argument("-o", dest="output_css", required=True)
    parser.add_argument("--config", dest="config", default=None)
    parser.add_argument("--minify", action="store_true")
    parser.add_argument("-h", "--help", action="store_true")
    args, _ = parser.parse_known_args()
    if args.help:
        parser.print_help()
        raise SystemExit(0)
    return args


def extract_js_array_literal(text: str, key: str) -> str | None:
    key_idx = text.find(key)
    if key_idx == -1:
        return None
    start = text.find("[", key_idx)
    if start == -1:
        return None

    depth = 0
    in_str: str | None = None
    escaped = False
    for idx in range(start, len(text)):
        ch = text[idx]
        if in_str is not None:
            if escaped:
                escaped = False
                continue
            if ch == "\\":
                escaped = True
                continue
            if ch == in_str:
                in_str = None
            continue

        if ch in {"'", '"'}:
            in_str = ch
            continue
        if ch == "[":
            depth += 1
            continue
        if ch == "]":
            depth -= 1
            if depth == 0:
                return text[start + 1 : idx]

    return None


def parse_js_string_list(array_src: str | None) -> list[str]:
    if not array_src:
        return []

    values: list[str] = []
    pattern = re.compile(r"'([^'\\]*(?:\\.[^'\\]*)*)'|\"([^\"\\]*(?:\\.[^\"\\]*)*)\"")
    for match in pattern.finditer(array_src):
        value = match.group(1) if match.group(1) is not None else match.group(2)
        if value is None:
            continue
        values.append(bytes(value, "utf-8").decode("unicode_escape"))
    return values


def load_tailwind_config(config_path: Path | None) -> tuple[list[Path], list[str]]:
    if config_path is None or not config_path.exists():
        return [], []

    raw = config_path.read_text(encoding="utf-8")
    content_items = parse_js_string_list(extract_js_array_literal(raw, "content"))
    safelist_items = parse_js_string_list(extract_js_array_literal(raw, "safelist"))

    config_dir = config_path.parent
    files: list[Path] = []
    for item in content_items:
        pattern = item.strip()
        if not pattern:
            continue
        path_pattern = (config_dir / pattern).resolve()
        if "*" in pattern or "?" in pattern or "[" in pattern:
            files.extend(sorted(p for p in config_dir.glob(pattern) if p.is_file()))
        elif path_pattern.is_file():
            files.append(path_pattern)

    deduped_files = sorted({f.resolve() for f in files})
    return deduped_files, safelist_items


def normalize_token(token: str) -> str:
    cleaned = token.strip().strip(",;")
    cleaned = cleaned.strip("'\"")
    return cleaned


def is_class_token(token: str) -> bool:
    if not token:
        return False
    if token.startswith("ph-") or token == "ph":
        return False
    return CLASS_TOKEN_RE.match(token) is not None


def looks_like_utility(token: str) -> bool:
    if token in DISPLAY or token in POSITION or token in SIMPLE:
        return True
    if re.match(r"^-?m[trblxy]?-.+$", token):
        return True
    if re.match(r"^p[trblxy]?-.+$", token):
        return True

    prefixes = (
        "text-",
        "bg-",
        "border",
        "ring-",
        "shadow",
        "rounded",
        "w-",
        "h-",
        "min-",
        "max-",
        "gap",
        "space-",
        "grid",
        "col-",
        "top-",
        "left-",
        "right-",
        "bottom-",
        "inset-",
        "z-",
        "hover:",
        "focus:",
        "active:",
        "disabled:",
        "group-hover:",
        "translate",
        "-translate",
        "scale-",
        "transition",
        "duration-",
        "cursor-",
        "overflow",
        "items-",
        "justify-",
        "font-",
        "leading-",
        "tracking-",
        "whitespace-",
        "break-",
        "select-",
        "pointer-events-",
        "backdrop-",
        "animate-",
        "origin-",
        "outline-",
    )
    return token.startswith(prefixes)


def extract_class_tokens_from_content(text: str) -> set[str]:
    tokens: set[str] = set()

    for match in re.finditer(r"\bclass\s*=\s*(['\"])(.*?)\1", text, flags=re.S):
        for token in match.group(2).split():
            token = normalize_token(token)
            if is_class_token(token):
                tokens.add(token)

    for match in re.finditer(r"\b:class\s*=\s*(['\"])(.*?)\1", text, flags=re.S):
        expression = match.group(2)
        for sm in re.finditer(r"'([^'\\]*(?:\\.[^'\\]*)*)'|\"([^\"\\]*(?:\\.[^\"\\]*)*)\"", expression):
            value = sm.group(1) if sm.group(1) is not None else sm.group(2)
            if value is None:
                continue
            for token in value.split():
                token = normalize_token(token)
                if is_class_token(token):
                    tokens.add(token)

    for match in re.finditer(r"classList\.(?:add|remove|toggle)\((.*?)\)", text, flags=re.S):
        args = match.group(1)
        for sm in re.finditer(r"'([^'\\]*(?:\\.[^'\\]*)*)'|\"([^\"\\]*(?:\\.[^\"\\]*)*)\"", args):
            value = sm.group(1) if sm.group(1) is not None else sm.group(2)
            if value is None:
                continue
            for token in value.split():
                token = normalize_token(token)
                if is_class_token(token):
                    tokens.add(token)

    for match in re.finditer(r"setAttribute\(\s*['\"]class['\"]\s*,\s*(['\"])(.*?)\1\s*\)", text, flags=re.S):
        for token in match.group(2).split():
            token = normalize_token(token)
            if is_class_token(token):
                tokens.add(token)

    return tokens


def parse_opacity(value: str) -> float | None:
    if value.endswith("%"):
        value = value[:-1]
    if not value:
        return None
    if re.fullmatch(r"\d+", value):
        number = int(value)
        if 0 <= number <= 100:
            return number / 100.0
        return None
    try:
        number = float(value)
        if 0.0 <= number <= 1.0:
            return number
    except ValueError:
        return None
    return None


def color_to_rgba(hex_color: str, opacity: float) -> str:
    hex_color = hex_color.lstrip("#")
    if len(hex_color) != 6:
        return f"rgb(0 0 0 / {opacity})"
    r = int(hex_color[0:2], 16)
    g = int(hex_color[2:4], 16)
    b = int(hex_color[4:6], 16)
    return f"rgb({r} {g} {b} / {opacity:g})"


def resolve_color(token: str) -> str | None:
    if "/" in token:
        color_key, alpha_key = token.split("/", 1)
        base = COLORS.get(color_key)
        alpha = parse_opacity(alpha_key)
        if base is None or alpha is None:
            return None
        if base == "transparent":
            return "transparent"
        if base.startswith("#"):
            return color_to_rgba(base, alpha)
        return base
    return COLORS.get(token)


def resolve_size(value: str) -> str | None:
    if value == "auto":
        return "auto"
    if value in FRACTIONS:
        return FRACTIONS[value]
    if value in SPACING:
        return SPACING[value]
    if value.startswith("[") and value.endswith("]"):
        return value[1:-1]
    return None


def negate_value(value: str) -> str:
    if value.startswith("-"):
        return value[1:]
    if value == "0":
        return value
    return f"-{value}"


def resolve_text(token: str) -> str | None:
    rest = token[len("text-") :]
    if rest in TEXT_SIZES:
        return TEXT_SIZES[rest]
    if rest.startswith("[") and rest.endswith("]"):
        value = rest[1:-1]
        return f"font-size:{value}"
    color = resolve_color(rest)
    if color:
        return f"color:{color}"
    return None


def resolve_background(token: str) -> str | None:
    rest = token[len("bg-") :]
    color = resolve_color(rest)
    if color:
        return f"background-color:{color}"
    return None


def resolve_border(token: str) -> str | None:
    if token == "border":
        return SIMPLE[token]

    rest = token[len("border-") :]
    if rest in {"0", "2", "4"}:
        return f"border-width:{rest}px" if rest != "0" else "border-width:0"
    if rest == "b":
        return "border-bottom-width:1px"
    if rest == "b-0":
        return "border-bottom-width:0"
    if rest == "b-4":
        return "border-bottom-width:4px"
    if rest == "l":
        return "border-left-width:1px"
    if rest == "l-4":
        return "border-left-width:4px"
    if rest == "t":
        return "border-top-width:1px"
    if rest == "none":
        return "border-style:none"
    if rest == "dashed":
        return "border-style:dashed"
    if rest == "t-transparent":
        return "border-top-color:transparent"

    if rest.startswith("l-"):
        color = resolve_color(rest[2:])
        if color:
            return f"border-left-color:{color}"
    if rest.startswith("t-"):
        color = resolve_color(rest[2:])
        if color:
            return f"border-top-color:{color}"

    color = resolve_color(rest)
    if color:
        return f"border-color:{color}"
    return None


def resolve_spacing_utility(token: str) -> str | None:
    match = re.match(r"^(-?)([mp])([trblxy]?)-(.+)$", token)
    if not match:
        return None

    negative = match.group(1) == "-"
    mode = match.group(2)
    axis = match.group(3)
    raw_value = match.group(4)
    value = resolve_size(raw_value)
    if value is None:
        return None

    if mode == "p" and negative:
        return None
    if negative:
        value = negate_value(value)

    prop = "margin" if mode == "m" else "padding"
    if axis == "":
        return f"{prop}:{value}"
    if axis == "x":
        return f"{prop}-left:{value};{prop}-right:{value}"
    if axis == "y":
        return f"{prop}-top:{value};{prop}-bottom:{value}"
    axis_map = {"t": "top", "r": "right", "b": "bottom", "l": "left"}
    return f"{prop}-{axis_map[axis]}:{value}"


def resolve_size_utility(token: str) -> str | None:
    for prefix, prop in (
        ("w-", "width"),
        ("h-", "height"),
        ("min-w-", "min-width"),
        ("max-w-", "max-width"),
        ("min-h-", "min-height"),
        ("max-h-", "max-height"),
    ):
        if not token.startswith(prefix):
            continue
        raw = token[len(prefix) :]
        if raw == "full":
            return f"{prop}:100%"
        if raw == "screen" and prop == "height":
            return "height:100vh"
        if raw == "md" and prop == "max-width":
            return "max-width:28rem"
        value = resolve_size(raw)
        if value is not None:
            return f"{prop}:{value}"
    return None


def resolve_position_utility(token: str) -> str | None:
    if token == "inset-0":
        return "top:0;right:0;bottom:0;left:0"

    for prefix, prop in (
        ("top-", "top"),
        ("right-", "right"),
        ("bottom-", "bottom"),
        ("left-", "left"),
    ):
        if token.startswith(prefix):
            value = resolve_size(token[len(prefix) :])
            if value is not None:
                return f"{prop}:{value}"
    return None


def resolve_grid_utility(token: str) -> str | None:
    if token.startswith("grid-cols-"):
        count = token[len("grid-cols-") :]
        if count.isdigit():
            return f"grid-template-columns:repeat({int(count)},minmax(0,1fr))"
    if token.startswith("col-span-"):
        count = token[len("col-span-") :]
        if count.isdigit():
            return f"grid-column:span {int(count)} / span {int(count)}"
    return None


def resolve_gap_utility(token: str) -> str | None:
    match = re.match(r"^gap(?:-([xy]))?-(.+)$", token)
    if not match:
        return None
    axis = match.group(1)
    value = resolve_size(match.group(2))
    if value is None:
        return None
    if axis is None:
        return f"gap:{value}"
    if axis == "x":
        return f"column-gap:{value}"
    return f"row-gap:{value}"


def resolve_space_y_rule(class_name: str) -> str | None:
    if not class_name.startswith("space-y-"):
        return None
    value = resolve_size(class_name[len("space-y-") :])
    if value is None:
        return None
    escaped = escape_class_name(class_name)
    return f".{escaped}>:not([hidden])~:not([hidden]){{margin-top:{value}}}"


def resolve_z_utility(token: str) -> str | None:
    if token.startswith("z-"):
        raw = token[2:]
        if raw.startswith("[") and raw.endswith("]"):
            return f"z-index:{raw[1:-1]}"
        if re.fullmatch(r"-?\d+", raw):
            return f"z-index:{raw}"
    return None


def resolve_transform_utility(token: str) -> str | None:
    translate = re.match(r"^(-?)translate-([xy])-(.+)$", token)
    if translate:
        negative = translate.group(1) == "-"
        axis = translate.group(2)
        raw = translate.group(3)
        value = resolve_size(raw)
        if value is None:
            return None
        if negative:
            value = negate_value(value)
        if axis == "x":
            return (
                f"--tw-translate-x:{value};"
                "transform:translate(var(--tw-translate-x,0),var(--tw-translate-y,0)) "
                "scaleX(var(--tw-scale-x,1)) scaleY(var(--tw-scale-y,1))"
            )
        return (
            f"--tw-translate-y:{value};"
            "transform:translate(var(--tw-translate-x,0),var(--tw-translate-y,0)) "
            "scaleX(var(--tw-scale-x,1)) scaleY(var(--tw-scale-y,1))"
        )

    scale = re.match(r"^scale-(\d+)$", token)
    if scale:
        value = int(scale.group(1)) / 100.0
        return (
            f"--tw-scale-x:{value:g};--tw-scale-y:{value:g};"
            "transform:translate(var(--tw-translate-x,0),var(--tw-translate-y,0)) "
            "scaleX(var(--tw-scale-x,1)) scaleY(var(--tw-scale-y,1))"
        )

    return None


def resolve_ring_utility(token: str) -> str | None:
    if token.startswith("ring-"):
        raw = token[len("ring-") :]
        if raw == "0":
            return SIMPLE["ring-0"]
        if raw.isdigit():
            return f"box-shadow:0 0 0 {int(raw)}px var(--tw-ring-color,rgba(59,130,246,0.5))"
        color = resolve_color(raw)
        if color:
            return f"--tw-ring-color:{color}"
    return None


def resolve_shadow_utility(token: str) -> str | None:
    if token in SHADOWS:
        return SHADOWS[token]
    if token.startswith("shadow-"):
        color = resolve_color(token[len("shadow-") :])
        if color:
            return f"--tw-shadow-color:{color}"
    return None


def resolve_base_utility(token: str) -> str | None:
    if token in DISPLAY:
        return DISPLAY[token]
    if token in POSITION:
        return POSITION[token]
    if token in SIMPLE:
        return SIMPLE[token]
    if token in TRANSITIONS:
        return TRANSITIONS[token]
    if token in ROUNDED:
        return ROUNDED[token]

    if token.startswith("duration-"):
        raw = token[len("duration-") :]
        if raw.isdigit():
            return f"transition-duration:{int(raw) / 1000:g}s"

    if token.startswith("text-"):
        resolved = resolve_text(token)
        if resolved:
            return resolved

    if token.startswith("bg-"):
        resolved = resolve_background(token)
        if resolved:
            return resolved

    if token.startswith("font-"):
        raw = token[len("font-") :]
        if raw in FONT_WEIGHTS:
            return f"font-weight:{FONT_WEIGHTS[raw]}"

    if token.startswith("leading-"):
        raw = token[len("leading-") :]
        if raw in LEADING:
            return f"line-height:{LEADING[raw]}"

    if token.startswith("tracking-"):
        raw = token[len("tracking-") :]
        if raw in TRACKING:
            return f"letter-spacing:{TRACKING[raw]}"

    if token.startswith("border"):
        resolved = resolve_border(token)
        if resolved:
            return resolved

    resolved = resolve_spacing_utility(token)
    if resolved:
        return resolved

    resolved = resolve_size_utility(token)
    if resolved:
        return resolved

    resolved = resolve_position_utility(token)
    if resolved:
        return resolved

    resolved = resolve_grid_utility(token)
    if resolved:
        return resolved

    resolved = resolve_gap_utility(token)
    if resolved:
        return resolved

    resolved = resolve_z_utility(token)
    if resolved:
        return resolved

    resolved = resolve_transform_utility(token)
    if resolved:
        return resolved

    resolved = resolve_ring_utility(token)
    if resolved:
        return resolved

    resolved = resolve_shadow_utility(token)
    if resolved:
        return resolved

    if token.startswith("opacity-"):
        raw = token[len("opacity-") :]
        opacity = parse_opacity(raw)
        if opacity is not None:
            return f"opacity:{opacity:g}"

    return None


def split_variants(token: str) -> tuple[list[str], str]:
    variants: list[str] = []
    rest = token
    while ":" in rest:
        prefix, remainder = rest.split(":", 1)
        if prefix not in VARIANTS:
            break
        variants.append(prefix)
        rest = remainder
    return variants, rest


def escape_class_name(name: str) -> str:
    escaped = name
    for ch in ["\\", ":", "/", ".", "%", "[", "]", "(", ")", "#", "!"]:
        escaped = escaped.replace(ch, f"\\{ch}")
    if escaped.startswith("-"):
        escaped = "\\-" + escaped[1:]
    return escaped


def apply_variants_to_class(class_name: str, variants: list[str]) -> list[str]:
    selectors = [f".{escape_class_name(class_name)}"]
    for variant in variants:
        if variant == "group-hover":
            selectors = [f".group:hover {selector}" for selector in selectors]
            continue
        if variant == "hover":
            selectors = [f"{selector}:hover" for selector in selectors]
            continue
        if variant == "focus":
            selectors = [f"{selector}:focus" for selector in selectors]
            continue
        if variant == "active":
            selectors = [f"{selector}:active" for selector in selectors]
            continue
        if variant == "disabled":
            selectors = [
                variant_selector
                for selector in selectors
                for variant_selector in (f"{selector}:disabled", f"{selector}[disabled]")
            ]
            continue
    return selectors


def apply_variants_to_selector(selector: str, variants: list[str]) -> list[str]:
    selectors = [selector]
    for variant in variants:
        if variant == "group-hover":
            selectors = [f".group:hover {sel}" for sel in selectors]
            continue
        if variant == "hover":
            selectors = [f"{sel}:hover" for sel in selectors]
            continue
        if variant == "focus":
            selectors = [f"{sel}:focus" for sel in selectors]
            continue
        if variant == "active":
            selectors = [f"{sel}:active" for sel in selectors]
            continue
        if variant == "disabled":
            selectors = [
                disabled_sel
                for sel in selectors
                for disabled_sel in (f"{sel}:disabled", f"{sel}[disabled]")
            ]
            continue
    return selectors


def compile_utility_class(token: str) -> str | None:
    variants, base = split_variants(token)
    if base.startswith("space-y-") and not variants:
        return resolve_space_y_rule(base)

    declarations = resolve_base_utility(base)
    if declarations is None:
        return None

    selectors = apply_variants_to_class(token, variants)
    return f"{','.join(selectors)}{{{declarations}}}"


def expand_apply_rule(selector: str, apply_tokens: list[str]) -> tuple[list[str], list[str]]:
    base_declarations: list[str] = []
    variant_rules: list[str] = []

    for token in apply_tokens:
        variants, base = split_variants(token)
        declaration = resolve_base_utility(base)
        if declaration is None:
            raise BuildError(f"Unsupported utility in @apply: {token}")

        if not variants:
            base_declarations.append(declaration)
            continue

        selectors = apply_variants_to_selector(selector, variants)
        variant_rules.append(f"{','.join(selectors)}{{{declaration}}}")

    return base_declarations, variant_rules


def transform_source_css(css_text: str) -> str:
    lines: list[str] = []
    rule_pattern = re.compile(r"^(?P<indent>\s*)(?P<selector>[^{}]+?)\s*\{(?P<body>[^{}]*)\}\s*$")

    for raw_line in css_text.splitlines():
        if re.match(r"^\s*@tailwind\s+\w+\s*;\s*$", raw_line):
            continue

        match = rule_pattern.match(raw_line)
        if not match or "@apply" not in raw_line:
            lines.append(raw_line)
            continue

        selector = match.group("selector").strip()
        body = match.group("body")

        apply_parts = re.findall(r"@apply\s+([^;]+);", body)
        tokens: list[str] = []
        for part in apply_parts:
            tokens.extend(part.split())

        body_no_apply = re.sub(r"@apply\s+[^;]+;", "", body)
        existing_decls = [decl.strip() for decl in body_no_apply.split(";") if decl.strip()]

        apply_decls, variant_rules = expand_apply_rule(selector, tokens)
        merged = existing_decls + apply_decls
        if merged:
            lines.append(f"{selector}{{{';'.join(merged)}}}")
        for variant_rule in variant_rules:
            lines.append(variant_rule)

    return "\n".join(lines) + "\n"


def compile_utility_layer(content_files: list[Path], safelist: list[str]) -> str:
    class_tokens: set[str] = set()

    for path in content_files:
        text = path.read_text(encoding="utf-8", errors="ignore")
        class_tokens.update(extract_class_tokens_from_content(text))

    for entry in safelist:
        for token in entry.split():
            token = normalize_token(token)
            if is_class_token(token):
                class_tokens.add(token)

    rules: list[str] = []
    unknown: list[str] = []

    for token in sorted(class_tokens):
        rule = compile_utility_class(token)
        if rule is not None:
            rules.append(rule)
        elif looks_like_utility(token):
            unknown.append(token)

    if unknown:
        preview = "\n".join(f"- {item}" for item in sorted(set(unknown))[:40])
        raise BuildError(
            "Unsupported utility classes in Tailwind stub compiler.\n"
            "Add support in tools/tailwind/bin/tailwindcss_stub.py:\n"
            f"{preview}"
        )

    keyframes = []
    if "animate-spin" in class_tokens:
        keyframes.append("@keyframes tw-spin{to{transform:rotate(360deg)}}")
    if "animate-pulse" in class_tokens:
        keyframes.append("@keyframes tw-pulse{0%,100%{opacity:1}50%{opacity:.5}}")

    if not rules and not keyframes:
        return ""

    return "@layer utilities{" + "".join(keyframes + rules) + "}"


def minify_css(css_text: str) -> str:
    no_comments = re.sub(r"/\*.*?\*/", "", css_text, flags=re.S)
    squeezed = re.sub(r"\s+", " ", no_comments).strip()
    squeezed = re.sub(r"\s*([{}:;,>~])\s*", r"\1", squeezed)
    squeezed = squeezed.replace(";}", "}")
    return squeezed + "\n"


def main() -> int:
    args = parse_args()
    input_path = Path(args.input_css)
    output_path = Path(args.output_css)
    config_path = Path(args.config) if args.config else None

    source_css = input_path.read_text(encoding="utf-8")
    transformed = transform_source_css(source_css)

    content_files, safelist = load_tailwind_config(config_path)
    utility_layer = compile_utility_layer(content_files, safelist)
    final_css = transformed + utility_layer + "\n"

    if args.minify:
        final_css = minify_css(final_css)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(final_css, encoding="utf-8")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except BuildError as exc:
        raise SystemExit(str(exc))
