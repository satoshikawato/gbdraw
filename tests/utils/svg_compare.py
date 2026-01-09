"""
SVG comparison utilities for regression testing.

This module provides functions to compare two SVG files for semantic equivalence,
handling common sources of non-significant differences such as:
- Floating point representation differences
- Attribute order variations
- Whitespace normalization
- XML declaration differences
"""

import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Union
import math


@dataclass
class SVGComparisonResult:
    """Result of comparing two SVG files."""

    equal: bool
    message: str
    differences: list[str]

    def __bool__(self) -> bool:
        return self.equal


def normalize_number(value: str, tolerance: float = 1e-4) -> str:
    """
    Normalize a numeric string for comparison.

    Rounds floating point numbers to avoid precision differences.

    Args:
        value: String that may contain a number
        tolerance: Relative tolerance for float comparison

    Returns:
        Normalized string representation
    """
    try:
        num = float(value)
        if math.isnan(num) or math.isinf(num):
            return value
        # Round to 4 decimal places for comparison
        if abs(num) < 1e-10:
            return "0"
        rounded = round(num, 4)
        # Format consistently
        if rounded == int(rounded):
            return str(int(rounded))
        return f"{rounded:.4f}".rstrip('0').rstrip('.')
    except ValueError:
        return value


def normalize_path_data(path_data: str) -> str:
    """
    Normalize SVG path data for comparison.

    Handles floating point precision and whitespace variations in path commands.

    Args:
        path_data: SVG path d attribute value

    Returns:
        Normalized path data string
    """
    # Split into commands and coordinates
    # Pattern matches path commands followed by numbers
    pattern = r'([MmZzLlHhVvCcSsQqTtAa])|(-?[\d.]+(?:[eE][+-]?\d+)?)'
    tokens = re.findall(pattern, path_data)

    result = []
    for cmd, num in tokens:
        if cmd:
            result.append(cmd)
        elif num:
            result.append(normalize_number(num))

    return ' '.join(result)


def normalize_transform(transform: str) -> str:
    """
    Normalize SVG transform attribute for comparison.

    Args:
        transform: SVG transform attribute value

    Returns:
        Normalized transform string
    """
    # Pattern to match transform functions
    pattern = r'(translate|rotate|scale|matrix|skewX|skewY)\s*\(([^)]+)\)'

    def normalize_func(match):
        func_name = match.group(1)
        params = match.group(2)
        # Normalize numbers in parameters
        normalized_params = ','.join(
            normalize_number(p.strip())
            for p in params.split(',')
        )
        return f"{func_name}({normalized_params})"

    return re.sub(pattern, normalize_func, transform)


def normalize_style(style: str) -> str:
    """
    Normalize SVG style attribute for comparison.

    Sorts style properties and normalizes values.

    Args:
        style: CSS style string

    Returns:
        Normalized style string
    """
    if not style:
        return ""

    # Parse style properties
    props = {}
    for prop in style.split(';'):
        prop = prop.strip()
        if ':' in prop:
            key, value = prop.split(':', 1)
            key = key.strip().lower()
            value = value.strip()

            # Normalize numeric values in style
            if re.match(r'^-?[\d.]+(?:px|pt|em|%)?$', value):
                num_match = re.match(r'^(-?[\d.]+)', value)
                if num_match:
                    unit = value[len(num_match.group(1)):]
                    value = normalize_number(num_match.group(1)) + unit

            props[key] = value

    # Sort and reconstruct
    return ';'.join(f"{k}:{v}" for k, v in sorted(props.items()))


def normalize_attribute_value(name: str, value: str, id_map: dict = None) -> str:
    """
    Normalize an SVG attribute value based on attribute type.

    Args:
        name: Attribute name
        value: Attribute value
        id_map: Optional dictionary for normalizing dynamic IDs

    Returns:
        Normalized value
    """
    name_lower = name.lower()

    if name_lower == 'd':
        return normalize_path_data(value)
    elif name_lower == 'transform':
        return normalize_transform(value)
    elif name_lower == 'style':
        return normalize_style(value)
    elif name_lower in ('x', 'y', 'x1', 'y1', 'x2', 'y2', 'cx', 'cy',
                         'r', 'rx', 'ry', 'width', 'height',
                         'stroke-width', 'font-size', 'opacity'):
        return normalize_number(value)
    elif name_lower == 'points':
        # Polygon/polyline points
        nums = re.findall(r'-?[\d.]+(?:[eE][+-]?\d+)?', value)
        return ' '.join(normalize_number(n) for n in nums)
    elif name_lower == 'id' and id_map is not None:
        # Normalize dynamic IDs to sequential IDs
        if value not in id_map:
            id_map[value] = f"id_{len(id_map)}"
        return id_map[value]
    elif name_lower == 'clip-path' and id_map is not None:
        # Normalize clip-path references
        match = re.match(r'url\(#(.+)\)', value)
        if match:
            ref_id = match.group(1)
            if ref_id not in id_map:
                id_map[ref_id] = f"id_{len(id_map)}"
            return f"url(#{id_map[ref_id]})"
        return value
    elif name_lower in ('xlink:href', 'href') and value.startswith('#') and id_map is not None:
        # Normalize href references
        ref_id = value[1:]
        if ref_id not in id_map:
            id_map[ref_id] = f"id_{len(id_map)}"
        return f"#{id_map[ref_id]}"
    else:
        return value.strip()


def element_to_canonical(elem: ET.Element, depth: int = 0, id_map: dict = None) -> list[str]:
    """
    Convert an element to a canonical string representation for comparison.

    Args:
        elem: XML element
        depth: Current nesting depth
        id_map: Dictionary for normalizing dynamic IDs

    Returns:
        List of canonical string lines
    """
    if id_map is None:
        id_map = {}

    lines = []
    indent = "  " * depth

    # Tag name
    tag = elem.tag
    # Remove namespace if present
    if '}' in tag:
        tag = tag.split('}')[1]

    # Normalized attributes (sorted by name)
    attrs = []
    for name, value in sorted(elem.attrib.items()):
        # Remove namespace from attribute names
        if '}' in name:
            name = name.split('}')[1]
        norm_value = normalize_attribute_value(name, value, id_map)
        attrs.append(f'{name}="{norm_value}"')

    attr_str = ' '.join(attrs)

    # Element text
    text = (elem.text or '').strip()

    # Build element representation
    if attr_str:
        lines.append(f"{indent}<{tag} {attr_str}>")
    else:
        lines.append(f"{indent}<{tag}>")

    if text:
        lines.append(f"{indent}  {text}")

    # Recursively process children
    for child in elem:
        lines.extend(element_to_canonical(child, depth + 1, id_map))

    lines.append(f"{indent}</{tag}>")

    return lines


def parse_svg(content: str) -> ET.Element:
    """
    Parse SVG content, handling namespace issues.

    Args:
        content: SVG file content as string

    Returns:
        Parsed XML element tree root
    """
    # Remove XML declaration to avoid parsing issues
    content = re.sub(r'<\?xml[^?]*\?>', '', content)

    # Register SVG namespace
    namespaces = {
        'svg': 'http://www.w3.org/2000/svg',
        'xlink': 'http://www.w3.org/1999/xlink',
    }
    for prefix, uri in namespaces.items():
        ET.register_namespace(prefix, uri)

    return ET.fromstring(content)


def compare_svgs(
    expected: Union[str, Path],
    actual: Union[str, Path],
    tolerance: float = 1e-4,
    ignore_comments: bool = True,
    max_differences: int = 10
) -> SVGComparisonResult:
    """
    Compare two SVG files for semantic equivalence.

    This function compares SVGs by parsing them and comparing their structure
    and attributes, normalizing floating point values and handling common
    sources of non-significant differences.

    Args:
        expected: Path to expected SVG file or SVG content string
        actual: Path to actual SVG file or SVG content string
        tolerance: Floating point comparison tolerance
        ignore_comments: Whether to ignore XML comments
        max_differences: Maximum number of differences to report

    Returns:
        SVGComparisonResult with comparison outcome and any differences
    """
    # Load content
    if isinstance(expected, Path) or (isinstance(expected, str) and not expected.strip().startswith('<')):
        expected_path = Path(expected)
        if not expected_path.exists():
            return SVGComparisonResult(
                equal=False,
                message=f"Expected file not found: {expected_path}",
                differences=[]
            )
        expected_content = expected_path.read_text(encoding='utf-8')
    else:
        expected_content = expected

    if isinstance(actual, Path) or (isinstance(actual, str) and not actual.strip().startswith('<')):
        actual_path = Path(actual)
        if not actual_path.exists():
            return SVGComparisonResult(
                equal=False,
                message=f"Actual file not found: {actual_path}",
                differences=[]
            )
        actual_content = actual_path.read_text(encoding='utf-8')
    else:
        actual_content = actual

    # Parse SVGs
    try:
        expected_root = parse_svg(expected_content)
    except ET.ParseError as e:
        return SVGComparisonResult(
            equal=False,
            message=f"Failed to parse expected SVG: {e}",
            differences=[]
        )

    try:
        actual_root = parse_svg(actual_content)
    except ET.ParseError as e:
        return SVGComparisonResult(
            equal=False,
            message=f"Failed to parse actual SVG: {e}",
            differences=[]
        )

    # Convert to canonical form
    expected_lines = element_to_canonical(expected_root)
    actual_lines = element_to_canonical(actual_root)

    # Compare line by line
    differences = []
    max_len = max(len(expected_lines), len(actual_lines))

    for i in range(max_len):
        if i >= len(expected_lines):
            differences.append(f"Line {i+1}: Extra in actual: {actual_lines[i][:100]}")
        elif i >= len(actual_lines):
            differences.append(f"Line {i+1}: Missing in actual: {expected_lines[i][:100]}")
        elif expected_lines[i] != actual_lines[i]:
            differences.append(
                f"Line {i+1}:\n"
                f"  Expected: {expected_lines[i][:100]}\n"
                f"  Actual:   {actual_lines[i][:100]}"
            )

        if len(differences) >= max_differences:
            differences.append(f"... (truncated, {max_len - i - 1} more lines to compare)")
            break

    if differences:
        return SVGComparisonResult(
            equal=False,
            message=f"Found {len(differences)} difference(s)",
            differences=differences
        )

    return SVGComparisonResult(
        equal=True,
        message="SVGs are semantically equivalent",
        differences=[]
    )


def compare_svg_files(expected_path: Path, actual_path: Path) -> SVGComparisonResult:
    """
    Convenience function to compare two SVG files.

    Args:
        expected_path: Path to expected SVG file
        actual_path: Path to actual SVG file

    Returns:
        SVGComparisonResult with comparison outcome
    """
    return compare_svgs(expected_path, actual_path)


def quick_hash_compare(file1: Path, file2: Path) -> bool:
    """
    Quick comparison using file hash.

    If hashes match, files are identical. If not, need deeper comparison.

    Args:
        file1: First file path
        file2: Second file path

    Returns:
        True if files are byte-identical
    """
    import hashlib

    def file_hash(path: Path) -> str:
        h = hashlib.sha256()
        h.update(path.read_bytes())
        return h.hexdigest()

    return file_hash(file1) == file_hash(file2)
