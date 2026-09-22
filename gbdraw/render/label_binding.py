"""Internal SVG metadata for complete feature-label visual units."""

from svgwrite.params import Parameter  # type: ignore[reportMissingImports]


LABEL_FEATURE_ID_ATTRIBUTE = "data-label-feature-id"
LABEL_BINDING_SCHEMA_ATTRIBUTE = "data-gbdraw-label-binding-schema"
LABEL_BINDING_SCHEMA = "1"


def bind_label_part(element, feature_id: str, *, complete: bool = False):
    """Bind one label text or leader segment to its rendered feature."""

    identity = str(feature_id or "").strip()
    if not identity:
        raise ValueError("Rendered label parts require a feature identity.")
    element.set_parameter(Parameter(debug=False, profile=element.profile))
    element.attribs[LABEL_FEATURE_ID_ATTRIBUTE] = identity
    if complete:
        element.attribs[LABEL_BINDING_SCHEMA_ATTRIBUTE] = LABEL_BINDING_SCHEMA
    return element


__all__ = [
    "LABEL_BINDING_SCHEMA",
    "LABEL_BINDING_SCHEMA_ATTRIBUTE",
    "LABEL_FEATURE_ID_ATTRIBUTE",
    "bind_label_part",
]
