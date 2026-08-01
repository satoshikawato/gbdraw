from __future__ import annotations

import json

from gbdraw.api import get_web_runtime_capabilities
from gbdraw.features.shapes import FEATURE_RENDERING_VALUES
from gbdraw.session_io import (
    DEPTH_FILE_ENCODING,
    LOSAT_DERIVED_CACHE_SCHEMA,
    NUCLEOTIDE_LOSAT_CACHE_SCHEMA,
    PROTEIN_IDENTITY_MANIFEST_SCHEMA,
    PROTEIN_LOSAT_CACHE_SCHEMA,
)
from gbdraw.session_request_codec import (
    CANONICAL_REQUEST_SCHEMA,
    SUPPORTED_CANONICAL_REQUEST_SCHEMAS,
    UNKNOWN_FIELD_POLICY,
)
from gbdraw.tracks import (
    SUPPORTED_CIRCULAR_TRACK_RENDERERS,
    SUPPORTED_LINEAR_TRACK_RENDERERS,
)
from gbdraw.web_support.capabilities import (
    WEB_RENDER_OPTIONS_SCHEMA,
    WEB_RENDER_PROTOCOL,
    WEB_RUNTIME_CAPABILITY_SCHEMA,
)


def test_web_runtime_capability_manifest_is_owner_backed_and_json_safe() -> None:
    capabilities = get_web_runtime_capabilities()

    assert WEB_RENDER_OPTIONS_SCHEMA == 2
    assert capabilities["rendering"]["featureRenderings"] == [
        "arrow",
        "arrowhead",
        "rectangle",
        "underlay",
    ]
    assert capabilities == {
        "schema": WEB_RUNTIME_CAPABILITY_SCHEMA,
        "renderProtocol": WEB_RENDER_PROTOCOL,
        "request": {
            "currentSchema": CANONICAL_REQUEST_SCHEMA,
            "supportedSchemas": sorted(SUPPORTED_CANONICAL_REQUEST_SCHEMAS),
            "unknownFieldPolicy": UNKNOWN_FIELD_POLICY,
        },
        "resources": {
            "encodings": ["base64", DEPTH_FILE_ENCODING],
        },
        "rendering": {
            "optionSchema": WEB_RENDER_OPTIONS_SCHEMA,
            "featureRenderings": sorted(FEATURE_RENDERING_VALUES),
            "circularTrackRenderers": sorted(
                SUPPORTED_CIRCULAR_TRACK_RENDERERS
            ),
            "linearTrackRenderers": sorted(SUPPORTED_LINEAR_TRACK_RENDERERS),
        },
        "analysis": {
            "proteinLosatCacheSchema": PROTEIN_LOSAT_CACHE_SCHEMA,
            "nucleotideLosatCacheSchema": NUCLEOTIDE_LOSAT_CACHE_SCHEMA,
            "losatDerivedCacheSchema": LOSAT_DERIVED_CACHE_SCHEMA,
            "proteinIdentityManifestSchema": PROTEIN_IDENTITY_MANIFEST_SCHEMA,
        },
    }
    assert json.loads(json.dumps(capabilities)) == capabilities


def test_web_runtime_capability_manifest_returns_fresh_state() -> None:
    first = get_web_runtime_capabilities()
    first["request"]["supportedSchemas"].append(999)

    second = get_web_runtime_capabilities()
    assert 999 not in second["request"]["supportedSchemas"]
