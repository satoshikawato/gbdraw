"""Migrations that bring a saved Session's Web-owned fields to the current writer.

The CLI re-save and :func:`gbdraw.session.upgrade_session_document` share
them: the canonical request is migrated by the typed bridge, and these
functions migrate the fields around it.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Mapping, NamedTuple, Sequence, cast

from gbdraw.exceptions import ValidationError
from gbdraw.session_io import (
    CURRENT_AUTHORITY_SESSION_MIN_VERSION,
    MODE_SCOPED_SESSION_MIN_VERSION,
    RETIRED_RENDERED_ID_FEATURE_FIELDS,
    SessionDraftMigration,
    _project_web_file_binding,
    empty_protein_identity_manifest,
    migrate_legacy_linear_comparison_draft_for_current_writer,
    migrate_session_flat_draft,
    mode_split_palette_colors,
    session_depth_source_widths,
    session_mode,
    split_draft_into_modes,
    unmapped_hash_selector_warning,
)

logger = logging.getLogger(__name__)

# The Web reader names a Session it cannot migrate as a Session field error.
_SESSION_INVALID = {"code": "INPUT_INVALID", "field": "schema", "reason": "FIELDS"}


def replace_current_derived_feature_state(
    payload: dict[str, Any],
    feature_catalog: Mapping[str, object] | None,
) -> None:
    """Install the Results' catalog and drop the payloads the writer derives."""

    # The catalog stays at the top level, with the committed set; the draft is
    # in the mode slices (Session 46).
    editor_state = payload.get("editorState")
    editor_state = (
        dict(editor_state) if isinstance(editor_state, Mapping) else {}
    )
    editor_state["featureCatalog"] = (
        dict(feature_catalog) if feature_catalog is not None else None
    )
    payload["editorState"] = editor_state

    orthogroup_state = payload.get("orthogroupState")
    orthogroup_state = (
        dict(orthogroup_state)
        if isinstance(orthogroup_state, Mapping)
        else {}
    )
    orthogroup_state.pop("groups", None)
    payload["orthogroupState"] = orthogroup_state


def _project_web_file_inventory(
    session: Mapping[str, Any],
) -> dict[str, Any] | None:
    web_files = session.get("webFiles")
    resources = session.get("resources")
    if not isinstance(web_files, Mapping) or not isinstance(resources, Mapping):
        return None
    bindings_value = web_files.get("bindings")
    has_current_bindings = (
        isinstance(bindings_value, Mapping) and bindings_value.get("schema") in (1, 2)
    )
    # has_current_bindings is true only when bindings_value is a Mapping.
    bindings = cast(
        "Mapping[str, Any]", bindings_value if has_current_bindings else {}
    )
    direct_source_fields = {
        "conservationLosatFastaSources": "c_conservation_fastas",
        "conservationSequenceSources": "c_conservation_sequence_sources",
    }
    has_direct_sources = any(
        isinstance(web_files.get(field), list) for field in direct_source_fields
    )
    if not has_current_bindings and not has_direct_sources:
        return None

    original_names_value = web_files.get("resourceOriginalNames")
    original_names = (
        original_names_value if isinstance(original_names_value, Mapping) else {}
    )

    def restore(value: Any) -> Any:
        return _project_web_file_binding(resources, value, schema=bindings["schema"])

    def restore_resource_id(value: Any) -> Any:
        if isinstance(value, list):
            return [restore_resource_id(item) for item in value]
        resource_id = str(value or "").strip()
        if not resource_id:
            return None
        return _project_web_file_binding(
            resources,
            {
                "resourceId": resource_id,
                "name": original_names.get(resource_id),
            },
            schema=None,
        )

    files: dict[str, Any] = {}
    for slot in (
        "c_gb",
        "c_gff",
        "c_fasta",
        "c_depth",
        "c_conservation_blasts",
        "c_conservation_fastas",
        "c_conservation_sequence_sources",
        "d_color",
        "t_color",
        "blacklist",
        "whitelist",
        "qualifier_priority",
    ):
        if slot in bindings:
            files[slot] = restore(bindings[slot])
    files["c_conservation_blasts_source"] = (
        "losat-cache"
        if bindings.get("c_conservation_blasts_source") == "losat-cache"
        else None
    )

    linear_sequences = bindings.get("linearSeqs")
    if isinstance(linear_sequences, list):
        files["linearSeqs"] = [
            {
                **dict(sequence),
                "gb": restore(sequence.get("gb")),
                "gff": restore(sequence.get("gff")),
                "fasta": restore(sequence.get("fasta")),
                "depth": restore(sequence.get("depth")),
                "blast": restore(sequence.get("blast")),
            }
            for sequence in linear_sequences
            if isinstance(sequence, Mapping)
        ]
    linear_comparisons = bindings.get("linearComparisons")
    if isinstance(linear_comparisons, list):
        files["linearComparisons"] = [
            {**dict(comparison), "file": restore(comparison.get("file"))}
            for comparison in linear_comparisons
            if isinstance(comparison, Mapping)
        ]
    for source_field, slot in direct_source_fields.items():
        source_ids = web_files.get(source_field)
        if isinstance(source_ids, list) and slot not in files:
            files[slot] = restore_resource_id(source_ids)
    return files


# Similarity-alignment flags a projected source session drops from its argv:
# the current flag and the spellings that legacy sessions carry.
_SIMILARITY_ALIGNMENT_FLAGS = frozenset(
    {
        "--similarity_alignment_feature",
        "--align_orthogroup_feature",
        "--align-orthogroup-feature",
    }
)


class LegacySourceRead(NamedTuple):
    """One source read of a Session without a feature catalog, as the Web reads it."""

    resource_id: str
    region_spec: str | None
    record_selector: str | None
    reverse: bool


def legacy_source_reads(request: Mapping[str, Any]) -> tuple[LegacySourceRead, ...] | None:
    """The source reads that name a catalog-less Session's features again.

    Web Load reads each GenBank source with the crop and orientation of its
    request record (``extractSessionSourceFeatures``): a Circular source whole,
    a Linear one with ``<start>-<end>[:rc]`` and no reverse flag when it is
    cropped, otherwise with the record's reverse flag, and with its record
    selector (``recordId`` value or ``#<index + 1>``). None when a source is
    not GenBank or a crop has one endpoint, where the Web reads nothing, or a
    record index is not an integer, which the request decode rejects.
    """

    records = request.get("records")
    if not isinstance(records, list) or not records:
        return None
    linear = request.get("mode") == "linear"
    reads: list[LegacySourceRead] = []
    for record in records if linear else records[:1]:
        record = record if isinstance(record, Mapping) else {}
        source = record.get("source")
        if not isinstance(source, Mapping) or source.get("kind") != "genbank":
            return None
        if not linear:
            reads.append(LegacySourceRead(str(source.get("resourceId")), None, None, False))
            continue
        region_value = record.get("region")
        region: Mapping[str, Any] = region_value if isinstance(region_value, Mapping) else {}
        presentation = record.get("presentation")
        presentation = presentation if isinstance(presentation, Mapping) else {}
        selector = region.get("selector") or record.get("selector")
        selector = selector if isinstance(selector, Mapping) else {}
        record_selector = ""
        if selector.get("kind") == "recordId":
            record_selector = str(selector.get("value") or "")
        elif selector.get("kind") == "recordIndex":
            index = selector.get("index")
            if isinstance(index, bool) or not isinstance(index, int):
                return None  # Invalid; the request decode names it.
            record_selector = f"#{index + 1}"
        reverse = bool(region.get("reverseComplement") or presentation.get("reverseComplement"))
        start, end = region.get("start"), region.get("end")
        if (start is None) != (end is None):
            return None
        reads.append(
            LegacySourceRead(
                str(source.get("resourceId")),
                f"{start}-{end}{':rc' if reverse else ''}" if start is not None else None,
                record_selector or None,
                reverse and start is None,
            )
        )
    return tuple(reads)


def read_legacy_source_features(
    session: Mapping[str, Any], resource_paths: Mapping[str, Path]
) -> dict[str, list[dict[str, Any]]] | None:
    """Read a catalog-less Session's sources again for its rendered-ID edits.

    The twin of Web Load's source read (``extractSessionSourceFeatures``) with
    the reads of :func:`legacy_source_reads`: ``extractedFeatures`` and
    ``biologicalFeatures`` of every source, each named by its source hash and,
    for Linear, its input (``fileIdx``). None when the Session saved a catalog
    or no rendered-ID edit, or a source cannot be read; the edit migration then
    uses the saved feature metadata. The Web's last fallback, the features of
    the saved SVG (its recovery plan), is not read.
    """

    from gbdraw.web_support.feature_metadata import extract_features_from_genbank_payload

    features = session.get("features")
    editor_state = session.get("editorState")
    if not isinstance(features, Mapping) or not any(
        isinstance(features.get(field), Mapping) and features[field]
        for field in RETIRED_RENDERED_ID_FEATURE_FIELDS
    ):
        return None
    if isinstance(editor_state, Mapping) and isinstance(editor_state.get("featureCatalog"), Mapping):
        return None
    request = session.get("renderRequest")
    reads = legacy_source_reads(request) if isinstance(request, Mapping) else None
    if reads is None:
        return None
    linear = cast(Mapping[str, Any], request).get("mode") == "linear"
    read_again: dict[str, list[dict[str, Any]]] = {"extractedFeatures": [], "biologicalFeatures": []}
    for input_index, read in enumerate(reads):
        path = resource_paths.get(read.resource_id)
        if path is None:
            return None
        try:
            payload = extract_features_from_genbank_payload(
                path,
                read.region_spec,
                read.record_selector,
                read.reverse,
                include_biological_features=True,
            )
        except Exception as exc:  # The Web reads no source that fails.
            logger.debug("Session source %s could not be read again: %s", read.resource_id, exc)
            return None
        listed = payload.get("features")
        biological = payload.get("biological_features")
        if not isinstance(listed, list):
            return None
        if not isinstance(biological, list) or not (biological or not listed):
            biological = listed
        for key, entries in (("extractedFeatures", listed), ("biologicalFeatures", biological)):
            for entry in entries:
                feature = dict(entry)
                stable_id = str(
                    next(
                        (
                            feature[field]
                            for field in ("stable_feature_id", "stableFeatureId", "stable_svg_id", "svg_id")
                            if feature.get(field)
                        ),
                        "",
                    )
                ).strip()
                if linear:
                    feature["fileIdx"] = input_index
                feature.update(svg_id=stable_id, stable_svg_id=stable_id, stable_feature_id=stable_id)
                read_again[key].append(feature)
    return read_again


def log_session_draft_migration(migration: SessionDraftMigration, source_version: int) -> None:
    """Report what the draft migrations of a Session 27-44 changed or could not move."""

    if migration.dropped_feature_edit_count:
        logger.warning(
            "WARNING: %d feature edit(s) from Session version %d could not "
            "be matched to a feature of its saved diagram and were dropped "
            "from the written Session.",
            migration.dropped_feature_edit_count,
            source_version,
        )
    if migration.narrowed_visibility_count:
        logger.warning(
            "WARNING: %d Feature visibility edit(s) from Session version %d "
            "hid every feature with the same hash; in the written Session "
            "each applies only to the feature that was edited.",
            migration.narrowed_visibility_count,
            source_version,
        )
    if migration.unmapped_hash_selector_count:
        logger.warning(
            "WARNING: %s",
            unmapped_hash_selector_warning(
                migration.unmapped_hash_selector_count, source_version
            ),
        )
    if migration.migrated_annotation_count:
        logger.info(
            "INFO: %d annotation(s) from Session version %d named a feature "
            "by hash=; in the written Session each names that feature by "
            "its source.",
            migration.migrated_annotation_count,
            source_version,
        )


def project_session_adjunct_for_current_write(
    session: Mapping[str, Any],
    *,
    source_version: int,
    source_features: Mapping[str, Any] | None = None,
) -> tuple[dict[str, Any], dict[str, Any] | None]:
    """Detach non-canonical state and migrate released Web-owned field names.

    ``source_features`` are the Session's sources read again
    (:func:`read_legacy_source_features`) for the rendered-ID edits of a
    Session without a feature catalog.
    """

    adjunct = {
        key: value
        for key, value in session.items()
        if key
        not in {
            "format",
            "version",
            "createdAt",
            "renderRequest",
            "resources",
            "files",
        }
    }
    if source_version < MODE_SCOPED_SESSION_MIN_VERSION:
        # The older draft migrations, in Web Load's order: the rendered-ID edit
        # maps become identity drafts through the Session's saved catalog (or
        # its sources read again, or its saved feature metadata), and
        # a hash= annotation target moves to its source feature where that
        # catalog makes the figure certain (R-7). The request keeps the targets
        # that drew the figure.
        migration = migrate_session_flat_draft(session, source_features=source_features)
        for key in ("config", "features"):
            if key in migration.session:
                adjunct[key] = migration.session[key]
            else:
                adjunct.pop(key, None)
        log_session_draft_migration(migration, source_version)
    orthogroup_state = adjunct.get("orthogroupState")
    if isinstance(orthogroup_state, Mapping):
        projected_orthogroup_state = dict(orthogroup_state)
        projected_orthogroup_state.pop(
            "selectedOrthogroupAlignmentFeature",
            None,
        )
        adjunct["orthogroupState"] = projected_orthogroup_state
    cli_invocation = adjunct.get("cliInvocation")
    if isinstance(cli_invocation, Mapping):
        projected_invocation = dict(cli_invocation)
        args = cli_invocation.get("args")
        bindings = cli_invocation.get("fileBindings")
        if isinstance(args, list):
            projected_args: list[str] = []
            retained_indexes: dict[int, int] = {}
            index = 0
            while index < len(args):
                token = str(args[index])
                if token in _SIMILARITY_ALIGNMENT_FLAGS:
                    index += 2
                    continue
                if token.startswith(
                    tuple(f"{flag}=" for flag in _SIMILARITY_ALIGNMENT_FLAGS)
                ):
                    index += 1
                    continue
                retained_indexes[index] = len(projected_args)
                projected_args.append(token)
                index += 1
            projected_invocation["args"] = projected_args
            if isinstance(bindings, list):
                projected_bindings = []
                for binding in bindings:
                    if not isinstance(binding, Mapping):
                        projected_bindings.append(binding)
                        continue
                    arg_index = binding.get("argIndex")
                    if arg_index not in retained_indexes:
                        raise ValidationError(
                            "Legacy similarity alignment cannot own a CLI file binding.",
                            diagnostic=_SESSION_INVALID,
                        )
                    projected_bindings.append(
                        {
                            **dict(binding),
                            "argIndex": retained_indexes[arg_index],
                        }
                    )
                projected_invocation["fileBindings"] = projected_bindings
        adjunct["cliInvocation"] = projected_invocation
    editor_state_value = adjunct.get("editorState")
    if isinstance(editor_state_value, Mapping):
        editor_state = dict(editor_state_value)
        catalog = editor_state.get("featureCatalog")
        if isinstance(catalog, Mapping) and catalog.get("schema") == 3:
            from gbdraw.web_support.feature_catalog import (
                promote_legacy_feature_catalog,
            )

            editor_state["featureCatalog"] = promote_legacy_feature_catalog(catalog)
            adjunct["editorState"] = editor_state
    web_file_inventory = _project_web_file_inventory(session)
    if source_version >= CURRENT_AUTHORITY_SESSION_MIN_VERSION:
        return _split_session_adjunct(adjunct, session, source_version), web_file_inventory
    config = adjunct.get("config")

    if isinstance(config, Mapping):
        source_files = session.get("files")
        has_source_file_inventory = (
            isinstance(source_files, Mapping) and bool(source_files)
        ) or web_file_inventory is not None
        migrated_config, migrated_files = (
            migrate_legacy_linear_comparison_draft_for_current_writer(
                config,
                source_files
                if isinstance(source_files, Mapping)
                else (web_file_inventory or {}),
                force_web_draft=(
                    isinstance(config.get("linearRecordLayout"), Mapping)
                    or not isinstance(config.get("cliOptions"), Mapping)
                ),
            )
        )
        adjunct["config"] = migrated_config
        web_file_inventory = migrated_files if has_source_file_inventory else None
    else:
        adjunct.pop("config", None)
    if isinstance(adjunct.get("ui"), Mapping):
        adjunct["ui"] = dict(adjunct["ui"])
        adjunct["ui"].pop("blastSource", None)
    else:
        adjunct.pop("ui", None)
    web_files_value = adjunct.get("webFiles")
    if isinstance(web_files_value, Mapping):
        web_files = dict(web_files_value)
        bindings_value = web_files.get("bindings")
        if isinstance(bindings_value, Mapping):
            bindings = dict(bindings_value)
            bindings.pop("linearCanonicalComparisons", None)
            if web_file_inventory is not None:
                bindings = {}
            else:
                linear_sequences = bindings.get("linearSeqs")
                if isinstance(linear_sequences, list):
                    bindings["linearSeqs"] = [
                        {
                            key: value
                            for key, value in sequence.items()
                            if key not in {"blast", "losat_filename"}
                        }
                        if isinstance(sequence, Mapping)
                        else sequence
                        for sequence in linear_sequences
                    ]
                bindings["linearComparisons"] = []
            web_files["bindings"] = bindings
        metadata_value = web_files.get("linearRecordMetadata")
        if isinstance(metadata_value, list):
            web_files["linearRecordMetadata"] = [
                {
                    key: value
                    for key, value in metadata.items()
                    if key not in {"losatFilename", "losat_filename"}
                }
                if isinstance(metadata, Mapping)
                else metadata
                for metadata in metadata_value
            ]
        adjunct["webFiles"] = web_files
    return _split_session_adjunct(adjunct, session, source_version), web_file_inventory


def _split_session_adjunct(
    adjunct: dict[str, Any],
    session: Mapping[str, Any],
    source_version: int,
) -> dict[str, Any]:
    """Move an older source's flat draft into the Session 46 mode slices.

    A Session 46 source keeps its ``modes`` (and ``otherModeResult``) as they
    are. Each mode's Depth sources are counted in the source's bindings.
    """

    if source_version >= MODE_SCOPED_SESSION_MIN_VERSION:
        return adjunct
    web_files = session.get("webFiles")
    bindings = web_files.get("bindings") if isinstance(web_files, Mapping) else None
    return split_draft_into_modes(
        adjunct,
        committed_mode=session_mode(session),
        depth_sources=session_depth_source_widths(
            bindings if isinstance(bindings, Mapping) else session.get("files")
        ),
        palette_colors=mode_split_palette_colors(adjunct.get("config")),
    )


def with_current_artifacts(
    fields: Mapping[str, Any],
    *,
    losat_cache_entries: Sequence[Mapping[str, Any]],
    protein_identity_manifest: Mapping[str, Any] | None,
    legacy_protein_raw_candidates: Sequence[Mapping[str, Any]] = (),
    legacy_protein_derived_evidence: Sequence[Mapping[str, Any]] = (),
    protein_id_map: Mapping[str, str] | None = None,
) -> dict[str, Any]:
    """``fields`` with the LOSAT artifacts of a fresh render or migration.

    The legacy protein IDs that the migration resolved are rewritten first.
    The derived cache is written empty, as every current writer writes it.
    """

    updated = dict(fields)
    if protein_id_map:
        from gbdraw.api.session_compat import rewrite_protein_artifact_references

        updated = rewrite_protein_artifact_references(updated, protein_id_map)
    for artifact_key in (
        "losatCache",
        "losatDerivedCache",
        "proteinIdentityManifest",
        "legacyArtifacts",
    ):
        updated.pop(artifact_key, None)
    updated["losatCache"] = {"entries": [dict(entry) for entry in losat_cache_entries]}
    updated["losatDerivedCache"] = {"entries": []}
    updated["proteinIdentityManifest"] = dict(
        protein_identity_manifest or empty_protein_identity_manifest()
    )
    legacy_artifacts: dict[str, Any] = {}
    if legacy_protein_raw_candidates:
        legacy_artifacts["proteinRawCandidates"] = {
            "schema": 1,
            "entries": [dict(entry) for entry in legacy_protein_raw_candidates],
        }
    if legacy_protein_derived_evidence:
        legacy_artifacts["proteinDerivedEvidence"] = {
            "schema": 1,
            "entries": [dict(entry) for entry in legacy_protein_derived_evidence],
        }
    if legacy_artifacts:
        updated["legacyArtifacts"] = legacy_artifacts
    return updated


__all__ = [
    "LegacySourceRead",
    "legacy_source_reads",
    "log_session_draft_migration",
    "project_session_adjunct_for_current_write",
    "read_legacy_source_features",
    "replace_current_derived_feature_state",
    "with_current_artifacts",
]
