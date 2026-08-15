"""CLI-independent values shared after request input normalization."""

from __future__ import annotations

from collections import OrderedDict
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Callable, Hashable, Iterator, Mapping, MutableMapping

from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]
from pandas import DataFrame  # type: ignore[reportMissingImports]

from gbdraw.exceptions import ValidationError
from gbdraw.features.colors import preprocess_color_tables
from gbdraw.features.visibility import compile_feature_visibility_rules


@dataclass(frozen=True)
class ResolvedFeatureInputs:
    """Feature tables and compiled rules resolved once for one diagram request."""

    color_table: DataFrame | None
    default_colors: DataFrame
    feature_visibility_table: DataFrame | None
    feature_visibility_rules: list[dict[str, Any]]
    specific_color_rules: Mapping[str, Any]
    default_color_map: Mapping[str, str]


@dataclass(frozen=True, order=True)
class PreparedResourceIdentity:
    """Compact identity for one Worker-owned canonical resource."""

    resource_id: str
    cache_token: str
    size: int


@dataclass(frozen=True)
class _PreparedCacheEntry:
    value: Any
    resource_identities: frozenset[PreparedResourceIdentity]
    owner_stamp: Hashable
    retained_bytes: int


@dataclass
class _PreparedCacheTransaction:
    cache: "PreparedBiologicalInputCache"
    resource_paths: Mapping[Path, PreparedResourceIdentity]
    active_identities: frozenset[PreparedResourceIdentity]
    diagnostics: MutableMapping[str, Any] | None
    pending: dict[str, OrderedDict[Hashable, _PreparedCacheEntry]] = field(
        default_factory=lambda: {
            "parsed": OrderedDict(),
            "resolved": OrderedDict(),
            "interactive": OrderedDict(),
        }
    )
    used: dict[str, set[Hashable]] = field(
        default_factory=lambda: {
            "parsed": set(),
            "resolved": set(),
            "interactive": set(),
        }
    )
    accessed: dict[str, dict[Hashable, Any]] = field(
        default_factory=lambda: {
            "parsed": {},
            "resolved": {},
            "interactive": {},
        }
    )
    decoded_resource_values: dict[int, tuple[Any, PreparedResourceIdentity]] = field(
        default_factory=dict
    )
    resolved_record_membership: dict[int, tuple[Hashable, int]] = field(
        default_factory=dict
    )


_ACTIVE_PREPARED_TRANSACTION: ContextVar[_PreparedCacheTransaction | None] = (
    ContextVar("gbdraw_active_prepared_transaction", default=None)
)

_CACHE_METRICS = (
    "parsedSourceCacheHitCount",
    "parsedSourceCacheMissCount",
    "parsedSourceParseCount",
    "resolvedRecordCacheHitCount",
    "resolvedRecordCacheMissCount",
    "resolvedRecordBuildCount",
    "interactiveContextCacheHitCount",
    "interactiveContextCacheMissCount",
    "interactiveContextBuildCount",
    "interactiveFeatureTraversalCount",
    "selectorSafetyScopeBuildCount",
    "preparedInputCacheEvictionCount",
    "preparedInputCacheRetainedBytes",
    "preparedInputCacheMutationViolationCount",
)

_LAYER_METRICS = {
    "parsed": (
        "parsedSourceCacheHitCount",
        "parsedSourceCacheMissCount",
        "parsedSourceParseCount",
    ),
    "resolved": (
        "resolvedRecordCacheHitCount",
        "resolvedRecordCacheMissCount",
        "resolvedRecordBuildCount",
    ),
    "interactive": (
        "interactiveContextCacheHitCount",
        "interactiveContextCacheMissCount",
        "interactiveContextBuildCount",
    ),
}


def _record_metric(
    transaction: _PreparedCacheTransaction | None,
    name: str,
    value: int = 1,
) -> None:
    if transaction is None or transaction.diagnostics is None:
        return
    metrics = transaction.diagnostics.setdefault("metrics", {})
    metrics[name] = int(metrics.get(name, 0)) + int(value)


def record_prepared_input_metric(name: str, value: int = 1) -> None:
    """Record a private Web preparation metric when its cache is active."""

    _record_metric(_ACTIVE_PREPARED_TRANSACTION.get(), name, value)


def prepared_input_cache_active() -> bool:
    """Return whether the current call is inside a Worker cache transaction."""

    return _ACTIVE_PREPARED_TRANSACTION.get() is not None


def prepared_resource_identity(
    value: str | Path,
) -> PreparedResourceIdentity | None:
    """Resolve an ephemeral staged path to its compact Worker identity."""

    transaction = _ACTIVE_PREPARED_TRANSACTION.get()
    if transaction is None:
        return None
    try:
        path = Path(value)
        candidates = (path.absolute(), path.resolve(strict=True))
    except (OSError, RuntimeError, ValueError):
        return None
    for candidate in candidates:
        identity = transaction.resource_paths.get(candidate)
        if identity is not None:
            return identity
    return None


def register_prepared_resource_value(value: Any, source: str | Path) -> Any:
    """Associate one decoded value with the resource that supplied it."""

    transaction = _ACTIVE_PREPARED_TRANSACTION.get()
    identity = prepared_resource_identity(source)
    if transaction is not None and identity is not None and value is not None:
        transaction.decoded_resource_values[id(value)] = (value, identity)
    return value


def prepared_resource_value_identity(value: Any) -> PreparedResourceIdentity | None:
    """Return the compact resource identity registered for a decoded value."""

    transaction = _ACTIVE_PREPARED_TRANSACTION.get()
    if transaction is None or value is None:
        return None
    registered = transaction.decoded_resource_values.get(id(value))
    if registered is None or registered[0] is not value:
        return None
    return registered[1]


def _shallow_value_stamp(value: Any) -> Hashable:
    if value is None or isinstance(value, (bool, int, float, str, bytes, Enum)):
        return value
    if isinstance(value, Mapping):
        return (
            "mapping",
            id(value),
            len(value),
            tuple(
                sorted(
                    (
                        str(key),
                        (
                            entry
                            if entry is None
                            or isinstance(entry, (bool, int, float, str, bytes, Enum))
                            else (type(entry).__name__, id(entry), len(entry))
                            if hasattr(entry, "__len__")
                            else (type(entry).__name__, id(entry))
                        ),
                    )
                    for key, entry in value.items()
                )
            ),
        )
    if isinstance(value, (list, tuple)):
        return (
            type(value).__name__,
            id(value),
            len(value),
            tuple(
                entry
                if entry is None
                or isinstance(entry, (bool, int, float, str, bytes, Enum))
                else (type(entry).__name__, id(entry), len(entry))
                if hasattr(entry, "__len__")
                else (type(entry).__name__, id(entry))
                for entry in value
            ),
        )
    return (type(value).__name__, id(value))


def _feature_stamp(feature: Any) -> Hashable:
    location = getattr(feature, "location", None)
    qualifiers = getattr(feature, "qualifiers", None)
    return (
        id(feature),
        str(getattr(feature, "type", "") or ""),
        id(location),
        int(location.start) if location is not None else None,
        int(location.end) if location is not None else None,
        getattr(location, "strand", None),
        _shallow_value_stamp(qualifiers),
    )


def _record_stamp(record: SeqRecord) -> Hashable:
    features = getattr(record, "features", ())
    annotations = getattr(record, "annotations", None)
    return (
        id(record),
        str(record.id),
        str(record.name),
        id(record.seq),
        len(record.seq),
        id(features),
        len(features),
        tuple(_feature_stamp(feature) for feature in features),
        _shallow_value_stamp(annotations),
    )


def _records_from_value(value: Any) -> tuple[SeqRecord, ...] | None:
    if isinstance(value, tuple) and value and all(
        isinstance(record, SeqRecord) for record in value
    ):
        return value
    records = getattr(value, "records", None)
    if isinstance(records, tuple) and records and all(
        isinstance(record, SeqRecord) for record in records
    ):
        return records
    return None


def _owner_stamp(value: Any) -> Hashable:
    records = _records_from_value(value)
    if records is not None:
        provenance = getattr(value, "provenance", None)
        return (
            "records",
            tuple(_record_stamp(record) for record in records),
            _shallow_value_stamp(provenance),
        )
    if all(
        hasattr(value, field_name)
        for field_name in (
            "features",
            "biological_features",
            "orthogroups",
            "annotations",
            "sequence_sources",
            "record_keys",
        )
    ):
        return (
            "interactive-context",
            tuple(
                (field_name, _shallow_value_stamp(getattr(value, field_name)))
                for field_name in (
                    "features",
                    "biological_features",
                    "orthogroups",
                    "annotations",
                    "sequence_sources",
                    "record_keys",
                    "collinearity_search_scope",
                )
            ),
        )
    return (type(value).__name__, id(value))


def _owner_mutation_detail(before: Hashable, after: Hashable) -> str:
    if not (
        isinstance(before, tuple)
        and isinstance(after, tuple)
        and before
        and after
        and before[0] == after[0] == "records"
    ):
        return "owner structure"
    before_records = before[1]
    after_records = after[1]
    if not isinstance(before_records, tuple) or not isinstance(after_records, tuple):
        return "record collection"
    if len(before_records) != len(after_records):
        return "record count"
    record_fields = (
        "owner",
        "id",
        "name",
        "sequence owner",
        "sequence length",
        "feature-list owner",
        "feature count",
        "features",
        "annotations",
    )
    for record_index, (before_record, after_record) in enumerate(
        zip(before_records, after_records, strict=True)
    ):
        if before_record == after_record:
            continue
        if not isinstance(before_record, tuple) or not isinstance(after_record, tuple):
            return f"record {record_index}"
        for field_index, (before_field, after_field) in enumerate(
            zip(before_record, after_record, strict=False)
        ):
            if before_field != after_field:
                if (
                    field_index == 7
                    and isinstance(before_field, tuple)
                    and isinstance(after_field, tuple)
                ):
                    feature_fields = (
                        "owner",
                        "type",
                        "location owner",
                        "start",
                        "end",
                        "strand",
                        "qualifiers",
                    )
                    for feature_index, (before_feature, after_feature) in enumerate(
                        zip(before_field, after_field, strict=False)
                    ):
                        if before_feature == after_feature:
                            continue
                        if isinstance(before_feature, tuple) and isinstance(
                            after_feature, tuple
                        ):
                            for feature_field_index, (
                                before_feature_field,
                                after_feature_field,
                            ) in enumerate(
                                zip(before_feature, after_feature, strict=False)
                            ):
                                if before_feature_field != after_feature_field:
                                    if (
                                        feature_field_index == 6
                                        and isinstance(before_feature_field, tuple)
                                        and isinstance(after_feature_field, tuple)
                                        and len(before_feature_field) >= 4
                                        and len(after_feature_field) >= 4
                                        and isinstance(before_feature_field[3], tuple)
                                        and isinstance(after_feature_field[3], tuple)
                                    ):
                                        before_qualifiers = dict(before_feature_field[3])
                                        after_qualifiers = dict(after_feature_field[3])
                                        changed_keys = sorted(
                                            key
                                            for key in (
                                                set(before_qualifiers)
                                                | set(after_qualifiers)
                                            )
                                            if before_qualifiers.get(key)
                                            != after_qualifiers.get(key)
                                        )
                                        if changed_keys:
                                            return (
                                                f"record {record_index} feature "
                                                f"{feature_index} qualifier owner "
                                                f"{changed_keys[0]}"
                                            )
                                    feature_field = (
                                        feature_fields[feature_field_index]
                                        if feature_field_index < len(feature_fields)
                                        else f"field {feature_field_index}"
                                    )
                                    return (
                                        f"record {record_index} feature "
                                        f"{feature_index} {feature_field}"
                                    )
                        return f"record {record_index} feature {feature_index}"
                field_name = (
                    record_fields[field_index]
                    if field_index < len(record_fields)
                    else f"field {field_index}"
                )
                return f"record {record_index} {field_name}"
    return "record provenance"


def _retained_bytes(value: Any) -> int:
    records = _records_from_value(value)
    if records is not None:
        total = 0
        for record in records:
            features = getattr(record, "features", ())
            total += len(record.seq) + len(features) * 512
            total += sum(
                len(getattr(feature, "qualifiers", None) or {}) * 64
                for feature in features
            )
        return total
    if all(
        hasattr(value, field_name)
        for field_name in ("features", "biological_features", "sequence_sources")
    ):
        total = 0
        for field_name in (
            "features",
            "biological_features",
            "orthogroups",
            "annotations",
        ):
            total += len(getattr(value, field_name, ()) or ()) * 512
        for source in getattr(value, "sequence_sources", ()) or ():
            if isinstance(source, Mapping):
                total += len(str(source.get("sequence") or "")) + 256
        return total
    return 0


class PreparedBiologicalInputCache:
    """One-project Worker cache for parsed, resolved, and interactive inputs.

    The Pyodide Worker owns one instance. Keys contain resource ID, cache token,
    byte size, and the semantic fields for the corresponding preparation layer.
    Values are parsed source records, one resolved record collection, and the
    interactive contexts used by the current request, including batch items.
    Publication waits until SVG generation and catalog construction succeed.

    A successful transaction is the one-project bound: it discards entries that
    the request did not use. A changed semantic key rebuilds its dependent layer;
    a changed or removed resource identity evicts every entry that refers to it.
    The retained-byte diagnostic estimates source bytes, sequence lengths,
    feature counts, qualifier counts, and interactive-context payloads without
    serializing them. Worker termination releases the cache with the Python
    runtime.
    """

    def __init__(self) -> None:
        self._layers: dict[str, OrderedDict[Hashable, _PreparedCacheEntry]] = {
            "parsed": OrderedDict(),
            "resolved": OrderedDict(),
            "interactive": OrderedDict(),
        }

    def _evict(
        self,
        layer: str,
        key: Hashable,
        transaction: _PreparedCacheTransaction,
    ) -> None:
        if self._layers[layer].pop(key, None) is not None:
            _record_metric(transaction, "preparedInputCacheEvictionCount")

    def _synchronize_resources(
        self,
        transaction: _PreparedCacheTransaction,
    ) -> None:
        for layer, entries in self._layers.items():
            for key, entry in tuple(entries.items()):
                if not entry.resource_identities.issubset(
                    transaction.active_identities
                ):
                    self._evict(layer, key, transaction)

    def _set_retained_metric(
        self,
        transaction: _PreparedCacheTransaction,
    ) -> None:
        retained = sum(
            entry.retained_bytes
            for entries in self._layers.values()
            for entry in entries.values()
        )
        if transaction.diagnostics is not None:
            metrics = transaction.diagnostics.setdefault("metrics", {})
            metrics["preparedInputCacheRetainedBytes"] = int(retained)

    @contextmanager
    def transaction(
        self,
        *,
        resource_paths: Mapping[str | Path, PreparedResourceIdentity],
        diagnostics: MutableMapping[str, Any] | None,
    ) -> Iterator[None]:
        """Activate one transactional Web render against the active manifest."""

        if _ACTIVE_PREPARED_TRANSACTION.get() is not None:
            raise ValidationError("Prepared input cache transactions cannot be nested.")
        normalized_paths: dict[Path, PreparedResourceIdentity] = {}
        for raw_path, identity in resource_paths.items():
            path = Path(raw_path)
            normalized_paths[path.absolute()] = identity
            normalized_paths[path.resolve(strict=True)] = identity
        transaction = _PreparedCacheTransaction(
            cache=self,
            resource_paths=normalized_paths,
            active_identities=frozenset(resource_paths.values()),
            diagnostics=diagnostics,
        )
        if diagnostics is not None:
            metrics = diagnostics.setdefault("metrics", {})
            for name in _CACHE_METRICS:
                metrics.setdefault(name, 0)
        self._synchronize_resources(transaction)
        self._set_retained_metric(transaction)
        token = _ACTIVE_PREPARED_TRANSACTION.set(transaction)
        try:
            yield
            self._commit(transaction)
        finally:
            self._set_retained_metric(transaction)
            _ACTIVE_PREPARED_TRANSACTION.reset(token)

    def _entry(
        self,
        layer: str,
        key: Hashable,
        transaction: _PreparedCacheTransaction,
    ) -> _PreparedCacheEntry | None:
        pending = transaction.pending[layer].get(key)
        if pending is not None:
            return pending
        entry = self._layers[layer].get(key)
        if entry is None:
            return None
        current_stamp = _owner_stamp(entry.value)
        if current_stamp != entry.owner_stamp:
            _record_metric(
                transaction,
                "preparedInputCacheMutationViolationCount",
            )
            self._evict(layer, key, transaction)
            raise ValidationError(
                f"Cached {layer} biological input ownership was mutated "
                f"({_owner_mutation_detail(entry.owner_stamp, current_stamp)})."
            )
        self._layers[layer].move_to_end(key)
        return entry

    def get_or_build(
        self,
        layer: str,
        key: Hashable,
        resource_identities: frozenset[PreparedResourceIdentity],
        builder: Callable[[], Any],
        *,
        publish: Callable[[Any], bool] | None = None,
    ) -> Any:
        transaction = _ACTIVE_PREPARED_TRANSACTION.get()
        if transaction is None or transaction.cache is not self:
            return builder()
        if key in transaction.accessed[layer]:
            return transaction.accessed[layer][key]
        hit_metric, miss_metric, build_metric = _LAYER_METRICS[layer]
        entry = self._entry(layer, key, transaction)
        if entry is not None:
            _record_metric(transaction, hit_metric)
            transaction.used[layer].add(key)
            transaction.accessed[layer][key] = entry.value
            return entry.value
        _record_metric(transaction, miss_metric)
        _record_metric(transaction, build_metric)
        value = builder()
        transaction.accessed[layer][key] = value
        if publish is not None and not publish(value):
            return value
        entry = _PreparedCacheEntry(
            value=value,
            resource_identities=resource_identities,
            owner_stamp=_owner_stamp(value),
            retained_bytes=(
                _retained_bytes(value)
                + sum(identity.size for identity in resource_identities)
            ),
        )
        transaction.pending[layer][key] = entry
        transaction.used[layer].add(key)
        return value

    def register_resolved_records(self, key: Hashable, value: Any) -> None:
        transaction = _ACTIVE_PREPARED_TRANSACTION.get()
        records = _records_from_value(value)
        if transaction is None or records is None:
            return
        for index, record in enumerate(records):
            transaction.resolved_record_membership[id(record)] = (key, index)

    def _commit(self, transaction: _PreparedCacheTransaction) -> None:
        replacements: dict[str, OrderedDict[Hashable, _PreparedCacheEntry]] = {}
        for layer in self._layers:
            retained: OrderedDict[Hashable, _PreparedCacheEntry] = OrderedDict()
            for key in transaction.used[layer]:
                entry = transaction.pending[layer].get(key) or self._layers[layer].get(
                    key
                )
                if entry is None:
                    continue
                current_stamp = _owner_stamp(entry.value)
                if current_stamp != entry.owner_stamp:
                    _record_metric(
                        transaction,
                        "preparedInputCacheMutationViolationCount",
                    )
                    self._evict(layer, key, transaction)
                    raise ValidationError(
                        f"Cached {layer} biological input ownership was mutated "
                        f"({_owner_mutation_detail(entry.owner_stamp, current_stamp)})."
                    )
                retained[key] = entry
            replacements[layer] = retained
        for layer, replacement in replacements.items():
            evicted = set(self._layers[layer]) - set(replacement)
            if evicted:
                _record_metric(
                    transaction,
                    "preparedInputCacheEvictionCount",
                    len(evicted),
                )
            self._layers[layer] = replacement
        self._set_retained_metric(transaction)


def get_or_build_parsed_source(
    key: Hashable,
    resource_identities: frozenset[PreparedResourceIdentity],
    builder: Callable[[], Any],
    *,
    publish: Callable[[Any], bool] | None = None,
) -> Any:
    """Get or transactionally stage one parsed biological source."""

    transaction = _ACTIVE_PREPARED_TRANSACTION.get()
    if transaction is None:
        return builder()
    return transaction.cache.get_or_build(
        "parsed",
        key,
        resource_identities,
        builder,
        publish=publish,
    )


def get_or_build_resolved_records(
    key: Hashable,
    resource_identities: frozenset[PreparedResourceIdentity],
    builder: Callable[[], Any],
) -> Any:
    """Get or transactionally stage one resolved record collection."""

    transaction = _ACTIVE_PREPARED_TRANSACTION.get()
    if transaction is None:
        return builder()
    value = transaction.cache.get_or_build(
        "resolved",
        key,
        resource_identities,
        builder,
    )
    transaction.cache.register_resolved_records(key, value)
    return value


def get_or_build_interactive_context(
    key: Hashable,
    resource_identities: frozenset[PreparedResourceIdentity],
    builder: Callable[[], Any],
) -> Any:
    """Get or transactionally stage one interactive metadata context."""

    transaction = _ACTIVE_PREPARED_TRANSACTION.get()
    if transaction is None:
        return builder()
    return transaction.cache.get_or_build(
        "interactive",
        key,
        resource_identities,
        builder,
    )


def prepared_record_membership(
    records: tuple[SeqRecord, ...],
) -> tuple[tuple[Hashable, int], ...] | None:
    """Return cached resolved-record membership without inspecting records."""

    transaction = _ACTIVE_PREPARED_TRANSACTION.get()
    if transaction is None:
        return None
    membership: list[tuple[Hashable, int]] = []
    for record in records:
        item = transaction.resolved_record_membership.get(id(record))
        if item is None:
            return None
        membership.append(item)
    return tuple(membership)


def resolve_feature_inputs(
    *,
    color_table: DataFrame | None,
    default_colors: DataFrame,
    feature_visibility_table: DataFrame | None,
) -> ResolvedFeatureInputs:
    """Compile already-loaded feature inputs into one reusable value."""

    specific_color_rules, default_color_map = preprocess_color_tables(
        color_table,
        default_colors,
    )
    return ResolvedFeatureInputs(
        color_table=color_table,
        default_colors=default_colors,
        feature_visibility_table=feature_visibility_table,
        feature_visibility_rules=compile_feature_visibility_rules(
            feature_visibility_table
        ),
        specific_color_rules=specific_color_rules,
        default_color_map=default_color_map,
    )


__all__ = ["ResolvedFeatureInputs", "resolve_feature_inputs"]
