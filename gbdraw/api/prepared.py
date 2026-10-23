"""CLI-independent values shared after request input normalization."""

from __future__ import annotations

import dataclasses
import itertools
import secrets
import sys

from collections import OrderedDict
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Callable, Hashable, Iterator, Mapping, MutableMapping

from Bio.SeqRecord import SeqRecord
from pandas import DataFrame

from gbdraw.core.fonts import registered_fonts
from gbdraw.exceptions import ValidationError
from gbdraw.features.source import SourceFeatureIdentity


@dataclass(frozen=True)
class ParsedRecordInputs:
    """One cached source parse and its immutable pre-transform feature catalogs."""

    records: tuple[SeqRecord, ...]
    source_feature_catalogs: tuple[tuple[SourceFeatureIdentity, ...], ...]


@dataclass(frozen=True)
class ResolvedFeatureInputs:
    """Feature tables and compiled rules resolved once for one diagram request."""

    color_table: DataFrame | None
    default_colors: DataFrame
    feature_visibility_table: DataFrame | None
    feature_visibility_rules: list[dict[str, Any]] | None
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
    # The full stamp walks every feature; the structure stamp stops at the
    # feature list (D-05). Both are taken when the entry is built.
    owner_stamp: Hashable
    structure_stamp: Hashable
    retained_bytes: int
    # Render stages only: the keys of the stages its build called, kept when
    # this entry is reused (design 2.3, M10).
    inner_keys: frozenset[Hashable] = frozenset()
    # Render stages only: the build number in its cache (step 3a). A hit
    # keeps it; a rebuild gets a new one, never reused within a cache.
    build_id: int = 0


# The fifth layer holds render-stage results (gbdraw.render_stages). Its stamp
# is a constant: the layer decides that, not the shape of the value.
_RENDER_STAGE_LAYER = "render_stage"
_RENDER_STAGE_STAMP = ("render-stage",)
_LAYERS = ("decoded", "parsed", "resolved", _RENDER_STAGE_LAYER)
# The input layers: their values are functions of the validated resource bytes
# alone, so a write transaction keeps them also when its render fails (OV-1575).
_INPUT_LAYERS = ("decoded", "parsed")


class _NoIdentity:
    """Marks a value that must not be keyed: an unstored render-stage result."""

    __slots__ = ()

    def __repr__(self) -> str:
        return "NO_IDENTITY"


NO_IDENTITY = _NoIdentity()


@dataclass
class _PreparedCacheTransaction:
    cache: "PreparedBiologicalInputCache"
    resource_paths: Mapping[Path, PreparedResourceIdentity]
    active_identities: frozenset[PreparedResourceIdentity]
    diagnostics: MutableMapping[str, Any] | None
    pending: dict[str, OrderedDict[Hashable, _PreparedCacheEntry]] = field(
        default_factory=lambda: {layer: OrderedDict() for layer in _LAYERS}
    )
    used: dict[str, set[Hashable]] = field(
        default_factory=lambda: {layer: set() for layer in _LAYERS}
    )
    accessed: dict[str, dict[Hashable, Any]] = field(
        default_factory=lambda: {layer: {} for layer in _LAYERS}
    )
    # One table maps a value to its identity for this transaction: a decoded
    # resource to its PreparedResourceIdentity, a resolved record to its
    # per-record identity, a render-stage result or a direct field of one to
    # the stage key (NO_IDENTITY when the result was not stored).
    value_identities: dict[int, tuple[object, Hashable]] = field(default_factory=dict)
    # The fonts registered when the transaction starts (D-11): text measured
    # inside a render stage depends on them, so every stage key includes them.
    font_environment: tuple[tuple[str, int, bool, str, int], ...] = ()
    # O-2 (D-07): a layer built an entry that holds records, or a render stage
    # that takes a record built (or ran uncacheable). Only then does the commit
    # check every used owner in full.
    record_touched: bool = False
    # A helper's lookup (``transaction(read_only=True)``): entries are read,
    # never published, pruned, evicted or reordered; misses build for the call.
    read_only: bool = False
    # The protein extraction's lookup also publishes the input layers it built
    # for the resources it staged, as a Generate does (OV-1595).
    publish_inputs: bool = False


# The first element of a resolved record's identity (``request_render``).
RECORD_IDENTITY_TAG = "record-v1"

_ACTIVE_PREPARED_TRANSACTION: ContextVar[_PreparedCacheTransaction | None] = (
    ContextVar("gbdraw_active_prepared_transaction", default=None)
)

_CACHE_METRICS = (
    "decodedResourceCacheHitCount",
    "decodedResourceCacheMissCount",
    "decodedResourceBuildCount",
    "parsedSourceCacheHitCount",
    "parsedSourceCacheMissCount",
    "parsedSourceParseCount",
    "resolvedRecordCacheHitCount",
    "resolvedRecordCacheMissCount",
    "resolvedRecordBuildCount",
    "interactiveFeatureTraversalCount",
    "preparedInputCacheEvictionCount",
    "preparedInputCacheRetainedBytes",
    "preparedInputCacheMutationViolationCount",
    "preparedInputCacheFullCheckCount",
    "renderStageHitCount",
    "renderStageBuildCount",
    "renderStageUncacheableCount",
    "renderStageUndeclaredReadCount",
)

_LAYER_METRICS = {
    "decoded": (
        "decodedResourceCacheHitCount",
        "decodedResourceCacheMissCount",
        "decodedResourceBuildCount",
    ),
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


def register_prepared_value_identity(value: object, identity: Hashable) -> None:
    """Record the identity of one value for the active transaction.

    The first registration of a value wins, so a decoded resource keeps its
    resource identity when a render stage returns it unchanged.
    """

    transaction = _ACTIVE_PREPARED_TRANSACTION.get()
    if transaction is None or value is None:
        return
    transaction.value_identities.setdefault(id(value), (value, identity))


def prepared_value_identity(value: object) -> Hashable | None:
    """Return the identity registered for ``value`` in the active transaction."""

    transaction = _ACTIVE_PREPARED_TRANSACTION.get()
    if transaction is None or value is None:
        return None
    registered = transaction.value_identities.get(id(value))
    if registered is None or registered[0] is not value:
        return None
    return registered[1]


def register_prepared_resource_value(value: Any, source: str | Path) -> Any:
    """Associate one decoded value with the resource that supplied it."""

    identity = prepared_resource_identity(source)
    if identity is not None:
        register_prepared_value_identity(value, identity)
        # The fields of a decoded dataclass (an orthogroup or collinearity
        # result) come from the same resource, so a stage keys them by it
        # instead of walking them.
        if dataclasses.is_dataclass(value) and not isinstance(value, type):
            for item in dataclasses.fields(value):
                child = getattr(value, item.name)
                if child is None or isinstance(child, (bool, int, float, str, bytes, Enum)) or (
                    type(child) in (tuple, frozenset) and not child
                ):
                    continue
                register_prepared_value_identity(child, ("field", identity, item.name))
    return value


def prepared_resource_value_identity(value: Any) -> PreparedResourceIdentity | None:
    """Return the compact resource identity registered for a decoded value."""

    identity = prepared_value_identity(value)
    return identity if isinstance(identity, PreparedResourceIdentity) else None


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


def _record_stamp(record: SeqRecord, *, walk_features: bool) -> Hashable:
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
        (
            tuple(_feature_stamp(feature) for feature in features)
            if walk_features
            else None
        ),
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


def _owner_stamp(value: Any, *, walk_features: bool) -> Hashable:
    """Stamp a cached value; ``walk_features=False`` gives the structure stamp.

    The structure stamp keeps every record field but the per-feature walk:
    record object, id, name, sequence and feature-list identity and length,
    and the annotations one level deep. Values without records have one stamp.
    """

    records = _records_from_value(value)
    if records is not None:
        provenance = getattr(value, "provenance", None)
        return (
            "records",
            tuple(
                _record_stamp(record, walk_features=walk_features)
                for record in records
            ),
            _shallow_value_stamp(provenance),
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


# The smallest str object: each object cell of a frame points to one beside its characters.
_STR_OBJECT_BYTES = sys.getsizeof("")


def _retained_bytes(value: object, source_bytes: int) -> int:
    """Estimate the bytes ``value`` holds without walking its rows.

    ``source_bytes`` is the size of the resources it was built from: the text
    in a decoded frame's cells, and the estimate of a decoded value of another
    kind (an orthogroup or collinearity result).
    """

    if isinstance(value, DataFrame):
        object_cells = len(value) * sum(1 for dtype in value.dtypes if dtype.kind == "O")
        return (
            int(value.memory_usage(index=True, deep=False).sum())
            + object_cells * _STR_OBJECT_BYTES
            + source_bytes
        )
    if isinstance(value, str):
        return len(value)
    records = _records_from_value(value)
    if isinstance(value, tuple) and records is None:
        # A nucleotide edge holds its text beside its frame.
        return sum(_retained_bytes(item, source_bytes) for item in value)
    if records is not None:
        total = 0
        for record in records:
            features = getattr(record, "features", ())
            total += len(record.seq) + len(features) * 512
            total += sum(
                len(getattr(feature, "qualifiers", None) or {}) * 64
                for feature in features
            )
        catalogs = getattr(value, "source_feature_catalogs", ())
        if not catalogs:
            catalogs = tuple(getattr(item, "source_feature_catalog", ()) or ()
                             for item in getattr(value, "provenance", ()))
        for catalog in catalogs:
            total += sum(256 + len(entry.location_parts) * 96
                         + sum(len(key) + sum(len(text) for text in values)
                               for key, values in entry.qualifiers)
                         for entry in catalog)
        return total
    return source_bytes


class PreparedBiologicalInputCache:
    """One-project Worker cache for decoded, parsed and resolved inputs and render stages.

    The Pyodide Worker owns one instance. Keys contain resource ID, cache token,
    byte size, and the semantic fields for the corresponding preparation layer.
    Values are decoded resources, parsed source records, one resolved record
    collection, and the render-stage results used by the current request,
    including batch items.
    Resolved records and render-stage results are published only when SVG
    generation and catalog construction succeed. Decoded resources and parsed
    sources are functions of the resource bytes, so a write transaction that
    fails publishes the ones it built, unless it changed them; it prunes
    nothing. A build that failed leaves no entry. The Generate after a failed
    one (the Label fit choice, for example) parses no source again.

    A successful transaction is the one-project bound: it discards entries that
    the request did not use. A changed semantic key rebuilds its dependent layer;
    a changed or removed resource identity evicts every entry that refers to it.
    The retained-byte diagnostic estimates what each entry holds without
    serializing it: sequence lengths, feature and qualifier counts of records;
    the columns, string cells and text of a comparison frame; the text of a
    nucleotide edge; the source size of another decoded value. Render-stage
    results are not counted. The total is kept as entries are published and
    removed, so no transaction walks the entries it does not use. Worker
    termination releases the cache with the Python runtime.

    Ownership guard (D-05): a reused entry is checked by structure when it is
    read and in full, every feature, when the transaction commits, so a change
    to a cached record fails the same request. The commit checks in full only
    when the request built an entry that holds records or ran a render stage
    that takes a record (O-2, D-07); otherwise a change is found by the next
    request that does. After a transaction that did not commit, the next one
    checks in full when it reads, too.

    A Worker helper that plans the request the Worker staged (alignment
    projection, Feature Edits TSV, stored protein frames) reads the cache in a
    read-only transaction: a lookup hits under the same keys, a miss builds a
    value that serves only that call, and nothing is published, pruned,
    evicted or reordered, so the render-stage entries stay for the next
    Generate. The owners it read are checked when it ends, as a commit checks
    them. The blastp protein extraction's lookup (``publish_inputs``) also
    publishes the decoded resources and parsed sources it built, which the
    page staged under the ids the next Generate names (OV-1595); it prunes
    and evicts nothing.
    """

    def __init__(self) -> None:
        self._layers: dict[str, OrderedDict[Hashable, _PreparedCacheEntry]] = {
            layer: OrderedDict() for layer in _LAYERS
        }
        # A transaction that did not commit may have changed what it read, so
        # the next one checks owners in full at entry, not only at commit.
        self._previous_transaction_failed = False
        # A render-stage build is named by this cache's random epoch and its
        # build number (gbdraw/web/CLAUDE.md, Runtime data flow). Without a
        # random source there is no epoch and no build is named (D-38 OQ-5).
        self._build_ids = itertools.count(1)
        self._epoch = _new_cache_epoch()
        # The retained bytes of every entry, kept as entries are published
        # and removed (OV-1583).
        self._retained_total = 0
        # Every resource identity an entry names (a superset after an
        # eviction): a manifest that holds them all removes no entry.
        self._held_identities: frozenset[PreparedResourceIdentity] = frozenset()

    def _evict(
        self,
        layer: str,
        key: Hashable,
        transaction: _PreparedCacheTransaction,
    ) -> None:
        # A rejected value this request built is not kept when it fails either.
        transaction.pending[layer].pop(key, None)
        entry = self._layers[layer].pop(key, None)
        if entry is not None:
            self._retained_total -= entry.retained_bytes
            _record_metric(transaction, "preparedInputCacheEvictionCount")

    def _synchronize_resources(
        self,
        transaction: _PreparedCacheTransaction,
    ) -> None:
        if self._held_identities <= transaction.active_identities:
            return
        held: set[PreparedResourceIdentity] = set()
        for layer, entries in self._layers.items():
            for key, entry in tuple(entries.items()):
                if entry.resource_identities <= transaction.active_identities:
                    held.update(entry.resource_identities)
                else:
                    self._evict(layer, key, transaction)
        self._held_identities = frozenset(held)

    def _set_retained_metric(
        self,
        transaction: _PreparedCacheTransaction,
    ) -> None:
        if transaction.diagnostics is not None:
            metrics = transaction.diagnostics.setdefault("metrics", {})
            metrics["preparedInputCacheRetainedBytes"] = self._retained_total

    @contextmanager
    def transaction(
        self,
        *,
        resource_paths: Mapping[str | Path, PreparedResourceIdentity],
        diagnostics: MutableMapping[str, Any] | None,
        read_only: bool = False,
        publish_inputs: bool = False,
    ) -> Iterator[None]:
        """Activate one transactional Web render against the active manifest.

        ``read_only`` opens a helper's lookup instead; with it, ``publish_inputs``
        lets that lookup publish the input layers it built (see the class
        docstring). A write transaction publishes them anyway.
        """

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
            read_only=read_only,
            publish_inputs=publish_inputs,
            font_environment=tuple(
                sorted(
                    (font.family, font.weight, font.italic, font.face.content_id, font.face.index)
                    for font in registered_fonts()
                )
            ),
        )
        if diagnostics is not None:
            metrics = diagnostics.setdefault("metrics", {})
            for name in _CACHE_METRICS:
                metrics.setdefault(name, 0)
        if not read_only:
            self._synchronize_resources(transaction)
        self._set_retained_metric(transaction)
        token = _ACTIVE_PREPARED_TRANSACTION.set(transaction)
        try:
            yield
            if read_only:
                self._verify_used_owners(transaction)
                if publish_inputs:
                    self._publish_inputs(transaction, verified=True)
            else:
                self._commit(transaction)
        except BaseException:
            if not read_only or publish_inputs:
                self._publish_inputs(transaction)
            self._previous_transaction_failed = True
            raise
        else:
            # A lookup checked only the owners it read.
            if not read_only:
                self._previous_transaction_failed = False
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
        # D-05: the structure here, every feature at commit; in full here too
        # after a transaction that did not commit. A render-stage entry has no
        # owner to check (one constant stamp).
        if layer != _RENDER_STAGE_LAYER:
            self._verify_owner(
                layer,
                key,
                entry,
                transaction,
                walk_features=self._previous_transaction_failed,
            )
        if not transaction.read_only:
            self._layers[layer].move_to_end(key)
        return entry

    def _verify_owner(
        self,
        layer: str,
        key: Hashable,
        entry: _PreparedCacheEntry,
        transaction: _PreparedCacheTransaction,
        *,
        walk_features: bool,
    ) -> None:
        expected = entry.owner_stamp if walk_features else entry.structure_stamp
        if walk_features and _records_from_value(entry.value) is not None:
            _record_metric(transaction, "preparedInputCacheFullCheckCount")
        current_stamp = _owner_stamp(entry.value, walk_features=walk_features)
        if current_stamp == expected:
            return
        _record_metric(
            transaction,
            "preparedInputCacheMutationViolationCount",
        )
        self._evict(layer, key, transaction)
        raise ValidationError(
            f"Cached {layer} biological input ownership was mutated "
            f"({_owner_mutation_detail(expected, current_stamp)})."
        )

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
        if _records_from_value(value) is not None:
            transaction.record_touched = True
        kept = not transaction.read_only or (transaction.publish_inputs and layer in _INPUT_LAYERS)
        if not kept or (publish is not None and not publish(value)):
            return value
        entry = _PreparedCacheEntry(
            value=value,
            resource_identities=resource_identities,
            owner_stamp=_owner_stamp(value, walk_features=True),
            structure_stamp=_owner_stamp(value, walk_features=False),
            retained_bytes=_retained_bytes(
                value, sum(identity.size for identity in resource_identities)
            ),
        )
        transaction.pending[layer][key] = entry
        transaction.used[layer].add(key)
        return value

    def lookup_render_stage(
        self,
        key: Hashable,
        transaction: _PreparedCacheTransaction,
    ) -> _PreparedCacheEntry | None:
        """Return a stored render-stage entry and mark it and its inner stages used."""

        entry = self._entry(_RENDER_STAGE_LAYER, key, transaction)
        if entry is None:
            return None
        used = transaction.used[_RENDER_STAGE_LAYER]
        used.add(key)
        used.update(entry.inner_keys)
        return entry

    def store_render_stage(
        self,
        key: Hashable,
        value: object,
        inner_keys: frozenset[Hashable],
        transaction: _PreparedCacheTransaction,
    ) -> None:
        """Stage one render-stage result; it is published if the render succeeds."""

        # A helper's lookup publishes nothing, so it does not walk the key either.
        if transaction.read_only:
            return
        transaction.pending[_RENDER_STAGE_LAYER][key] = _PreparedCacheEntry(
            value=value,
            resource_identities=_key_resource_identities(key),
            owner_stamp=_RENDER_STAGE_STAMP,
            structure_stamp=_RENDER_STAGE_STAMP,
            retained_bytes=0,
            inner_keys=inner_keys,
            build_id=next(self._build_ids),
        )
        transaction.used[_RENDER_STAGE_LAYER].add(key)

    def render_stage_build(
        self,
        key: Hashable,
        transaction: _PreparedCacheTransaction,
    ) -> tuple[str, int] | None:
        """Name the build stored under ``key``: this cache's epoch and its build number.

        Read without verifying or moving the entry: the stage that returned it
        already did.
        """

        if self._epoch is None:
            return None
        entry = transaction.pending[_RENDER_STAGE_LAYER].get(key)
        if entry is None:
            entry = self._layers[_RENDER_STAGE_LAYER].get(key)
        if entry is None:
            return None
        return (self._epoch, entry.build_id)

    def _used_entries(
        self,
        layer: str,
        transaction: _PreparedCacheTransaction,
    ) -> Iterator[tuple[Hashable, _PreparedCacheEntry]]:
        current = self._layers[layer]
        pending = transaction.pending[layer]
        for key in transaction.used[layer]:
            entry = pending.get(key)
            if entry is None:
                entry = current.get(key)
            if entry is not None:
                yield key, entry

    def _verify_used_owners(self, transaction: _PreparedCacheTransaction) -> None:
        # O-2: without a record built or read by a stage, the structure check
        # (made when each entry was read) is repeated, not the feature walk.
        walk_features = transaction.record_touched
        for layer in self._layers:
            # A render-stage entry has no owner to check (one constant stamp).
            if layer == _RENDER_STAGE_LAYER:
                continue
            for key, entry in tuple(self._used_entries(layer, transaction)):
                self._verify_owner(
                    layer,
                    key,
                    entry,
                    transaction,
                    walk_features=walk_features,
                )

    def _publish_inputs(
        self, transaction: _PreparedCacheTransaction, *, verified: bool = False,
    ) -> None:
        """Publish the decoded and parsed inputs a transaction built.

        A failed write transaction, and the protein extraction's lookup, call
        it. An entry whose value the request changed is not kept (checked in
        full, as a commit checks; ``verified`` when the lookup's end has just
        checked it). Nothing is pruned: neither is the one-project bound.
        """

        for layer in _INPUT_LAYERS:
            current = self._layers[layer]
            for key, entry in transaction.pending[layer].items():
                if not verified:
                    try:
                        self._verify_owner(layer, key, entry, transaction, walk_features=True)
                    except Exception:
                        # Changed, or its stamp cannot be built: the request's own
                        # error stays the one it raises.
                        continue
                previous = current.pop(key, None)
                self._retained_total += entry.retained_bytes - (
                    previous.retained_bytes if previous is not None else 0
                )
                current[key] = entry
                self._held_identities |= entry.resource_identities

    def _commit(self, transaction: _PreparedCacheTransaction) -> None:
        self._verify_used_owners(transaction)
        replacements: dict[str, tuple[OrderedDict[Hashable, _PreparedCacheEntry], int]] = {}
        retained_total = 0
        held: set[PreparedResourceIdentity] = set()
        for layer in self._layers:
            current = self._layers[layer]
            retained: OrderedDict[Hashable, _PreparedCacheEntry] = OrderedDict()
            kept = 0
            for key, entry in self._used_entries(layer, transaction):
                retained[key] = entry
                retained_total += entry.retained_bytes
                held.update(entry.resource_identities)
                # A pending key is a new entry, not one kept from the layer.
                if key in current:
                    kept += 1
            replacements[layer] = (retained, kept)
        for layer, (replacement, kept) in replacements.items():
            evicted = len(self._layers[layer]) - kept
            if evicted:
                _record_metric(
                    transaction,
                    "preparedInputCacheEvictionCount",
                    evicted,
                )
            self._layers[layer] = replacement
        self._retained_total = retained_total
        self._held_identities = frozenset(held)
        self._set_retained_metric(transaction)


def _new_cache_epoch() -> str | None:
    """A random epoch per cache, or ``None`` without a random source (never a time or id seed)."""

    try:
        return secrets.token_hex(8)
    except (NotImplementedError, OSError):
        return None


def _key_resource_identities(key: Hashable) -> frozenset[PreparedResourceIdentity]:
    """Collect the resource identities a key names, so their removal evicts it."""

    identities: set[PreparedResourceIdentity] = set()
    pending: list[object] = [key]
    while pending:
        item = pending.pop()
        if isinstance(item, PreparedResourceIdentity):
            identities.add(item)
        elif isinstance(item, (tuple, list, frozenset)):
            pending.extend(item)
    return frozenset(identities)


def key_names_record(key: Hashable) -> bool:
    """Whether a render-stage key takes a resolved record itself.

    A stage result a key takes (another stage's key, or a field of one) does
    not count: that stage read the record, this one reads its result.
    """

    pending: list[object] = [key]
    while pending:
        item = pending.pop()
        if not isinstance(item, tuple):
            continue
        if len(item) == 2 and item[0] == "identity":
            identity = item[1]
            if isinstance(identity, tuple) and identity and identity[0] == RECORD_IDENTITY_TAG:
                return True
            continue
        pending.extend(item)
    return False


def mark_record_touched(transaction: _PreparedCacheTransaction) -> None:
    """A render stage that may read a record built: the commit checks owners in full (O-2)."""

    transaction.record_touched = True


def active_render_stage_transaction() -> _PreparedCacheTransaction | None:
    """Return the active transaction for ``gbdraw.render_stages``, if any."""

    return _ACTIVE_PREPARED_TRANSACTION.get()


def lookup_render_stage(
    transaction: _PreparedCacheTransaction,
    key: Hashable,
) -> _PreparedCacheEntry | None:
    """Return the render-stage entry stored under ``key``, if any."""

    return transaction.cache.lookup_render_stage(key, transaction)


def render_stage_build(
    transaction: _PreparedCacheTransaction,
    key: Hashable,
) -> tuple[str, int] | None:
    """Return the cache epoch and build number of the render stage stored under ``key``."""

    return transaction.cache.render_stage_build(key, transaction)


def store_render_stage(
    transaction: _PreparedCacheTransaction,
    key: Hashable,
    value: object,
    inner_keys: frozenset[Hashable],
) -> None:
    """Stage one render-stage result for publication at commit."""

    transaction.cache.store_render_stage(key, value, inner_keys, transaction)


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


def get_or_build_decoded_resource(
    source: str | Path,
    key: Hashable,
    builder: Callable[[], Any],
) -> Any:
    """Reuse one decoded immutable value while its Worker resource is unchanged."""

    transaction = _ACTIVE_PREPARED_TRANSACTION.get()
    identity = prepared_resource_identity(source)
    if transaction is None or identity is None:
        return builder()
    value = transaction.cache.get_or_build(
        "decoded",
        (identity, key),
        frozenset({identity}),
        builder,
    )
    return register_prepared_resource_value(value, source)


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
    return value


__all__ = ["ResolvedFeatureInputs"]
