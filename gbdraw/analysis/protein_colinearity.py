#!/usr/bin/env python
# coding: utf-8

"""Protein-level colinearity helpers for linear pairwise comparisons."""

from __future__ import annotations

import copy
from collections import namedtuple
from dataclasses import dataclass, field, replace
from functools import lru_cache
import base64
import hashlib
from io import StringIO
import itertools
import json
import logging
import math
import re
import unicodedata
from typing import (
    AbstractSet,
    Callable,
    Iterable,
    Literal,
    Mapping,
    NamedTuple,
    NoReturn,
    Sequence,
    SupportsFloat,
    TypedDict,
    TypeGuard,
    cast,
)

import numpy as np
import pandas as pd
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from pandas import DataFrame

from gbdraw.analysis.ortholog_paths import OrthologPath, OrthologPathCollection, ortholog_edge_id
from gbdraw.contract_shapes import Spellings, read_spelling
from gbdraw.comparisons.losat_jobs import (
    LosatBatch,
    RecordPair,
    losat_entry_filename,
    losat_record_uids,
    losat_source_ids,
    losatp_job_specs,
    plan_losat_jobs,
    prepare_losat_batches,
    split_losat_batch_result,
)
from gbdraw.comparisons.losat_runtime import (
    LosatCacheWrite,
    LosatRawCache,
    LosatRuntimeCallback,
    LosatRuntimeRecord,
    LosatSearchArgs,
    losat_cache_args,
    run_losat_search,
)
from gbdraw.core.numeric import compensated_sum
from gbdraw.core.record_metadata import _source_feature_index
from gbdraw.core.sequence import translate_cds
from gbdraw.exceptions import ParseError, ValidationError
from gbdraw.features.ids import compute_feature_hash_from_location_parts, source_feature_hash_and_parts
from gbdraw.features.visibility import should_include_feature_in_analysis
from gbdraw.io.comparisons import COMPARISON_COLUMNS, read_comparison_table

# Historical internal import; edge identity is owned by ortholog_paths.
_edge_id = ortholog_edge_id

logger = logging.getLogger(__name__)

PROTEIN_IDENTITY_MANIFEST_SCHEMA = 2
PROTEIN_LOSAT_CACHE_SCHEMA = 4
LEGACY_PROTEIN_LOSAT_CACHE_SCHEMA = 2
LEGACY_PROTEIN_RAW_CANDIDATE_SCHEMA = 1
# Backwards-compatible import name. Protein raw-cache writers now emit schema 4.
LOSAT_CACHE_SCHEMA = PROTEIN_LOSAT_CACHE_SCHEMA

_VALID_PROTEIN_RE = re.compile(r"^[A-Z]+$")
_VALID_FASTA_ID_RE = re.compile(r"^\S+$")
_VALID_RUNTIME_HANDLE_RE = re.compile(r"^h_[a-z2-7]{26}$")
_RUNTIME_HANDLE_SEARCH_RE = re.compile(r"(?<![a-z2-7_])h_[a-z2-7]{26}(?![a-z2-7])")
_MANIFEST_PRESENTATION_IDENTITY_FRAGMENTS = frozenset(
    {
        "viewfeaturesvgid",
        "viewfeaturehashparts",
        "renderedfeaturesvgid",
        "renderedsvgid",
    }
)
_MANIFEST_SCALAR_TYPES = frozenset({str, int, float, bool, type(None)})
_NUMERIC_COMPARISON_COLUMNS = COMPARISON_COLUMNS[2:]
_INTEGER_COMPARISON_COLUMNS = frozenset(
    {
        "alignment_length",
        "mismatches",
        "gap_opens",
        "qstart",
        "qend",
        "sstart",
        "send",
    }
)
_STRICT_DECIMAL_INTEGER = r"[+-]?[0-9]+"
# Digits with an optional point, or a point and digits; then an optional exponent.
# "[0-9]+\.?[0-9]*" spells "[0-9]+(?:\.[0-9]*)?" without a group, which re matches faster.
_STRICT_DECIMAL_NUMBER = r"[+-]?(?:[0-9]+\.?[0-9]*|\.[0-9]+)(?:[eE][+-]?[0-9]+)?"
_STRICT_DECIMAL_INTEGER_RE = re.compile(rf"^{_STRICT_DECIMAL_INTEGER}$")
_STRICT_DECIMAL_NUMBER_RE = re.compile(rf"^{_STRICT_DECIMAL_NUMBER}$")
# The start of a data line that is not 12 tab-separated fields of the strict
# grammar above, searched in "\n" + data: a literal "\n" lets re skip to each line
# start, where "^" with re.MULTILINE would try every character.
_NON_STRICT_OUTFMT6_LINE_RE = re.compile(
    "\n(?!"
    + "\t".join(
        [r"[^\t\n]*", r"[^\t\n]*"]
        + [
            _STRICT_DECIMAL_INTEGER if column in _INTEGER_COMPARISON_COLUMNS else _STRICT_DECIMAL_NUMBER
            for column in COMPARISON_COLUMNS[2:]
        ]
    )
    + "(?:\n|\\Z))"
)
LOSATP_METADATA_COLUMNS = (
    "query_protein_id",
    "subject_protein_id",
    "query_source_protein_id",
    "subject_source_protein_id",
    "query_record_index",
    "subject_record_index",
    "query_feature_index",
    "subject_feature_index",
    "query_feature_svg_id",
    "subject_feature_svg_id",
    "query_view_feature_svg_id",
    "subject_view_feature_svg_id",
    "orthogroup_id",
    "rbh_orthogroup_id",
    "ortholog_path_id",
    "edge_kind",
    "render_role",
    "query_orthogroup_representative",
    "subject_orthogroup_representative",
    "query_orthogroup_member_count",
    "subject_orthogroup_member_count",
    "query_orthogroup_role",
    "subject_orthogroup_role",
    "query_orthogroup_confidence",
    "subject_orthogroup_confidence",
    "query_orthogroup_assignment_reason",
    "subject_orthogroup_assignment_reason",
)
LOSATP_COMPARISON_COLUMNS = tuple(COMPARISON_COLUMNS) + LOSATP_METADATA_COLUMNS
PROTEIN_BLASTP_MODES = ("none", "pairwise", "orthogroup", "collinear")
ProteinBlastpMode = Literal["none", "pairwise", "orthogroup", "collinear"]
ORTHOGROUP_INFERENCE_VERSION = "anchor_core_v1"
ORTHOGROUP_MEMBERSHIP_MODES = (ORTHOGROUP_INFERENCE_VERSION,)
_LEGACY_ORTHOGROUP_MEMBERSHIP_MODES = ("rbh", "family_merge", "distribution_split")
# Legacy user and Session values that normalize_orthogroup_membership_mode reads
# as the active model (gbdraw.contract_shapes.read_spelling).
ORTHOGROUP_MEMBERSHIP_MODE_SPELLINGS: Spellings = {
    "aliases": {
        alias: ORTHOGROUP_INFERENCE_VERSION
        for alias in (
            "legacy", "rbh", "rbh_only", "merge", "family", "family_merge", "local_split",
            "density_split", "outparalog_split", "distribution_split", "orthogroups", "anchor_core",
        )
    },
    "foldHyphens": True,
    "blank": ORTHOGROUP_INFERENCE_VERSION,
}
OrthogroupMembershipMode = Literal["anchor_core_v1"]
OrthogroupScope = Literal["cross_record", "record_local"]
OrthogroupMemberRole = Literal["anchor", "coortholog", "inparalog", "low_confidence", "local_paralog"]
OrthogroupMemberConfidence = Literal["high", "medium", "low"]
OrthologEdgeKind = Literal[
    "rbh",
    "coortholog",
    "same_record_inparalog",
    "record_local_paralog",
    "related_homolog",
    "ambiguous_paralog",
    "weak_bridge",
    "domain_only",
    "related_paralog",
]
OrthologRenderRole = Literal["block_anchor", "display_edge"]
ORTHOLOG_DISPLAY_EDGE_KIND_RANK = {
    "rbh": 0,
    "coortholog": 1,
    "same_record_inparalog": 2,
    "record_local_paralog": 3,
    "related_homolog": 4,
    "ambiguous_paralog": 5,
    "weak_bridge": 6,
    "domain_only": 7,
    "related_paralog": 8,
}
_NEAR_RECIPROCAL_MIN_RATIO = 0.85
_INPARALOG_MIN_SAME_RECORD_RATIO_TO_LOCAL_ANCHOR = 0.75
_MIN_BEST_SECOND_CORE_RATIO = 1.25
_MIN_MEMBERSHIP_MIN_COVERAGE = 0.30
_DOMAIN_ONLY_MAX_MIN_COVERAGE = 0.25
_LOW_CONFIDENCE_MIN_BEST_SECOND_CORE_RATIO = 1.10
_ORTHOGROUP_SPLIT_TOP_FRACTION = 0.05
_ORTHOGROUP_SPLIT_MIN_COVERAGE = 0.30
_ORTHOGROUP_SPLIT_MIN_FIT_HITS = 12
_ORTHOGROUP_SPLIT_EPSILON = 1e-12
_RECORD_LOCAL_NEAR_RECIPROCAL_MIN_RATIO = 0.85
_RECORD_LOCAL_MIN_MEMBERS = 2
_RECORD_LOCAL_MIN_MEMBERSHIP_MIN_COVERAGE = _MIN_MEMBERSHIP_MIN_COVERAGE
_RECORD_LOCAL_CORE_COMPETITION_RATIO = 1.10


@dataclass(frozen=True)
class CdsProtein:
    """A protein sequence and its genomic CDS coordinates."""

    protein_id: str
    record_index: int
    feature_index: int
    record_id: str
    start: int
    end: int
    strand: int | None
    label: str
    protein_length: int
    sequence: str
    source_protein_id: str | None = None
    feature_svg_id: str | None = None
    gene: str | None = None
    product: str | None = None
    note: str | None = None
    locus_tag: str | None = None
    gene_id: str | None = None
    old_locus_tag: str | None = None
    db_xref: tuple[str, ...] = ()
    gff_id: str | None = None
    parent_ids: tuple[str, ...] = ()
    gene_parent_id: str | None = None
    feature_type: str = "CDS"
    feature_hash_start: int | None = None
    feature_hash_end: int | None = None
    feature_hash_strand: int | None = None
    feature_hash_parts: tuple[tuple[int, int, int | None], ...] = ()
    location_operator: str = ""
    source_feature_position: int | None = None
    same_location_ordinal: int | None = None
    feature_analysis_id: str | None = None
    display_alias: str | None = None
    runtime_handle: str | None = None
    aa_sha256: str | None = None
    record_instance_key: str | None = None
    record_analysis_id: str | None = None
    protein_set_hash: str | None = None
    runtime_binding_hash: str | None = None
    display_binding_hash: str | None = None
    view_feature_svg_id: str | None = None


def _refuse_change(self: object, *_args: object, **_kwargs: object) -> NoReturn:
    raise TypeError(
        f"{type(self).__name__} is read-only: a validated protein identity manifest "
        "cannot change. Change a copy from ProteinIdentityManifest.to_dict()."
    )


class _FrozenDict(dict[str, object]):
    """A read-only dict inside a protein identity manifest proof.

    It stays a ``dict``, so equality with a plain manifest and ``json.dumps``
    bytes are those of the same plain data. It is built from a mapping or
    from key-value pairs (as ``dataclasses.asdict`` rebuilds one) and freezes
    what it holds as the seal does; calling ``__init__`` again does not
    change it.
    """

    __slots__ = ()

    def __new__(
        cls,
        items: Mapping[str, object] | Iterable[tuple[str, object]] = (),
    ) -> _FrozenDict:
        flag = [False]
        frozen = dict.__new__(cls)
        dict.update(frozen, {key: _frozen(item, flag) for key, item in dict(items).items()})
        return frozen

    def __init__(self, *_args: object, **_kwargs: object) -> None:
        # dict.__init__ would fill the dict again; __new__ already did.
        pass

    __setitem__ = __delitem__ = __ior__ = _refuse_change
    clear = pop = popitem = setdefault = update = _refuse_change

    def __copy__(self) -> _FrozenDict:
        return self

    def __deepcopy__(self, _memo: object) -> _FrozenDict:
        return self

    def __reduce__(self) -> tuple[object, ...]:
        return (_FrozenDict, (dict(self),))


class _FrozenList(list[object]):
    """A read-only list inside a protein identity manifest proof.

    Like :class:`_FrozenDict`, it freezes what it holds.
    """

    __slots__ = ()

    def __new__(cls, items: Iterable[object] = ()) -> _FrozenList:
        flag = [False]
        frozen = list.__new__(cls)
        list.extend(frozen, [_frozen(item, flag) for item in items])
        return frozen

    def __init__(self, *_args: object, **_kwargs: object) -> None:
        # list.__init__ would empty the list; __new__ already filled it.
        pass

    __setitem__ = __delitem__ = __iadd__ = __imul__ = _refuse_change
    append = extend = insert = pop = remove = clear = sort = reverse = _refuse_change

    def __copy__(self) -> _FrozenList:
        return self

    def __deepcopy__(self, _memo: object) -> _FrozenList:
        return self

    def __reduce__(self) -> tuple[object, ...]:
        return (_FrozenList, (list(self),))


class ProteinIdentityManifest(_FrozenDict):
    """Schema-2 authority for protein sets, analyses, and runtime bindings.

    A proof: a read-only mapping whose items are the manifest. Calling the
    class checks its argument like :func:`validate_protein_identity_manifest`
    (which returns a proof unchanged); only that check and
    :func:`extract_protein_identity_manifest` seal one. It equals, and
    ``json.dumps`` writes, the same plain manifest; every mapping, sequence
    and set inside it is read-only, and ``copy.copy`` and ``copy.deepcopy``
    return it. Set items, dict keys and other objects, such as custom
    objects, are kept as given.
    """

    __slots__ = ("_json_exact",)
    _json_exact: bool

    schema = PROTEIN_IDENTITY_MANIFEST_SCHEMA

    def __new__(
        cls,
        items: Mapping[str, object] | Iterable[tuple[str, object]] = (),
        **fields: object,
    ) -> ProteinIdentityManifest:
        return validate_protein_identity_manifest({**dict(items), **fields})

    __setattr__ = __delattr__ = _refuse_change

    def __reduce__(self) -> tuple[object, ...]:
        return (validate_protein_identity_manifest, (self.to_dict(),))

    @property
    def protein_sets(self) -> Mapping[str, Mapping[str, object]]:
        return cast("Mapping[str, Mapping[str, object]]", self["proteinSets"])

    @property
    def record_analyses(self) -> Mapping[str, Mapping[str, object]]:
        return cast("Mapping[str, Mapping[str, object]]", self["recordAnalyses"])

    @property
    def record_instances(self) -> Mapping[str, Mapping[str, object]]:
        return cast("Mapping[str, Mapping[str, object]]", self["recordInstances"])

    def to_dict(self) -> dict[str, object]:
        """A detached, mutable plain copy of the manifest."""

        return cast("dict[str, object]", _thawed(self))

    @classmethod
    def from_dict(cls, value: Mapping[str, object]) -> ProteinIdentityManifest:
        return validate_protein_identity_manifest(value)

    def binding_for(self, record_instance_key: str) -> Mapping[str, object]:
        try:
            return self.record_instances[str(record_instance_key)]
        except KeyError as exc:
            raise ValidationError(
                f"Protein manifest has no record instance {record_instance_key!r}."
            ) from exc


def is_json_exact_protein_identity_manifest(
    value: object,
) -> TypeGuard[ProteinIdentityManifest]:
    """Whether ``value`` is a proof that holds only plain JSON data.

    Plain JSON data is data that equals its JSON clone: dicts with ``str``
    keys, lists, and exactly ``str``, ``int``, finite ``float``, ``bool`` or
    ``None``. Such a proof can stand where a deep copy or a JSON clone of its
    input stood; a proof of other Python values (tuples, sets, other
    mappings, other keys or scalars, NaN or infinities) is used only as the
    checked authority.
    """

    return isinstance(value, ProteinIdentityManifest) and value._json_exact


def _new_frozen_dict(items: Mapping[object, object]) -> _FrozenDict:
    """A read-only dict of ``items``, which ``_frozen`` has already frozen."""

    frozen = dict.__new__(_FrozenDict)
    dict.update(frozen, cast("Mapping[str, object]", items))
    return frozen


def _new_frozen_list(items: list[object]) -> _FrozenList:
    """A read-only list of ``items``, which ``_frozen`` has already frozen."""

    frozen = list.__new__(_FrozenList)
    list.extend(frozen, items)
    return frozen


def _frozen(value: object, inexact: list[bool]) -> object:
    """``value`` with every mapping and sequence in it read-only.

    A mapping becomes a read-only dict, a list a read-only list, any other
    sequence but ``str`` and ``bytes`` a tuple of frozen items, and a set a
    frozenset. Set items, dict keys and other objects, such as custom
    objects, are kept as they are.
    Sets ``inexact[0]`` when ``value`` holds anything but plain JSON data. A
    read-only container's constructor froze what it holds, so it is kept,
    and walked only to tell whether it is plain JSON data.
    """

    kind = type(value)
    if kind in _MANIFEST_SCALAR_TYPES:
        if kind is float and not math.isfinite(cast("float", value)):
            # A JSON clone of NaN or an infinity is a new float that does
            # not equal it, so it cannot stand for the input.
            inexact[0] = True
        return value
    if kind is _FrozenDict or kind is _FrozenList:
        if not inexact[0]:
            container = cast("_FrozenDict | _FrozenList", value)
            if kind is _FrozenDict and any(type(key) is not str for key in container):
                inexact[0] = True
            for item in container.values() if isinstance(container, dict) else container:
                if inexact[0]:
                    break
                _frozen(item, inexact)
        return value
    if kind is dict:
        mapping = cast("dict[object, object]", value)
        if not inexact[0] and any(type(key) is not str for key in mapping):
            inexact[0] = True
        return _new_frozen_dict({key: _frozen(item, inexact) for key, item in mapping.items()})
    if kind is list:
        return _new_frozen_list([_frozen(item, inexact) for item in cast("list[object]", value)])
    inexact[0] = True
    if isinstance(value, Mapping):
        return _new_frozen_dict({key: _frozen(item, inexact) for key, item in value.items()})
    if isinstance(value, list):
        return _new_frozen_list([_frozen(item, inexact) for item in value])
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        return tuple(_frozen(item, inexact) for item in value)
    if isinstance(value, AbstractSet):
        # Set items are hashable, and are kept as given, like dict keys.
        return frozenset(value)
    return value


def _thawed(value: object) -> object:
    """A detached copy of ``value`` with plain, mutable dicts and lists."""

    if type(value) in _MANIFEST_SCALAR_TYPES:
        return value
    if isinstance(value, dict):
        return {key: _thawed(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_thawed(item) for item in value]
    if isinstance(value, tuple):
        return tuple(_thawed(item) for item in value)
    return copy.deepcopy(value)


def _seal_protein_identity_manifest(
    items: Iterable[tuple[object, object]],
    *,
    inexact: bool = False,
) -> ProteinIdentityManifest:
    """The unchecked proof maker: only the check and the extraction call it."""

    flag = [inexact]
    frozen: dict[object, object] = {}
    for key, item in items:
        if type(key) is not str:
            flag[0] = True
        frozen[key] = _frozen(item, flag)
    proof = dict.__new__(ProteinIdentityManifest)
    # Keys that are not str make the proof inexact; they are kept as given.
    dict.update(proof, cast("dict[str, object]", frozen))
    object.__setattr__(proof, "_json_exact", not flag[0])
    return proof


@dataclass(frozen=True)
class ProteinLosatPairIdentity:
    """Directional raw-cache identity for one query/subject protein pair."""

    query_protein_set_hash: str
    subject_protein_set_hash: str
    query_runtime_binding_hash: str
    subject_runtime_binding_hash: str
    query_record_instance_key: str
    subject_record_instance_key: str

    def cache_payload(
        self,
        *,
        args: Sequence[str],
        program: str = "blastp",
        outfmt: str = "6",
    ) -> dict[str, object]:
        return {
            "cacheSchema": PROTEIN_LOSAT_CACHE_SCHEMA,
            "identityKind": "protein",
            "idEncoding": "runtime-handle-v1",
            "program": str(program),
            "outfmt": str(outfmt),
            "args": [str(arg) for arg in args],
            "queryProteinSetHash": self.query_protein_set_hash,
            "subjectProteinSetHash": self.subject_protein_set_hash,
            "queryRuntimeBindingHash": self.query_runtime_binding_hash,
            "subjectRuntimeBindingHash": self.subject_runtime_binding_hash,
            "queryRecordInstanceKey": self.query_record_instance_key,
            "subjectRecordInstanceKey": self.subject_record_instance_key,
        }

@dataclass(frozen=True)
class LegacyProteinRawCachePromotion:
    """A verified, copy-on-write schema-2 to schema-4 promotion."""

    candidate_index: int
    entry: dict[str, object]
    rewritten_tsv: str
    legacy_query_fasta: str
    legacy_subject_fasta: str


@dataclass(frozen=True)
class LegacyProteinRawCacheRejection:
    """Why a legacy candidate could not be promoted for the current pair."""

    candidate_index: int
    reason: str


@dataclass(frozen=True)
class LegacyProteinRawCachePromotionScan:
    """Promotion result plus pair-local rejection diagnostics."""

    promotion: LegacyProteinRawCachePromotion | None
    rejections: tuple[LegacyProteinRawCacheRejection, ...] = ()


@dataclass(frozen=True)
class ProteinExtractionResult:
    """CDS protein extraction output grouped for pairwise LOSATP runs."""

    proteins_by_record: list[list[CdsProtein]]
    protein_map: dict[str, CdsProtein]
    identity_manifest: ProteinIdentityManifest | None = None
    record_instance_keys: tuple[str, ...] = ()
    protein_set_hashes: tuple[str, ...] = ()
    record_analysis_ids: tuple[str, ...] = ()
    runtime_binding_hashes: tuple[str, ...] = ()
    display_binding_hashes: tuple[str, ...] = ()


@dataclass(frozen=True)
class OrthogroupMember:
    """A CDS/protein member assigned to an orthogroup."""

    orthogroup_id: str
    protein_id: str
    record_index: int
    feature_index: int
    record_id: str
    label: str
    start: int
    end: int
    strand: int | None
    feature_svg_id: str | None
    source_protein_id: str | None
    gene: str | None = None
    product: str | None = None
    note: str | None = None
    representative: bool = False
    role: OrthogroupMemberRole = "anchor"
    confidence: OrthogroupMemberConfidence = "high"
    assignment_reason: str = ""
    supporting_edges: tuple[str, ...] = ()
    best_core_support: float = 0.0
    second_best_core_support: float = 0.0


@dataclass(frozen=True)
class OrthologEdge:
    """Selected evidence edge describing relationships inside or near an orthogroup."""

    orthogroup_id: str
    source_rbh_orthogroup_id: str | None
    target_rbh_orthogroup_id: str | None
    query_protein_id: str
    subject_protein_id: str
    query_record_index: int
    subject_record_index: int
    edge_kind: OrthologEdgeKind
    render_role: OrthologRenderRole
    path_id: str | None
    identity: float
    evalue: float
    bitscore: float
    alignment_length: int


@dataclass(frozen=True)
class OrthogroupNameCandidate:
    """Annotation-derived display-name candidate for an orthogroup."""

    text: str
    source: str
    member_count: int
    record_coverage_count: int
    representative_count: int
    score: float


@dataclass(frozen=True)
class OrthogroupResult:
    """Connected-component orthogroups and per-protein lookup metadata."""

    orthogroups: dict[str, list[OrthogroupMember]]
    member_by_protein_id: dict[str, OrthogroupMember]
    names_by_orthogroup_id: dict[str, str] = field(default_factory=dict)
    descriptions_by_orthogroup_id: dict[str, str] = field(default_factory=dict)
    name_candidates_by_orthogroup_id: dict[str, list[OrthogroupNameCandidate]] = field(default_factory=dict)
    confidence_by_orthogroup_id: dict[str, str] = field(default_factory=dict)
    rbh_orthogroups: dict[str, tuple[str, ...]] = field(default_factory=dict)
    ortholog_edges_by_orthogroup_id: dict[str, tuple[OrthologEdge, ...]] = field(default_factory=dict)
    ortholog_paths_by_orthogroup_id: dict[str, tuple[OrthologPath, ...]] = field(default_factory=dict)
    related_edges_by_orthogroup_id: dict[str, tuple[OrthologEdge, ...]] = field(default_factory=dict)
    scope_by_orthogroup_id: dict[str, OrthogroupScope] = field(default_factory=dict)
    source_record_index_by_orthogroup_id: dict[str, int] = field(default_factory=dict)


@dataclass(frozen=True)
class OrthogroupGraphResult:
    """Orthogroups with compact lossless path indexes (the default result)."""

    orthogroups: dict[str, list[OrthogroupMember]]
    member_by_protein_id: dict[str, OrthogroupMember]
    names_by_orthogroup_id: dict[str, str] = field(default_factory=dict)
    descriptions_by_orthogroup_id: dict[str, str] = field(default_factory=dict)
    name_candidates_by_orthogroup_id: dict[str, list[OrthogroupNameCandidate]] = field(default_factory=dict)
    confidence_by_orthogroup_id: dict[str, str] = field(default_factory=dict)
    rbh_orthogroups: dict[str, tuple[str, ...]] = field(default_factory=dict)
    ortholog_edges_by_orthogroup_id: dict[str, tuple[OrthologEdge, ...]] = field(default_factory=dict)
    path_indexes_by_orthogroup_id: dict[str, OrthologPathCollection] = field(default_factory=dict)
    related_edges_by_orthogroup_id: dict[str, tuple[OrthologEdge, ...]] = field(default_factory=dict)
    scope_by_orthogroup_id: dict[str, OrthogroupScope] = field(default_factory=dict)
    source_record_index_by_orthogroup_id: dict[str, int] = field(default_factory=dict)

    def __post_init__(self):
        for group_id, index in self.path_indexes_by_orthogroup_id.items():
            if not isinstance(index, OrthologPathCollection) or group_id != index.orthogroup_id:
                raise ValidationError("Orthogroup path index has an inconsistent group ID")
            index.validate_edges(self.ortholog_edges_by_orthogroup_id.get(group_id, ()))


def _orthogroup_metadata_fields(result):
    return {name: getattr(result, name) for name in OrthogroupResult.__dataclass_fields__
            if name != "ortholog_paths_by_orthogroup_id"}


def compact_ortholog_paths(result: OrthogroupResult | OrthogroupGraphResult) -> OrthogroupGraphResult:
    """Adapt a legacy corpus once, without reconstructing its edge closure."""
    if isinstance(result, OrthogroupGraphResult):
        return result
    return OrthogroupGraphResult(
        **_orthogroup_metadata_fields(result),
        path_indexes_by_orthogroup_id={
            group_id: OrthologPathCollection(group_id, "explicit", paths=paths)
            for group_id, paths in result.ortholog_paths_by_orthogroup_id.items()
        },
    )


def materialize_ortholog_paths(result: OrthogroupResult | OrthogroupGraphResult) -> OrthogroupResult:
    """Explicit exhaustive legacy output; time and memory grow with all paths."""
    if isinstance(result, OrthogroupResult):
        return result
    return OrthogroupResult(
        **_orthogroup_metadata_fields(result),
        ortholog_paths_by_orthogroup_id={
            group_id: tuple(index.iter_paths())
            for group_id, index in result.path_indexes_by_orthogroup_id.items()
        },
    )


@dataclass(frozen=True)
class ProteinBlastpResult:
    """LOSATP blastp display comparisons plus optional orthogroup metadata."""

    comparisons: list[DataFrame]
    orthogroups: OrthogroupResult | OrthogroupGraphResult | None = None


@dataclass(frozen=True)
class OrthogroupEdgeSelectionResult:
    """Orthogroup selection plus separate adjacent anchor and display edges."""

    orthogroups: OrthogroupResult | OrthogroupGraphResult
    all_edges_by_pair: dict[tuple[int, int], DataFrame]
    adjacent_anchor_edges_by_pair: dict[tuple[int, int], DataFrame]
    adjacent_display_edges_by_pair: dict[tuple[int, int], DataFrame]


class _ComparisonColumns:
    """One normalized table's values that only output edges read, row by row.

    The comparison columns are copies, so the table itself can be freed once
    the hits are built; a column view would keep its whole block alive.
    """

    __slots__ = ("arrays", "identity", "bitscore", "alignment_length")

    def __init__(self, table: DataFrame, identity: np.ndarray, bitscore: np.ndarray, alignment_length: np.ndarray) -> None:
        self.arrays = tuple(table[column].to_numpy().copy() for column in COMPARISON_COLUMNS[2:])
        self.identity = identity
        self.bitscore = bitscore
        self.alignment_length = alignment_length

    def values(self, row: int) -> tuple[object, ...]:
        # ndarray.item gives the Python scalar that itertuples gives.
        return tuple(array.item(row) for array in self.arrays)


class _DirectedHit(NamedTuple):
    """The best hit from one protein to another, with every value inference reads.

    ``rank`` is the hit's position in the anchor-core order: higher normalized
    score, lower e-value, higher minimum coverage, higher identity, longer
    alignment, then query record, query ID, subject record and subject ID.
    ``comparison_row`` is the hit's row in ``comparison_source``, which also
    holds the values only output edges read.
    """

    query: str
    subject: str
    query_record_index: int
    subject_record_index: int
    score: float
    evalue: float
    min_coverage: float
    domain_only: bool
    supports_membership: bool
    rank: int
    comparison_source: _ComparisonColumns
    comparison_row: int

    @property
    def identity(self) -> float:
        return float(self.comparison_source.identity.item(self.comparison_row))

    @property
    def bitscore(self) -> float:
        return float(self.comparison_source.bitscore.item(self.comparison_row))

    @property
    def alignment_length(self) -> float:
        return float(self.comparison_source.alignment_length.item(self.comparison_row))


@dataclass(frozen=True, eq=False)
class _DirectedHitTable:
    """The best normalized hit per (query, subject) direction, one row per direction.

    Rows are in rank order (the row is the hit's ``rank``). Proteins are integer
    codes in (record index, ID) order, so the lower code of a pair is its
    canonical query: the earlier record, then the smaller ID. ``reverse`` is the
    row of the opposite direction, or -1; ``source`` and ``source_row`` locate the
    row in ``sources``. ``hit`` builds the record of one row, for the few hits
    that later steps keep.
    """

    protein_ids: list[str]
    protein_index: pd.Index
    protein_record: np.ndarray
    query: np.ndarray
    subject: np.ndarray
    score: np.ndarray
    evalue: np.ndarray
    min_coverage: np.ndarray
    domain_only: np.ndarray
    supports_membership: np.ndarray
    reverse: np.ndarray
    sources: list[_ComparisonColumns]
    source: np.ndarray
    source_row: np.ndarray

    def __len__(self) -> int:
        return len(self.query)

    def codes(self, protein_ids: Iterable[str]) -> np.ndarray:
        """Codes of the given IDs; -1 for an ID outside the protein map."""
        return np.asarray(self.protein_index.get_indexer(list(protein_ids)), dtype=np.int64)

    def hit(self, row: int) -> _DirectedHit:
        # ndarray.item gives the Python scalar that tolist gives.
        query, subject = self.query.item(row), self.subject.item(row)
        return _DirectedHit(
            self.protein_ids[query],
            self.protein_ids[subject],
            self.protein_record.item(query),
            self.protein_record.item(subject),
            self.score.item(row),
            self.evalue.item(row),
            self.min_coverage.item(row),
            self.domain_only.item(row),
            self.supports_membership.item(row),
            row,
            self.sources[self.source.item(row)],
            self.source_row.item(row),
        )


@dataclass(frozen=True)
class _AnchorCoreEvidenceEdge:
    query_id: str
    subject_id: str
    hit: _DirectedHit
    score: float
    reverse_score: float
    edge_kind: Literal["rbh", "coortholog", "record_local_paralog"]


@dataclass(frozen=True)
class _CoreSupportCandidate:
    group_id: str
    support: float
    diagnostic_score: float
    same_record_score: float
    cross_record_score: float
    evidence: _DirectedHit
    relation_kind: OrthologEdgeKind
    min_coverage: float
    domain_only: bool
    high_confidence_pass: bool
    low_confidence_pass: bool
    reason: str


@dataclass(frozen=True)
class _BestCoreEvidence:
    """The best membership-supporting hit and the best hit of any coverage."""

    support: _DirectedHit | None
    diagnostic: _DirectedHit | None


@dataclass(frozen=True)
class _CdsProteinCandidate:
    protein_id: str
    source_protein_id: str | None
    record_index: int
    feature_index: int
    record_id: str
    start: int
    end: int
    strand: int | None
    label: str
    protein_length: int
    sequence: str
    feature_svg_id: str | None
    view_feature_svg_id: str | None
    gene: str | None
    product: str | None
    note: str | None
    locus_tag: str | None
    gene_id: str | None
    old_locus_tag: str | None
    db_xref: tuple[str, ...]
    gff_id: str | None
    parent_ids: tuple[str, ...]
    gene_parent_id: str | None
    feature_type: str
    feature_hash_start: int | None
    feature_hash_end: int | None
    feature_hash_strand: int | None
    feature_hash_parts: tuple[tuple[int, int, int | None], ...]
    location_operator: str
    source_feature_position: int


LosatpRunner = Callable[[str, str], DataFrame]


def _hash_text_sha256(text: str) -> str:
    return hashlib.sha256(str(text).encode("utf-8")).hexdigest()


def _web_json_dumps(value: object) -> str:
    return json.dumps(value, ensure_ascii=False, separators=(",", ":"))


def canonical_protein_identity_json(value: object) -> str:
    """Serialize a protein identity value using the normative JSON contract."""

    try:
        return json.dumps(
            value,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        )
    except (TypeError, ValueError) as exc:
        raise ValidationError(f"Protein identity value is not canonical JSON: {exc}") from exc


def _identity_sha256(value: object) -> str:
    digest = hashlib.sha256(canonical_protein_identity_json(value).encode("utf-8")).hexdigest()
    return f"sha256:{digest}"


def _normalize_identity_text(value: object) -> str:
    return unicodedata.normalize("NFC", str(value or ""))


def _canonical_identity_context_value(value: object) -> object:
    if value is None or isinstance(value, (bool, int, float)):
        return value
    if isinstance(value, str):
        return _normalize_identity_text(value)
    if isinstance(value, Mapping):
        return {
            _normalize_identity_text(key): _canonical_identity_context_value(item)
            for key, item in value.items()
        }
    if isinstance(value, (list, tuple)):
        return [_canonical_identity_context_value(item) for item in value]
    raise ValidationError(
        f"Protein record analysis context contains unsupported {type(value).__name__}."
    )


def percent_encode_losat_transport_field(value: object) -> str:
    """Encode one transport field with NFC/UTF-8 and the fixed safe alphabet."""

    safe = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789._-"
    safe_bytes = set(safe)
    encoded: list[str] = []
    for byte in _normalize_identity_text(value).encode("utf-8"):
        encoded.append(chr(byte) if byte in safe_bytes else f"%{byte:02X}")
    return "".join(encoded)


def build_protein_runtime_handle(
    *,
    feature_analysis_id: str,
    record_instance_key: object,
) -> str:
    """Build the deterministic session-global runtime handle for one CDS."""

    feature_id = str(feature_analysis_id)
    instance_key = _normalize_identity_text(record_instance_key).strip()
    if not re.fullmatch(r"f_[0-9a-f]{64}", feature_id):
        raise ValidationError(
            "feature_analysis_id must be f_ followed by 64 lowercase hex digits."
        )
    if not instance_key:
        raise ValidationError("record_instance_key must be non-empty.")
    payload = canonical_protein_identity_json(
        {
            "featureAnalysisId": feature_id,
            "recordInstanceKey": instance_key,
        }
    ).encode("utf-8")
    digest = hashlib.sha256(b"gbdraw-runtime-handle-v1\0" + payload).digest()[:16]
    token = base64.b32encode(digest).decode("ascii").lower().rstrip("=")
    handle = f"h_{token}"
    if not _VALID_RUNTIME_HANDLE_RE.fullmatch(handle):  # pragma: no cover
        raise ValidationError("Generated protein runtime handle is invalid.")
    return handle


def canonical_feature_analysis_id(
    *,
    feature_type: str,
    location_operator: str | None,
    location_parts: Sequence[Sequence[int | None]],
    strand: int | None,
    same_location_ordinal: int,
) -> str:
    """Return the canonical, qualifier-independent ID for one CDS feature."""

    ordinal = int(same_location_ordinal)
    if ordinal <= 0:
        raise ValidationError("same_location_ordinal must be a positive integer.")
    normalized_parts: list[list[int | None]] = []
    for raw_part in location_parts:
        if len(raw_part) != 3:
            raise ValidationError("Each location part must contain start, end, and strand.")
        # Only the strand part may be None.
        start = int(cast(int, raw_part[0]))
        end = int(cast(int, raw_part[1]))
        part_strand = int(raw_part[2]) if raw_part[2] in {-1, 1} else None
        if start < 0 or end <= start:
            raise ValidationError("Feature identity location parts must be 0-based, end-exclusive intervals.")
        normalized_parts.append([start, end, part_strand])
    if not normalized_parts:
        raise ValidationError("Feature identity requires at least one location part.")
    normalized_strand = int(strand) if strand in {-1, 1} else None
    payload = {
        "featureType": _normalize_identity_text(feature_type),
        "locationOperator": _normalize_identity_text(location_operator),
        "locationParts": normalized_parts,
        "sameLocationOrdinal": ordinal,
        "strand": normalized_strand,
    }
    digest = hashlib.sha256(canonical_protein_identity_json(payload).encode("utf-8")).hexdigest()
    return f"f_{digest}"


def _strand_display(value: int | None) -> str:
    if value == 1:
        return "+"
    if value == -1:
        return "-"
    return "0"


def select_protein_display_alias(protein: CdsProtein) -> str:
    """Choose the readable alias without using it as machine identity."""

    for value in (protein.source_protein_id, protein.locus_tag, protein.gff_id):
        normalized = _normalize_identity_text(value).strip()
        if normalized:
            return normalized
    parts = protein.feature_hash_parts or (
        (protein.start, protein.end, protein.strand),
    )
    location = ",".join(f"{int(start)}-{int(end)}" for start, end, _ in parts)
    operator = _normalize_identity_text(protein.location_operator).strip()
    if operator and len(parts) > 1:
        location = f"{operator}({location})"
    return f"CDS:{location}:{_strand_display(protein.feature_hash_strand)}"


def _validate_protein_sets_and_record_analyses(
    raw_sets: Mapping[object, object],
    raw_analyses: Mapping[object, object],
) -> tuple[dict[str, dict[str, object]], dict[str, dict[str, object]]]:
    """Validate manifest authorities shared independently of top-level schema."""

    protein_sets: dict[str, dict[str, object]] = {}
    for key, raw_set in raw_sets.items():
        set_key = str(key)
        if not re.fullmatch(r"sha256:[0-9a-f]{64}", set_key) or not isinstance(raw_set, Mapping):
            raise ValidationError("Protein set keys must be sha256 digests.")
        raw_proteins = raw_set.get("proteins")
        if raw_set.get("schema") != 1 or not isinstance(raw_proteins, list):
            raise ValidationError(f"Protein set {set_key!r} is invalid.")
        identity_pairs: list[tuple[str, str]] = []
        feature_ids: set[str] = set()
        for protein in raw_proteins:
            if not isinstance(protein, Mapping):
                raise ValidationError(f"Protein set {set_key!r} contains an invalid protein.")
            feature_id = str(protein.get("featureAnalysisId") or "")
            aa_sha256 = str(protein.get("aaSha256") or "")
            if not re.fullmatch(r"f_[0-9a-f]{64}", feature_id):
                raise ValidationError(f"Protein set {set_key!r} contains an invalid feature ID.")
            if not re.fullmatch(r"[0-9a-f]{64}", aa_sha256) or feature_id in feature_ids:
                raise ValidationError(f"Protein set {set_key!r} has invalid or duplicate protein identity.")
            try:
                expected_feature_id = canonical_feature_analysis_id(
                    feature_type=str(protein.get("featureType") or ""),
                    location_operator=str(protein.get("locationOperator") or ""),
                    location_parts=protein.get("locationParts") or (),
                    strand=protein.get("strand"),
                    same_location_ordinal=int(protein.get("sameLocationOrdinal") or 0),
                )
            except (TypeError, ValueError) as exc:
                raise ValidationError(f"Protein set {set_key!r} contains invalid feature identity fields.") from exc
            if expected_feature_id != feature_id:
                raise ValidationError(f"Protein set {set_key!r} feature ID does not match its content.")
            feature_ids.add(feature_id)
            identity_pairs.append((feature_id, aa_sha256))
        identity_pairs.sort(key=lambda pair: pair[0])
        if [
            str(protein.get("featureAnalysisId"))
            for protein in raw_proteins
        ] != [feature_id for feature_id, _aa_sha256 in identity_pairs]:
            raise ValidationError(f"Protein set {set_key!r} is not in canonical feature order.")
        if _identity_sha256(_protein_set_hash_payload(identity_pairs)) != set_key:
            raise ValidationError(f"Protein set {set_key!r} hash does not match its content.")
        protein_sets[set_key] = dict(raw_set)

    record_analyses: dict[str, dict[str, object]] = {}
    for key, raw_analysis in raw_analyses.items():
        analysis_key = str(key)
        if not re.fullmatch(r"sha256:[0-9a-f]{64}", analysis_key) or not isinstance(raw_analysis, Mapping):
            raise ValidationError("Record analysis keys must be sha256 digests.")
        protein_set_hash = str(raw_analysis.get("proteinSetHash") or "")
        record_source_id = raw_analysis.get("recordSourceId")
        if (
            raw_analysis.get("schema") != 1
            or protein_set_hash not in protein_sets
            or not isinstance(record_source_id, str)
            or not record_source_id
        ):
            raise ValidationError(f"Record analysis {analysis_key!r} has an invalid protein set reference.")
        analysis_hash_payload = _record_analysis_hash_payload(
            protein_set_hash=protein_set_hash,
            record_source_id=raw_analysis.get("recordSourceId"),
            region=raw_analysis.get("region"),
            selector=raw_analysis.get("selector"),
        )
        if _identity_sha256(analysis_hash_payload) != analysis_key:
            raise ValidationError(f"Record analysis {analysis_key!r} hash does not match its content.")
        record_analyses[analysis_key] = dict(raw_analysis)

    return protein_sets, record_analyses


def _protein_set_hash_payload(
    identity_pairs: Sequence[tuple[str, str]],
) -> dict[str, object]:
    return {
        "proteins": [
            {"featureAnalysisId": feature_id, "aaSha256": aa_sha256}
            for feature_id, aa_sha256 in identity_pairs
        ]
    }


def _record_analysis_hash_payload(
    *,
    protein_set_hash: str,
    record_source_id: object,
    region: object,
    selector: object,
) -> dict[str, object]:
    return {
        "proteinSetHash": protein_set_hash,
        "recordSourceId": record_source_id,
        "region": region,
        "selector": selector,
    }


def _runtime_binding_payload(
    *,
    protein_set_hash: str,
    record_instance_key: str,
    runtime_ids: Mapping[str, str],
) -> dict[str, object]:
    return {
        "encoding": "runtime-handle-v1",
        "proteinSetHash": str(protein_set_hash),
        "recordInstanceKey": str(record_instance_key),
        "runtimeIds": dict(runtime_ids),
    }


def _display_binding_payload(
    *,
    record_analysis_id: str,
    record_source_id: str,
    record_instance_key: str,
    feature_metadata: Mapping[str, object],
) -> dict[str, object]:
    return {
        "recordAnalysisId": str(record_analysis_id),
        "recordSourceId": str(record_source_id),
        "recordInstanceKey": str(record_instance_key),
        "featureMetadata": dict(feature_metadata),
    }


def _expected_export_ordinals(
    display_aliases: Mapping[str, object],
) -> dict[str, int | None]:
    feature_ids_by_alias: dict[str, list[str]] = {}
    for feature_id, display_alias in display_aliases.items():
        alias = _normalize_identity_text(display_alias).strip()
        if not alias:
            raise ValidationError("Protein feature metadata requires a display alias.")
        feature_ids_by_alias.setdefault(alias, []).append(str(feature_id))
    expected: dict[str, int | None] = {}
    for feature_ids in feature_ids_by_alias.values():
        ordered = sorted(feature_ids)
        for ordinal, feature_id in enumerate(ordered, start=1):
            expected[feature_id] = ordinal if len(ordered) > 1 else None
    return expected


@lru_cache(maxsize=4096)
def _is_presentation_identity_key(raw_key: str) -> bool:
    normalized_key = re.sub(
        r"[^a-z0-9]+",
        "",
        unicodedata.normalize("NFKC", raw_key).casefold(),
    )
    return any(
        fragment in normalized_key
        for fragment in _MANIFEST_PRESENTATION_IDENTITY_FRAGMENTS
    )


def _reject_manifest_presentation_identity(value: Mapping[str, object]) -> None:
    """Keep protein identity manifests independent of rendered view coordinates."""

    pending: list[object] = [value]
    visited: set[int] = set()
    while pending:
        current = pending.pop()
        if type(current) is dict or isinstance(current, Mapping):
            object_id = id(current)
            if object_id in visited:
                continue
            visited.add(object_id)
            for raw_key, item in current.items():
                if _is_presentation_identity_key(str(raw_key)):
                    raise ValidationError(
                        "Protein identity manifests cannot contain rendered-view identity fields."
                    )
                # Only containers can hold more keys.
                if type(item) not in _MANIFEST_SCALAR_TYPES:
                    pending.append(item)
        elif type(current) is list or isinstance(current, (list, tuple)):
            pending.extend(
                item for item in current if type(item) not in _MANIFEST_SCALAR_TYPES
            )


def validate_protein_identity_manifest(
    value: Mapping[str, object],
) -> ProteinIdentityManifest:
    """Check a schema-2 compact protein manifest and return its proof.

    A :class:`ProteinIdentityManifest` is returned unchanged: it passed this
    check when it was made and cannot change. Any other value gets the full
    check, and its proof keeps the value's items in their order, with the
    keys and runtime handles the check reads as text made ``str``. For JSON
    data the proof equals the value and ``json.dumps`` writes the same bytes.
    """

    if isinstance(value, ProteinIdentityManifest):
        return value
    checked = dict(
        zip(
            ("proteinSets", "recordAnalyses", "recordInstances"),
            _check_protein_identity_manifest(value),
        )
    )
    # Plain JSON data is sealed as it is: the check's maps hold the same
    # content then. Other Python values get the check's str() coercions.
    proof = _seal_protein_identity_manifest(value.items(), inexact=type(value) is not dict)
    if proof._json_exact:
        return proof
    return _seal_protein_identity_manifest(
        ((key, checked[key] if key in checked else item) for key, item in value.items()),
        inexact=True,
    )


def _check_protein_identity_manifest(
    value: Mapping[str, object],
) -> tuple[
    dict[str, dict[str, object]],
    dict[str, dict[str, object]],
    dict[str, dict[str, object]],
]:
    """The full check: the protein sets, analyses and bindings it read, in input order."""

    if not isinstance(value, Mapping) or value.get("schema") != PROTEIN_IDENTITY_MANIFEST_SCHEMA:
        raise ValidationError("Protein identity manifest must use schema 2.")
    _reject_manifest_presentation_identity(value)
    raw_sets = value.get("proteinSets")
    raw_analyses = value.get("recordAnalyses")
    raw_instances = value.get("recordInstances")
    if (
        not isinstance(raw_sets, Mapping)
        or not isinstance(raw_analyses, Mapping)
        or not isinstance(raw_instances, Mapping)
    ):
        raise ValidationError("Protein identity manifest maps are missing or invalid.")

    protein_sets, record_analyses = _validate_protein_sets_and_record_analyses(
        raw_sets,
        raw_analyses,
    )

    record_instances: dict[str, dict[str, object]] = {}
    all_runtime_handles: set[str] = set()
    for key, raw_binding in raw_instances.items():
        instance_key = str(key)
        if not instance_key or not isinstance(raw_binding, Mapping) or raw_binding.get("schema") != 2:
            raise ValidationError("Record instance bindings must be non-empty schema-2 objects.")
        analysis_id = str(raw_binding.get("recordAnalysisId") or "")
        runtime_binding_hash = str(raw_binding.get("runtimeBindingHash") or "")
        display_binding_hash = str(raw_binding.get("displayBindingHash") or "")
        runtime_ids = raw_binding.get("runtimeIds")
        feature_metadata = raw_binding.get("featureMetadata")
        if analysis_id not in record_analyses:
            raise ValidationError(
                f"Record instance {instance_key!r} has an invalid analysis reference."
            )
        if not re.fullmatch(r"sha256:[0-9a-f]{64}", runtime_binding_hash):
            raise ValidationError(
                f"Record instance {instance_key!r} has an invalid runtime binding hash."
            )
        if not re.fullmatch(r"sha256:[0-9a-f]{64}", display_binding_hash):
            raise ValidationError(
                f"Record instance {instance_key!r} has an invalid display binding hash."
            )
        if not isinstance(runtime_ids, Mapping) or not isinstance(feature_metadata, Mapping):
            raise ValidationError(
                f"Record instance {instance_key!r} has invalid runtime or metadata maps."
            )

        protein_set_hash = str(record_analyses[analysis_id]["proteinSetHash"])
        expected_feature_ids = {
            str(protein["featureAnalysisId"])
            for protein in cast(
                "Sequence[Mapping[str, object]]", protein_sets[protein_set_hash]["proteins"]
            )
        }
        normalized_runtime_ids: dict[str, str] = {}
        for feature_id, runtime_handle in runtime_ids.items():
            feature_key = str(feature_id)
            handle = str(runtime_handle)
            if feature_key not in expected_feature_ids:
                raise ValidationError("Runtime map contains an unknown feature analysis ID.")
            if (
                not _VALID_RUNTIME_HANDLE_RE.fullmatch(handle)
                or handle
                != build_protein_runtime_handle(
                    feature_analysis_id=feature_key,
                    record_instance_key=instance_key,
                )
                or handle in all_runtime_handles
            ):
                raise ValidationError(
                    "Runtime handles must be deterministic and globally unique."
                )
            normalized_runtime_ids[feature_key] = handle
            all_runtime_handles.add(handle)
        if set(normalized_runtime_ids) != expected_feature_ids:
            raise ValidationError(
                f"Record instance {instance_key!r} runtime map does not match its protein set."
            )

        normalized_metadata: dict[str, dict[str, object]] = {}
        for feature_id, metadata in feature_metadata.items():
            feature_key = str(feature_id)
            if feature_key not in expected_feature_ids or not isinstance(metadata, Mapping):
                raise ValidationError(
                    f"Record instance {instance_key!r} contains invalid feature metadata."
                )
            normalized_metadata[feature_key] = dict(metadata)
        if set(normalized_metadata) != expected_feature_ids:
            raise ValidationError(
                f"Record instance {instance_key!r} metadata does not match its protein set."
            )
        expected_ordinals = _expected_export_ordinals(
            {
                feature_id: metadata.get("displayAlias")
                for feature_id, metadata in normalized_metadata.items()
            }
        )
        for feature_id, expected_ordinal in expected_ordinals.items():
            actual = normalized_metadata[feature_id].get("exportOrdinal")
            if actual != expected_ordinal:
                raise ValidationError(
                    f"Record instance {instance_key!r} has a non-canonical export ordinal."
                )

        runtime_payload = _runtime_binding_payload(
            protein_set_hash=protein_set_hash,
            record_instance_key=instance_key,
            runtime_ids=normalized_runtime_ids,
        )
        if _identity_sha256(runtime_payload) != runtime_binding_hash:
            raise ValidationError(
                f"Record instance {instance_key!r} runtime binding hash does not match."
            )
        display_payload = _display_binding_payload(
            record_analysis_id=analysis_id,
            record_source_id=str(record_analyses[analysis_id]["recordSourceId"]),
            record_instance_key=instance_key,
            feature_metadata=normalized_metadata,
        )
        if _identity_sha256(display_payload) != display_binding_hash:
            raise ValidationError(
                f"Record instance {instance_key!r} display binding hash does not match."
            )

        normalized_binding = dict(raw_binding)
        normalized_binding["runtimeIds"] = normalized_runtime_ids
        normalized_binding["featureMetadata"] = normalized_metadata
        record_instances[instance_key] = normalized_binding

    return protein_sets, record_analyses, record_instances


def _safe_web_protein_id_token(value: object) -> str:
    """Legacy schema-2 token sanitizer used only for verified migration."""

    text = str(value or "").strip()
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", text)
    text = text.strip("._-")
    return text or "record"


def _with_stable_web_protein_ids(
    proteins: Sequence[CdsProtein],
    record_instance_key: object,
) -> list[CdsProtein]:
    """Reconstruct legacy schema-2 Web IDs for verified migration only."""

    record_key = _safe_web_protein_id_token(record_instance_key)
    remapped: list[CdsProtein] = []
    used: set[str] = set()
    for protein in proteins:
        strand = protein.strand if protein.strand in {-1, 1} else 0
        aa_hash = hashlib.sha256(str(protein.sequence or "").encode("utf-8")).hexdigest()[:12]
        base_id = f"p_{record_key}_{int(protein.start)}_{int(protein.end)}_{strand}_{aa_hash}"
        protein_id = base_id
        suffix = 2
        while protein_id in used:
            protein_id = f"{base_id}_{suffix}"
            suffix += 1
        used.add(protein_id)
        remapped.append(replace(protein, protein_id=protein_id))
    return remapped


def extract_web_stable_cds_proteins(
    records: Sequence[SeqRecord],
    *,
    record_instance_keys: Sequence[str] | None = None,
    record_source_ids: Sequence[str] | None = None,
    record_selectors: Sequence[object] | None = None,
    regions: Sequence[object] | None = None,
    record_index_offset: int = 0,
    feature_visibility_rules: list[dict[str, object]] | None = None,
    feature_overrides: Sequence[Mapping[int, object]] | None = None,
) -> ProteinExtractionResult:
    """Extract CDS proteins with the current manifest and runtime handles."""

    return extract_protein_identity_manifest(
        records,
        record_instance_keys=record_instance_keys,
        record_source_ids=record_source_ids,
        record_selectors=record_selectors,
        regions=regions,
        record_index_offset=record_index_offset,
        feature_visibility_rules=feature_visibility_rules,
        feature_overrides=feature_overrides,
    )


def _per_record_values(
    values: Sequence[object] | None,
    records: Sequence[SeqRecord],
    *,
    name: str,
    default: Callable[[int, SeqRecord], object],
) -> list[object]:
    if values is not None and len(values) != len(records):
        raise ValidationError(f"{name} must have one value per record.")
    if values is None:
        return [default(index, record) for index, record in enumerate(records)]
    return list(values)


def _protein_location_identity_key(protein: CdsProtein) -> tuple[object, ...]:
    return (
        _normalize_identity_text(protein.feature_type),
        _normalize_identity_text(protein.location_operator),
        tuple(
            (int(start), int(end), int(strand) if strand in {-1, 1} else None)
            for start, end, strand in protein.feature_hash_parts
        ),
        (
            int(protein.feature_hash_strand)
            if protein.feature_hash_strand in {-1, 1}
            else None
        ),
    )


def _same_location_sort_key(protein: CdsProtein) -> tuple[object, ...]:
    aa_sha256 = hashlib.sha256(protein.sequence.encode("utf-8")).hexdigest()
    # Display qualifiers are intentionally excluded; source position only breaks
    # ties between otherwise identical same-location protein sequences.
    source_position = (
        int(protein.source_feature_position)
        if protein.source_feature_position is not None
        else int(protein.feature_index)
    )
    return aa_sha256, source_position


def _protein_identity_core(protein: CdsProtein) -> dict[str, object]:
    return {
        "featureAnalysisId": str(protein.feature_analysis_id),
        "sameLocationOrdinal": int(protein.same_location_ordinal or 0),
        "featureType": protein.feature_type,
        "locationOperator": protein.location_operator,
        "locationParts": [list(part) for part in protein.feature_hash_parts],
        "strand": protein.feature_hash_strand,
        "aaSha256": str(protein.aa_sha256),
    }


def _protein_binding_metadata(
    protein: CdsProtein,
    *,
    export_ordinal: int | None,
) -> dict[str, object]:
    return {
        "displayAlias": protein.display_alias,
        "exportOrdinal": export_ordinal,
        "sourceFeaturePosition": protein.source_feature_position,
        "sourceProteinId": protein.source_protein_id,
        "locusTag": protein.locus_tag,
        "gffId": protein.gff_id,
        "featureSvgId": protein.feature_svg_id,
        "gene": protein.gene,
        "product": protein.product,
        "note": protein.note,
        "oldLocusTag": protein.old_locus_tag,
        "geneId": protein.gene_id,
        "dbXref": list(protein.db_xref),
        "parentIds": list(protein.parent_ids),
        "geneParentId": protein.gene_parent_id,
        "start": protein.start,
        "end": protein.end,
    }


def extract_protein_identity_manifest(
    records: Sequence[SeqRecord],
    *,
    record_instance_keys: Sequence[str] | None = None,
    record_source_ids: Sequence[str] | None = None,
    record_selectors: Sequence[object] | None = None,
    regions: Sequence[object] | None = None,
    record_index_offset: int = 0,
    feature_visibility_rules: list[dict[str, object]] | None = None,
    feature_overrides: Sequence[Mapping[int, object]] | None = None,
) -> ProteinExtractionResult:
    """Build the canonical protein manifest and runtime-bound extraction."""

    base = extract_cds_proteins(
        records,
        record_index_offset=record_index_offset,
        prefer_source_ids=False,
        feature_visibility_rules=feature_visibility_rules,
        feature_overrides=feature_overrides,
    )
    instance_values = _per_record_values(
        record_instance_keys,
        records,
        name="record_instance_keys",
        default=lambda index, record: (
            (getattr(record, "annotations", {}) or {}).get("gbdraw_record_key")
            or f"record-{index + 1}"
        ),
    )
    source_values = _per_record_values(
        record_source_ids,
        records,
        name="record_source_ids",
        default=lambda _index, record: record.id,
    )
    selector_values = _per_record_values(
        record_selectors,
        records,
        name="record_selectors",
        default=lambda _index, _record: None,
    )
    region_values = _per_record_values(
        regions,
        records,
        name="regions",
        default=lambda _index, _record: None,
    )
    instance_keys = [_normalize_identity_text(value).strip() for value in instance_values]
    if any(not key for key in instance_keys):
        raise ValidationError("record_instance_keys must be non-empty.")
    if len(set(instance_keys)) != len(instance_keys):
        raise ValidationError("record_instance_keys must be unique across the extraction.")

    protein_sets: dict[str, dict[str, object]] = {}
    record_analyses: dict[str, dict[str, object]] = {}
    record_instances: dict[str, dict[str, object]] = {}
    proteins_by_record: list[list[CdsProtein]] = []
    protein_map: dict[str, CdsProtein] = {}
    protein_set_hashes: list[str] = []
    record_analysis_ids: list[str] = []
    runtime_binding_hashes: list[str] = []
    display_binding_hashes: list[str] = []

    for record_index, base_proteins in enumerate(base.proteins_by_record):
        grouped: dict[tuple[object, ...], list[CdsProtein]] = {}
        for protein in base_proteins:
            grouped.setdefault(_protein_location_identity_key(protein), []).append(protein)

        identified_by_position: dict[int, CdsProtein] = {}
        for group in grouped.values():
            for ordinal, protein in enumerate(sorted(group, key=_same_location_sort_key), start=1):
                aa_sha256 = hashlib.sha256(protein.sequence.encode("utf-8")).hexdigest()
                feature_analysis_id = canonical_feature_analysis_id(
                    feature_type=protein.feature_type,
                    location_operator=protein.location_operator,
                    location_parts=protein.feature_hash_parts,
                    strand=protein.feature_hash_strand,
                    same_location_ordinal=ordinal,
                )
                display_alias = select_protein_display_alias(protein)
                identified_by_position[id(protein)] = replace(
                    protein,
                    source_feature_position=(
                        protein.source_feature_position
                        if protein.source_feature_position is not None
                        else protein.feature_index
                    ),
                    same_location_ordinal=ordinal,
                    feature_analysis_id=feature_analysis_id,
                    display_alias=display_alias,
                    aa_sha256=aa_sha256,
                )

        identified = [identified_by_position[id(protein)] for protein in base_proteins]
        sorted_core = sorted(
            (_protein_identity_core(protein) for protein in identified),
            key=lambda entry: str(entry["featureAnalysisId"]),
        )
        protein_set_hash = _identity_sha256(
            _protein_set_hash_payload(
                [
                    (str(entry["featureAnalysisId"]), str(entry["aaSha256"]))
                    for entry in sorted_core
                ]
            )
        )
        protein_set = {"schema": 1, "proteins": sorted_core}
        existing_set = protein_sets.get(protein_set_hash)
        if (
            existing_set is not None
            and canonical_protein_identity_json(existing_set)
            != canonical_protein_identity_json(protein_set)
        ):
            raise ValidationError("Protein-set hash resolved to conflicting canonical content.")
        protein_sets.setdefault(protein_set_hash, protein_set)

        record_source_id = _normalize_identity_text(source_values[record_index]).strip()
        if not record_source_id:
            raise ValidationError("record_source_ids must be non-empty.")
        selector = _canonical_identity_context_value(selector_values[record_index])
        region = _canonical_identity_context_value(region_values[record_index])
        record_analysis_id = _identity_sha256(
            _record_analysis_hash_payload(
                protein_set_hash=protein_set_hash,
                record_source_id=record_source_id,
                region=region,
                selector=selector,
            )
        )
        record_analysis = {
            "schema": 1,
            "recordSourceId": record_source_id,
            "selector": selector,
            "region": region,
            "proteinSetHash": protein_set_hash,
        }
        existing_analysis = record_analyses.get(record_analysis_id)
        if (
            existing_analysis is not None
            and canonical_protein_identity_json(existing_analysis)
            != canonical_protein_identity_json(record_analysis)
        ):
            raise ValidationError("Record-analysis hash resolved to conflicting canonical content.")
        record_analyses.setdefault(record_analysis_id, record_analysis)

        instance_key = instance_keys[record_index]
        runtime_ids: dict[str, str] = {}
        for protein in sorted(
            identified,
            key=lambda item: str(item.feature_analysis_id),
        ):
            feature_id = str(protein.feature_analysis_id)
            runtime_ids[feature_id] = build_protein_runtime_handle(
                feature_analysis_id=feature_id,
                record_instance_key=instance_key,
            )

        export_ordinals = _expected_export_ordinals(
            {str(protein.feature_analysis_id): protein.display_alias for protein in identified}
        )
        feature_metadata = {
            str(protein.feature_analysis_id): _protein_binding_metadata(
                protein,
                export_ordinal=export_ordinals[str(protein.feature_analysis_id)],
            )
            for protein in identified
        }
        feature_metadata = dict(sorted(feature_metadata.items()))
        runtime_binding_hash = _identity_sha256(
            _runtime_binding_payload(
                protein_set_hash=protein_set_hash,
                record_instance_key=instance_key,
                runtime_ids=runtime_ids,
            )
        )
        display_binding_hash = _identity_sha256(
            _display_binding_payload(
                record_analysis_id=record_analysis_id,
                record_source_id=record_source_id,
                record_instance_key=instance_key,
                feature_metadata=feature_metadata,
            )
        )

        bound: list[CdsProtein] = []
        for protein in identified:
            feature_id = str(protein.feature_analysis_id)
            runtime_handle = runtime_ids[feature_id]
            bound_protein = replace(
                protein,
                protein_id=runtime_handle,
                runtime_handle=runtime_handle,
                record_instance_key=instance_key,
                record_analysis_id=record_analysis_id,
                protein_set_hash=protein_set_hash,
                runtime_binding_hash=runtime_binding_hash,
                display_binding_hash=display_binding_hash,
            )
            if runtime_handle in protein_map:
                raise ValidationError(
                    f"Duplicate protein runtime handle {runtime_handle!r}."
                )
            bound.append(bound_protein)
            protein_map[runtime_handle] = bound_protein

        binding = {
            "schema": 2,
            "recordAnalysisId": record_analysis_id,
            "runtimeBindingHash": runtime_binding_hash,
            "displayBindingHash": display_binding_hash,
            "runtimeIds": runtime_ids,
            "featureMetadata": feature_metadata,
        }
        record_instances[instance_key] = binding
        proteins_by_record.append(bound)
        protein_set_hashes.append(protein_set_hash)
        record_analysis_ids.append(record_analysis_id)
        runtime_binding_hashes.append(runtime_binding_hash)
        display_binding_hashes.append(display_binding_hash)

    manifest = {
        "schema": PROTEIN_IDENTITY_MANIFEST_SCHEMA,
        "proteinSets": protein_sets,
        "recordAnalyses": record_analyses,
        "recordInstances": record_instances,
    }
    # The validator's other checks re-derive what this loop computed with the same helpers;
    # only the key walk can fail on caller input (record keys, selectors, regions).
    _reject_manifest_presentation_identity(manifest)
    return ProteinExtractionResult(
        proteins_by_record=proteins_by_record,
        protein_map=protein_map,
        identity_manifest=_seal_protein_identity_manifest(manifest.items()),
        record_instance_keys=tuple(instance_keys),
        protein_set_hashes=tuple(protein_set_hashes),
        record_analysis_ids=tuple(record_analysis_ids),
        runtime_binding_hashes=tuple(runtime_binding_hashes),
        display_binding_hashes=tuple(display_binding_hashes),
    )


def build_web_losat_cache_key(
    *,
    query_fasta: str,
    subject_fasta: str,
    args: Sequence[str],
    program: str = "blastp",
    outfmt: str = "6",
    runtime_compatibility: str = "threaded-compatible-v1",
    threads_per_job: str | int = "auto",
) -> tuple[str, str, str]:
    """Rebuild a legacy schema-2 full-FASTA cache key for verification.

    Thread/runtime options are intentionally excluded: they control execution
    resources, not the raw LOSAT result identity.
    """

    query_hash = _hash_text_sha256(query_fasta)
    subject_hash = _hash_text_sha256(subject_fasta)
    payload = {
        "cacheSchema": LEGACY_PROTEIN_LOSAT_CACHE_SCHEMA,
        "program": program,
        "outfmt": str(outfmt),
        "args": [str(arg) for arg in args],
        "queryCanonicalHash": query_hash,
        "subjectCanonicalHash": subject_hash,
    }
    return _hash_text_sha256(_web_json_dumps(payload)), query_hash, subject_hash


def build_protein_losat_pair_identity(
    manifest: ProteinIdentityManifest | Mapping[str, object],
    *,
    query_record_instance_key: str,
    subject_record_instance_key: str,
) -> ProteinLosatPairIdentity:
    """Resolve one directional pair from a validated identity manifest."""

    authority = validate_protein_identity_manifest(manifest)

    def resolve(instance_key: str) -> tuple[str, str]:
        binding = authority.binding_for(instance_key)
        analysis_id = str(binding["recordAnalysisId"])
        analysis = authority.record_analyses[analysis_id]
        return str(analysis["proteinSetHash"]), str(binding["runtimeBindingHash"])

    query_set_hash, query_runtime_binding_hash = resolve(
        str(query_record_instance_key)
    )
    subject_set_hash, subject_runtime_binding_hash = resolve(
        str(subject_record_instance_key)
    )
    return ProteinLosatPairIdentity(
        query_protein_set_hash=query_set_hash,
        subject_protein_set_hash=subject_set_hash,
        query_runtime_binding_hash=query_runtime_binding_hash,
        subject_runtime_binding_hash=subject_runtime_binding_hash,
        query_record_instance_key=str(query_record_instance_key),
        subject_record_instance_key=str(subject_record_instance_key),
    )


def build_protein_losat_cache_key(
    pair_identity: ProteinLosatPairIdentity,
    *,
    args: Sequence[str],
    program: str = "blastp",
    outfmt: str = "6",
    search_context: str | None = None,
) -> str:
    """Return the Web-compatible directional schema-4 protein raw key."""

    payload = pair_identity.cache_payload(args=args, program=program, outfmt=outfmt)
    if search_context is not None:
        if not re.fullmatch(r"[0-9a-f]{64}", search_context):
            raise ValidationError("LOSAT search context must be a SHA-256 digest")
        payload["searchContext"] = search_context
    return _hash_text_sha256(_web_json_dumps(payload))


def parse_losat_fasta_ids(fasta_text: str) -> tuple[str, ...]:
    """Return FASTA identifier tokens and reject malformed/duplicate headers."""

    ids: list[str] = []
    seen: set[str] = set()
    for line_number, raw_line in enumerate(str(fasta_text).splitlines(), start=1):
        if not raw_line.startswith(">"):
            continue
        token = raw_line[1:].strip().split(maxsplit=1)[0] if raw_line[1:].strip() else ""
        if not token or not _VALID_FASTA_ID_RE.fullmatch(token):
            raise ValidationError(f"Invalid FASTA identifier on line {line_number}.")
        if token in seen:
            raise ValidationError(f"Duplicate FASTA identifier {token!r}.")
        seen.add(token)
        ids.append(token)
    if not ids and str(fasta_text).strip():
        raise ValidationError("Protein FASTA contains no identifier headers.")
    return tuple(ids)


def _parse_losat_fasta_sequences(fasta_text: str) -> dict[str, str]:
    ids = parse_losat_fasta_ids(fasta_text)
    if not ids:
        raise ValidationError("Protein FASTA contains no sequences.")
    sequences: dict[str, list[str]] = {protein_id: [] for protein_id in ids}
    current_id: str | None = None
    for line_number, raw_line in enumerate(str(fasta_text).splitlines(), start=1):
        if raw_line.startswith(">"):
            current_id = raw_line[1:].strip().split(maxsplit=1)[0]
            continue
        sequence_part = "".join(raw_line.split()).upper()
        if not sequence_part:
            continue
        if current_id is None:
            raise ValidationError(
                f"Protein FASTA sequence appears before its header on line {line_number}."
            )
        sequences[current_id].append(sequence_part)
    normalized = {
        protein_id: "".join(parts)
        for protein_id, parts in sequences.items()
    }
    if any(not sequence or not _VALID_PROTEIN_RE.fullmatch(sequence) for sequence in normalized.values()):
        raise ValidationError("Protein FASTA contains an empty or invalid amino-acid sequence.")
    return normalized


def _binding_runtime_ids(
    manifest: ProteinIdentityManifest,
    record_instance_key: str,
) -> set[str]:
    binding = manifest.binding_for(record_instance_key)
    runtime_ids = binding.get("runtimeIds")
    if not isinstance(runtime_ids, Mapping):
        raise ValidationError("Protein record binding has no runtime handles.")
    return {str(value) for value in runtime_ids.values()}


def _fasta_record_instance_key(
    fasta_text: str,
    manifest: ProteinIdentityManifest,
) -> str:
    fasta_sequences = _parse_losat_fasta_sequences(fasta_text)
    fasta_ids = set(fasta_sequences)
    matches: list[str] = []
    for instance_key, binding in manifest.record_instances.items():
        if fasta_ids != _binding_runtime_ids(manifest, instance_key):
            continue
        analysis = manifest.record_analyses[str(binding["recordAnalysisId"])]
        protein_set = manifest.protein_sets[str(analysis["proteinSetHash"])]
        aa_by_feature_id = {
            str(protein["featureAnalysisId"]): str(protein["aaSha256"])
            for protein in cast("Sequence[Mapping[str, object]]", protein_set["proteins"])
        }
        runtime_ids = binding["runtimeIds"]
        if not isinstance(runtime_ids, Mapping):
            continue
        if all(
            hashlib.sha256(fasta_sequences[str(runtime_handle)].encode("utf-8")).hexdigest()
            == aa_by_feature_id[str(feature_id)]
            for feature_id, runtime_handle in runtime_ids.items()
        ):
            matches.append(instance_key)
    if len(matches) != 1:
        raise ValidationError(
            "Protein FASTA IDs and sequences must match exactly one complete record-instance binding."
        )
    return matches[0]


# The verdict depends only on the text and the runtime handles of the two
# bindings, which a validated manifest's runtime binding hash digests. Key by the
# text's digest and both hashes, so the memo does not keep the (possibly
# multi-megabyte) texts alive and an entry that several layers check reads its
# rows once; bound it.
_RAW_TSV_VERDICTS: dict[tuple[bytes, str, str], bool] = {}
_RAW_TSV_VERDICTS_MAX = 1024


def _raw_tsv_verdict(text: str, query_ids: set[str], subject_ids: set[str]) -> bool:
    for raw_line in text.splitlines():
        if not raw_line.strip() or raw_line.lstrip().startswith("#"):
            continue
        columns = raw_line.split("\t")
        if (
            len(columns) != len(COMPARISON_COLUMNS)
            or columns[0] not in query_ids
            or columns[1] not in subject_ids
        ):
            return False
    try:
        parse_losatp_outfmt6(text)
    except ParseError:
        return False
    return True


def raw_protein_tsv_matches_bindings(
    text: str,
    manifest: ProteinIdentityManifest,
    *,
    query_record_instance_key: str,
    subject_record_instance_key: str,
) -> bool:
    """Return whether the TSV parses and every non-comment row belongs to the directional pair."""

    text = str(text)
    key = (
        hashlib.sha256(text.encode("utf-8", "surrogatepass")).digest(),
        str(manifest.binding_for(query_record_instance_key)["runtimeBindingHash"]),
        str(manifest.binding_for(subject_record_instance_key)["runtimeBindingHash"]),
    )
    verdict = _RAW_TSV_VERDICTS.get(key)
    if verdict is None:
        verdict = _raw_tsv_verdict(
            text,
            _binding_runtime_ids(manifest, query_record_instance_key),
            _binding_runtime_ids(manifest, subject_record_instance_key),
        )
        if len(_RAW_TSV_VERDICTS) >= _RAW_TSV_VERDICTS_MAX:
            _RAW_TSV_VERDICTS.clear()
        _RAW_TSV_VERDICTS[key] = verdict
    return verdict


def is_protein_losat_cache_entry(entry: Mapping[str, object] | object) -> TypeGuard[Mapping[str, object]]:
    """Discriminate a current schema-4 protein raw entry."""

    if not isinstance(entry, Mapping):
        return False
    required_text = (
        "key",
        "queryProteinSetHash",
        "subjectProteinSetHash",
        "queryRuntimeBindingHash",
        "subjectRuntimeBindingHash",
        "queryRecordInstanceKey",
        "subjectRecordInstanceKey",
    )
    return (
        entry.get("schema") == PROTEIN_LOSAT_CACHE_SCHEMA
        and entry.get("kind") == "raw-losat"
        and entry.get("identityKind") == "protein"
        and entry.get("idEncoding") == "runtime-handle-v1"
        and str(entry.get("program") or "").lower() == "blastp"
        and isinstance(entry.get("text"), str)
        and all(isinstance(entry.get(key), str) and bool(entry.get(key)) for key in required_text)
    )


def is_legacy_protein_losat_cache_entry(
    entry: Mapping[str, object] | object,
) -> TypeGuard[Mapping[str, object]]:
    """Discriminate an imported schema-2 blastp candidate."""

    return (
        isinstance(entry, Mapping)
        and entry.get("schema") == LEGACY_PROTEIN_LOSAT_CACHE_SCHEMA
        and entry.get("kind") == "raw-losat"
        and isinstance(entry.get("text"), str)
        and str(entry.get("program") or "blastp").lower() == "blastp"
        and str(entry.get("flow") or "").lower() not in {"blastn", "tblastx", "nucleotide"}
    )


def validate_protein_raw_entry_references(
    entry: Mapping[str, object] | object,
    manifest: ProteinIdentityManifest | Mapping[str, object],
) -> bool:
    """Validate a schema-4 key, manifest references, and directional TSV IDs."""

    if not is_protein_losat_cache_entry(entry):
        return False
    try:
        authority = validate_protein_identity_manifest(manifest)
        pair_identity = build_protein_losat_pair_identity(
            authority,
            query_record_instance_key=str(entry["queryRecordInstanceKey"]),
            subject_record_instance_key=str(entry["subjectRecordInstanceKey"]),
        )
        if (
            str(entry["queryProteinSetHash"]) != pair_identity.query_protein_set_hash
            or str(entry["subjectProteinSetHash"]) != pair_identity.subject_protein_set_hash
            or str(entry["queryRuntimeBindingHash"])
            != pair_identity.query_runtime_binding_hash
            or str(entry["subjectRuntimeBindingHash"])
            != pair_identity.subject_runtime_binding_hash
        ):
            return False
        args = entry.get("args")
        if not isinstance(args, list):
            return False
        expected_key = build_protein_losat_cache_key(
            pair_identity,
            args=args,
            program=str(entry.get("program") or "blastp"),
            outfmt=str(entry.get("outfmt") or "6"),
            # A non-string context raises TypeError in the key builder, handled below.
            search_context=cast("str | None", entry.get("searchContext")),
        )
        if str(entry.get("key") or "") != expected_key:
            return False
        return raw_protein_tsv_matches_bindings(
            str(entry.get("text") or ""),
            authority,
            query_record_instance_key=pair_identity.query_record_instance_key,
            subject_record_instance_key=pair_identity.subject_record_instance_key,
        )
    except (KeyError, TypeError, ValueError, ValidationError):
        return False


def build_protein_export_id_map(
    manifest: ProteinIdentityManifest | Mapping[str, object],
    record_instance_key: str,
) -> dict[str, str]:
    """Resolve runtime handles to deterministic user-visible export IDs."""

    authority = validate_protein_identity_manifest(manifest)
    binding = authority.binding_for(str(record_instance_key))
    runtime_ids = binding.get("runtimeIds")
    feature_metadata = binding.get("featureMetadata")
    if not isinstance(runtime_ids, Mapping) or not isinstance(feature_metadata, Mapping):
        raise ValidationError("Protein record binding cannot resolve export IDs.")

    result: dict[str, str] = {}
    for feature_id, runtime_handle in runtime_ids.items():
        metadata = feature_metadata.get(feature_id)
        if not isinstance(metadata, Mapping):
            raise ValidationError("Protein export metadata is incomplete.")
        alias = _normalize_identity_text(metadata.get("displayAlias")).strip()
        if not alias:
            raise ValidationError("Protein export alias is empty.")
        export_id = percent_encode_losat_transport_field(alias)
        ordinal = metadata.get("exportOrdinal")
        if ordinal is not None:
            if isinstance(ordinal, bool) or not isinstance(ordinal, int) or ordinal <= 0:
                raise ValidationError("Protein export ordinal is invalid.")
            export_id = f"{export_id}~{ordinal}"
        handle = str(runtime_handle)
        if handle in result:
            raise ValidationError("Protein runtime handle has multiple export IDs.")
        result[handle] = export_id
    if len(set(result.values())) != len(result):
        raise ValidationError("Protein export IDs are not unique within the record instance.")
    return result


def hydrate_protein_losat_tsv(
    entry: Mapping[str, object],
    manifest: ProteinIdentityManifest | Mapping[str, object],
) -> str:
    """Replace only QUERY/SUBJECT runtime handles for a raw TSV download."""

    authority = validate_protein_identity_manifest(manifest)
    if not validate_protein_raw_entry_references(entry, authority):
        raise ValidationError(
            "Protein raw LOSAT entry does not resolve through the current manifest."
        )
    query_ids = build_protein_export_id_map(
        authority,
        str(entry["queryRecordInstanceKey"]),
    )
    subject_ids = build_protein_export_id_map(
        authority,
        str(entry["subjectRecordInstanceKey"]),
    )
    hydrated: list[str] = []
    for raw_line in str(entry.get("text") or "").splitlines(keepends=True):
        body = raw_line.rstrip("\r\n")
        ending = raw_line[len(body) :]
        if not body.strip() or body.lstrip().startswith("#"):
            hydrated.append(raw_line)
            continue
        columns = body.split("\t")
        if len(columns) != len(COMPARISON_COLUMNS):
            raise ValidationError("Protein LOSAT TSV must contain exactly 12 columns.")
        try:
            columns[0] = query_ids[columns[0]]
            columns[1] = subject_ids[columns[1]]
        except KeyError as exc:
            raise ValidationError(
                "Protein LOSAT TSV contains an unresolved or wrong-binding runtime handle."
            ) from exc
        hydrated.append("\t".join(columns) + ending)
    result = "".join(hydrated)
    if _RUNTIME_HANDLE_SEARCH_RE.search(result):
        raise ValidationError("Protein LOSAT export still contains an internal runtime handle.")
    return result


def make_legacy_protein_raw_candidate(
    entry: Mapping[str, object],
    *,
    state: Literal["pending", "promoted", "rejected"] = "pending",
    rejection_reason: str | None = None,
) -> dict[str, object]:
    """Wrap an imported legacy entry without modifying its original payload."""

    if not is_legacy_protein_losat_cache_entry(entry):
        raise ValidationError("Legacy protein candidate must contain a schema-2 blastp entry.")
    if state not in {"pending", "promoted", "rejected"}:
        raise ValidationError(f"Unsupported legacy protein candidate state {state!r}.")
    return {
        "state": state,
        "originalEntry": copy.deepcopy(dict(entry)),
        "rejectionReason": str(rejection_reason) if rejection_reason is not None else None,
    }


def validate_legacy_protein_raw_candidate_envelope(
    value: Mapping[str, object],
) -> dict[str, object]:
    """Validate the lossless schema-1 legacy candidate envelope."""

    if value.get("schema") != LEGACY_PROTEIN_RAW_CANDIDATE_SCHEMA or not isinstance(value.get("entries"), list):
        raise ValidationError("Legacy protein raw candidate envelope must use schema 1.")
    entries: list[dict[str, object]] = []
    for candidate in cast("Sequence[object]", value["entries"]):
        if not isinstance(candidate, Mapping):
            raise ValidationError("Legacy protein candidate must be an object.")
        state = str(candidate.get("state") or "")
        original = candidate.get("originalEntry")
        if state not in {"pending", "promoted", "rejected"} or not isinstance(original, Mapping):
            raise ValidationError("Legacy protein candidate state or original entry is invalid.")
        if not is_legacy_protein_losat_cache_entry(original):
            raise ValidationError("Legacy protein candidate contains a non-protein entry.")
        entries.append(
            make_legacy_protein_raw_candidate(
                original,
                state=cast(Literal["pending", "promoted", "rejected"], state),
                rejection_reason=(
                    str(candidate.get("rejectionReason"))
                    if candidate.get("rejectionReason") is not None
                    else None
                ),
            )
        )
    return {"schema": LEGACY_PROTEIN_RAW_CANDIDATE_SCHEMA, "entries": entries}


_LEGACY_WEB_PROTEIN_ID_RE = re.compile(
    r"^p_(?P<record_token>.+)_(?P<start>\d+)_(?P<end>\d+)_"
    r"(?P<strand>-1|0|1)_(?P<aa_sha12>[0-9a-f]{12})(?:_(?P<duplicate>[2-9][0-9]*))?$"
)


def _legacy_record_tokens_from_tsv(text: str) -> tuple[str, str]:
    query_tokens: set[str] = set()
    subject_tokens: set[str] = set()
    found_row = False
    for raw_line in str(text).splitlines():
        if not raw_line.strip() or raw_line.lstrip().startswith("#"):
            continue
        columns = raw_line.split("\t")
        if len(columns) != len(COMPARISON_COLUMNS):
            raise ValidationError("Legacy LOSAT TSV does not contain 12 outfmt-6 columns.")
        query_match = _LEGACY_WEB_PROTEIN_ID_RE.fullmatch(columns[0])
        subject_match = _LEGACY_WEB_PROTEIN_ID_RE.fullmatch(columns[1])
        if query_match is None or subject_match is None:
            raise ValidationError("Legacy LOSAT TSV contains an unrecognized protein ID.")
        query_tokens.add(query_match.group("record_token"))
        subject_tokens.add(subject_match.group("record_token"))
        found_row = True
    if not found_row:
        raise ValidationError("Empty legacy LOSAT output does not prove record tokens.")
    if len(query_tokens) != 1 or len(subject_tokens) != 1:
        raise ValidationError("Legacy LOSAT TSV does not resolve to one directional record pair.")
    return next(iter(query_tokens)), next(iter(subject_tokens))


def _legacy_record_token_evidence(
    entries: Sequence[Mapping[str, object]],
) -> dict[str, set[str]]:
    tokens_by_fasta_hash: dict[str, set[str]] = {}
    for candidate in entries:
        entry = _legacy_original_entry(candidate)
        if not is_legacy_protein_losat_cache_entry(entry):
            continue
        try:
            query_token, subject_token = _legacy_record_tokens_from_tsv(
                str(entry.get("text") or "")
            )
        except ValidationError:
            continue
        for hash_field, token in (
            ("queryCanonicalHash", query_token),
            ("subjectCanonicalHash", subject_token),
        ):
            fasta_hash = str(entry.get(hash_field) or "")
            if fasta_hash:
                tokens_by_fasta_hash.setdefault(fasta_hash, set()).add(token)
    return tokens_by_fasta_hash


def _legacy_fasta_and_id_map(
    proteins: Sequence[CdsProtein],
    record_token: str,
) -> tuple[str, dict[str, str]]:
    legacy = _with_stable_web_protein_ids(proteins, record_token)
    mapping: dict[str, str] = {}
    for legacy_protein, current_protein in zip(legacy, proteins, strict=True):
        if legacy_protein.protein_id in mapping:
            raise ValidationError("Legacy protein IDs are not one-to-one.")
        mapping[legacy_protein.protein_id] = current_protein.protein_id
    if len(set(mapping.values())) != len(mapping):
        raise ValidationError("Current runtime handles are not one-to-one.")
    return proteins_to_fasta(legacy), mapping


def _legacy_reverse_record_lengths(
    proteins: Sequence[CdsProtein],
    references: set[str],
) -> set[int]:
    """Infer lengths that make references exact reverse-complement spans."""

    lengths: set[int] = set()
    for reference in references:
        match = _LEGACY_WEB_PROTEIN_ID_RE.fullmatch(reference)
        if match is None:  # Already validated by the caller.
            continue
        reference_start = int(match.group("start"))
        reference_end = int(match.group("end"))
        reference_strand = int(match.group("strand"))
        reference_aa_sha12 = match.group("aa_sha12")
        for protein in proteins:
            protein_start = (
                int(protein.feature_hash_start)
                if protein.feature_hash_start is not None
                else int(protein.start)
            )
            protein_end = (
                int(protein.feature_hash_end)
                if protein.feature_hash_end is not None
                else int(protein.end)
            )
            identity_strand = (
                protein.feature_hash_strand
                if protein.feature_hash_strand is not None
                else protein.strand
            )
            protein_strand = identity_strand if identity_strand in {-1, 1} else 0
            reversed_strand = -protein_strand if protein_strand else 0
            if (
                reference_aa_sha12
                != hashlib.sha256(protein.sequence.encode("utf-8")).hexdigest()[:12]
                or reference_end - reference_start != protein_end - protein_start
                or reference_strand != reversed_strand
            ):
                continue
            length_from_start = reference_start + protein_end
            length_from_end = reference_end + protein_start
            if (
                length_from_start == length_from_end
                and length_from_start >= max(reference_end, protein_end)
            ):
                lengths.add(length_from_start)
    return lengths


def _legacy_reverse_id_map(
    proteins: Sequence[CdsProtein],
    record_token: str,
    record_length: int,
) -> dict[str, str]:
    """Reconstruct IDs emitted for a legacy reverse-complemented record."""

    reversed_proteins: list[CdsProtein] = []
    for protein in reversed(proteins):
        protein_start = (
            int(protein.feature_hash_start)
            if protein.feature_hash_start is not None
            else int(protein.start)
        )
        protein_end = (
            int(protein.feature_hash_end)
            if protein.feature_hash_end is not None
            else int(protein.end)
        )
        start = int(record_length) - protein_end
        end = int(record_length) - protein_start
        if start < 0 or end <= start:
            return {}
        identity_strand = (
            protein.feature_hash_strand
            if protein.feature_hash_strand is not None
            else protein.strand
        )
        protein_strand = identity_strand if identity_strand in {-1, 1} else 0
        reversed_proteins.append(
            replace(
                protein,
                start=start,
                end=end,
                strand=-protein_strand if protein_strand else None,
                feature_hash_start=start,
                feature_hash_end=end,
                feature_hash_strand=-protein_strand if protein_strand else None,
            )
        )
    _, mapping = _legacy_fasta_and_id_map(reversed_proteins, record_token)
    return mapping


def build_legacy_protein_reference_map(
    extraction: ProteinExtractionResult,
    reference_ids: Sequence[str],
) -> dict[str, str]:
    """Resolve legacy ``p_r_`` artifact IDs through current protein identity."""

    if extraction.identity_manifest is None:
        raise ValidationError(
            "Legacy protein reference migration requires an identity manifest."
        )
    references_by_token: dict[str, set[str]] = {}
    for raw_reference in reference_ids:
        reference = str(raw_reference)
        match = _LEGACY_WEB_PROTEIN_ID_RE.fullmatch(reference)
        if match is None:
            raise ValidationError(
                f"Unrecognized legacy protein reference {reference!r}."
            )
        references_by_token.setdefault(
            match.group("record_token"), set()
        ).add(reference)

    resolved: dict[str, str] = {}
    for record_token, references in references_by_token.items():
        candidate_maps: list[dict[str, str]] = []
        candidate_signatures: set[tuple[tuple[str, str], ...]] = set()
        for proteins in extraction.proteins_by_record:
            _, direct_map = _legacy_fasta_and_id_map(proteins, record_token)
            unresolved = references.difference(direct_map)
            combined_maps = [direct_map] if not unresolved else []
            for record_length in _legacy_reverse_record_lengths(
                proteins,
                unresolved,
            ):
                reverse_map = _legacy_reverse_id_map(
                    proteins,
                    record_token,
                    record_length,
                )
                if any(
                    reference in direct_map
                    and direct_map[reference] != runtime_handle
                    for reference, runtime_handle in reverse_map.items()
                ):
                    continue
                combined_maps.append({**direct_map, **reverse_map})
            for candidate_map in combined_maps:
                if not references.issubset(candidate_map):
                    continue
                projected_map = {
                    reference: candidate_map[reference]
                    for reference in references
                }
                signature = tuple(sorted(projected_map.items()))
                if signature in candidate_signatures:
                    continue
                candidate_signatures.add(signature)
                candidate_maps.append(projected_map)
        if len(candidate_maps) != 1:
            raise ValidationError(
                f"Legacy protein record token {record_token!r} resolves to "
                f"{len(candidate_maps)} current record instances."
            )
        candidate_map = candidate_maps[0]
        for reference in references:
            runtime_handle = candidate_map[reference]
            previous = resolved.get(reference)
            if previous is not None and previous != runtime_handle:
                raise ValidationError(
                    f"Legacy protein reference {reference!r} is ambiguous."
                )
            resolved[reference] = runtime_handle
    return resolved


def _rewrite_legacy_protein_tsv(
    text: str,
    query_id_map: Mapping[str, str],
    subject_id_map: Mapping[str, str],
) -> str:
    rewritten: list[str] = []
    for raw_line in str(text).splitlines(keepends=True):
        body = raw_line.rstrip("\r\n")
        ending = raw_line[len(body) :]
        if not body.strip() or body.lstrip().startswith("#"):
            rewritten.append(raw_line)
            continue
        columns = body.split("\t")
        if len(columns) != len(COMPARISON_COLUMNS):
            raise ValidationError("Legacy LOSAT TSV does not contain 12 outfmt-6 columns.")
        if columns[0] not in query_id_map or columns[1] not in subject_id_map:
            raise ValidationError("Legacy LOSAT TSV references a protein outside the verified pair.")
        columns[0] = query_id_map[columns[0]]
        columns[1] = subject_id_map[columns[1]]
        rewritten.append("\t".join(columns) + ending)
    result = "".join(rewritten)
    # Parsing also proves that all numeric fields survived the rewrite.
    parse_losatp_outfmt6(result)
    return result


def _legacy_original_entry(candidate: Mapping[str, object]) -> Mapping[str, object]:
    original = candidate.get("originalEntry")
    return original if isinstance(original, Mapping) else candidate


def promote_legacy_protein_raw_cache_entries(
    entries: Sequence[Mapping[str, object]],
    *,
    query_proteins: Sequence[CdsProtein],
    subject_proteins: Sequence[CdsProtein],
    query_fasta: str,
    subject_fasta: str,
    identity_manifest: ProteinIdentityManifest | Mapping[str, object],
    expected_args: Sequence[str],
    expected_program: str = "blastp",
    expected_outfmt: str = "6",
) -> LegacyProteinRawCachePromotionScan:
    """Find and verify a legacy pair, then return a schema-4 copy.

    No input object is mutated. A failed candidate contributes a pair-local
    rejection diagnostic and never enters the current cache.
    """

    manifest = validate_protein_identity_manifest(identity_manifest)
    if proteins_to_fasta(query_proteins) != str(query_fasta):
        raise ValidationError("query_fasta does not match the current query extraction.")
    if proteins_to_fasta(subject_proteins) != str(subject_fasta):
        raise ValidationError("subject_fasta does not match the current subject extraction.")
    query_instance_key = _fasta_record_instance_key(query_fasta, manifest)
    subject_instance_key = _fasta_record_instance_key(subject_fasta, manifest)
    pair_identity = build_protein_losat_pair_identity(
        manifest,
        query_record_instance_key=query_instance_key,
        subject_record_instance_key=subject_instance_key,
    )
    expected_arg_tuple = tuple(str(arg) for arg in expected_args)
    rejections: list[LegacyProteinRawCacheRejection] = []
    token_evidence: dict[str, set[str]] | None = None

    for candidate_index, candidate in enumerate(entries):
        entry = _legacy_original_entry(candidate)
        try:
            if not is_legacy_protein_losat_cache_entry(entry):
                raise ValidationError("Candidate is not a schema-2 protein raw entry.")
            if str(entry.get("program") or "blastp").lower() != str(expected_program).lower():
                raise ValidationError("Legacy program does not match the requested program.")
            if str(entry.get("outfmt") or "") != str(expected_outfmt):
                raise ValidationError("Legacy outfmt is missing or does not match.")
            raw_args = entry.get("args")
            if not isinstance(raw_args, list) or tuple(str(arg) for arg in raw_args) != expected_arg_tuple:
                raise ValidationError("Legacy search args are missing or do not match.")

            raw_text = str(entry.get("text") or "")
            if any(
                line.strip() and not line.lstrip().startswith("#")
                for line in raw_text.splitlines()
            ):
                query_token, subject_token = _legacy_record_tokens_from_tsv(raw_text)
            else:
                if token_evidence is None:
                    token_evidence = _legacy_record_token_evidence(entries)
                query_candidates = token_evidence.get(
                    str(entry.get("queryCanonicalHash") or ""),
                    set(),
                )
                subject_candidates = token_evidence.get(
                    str(entry.get("subjectCanonicalHash") or ""),
                    set(),
                )
                if len(query_candidates) != 1 or len(subject_candidates) != 1:
                    raise ValidationError(
                        "Empty legacy LOSAT output does not prove record tokens."
                    )
                query_token = next(iter(query_candidates))
                subject_token = next(iter(subject_candidates))
            legacy_query_fasta, query_id_map = _legacy_fasta_and_id_map(
                query_proteins,
                query_token,
            )
            legacy_subject_fasta, subject_id_map = _legacy_fasta_and_id_map(
                subject_proteins,
                subject_token,
            )
            if _hash_text_sha256(legacy_query_fasta) != str(entry.get("queryCanonicalHash") or ""):
                raise ValidationError("Legacy query FASTA hash does not match.")
            if _hash_text_sha256(legacy_subject_fasta) != str(entry.get("subjectCanonicalHash") or ""):
                raise ValidationError("Legacy subject FASTA hash does not match.")
            rewritten_tsv = _rewrite_legacy_protein_tsv(
                raw_text,
                query_id_map,
                subject_id_map,
            )
            if not raw_protein_tsv_matches_bindings(
                rewritten_tsv,
                manifest,
                query_record_instance_key=query_instance_key,
                subject_record_instance_key=subject_instance_key,
            ):
                raise ValidationError("Rewritten TSV does not match the current bindings.")

            cache_key = build_protein_losat_cache_key(
                pair_identity,
                args=expected_arg_tuple,
                program=expected_program,
                outfmt=expected_outfmt,
            )
            promoted: dict[str, object] = {
                "schema": PROTEIN_LOSAT_CACHE_SCHEMA,
                "kind": "raw-losat",
                "identityKind": "protein",
                "idEncoding": "runtime-handle-v1",
                "key": cache_key,
                "filename": losat_entry_filename(entry),
                "display": bool(entry.get("display")) if "display" in entry else False,
                "text": rewritten_tsv,
                "program": str(expected_program),
                "outfmt": str(expected_outfmt),
                "args": list(expected_arg_tuple),
                "queryProteinSetHash": pair_identity.query_protein_set_hash,
                "subjectProteinSetHash": pair_identity.subject_protein_set_hash,
                "queryRuntimeBindingHash": pair_identity.query_runtime_binding_hash,
                "subjectRuntimeBindingHash": pair_identity.subject_runtime_binding_hash,
                "queryRecordInstanceKey": pair_identity.query_record_instance_key,
                "subjectRecordInstanceKey": pair_identity.subject_record_instance_key,
            }
            if not is_protein_losat_cache_entry(promoted):
                raise ValidationError("Promoted schema-4 entry failed validation.")
            return LegacyProteinRawCachePromotionScan(
                promotion=LegacyProteinRawCachePromotion(
                    candidate_index=candidate_index,
                    entry=promoted,
                    rewritten_tsv=rewritten_tsv,
                    legacy_query_fasta=legacy_query_fasta,
                    legacy_subject_fasta=legacy_subject_fasta,
                ),
                rejections=tuple(rejections),
            )
        except (ParseError, ValidationError) as exc:
            rejections.append(
                LegacyProteinRawCacheRejection(
                    candidate_index=candidate_index,
                    reason=str(exc),
                )
            )

    return LegacyProteinRawCachePromotionScan(
        promotion=None,
        rejections=tuple(rejections),
    )


@dataclass(frozen=True)
class LosatpCacheState:
    """What a protein raw cache holds; the protein analysis is keyed by it."""

    entries: tuple[tuple[str, Mapping[str, object]], ...]
    display: tuple[tuple[str, str], ...]
    legacy_entries: tuple[Mapping[str, object], ...]
    identity_manifest: ProteinIdentityManifest | None
    proteins_by_instance: tuple[tuple[str, tuple[CdsProtein, ...]], ...]
    threads_per_job: str | int
    runtime_compatibility: str


class LosatpCacheManager(LosatRawCache):
    """Validated schema-4 protein raw-cache lookup and collection."""

    def __init__(
        self,
        entries: Sequence[Mapping[str, object]] | None = None,
        *,
        identity_manifest: ProteinIdentityManifest | Mapping[str, object] | None = None,
        threads_per_job: str | int = "auto",
        runtime_compatibility: str = "threaded-compatible-v1",
    ) -> None:
        super().__init__()
        self.threads_per_job = threads_per_job
        self.runtime_compatibility = runtime_compatibility
        self.identity_manifest = (
            validate_protein_identity_manifest(identity_manifest)
            if identity_manifest is not None
            else None
        )
        self._legacy_entries: list[dict[str, object]] = []
        self._proteins_by_instance: dict[str, tuple[CdsProtein, ...]] = {}
        for entry in entries or ():
            if is_legacy_protein_losat_cache_entry(entry):
                self._legacy_entries.append(copy.deepcopy(dict(entry)))
                continue
            if not is_protein_losat_cache_entry(entry):
                continue
            key = str(entry.get("key") or "").strip()
            if not key:
                continue
            normalized = {
                "schema": PROTEIN_LOSAT_CACHE_SCHEMA,
                "kind": "raw-losat",
                "identityKind": "protein",
                "idEncoding": "runtime-handle-v1",
                "key": key,
                "filename": losat_entry_filename(entry),
                "display": bool(entry.get("display")) if "display" in entry else False,
                "text": str(entry.get("text") or ""),
                "program": str(entry.get("program") or "blastp"),
                "outfmt": str(entry.get("outfmt") or "6"),
                "args": [str(arg) for arg in cast("Sequence[object]", entry.get("args") or ())],
                "queryProteinSetHash": str(entry.get("queryProteinSetHash") or ""),
                "subjectProteinSetHash": str(entry.get("subjectProteinSetHash") or ""),
                "queryRuntimeBindingHash": str(
                    entry.get("queryRuntimeBindingHash") or ""
                ),
                "subjectRuntimeBindingHash": str(
                    entry.get("subjectRuntimeBindingHash") or ""
                ),
                "queryRecordInstanceKey": str(entry.get("queryRecordInstanceKey") or ""),
                "subjectRecordInstanceKey": str(entry.get("subjectRecordInstanceKey") or ""),
            }
            if "searchContext" in entry:
                normalized["searchContext"] = entry["searchContext"]
            if self.identity_manifest is not None and not validate_protein_raw_entry_references(
                normalized,
                self.identity_manifest,
            ):
                continue
            self._add_loaded_entry(key, normalized, entry)

    def state(self) -> LosatpCacheState:
        entries, display = self._raw_state()
        return LosatpCacheState(
            entries=entries,
            display=display,
            legacy_entries=tuple(self._legacy_entries),
            identity_manifest=self.identity_manifest,
            proteins_by_instance=tuple(self._proteins_by_instance.items()),
            threads_per_job=self.threads_per_job,
            runtime_compatibility=self.runtime_compatibility,
        )

    @classmethod
    def from_state(cls, state: LosatpCacheState) -> LosatpCacheManager:
        """A cache in ``state`` whose writes do not reach the cache it came from."""

        cache = cls(
            threads_per_job=state.threads_per_job,
            runtime_compatibility=state.runtime_compatibility,
        )
        cache.identity_manifest = state.identity_manifest
        cache._restore_raw_state(state.entries, state.display)
        cache._legacy_entries = [copy.deepcopy(dict(entry)) for entry in state.legacy_entries]
        cache._proteins_by_instance = dict(state.proteins_by_instance)
        return cache

    def _apply_write(self, write: LosatCacheWrite) -> None:
        if write.kind == "retire-legacy":
            self._retire_legacy(write.index)
        else:
            super()._apply_write(write)

    def _retire_legacy(self, index: int) -> None:
        del self._legacy_entries[index]
        self._writes.append(LosatCacheWrite("retire-legacy", index=index))

    @property
    def has_legacy_candidates(self) -> bool:
        return bool(self._legacy_entries)

    def set_identity_manifest(
        self,
        manifest: ProteinIdentityManifest | Mapping[str, object],
    ) -> None:
        self.identity_manifest = validate_protein_identity_manifest(manifest)
        authority = self.identity_manifest
        self._retain_entries(
            lambda entry: validate_protein_raw_entry_references(entry, authority)
        )

    def set_protein_extraction(self, extraction: ProteinExtractionResult) -> None:
        """Attach the current extraction needed to verify legacy cache candidates."""

        if extraction.identity_manifest is None:
            raise ValidationError(
                "Stable protein extraction must include an identity manifest."
            )
        if len(extraction.record_instance_keys) != len(extraction.proteins_by_record):
            raise ValidationError(
                "Stable protein extraction record bindings are incomplete."
            )

        self.set_identity_manifest(extraction.identity_manifest)
        proteins_by_instance: dict[str, tuple[CdsProtein, ...]] = {}
        for instance_key, proteins in zip(
            extraction.record_instance_keys,
            extraction.proteins_by_record,
            strict=True,
        ):
            normalized_key = str(instance_key).strip()
            if not normalized_key or normalized_key in proteins_by_instance:
                raise ValidationError(
                    "Stable protein extraction record-instance keys must be non-empty and unique."
                )
            protein_tuple = tuple(proteins)
            fasta_text = proteins_to_fasta(protein_tuple)
            # set_identity_manifest above stored the (non-None) extraction manifest.
            if (
                _fasta_record_instance_key(
                    fasta_text, cast(ProteinIdentityManifest, self.identity_manifest)
                )
                != normalized_key
            ):
                raise ValidationError(
                    "Stable protein extraction does not match its record-instance binding."
                )
            proteins_by_instance[normalized_key] = protein_tuple
        self._proteins_by_instance = proteins_by_instance

    def legacy_candidate_envelope(self) -> dict[str, object]:
        return {
            "schema": LEGACY_PROTEIN_RAW_CANDIDATE_SCHEMA,
            "entries": [make_legacy_protein_raw_candidate(entry) for entry in self._legacy_entries],
        }

    def add_promoted_entry(self, entry: Mapping[str, object]) -> None:
        if not is_protein_losat_cache_entry(entry):
            raise ValidationError(
                "Only current schema-4 protein entries can be promoted into the cache."
            )
        if self.identity_manifest is None or not validate_protein_raw_entry_references(
            entry,
            self.identity_manifest,
        ):
            raise ValidationError("Promoted protein entry does not match the current manifest.")
        self._add_search_entry(str(entry["key"]), dict(entry))

    def _pair_identity_from_fasta(
        self,
        query_fasta: str,
        subject_fasta: str,
    ) -> ProteinLosatPairIdentity:
        if self.identity_manifest is None:
            raise ValidationError(
                "Protein raw cache schema 4 requires a protein identity manifest."
            )
        query_instance_key = _fasta_record_instance_key(query_fasta, self.identity_manifest)
        subject_instance_key = _fasta_record_instance_key(subject_fasta, self.identity_manifest)
        return build_protein_losat_pair_identity(
            self.identity_manifest,
            query_record_instance_key=query_instance_key,
            subject_record_instance_key=subject_instance_key,
        )

    def _find_cached_entry(
        self,
        *,
        cache_key: str,
        pair_identity: ProteinLosatPairIdentity,
        program: str,
        outfmt: str,
        args: Sequence[str],
    ) -> dict[str, object] | None:
        cached = self._cached_entry(cache_key)
        if cached is None or self.identity_manifest is None:
            return None
        expected = pair_identity.cache_payload(args=args, program=program, outfmt=outfmt)
        for payload_key, entry_key in (
            ("queryProteinSetHash", "queryProteinSetHash"),
            ("subjectProteinSetHash", "subjectProteinSetHash"),
            ("queryRuntimeBindingHash", "queryRuntimeBindingHash"),
            ("subjectRuntimeBindingHash", "subjectRuntimeBindingHash"),
            ("queryRecordInstanceKey", "queryRecordInstanceKey"),
            ("subjectRecordInstanceKey", "subjectRecordInstanceKey"),
        ):
            if str(cached.get(entry_key) or "") != str(expected[payload_key]):
                return None
        if str(cached.get("program") or "") != program:
            return None
        if str(cached.get("outfmt") or "") != outfmt:
            return None
        if tuple(str(arg) for arg in cast("Sequence[object]", cached.get("args") or ())) != tuple(str(arg) for arg in args):
            return None
        if not raw_protein_tsv_matches_bindings(
            str(cached.get("text") or ""),
            self.identity_manifest,
            query_record_instance_key=pair_identity.query_record_instance_key,
            subject_record_instance_key=pair_identity.subject_record_instance_key,
        ):
            return None
        return cached

    def cached_text(
        self,
        query_fasta: str,
        subject_fasta: str,
        *,
        args: Sequence[str],
        search_context: str | None,
        display: bool = False,
        filename: str = "",
    ) -> str | None:
        """Raw TSV of one record pair from a current or promoted entry, if any."""

        pair_identity = self._pair_identity_from_fasta(query_fasta, subject_fasta)
        cache_key = build_protein_losat_cache_key(
            pair_identity,
            args=args,
            program="blastp",
            outfmt="6",
            search_context=search_context,
        )
        cached = self._find_cached_entry(
            cache_key=cache_key,
            pair_identity=pair_identity,
            program="blastp",
            outfmt="6",
            args=args,
        )
        if cached is not None:
            if display:
                self._mark_display(cache_key, filename)
            return str(cached.get("text") or "")
        # Legacy entries predate source batches; they answer only a search
        # whose sides are the two records themselves (as in the Web).
        if search_context is not None or not self._legacy_entries:
            return None
        query_proteins = self._proteins_by_instance.get(
            pair_identity.query_record_instance_key
        )
        subject_proteins = self._proteins_by_instance.get(
            pair_identity.subject_record_instance_key
        )
        if query_proteins is None or subject_proteins is None:
            return None
        scan = promote_legacy_protein_raw_cache_entries(
            self._legacy_entries,
            query_proteins=query_proteins,
            subject_proteins=subject_proteins,
            query_fasta=query_fasta,
            subject_fasta=subject_fasta,
            # _pair_identity_from_fasta above raised when the manifest is None.
            identity_manifest=cast(ProteinIdentityManifest, self.identity_manifest),
            expected_args=args,
            expected_program="blastp",
            expected_outfmt="6",
        )
        if scan.promotion is None:
            return None
        promotion = scan.promotion
        self.add_promoted_entry(promotion.entry)
        self._retire_legacy(promotion.candidate_index)
        if display:
            self._mark_display(cache_key, filename)
        return promotion.rewritten_tsv

    def store_search(
        self,
        query_fasta: str,
        subject_fasta: str,
        *,
        args: Sequence[str],
        search_context: str | None,
        text: str,
        runtime: Mapping[str, object] | None,
        display: bool = False,
        filename: str = "",
    ) -> None:
        """Store the raw TSV that a source job produced for one record pair."""

        pair_identity = self._pair_identity_from_fasta(query_fasta, subject_fasta)
        entry: dict[str, object] = {
            "schema": PROTEIN_LOSAT_CACHE_SCHEMA,
            "kind": "raw-losat",
            "identityKind": "protein",
            "idEncoding": "runtime-handle-v1",
            "key": build_protein_losat_cache_key(
                pair_identity,
                args=args,
                program="blastp",
                outfmt="6",
                search_context=search_context,
            ),
            "filename": filename if display else "",
            "display": bool(display),
            "text": str(text),
            "program": "blastp",
            "outfmt": "6",
            "args": [str(arg) for arg in args],
            **({"searchContext": search_context} if search_context else {}),
            "queryProteinSetHash": pair_identity.query_protein_set_hash,
            "subjectProteinSetHash": pair_identity.subject_protein_set_hash,
            "queryRuntimeBindingHash": pair_identity.query_runtime_binding_hash,
            "subjectRuntimeBindingHash": pair_identity.subject_runtime_binding_hash,
            "queryRecordInstanceKey": pair_identity.query_record_instance_key,
            "subjectRecordInstanceKey": pair_identity.subject_record_instance_key,
        }
        if not raw_protein_tsv_matches_bindings(
            str(entry["text"]),
            cast(ProteinIdentityManifest, self.identity_manifest),
            query_record_instance_key=pair_identity.query_record_instance_key,
            subject_record_instance_key=pair_identity.subject_record_instance_key,
        ):
            raise ValidationError("LOSAT output references IDs outside the current protein bindings.")
        self._add_search_entry(
            str(entry["key"]),
            entry,
            runtime=runtime,
            display=display,
            filename=filename,
        )


def _first_qualifier(feature: SeqFeature, key: str) -> str | None:
    value = feature.qualifiers.get(key)
    if value is None:
        return None
    if isinstance(value, (list, tuple)):
        if not value:
            return None
        value = value[0]
    text = str(value).strip()
    return text or None


def _qualifier_values(feature: SeqFeature, key: str) -> tuple[str, ...]:
    value = feature.qualifiers.get(key)
    if value is None:
        return ()
    if not isinstance(value, (list, tuple)):
        value = (value,)
    values: list[str] = []
    for item in value:
        text = str(item).strip()
        if text:
            values.append(text)
    return tuple(values)


def _first_qualifier_from_any(feature: SeqFeature, keys: Sequence[str]) -> str | None:
    for key in keys:
        value = _first_qualifier(feature, key)
        if value:
            return value
    return None


def _fasta_safe_protein_id(protein_id: str | None) -> str | None:
    if protein_id is None:
        return None
    return protein_id if _VALID_FASTA_ID_RE.fullmatch(protein_id) else None


def _unique_protein_id(
    preferred_id: str | None,
    fallback_id: str,
    used_ids: set[str],
) -> str:
    if preferred_id is not None and preferred_id not in used_ids:
        return preferred_id
    if fallback_id not in used_ids:
        return fallback_id
    suffix = 2
    while f"{fallback_id}_{suffix}" in used_ids:
        suffix += 1
    return f"{fallback_id}_{suffix}"


def _clean_protein_sequence(sequence: object | None) -> str | None:
    if sequence is None:
        return None
    protein = "".join(str(sequence).split()).upper()
    while protein.endswith("*"):
        protein = protein[:-1]
    if not protein or "*" in protein:
        return None
    if not _VALID_PROTEIN_RE.fullmatch(protein):
        return None
    return protein


def _feature_identifier(feature: SeqFeature) -> str | None:
    feature_id = _first_qualifier(feature, "ID")
    if feature_id:
        return feature_id
    raw_id = str(getattr(feature, "id", "") or "").strip()
    if raw_id and raw_id != "<unknown id>":
        return raw_id
    return None


def _iter_record_features(
    record: SeqRecord,
) -> list[tuple[int, SeqFeature, tuple[SeqFeature, ...]]]:
    items: list[tuple[int, SeqFeature, tuple[SeqFeature, ...]]] = []

    def walk(feature: SeqFeature, ancestors: tuple[SeqFeature, ...]) -> None:
        source_index = _source_feature_index(feature)
        items.append(
            (
                len(items) if source_index is None else source_index,
                feature,
                ancestors,
            )
        )
        sub_features = getattr(feature, "sub_features", None) or []
        for sub_feature in sub_features:
            walk(sub_feature, ancestors + (feature,))

    for feature in record.features:
        walk(feature, ())
    return items


def _build_parent_graph(
    feature_items: Sequence[tuple[int, SeqFeature, tuple[SeqFeature, ...]]],
) -> dict[str, tuple[str, tuple[str, ...]]]:
    graph: dict[str, tuple[str, tuple[str, ...]]] = {}
    for _feature_index, feature, ancestors in feature_items:
        feature_id = _feature_identifier(feature)
        if not feature_id:
            continue
        parent_ids = _qualifier_values(feature, "Parent")
        if not parent_ids and ancestors:
            ancestor_id = _feature_identifier(ancestors[-1])
            if ancestor_id:
                parent_ids = (ancestor_id,)
        graph[feature_id] = (str(feature.type), tuple(parent_ids))
    return graph


def _resolve_gene_parent_id(
    feature: SeqFeature,
    ancestors: Sequence[SeqFeature],
    parent_graph: Mapping[str, tuple[str, tuple[str, ...]]],
) -> str | None:
    for ancestor in reversed(tuple(ancestors)):
        if str(ancestor.type).lower() == "gene":
            ancestor_id = _feature_identifier(ancestor)
            if ancestor_id:
                return ancestor_id

    start_ids = list(_qualifier_values(feature, "Parent"))
    if not start_ids and ancestors:
        ancestor_id = _feature_identifier(ancestors[-1])
        if ancestor_id:
            start_ids.append(ancestor_id)

    seen: set[str] = set()

    def resolve(feature_id: str) -> str | None:
        if feature_id in seen:
            return None
        seen.add(feature_id)
        feature_type, parent_ids = parent_graph.get(feature_id, ("", ()))
        if str(feature_type).lower() == "gene":
            return feature_id
        for parent_id in parent_ids:
            resolved = resolve(parent_id)
            if resolved:
                return resolved
        return None

    for parent_id in start_ids:
        resolved = resolve(parent_id)
        if resolved:
            return resolved
    return None


def _translate_cds_feature(record: SeqRecord, feature: SeqFeature) -> str | None:
    try:
        protein = translate_cds(feature, feature.extract(record.seq))
    except Exception as exc:
        logger.debug(
            "Skipping CDS feature %s on %s: translation failed: %s",
            getattr(feature, "id", ""),
            record.id,
            exc,
        )
        return None
    return _clean_protein_sequence(protein)


def _cds_span(feature: SeqFeature) -> tuple[int, int, int | None] | None:
    location = feature.location
    if location is None:
        return None
    start = int(location.start)
    end = int(location.end)
    if start < 0 or end <= start:
        return None
    strand = location.strand if location.strand in {-1, 1} else None
    return start, end, strand


def _view_feature_hash_parts(
    feature: SeqFeature,
) -> tuple[tuple[int, int, int | None], ...]:
    """Return location parts in the currently rendered record view."""

    location = feature.location
    if location is None:
        return ()
    raw_parts = location.parts if hasattr(location, "parts") and location.parts else [location]
    return tuple(
        (
            int(part.start),
            int(part.end),
            int(part.strand) if part.strand in {-1, 1} else None,
        )
        for part in raw_parts
    )


def extract_cds_proteins(
    records: Sequence[SeqRecord],
    *,
    record_index_offset: int = 0,
    prefer_source_ids: bool = True,
    feature_visibility_rules: list[dict[str, object]] | None = None,
    feature_overrides: Sequence[Mapping[int, object]] | None = None,
) -> ProteinExtractionResult:
    """Extract CDS proteins using protein IDs where possible and genomic spans.

    Coordinates in the returned metadata are Python/Biopython-style 0-based start
    and 0-based exclusive end values. ``feature_overrides`` holds, per record,
    the identity-addressed edits by source feature index, aligned with
    ``records``; they decide before the visibility rules.
    """

    candidate_records: list[list[_CdsProteinCandidate]] = []
    source_id_counts: dict[str, int] = {}
    protein_map: dict[str, CdsProtein] = {}

    record_index_offset = int(record_index_offset)

    for record_index, record in enumerate(records):
        global_record_index = record_index + record_index_offset
        record_candidates: list[_CdsProteinCandidate] = []
        cds_count = 0
        feature_items = _iter_record_features(record)
        parent_graph = _build_parent_graph(feature_items)
        record_overrides = feature_overrides[record_index] if feature_overrides else {}
        for feature_index, feature, ancestors in feature_items:
            if feature.type != "CDS":
                continue
            if not should_include_feature_in_analysis(
                feature,
                feature_visibility_rules,
                record_id=record.id,
                feature_override=record_overrides.get(feature_index),
            ):
                continue

            span = _cds_span(feature)
            if span is None:
                continue
            start, end, strand = span

            protein_sequence = _clean_protein_sequence(
                _first_qualifier(feature, "translation")
            )
            if protein_sequence is None:
                protein_sequence = _translate_cds_feature(record, feature)
            if protein_sequence is None:
                continue

            cds_count += 1
            feature_svg_id, hash_parts = source_feature_hash_and_parts(feature, record)
            view_hash_parts = _view_feature_hash_parts(feature)
            view_feature_svg_id = (
                compute_feature_hash_from_location_parts(
                    str(feature.type or ""),
                    view_hash_parts,
                    record_id=record.id,
                )
                if view_hash_parts
                else None
            )
            hash_start, hash_end, hash_strand = hash_parts[0] if hash_parts else (None, None, None)
            synthetic_protein_id = f"gbd_r{global_record_index + 1:04d}_cds{cds_count:06d}"
            source_protein_id = _first_qualifier(feature, "protein_id")
            fasta_safe_source_id = _fasta_safe_protein_id(source_protein_id)
            if fasta_safe_source_id is not None:
                source_id_counts[fasta_safe_source_id] = source_id_counts.get(fasta_safe_source_id, 0) + 1
            gene = _first_qualifier(feature, "gene")
            product = _first_qualifier(feature, "product")
            note = _first_qualifier(feature, "note")
            locus_tag = _first_qualifier(feature, "locus_tag")
            gene_id = _first_qualifier(feature, "gene_id")
            old_locus_tag = _first_qualifier(feature, "old_locus_tag")
            db_xref = _qualifier_values(feature, "db_xref")
            gff_id = _feature_identifier(feature)
            parent_ids = _qualifier_values(feature, "Parent")
            if not parent_ids and ancestors:
                ancestor_id = _feature_identifier(ancestors[-1])
                if ancestor_id:
                    parent_ids = (ancestor_id,)
            gene_parent_id = _resolve_gene_parent_id(feature, ancestors, parent_graph)
            label = (
                _first_qualifier_from_any(
                    feature, ("locus_tag", "protein_id", "gene", "product")
                )
                or synthetic_protein_id
            )
            record_candidates.append(
                _CdsProteinCandidate(
                    protein_id=synthetic_protein_id,
                    source_protein_id=source_protein_id,
                    record_index=global_record_index,
                    feature_index=feature_index,
                    record_id=str(record.id),
                    start=start,
                    end=end,
                    strand=strand,
                    label=label,
                    protein_length=len(protein_sequence),
                    sequence=protein_sequence,
                    feature_svg_id=feature_svg_id,
                    view_feature_svg_id=view_feature_svg_id,
                    gene=gene,
                    product=product,
                    note=note,
                    locus_tag=locus_tag,
                    gene_id=gene_id,
                    old_locus_tag=old_locus_tag,
                    db_xref=db_xref,
                    gff_id=gff_id,
                    parent_ids=parent_ids,
                    gene_parent_id=gene_parent_id,
                    feature_type=str(feature.type),
                    feature_hash_start=hash_start,
                    feature_hash_end=hash_end,
                    feature_hash_strand=hash_strand,
                    feature_hash_parts=hash_parts,
                    location_operator=str(getattr(feature.location, "operator", None) or ""),
                    source_feature_position=feature_index,
                )
            )

        candidate_records.append(record_candidates)

    proteins_by_record: list[list[CdsProtein]] = []
    used_ids: set[str] = set()
    for record_candidates in candidate_records:
        record_proteins: list[CdsProtein] = []
        for candidate in record_candidates:
            fasta_safe_source_id = _fasta_safe_protein_id(candidate.source_protein_id)
            preferred_id = (
                fasta_safe_source_id
                if prefer_source_ids
                and fasta_safe_source_id is not None
                and source_id_counts.get(fasta_safe_source_id, 0) == 1
                else candidate.protein_id
            )
            protein_id = _unique_protein_id(preferred_id, candidate.protein_id, used_ids)
            cds_protein = CdsProtein(
                protein_id=protein_id,
                record_index=candidate.record_index,
                feature_index=candidate.feature_index,
                record_id=candidate.record_id,
                start=candidate.start,
                end=candidate.end,
                strand=candidate.strand,
                label=candidate.label,
                protein_length=candidate.protein_length,
                sequence=candidate.sequence,
                source_protein_id=candidate.source_protein_id,
                feature_svg_id=candidate.feature_svg_id,
                view_feature_svg_id=candidate.view_feature_svg_id,
                gene=candidate.gene,
                product=candidate.product,
                note=candidate.note,
                locus_tag=candidate.locus_tag,
                gene_id=candidate.gene_id,
                old_locus_tag=candidate.old_locus_tag,
                db_xref=candidate.db_xref,
                gff_id=candidate.gff_id,
                parent_ids=candidate.parent_ids,
                gene_parent_id=candidate.gene_parent_id,
                feature_type=candidate.feature_type,
                feature_hash_start=candidate.feature_hash_start,
                feature_hash_end=candidate.feature_hash_end,
                feature_hash_strand=candidate.feature_hash_strand,
                feature_hash_parts=candidate.feature_hash_parts,
                location_operator=candidate.location_operator,
                source_feature_position=candidate.source_feature_position,
            )
            record_proteins.append(cds_protein)
            protein_map[protein_id] = cds_protein
            used_ids.add(protein_id)

        proteins_by_record.append(record_proteins)

    return ProteinExtractionResult(
        proteins_by_record=proteins_by_record,
        protein_map=protein_map,
    )


def proteins_to_fasta(proteins: Sequence[CdsProtein]) -> str:
    """Return FASTA text for CDS proteins."""

    lines: list[str] = []
    for protein in proteins:
        lines.append(f">{protein.protein_id} {protein.label}")
        sequence = protein.sequence
        lines.extend(sequence[index : index + 80] for index in range(0, len(sequence), 80))
    return "\n".join(lines) + ("\n" if lines else "")


def _finite_numeric_values(df: DataFrame) -> np.ndarray | None:
    """The numeric comparison columns as one float matrix, or None unless every value is finite."""

    try:
        values = df[list(_NUMERIC_COMPARISON_COLUMNS)].to_numpy(dtype=float)
    except (TypeError, ValueError, OverflowError):
        return None
    return values if bool(np.isfinite(values).all()) else None


def _coerce_outfmt6_numeric_columns(df: DataFrame) -> DataFrame:
    coerced = df.copy()
    for column in _NUMERIC_COMPARISON_COLUMNS:
        coerced[column] = pd.to_numeric(coerced[column], errors="coerce")
    if _finite_numeric_values(coerced) is None:
        raise ParseError("LOSATP blastp output contains non-numeric outfmt 6 fields.")
    return coerced


def _is_strict_outfmt6_numeric_field(column: str, value: str) -> bool:
    pattern = (
        _STRICT_DECIMAL_INTEGER_RE
        if column in _INTEGER_COMPARISON_COLUMNS
        else _STRICT_DECIMAL_NUMBER_RE
    )
    if pattern.fullmatch(value) is None:
        return False
    try:
        return math.isfinite(float(value))
    except ValueError:  # pragma: no cover - the grammar excludes float parse failures
        return False


def _validate_losatp_outfmt6_lines(data_lines: Sequence[str]) -> None:
    for line in data_lines:
        columns = line.split("\t")
        if len(columns) != len(COMPARISON_COLUMNS):
            raise ParseError("LOSATP blastp output does not match outfmt 6 columns.")
        if any(
            not _is_strict_outfmt6_numeric_field(column, columns[index])
            for index, column in enumerate(COMPARISON_COLUMNS[2:], start=2)
        ):
            raise ParseError(
                "LOSATP blastp output contains non-numeric outfmt 6 fields."
            )


def parse_losatp_outfmt6(text: str) -> DataFrame:
    """Parse LOSATP blastp outfmt 6 text into the standard comparison columns."""

    data_lines = [
        line
        for line in text.splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    if not data_lines:
        return pd.DataFrame(columns=COMPARISON_COLUMNS)
    data = "\n".join(data_lines)
    # One regex pass accepts the strict grammar. The line-by-line check runs
    # only to choose the error that the first failing line raises, and to
    # settle a value whose float is not finite, for which read_comparison_table
    # raises.
    if _NON_STRICT_OUTFMT6_LINE_RE.search("\n" + data) is not None:
        _validate_losatp_outfmt6_lines(data_lines)
    try:
        table = read_comparison_table(StringIO(data), label="LOSATP blastp output")
    except ValidationError as exc:
        _validate_losatp_outfmt6_lines(data_lines)
        raise ParseError(f"Failed to parse LOSATP blastp output: {exc}") from exc
    return table


def _validate_max_hits(max_hits: int, *, option_name: str = "protein_blastp_max_hits") -> None:
    if int(max_hits) <= 0:
        raise ValidationError(f"{option_name} must be > 0")


def _validate_candidate_limit(candidate_limit: int | None) -> None:
    if candidate_limit is not None and int(candidate_limit) <= 0:
        raise ValidationError("protein_blastp_candidate_limit must be > 0 or None")


def _validate_losatp_threads(threads: int | None) -> None:
    if threads is not None and int(threads) <= 0:
        raise ValidationError("losatp_threads must be > 0 or None")


def normalize_protein_blastp_mode(mode: str | None) -> ProteinBlastpMode:
    """Return a validated LOSATP blastp mode."""

    normalized = str(mode or "none").strip().lower()
    if normalized not in PROTEIN_BLASTP_MODES:
        raise ValidationError(
            "protein_blastp_mode must be one of: "
            + ", ".join(PROTEIN_BLASTP_MODES)
        )
    return cast(ProteinBlastpMode, normalized)


def normalize_orthogroup_membership_mode(mode: str | None) -> OrthogroupMembershipMode:
    """Return the active orthogroup inference model.

    Legacy user/session values are accepted as aliases so old configurations
    load without reactivating the removed inference modes.
    """

    normalized = read_spelling(mode, ORTHOGROUP_MEMBERSHIP_MODES, ORTHOGROUP_MEMBERSHIP_MODE_SPELLINGS)
    if normalized is None:
        raise ValidationError(
            "orthogroup_membership_mode must be "
            f"{ORTHOGROUP_INFERENCE_VERSION} or a legacy alias: "
            + ", ".join(_LEGACY_ORTHOGROUP_MEMBERSHIP_MODES)
        )
    return cast(OrthogroupMembershipMode, normalized)


def _validate_comparison_columns(hits: DataFrame) -> None:
    missing_columns = set(COMPARISON_COLUMNS).difference(hits.columns)
    if missing_columns:
        raise ParseError(
            "LOSATP blastp output is missing required columns: "
            + ", ".join(sorted(missing_columns))
        )


def filter_protein_hits_by_thresholds(
    hits: DataFrame,
    *,
    evalue: float,
    bitscore: float,
    identity: float,
    alignment_length: int,
) -> DataFrame:
    """Keep LOSATP blastp hits that pass visible pairwise-match thresholds."""

    if hits.empty:
        return hits.copy()
    _validate_comparison_columns(hits)

    try:
        evalue_threshold = float(evalue)
        bitscore_threshold = float(bitscore)
        identity_threshold = float(identity)
        alignment_length_threshold = int(alignment_length)
    except (TypeError, ValueError) as exc:
        raise ValidationError("protein colinearity thresholds must be numeric") from exc
    if alignment_length_threshold < 0:
        raise ValidationError("alignment_length must be >= 0")

    coerced_hits = _coerce_outfmt6_numeric_columns(hits)
    return coerced_hits.loc[
        (coerced_hits["evalue"] <= evalue_threshold)
        & (coerced_hits["bitscore"] >= bitscore_threshold)
        & (coerced_hits["identity"] >= identity_threshold)
        & (coerced_hits["alignment_length"] >= alignment_length_threshold)
    ].reset_index(drop=True)


def _sort_hits_for_query_best(hits: DataFrame) -> DataFrame:
    return hits.sort_values(
        ["query", "bitscore", "evalue", "identity", "alignment_length", "subject"],
        ascending=[True, False, True, False, False, True],
        kind="mergesort",
    )


def _sort_hits_for_subject_best(hits: DataFrame) -> DataFrame:
    return hits.sort_values(
        ["subject", "bitscore", "evalue", "identity", "alignment_length", "query"],
        ascending=[True, False, True, False, False, True],
        kind="mergesort",
    )


def select_top_hits_per_query(
    hits: DataFrame,
    *,
    max_hits: int,
) -> DataFrame:
    """Keep the strongest distinct subject protein hits per query protein."""

    _validate_max_hits(max_hits)
    if hits.empty:
        return hits.copy()
    _validate_comparison_columns(hits)

    sorted_hits = _sort_hits_for_query_best(_coerce_outfmt6_numeric_columns(hits))
    sorted_hits = sorted_hits.drop_duplicates(["query", "subject"], keep="first")
    return (
        sorted_hits.groupby("query", group_keys=False, sort=False)
        .head(int(max_hits))
        .reset_index(drop=True)
    )


def _select_member_candidate_hits_per_query(
    hits: DataFrame,
    *,
    max_hits: int | None,
) -> DataFrame:
    """Keep every HSP for the strongest distinct member candidates.

    Without a limit every HSP is kept, and the result shares the columns of
    ``hits`` instead of copying them; callers only read it.
    """

    if max_hits is None:
        _validate_comparison_columns(hits)
        shared = hits.copy(deep=False)
        shared.index = pd.RangeIndex(len(shared))
        return shared
    selected_pairs = select_top_hits_per_query(hits, max_hits=max_hits)
    if selected_pairs.empty:
        return hits.iloc[0:0].copy()
    pair_index = pd.MultiIndex.from_arrays([selected_pairs["query"], selected_pairs["subject"]])
    hit_index = pd.MultiIndex.from_arrays([hits["query"], hits["subject"]])
    return hits.loc[hit_index.isin(pair_index)].reset_index(drop=True)


def select_reciprocal_best_hits(hits: DataFrame) -> DataFrame:
    """Keep query-subject pairs that are mutual best hits in one adjacent-pair table."""

    if hits.empty:
        return hits.copy()
    _validate_comparison_columns(hits)

    coerced_hits = _coerce_outfmt6_numeric_columns(hits)
    query_best = (
        _sort_hits_for_query_best(coerced_hits)
        .drop_duplicates("query", keep="first")
    )
    subject_best = (
        _sort_hits_for_subject_best(coerced_hits)
        .drop_duplicates("subject", keep="first")
    )
    reciprocal_pairs = {
        (str(query), str(subject))
        for query, subject in zip(subject_best["query"], subject_best["subject"])
    }
    return query_best.loc[
        [
            (str(query), str(subject)) in reciprocal_pairs
            for query, subject in zip(query_best["query"], query_best["subject"])
        ]
    ].reset_index(drop=True)


def select_best_hits_per_query(
    hits: DataFrame,
) -> DataFrame:
    """Choose the deterministic best subject for each query protein."""

    if hits.empty:
        return hits.copy()
    _validate_comparison_columns(hits)

    return (
        _sort_hits_for_query_best(_coerce_outfmt6_numeric_columns(hits))
        .drop_duplicates("query", keep="first")
        .reset_index(drop=True)
    )


def select_reciprocal_best_hit_edges(
    forward_hits: DataFrame,
    reverse_hits: DataFrame,
) -> DataFrame:
    """Keep RBH edges from directional best-hit tables.

    Returned rows use the forward orientation: query record i, subject record j.
    """

    if forward_hits.empty or reverse_hits.empty:
        return forward_hits.iloc[0:0].copy()

    forward_best = select_best_hits_per_query(forward_hits)
    reverse_best = select_best_hits_per_query(reverse_hits)
    reverse_best_by_query = {
        str(query): str(subject)
        for query, subject in zip(reverse_best["query"], reverse_best["subject"])
    }
    keep_mask = [
        reverse_best_by_query.get(str(subject)) == str(query)
        for query, subject in zip(forward_best["query"], forward_best["subject"])
    ]
    return forward_best.loc[keep_mask].reset_index(drop=True)


def _normalized_fit_row_positions(
    normalized_hits: DataFrame,
    *,
    top_fraction: float,
) -> list[int]:
    """Positions of the fit rows, in fit order.

    The length-sorted table is cut into at most eight bins; each bin keeps the
    top fraction of its rows by bitscore. One stable sort with the bin as its
    first key orders every bin as a separate sort of that bin would.
    """

    projection = normalized_hits[["length_product", "query", "subject", "bitscore"]].reset_index(drop=True)
    sorted_hits = projection.sort_values(
        ["length_product", "query", "subject"],
        ascending=[True, True, True],
        kind="mergesort",
    )
    row_count = len(sorted_hits)
    bin_size = max(4, int(math.ceil(row_count / 8)))
    ranked = sorted_hits.assign(_bin=np.arange(row_count) // bin_size).sort_values(
        ["_bin", "bitscore", "query", "subject"],
        ascending=[True, False, True, True],
        kind="mergesort",
    )
    top_counts = np.array(
        [
            max(1, int(math.ceil(min(bin_size, row_count - start) * float(top_fraction))))
            for start in range(0, row_count, bin_size)
        ],
        dtype=np.int64,
    )
    ranked_bins = ranked["_bin"].to_numpy()
    top = (np.arange(row_count) - ranked_bins * bin_size) < top_counts[ranked_bins]
    return [int(position) for position in ranked.index[top]]


def _select_normalized_fit_rows(
    normalized_hits: DataFrame,
    *,
    top_fraction: float,
) -> list[tuple[float, float]]:
    if normalized_hits.empty:
        return []
    positions = _normalized_fit_row_positions(normalized_hits, top_fraction=top_fraction)
    fit_points: list[tuple[float, float]] = []
    for length_product, bitscore in zip(
        normalized_hits["length_product"].to_numpy(dtype=float)[positions].tolist(),
        normalized_hits["bitscore"].to_numpy(dtype=float)[positions].tolist(),
    ):
        if length_product <= 0.0 or bitscore <= 0.0:
            continue
        fit_points.append((math.log10(length_product), math.log10(bitscore)))
    return fit_points


def _fit_expected_bitscore_model(
    normalized_hits: DataFrame,
    *,
    top_fraction: float,
) -> tuple[float, float] | None:
    if len(normalized_hits) < _ORTHOGROUP_SPLIT_MIN_FIT_HITS:
        return None
    fit_points = _select_normalized_fit_rows(
        normalized_hits,
        top_fraction=top_fraction,
    )
    if len(fit_points) < 3:
        return None
    x_values = [point[0] for point in fit_points]
    y_values = [point[1] for point in fit_points]
    if len(set(round(value, 9) for value in x_values)) < 2:
        return None
    # compensated_sum: CPython 3.12+ sum() on every Python (OV-1237).
    x_mean = compensated_sum(x_values) / len(x_values)
    y_mean = compensated_sum(y_values) / len(y_values)
    denominator = compensated_sum((x_value - x_mean) ** 2 for x_value in x_values)
    if denominator <= _ORTHOGROUP_SPLIT_EPSILON:
        return None
    slope = compensated_sum(
        (x_value - x_mean) * (y_value - y_mean)
        for x_value, y_value in zip(x_values, y_values)
    ) / denominator
    intercept = y_mean - slope * x_mean
    if not (math.isfinite(slope) and math.isfinite(intercept)):
        return None
    return slope, intercept


def _empty_normalized_hit_table(columns: Sequence[str]) -> DataFrame:
    extra_columns = (
        "query_length",
        "subject_length",
        "length_product",
        "hsp_count",
        "query_covered_length",
        "subject_covered_length",
        "query_coverage",
        "subject_coverage",
        "min_coverage",
        "representative_alignment_length",
        "total_hsp_alignment_length",
        "coverage_source",
        "normalized_score",
        "normalization_fallback",
        "record_pair_normalization_source",
        "domain_only",
    )
    return pd.DataFrame(columns=tuple(dict.fromkeys([*columns, *extra_columns])))


def _python_scalar_columns(rows: DataFrame) -> dict[str, object]:
    """The columns ``DataFrame.from_records`` builds from these rows' ``itertuples``.

    Names follow the namedtuple renaming of ``itertuples``. Numpy numeric
    columns become int64, float64 or bool, as their Python scalars would;
    every other column is rebuilt from its Python values.
    """

    names = cast(type[NamedTuple], namedtuple("Pandas", list(rows.columns), rename=True))._fields
    numeric: dict[str, np.ndarray] = {}
    other_names: list[str] = []
    other_values: list[list[object]] = []
    for position, name in enumerate(names):
        column = rows.iloc[:, position]
        dtype = column.dtype
        if isinstance(dtype, np.dtype) and dtype.kind == "f":
            numeric[name] = column.to_numpy(dtype=np.float64)
        elif isinstance(dtype, np.dtype) and dtype.kind == "b":
            numeric[name] = column.to_numpy(dtype=bool)
        elif isinstance(dtype, np.dtype) and (dtype.kind == "i" or (dtype.kind == "u" and dtype.itemsize < 8)):
            numeric[name] = column.to_numpy(dtype=np.int64)
        else:
            other_names.append(name)
            other_values.append(column.tolist())
    # Row tuples take the same per-column type inference as the records of itertuples.
    other = pd.DataFrame.from_records(list(zip(*other_values)), columns=other_names)
    return {name: numeric[name] if name in numeric else other[name] for name in names}


def _first_minimum_per_group(
    groups: np.ndarray,
    keys: Sequence[np.ndarray],
    group_count: int,
) -> np.ndarray:
    """Position of each group's smallest row, comparing keys in order; the earliest row wins a tie.

    Groups come in ascending order, and a group without rows is left out.
    Only rows that share their group are sorted.
    """

    best = np.full(group_count, -1, dtype=np.int64)
    shared = np.bincount(groups, minlength=group_count)[groups] > 1
    single = np.flatnonzero(~shared)
    best[groups[single]] = single
    tied = np.flatnonzero(shared)
    if len(tied):
        # lexsort is stable, so rows with equal keys keep their order.
        order = tied[np.lexsort((*(key[tied] for key in reversed(keys)), groups[tied]))]
        ordered_groups = groups[order]
        first = order[np.r_[True, ordered_groups[1:] != ordered_groups[:-1]]]
        best[groups[first]] = first
    return best[best >= 0]


def _first_occurrences(keys: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Each key's code, numbered in order of first appearance, and each code's first position.

    The positions ascend, since the codes follow first appearance.
    """

    codes = pd.factorize(keys)[0]
    return codes, _first_minimum_per_group(codes, (), int(codes.max(initial=-1)) + 1)


def _covered_lengths(
    groups: np.ndarray,
    start: np.ndarray,
    end: np.ndarray,
    lengths: np.ndarray,
    group_count: int,
) -> np.ndarray:
    """Length of the union of each group's HSP intervals, clamped to 1..length.

    Intervals that touch (end + 1 == next start) merge. Within each group,
    sorted by start, a running maximum of the ends marks where a new run begins.
    """

    low = np.maximum(1.0, np.minimum(start, end))
    high = np.minimum(lengths[groups].astype(float), np.maximum(start, end))
    inside = low <= high
    covered = np.zeros(group_count, dtype=np.int64)
    if not inside.any():
        return covered
    group = groups[inside]
    low_int = low[inside].astype(np.int64)
    high_int = high[inside].astype(np.int64)
    order = np.lexsort((high_int, low_int, group))
    group, low_int, high_int = group[order], low_int[order], high_int[order]
    # Offsetting each group by more than any end keeps one running maximum per group.
    offset = group.astype(np.int64) * (int(high_int.max()) + 2)
    low_int, high_int = low_int + offset, high_int + offset
    running_end = np.maximum.accumulate(high_int)
    run_starts = np.flatnonzero(
        np.r_[True, (group[1:] != group[:-1]) | (low_int[1:] > running_end[:-1] + 1)]
    )
    run_ends = np.r_[run_starts[1:] - 1, len(group) - 1]
    np.add.at(covered, group[run_starts], running_end[run_ends] - low_int[run_starts] + 1)
    return covered


def _aggregate_hsps_by_protein_pair(
    hits: DataFrame,
    protein_map: Mapping[str, CdsProtein],
) -> DataFrame:
    """One row per protein pair: its best HSP plus the union coverage of all its HSPs.

    The numeric columns must be finite numbers, as ``_coerce_outfmt6_numeric_columns``
    leaves them. Pairs keep their order of first appearance; the best HSP is the
    smallest (-bitscore, evalue, -identity, -alignment length, qstart, qend,
    sstart, send, row position), with coordinates truncated to integers.
    """

    if hits is None or hits.empty:
        columns = tuple(hits.columns) if hits is not None else tuple(COMPARISON_COLUMNS)
        return _empty_normalized_hit_table(columns)
    values = _finite_numeric_values(hits)
    if values is None:
        raise ParseError("LOSATP blastp output contains non-numeric outfmt 6 fields.")

    # Pandas groups the ID pairs (its key coercion and missing-ID handling) in
    # order of appearance; every later step works on whole columns.
    grouped = hits.groupby(["query", "subject"], sort=False, observed=True)
    group_count = int(grouped.ngroups)
    pair_keys = cast(pd.MultiIndex, grouped.size().index)

    def protein_lengths(level: int) -> np.ndarray:
        # One lookup per distinct ID; an unknown protein has length 0.
        lengths = [
            max(0, int(protein.protein_length)) if (protein := protein_map.get(str(value))) is not None else 0
            for value in pair_keys.levels[level]
        ]
        return np.asarray(lengths, dtype=np.int64)[pair_keys.codes[level]]

    query_lengths, subject_lengths = protein_lengths(0), protein_lengths(1)
    row_groups = np.asarray(grouped.ngroup(), dtype=float)
    rows = np.flatnonzero(np.isfinite(row_groups) & (row_groups >= 0))
    groups = row_groups[rows].astype(np.int64)
    in_pair = (query_lengths[groups] > 0) & (subject_lengths[groups] > 0)
    rows, groups = rows[in_pair], groups[in_pair]
    if not len(rows):
        return _empty_normalized_hit_table(tuple(hits.columns))

    column = {name: values[rows, index] for index, name in enumerate(_NUMERIC_COMPARISON_COLUMNS)}
    bitscore = column["bitscore"]
    alignment_length = np.trunc(column["alignment_length"])
    qstart, qend, sstart, send = (np.trunc(column[name]) for name in ("qstart", "qend", "sstart", "send"))
    representative = _first_minimum_per_group(
        groups,
        (-bitscore, column["evalue"], -column["identity"], -alignment_length, qstart, qend, sstart, send),
        group_count,
    )
    pair_groups = groups[representative]
    keep = (bitscore[representative] > 0.0) & (alignment_length[representative] > 0)
    if not keep.any():
        return _empty_normalized_hit_table(tuple(hits.columns))
    pair_groups, representative = pair_groups[keep], representative[keep]

    query_length = query_lengths[pair_groups]
    subject_length = subject_lengths[pair_groups]
    query_covered = _covered_lengths(groups, qstart, qend, query_lengths, group_count)[pair_groups]
    subject_covered = _covered_lengths(groups, sstart, send, subject_lengths, group_count)[pair_groups]
    query_coverage = np.minimum(1.0, query_covered / query_length)
    subject_coverage = np.minimum(1.0, subject_covered / subject_length)
    total_alignment = np.bincount(
        groups, weights=np.where(alignment_length > 0, alignment_length, 0.0), minlength=group_count
    )[pair_groups]

    table = _python_scalar_columns(hits.iloc[rows[representative]])
    # A summary column that the hits already have keeps its position.
    table.update(
        {
            "query_length": query_length,
            "subject_length": subject_length,
            "length_product": (query_length * subject_length).astype(float),
            "hsp_count": np.bincount(groups, minlength=group_count)[pair_groups].astype(np.int64),
            "query_covered_length": query_covered,
            "subject_covered_length": subject_covered,
            "query_coverage": query_coverage,
            "subject_coverage": subject_coverage,
            "min_coverage": np.minimum(query_coverage, subject_coverage),
            "representative_alignment_length": alignment_length[representative].astype(np.int64),
            "total_hsp_alignment_length": total_alignment.astype(np.int64),
            "coverage_source": "hsp_union",
        }
    )
    return pd.DataFrame(table)


def _normalize_directional_hit_table(
    hits: DataFrame,
    protein_map: Mapping[str, CdsProtein],
    *,
    top_fraction: float = _ORTHOGROUP_SPLIT_TOP_FRACTION,
    min_coverage: float = _ORTHOGROUP_SPLIT_MIN_COVERAGE,
) -> DataFrame:
    if hits is None or hits.empty:
        columns = tuple(hits.columns) if hits is not None else tuple(COMPARISON_COLUMNS)
        return _empty_normalized_hit_table(columns)
    _validate_comparison_columns(hits)
    coerced_hits = _coerce_outfmt6_numeric_columns(hits)
    aggregated_hits = _aggregate_hsps_by_protein_pair(coerced_hits, protein_map)
    if aggregated_hits.empty:
        return _empty_normalized_hit_table(tuple(coerced_hits.columns))

    normalized_hits = aggregated_hits.loc[
        aggregated_hits["min_coverage"].astype(float) >= float(min_coverage)
    ].reset_index(drop=True)
    if normalized_hits.empty:
        return _empty_normalized_hit_table(tuple(coerced_hits.columns))

    model = _fit_expected_bitscore_model(
        normalized_hits,
        top_fraction=top_fraction,
    )
    fallback = model is None
    normalized_scores: list[float] = []
    # Scalar math keeps math.log10 and ** exactly as before; numpy's can differ in the last bit.
    for bitscore, length_product in zip(
        normalized_hits["bitscore"].to_numpy(dtype=float).tolist(),
        normalized_hits["length_product"].to_numpy(dtype=float).tolist(),
    ):
        if model is None:
            denominator = math.sqrt(max(length_product, _ORTHOGROUP_SPLIT_EPSILON))
            normalized_scores.append(bitscore / denominator)
            continue
        slope, intercept = model
        expected_bitscore = 10 ** (slope * math.log10(length_product) + intercept)
        if not math.isfinite(expected_bitscore) or expected_bitscore <= 0.0:
            denominator = math.sqrt(max(length_product, _ORTHOGROUP_SPLIT_EPSILON))
            normalized_scores.append(bitscore / denominator)
        else:
            normalized_scores.append(bitscore / expected_bitscore)
    normalized_hits["normalized_score"] = normalized_scores
    normalized_hits["normalization_fallback"] = bool(fallback)
    normalized_hits["record_pair_normalization_source"] = (
        "sqrt_length_fallback" if fallback else "fit"
    )
    normalized_hits["domain_only"] = normalized_hits["min_coverage"].astype(float) <= _DOMAIN_ONLY_MAX_MIN_COVERAGE
    return normalized_hits.reset_index(drop=True)


def _normalize_directional_hit_tables(
    directional_tables: Mapping[tuple[int, int], DataFrame],
    protein_map: Mapping[str, CdsProtein],
    *,
    min_coverage: float = _ORTHOGROUP_SPLIT_MIN_COVERAGE,
) -> dict[tuple[int, int], DataFrame]:
    return {
        (int(pair[0]), int(pair[1])): _normalize_directional_hit_table(
            hits,
            protein_map,
            min_coverage=float(min_coverage),
        )
        for pair, hits in directional_tables.items()
    }


def _comparison_columns_only(hits: DataFrame) -> DataFrame:
    if hits.empty:
        return pd.DataFrame(columns=COMPARISON_COLUMNS)
    return hits.loc[:, list(COMPARISON_COLUMNS)].copy()


class _UnionFind:
    def __init__(self) -> None:
        self._parent: dict[str, str] = {}

    def add(self, item: str) -> None:
        if item not in self._parent:
            self._parent[item] = item

    def find(self, item: str) -> str:
        self.add(item)
        parent = self._parent[item]
        if parent != item:
            parent = self.find(parent)
            self._parent[item] = parent
        return parent

    def union(self, left: str, right: str) -> None:
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root == right_root:
            return
        if right_root < left_root:
            left_root, right_root = right_root, left_root
        self._parent[right_root] = left_root


def _is_better_member_rank(
    candidate: tuple[float, float, float, float],
    current: tuple[float, float, float, float] | None,
) -> bool:
    if current is None:
        return True
    return (
        candidate[0] > current[0]
        or (candidate[0] == current[0] and candidate[1] < current[1])
        or (
            candidate[0] == current[0]
            and candidate[1] == current[1]
            and candidate[2] > current[2]
        )
        or (
            candidate[0] == current[0]
            and candidate[1] == current[1]
            and candidate[2] == current[2]
            and candidate[3] > current[3]
        )
    )


def _protein_sort_key(protein: CdsProtein) -> tuple[int, int, int, str]:
    return (
        int(protein.record_index),
        int(protein.start),
        int(protein.end),
        str(protein.protein_id),
    )


_ORTHOGROUP_NAME_SOURCE_WEIGHTS = {
    "gene": 80,
    "product": 55,
    "note": 35,
    "label": 10,
}
_ORTHOGROUP_NAME_SOURCE_ORDER = ("gene", "product", "note", "label")
_ORTHOGROUP_NAME_SOURCE_RANK = {
    source: index for index, source in enumerate(_ORTHOGROUP_NAME_SOURCE_ORDER)
}
_NOTE_PREFIX_RE = re.compile(r"^(?:product|gene|note)\s*:\s*", re.IGNORECASE)
_DUF_CANDIDATE_RE = re.compile(r"\bDUF\d+\b.*\bdomain-containing protein\b", re.IGNORECASE)
_ANNOTATION_PROVENANCE_NOTE_RE = re.compile(
    r"^derived by automated computational analysis\b",
    re.IGNORECASE,
)
_GENERIC_ORTHOGROUP_NAME_CANDIDATES = {
    "hypothetical protein",
    "uncharacterized protein",
    "unknown protein",
    "predicted protein",
    "putative protein",
}
_WEAK_ORTHOGROUP_NAME_CANDIDATES = {
    "hypothetical protein",
    "uncharacterized protein",
    "unknown protein",
    "predicted protein",
    "putative protein",
}


# Annotation texts repeat across members and groups; both checks are pure.
@lru_cache(maxsize=8192)
def _normalize_orthogroup_name_text(
    text: str | None,
    source: str,
) -> tuple[str, str] | None:
    if text is None:
        return None
    normalized = re.sub(r"\s+", " ", str(text)).strip()
    normalized = _NOTE_PREFIX_RE.sub("", normalized).strip()
    normalized = normalized.strip(" \t\r\n.,;:()[]{}")
    if not normalized:
        return None
    if source == "note" and _ANNOTATION_PROVENANCE_NOTE_RE.search(normalized):
        return None
    compact_length = len(re.sub(r"[^0-9A-Za-z]+", "", normalized))
    if source != "gene" and compact_length < 4:
        return None
    key = normalized.casefold()
    if key in _GENERIC_ORTHOGROUP_NAME_CANDIDATES:
        return None
    return key, normalized


@lru_cache(maxsize=8192)
def _is_weak_orthogroup_name_candidate(text: str, source: str) -> bool:
    normalized = re.sub(r"\s+", " ", str(text)).strip().casefold()
    if source == "label":
        return True
    if normalized in _WEAK_ORTHOGROUP_NAME_CANDIDATES:
        return True
    return bool(_DUF_CANDIDATE_RE.search(text))


def _member_annotation_value(member: OrthogroupMember, source: str) -> str | None:
    if source == "product":
        return member.product
    if source == "gene":
        return member.gene
    if source == "note":
        return member.note
    if source == "label":
        return member.label
    return None


class _OrthogroupNameAccumulator(TypedDict):
    text: str
    source: str
    source_weight: int
    source_rank: int
    members: set[str]
    records: set[int]
    representatives: set[str]
    weak: bool


def _iter_orthogroup_name_candidates(
    members: Sequence[OrthogroupMember],
    *,
    include_label: bool,
) -> list[OrthogroupNameCandidate]:
    candidate_sources = _ORTHOGROUP_NAME_SOURCE_ORDER if include_label else _ORTHOGROUP_NAME_SOURCE_ORDER[:-1]
    accumulators: dict[str, _OrthogroupNameAccumulator] = {}

    for member in members:
        member_id = str(member.protein_id)
        for source in candidate_sources:
            normalized = _normalize_orthogroup_name_text(
                _member_annotation_value(member, source),
                source,
            )
            if normalized is None:
                continue
            key, display_text = normalized
            source_weight = _ORTHOGROUP_NAME_SOURCE_WEIGHTS[source]
            source_rank = _ORTHOGROUP_NAME_SOURCE_RANK[source]
            weak = _is_weak_orthogroup_name_candidate(display_text, source)
            accumulator = accumulators.get(key)
            if accumulator is None:
                accumulator = accumulators[key] = {
                    "text": display_text,
                    "source": source,
                    "source_weight": source_weight,
                    "source_rank": source_rank,
                    "members": set(),
                    "records": set(),
                    "representatives": set(),
                    "weak": weak,
                }
            if (
                source_weight > int(accumulator["source_weight"])
                or (
                    source_weight == int(accumulator["source_weight"])
                    and source_rank < int(accumulator["source_rank"])
                )
            ):
                accumulator["text"] = display_text
                accumulator["source"] = source
                accumulator["source_weight"] = source_weight
                accumulator["source_rank"] = source_rank
            if not weak:
                accumulator["weak"] = False
            members_set = accumulator["members"]
            records_set = accumulator["records"]
            representatives_set = accumulator["representatives"]
            if isinstance(members_set, set):
                members_set.add(member_id)
            if isinstance(records_set, set):
                records_set.add(int(member.record_index))
            if member.representative and isinstance(representatives_set, set):
                representatives_set.add(member_id)

    candidates: list[OrthogroupNameCandidate] = []
    for accumulator in accumulators.values():
        members_set = accumulator["members"]
        records_set = accumulator["records"]
        representatives_set = accumulator["representatives"]
        member_count = len(members_set) if isinstance(members_set, set) else 0
        record_coverage_count = len(records_set) if isinstance(records_set, set) else 0
        representative_count = len(representatives_set) if isinstance(representatives_set, set) else 0
        source_weight = int(accumulator["source_weight"])
        weak_penalty = 45 if bool(accumulator["weak"]) else 0
        score = (
            source_weight
            + member_count * 10
            + record_coverage_count * 20
            + representative_count * 5
            - weak_penalty
        )
        candidates.append(
            OrthogroupNameCandidate(
                text=str(accumulator["text"]),
                source=str(accumulator["source"]),
                member_count=member_count,
                record_coverage_count=record_coverage_count,
                representative_count=representative_count,
                score=float(score),
            )
        )

    return sorted(candidates, key=_orthogroup_name_candidate_sort_key)


def _orthogroup_name_candidate_sort_key(
    candidate: OrthogroupNameCandidate,
) -> tuple[int, float, int, int, int, str]:
    gene_priority = (
        candidate.source == "gene"
        and not _is_weak_orthogroup_name_candidate(candidate.text, candidate.source)
        and (
            candidate.record_coverage_count >= 2
            or candidate.representative_count > 0
        )
    )
    return (
        0 if gene_priority else 1,
        -candidate.score,
        _ORTHOGROUP_NAME_SOURCE_RANK.get(candidate.source, 99),
        -candidate.record_coverage_count,
        -candidate.member_count,
        candidate.text.casefold(),
    )


def _orthogroup_name_confidence(
    candidate: OrthogroupNameCandidate | None,
) -> str:
    if candidate is None:
        return "none"
    if _is_weak_orthogroup_name_candidate(candidate.text, candidate.source):
        return "low"
    if candidate.member_count <= 1:
        return "low"
    if candidate.record_coverage_count >= 2:
        return "high"
    return "medium"


def _orthogroup_description(
    candidate: OrthogroupNameCandidate | None,
    confidence: str,
    total_record_coverage: int,
) -> str:
    if candidate is None or confidence in {"none", "low"}:
        return "No informative product/gene/note consensus was found."
    record_word = "record" if total_record_coverage == 1 else "records"
    return (
        f"Suggested from {candidate.source} annotations in "
        f"{candidate.record_coverage_count} of {total_record_coverage} {record_word}."
    )


def _build_orthogroup_name_metadata(
    orthogroups: Mapping[str, Sequence[OrthogroupMember]],
) -> tuple[
    dict[str, str],
    dict[str, str],
    dict[str, list[OrthogroupNameCandidate]],
    dict[str, str],
]:
    names_by_orthogroup_id: dict[str, str] = {}
    descriptions_by_orthogroup_id: dict[str, str] = {}
    name_candidates_by_orthogroup_id: dict[str, list[OrthogroupNameCandidate]] = {}
    confidence_by_orthogroup_id: dict[str, str] = {}

    for orthogroup_id, members in orthogroups.items():
        annotation_candidates = _iter_orthogroup_name_candidates(members, include_label=False)
        candidates = annotation_candidates or _iter_orthogroup_name_candidates(members, include_label=True)
        top_candidate = candidates[0] if candidates else None
        total_record_coverage = len({int(member.record_index) for member in members})
        confidence = _orthogroup_name_confidence(top_candidate)

        if top_candidate is not None and confidence not in {"none", "low"}:
            names_by_orthogroup_id[orthogroup_id] = top_candidate.text
        descriptions_by_orthogroup_id[orthogroup_id] = _orthogroup_description(
            top_candidate,
            confidence,
            total_record_coverage,
        )
        name_candidates_by_orthogroup_id[orthogroup_id] = candidates
        confidence_by_orthogroup_id[orthogroup_id] = confidence

    return (
        names_by_orthogroup_id,
        descriptions_by_orthogroup_id,
        name_candidates_by_orthogroup_id,
        confidence_by_orthogroup_id,
    )



def _make_orthogroup_member(
    orthogroup_id: str,
    protein: CdsProtein,
    *,
    representative: bool,
    role: OrthogroupMemberRole = "anchor",
    confidence: OrthogroupMemberConfidence = "high",
    assignment_reason: str = "",
    supporting_edges: Sequence[str] = (),
    best_core_support: float = 0.0,
    second_best_core_support: float = 0.0,
) -> OrthogroupMember:
    return OrthogroupMember(
        orthogroup_id=orthogroup_id,
        protein_id=protein.protein_id,
        record_index=protein.record_index,
        feature_index=protein.feature_index,
        record_id=protein.record_id,
        label=protein.label,
        start=protein.start,
        end=protein.end,
        strand=protein.strand,
        feature_svg_id=protein.feature_svg_id,
        source_protein_id=protein.source_protein_id,
        gene=protein.gene,
        product=protein.product,
        note=protein.note,
        representative=representative,
        role=role,
        confidence=confidence,
        assignment_reason=assignment_reason,
        supporting_edges=tuple(str(edge_id) for edge_id in supporting_edges),
        best_core_support=float(best_core_support),
        second_best_core_support=float(second_best_core_support),
    )


def _orthogroup_result_from_member_ids(
    group_member_ids: Mapping[str, set[str]],
    representative_ids_by_group: Mapping[str, set[str]],
    protein_map: Mapping[str, CdsProtein],
    *,
    group_order: Sequence[str],
    rbh_orthogroups: Mapping[str, Sequence[str]],
    ortholog_edges_by_orthogroup_id: Mapping[str, Sequence[OrthologEdge]],
    path_indexes_by_orthogroup_id: Mapping[str, OrthologPathCollection],
    related_edges_by_orthogroup_id: Mapping[str, Sequence[OrthologEdge]],
    member_roles_by_protein_id: Mapping[str, OrthogroupMemberRole] | None = None,
    member_confidence_by_protein_id: Mapping[str, OrthogroupMemberConfidence] | None = None,
    assignment_reason_by_protein_id: Mapping[str, str] | None = None,
    supporting_edges_by_protein_id: Mapping[str, Sequence[str]] | None = None,
    best_core_support_by_protein_id: Mapping[str, float] | None = None,
    second_best_core_support_by_protein_id: Mapping[str, float] | None = None,
    scope_by_orthogroup_id: Mapping[str, OrthogroupScope] | None = None,
    source_record_index_by_orthogroup_id: Mapping[str, int] | None = None,
) -> OrthogroupGraphResult:
    orthogroups: dict[str, list[OrthogroupMember]] = {}
    member_by_protein_id: dict[str, OrthogroupMember] = {}
    ordered_group_ids = [
        group_id
        for group_id in group_order
        if group_id in group_member_ids and group_member_ids[group_id]
    ]
    remaining_group_ids = sorted(
        set(group_member_ids).difference(ordered_group_ids),
        key=lambda group_id: min(
            _protein_sort_key(protein_map[member_id])
            for member_id in group_member_ids[group_id]
            if member_id in protein_map
        ),
    )
    for group_id in [*ordered_group_ids, *remaining_group_ids]:
        member_ids = sorted(
            group_member_ids[group_id],
            key=lambda member_id: _protein_sort_key(protein_map[member_id]),
        )
        representative_ids = representative_ids_by_group.get(group_id, set())
        members: list[OrthogroupMember] = []
        for member_id in member_ids:
            protein = protein_map[member_id]
            member = _make_orthogroup_member(
                group_id,
                protein,
                representative=member_id in representative_ids,
                role=(member_roles_by_protein_id or {}).get(member_id, "anchor"),
                confidence=(member_confidence_by_protein_id or {}).get(member_id, "high"),
                assignment_reason=(assignment_reason_by_protein_id or {}).get(member_id, ""),
                supporting_edges=(supporting_edges_by_protein_id or {}).get(member_id, ()),
                best_core_support=(best_core_support_by_protein_id or {}).get(member_id, 0.0),
                second_best_core_support=(second_best_core_support_by_protein_id or {}).get(member_id, 0.0),
            )
            members.append(member)
            member_by_protein_id[member_id] = member
        orthogroups[group_id] = members

    (
        names_by_orthogroup_id,
        descriptions_by_orthogroup_id,
        name_candidates_by_orthogroup_id,
        confidence_by_orthogroup_id,
    ) = _build_orthogroup_name_metadata(orthogroups)

    return OrthogroupGraphResult(
        orthogroups=orthogroups,
        member_by_protein_id=member_by_protein_id,
        names_by_orthogroup_id=names_by_orthogroup_id,
        descriptions_by_orthogroup_id=descriptions_by_orthogroup_id,
        name_candidates_by_orthogroup_id=name_candidates_by_orthogroup_id,
        confidence_by_orthogroup_id=confidence_by_orthogroup_id,
        rbh_orthogroups={
            str(key): tuple(str(protein_id) for protein_id in value)
            for key, value in rbh_orthogroups.items()
        },
        ortholog_edges_by_orthogroup_id={
            str(key): tuple(value)
            for key, value in ortholog_edges_by_orthogroup_id.items()
        },
        path_indexes_by_orthogroup_id=dict(path_indexes_by_orthogroup_id),
        related_edges_by_orthogroup_id={
            str(key): tuple(value)
            for key, value in related_edges_by_orthogroup_id.items()
        },
        scope_by_orthogroup_id={
            group_id: (scope_by_orthogroup_id or {}).get(group_id, "cross_record")
            for group_id in orthogroups
        },
        source_record_index_by_orthogroup_id={
            str(group_id): int(source_record_index_by_orthogroup_id[group_id])
            for group_id in orthogroups
            if source_record_index_by_orthogroup_id is not None
            and group_id in source_record_index_by_orthogroup_id
        },
    )


def _make_ortholog_edge(
    *,
    orthogroup_id: str,
    source_rbh_orthogroup_id: str | None,
    target_rbh_orthogroup_id: str | None,
    query_id: str,
    subject_id: str,
    hit: _DirectedHit,
    protein_map: Mapping[str, CdsProtein],
    edge_kind: OrthologEdgeKind,
    render_role: OrthologRenderRole,
    path_id: str | None = None,
) -> OrthologEdge:
    query_protein = protein_map[query_id]
    subject_protein = protein_map[subject_id]
    return OrthologEdge(
        orthogroup_id=str(orthogroup_id),
        source_rbh_orthogroup_id=source_rbh_orthogroup_id,
        target_rbh_orthogroup_id=target_rbh_orthogroup_id,
        query_protein_id=str(query_id),
        subject_protein_id=str(subject_id),
        query_record_index=int(query_protein.record_index),
        subject_record_index=int(subject_protein.record_index),
        edge_kind=edge_kind,
        render_role=render_role,
        path_id=path_id,
        identity=hit.identity,
        evalue=hit.evalue,
        bitscore=hit.bitscore,
        alignment_length=int(hit.alignment_length),
    )


def _build_ortholog_path_indexes(edges_by_group, protein_map):
    updated, indexes = {}, {}
    for group_id, edges in edges_by_group.items():
        index = OrthologPathCollection.from_edges(group_id, edges, protein_map)
        updated[group_id] = tuple(replace(edge, path_id=index.first_path_id(edge)) for edge in edges)
        indexes[group_id] = index
    return updated, indexes


_DIRECTED_HIT_FLOAT_COLUMNS = ("normalized_score", "evalue", "min_coverage", "identity", "bitscore", "alignment_length")


def _directed_hits_from_tables(
    normalized_tables: Mapping[tuple[int, int], DataFrame],
    protein_map: Mapping[str, CdsProtein],
) -> _DirectedHitTable:
    """The best normalized hit per (query, subject) direction, in rank order.

    Self hits, hits to a protein outside ``protein_map`` and hits without a
    positive finite normalized score are left out. Of the hits in one direction
    the highest score wins, then the lowest e-value, highest minimum coverage,
    highest identity and longest alignment; the earliest wins a complete tie.
    """

    record_index_by_id = {protein_id: int(protein.record_index) for protein_id, protein in protein_map.items()}
    protein_ids = sorted(record_index_by_id, key=lambda protein_id: (record_index_by_id[protein_id], protein_id))
    protein_index = pd.Index(protein_ids, dtype=object)
    tables = [table for table in normalized_tables.values() if table is not None and not table.empty]
    query_codes = protein_index.get_indexer([str(value) for table in tables for value in table["query"].tolist()])
    subject_codes = protein_index.get_indexer([str(value) for table in tables for value in table["subject"].tolist()])
    # The empty first array gives the dtype when there are no tables.
    columns = {
        column: np.concatenate([np.empty(0), *(table[column].to_numpy(dtype=float) for table in tables)])
        for column in _DIRECTED_HIT_FLOAT_COLUMNS
    }
    domain_only_column = np.concatenate(
        [np.empty(0, dtype=bool), *(table["domain_only"].to_numpy(dtype=bool) for table in tables)]
    )
    offsets = np.cumsum([0, *(len(table) for table in tables)])
    sources = [
        _ComparisonColumns(table, *(columns[column][start:end] for column in ("identity", "bitscore", "alignment_length")))
        for table, start, end in zip(tables, offsets[:-1].tolist(), offsets[1:].tolist())
    ]
    # Only the arrays above are read from here on. On CPython 3.11 and later this
    # frame holds the only reference to a temporary argument, so dropping it frees
    # the tables before the rows are selected; older interpreters free them on return.
    del normalized_tables, tables
    raw_score = columns["normalized_score"]
    columns["normalized_score"] = np.where(np.isfinite(raw_score) & (raw_score > 0.0), raw_score, 0.0)
    known = (query_codes >= 0) & (subject_codes >= 0) & (query_codes != subject_codes)
    rows = np.flatnonzero(known & (columns["normalized_score"] > 0.0))

    count = len(protein_ids)
    query, subject = query_codes[rows].astype(np.int64), subject_codes[rows].astype(np.int64)
    # One group per direction, numbered in order of first appearance.
    directions = pd.factorize(query * count + subject)[0]
    score, evalue, min_coverage, identity, alignment_length = (
        columns[column][rows] for column in ("normalized_score", "evalue", "min_coverage", "identity", "alignment_length")
    )
    truncated_length = np.trunc(alignment_length)
    best = _first_minimum_per_group(
        directions,
        (-score, evalue, -min_coverage, -identity, -truncated_length),
        int(directions.max(initial=-1)) + 1,
    )
    # Record index then ID orders the endpoints; the codes keep that order.
    ranked = best[np.lexsort((
        query[best] * count + subject[best],
        -truncated_length[best],
        -identity[best],
        -min_coverage[best],
        evalue[best],
        -score[best],
    ))]
    ranked_rows = rows[ranked]
    query, subject, min_coverage = query[ranked], subject[ranked], min_coverage[ranked]
    domain_only = domain_only_column[ranked_rows] | (min_coverage <= _DOMAIN_ONLY_MAX_MIN_COVERAGE)
    table_of_row = np.searchsorted(offsets, ranked_rows, side="right") - 1
    return _DirectedHitTable(
        protein_ids=protein_ids,
        protein_index=protein_index,
        protein_record=np.asarray([record_index_by_id[protein_id] for protein_id in protein_ids], dtype=np.int64),
        query=query,
        subject=subject,
        score=score[ranked],
        evalue=evalue[ranked],
        min_coverage=min_coverage,
        domain_only=domain_only,
        supports_membership=(min_coverage >= _MIN_MEMBERSHIP_MIN_COVERAGE) & ~domain_only,
        reverse=np.asarray(pd.Index(query * count + subject).get_indexer(subject * count + query), dtype=np.int64),
        sources=sources,
        source=table_of_row.astype(np.int64),
        source_row=(ranked_rows - offsets[table_of_row]).astype(np.int64),
    )


def _comparison_values(hit: _DirectedHit) -> tuple[object, ...]:
    """The hit's comparison column values, in COMPARISON_COLUMNS order with the IDs as strings."""

    return (hit.query, hit.subject, *hit.comparison_source.values(hit.comparison_row))


def _comparison_record_for_ids(values: Sequence[object], query_id: str, subject_id: str) -> dict[str, object]:
    """A comparison row for query_id -> subject_id from values in COMPARISON_COLUMNS order."""

    record = dict(zip(COMPARISON_COLUMNS, values))
    if str(values[0]) == subject_id and str(values[1]) == query_id:
        record["qstart"], record["sstart"] = record["sstart"], record["qstart"]
        record["qend"], record["send"] = record["send"], record["qend"]
    record["query"] = query_id
    record["subject"] = subject_id
    return record


def _anchor_core_edge_sort_key(
    edge: _AnchorCoreEvidenceEdge,
    protein_map: Mapping[str, CdsProtein],
) -> tuple[tuple[int, int, int, str], tuple[int, int, int, str], int, float, float, str, str]:
    return (
        _protein_sort_key(protein_map[edge.query_id]),
        _protein_sort_key(protein_map[edge.subject_id]),
        0 if edge.edge_kind == "rbh" else 1,
        -max(float(edge.score), float(edge.reverse_score)),
        edge.hit.evalue,
        edge.query_id,
        edge.subject_id,
    )


def _group_codes(hits: _DirectedHitTable, group_by_protein: Mapping[str, str]) -> tuple[list[str], np.ndarray]:
    """Group IDs in order of first appearance, and each protein code's group index (-1 if none)."""

    group_ids = list(dict.fromkeys(group_by_protein.values()))
    index_by_group = {group_id: index for index, group_id in enumerate(group_ids)}
    codes = hits.codes(group_by_protein)
    groups = np.asarray([index_by_group[group_id] for group_id in group_by_protein.values()], dtype=np.int64)
    group_of = np.full(len(hits.protein_ids), -1, dtype=np.int64)
    group_of[codes[codes >= 0]] = groups[codes >= 0]
    return group_ids, group_of


def _first_row_per_pair(hits: _DirectedHitTable, rows: np.ndarray) -> np.ndarray:
    """Of the given ascending rows, the first of each unordered protein pair."""

    count = len(hits.protein_ids)
    low = np.minimum(hits.query[rows], hits.subject[rows])
    high = np.maximum(hits.query[rows], hits.subject[rows])
    return rows[_first_occurrences(low * count + high)[1]]


def _evidence_edges(
    hits: _DirectedHitTable,
    protein_map: Mapping[str, CdsProtein],
    rows: np.ndarray,
    kinds: Iterable[Literal["rbh", "coortholog", "record_local_paralog"]],
) -> tuple[_AnchorCoreEvidenceEdge, ...]:
    """Edges for the rows' pairs between canonical endpoints (the lower code first), in protein order."""

    reverse = hits.reverse[rows]
    forward = hits.query[rows] < hits.subject[rows]
    edge_rows = np.where(forward, rows, reverse)
    reverse_rows = np.where(forward, reverse, rows)
    edges = [
        _AnchorCoreEvidenceEdge(
            query_id=hits.protein_ids[query],
            subject_id=hits.protein_ids[subject],
            hit=hits.hit(row),
            score=score,
            reverse_score=reverse_score,
            edge_kind=kind,
        )
        for query, subject, row, score, reverse_score, kind in zip(
            hits.query[edge_rows].tolist(),
            hits.subject[edge_rows].tolist(),
            edge_rows.tolist(),
            hits.score[edge_rows].tolist(),
            hits.score[reverse_rows].tolist(),
            kinds,
        )
    ]
    return tuple(sorted(edges, key=lambda edge: _anchor_core_edge_sort_key(edge, protein_map)))


def _select_anchor_core_edges(
    hits: _DirectedHitTable,
    protein_map: Mapping[str, CdsProtein],
) -> tuple[_AnchorCoreEvidenceEdge, ...]:
    query, subject, score = hits.query, hits.subject, hits.score
    query_record, subject_record = hits.protein_record[query], hits.protein_record[subject]
    cross = np.flatnonzero(query_record != subject_record)
    # A query's best hit into another record is its first there in rank order.
    # Dense record codes keep the key distinct for any record index.
    record_code = pd.factorize(hits.protein_record)[0]
    record_count = int(record_code.max(initial=-1)) + 1
    target, first = _first_occurrences(query[cross] * record_count + record_code[subject[cross]])
    best_in_target = np.full(len(hits), -1, dtype=np.int64)
    best_in_target[cross] = cross[first][target]
    # Every test below is symmetric in the two directions, and both give the
    # same candidate, so the earlier direction decides the pair.
    rows = _first_row_per_pair(hits, cross)
    rows = rows[hits.supports_membership[rows] & (hits.reverse[rows] >= 0)]
    reverse = hits.reverse[rows]
    rows, reverse = rows[hits.supports_membership[reverse]], reverse[hits.supports_membership[reverse]]
    # Both directions are cross-record hits, so each direction has a best hit.
    query_best, subject_best = best_in_target[rows], best_in_target[reverse]
    keep = ~(score[rows] < _NEAR_RECIPROCAL_MIN_RATIO * score[query_best]) & ~(
        score[reverse] < _NEAR_RECIPROCAL_MIN_RATIO * score[subject_best]
    )
    rows, query_best, subject_best = rows[keep], query_best[keep], subject_best[keep]
    strict_rbh = (subject[query_best] == subject[rows]) & (subject[subject_best] == query[rows])
    return _evidence_edges(
        hits, protein_map, rows, ("rbh" if strict else "coortholog" for strict in strict_rbh.tolist())
    )


def _select_record_local_paralog_edges(
    hits: _DirectedHitTable,
    protein_map: Mapping[str, CdsProtein],
    thresholds: Mapping[str, float],
    group_by_protein: Mapping[str, str],
) -> tuple[_AnchorCoreEvidenceEdge, ...]:
    query, subject, score = hits.query, hits.subject, hits.score
    same = np.flatnonzero(hits.protein_record[query] == hits.protein_record[subject])
    _, group_of = _group_codes(hits, group_by_protein)
    unassigned = (group_of[query] < 0) & (group_of[subject] < 0)
    # Each query's best same-record score among hits between unassigned proteins.
    best_same_score = np.zeros(len(hits.protein_ids))
    free = same[unassigned[same]]
    np.maximum.at(best_same_score, query[free], score[free])
    threshold_score = np.zeros(len(hits.protein_ids))
    has_threshold = np.zeros(len(hits.protein_ids), dtype=bool)
    threshold_codes = hits.codes(thresholds)
    mapped = threshold_codes >= 0
    threshold_score[threshold_codes[mapped]] = np.fromiter(thresholds.values(), dtype=float, count=len(thresholds))[mapped]
    has_threshold[threshold_codes[mapped]] = True

    def passes_score_floor(protein: np.ndarray, value: np.ndarray) -> np.ndarray:
        best = best_same_score[protein]
        fallback = (best > 0.0) & (value >= _RECORD_LOCAL_NEAR_RECIPROCAL_MIN_RATIO * best)
        return np.where(has_threshold[protein], value >= threshold_score[protein], fallback)

    supports_local = hits.supports_membership & (hits.min_coverage >= _RECORD_LOCAL_MIN_MEMBERSHIP_MIN_COVERAGE)
    # As for anchors: symmetric tests and one candidate per pair.
    rows = _first_row_per_pair(hits, same)
    rows = rows[unassigned[rows] & supports_local[rows] & (hits.reverse[rows] >= 0)]
    reverse = hits.reverse[rows]
    keep = (
        supports_local[reverse]
        & passes_score_floor(query[rows], score[rows])
        & passes_score_floor(subject[rows], score[reverse])
    )
    return _evidence_edges(hits, protein_map, rows[keep], itertools.repeat("record_local_paralog"))


def _record_local_components_from_edges(
    edges: Sequence[_AnchorCoreEvidenceEdge],
    protein_map: Mapping[str, CdsProtein],
) -> tuple[tuple[str, ...], ...]:
    union_find = _UnionFind()
    edge_nodes: set[str] = set()
    for edge in edges:
        if edge.query_id not in protein_map or edge.subject_id not in protein_map:
            continue
        union_find.union(edge.query_id, edge.subject_id)
        edge_nodes.update((edge.query_id, edge.subject_id))
    components_by_root: dict[str, set[str]] = {}
    for protein_id in edge_nodes:
        components_by_root.setdefault(union_find.find(protein_id), set()).add(protein_id)

    components: list[tuple[str, ...]] = []
    for member_ids in components_by_root.values():
        if len(member_ids) < _RECORD_LOCAL_MIN_MEMBERS:
            continue
        record_indices = {int(protein_map[member_id].record_index) for member_id in member_ids}
        if len(record_indices) != 1:
            continue
        components.append(
            tuple(
                sorted(
                    member_ids,
                    key=lambda member_id: _protein_sort_key(protein_map[member_id]),
                )
            )
        )
    return tuple(
        sorted(
            components,
            key=lambda component: min(_protein_sort_key(protein_map[member_id]) for member_id in component),
        )
    )


def _record_local_support_by_member(
    edges: Sequence[_AnchorCoreEvidenceEdge],
) -> dict[str, float]:
    support_by_member: dict[str, float] = {}
    for edge in edges:
        support_by_member[edge.query_id] = max(
            support_by_member.get(edge.query_id, 0.0),
            float(edge.score),
        )
        support_by_member[edge.subject_id] = max(
            support_by_member.get(edge.subject_id, 0.0),
            float(edge.reverse_score),
        )
    return support_by_member


def _record_local_component_has_competing_core_support(
    member_ids: Sequence[str],
    has_core_groups: bool,
    evidence_by_protein: Mapping[str, Mapping[str, Sequence[_DirectedHit]]],
    thresholds: Mapping[str, float],
    local_support_by_member: Mapping[str, float],
) -> bool:
    if not has_core_groups:
        return False
    credible_core_ids: set[str] = set()
    for member_id in member_ids:
        member_local_support = float(local_support_by_member.get(member_id, 0.0))
        if member_local_support <= 0.0:
            return True
        candidates = [
            candidate
            for group_id, evidence in evidence_by_protein.get(member_id, {}).items()
            if (
                candidate := _build_core_support_candidate(
                    member_id,
                    group_id,
                    evidence,
                    thresholds,
                )
            )
            is not None
        ]
        if not candidates:
            continue
        strongest = max(
            candidates,
            key=lambda candidate: (
                max(
                    float(candidate.same_record_score),
                    float(candidate.cross_record_score),
                    float(candidate.diagnostic_score),
                ),
                float(candidate.diagnostic_score),
                candidate.group_id,
            ),
        )
        strongest_score = max(
            float(strongest.same_record_score),
            float(strongest.cross_record_score),
            float(strongest.diagnostic_score),
        )
        competition_floor = member_local_support / _RECORD_LOCAL_CORE_COMPETITION_RATIO
        if strongest.domain_only and strongest_score >= competition_floor:
            return True
        for candidate in candidates:
            core_support = max(
                float(candidate.cross_record_score),
                float(candidate.same_record_score),
            )
            if core_support <= 0.0:
                continue
            if candidate.high_confidence_pass:
                return True
            if candidate.domain_only:
                if core_support >= competition_floor:
                    return True
                continue
            if candidate.low_confidence_pass and core_support >= competition_floor:
                credible_core_ids.add(candidate.group_id)
    return len(credible_core_ids) >= 2


def _next_orthogroup_index(group_ids: Sequence[str]) -> int:
    max_index = 0
    for group_id in group_ids:
        match = re.fullmatch(r"og_(\d+)", str(group_id))
        if match is not None:
            max_index = max(max_index, int(match.group(1)))
    return max_index + 1


def _append_record_local_orthogroups(
    *,
    group_member_ids: dict[str, set[str]],
    group_by_protein: dict[str, str],
    representative_ids_by_group: dict[str, set[str]],
    member_roles: dict[str, OrthogroupMemberRole],
    member_confidence: dict[str, OrthogroupMemberConfidence],
    assignment_reasons: dict[str, str],
    supporting_edges: dict[str, tuple[str, ...]],
    best_core_support: dict[str, float],
    second_core_support: dict[str, float],
    ortholog_edges_by_group: dict[str, list[OrthologEdge]],
    scope_by_group: dict[str, OrthogroupScope],
    source_record_index_by_group: dict[str, int],
    group_order: list[str],
    components: Sequence[Sequence[str]],
    accepted_edges: Sequence[_AnchorCoreEvidenceEdge],
    local_support_by_member: Mapping[str, float],
    protein_map: Mapping[str, CdsProtein],
) -> None:
    next_group_index = _next_orthogroup_index(group_order or list(group_member_ids))
    for component in components:
        component_members = tuple(
            sorted(
                (member_id for member_id in component if member_id in protein_map),
                key=lambda member_id: _protein_sort_key(protein_map[member_id]),
            )
        )
        if len(component_members) < _RECORD_LOCAL_MIN_MEMBERS:
            continue
        component_set = set(component_members)
        record_indices = {int(protein_map[member_id].record_index) for member_id in component_members}
        if len(record_indices) != 1:
            continue
        group_id = f"og_{next_group_index}"
        next_group_index += 1
        group_order.append(group_id)
        group_member_ids[group_id] = set(component_members)
        scope_by_group[group_id] = "record_local"
        source_record_index_by_group[group_id] = next(iter(record_indices))
        for member_id in component_members:
            group_by_protein[member_id] = group_id
            member_roles[member_id] = "local_paralog"
            member_confidence[member_id] = "high"
            assignment_reasons[member_id] = "record-local reciprocal paralog cluster"
            best_core_support[member_id] = min(
                float(local_support_by_member.get(component_member_id, 0.0))
                for component_member_id in component_members
            )
            second_core_support[member_id] = 0.0

        representative_id = min(
            component_members,
            key=lambda member_id: (
                -float(local_support_by_member.get(member_id, 0.0)),
                _protein_sort_key(protein_map[member_id]),
            ),
        )
        representative_ids_by_group[group_id] = {representative_id}

        component_edges = [
            edge
            for edge in accepted_edges
            if edge.query_id in component_set and edge.subject_id in component_set
        ]
        for edge in component_edges:
            edge_id = f"{edge.query_id}->{edge.subject_id}:record_local_paralog"
            for member_id in (edge.query_id, edge.subject_id):
                supporting_edges[member_id] = (*supporting_edges.get(member_id, ()), edge_id)
            ortholog_edges_by_group.setdefault(group_id, []).append(
                _make_ortholog_edge(
                    orthogroup_id=group_id,
                    source_rbh_orthogroup_id=group_id,
                    target_rbh_orthogroup_id=group_id,
                    query_id=edge.query_id,
                    subject_id=edge.subject_id,
                    hit=edge.hit,
                    protein_map=protein_map,
                    edge_kind="record_local_paralog",
                    render_role="display_edge",
                )
            )


def _derive_anchor_core_thresholds(
    hits: _DirectedHitTable,
    anchor_edges: Sequence[_AnchorCoreEvidenceEdge],
) -> dict[str, float]:
    """Each protein's same-record score threshold, by protein ID.

    A fixed share of the protein's lowest positive, finite anchor score, or of
    its best cross-record score when it has none; absent without either.
    """
    cross = hits.protein_record[hits.query] != hits.protein_record[hits.subject]
    cross_best = np.full(len(hits.protein_ids), -np.inf)
    np.maximum.at(cross_best, hits.query[cross], hits.score[cross])
    has_cross = np.bincount(hits.query[cross], minlength=len(hits.protein_ids)) > 0

    endpoint = hits.codes([*(edge.query_id for edge in anchor_edges), *(edge.subject_id for edge in anchor_edges)])
    endpoint_score = np.asarray(
        [*(edge.score for edge in anchor_edges), *(edge.reverse_score for edge in anchor_edges)], dtype=float
    )
    usable = (endpoint >= 0) & (endpoint_score > 0.0) & np.isfinite(endpoint_score)
    anchor_floor = np.full(len(hits.protein_ids), np.inf)
    np.minimum.at(anchor_floor, endpoint[usable], endpoint_score[usable])
    has_anchor = anchor_floor < np.inf
    floor = np.where(has_anchor, anchor_floor, cross_best) * _INPARALOG_MIN_SAME_RECORD_RATIO_TO_LOCAL_ANCHOR
    codes = np.flatnonzero(has_anchor | has_cross)
    return dict(zip((hits.protein_ids[code] for code in codes.tolist()), floor[codes].tolist()))


def _index_core_support_evidence(
    hits: _DirectedHitTable,
    group_by_protein: Mapping[str, str],
) -> dict[str, dict[str, list[_DirectedHit]]]:
    """Snapshot incoming/outgoing evidence from unassigned proteins to members.

    Only hits with exactly one assigned endpoint become records. Buckets never
    alias the mutable membership map; rebuild after membership expansion.
    """
    group_ids, group_of = _group_codes(hits, group_by_protein)
    query_group, subject_group = group_of[hits.query], group_of[hits.subject]
    rows = np.flatnonzero((query_group < 0) != (subject_group < 0))
    evidence_by_protein: dict[str, dict[str, list[_DirectedHit]]] = {}
    for row, query_index, subject_index in zip(rows.tolist(), query_group[rows].tolist(), subject_group[rows].tolist()):
        hit = hits.hit(row)
        if query_index < 0:
            evidence_by_protein.setdefault(hit.query, {}).setdefault(group_ids[subject_index], []).append(hit)
        else:
            evidence_by_protein.setdefault(hit.subject, {}).setdefault(group_ids[query_index], []).append(hit)
    return evidence_by_protein


def _best_evidence_between_protein_and_members(
    evidence: Sequence[_DirectedHit],
) -> tuple[_BestCoreEvidence, _BestCoreEvidence]:
    """Reduce connected evidence once, returning same-record then cross-record.

    The higher score wins, and the earlier hit in rank order on equal scores.
    """
    supports: list[_DirectedHit | None] = [None, None]
    diagnostics: list[_DirectedHit | None] = [None, None]
    for hit in evidence:
        record_kind = int(hit.query_record_index != hit.subject_record_index)
        diagnostic = diagnostics[record_kind]
        if diagnostic is None or hit.score > diagnostic.score or (hit.score == diagnostic.score and hit.rank < diagnostic.rank):
            diagnostics[record_kind] = hit
        if not hit.supports_membership:
            continue
        support = supports[record_kind]
        if support is None or hit.score > support.score or (hit.score == support.score and hit.rank < support.rank):
            supports[record_kind] = hit
    return (
        _BestCoreEvidence(support=supports[0], diagnostic=diagnostics[0]),
        _BestCoreEvidence(support=supports[1], diagnostic=diagnostics[1]),
    )


def _support_gap_ratio(best_support: float, second_support: float) -> float:
    if best_support <= 0.0:
        return 0.0
    if second_support <= _ORTHOGROUP_SPLIT_EPSILON:
        return float("inf")
    return float(best_support) / max(float(second_support), _ORTHOGROUP_SPLIT_EPSILON)


def _hit_score(hit: _DirectedHit | None) -> float:
    return hit.score if hit is not None else 0.0


def _build_core_support_candidate(
    protein_id: str,
    group_id: str,
    evidence: Sequence[_DirectedHit],
    thresholds: Mapping[str, float],
) -> _CoreSupportCandidate | None:
    same_evidence, cross_evidence = _best_evidence_between_protein_and_members(evidence)
    if same_evidence.diagnostic is None and cross_evidence.diagnostic is None:
        return None

    same_score = _hit_score(same_evidence.support)
    cross_score = _hit_score(cross_evidence.support)
    support = float(cross_score) + 0.5 * float(same_score)
    same_diagnostic_score = _hit_score(same_evidence.diagnostic)
    cross_diagnostic_score = _hit_score(cross_evidence.diagnostic)
    diagnostic_score = max(same_diagnostic_score, cross_diagnostic_score)

    candidate_threshold = thresholds.get(protein_id)
    same_member_id = ""
    if same_evidence.support is not None:
        same_member_id = (
            same_evidence.support.subject
            if same_evidence.support.query == protein_id
            else same_evidence.support.query
        )
    same_member_threshold = thresholds.get(same_member_id)
    same_pass = (
        same_evidence.support is not None
        and (
            (
                candidate_threshold is not None
                and same_score >= candidate_threshold
            )
            or (
                same_member_threshold is not None
                and same_score >= same_member_threshold
            )
        )
    )
    cross_pass = (
        cross_evidence.support is not None
        and candidate_threshold is not None
        and cross_score >= candidate_threshold
    )

    evidence_hit: _DirectedHit | None
    relation_kind: OrthologEdgeKind
    if same_pass:
        evidence_hit, relation_kind = same_evidence.support, "same_record_inparalog"
    elif cross_pass:
        evidence_hit, relation_kind = cross_evidence.support, "coortholog"
    elif same_evidence.support is not None and (
        cross_evidence.support is None or same_score >= cross_score
    ):
        evidence_hit, relation_kind = same_evidence.support, "same_record_inparalog"
    elif cross_evidence.support is not None:
        evidence_hit, relation_kind = cross_evidence.support, "coortholog"
    elif same_evidence.diagnostic is not None and (
        cross_evidence.diagnostic is None
        or same_diagnostic_score >= cross_diagnostic_score
    ):
        evidence_hit, relation_kind = same_evidence.diagnostic, "same_record_inparalog"
    else:
        evidence_hit, relation_kind = cross_evidence.diagnostic, "coortholog"
    if evidence_hit is None:
        return None

    min_coverage = evidence_hit.min_coverage
    domain_only = evidence_hit.domain_only
    has_support_row = same_evidence.support is not None or cross_evidence.support is not None
    high_confidence_pass = bool(has_support_row and (same_pass or cross_pass)) and min_coverage >= _MIN_MEMBERSHIP_MIN_COVERAGE and not domain_only
    low_confidence_pass = bool(has_support_row) and min_coverage >= _MIN_MEMBERSHIP_MIN_COVERAGE and not domain_only
    if same_pass:
        reason = "same-record support passed local anchor threshold"
    elif cross_pass:
        reason = "cross-record support passed local fallback threshold"
    elif domain_only:
        reason = "support is domain-only"
    elif min_coverage < _MIN_MEMBERSHIP_MIN_COVERAGE:
        reason = "support is below membership coverage"
    else:
        reason = "support is separated but below high-confidence threshold"
    return _CoreSupportCandidate(
        group_id=group_id,
        support=support,
        diagnostic_score=diagnostic_score,
        same_record_score=float(same_score),
        cross_record_score=float(cross_score),
        evidence=evidence_hit,
        relation_kind=relation_kind,
        min_coverage=float(min_coverage),
        domain_only=bool(domain_only),
        high_confidence_pass=high_confidence_pass,
        low_confidence_pass=low_confidence_pass,
        reason=reason,
    )


def _core_support_sort_key(candidate: _CoreSupportCandidate) -> tuple[float, float, str]:
    return (-float(candidate.support), -float(candidate.diagnostic_score), candidate.group_id)


def _anchor_core_edge_tables(
    anchor_edges: Sequence[_AnchorCoreEvidenceEdge],
    protein_map: Mapping[str, CdsProtein],
) -> dict[tuple[int, int], DataFrame]:
    rows_by_pair: dict[tuple[int, int], list[dict[str, object]]] = {}
    for edge in anchor_edges:
        query_record_index = int(protein_map[edge.query_id].record_index)
        subject_record_index = int(protein_map[edge.subject_id].record_index)
        if query_record_index == subject_record_index:
            continue
        pair = (min(query_record_index, subject_record_index), max(query_record_index, subject_record_index))
        query_id, subject_id = edge.query_id, edge.subject_id
        if query_record_index > subject_record_index:
            query_id, subject_id = subject_id, query_id
        rows_by_pair.setdefault(pair, []).append(_comparison_record_for_ids(_comparison_values(edge.hit), query_id, subject_id))
    return {
        pair: pd.DataFrame.from_records(rows, columns=COMPARISON_COLUMNS)
        for pair, rows in rows_by_pair.items()
    }


_RELATED_EDGE_KINDS: tuple[OrthologEdgeKind, ...] = ("related_homolog", "same_record_inparalog", "weak_bridge", "domain_only")


def _append_related_edges(
    hits: _DirectedHitTable,
    protein_map: Mapping[str, CdsProtein],
    *,
    group_by_protein: Mapping[str, str],
    member_roles: Mapping[str, OrthogroupMemberRole],
    ortholog_edges_by_group: Mapping[str, Sequence[OrthologEdge]],
    related_edges_by_group: dict[str, list[OrthologEdge]],
    max_related_edges_per_orthogroup: int,
) -> None:
    """Fill each group's related edges from hits in rank order, up to the cap.

    The first hit of a protein pair and relation kind counts, unless a membership
    edge joins the pair. The group of the canonical query, then the subject's,
    takes its edge unless the group holds an equal edge (in either direction) or
    is full. Groups never share room, so each takes its first fresh edges. Every
    group in ``related_edges_by_group`` has members in ``group_by_protein``.
    """

    count, kinds = len(hits.protein_ids), len(_RELATED_EDGE_KINDS)
    group_ids, group_of = _group_codes(hits, group_by_protein)
    low, high = np.minimum(hits.query, hits.subject), np.maximum(hits.query, hits.subject)
    low_group, high_group = group_of[low], group_of[high]
    inparalog = np.zeros(count, dtype=bool)
    inparalog_codes = hits.codes(protein_id for protein_id, role in member_roles.items() if role == "inparalog")
    inparalog[inparalog_codes[inparalog_codes >= 0]] = True
    # Positions in _RELATED_EDGE_KINDS; the first condition that holds decides.
    kind = np.select(
        [
            hits.domain_only,
            (low_group >= 0) & (high_group >= 0) & (low_group != high_group),
            (low_group == high_group) & (inparalog[low] | inparalog[high]),
        ],
        [3, 2, 1],
        0,
    )
    pair = low * count + high
    membership_query = hits.codes(edge.query_protein_id for edges in ortholog_edges_by_group.values() for edge in edges)
    membership_subject = hits.codes(edge.subject_protein_id for edges in ortholog_edges_by_group.values() for edge in edges)
    membership_pairs = np.minimum(membership_query, membership_subject) * count + np.maximum(membership_query, membership_subject)
    _, rows = _first_occurrences(pair * kinds + kind)
    rows = rows[~np.isin(pair[rows], membership_pairs)]
    # One incidence per row and group, keyed by the row's position in rows.
    row_low_group, row_high_group = low_group[rows], high_group[rows]
    on_low = row_low_group >= 0
    on_high = (row_high_group >= 0) & (row_high_group != row_low_group)
    position = np.concatenate([np.flatnonzero(on_low), np.flatnonzero(on_high)])
    group = np.concatenate([row_low_group[on_low], row_high_group[on_high]])
    # The stable sort keeps the query's group before the subject's.
    order = np.argsort(position, kind="stable")
    position, group = position[order], group[order]
    # An edge a group holds uses its room; an equal incidence is not fresh.
    index_by_group = {group_id: index for index, group_id in enumerate(group_ids)}
    room = np.full(len(group_ids), max_related_edges_per_orthogroup, dtype=np.int64)
    held_group: list[int] = []
    held_key: list[int] = []
    for group_id, edges in related_edges_by_group.items():
        room[index_by_group[group_id]] -= len(edges)
        for edge in edges:
            if edge.edge_kind in _RELATED_EDGE_KINDS:
                low_code, high_code = sorted(hits.codes((edge.query_protein_id, edge.subject_protein_id)).tolist())
                held_group.append(index_by_group[group_id])
                held_key.append((low_code * count + high_code) * kinds + _RELATED_EDGE_KINDS.index(edge.edge_kind))
    held_position = pd.Index(pair[rows] * kinds + kind[rows]).get_indexer(held_key)
    found = held_position >= 0
    held = np.asarray(held_group, dtype=np.int64)[found] * len(rows) + held_position[found]
    fresh = ~np.isin(group * len(rows) + position, held)
    position, group = position[fresh], group[fresh]
    accepted = pd.Series(group).groupby(group).cumcount().to_numpy() < room[group]
    for row, group_index in zip(rows[position[accepted]].tolist(), group[accepted].tolist()):
        query_group, subject_group = low_group.item(row), high_group.item(row)
        group_id = group_ids[group_index]
        related_edges_by_group.setdefault(group_id, []).append(
            _make_ortholog_edge(
                orthogroup_id=group_id,
                source_rbh_orthogroup_id=group_ids[query_group] if query_group >= 0 else None,
                target_rbh_orthogroup_id=group_ids[subject_group] if subject_group >= 0 else None,
                query_id=hits.protein_ids[low.item(row)],
                subject_id=hits.protein_ids[high.item(row)],
                hit=hits.hit(row),
                protein_map=protein_map,
                edge_kind=_RELATED_EDGE_KINDS[kind.item(row)],
                render_role="display_edge",
            )
        )


def _build_anchor_core_orthogroups(
    hits: _DirectedHitTable,
    anchor_edges: Sequence[_AnchorCoreEvidenceEdge],
    protein_map: Mapping[str, CdsProtein],
    *,
    include_singletons: bool,
    max_related_edges_per_orthogroup: int,
) -> OrthogroupGraphResult:
    thresholds = _derive_anchor_core_thresholds(hits, anchor_edges)
    union_find = _UnionFind()
    member_ranks: dict[str, tuple[float, float, float, float]] = {}
    member_roles: dict[str, OrthogroupMemberRole] = {}
    member_confidence: dict[str, OrthogroupMemberConfidence] = {}
    assignment_reasons: dict[str, str] = {}
    supporting_edges: dict[str, tuple[str, ...]] = {}
    best_core_support: dict[str, float] = {}
    second_core_support: dict[str, float] = {}
    ortholog_edges_by_group: dict[str, list[OrthologEdge]] = {}
    scope_by_group: dict[str, OrthogroupScope] = {}
    source_record_index_by_group: dict[str, int] = {}

    for edge in anchor_edges:
        union_find.union(edge.query_id, edge.subject_id)
        rank = (edge.hit.bitscore, edge.hit.evalue, edge.hit.identity, edge.hit.alignment_length)
        if _is_better_member_rank(rank, member_ranks.get(edge.query_id)):
            member_ranks[edge.query_id] = rank
        if _is_better_member_rank(rank, member_ranks.get(edge.subject_id)):
            member_ranks[edge.subject_id] = rank
        for protein_id in (edge.query_id, edge.subject_id):
            current_role = member_roles.get(protein_id)
            if edge.edge_kind == "rbh" or current_role is None:
                member_roles[protein_id] = "anchor" if edge.edge_kind == "rbh" else "coortholog"
            member_confidence[protein_id] = "high"
            assignment_reasons[protein_id] = (
                "cross-record reciprocal best anchor"
                if edge.edge_kind == "rbh"
                else "cross-record near-reciprocal anchor"
            )
            supporting_edges[protein_id] = (*supporting_edges.get(protein_id, ()), f"{edge.query_id}->{edge.subject_id}:{edge.edge_kind}")

    components: dict[str, set[str]] = {}
    for edge in anchor_edges:
        root = union_find.find(edge.query_id)
        components.setdefault(root, set()).update((edge.query_id, edge.subject_id))
    if include_singletons:
        # anchor_core_v1 keeps isolated proteins unassigned; record-local
        # orthogroups are added only from accepted same-record evidence below.
        pass

    sorted_components = sorted(
        (member_ids for member_ids in components.values() if member_ids),
        key=lambda member_ids: min(_protein_sort_key(protein_map[member_id]) for member_id in member_ids),
    )
    group_member_ids: dict[str, set[str]] = {
        f"og_{index}": set(member_ids)
        for index, member_ids in enumerate(sorted_components, start=1)
    }
    group_order = list(group_member_ids)
    group_by_protein: dict[str, str] = {
        protein_id: group_id
        for group_id, member_ids in group_member_ids.items()
        for protein_id in member_ids
    }
    scope_by_group.update({group_id: "cross_record" for group_id in group_member_ids})

    def representative_member_key(member_id: str) -> tuple[float, float, float, float, str]:
        rank = member_ranks[member_id]
        return (-rank[0], rank[1], -rank[2], -rank[3], str(member_id))

    representative_ids_by_group: dict[str, set[str]] = {}
    for group_id, member_ids in group_member_ids.items():
        members_by_record: dict[int, list[str]] = {}
        for member_id in member_ids:
            protein = protein_map[member_id]
            members_by_record.setdefault(int(protein.record_index), []).append(member_id)
        for record_member_ids in members_by_record.values():
            representative_id = min(
                record_member_ids,
                key=representative_member_key,
            )
            representative_ids_by_group.setdefault(group_id, set()).add(representative_id)

    for edge in anchor_edges:
        edge_group_id = group_by_protein.get(edge.query_id)
        if edge_group_id is None:
            continue
        ortholog_edges_by_group.setdefault(edge_group_id, []).append(
            _make_ortholog_edge(
                orthogroup_id=edge_group_id,
                source_rbh_orthogroup_id=edge_group_id,
                target_rbh_orthogroup_id=edge_group_id,
                query_id=edge.query_id,
                subject_id=edge.subject_id,
                hit=edge.hit,
                protein_map=protein_map,
                edge_kind=edge.edge_kind,
                render_role="block_anchor",
            )
        )

    # Freeze core membership before any additions; evidence cannot chain through
    # proteins assigned later in this loop.
    core_member_snapshot = _index_core_support_evidence(hits, group_by_protein)
    related_edges_by_group: dict[str, list[OrthologEdge]] = {}
    for protein_id in sorted(protein_map, key=lambda item: _protein_sort_key(protein_map[item])):
        if protein_id in group_by_protein:
            continue
        candidates = [
            candidate
            for group_id, evidence in core_member_snapshot.get(protein_id, {}).items()
            if (
                candidate := _build_core_support_candidate(
                    protein_id,
                    group_id,
                    evidence,
                    thresholds,
                )
            )
            is not None
        ]
        if not candidates:
            continue
        candidates.sort(key=_core_support_sort_key)
        best = candidates[0]
        second = candidates[1] if len(candidates) > 1 else None
        second_support = float(second.support) if second is not None else 0.0
        gap_ratio = _support_gap_ratio(best.support, second_support)
        assigned_role: OrthogroupMemberRole | None = None
        assigned_confidence: OrthogroupMemberConfidence | None = None
        assigned_reason = best.reason
        if best.high_confidence_pass and gap_ratio >= _MIN_BEST_SECOND_CORE_RATIO:
            assigned_role = "inparalog" if best.relation_kind == "same_record_inparalog" else "coortholog"
            assigned_confidence = "high"
        elif best.low_confidence_pass and gap_ratio >= _LOW_CONFIDENCE_MIN_BEST_SECOND_CORE_RATIO:
            assigned_role = "low_confidence"
            assigned_confidence = "low"
            assigned_reason = f"low-confidence assignment: {best.reason}"

        if assigned_role is None or assigned_confidence is None:
            if best.domain_only:
                unassigned_relation_kind: OrthologEdgeKind = "domain_only"
            elif second is not None and best.low_confidence_pass:
                unassigned_relation_kind = "ambiguous_paralog"
            else:
                unassigned_relation_kind = best.relation_kind
            if unassigned_relation_kind in {"ambiguous_paralog", "domain_only"}:
                related_edges = related_edges_by_group.setdefault(best.group_id, [])
                # Each unassigned protein adds at most one edge, to a core member,
                # so no edge here repeats another.
                if len(related_edges) < max_related_edges_per_orthogroup:
                    related_edges.append(
                        _make_ortholog_edge(
                            orthogroup_id=best.group_id,
                            source_rbh_orthogroup_id=best.group_id,
                            target_rbh_orthogroup_id=None,
                            query_id=best.evidence.query,
                            subject_id=best.evidence.subject,
                            hit=best.evidence,
                            protein_map=protein_map,
                            edge_kind=unassigned_relation_kind,
                            render_role="display_edge",
                        )
                    )
            continue

        group_member_ids.setdefault(best.group_id, set()).add(protein_id)
        group_by_protein[protein_id] = best.group_id
        member_roles[protein_id] = assigned_role
        member_confidence[protein_id] = assigned_confidence
        assignment_reasons[protein_id] = assigned_reason
        best_core_support[protein_id] = float(best.support)
        second_core_support[protein_id] = float(second_support)
        support_edge_id = f"{best.evidence.query}->{best.evidence.subject}:{best.relation_kind}"
        supporting_edges[protein_id] = (support_edge_id,)
        ortholog_edges_by_group.setdefault(best.group_id, []).append(
            _make_ortholog_edge(
                orthogroup_id=best.group_id,
                source_rbh_orthogroup_id=best.group_id,
                target_rbh_orthogroup_id=best.group_id,
                query_id=best.evidence.query,
                subject_id=best.evidence.subject,
                hit=best.evidence,
                protein_map=protein_map,
                edge_kind=best.relation_kind,
                render_role="display_edge",
            )
        )

    del core_member_snapshot
    record_local_edges = _select_record_local_paralog_edges(hits, protein_map, thresholds, group_by_protein)
    record_local_components = _record_local_components_from_edges(record_local_edges, protein_map)
    local_support_by_member = _record_local_support_by_member(record_local_edges)
    # Expanded cross-record membership is stable throughout component screening.
    # All components are screened before any record-local groups are appended.
    expanded_core_evidence = (
        _index_core_support_evidence(hits, group_by_protein)
        if record_local_components and group_member_ids else {}
    )
    accepted_record_local_components = tuple(
        component
        for component in record_local_components
        if not _record_local_component_has_competing_core_support(
            component,
            bool(group_member_ids),
            expanded_core_evidence,
            thresholds,
            local_support_by_member,
        )
    )
    del expanded_core_evidence, thresholds
    _append_record_local_orthogroups(
        group_member_ids=group_member_ids,
        group_by_protein=group_by_protein,
        representative_ids_by_group=representative_ids_by_group,
        member_roles=member_roles,
        member_confidence=member_confidence,
        assignment_reasons=assignment_reasons,
        supporting_edges=supporting_edges,
        best_core_support=best_core_support,
        second_core_support=second_core_support,
        ortholog_edges_by_group=ortholog_edges_by_group,
        scope_by_group=scope_by_group,
        source_record_index_by_group=source_record_index_by_group,
        group_order=group_order,
        components=accepted_record_local_components,
        accepted_edges=record_local_edges,
        local_support_by_member=local_support_by_member,
        protein_map=protein_map,
    )

    _append_related_edges(
        hits,
        protein_map,
        group_by_protein=group_by_protein,
        member_roles=member_roles,
        ortholog_edges_by_group=ortholog_edges_by_group,
        related_edges_by_group=related_edges_by_group,
        max_related_edges_per_orthogroup=max_related_edges_per_orthogroup,
    )

    updated_edges_by_group, paths_by_group = _build_ortholog_path_indexes(
        {
            group_id: sorted(
                edges,
                key=lambda edge: (
                    edge.query_record_index,
                    edge.subject_record_index,
                    edge.query_protein_id,
                    edge.subject_protein_id,
                    edge.edge_kind,
                ),
            )
            for group_id, edges in ortholog_edges_by_group.items()
        },
        protein_map,
    )
    rbh_orthogroups: dict[str, tuple[str, ...]] = {}
    for group_id, member_ids in group_member_ids.items():
        anchor_member_ids = [
            member_id
            for member_id in member_ids
            if member_roles.get(member_id) == "anchor"
        ]
        rbh_orthogroups[group_id] = tuple(
            sorted(
                anchor_member_ids or member_ids,
                key=lambda member_id: _protein_sort_key(protein_map[member_id]),
            )
        )
    return _orthogroup_result_from_member_ids(
        group_member_ids,
        representative_ids_by_group,
        protein_map,
        group_order=group_order,
        rbh_orthogroups=rbh_orthogroups,
        ortholog_edges_by_orthogroup_id=updated_edges_by_group,
        path_indexes_by_orthogroup_id=paths_by_group,
        related_edges_by_orthogroup_id={
            group_id: tuple(
                sorted(
                    edges,
                    key=lambda edge: (
                        ORTHOLOG_DISPLAY_EDGE_KIND_RANK.get(str(edge.edge_kind), 99),
                        edge.query_record_index,
                        edge.subject_record_index,
                        edge.query_protein_id,
                        edge.subject_protein_id,
                    ),
                )
            )
            for group_id, edges in related_edges_by_group.items()
        },
        member_roles_by_protein_id=member_roles,
        member_confidence_by_protein_id=member_confidence,
        assignment_reason_by_protein_id=assignment_reasons,
        supporting_edges_by_protein_id=supporting_edges,
        best_core_support_by_protein_id=best_core_support,
        second_best_core_support_by_protein_id=second_core_support,
        scope_by_orthogroup_id=scope_by_group,
        source_record_index_by_orthogroup_id=source_record_index_by_group,
    )


def _select_anchor_core_orthogroup_edges_from_directional_hits(
    directional_hits_by_pair: Mapping[tuple[int, int], DataFrame],
    protein_map: Mapping[str, CdsProtein],
    *,
    record_count: int | None,
    include_singletons: bool,
    max_related_edges_per_orthogroup: int,
    comparison_pairs: Sequence[tuple[int, int]] | None,
) -> OrthogroupEdgeSelectionResult:
    max_related_edges_per_orthogroup = int(max_related_edges_per_orthogroup)
    if max_related_edges_per_orthogroup <= 0:
        raise ValidationError("collinear_max_paralog_links_per_orthogroup must be > 0")
    # The normalized tables are freed once the hits are built.
    hits = _directed_hits_from_tables(
        _normalize_directional_hit_tables(directional_hits_by_pair, protein_map, min_coverage=0.0),
        protein_map,
    )
    anchor_edges = _select_anchor_core_edges(hits, protein_map)
    anchor_edge_tables = _anchor_core_edge_tables(anchor_edges, protein_map)
    if record_count is None:
        record_count = 0
        for query_index, subject_index in directional_hits_by_pair:
            record_count = max(record_count, int(query_index) + 1, int(subject_index) + 1)

    all_edges_by_pair = {
        pair: _comparison_columns_only(table)
        for pair, table in anchor_edge_tables.items()
    }
    display_pairs = (
        tuple(comparison_pairs) if comparison_pairs is not None else
        tuple((index, index + 1) for index in range(max(0, int(record_count) - 1)))
    )
    if any(
        not (0 <= query < int(record_count) and 0 <= subject < int(record_count))
        or query == subject for query, subject in display_pairs
    ):
        raise ValidationError("comparison_pairs contains an invalid record-index pair.")
    if len(set(display_pairs)) != len(display_pairs):
        raise ValidationError("comparison_pairs must not contain duplicates.")
    adjacent_anchor_edges_by_pair = {}
    for query, subject in display_pairs:
        table = anchor_edge_tables.get(cast(tuple[int, int], tuple(sorted((query, subject)))), _empty_comparison_hits())
        if query > subject:
            table = pd.DataFrame.from_records(
                [_comparison_record_for_ids(row, str(row.subject), str(row.query))
                 for row in table.itertuples(index=False)],
                columns=COMPARISON_COLUMNS,
            )
        adjacent_anchor_edges_by_pair[(query, subject)] = _comparison_columns_only(table)

    orthogroups = _build_anchor_core_orthogroups(
        hits,
        anchor_edges,
        protein_map,
        include_singletons=include_singletons,
        max_related_edges_per_orthogroup=max_related_edges_per_orthogroup,
    )
    adjacent_candidate_edges_by_pair = {
        pair: _comparison_columns_only(directional_hits_by_pair.get(pair, _empty_comparison_hits()))
        for pair in display_pairs
    }
    adjacent_display_edges_by_pair = _build_adjacent_display_edges_by_pair(
        adjacent_anchor_edges_by_pair,
        orthogroups,
        max_display_edges_per_orthogroup=max_related_edges_per_orthogroup,
        adjacent_candidate_edges_by_pair=adjacent_candidate_edges_by_pair,
    )
    return OrthogroupEdgeSelectionResult(
        orthogroups=orthogroups,
        all_edges_by_pair=all_edges_by_pair,
        adjacent_anchor_edges_by_pair=adjacent_anchor_edges_by_pair,
        adjacent_display_edges_by_pair=adjacent_display_edges_by_pair,
    )


def _genomic_link_coordinates(protein: CdsProtein) -> tuple[int, int]:
    if protein.strand == -1:
        return protein.end, protein.start + 1
    return protein.start + 1, protein.end


def convert_protein_hits_to_genomic_links(
    hits: DataFrame,
    protein_map: Mapping[str, CdsProtein],
    orthogroups: OrthogroupResult | OrthogroupGraphResult | None = None,
) -> DataFrame:
    """Convert protein hit rows to genomic-coordinate comparison rows."""

    return convert_pair_protein_hits_to_genomic_links(
        hits,
        protein_map,
        protein_map,
        orthogroups=orthogroups,
    )


def _view_feature_svg_id(protein: CdsProtein) -> str:
    return str(
        protein.view_feature_svg_id
        or protein.feature_svg_id
        or ""
    )


def _protein_metadata_value(value: object | None) -> object:
    return "" if value is None else value


def _orthogroup_member_counts(members: Sequence[OrthogroupMember]) -> dict[int, int]:
    counts: dict[int, int] = {}
    for member in members:
        record_index = int(member.record_index)
        counts[record_index] = counts.get(record_index, 0) + 1
    return counts


_OrthogroupEdgeIndex = tuple[
    dict[tuple[str, str], OrthologEdge],
    int,
]


def _index_orthogroup_edges(
    orthogroups: OrthogroupResult | OrthogroupGraphResult,
    orthogroup_id: str,
    requested: tuple[str, str],
    index: _OrthogroupEdgeIndex,
) -> _OrthogroupEdgeIndex:
    """Extend an endpoint index to the requested edge without rereading a prefix."""
    indexed, offset = index
    membership = orthogroups.ortholog_edges_by_orthogroup_id.get(orthogroup_id, ())
    related = orthogroups.related_edges_by_orthogroup_id.get(orthogroup_id, ())
    boundary = len(membership)
    total = boundary + len(related)
    for position in range(offset, total):
        edge = membership[position] if position < boundary else related[position - boundary]
        query_id, subject_id = edge.query_protein_id, edge.subject_protein_id
        key = (query_id, subject_id) if query_id <= subject_id else (subject_id, query_id)
        indexed.setdefault(key, edge)
        if key == requested:
            return indexed, position + 1
    return indexed, total


def _edge_metadata_for_protein_pair(
    orthogroups: OrthogroupResult | OrthogroupGraphResult | None,
    orthogroup_id: str,
    query_id: str,
    subject_id: str,
    edge_indexes: dict[str, _OrthogroupEdgeIndex],
) -> dict[str, object]:
    metadata: dict[str, object] = {
        "rbh_orthogroup_id": "",
        "ortholog_path_id": "",
        "edge_kind": "",
        "render_role": "",
    }
    if orthogroups is None or not orthogroup_id:
        return metadata
    index = edge_indexes.get(orthogroup_id)
    if index is None:
        index = ({}, 0)
    key = (query_id, subject_id) if query_id <= subject_id else (subject_id, query_id)
    if key not in index[0]:
        index = _index_orthogroup_edges(orthogroups, orthogroup_id, key, index)
    edge_indexes[orthogroup_id] = index
    edge = index[0].get(key)
    if edge is not None:
        source_group = str(edge.source_rbh_orthogroup_id or "")
        target_group = str(edge.target_rbh_orthogroup_id or "")
        if source_group and target_group and source_group != target_group:
            rbh_group = f"{source_group};{target_group}"
        else:
            rbh_group = source_group or target_group
        metadata.update(
            {
                "rbh_orthogroup_id": rbh_group,
                "ortholog_path_id": str(edge.path_id or ""),
                "edge_kind": edge.edge_kind,
                "render_role": edge.render_role,
            }
        )
    return metadata


def convert_pair_protein_hits_to_genomic_links(
    hits: DataFrame,
    query_protein_map: Mapping[str, CdsProtein],
    subject_protein_map: Mapping[str, CdsProtein],
    orthogroups: OrthogroupResult | OrthogroupGraphResult | None = None,
) -> DataFrame:
    """Convert pairwise protein hit rows using separate query and subject maps."""

    if hits.empty:
        return pd.DataFrame(columns=LOSATP_COMPARISON_COLUMNS)

    # Build each touched group once, preserving lazy handling of unused groups.
    edge_indexes: dict[str, _OrthogroupEdgeIndex] = {}
    member_counts: dict[str, dict[int, int]] = {}
    rows: list[dict[str, object]] = []
    missing_ids: set[str] = set()
    for row in hits.itertuples(index=False):
        query_id = str(row.query)
        subject_id = str(row.subject)
        query_protein = query_protein_map.get(query_id)
        subject_protein = subject_protein_map.get(subject_id)
        if query_protein is None:
            missing_ids.add(query_id)
        if subject_protein is None:
            missing_ids.add(subject_id)
        if query_protein is None or subject_protein is None:
            continue

        qstart, qend = _genomic_link_coordinates(query_protein)
        sstart, send = _genomic_link_coordinates(subject_protein)
        query_member = (
            orthogroups.member_by_protein_id.get(query_id)
            if orthogroups is not None
            else None
        )
        subject_member = (
            orthogroups.member_by_protein_id.get(subject_id)
            if orthogroups is not None
            else None
        )
        orthogroup_id = ""
        if query_member is not None and subject_member is not None:
            if query_member.orthogroup_id == subject_member.orthogroup_id:
                orthogroup_id = query_member.orthogroup_id
        edge_metadata = _edge_metadata_for_protein_pair(
            orthogroups,
            orthogroup_id,
            query_id,
            subject_id,
            edge_indexes,
        )
        if orthogroup_id and orthogroup_id not in member_counts:
            # A non-empty orthogroup_id comes from orthogroup members, so orthogroups is set.
            member_counts[orthogroup_id] = _orthogroup_member_counts(
                cast("OrthogroupResult | OrthogroupGraphResult", orthogroups).orthogroups.get(orthogroup_id, [])
            )
        counts = member_counts.get(orthogroup_id, {})
        rows.append(
            {
                "query": query_protein.record_id,
                "subject": subject_protein.record_id,
                "identity": row.identity,
                "alignment_length": row.alignment_length,
                "mismatches": row.mismatches,
                "gap_opens": row.gap_opens,
                "qstart": qstart,
                "qend": qend,
                "sstart": sstart,
                "send": send,
                "evalue": row.evalue,
                "bitscore": row.bitscore,
                "query_protein_id": query_protein.protein_id,
                "subject_protein_id": subject_protein.protein_id,
                "query_source_protein_id": _protein_metadata_value(query_protein.source_protein_id),
                "subject_source_protein_id": _protein_metadata_value(subject_protein.source_protein_id),
                "query_record_index": query_protein.record_index,
                "subject_record_index": subject_protein.record_index,
                "query_feature_index": query_protein.feature_index,
                "subject_feature_index": subject_protein.feature_index,
                "query_feature_svg_id": str(
                    query_protein.feature_svg_id or ""
                ),
                "subject_feature_svg_id": str(
                    subject_protein.feature_svg_id or ""
                ),
                "query_view_feature_svg_id": _view_feature_svg_id(
                    query_protein
                ),
                "subject_view_feature_svg_id": _view_feature_svg_id(
                    subject_protein
                ),
                "orthogroup_id": orthogroup_id,
                "rbh_orthogroup_id": edge_metadata["rbh_orthogroup_id"],
                "ortholog_path_id": edge_metadata["ortholog_path_id"],
                "edge_kind": edge_metadata["edge_kind"],
                "render_role": edge_metadata["render_role"],
                "query_orthogroup_representative": (
                    bool(query_member.representative) if query_member is not None else False
                ),
                "subject_orthogroup_representative": (
                    bool(subject_member.representative) if subject_member is not None else False
                ),
                "query_orthogroup_member_count": counts.get(int(query_protein.record_index), 0) if counts else 0,
                "subject_orthogroup_member_count": counts.get(int(subject_protein.record_index), 0) if counts else 0,
                "query_orthogroup_role": (
                    str(query_member.role) if query_member is not None else ""
                ),
                "subject_orthogroup_role": (
                    str(subject_member.role) if subject_member is not None else ""
                ),
                "query_orthogroup_confidence": (
                    str(query_member.confidence) if query_member is not None else ""
                ),
                "subject_orthogroup_confidence": (
                    str(subject_member.confidence) if subject_member is not None else ""
                ),
                "query_orthogroup_assignment_reason": (
                    str(query_member.assignment_reason) if query_member is not None else ""
                ),
                "subject_orthogroup_assignment_reason": (
                    str(subject_member.assignment_reason) if subject_member is not None else ""
                ),
            }
        )

    if missing_ids:
        raise ParseError(
            "LOSATP blastp output contains unknown protein IDs: "
            + ", ".join(sorted(missing_ids))
        )

    return pd.DataFrame.from_records(rows, columns=LOSATP_COMPARISON_COLUMNS)


def run_losatp_blastp(
    query_fasta: str,
    subject_fasta: str,
    *,
    losatp_bin: str = "losat",
    ncbi_blastp_bin: str | None = None,
    max_hits: int | None = None,
    max_hsps_per_subject: int | None = 1,
    threads: int | None = None,
    raw_output_callback: Callable[[str], None] | None = None,
    runtime_callback: LosatRuntimeCallback | None = None,
) -> DataFrame:
    """Run LOSATP through the shared LOSAT runtime and parse outfmt 6 output."""

    if max_hits is not None:
        _validate_max_hits(max_hits, option_name="protein_blastp_candidate_limit")
    if max_hsps_per_subject is not None:
        _validate_max_hits(max_hsps_per_subject, option_name="max_hsps_per_subject")
    _validate_losatp_threads(threads)
    raw_text = run_losat_search(
        "losatp",
        query_fasta,
        subject_fasta,
        options=LosatSearchArgs(
            max_hsps=max_hsps_per_subject,
            max_target_seqs=max_hits,
        ),
        losat_bin=losatp_bin,
        ncbi_blast_bin=ncbi_blastp_bin,
        threads=threads,
        runtime_callback=runtime_callback,
    )
    if raw_output_callback is not None:
        raw_output_callback(raw_text)
    return parse_losatp_outfmt6(raw_text)


def _validate_extraction_has_proteins(
    records: Sequence[SeqRecord],
    extraction: ProteinExtractionResult,
    *,
    option_name: str,
) -> None:
    empty_record_ids = [
        str(records[index].id)
        for index, proteins in enumerate(extraction.proteins_by_record)
        if not proteins
    ]
    if empty_record_ids:
        raise ValidationError(
            f"{option_name} requires at least one CDS protein in each record; "
            "no CDS proteins were found in: "
            + ", ".join(empty_record_ids)
        )


def _split_protein_hits(
    hits: DataFrame,
    batch: LosatBatch,
    pairs: Sequence[RecordPair],
) -> dict[RecordPair, DataFrame]:
    """Split the hits of one source job into record-pair hits, in runtime order."""

    query_records = [batch.query.ids.get(str(value), (None, ""))[0] for value in hits["query"]]
    subject_records = [
        batch.subject.ids.get(str(value), (None, ""))[0] for value in hits["subject"]
    ]
    if None in query_records or None in subject_records:
        raise ValidationError(
            "LOSAT output references IDs outside the current protein bindings.",
            diagnostic={"code": "LOSAT_RUNTIME", "reason": "OUTPUT"},
        )
    return {
        pair: hits.loc[
            [
                query == pair[0] and subject == pair[1]
                for query, subject in zip(query_records, subject_records)
            ]
        ].reset_index(drop=True)
        for pair in pairs
    }


def _consecutive_display_filenames(
    record_count: int,
    cache_filenames: Sequence[str] | None,
) -> dict[RecordPair, str]:
    """Displayed pairs are the consecutive records; ``cache_filenames[i]`` names (i, i+1)."""

    return {
        (index, index + 1): (
            str(cache_filenames[index])
            if cache_filenames is not None and index < len(cache_filenames)
            else ""
        )
        for index in range(max(0, record_count - 1))
    }


def _search_losatp_record_pairs(
    records: Sequence[SeqRecord],
    extraction: ProteinExtractionResult,
    specs: Sequence[RecordPair],
    *,
    losatp_bin: str,
    ncbi_blastp_bin: str | None,
    losatp_threads: int | None,
    candidate_limit: int | None,
    max_hsps_per_subject: int | None,
    runner: LosatpRunner | None,
    losatp_cache: LosatpCacheManager | None,
    display_filenames: Mapping[RecordPair, str] | None = None,
) -> dict[RecordPair, DataFrame]:
    """Answer directed record-pair searches with the LOSAT job plan (design D7).

    A source file is one database: each source job of
    :func:`gbdraw.comparisons.losat_jobs.plan_losat_jobs` runs once, and its
    result is split into record-pair hits and raw cache entries under the Web
    raw keys. ``display_filenames`` names the displayed pairs.
    """

    args = tuple(
        losat_cache_args(
            "losatp",
            LosatSearchArgs(max_hsps=max_hsps_per_subject, max_target_seqs=candidate_limit),
        )
    )
    fasta_by_record: dict[int, str] = {}

    def record_fasta(index: int) -> str:
        if index not in fasta_by_record:
            fasta_by_record[index] = proteins_to_fasta(extraction.proteins_by_record[index])
        return fasta_by_record[index]

    batches = prepare_losat_batches(
        plan_losat_jobs(
            source_ids=losat_source_ids(records),
            specs=specs,
            build_args=lambda _query, _subject: args,
        ),
        uids=losat_record_uids(records),
        record_fasta=record_fasta,
        protein=True,
    )
    displayed = dict(display_filenames or {})
    cache = losatp_cache if runner is None else None
    hits: dict[RecordPair, DataFrame] = {}
    for batch in batches:
        texts: dict[RecordPair, str] = {}
        if cache is not None:
            for pair in batch.job.specs:
                text = cache.cached_text(
                    record_fasta(pair[0]),
                    record_fasta(pair[1]),
                    args=args,
                    search_context=batch.search_context,
                    display=pair in displayed,
                    filename=displayed.get(pair, ""),
                )
                if text is not None:
                    texts[pair] = text
        for pair, text in texts.items():
            hits[pair] = parse_losatp_outfmt6(text)
        missing = [pair for pair in batch.job.specs if pair not in texts]
        if not missing:
            continue
        raw_text: list[str] = []
        runtime_records: list[LosatRuntimeRecord] = []
        if runner is not None:
            batch_hits = runner(batch.query.fasta, batch.subject.fasta)
        else:
            batch_hits = run_losatp_blastp(
                batch.query.fasta,
                batch.subject.fasta,
                losatp_bin=losatp_bin,
                ncbi_blastp_bin=ncbi_blastp_bin,
                max_hits=candidate_limit,
                max_hsps_per_subject=max_hsps_per_subject,
                threads=losatp_threads,
                raw_output_callback=raw_text.append,
                runtime_callback=runtime_records.append,
            )
        hits.update(_split_protein_hits(batch_hits, batch, missing))
        if cache is None:
            continue
        split = split_losat_batch_result("".join(raw_text), batch, missing)
        for pair in missing:
            cache.store_search(
                record_fasta(pair[0]),
                record_fasta(pair[1]),
                args=args,
                search_context=batch.search_context,
                text=split[pair],
                runtime=runtime_records[0] if runtime_records else None,
                display=pair in displayed,
                filename=displayed.get(pair, ""),
            )
    return hits


def _empty_comparison_hits() -> DataFrame:
    return pd.DataFrame(columns=COMPARISON_COLUMNS)


def _orthogroup_scope(
    orthogroups: OrthogroupResult | OrthogroupGraphResult,
    orthogroup_id: str,
) -> OrthogroupScope:
    return orthogroups.scope_by_orthogroup_id.get(str(orthogroup_id), "cross_record")


def _ortholog_edge_display_rank(
    edge: OrthologEdge,
    orthogroups: OrthogroupResult | OrthogroupGraphResult,
) -> tuple[object, ...]:
    query_member = orthogroups.member_by_protein_id.get(edge.query_protein_id)
    subject_member = orthogroups.member_by_protein_id.get(edge.subject_protein_id)
    both_representative = (
        query_member is not None
        and subject_member is not None
        and bool(query_member.representative)
        and bool(subject_member.representative)
    )
    return (
        ORTHOLOG_DISPLAY_EDGE_KIND_RANK.get(str(edge.edge_kind), 99),
        0 if both_representative else 1,
        float(edge.evalue),
        -float(edge.bitscore),
        -float(edge.identity),
        -int(edge.alignment_length),
        int(edge.query_record_index),
        int(edge.subject_record_index),
        str(edge.query_protein_id),
        str(edge.subject_protein_id),
        str(edge.edge_kind),
    )


def _display_pair_orthogroup_id(
    query_id: str,
    subject_id: str,
    orthogroups: OrthogroupResult | OrthogroupGraphResult,
) -> str | None:
    query_member = orthogroups.member_by_protein_id.get(str(query_id))
    subject_member = orthogroups.member_by_protein_id.get(str(subject_id))
    if query_member is None or subject_member is None:
        return None
    if query_member.orthogroup_id != subject_member.orthogroup_id:
        return None
    return str(query_member.orthogroup_id)


def _comparison_row_from_ortholog_edge_for_pair(
    edge: OrthologEdge,
    pair: tuple[int, int],
) -> dict[str, object] | None:
    query_index, subject_index = int(pair[0]), int(pair[1])
    if (
        int(edge.query_record_index) == query_index
        and int(edge.subject_record_index) == subject_index
    ):
        query_id = edge.query_protein_id
        subject_id = edge.subject_protein_id
    elif (
        int(edge.query_record_index) == subject_index
        and int(edge.subject_record_index) == query_index
    ):
        query_id = edge.subject_protein_id
        subject_id = edge.query_protein_id
    else:
        return None

    alignment_length = max(0, int(edge.alignment_length))
    return {
        "query": str(query_id),
        "subject": str(subject_id),
        "identity": float(edge.identity),
        "alignment_length": alignment_length,
        "mismatches": 0,
        "gap_opens": 0,
        "qstart": 1,
        "qend": alignment_length,
        "sstart": 1,
        "send": alignment_length,
        "evalue": float(edge.evalue),
        "bitscore": float(edge.bitscore),
    }


def _comparison_row_key(row: Mapping[str, object]) -> tuple[str, str]:
    return str(row["query"]), str(row["subject"])


def _comparison_row_display_rank(
    row: Mapping[str, object],
    orthogroups: OrthogroupResult | OrthogroupGraphResult,
    *,
    edge_kind_rank: int,
) -> tuple[object, ...]:
    query_id, subject_id = _comparison_row_key(row)
    query_member = orthogroups.member_by_protein_id.get(query_id)
    subject_member = orthogroups.member_by_protein_id.get(subject_id)
    both_representative = (
        query_member is not None
        and subject_member is not None
        and bool(query_member.representative)
        and bool(subject_member.representative)
    )
    return (
        int(edge_kind_rank),
        0 if both_representative else 1,
        float(cast(SupportsFloat, row.get("evalue", float("inf")))),
        -float(cast(SupportsFloat, row.get("bitscore", 0.0))),
        -float(cast(SupportsFloat, row.get("identity", 0.0))),
        -int(float(cast(SupportsFloat, row.get("alignment_length", 0.0)))),
        int(query_member.record_index) if query_member is not None else 0,
        int(subject_member.record_index) if subject_member is not None else 0,
        query_id,
        subject_id,
        "direct_orthogroup",
    )


def _comparison_row_from_hit_row(row: object) -> dict[str, object]:
    return {
        column: getattr(row, column)
        for column in COMPARISON_COLUMNS
    }


def _has_uncovered_display_alternative(
    row_key: tuple[str, str],
    *,
    current_kind_rank: int,
    candidate_subjects_by_query: Mapping[str, set[str]],
    candidate_queries_by_subject: Mapping[str, set[str]],
    candidate_kind_rank_by_pair: Mapping[tuple[str, str], int],
    covered_queries: set[str],
    covered_subjects: set[str],
    existing_pairs: set[tuple[str, str]],
) -> bool:
    query_id, subject_id = row_key
    if query_id not in covered_queries and subject_id in covered_subjects:
        for alternative_subject_id in candidate_subjects_by_query.get(query_id, set()):
            alternative_key = (query_id, alternative_subject_id)
            if alternative_subject_id == subject_id:
                continue
            if alternative_subject_id in covered_subjects:
                continue
            if alternative_key in existing_pairs or (alternative_key[1], alternative_key[0]) in existing_pairs:
                continue
            if candidate_kind_rank_by_pair.get(alternative_key, 99) > current_kind_rank:
                continue
            return True
    if query_id in covered_queries and subject_id not in covered_subjects:
        for alternative_query_id in candidate_queries_by_subject.get(subject_id, set()):
            alternative_key = (alternative_query_id, subject_id)
            if alternative_query_id == query_id:
                continue
            if alternative_query_id in covered_queries:
                continue
            if alternative_key in existing_pairs or (alternative_key[1], alternative_key[0]) in existing_pairs:
                continue
            if candidate_kind_rank_by_pair.get(alternative_key, 99) > current_kind_rank:
                continue
            return True
    return False


def _add_display_row(
    row: dict[str, object],
    *,
    existing_pairs: set[tuple[str, str]],
    covered_queries: set[str],
    covered_subjects: set[str],
    added_rows: list[dict[str, object]],
) -> None:
    row_key = _comparison_row_key(row)
    existing_pairs.add(row_key)
    existing_pairs.add((row_key[1], row_key[0]))
    covered_queries.add(row_key[0])
    covered_subjects.add(row_key[1])
    added_rows.append(row)


def _orthogroup_display_table(
    table: DataFrame,
    orthogroups: OrthogroupResult | OrthogroupGraphResult,
) -> DataFrame:
    if table is None or table.empty:
        return _empty_comparison_hits()

    rows: list[dict[str, object]] = []
    for row in table.itertuples(index=False):
        query_id = str(row.query)
        subject_id = str(row.subject)
        orthogroup_id = _display_pair_orthogroup_id(
            query_id,
            subject_id,
            orthogroups,
        )
        if orthogroup_id is None:
            continue
        if _orthogroup_scope(orthogroups, orthogroup_id) == "record_local":
            continue
        rows.append(_comparison_row_from_hit_row(row))
    if not rows:
        return _empty_comparison_hits()
    return pd.DataFrame.from_records(rows, columns=COMPARISON_COLUMNS)


def _build_adjacent_display_edges_by_pair(
    adjacent_anchor_edges_by_pair: Mapping[tuple[int, int], DataFrame],
    orthogroups: OrthogroupResult | OrthogroupGraphResult,
    *,
    max_display_edges_per_orthogroup: int,
    adjacent_candidate_edges_by_pair: Mapping[tuple[int, int], DataFrame] | None = None,
) -> dict[tuple[int, int], DataFrame]:
    """Overlay bounded non-anchor orthogroup edges for adjacent display."""

    cap = int(max_display_edges_per_orthogroup)
    if cap <= 0:
        raise ValidationError("collinear_max_paralog_links_per_orthogroup must be > 0")

    display_edges_by_pair: dict[tuple[int, int], DataFrame] = {
        (int(pair[0]), int(pair[1])): _orthogroup_display_table(table, orthogroups)
        for pair, table in adjacent_anchor_edges_by_pair.items()
    }

    all_secondary_edges_by_group = {
        orthogroup_id: tuple(
            edge
            for edge in (
                *orthogroups.ortholog_edges_by_orthogroup_id.get(orthogroup_id, ()),
                *orthogroups.related_edges_by_orthogroup_id.get(orthogroup_id, ()),
            )
            if edge.render_role == "display_edge"
        )
        for orthogroup_id in orthogroups.orthogroups
        if _orthogroup_scope(orthogroups, orthogroup_id) != "record_local"
    }

    for pair, anchor_table in list(display_edges_by_pair.items()):
        existing_pairs = {
            (str(row.query), str(row.subject))
            for row in anchor_table.itertuples(index=False)
        }
        existing_pairs.update((subject_id, query_id) for query_id, subject_id in list(existing_pairs))
        # Secondary edges should reveal extra paralog endpoints, not redraw
        # cross-links between endpoints already covered by stronger edges.
        covered_query_by_group: dict[str, set[str]] = {}
        covered_subject_by_group: dict[str, set[str]] = {}
        for row in anchor_table.itertuples(index=False):
            query_id = str(row.query)
            subject_id = str(row.subject)
            orthogroup_id = _display_pair_orthogroup_id(
                query_id,
                subject_id,
                orthogroups,
            )
            if orthogroup_id is None:
                continue
            if _orthogroup_scope(orthogroups, orthogroup_id) == "record_local":
                continue
            covered_query_by_group.setdefault(orthogroup_id, set()).add(query_id)
            covered_subject_by_group.setdefault(orthogroup_id, set()).add(subject_id)
        added_rows: list[dict[str, object]] = []
        direct_candidates_by_group: dict[str, list[tuple[tuple[object, ...], dict[str, object]]]] = {}
        direct_candidate_table = (
            adjacent_candidate_edges_by_pair.get(pair)
            if adjacent_candidate_edges_by_pair is not None
            else None
        )
        if direct_candidate_table is not None and not direct_candidate_table.empty:
            direct_rows = _sort_hits_for_query_best(
                _coerce_outfmt6_numeric_columns(
                    direct_candidate_table.loc[:, list(COMPARISON_COLUMNS)]
                )
            ).drop_duplicates(["query", "subject"], keep="first")
            for row in direct_rows.itertuples(index=False):
                query_id = str(row.query)
                subject_id = str(row.subject)
                orthogroup_id = _display_pair_orthogroup_id(
                    query_id,
                    subject_id,
                    orthogroups,
                )
                if orthogroup_id is None:
                    continue
                if _orthogroup_scope(orthogroups, orthogroup_id) == "record_local":
                    continue
                row_record = _comparison_row_from_hit_row(row)
                direct_candidates_by_group.setdefault(orthogroup_id, []).append(
                    (
                        _comparison_row_display_rank(
                            row_record,
                            orthogroups,
                            edge_kind_rank=ORTHOLOG_DISPLAY_EDGE_KIND_RANK["coortholog"],
                        ),
                        row_record,
                    )
                )

        for orthogroup_id in sorted(set(all_secondary_edges_by_group).union(direct_candidates_by_group)):
            if _orthogroup_scope(orthogroups, orthogroup_id) == "record_local":
                continue
            candidates: list[tuple[tuple[object, ...], dict[str, object]]] = []
            candidates.extend(direct_candidates_by_group.get(orthogroup_id, []))
            for edge in all_secondary_edges_by_group[orthogroup_id]:
                row = _comparison_row_from_ortholog_edge_for_pair(edge, pair)
                if row is None:
                    continue
                row_key = (str(row["query"]), str(row["subject"]))
                row_orthogroup_id = _display_pair_orthogroup_id(
                    row_key[0],
                    row_key[1],
                    orthogroups,
                )
                if row_orthogroup_id != orthogroup_id:
                    continue
                if row_key in existing_pairs or (row_key[1], row_key[0]) in existing_pairs:
                    continue
                candidates.append((_ortholog_edge_display_rank(edge, orthogroups), row))
            added_for_group = 0
            covered_queries = covered_query_by_group.setdefault(orthogroup_id, set())
            covered_subjects = covered_subject_by_group.setdefault(orthogroup_id, set())
            sorted_candidates = sorted(candidates, key=lambda item: item[0])
            candidate_subjects_by_query: dict[str, set[str]] = {}
            candidate_queries_by_subject: dict[str, set[str]] = {}
            candidate_kind_rank_by_pair: dict[tuple[str, str], int] = {}
            for rank, row in sorted_candidates:
                query_id, subject_id = _comparison_row_key(row)
                candidate_subjects_by_query.setdefault(query_id, set()).add(subject_id)
                candidate_queries_by_subject.setdefault(subject_id, set()).add(query_id)
                candidate_kind_rank_by_pair[(query_id, subject_id)] = int(cast(int, rank[0]))

            # Display selection is stricter than membership expansion: prefer
            # same-quality links that introduce new endpoints on both records.
            for require_both_uncovered in (True, False):
                for rank, row in sorted_candidates:
                    if added_for_group >= cap:
                        break
                    row_key = _comparison_row_key(row)
                    if row_key in existing_pairs or (row_key[1], row_key[0]) in existing_pairs:
                        continue
                    query_covered = row_key[0] in covered_queries
                    subject_covered = row_key[1] in covered_subjects
                    if query_covered and subject_covered:
                        continue
                    if require_both_uncovered:
                        if query_covered or subject_covered:
                            continue
                    else:
                        if not (query_covered or subject_covered):
                            continue
                        if _has_uncovered_display_alternative(
                            row_key,
                            current_kind_rank=int(cast(int, rank[0])),
                            candidate_subjects_by_query=candidate_subjects_by_query,
                            candidate_queries_by_subject=candidate_queries_by_subject,
                            candidate_kind_rank_by_pair=candidate_kind_rank_by_pair,
                            covered_queries=covered_queries,
                            covered_subjects=covered_subjects,
                            existing_pairs=existing_pairs,
                        ):
                            continue
                    _add_display_row(
                        row,
                        existing_pairs=existing_pairs,
                        covered_queries=covered_queries,
                        covered_subjects=covered_subjects,
                        added_rows=added_rows,
                    )
                    added_for_group += 1
                if added_for_group >= cap:
                    break
        if added_rows:
            added_frame = pd.DataFrame.from_records(added_rows, columns=COMPARISON_COLUMNS)
            if anchor_table.empty:
                display_edges_by_pair[pair] = added_frame.reset_index(drop=True)
            else:
                display_edges_by_pair[pair] = pd.concat(
                    [
                        anchor_table.loc[:, list(COMPARISON_COLUMNS)],
                        added_frame,
                    ],
                    ignore_index=True,
                )
    return display_edges_by_pair


def select_rbh_orthogroup_edges_from_directional_hits(
    directional_hits_by_pair: Mapping[tuple[int, int], DataFrame],
    protein_map: Mapping[str, CdsProtein],
    *,
    record_count: int | None = None,
    include_singletons: bool = False,
    orthogroup_membership_mode: OrthogroupMembershipMode | str = ORTHOGROUP_INFERENCE_VERSION,
    orthogroup_member_max_hits: int | None = None,
    max_related_edges_per_orthogroup: int = 2,
    path_representation: Literal["graph", "exhaustive"] = "graph",
    comparison_pairs: Sequence[tuple[int, int]] | None = None,
) -> OrthogroupEdgeSelectionResult:
    """Infer groups from all evidence and project links onto requested display pairs.

    Omitted comparison_pairs preserves consecutive-record display. An empty
    sequence produces membership without display links.
    """

    if path_representation not in {"graph", "exhaustive"}:
        raise ValidationError("path_representation must be graph or exhaustive")
    normalize_orthogroup_membership_mode(str(orthogroup_membership_mode))
    if orthogroup_member_max_hits is not None:
        _validate_max_hits(
            orthogroup_member_max_hits,
            option_name="orthogroup_member_max_hits",
        )
    member_hits = {
        pair: _select_member_candidate_hits_per_query(
            hits,
            max_hits=orthogroup_member_max_hits,
        )
        for pair, hits in directional_hits_by_pair.items()
    }
    result = _select_anchor_core_orthogroup_edges_from_directional_hits(
        member_hits,
        protein_map,
        record_count=record_count,
        include_singletons=include_singletons,
        max_related_edges_per_orthogroup=max_related_edges_per_orthogroup,
        comparison_pairs=comparison_pairs,
    )

    if path_representation == "exhaustive":
        return replace(result, orthogroups=materialize_ortholog_paths(result.orthogroups))
    return result


def build_pairwise_protein_blastp_comparisons(
    records: Sequence[SeqRecord],
    *,
    losatp_bin: str = "losat",
    ncbi_blastp_bin: str | None = None,
    losatp_threads: int | None = None,
    max_hits: int = 5,
    candidate_limit: int | None = None,
    evalue: float = 1e-5,
    bitscore: float = 50.0,
    identity: float = 70.0,
    alignment_length: int = 0,
    runner: LosatpRunner | None = None,
    losatp_cache: LosatpCacheManager | None = None,
    protein_extraction: ProteinExtractionResult | None = None,
    feature_visibility_rules: list[dict[str, object]] | None = None,
    cache_filenames: Sequence[str] | None = None,
    pairs: Sequence[RecordPair] | None = None,
) -> ProteinBlastpResult:
    """Generate LOSATP blastp display comparisons for record pairs.

    ``pairs`` are directed (query, subject) record indexes; omitted, they are
    the consecutive records. ``comparisons`` follows ``pairs`` and
    ``cache_filenames`` names the raw TSV of each pair.
    """

    if len(records) < 2:
        raise ValidationError("protein_blastp_mode='pairwise' requires at least two records")
    _validate_max_hits(max_hits)
    _validate_losatp_threads(losatp_threads)
    _validate_candidate_limit(candidate_limit)
    if int(alignment_length) < 0:
        raise ValidationError("alignment_length must be >= 0")

    extraction = protein_extraction or extract_cds_proteins(
        records,
        feature_visibility_rules=feature_visibility_rules,
    )
    specs = losatp_job_specs("pairwise", record_count=len(records), pairs=pairs)
    searched = sorted({index for pair in specs for index in pair})
    _validate_extraction_has_proteins(
        [records[index] for index in searched],
        replace(
            extraction,
            proteins_by_record=[extraction.proteins_by_record[index] for index in searched],
        ),
        option_name="protein_blastp_mode='pairwise'",
    )
    hits_by_pair = _search_losatp_record_pairs(
        records,
        extraction,
        specs,
        losatp_bin=losatp_bin,
        ncbi_blastp_bin=ncbi_blastp_bin,
        losatp_threads=losatp_threads,
        candidate_limit=candidate_limit,
        max_hsps_per_subject=1,
        runner=runner,
        losatp_cache=losatp_cache,
        display_filenames={
            pair: str(cache_filenames[index])
            if cache_filenames is not None and index < len(cache_filenames)
            else ""
            for index, pair in enumerate(specs)
        },
    )
    comparisons: list[DataFrame] = []
    for pair in specs:
        filtered_hits = filter_protein_hits_by_thresholds(
            hits_by_pair[pair],
            evalue=evalue,
            bitscore=bitscore,
            identity=identity,
            alignment_length=alignment_length,
        )
        display_hits = select_top_hits_per_query(filtered_hits, max_hits=max_hits)
        comparisons.append(
            convert_protein_hits_to_genomic_links(
                display_hits,
                extraction.protein_map,
                orthogroups=None,
            )
        )
    return ProteinBlastpResult(comparisons=comparisons, orthogroups=None)


def build_rbh_orthogroup_protein_blastp_comparisons(
    records: Sequence[SeqRecord],
    *,
    path_representation: Literal["graph", "exhaustive"] = "graph",
    losatp_bin: str = "losat",
    ncbi_blastp_bin: str | None = None,
    losatp_threads: int | None = None,
    candidate_limit: int | None = None,
    orthogroup_membership_mode: OrthogroupMembershipMode | str = ORTHOGROUP_INFERENCE_VERSION,
    orthogroup_member_max_hits: int | None = None,
    max_related_edges_per_orthogroup: int = 2,
    evalue: float = 1e-5,
    bitscore: float = 50.0,
    identity: float = 70.0,
    alignment_length: int = 0,
    runner: LosatpRunner | None = None,
    losatp_cache: LosatpCacheManager | None = None,
    protein_extraction: ProteinExtractionResult | None = None,
    feature_visibility_rules: list[dict[str, object]] | None = None,
    cache_filenames: Sequence[str] | None = None,
) -> ProteinBlastpResult:
    """Infer all-vs-all RBH-seeded orthogroups and return adjacent display links."""

    if path_representation not in {"graph", "exhaustive"}:
        raise ValidationError("path_representation must be graph or exhaustive")
    if len(records) < 2:
        raise ValidationError("protein_blastp_mode='orthogroup' requires at least two records")
    _validate_losatp_threads(losatp_threads)
    _validate_candidate_limit(candidate_limit)
    normalized_membership_mode = normalize_orthogroup_membership_mode(str(orthogroup_membership_mode))
    if orthogroup_member_max_hits is not None:
        _validate_max_hits(orthogroup_member_max_hits, option_name="orthogroup_member_max_hits")
    if int(max_related_edges_per_orthogroup) <= 0:
        raise ValidationError("collinear_max_paralog_links_per_orthogroup must be > 0")
    if int(alignment_length) < 0:
        raise ValidationError("alignment_length must be >= 0")

    extraction = protein_extraction or extract_cds_proteins(
        records,
        feature_visibility_rules=feature_visibility_rules,
    )
    _validate_extraction_has_proteins(
        records,
        extraction,
        option_name="protein_blastp_mode='orthogroup'",
    )

    specs = losatp_job_specs("orthogroup", record_count=len(records))
    hits_by_pair = _search_losatp_record_pairs(
        records,
        extraction,
        specs,
        losatp_bin=losatp_bin,
        ncbi_blastp_bin=ncbi_blastp_bin,
        losatp_threads=losatp_threads,
        candidate_limit=candidate_limit,
        max_hsps_per_subject=None,
        runner=runner,
        losatp_cache=losatp_cache,
        display_filenames=_consecutive_display_filenames(len(records), cache_filenames),
    )
    directional_hits_by_pair: dict[tuple[int, int], DataFrame] = {
        pair: filter_protein_hits_by_thresholds(
            hits_by_pair[pair],
            evalue=evalue,
            bitscore=bitscore,
            identity=identity,
            alignment_length=alignment_length,
        )
        for pair in specs
    }

    edge_selection = select_rbh_orthogroup_edges_from_directional_hits(
        directional_hits_by_pair,
        extraction.protein_map,
        path_representation=path_representation,
        record_count=len(records),
        orthogroup_membership_mode=normalized_membership_mode,
        orthogroup_member_max_hits=orthogroup_member_max_hits,
        max_related_edges_per_orthogroup=max_related_edges_per_orthogroup,
    )
    comparisons = [
        convert_protein_hits_to_genomic_links(
            edge_selection.adjacent_display_edges_by_pair.get(
                (record_index, record_index + 1),
                _empty_comparison_hits(),
            ),
            extraction.protein_map,
            orthogroups=edge_selection.orthogroups,
        )
        for record_index in range(len(records) - 1)
    ]
    return ProteinBlastpResult(comparisons=comparisons, orthogroups=edge_selection.orthogroups)


__all__ = [
    "CdsProtein",
    "LEGACY_PROTEIN_LOSAT_CACHE_SCHEMA",
    "LEGACY_PROTEIN_RAW_CANDIDATE_SCHEMA",
    "LOSAT_CACHE_SCHEMA",
    "LOSATP_COMPARISON_COLUMNS",
    "LOSATP_METADATA_COLUMNS",
    "LegacyProteinRawCachePromotion",
    "LegacyProteinRawCachePromotionScan",
    "LegacyProteinRawCacheRejection",
    "LosatpCacheManager",
    "LosatpCacheState",
    "ORTHOGROUP_INFERENCE_VERSION",
    "ORTHOGROUP_MEMBERSHIP_MODES",
    "PROTEIN_BLASTP_MODES",
    "OrthogroupEdgeSelectionResult",
    "OrthogroupMember",
    "OrthogroupMemberConfidence",
    "OrthogroupMemberRole",
    "OrthogroupMembershipMode",
    "OrthogroupNameCandidate",
    "OrthogroupResult",
    "OrthogroupGraphResult",
    "OrthologPathCollection",
    "materialize_ortholog_paths",
    "OrthogroupScope",
    "OrthologEdge",
    "OrthologEdgeKind",
    "OrthologPath",
    "OrthologRenderRole",
    "ProteinBlastpMode",
    "ProteinBlastpResult",
    "ProteinExtractionResult",
    "ProteinIdentityManifest",
    "ProteinLosatPairIdentity",
    "PROTEIN_IDENTITY_MANIFEST_SCHEMA",
    "PROTEIN_LOSAT_CACHE_SCHEMA",
    "build_legacy_protein_reference_map",
    "build_protein_export_id_map",
    "build_protein_losat_cache_key",
    "build_protein_losat_pair_identity",
    "build_protein_runtime_handle",
    "build_web_losat_cache_key",
    "build_pairwise_protein_blastp_comparisons",
    "build_rbh_orthogroup_protein_blastp_comparisons",
    "convert_pair_protein_hits_to_genomic_links",
    "convert_protein_hits_to_genomic_links",
    "canonical_feature_analysis_id",
    "canonical_protein_identity_json",
    "extract_cds_proteins",
    "extract_protein_identity_manifest",
    "extract_web_stable_cds_proteins",
    "filter_protein_hits_by_thresholds",
    "hydrate_protein_losat_tsv",
    "normalize_orthogroup_membership_mode",
    "normalize_protein_blastp_mode",
    "parse_losatp_outfmt6",
    "parse_losat_fasta_ids",
    "percent_encode_losat_transport_field",
    "proteins_to_fasta",
    "promote_legacy_protein_raw_cache_entries",
    "raw_protein_tsv_matches_bindings",
    "run_losatp_blastp",
    "select_protein_display_alias",
    "select_best_hits_per_query",
    "select_rbh_orthogroup_edges_from_directional_hits",
    "select_reciprocal_best_hit_edges",
    "select_reciprocal_best_hits",
    "select_top_hits_per_query",
    "is_legacy_protein_losat_cache_entry",
    "is_protein_losat_cache_entry",
    "make_legacy_protein_raw_candidate",
    "validate_legacy_protein_raw_candidate_envelope",
    "validate_protein_identity_manifest",
    "is_json_exact_protein_identity_manifest",
    "validate_protein_raw_entry_references",
]
