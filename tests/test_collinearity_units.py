from __future__ import annotations

from dataclasses import fields, is_dataclass, replace
import logging
import random
from unittest.mock import patch

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from gbdraw.analysis import collinearity_units as units_module
from gbdraw.analysis.collinearity_units import build_collinearity_unit_index
from gbdraw.analysis.protein_colinearity import CdsProtein, ProteinExtractionResult, extract_cds_proteins
from gbdraw.exceptions import ValidationError
from tests.prototypes import collinearity_units_frozen as frozen
from tests.prototypes import collinearity_units_oracle as oracle


def _record(record_id: str, features: list[SeqFeature]) -> SeqRecord:
    record = SeqRecord(Seq("ATGAAATAG" * 20), id=record_id)
    record.features = features
    return record


def _cds(start: int, end: int, qualifiers: dict[str, list[str]]) -> SeqFeature:
    merged = {"translation": ["MK"]}
    merged.update(qualifiers)
    return SeqFeature(FeatureLocation(start, end, strand=1), type="CDS", qualifiers=merged)


@pytest.mark.linear
def test_auto_unit_mode_collapses_strong_locus_ids() -> None:
    record = _record(
        "record_a",
        [
            _cds(0, 9, {"locus_tag": ["locus_a"], "protein_id": ["p1"]}),
            _cds(12, 24, {"locus_tag": ["locus_a"], "protein_id": ["p2"]}),
            _cds(30, 39, {"locus_tag": ["locus_b"], "protein_id": ["p3"]}),
        ],
    )
    extraction = extract_cds_proteins([record])

    unit_index = build_collinearity_unit_index(extraction, records=[record], mode="auto")

    units = unit_index.units_by_record[0]
    assert [unit.unit_kind for unit in units] == ["locus", "locus"]
    assert units[0].locus_id == "locus_a"
    assert units[0].cds_members == ("p1", "p2")
    assert units[0].start == 0
    assert units[0].end == 24


@pytest.mark.linear
def test_gene_labels_do_not_drive_unit_collapse() -> None:
    record = _record(
        "record_a",
        [
            _cds(0, 9, {"gene": ["abc"], "protein_id": ["p1"]}),
            _cds(12, 24, {"gene": ["abc"], "protein_id": ["p2"]}),
        ],
    )
    extraction = extract_cds_proteins([record])

    unit_index = build_collinearity_unit_index(extraction, records=[record], mode="auto")

    units = unit_index.units_by_record[0]
    assert [unit.unit_kind for unit in units] == ["cds", "cds"]
    assert units[0].unit_id != units[1].unit_id
    assert "abc" in unit_index.ambiguous_aliases_by_record[0]


@pytest.mark.linear
def test_locus_unit_mode_rejects_product_only_cds() -> None:
    record = _record(
        "record_a",
        [_cds(0, 9, {"product": ["hypothetical protein"]})],
    )
    extraction = extract_cds_proteins([record])

    with pytest.raises(ValidationError, match="requires stable locus identifiers"):
        build_collinearity_unit_index(extraction, records=[record], mode="locus")


@pytest.mark.linear
def test_cds_unit_mode_never_collapses_same_locus() -> None:
    record = _record(
        "record_a",
        [
            _cds(0, 9, {"locus_tag": ["locus_a"], "protein_id": ["p1"]}),
            _cds(12, 24, {"locus_tag": ["locus_a"], "protein_id": ["p2"]}),
        ],
    )
    extraction = extract_cds_proteins([record])

    unit_index = build_collinearity_unit_index(extraction, records=[record], mode="cds")

    units = unit_index.units_by_record[0]
    assert [unit.unit_kind for unit in units] == ["cds", "cds"]
    assert [unit.representative_protein_id for unit in units] == ["p1", "p2"]

# The frozen builder and independent oracle protect different failure modes.
def _ordered(value):
    if is_dataclass(value):
        return type(value).__name__, tuple((f.name, _ordered(getattr(value, f.name))) for f in fields(value))
    if isinstance(value, dict):
        return 'dict', tuple((_ordered(k), _ordered(v)) for k, v in value.items())
    if isinstance(value, (list, tuple)):
        return type(value).__name__, tuple(map(_ordered, value))
    if isinstance(value, set):
        return 'set', frozenset(value)
    return value


def _protein(pid, **kwargs):
    defaults = dict(protein_id=pid, record_index=0, feature_index=0, record_id='r0',
                    start=10, end=40, strand=1, label=pid, protein_length=10, sequence='M'*10,
                    feature_svg_id=f'svg-{pid}')
    return CdsProtein(**(defaults | kwargs))


def _extraction(rows):
    return ProteinExtractionResult(rows, {p.protein_id: p for row in rows for p in row})


def _observe(builder, extraction, caplog, **kwargs):
    before = _ordered(extraction)
    records_before = repr(kwargs.get('records'))
    caplog.clear()
    with caplog.at_level(logging.INFO):
        try:
            value = builder(extraction, **kwargs)
            outcome = ('ok', _ordered(value))
        except Exception as exc:
            outcome = ('error', type(exc), str(exc))
    assert _ordered(extraction) == before
    assert repr(kwargs.get('records')) == records_before
    return outcome, [(r.levelname, r.getMessage()) for r in caplog.records]


def _parity(extraction, caplog, **kwargs):
    expected = _observe(oracle.build, extraction, caplog, **kwargs)
    assert _observe(frozen.build_collinearity_unit_index, extraction, caplog, **kwargs) == expected
    assert _observe(build_collinearity_unit_index, extraction, caplog, **kwargs) == expected
    return expected


def _boundaries():
    # Independent minima differ from the first member's tuple: B precedes A.
    a = [_protein('z', start=0, end=100, feature_index=9, locus_tag='A'),
         _protein('a', start=20, end=30, feature_index=0, locus_tag='A', strand=-1)]
    b = [_protein('b', start=0, end=90, feature_index=8, locus_tag='B'),
         _protein('c', start=20, end=25, feature_index=1, locus_tag='B', strand=None)]
    # All locus sources are present on the first member; only its parent wins.
    priority = [_protein('parent', gene_parent_id=' parent ', locus_tag='tag', gene_id='id', db_xref=('GeneID:9',)),
                _protein('tag', gene_parent_id=' ', locus_tag=' tag ', gene_id='id', db_xref=('GeneID:9',)),
                _protein('id', locus_tag=' ', gene_id=' id ', db_xref=('GeneID:9',)),
                _protein('xref', db_xref=('GeneID:', 'Other:1', ' GeneID:9 ', 'GeneID:10'))]
    aliases = [_protein(f'p{i}', start=200+i*10, locus_tag=f'L{i}', gene='shared',
                        old_locus_tag='twice' if i < 2 else 'unique', source_protein_id=f'p{i}') for i in range(4)]
    return [a+b+priority+aliases, [], [_protein('other', record_index=2, record_id='r2', gene='shared')]]


@pytest.mark.parametrize('mode', ['auto', 'cds', 'locus', ' AUTO ', '', None, 'invalid'])
@pytest.mark.parametrize('records', [None, [], [_record('override', [])],
                                     [_record('first', []), _record('empty', []), _record('last', [])]])
def test_full_index_oracles_at_boundaries(mode, records, caplog):
    rows = _boundaries()
    random.Random(73).shuffle(rows[0])
    _parity(_extraction(rows), caplog, records=records, mode=mode)


@pytest.mark.parametrize('mode', ['auto', 'cds', 'locus'])
@pytest.mark.parametrize('rows', [[], [[]], [[], []]])
def test_empty_index_oracles(rows, mode, caplog):
    _parity(_extraction(rows), caplog, mode=mode)


@pytest.mark.parametrize('mode', ['auto', 'cds', 'locus'])
@pytest.mark.parametrize('seed', range(20))
def test_shuffled_multicds_full_index_differential(seed, mode, caplog):
    rng = random.Random(seed)
    rows = []
    for record_index in range(3):
        row = [_protein(f'p{record_index}-{i}', record_index=record_index, record_id=f'r{record_index}',
                        start=rng.randrange(5)*10, end=rng.randrange(5)*10+50,
                        feature_index=rng.randrange(4), protein_length=rng.randrange(3)+10,
                        locus_tag=rng.choice([None, 'A', 'B', 'C']),
                        gene=rng.choice(['same', None, 'another']),
                        source_protein_id=rng.choice([None, 'source', ' ', f'p{record_index}-{i}']),
                        strand=rng.choice([-1, 1, None, 0]), old_locus_tag=rng.choice([None, 'old']))
               for i in range(rng.randrange(1, 15))]
        rng.shuffle(row)
        rows.append(row)
    _parity(_extraction(rows), caplog, mode=mode)


@pytest.mark.parametrize('changes,winner', [
    ({'protein_length': 11}, 'b'),
    ({'source_protein_id': 'source'}, 'b'),
    ({'end': 41}, 'b'),
    ({'feature_index': 1}, 'a'),
    ({}, 'a'),
])
def test_each_representative_tiebreak(changes, winner, caplog):
    row = [_protein('b', locus_tag='L', **changes), _protein('a', locus_tag='L')]
    extraction = _extraction([row])
    _parity(extraction, caplog)
    unit, = build_collinearity_unit_index(extraction).units_by_record[0]
    assert unit.representative_protein_id == winner
    assert unit.representative_feature_svg_id == f'svg-{winner}'


@pytest.mark.parametrize('strands,expected', [([-1, -1], -1), ([-1, None], -1),
                                            ([1, -1], None), ([None, 0], None), ([1, None], 1)])
def test_strand_reduction(strands, expected, caplog):
    ex = _extraction([[_protein(str(i), locus_tag='L', strand=s) for i, s in enumerate(strands)]])
    _parity(ex, caplog)
    assert build_collinearity_unit_index(ex).units_by_record[0][0].strand == expected


def test_independent_minima_and_alias_collision_order(caplog):
    ex = _extraction(_boundaries())
    _parity(ex, caplog)
    result = build_collinearity_unit_index(ex)
    assert [u.locus_id for u in result.units_by_record[0]][:2] == ['B', 'A']
    assert {'shared', 'twice', 'unique'} <= result.ambiguous_aliases_by_record[0]
    assert result.aliases_by_record[2]['shared'] == result.units_by_record[2][0].unit_id
    assert result.units_by_record[2][0].unit_id.endswith(f'{len(result.unit_by_id):06d}')


def test_duplicate_protein_mapping_keeps_first_key_position_last_value(caplog):
    ex = _extraction([[_protein('dup', start=0, locus_tag='A'), _protein('middle', start=2),
                       _protein('dup', start=4, locus_tag='B')]])
    _parity(ex, caplog)
    result = build_collinearity_unit_index(ex)
    assert list(result.unit_by_protein_id) == ['dup', 'middle']
    assert result.unit_by_protein_id['dup'].locus_id == 'B'
    assert 'dup' in result.ambiguous_aliases_by_record[0]


def test_repeated_alias_within_one_unit_is_not_ambiguous(monkeypatch, caplog):
    ex = _extraction([[_protein('same', locus_tag='same', source_protein_id='same', gene='same')]])
    _parity(ex, caplog)
    original = units_module._unit_aliases
    monkeypatch.setattr(units_module, '_unit_aliases', lambda **kw: original(**kw)*3)
    result = build_collinearity_unit_index(ex)
    assert result.ambiguous_aliases_by_record == [set()]
    assert result.aliases_by_record[0]['same'] == result.units_by_record[0][0].unit_id


@pytest.mark.parametrize('field,value', [('start', 'bad'), ('end', None), ('feature_index', 'bad'),
                                        ('protein_length', float('inf'))])
def test_invalid_field_exception_parity(field, value, caplog):
    ex = _extraction([[replace(_protein('p'), **{field: value})]])
    expected = _observe(frozen.build_collinearity_unit_index, ex, caplog)
    assert expected[0][0] == 'error'
    assert _observe(build_collinearity_unit_index, ex, caplog) == expected


def test_warning_and_exception_preview_limits(caplog):
    ex = _extraction([[_protein(f'p{r}-{i}', record_id=f'r{r}', start=i) for i in range(8)] for r in range(7)])
    for mode in ('auto', 'locus'):
        _parity(ex, caplog, mode=mode, records=[])


@pytest.mark.parametrize('mode', ['auto', 'cds', 'locus'])
def test_unit_allocation_and_removed_work(mode):
    ex = _extraction([[_protein(str(i), start=20-i,
                     locus_tag=None if mode != 'locus' and i % 5 == 0 else f'L{i%2}')
                     for i in range(20)], []])
    with patch.object(units_module, 'CollinearityUnit', wraps=units_module.CollinearityUnit) as construct, \
         patch.object(units_module, 'sorted', wraps=sorted, create=True) as sort, \
         patch.object(units_module, 'set', wraps=set, create=True) as sets, \
         patch.object(units_module, 'strong_locus_id', wraps=units_module.strong_locus_id) as resolve_locus, \
         patch.object(units_module, '_add_alias', wraps=units_module._add_alias) as add_alias:
        result = build_collinearity_unit_index(ex, mode=mode)
    assert construct.call_count == len(result.unit_by_id)
    protein_count = sum(map(len, ex.proteins_by_record))
    assert resolve_locus.call_count == protein_count
    assert add_alias.call_count == 7*protein_count + 3*len(result.unit_by_id)
    # Only the original protein and group sorts, even for a large collapsed locus.
    assert sort.call_count == 2*len(ex.proteins_by_record)
    # One returned ambiguity set per record; no alias-target temporary sets.
    assert sets.call_count == len(ex.proteins_by_record)
    assert all(call.kwargs['aliases'] for call in construct.call_args_list)


@pytest.mark.parametrize('scope', ['adjacent', 'all'])
@pytest.mark.parametrize('infer', [False, True])
@pytest.mark.parametrize('limit', [1, 5, None])
@pytest.mark.parametrize('mode', ['auto', 'cds', 'locus'])
@pytest.mark.parametrize('anchor', ['rbh', 'one_to_one', 'all'])
@pytest.mark.parametrize('scenario,reverse', [('sparse-2', False), ('dense-4', True), ('support-edges', False)])
def test_complete_collinear_outputs_with_frozen_units(scope, infer, limit, mode, anchor, scenario, reverse):
    from gbdraw.analysis import collinearity as cc
    from gbdraw.analysis import protein_colinearity as pc
    from gbdraw.session_request_codec import encode_canonical_typed_resource
    from tests.test_protein_comparison_benchmark import benchmark

    pm, tables = benchmark.synthetic(pc, scenario, 20260916)
    pm = {pid: replace(p, locus_tag=f'L{p.feature_index//2}',
                      start=4000-p.end if reverse else p.start,
                      end=4000-p.start if reverse else p.end,
                      strand=-1 if reverse else 1) for pid, p in pm.items()}
    ex = _extraction([[p for p in pm.values() if p.record_index == r] for r in range(3)])
    records = [SeqRecord(Seq('A'*5000), id=f'record_{r}') for r in range(3)]
    before = benchmark.canonical((ex, tables))
    def run():
        result = cc.build_orthogroup_collinearity_blocks_from_hits(tables, ex, records=records,
                    search_scope=scope, infer_orthogroups=infer, orthogroup_member_max_hits=limit,
                    unit_mode=mode, edge_mode=anchor)
        return benchmark.canonical((result,
                   cc.convert_collinearity_blocks_to_pair_comparisons(result, records=records),
                   encode_canonical_typed_resource('result', result)))
    with patch.object(cc, 'build_collinearity_unit_index', frozen.build_collinearity_unit_index):
        expected = run()
    assert run() == expected
    assert benchmark.canonical((ex, tables)) == before
