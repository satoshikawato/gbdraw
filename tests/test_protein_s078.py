"""S07.8 exact baseline and independent boundary/structural checks."""
import ast
import gzip
import hashlib
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest

from gbdraw.analysis import collinearity_units as cu, protein_colinearity as pc
from tests.prototypes.protein_s078 import frozen, frozen_units
from tests.test_protein_colinearity import _hit_row, _protein_map_for_lengths
from tests.test_collinearity_units import _boundaries, _extraction, _observe


def test_frozen_source_identity_and_isolation():
    for module, live, sha in [
        (frozen, pc, '608abd7f9e0ff266e0746c7a961712c74d358ae85c27653355aa8ccce4fbc17d'),
        (frozen_units, cu, '975873bfd0c4cf685df5790c0bcbcea85c71a1498f27a511e2ce07b069c564da'),
    ]:
        source = gzip.decompress((Path(__file__).parent/'prototypes'/f'{live.__name__.rsplit(".", 1)[-1]}_s077.py.gz').read_bytes())
        assert hashlib.sha256(source).hexdigest() == sha
        for node in ast.parse(source).body:
            if isinstance(node, ast.FunctionDef):
                fn = getattr(module, node.name)
                assert fn is not getattr(live, node.name)
                assert fn.__globals__ is vars(module)
            elif isinstance(node, ast.ClassDef) and node.name.startswith('_'):
                assert getattr(module, node.name) is not getattr(live, node.name)


@pytest.mark.parametrize('interleaved_error', [False, True])
@pytest.mark.parametrize('column', ['bitscore', 'evalue', 'identity', 'alignment_length', 'qstart', 'qend', 'sstart', 'send'])
@pytest.mark.parametrize('value', [float('nan'), float('inf'), -float('inf'), 'bad', None, pd.NA, -0.5, 0.5, '2.5'])
def test_hsp_exact_private_errors_values_and_all_columns(column, value, interleaved_error):
    pm = _protein_map_for_lengths({'q': 100, 's': 100, 't': 100})
    hits = pd.DataFrame([_hit_row('q', 's'), _hit_row('q', 't', alignment_length=float('inf') if interleaved_error else 10),
                         _hit_row('q', 's') | {column: value}])
    hits['extra column'] = ['first', 'second', 'third']
    hits.index = [9, 9, 1]
    try:
        expected = frozen._aggregate_hsps_by_protein_pair(hits, pm)
    except Exception as error:
        with pytest.raises(type(error)) as caught:
            pc._aggregate_hsps_by_protein_pair(hits, pm)
        assert str(caught.value) == str(error)
    else:
        pd.testing.assert_frame_equal(pc._aggregate_hsps_by_protein_pair(hits, pm), expected, check_exact=True)


@pytest.mark.parametrize('mode', ['auto', 'cds', 'locus', 'invalid'])
def test_complete_unit_index_s077(mode, caplog):
    ex = _extraction(_boundaries())
    assert _observe(cu.build_collinearity_unit_index, ex, caplog, mode=mode) == \
           _observe(frozen_units.build_collinearity_unit_index, ex, caplog, mode=mode)


def test_stateful_scalar_reuse_counterexample():
    class Scalar:
        def __init__(self):
            self.calls = 0
        def __float__(self):
            self.calls += 1
            return float(10 * self.calls)
    pm = _protein_map_for_lengths({'q': 100, 's': 100})
    results = []
    for module in (frozen, pc):
        value = Scalar()
        hits = pd.DataFrame([_hit_row('q', 's', alignment_length=value)])
        result = module._aggregate_hsps_by_protein_pair(hits, pm)
        assert value.calls == 3
        assert result.iloc[0].total_hsp_alignment_length == 20
        assert result.iloc[0].representative_alignment_length == 30
        results.append(result.drop(columns='alignment_length'))
    pd.testing.assert_frame_equal(*results, check_exact=True)


def test_hsp_avoids_copying_asdict():
    pm = _protein_map_for_lengths({'q': 100, 's': 100})
    hits = pd.DataFrame([_hit_row('q', 's')])
    calls = []
    for module in (frozen, pc):
        with patch.object(module, 'dict', wraps=dict, create=True) as copied:
            result = module._aggregate_hsps_by_protein_pair(hits, pm)
        calls.append(copied.call_count)
        assert result.hsp_count.tolist() == [1]
    assert calls == [1, 0]


@pytest.mark.parametrize('mode', ['auto', 'cds', 'locus'])
def test_aliases_remain_ordered_for_large_collapsed_locus(mode, caplog):
    from tests.test_collinearity_units import _protein
    from tests.prototypes.collinearity_units_oracle import build
    rows = [[_protein(f'p{i}', start=300-i, locus_tag='one', source_protein_id=f'p{i}',
                      gene=f'gene{i%3}', old_locus_tag='repeat', strand=(-1)**i) for i in range(100)]]
    ex = _extraction(rows)
    expected = _observe(build, ex, caplog, mode=mode)
    assert _observe(frozen_units.build_collinearity_unit_index, ex, caplog, mode=mode) == expected
    assert _observe(cu.build_collinearity_unit_index, ex, caplog, mode=mode) == expected


@pytest.mark.parametrize('scenario', ['sparse-2', 'dense-4', 'support-edges'])
@pytest.mark.parametrize('limit', [0, 1, 200])
def test_related_limits_complete_s077_result(scenario, limit):
    from gbdraw.exceptions import ValidationError
    from tests.test_protein_comparison_benchmark import benchmark
    pm, tables = benchmark.synthetic(pc, scenario, benchmark.SEED)
    def run(module):
        try:
            result = module.select_rbh_orthogroup_edges_from_directional_hits(
                tables, pm, record_count=3, max_related_edges_per_orthogroup=limit)
        except ValidationError as error:
            assert limit == 0
            return type(error), str(error)
        assert limit > 0
        return benchmark.canonical(result)
    assert run(pc) == run(frozen)


def test_unit_key_independent_minima_conversion_order_and_empty():
    from types import SimpleNamespace
    # Error in the second start must precede the first end conversion.
    cases = [[], [SimpleNamespace(start=0, end=None, feature_index=0, protein_id='a'),
                  SimpleNamespace(start='bad', end=2, feature_index=1, protein_id='b')]]
    for proteins in cases:
        with pytest.raises(Exception) as expected:
            frozen_units._unit_sort_key(proteins)
        with pytest.raises(type(expected.value)) as actual:
            cu._unit_sort_key(proteins)
        assert str(actual.value) == str(expected.value)
    proteins = _boundaries()[0][:2]
    assert cu._unit_sort_key(proteins) == frozen_units._unit_sort_key(proteins) == (0, 30, 0, 'a')


@pytest.mark.parametrize('dtype', ['object', 'string', 'category'])
def test_hsp_s077_pair_order_missing_ids_and_extra_columns(dtype):
    pm = _protein_map_for_lengths({'q': 100, 's': 100, 'nan': 100})
    rows = [_hit_row(q, 's', qstart=start, qend=end) | {'extra data': [ordinal]}
            for ordinal, (q, start, end) in enumerate([
                ('q', 20, -5), (None, 1, 20), ('nan', 1, 3), ('q', 21, 40),
                (pd.NA, 1, 2), ('unknown', 1, 5), ('q', 90, 120), ('q', 21, 40)])]
    hits = pd.DataFrame(rows).astype({'query': dtype, 'subject': dtype})
    hits.index = [3]*len(hits)
    expected = frozen._aggregate_hsps_by_protein_pair(hits, pm)
    actual = pc._aggregate_hsps_by_protein_pair(hits, pm)
    pd.testing.assert_frame_equal(actual, expected, check_exact=True)
    assert actual['query'].tolist() == ['q', 'nan']
    assert actual['hsp_count'].tolist() == [4, 1]
    assert actual['query_covered_length'].tolist() == [51, 3]
