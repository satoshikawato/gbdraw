"""Exact S07.6 differential plus independent oracles for local reductions."""
import math
from dataclasses import replace
import random
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pandas.testing import assert_frame_equal

from gbdraw.analysis import collinearity as cc, protein_colinearity as pc
from tests.prototypes.protein_local_reductions import frozen, fit_oracle, hit_ordinals, frozen_collinear_callers
from tests.test_protein_comparison_benchmark import benchmark
from gbdraw.session_request_codec import encode_canonical_typed_resource


def fit_frame(size, ids='string', duplicate_index=False):
    frame = pd.DataFrame({
        'query': [f'q{i % 3}' for i in range(size)],
        'subject': [f's{i % 2}' for i in range(size)],
        'length_product': [float(100 + (i % 7) * 17) for i in range(size)],
        'bitscore': [float(10 + (i % 4) * 3) for i in range(size)],
        **{f'extra_{j}': [f'{i}-{j}' for i in range(size)] for j in range(20)},
    }, index=[i // 2 if duplicate_index else i * 3 + 7 for i in range(size)])
    for column in ('query', 'subject'):
        if ids == 'numeric':
            frame[column] = [int(v[1:]) for v in frame[column]]
        elif ids == 'categorical':
            frame[column] = pd.Categorical(frame[column], categories=sorted(set(frame[column]), reverse=True))
    return frame


@pytest.mark.parametrize('size', [0, 1, 3, 4, 7, 8, 11, 12, 15, 16, 31, 32, 33, 65])
@pytest.mark.parametrize('fraction', [.05, .25, .26, .5, 1.])
@pytest.mark.parametrize('ids', ['string', 'numeric', 'categorical'])
def test_fit_points_model_ordinals_and_input(size, fraction, ids):
    frame = fit_frame(size, ids)
    before = frame.copy(deep=True)
    ordinals, points, model = fit_oracle(frame, fraction)
    assert pc._normalized_fit_row_positions(frame, top_fraction=fraction) == ordinals
    actual = pc._select_normalized_fit_rows(frame, top_fraction=fraction)
    assert actual == points == frozen._select_normalized_fit_rows(frame, top_fraction=fraction)
    assert pc._fit_expected_bitscore_model(frame, top_fraction=fraction) == model
    assert model == frozen._fit_expected_bitscore_model(frame, top_fraction=fraction)
    assert_frame_equal(frame, before, check_exact=True)
    frame.index = [i // 2 for i in range(size)]
    for col in frame.columns[4:]:
        frame[col] = None
    assert pc._select_normalized_fit_rows(frame, top_fraction=fraction) == points


@pytest.mark.parametrize('kind', ['full_tie', 'zero', 'round_tie', 'epsilon', 'few_points', 'nonfinite'])
def test_fit_degenerate_models(kind):
    frame = fit_frame(16)
    if kind == 'full_tie':
        frame.loc[:, ['query', 'subject']] = 'tie'
        frame['length_product'] = 100.
        frame['bitscore'] = 20.
    elif kind == 'zero':
        frame['length_product'] = 0.
    elif kind == 'round_tie':
        frame['length_product'] = [100. + i * 1e-10 for i in range(16)]
    elif kind == 'epsilon':
        frame['length_product'] = [100. + i * 1e-6 for i in range(16)]
    elif kind == 'few_points':
        frame['bitscore'] = [1.] + [0.] * 15
    else:
        frame['bitscore'] = float('inf')
    ordinals, expected, model = fit_oracle(frame, .25)
    assert model is None
    assert pc._normalized_fit_row_positions(frame, top_fraction=.25) == ordinals
    assert pc._select_normalized_fit_rows(frame, top_fraction=.25) == expected
    if kind == 'full_tie':
        assert ordinals == [0, 4, 8, 12]
    assert pc._fit_expected_bitscore_model(frame, top_fraction=.25) is None


def test_fit_sorts_the_four_fit_columns_in_two_passes():
    frame = fit_frame(65)
    widths = []
    original = pd.DataFrame.sort_values
    def sort(df, *args, **kwargs):
        widths.append(tuple(df.columns))
        return original(df, *args, **kwargs)
    with patch.object(pd.DataFrame, 'sort_values', sort):
        pc._select_normalized_fit_rows(frame, top_fraction=.05)
    fit_columns = ('length_product', 'query', 'subject', 'bitscore')
    assert widths == [fit_columns, (*fit_columns, '_bin')]


@pytest.mark.parametrize('size', [0, 1, 11, 12, 33])
@pytest.mark.parametrize('fraction', [.05, .26, 1.])
def test_normalized_table_all_columns_and_exact_scores(size, fraction):
    pm = {}
    rows = []
    for i in range(size):
        q, s = f'q{i}', f's{i}'
        pm[q] = benchmark.protein(pc, q, 0, i, length=100 + i * 3)
        pm[s] = benchmark.protein(pc, s, 1, i, length=110 + i * 5)
        rows.append(benchmark.hit(q, s, bitscore=100. + i * 7, extra=f'keep-{i}'))
    hits = pd.DataFrame(rows, columns=[*pc.COMPARISON_COLUMNS, 'extra'])
    hits.index = [i // 2 for i in range(size)]
    before = hits.copy(deep=True)
    actual = pc._normalize_directional_hit_table(hits, pm, top_fraction=fraction)
    expected = frozen._normalize_directional_hit_table(hits, pm, top_fraction=fraction)
    assert_frame_equal(actual, expected, check_exact=True)
    assert_frame_equal(hits, before, check_exact=True)
    if size:
        assert actual['extra'].tolist() == hits['extra'].tolist()
        _, _, model = fit_oracle(actual, fraction)
        scores = []
        for length, score in zip(actual['length_product'], actual['bitscore']):
            denominator = (math.sqrt(max(length, 1e-12)) if model is None else
                           10 ** (model[0] * math.log10(length) + model[1]))
            scores.append(score / denominator)
        assert actual['normalized_score'].tolist() == scores
        assert actual['normalization_fallback'].tolist() == [model is None] * size


def member_frame(ids='string'):
    rows = [benchmark.hit(q, s, bitscore=str(score), extra=f'row-{i}')
            for i, (q, s, score) in enumerate([
                ('q1', 's1', 9), ('q0', 's2', 100), ('q0', 's1', 100),
                ('q1', 's0', 200), ('q0', 's1', 80), ('q1', 's0', 200),
                ('q0', 's0', 50), ('q1', 's1', 90)])]
    frame = pd.DataFrame(rows, index=[8, 3, 3, 1, 9, 4, 0, 0])
    frame['mismatches'] = frame['mismatches'].astype('Int64')
    for column in ('query', 'subject'):
        if ids == 'numeric':
            frame[column] = [int(v[1:]) for v in frame[column]]
        elif ids == 'categorical':
            frame[column] = pd.Categorical(frame[column], categories=sorted(set(frame[column]), reverse=True))
    return frame


def numeric_frame(frame):
    result = frame.copy()
    for column in pc.COMPARISON_COLUMNS[2:]:
        result[column] = pd.to_numeric(result[column])
    return result


@pytest.mark.parametrize('ids', ['string', 'numeric', 'categorical'])
@pytest.mark.parametrize('limit', [1, 5, None])
def test_member_top_best_and_reciprocity_independent_oracle(ids, limit):
    hits = member_frame(ids)
    before = hits.copy(deep=True)
    if limit is None:
        positions = list(range(len(hits)))
    else:
        chosen = hit_ordinals(hits, limit)
        pairs = {(hits['query'].iloc[i], hits['subject'].iloc[i]) for i in chosen}
        positions = [i for i, pair in enumerate(zip(hits['query'], hits['subject'])) if pair in pairs]
        top = pc.select_top_hits_per_query(hits, max_hits=limit)
        assert_frame_equal(top, numeric_frame(hits).iloc[chosen].reset_index(drop=True), check_exact=True)
        assert_frame_equal(top, frozen.select_top_hits_per_query(hits, max_hits=limit), check_exact=True)
    members = pc._select_member_candidate_hits_per_query(hits, max_hits=limit)
    assert_frame_equal(members, hits.iloc[positions].reset_index(drop=True), check_exact=True)
    assert_frame_equal(members, frozen._select_member_candidate_hits_per_query(hits, max_hits=limit), check_exact=True)
    query_best = hit_ordinals(members)
    actual = pc.select_best_hits_per_query(members)
    assert_frame_equal(actual, numeric_frame(members).iloc[query_best].reset_index(drop=True), check_exact=True)
    assert_frame_equal(actual, frozen.select_best_hits_per_query(members), check_exact=True)
    subject_best = hit_ordinals(members, subject_first=True)
    reciprocal = {(str(members['query'].iloc[i]), str(members['subject'].iloc[i])) for i in subject_best}
    mutual = [i for i in query_best if (str(members['query'].iloc[i]), str(members['subject'].iloc[i])) in reciprocal]
    assert_frame_equal(pc.select_reciprocal_best_hits(members),
                       numeric_frame(members).iloc[mutual].reset_index(drop=True), check_exact=True)
    reverse = members.rename(columns={'query': 'subject', 'subject': 'query'}).iloc[::-1].copy()
    reverse_best = hit_ordinals(reverse)
    back = {str(reverse['query'].iloc[i]): str(reverse['subject'].iloc[i]) for i in reverse_best}
    mutual = [i for i in query_best if back.get(str(members['subject'].iloc[i])) == str(members['query'].iloc[i])]
    assert_frame_equal(pc.select_reciprocal_best_hit_edges(members, reverse),
                       numeric_frame(members).iloc[mutual].reset_index(drop=True), check_exact=True)
    assert_frame_equal(pc.select_reciprocal_best_hit_edges(members, reverse),
                       frozen.select_reciprocal_best_hit_edges(members, reverse), check_exact=True)
    assert_frame_equal(hits, before, check_exact=True)
    if limit is None:
        # Keeping every HSP copies no column data; callers only read the result.
        assert np.shares_memory(members['extra'].to_numpy(), hits['extra'].to_numpy())
    else:
        members.iloc[0, members.columns.get_loc('extra')] = 'changed'
        assert_frame_equal(hits, before, check_exact=True)


@pytest.mark.parametrize('mode', ['rbh', 'one_to_one', 'all'])
@pytest.mark.parametrize('missing', [False, True])
@pytest.mark.parametrize('limit', [1, 5, None])
def test_edge_caller_missing_direction_and_member_limit(mode, missing, limit):
    hits = member_frame()
    tables = {(0, 1): hits}
    if not missing:
        tables[1, 0] = hits.rename(columns={'query': 'subject', 'subject': 'query'}).iloc[::-1]
    def run():
        selected = {pair: cc._select_member_candidate_hits_per_query(h, max_hits=limit) for pair, h in tables.items()}
        return cc._select_orthogroup_edges(selected, edge_mode=mode)
    with frozen_collinear_callers():
        expected = run()
    actual = run()
    assert list(actual) == list(expected)
    for pair in actual:
        assert_frame_equal(actual[pair], expected[pair], check_exact=True)


@pytest.mark.parametrize('column', list(pc.COMPARISON_COLUMNS[2:]))
@pytest.mark.parametrize('bad', ['bad', None, float('nan'), float('inf'), -float('inf')])
def test_public_numeric_rejection_exact_type_message(column, bad):
    hits = member_frame()
    hits[column] = hits[column].astype(object)
    hits.loc[8, column] = bad
    for name, args, kwargs in [
        ('select_top_hits_per_query', (hits,), {'max_hits': 5}),
        ('_select_member_candidate_hits_per_query', (hits,), {'max_hits': 1}),
        ('select_best_hits_per_query', (hits,), {}),
        ('select_reciprocal_best_hits', (hits,), {}),
        ('select_reciprocal_best_hit_edges', (hits, member_frame()), {}),
        ('select_reciprocal_best_hit_edges', (member_frame(), hits), {}),
    ]:
        with pytest.raises(Exception) as old_error:
            getattr(frozen, name)(*args, **kwargs)
        with pytest.raises(type(old_error.value)) as new_error:
            getattr(pc, name)(*args, **kwargs)
        assert str(new_error.value) == str(old_error.value)
    # Unlimited member selection returns original rows, including raw numeric fields.
    assert_frame_equal(pc._select_member_candidate_hits_per_query(hits, max_hits=None),
                       frozen._select_member_candidate_hits_per_query(hits, max_hits=None), check_exact=True)


@pytest.mark.parametrize('shape', ['empty', 'missing_columns', 'empty_missing_columns'])
def test_selector_empty_and_validation_order(shape):
    hits = member_frame()
    if 'empty' in shape:
        hits = hits.iloc[:0]
    if 'missing_columns' in shape:
        hits = hits.drop(columns='bitscore')
    for name in ['select_top_hits_per_query', '_select_member_candidate_hits_per_query',
                 'select_best_hits_per_query', 'select_reciprocal_best_hits', 'select_reciprocal_best_hit_edges']:
        limits = [1, 5, None, 0, -1, 'bad'] if 'per_query' in name and name != 'select_best_hits_per_query' else [None]
        for limit in limits:
            kwargs = {'max_hits': limit} if name in ['select_top_hits_per_query', '_select_member_candidate_hits_per_query'] else {}
            args = (hits, hits) if name == 'select_reciprocal_best_hit_edges' else (hits,)
            try:
                expected = getattr(frozen, name)(*args, **kwargs)
            except Exception as old_error:
                with pytest.raises(type(old_error)) as new_error:
                    getattr(pc, name)(*args, **kwargs)
                assert str(new_error.value) == str(old_error)
            else:
                assert_frame_equal(getattr(pc, name)(*args, **kwargs), expected, check_exact=True)


def test_removed_member_copies_pair_dedup_and_full_row_conversion():
    hits = member_frame()
    original_copy, original_getitem = pd.DataFrame.copy, pd.DataFrame.__getitem__
    def counts(module, limit):
        copied, projections = [], []
        def copy(df, *args, **kwargs):
            copied.append(tuple(df.columns))
            return original_copy(df, *args, **kwargs)
        def getitem(df, key):
            if isinstance(key, list):
                projections.append(tuple(key))
            return original_getitem(df, key)
        with patch.object(pd.DataFrame, 'copy', copy), patch.object(pd.DataFrame, '__getitem__', getitem):
            module._select_member_candidate_hits_per_query(hits, max_hits=limit)
        return copied, projections
    assert len(counts(frozen, None)[0]) - len(counts(pc, None)[0]) == 1
    assert counts(frozen, 1)[1].count(('query', 'subject')) == 2
    assert counts(pc, 1)[1].count(('query', 'subject')) == 0
    original_dedupe = pd.DataFrame.drop_duplicates
    for name, expected_removed in [('select_best_hits_per_query', 1), ('select_reciprocal_best_hits', 2),
                                    ('select_reciprocal_best_hit_edges', 2)]:
        observed = []
        for module in [frozen, pc]:
            calls = []
            def dedupe(df, *args, **kwargs):
                calls.append(args)
                return original_dedupe(df, *args, **kwargs)
            with patch.object(pd.DataFrame, 'drop_duplicates', dedupe):
                args = (hits, hits) if name == 'select_reciprocal_best_hit_edges' else (hits,)
                getattr(module, name)(*args)
            observed.append(len(calls))
        assert observed[0] - observed[1] == expected_removed
    with patch.object(pd.DataFrame, 'itertuples', side_effect=AssertionError('full row conversion')):
        pc.select_reciprocal_best_hits(hits)
        pc.select_reciprocal_best_hit_edges(hits, hits)


def independent_evidence_rank(row, proteins):
    return (-float(row.normalized_score), float(row.evalue), -float(row.min_coverage),
            -float(row.identity), -int(row.alignment_length),
            proteins[row.query].record_index, row.query,
            proteins[row.subject].record_index, row.subject)


def test_stable_evidence_rank_and_first_occurrence_oracle():
    proteins = {pid: benchmark.protein(pc, pid, record)
                for pid, record in [('b', 1), ('a', 0), ('c', 1), ('d', 0), ('e', 2)]}
    pairs = [('a', 'c'), ('a', 'b'), ('a', 'd'), ('b', 'a'), ('a', 'e'), ('d', 'a')]
    rows = [(pair, SimpleNamespace(**benchmark.hit(*pair), normalized_score=1., min_coverage=1.))
            for pair in pairs]
    hits = pc._directed_hits_from_tables(
        {(0, 1): pd.DataFrame([vars(row) | {'domain_only': False} for _, row in rows])}, proteins)
    # The table holds the rows in rank order, and each row is its record's rank.
    ranked = [hits.hit(row) for row in range(len(hits))]
    assert [(hit.query, hit.subject) for hit in ranked] == [
        pair for pair, _ in sorted(rows, key=lambda item: independent_evidence_rank(item[1], proteins))]
    assert [hit.rank for hit in ranked] == list(range(len(pairs)))
    for hit, reverse in zip(ranked, hits.reverse.tolist()):
        assert (reverse >= 0) == ((hit.subject, hit.query) in pairs)
        if reverse >= 0:
            assert (ranked[reverse].query, ranked[reverse].subject) == (hit.subject, hit.query)
    # Every per-bucket first pick (best hit into a record, first row of a pair)
    # is the earliest position of its key; brute-force it with heavy ties.
    rng = random.Random(20261025)
    for size in [0, 1, 2, 7, 50, 300]:
        keys = np.array([rng.choice([-3, 0, 5, 2**40, 7]) for _ in range(size)], dtype=np.int64)
        codes, first = pc._first_occurrences(keys)
        expected = sorted({key: position for position, key in reversed(list(enumerate(keys.tolist())))}.values())
        assert first.tolist() == expected
        assert keys[first].tolist() == list(dict.fromkeys(keys.tolist()))
        assert (keys[first][codes] == keys).all()


@pytest.mark.parametrize('scenario', ['sparse-2', 'dense-4', 'support-edges'])
def test_evidence_is_ranked_once_from_precomputed_ranks(scenario):
    pm, tables = benchmark.synthetic(pc, scenario, 20260916)
    expected = benchmark.canonical(frozen.select_rbh_orthogroup_edges_from_directional_hits(tables, pm, record_count=3))
    sorts = []
    def sort(values, *args, **kwargs):
        values = list(values)
        if values and all(isinstance(item, pc._DirectedHit) for item in values):
            sorts.append(len(values))
        return sorted(values, *args, **kwargs)
    with patch.object(pc, 'sorted', sort, create=True):
        result = pc.select_rbh_orthogroup_edges_from_directional_hits(tables, pm, record_count=3)
    assert benchmark.canonical(result) == expected
    # The table is built in rank order once; no membership phase sorts hit records.
    assert sorts == []


@pytest.mark.parametrize('baseline', ['s076', 's077'])
@pytest.mark.parametrize('scope', ['adjacent', 'all'])
@pytest.mark.parametrize('infer', [False, True])
@pytest.mark.parametrize('limit', [1, 5, None])
@pytest.mark.parametrize('mode', ['auto', 'cds', 'locus'])
@pytest.mark.parametrize('anchor', ['rbh', 'one_to_one', 'all'])
@pytest.mark.parametrize('scenario,reverse', [('sparse-2', False), ('dense-4', True), ('support-edges', False)])
def test_complete_collinear_result_differential(scope, infer, limit, mode, anchor, scenario, reverse, baseline):
    from tests.prototypes import protein_s078
    baseline_callers = frozen_collinear_callers if baseline == "s076" else protein_s078.frozen_collinear_callers
    pm, tables = benchmark.synthetic(pc, scenario, 20260916)
    pm = {pid: replace(p, locus_tag=f'L{p.feature_index // 2}',
                      start=4000 - p.end if reverse else p.start,
                      end=4000 - p.start if reverse else p.end,
                      strand=-1 if reverse else 1) for pid, p in pm.items()}
    ex = pc.ProteinExtractionResult([[p for p in pm.values() if p.record_index == r] for r in range(3)], pm)
    records = [SeqRecord(Seq('A' * 5000), id=f'record_{r}') for r in range(3)]
    before = benchmark.canonical((ex, tables))
    def run():
        result = cc.build_orthogroup_collinearity_blocks_from_hits(tables, ex, records=records,
            search_scope=scope, infer_orthogroups=infer, orthogroup_member_max_hits=limit,
            unit_mode=mode, edge_mode=anchor)
        return benchmark.canonical((result, cc.convert_collinearity_blocks_to_pair_comparisons(result, records=records),
                                    encode_canonical_typed_resource('result', result)))
    with baseline_callers():
        expected = run()
    assert run() == expected
    assert benchmark.canonical((ex, tables)) == before


@pytest.mark.parametrize('baseline', ['s076', 's077'])
@pytest.mark.parametrize('limit', [1, 5, None])
@pytest.mark.parametrize('scenario', ['sparse-2', 'dense-4', 'support-edges'])
@pytest.mark.parametrize('reverse', [False, True])
def test_complete_similarity_result_differential(limit, scenario, reverse, baseline):
    from tests.prototypes import protein_s078
    baseline_module = frozen if baseline == "s076" else protein_s078.frozen
    pm, tables = benchmark.synthetic(pc, scenario, 20260916)
    if reverse:
        pm = {pid: replace(p, start=4000 - p.end, end=4000 - p.start, strand=-1) for pid, p in pm.items()}
    before = benchmark.canonical((pm, tables))
    def run(module):
        result = module.select_rbh_orthogroup_edges_from_directional_hits(tables, pm, record_count=3,
            orthogroup_member_max_hits=limit, include_singletons=True)
        return benchmark.canonical((result, encode_canonical_typed_resource('result', result.orthogroups)))
    assert run(pc) == run(baseline_module)
    assert benchmark.canonical((pm, tables)) == before
