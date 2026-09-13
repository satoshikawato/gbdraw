"""Portable regressions for the retained runner's operation adapters."""
import runpy
from pathlib import Path
from types import SimpleNamespace

import pytest


RUNNER = Path(__file__).parent / 'web/helpers/retained-visual-journeys.py'
install = runpy.run_path(str(RUNNER))['install_operation_adapters']


@pytest.fixture
def journey():
    class Journey:
        jid = 'J42'

        def __init__(self):
            self.calls = []
            self.events = []
            self.pending = True
            self.wait_error = None
            self.observation_error = None
            self.comparison_bytes = 100
            self.failures = []
            self.page = SimpleNamespace(wait_for_function=self.wait)

        def wait(self, predicate, **options):
            self.calls.append(('wait', predicate, options))
            if self.wait_error:
                raise self.wait_error
            self.pending = False

        def load(self, seed, valid=True):
            self.calls.append(('load', seed, valid))
            return 'accepted' if valid else 'rejected'

        def edit(self, kind, *args, **kwargs):
            self.calls.append(('pointer edit', kind, args, kwargs))
            self.observe('before drawer/popup close')
            self.calls.append(('close',))
            self.pending = False
            self.observe('after completed ' + kind + ' edit')
            return {'svg_id': 'target'}

        def observe(self, reason):
            self.calls.append(('observe', reason))
            if self.observation_error:
                raise self.observation_error
            state = {'resultCount': 1, 'failures': self.failures,
                     'observation_phase': 'IN_PROGRESS' if self.pending else 'SETTLED'}
            if not self.pending:
                state.update({key: {'uncompressed_bytes': self.comparison_bytes} for key in [
                    'selected_semantics', 'mounted_semantics', 'export_semantics',
                    'export_boundary_selected_semantics', 'export_boundary_mounted_semantics']})
            self.record('authority', 'PASS', phase=state['observation_phase'])
            return state

        def record(self, kind, status, **data):
            event = {'kind': kind, 'status': status, **data}
            self.events.append(event)
            return event

    install(Journey)
    return Journey()


@pytest.mark.parametrize('args,kwargs,valid', [
    (('seed',), {}, True),
    (('seed', True), {}, True),
    (('seed', False), {}, False),
    (('seed',), {'valid': False}, False),
])
def test_load_preserves_original_valid_argument(journey, args, kwargs, valid):
    assert journey.load(*args, **kwargs) == ('accepted' if valid else 'rejected')
    assert journey.calls == [('load', 'seed', valid)]


def test_c01_waits_for_definition_callback_after_valid_load(journey):
    journey.jid = 'C01'
    assert journey.load('C01-03-save.json') == 'accepted'
    assert [call[0] for call in journey.calls] == ['load', 'wait', 'observe']
    assert journey.calls[1][1] == 'window.__DEFINITION_COMPLETED__ > 0'


def test_color_comparison_precedes_close_and_preserves_pointer_arguments(journey):
    target = journey.edit('color', '#123456', record=2)
    assert target == {'svg_id': 'target'}
    assert [call[0] for call in journey.calls] == ['pointer edit', 'wait', 'observe', 'close', 'observe']
    assert journey.calls[0] == ('pointer edit', 'color', ('#123456',), {'record': 2})
    assert any(row['kind'] == 'completion_boundary' and row['status'] == 'PASS' for row in journey.events)


@pytest.mark.parametrize('failure', ['timeout', 'observer', 'mismatch', 'empty', 'missing'])
def test_unverified_color_boundary_cannot_close_or_pass(journey, failure):
    if failure == 'timeout':
        journey.wait_error = TimeoutError('color owner still pending')
    elif failure == 'observer':
        journey.observation_error = RuntimeError('capture failed')
    elif failure == 'mismatch':
        journey.failures = ['selected / mounted differ']
    elif failure == 'empty':
        journey.comparison_bytes = 2  # [] is not a visual comparison.
    else:
        journey.wait = lambda *args, **kwargs: None
        journey.page.wait_for_function = journey.wait  # Observation is still IN_PROGRESS.
    with pytest.raises((AssertionError, TimeoutError, RuntimeError)):
        journey.edit('color')
    assert ('close',) not in journey.calls
    assert not any(row['kind'] == 'completion_boundary' and row['status'] == 'PASS' for row in journey.events)
    # A failed edit must not impose its wait on a later independent operation.
    journey.wait_error = journey.observation_error = None
    journey.calls.clear()
    journey.observe('before drawer/popup close')
    assert not any(call[0] == 'wait' for call in journey.calls)


def test_in_progress_observation_is_explicit_without_blocking_async_generate(journey):
    journey.observe('after original statement starting Generate')
    assert journey.events[-1]['status'] == 'IN_PROGRESS'
    assert not any(call[0] == 'wait' for call in journey.calls)
    assert journey.record('authority', 'OBSERVATION', phase='IN_PROGRESS')['status'] == 'OBSERVATION'


def test_mobile_label_pointer_path_does_not_gain_color_wait(journey):
    journey.jid = 'J43'
    journey.edit('label', 'label text', keyboard=False)
    assert not any(call[0] == 'wait' for call in journey.calls)
    assert journey.calls[0] == ('pointer edit', 'label', ('label text',), {'keyboard': False})
