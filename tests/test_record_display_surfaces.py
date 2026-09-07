"""Public rotation adapters converge on the accepted typed planner."""

from dataclasses import replace

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

import gbdraw
from gbdraw.api import RecordDisplayOptions
from gbdraw.api.record_planning import record_input_manifest_from_table
from gbdraw.api.requests import CircularDiagramRequest, InMemoryRecordSource, RecordInput
from gbdraw.exceptions import ValidationError
from gbdraw.session_request_codec import encode_canonical_request
from gbdraw.session_io import load_session


def record():
    result = SeqRecord(Seq("ACGT" * 30), id="duplicate", annotations={
        "molecule_type": "DNA", "topology": "circular",
    })
    result.features = [SeqFeature(SimpleLocation(10, 35, strand=1), type="CDS")]
    return result


def test_root_exports_the_typed_display_class():
    assert gbdraw.RecordDisplayOptions is RecordDisplayOptions


@pytest.mark.parametrize("mode", ["circular", "linear"])
@pytest.mark.parametrize("start", [None, 1, 120, 25])
def test_root_display_reaches_the_typed_renderer(mode, start):
    source = record()
    before = str(source.seq), str(source.features[0].location)
    diagram = getattr(gbdraw, f"draw_{mode}")(
        source, record_displays=[RecordDisplayOptions(start_coordinate=start)],
    )
    assert "<svg" in diagram._drawing.tostring()
    assert (str(source.seq), str(source.features[0].location)) == before


def test_records_table_preserves_order_and_per_instance_display(tmp_path):
    path = tmp_path / "input.gbk"
    SeqIO.write([record(), record()], path, "genbank")
    table = tmp_path / "records.tsv"
    table.write_text("gbk\trecord_id\torder\ttopology\tdisplay_start\n"
                     "input.gbk\t#1\t2\tauto\t1\n"
                     "input.gbk\t#2\t1\tcircular\t25\n")
    manifest = record_input_manifest_from_table(str(table))
    assert [r.display for r in manifest.records] == [
        RecordDisplayOptions(True, 25), RecordDisplayOptions(None, 1),
    ]
    assert [r.selector.record_index for r in manifest.records] == [1, 0]


@pytest.mark.parametrize("display", [RecordDisplayOptions(True), RecordDisplayOptions(False),
                                     RecordDisplayOptions(start_coordinate=1)])
def test_current_writer_preserves_raw_display(display):
    request = CircularDiagramRequest(records=(RecordInput(InMemoryRecordSource(record())),))
    request = replace(request, records=(replace(request.records[0], display=display),))
    payload = encode_canonical_request(request).payload
    assert payload["records"][0]["display"] == {
        "isCircular": display.is_circular, "startCoordinate": display.start_coordinate,
    }


@pytest.mark.parametrize('mode', ['circular', 'linear'])
@pytest.mark.parametrize('values', ['bad', b'bad', {}, [], [None], [RecordDisplayOptions()] * 2])
def test_root_validates_display_sequence(mode, values):
    with pytest.raises(ValidationError, match='record_displays'):
        getattr(gbdraw, f'draw_{mode}')(record(), record_displays=values)


@pytest.mark.parametrize('mode', ['circular', 'linear'])
@pytest.mark.parametrize('topology,override,start', [
    ('linear', None, 1), ('unknown', None, 25), ('circular', None, 121),
])
def test_root_uses_shared_biological_validation(mode, topology, override, start):
    source = record()
    source.annotations['topology'] = topology
    with pytest.raises(ValidationError):
        getattr(gbdraw, f'draw_{mode}')(source, record_displays=[RecordDisplayOptions(override, start)])


@pytest.mark.parametrize('mode', ['circular', 'linear'])
@pytest.mark.parametrize('start', [None, 1, 120, 25])
def test_root_and_typed_svg_parity(mode, start):
    from gbdraw.api.request_render import build_request_diagram
    from gbdraw.api.requests import LinearDiagramRequest
    import gbdraw.interface as interface
    source = record()
    display = RecordDisplayOptions(start_coordinate=start)
    options = getattr(gbdraw, f'{mode.title()}Options')()
    typed_options = getattr(interface, f'_{mode}_options')(options, record_count=1)
    request_type = CircularDiagramRequest if mode == 'circular' else LinearDiagramRequest
    typed = build_request_diagram(request_type(
        records=(RecordInput(InMemoryRecordSource(source), display=display),), options=typed_options,
    ))
    public = getattr(gbdraw, f'draw_{mode}')(source, options=options, record_displays=[display])
    assert public._drawing.tostring() == typed.drawing.tostring()
    if start is None:
        assert public._drawing.tostring() == getattr(gbdraw, f'draw_{mode}')(source)._drawing.tostring()


@pytest.mark.parametrize('mode', ['circular', 'linear'])
@pytest.mark.parametrize('start', [None, 1, 120, 25])
@pytest.mark.parametrize('reverse', [False, True])
def test_cli_records_table_resolves_and_renders_like_typed(mode, start, reverse, tmp_path):
    import importlib
    from gbdraw.api.request_render import plan_request, build_request_plan_diagram
    from gbdraw.api.requests import LinearDiagramRequest, RecordPresentation
    cli = importlib.import_module(f'gbdraw.{mode}')
    source = record()
    path = tmp_path / 'source.gbk'
    SeqIO.write(source, path, 'genbank')
    source = SeqIO.read(path, 'genbank')
    table = tmp_path / 'records.tsv'
    table.write_text('gbk\trecord_id\treverse_complement\ttopology\tdisplay_start\n'
                     f'source.gbk\t#1\t{str(reverse).lower()}\tauto\t{start or ""}\n')
    flags = ['--records_table', str(table), '-o', str(tmp_path / 'result')]
    if mode == 'circular':
        flags += ['--multi_record_canvas']
    result = getattr(cli, f'run_{mode}_from_namespace')(cli._get_args(flags))
    request = result.canonical_request
    assert request.records[0].display == RecordDisplayOptions(None, start)
    plan = plan_request(request)
    assert plan.displays[0].detected_topology == 'circular'
    assert plan.transforms[0].source_step == (-1 if reverse else 1)
    if start is not None:
        assert plan.transforms[0].source_base_to_display_index(start) == 0
    request_type = CircularDiagramRequest if mode == 'circular' else LinearDiagramRequest
    typed = request_type(records=(RecordInput(InMemoryRecordSource(source),
        presentation=RecordPresentation(reverse_complement=reverse), display=RecordDisplayOptions(None, start)),),
        options=request.options, layout=request.layout, output=request.output)
    typed_plan = plan_request(typed)
    assert typed_plan.displays == plan.displays
    # Whole generated SVG equality after parsing removes only serialization differences.
    import xml.etree.ElementTree as ET
    actual = ET.parse(result.outputs[0].svg_path).getroot()
    expected = ET.fromstring(build_request_plan_diagram(typed_plan).drawing.tostring())
    assert ET.tostring(actual) == ET.tostring(expected)


@pytest.mark.parametrize('mode', ['circular', 'linear'])
def test_cli_direct_flags_validate_cardinality_and_table_conflict(mode, tmp_path):
    import importlib
    cli = importlib.import_module(f'gbdraw.{mode}')
    path = tmp_path / 'source.gbk'
    SeqIO.write(record(), path, 'genbank')
    flags = ['--gbk', str(path), '--record_topology', 'circular', '--display_start_coordinate', '25', '-o', str(tmp_path / 'direct')]
    result = getattr(cli, f'run_{mode}_from_namespace')(cli._get_args(flags))
    assert result.canonical_request.records[0].display == RecordDisplayOptions(True, 25)
    SeqIO.write([record(), record()], path, 'genbank')
    with pytest.raises(ValidationError, match='--records_table'):
        getattr(cli, f'run_{mode}_from_namespace')(cli._get_args(flags))
    with pytest.raises(SystemExit):
        cli._get_args(['--records_table', 'records.tsv', '--record_topology', 'auto'])


@pytest.mark.parametrize('topology,start,region,column', [
    ('wrong', '', '', 'topology'), ('auto', '0', '', 'display_start'),
    ('auto', '1.5', '', 'display_start'), ('linear', '1', '', 'display_start'),
    ('circular', '1', '1-20', 'display_start'),
])
def test_table_reports_invalid_cell(tmp_path, topology, start, region, column):
    path = tmp_path / 'records.tsv'
    path.write_text(f'gbk\ttopology\tdisplay_start\tregion\nsource.gbk\t{topology}\t{start}\t{region}\n')
    with pytest.raises(ValidationError, match=rf'row 2.*{column}'):
        record_input_manifest_from_table(str(path))


@pytest.mark.parametrize('mode', ['circular', 'linear'])
def test_api_session_save_preserves_display(mode, tmp_path):
    from gbdraw.api import save_session_document
    from gbdraw.api.requests import LinearDiagramRequest
    request_type = CircularDiagramRequest if mode == 'circular' else LinearDiagramRequest
    request = request_type(records=(RecordInput(InMemoryRecordSource(record()),
                                               display=RecordDisplayOptions(None, 1)),))
    path = tmp_path / 'session.json'
    save_session_document(path, request)
    assert load_session(path)['renderRequest']['records'][0]['display']['startCoordinate'] == 1


@pytest.mark.parametrize('mode', ['circular', 'linear'])
def test_cli_session_save_preserves_display(mode, tmp_path):
    import importlib
    cli = importlib.import_module(f'gbdraw.{mode}')
    path = tmp_path / 'source.gbk'
    SeqIO.write(record(), path, 'genbank')
    session = tmp_path / 'session.json'
    getattr(cli, f'{mode}_main')(['--gbk', str(path), '--display_start_coordinate', '25',
                                    '-o', str(tmp_path / 'diagram'), '--session_output', str(session)])
    assert load_session(session)['renderRequest']['records'][0]['display']['startCoordinate'] == 25


@pytest.mark.parametrize('mode', ['circular', 'linear'])
def test_direct_cli_command_executes_approved_flags(mode, tmp_path):
    import subprocess
    import sys
    from pathlib import Path
    path = tmp_path / 'source.gbk'
    SeqIO.write(record(), path, 'genbank')
    result = subprocess.run([sys.executable, '-m', 'gbdraw.cli', mode, '--gbk', str(path),
                             '--record_topology', 'circular', '--display_start_coordinate', '25',
                             '-o', str(tmp_path / 'command')], capture_output=True, text=True,
                            cwd=Path(__file__).resolve().parents[1])
    assert result.returncode == 0, result.stderr
    assert (tmp_path / 'command.svg').is_file()
