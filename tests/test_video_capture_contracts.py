"""Timing and edit contracts for the recorded gbdraw walkthrough and tutorial videos."""

from __future__ import annotations

import ast
import dataclasses
import re
import sys
from pathlib import Path

import pytest

CAPTURE_ROOT = Path(__file__).resolve().parents[1] / 'docs' / 'capture'
sys.path.insert(0, str(CAPTURE_ROOT))
from video import walkthrough_render  # noqa: E402
from video.walkthrough import PORTRAIT_TAB, Recorder, journeys  # noqa: E402
from video.walkthrough_render import (  # noqa: E402
    FPS, HIGHLIGHTS, LANDSCAPE, PORTRAIT, STORIES, TUTORIAL_IDS, VERTICAL, Camera, Language, Overlay, TimeMap,
    Type, _cjk_fonts, _clean_breaks, _fonts, _moments, _runs, wrap_clauses)

# The tall window without its overview inset, for the framing rules that predate it.
PORTRAIT_PLAIN = dataclasses.replace(PORTRAIT, inset=0)


def test_highlights_cut_only_at_marks_the_journey_records() -> None:
    source = (CAPTURE_ROOT / 'video' / 'walkthrough.py').read_text(encoding='utf-8')
    marks = set(re.findall(r'(?:rec\.mark\(|mark=)"([^"]+)"', source))
    if 'rec.mark(f"rule-{index + 1}")' in source:
        marks |= {f'rule-{index}' for index in range(1, 5)}
    wanted = {name for first, _, last, _ in HIGHLIGHTS for name in (first, last)}
    assert wanted <= marks, sorted(wanted - marks)


def test_walkthrough_time_map_compresses_only_marked_waits() -> None:
    events = [
        {'t': 2.0, 'kind': 'speed_start', 'factor': None, 'target': 1.0},
        {'t': 10.0, 'kind': 'speed_end'},
        {'t': 12.0, 'kind': 'speed_start', 'factor': 3.0, 'target': None},
        {'t': 15.0, 'kind': 'speed_end'},
    ]
    mapping = TimeMap.build(events, 0.0, 20.0)
    assert mapping.out(2.0) == pytest.approx(2.0)
    assert mapping.out(10.0) == pytest.approx(3.0)  # 8 s wait shown in 1 s
    assert mapping.out(15.0) == pytest.approx(6.0)  # 3 s at 3x
    assert mapping.duration == pytest.approx(11.0)
    for moment in (0.5, 4.0, 11.0, 13.5, 19.0):
        assert mapping.source(mapping.out(moment)) == pytest.approx(moment)
    assert mapping.factor(mapping.out(4.0)) == pytest.approx(8.0)


def test_walkthrough_camera_eases_and_stays_inside_the_page() -> None:
    events = [
        {'t': 0.0, 'kind': 'camera', 'cx': 960, 'cy': 540, 'zoom': 1.0, 'duration': 0},
        {'t': 1.0, 'kind': 'camera', 'cx': 10, 'cy': 10, 'zoom': 2.0, 'duration': 1.0},
    ]
    camera = Camera(events, TimeMap.build([], 0.0, 5.0), (1920, 1080), LANDSCAPE)
    assert camera.at(0.5) == pytest.approx((960, 540, 1080))
    assert camera.at(3.0) == pytest.approx((480, 270, 540))  # clamped to the top-left corner
    x, y, height = camera.at(1.5)
    assert 480 < x < 960 and 270 < y < 540 and 540 < height < 1080
    assert camera.inset_at(0.5) == camera.inset_at(3.0) == (None, 0.0)   # the wide window has no inset


def test_portrait_camera_frames_the_diagram_for_whole_page_shots() -> None:
    page = {'x': 0, 'y': 0, 'width': 1920, 'height': 1080}
    diagram = {'x': 900, 'y': 150, 'width': 800, 'height': 800}
    button = {'x': 1500, 'y': 70, 'width': 60, 'height': 20}
    events = [
        {'t': 0.0, 'kind': 'camera', 'cx': 960, 'cy': 540, 'zoom': 1.0, 'duration': 0,
         'rect': page, 'subject': diagram, 'portrait': None},
        {'t': 2.0, 'kind': 'camera', 'cx': 1530, 'cy': 80, 'zoom': 2.1, 'duration': 0,
         'rect': button, 'subject': diagram, 'portrait': None},
    ]
    camera = Camera(events, TimeMap.build([], 0.0, 5.0), (1920, 1080), PORTRAIT_PLAIN)
    x, y, height = camera.at(1.0)
    assert (x, y) == pytest.approx((1300, 540)) and height == pytest.approx(1080)
    x, y, height = camera.at(3.0)
    assert height == pytest.approx(1080 / 2.6)  # small targets stop at the maximum zoom
    assert y == pytest.approx(height / 2)       # and stay inside the page


def test_portrait_camera_keeps_partial_targets_and_the_logged_zoom() -> None:
    diagram = {'x': 700, 'y': 150, 'width': 800, 'height': 800}
    # A clicked feature and its popup: wider than the tall window, but not the whole diagram.
    popup = {'x': 1000, 'y': 300, 'width': 900, 'height': 500}
    column = {'x': 0, 'y': 120, 'width': 450, 'height': 900}
    events = [
        {'t': 0.0, 'kind': 'camera', 'cx': 1450, 'cy': 550, 'zoom': 1.55, 'duration': 0,
         'rect': popup, 'subject': diagram, 'portrait': None},
        {'t': 2.0, 'kind': 'camera', 'cx': 225, 'cy': 570, 'zoom': 2.0, 'duration': 0,
         'rect': column, 'subject': diagram, 'portrait': None},
    ]
    camera = Camera(events, TimeMap.build([], 0.0, 5.0), (1920, 1080), PORTRAIT_PLAIN)
    x, _, height = camera.at(1.0)
    assert x == pytest.approx(1920 - height * 936 / 1060 / 2)  # centered on the popup, clamped to the page
    _, _, height = camera.at(3.0)
    assert height == pytest.approx(1080 / 2.0 * 1.35)  # not the column's full height


def _window_box(camera: Camera, output: float, rect: dict) -> tuple[float, float, float, float]:
    """``rect`` (page px) in window px at ``output``: left, top, right, bottom."""

    cx, cy, height = camera.at(output)
    scale = PORTRAIT.card[3] / height
    left, top = cx - height * PORTRAIT.aspect / 2, cy - height / 2
    return ((rect['x'] - left) * scale, (rect['y'] - top) * scale,
            (rect['x'] + rect['width'] - left) * scale, (rect['y'] + rect['height'] - top) * scale)


def _overlaps(a: tuple[float, float, float, float], b: tuple[float, float, float, float]) -> bool:
    return a[0] < b[2] and a[2] > b[0] and a[1] < b[3] and a[3] > b[1]


def _inset(corner: str) -> tuple[int, int, int, int]:
    x, y, width, height = PORTRAIT.inset_box((1920, 1080), corner)
    return x, y, x + width, y + height


def test_portrait_camera_centers_targets_above_the_overview_inset() -> None:
    field = {'x': 900, 'y': 500, 'width': 300, 'height': 30}
    events = [{'t': 0.0, 'kind': 'camera', 'cx': 1050, 'cy': 515, 'zoom': 2.0, 'duration': 0,
               'rect': field, 'subject': None, 'portrait': None}]
    camera = Camera(events, TimeMap.build([], 0.0, 5.0), (1920, 1080), PORTRAIT)
    assert camera.inset_at(1.0) == ('bottom', 1.0)
    box = _window_box(camera, 1.0, field)
    top, bottom = PORTRAIT.free_band((1920, 1080), 'bottom')
    assert (box[1] + box[3]) / 2 == pytest.approx((top + bottom) / 2 * PORTRAIT.card[3])
    assert (box[0] + box[2]) / 2 == pytest.approx(PORTRAIT.card[2] / 2)
    assert not _overlaps(box, _inset('bottom'))


def test_portrait_inset_moves_to_the_top_for_a_target_at_the_page_foot() -> None:
    middle = {'x': 900, 'y': 500, 'width': 300, 'height': 30}
    button = {'x': 20, 'y': 963, 'width': 426, 'height': 68}   # Generate, at the bottom of the settings column
    events = [{'t': 0.0, 'kind': 'camera', 'cx': 0, 'cy': 0, 'zoom': 2.0, 'duration': 0,
               'rect': middle, 'subject': None, 'portrait': None},
              {'t': 2.0, 'kind': 'camera', 'cx': 0, 'cy': 0, 'zoom': 1.9, 'duration': 0.9,
               'rect': button, 'subject': None, 'portrait': None}]
    camera = Camera(events, TimeMap.build([], 0.0, 6.0), (1920, 1080), PORTRAIT)
    # The page edge keeps the button low in the window, so the inset fades to the top corner.
    corner, alpha = camera.inset_at(2.1)
    assert corner == 'bottom' and 0 < alpha < 1
    corner, alpha = camera.inset_at(2.8)
    assert corner == 'top' and 0 < alpha < 1
    assert camera.inset_at(3.5) == ('top', 1.0)
    box = _window_box(camera, 3.5, button)
    assert not _overlaps(box, _inset('top')) and _overlaps(box, _inset('bottom'))


def test_portrait_inset_hides_while_the_stage_is_a_smaller_tab() -> None:
    tab = {'x': 483.0, 'y': 0.0, 'width': 954, 'height': 1080}
    events = [{'t': 0.0, 'kind': 'camera', 'cx': 960, 'cy': 540, 'zoom': 1.0, 'duration': 0,
               'rect': {'x': 0, 'y': 0, 'width': 1920, 'height': 1080}, 'subject': None, 'portrait': None},
              {'t': 2.0, 'kind': 'stage', 'rect': tab},
              {'t': 2.0, 'kind': 'camera', 'cx': 960, 'cy': 540, 'zoom': 1.0, 'duration': 0,
               'rect': tab, 'subject': None, 'portrait': None}]
    camera = Camera(events, TimeMap.build([], 0.0, 6.0), (1920, 1080), PORTRAIT)
    assert camera.inset_at(1.0) == ('bottom', 1.0)
    assert camera.inset_at(2.25) == ('bottom', pytest.approx(0.5))   # fades out
    assert camera.inset_at(3.0) == (None, 0.0)
    # The tab fills the window: no rows are kept free for an inset.
    assert camera.at(3.0) == pytest.approx((960, 540, 1080))


def test_toasts_and_the_badge_stay_clear_of_the_inset() -> None:
    viewport = (1920, 1080)
    x, y = PORTRAIT.toast_origin(viewport, (PORTRAIT.card[2], 140), 'bottom', 1.0)
    toast = (x - PORTRAIT.card[0], y - PORTRAIT.card[1], x - PORTRAIT.card[0] + PORTRAIT.card[2],
             y - PORTRAIT.card[1] + 140)
    assert not _overlaps(toast, _inset('bottom'))
    assert PORTRAIT.toast_origin(viewport, (500, 140), 'top', 1.0)[1] == PORTRAIT.toast_origin(
        viewport, (500, 140), None, 0.0)[1]   # at the window's foot when the inset is elsewhere
    width, height = 320, 62   # wider than any badge, "fast" in either language
    x, y = PORTRAIT.badge_origin((width, height))
    badge = (x - PORTRAIT.card[0], y - PORTRAIT.card[1], x - PORTRAIT.card[0] + width, y - PORTRAIT.card[1] + height)
    assert not _overlaps(badge, _inset('top'))


def _chars(text: str) -> float:
    return float(len(text))


def test_caption_lines_break_only_after_a_clause() -> None:
    assert wrap_clauses('Set the match filters', _chars, 40) == ['Set the match filters']
    assert wrap_clauses('Ready in seconds; the first run starts Python', _chars, 30) == [
        'Ready in seconds;', 'the first run starts Python']
    assert wrap_clauses('Layout › Track Layout: Above · Separate Strands', _chars, 30) == [
        'Layout › Track Layout: Above ·', 'Separate Strands']
    # No break inside a number, a file name, or a word.
    with pytest.raises(ValueError, match='Draw a 750,000 bp ruler on each axis'):
        wrap_clauses('Draw a 750,000 bp ruler on each axis', _chars, 20)
    with pytest.raises(ValueError):
        wrap_clauses('GCF_000196095.1 and GCF_000354175.2', _chars, 20)
    # Chinese breaks after its own punctuation, never at 的 or between two characters.
    assert wrap_clauses('两个基因组，各有两条染色体', _chars, 8) == ['两个基因组，', '各有两条染色体']
    with pytest.raises(ValueError, match='之间的蛋白质相似性'):
        wrap_clauses('生物合成基因簇之间的蛋白质相似性', _chars, 10)


@pytest.fixture(scope='module')
def caption_types(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Type]:
    types = {'en': Type(_fonts(tmp_path_factory.mktemp('fonts')))}
    try:
        types['zh-Hans'] = Type(_cjk_fonts())
    except FileNotFoundError:
        pass
    return types


@pytest.mark.parametrize(('code', 'layout', 'step', 'title', 'subtitle', 'lines'), [
    ('en', LANDSCAPE, 4, 'Generate a first preview',
     'Without a comparison it is ready in seconds; the first run starts Python in the browser',
     (['Generate a first preview'],
      ['Without a comparison it is ready in seconds; the first run starts Python in the browser'])),
    ('en', PORTRAIT, 4, 'Generate a first preview',
     'Without a comparison it is ready in seconds; the first run starts Python in the browser',
     (['Generate a first preview'],
      ['Without a comparison it is ready in seconds;', 'the first run starts Python in the browser'])),
    ('en', PORTRAIT, 9, 'Show each definition as species and strain', '', None),
    ('zh-Hans', PORTRAIT, 7, '与 Gallery 会话相同，对蛋白进行分组', '',
     (['与 Gallery 会话相同，', '对蛋白进行分组'], [])),
    ('zh-Hans', PORTRAIT, 3, '氨基糖苷类生物合成基因簇之间的蛋白质相似性', '', None),
    ('zh-Hans', LANDSCAPE, 3, '氨基糖苷类生物合成基因簇之间的蛋白质相似性', '',
     (['氨基糖苷类生物合成基因簇之间的蛋白质相似性'], [])),
])
def test_captions_fit_their_window_or_fail_the_render(caption_types, code, layout, step, title, subtitle,
                                                      lines) -> None:
    if code not in caption_types:
        pytest.skip('Noto Sans SC is not installed')
    overlay = Overlay(caption_types[code], layout, Language(code))
    if lines is None:
        with pytest.raises(ValueError, match=re.escape(title)):
            overlay.caption_lines(step, title, subtitle)
    else:
        assert overlay.caption_lines(step, title, subtitle) == lines


def test_every_tutorial_has_a_portrait_take_for_phone_feeds() -> None:
    takes = journeys()
    for name in TUTORIAL_IDS:
        assert STORIES[name].outputs == (('walkthrough', LANDSCAPE, f'gbdraw-tutorial-{name}'),)
        assert STORIES[name + VERTICAL].outputs == (('walkthrough', PORTRAIT, f'gbdraw-tutorial-{name}-vertical'),)
        assert not takes[name].portrait and takes[name + VERTICAL].portrait
        assert takes[name + VERTICAL].run is takes[name].run  # the same steps, shaped for the tall window
    # The tab for the opened file has the tall window's shape at full page height.
    assert PORTRAIT_TAB[0] / PORTRAIT_TAB[1] == pytest.approx(PORTRAIT.aspect, abs=0.002)


class _StubPage:
    """Just enough of a Playwright page for the recorder's pointer and camera."""

    def __init__(self, width: int, height: int) -> None:
        self.viewport_size = {'width': width, 'height': height}
        self.moves: list[tuple[float, float]] = []
        self.mouse = type('Mouse', (), {'move': lambda _, x, y: self.moves.append((x, y))})()

    def wait_for_timeout(self, _: float) -> None:
        pass

    def evaluate(self, _: str) -> None:
        return None


def test_recorder_logs_a_smaller_tab_on_the_stage(tmp_path: Path) -> None:
    rec = Recorder(_StubPage(*PORTRAIT_TAB), tmp_path, portrait=True)
    rec._origin = 0.0
    rec.origin = (483.0, 0.0)
    rec.x, rec.y = 583.0, 100.0
    rec.move_to(300, 200)
    assert rec.page.moves[-1] == pytest.approx((300, 200))   # the tab's own coordinates
    assert (rec.events[-1]['x'], rec.events[-1]['y']) == pytest.approx((783, 200))
    rec.camera(duration=0)
    assert rec.events[-1]['rect'] == {'x': 483.0, 'y': 0.0, 'width': PORTRAIT_TAB[0], 'height': PORTRAIT_TAB[1]}


def test_intro_titles_break_chinese_only_where_no_word_splits() -> None:
    title = 'Hepatoplasmataceae 基因组的共线性蛋白匹配区块'
    assert _clean_breaks(title, ['Hepatoplasmataceae 基因组的', '共线性蛋白匹配区块'])
    assert not _clean_breaks(title, ['Hepatoplasmataceae 基', '因组的共线性蛋白匹配区块'])
    tobacco = 'Nicotiana tabacum（烟草）叶绿体基因组区域'
    assert _clean_breaks(tobacco, ['Nicotiana tabacum', '（烟草）叶绿体基因组区域'])


def test_highlight_moments_skim_waits_but_not_fast_forwards() -> None:
    events = [
        {'t': 1.0, 'kind': 'speed_start', 'factor': None, 'target': 1.0},
        {'t': 5.0, 'kind': 'speed_end'},
        {'t': 6.0, 'kind': 'speed_start', 'factor': 3.0, 'target': None},
        {'t': 9.0, 'kind': 'speed_end'},
    ]
    mapping = TimeMap.build(events, 0.0, 12.0)
    assert [kind for *_, kind in mapping.segments] == ['wait', 'fast']
    moments = _moments(mapping, 0.0, mapping.duration, 1.0, 2.0)
    # 1 s + 1 s wait at 2x + 1 s + 1 s fast-forward + 3 s at normal pace.
    assert len(moments) / FPS == pytest.approx(0.5 + 6.0, abs=0.05)
    assert moments == tuple(sorted(moments))


def test_every_journey_has_a_story_and_the_reverse() -> None:
    assert set(journeys()) == set(STORIES)


def test_the_video_workflow_builds_every_journey() -> None:
    workflow = (CAPTURE_ROOT.parents[1] / '.github' / 'workflows' / 'walkthrough-video.yml').read_text(encoding='utf-8')
    matrix = re.search(r'journey: \[([^\]]*)\]', workflow).group(1)
    assert sorted(name.strip() for name in matrix.split(',')) == sorted(STORIES)


def test_intro_titles_split_italic_markup_into_runs() -> None:
    assert _runs('<i>Nicotiana tabacum</i> chloroplast genome regions') == [
        ('Nicotiana tabacum', True), (' chloroplast genome regions', False)]
    assert _runs('Aminoglycoside clusters from <i>Streptomyces</i> spp.') == [
        ('Aminoglycoside clusters from ', False), ('Streptomyces', True), (' spp.', False)]
    assert _runs('Human mitochondrial genome (AT skew)') == [('Human mitochondrial genome (AT skew)', False)]


def _screen_text() -> set[str]:
    """Every English string the recorder or renderer can put on screen."""

    strings: set[str] = {'NC_012920.1.svg downloaded (89 KB)', 'fast'}
    strings |= {source.name for journey in journeys().values() for source in journey.sources}
    for module in ('walkthrough.py', 'tutorials.py'):
        tree = ast.parse((CAPTURE_ROOT / 'video' / module).read_text(encoding='utf-8'))
        for node in ast.walk(tree):
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) and node.func.attr == 'caption':
                strings |= {arg.value for arg in node.args
                            if isinstance(arg, ast.Constant) and isinstance(arg.value, str)}
            if isinstance(node, ast.Assign) and any(getattr(t, 'id', '') == 'GALLERY' for t in node.targets):
                strings |= {row.elts[1].value for row in node.value.elts}
    render = walkthrough_render
    strings |= {text for text in (*render.FINALE_TEXT, *render.MONTAGE_CAPTION[1:])}
    strings |= {text for story in STORIES.values() for text in story.intro}
    strings |= {row[0] for row in render.Composer.OUTRO_ROWS}
    return {text for text in strings if text}


def test_every_screen_text_has_a_chinese_version() -> None:
    lang = Language('zh-Hans')
    missing = sorted(text for text in _screen_text() if _translate(lang, text) is None)
    assert not missing, missing
    unused = sorted(set(lang.strings) - _screen_text())
    assert not unused, unused


def _translate(lang: Language, text: str) -> str | None:
    try:
        return lang(text)
    except KeyError:
        return None


def test_collect_final_copies_every_file_and_never_overwrites(tmp_path: Path) -> None:
    import build_video

    final = tmp_path / 'run' / 'final'
    final.mkdir(parents=True)
    (final / 'video.mp4').write_bytes(b'mp4')
    (final / 'video.en.srt').write_text('1\n', encoding='utf-8')
    target = tmp_path / 'collected' / 'HmmtDNA_ATskew'
    copied = build_video.collect_final(final, target)
    assert [path.name for path in copied] == ['video.en.srt', 'video.mp4']
    assert (target / 'video.mp4').read_bytes() == b'mp4'
    with pytest.raises(FileExistsError):
        build_video.collect_final(final, target)


def test_build_all_refuses_used_directories_before_building(tmp_path: Path) -> None:
    import build_video

    collect = tmp_path / 'videos'
    (collect / 'tobacco-chloroplast').mkdir(parents=True)
    (collect / 'tobacco-chloroplast' / 'poster.png').write_bytes(b'png')
    with pytest.raises(FileExistsError, match='tobacco-chloroplast'):
        build_video.build_all(['tobacco-chloroplast'], tmp_path / 'out', collect=collect)
    assert not (tmp_path / 'out').exists()
