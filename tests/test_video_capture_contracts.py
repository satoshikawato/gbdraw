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
    FPS, HIGHLIGHTS, LANDSCAPE, PORTRAIT, STORIES, TUTORIAL_IDS, TUTORIAL_LANGUAGES, Camera, Language,
    Overlay, TimeMap, Type, _cjk_fonts, _clean_breaks, _fonts, _moments, _outputs, _runs, take_name, wrap_clauses)

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
    # The landscape band below the window holds one title line and two subtitle lines.
    ('en', LANDSCAPE, 5, 'Three subtitle lines',
     'The first clause of this subtitle runs long enough to fill most of one line; '
     'the second clause of this subtitle runs long enough to fill another line; '
     'and the third clause of this subtitle would fall below the caption band', None),
    ('en', LANDSCAPE, 5, 'A title clause long enough for one line of the band, '
     'and a second clause that wraps the title', 'Any subtitle', None),
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
        for language in TUTORIAL_LANGUAGES[name]:
            wide, tall = take_name(name, language=language), take_name(name, vertical=True, language=language)
            assert STORIES[wide].outputs == (('walkthrough', LANDSCAPE, f'gbdraw-tutorial-{name}'),)
            assert STORIES[tall].outputs == (('walkthrough', PORTRAIT, f'gbdraw-tutorial-{name}-vertical'),)
            assert STORIES[wide].language == STORIES[tall].language == language
            assert not takes[wide].portrait and takes[tall].portrait
            # The same steps, shaped for the tall window or recorded in another UI language.
            assert takes[wide].run is takes[tall].run is takes[name].run
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

    def bring_to_front(self) -> None:
        pass

    def evaluate(self, _: str) -> None:
        return None


def test_a_take_renders_only_the_language_it_recorded() -> None:
    assert take_name('BGC0000708-BGC0000713', vertical=True, language='zh-Hans') \
        == 'BGC0000708-BGC0000713-vertical.zh-Hans'
    assert 'zh-Hans' in TUTORIAL_LANGUAGES['vibrio-harveyi-group-collinear']
    for story in STORIES.values():
        outputs = list(_outputs(story))
        assert {lang.code for lang, *_ in outputs} == {story.language}
        suffix = '' if story.language == 'en' else f'.{story.language}'
        assert [stem for *_, stem in outputs] == [stem + suffix for *_, stem in story.outputs]


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


def test_switching_to_a_tab_logs_the_stage_it_fills(tmp_path: Path) -> None:
    rec = Recorder(_StubPage(1920, 1080), tmp_path, portrait=True)
    rec._origin = 0.0
    rec._cdp = type('Cdp', (), {'send': lambda *_: None})()
    rec._screencast = lambda: None
    rec._first_frame = lambda _: None
    rec.switch_page(_StubPage(*PORTRAIT_TAB), origin=(483.0, 0.0))
    stages = [event for event in rec.events if event['kind'] == 'stage']
    assert [event['rect'] for event in stages] == [
        {'x': 483.0, 'y': 0.0, 'width': PORTRAIT_TAB[0], 'height': PORTRAIT_TAB[1]}]


def test_portrait_camera_stays_inside_a_smaller_tab() -> None:
    tab = {'x': 483.0, 'y': 0.0, 'width': 954, 'height': 1080}
    popup = {'x': 1180.0, 'y': 300.0, 'width': 240, 'height': 300}   # at the tab's right edge
    events = [{'t': 0.0, 'kind': 'stage', 'rect': tab},
              {'t': 0.0, 'kind': 'camera', 'cx': 1300, 'cy': 450, 'zoom': 2.0, 'duration': 0,
               'rect': popup, 'subject': None, 'portrait': None}]
    camera = Camera(events, TimeMap.build([], 0.0, 4.0), (1920, 1080), PORTRAIT)
    cx, cy, height = camera.at(2.0)
    width = height * PORTRAIT.aspect
    assert tab['x'] <= cx - width / 2 and cx + width / 2 <= tab['x'] + tab['width'] + 1e-6
    assert 0 <= cy - height / 2 and cy + height / 2 <= tab['height'] + 1e-6


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
    for module in ('walkthrough.py', 'tutorials.py', 'promo.py'):
        tree = ast.parse((CAPTURE_ROOT / 'video' / module).read_text(encoding='utf-8'))
        for node in ast.walk(tree):
            # Captions come from ``rec.caption`` and from helpers that caption a view, such as ``_show_figure``.
            name = getattr(node.func, 'attr', None) or getattr(node.func, 'id', None) if isinstance(node, ast.Call) else None
            if name in ('caption', '_show_figure'):
                strings |= {arg.value for arg in node.args
                            if isinstance(arg, ast.Constant) and isinstance(arg.value, str)}
            if isinstance(node, ast.Assign) and any(getattr(t, 'id', '') == 'GALLERY' for t in node.targets):
                strings |= {row.elts[1].value for row in node.value.elts}
    render = walkthrough_render
    strings |= {text for text in (*render.FINALE_TEXT, *render.MONTAGE_CAPTION[1:])}
    strings |= {text for story in STORIES.values() for text in story.intro}
    strings |= {text for story in STORIES.values() for text in story.outro}
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


class _PreviewPage(_StubPage):
    """The Result Preview as the detail holds see it: a drawing scaled from the canvas's top center.

    "Zoom in" and "Zoom out" change the zoom by 0.1 like the app's buttons;
    a press that moves before its release pans the drawing by the same amount.
    """

    CANVAS = {'x': 461.0, 'y': 211.0, 'width': 1442.0, 'height': 737.0}
    DRAWING = (2680.0, 1020.0)          # the drawing at zoom 1, in CSS px
    ROWS = (0.0, 0.0, 1.0, 0.8)         # the detail part, as shares of the drawing
    LABEL = 24.0                        # the record label font at zoom 1
    BUTTONS = {'Zoom in': {'x': 1500.0, 'y': 960.0, 'width': 40.0, 'height': 40.0},
               'Zoom out': {'x': 1590.0, 'y': 960.0, 'width': 40.0, 'height': 40.0}}

    def __init__(self) -> None:
        super().__init__(1920, 1080)
        self.zoom, self.pan, self.clicks = 0.2, [0.0, 0.0], []
        self._at = (0.0, 0.0)
        self._down = None

        def move(x: float, y: float) -> None:
            self._at = (x, y)
        self.mouse = type('Mouse', (), {'move': lambda _, x, y: move(x, y),
                                        'down': lambda _: setattr(self, '_down', self._at),
                                        'up': lambda _: self._release()})()

    def _release(self) -> None:
        (sx, sy), (ex, ey) = self._down, self._at
        if (sx, sy) != (ex, ey):
            self.pan[0] += ex - sx
            self.pan[1] += ey - sy
            return
        for name, box in self.BUTTONS.items():
            if box['x'] <= ex <= box['x'] + box['width'] and box['y'] <= ey <= box['y'] + box['height']:
                self.clicks.append(name)
                self.zoom = round(self.zoom + (0.1 if name == 'Zoom in' else -0.1), 1)

    def svg(self) -> dict:
        canvas, (width, height) = self.CANVAS, self.DRAWING
        return {'x': canvas['x'] + canvas['width'] / 2 - width * self.zoom / 2 + self.pan[0],
                'y': canvas['y'] + 8 + self.pan[1], 'width': width * self.zoom, 'height': height * self.zoom}

    def rows(self) -> dict:
        svg, (x, y, width, height) = self.svg(), self.ROWS
        return {'x': svg['x'] + x * svg['width'], 'y': svg['y'] + y * svg['height'],
                'width': width * svg['width'], 'height': height * svg['height']}

    def get_by_role(self, _: str, *, name: str, exact: bool) -> object:
        box = self.BUTTONS[name]
        return type('Button', (), {'bounding_box': lambda _: dict(box)})()

    def evaluate(self, script: str, arg: object = None) -> object:
        from video import tutorials
        from video.walkthrough import _FEATURE_BOUNDS, _PREVIEW
        canvas = self.CANVAS
        if script == _PREVIEW:
            blank = [(canvas['x'] + canvas['width'] * fx, canvas['y'] + canvas['height'] * fy)
                     for fx in (0.04, 0.96, 0.5, 0.25, 0.75) for fy in (0.5, 0.3, 0.7)]
            return {'canvas': dict(canvas), 'svg': self.svg(), 'blank': blank, 'zoom': self.zoom}
        if script == _FEATURE_BOUNDS:
            return self.svg()
        if script == tutorials._DETAIL:
            return self.rows()
        if script == getattr(tutorials, '_RECORD_LABEL', None):
            return self.LABEL * self.zoom
        if 'window.__GBDRAW_APP__.t(' in script:
            return arg[0]
        return None


def _detail_hold(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, *, portrait: bool) -> tuple[Recorder, _PreviewPage]:
    from video import tutorials, walkthrough

    expect = type('Expect', (), {'to_be_visible': lambda _: None})()
    monkeypatch.setattr(walkthrough, 'expect', lambda _: expect)
    page = _PreviewPage()
    rec = Recorder(page, tmp_path, portrait=portrait)
    rec._origin = 0.0
    tutorials._show_figure(rec, 4, 'Five clusters, one row each', detail=tutorials._preview_part('[data-x]'))
    return rec, page


def test_portrait_detail_holds_zoom_the_preview_until_labels_read_and_fit_again(
        tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    rec, page = _detail_hold(tmp_path, monkeypatch, portrait=True)
    zoom_in = page.clicks.count('Zoom in')
    assert zoom_in >= 1 and page.clicks == ['Zoom in'] * zoom_in + ['Zoom out'] * zoom_in   # the app's control, then Fit
    assert page.zoom == pytest.approx(0.2)
    kinds = [event['kind'] for event in rec.events]
    # The clicks are skimmed under the fast-forward badge.
    first_down = kinds.index('down')
    assert 'speed_start' in kinds[:first_down]
    detail = next(event for event in rec.events if event['kind'] == 'detail')
    assert detail['cap_px'] >= tutorials_cap() and detail['stops'] == 2
    stop = next(event for event in rec.events if event['kind'] == 'camera' and event['t'] >= detail['t'])['portrait']
    # The view of the first stop starts at the part's left end and lies inside the canvas.
    cx, cy, height = PORTRAIT.portrait_view((1920, 1080), stop, 'bottom')
    canvas = page.CANVAS
    assert canvas['x'] <= cx - height * PORTRAIT.aspect / 2 and cx + height * PORTRAIT.aspect / 2 <= canvas['x'] + canvas['width']
    assert canvas['y'] <= cy - height / 2 and cy + height / 2 <= canvas['y'] + canvas['height']
    # The pan to the second stop moves the drawing left, before the fit zooms out.
    after = rec.events[rec.events.index(next(e for e in rec.events if e['kind'] == 'camera' and e['t'] >= detail['t'])):]
    downs = [e for e in after if e['kind'] == 'down']
    ups = [e for e in after if e['kind'] == 'up']
    assert ups[0]['x'] < downs[0]['x']
    # The hold ends on the whole figure, fitted as before.
    whole = next(event for event in rec.events if event['kind'] == 'camera')
    assert rec.events[-1]['kind'] == 'camera'
    assert (rec.events[-1]['rect']['width'], rec.events[-1]['rect']['height']) \
        == (whole['rect']['width'], whole['rect']['height'])


def tutorials_cap() -> float:
    from video import tutorials

    return getattr(tutorials, 'DETAIL_CAP_PX', 11) - 0.5


def test_landscape_detail_holds_only_move_the_camera(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    rec, page = _detail_hold(tmp_path, monkeypatch, portrait=False)
    assert page.clicks == [] and page.pan == [0.0, 0.0] and page.zoom == 0.2
    cameras = [event for event in rec.events if event['kind'] == 'camera']
    assert len(cameras) == 3                                   # whole figure, the part, whole figure
    part = page.rows()
    zoom = min(1920 / part['width'], 1080 / part['height']) / 1.08
    assert cameras[1]['rect'] == {key: round(part[key], 1) for key in ('x', 'y', 'width', 'height')}
    assert cameras[1]['zoom'] == pytest.approx(max(1.0, min(zoom, 3.0)))


def test_title_card_curl_reads_the_loader_geometry_from_the_app() -> None:
    mark = walkthrough_render.loader_mark()
    first, second = mark.pens
    # Pen A starts the stroke; pen B starts before pen A ends, as in the app's loader.
    assert first.offset == 0 and 0 < second.offset < first.path_length
    assert mark.total == pytest.approx(second.offset + second.path_length)
    assert mark.drawn(0) == 0 and mark.drawn(walkthrough_render.CURL_SECONDS) == 1
    assert 0 < mark.drawn(walkthrough_render.CURL_SECONDS / 2) < 1
    index = walkthrough_render.INDEX_HTML.read_text(encoding='utf-8')
    with pytest.raises(ValueError):
        walkthrough_render.loader_mark(index.replace('id="gbdraw-loader-template"', 'id="moved"'))
    with pytest.raises(ValueError):
        walkthrough_render.loader_mark(index.replace('animation-timing-function: cubic-bezier', 'timing: cubic-bezier'))


@pytest.mark.parametrize('width', [820, 760])
def test_the_finished_curl_is_the_static_logo_mark(width: int) -> None:
    curl = walkthrough_render.logo_curl(width)
    logo = walkthrough_render._white_logo(width)
    split = curl.split
    static = logo.getchannel('A').crop((0, 0, split, logo.height))
    drawn = curl.mark(1.0).getchannel('A')
    assert drawn.size == logo.size and drawn.crop((split, 0, logo.width, logo.height)).getbbox() is None
    drawn = drawn.crop((0, 0, split, logo.height))
    # Only the antialiased rim may differ: under 4% of the mark's coverage at the intro widths.
    difference = sum(abs(a - b) for a, b in zip(static.getdata(), drawn.getdata()))
    assert difference < 0.04 * sum(static.getdata())
    # The word keeps everything right of the mark, and nothing of the mark.
    word = curl.word.getchannel('A')
    assert word.crop((0, 0, split, logo.height)).getbbox() is None
    assert list(word.crop((split, 0, logo.width, logo.height)).getdata()) == \
        list(logo.getchannel('A').crop((split, 0, logo.width, logo.height)).getdata())
    assert curl.mark(0.0).getchannel('A').getbbox() is None
    half = sum(curl.mark(0.5).getchannel('A').getdata())
    assert 0.2 * sum(drawn.getdata()) < half < 0.8 * sum(drawn.getdata())
