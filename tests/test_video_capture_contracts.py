"""Timing and edit contracts for the recorded gbdraw walkthrough and tutorial videos."""

from __future__ import annotations

import ast
import re
import sys
from pathlib import Path

import pytest

CAPTURE_ROOT = Path(__file__).resolve().parents[1] / 'docs' / 'capture'
sys.path.insert(0, str(CAPTURE_ROOT))
from video import walkthrough_render  # noqa: E402
from video.walkthrough import PORTRAIT_TAB, Recorder, journeys  # noqa: E402
from video.walkthrough_render import (  # noqa: E402
    FPS, HIGHLIGHTS, LANDSCAPE, PORTRAIT, STORIES, TUTORIAL_IDS, TUTORIAL_LANGUAGES, Camera, Language, TimeMap,
    _clean_breaks, _moments, _outputs, _runs, take_name)


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
    camera = Camera(events, TimeMap.build([], 0.0, 5.0), (1920, 1080), 16 / 9, False)
    assert camera.at(0.5) == pytest.approx((960, 540, 1080))
    assert camera.at(3.0) == pytest.approx((480, 270, 540))  # clamped to the top-left corner
    x, y, height = camera.at(1.5)
    assert 480 < x < 960 and 270 < y < 540 and 540 < height < 1080


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
    camera = Camera(events, TimeMap.build([], 0.0, 5.0), (1920, 1080), 936 / 1060, True)
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
    camera = Camera(events, TimeMap.build([], 0.0, 5.0), (1920, 1080), 936 / 1060, True)
    x, _, height = camera.at(1.0)
    assert x == pytest.approx(1920 - height * 936 / 1060 / 2)  # centered on the popup, clamped to the page
    _, _, height = camera.at(3.0)
    assert height == pytest.approx(1080 / 2.0 * 1.35)  # not the column's full height


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
