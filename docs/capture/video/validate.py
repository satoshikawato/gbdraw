"""Read back finished media and compare representative frames with FFmpeg SSIM."""

from __future__ import annotations

import json
import re
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

from PIL import Image

from video.model import load_assets, sha256
from video.render import CHAPTERS, _probe_video


def _assert_video(path: Path, frames: int) -> dict:
    probe = _probe_video(path)
    streams = probe.get('streams', [])
    if len(streams) != 1 or streams[0].get('codec_type') != 'video' or \
            streams[0].get('codec_name') != 'h264' or streams[0].get('width') != 1920 or \
            streams[0].get('height') != 1080 or streams[0].get('r_frame_rate') != '30/1' or streams[0].get('avg_frame_rate') != '30/1' or \
            streams[0].get('sample_aspect_ratio') != '1:1' or \
            streams[0].get('pix_fmt') != 'yuv420p' or int(streams[0].get('nb_read_frames', -1)) != frames:
        raise AssertionError(f'Invalid {path.name}: {probe}')
    subprocess.run(['ffmpeg', '-v', 'error', '-i', str(path), '-f', 'null', '-'], check=True)
    return probe


def _check_human_evidence(source_root: Path) -> None:
    evidence = json.loads(
        (source_root / 'evidence/human-edits/transitions.json').read_text(encoding='utf-8')
    )
    source = Path(__file__).resolve().parents[3] / 'gbdraw/web/tutorial-data/human-mitochondrion/HmmtDNA.gbk'
    if sha256(source) != evidence['source_sha256']:
        raise AssertionError('Human mitochondrial source changed since capture')
    states = evidence['states']
    names = ('human.labels-product', 'human.labels-gene',
             'human.functional-colors', 'human.dloop-bracket')
    identifiers = [
        {row['id'] for row in states[name]['semantics']['features']}
        for name in names
    ]
    if any(len(rows) != 37 or rows != identifiers[0] for rows in identifiers):
        raise AssertionError('Mitochondrial source features changed across P0–P3')
    if states['roundtrip']['p3_svg_sha256'] != states[names[-1]]['svg_sha256']:
        raise AssertionError('P3 did not survive the annotation-slot round trip')
    regions = source_root / 'evidence/human-edits/tables/mitochondrial_regions.tsv'
    if sha256(regions) != evidence['region_table_sha256'] or \
            '\td_loop\tbracket\tNC_012920.1\t16024\t576\tsource\ttrue\tD-loop\t' not in regions.read_text():
        raise AssertionError('The origin-spanning D-loop table changed')
    svg = source_root / 'figures/human.dloop-bracket.svg'
    root = ET.parse(svg).getroot()
    annotations = [element for element in root.iter()
                   if element.get('data-gbdraw-annotation-id') == 'd_loop']
    if len(annotations) != 1:
        raise AssertionError('Expected one D-loop annotation group')
    annotation = annotations[0]
    expected = {
        'data-gbdraw-annotation-mark': 'bracket',
        'data-gbdraw-annotation-set-id': 'mitochondrial_regions',
        'data-gbdraw-annotation-track-id': 'mitochondrial_regions',
        'data-gbdraw-record-id': 'NC_012920.1',
    }
    if any(annotation.get(key) != value for key, value in expected.items()):
        raise AssertionError('D-loop bracket has the wrong source or track')
    arcs = [element for element in annotation.iter()
            if element.tag.rsplit('}', 1)[-1] == 'path' and ' A ' in element.get('d', '')]
    if len(arcs) != 2:
        raise AssertionError('Origin-spanning D-loop must render as two arcs in one annotation')
    final_svg = source_root / 'raw/final-export.svg'
    if sha256(final_svg) != states[names[-1]]['svg_sha256']:
        raise AssertionError('Recorded final export differs from P3')
    legend = {row['key']: row['fill'] for row in states['human.functional-colors']['semantics']['legend']}
    expected_legend = {
        'NADH dehydrogenase': '#3b82f6',
        'Cytochrome c oxidase': '#ef4444',
        'ATP synthase': '#f59e0b',
        'Cytochrome b': '#8b5cf6',
    }
    if any(legend.get(key) != color for key, color in expected_legend.items()):
        raise AssertionError('Functional color legend differs from the four CDS rules')


def check_video(out: Path, *, visual: bool = False, baseline: Path | None = None) -> dict:
    report_path = out / 'reports' / 'build-report.json'
    report = json.loads(report_path.read_text(encoding='utf-8'))
    assets_path = (out / report['asset_manifest']).resolve()
    if sha256(assets_path) != report['asset_manifest_sha256']:
        raise AssertionError('Captured asset manifest changed after editing')
    load_assets(assets_path)
    _check_human_evidence(assets_path.parent)
    final = out / 'final'
    video = final / 'meet-gbdraw.mp4'
    if sha256(video) != report['video_sha256']:
        raise AssertionError('Final video checksum changed')
    probe = _assert_video(video, 1140)
    for chapter, _, frames in CHAPTERS:
        path = final / 'chapters' / chapter
        if sha256(path) != report['chapters'][chapter]:
            raise AssertionError(f'Chapter checksum changed: {chapter}')
        _assert_video(path, frames)
    if not (final / 'meet-gbdraw.en.srt').is_file() or not (final / 'meet-gbdraw.en.ass').is_file():
        raise AssertionError('Caption files are missing')
    with Image.open(final / 'poster.png') as poster:
        if poster.size != (1920, 1080):
            raise AssertionError('Poster resolution changed')
    webp = final / 'meet-gbdraw.webp'
    with Image.open(webp) as animation:
        if animation.size != (960, 540) or not animation.is_animated or animation.info.get('loop') != 0:
            raise AssertionError('README WebP is not animated at 960×540')
    if webp.stat().st_size > 2_000_000:
        raise AssertionError(f'README WebP is too large: {webp.stat().st_size}')
    result = {'status': 'PASS', 'frames': 1140, 'seconds': 38,
              'video_sha256': report['video_sha256'], 'webp_bytes': webp.stat().st_size,
              'source_bundle_verified': True, 'probe': probe}
    if visual:
        if baseline is None or not baseline.is_dir():
            raise ValueError('A reviewed reference-frame directory is required for SSIM')
        scores = {}
        excluded = {
            'export-svg-first.png': 'Live browser capture begins at a non-deterministic UI frame',
            'export-svg-last.png': 'Live browser capture ends at a non-deterministic UI frame',
        }
        for current in sorted((out / 'reports' / 'review-frames').glob('*.png')):
            if current.name in excluded:
                continue
            reference = baseline / current.name
            if not reference.is_file():
                raise FileNotFoundError(reference)
            completed = subprocess.run(
                ['ffmpeg', '-hide_banner', '-i', str(reference), '-i', str(current),
                 '-lavfi', 'ssim', '-f', 'null', '-'], text=True, capture_output=True, check=True,
            )
            matches = re.findall(r'All:([0-9.]+)', completed.stderr)
            if not matches:
                raise AssertionError(f'FFmpeg did not report SSIM for {current.name}')
            scores[current.name] = min(map(float, matches))
        if not scores or min(scores.values()) < 0.98:
            raise AssertionError(f'Visual regression: minimum SSIM {min(scores.values()) if scores else None}')
        result['ssim'] = scores
        result['ssim_excluded'] = excluded
        (out / 'reports' / 'ssim-report.json').write_text(
            json.dumps({'scores': scores, 'excluded': excluded}, indent=2), encoding='utf-8'
        )
    (out / 'reports' / 'validation.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    return result
