"""Edit verified figure assets into a fixed 38-second, silent video."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
from pathlib import Path

from PIL import Image, ImageChops, ImageDraw, ImageFont, ImageOps

from video.model import load_assets, load_storyboard, sha256

FONT = Path('/usr/share/fonts/truetype/lato/Lato-Bold.ttf')
CHAPTERS = (
    ('01_introduction.mp4', 0, 60),
    ('02_single_genome.mp4', 60, 180),
    ('03_synteny_comparison.mp4', 240, 300),
    ('04_biological_styling.mp4', 540, 360),
    ('05_export_and_web.mp4', 900, 240),
)


def _run(*args: str, cwd: Path | None = None, capture: bool = False) -> subprocess.CompletedProcess:
    return subprocess.run(list(args), cwd=cwd, check=True, text=True,
                          capture_output=capture)


def _timestamp(frame: int, *, srt: bool) -> str:
    seconds, remainder = divmod(frame, 30)
    hours, seconds = divmod(seconds, 3600)
    minutes, seconds = divmod(seconds, 60)
    if srt:
        return f'{hours:02}:{minutes:02}:{seconds:02},{round(remainder * 1000 / 30):03}'
    return f'{hours}:{minutes:02}:{seconds:02}.{round(remainder * 100 / 30):02}'


def _subtitles(scenes: list[dict], final: Path) -> None:
    srt = []
    events = []
    for number, scene in enumerate(scenes, 1):
        start, end = scene['start_frame'], scene['start_frame'] + scene['frames']
        caption = scene['caption']
        srt.append(f'{number}\n{_timestamp(start, srt=True)} --> {_timestamp(end, srt=True)}\n{caption}\n')
        escaped = caption.replace('\\', r'\\').replace('{', r'\{').replace('}', r'\}')
        events.append(f'Dialogue: 0,{_timestamp(start, srt=False)},{_timestamp(end, srt=False)},Caption,,0,0,0,,{escaped}')
    (final / 'meet-gbdraw.en.srt').write_text('\n'.join(srt), encoding='utf-8')
    ass = '''[Script Info]
ScriptType: v4.00+
PlayResX: 1920
PlayResY: 1080
WrapStyle: 2
ScaledBorderAndShadow: yes

[V4+ Styles]
Format: Name, Fontname, Fontsize, PrimaryColour, SecondaryColour, OutlineColour, BackColour, Bold, Italic, Underline, StrikeOut, ScaleX, ScaleY, Spacing, Angle, BorderStyle, Outline, Shadow, Alignment, MarginL, MarginR, MarginV, Encoding
Style: Caption,Lato,46,&H00FFFFFF,&H00FFFFFF,&H00283344,&H80000000,-1,0,0,0,100,100,0,0,1,3,0,2,100,100,65,1

[Events]
Format: Layer, Start, End, Style, Name, MarginL, MarginR, MarginV, Effect, Text
'''
    (final / 'meet-gbdraw.en.ass').write_text(ass + '\n'.join(events) + '\n', encoding='utf-8')


def _content_crop(image: Image.Image) -> Image.Image:
    white = Image.new('RGB', image.size, 'white')
    diff = ImageChops.difference(image.convert('RGB'), white)
    bbox = diff.point(lambda value: 255 if value > 18 else 0).getbbox()
    if not bbox:
        raise AssertionError('Figure source is blank')
    pad = 20
    return image.crop((max(0, bbox[0]-pad), max(0, bbox[1]-pad),
                       min(image.width, bbox[2]+pad), min(image.height, bbox[3]+pad)))


def _intro(assets: dict, root: Path, target: Path) -> None:
    canvas = Image.new('RGB', (1920, 1080), '#f8fafc')
    draw = ImageDraw.Draw(canvas)
    font = ImageFont.truetype(FONT, 27)
    labels = ('Circular genome', 'Plastome', 'Whole genomes', 'Biosynthetic clusters')
    for index, asset_id in enumerate(('human.circular', 'tobacco.plastome',
                                       'lambda-de3.comparison', 'bgc.comparison')):
        x = 64 + (index % 2) * 912
        y = 45 + (index // 2) * 450
        draw.rounded_rectangle((x, y, x+880, y+420), radius=18, fill='white', outline='#d8e2ec', width=2)
        image = Image.open(root / assets[asset_id]['path']).convert('RGB')
        figure = ImageOps.contain(_content_crop(image), (840, 345), Image.Resampling.LANCZOS)
        canvas.paste(figure, (x+(880-figure.width)//2, y+20+(345-figure.height)//2))
        draw.text((x+26, y+376), labels[index], font=font, fill='#24364b')
    canvas.save(target)


def _spotlight(source: Path, target: Path, kind: str) -> None:
    image = Image.open(source).convert('RGBA')
    overlay = Image.new('RGBA', image.size, (0, 0, 0, 0))
    draw = ImageDraw.Draw(overlay)
    boxes = {
        'labels': (1060, 290, 1430, 605),
        'legend': (1275, 300, 1725, 760),
        'dloop': (825, 75, 1130, 260),
    }
    box = boxes[kind]
    draw.rounded_rectangle(box, radius=28, outline=(245, 158, 11, 150), width=7)
    draw.rounded_rectangle(tuple(v+(-10 if i < 2 else 10) for i,v in enumerate(box)),
                           radius=36, outline=(245, 158, 11, 45), width=8)
    Image.alpha_composite(image, overlay).convert('RGB').save(target)


def _probe_video(path: Path) -> dict:
    result = _run('ffprobe', '-v', 'error', '-count_frames', '-show_entries',
                  'stream=codec_name,width,height,r_frame_rate,avg_frame_rate,sample_aspect_ratio,pix_fmt,nb_read_frames,codec_type',
                  '-show_entries', 'format=duration', '-of', 'json', str(path), capture=True)
    return json.loads(result.stdout)


def extract_review_frames(video: Path, scenes: list[dict], review: Path) -> None:
    """Save scene midpoints and both sides of every cut by exact frame index."""

    review.mkdir(parents=True, exist_ok=True)
    selections = {}
    for scene in scenes:
        start = scene['start_frame']
        end = start + scene['frames'] - 1
        selections[start] = f"{scene['id']}-first.png"
        selections[start + scene['frames'] // 2] = f"{scene['id']}.png"
        selections[end] = f"{scene['id']}-last.png"
    frames = sorted(selections)
    expression = "select='" + '+'.join(f'eq(n,{frame})' for frame in frames) + "'"
    _run('ffmpeg', '-hide_banner', '-loglevel', 'error', '-y', '-i', str(video),
         '-vf', expression, '-vsync', '0', str(review / 'frame-%03d.png'))
    captures = sorted(review.glob('frame-*.png'))
    if len(captures) != len(frames):
        raise AssertionError(f'Expected {len(frames)} review frames, found {len(captures)}')
    for frame, capture in zip(frames, captures):
        capture.rename(review / selections[frame])


def render_video(assets_path: Path, storyboard_path: Path, out: Path) -> Path:
    storyboard = load_storyboard(storyboard_path)
    manifest = load_assets(assets_path)
    if not FONT.is_file():
        raise FileNotFoundError(f'Pinned editing font is missing: {FONT}')
    if out.exists() and any(out.iterdir()) and out.resolve() != assets_path.parent.resolve():
        raise FileExistsError(f'Render output is not empty: {out}')
    out.mkdir(parents=True, exist_ok=True)
    final = out / 'final'
    reports = out / 'reports'
    frames = out / 'frames'
    final.mkdir(exist_ok=True)
    reports.mkdir(exist_ok=True)
    frames.mkdir(exist_ok=True)
    source_root = assets_path.parent
    assets = manifest['assets']
    _subtitles(storyboard['scenes'], final)
    scene_images = out / 'scene-images'
    scene_images.mkdir(exist_ok=True)
    intro = scene_images / 'intro.png'
    _intro(assets, source_root, intro)
    poster = Image.open(intro).convert('RGB')
    poster_draw = ImageDraw.Draw(poster)
    poster_draw.text((84, 946), 'Meet gbdraw', font=ImageFont.truetype(FONT, 58), fill='#17324d')
    poster_draw.text((84, 1011), 'Publication-ready genome figures',
                     font=ImageFont.truetype(FONT, 30), fill='#40566e')
    poster.save(final / 'poster.png')

    export_frames = scene_images / 'export'
    export_frames.mkdir(exist_ok=True)
    _run('ffmpeg', '-hide_banner', '-loglevel', 'error', '-y', '-i',
         str(source_root / assets['human.svg-export']['path']), '-vsync', '0',
         '-vf', 'scale=1600:900,pad=1920:1080:160:20:white',
         '-frames:v', '120', str(export_frames / '%04d.png'))
    export_pngs = sorted(export_frames.glob('*.png'))
    if len(export_pngs) != 120:
        raise AssertionError(f'Export recording has {len(export_pngs)} decoded frames, expected 120')

    frame_number = 1
    for scene in storyboard['scenes']:
        if scene['id'] == 'export-svg':
            sources = export_pngs
        else:
            source = intro if scene['id'] == 'intro' else source_root / assets[scene['assets'][0]]['path']
            if scene.get('spotlight'):
                focused = scene_images / f"{scene['id']}.png"
                _spotlight(source, focused, scene['spotlight'])
                source = focused
            sources = [source] * scene['frames']
        for source in sources:
            (frames / f'{frame_number:04d}.png').symlink_to(source.resolve())
            frame_number += 1
    if frame_number != 1141:
        raise AssertionError(f'Timeline has {frame_number-1} frames')

    video = final / 'meet-gbdraw.mp4'
    _run('ffmpeg', '-hide_banner', '-loglevel', 'error', '-y', '-framerate', '30',
         '-i', 'frames/%04d.png', '-vf', 'ass=final/meet-gbdraw.en.ass:fontsdir=/usr/share/fonts/truetype/lato,setsar=1',
         '-frames:v', '1140', '-an', '-c:v', 'libx264', '-preset', 'veryfast', '-crf', '20',
         '-pix_fmt', 'yuv420p', '-movflags', '+faststart', 'final/meet-gbdraw.mp4', cwd=out)
    _run('ffmpeg', '-v', 'error', '-i', str(video), '-f', 'null', '-')
    probe = _probe_video(video)
    streams = probe['streams']
    if len(streams) != 1 or streams[0]['codec_name'] != 'h264' or streams[0]['pix_fmt'] != 'yuv420p' \
            or (streams[0]['width'], streams[0]['height']) != (1920,1080) \
            or streams[0]['r_frame_rate'] != '30/1' or streams[0]['avg_frame_rate'] != '30/1' \
            or streams[0].get('sample_aspect_ratio') != '1:1' or int(streams[0]['nb_read_frames']) != 1140:
        raise AssertionError(f'Final MP4 contract failed: {probe}')

    for name, start, count in CHAPTERS:
        chapter = final / 'chapters' / name
        chapter.parent.mkdir(exist_ok=True)
        _run('ffmpeg', '-hide_banner', '-loglevel', 'error', '-y', '-ss', str(start/30),
             '-i', str(video), '-frames:v', str(count), '-vf', 'setsar=1', '-an', '-c:v', 'libx264',
             '-preset', 'veryfast', '-crf', '20', '-pix_fmt', 'yuv420p', '-movflags', '+faststart', str(chapter))
    webp = final / 'meet-gbdraw.webp'
    _run('ffmpeg', '-hide_banner', '-loglevel', 'error', '-y', '-i', str(video),
         '-vf', 'fps=15,scale=960:540:flags=lanczos', '-an', '-c:v', 'libwebp_anim',
         '-quality', '95', '-loop', '0', str(webp))

    review = reports / 'review-frames'
    extract_review_frames(video, storyboard['scenes'], review)
    report = {
        'status': 'PASS', 'video': 'final/meet-gbdraw.mp4', 'video_sha256': sha256(video),
        'asset_manifest': os.path.relpath(assets_path.resolve(), out.resolve()),
        'asset_manifest_sha256': sha256(assets_path),
        'storyboard_sha256': sha256(storyboard_path), 'font_sha256': sha256(FONT),
        'probe': probe, 'webp_bytes': webp.stat().st_size,
        'chapters': {name: sha256(final / 'chapters' / name) for name,_,_ in CHAPTERS},
    }
    (reports / 'build-report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    shutil.rmtree(frames)
    shutil.rmtree(scene_images)
    return video
