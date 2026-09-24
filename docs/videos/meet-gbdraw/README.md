# Meet gbdraw video

This directory owns the 38-second, silent Web introduction to gbdraw. The English captions and exact 13-scene timeline are in [`storyboard.json`](storyboard.json). A build captures nine source assets from the local Web app, verifies their biological content and file hashes, then renders the video. The capture and editor are separate: caption changes can be rendered from an existing verified asset bundle without starting Chromium, the app server, or LOSAT.

## Requirements

- Python with the repository's development dependencies and Python Playwright 1.61.0
- Playwright Chromium 149.0.7827.55
- FFmpeg/ffprobe with libx264, libass, libwebp, and SSIM filters
- `/usr/share/fonts/truetype/lato/Lato-Bold.ttf` (OFL 1.1)

The measured local versions and font/wheel hashes are in [`environment.json`](environment.json). Browser captures use the repository's offline, same-origin server and generated browser wheel.

## Build and check

From the repository root:

```bash
python tools/prepare_browser_wheel.py
python docs/capture/build_video.py build \
  --storyboard docs/videos/meet-gbdraw/storyboard.json \
  --out build/videos/meet-gbdraw/run-001
python docs/capture/build_video.py check --out build/videos/meet-gbdraw/run-001
```

Choose a new or empty `--out` directory for every build. The command fails on a capture, semantic, or media-contract error; it never borrows assets from an older run. Review `reports/review-frames/*.png` at full size and play the MP4 before accepting a result.

To revise captions or visual cues from the same source bundle:

```bash
python docs/capture/build_video.py render \
  --assets build/videos/meet-gbdraw/run-001/assets.json \
  --storyboard docs/videos/meet-gbdraw/storyboard.json \
  --out build/videos/meet-gbdraw/render-002
```

The source PNG, SVG, and raw recording hashes are checked before render. If the storyboard changes, use a new output directory. `render` imports no browser capture module.

To compare representative frames from two independent builds:

```bash
python docs/capture/build_video.py check \
  --out build/videos/meet-gbdraw/run-002 \
  --visual \
  --baseline build/videos/meet-gbdraw/run-001/reports/review-frames
```

The SSIM threshold is 0.98 for 37 deterministic review frames: each scene midpoint and both sides of static scene cuts. The first and last frames of the live SVG-export recording are reported as excluded because browser recording starts and stops at variable UI frames; its midpoint is compared. `check` also fully decodes the MP4 and five independent chapter clips and verifies stream format, dimensions, frame counts, poster size, animation size, asset integrity, and output hashes.

## Outputs

`final/meet-gbdraw.mp4` is 1920×1080, H.264/yuv420p, 30 fps, 1,140 frames, 38 seconds, with no audio. Captions are burnt from the same scene text as `final/meet-gbdraw.en.srt`; `final/meet-gbdraw.en.ass` records the exact Lato styling. `final/poster.png`, `final/meet-gbdraw.webp`, and `final/chapters/*.mp4` are generated from the same timeline. Large run files live under gitignored `build/` and should not be added to the source repository.

`assets.json` binds each source to a SHA-256 digest, its downloaded SVG or browser recording, semantic evidence, source revision, capture code and wheel hashes, and environment. P0–P3 evidence under `evidence/human-edits/` proves all 13 CDS labels, the 7/3/2/1 functional color groups, the origin-spanning D-loop annotation, and P3→P2→P3 restoration. The BGC capture uses the exact alignment anchors already recorded in the approved Gallery session.

The video source shows these fixed examples: human mitochondrial genome, annotated tobacco plastome, Lambda–DE3 comparison, and five BGC records. The editing sequence changes only the human mitochondrial presentation: CDS product names to gene symbols, CDS colors by respiratory complex, then a D-loop bracket imported through Region Annotations and placed in an inner Custom Track Slot. The export clip is a real SVG download from the final P3 result.

## Editing the story and source

`storyboard.json` is the caption, scene order, duration, and spotlight source. The 13 scene durations must add to exactly 1,140 frames. After a caption or cue change, use `render` against an existing verified `assets.json`; after a diagram or Web operation change, use `build` in a fresh directory. The three `spotlight` values in scenes 7, 9, and 11 are `labels`, `legend`, and `dloop`. Their inset rounded-rectangle bounds and colors are in `docs/capture/video/render.py` (`_spotlight`). They are drawn over the video frame, leaving the underlying exported SVG untouched.

The D-loop is an actual Region Annotations import. `docs/capture/video/human_edits.py` creates `mitochondrial_regions.tsv` with record `NC_012920.1`, source coordinates `16024` to `576`, `wraps_origin=true`, mark `bracket`, and label `D-loop`. The GUI flow imports that TSV and assigns the `mitochondrial_regions` track to an inner Custom Track Slot. The final SVG has one annotation group containing two circular arc paths across the origin. The build saves both the TSV and the exported SVG under `evidence/human-edits/` and `raw/`.

## Troubleshooting

- If the browser wheel is absent or stale, rerun `python tools/prepare_browser_wheel.py` before `build`. Capture uses the local app and fails instead of substituting old figures.
- If Chromium fails to launch in a restricted Linux sandbox, allow its local browser process and rerun the same command. Python Playwright and its installed Chromium are required.
- If `render` reports a hash mismatch, do not edit the captured PNG/SVG/WebM in place. Capture a fresh source bundle or restore the verified file.
- If `check --visual` fails, open the named current and baseline frames at full size, then inspect the changed diagram or caption. Use another independently built and reviewed `reports/review-frames` directory as the baseline.
- Every `build` and `render` needs a new output directory; this keeps evidence from separate runs distinct. MP4 and WebP files under `build/` are gitignored.

## Hands-on walkthrough

One recording of the real Web app produces two videos. The walkthrough (about 2 minutes) shows the app used at a human pace: load `HmmtDNA.gbk`, generate the map, turn on labels, switch CDS labels to gene symbols, add four functional color rules, mark the origin-spanning D-loop, inspect a feature, and export the SVG. The highlights (about 45 seconds) open with six Gallery figures and then cut the same recording down to the key moments at 1.4× speed. Both are rebuilt from the current code with one command:

```bash
python tools/prepare_browser_wheel.py
python docs/capture/build_video.py walkthrough --out build/videos/walkthrough/run-001
```

The command records the journey, renders `final/gbdraw-walkthrough.mp4` and `final/gbdraw-highlights.mp4` (1920×1080, 30 fps, H.264, silent), a subtitle file for each, and `final/poster.png`, and then decodes and checks both MP4s. It fails if the exported SVG lacks the 13 gene symbols, the four functional colors and legend entries, or the D-loop. The **Walkthrough video** GitHub Actions workflow runs the same command for every published release and on manual dispatch, and uploads the outputs as a workflow artifact.

How it works:

- `docs/capture/video/walkthrough.py` drives the local, network-isolated app with an eased pointer, typed input, and wheel scrolling. Chromium runs with `--force-device-scale-factor=2`, so its screencast delivers 3840×2160 frames. The settings column is widened to 450 px before recording so short form fields show their values. The pointer path, clicks, camera cues, captions, toasts, waits, and named marks are saved in `raw/walkthrough/events.json`. The downloaded SVG is then re-rendered at increasing magnification for the closing vector zoom, and the six Gallery source SVGs in `gbdraw/web/gallery/sources/` are rendered for the highlights montage.
- `docs/capture/video/walkthrough_render.py` composes the saved bundle without a browser. It places the app in a window over the gbdraw background, eases the camera between logged targets, draws the pointer and click ripples, and writes captions below the window with the app's vendored Inter font. Generation waits and the three repeated color rules are fast-forwarded. A `▶▶` badge marks every fast-forwarded span, and `reports/walkthrough-report.json` lists them.
- The highlights are defined in `HIGHLIGHTS` in `walkthrough_render.py` as spans between named marks (for example `labels` to `crowded`), so a new recording with different timings yields the same story. Generation waits play 2.5× faster again in the highlights, and the badge appears only at 2× or more.
- The recording hides only the floating feature-search palette, which would otherwise cover the result. Headless Chromium does not paint native `<select>` menus or color pickers, so the recorder logs each select's real option list and each chosen color, and the editor draws the menu or color popover while the pointer picks it.

To change captions, camera framing, or pacing without recording again, edit `walkthrough_render.py` (or the cues in `events.json` of a copy) and run:

```bash
python docs/capture/build_video.py walkthrough-render \
  --recording build/videos/walkthrough/run-001 --out build/videos/walkthrough/render-002
```

Changing the journey itself, including caption text, means editing `_journey` in `walkthrough.py` and running `walkthrough` again in a new directory. The rendering font is converted from `gbdraw/web/vendor/fonts/inter/*.woff2`, which requires the `brotli` Python package.
