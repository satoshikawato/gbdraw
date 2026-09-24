# Meet gbdraw videos

One recording of the real Web app produces three videos:

- **Walkthrough** (1920×1080, about 2 minutes): load `HmmtDNA.gbk`, generate the map, turn on labels, switch CDS labels to gene symbols, add four functional color rules, mark the origin-spanning D-loop, inspect a feature, and export the SVG, all at a human pace.
- **Highlights** (1920×1080, about 45 seconds): six Gallery figures, then the key moments of the same recording at 1.4× speed.
- **Vertical highlights** (1080×1920, same cut) for TikTok, RedNote, and other phone feeds. Captions sit above the app window and the brand below it, clear of the top 400 px and bottom 430 px where those apps draw their own controls. The camera re-frames each shot for the tall window: whole-page shots show the diagram, and close-ups follow the control being used.

All three end with the exported SVG magnified 7× and a gbdraw.app card, and each is rendered in English and in Simplified Chinese (`*.zh-Hans.mp4`).

## Build

From the repository root:

```bash
python tools/prepare_browser_wheel.py
python docs/capture/build_video.py build --out build/videos/walkthrough/run-001
```

The command records the journey, renders `final/gbdraw-walkthrough.mp4`, `final/gbdraw-highlights.mp4`, and `final/gbdraw-highlights-vertical.mp4` (30 fps, H.264), a subtitle file for each, `final/poster.png`, and two covers for phone feeds (`cover-vertical.png` at 9:16 and `cover-3x4.png` for RedNote), and then decodes and checks every MP4. The landscape videos have no audio track; the vertical one carries a silent AAC track because some social upload paths reject video-only files. It fails if the exported SVG lacks the 13 gene symbols, the four functional colors and legend entries, or the D-loop. Use a new or empty `--out` directory for every build; `build/` is gitignored.

The **Walkthrough video** GitHub Actions workflow runs the same command for every published release and on manual dispatch, and uploads the outputs as a workflow artifact for 30 days.

Requirements: the repository's `dev` extra (Python Playwright 1.61.0 and its Chromium), FFmpeg with libx264, the `brotli` Python package (it converts the app's vendored Inter WOFF2 fonts for English text), and Noto Sans SC for Chinese text. The renderer looks for `NotoSansSC-{Regular,Medium,Bold}.otf` in `$GBDRAW_VIDEO_FONT_DIR` or `~/.cache/gbdraw-video-fonts/`, and then for the `fonts-noto-cjk` package (`/usr/share/fonts/opentype/noto/NotoSansCJK-*.ttc`), which the workflow installs.

## How it works

- `docs/capture/video/walkthrough.py` drives the local, network-isolated app with an eased pointer, typed input, and wheel scrolling. Chromium runs with `--force-device-scale-factor=2`, so its screencast delivers 3840×2160 frames and camera zooms stay sharp. The settings column is widened to 450 px before recording so short form fields show their values. The pointer path, clicks, camera cues (with the target and the drawn diagram's bounds), captions, toasts, waits, and named marks are saved in `raw/walkthrough/events.json`. The downloaded SVG is then re-rendered at increasing magnification for the vector zoom, in both window shapes, and six Gallery source SVGs from `gbdraw/web/gallery/sources/` are rendered for the highlights montage.
- `docs/capture/video/walkthrough_render.py` composes the saved bundle without a browser. It places the app in a window over the gbdraw background, eases the camera between logged targets, draws the pointer and click ripples, and writes captions below the window. Generation waits and the three repeated color rules are fast-forwarded, and a `▶▶` badge marks each fast-forwarded span. `reports/walkthrough-report.json` lists the spans and each video's timeline.
- The highlights are defined in `HIGHLIGHTS` as spans between named marks (for example `labels` to `crowded`), so a new recording with different timings still tells the same story. Generation waits play 2.5× faster again in the highlights, and the badge appears there only at 2× or more.
- The recording hides only the floating feature-search palette, which would otherwise cover the result. Headless Chromium does not paint native `<select>` menus or color pickers, so the recorder logs each select's real option list and each chosen color, and the editor draws the menu or color popover while the pointer picks it.

## Chinese text

`captions.zh-Hans.json` maps each English on-screen string to its Simplified Chinese text, plus a pattern for the download toast. Web UI names such as `Labels › Label Mode` stay in English because the app itself is English. Rendering stops with the missing string if the journey or renderer shows text the file lacks, and `tests/test_video_capture_contracts.py` checks that every string has a translation and that no translation is unused. Have a native speaker review changes to this file.

## Editing

To change camera framing, pacing, overlays, or the highlights cut without recording again, edit `walkthrough_render.py` and run:

```bash
python docs/capture/build_video.py render \
  --recording build/videos/walkthrough/run-001 --out build/videos/walkthrough/render-002
```

Changing the journey itself, including caption text and marks, means editing `_journey` in `walkthrough.py` and running `build` again in a new directory. `tests/test_video_capture_contracts.py` checks that every mark used by `HIGHLIGHTS` is still recorded by the journey.

## Troubleshooting

- If the browser wheel is missing or stale, rerun `python tools/prepare_browser_wheel.py`.
- If Chromium cannot start in a restricted Linux sandbox, allow its local browser process and rerun the same command.
- If the recorder stops with `Could not scroll ... into the settings view` or `did not focus it`, the Web UI layout changed under the journey. Update the locator or framing in `_journey` rather than retrying.
- If rendering fails with a missing Inter glyph, change the caption text; the vendored font covers the Latin subset only.
