# Meet gbdraw videos

One recording of the real Web app produces three videos:

- **Walkthrough** (1920×1080, about 2 minutes): load `HmmtDNA.gbk`, generate the map, switch the track preset to Middle and put both strands on one track, turn on labels, switch CDS labels to gene symbols, add four functional color rules, mark the origin-spanning D-loop, inspect a feature, and export the SVG, all at a human pace.
- **Highlights** (1920×1080, about 45 seconds): six Gallery figures, then the key moments of the same recording at 1.4× speed.
- **Vertical highlights** (1080×1920, same cut) for TikTok, RedNote, and other phone feeds. Captions sit above the app window and the brand below it, clear of the top 400 px and bottom 430 px where those apps draw their own controls. The camera re-frames each shot for the tall window: whole-page shots show the diagram, close-ups follow the control being used, a tall target such as the settings column stays near the recorded zoom, and a clicked feature stays framed with its popup.

All three end with the exported SVG magnified 7× and a gbdraw.app card. They are in English only.

## Gallery tutorial videos

Six more journeys record a Gallery tutorial from an empty app to its exported Interactive SVG, as a landscape walkthrough (`gbdraw-tutorial-<id>.mp4`, 1920×1080). Each ends by opening the downloaded `*.interactive.svg` in a new browser tab and using it there, then the gbdraw.app card:

| Journey (`--journey`) | Gallery example | What the video shows |
|---|---|---|
| `HmmtDNA_ATskew` | Human mitochondrial genome (AT skew) | Middle layout on one strand track, 500 nt window, Custom Track Slots with an added AT skew row and inner ticks, gene labels, the RNR2 feature details; in the Interactive SVG, the CYTB details and sequences and a Label search for COX |
| `tobacco-chloroplast` | *Nicotiana tabacum* chloroplast genome regions | Radial labels on both sides, outlined features and a thicker axis, the chloroplast color and label tables, the imported LSC/IRb/SSC/IRa table bound to an annotation track between the genes and GC content; in the Interactive SVG, the rbcL qualifiers and sequences and a Label search for psb |
| `BGC0000708-BGC0000713` | Aminoglycoside biosynthetic gene clusters from *Streptomyces* spp. | The tutorial's 12 steps: five MIBiG files added one at a time with their names, a first Generate with the defaults, LOSATP Similarity groups (LOSAT runs once), BGC0000713 reverse-complemented, the rows aligned on og_6 through Review alignment options…, the orange palette with four gene_kind color rules, labels on the first cluster, the finished figure, the og_18 ribbon, and the livE qualifiers; then, unnumbered, the Interactive SVG with the livE popup |
| `Vnig_TUMSAT-TG-2018` | *Vibrio nigripulchritudo* TUMSAT-TG-2018 | One RefSeq GBFF with two chromosomes and four plasmids, a first map, Tuckin, orchid, the left legend, Multi-Record Canvas with the chromosomes in row 1 and the plasmids in row 2, the plot title, and the topB details, qualifiers, and sequence on pVNTG3; in the Interactive SVG, a Label search for dnaA whose Open zooms to it on chromosome 1, and its sequences |
| `hepatoplasmataceae_collinear` | Hepatoplasmataceae collinear protein-match blocks | Five genomes, LOSATP Collinear blocks in the browser with the tutorial's block rules, filters, and Orientation + identity colors, a centered Middle layout with both GC tracks, a ruler, and ajisai, the block_0024 ribbon, and the DnaA details; in the Interactive SVG, DnaA with its og_1 members lit and its sequences, and a Label search for ribosomal protein |
| `vibrio-harveyi-group-collinear` | Collinearity analysis of multi-replicon bacterial genomes (*Vibrio* spp.) | The tutorial's 10 steps: two RefSeq GBFFs added one at a time with their Names for all records in this file, a first Generate with the defaults (no Accession, Replicon, or Length), LOSATP Collinear blocks (LOSAT runs once), the feature search for dnaA and RctB with each record rotated to start at its initiator, Generate again, then rectangle features, the Above layout, a 750,000 bp ruler on each axis, and the record-label fonts; then, unnumbered, the Interactive SVG with the dnaA popup |

Each tutorial is recorded a second time for RedNote and other phone feeds, as the journey `<id>-vertical` (`gbdraw-tutorial-<id>-vertical.mp4`, 1080×1920).

A recording shows one UI language, and its videos show text only in that language. The BGC and *Vibrio* tutorials are also recorded with the app in Simplified Chinese, as the journeys `<id>.zh-Hans` and `<id>-vertical.zh-Hans` (`gbdraw-tutorial-<id>[-vertical].zh-Hans.mp4`): such a take first switches the header's Language menu to 简体中文 on camera, and the journey finds every control by the app's own text for it (`Recorder.t`, the app's `t()`), so the same steps serve both languages. `walkthrough_render.TUTORIAL_LANGUAGES` lists each tutorial's languages, and `Story.language` gives each journey its one language. The exported Interactive SVG's own controls are English in every take. This portrait take runs the same steps and checks, laid out like the vertical highlights with the step bar just below the app window, but shapes the page for the tall window: the preview figure is zoomed out until it is at most 795 CSS px wide, so every overview frames the whole figure (the extra zoom-out clicks run at 3× under the `▶▶` badge while the camera holds the canvas); a clicked feature's popup is framed alone; and the downloaded Interactive SVG opens in a 954×1080 tab, the window's shape at full page height, where the file fits its drawing to the tab as it does to any window. The recorder logs that tab's pointer and camera at its place on the 1920×1080 stage and pads its frames to the stage after recording, so the renderer treats both tabs alike.

Captions carry the tutorial's own step numbers and UI labels (`gbdraw/web/gallery/tutorials/<id>.json`), so the video can be followed beside the written steps; the export is the last, extra step. The circular tutorials also press Generate at one or two earlier checkpoints so a change is seen as it lands. Repeated entries (records 2–5, rules 2–4, definition-line styles, annotation gaps, record rows) run at 2–3× under the `▶▶` badge, and the wait while a large Interactive SVG opens in its tab is compressed like a generation wait. The videos reproduce the Gallery session, not only the written steps: where the session sets values the tutorial text does not list, the video sets them too and says so in the caption (for example Separate Strands off and misc_RNA removed for AT skew, Multi-Record Canvas off, outlines, the axis width, and rep_origin for the chloroplast, the default palette with CDS #d3d3d3, feature height, outlines, axis, scale, and legend sizes for the BGCs, and Infer orthogroups with self-comparisons for the Hepatoplasmataceae blocks, which the Web app leaves off and the Gallery command's CLI default turns on; the tutorial's group IDs and group counts come from it). The AT skew video turns off Layout › Multi-Record Canvas, as the Gallery session does, because the multi-record canvas drops the legend entries of added skew rows. The exported figures then match the Gallery images, except that the two circular sessions space the center label lines 30 px apart (`objects.definition.circular.interval`), which the Web app has no control for; the videos keep the default of 20. The Hepatoplasmataceae records also show each genome's organism and strain, which the Web app takes from the GenBank file (Record options › Organism / strain, "Using file default"); the CLI that made the Gallery figure draws no name line, and the Web app has no control to hide it. The *Vibrio* figure also draws its three ncRNA features and their legend swatch in the orchid palette's `default` gray, #dddddd: the orchid palette has no ncRNA color, and the CLI that made the Gallery figure falls back to #d3d3d3 instead of the palette's `default`. The intro shows the Gallery example's title from `examples.json`, with taxon names slanted for italics. Each recording checks the frozen tutorial inputs against `gbdraw/web/tutorial-data/manifest.json` and fails if the exported SVG is not interactive or lacks the tutorial's result (for example the AT skew legend and colors, the four region brackets, the four gene_kind colors and similarity links, the six records placed as in the Gallery figure, or the Gallery figure's 500 collinear blocks). The *Vibrio* GBFF is not in the tutorial data package, so that journey pins its size and SHA-256 instead. The opened file is served at a same-origin address, so the recording stays network-isolated, and the journey fails if the file logs an error, a search finds nothing, or the og_1 FASTA does not hold every member.

## 0.14.0 promotion

The journey `promo-zh-0140-vertical` (`docs/capture/video/promo.py`) records one portrait video of what is new in 0.14.0, with the app in Simplified Chinese only (`final/gbdraw-0.14.0-vertical.zh-Hans.mp4`, 1080×1920, about 110 s, no audio track, with its subtitles and the two covers). It switches the header's Language menu to 简体中文 with all 11 languages shown, loads `HmmtDNA_basic_circular` with Load example, loads Noto Sans SC as the figure font with the plot title 人类线粒体基因组, rotates the record to start at COX1 with Apply now, puts ATP6 on Outer lane 1 and shows Generate's stage line, saves the Session and shows Run info's Source recipe and Exact replay commands for that Generate, then loads `BGC0000708-BGC0000713`, resets its og_6 alignment (skimmed) and aligns again from the livU popup's Align…, and ends over the aligned figure. The font comes from `GBDRAW_VIDEO_FONT_DIR` or `~/.cache/gbdraw-video-fonts` (`NotoSansSC-Regular.otf`). Its pointer moves are paced (`Journey.paced`): on a busy page each move keeps its planned duration instead of stretching with the page's latency.

## Build

From the repository root:

```bash
python tools/prepare_browser_wheel.py
python docs/capture/build_video.py build --out build/videos/walkthrough/run-001
python docs/capture/build_video.py build --journey tobacco-chloroplast --out build/videos/tobacco-chloroplast/run-001
```

Without `--journey`, the command records the Meet gbdraw journey, renders `final/gbdraw-walkthrough.mp4`, `final/gbdraw-highlights.mp4`, and `final/gbdraw-highlights-vertical.mp4` (30 fps, H.264), a subtitle file for each, `final/poster.png`, and two covers for phone feeds (`cover-vertical.png` at 9:16 and `cover-3x4.png` for RedNote), and then decodes and checks every MP4. The landscape videos have no audio track; the vertical one carries a silent AAC track because some social upload paths reject video-only files. It fails if the exported SVG lacks the 13 gene symbols, the four functional colors and legend entries, or the D-loop. Use a new or empty `--out` directory for every build; `build/` is gitignored.

To make every video in one command, use `build-all`. It builds each journey in its own process into `OUT/<journey>/`, keeps each build log in `OUT/logs/`, prints PASS or FAIL per journey, and, with `--collect`, copies each journey's `final/` files to `COLLECT/<journey>/` after all builds pass:

```bash
python docs/capture/build_video.py build-all --prepare-wheel --jobs 4 \
  --out build/videos/all/run-001 --collect /mnt/c/Users/<you>/Videos/gbdraw
```

`--journey` (repeatable) limits the run to some journeys. `--jobs N` records N journeys at once (four at once took about half an hour on a 32-thread workstation); the default of 1 is safer on smaller machines. The output and every collect folder must be new or empty, and this is checked before anything runs; collected files are never overwritten and are checked by SHA-256 after copying. `--prepare-wheel` runs `tools/prepare_browser_wheel.py` first, which deletes the whole `build/` directory, including earlier video runs, so move any runs you want to keep out of `build/` before using it.

A tutorial build writes `final/gbdraw-tutorial-<id>.mp4`, its subtitle file, and `final/poster.png`; a `<id>-vertical` build writes `final/gbdraw-tutorial-<id>-vertical.mp4`, its subtitles, and the two phone-feed covers. A `.zh-Hans` journey writes the same files with `.zh-Hans` before the extension. All are checked the same way. For RedNote, upload `gbdraw-tutorial-<id>-vertical.zh-Hans.mp4` with `cover-3x4.zh-Hans.png`. The BGC journey runs LOSATP for all 25 record pairs in the browser, and the Hepatoplasmataceae journey for each genome against itself and its neighbor; those waits are fast-forwarded.

The **Walkthrough video** GitHub Actions workflow runs the same command for every journey for every published release and on manual dispatch, and uploads each journey's outputs as a workflow artifact for 30 days.

Requirements: the repository's `dev` extra (Python Playwright 1.61.0 and its Chromium), FFmpeg with libx264, the `brotli` Python package (it converts the app's vendored Inter WOFF2 fonts for English text), and Noto Sans SC for Chinese text. The renderer looks for `NotoSansSC-{Regular,Medium,Bold}.otf` in `$GBDRAW_VIDEO_FONT_DIR` or `~/.cache/gbdraw-video-fonts/`, and then for the `fonts-noto-cjk` package (`/usr/share/fonts/opentype/noto/NotoSansCJK-*.ttc`), which the workflow installs.

## How it works

- `docs/capture/video/walkthrough.py` drives the local, network-isolated app; `docs/capture/video/tutorials.py` holds the five Gallery tutorial journeys and their portrait takes, and `walkthrough_render.STORIES` names each journey's intro and outputs with an eased pointer, typed input, and wheel scrolling. Chromium runs with `--force-device-scale-factor=2`, so its screencast delivers 3840×2160 frames and camera zooms stay sharp. The settings column is widened to 450 px before recording so short form fields show their values. The pointer path, clicks, camera cues (with the target and the drawn diagram's bounds), captions, toasts, waits, and named marks are saved in `raw/walkthrough/events.json`. A tutorial opens its downloaded Interactive SVG in a second tab and the recorder moves its screencast there (`Recorder.switch_page`), so the file's own popups, search, and zoom are part of the same timeline. For the Meet gbdraw journey, the downloaded SVG is then re-rendered at increasing magnification for the vector zoom, in both window shapes, and six Gallery source SVGs from `gbdraw/web/gallery/sources/` are rendered for the highlights montage.
- `docs/capture/video/walkthrough_render.py` composes the saved bundle without a browser. It places the app in a window over the gbdraw background, eases the camera between logged targets, draws the pointer and click ripples, and writes captions below the window. Generation waits and the three repeated color rules are fast-forwarded, and a `▶▶` badge marks each fast-forwarded span. `reports/walkthrough-report.json` lists the spans and each video's timeline.
- The highlights are defined in `HIGHLIGHTS` as spans between named marks (for example `labels` to `crowded`), so a new recording with different timings still tells the same story. Generation waits play 2.5× faster again in the highlights, and the badge appears there only at 2× or more.
- The recording hides nothing: the app's feature search bar sits above the preview. Headless Chromium does not paint native `<select>` menus or color pickers, so the recorder logs each select's real option list and each chosen color, and the editor draws the menu or color popover while the pointer picks it.

## Chinese text

`captions.zh-Hans.json` maps each English on-screen string of every journey to its Simplified Chinese text, plus patterns for the download toast and for picked file names, which stay unchanged. A zh-Hans take shows the app in 简体中文, so a Web UI name in its captions is the app's own zh-Hans text (`gbdraw/web/js/utils/i18n-catalogs/zh-Hans.js`); look each name up there instead of translating it. Rendering stops with the missing string if the journey or renderer shows text the file lacks, and `tests/test_video_capture_contracts.py` checks that every string has a translation and that no translation is unused. Have a native speaker review changes to this file.

## Editing

To change camera framing, pacing, overlays, or the highlights cut without recording again, edit `walkthrough_render.py` and run:

```bash
python docs/capture/build_video.py render \
  --recording build/videos/walkthrough/run-001 --out build/videos/walkthrough/render-002
```

Changing a journey itself, including caption text and marks, means editing `_meet_gbdraw` in `walkthrough.py` or the tutorial in `tutorials.py` and running `build` again in a new directory. `tests/test_video_capture_contracts.py` checks that every mark used by `HIGHLIGHTS` is still recorded by the journey.

## Troubleshooting

- If the browser wheel is missing or stale, rerun `python tools/prepare_browser_wheel.py`.
- If Chromium cannot start in a restricted Linux sandbox, allow its local browser process and rerun the same command.
- If the recorder stops with `Could not scroll ... into the settings view` or `did not focus it`, the Web UI layout changed under the journey. Update the locator or framing in that journey rather than retrying.
- If rendering fails with a missing Inter glyph, change the caption text; the vendored font covers the Latin subset only.
