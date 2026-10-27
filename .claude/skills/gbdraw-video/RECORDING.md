# Recording a gbdraw video reproducibly

How to get the best take, and how to get it again. `SKILL.md` decides what a
video shows. This file covers how to record it, check it, and deliver it. The
lessons come from the 0.14.0 promotion and tutorial takes (2026-10-23).

## 1. Set up

- **One recorder.** The recorder (`docs/capture/build_video.py`,
  `docs/capture/video/`) lives only on local `video/*` branches.
  - Start from the newest branch that the other video work builds on
    (`git branch --list 'video/*' --sort=-committerdate`), and fast-forward or
    merge it. Never copy recorder code.
  - Promotion and tutorial journeys share the helpers in
    `docs/capture/video/walkthrough.py`, for example `switch_language` and
    `_skimmed`. A step that two journeys need is one helper with parameters,
    not a copy in one journey.
- **Other video sessions.**
  - When another session works on a video branch, send it each commit hash
    you make: a message, plus one line in its task folder's `LOG.md`.
  - Merge its hashes before you record. Say which of your changes alter every
    render or every caller of a helper.
- **Environment.**
  - Use your own worktree with the prepared browser wheel.
  - Put the Noto Sans SC fonts in `~/.cache/gbdraw-video-fonts` (or set
    `GBDRAW_VIDEO_FONT_DIR`).
  - Set `PYTHONPATH=<worktree>`, and point `TMPDIR` at the task folder.

## 2. Make the take deterministic

Journey rules:

- **Wait for state, then hold.** Wait for the app's state, never for a fixed
  time: `uiLocale.value === code`, `document.fonts.status === 'loaded'`, the
  result drawn. Then hold.
- **Camera cue first.** Start each step with a camera cue. Until the first cue,
  the recording shows the whole page.
- **Clear the figure.** Move the pointer off the figure before each result
  hold.
- **One framing per sequence.** Measure a framing once and reuse it for every
  version of a sequence. A framing measured again after each Generate drifts.
- **Edges.** Keep every label and leader line of the change inside the frame;
  check the outermost one. A toolbar or search box at the frame edge is either
  whole or out of the frame.
- **Captions.** A caption set for a sub-step (such as the Generate stages) is
  set back to the step's caption after it. A caption must be true for its
  whole span.
- **Marks and opening clips.** Highlights are data:
  - Set `rec.mark("hl-…")` at the end of each result's hold, and list the
    clips in `Story.opening` as `Clip(mark, start, mark, 0.0, label)`.
  - End a clip at its mark. A clip that ends 0.5 s before its mark showed a
    switched UI for only 0.3 s.
  - Consecutive clips with the same label play as one highlight: the label
    once, hard cuts, and one SRT cue.
- **Version sequences** (SKILL.md, "Show a change as it happens"):
  - Record the first version at normal speed in its step. Record the others as
    fast-forwarded repeats, each ending with a short hold and a mark.
  - Make the moves in one direction and cumulative, so the change builds up
    visibly. For lanes, each moved gene stays out.
  - Frame tightly on what changes, with the Legend and center text outside the
    frame.
- **Language switch.** Use `walkthrough.switch_language`. It shows the settings
  pane in English, the menu pick, then the same pane in the target language.
  It sets `mark("language")` before the menu opens.
- **Drawn menus.** `Recorder.choose` keeps its drawn menu open until the click,
  because on a busy page the pointer can take seconds to reach the row.

## 3. Record and re-render

- **Record.** Each take goes to a new folder; never reuse a run folder. Start
  it in the background and wait once for it to finish:

  ```bash
  PYTHONPATH=$PWD /home/kawato/gbdraw-baselines/bin/heavy \
    python docs/capture/build_video.py build --journey <id> --out <task>/videos/runNN
  ```

- **Edit; do not record again (Owner 2026-10-23, repeated 2026-10-26).**
  Fix what a review finds by editing the recorded material: re-render it. It
  takes minutes and keeps the same frames:

  ```bash
  python docs/capture/build_video.py render --recording <run> --out <new folder>
  ```

  - Edits include clip offsets, holds, speed, cuts, skims, opening clips, the
    outro text, labels, captions, fades, the
    badge, the camera framing (crop, zoom, letterbox), and where the pointer
    is drawn.
  - When the renderer cannot make such a fix from what the take logged,
    change the renderer. If the take lacks the geometry for it, make the
    recorder log that geometry for the next video. Do not record again for it.
  - **Take edits.** A take is read-only. Add what it lacks to `TAKE_EDITS` in
    `docs/capture/video/walkthrough_render.py`: one `TakeEdits` entry per
    take, keyed by the journey and the SHA-256 of its `walkthrough.json`, with
    each edit anchored on the take's own events (`At(mark=...)`, or
    `At(title=...)` for a caption, `press=n` for the n-th press after it).
    - `AddMark` adds a mark (with `view`, the rows a clip at it shows).
    - `AddSkim` and `Retime` add or change a fast-forward.
    - `CloseRead` reads a still result hold closely: it plays the hold longer
      under camera stops.
    - `ShowOnly` makes recorded camera stops show only their rect, cut to
      `within` when an app control overlaps it.
    - `render` applies them and counts them in its report. It fails when an
      anchor is missing or ambiguous, a span overlaps another, or a still's
      frames change.
    - Give the journey the same marks and skims in the same change, so the
      next take needs no entry.
  - Record a new take only for app operations that the take never captured.
    First collect every such fix from one review, and tell the Owner why an
    edit cannot make it. Then record one take for all of them.
  - Put this rule in every brief that starts recordings.
- **Machine load.** Check `heavy --status` first.
  - A take takes 8-10 minutes when the load is about 40.
  - Its length varies by several seconds from take to take; three takes of one
    journey ran 109.5-121.8 s.
  - "Chromium screencast produced no frames" at the start is load: retry once.
- **Keep every take.** Never delete a take before delivery. If there is more
  than one, say in the handoff why you dropped each of the others.

## 4. Review a take

Run the review script and read both sheets:

```bash
python3 .claude/skills/gbdraw-video/scripts/review_take.py <mp4> <srt> <out-dir> \
  --opening <seconds> --min <s> --max <s> --no-audio
python docs/capture/build_video.py check --out <run>   # the recorder's own decode check
```

- **`cues.png`** has the start, middle and end of every caption cue, plus a
  frame every 3 s in a long cue, each with its cue text. Every frame must match
  its text. A long cue is where a caption spans actions it does not describe.
  The script warns about cues longer than `--max-cue` (12 s).
- **`opening.png`** samples the opening every 0.25 s. Each version of a change
  highlight must differ visibly at phone size.
- **Failures.** The script exits 1 on:
  - SRT overlaps, or a cue that ends after the video;
  - a length outside `--min`/`--max`;
  - a size other than 1080x1920 or 1920x1080, or a codec other than h264;
  - an audio stream with `--no-audio`.
- **Full size.** Then look at these at full size:
  - each key result;
  - every frame edge (cut labels, toolbars);
  - in a localized take, the frames after the switch;
  - the fast-forward badge, which must appear only on skims, while first uses
    play at normal speed.

## 5. Deliver

- Copy the finals to the task folder's `final/`.
- Then copy them to a new folder for each round
  (`<Videos>/<name>-r<N>/`), never over an earlier delivery, and `cmp` every
  file.
- In the report, list:
  - what changed;
  - the takes you dropped and why;
  - the decisions you made under the Owner's standing instruction.
