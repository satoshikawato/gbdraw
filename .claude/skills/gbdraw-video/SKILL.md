---
name: gbdraw-video
description: Plan, record, edit, and check gbdraw videos - tutorial takes, feature promotions, landscape and portrait (phone) cuts, and localized takes - by the Owner's video rules - example choice, opening highlights, fast-forward of redundant scenes, progressive tutorials, framing, captions, and the checks before delivery. Use before planning or recording any gbdraw video, when changing a recording journey or the renderer, and when reviewing a finished take.
---

# gbdraw videos

These are the Owner's rules for every gbdraw video, with the date each was
given. This skill decides what a video shows and how it is paced. The recorder
(`docs/capture/build_video.py` and `docs/capture/video/`) lives on local `video/*`
branches, which are not pushed. Find the newest one with
`git branch --list 'video/*' --sort=-committerdate`. A video never changes the
app: a bug found while recording gets an OV ID and a normal PR to `dev`.

## Before recording

- **Purpose and kind.** A tutorial teaches one Gallery tutorial step by step. A
  promotion shows a release's features in about 70-100 s. Decide which one you
  are making and who watches it before you write a journey.
- **Build.** Record the released app: the newest release tag, which is what
  users install and what gbdraw.app serves. Do not record `dev`.
- **Examples: show the richest realistic figure.** Use a Gallery example that
  shows tracks and context, not a minimal map. For example, use
  `HmmtDNA_ATskew` (AT content and AT skew), not `HmmtDNA_basic_circular`
  (Owner 2026-10-23). Keep the example's labels, Legend, colors, and tracks, and
  change only the setting being shown (AGENTS.md "Showcase Figure Quality").
- **Language.** A localized take switches the app with the header **Language**
  menu on camera before the first step. Captions name controls exactly as that
  locale's catalog names them. Never call a translation reviewed.

## Structure

1. **Highlights first (Owner 2026-10-23).** Open with the strongest results
   before the title card. A promotion opens with a 5-8 s montage of the
   features' results, one short label per clip. A tutorial opens with its
   finished figure for 2-3 s. Cut the highlights from the same recording; they
   need no second recording.
2. A short title card (about 2 s).
3. The steps, numbered as in the tutorial or the caption sheet.
4. An outro over the final figure, with where to get gbdraw.

## Pacing: fast-forward what is redundant (Owner 2026-10-23, every video)

Fast-forward these, with the fast-forward badge on screen:

- Waiting: runtime start, example and Session loads, uploads, LOSAT runs, every
  Generate except the one that shows the Generate stages, and rerenders.
- Navigation: scrolling, opening sections and popups, and long pointer travel.
- The file picker, and typing after the first few characters.
- Any repeat of an operation the viewer has already seen once, such as a second
  example load, the same control in another record, or an undo before the
  real step.

Keep these at normal speed:

- The first use of each feature's key control: the menu pick, and the click
  that applies the change.
- Each result, held at least 2-3 s so a phone viewer can read it.

## Tutorials (Owner 2026-10-11)

- **Progressive.** Start from the defaults and Generate. Name what looks
  wrong, change it, and Generate again. Show each intermediate result. The
  tutorial does not dictate the final settings up front.
- **Files.** Add input files one at a time, and fill each file's defaults before
  adding the next.
- **LOSAT.** Run LOSAT once, and show that a later Generate reuses its cached
  results.
- **Vibrio.** Vibrio records carry no Accession, Replicon, or Length labels.

## Framing

- **Desktop UI only (Owner 2026-10-23).** Record the desktop UI even for a phone
  video (a 1920x1080 window). Never use a mobile viewport.
- **Portrait cuts (Owner 2026-10-11).** A portrait cut shows a picture-in-picture
  overview of the whole app, with a rectangle that marks the shown region. Keep
  the demonstrated control and the changed part of the figure centered and
  fully in frame.
- **English blocks in a localized take.** Do not let a block of English text
  fill a localized take. Frame around it, and report it as a translation bug.
  Do not hide it with CSS.

## Captions (Owner 2026-10-11, 2026-10-22)

- Break a caption line only after a clause or sentence mark. Never split a
  sentence mid-line, and never split a word over two lines. If a line does not
  fit, shorten it.
- No arrow characters; use `›` for paths.
- Every caption must be true for the frames it covers.

## Checks before delivery

- Extract a frame sheet with one frame per caption and per highlight clip,
  plus evenly spaced frames. Look at it, and at the key moments at full size.
- Check each item:
  - every intermediate result is readable;
  - nothing is cut off;
  - captions follow the break rule;
  - a localized take shows that UI in every frame after the switch;
  - the length is within its target.
- Deliver the MP4 and its `.srt` to a new dated folder, and never overwrite an
  earlier folder. In the report, list the steps, any differences from the
  tutorial or the caption sheet, and the decisions you made.
