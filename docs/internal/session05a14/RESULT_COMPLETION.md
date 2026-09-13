# Result completion and layout corrections

The candidate fixes B-01/A-01, B-04, B-05 and B-06 in the existing style,
Result and composition owners. Palette fallback (B-02) and mobile popup access
(B-03) are separate dependent review packages.

The failing production witness is `550540ef7f1bc9132c22162aa314617dd48096ca`
(tree `62643107805b61586ea55affa179356ee664ce5a`). Fetching origin before
implementation and again before review preparation found the same dev SHA.
All work took place in `/tmp/gbdraw-session05a14`; the user's dirty workspace
was not edited. A/B evidence archives remain immutable.

## Cause, correction and authority

- **B-01/A-01 and B-04:** `svg-styles.js` returned after marking the mounted
  Result dirty. History replay and depth visibility could finish with stale
  `Result.content`. Its existing `persistSvgEdit` now always serializes and
  replaces the selected Result synchronously. The dirty-only branch and its
  caller parameter are removed. No request or draft History is rewritten.
- **B-05:** after Web Load's delayed definition callback, composition correctly
  changed the legend by 1 px / 0.5 px and canvas width by 1 px. A later Generate
  could repeat the previous initial SVG bytes; the old Vue key and cached
  mounted-content string then reused the edited DOM for a new Result. The
  template now keys the mount by the existing committed artifact identity,
  projected through `state.js`. Live edits retain that identity. Delayed
  definition responses also check their scheduled revision, Result identity,
  selected index and mounted root before applying changes.
- **B-06:** incremental drag binding reused the previous SVG's scale baseline
  (100) for a replacement SVG whose renderer-owned translation was 144.5.
  `diagram-drag.js` now adopts the replacement element's geometry while
  preserving same-element user offsets. It does not copy stale geometry into
  the new preview.

The base Web guidance's [live-edit boundary](../../../gbdraw/web/CLAUDE.md)
requires the mounted target and current Result to update before an action
completes, without using Save, Generate or drawer close as edit triggers.
The existing Result ingestion identity and composition/drag ownership determine
mounting and generated geometry. These corrections implement that authority;
they do not select a new Product outcome or cite a candidate `BD-###`.

Architecture review uses the ordinary non-increasing path: the same style,
ingestion, definition and composition owners remain; no canonical store,
compatibility reader, Worker, watcher, polling loop or privileged importer is
added. The new computed value only projects existing ingestion identity through
an already permitted import. The style owner keeps its existing serialization
path and removes the alternate dirty-only completion. Architecture tests pass;
the policy Gate passes and Review remains REQUIRED for the reactive declaration.
A maintainer must supply that review.

## Permanent evidence

Before production edits, all ten browser reproductions failed at their intended
actions in `red-confirmed.log` under `/tmp/gbdraw-session05a14-evidence`.
The corresponding Result assertions now pass:

| Family | Discoverable test | Original program | Original failed property |
| --- | --- | --- | --- |
| B-01/A-01 | `visual-state-regressions.playwright.spec.js`, Circular and Linear color Redo | J09, J10 | selected fill stayed stale after Redo |
| B-04 | same file, Circular and Linear depth OFF/ON | J27, J28 | selected ancestor lacked mounted/exported `display="none"` |
| B-05 | `definition-replay-visual-state.playwright.spec.js`, single and two-source composite | C01 | legend ancestor transform and root width differed after completed layout |
| B-06 | `visual-state-regressions.playwright.spec.js`, rejected FASTA recovery and label regeneration | J36 | scale line and `5 kbp` parent translation differed by 44.5 px |

`definition-layout-completion.test.mjs` controls the timer and Worker response,
proving both legitimate delayed completion and stale-response rejection. Browser
B-05 waits for the existing definition-completed event, not equality or a fixed
700 ms delay. The fast CLI path is checked separately. The existing layout
delay and test routing/timeouts are unchanged.

`helpers/svg-visual-semantics.mjs` promotes B's visual projection into normal
tests. It compares biological/rendered identities, multipart multiplicity,
geometry/paint/text, depth/legend/scale nodes, ancestor transforms/style and
root dimensions/viewBox. Deterministic controls reject the retained fill,
display, 1 px / 0.5 px / 44.5 px, and missing/duplicate-part mutations. They also
reject a coherently wrong non-target color independently of representation
agreement. Numeric notation normalizes without rounding; only existing
transient UI styling, default opacity 1, and missing-to-valid-catalog label
binding enrichment receive narrow equivalence treatment.

`helpers/visual-state.cjs` records completed selected/mounted SVG before export,
then the actual public download and adjacent state separately. It never uses a
repair action to observe completion. `helpers/retained-visual-journeys.py`
imports this same projection and comparator into the immutable B statement
runner. Original operation programs, failure retention, edit postconditions,
per-Undo/Redo observations and export-boundary capture remain in use.

Run the affected programs with a local server on port 4194:

```bash
python tests/web/helpers/retained-visual-journeys.py \
  --archive /path/to/gbdraw_v014_session05a13_b_evidence_2026-09-12 \
  --output /path/to/new-evidence --root "$PWD" \
  --journeys J09 J10 J27 J28 J36 C01
```

The adapter disables bytecode cache writes and verifies original source and seed
hashes before execution. It can
serve the retained runner's other explicit journey IDs later; this remediation
runs only the eight affected original journeys and supplemental C01.
