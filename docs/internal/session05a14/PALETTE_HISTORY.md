# Preserve unrelated feature colors during placement History

B-02/A-02 is an output-correctness defect separate from Result synchronization.
On base `550540ef7f1bc9132c22162aa314617dd48096ca`, the permanent placement
History test observed unrelated `rep_origin` fills change from the active
palette default `#d3d3d3` to `#cccccc`. All three SVG representations could
agree on that wrong color; coherence alone could not reject it.

`svg-styles.js` now resolves an unmatched type through the existing applied
palette's `default` entry in palette application and specific-rule replay.
The existing palette normalization/loader remains the color authority. Explicit
type colors, specific rules and hash-rule precedence remain intact; there is
no species, feature-ID or color-value special case and no additional owner.

`palette-history-preservation.test.mjs` covers the palette default, another
custom default, explicit type color, matching specific rule and hash precedence.
`palette-history-preservation.playwright.spec.js` compares every non-target
feature after each placement Undo/Redo, including identity, paint, visibility
and geometry. It reuses the first review package's visual comparator. Its red
proof predates production edits and its unchanged preservation assertion passes
with this correction.

The full original J13 program also passes placement, label edit, mode changes,
each History step and regeneration. The retained adapter additionally applies
the same non-target assertion at every individual J13 Undo and Redo, independent
of selected/mounted/download agreement. Immutable source/seed references and
common reproduction instructions are in [RESULT_COMPLETION.md](RESULT_COMPLETION.md).
