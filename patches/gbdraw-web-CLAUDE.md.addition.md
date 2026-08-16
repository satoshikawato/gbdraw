# Addition for gbdraw/web/CLAUDE.md

Insert this after the live-edit ownership section:

```markdown
## Behavior-preserving refactor invariants

Renderer-owned state and editor-owned overrides are separate provenance domains.
Clearing an editor override restores the renderer baseline; it must not infer or
overwrite that baseline from the current DOM.

Validated large owners such as Feature catalogs, extracted Features, biological
Features, orthogroups, sequences, SVG, Results, and LOSAT payloads are borrowed
read-only unless a boundary explicitly documents ownership transfer. Do not
JSON-clone or JSON-serialize them as defensive copying.

Mounted SVG and selected Result ownership is synchronous by default. It must not
cross an `await` unless an explicit Result/SVG/document lease is revalidated
before mutation and commit.

An empty editor projection is an identity operation on renderer output.

Regression tests observe production state, production counters, or independently
derived output. Locally initialized dummy counters and comparison of two paths
using the same new helper are not behavioral proof.

Behavior-preserving refactors establish characterization against the base before
production changes, record base/head differential evidence, and receive a
separate adversarial review before merge.
```
