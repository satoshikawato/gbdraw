# Product Decision Pack — <title>

Status: non-normative working template

Use this template for either a registered Product Impact change or developer
preflight that finds unmapped material user-effect ambiguity. The normative
lifecycle, classifications, authority rules, and routes are defined by the
[Product Impact Ratchet](./PRODUCT_IMPACT_RATCHET.md); this working document is
not Product authority.

## Identity

- Concern key: <stable dotted key>
- Scenario revision: <positive integer>
- Discovery lane: <Product Impact Ratchet | developer preflight>
- Prepared from base SHA:
- Prepared for proposed branch/head:
- Prepared by:
- Related issue/PR:

## Trigger

<What proposed change exposed the unresolved Product outcome?>

## Authority search

| Source inspected | Result | Conflict or gap |
| --- | --- | --- |
| Product Impact concern/map | | |
| Active durable `BD-###` | | |
| Static Product Contract, when available | | |
| Domain/scientific/integrity rule | | |
| Released compatibility evidence | | |
| Eligible exact-head current decision | | |
| Current code/tests as evidence only | | |

Result: <UNRESOLVED | CONFLICT | INSUFFICIENT_EVIDENCE>

Procedural classification:
<IMPLEMENT_EXISTING_AUTHORITY | EVIDENCE_REQUIRED | PRODUCT_DECISION_REQUIRED | NOT_ALLOWED>

## User journey

- Actor:
- Context:
- Goal:
- Journey ID/title:
- Entry point:
- Preconditions:
- Major steps:
- Affected checkpoints:
- Failure/recovery path:
- Persistence/export implications:
- Next available action after each candidate outcome:

## Current observed behavior

<Describe observable behavior and evidence. Do not call it correct merely
because it exists.>

## Non-waivable constraints

- Architecture:
- Security/privacy:
- Scientific correctness:
- Persisted compatibility:
- Performance/resource safety:
- Required contracts/evidence:

## Choice A — <stable outcome ID and title>

Describe the complete product outcome without using implementation filenames as
the choice.

- Complete normative outcome:
- Preserved effects:
- Added effects:
- Lost effects:
- Retired effects:
- Discoverability/accessibility:
- Canonical state update:
- Undo/Redo consequence:
- Session/regeneration consequence:
- Export/artifact consequence:
- Validation/error consequence:
- Failure/recovery consequence:
- Scientific-output consequence:
- Cache/provenance consequence:
- Performance consequence:
- Compatibility consequence:
- Architecture consequence:
- Evidence available/missing:
- Residual risk:
- Route: <PR_LOCAL_ALLOWED | DURABLE_AUTHORITY_REQUIRED | EVIDENCE_REQUIRED | NOT_ALLOWED>
- Next action if selected:

## Choice B — <stable outcome ID and title>

<Repeat every Choice A field. Add Choice C only when it is materially
distinct.>

## Comparison matrix

| Dimension | Choice A | Choice B | Choice C |
| --- | --- | --- | --- |
| Entry point / discoverability | | | |
| Immediate feedback | | | |
| Canonical state update | | | |
| Undo / Redo | | | |
| Session round trip | | | |
| Regeneration | | | |
| Export / output | | | |
| Validation / error | | | |
| Failure recovery | | | |
| Accessibility | | | |
| Performance | | | |
| Compatibility | | | |
| Scientific-output consequence | | | |
| Cache / provenance consequence | | | |
| Preserved / added / lost / retired effects | | | |
| Residual risk | | | |
| Route | | | |
| Next available action | | | |

## Evidence-first option, when applicable

- Question evidence must answer:
- Representative deterministic inputs:
- Compared explicit values/outcomes:
- Required measurements and outputs:
- Limitations:
- Runtime/default/authority changes prohibited in the evidence PR:

## Engineering recommendation

- Recommended option, if any:
- Engineering reason:
- Explicit warning: This recommendation is not Product authority.

## Product Decision Owner response

```text
PRODUCT_DECISION
Concern: <concern key/title>
Scenario revision: <revision>
Choice: <stable choice code and outcome ID>
Rationale: <product-level reason>
Must preserve: <effects and affordances>
May retire: <none or explicit scope>
Accepted residual risk: <bounded risk or none>
Owner: <maintainer identity>
Decision date: <YYYY-MM-DD>
```

A choice letter or value alone is insufficient. Complete wording controls. A
coding agent may serialize this explicit receipt but must not infer omitted
rationale, preservation intent, retirement scope, risk, owner, or date.
