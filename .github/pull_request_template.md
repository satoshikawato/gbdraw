<!--
Use a concrete action and object in the PR title. Assume the reader has not read the implementation plan or earlier PRs.
-->

## Plain-language summary

<!--
In 2-4 sentences, explain what changes, why it changes, and what will be different after merge.
Keep exact check names, paths, and identifiers in backticks, then explain their role in ordinary language.
-->

## Change class

Select exactly one:

- [ ] STANDARD
- [ ] GOVERNANCE
- [ ] ARCHITECTURE_EXCEPTION
- [ ] PROMOTION

## Purpose

## User-visible impact

## Changes

## Verification

## Risk and review notes

- Gate result and any blocker:
- Review result and why it is `CLEAR` or `REQUIRED`:
- Architecture impact: state whether **This is not architecture-bearing** or
  **This is architecture-bearing**, and identify the affected owner, path, or
  responsibility.
- `architecture-change` label, if relevant:

## Rollback

## Conditional Product Impact

<!--
Use `N/A` and leave `decisions` empty for the normal path.
Nonempty decisions are Product Decision Owner-only and bind to the exact head SHA.
The human selects an outcome with the documented `PRODUCT_DECISION` response; Codex fills the machine-only fields and shows the generated block for review.
Every nonempty current decision requires a product-level `rationale`.
This block does not waive any other gate. Product change or retirement requires durable base authority.
-->

- Product-impact role: N/A | EVIDENCE_ONLY | DECISION_ONLY | IMPLEMENTATION
- Product preflight classification: N/A | IMPLEMENT_EXISTING_AUTHORITY | EVIDENCE_REQUIRED | PRODUCT_DECISION_REQUIRED | NOT_ALLOWED
- Affected concern(s), if known:
- Authority and decision references:
- Decision Pack or evidence reference, if required:
- Behavior contracts and residual risk:

<!-- gbdraw-product-impact-decision:start -->
{"schemaVersion":1,"headSha":"","decisions":[]}
<!-- gbdraw-product-impact-decision:end -->

## Conditional evidence

Complete only the section for the selected change class.

### GOVERNANCE

Complete only when `GOVERNANCE` is selected.

- Authority or evidence-producer files touched:
- Checker implementation files touched:
- Self-authorization separation evidence:
- Branch-protection or ruleset impact, including unchanged settings:
- Governance-specific rollback or protection restore point:

### ARCHITECTURE_EXCEPTION

Complete only when `ARCHITECTURE_EXCEPTION` is selected.

- Exception trigger and why ordinary non-increasing evidence is insufficient:
- Why this changed-capability scope is complete:

Owner-excess rows (repeat for every changed capability):

| Row ID | Capability | Semantic owners before | O before | T before | OE before | Semantic owners after | O after | T after | OE after | delta(OE) |
| --- | --- | --- | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: |
| OE-1 |  |  |  |  |  |  |  |  |  |  |

Path-excess rows (repeat for every changed behavior):

| Row ID | Behavior | Canonical production paths before | P before | PE before | Canonical production paths after | P after | PE after | delta(PE) |
| --- | --- | --- | ---: | ---: | --- | ---: | ---: | ---: |
| PE-1 |  |  |  |  |  |  |  |  |

Compatibility-burden rows (repeat for every changed compatibility namespace):

| Row ID | Namespace | Compatibility path stable IDs before | CB before | Compatibility path stable IDs after | CB after | delta(CB) |
| --- | --- | --- | ---: | --- | ---: | ---: |
| CB-1 |  |  |  |  |  |  |

Changed-scope totals:

- `OE before -> OE after; delta(OE) = <integer>`:
- `PE before -> PE after; delta(PE) = <integer>`:
- `CB before -> CB after; delta(CB) = <integer>`:

List superseded sets separately:

- Superseded semantic owners:
- Superseded canonical production paths:
- Superseded compatibility paths, with stable IDs:

The sets and O, T, and P inputs must reproduce every row and total. `<= 0`, `non-positive`, or `yes` alone is insufficient evidence.

- Best non-exception alternative and why it was not selected:
- Expiry or measurable removal condition:
- Persisted-compatibility, deterministic-check, performance, and scientific-output evidence, as applicable:
- Maintainer architecture decision comment permalink:
- Reviewed exact head SHA:

### PROMOTION

Complete only when `PROMOTION` is selected. A promotion contains no
implementation changes.

- Base (`main`) SHA:
- Head (`dev`) SHA:
- Exact-head integrated staging evidence and URLs:
- Exact-head promotion-readiness evidence and URLs:
- Source-coverage and result-tree identity proof:
- Required-check/protection notes:
- Release verification and deploy notes:
