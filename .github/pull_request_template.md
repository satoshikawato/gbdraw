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
