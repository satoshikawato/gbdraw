## Purpose

## Verification

## Product behavior semantics

- **Semantics:** `NONE` / `DEFINED` / `UNCERTAIN`
- **PR role:** `EVIDENCE_ONLY` / `DECISION_ONLY` / `IMPLEMENTATION`
- **Change intent:** `NOT_APPLICABLE` / `PRESERVE` / `RESTORE` / `IMPROVE`
- **Basis:** `none` / `BD-###` / `SD-###` / `explicit promise` /
  `compatibility commitment` / `authoritative domain rule` /
  `precise unambiguous valid-use invariant`
- **Rationale:**

Human reviewer comment:

```text
SEMANTIC-DISPOSITION:
  NO_USER_VISIBLE_DIFFERENCE
  COVERED_BY_AUTHORITY
  ORDINARY_REGRESSION
  DOMAIN_DERIVED
  SCENARIO_CANDIDATE
  INSUFFICIENT_EVIDENCE

BASIS:
  <source or explanation>
```

## Architecture impact

Select exactly one:

- [ ] This is not architecture-bearing. Rationale:
- [ ] This is architecture-bearing; the evidence below is complete.

For an architecture-bearing change:

- Why this changed-capability scope is complete:

Humans must declare the complete changed-scope sets and arithmetic below. Add a
row for every changed capability, behavior, or compatibility namespace. Use
`none` for an empty set.

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

The sets and O, T, and P inputs must reproduce every row and total.
`<= 0`, `non-positive`, or `yes` alone is insufficient evidence.

- New module classification:
- Persisted-compatibility evidence and removal condition:
- Performance/scientific-output evidence, when material:
- Deterministic checker evidence:
- Maintainer architecture decision comment permalink, required before merge:

- [ ] No second parser, normalizer, state mirror, lifecycle owner, or canonical pipeline was added.
- [ ] Every superseded implementation was removed in this change.
- [ ] No accepted deterministic violation was added in this implementation PR.
