# Product Impact Ratchet MVP acceptance record

Status: accepted on the audited `dev` baseline

Date: 2026-08-27

This record applies the final acceptance checklist in
`docs/internal/_local_v014_pr_instructions/07_FINAL_ACCEPTANCE_AND_FUTURE_CONCERN_TEMPLATE.md`
to the six-PR Product Impact Ratchet implementation. The decision is `GO` for
the two-rule pilot on `dev`.

## Repository state

| Field | Audited value |
| --- | --- |
| `dev` SHA | `45c0beb034018d9eb868a7f9bbe2c61586a8cd5a` |
| `main` SHA | `5ff8c95c3fbfc26d7cdf175e6c20b5cfeb45c4af` |
| Required `dev` checks | `Web base policy (trusted base)` and `PR / gate` |
| Product Impact map schema | `1` |
| Decision authority schema | `1` |
| Pilot concern | `product.canonical-render-request-boundary` |
| Pilot rules | `canonical-path.render-request`; `semantic-owner.render-request` |
| Pilot enforcement | `hard`; coverage is `complete` for both rules |

GitHub branch protection reported exactly the two required checks above with
strict up-to-date enforcement. No repository ruleset adds another status.

## Implementation and separation history

| Stage | Pull request and merge | File boundary confirmed |
| --- | --- | --- |
| Normative policy | PR #398, `413dcbfd48bac11884db7ef278a3aec2db407d28` | Policy, guidance, and PR template only |
| Pure mechanics | PR #400, `fa635381c17a36802e105ab13be48ecce5e2e9b0` | Pure evaluators, parser, guard classification, and fixtures; no authority or workflow |
| Inert pilot authority | PR #401, `4ece3e972921edea77f644f3745266cf83f33325` | `tools/web-product-impact-map.json` and `tools/web-product-decisions.json` only |
| Trusted report-only integration | PR #402, `04689eabbf84d80b6288ae0d9ed822a8bedfe594` | Existing checker and integration contracts; no authority or workflow |
| PR-body edit rerun | PR #403, `42e4dac7afabf7cd2000ffbd37e8482693e92a43` | Trusted workflow and its contract only |
| Hard evidence sequence | PRs #404, #405, #406, and #407 | Direct PR-gate evidence, fixture registration, map registration, and enforcement-independent fixtures remained separate |
| Pilot hard enforcement | PR #409, `0d73ff30ea6c51ea70b6666ccf2179836835994f` | One map file; exactly two `report-only` to `hard` value changes |

The history preserves the required order. Candidate authority did not admit the
same candidate runtime, and no mapped contract replacement proved a same-PR
hard runtime change.

## Ownership audit

| Responsibility | Canonical owner | Result |
| --- | --- | --- |
| Architecture source facts | `tools/web-architecture-detectors.mjs` | Confirmed |
| Structured architecture delta | `tools/web-architecture-evaluation.mjs` | Confirmed |
| Architecture authority | `tools/web-architecture-rules.json` | Confirmed |
| Product mapping authority | `tools/web-product-impact-map.json` | Confirmed |
| Durable decisions and maintainer allowlist | `tools/web-product-decisions.json` | Confirmed |
| Product evaluation | `tools/web-product-impact-evaluation.mjs` | Confirmed |
| PR-body parsing | `tools/web-product-impact-decision-source.mjs` | Confirmed |
| Git, file I/O, reporting, and Gate aggregation | `tools/check-web-change-budget.mjs` | Confirmed |
| Trusted admission | `.github/workflows/web-base-policy.yml` | Confirmed |
| Normative Product Impact policy | `docs/internal/PRODUCT_IMPACT_RATCHET.md` | Confirmed |

The pure modules have no Git, filesystem, environment, reporting, or CLI
ownership. `tools/check-web-change-budget.mjs` remains the only Web policy CLI.

## Schema, model, and authority audit

Direct validation returned `{ "valid": true, "errors": [] }` for both the
Product Impact map and decision authority. The active map has these properties:

- Both hard rules use complete coverage, and their allowed subjects are mapped.
- `REQ-RENDER-REQUEST-ENTRY` and `REQ-RENDER-REQUEST-OWNER` remain separate.
  The option requires both; alternatives are permitted only within one
  requirement's `anyOf` list.
- The stable effects refer to `generate-submit` and
  `roundtrip-regeneration`. They describe supported user continuations rather
  than function calls.
- Static authority selects `canonical-typed-request-boundary` from the runtime
  data-flow and request/session promises in `gbdraw/web/CLAUDE.md`.
- `tools/web-product-decisions.json` contains only the human maintainer login
  `satoshikawato` and no durable decision.
- The two blocking requirement checkpoints have direct `PR_GATE` evidence.
  The additional `fresh-result` checkpoint has `DEV_STAGING` browser evidence.

Strict validators reject unknown fields, invalid cross-references, unsorted or
duplicate authority arrays, `hard` with partial coverage, unmapped complete
subjects, manual-only hard evidence, stale revisions, and incompatible decision
scope.

## Trust and workflow audit

The trusted workflow checks out the pull request base SHA with credentials
disabled, fetches the head only as Git data, and runs the base checker. It does
not check out or import candidate code. The PR body parser accepts one bounded
JSON block and rejects duplicate markers, malformed JSON, stale heads, unknown
fields, oversized inputs, and non-`AFFORDANCE_PRESERVED` current decisions.
Errors are sanitized and do not echo the full body.

Only `.github/workflows/web-base-policy.yml` declares the `edited` activity.
`.github/workflows/test.yml` omits it, so a body edit does not schedule the
candidate test graph after the workflow definition reaches the default branch.
The audited `main` SHA predates PR #403. GitHub therefore continues to use the
older default-branch event list until the next normal `dev` to `main` promotion.
The source contract and decision rerun path are fixture-proven on `dev`; this
deployment timing does not change the exact-head staleness rule.

## Gate and smoke evidence

`node --test tests/web/architecture-contracts.test.mjs` passed 131 tests in
103.160 seconds after the documentation update. The suite covers the required
acceptance scenarios:

| Scenario | Executable evidence | Observed result |
| --- | --- | --- |
| A: unrelated PR | `trusted Product Impact stays absent for unchanged source and non-authoritative pull requests` | No Product Impact section or Gate change |
| B: provider substitution | `trusted Product Impact reports a mapped provider substitution as one conforming concern` | `NO_USER_VISIBLE_DIFFERENCE`, `CONFORMING` |
| C: lost separate requirement | `separate jointly required contributions cannot be hidden by a shared option ID` | The option becomes unrealized |
| D: effect-equivalent transition | `Decision Pack routing is outcome-oriented and accepts only an eligible exact-head decision` | Unresolved without authority; conforming with an eligible current decision |
| E: stale decision | `stale and non-maintainer Product Impact decisions remain unresolved and routed` | The old head is ineligible |
| F: product change or retirement | `a stable effect change or retirement is not eligible for a current decision` | Durable base authority is required |
| G: authority conflict | `conflicting static and current Product Impact authority is report-only but explicit` and the hard matrix | `AUTHORITY_CONFLICT`; hard Gate fails |
| H: candidate self-authorization | `candidate Product Impact authority is validation-only and cannot authorize candidate runtime` | Base authority remains decisive |
| I: mapped contract modification | `mapped contract changes are surfaced while unrelated tests add no Product Impact burden` and the hard matrix | `CANDIDATE_MODIFIED`; hard Gate fails when it is sole evidence |
| J: human serialization and body edit | Decision Pack/current-decision fixtures plus `workflow triggers separate dev admission, dev staging, promotion, and deployment` | Exact-head JSON is validated; only trusted policy has `edited` in the `dev` source |

The prospective hard matrix passed all four blocking observations:
`ORDINARY_REGRESSION`, `UNRESOLVED_DECISION`, `AUTHORITY_CONFLICT`, and
`INSUFFICIENT_EVIDENCE`. Existing Architecture Ratchet failures remain
independent blockers. Report-only fixtures remain explicit and do not inherit
the active hard values.

## Verification and performance

The acceptance run produced these results:

| Command or measurement | Result |
| --- | --- |
| `node --test tests/web/product-impact-ratchet-fixtures.test.mjs` | Passed; 253 ms |
| `node --test tests/web/architecture-ratchet-fixtures.test.mjs` | Passed; 124 ms |
| `node --test tests/web/architecture-contracts.test.mjs` | Passed; 131 tests; 103.160 s |
| Direct map and decision validation | Both valid; zero errors |
| `node tools/check-web-change-budget.mjs --base HEAD --head HEAD` | `Gate: PASS`; `Review: CLEAR` |

The unchanged checker took 11.87 seconds at the pre-integration PR #401 merge
and 12.08 seconds at the audited `dev` baseline on the same machine. The 0.21
second wall-time difference is below the one-second pilot budget. The checker
does not use a network or browser, and it does not add another repository scan.

Ordinary PR authors add no Product Impact fields beyond leaving the existing
decision array empty. No label, workflow, status context, service, database,
browser run, or external decision UI was added.

## Acceptance and rollback

Decision: `GO` for the two-rule `dev` pilot.

All map, realization, authority, evidence, trust, Gate, and low-friction checks
required for hard enforcement passed on the audited baseline. The default
branch promotion timing above remains visible so the body-only rerun behavior
is not overstated.

Rollback changes only these two values in
`tools/web-product-impact-map.json`:

```diff
-      "enforcement": "hard"
+      "enforcement": "report-only"
```

Apply that replacement once for `canonical-path.render-request` and once for
`semantic-owner.render-request`. Architecture Ratchet enforcement and all
traceability data remain active after the rollback.
