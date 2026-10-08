---
name: gbdraw-promotion
description: Prepare the evidence for a gbdraw dev-to-main promotion - the promotion checklist, carried CI evidence, deploy-only suites, and version literals. Run only when the Owner asks for promotion preparation.
disable-model-invocation: true
---

# Prepare a dev-to-main promotion

The Owner decides when to promote and opens or approves the promotion PR. This
skill gathers the evidence. The authority is the "Promotion PR checklist" and
"Carrying evidence forward" in
[`docs/internal/WEB_PERIODIC_AUDIT.md`](../../../docs/internal/WEB_PERIODIC_AUDIT.md),
and "Carrying evidence to a later commit" in
[`docs/internal/SELECTIVE_CI.md`](../../../docs/internal/SELECTIVE_CI.md).

1. Fix the candidate head `H` (the `dev` tip) and the newest ancestor `E` that
   already has release-tier evidence.
2. Run `node tools/ci-impact.mjs classify --base <E> --head <H>`. Redo only the
   evidence whose verdict (`releaseEvidenceCarries`,
   `generatedArtifactChecksCarry`, `localTestEvidenceCarries`) is false; link
   E's evidence and the classify output.
3. Run the suites that run only on `main` or on deploy before promoting:
   dispatch `deploy_web.yml` on `dev`, or run its suites locally (for example
   `npm run test:web:vibrio-generate`).
4. Grep `tests/`, `docs/capture/`, `docs/recipes/`, and tool fixtures for
   schema or version literals that changed since `main`.
5. Fill every item of the checklist with links. Report each item that is not
   met, and the evidence for it, to the Owner instead of working around it.
