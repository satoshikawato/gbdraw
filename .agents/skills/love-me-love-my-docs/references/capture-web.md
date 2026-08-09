# Reproducible Web capture

Use Playwright against a local or explicitly approved test host. Keep the base
URL, viewport, device scale factor, locale, timezone, theme, and network policy
in one configuration owner.

## Flow contract

1. Start from a fresh browser context. A restore chapter may reload only the
   session that its own earlier steps created from original source inputs; do
   not preload a bundled Gallery session or finished project.
2. Reach the UI through normal navigation.
3. Upload files and change values through visible, accessible controls.
4. Prefer `get_by_role`, then `get_by_label` for a form control with an
   associated accessible label, then a stable `get_by_test_id`. A label locator
   is compliant and does not imply that the app needs a test ID. Record any
   surviving CSS or bare-text selector as an application finding for missing
   accessible semantics or a stable test ID.
5. Wait for a semantic ready/result condition rather than an arbitrary delay.
6. Assert the state named by the step before capturing.
7. Capture the smallest truthful region that keeps the operated control,
   selected value, and locating context readable.
8. Capture downloads through the browser event and validate the saved file.

Do not inject a completed application state, call private generation methods,
or hand-crop a one-off image. Capture-only highlights or overlays may improve
readability when they do not change application state.

Block unapproved external requests. Treat every frame as public: no real user
data, credentials, local absolute paths, tokens, or transient notifications.

For a new screenshot system, use device scale factor 2 or higher for UI text.
An established repository may retain a lower factor only when one documented,
tested capture contract owns it and the final text remains readable at its
published size; treat a density change as a repository-wide artifact migration,
not a per-chapter override. Recapture at the selected density and never upscale
an existing bitmap. Compare replacement captures with the old image at the same
rendered size before accepting them.
