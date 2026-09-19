# Reproducible Web capture

Use Playwright against the requested application; apply the external-mutation
boundary in `SKILL.md` when the workflow writes remote state. Keep the base
URL, viewport, device scale factor, locale, timezone, theme, and network policy
in one configuration owner.

## Flow contract

1. Start from a fresh browser context. A restore scenario may reload only the
   session that its own earlier steps created from original source inputs; do
   not preload a bundled Gallery session or finished project.
2. Reach the UI through normal navigation.
3. Upload files and change values through visible, accessible controls.
4. Prefer `get_by_role`, then `get_by_label` for a form control with an
   associated accessible label, then a stable `get_by_test_id`. A label locator
   is compliant and does not imply that the app needs a test ID. Record any
   inaccessible interactive control as an application finding. Stable CSS
   selectors for SVG geometry or capture regions are acceptable; they are not
   by themselves evidence of an accessibility defect.
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
published size. For a requested local resolution change, use a supported scoped
capture setting when available. If changing a shared capture configuration would
affect other outputs, identify that affected contract before broad regeneration;
a one-image request does not authorize a repository-wide artifact migration.
Recapture at the selected density and never upscale an existing bitmap. Compare replacement captures with the old image at the same
rendered size before accepting them.
