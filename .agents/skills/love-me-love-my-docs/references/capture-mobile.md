# Reproducible mobile capture

Use this path only when the user explicitly requests a mobile manual.

Prefer committed Maestro flows with stable accessibility identifiers and
`takeScreenshot`. Pin the device or emulator model, OS version, locale, theme,
font scale, orientation, and seeded demo state. Install tools only from an
auditable package manager or official release; never use `curl | bash`.

Each flow starts from a declared reset state, performs the documented taps and
text entry, asserts the resulting screen, and saves screenshots in step order.
Do not include notifications, device-owner data, accounts, tokens, or other
private content.

If no simulator or emulator is available, commit the flow and an exact manual
capture checklist, mark every missing image with a visible TODO, and report the
degraded state. Do not silently substitute Web screenshots.
