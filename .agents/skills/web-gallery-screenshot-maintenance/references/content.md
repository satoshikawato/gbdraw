# Gallery tutorial content

Read for captions, instructions, tables, or tutorial structure changes. Apply
the rules to the requested scope; a local edit does not require a full audit.

## Writing Rules

Keep tutorial text aligned with the actual workflow.

- Use user-facing labels in tutorial text. Avoid implementation/UI-pattern jargon such as `segmented control`; write the action the reader should take, such as `Select Circular`.
- Omit generic Requirements sections when they only restate that the reader needs a browser or the listed input file.
- Do not add a Files tab step unless checking Gallery artifacts is part of the workflow. Input filenames can usually live in `downloads` or the upload step.
- Do not add Color rules or Post-generation edits sections when the example does not require color-rule setup or real post-generation editing.
- Prefer screenshots over dense setup lists when the UI can show the values clearly. Keep any remaining text to the values the reader must type or choose.
- Do not describe automatically generated output as a manual editor task. If `Generate Diagram` creates legend entries, tracks, labels, or previews, say they are generated and reserve drawers/editors for review or optional tweaks.
- Use action text for operations and state text for generated results. Avoid vague instructions such as `keep visible` when the UI already produced the state.
- Omit operation `title` and `body` when a single media item already sits under a clear step title/body and the extra operation text would only create a redundant bold subheading.
- Keep pre-generation setup, generated-result inspection, and post-generation edits as separate concepts.
- When a step already has a `table`, do not repeat the same fields as slash-delimited operation text or captions. Refer to the table row or category instead.

Write captions and alt text as action/state descriptions:

- Good: `Select Circular.`
- Good: `Set Dinucleotide to GC, Window to 500, and Step to 50.`
- Good: `Enter CDS / product / tyrosine recombinase in the SPECIFIC RULES (-t) row.`
- Bad: `The rendered web app crop shows ...`
- Bad: `Confirm ...`
- Bad: `Web app session view.`

When replacing a screenshot, update the caption and alt text in the same change unless they are already concrete and accurate.

## Structured Content Rules

Use structured display when tutorial content has repeated fields.

- Prefer a table over slash-delimited bullets for color rules, track-slot recipes, file mappings, record metadata, or any list where each row shares the same fields.
- Convert true repeated-field content within the requested tutorial inventory; keep simple one-column checklists as bullets. Report patterns noticed outside that scope without expanding a local edit into a sibling-tutorial audit.
- Use clear column headers that match the UI or data model, such as `Feature`, `Qualifier`, `Value pattern`, `Color`, and `Legend caption`.
- Keep caveats such as regex specificity in a separate `note`; do not hide them inside a dense row.
- Use the renderer's existing structured fields, such as `table: { columns, rows }`, rather than embedding table markup in JSON strings. If support is missing, change the renderer only when implementation is within the request; otherwise report that limitation.
- When renderer/layout work is in scope, use compact CSS with horizontal overflow for narrow viewports and update focused browser tests for headers and representative rows.
- Do not allow table text to break inside short file names, accessions, extensions, or other atomic tokens such as `BGC0000708.gbk`. If a table is too narrow, preserve the token and rely on the table wrapper's horizontal overflow instead of using character-level wrapping.
- Keep existing `items` bullets for simple one-column checklists where a table would add noise.
