- Fixed: In the Web app, switching drawings, including to the copy a **Duplicate** shows, shows the arriving
  figure fitted to the preview. Before, a wide Linear figure could sit half past the right edge, and while a
  **Review alignment** panel was open it could lie entirely out of view (OV-3620, OV-3621).
- Fixed: In the Web app, while an alignment review resolves or applies, **Undo**, **Redo** and the History list's
  rows are disabled, with the reason, when the step they would take creates, copies or deletes a drawing. Before,
  they looked available and a press did nothing.
- Fixed: In the Web app, while alignment work runs, **Retry Generate**, the label fit dialog's **Apply and
  Generate**, and the Generate choices of the missing BLAST TSV dialog are disabled with the reason, as **Generate
  Diagram** is. Before, they started nothing without saying why, and the two dialog choices still changed their
  setting.
- Fixed: In the Web app, a Session Load that fails after it began replacing the drawings leaves an open alignment
  review open on the restored drawing; before, the review closed.
