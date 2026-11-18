- Fixed: Web app: after you load a Session saved while a deleted Legend row's color rule had been removed or
  recaptioned (or its feature given another Legend name), **Restore** no longer lists that row, and the next
  **Generate Diagram** no longer fails with "The generated result could not be accepted" (OV-2412). The first
  Restore after a Session Load now updates the diagram before it returns the rows, except for a Session saved by
  gbdraw versions that stored no render request, which keeps its Restore as before until the next Generate.
- Fixed: Web app: when the diagram update that follows an edit fails or is canceled, a **Restore** of deleted Legend
  rows pressed meanwhile asks for the update once more and, if that fails too, is refused with "Updating diagram.
  Retry after the update finishes." instead of listing rows the next **Generate Diagram** may refuse (OV-2412).
- Fixed: Web app: a **Restore** of a deleted Legend row that waits for a diagram update no longer restores nothing
  when the deleted rows are rewritten during the wait (an Undo, a color-rule change).
