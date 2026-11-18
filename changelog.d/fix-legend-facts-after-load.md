- Fixed: Web app: after you load a Session saved while a deleted Legend row's color rule had been removed or
  recaptioned (or its feature given another Legend name), **Restore** no longer lists that row, and the next
  **Generate Diagram** no longer fails with "The generated result could not be accepted" (OV-2412). The first
  Restore after a Session Load now updates the diagram before it returns the rows, except for a Session saved by
  gbdraw versions that stored no render request, which keeps its Restore as before until the next Generate.
