- Fixed: Web app: after opening a saved Session whose drawing does not use a custom track stack,
  Track Preset, Add feature density track and Use custom stack work as on a new drawing: Use custom
  stack builds the stack from the Track Preset chosen after opening, and Add feature density track no
  longer changes the hidden stack. Before, the opened drawing kept a stack placed for its old preset,
  so the figure after Use custom stack could draw the features on the wrong side of the rings. The
  figure of an opened Session does not change.
- Fixed: Web app: in a Circular drawing whose custom track stack is off, Use custom stack after Hide
  GC Content or a Depth file now draws the features where the Track Preset draws them (outside the
  Axis for Spreadout). Before, the stack placed every row inside the Axis.
- Fixed: Web app: Undo of choosing Record layout "One circle" now restores the drawing exactly as it
  was, without keeping a note of the settings One circle replaced.
- Fixed: Web app: a Legend row or feature stroke edit on a figure part drawn without a stroke is
  saved and opened again unchanged.
