- Added: `paste_region_drawing()` in `gbdraw.api` appends records or regions
  of one drawing of a materialized Session to another drawing of it. The
  drawing keeps its own settings, colors, and alignment; the records'
  per-feature color, label, and visibility edits come along under their new
  record keys (`edits="records"` brings the records only). A span the drawing already
  draws (the same record of the same file, and the same span) is listed in
  `already_present` and left out unless `duplicates="add"`. See
  [Typed requests](./docs/REFERENCE/typed-requests.md).
- Changed: `derive_region_drawing()` takes `settings="defaults"`, which draws
  the chosen records with the mode's defaults instead of the source drawing's
  settings and edits. Its mode now defaults to the source drawing's mode
  whenever every region is a whole record, so a Circular drawing of several
  whole records gives a Circular Multi-Record Canvas instead of a Linear
  drawing; a Circular drawing that would mix records and ranges is refused
  with `CopyRefused`.
