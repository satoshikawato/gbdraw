- Added: `save_session_document(path, document=...)` writes a Session
  document that is built already, such as the `document` of
  `paste_region_drawing()` or of `upgrade_session_document()`, without
  building it again. See [Typed requests](./docs/REFERENCE/typed-requests.md).
- Changed: in the Web app, **New drawing from a region…** is now **Copy to
  another drawing…**, and it also opens from a record's popup and from the
  records a Shift+drag selects ("N records: <names>" in the selection
  toolbar). The dialog chooses the new drawing's mode and whether it keeps
  the drawing's settings and edits or starts from the mode's defaults;
  several whole Circular records make one Multi-Record Canvas. It also adds
  records or regions below the records of an existing drawing, with or
  without the edits of those records, says when the drawing has one of them
  already (**Skip it** or **Add it again**), and is one Undo step of that
  drawing. See [Web app](./docs/REFERENCE/web-app.md).
- Fixed: the Web app drew every record of a file when a Session written by
  `gbdraw circular --records_table --multi_record_canvas` drew only some
  records of that file on one canvas; it now draws the records the Session
  names.
- Fixed: the Web app lost the rows (`--multi_record_position`) of such a
  Session's records, and of the records pasted into a Circular canvas, and
  drew them on its own grid; it now keeps the rows the request gives.
- Fixed: **Copy to another drawing…** with the drawing's settings kept
  drew the new drawing without the plot title's position (and so, on a
  Circular canvas, without the title); a copy into the same mode now keeps
  the title's text, size and position, and a copy into the other mode keeps
  its text and size.
- Changed: `derive_region_drawing()` of several whole Circular records
  keeps each record on the row the source canvas draws it on, without the
  rows that hold none of them, instead of placing the records in order; and
  `paste_region_drawing()` into a Circular canvas that names no places keeps
  its records where they are drawn and adds the pasted records on a new row
  below them.
- Fixed: in a Linear row of region windows, Shift+drag over the row's list
  of windows in the Web app selected no record; it now selects the records of
  those windows. Each list line of the SVG names its record
  (`data-gbdraw-record-index`, `data-gbdraw-record-id`); the figure looks the
  same.
- Fixed: in the Web app, a Save that fails because the diagram runtime cannot
  start no longer tells you to save a Session; it offers Retry and Reload.
