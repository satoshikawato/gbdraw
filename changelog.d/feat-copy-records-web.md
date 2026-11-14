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
- Changed: `derive_region_drawing()` of several whole Circular records
  keeps each record on the row the source canvas draws it on, without the
  rows that hold none of them, instead of placing the records in order; and
  `paste_region_drawing()` into a Circular canvas that names no places keeps
  its records where they are drawn and adds the pasted records on a new row
  below them.
