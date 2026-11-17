- Fixed: the Web app reads the row of a Record Order entry (`#1@2`, the number
  after the last `@`) as the command line does: a whole number above zero, with
  an optional `+`, leading zeros and spaces around it. A hand-edited Session
  whose entry has a row such as `1.0`, `0x1` or `1e0` was drawn on the Web with
  that row as 1, while `gbdraw render` refused it; the Web now loads such a
  drawing without its saved diagram and names the entry's place, where the
  command line refuses it. A one-circle drawing, which the command line draws
  without reading its entries, loads as before. An older Linear Session
  (render request schema 2 to 7) reads and refuses its entries as the command
  line does, so a row of `+2` there is row 2 on both. Sessions that gbdraw
  writes are unchanged (OV-3313).
- Fixed: Web app: when Load cannot read the saved diagram of a drawing, its
  dialog now says that the drawing keeps its settings and the input files they
  name, instead of saying that the inputs are left out (OV-3472).
- Fixed: Web app: a Linear Session written by the command line with a Record
  Order (`--multi_record_position`) places each record in the row of the entry
  that names it (`#2@1` puts record 2 in row 1), as `gbdraw render` draws it;
  before, the Web gave each record the row of the entry at its own position.
  Load now also refuses, as the command line does, a Linear drawing whose
  records have grid rows on some records only, grid rows together with a
  Record Order, or a Record Order that names a record twice, names a record
  the drawing does not have, or lists another number of records: the drawing
  loads without its saved diagram (OV-3473).
- Fixed: Web app: Save Session keeps a drawing whose only change is a color
  table or another card file, or its input type; before, such a drawing was
  left out of the file (OV-3701).
- Fixed: Web app: the Load dialog fits a phone screen when it names a long
  place in the file; the place may now wrap after a `.` or `[` (OV-3761).
