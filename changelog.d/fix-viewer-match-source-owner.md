- Fixed: An exported Interactive SVG now gives each match span the same
  sequence, or the same reason, as the Web app. Before, a homology band that
  names no comparison source said its source identity was invalid instead of
  giving its span, a reference named by an alias that two records share
  copied one of them instead of saying the alias is ambiguous, and spans past
  the record's end, on a record the file lacks, with no record named or with
  an index that is not a whole number gave other reasons than the app
  (OV-3542).
- Fixed: In an exported Interactive SVG, the coordinates of a match span, the
  rows of a collinear block or similarity group, and their notes no longer
  touch the border of their box (OV-3543).
