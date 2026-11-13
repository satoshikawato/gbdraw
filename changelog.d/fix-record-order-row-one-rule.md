- Fixed: the Web app reads the row of a Record Order entry (`#1@2`, the number
  after the last `@`) as the command line does: a whole number above zero, with
  an optional `+`, leading zeros and spaces around it. A hand-edited Session
  whose entry has a row such as `1.0`, `0x1` or `1e0` was drawn on the Web with
  that row as 1, while `gbdraw render` refused it; the Web now loads such a
  Session as settings only, as it does for a row of `0`, and names the entry's
  place. Sessions that gbdraw writes are unchanged (OV-3313).
