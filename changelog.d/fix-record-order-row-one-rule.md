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
