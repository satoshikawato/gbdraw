- Fixed: the Web app reads the row of a Record Order entry (`#1@2`, the number
  after the last `@`) as the command line does: a whole number above zero, with
  an optional `+`, leading zeros and spaces around it. A hand-edited Session
  whose entry has a row such as `1.0`, `0x1` or `1e0` was drawn on the Web with
  that row as 1, while `gbdraw render` refused it; the Web now loads such a
  drawing without its saved diagram and names the entry's place, and a drawing
  saved with its input files keeps them and its settings. An older Linear
  Session (render request schema 2 to 5) reads each row the same way on both,
  so the command line now reads a hand-edited `+2` there as row 2. Sessions
  that gbdraw writes are unchanged (OV-3313).
