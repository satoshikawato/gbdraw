- Changed: a Session no longer saves, with each Legend row, a copy of the
  color Python drew for it or of the Legend color you set; it keeps a row's
  color only for a row that Python does not draw. Sessions saved before
  show the same Legend colors after Load.
- Fixed: a Legend row that you gave a Legend color of your own keeps
  showing that color in the Legend panel after you recolor the color rule
  that draws it, as the figure draws it.
