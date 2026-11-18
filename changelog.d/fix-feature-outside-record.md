- Fixed: A feature that does not fit in its record (one that starts after the
  record ends, or is longer than the record) is now left out of the figure with
  a warning that names it, in Circular and Linear diagrams, on the command line,
  in the Python API and in the Web app. `gbdraw circular` stopped with "A circle
  arc spans at most 360 degrees" on such a feature, and `gbdraw linear` drew it
  past the end of the record without a word. Features that fit are drawn as
  before, including those that cross a circular record's origin.
