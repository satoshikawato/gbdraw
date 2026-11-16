- Fixed: A feature that starts at the record's first base is drawn from that
  base in Circular and Linear diagrams, instead of one base short at its left
  end, and the line between two parts of a feature reaches both parts instead
  of stopping one base short at each end. The line is drawn under the
  blocks, so each block's outline stays whole where the line meets it. Lanes
  and labels read the same bounds, so a feature that shares only the first
  base with another takes its own lane.
- Fixed: The Web app reads a feature's block stroke where it is drawn: on a
  multi-part feature from one of its blocks, not from its connecting line,
  and on a partial feature from its outline instead of "none". This applies
  to the feature popup's stroke, the stroke it shows after a stroke edit is
  reset, a new Legend row's swatch outline, and the Result's block stroke
  that a Session keeps.
- Fixed: On a circular record, a feature whose parts cross the origin with a
  gap there draws its connecting line over the origin in every view, not only
  from a display start, and keeps one lane with the features it does not
  overlap. A feature that covers the whole record is drawn as a full ring in
  Circular diagrams.
