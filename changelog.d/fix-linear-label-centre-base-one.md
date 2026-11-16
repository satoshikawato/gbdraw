- Fixed: A feature that starts at the record's first base is drawn from that
  base in Circular and Linear diagrams, instead of one base short at its left
  end, and the line between two parts of a feature reaches both parts instead
  of stopping one base short at each end. The line is drawn under the
  blocks, so each block's outline stays whole where the line meets it. Lanes
  and labels read the same bounds, so a feature that shares only the first
  base with another takes its own lane.
- Fixed: The Web app reads a multi-part feature's block stroke from one of
  its blocks, not from its connecting line, for the Result's block stroke
  that a Session keeps, a new Legend row's swatch outline, and the stroke the
  feature popup shows after a stroke edit is reset.
- Fixed: On a circular record, a feature whose parts cross the origin with a
  gap there draws its connecting line over the origin in every view, not only
  from a display start, and keeps one lane with the features it does not
  overlap. A feature that covers the whole record is drawn as a full ring in
  Circular diagrams.
