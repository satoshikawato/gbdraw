- Fixed: A captioned Specific color rule for every feature type (`*`) now draws
  its Legend row, read as the same rule written for each type; before, every
  feature took its color while the Legend kept the types' palette rows.
- Fixed: In the Web app, the feature popup, its row-scope dialogs, and the
  Legend editor put a feature that no captioned rule colors in the row the
  figure draws for it: its type's `other <feature type>s` row while a
  captioned rule of its type is used, as the figure shows, not its type's row.
  The Legend item popup of that row now counts its features, and its Color
  sets the type's default color, so the features change with the swatch.
- Fixed: A Legend no longer shows `other <feature type>s` for a feature type
  that a captioned rule colors none of. A rule now counts as used only for the
  types whose features take it, so such a type keeps its own row, both for a
  table that writes one caption and color for several types and for a rule
  for every type (`*`).
