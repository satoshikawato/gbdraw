- Fixed: A captioned Specific color rule for every feature type (`*`) now draws
  its Legend row, read as the same rule written for each type; before, every
  feature took its color while the Legend kept the types' palette rows.
- Fixed: In the Web app, the feature popup, its row-scope dialogs, and the
  Legend editor put a feature that no captioned rule colors in the row the
  figure draws for it: its type's `other <feature type>s` row while a
  captioned rule of its type is used, as the figure shows, not its type's row.
