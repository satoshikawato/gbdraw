- Fixed: Renaming a feature type's Legend row in the Web app, such as `CDS`
  renamed to `Coding`, now changes only the row's name, as `--legend_rows`
  does on the command line. The features keep their default color, the
  Specific color rules table stays as it was, and **Reset** returns the name.
  This holds from the Legend editor and the popup, on a row of one feature,
  and for a type's "other" row such as `other proteins`. The name is the
  type's: while a color rule with a Legend name takes some of the type's
  features, the rest are drawn as `Coding` in the type's "other" row, and the
  name stays the same when the rule goes, whichever of the two rows was
  renamed. Typing the row's default name gives the default name back. Before,
  the rename added one color rule per product to the table. Merging into
  another row, picking a row from the popup's list, typing the name of a color
  rule group in the row's color, "This feature only" on a row of several
  features, and renaming a row that color rules or a per-feature color edit
  draw still write rules.
  A Legend name of a row that a Generate does not draw waits for the row,
  also when a File is replaced; before, replacing a File dropped it.
