- Fixed: Renaming a feature type's Legend row in the Web app, such as `CDS`
  renamed to `Coding`, now changes only the row's name, as `--legend_rows`
  does on the command line. The features keep their default color, the
  Specific color rules table stays as it was, and **Reset** returns the name.
  This holds from the Legend editor and the popup, on a row of one feature,
  and for a type's "other" row such as `other proteins`. The name follows the
  type's row: while a color rule with a Legend name takes some of the type's
  features, the rest are drawn as `Coding` in the type's "other" row, unless
  that row has a name of its own. Before, the rename added one color rule per
  product to the table. Merging into another row, picking a row from the
  popup's list, typing the name of a color rule group in the row's color,
  "This feature only" on a row of several features, and renaming a row that
  color rules or a per-feature color edit draw still write rules.
