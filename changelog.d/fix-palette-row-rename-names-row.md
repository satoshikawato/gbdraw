- Fixed: Renaming a feature type's Legend row in the Web app, such as `CDS`
  renamed to `Coding`, now changes only the row's name, as `--legend_rows`
  does on the command line. The features keep their default color, the
  Specific color rules table stays as it was, and **Reset** returns the name.
  Before, the rename added one color rule per product to the table.
