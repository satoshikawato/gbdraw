- Fixed: after Generate in the web app, and when a saved Session is drawn
  again, a LOSATP record pair without a match is reported as "Protein
  comparison of ... has no protein match to draw.", as `gbdraw linear` reports
  it, instead of as a comparison table with no rows (OV-3512).
- Changed: a Session now saves a BLAST table given in `--comparisons_table`,
  or as a Python `LinearComparison`, as a nucleotide comparison table, and
  stores with each comparison the name of the file the user gave. With
  `gbdraw render --session` and the Python API, a Session saved by gbdraw
  0.14.0 to 0.14.5 still reads such a table as a nucleotide table.
- Fixed: after Generate in the web app, and when a saved Session is drawn
  again with `gbdraw render --session` or the Python API, a comparison is named
  by the file the user gave, as the first drawing named it, spaces and all: a
  LOSAT nucleotide comparison reads "Comparison of record #1 ... and record #2
  ..." instead of naming a table file gbdraw made, an empty `--comparisons_table` table keeps
  its file name, two files of one name keep it, and a BLAST table that the web
  app saved is named without its resource prefix (OV-3557).
- Fixed: a Python comparison table whose rows are bound to features but that
  holds no protein matches is reported as a comparison table when its Session
  is drawn again, as the first drawing reported it, not as a protein
  comparison (OV-3574).
