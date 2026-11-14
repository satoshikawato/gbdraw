- Fixed: after Generate in the web app, and when a saved Session is drawn
  again, a LOSATP record pair without a match is reported as "Protein
  comparison of ... has no protein match to draw.", as `gbdraw linear` reports
  it, instead of as a comparison table with no rows (OV-3512).
- Changed: a Session now saves a BLAST table given in `--comparisons_table`,
  or as a Python `LinearComparison`, as a nucleotide comparison table under its
  file name, so drawing the Session again names the file as the first drawing
  did. A Session saved by gbdraw 0.14.0 to 0.14.5 still reads such a table as
  a nucleotide table.
- Fixed: when `gbdraw render --session` or the Python API draws a saved
  Session again, a comparison is named by the file the user gave, as the
  first drawing named it: a LOSAT nucleotide comparison reads "Comparison of
  record #1 ... and record #2 ..." instead of naming a table file gbdraw made,
  an empty `--comparisons_table` table keeps its file name, and a BLAST table
  that the web app saved is named without the web app's resource prefix
  (OV-3557).
