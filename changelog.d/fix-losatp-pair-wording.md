- Fixed: after Generate in the web app, and when a saved Session is drawn
  again, a LOSATP record pair without a match is reported as "Protein
  comparison of ... has no protein match to draw.", as `gbdraw linear` reports
  it, instead of as a comparison table with no rows (OV-3512).
- Changed: a Session now saves a BLAST table given in `--comparisons_table`,
  or as a Python `LinearComparison`, as a nucleotide comparison table under its
  file name, so drawing the Session again names the file as the first drawing
  did. A Session saved by gbdraw 0.14.0 to 0.14.5 still reads such a table as
  a nucleotide table.
