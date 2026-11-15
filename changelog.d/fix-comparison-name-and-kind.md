- Fixed: the Web app loads the BLAST table of a Session saved by
  `gbdraw linear --comparisons_table` in gbdraw 0.14.0 to 0.14.5 as a BLAST
  comparison table, and reports it as "Comparison of record #1 ... and record
  #2 ...", as `gbdraw render --session` does; it had opened the table as a
  LOSATP Pairwise comparison and, after Generate, reported it as a protein
  comparison. A LOSATP Session of those releases still loads as LOSATP.
- Fixed: after Generate in the Web app, Run Info names a comparison as
  `gbdraw render --session` names it: a LOSAT nucleotide comparison reads
  "Comparison of record #1 ... and record #2 ..." instead of "Comparison table
  comparison-resolved-1-comparison-resolved-1.tsv", and an uploaded BLAST table
  keeps the name of the file the user chose, spaces and all. A Session saved
  after Load keeps these names, so its comparisons are named the same when it
  is drawn again in the Web app or with `gbdraw render --session`.
- Fixed: a Python comparison table whose rows are bound to features, or whose
  BLAST columns are in another order, is reported as a comparison table, not
  as a protein comparison, when its Session is drawn again with
  `gbdraw render --session`, the Python API, or the Web app, as the first
  drawing reported it. This holds for Sessions of gbdraw 0.14.0 to 0.14.5 and
  for Sessions saved from now on, and a Session saved again in the Web app
  keeps it.
- Fixed: a LOSATP Pairwise drawing whose records all share one row keeps LOSATP
  Pairwise in its Session, from `gbdraw linear` and from the Web app alike, so
  the Web app loads it with LOSATP Pairwise selected and a Generate after
  moving a record to another row runs the comparison. The figure is unchanged,
  and redrawing it still needs no LOSAT.
- Fixed: when the Web app loads a Session whose BLAST tables its comparison
  controls cannot edit, the message says "The saved comparison tables can run,
  but the current controls cannot represent them losslessly." instead of
  calling them a protein comparison.
