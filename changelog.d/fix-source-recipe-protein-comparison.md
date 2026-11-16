- Fixed: Run info in the Web app gives the CLI rerun command for a Linear drawing with a protein comparison (LOSATP
  Pairwise, Similarity groups, or Collinear) and for a Session that `gbdraw linear` wrote with `-b` or
  `--comparisons_table`, where it said that the source recipe was unavailable. The command runs the protein search
  with the settings the Web used and lists selected record pairs as `source=losat` rows of its comparisons table.
