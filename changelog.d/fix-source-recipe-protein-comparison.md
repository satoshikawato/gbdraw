- Fixed: Run info in the Web app gives the CLI rerun command for a Linear drawing with a protein comparison (LOSATP
  Pairwise, Similarity groups, or Collinear), and for a Session that `gbdraw linear` wrote with `-b` or
  `--comparisons_table` once it is generated in the Web, where it said that the source recipe was unavailable. The
  command runs the protein search with the settings the Web used; for Pairwise it lists the compared record pairs as
  `source=losat` rows of its comparisons table. A drawing with a BLAST table beside a protein comparison still has no
  rerun command, with a reason in words.
