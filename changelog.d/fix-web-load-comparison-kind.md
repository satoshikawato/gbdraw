- Fixed: the Web app loads a 0.14.0-0.14.5 `gbdraw linear --comparisons_table`
  Session's BLAST table as a BLAST comparison; it had opened it as a LOSATP
  Pairwise comparison.
- Fixed: a LOSATP Pairwise drawing whose records all share one row keeps LOSATP
  Pairwise in its Session, from `gbdraw linear` and from the Web app alike, so the
  Web app loads it with LOSATP Pairwise selected and a Generate after moving a
  record to another row runs the comparison. The figure is unchanged, and
  redrawing it still needs no LOSAT.
