- Fixed: `gbdraw linear --losat losatp --losatp_mode pairwise` draws a comparisons table's `source=table` BLAST rows
  beside the record pairs its `source=losat` rows search, as the Web draws that plan; before, it stopped with
  "Pass either losat_search pairs or linear_comparisons, not both". A pair that both kinds of row name stops the run
  with a reason. The Web's Source recipe of such a drawing is now available
  ([web app reference](./docs/REFERENCE/web-app.md)).
  The Python API likewise refuses a given comparison of an adjacent pair beside a LOSATP Pairwise search that lists no
  pairs, instead of merging it into the protein matches.
