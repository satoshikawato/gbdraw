- Fixed: A Linear protein comparison (LOSATP collinear, pairwise, or Similarity groups) over a circular
  record whose CDS is written as GFF3 writes a span over the origin (one part past the record's end,
  `68680..69018` on a 68,848 bp record) no longer stops with "Comparison source feature index conflicts
  with its view feature ID". The CDS reads as the join over the origin it is, as a GenBank file writes
  it, so its matches cover its two parts at both ends of the record, in the Python API, the CLI and the
  Web app.
