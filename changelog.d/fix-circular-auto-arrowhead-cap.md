- Fixed: on a short Circular record the automatic arrowhead was a long, thin
  spike: it spanned at least 30 bp, which on a 240 bp record is 45 degrees of
  the circle, about eight times the feature's thickness. The automatic head is
  now at most half the feature's thickness, so short records get short arrow
  tips. A head that was already shorter is drawn as before; with the default
  track layout that is every record above about 6 kb except those of 75-85 kb,
  whose heads become up to 0.3 px shorter. With other track layouts, and on
  records drawn with thinner features than usual, such as a mitochondrion
  beside a large genome on a multi-record canvas, the cap can also shorten the
  heads of somewhat longer records. An explicit `--arrow_head_length_ratio` is
  unchanged.
