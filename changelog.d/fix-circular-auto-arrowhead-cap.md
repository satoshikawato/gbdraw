- Fixed: on a short Circular record the automatic arrowhead was a long, thin
  spike: it spanned at least 30 bp, which on a 240 bp record is 45 degrees of
  the circle, about eight times the feature's thickness. The automatic head is
  now at most half the feature's thickness, the proportion a 5-10 kb record
  already had, so short records get short arrow tips. At default sizes, records
  above about 6 kb keep their heads, except those of 75-85 kb, whose heads are
  up to 0.3 px shorter. A record drawn with thinner features than usual, such
  as a mitochondrion beside a large genome on a multi-record canvas, also gets
  slightly shorter heads. An explicit `--arrow_head_length_ratio` is
  unchanged.
