- Fixed: on a short Circular record (below about 6 kb at default sizes) the
  automatic arrowhead was a long, thin spike: it spanned at least 30 bp, which
  on a 240 bp record is 45 degrees of the circle, about eight times the
  feature's thickness. The automatic head is now at most half the feature's
  thickness, the proportion a 5-10 kb record already had, so short records get
  short arrow tips. Longer records draw as before, and an explicit
  `--arrow_head_length_ratio` is unchanged.
