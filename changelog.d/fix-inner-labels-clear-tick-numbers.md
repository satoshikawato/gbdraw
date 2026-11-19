- Fixed: on a Circular figure whose scale stands between the features and the
  inner labels (for example `--separate_strands --track_type spreadout
  --labels both`), inner feature labels are no longer drawn over the scale's
  numbers. The inner labels move just inside the numbers when that adds no
  overlap, crossing leader or region-label warning, also when that puts them
  on GC content, GC skew, depth, conservation ring or annotation paint;
  otherwise they are placed as before.
- Changed: Circular `--label_protection` has a new default, `auto` (web app
  **Label Protection** **Auto**): a white box, as `box` draws it, behind each
  inner label whose letters lie on a data track (GC content, GC or AT skew,
  depth, feature density, a conservation ring or an annotation track), so the
  letters stay readable there. No other label gets a box, the boxes move
  no label or track, and a box leaves a hole where a region label or a scale
  number lies under it. `--label_protection none` draws no box. A Session or
  request saved before this release, which does not state the protection,
  replays with `none` as it was drawn.
