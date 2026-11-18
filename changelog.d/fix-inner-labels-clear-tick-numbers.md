- Fixed: on a Circular figure whose scale stands between the features and the
  inner labels (for example `--separate_strands --track_type spreadout
  --labels both`), inner feature labels are no longer drawn over the scale's
  numbers. The inner labels move just inside the numbers when that adds no
  overlap or crossing leader and puts no label onto GC content, GC skew,
  depth, conservation ring or annotation paint it was not already on;
  otherwise they are placed as before. When such a track sits directly inside
  the scale, there is no room between them, and labels can still be drawn over
  the numbers.
