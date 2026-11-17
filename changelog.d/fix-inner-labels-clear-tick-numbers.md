- Fixed: on a Circular figure whose scale stands between the features and the
  inner labels (for example `--separate_strands --track_type spreadout
  --labels both`), inner feature labels are no longer drawn over the scale's
  numbers. The inner labels move just inside the numbers when that adds no
  overlap or crossing leader; otherwise they are placed as before.
