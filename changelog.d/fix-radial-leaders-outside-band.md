- Fixed: Circular radial labels: a leader now ends on the edge of its feature's outline that faces the label, at the
  point nearest the label, and its slanted part runs outside the feature band, so features no longer hide it.
  In a crowded cluster the inner labels keep enough room for their leaders to rise clear of the band; when the
  circle does not have that room, the figure grows instead of refusing the labels.
- Fixed: Circular labels drawn from `--display_start_coordinate` (or a cropped record) kept clear of feature
  blocks read where the record's coordinates put them, not where the blocks are drawn.
