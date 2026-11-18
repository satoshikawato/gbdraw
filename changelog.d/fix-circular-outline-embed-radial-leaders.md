- Fixed: Circular labels inside a feature: a label is drawn inside its feature exactly when its letters fit
  inside the feature's drawn shape, clear of its outline. A label near the end of its feature no longer runs
  into the arrowhead or touches its outline: it moves along its feature until its letters clear the arrowhead, and
  goes outside only when no position fits them (for example ndhD in the tobacco plastome of the Gallery), and a
  label whose letters fit is drawn inside even when its laid-out box with its spacing is longer than the
  feature (for example ND4L, and "NADH dehydrogenase subunit 4" in the human mitochondrial genome).
- Fixed: Circular radial labels: a leader now ends on the edge of its feature's outline that faces the label, at the
  point nearest the label, and its slanted part runs outside the feature band, so features no longer hide it.
  In a crowded cluster the inner labels keep enough room for their leaders to rise clear of the band; when the
  circle does not have that room, the figure grows instead of refusing the labels.
- Fixed: Circular labels drawn from `--display_start_coordinate` (or a cropped record) kept clear of feature
  blocks read where the record's coordinates put them, not where the blocks are drawn.
- Fixed: Circular radial labels: when the features of another lane cover the whole edge of a labelled feature that
  faces its label, the label moves along the circle just far enough for its leader to pass them, instead of
  running under them (for example NitaCp079 in the rotated tobacco plastome of the command-line reference).
