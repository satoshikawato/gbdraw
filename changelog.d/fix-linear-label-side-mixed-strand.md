- Fixed: in Linear diagrams with `--label_placement above_feature` and
  `--separate_strands`, the label of a feature whose parts lie on both strands
  (a trans-spliced gene such as plastid rps12) now sits on the same side as the
  other labels of the lane its feature is drawn in. Before, it could sit on the
  other side of its feature, with its text across the axis.
