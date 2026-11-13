- Fixed: A Linear label drawn inside its feature no longer runs over the
  feature's outline. Its letters, with 0.125 em more at each end, must lie
  inside the outline the feature is drawn with: the shaft (lower when the
  shaft is narrowed) and, in the arrowhead, a height that falls to 0 at the
  tip. A label that fits at the middle stays there; one that does not moves
  sideways to the nearest place where it fits; one that fits nowhere is drawn
  outside, with a leader (with Label Placement inside only, it is not drawn).
  A label also keeps off the outline of a feature drawn after its own on the
  same track, which used to cross its letters. Labels drawn outside before
  stay outside.
