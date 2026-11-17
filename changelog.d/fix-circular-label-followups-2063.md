- Fixed: Circular radial labels: a crowded group of labels spreads to both sides of its features in every
  figure. Before, the group was pushed one way unless its leaders crossed, so a one-pixel change of a leader
  end could turn a whole group by several degrees and slant its leaders (for example the rotated tobacco
  plastome in the command-line reference).
- Fixed: Circular labels inside a feature: whether a label fits inside its feature is measured at the radius
  the label is drawn at. A features lane inside the axis circle (the inner copy of two features slots, or
  the minus-strand lane with `--separate_strands`) no longer draws a label longer than its feature inside
  it; the label goes outside instead.
