- Fixed: a Linear figure with a GC content track in deviation mode (the default) kept 0.4 px of empty room above
  and below that track for a band border that only percent mode draws, so its rows were up to 0.8 px taller than
  their tracks; the room is now kept only where the border is drawn.
- Fixed: a GC content, GC skew or Depth track whose outline colour is `none` kept room for an outline width it does
  not draw, in its Linear row and at the figure's edge in both Circular and Linear figures, so a wide invisible
  outline made the figure larger; it no longer does.
- Fixed: in a Linear figure, a Depth track outline wider than the axis line was cut off at the figure's right edge;
  the figure now includes it.
