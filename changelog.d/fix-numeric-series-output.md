- Fixed: a Linear figure with a GC content track in deviation mode (the default) kept 0.4 px of empty room above
  and below that track for a band border that only percent mode draws, so its rows were up to 0.8 px taller than
  their tracks; the room is now kept only where the border is drawn.
- Fixed: a GC content, GC skew or Depth track whose outline colour is `none` kept room for an outline width it does
  not draw, in its Linear row, and a GC content or GC skew track also at the figure's edge (and a Depth track at
  the edge of a Circular figure), so a wide invisible outline made the figure larger; it no longer does.
- Fixed: in a Linear figure, a Depth track outline wider than the other lines of its record ran into the figure's
  right margin, and a wide one past the figure's edge; the figure now keeps its full margin beyond it.
- Fixed: in a Linear figure, a GC content or GC skew outline kept its room above and below the record a second
  time, and kept room on records that do not draw that track; the figure now keeps that room once, where the track
  is drawn.
