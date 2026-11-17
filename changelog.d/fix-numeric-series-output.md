- Fixed: a Linear figure with a GC content track in deviation mode (the default) kept 0.4 px of empty room above
  and below that track for a band border that only percent mode draws, so each record row was 0.8 px taller than
  its tracks; the room is now kept only where the border is drawn. Likewise a Linear GC content, GC skew or Depth
  track whose outline colour is `none` no longer keeps room for an outline width it does not draw.
