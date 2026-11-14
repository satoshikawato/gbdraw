- Fixed: In the Web app, each Linear drawing keeps its own Similarity
  alignment and record positions. Before, drawings shared one: aligning a copy
  showed the alignment on the drawing it was copied from and Generate drew that
  drawing aligned too, Undo in one drawing cleared the other's, a Session saved
  with one drawing aligned loaded with every drawing aligned, adding a file to
  a new drawing cleared another drawing's alignment, and a record moved in a
  copy moved in the drawing it was copied from. Duplicate copies the alignment;
  New drawing, New drawing with this style, and a Circular-only Session start
  without one.
