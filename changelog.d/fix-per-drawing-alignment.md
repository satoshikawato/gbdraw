- Fixed: In the Web app, each Linear drawing keeps its own Similarity
  alignment and record positions. Before, drawings shared one: aligning a copy
  showed the alignment on the drawing it was copied from and Generate drew that
  drawing aligned too, Undo in one drawing cleared the other's, a Session saved
  with one drawing aligned loaded with every drawing aligned, adding a file to
  a new drawing cleared another drawing's alignment, and a record moved in a
  copy moved in the drawing it was copied from. Duplicate copies the alignment;
  New drawing, New drawing with this style, and a Circular-only Session start
  without one.
- Fixed: In the Web app, an open Similarity alignment review stays with the
  drawing it was opened in. While another drawing is shown, the review stays
  open with its choices but does not mark them on that drawing's figure, and
  Apply waits with "Show <drawing> to apply this alignment."; before, Apply
  there dropped the review and its choices with "The diagram, source, crop, or
  Similarity Group changed." although nothing had changed. Showing the review's
  drawing again lets you Apply. A Generate of another drawing leaves the review
  open, and deleting the review's drawing ends the review as Cancel does.
- Changed: In the Web app, switching drawings waits while a Similarity
  alignment is resolved, applied, or reset.
- Changed: In the Web app, switching drawings clears the "Similarity alignment
  applied" summary, which described the drawing shown before.
