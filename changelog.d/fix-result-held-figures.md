- Fixed: Web app: with Auto Reflow on, the Undo of a feature hide or a Label
  Off now shows the figure from before the edit at once, instead of waiting
  for the labels to be placed again; its Redo shows the edited figure at once
  too.
- Fixed: Web app: in a duplicated drawing, Undo, Redo and Legend moves no
  longer change the drawing it was copied from, and edits in that drawing no
  longer change the copy.
- Changed: Web app: on a figure drawn to a requested size, the Undo of a
  Legend rename now shows the fitted figure from before the rename at once,
  without fitting the figure again.
- Performance: Web app: an Undo or Redo back to Label visibility On shows the
  label at once, from the figure drawn for that step, without drawing it
  again.
