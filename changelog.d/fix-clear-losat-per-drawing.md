- Changed: in a Web project with several drawings, **Clear LOSAT results**
  clears only the LOSAT results of the drawing you are looking at, and its
  count ("Clear 3 LOSAT results") is the number of results that drawing lists.
  The other drawings keep their results and their **Save Raw LOSAT TSV**. A
  search another drawing still holds is reused by the cleared drawing's next
  Generate instead of running LOSAT again; the result is the same. In a project
  with one drawing, Clear removes every LOSAT result, as before.
- Changed: **Clear LOSAT results** is a step in the drawing's History. **Undo**
  brings back the results it cleared, and **Redo** clears them again.
- Fixed: each drawing keeps its own LOSAT results through Delete, Duplicate,
  Undo and Redo. A deleted drawing that **Undo** brings back has the results it
  had, even after a Clear in another drawing; a Linear comparison edit in one
  drawing no longer removes another drawing's results; and **Save Session**
  writes into each drawing only the LOSAT results that drawing holds.
- Changed: after **Undo** of a Generate, the next Generate runs the LOSAT
  searches the undone Generate had added again; it no longer takes them from
  the undone step.
