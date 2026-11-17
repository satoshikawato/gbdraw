- Changed: in a Web project with several drawings, **Clear LOSAT results**
  clears only the LOSAT results of the drawing you are looking at, and its
  count ("Clear 3 LOSAT results") is the number of results that drawing lists.
  The other drawings keep their results and their **Save Raw LOSAT TSV**. A
  search another drawing's diagram uses is reused by the cleared drawing's next
  Generate instead of running LOSAT again; the result is the same. A search
  that another drawing only keeps from an earlier Generate is not reused, so
  LOSAT runs it again; this also holds after **Save Session** and Load for a
  file this version saved. In a project with one drawing, Clear removes every
  LOSAT result, as before.
- Changed: **Clear LOSAT results** is a step in the drawing's History. **Undo**
  brings back the results it cleared, and **Redo** clears them again. It is
  unavailable while a diagram is being generated, as **Reset Settings** is;
  before, a Generate that finished after a Clear brought the cleared results
  back.
- Fixed: each drawing keeps its own LOSAT results through Delete, Duplicate,
  Undo and Redo. A deleted drawing that **Undo** brings back has the results it
  had, even after a Clear in another drawing; a Linear comparison edit in one
  drawing no longer removes another drawing's results; and **Save Session**
  writes into each drawing only the LOSAT results that drawing holds.
- Changed: **Save Session** writes into each drawing only the LOSAT searches
  its diagram uses, not one it only keeps from an earlier Generate under other
  settings. After the file is opened again, a Generate under those earlier
  settings runs LOSAT once more; the result is the same. A file an earlier
  version saved may still hold such searches; after Load, each drawing counts
  them as used until its next Generate.
- Fixed: after **Clear LOSAT results** in a Linear drawing with a LOSATP
  comparison, a Generate with nothing changed ran no LOSAT search and drew the
  comparison from the stored results, so the drawing listed no results. It now
  runs LOSAT again, or reuses a search another drawing's diagram uses. A file
  saved after the Clear and opened again still draws the stored comparison
  without a search, as a file saved without LOSAT searches does.
- Changed: after **Undo** of a Generate, the next Generate runs the LOSAT
  searches the undone Generate had added again; it no longer takes them from
  the undone step.
