- Changed: in a Web project with several drawings, **Clear LOSAT results**
  clears only the LOSAT results of the drawing you are looking at, and its
  count ("Clear 3 LOSAT results") is the number of results that drawing lists.
  The other drawings keep their results and their **Save Raw LOSAT TSV**. A
  search another drawing also uses stays saved, so the cleared drawing's next
  Generate reuses it instead of running LOSAT again; the result is the same.
  In a project with one drawing, Clear removes every LOSAT result, as before.
