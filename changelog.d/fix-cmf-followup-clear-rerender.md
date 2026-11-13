- Changed: the Web's **Clear LOSAT results** button now says the count is the
  project's ("Clear the project's 5 LOSAT results"). It clears the LOSAT
  results of every drawing in the project, so its count can be larger than the
  list of results the shown drawing's panel shows.
- Fixed: in a Linear drawing of the Web app, a rerender of the shown figure
  (rotating a record to a feature, Align, Align at a feature) no longer writes
  the comparisons of the last Generate into the drawing's settings; the next
  Generate uses the comparisons the settings held, as a Circular drawing already
  does for its comparison rings (OV-3314).
- Fixed: a Session 31–39 whose LOSAT comparison rings were moved after the last
  Generate opened with the genetic codes of the rings swapped; each ring now
  opens with its own code (OV-3179).
