- Fixed: In the Web app, a Session whose protein comparison has no saved
  results (for example one written with the Python API before any comparison
  ran) no longer answers a label edit, Reset record positions, Align at
  feature, or a record rotation before the first Generate with "Install LOSAT
  or NCBI BLAST+". The diagram is not redrawn, and a notice says that the
  Session has no saved comparison results and asks for one Generate, which
  runs the comparison in the browser; after it, edits redraw as before.
