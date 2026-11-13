- Fixed: In the Web app, a number field you clear stays empty after Undo,
  Redo, Save Session and Load. Before, a cleared LOSAT Gencode on a Linear
  file card came back as 1 (and lost its "Gencode required" mark), a cleared
  Depth series tick came back as the drawing's tick after Undo, and a Session
  could hold a cleared field in two different forms. Sessions saved by earlier
  versions load with the same fields and figures as before; a cleared Depth
  series tick in such a Session now stays empty instead of showing the
  drawing's tick.
