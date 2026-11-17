- Fixed: In the Web app, a number field you clear stays empty after Undo,
  Redo, Save Session and Load. Before, a cleared LOSAT Gencode on a Linear
  file card came back as 1 (and lost its "Gencode required" mark), a cleared
  Depth series tick came back as the drawing's tick after Undo, and a Session
  could hold a cleared field in two different forms. Sessions saved by earlier
  versions load with the same figures as before.
- Fixed: When a drawing has its own Depth tick (set by a configuration or a
  command-line Session), an empty Large Tick, Small Tick or Tick Font Size
  field of a Depth series now shows that tick as its Auto value, which is the
  tick the figure draws. Before, it showed "endpoints only", "None" or the
  track's own size.
