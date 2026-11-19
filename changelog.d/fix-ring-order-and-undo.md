- Fixed: Web app: with **Use custom stack** on, an Undo back past a Track-stack edit made before a
  Generate that ran the ring searches no longer stops each later **Generate Diagram**; the Undo and Redo
  of that Generate give back its rings and their stack rows together, also after a Save and Load, and keep the
  Track stack and layout edits made while that Generate ran (OV-2445).
- Fixed: Web app: while **Generate Diagram** runs, the comparison-ring controls (the ring source,
  **Add comparison sequence**, the ring files, **Move series**, **Remove series**, and the colour
  re-assignment) and the moves of a ring's Track stack row refuse with the run's message instead of changing the
  rings under the run (OV-2317).
- Fixed: Web app: a Linear protein comparison (LOSATP) that the Web app of another LOSAT build searched, such
  as gbdraw 0.14.5's, is searched again at the first **Generate Diagram** after a Session load, as it is when
  any setting changes, instead of being drawn from the stored result with nothing changed.
- Fixed: Web app: with **Use custom stack** on, the comparison rings have one order: **Move series** moves the
  ring's Track stack row with it, a move of a ring's Track stack row moves the ring in the list, and the figure
  draws the rings in that order; each move is one Undo step. **Move series** never moves a ring across the Axis,
  as the stack's own arrows do not. A Session whose stack drew its rings in another order than its list opens with
  the list in the drawn order (OV-3852).
