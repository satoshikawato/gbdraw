- Fixed: Web app: with **Use custom stack** on, an Undo back past a Track-stack edit made before a
  Generate that ran the ring searches no longer stops each later **Generate Diagram**; the Undo and Redo
  of that Generate give back its rings and their stack rows together, also after a Save and Load (OV-2445).
- Fixed: Web app: while **Generate Diagram** runs, the comparison-ring controls (the ring source,
  **Add comparison sequence**, the ring files, **Move series**, **Remove series**, and the colour
  re-assignment) refuse with the run's message instead of changing the rings under the run (OV-2317).
- Fixed: Web app: a Linear protein comparison (LOSATP) that the Web app of another LOSAT build searched, such
  as gbdraw 0.14.5's, is searched again at the first **Generate Diagram** after a Session load, as it is when
  any setting changes, instead of being drawn from the stored result with nothing changed.
