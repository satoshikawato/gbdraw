- Fixed: After you load a Linear Session saved by an older Web app (Session 31 to 44) whose LOSATP
  comparison was stored as a search with its saved LOSAT rows, Reset record positions, Align at feature and
  label edits before the first Generate redraw the figure from those rows. They stopped before with "The
  comparison search could not run or returned unusable output", because the browser has no LOSAT to search
  with. A search whose rows the Session did not save still runs at Generate.
