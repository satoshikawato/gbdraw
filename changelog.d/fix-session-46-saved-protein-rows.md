- Fixed: After you load a Linear Session 46 or 47 that `gbdraw linear --losat losatp --session_output` wrote
  between 2026-10-11 and 2026-10-23, a label edit or Reset record positions before the first Generate redraws
  the figure from the LOSAT rows the Session saved. It stopped before with "The comparison search could not run
  or returned unusable output", because the browser has no LOSAT to search with.
