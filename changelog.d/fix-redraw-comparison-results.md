- Fixed: After you load a Linear Session saved by an older Web app (Session 31 to 44) whose LOSATP
  comparison was stored as a search with its saved LOSAT rows, a label edit before the first Generate redraws
  the figure from those rows. It stopped before with "The comparison search could not run or returned unusable
  output", because the browser has no LOSAT to search with; Reset record positions and Align at feature no
  longer stop with that message either. A search whose rows the Session did not save still runs at Generate.
