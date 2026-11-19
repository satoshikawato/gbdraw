- Fixed: `gbdraw linear --losat losatp --losatp_mode pairwise` without listed pairs, and the Python API, now take
  the upper-row record as the query of each adjacent-row pair when a later record stands in an upper row, as
  LOSATN, TLOSATX and the Web app do. Before, the record listed first was the query, so such a layout drew slightly
  different protein matches on the command line than in the Web app.
