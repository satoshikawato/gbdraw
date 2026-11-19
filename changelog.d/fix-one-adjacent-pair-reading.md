- Fixed: `gbdraw linear --losat losatp --losatp_mode pairwise` without listed pairs, and the Python API, now search
  each record of a row against each record of the next row when a later record stands in an upper row, as LOSATN,
  TLOSATX and the Web app do. Before, the record listed first was the query, so such a layout drew slightly
  different protein matches on the command line than in the Web app.
