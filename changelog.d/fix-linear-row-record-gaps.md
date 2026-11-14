- Fixed: A Linear row of many records, such as the 97 contigs of a draft assembly, draws on both the command line
  and the Web app instead of failing with "Linear row cannot fit a positive bp scale", and a row of about 70 records
  no longer spends most of the figure's width on gaps. The default record gap (24 px) now narrows in a row whose
  gaps would take more than a quarter of its width; at the default width of 2000 px, rows of up to 21 records keep
  24 px. When the gaps or the marks between the windows of one row need more width than the figure has, the Web app
  now says so instead of naming only a Python exception.
