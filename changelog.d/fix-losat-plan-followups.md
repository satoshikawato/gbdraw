- Fixed: A TLOSATX genetic code that is not one of the NCBI tables LOSAT reads (1–6, 9–16, 21–33) is refused before
  any search, on the command line (`--losat_gencode`, `--conservation_losat_gencode`, and the `losat_gencode` table
  columns) and in the Web app (a Linear card's Gencode, a comparison ring's Comparison gencode, the Circular
  Reference gencode), with a message that names the field and the tables. Before, the Web app's LOSAT stopped
  unexpectedly with such a code (OV-3326).
- Fixed: Downloading the LOSAT table of a window pair whose card pair has a typed file name no longer reads every
  record card first.
