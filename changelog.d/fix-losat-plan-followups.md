- Fixed: A TLOSATX genetic code that is not one of the NCBI tables LOSAT reads (1–6, 9–16, 21–33) is refused before
  any search, on the command line (`--losat_gencode`, `--conservation_losat_gencode`, and the `losat_gencode` table
  columns, whose help now lists the tables) and in the Web app (a Linear card's LOSAT Gencode, a comparison ring's
  Comparison gencode, the Circular Reference gencode), with a message that names the field and the tables. Show in
  Comparison opens the field (in Linear, the card's LOSAT Gencode), and before Generate the card lists the tables
  LOSAT reads beside such a code. LOSATN ignores a TLOSATX Gencode left on a card, and a Circular ring code that
  LOSATN would not use is reported as not applying before its table is checked. Before, the Web app's LOSAT stopped unexpectedly with
  such a code (OV-3326).
- Fixed: Downloading the LOSAT table of a window pair whose card pair has a typed file name no longer reads every
  record card first.
