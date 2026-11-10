- Changed: LOSATN and TLOSATX E-values in Linear diagrams use one scale for
  each pair of genomes, in the Web app, the CLI, and the Python API alike
  (OV-3143). Files whose records share a row are one genome, with every record
  of those files counted, drawn or not; a GFF3 and FASTA pair counts every
  FASTA sequence. A crop, a window, or a genome split over several files no
  longer changes a hit's E-value. Figures in which every drawn record is the
  whole of a one-record file are unchanged: LOSAT keeps its own scale for
  them, which it adjusts for the length of a hit, so their E-values differ
  slightly from those of the same genomes cropped. Other figures with LOSAT comparisons
  search again at their next Generate, and links near the E-value threshold
  may appear or disappear. LOSATP E-values (Pairwise, Similarity groups,
  Collinear) use one scale for each subject genome: the database size is
  every protein of every record of the subject genome's files, drawn or not.
  With LOSAT v0.1.0, which has no `-searchsp` or `-dbsize`, the CLI keeps the
  previous scale and says so in a warning; one scale per genome needs a newer
  LOSAT or NCBI BLAST+. The Web app's LOSAT is now built from LOSAT
  `5234c754f`, which reads these scales (`-searchsp`, `-dbsize`) in every
  execution mode, Serial included; searches without them give the same
  results as before.
