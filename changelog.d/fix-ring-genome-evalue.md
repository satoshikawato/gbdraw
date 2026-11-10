- Changed: LOSATN and TLOSATX E-values of Circular similarity rings use one
  scale for each ring, in the Web app, the CLI, and the Python API alike
  (OV-3319). The reference genome is the files of the drawn records and the
  ring's genome is the files of that ring, with every record of those files
  counted, drawn or not; a GFF3 and FASTA pair counts every FASTA sequence. A
  crop, the records drawn, or a genome split over several files no longer
  changes a ring's E-values. Rings whose genome and reference are each one
  file of one record are unchanged; other LOSAT rings search again at their
  next Generate, and hits near the E-value threshold may appear or disappear.
  With LOSAT v0.1.0, which has no `-searchsp`, the CLI keeps the previous
  scale and says so in a warning.
