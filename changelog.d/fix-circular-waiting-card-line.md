- Fixed: Web Circular input: while a File card still waits for a file (a GFF3
  file without its FASTA, for example), the Source records line says which
  card needs which file, such as "File 2 needs its FASTA file before the
  records can be inspected." It used to show nothing, so it was unclear why no
  record of any card was listed.
