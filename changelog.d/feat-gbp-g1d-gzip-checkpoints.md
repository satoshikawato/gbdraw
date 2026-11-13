- Performance: Large genomes draw several times faster and with less memory.
  A record's sequence is read where it lies in its file: crops, reverse
  complements and origin crops no longer copy it, and its bases are counted
  once, as the file is parsed, for the letter check, the GC % in the
  definition and the GC content and skew tracks. A gzip input is read
  without writing an expanded copy to disk. The whole human genome
  (T2T-CHM13v2.0, 3.1 Gbp) draws in the CLI in about 8 s from its FASTA and
  11 s from its `.fna.gz`, against about 65 s before; every figure is the same.
- Changed: `crop_genbank` reads its input as the drawing does. It writes
  feature keys in their INSDC spelling (`cds` as `CDS`), reads EMBL files,
  tab-indented features and spaces in a LOCUS name, and refuses an empty or
  unreadable file with the reader's message.
- Changed: an input file that changes while gbdraw reads it twice is refused
  with a message naming it, for plain files as already for gzip files.
- Fixed: a FASTA entry holding a letter of two UTF-8 bytes, such as `é`, is
  refused as not a nucleotide sequence, naming the letter, instead of
  stopping with an uncaught decoding error (OV-3333).
- Fixed: a reverse-complemented crop of an RNA record that was edited is
  saved in the Session as its source and window, so the replay's GC tracks
  read the flanks again (OV-3332).
