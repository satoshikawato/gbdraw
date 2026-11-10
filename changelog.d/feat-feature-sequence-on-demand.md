- Changed: a Session no longer stores the nucleotide and amino acid sequence of
  each feature; the Web app reads it from the Session's input files when you
  open a feature's Sequence tab, copy a match span or a Similarity group member,
  search the Nucleotide or Amino acid field, or export an Interactive SVG. The
  Gallery Sessions are 17.5 % smaller in all (78.7 to 64.9 MiB), and the
  figures and exported files stay the same. The first such read after a Load
  starts Python, which takes a few seconds on large Sessions; a Session that
  holds only its Results says that its sequences are read after the next
  Generate. Sessions saved before this change load as before.
