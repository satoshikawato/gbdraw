- Fixed: the **File name** field of a pair of Linear records now shows the
  raw LOSAT file name that Generate and **Save Raw LOSAT TSV** give, built
  from the records' IDs when the cards have no Definition (for example
  `NC_000913.3.AP009048.1.losatn.tsv`); before, it showed
  `seq_1.seq_2.losatn.tsv`. After Load, opening the field reads the
  records, as the record list does, and it then shows the name; **Save Raw
  LOSAT TSV** of a pair saved without a name reads them too. A card that
  draws every record of a file of several records (a loaded Session's card
  with no record chosen) is named `seq_<n>`, its number among the drawn
  cards, where Generate used its first record's ID (OV-3227).
