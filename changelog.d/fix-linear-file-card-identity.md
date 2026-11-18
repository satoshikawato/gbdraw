- Fixed: two Linear File cards whose files have the same bytes stay two File
  cards after Save Session and Load, and Generate reads them as two files, as
  before Save; before, Load showed them as one File with every record,
  Replace, Remove and the LOSAT search of one card then acted on both, and
  Generate read both cards from one file. A Session saved before this fix
  still loads its rows as before.
