- Fixed: two Linear File cards whose files have the same bytes stay two File
  cards after Save Session and Load, and Generate reads them as two files, as
  before Save; before, Load showed them as one File with every record,
  Replace, Remove and the LOSAT search of one card then acted on both, and
  Generate read both cards from one file. A Session saved before this fix
  still loads its rows as before.
- Fixed: a Linear Result of two files with the same bytes, saved after a
  Generate or written by the CLI (two such `--gbk` files, or one file given
  twice), loads matched to both File cards, also after the next Generate;
  before, the second card's record matched no record of the Result, so its
  record name edits waited for Generate.
