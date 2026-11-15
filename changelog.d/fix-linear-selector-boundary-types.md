- Fixed: replaying an older Web Linear Session with `gbdraw linear --session`
  crops the record of the card that set the crop or window. When an earlier
  card drew every record of a file of several records, the command line named
  the cropped record by card position, so the CLI cropped a record of the
  earlier file and drew the intended record whole. The region now names the
  first record of its card's file (`<file>:#1:<start>-<end>`); Sessions whose
  cards each hold one record draw as before. The command line the replay
  records names two cards of the same file name apart (`HmmtDNA.gbk`,
  `HmmtDNA.2.gbk`), so rerunning it crops the same card.
