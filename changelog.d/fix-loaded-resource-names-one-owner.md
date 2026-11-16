- Fixed: `gbdraw render --session` now names a Session's files as the Web app
  does when it opens that Session. For an older Session that stored no file
  names, the command line printed names such as `resource-0001-record-1.gbk`
  and the Web app printed `record-1.gbk`. The second of two files saved with
  the same name is no longer called `record.2.gb`.
- Fixed: in the Web app, saving or generating again after opening a Session
  keeps the name of each comparison's BLAST table. When two comparisons of the
  Session used one table, or a saved table kept its stored name, the new Session
  could name a table by its internal resource name or lose the name.
