- Fixed: `gbdraw render --session` now names a Session's files as the Web
  app does when it opens that Session (Session version 31 and later). For an
  older Session that stored no file names, the command line printed names such
  as `resource-0001-record-1.gbk` or `record-1.gbk` where the Web app showed
  `record-1.gbk` or `SARS-CoV-1.gbk`.
  The second of two files saved with the same name is no longer called
  `record.2.gb`.
- Fixed: in the Web app, generating or saving again after opening a Session
  keeps the name each comparison's BLAST table had, also after Undo and Redo.
  When two comparisons of the Session used one table, or a table gbdraw made
  was saved with the Session, the new request could name a table by another
  comparison's name or by its file name.
