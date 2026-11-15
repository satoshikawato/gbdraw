- Fixed: After you load a saved Session in the Web app and generate again,
  Run info names each record by the file you chose (for example
  `BGC0000708.gbk`), not by the file's name inside the Session
  (`record-1-genbank-BGC0000708.gbk`). `gbdraw render --session` and the
  Python Session API name a record of a Session the Web app saved the same way.
