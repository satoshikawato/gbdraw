- Fixed: `gbdraw linear` and `gbdraw render` list the proteins of a record or
  window drawn reversed in the order the Web app searches them (the record
  before its reverse complement), so a LOSATP comparison whose side holds
  several records of one file, one of them drawn reversed, reuses the searches
  a Web Session saved instead of reporting them missing from its LOSAT cache.
