- Fixed: Opening a Session saved by an older gbdraw (Session version 39 or
  earlier) whose Linear comparisons name a pair of records by positions that
  do not exist (a hand-edited file) no longer stops with an internal error.
  The comparison is kept without a LOSAT result file name. A position before
  the first record no longer takes the LOSAT result file name of the last
  record.
