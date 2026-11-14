- Fixed: after you load a Session saved by an older gbdraw (Session version
  before 40) with an active Similarity alignment, moving a record by drag,
  starting or reviewing an alignment, and Align at feature now work before the
  first Generate. Before, they stopped on a record the older figure draws
  reverse complemented, because its features are drawn under other IDs than
  the alignment names. Feature rows such a Session saved with a conflicting
  record number from the second input on are now read again from the source
  files at Load, so their Similarity group members show their source strand.
