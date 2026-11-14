- Fixed: `gbdraw render --session` and `gbdraw circular|linear --session` now
  log the wrap line of a wrapped Linear figure ("Wrapping at ... per row"), as
  `gbdraw linear` does, and a Session with several drawings logs each
  drawing's warnings and input summary right after that drawing's figures.
  Before, the wrap line was missing, all drawings' lines came after the last
  figure, and a drawing that failed to render kept the drawings before it
  from logging theirs.
