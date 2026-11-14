- Fixed: `gbdraw render --session` and `gbdraw circular|linear --session` now
  log the same warnings and input summary as the `gbdraw circular` or
  `gbdraw linear` command that saved the Session: region annotations clipped
  to a record or skipped because a feature selector matched nothing, the line
  for each drawn record and comparison table, and the warnings for a record or
  comparison table that draws nothing. Before, they logged only the notices
  about saved feature edits.
