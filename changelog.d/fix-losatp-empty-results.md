- Fixed: `gbdraw linear --losat losatp --losatp_mode pairwise --save_session`
  no longer stops with "computed no pairwise result" when every record shares
  one row, so there is no record pair to compare; the Session stores the empty
  result and draws the same figure again (OV-3511).
- Fixed: `gbdraw render --session` and `gbdraw linear|circular --session` log
  what the replay read and drew (the per-record and per-comparison lines and
  their warnings), as a fresh render does; before, a replay logged none of
  them (OV-3528).
