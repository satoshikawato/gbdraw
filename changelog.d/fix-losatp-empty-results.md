- Fixed: `gbdraw linear --losat losatp --losatp_mode pairwise --save_session`
  no longer stops with "computed no pairwise result" when every record shares
  one row, so there is no record pair to compare; the Session stores the empty
  result and draws the same figure again (OV-3511).
- Fixed: a LOSATP record pair without a match is reported as a protein
  comparison with no protein match to draw, also when the comparison is read
  back from a Session or redrawn in the Web app after Generate; it was reported
  as a comparison table with no rows (OV-3512).
