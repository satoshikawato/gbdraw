- Fixed: In the Web app, a Linear drawing that Find in other files creates, with
  one window per record, now opens with a comparison the Comparison controls can
  change, its window pairs as the Selected pairs; before, it kept the comparison
  read-only until "Replace with current controls" was clicked. Sessions that the
  Python API of 0.14.5 and earlier saved with unset Collinear block parameters
  open the same way. A drawing or Session without saved Web controls whose LOSATP
  Pairwise comparison lists pairs now opens with those pairs selected; before, it
  opened with all adjacent pairs, so Generate could compare pairs it was not
  given. Python now writes the default Collinear block parameters instead of
  `null` (OV-3382, OV-3383).
