- Fixed: A feature that does not fit in its record (one that starts after the
  record ends, or is longer than the record) is now left out of the figure with
  a warning that names it, in Circular and Linear diagrams, on the command line,
  in the Python API and in the Web app. `gbdraw circular` stopped with "A circle
  arc spans at most 360 degrees" on such a feature, and `gbdraw linear` drew it
  past the end of the record without a word. The Web app's feature list leaves
  it out too, so a cropped or reversed record no longer fails to list its
  features, and protein comparisons (LOSATP) search and count only the
  features that are drawn, on the command line and in the Web app.
- Fixed: A feature written past the end of a circular record, the GFF3 way to
  write a feature over the origin (`81..110` on a 100 bp record), is now drawn
  in Linear diagrams as the same feature written `join(81..100,1..10)`: in two
  blocks, at the end and at the start of the record. Linear drew it past the
  end of the axis, a reversed record drew only its first block, and a moved
  display start (`--display_start_coordinate`) stopped the run in both modes.
  On a record that is not circular such a feature does not fit and is left out
  with the warning. Circular diagrams draw it in one arc as before, also
  reversed, and color rules and feature selectors still match its location as
  written. A feature density track and a `region` crop now count and draw its
  part from the start of the record, and its protein translation and the Web
  popup's sequence read all of its bases; they stopped at the record's end.
- Fixed: A GFF3 record whose landmark (the `region` line covering the whole
  sequence) says `Is_circular=true` is now read as circular, as a GenBank
  record whose LOCUS line says `circular`: its GC content and skew windows wrap
  over the origin in Circular and Linear diagrams, which changes the GC content
  and skew tracks near both ends of a Linear diagram and their scale, and a
  window across its origin can be drawn.
