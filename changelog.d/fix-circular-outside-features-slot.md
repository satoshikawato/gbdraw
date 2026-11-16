- Fixed: a second Circular features slot (`--circular_track_slot` with the
  `features` renderer, or a second Features row in the Web app's custom track
  stack) now draws its features at its own radii; it drew them at the first
  slot's radii, one copy over the other. Outer labels name the features of the
  outermost features slot and inner labels those of the innermost, so no label
  reaches its feature across the other copy. The room radial labels keep
  around a scale now comes from those same lanes, not from the first features
  slot; with a scale between two features slots, a large label font may still
  need a smaller size to fit.
- Fixed: when the Circular features lanes of the outer labels lie outside the
  axis circle (a features slot with `side=outside` or
  `lane_direction=outside`, `--track_type spreadout`, and radial labels with
  `--labels both`), outer labels keep
  clear of the feature lane: each leader leaves its feature at the feature's
  outer edge, so it no longer runs at a shallow angle under the neighbouring
  features and seems to name one of them, and the label text stays off the
  features.
- Fixed: a Circular leader that starts at its feature's edge (the case above,
  `--resolve_overlaps`, and a feature placement that moves a feature to another
  lane) now starts on the feature's drawn outline; beside an arrowhead, which
  narrows toward its tip, it started up to 5 px short of the head, in the white
  space past it, and seemed to name nothing.
