- Fixed: a second Circular features slot (`--circular_track_slot` with the
  `features` renderer, or a second Features row in the Web app's custom track
  stack) now draws its features at its own radii; it drew them at the first
  slot's radii, one copy over the other. Outer labels name the features of the
  outermost features slot and inner labels those of the innermost, so no label
  reaches its feature across the other copy.
- Fixed: when the Circular features slot lies outside the axis circle
  (`side=outside`), outer labels keep clear of the feature lane: each leader
  leaves its feature at the feature's outer edge, so it no longer runs at a
  shallow angle under the neighbouring features and seems to name one of them,
  and the label text stays off the features.
