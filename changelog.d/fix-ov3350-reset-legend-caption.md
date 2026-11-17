- Fixed: after **Reset fill color** gives a feature its palette color back, the
  Legend keeps the feature type's own caption (for example `rRNA`) when no
  other color rule with a caption colors a drawn feature of that type, instead
  of `other rRNAs`; the live Legend, Generate, a CLI replay of the Session, and
  the feature popup's Legend name show the same row. While another rule still
  colors a feature of the type, the reset feature stays in the `other ...` row
  as before.
