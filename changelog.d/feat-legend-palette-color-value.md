- Added: a Specific-colors `color` can be `palette` (any letter case): the
  rule draws its features in their type's color as the figure is drawn, so a
  palette or Default colors change reaches them. Without a caption, its
  features are in their type's legend row. See
  [Palettes and feature rules](./docs/REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md).
- Fixed: in the Web app, a feature given its palette color back with
  **Reset fill color** kept the old color after a palette or Default colors
  change, and its legend row split in two; the reset now writes a `palette`
  rule without a caption, which the rule table shows with the tag
  **Palette color**. A Session saved by 0.14.5 reads its reset rules so.
- Fixed: a rule a user named `other rRNAs` (or `other proteins`) in the
  type's palette color was drawn in the type's own legend row; it draws its
  own row again, as in 0.14.5.
