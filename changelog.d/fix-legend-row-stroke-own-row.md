- Fixed: In the Web app, a stroke on a Legend row now outlines only the features drawn in that row, not every feature
  drawn in the row's color. When two rows share a color (two color rules of one color, or a rule in a feature type's
  color), a stroke on one row no longer outlines the other row's features. A Session saved with such a stroke opens as
  saved; the next Generate, or the next stroke or color edit, outlines only that row's features.
