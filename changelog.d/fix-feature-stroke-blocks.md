- Fixed: In the Web app, a stroke set on a feature or on its Legend row no
  longer closes the open end of a partial feature (for example `<101..300`):
  the outline takes the new stroke and the partial end stays open, as gbdraw
  draws it. The connector lines between the blocks keep their own line stroke.
