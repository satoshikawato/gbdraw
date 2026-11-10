- Fixed: A saved Session follows a later change of its palette instead of
  keeping the palette's old colors. `--save_session` stores the `-d`/`--default_colors`
  rows you gave, not every row of the palette; a Web Session whose drawing
  changed none of its palette's colors stores no colors; and the Gallery
  Sessions keep only the Default colors that differ from their palette.
