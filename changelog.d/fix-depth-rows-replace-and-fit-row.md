- Fixed: When a Circular track row's typed width is wider than the room inside the axis (a bare number is a multiple
  of the radius, so `w=8` asks for eight radii), the command line and the Web app now name that row and its width in
  the error, instead of the next row with an Auto width (OV-2314).
- Fixed: In the Web app, replacing the file of a Circular File card with a file that has fewer records now asks before
  it removes the Depth files of the records after the new file's last one, also after a Session load and in every
  other drawing a Replace of the file reaches: the dialog names the files (with their drawing), Replace file is one
  Undo step, and Cancel changes nothing. A Replace whose drawing is not shown when its new file has been read ends
  without writing, and a notice asks you to pick the file again; picking another file for the card, or an Undo that
  changes the card, ends it without a notice (OV-2318).
