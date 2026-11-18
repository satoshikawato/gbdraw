- Fixed: Web app: on a figure with thousands of labels, editing or hiding a label no longer blocks the page for many
  seconds; the editor now finds the parts of every label in one pass over the figure (OV-2438).
- Fixed: Web app: an exported SVG or Interactive SVG is now the same whether or not you clicked a label first, and it
  no longer keeps the original text of a label you renamed (OV-3848).
- Fixed: Web app: Reset Label Edits after renaming every label of one text now redraws a label the editor could not
  tie to a feature, so it shows its original text instead of the renamed one (OV-3880).
