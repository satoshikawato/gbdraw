- Fixed: Web app: Reset Label Edits after renaming every label of one text now redraws a label the editor could not
  tie to a feature, so it shows its original text instead of the renamed one (OV-3880).
- Fixed: Web app: the original text of a label you renamed (Undo needs it) is now in the saved Session whether or
  not you clicked a label first, and whether or not the shown figure was written again: Generate writes it, and Load
  writes it for each renamed label it can tie to its feature. A Session saved by an earlier version loads without the
  label editor's marks; for a renamed label that Load cannot tie to its feature, Generate writes it (OV-3848).
