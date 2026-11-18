- Fixed: Web app: Reset Label Edits after renaming every label of one text now redraws a label the editor could not
  tie to a feature, so it shows its original text instead of the renamed one (OV-3880).
- Fixed: Web app: a saved Session no longer depends on whether you clicked a label first, or on whether the shown
  figure was written again after Generate or Load. A label you renamed keeps its original text in the Session (Undo
  needs it), from the moment Generate draws it or Load reads it; a Session saved by an earlier version is brought to
  the same form when it loads (OV-3848).
