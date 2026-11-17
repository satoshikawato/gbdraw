- Fixed: Web app: after Duplicate, Align at feature and the popup's record start in the copy change the copy's
  records, not its source's.
- Fixed: Web app: a Session saved after an input file changed on disk keeps the file as it was chosen while a drawn
  figure still uses it; when the file cannot be read, Save names it.
- Fixed: Web app: a Session Load that fails and is rolled back keeps the previous drawings' input files, CLI options
  and Preview searches, and Undo and Redo of a Sequence search made before it still redraw the figure at once.
- Fixed: Web app: a new drawing made after a drawing was deleted no longer inherits the deleted drawing's Preview
  search or the CLI options of a loaded Session, and its first Sequence search set is named "Search 1".
