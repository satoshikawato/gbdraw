- Fixed: in the Web app, the drawing that Find in other files creates draws
  with its first Generate when a File draws one record in several windows;
  before, Generate stopped and asked to replace or clear the comparison
  (OV-3365).
- Fixed: in the Web app, a (−) window that covers a whole record, as Find in
  other files proposes it, is drawn reverse-complemented; before, it was drawn
  forward when the drawing it came from showed the same record forward
  (OV-3376).
- Fixed: in the Web app, Input Genomes no longer says that record rotation,
  feature placement, or tolerance has changes pending Generate right after a
  Generate of a drawing that Find in other files created (OV-3375), and
  Similarity alignment works in such a drawing and in the drawing it came from
  (OV-3377). In the drawing it came from, the feature popup's record start and
  orientation actions also work again on a File drawn whole.
- Fixed: in the Web app, Find in other files no longer drops a Search started
  right after an edit of Start or End and a tick of Across the origin; before,
  the dialog went back to the search without a proposal or a message.
