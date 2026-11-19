- Fixed: In the Web app, replacing the file of a Circular File card with a file that holds the same records in another
  order no longer draws one record's Depth files and feature edits on another record, also when the Session is saved
  and loaded, or was loaded, before the next Generate, and when the file is replaced more than once before it: each
  record keeps them with the record of the new file with its accession (its record ID without the version) or its
  record ID, and by its place in the file only when record IDs repeat or are missing, counted from the records the
  last Generate drew. A record the new file does not hold keeps neither: the question before the Replace names its
  Depth files, also when one record is replaced by a record with another ID (OV-3846). When those are the only files
  of a Depth series, the question also says that the series stays empty and that Generate will ask for a file or its
  removal (OV-3954). The question lists a record's Depth files in the shown language's list form.
- Fixed: In the Web app, the unit selector of an empty (Auto) Circular track Width now starts at px, the unit of the
  Auto width shown under it, so a number typed over the shown width means pixels; an empty Radius still starts at ×R
  (OV-2314).
- Fixed: In the Web app, the messages of an Align on a feature that matches no member of its Similarity group, or
  more than one, and of an Align on a group the shown diagram no longer has name the button to select instead: Align…
  in the Similarity groups tab, after choosing the exact reference record and feature.
