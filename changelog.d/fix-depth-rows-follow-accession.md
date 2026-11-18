- Fixed: In the Web app, replacing the file of a Circular File card with a file that holds the same records in another
  order no longer draws one record's Depth files and feature edits on another record, also when the Session is saved
  and loaded, or was loaded, before the next Generate: each record keeps them by its record ID when every record ID of
  both files is unique, and by its place in the file otherwise. The question before
  such a Replace names the records whose ID the new file lacks (OV-3846).
- Fixed: In the Web app, the unit selector of an empty (Auto) Circular track Width now starts at px, the unit of the
  Auto width shown under it, so a number typed over the shown width means pixels; an empty Radius still starts at ×R
  (OV-2314).
