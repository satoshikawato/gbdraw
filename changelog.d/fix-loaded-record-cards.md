- Fixed: In the Web app, a loaded Session whose Linear card draws every record of a multi-record file, with a saved
  comparison that names those records one by one, now draws after you choose Replace with current controls: Generate
  splits the card into one card per record, as it does for a new comparison, so LOSAT searches every record. Before,
  LOSAT searched only the first record and Generate failed (OV-3324).
- Fixed: In the Web app, a Linear card that names its record by number (`#1`) shows that record in its Record select
  instead of "#1 (not found in current file)" with a warning, as Generate already drew it (OV-3325).
