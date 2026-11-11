- Changed: When `gbdraw circular` draws only one record without `--multi_record_canvas`, its Session saves that drawing
  as a single-record drawing, as the Web app saves its own one-record drawings. In the Web app, the Record select then
  shows that record, and a records-table crop appears in the single-record crop fields, where it can be edited. The
  figure files, their names, and their content do not change. Sessions saved before this change still open as
  "All records (separate diagrams)" until you choose the record, and a typed Python API `CircularBatchRequest` of one
  record still saves as a batch.
