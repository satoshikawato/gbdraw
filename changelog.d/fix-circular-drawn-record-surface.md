- Fixed: In the Web app, the Circular status line under the file card counts
  the records the next Generate draws. With one record chosen in the Record
  select, a two-record file now reads "1 of 2 records drawn" instead of "2 of
  2". The Topology rows list the same records, so a Multi-Record Canvas lists
  every record even while the Record select keeps a chosen one.
- Fixed: In the Web app, the disabled Display start field of a cropped record
  shows where the drawn crop starts (the Circular 12 o'clock, the Linear left
  end): 3000 for 500..3000 drawn reverse-complemented, 500 drawn forward. It
  showed an end of the whole record before.
