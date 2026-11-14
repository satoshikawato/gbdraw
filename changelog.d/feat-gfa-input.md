- Added: a GFA 1 assembly graph as input, in both diagram modes: `--gfa` in
  the command line, a `gfa` column in `--records_table`, `read_gfa()` and
  `GfaInputSource` in Python, and the **GFA** source choice in the Web app.
  Each segment (S line) is a record with its sequence and no features, in file
  order, named by the segment name. A segment with a link from its end back to
  its own start is circular and drops the last n bases of an `nM` overlap, so
  each base of the circle is drawn once; in the Web app its **Circular record**
  box starts on. Other links, and path and walk lines, are not drawn, and one
  notice per file says so. A GFA 2 file, a segment without its sequence (`*`)
  and a self-loop overlap gbdraw cannot drop are refused with the reason.
