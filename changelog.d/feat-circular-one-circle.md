- Added: `gbdraw circular --one_circle` draws all records of the input one
  after another on one circle (the chromosomes of a eukaryote, the contigs of a
  draft assembly), each on its own arc with its own ticks numbered from its
  start, its GC content in percent, and a band that names it (chromosome,
  plasmid or organelle). `--one_circle_gap_degrees` sets the gap between
  records (default 1.0), `--one_circle_min_length` leaves shorter records off
  the circle and counts them in the centre text, also when the definition is
  hidden (default 1000 bp), and
  `--one_circle_no_record_band` leaves the band out. GC skew is off unless
  `--skew` is given. Depth, conservation and region tracks are refused on one
  circle. The Python API (`CircularMultiRecordOptions(record_layout=
  "one_circle")`) and the render request (schema 10) carry the same layout.
  The centre text names an organism only when every drawn record shares it.
  In the Web app, **One circle** sits beside **Multi-Record Canvas** with the
  same three settings; turning it off returns to the layout, GC content mode
  and Hide GC Skew it left. A one-circle Session from the command line opens
  and Generates the same figure. Feature density tracks are drawn per record on
  one circle (#200).
