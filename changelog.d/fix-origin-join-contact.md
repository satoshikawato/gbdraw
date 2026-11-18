- Fixed: A feature whose two parts reach both ends of a circular record, such as
  `join(1..300,2701..3000)` on a 3,000 bp record or the same gene written as two
  sorted GFF3 rows, is now drawn as one block over the origin from a display
  start too. A Circular diagram with `--display_start_coordinate` drew a
  connector line all round the ring between the two parts.
- Fixed: A join over the origin that leaves base 1 out, such as
  `join(2901..3000,2..100)`, now shows its 1-base gap and connector line in the
  whole Circular view as it does from a display start, and its label stands on
  its longest part. The whole view joined the two parts over base 1.
- Fixed: With `--resolve_overlaps`, a feature that lies under the connector line
  of a join written in order over both ends of a record
  (`join(1..300,2701..3000)` with a feature at 1201..1800) now takes its own
  lane in Linear diagrams, and in Circular diagrams where the join is drawn as
  written. It shared the join's lane and sat on its connector line, where it
  read as a third part of the join.
- Changed: In a Circular diagram of a linear record, a join of parts at both
  ends of the record is drawn as written, with its connector line, since a
  linear record has no origin to join over; and two parts on different strands
  at both ends of a circular record are drawn as two blocks, each with its own
  strand, in the whole view as from a display start. From a display start, the
  label of a two-part join that covers both ends of a circular record with
  parts that touch or overlap each other stands on the middle of its longer
  part.
