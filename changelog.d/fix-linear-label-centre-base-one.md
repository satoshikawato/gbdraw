- Fixed: A feature that starts at the record's first base is drawn from that
  base in Circular and Linear diagrams, instead of one base short at its left
  end, and the line between two parts of a feature reaches both parts instead
  of stopping one base short at each end. Lanes and labels read the same
  bounds, so a feature that shares only the first base with another takes
  its own lane.
