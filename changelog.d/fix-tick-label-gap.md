- Changed: Circular scale labels drawn inside the scale (the default figure,
  `label_in_tick_out`) now stand as close to their tick marks as labels drawn
  outside: 3.2-5.9 px of ink instead of 7.1-9.6 px in the default human
  mitochondrion figure (OV-3300). Their text box keeps the label margin from
  the ticks on both sides; before, labels inside kept a quarter of their
  height more. Every default Circular figure changes: the labels move toward
  the ticks and the tracks inside them move outward by the same few pixels.
  Figures with the labels outside (`label_out_tick_in`) do not change. The
  Circular tick tiers above 10 Mbp end 1.17 % later, since the labels now sit
  on a larger circle.
