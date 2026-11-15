- Fixed: a Linear row that its record gaps could not hold was told to "lower Record gap or leave it blank (Auto)",
  with no value, and with the Auto gap even when that gap was the one refused. The Web now asks for a smaller typed
  Record gap, such as 0, and names the Auto gaps when they are the cause. The command line and Python name the set
  gap below which every row of the figure draws (and the Auto gap in px when it was refused), and suggest Auto only
  when every row draws with it. A row that its window marks fill is named first, as no gap helps it.
