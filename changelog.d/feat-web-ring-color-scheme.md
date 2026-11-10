- Added: **Ring colors** in the Web's Circular **Pairwise Comparisons** card,
  the Web side of `--conservation_color_scheme`: **Same color** (the default),
  **10 colors**, or **One color per 5 rings**. A color you set on a ring wins on
  every scheme and stays when you switch; **Use the scheme color** next to such
  a ring drops it. Scheme colors stay with their rings when you reorder them;
  **Re-assign colors when tracks are reordered** (off by default) and
  **Re-assign colors in the current order** give them out again in the drawn
  order, and never change a color you set; under **Same color** both are
  unavailable and a reorder gives out no colors. A ring added after a removal takes
  the first scheme color no ring uses. Every ring color in a Session saved by
  gbdraw 0.14.5 or earlier is a color of its own, so those Sessions open
  unchanged (PD-OI-110 scenario revision 3).
- Fixed: a CLI or Python Session saved with `--conservation_color_scheme` and
  rings without their own colors opened in the Web with every ring in one
  blue (`#4e79a7`), not the colors the CLI drew, and Save wrote that blue as
  each ring's own color. The Web now opens such a Session with that scheme under **Ring
  colors**, draws the CLI's colors, and saves the rings without colors of their
  own (OV-3054). A Session that the gbdraw 0.14.4 or 0.14.5 Web saved from such
  a CLI Session opens with that scheme too. When the Circular request it saved
  names no ring colors, each blue ring draws the scheme color that request draws
  it with, as the CLI does, also when the rings were moved before the Save; when
  that request names ring colors (that Web generated the figure before the
  Save), the Session holds no Circular request, or the request's track stack
  cannot be read, the rings keep the colors they were saved with.
