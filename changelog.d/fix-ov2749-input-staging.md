- Fixed: in the Web app, Generate passes each input file to the drawing Worker as a
  file, no longer as a copy of its bytes in the message, which Chromium refuses
  for 2,047 MiB or more. Such an input no longer stops Generate with "Input
  resource preparation failed. Reselect the input and retry." (OV-2749). When
  the page runs out of memory while it prepares the inputs, the error now says
  so and offers Retry, Save Session, and reloading the page, instead of asking
  to reselect the input. The Undo history limit no longer counts chosen input
  files, plain or gzip, whose bytes Undo does not keep, so choosing an input
  larger than about 200 MiB no longer discards the older Undo steps (OV-2750).
