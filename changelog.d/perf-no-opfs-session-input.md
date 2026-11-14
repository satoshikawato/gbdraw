- Performance: In a browser that lets the Web app write no temporary input copies, Generate after
  loading a Session file no longer reads each input whole into the page and hashes it again: the
  Load hashes each input while it checks it, and Generate streams the input from the Session file.
  With a 220 MiB depth table, the page's peak memory during Generate fell by about a quarter.
