- Added: Web UI: العربية (Arabic) and فارسی (Persian) in the Language menu.
- Fixed: Ring names and the Label filtering Blacklist take the Arabic comma ، between items, in the Web UI and
  in the `--label_blacklist` option of the command line, which reads the same list. `--label_blacklist`
  now also starts a new keyword at a line break, as the Web UI does, and takes a byte order mark (U+FEFF),
  U+0085 and U+001C-U+001F off a keyword's ends, as it already did with spaces.
