- Fixed: An Interactive SVG exported from the Web app shows every text that
  contains "@" (a qualifier value such as an e-mail address in a note, a match
  popup row, a search entry); its viewer took every such text for an internal
  protein ID and left it out, while a file exported with the command line or
  Python showed it.
