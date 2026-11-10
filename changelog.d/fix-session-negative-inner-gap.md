- Fixed: loading a Session in the web app whose Circular track has an inner
  or outer gap below zero (a hand-edited file) no longer stops with an
  unrecognized error; the gap reads as Auto and the Load dialog lists it, as
  it does for a gap that is not a number. The CLI and the Python API still
  refuse such a gap (OV-3237).
