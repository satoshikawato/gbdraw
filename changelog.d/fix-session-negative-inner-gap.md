- Fixed: loading a Session in the web app whose Circular track has an inner
  or outer gap below zero (a hand-edited file) no longer stops with an
  unrecognized error; the gap reads as Auto and the Load dialog lists it, as
  it does for a gap that is not a number. A number too large to hold (such as
  `1e400`) in any saved number setting reads as its default the same way. The
  CLI and the Python API still refuse both, and now report an infinite
  annotation coordinate or lane as a Session error rather than an internal
  error (OV-3237, OV-3390).
