- Changed: a Session is one `.gbdraw` file (a ZIP file) that holds each input as
  it is instead of as text inside one JSON file, so a Session no longer stops at
  about 384 MiB of inputs. The Web app, `save_session_document`,
  `--session_output` and `--save_session` write only this file, and refuse a
  name that ends with `.json` or `.json.gz`; every reader opens it and the older
  JSON Sessions by content. gbdraw 0.14.5 and older cannot open it: Session 51
  files open in v0.15.0 and later. Convert an older Session with
  `gbdraw render --session old.json.gz --session_output new.gbdraw`.
- Changed: in the Python API, `SessionDocument.to_dict()` no longer holds
  input bytes; they stay in `document.resources`. To load an edited mapping,
  pass `load_session_document(mapping, resources=document.resources)` (the
  error without them says so), or save and load the `.gbdraw` file.
