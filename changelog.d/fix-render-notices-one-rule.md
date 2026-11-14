- Fixed: a warning about one input file is now given once per run, and once
  per drawing of a Session, on every route. With the same file given twice
  (for example `gbdraw linear --gbk a.gb a.gb -k tRNA`), `gbdraw linear`,
  the Python API's `draw_linear` and `draw_circular`, and the Web app's
  preview notice gave the file's warning (such as "has no feature to draw")
  twice, while `gbdraw render --session` and `gbdraw linear --session` gave
  it once; all of them now give it once, as `gbdraw circular` does.
