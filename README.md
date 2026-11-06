# gbdraw test data

This tag holds no gbdraw source code and is not a gbdraw release. Its release
assets are the large input files of gbdraw's test suite.

In a gbdraw checkout, `tests/test_data.json` lists each file with its path,
size and sha256, and `python tools/test_data.py fetch` puts the files at those
paths. Asset names start with the file's sha256, so a changed file is a new
asset and older commits keep finding theirs.
