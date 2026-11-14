- Fixed: A Linear LOSATP comparison (every mode) now draws when some records
  have no CDS protein, such as the short contigs of a draft assembly, on the
  command line and in the web app. Such a record draws without protein links,
  and a record turned off in Choose records adds nothing to the searched
  database, so no E-value changes. Generate stops only when no compared pair
  of records has CDS proteins in both, and then names a record without them.
- Changed: In Collinear LOSATP comparisons, every displayed pair of records
  now keeps a comparison, empty when the pair has no collinear block (such as
  a pair with a record that has no CDS protein), so the figure reports that
  the pair has no protein match to draw.
- Fixed: When a Session saved by gbdraw 0.11 to 0.13 holds two stored copies
  of one LOSATP search (one the web app wrote and one the command line wrote),
  the web app and the command line now drop the copy they do not use once the
  search is answered, with the older derived results kept for it, instead of
  carrying it in every saved Session.
