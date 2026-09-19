# Sequence input provenance

Read before adding or changing sequence inputs, mirrored fixtures, or visible
results in public procedural documentation.

- For sequence-based tools, make the reader obtain original sequence records
  from an authoritative public database by accession. A public procedural page
  must not use a repository-bundled sequence file or a prebuilt project/session
  as its reader-facing input. Give the database page or API, format choice,
  exact save name, and identity checks. A page may reload only a session the
  reader created earlier in that same page
  from the original inputs.
- Do not assume that an accession's nucleotide version also freezes its
  feature table. When rendered features or protein searches depend on an
  annotation revision, link an official database revision-history snapshot or
  another authoritative annotation release that readers can download, and
  record that exact request in the evidence.
- A frozen repository copy may support deterministic offline automation only
  after a mirror-verification record captures the authoritative database URL
  or exact API request, versioned accession, requested format, UTC retrieval
  date, source byte size and SHA-256, mirror byte size and SHA-256, and the
  direct comparison result. Byte-preserving copies require equal hashes. If a
  documented deterministic normalization is necessary, retain the original
  source hash and executable transformation, then compare the mirror hash with
  a freshly derived output. Accession, length, topology, or parser checks are
  useful additional guards but do not replace the byte comparison.
- Mark a mirror without that evidence `legacy-unverified`. Do not use it to
  regenerate public screenshots or other visible results until it is fetched
  again from the authoritative source and compared or rebuilt. If the exact
  source cannot be recovered, report the blocked evidence path or adopt a new
  authoritative input and regenerate the documentation; never grandfather the
  mirror. Once verified, routine offline runs need only check the pinned mirror
  hash and identity metadata and must still keep the mirror out of the reader's
  acquisition path.
- Distinguish files readers download, create, generate, and compare. Show the
  complete contents or reproducible derivation for support tables and other
  non-sequence inputs.
