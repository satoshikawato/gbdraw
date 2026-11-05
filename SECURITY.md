# Security policy

## Supported versions

| Version | Supported |
| --- | --- |
| 0.14.x | Yes |
| Older | No |

A version line appears here when its first release is tagged.

## Reporting a vulnerability

Report a vulnerability privately with GitHub Security Advisories: open
<https://github.com/satoshikawato/gbdraw/security/advisories/new> (the
"Report a vulnerability" button on the repository's Security tab). Do not open
a public issue or pull request for a vulnerability.

Include the gbdraw version, how you installed it, the input that triggers the
problem (a minimal file if possible), and what you expected to happen.

## Scope

In scope:

- the `gbdraw` command-line interface;
- the Python API (`import gbdraw`);
- the hosted Web app at <https://gbdraw.app/>.

Out of scope:

- bugs in third-party dependencies themselves (report them upstream; tell us
  if gbdraw uses one in an unsafe way);
- LOSAT, BLAST, or other binaries that you run yourself.

## Threat model

This section says which inputs an attacker can control and what gbdraw does
with them, so that you can tell whether a problem is a vulnerability.

### What an attacker controls

Any file a user opens may come from someone else: GenBank, EMBL, GFF3 and
FASTA records, BLAST and LOSAT tables, the TSV tables (colors, label lists,
feature visibility), configuration TOML, and Session files
(`.gbdraw-session.json`, `.gbdraw-session.json.gz`). The attacker controls
their content and size, including the text that ends up in a diagram:
qualifiers, labels, record definitions and Session fields.

### Trust boundaries

- **Input files.** gbdraw reads every input as data. It never evaluates code
  from an input: TOML is parsed with `tomllib`, and no input is unpickled or
  run.
- **Diagram text.** The Python renderer escapes text and attribute values when
  it writes SVG. The Web app passes every SVG through DOMPurify before it
  shows it, which removes scripts, `foreignObject`, embedded content and event
  handlers.
- **The browser.** In the hosted Web app at <https://gbdraw.app/>, parsing,
  LOSAT searches, rendering, Session save and load, and export run in the
  browser, in Web Workers (Pyodide and LOSAT WebAssembly). Input files are not
  uploaded: the app sends no request that carries file contents
  ([Execution, privacy, and offline use](docs/REFERENCE/web-app.md#execution-privacy-and-offline-use)).
  The site's Worker serves the packaged files and forwards large Gallery files
  from `raw.githubusercontent.com` at a fixed commit; those requests carry no
  user data.
- **The local GUI.** `gbdraw gui` serves the same app on `127.0.0.1` only.
- **Sessions on the command line and in Python.** Files embedded in a Session
  are written by base name into a new temporary directory. A Session's output
  prefix must be a bare file name, and the overwrite setting stored in a
  Session is ignored: only your own `--overwrite` replaces existing files.
- **External programs.** gbdraw runs LOSAT and BLAST+ with an argument list,
  never through a shell. It runs the binary you name (`--losat_bin`,
  `--ncbi_blast_bin`) or the one it finds in the active conda environment, the
  managed install, or `PATH`; that binary is trusted. `gbdraw setup-losat` installs a LOSAT
  release only when the archive and the binary match the size and SHA-256
  pinned in the gbdraw package.

Examples of vulnerabilities: script that runs in the Web app from a crafted
input or Session; input data that leaves the browser; a crafted Session that
makes the command line or Python read or write outside the directories the
user gave; `gbdraw setup-losat` installing a binary that does not match the
pinned hash.

### Known limits

- The Web app's Content Security Policy still allows `'unsafe-inline'` and
  `'unsafe-eval'` for scripts, so the SVG sanitizer, not the policy, is the
  main defense against script in a crafted input.
- The command line and the Python API read a whole input or Session into
  memory, gzip expanded, with no size limit. The Web app refuses a Session over
  200 MiB compressed or 512 MiB expanded.
- The LOSAT release is checked by hash, not by signature; the hashes are as
  trustworthy as the gbdraw package that carries them.

## Response

This is a small, volunteer-maintained project. We will acknowledge a report
within a reasonable time, say whether we consider it a vulnerability, and
publish a fix and an advisory for supported versions. We give no response-time
guarantee.
