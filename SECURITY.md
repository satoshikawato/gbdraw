# Security policy

## Supported versions

| Version | Supported |
| --- | --- |
| 0.15.x | Yes |
| 0.14.x | Yes |
| Older | No |

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

## Response

This is a small, volunteer-maintained project. We will acknowledge a report
within a reasonable time, say whether we consider it a vulnerability, and
publish a fix and an advisory for supported versions. We give no response-time
guarantee.
