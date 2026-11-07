[Documentation home](../DOCS.md) | [Tutorials](../TUTORIALS/README.md) | [Technical documentation](README.md) | [FAQ](../FAQ.md)

# Compatibility and deprecation

This page says which gbdraw releases receive fixes, and how long an old
command-line option, Python API name, saved Session, or Python version keeps
working. gbdraw is before 1.0: a `0.MINOR.0` release may include breaking
changes, and the rules below say how they are made and announced.

## Supported releases

Only the newest minor line receives fixes, including security fixes. When a
new `0.MINOR.0` release is published, the previous line receives no more
releases; upgrade to get a fix. [Security policy](../../SECURITY.md#supported-versions)
lists the line that is supported now. The hosted Web app at
<https://gbdraw.app/> always runs the newest release.

## Command-line options

A command-line option or option value is retired in a `0.MINOR.0` release,
without a warning period. A retired option keeps a message: using it exits with
status 2 and says what to use instead. [Retired inputs](session-and-request-compatibility.md#retired-inputs)
lists the retired names, and how a Session or request saved by an earlier
release that still holds one is read.

## Python API

The public Python API is the names that the `gbdraw` package exports
(`gbdraw.__all__`) and the `gbdraw.api` module. Other modules are internal and
may change in any release.

- A public function, class, or argument that is removed first keeps working for
  one minor line and raises a `DeprecationWarning` that names the release that
  removes it. The next `0.MINOR.0` release removes it.
- A field of an option or request dataclass is renamed or removed at once. The
  old field name raises `TypeError`; there is no alias.

## Saved Sessions

A Session saved by gbdraw 0.12.0 or later opens in every later release, with no
end date. If a release stops reading an old Session version, its release
notes say so.

A Session saved by a newer release may not open in an older one. The release
notes of each release that changes the Session version say which releases can
open its Sessions. [Session and request compatibility](session-and-request-compatibility.md)
lists the Session versions that the current release reads.

## Python versions

gbdraw supports the Python versions that [Installation](../INSTALL.md#supported-and-verified-environments)
lists. After the Python project ends its support for a version, the release
notes of the next minor line announce that gbdraw will drop it, and the
`0.MINOR.0` release after that drops it.

## How a deprecation is announced

- `CHANGELOG.md` has a `Deprecated:` entry that names the release that removes
  the item.
- The release notes say what to use instead.
- At run time, the command line prints a warning, Python raises a
  `DeprecationWarning`, and the Web app lists the change when it loads a
  Session that uses the item.
