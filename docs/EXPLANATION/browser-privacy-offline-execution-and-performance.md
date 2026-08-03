[Documentation home](../DOCS.md) | [Tutorials](../TUTORIALS/README.md) | [How-to guides](../HOW_TO/README.md) | [Reference](../REFERENCE/README.md)

# Understand browser privacy, offline execution, and performance

The hosted page must be downloaded from `gbdraw.app`, but genome inputs are parsed and rendered in the browser. Browser LOSAT searches run in bundled workers. gbdraw does not require an application-server upload for generation, search, session save, or export. The production site may load the aggregate page-usage analytics disclosed in the project README; uploaded genome files and generated diagrams are not analytics payloads.

`gbdraw gui` serves the same static application and bundled runtimes from the installed Python package. After installation, the normal drawing workflow can run without an internet connection. Browser security headers and the application content-security policy still apply in local mode.

This boundary does not make every browser extension, proxy, operating system, or downloaded file trustworthy. Use a controlled browser profile for sensitive data, inspect the deployed origin, and keep the application version with saved sessions. A session contains the resources needed to reproduce its state and should be handled like the input data it embeds.

Performance depends on record length, feature count, number of records, comparison-pair count, search mode, thresholds, browser memory, and CPU. There is no single safe genome-size limit. Begin with serial execution and one worker for a reproducible comparison, then increase concurrency only after measuring the same inputs. Large assemblies, dense all-pairs searches, and large depth tables are better suited to a local controlled environment or prepared CLI evidence.

Comparison caches avoid repeating identical work when the program, input identities, order, and settings match. Changing an input, threshold, mode, or relevant runtime identity invalidates that reuse. A cache is an optimization, not the evidence record; retain the exported raw results.

The committed documentation capture blocks non-local requests and starts every scenario in a fresh context. This verifies that public Tutorials do not depend on accounts, mutable remote data, or hidden browser state.
