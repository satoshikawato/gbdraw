[Documentation home](../DOCS.md) | [Tutorials](../TUTORIALS/README.md) | [How-to guides](../HOW_TO/README.md) | [Reference](../REFERENCE/README.md)

# Prepare publication and reproducible figures

Keep an editable vector master when the publication workflow allows it. SVG preserves groups, text, stable element identifiers, and metadata. PDF is suitable for vector submission and print. Use PNG only when a raster file is required, and derive its pixel dimensions from the final physical size and requested DPI.

Before submission, inspect the figure at its final placement size. Check label collisions, italic species names, thin strokes, legend order, grayscale contrast, and common color-vision deficiencies. Reopen converted PDF, EPS, or PS output in a second viewer to catch font substitution and changed line breaks.

The reproducibility package should contain:

- immutable input accessions or checksums;
- gbdraw and comparison-runtime versions;
- search program, arguments, thresholds, input order, and retained raw evidence;
- the CLI command, Python script, or current session;
- the editable SVG and exact submitted derivative;
- a note for any manual vector-editor changes.

Editing or optimizing an SVG can remove IDs, data attributes, embedded controls, and popup metadata while leaving the static drawing intact. Edit a copy and retain the original generated file. Validate interactive output in a browser after any downstream tool touches it.

See [Output format and export reference](../REFERENCE/output-formats-and-export.md) for format support and [Session and request compatibility](../REFERENCE/session-and-request-compatibility.md) for persisted-state versions.
