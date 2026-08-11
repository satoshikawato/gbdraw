[Documentation home](../DOCS.md) | [Tutorials](../TUTORIALS/README.md) | [Technical documentation](README.md) | [Publication FAQ](../FAQ.md#how-should-i-prepare-a-figure-for-publication)

# Output formats and export

| Format | Web app | Command line | Python | Main use |
|---|---:|---:|---:|---|
| SVG | Yes | Yes | Yes | Editable vector master and web display |
| Interactive SVG | Yes | Yes | Yes | Browser search, popups, and sequence actions |
| PNG | Yes | Yes | Yes | Fixed-pixel raster output |
| PDF | Yes | Yes | Yes | Vector print and submission |
| EPS | No | Yes | Yes | Requested legacy vector workflow |
| PS | No | Yes | Yes | Requested legacy PostScript workflow |

## Names, dependencies, and overwrite

SVG is the base render. Static SVG uses `<prefix>.svg`; Interactive SVG uses
`<prefix>.interactive.svg`. A session uses `.gbdraw-session.json` or
`.gbdraw-session.json.gz`. Command-line protein evidence uses the path supplied
to `--protein_blastp_output`.

Command-line `-f` or `--format` and Python format arguments use `svg`,
`interactive_svg`, `png`, `pdf`, `eps`, or `ps`.

Command-line and Python conversion to PNG, PDF, EPS, or PS requires CairoSVG
and its runtime dependencies. SVG and Interactive SVG do not. The web app uses
its packaged browser conversion path.

Command-line `-o` or `--output` accepts a path-like prefix, and
`Diagram.save(path)` accepts a path. Typed `RenderOutputRequest.output_prefix`
is one filename component; its directory belongs in `output_directory`.
Existing regular files are replaced only when overwrite is explicit.
Directories, special files, dangling symlinks, unsafe parents, invalid prefix
components, and output-target collisions are errors. A multi-format request
writes formats sequentially; a later conversion failure does not remove files
already written.

## Static, interactive, and raster output

Static SVG contains the diagram and supported semantic groups but no embedded
interactive application. Interactive SVG adds a self-contained style and
script, searchable metadata, popups, zoom and reset controls, and supported
sequence downloads. Open it in a modern browser; many image viewers display
only the static artwork or block its script.

Interactive output intentionally contains script. Input-derived text must not
become executable markup, event-handler attributes, or unsafe links. SVG
optimizers and vector editors can remove the IDs, data attributes, metadata,
or script needed for interaction while leaving the visible drawing intact.
Keep the original generated file and edit a copy.

DPI changes raster pixel dimensions, not SVG or PDF geometry. In the web app,
PNG dimensions scale from the SVG canvas by `DPI / 96`; **96 (Screen)** keeps
the canvas pixel dimensions and **300 (Print)** produces a larger raster.
Check the downloaded dimensions and file signature instead of trusting the
extension.

Mixed inline text formatting, such as italic markup inside a species label,
does not reliably survive conversion to PNG, PDF, EPS, or PS. Keep SVG when
that formatting must remain exact.

See [Interactive SVG and semantic
hooks](interactive-svg-and-semantic-hooks.md) for stable integration tokens.

## Evidence and FASTA downloads

Raw LOSAT and LOSATP exports contain the search rows before display filtering
or grouping. Generated protein exports use stable readable aliases instead of
session-only runtime handles; uploaded comparison files are left unchanged.
See [Comparison programs, thresholds, and result
semantics](comparison-programs-thresholds-and-results.md#raw-results-and-cache-identity).

When source sequence and endpoint metadata are available, feature and match
popups can download nucleotide or amino-acid FASTA. Pairwise and Circular
matches can export either side or both spans. Collinear block downloads use the
query and subject envelopes represented by the block. Similarity-group actions
operate on the documented member set, not on an inferred phylogeny.

## Publication and reproducible handoff

Keep SVG as the editable vector master when the publication workflow accepts
it. PDF is suitable for vector submission and print. Use PNG only when a raster
is required, and calculate its pixel dimensions from the final physical size
and requested DPI. Use EPS or PS only for a workflow that asks for them.

Inspect the figure at its final placement size. Check label collisions, italic
species names, thin strokes, legend order, grayscale contrast, and common
color-vision deficiencies. Reopen converted PDF, EPS, or PS in a second viewer
to catch font substitution, changed line breaks, transparency loss, or stroke
changes.

A reproducible handoff should contain:

- immutable input accessions or checksums;
- gbdraw and comparison-runtime versions;
- the search program, arguments, thresholds, input order, and retained raw
  comparison evidence;
- the command, Python script, or current session used to render the figure;
- the original SVG and the exact submitted derivative; and
- a note describing any manual edits.

Interactive SVG is for browser inspection, not a substitute for a static
publication file unless the receiving system explicitly accepts it. After any
downstream edit or optimization, validate the interactive copy in a browser and
proof the static derivative separately.
