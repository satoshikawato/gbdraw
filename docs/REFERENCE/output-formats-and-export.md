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
`.gbdraw-session.json.gz`. The command line writes the raw LOSAT rows of a run inside the
directory supplied to `--losat_output_dir`: `losatp.raw.tsv` for LOSATP, or one
TSV per record pair and `comparisons.tsv` for LOSATN and TLOSATX.

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
components, and output-target collisions are errors. Each diagram writes its
formats all or nothing: every format is generated before any output is
written, so when a conversion fails, no existing file changes and no new file
is left, with or without overwrite. If the filesystem fails while the outputs
are being written, the error names the files already written. A Circular batch
writes record by record: when one record fails, the records before it stay
written. Without CairoSVG, the command line skips PNG, PDF, EPS, and PS with a
warning and still writes the SVG formats; Python `save_figure_to()`,
`Diagram.save()`, and typed requests fail before writing anything.

## Static, interactive, and raster output

Static SVG contains the diagram and supported semantic groups but no embedded
interactive application. Interactive SVG adds a self-contained style and
script, searchable metadata, popups, zoom and reset controls, and supported
sequence downloads. Open it in a modern browser; many image viewers display
only the static artwork or block its script.

In the Web app, SVG, Interactive SVG, PNG, and PDF export the current Result.
Supported manual legend, title, and Linear scale positions in that Result are
included. Preview-only Layout edit hints, hover outlines, cursors, and the
search and toolbar rows are excluded. Export does not apply draft settings;
select **Generate Diagram** first if they should affect the file. A command or
raw Python recipe alone does not reproduce positions dragged in the browser;
keep the Web Session with the exported figure for that handoff.

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

SVG lengths are CSS pixels (96 px per inch). A web PDF page is the SVG canvas
at 0.75 pt per px, the same physical size as the CLI (CairoSVG) PDF; a
1491 × 984 px canvas becomes a 1118 × 738 pt page. The PDF text layer keeps
the spaces in curved and tick labels, so text such as `1 kbp` can be searched
and copied.

CairoSVG ignores `dominant-baseline` on curved text, such as circular tick
labels, and places `hanging` and `middle` text on other baselines than
browsers. It also centers or right-aligns each italic or roman part of a
mixed-style caption on that part's own width, so the parts overlap. Before
conversion, gbdraw computes the browser baseline offset for that text and the
browser start for each caption part from the packaged font metrics. As a result,
command-line and Python PNG, PDF, EPS, and PS text sits where browsers draw
the SVG. The SVG file is not changed.

Mixed inline text formatting, such as italic markup inside a species label,
does not reliably survive conversion to PNG, PDF, EPS, or PS. Keep SVG when
that formatting must remain exact.

See [Interactive SVG and semantic
hooks](interactive-svg-and-semantic-hooks.md) for stable integration tokens.

## Fonts

Diagram text uses Liberation Sans unless you choose another font: **Font** in
the web app, `--font_family` on the command line, or the
`objects.text.font_family` config override in Python. gbdraw measures labels,
the legend, and titles with the font itself, so they fit it. Liberation Sans,
Serif, and Mono, the names they stand for (such as Arial, Helvetica, Times New
Roman, and Courier), and the generic families (`serif`, `sans-serif`,
`monospace`) are built-in fonts. Other fonts, and a font you load under a
built-in name, are custom fonts.

SVG and Interactive SVG embed each custom font they use, by default, so they
show it on computers that do not have it:

- The file starts with `<defs id="gbdraw-embedded-fonts">`, whose `<style>`
  holds one `@font-face` per font face as a `data:` URL: `font/ttf` for
  TrueType outlines, `font/otf` for OpenType CFF outlines.
- Only the characters the diagram draws are embedded. A font whose OS/2
  `fsType` forbids subsetting (bit `0x0100`) is embedded whole.
- A font whose `fsType` forbids embedding (bit `0x0002` without `0x0004` or
  `0x0008`) or allows only bitmaps (bit `0x0200`) is not embedded, and a
  warning names it. The SVG shows it only where it is installed.
- Built-in fonts are never embedded.
- If embedding fails, the SVG export is an error that names the opt-out below.
- DejaVu Serif adds about 11 KB (+12.6%) to the Circular `HmmtDNA` example and
  about 14 KB (+5.5%) to the Linear `BGC0000708`–`713` comparison.

To keep only the font name in SVG and Interactive SVG, turn embedding off:
`--no_font_embedding` on the command line (also with `--session`),
`embed_fonts=False` in `Diagram.save()`, `to_svg()`, `to_bytes()`, or
`RenderOutputRequest`, and **Embed font in exported SVG** in the web app's
**Font** section. These settings do not change PNG, PDF, EPS, or PS. A Session and the web app's
preview never contain a font file; see [Session and request
compatibility](session-and-request-compatibility.md#fonts).

PNG always draws the chosen font, whatever the embedding setting and the
font's license: a raster image contains no font file. In the web app, when a
font cannot be subset for a PNG, the whole font is used; if that fails too,
the PNG export is an error that asks you to choose another font.

PDF, EPS, and PS carry the fonts they draw with, always as subsets. They never
contain a custom font whose `fsType` forbids embedding, allows only bitmaps, or
forbids subsetting (bit `0x0100`): that text is drawn in the nearest Liberation
font, and a warning (command line, Python) or a notice (web app) names both
fonts. The embedding setting does not change them.

In the web app, SVG and Interactive SVG exports follow **Embed font in exported
SVG**, and PNG export always draws with the loaded fonts. PDF export embeds the
TrueType outlines of the characters each custom font draws, converting an
OpenType CFF font; if the conversion fails, the PDF draws that text in the
nearest Liberation font. A notice after the export names each font that its
license keeps out of an SVG, and each font the PDF replaced together with its
replacement. Exports are not Undo steps.

PNG, PDF, EPS, and PS from the command line and Python go through CairoSVG.
gbdraw adds each `--font_file` font the export draws with to fontconfig's
application fonts, so these formats draw it even when it is not installed where
cairo finds fonts through fontconfig, as on Linux. Elsewhere they draw an
installed font of that name or a substitute; install the font or export SVG.
gbdraw warns once per file when it cannot add it.

## Raw results and FASTA downloads

Raw LOSAT and LOSATP exports contain the search rows before display filtering
or grouping. Generated protein exports use stable readable aliases instead of
session-only runtime handles; uploaded comparison files are left unchanged.
See [Comparison programs, thresholds, and result
semantics](comparison-programs-thresholds-and-results.md#raw-results-and-cache-identity).

When source sequence and endpoint metadata are available, feature and match
popups can download nucleotide or amino-acid FASTA. Pairwise and Circular
matches can export either side or both spans. Collinear block downloads use the
query and subject envelopes represented by the block. Similarity-group actions
operate on the documented member set, not on an inferred phylogeny. In an
Interactive SVG match popup, a similarity group shows its name and **N genes in
M records**, then the member table, with the group ID and counts under
**Technical details**, as in the web app.

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
- the search program, arguments, thresholds, input order, and the raw
  comparison results;
- the command, Python script, or current session used to render the figure;
- the original SVG and the exact submitted derivative; and
- a note describing any manual edits.

Interactive SVG is for browser inspection, not a substitute for a static
publication file unless the receiving system explicitly accepts it. After any
downstream edit or optimization, validate the interactive copy in a browser and
proof the static derivative separately.
