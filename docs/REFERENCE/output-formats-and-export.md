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
directory supplied to `--losat_output_dir`: `losatp.raw.tsv` for LOSATP; for Linear
LOSATN and TLOSATX, one TSV per record pair and `comparisons.tsv`; for Circular
rings, one TSV per ring and `conservation.tsv`.

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
written. Without CairoSVG, a requested PNG, PDF, EPS, or PS fails before
anything is written, with one error that names those formats and the install
command: the command line prints it and exits with status 1, and Python
`save_figure_to()`, `render_to_bytes()`, `Diagram.save()`, and typed requests
raise it as `ValidationError`.

## Figure size and PNG DPI

A figure can be drawn to a requested size: `--figure_width` and
`--figure_height` on the command line, `figure_width` and `figure_height` on
`LinearOptions` and `CircularOptions`, or the config keys `canvas.figure_width`
and `canvas.figure_height`. Each takes `auto` (the default) or a length with a
unit, `mm`, `cm`, `in`, or `px`, such as `180mm`, `18 cm`, or `7in`; a bare
number is refused. A px is a CSS px, 96 per inch. The size is the whole figure,
legend and title included.

gbdraw lays the figure out again to fit, while fonts, line widths, the legend,
record names, and padding keep their size. The figure fills the requested size
to within 3 px, and the rest is added evenly to both sides. A requested size
takes precedence over the length it sets: `canvas.linear.width` for Linear and
`canvas.circular.radius` for Circular, where a set value is only the starting
point.

- **Linear:** the sequence width (`canvas.linear.width`) changes to fill the
  requested width. A requested height adds room above and below; it does not
  stretch the figure. A width narrower than the figure drawn with a 1 px
  sequence, or a height shorter than the figure drawn at the requested width,
  is refused. The message names the record names and labels past the sequence
  ends, a legend beside the records, and the padding, or a title or a legend
  above or below the records when that is wider.
- **Circular:** the map radius changes, for one map or a multi-record canvas
  (where it is the largest record's radius), until the map fits both requested
  sides: the largest radius that fits both, so one side fills the request and
  the other may have room left, added evenly to both sides. The map stops
  shrinking where the centre definition keeps the tracks out, or where the
  tracks no longer fit; a size smaller than the map at its smallest is refused,
  with the map at its smallest, the legend, the title when wider, and the
  padding for a width, or the same parts top to bottom for a height.

`--canvas_padding TOP,RIGHT,BOTTOM,LEFT` (px; config key `canvas.padding`,
default `[0, 0, 0, 0]`) adds room on each side outside the edge padding. At
Auto the canvas grows by the padding; at a requested size the figure fits
inside that size with the padding, and a refusal counts it with the padding.
Legend row edits (`--legend_rows`, see
[`--legend_rows`](../CLI_Reference.md#--legend_rows)) are measured as drawn:
a renamed row takes the width of its new name, and a hidden row takes no room.

Each refusal states the smallest size that works and offers placing a side
legend at the top or bottom when the primary row and that legend set the width.
When labels that need much more room than the map has are not placed at a
requested size, the label error also offers a larger figure with the size the
labels need, for example `make the figure larger: --figure_width about 410.5mm;`.

The SVG, PDF, EPS, and PS files that the command line, `save_figure_to()`,
`render_to_bytes()`, and `Diagram.save()` write state the requested size in its
unit (`width="180mm"`; a side left auto is written in the other side's unit),
with the `viewBox` in px, so a PDF page is 180 mm wide. The in-memory drawing
keeps px.

In the Web app, **Figure size** in the **Layout** card sets the same two sides
(see [Figure size](web-app.md#figure-size)); it applies on Generate, and the
line under it shows the drawn size. The preview keeps px. Its SVG export states
the requested size as the command line writes it; its PDF page has the same
physical size (0.75 pt per px), its PNG has `dpi / 96` pixels per px with the
DPI in the file, and Interactive SVG still opens to fit the browser window. A
Result drawn at Auto exports as before.

`--dpi` on the command line (also with `--session` and on `gbdraw render`)
and `dpi` on `save_figure_to()`, `render_to_bytes()`, `Diagram.to_bytes()`,
`Diagram.save()`, and `RenderOutputRequest` scale a PNG to `dpi / 96` pixels
per SVG pixel. The default is 300, as in the web app, so a PNG is about three
times as wide and tall as its SVG in px; `--dpi 96` (`dpi=96`) draws one pixel
per SVG pixel. Other formats ignore it,
and the command line warns when `--dpi` is given without PNG output. The DPI
is the caller's export choice, like `--font_file`: a Session saved by the
command line does not record `--dpi`. Without `--dpi`, `--session` and
`gbdraw render` draw a PNG at the **PNG DPI** a web app Session saved (72, 96,
150, 300, or 600), as do `render_session()` and `render_session_drawings()` in
Python without `dpi=`; a Session without one draws at the default, 300, and
another stored value is noted and drawn at 300, as the web app reads it. When
the DPI makes a PNG too large, the message names the DPI, says whether it is
the default or the Session's, and names a lower DPI that fits; `--dpi` or
`dpi=` draws it at another. Saving the Session again keeps its PNG DPI. Every
PNG that the command line or Python writes carries a pHYs chunk with that DPI
(the default, 300, when none is given) in pixels per metre, as the web app's
PNG does.

## PNG and PDF size limits

On the command line and in Python, a PNG side cannot exceed 32,767 px, and a
wrapped figure's PDF page side cannot exceed 14,400 pt. A PNG there has
`dpi / 96` pixels per SVG pixel (3.125 at the default 300 DPI, so a side over
10,485 SVG px needs a lower `--dpi`) and a PDF page has
0.75 pt per SVG pixel, so 14,400 pt is 19,200 SVG px. A figure over a limit stops before anything is written, as a failed
conversion does. For a [wrapped figure](command-line.md#wrap-one-record-over-several-rows)
the error names the next longer 1-2-5 row length that should fit, or says that
none would. That length is estimated from the drawn height and the row count,
so the true height can differ when labels or the size class change with the row
length. SVG, Interactive SVG, EPS, and PS have no such limit.

An unwrapped PDF over 14,400 pt is still written, as before; some PDF viewers
may refuse such a page.

The Web app checks the same limits before it draws a PNG or a wrapped figure's
PDF, and downloads nothing when a figure is over. Its PNG has the SVG size times
the chosen DPI / 96 pixels per side, so the message names the highest DPI that
fits when one does, and for a wrapped figure a longer row length that should
fit. SVG downloads have no limit.

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

DPI changes raster pixel dimensions, not SVG or PDF geometry. On the command
line, in Python, and in the web app,
PNG dimensions scale from the SVG canvas by `DPI / 96`; **96 (Screen)** keeps
the canvas pixel dimensions and **300 (Print)** produces a larger raster.
Check the downloaded dimensions and file signature instead of trusting the
extension.

gbdraw draws no background. SVG and Interactive SVG have no background shape,
and command-line, Python, and web app PNG keep the background transparent. PDF,
EPS, and PS paint no page background, so a viewer shows its own page color,
usually white. The web app PNG export fills the background white when
**White background** is on next to the DPI choice; it is off by default, is a
page setting that is not saved in Sessions, and is not an Undo step. The web
app PDF has no background either way.

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
