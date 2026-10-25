const { test, expect } = require('@playwright/test');
const { execFileSync } = require('node:child_process');
const { createHash } = require('node:crypto');
const { resolve } = require('node:path');
const { readFileSync, writeFileSync } = require('node:fs');
const {
  openApp, reveal, generateAndWaitForResult, expectSessionLoadedWithoutWarnings, assertSessionLoadLeftWorkerIdle
} = require('./helpers/app-lifecycle.cjs');

// Custom fonts (D-01, D-06): a font file chosen in the Font section is
// registered in the page and in the diagram Worker, the drawing's font becomes
// its CSS list, and Generate draws and measures the diagram with it.

const SESSION = resolve('gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json');
const FAMILY = 'Gbdraw Test Serif';
const FONT_LIST = `'${FAMILY}', serif`;

// A bundled Liberation Serif face renamed, so the test needs no system font;
// `fsType` is its OS/2 embedding license.
const renamedFont = (target, family = FAMILY, fsType = 0) => {
  execFileSync('python', ['-c', `
import sys
from fontTools.ttLib import TTFont
font = TTFont("gbdraw/data/LiberationSerif-Regular.ttf")
names = font["name"]
names.names = [record for record in names.names if record.nameID not in {1, 2, 4, 6, 16, 17}]
for name_id, value in ((1, sys.argv[1]), (2, "Regular"), (4, sys.argv[1] + " Regular"), (6, sys.argv[1].replace(" ", "")), (16, sys.argv[1])):
    names.setName(value, name_id, 3, 1, 0x409)
font["OS/2"].fsType = int(sys.argv[3])
font.save(sys.argv[2])
`, family, target, String(fsType)], { cwd: resolve('.'), stdio: 'pipe' });
  return target;
};

// Fonts built from Liberation Serif with its metrics, so Python lays out the
// diagram exactly as with Liberation Serif: "box" draws every glyph as its
// bounding box (unmistakable in a raster), "cff" has the same outlines in CFF
// (an OpenType .otf that a PDF needs converted to TrueType).
const builtFont = (target, family, kind) => {
  execFileSync('python', ['-c', `
import sys
from fontTools.fontBuilder import FontBuilder
from fontTools.pens.boundsPen import BoundsPen
from fontTools.pens.t2CharStringPen import T2CharStringPen
from fontTools.pens.ttGlyphPen import TTGlyphPen
from fontTools.ttLib import TTFont
family, target, kind = sys.argv[1:4]
font = TTFont("gbdraw/data/LiberationSerif-Regular.ttf")
order = font.getGlyphOrder()
glyphs = font.getGlyphSet()
if kind == "box":
    bounds = {}
    for name in order:
        pen = BoundsPen(glyphs)
        glyphs[name].draw(pen)
        bounds[name] = pen.bounds
    for name in order:
        pen = TTGlyphPen(None)
        if bounds[name]:
            x0, y0, x1, y1 = bounds[name]
            pen.moveTo((x0, y0)); pen.lineTo((x0, y1)); pen.lineTo((x1, y1)); pen.lineTo((x1, y0)); pen.closePath()
        font["glyf"][name] = pen.glyph()
    names = font["name"]
    names.names = [record for record in names.names if record.nameID not in {1, 2, 4, 6, 16, 17}]
    for name_id, value in ((1, family), (2, "Regular"), (4, family + " Regular"), (6, family.replace(" ", "")), (16, family)):
        names.setName(value, name_id, 3, 1, 0x409)
    font.save(target)
else:
    charstrings = {}
    for name in order:
        pen = T2CharStringPen(font["hmtx"][name][0], glyphs)
        glyphs[name].draw(pen)
        charstrings[name] = pen.getCharString()
    builder = FontBuilder(font["head"].unitsPerEm, isTTF=False)
    builder.setupGlyphOrder(order)
    builder.setupCharacterMap(font.getBestCmap())
    builder.setupCFF(family.replace(" ", ""), {"FullName": family}, charstrings, {})
    builder.setupHorizontalMetrics({name: font["hmtx"][name] for name in order})
    builder.setupHorizontalHeader(ascent=font["hhea"].ascent, descent=font["hhea"].descent)
    builder.setupNameTable({"familyName": family, "styleName": "Regular"})
    os2 = font["OS/2"]
    builder.setupOS2(usWeightClass=400, fsType=0, panose=os2.panose, sTypoAscender=os2.sTypoAscender,
                     sTypoDescender=os2.sTypoDescender, usWinAscent=os2.usWinAscent, usWinDescent=os2.usWinDescent)
    builder.setupPost()
    builder.save(target)
`, family, target, kind], { cwd: resolve('.'), stdio: 'pipe' });
  return target;
};

// Load the Session, load a font file in Font (it becomes the drawing's font), and Generate.
const generateWithFontFile = async (page, fontFile, family) => {
  await openApp(page);
  await page.locator('input[type="file"][accept^=".json,"]').setInputFiles(SESSION);
  await expectSessionLoadedWithoutWarnings(page, { timeout: 180_000 });
  await loadFontFile(page, fontFile, family);
  await generateAndWaitForResult(page);
};

const loadFontFile = async (page, fontFile, family) => {
  const select = page.locator('[data-font-select]');
  await reveal(select);
  const chooser = page.waitForEvent('filechooser');
  await select.selectOption('load-file');
  await (await chooser).setFiles(fontFile);
  await expect(select.locator('option:checked')).toHaveText(family, { timeout: 180_000 });
};

/** Run an export action and return its outcome and the downloaded bytes. */
const exported = async (page, action) => {
  const pending = page.waitForEvent('download');
  const outcome = await page.evaluate((name) => window.__GBDRAW_APP__[name](), action);
  return { outcome, bytes: readFileSync(await (await pending).path()) };
};

const drawnFonts = (page) => page.evaluate(() => {
  const content = String(window.__GBDRAW_APP__.currentResults()[0]?.content || '');
  const label = document.querySelector('.gbdraw-preview-surface svg text[font-family]');
  return {
    families: [...new Set([...content.matchAll(/font-family="([^"]*)"/g)].map((match) => match[1]))],
    legendFace: (content.match(/fontFile(?:&quot;|")\s*:\s*(?:&quot;|")([^&"]*)/) || [])[1] || null,
    previewText: label ? getComputedStyle(label).fontFamily : null
  };
});

test('a font file chosen in Font is one History step and draws and measures the Generated diagram', async ({ page }, testInfo) => {
  test.setTimeout(300_000);
  const fontFile = renamedFont(testInfo.outputPath('gbdraw-test-serif.ttf'));
  await openApp(page);
  await page.locator('input[type="file"][accept^=".json,"]').setInputFiles(SESSION);
  await expectSessionLoadedWithoutWarnings(page, { timeout: 180_000 });

  const select = page.locator('[data-font-select]');
  await reveal(select);
  await expect(select).toHaveValue('default');
  await expect(page.locator('[data-font-custom-note]')).toHaveCount(0);
  // R8: the select's Load font file… opens the file picker.
  const chooser = page.waitForEvent('filechooser');
  await select.selectOption('load-file');
  await (await chooser).setFiles(fontFile);
  // Registration starts the diagram Worker's Python when it is not running yet.
  await expect(select).toHaveValue(FONT_LIST, { timeout: 180_000 });
  await expect(select.locator('option:checked')).toHaveText(FAMILY);
  await expect(page.locator('[data-font-custom-note]')).toBeVisible();
  await expect(page.locator('[data-font-load-error]')).toHaveCount(0);
  expect(await page.evaluate((family) => document.fonts.check(`12px "${family}"`), FAMILY)).toBe(true);

  // Loading and choosing the font is one History step; the font stays loaded.
  await page.getByRole('button', { name: 'Undo', exact: true }).click();
  await expect(select).toHaveValue('default');
  await expect(page.locator('[data-font-custom-note]')).toHaveCount(0);
  await expect(select.locator('option', { hasText: FAMILY })).toHaveCount(1);
  await page.getByRole('button', { name: 'Redo', exact: true }).click();
  await expect(select).toHaveValue(FONT_LIST);

  await generateAndWaitForResult(page);
  const drawn = await drawnFonts(page);
  expect(drawn.families).toEqual([FONT_LIST]);
  expect(drawn.legendFace).toMatch(/^user-[0-9a-f]{16}-Regular$/);
  expect(drawn.previewText).toContain(FAMILY);

  // R1: a new diagram Worker gets the page's fonts again, so a Generate after
  // the Worker is replaced still measures with the loaded face.
  await page.evaluate(async () => {
    const { disposeDiagramGenerationWorker } = await import('/gbdraw/web/js/services/diagram-generation.js');
    disposeDiagramGenerationWorker();
  });
  await generateAndWaitForResult(page);
  expect((await drawnFonts(page)).legendFace).toBe(drawn.legendFace);
  await expect(page.locator('[data-app-notice]', { hasText: 'could not be prepared again' })).toHaveCount(0);
});

// Receipt 〔読み込むとき〕 and 〔警告〕: a Session whose font is neither bundled
// nor loaded asks before it replaces the current work, without the diagram
// Worker. Cancel loads nothing; Load with <nearest> loads it, keeps the name in
// Font, draws the nearest font, and a Save says the Session keeps only the name.
test('a Session with a font that is not available asks first, and Load with the nearest font keeps its name', async ({ page }, testInfo) => {
  test.setTimeout(300_000);
  const MISSING = "'Garamond Premier Pro', serif";
  const session = JSON.parse(readFileSync(SESSION, 'utf8'));
  session.modes.circular.config.adv.font_family = MISSING;
  session.renderRequest.diagramOptions.configOverrides = {
    ...session.renderRequest.diagramOptions.configOverrides, 'objects.text.font_family': MISSING
  };
  const sessionFile = testInfo.outputPath('garamond.gbdraw-session.json');
  writeFileSync(sessionFile, JSON.stringify(session));
  await openApp(page);
  const sessionInput = page.locator('input[type="file"][accept^=".json,"]');
  const prompt = page.locator('[data-choice-prompt]');
  const select = page.locator('[data-font-select]');

  await sessionInput.setInputFiles(sessionFile);
  await expect(prompt.getByRole('heading', { name: 'This Session uses a font that is not available' })).toBeVisible({ timeout: 120_000 });
  await expect(prompt).toContainText('Its text is set in Garamond Premier Pro, which is not part of gbdraw and is not loaded in this page.');
  await expect(prompt).toContainText('Load with Liberation Serif: the Session opens and its text is drawn in Liberation Serif, the nearest available font. The Font setting keeps the name Garamond Premier Pro, so the font is used again once it is loaded. The saved preview keeps its look until you Generate.');
  await expect(prompt.locator('[data-choice-prompt-action="file"]')).toHaveText('Load font file…');
  await prompt.getByRole('button', { name: 'Cancel', exact: true }).click();
  await expect(prompt).toHaveCount(0);
  await expect(page.locator('[data-app-notice]', { hasText: 'Session loaded' })).toHaveCount(0);
  await reveal(select);
  await expect(select).toHaveValue('default');
  await assertSessionLoadLeftWorkerIdle(page);

  await sessionInput.setInputFiles(sessionFile);
  await prompt.locator('[data-choice-prompt-action="substitute"]').click();
  await expectSessionLoadedWithoutWarnings(page, { timeout: 180_000 });
  await expect(select).toHaveValue(MISSING);
  await expect(select.locator('option:checked')).toHaveText('Garamond Premier Pro — not available, drawn in Liberation Serif');
  await expect(page.locator('[data-font-custom-note]')).toBeVisible();
  await assertSessionLoadLeftWorkerIdle(page);

  await generateAndWaitForResult(page);
  expect((await drawnFonts(page)).legendFace).toBe('LiberationSerif-Regular');
  await expect(select).toHaveValue(MISSING);

  await page.getByRole('button', { name: 'Save Session', exact: true }).first().click();
  await expect(page.locator('[data-app-notice]')).toContainText('Session saved', { timeout: 120_000 });
  await expect(page.locator('[data-app-notice] li')).toHaveText([
    'The Session keeps the font name Garamond Premier Pro, not the font file. To open it elsewhere, the same font must be available there, or it is replaced with the nearest available font.'
  ]);
});

// Receipt 〔書き出し〕 (D-05, D-13): SVG and Interactive SVG exports carry a
// subset of the chosen font, which a page without the font shows; PNG export
// draws with it. Embed font in exported SVG off leaves the font name only in
// the SVG; the PNG still draws the font (OV-741).
test('SVG, Interactive SVG and PNG exports carry the chosen font, and the checkbox turns it off for SVG only', async ({ page, browser }, testInfo) => {
  test.setTimeout(300_000);
  const BOX = 'Gbdraw Box Serif';
  await generateWithFontFile(page, builtFont(testInfo.outputPath('box.ttf'), BOX, 'box'), BOX);
  const embed = page.getByRole('checkbox', { name: 'Embed font in exported SVG' });
  await expect(embed).toBeChecked();

  const svg = await exported(page, 'downloadSVG');
  expect(svg.outcome).toMatchObject({ fontWarnings: [] });
  const text = svg.bytes.toString('utf8');
  expect(text).toMatch(new RegExp(`^<svg[^>]*><defs id="gbdraw-embedded-fonts"><style type="text/css">@font-face \\{ font-family: '${BOX}'; font-weight: 400; font-style: normal; src: url\\(data:font/ttf;base64,`));
  expect(text.match(/@font-face/g)).toHaveLength(1);
  const fresh = await browser.newPage();
  await fresh.setContent(`<!doctype html><body>${text}</body>`);
  const faces = await fresh.evaluate(async (family) => {
    await document.fonts.ready;
    return [...document.fonts].filter((face) => face.family.replace(/["']/g, '') === family).map((face) => face.status);
  }, BOX);
  expect(faces).toEqual(['loaded']);
  await fresh.close();

  const ink = (bytes) => page.evaluate(async (base64) => {
    const bitmap = await createImageBitmap(new Blob([Uint8Array.from(atob(base64), (c) => c.charCodeAt(0))], { type: 'image/png' }));
    const canvas = new OffscreenCanvas(bitmap.width, bitmap.height);
    const context = canvas.getContext('2d');
    context.drawImage(bitmap, 0, 0);
    const { data } = context.getImageData(0, 0, bitmap.width, bitmap.height);
    let dark = 0;
    for (let index = 0; index < data.length; index += 4) if (data[index] + data[index + 1] + data[index + 2] < 192) dark += 1;
    return dark;
  }, bytes.toString('base64'));
  const embeddedPng = await exported(page, 'downloadPNG');
  expect(embeddedPng.outcome).toMatchObject({ fontWarnings: [] });

  // Interactive SVG carries the same block, first in the root, beside its own script.
  const interactive = await exported(page, 'downloadInteractiveSVG');
  expect(interactive.outcome).toMatchObject({ fontWarnings: [] });
  const interactiveText = interactive.bytes.toString('utf8');
  expect(interactiveText).toMatch(/^<svg[^>]*><defs id="gbdraw-embedded-fonts"><style type="text\/css">@font-face \{ font-family: 'Gbdraw Box Serif';/);
  expect(interactiveText).toContain('<script');
  expect(interactiveText.match(/gbdraw-embedded-fonts/g)).toHaveLength(1);

  await embed.uncheck();
  const plainSvg = await exported(page, 'downloadSVG');
  expect(plainSvg.bytes.toString('utf8')).not.toContain('gbdraw-embedded-fonts');
  const plainInteractive = await exported(page, 'downloadInteractiveSVG');
  expect(plainInteractive.bytes.toString('utf8')).not.toContain('gbdraw-embedded-fonts');
  // OV-741: the PNG draws the chosen font whatever the checkbox.
  const plainPng = await exported(page, 'downloadPNG');
  const embeddedInk = await ink(embeddedPng.bytes);
  expect(Math.abs(await ink(plainPng.bytes) - embeddedInk)).toBeLessThanOrEqual(embeddedInk * 0.01);
  // The checkbox is page state and records no History step: Undo leaves it off.
  await expect(page.getByRole('button', { name: 'Undo', exact: true })).toBeEnabled();
  await page.getByRole('button', { name: 'Undo', exact: true }).click();
  await expect(embed).not.toBeChecked();
});

// Receipt 〔書き出し〕 "PDF も選んだ font で描く": a custom TrueType face, and an
// OpenType (CFF) face converted to TrueType outlines in the Worker.
test('PDF export embeds a custom TrueType face and a converted OpenType CFF face', async ({ page }, testInfo) => {
  test.setTimeout(300_000);
  const BOX = 'Gbdraw Box Serif';
  const CFF = 'Gbdraw Cff Serif';
  // jsPDF names a face by the identity of its file: user-<16 hex of its SHA-256>.
  const pdfFontName = (file) => `/BaseFont /user-${createHash('sha256').update(readFileSync(file)).digest('hex').slice(0, 16)}`;
  const boxFile = builtFont(testInfo.outputPath('box.ttf'), BOX, 'box');
  await generateWithFontFile(page, boxFile, BOX);
  const ttf = await exported(page, 'downloadPDF');
  expect(ttf.outcome).toMatchObject({ fontWarnings: [] });
  const ttfText = ttf.bytes.toString('latin1');
  expect(ttfText).toContain(pdfFontName(boxFile));
  expect(ttfText).toContain('/FontFile2');
  expect(ttfText).not.toMatch(/\/BaseFont \/LiberationSerif/);

  const cffFile = builtFont(testInfo.outputPath('cff.otf'), CFF, 'cff');
  await loadFontFile(page, cffFile, CFF);
  await generateAndWaitForResult(page);
  const cff = await exported(page, 'downloadPDF');
  expect(cff.outcome).toMatchObject({ fontWarnings: [] });
  const cffText = cff.bytes.toString('latin1');
  expect(cffText).toContain(pdfFontName(cffFile));
  expect(cffText).toContain('/FontFile2');
  await expect(page.locator('[data-app-notice]', { hasText: 'Exported' })).toHaveCount(0);
});

// D-13: a font whose license forbids embedding is named in the export notice:
// an SVG keeps only its name, a PDF draws it in the nearest Liberation font, and
// a PNG draws it (rasterizing is not embedding; OV-741).
test('a font whose license forbids embedding is named in the SVG and PDF notices and drawn in the PNG', async ({ page }, testInfo) => {
  test.setTimeout(300_000);
  const LOCKED = 'Gbdraw Locked Serif';
  await generateWithFontFile(page, renamedFont(testInfo.outputPath('locked.ttf'), LOCKED, 0x0002), LOCKED);
  const notice = page.locator('[data-app-notice]');

  const svg = await exported(page, 'downloadSVG');
  const svgWarning = `Font ${LOCKED} is not embedded: its license (OS/2 fsType) does not allow embedding. The SVG shows it only where ${LOCKED} is installed.`;
  expect(svg.outcome.fontWarnings).toEqual([svgWarning]);
  expect(svg.bytes.toString('utf8')).not.toContain('gbdraw-embedded-fonts');
  await expect(notice).toContainText(`Exported ${svg.outcome.filename}`);
  await expect(notice.locator('li')).toHaveText([svgWarning]);

  const pdf = await exported(page, 'downloadPDF');
  const pdfWarning = `Font ${LOCKED} is not embedded in the PDF: its license (OS/2 fsType) does not allow embedding. Its text is drawn in Liberation Serif.`;
  expect(pdf.outcome.fontWarnings).toEqual([pdfWarning]);
  const pdfText = pdf.bytes.toString('latin1');
  expect(pdfText).toMatch(/\/BaseFont \/LiberationSerif/);
  expect(pdfText).not.toMatch(/\/BaseFont \/user-/);
  await expect(notice.locator('li')).toHaveText([pdfWarning]);

  const png = await exported(page, 'downloadPNG');
  expect(png.outcome.fontWarnings).toEqual([]);
});

// Choose installed font… (B, Local Font Access) with `queryLocalFonts`
// stubbed before the app loads: the picker lists the families, searches them,
// and loads the chosen family's faces; a refused permission is named.
test('Choose installed font… lists, searches, and loads an installed family, and names a refused permission', async ({ page }, testInfo) => {
  test.setTimeout(300_000);
  const INSTALLED = 'Gbdraw Installed Serif';
  const base64 = readFileSync(renamedFont(testInfo.outputPath('installed.ttf'), INSTALLED)).toString('base64');
  await page.addInitScript(({ family, data }) => {
    const bytes = Uint8Array.from(atob(data), (c) => c.charCodeAt(0));
    window.__gbdrawLocalFontsDenied = false;
    window.queryLocalFonts = async () => {
      if (window.__gbdrawLocalFontsDenied) throw new DOMException('Permission denied', 'NotAllowedError');
      return [
        { family, fullName: `${family} Regular`, style: 'Regular', blob: async () => new Blob([bytes]) },
        { family: 'Other Installed Sans', fullName: 'Other Installed Sans Regular', style: 'Regular', blob: async () => new Blob([bytes]) }
      ];
    };
  }, { family: INSTALLED, data: base64 });
  await openApp(page);
  await page.locator('input[type="file"][accept^=".json,"]').setInputFiles(SESSION);
  await expectSessionLoadedWithoutWarnings(page, { timeout: 180_000 });
  const select = page.locator('[data-font-select]');
  await reveal(select);
  await expect(page.locator('[data-font-installed-unavailable]')).toHaveCount(0);
  const picker = page.locator('[data-font-installed-picker]');

  await select.selectOption('installed');
  await expect(picker.getByRole('heading', { name: 'Choose installed font' })).toBeVisible();
  const families = picker.getByRole('list', { name: 'Installed font families' }).getByRole('button');
  await expect(families).toHaveText([`${INSTALLED} (1)`, 'Other Installed Sans (1)']);
  await picker.getByRole('searchbox', { name: 'Search installed fonts' }).fill('serif');
  await expect(families).toHaveText([`${INSTALLED} (1)`]);
  await families.first().click();
  await expect(picker).toHaveCount(0);
  await expect(select).toHaveValue(`'${INSTALLED}', serif`, { timeout: 180_000 });
  await expect(page.locator('[data-font-custom-note]')).toBeVisible();

  await page.evaluate(() => { window.__gbdrawLocalFontsDenied = true; });
  await select.selectOption('installed');
  await expect(picker.getByRole('alert')).toHaveText('Installed fonts are not available: the browser did not allow access to local fonts.');
  await picker.getByRole('button', { name: 'Cancel', exact: true }).click();
  await expect(picker).toHaveCount(0);
  await expect(select).toHaveValue(`'${INSTALLED}', serif`);
});

// The PDF draws each label once (a halo stays in the SVG, LP-D102), the
// default Linear box of each label at its opacity (LP-D106), and a centred or
// end-anchored text where the SVG draws it: svg2pdf measured the anchor with
// the browser, which knows the PDF face under no name, and drew long lambda
// labels up to 30 px right (OV-1303).
test('PDF export draws each label once on its box, centred on its SVG anchor', async ({ page }, testInfo) => {
  test.setTimeout(300_000);
  await openApp(page);
  await page.locator('input[type="file"][accept^=".json,"]').setInputFiles(resolve('gbdraw/web/gallery/sessions/lambda_basic_linear.gbdraw-session.json'));
  await expectSessionLoadedWithoutWarnings(page, { timeout: 180_000 });
  await generateAndWaitForResult(page);
  const svgFile = testInfo.outputPath('lambda.svg');
  writeFileSync(svgFile, await page.evaluate(() => String(window.__GBDRAW_APP__.currentResults()[0].content)));
  const pdf = await exported(page, 'downloadPDF');
  const pdfFile = testInfo.outputPath('lambda.pdf');
  writeFileSync(pdfFile, pdf.bytes);
  // Per anchored SVG text: the PDF runs of its characters, and the distance (px) between the PDF run's
  // anchor point (its left x plus its share of its width in the PDF face) and the SVG anchor x.
  const report = JSON.parse(execFileSync('python', ['-c', `
import json, re, sys, xml.etree.ElementTree as ET
from fontTools.ttLib import TTFont
svg, pdf = sys.argv[1:3]
font = TTFont("gbdraw/data/LiberationSans-Regular.ttf")
order, units, advance = font.getGlyphOrder(), font["head"].unitsPerEm, font["hmtx"]
by_glyph = {}
for code, name in font.getBestCmap().items():
    by_glyph.setdefault(name, chr(code))
runs = {}
for x, glyphs in re.findall(r"[-\\d.]+ [-\\d.]+ [-\\d.]+ [-\\d.]+ ([-\\d.]+) [-\\d.]+ Tm\\n<([0-9a-f]*)> Tj", open(pdf, "rb").read().decode("latin1")):
    names = [order[int(glyphs[i:i + 4], 16)] for i in range(0, len(glyphs), 4)]
    runs.setdefault("".join(by_glyph.get(name, "?") for name in names), []).append((float(x) / 0.75, sum(advance[name][0] for name in names) / units))
root = ET.parse(svg).getroot()
parents = {child: parent for parent in root.iter() for child in parent}
def left(element):
    x = 0.0
    while element is not None:
        x += sum(float(dx) for dx, _ in re.findall(r"translate\\(([^,)\\s]+)[,\\s]+([^)]+)\\)", element.get("transform", "")))
        element = parents.get(element)
    return x
texts = [text for text in root.iter("{http://www.w3.org/2000/svg}text") if text.get("data-label-feature-id")]
boxes = [rect for rect in root.iter("{http://www.w3.org/2000/svg}rect") if rect.get("data-label-feature-id")]
opacities = sorted({rect.get("fill-opacity") for rect in boxes})
data = open(pdf, "rb").read().decode("latin1")
# The graphics states of the box opacity, by resource name, and how often the page sets one of them.
states = {number for number, body in re.findall(r"(\\d+) 0 obj\\s*<<(.*?)>>", data, re.S) if re.search(r"/ca\\s+" + re.escape(opacities[0]) + r"\\b", body)} if len(opacities) == 1 else set()
state_names = {name for name, number in re.findall(r"/(\\w+)\\s+(\\d+) 0 R", data) if number in states}
same = {}
for text in texts:
    same[text.text] = same.get(text.text, 0) + 1
offsets = []
for text in texts:
    share = {"start": 0.0, "middle": 0.5, "end": 1.0}[text.get("text-anchor", "start")]
    anchor = left(text) + float(text.get("x", "0"))
    offsets.append(min((abs(x + share * width * float(text.get("font-size")) - anchor) for x, width in runs.get(text.text, [])), default=None))
print(json.dumps({
    "labels": len(texts),
    "boxes": len(boxes),
    "boxOpacities": opacities,
    "boxStateUses": sum(len(re.findall(r"/" + name + r"\\s+gs", data)) for name in state_names),
    "runsPerLabel": sorted({len(runs.get(text.text, [])) / same[text.text] for text in texts}),
    "maxOffset": max((offset for offset in offsets if offset is not None), default=None),
}))
`, svgFile, pdfFile], { cwd: resolve('.'), encoding: 'utf8' }));
  expect(report.labels).toBeGreaterThan(50);
  expect(report.boxes, 'labels with a box in the SVG').toBeGreaterThan(50);
  expect(report.boxOpacities, 'box opacity in the SVG').toEqual(['0.8']);
  expect(report.boxStateUses, 'PDF fills at the box opacity').toBeGreaterThanOrEqual(report.boxes);
  expect(report.runsPerLabel, 'PDF text runs per label').toEqual([1]);
  expect(report.maxOffset, 'PDF anchor offset (px)').toBeLessThan(0.05);
});
