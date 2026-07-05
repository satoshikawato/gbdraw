const { test, expect } = require('@playwright/test');
const { createReadStream, existsSync, readFileSync, statSync } = require('node:fs');
const { createServer } = require('node:http');
const { extname, join, normalize, resolve, sep } = require('node:path');

const repoRoot = resolve(process.env.GBDRAW_REPO || process.cwd());
const webRoot = join(repoRoot, 'gbdraw/web');

const contentTypes = {
  '.css': 'text/css; charset=utf-8',
  '.html': 'text/html; charset=utf-8',
  '.js': 'text/javascript; charset=utf-8',
  '.json': 'application/json; charset=utf-8',
  '.mp4': 'video/mp4',
  '.ogg': 'video/ogg',
  '.png': 'image/png',
  '.svg': 'image/svg+xml',
  '.webm': 'video/webm',
  '.webp': 'image/webp'
};

let server;
let baseUrl;

const readGalleryExamples = () =>
  JSON.parse(readFileSync(join(webRoot, 'gallery/examples.json'), 'utf8'));

const readTutorialJson = (filename) =>
  readFileSync(join(webRoot, 'gallery/tutorials', filename), 'utf8');

const readTutorialWithOperationMedia = () => {
  const tutorial = JSON.parse(readTutorialJson('hepatoplasmataceae_collinear.json'));
  ['manualSteps', 'colorRules', 'postGenerationEdits'].forEach((sectionName) => {
    (tutorial[sectionName] || []).forEach((step) => {
      if (step && typeof step === 'object') delete step.operations;
    });
  });
  tutorial.manualSteps[0].operations = [
    {
      title: 'Select Linear mode',
      body: 'Use Linear mode before adding the five GenBank inputs.',
      media: {
        src: './media/hepatoplasmataceae_collinear/manual-03-01-open-pairwise.webp',
        alt: 'Pairwise Comparisons settings cropped as an operation-level tutorial screenshot.',
        caption: 'Operation-level media uses the same image path as the current tutorial.'
      }
    }
  ];
  return JSON.stringify(tutorial);
};

test.beforeAll(async () => {
  await new Promise((resolveServer, rejectServer) => {
    server = createServer((request, response) => {
      const url = new URL(request.url || '/', 'http://127.0.0.1');
      const requestedPath = normalize(decodeURIComponent(url.pathname)).replace(/^(\.\.(?:\/|\\|$))+/, '');
      const filePath = resolve(webRoot, requestedPath.replace(/^[/\\]+/, ''));
      if (!filePath.startsWith(`${webRoot}${sep}`) && filePath !== webRoot) {
        response.writeHead(403);
        response.end('Forbidden');
        return;
      }

      if (!existsSync(filePath)) {
        response.writeHead(404);
        response.end('Not found');
        return;
      }

      const finalPath = statSync(filePath).isDirectory() ? join(filePath, 'index.html') : filePath;
      if (!existsSync(finalPath)) {
        response.writeHead(404);
        response.end('Not found');
        return;
      }

      response.writeHead(200, {
        'Content-Type': contentTypes[extname(finalPath)] || 'application/octet-stream'
      });
      createReadStream(finalPath).pipe(response);
    });
    server.once('error', rejectServer);
    server.listen(0, '127.0.0.1', () => {
      const address = server.address();
      baseUrl = `http://127.0.0.1:${address.port}`;
      resolveServer();
    });
  });
});

test.afterAll(async () => {
  await new Promise((resolveClose) => server.close(resolveClose));
});

test('Gallery renders the Hepatoplasmataceae tutorial and files panels', async ({ page }) => {
  const pageErrors = [];
  page.on('pageerror', (error) => pageErrors.push(error.message));

  await page.goto(`${baseUrl}/gallery/#hepatoplasmataceae_collinear`, { waitUntil: 'domcontentloaded' });
  await expect(page.getByRole('heading', { name: /Hepatoplasmataceae.*collinear analysis/i })).toBeVisible();

  await page.getByRole('tab', { name: 'Tutorial' }).click();
  await expect(
    page.getByRole('heading', { name: 'Reproduce the Hepatoplasmataceae collinear analysis in the web app' })
  ).toBeVisible();
  await expect(page.getByText('Web app steps')).toBeVisible();
  await expect(page.getByText('Use browser LOSAT')).toBeVisible();
  await expect(page.getByText('Use Orientation + identity for Color mode.')).toBeVisible();
  const tutorialPanel = page.getByRole('tabpanel', { name: 'Tutorial' });
  await expect(tutorialPanel.getByRole('heading', { name: 'Quick reproduce' })).toHaveCount(0);
  const mediaImages = tutorialPanel.getByRole('img');
  await expect(mediaImages).toHaveCount(12);
  await expect(mediaImages.first()).toHaveAttribute('src', /manual-01-01-linear-mode\.webp$/);
  await expect(tutorialPanel.locator('img[src$="manual-02-01-upload-row-context.webp"]')).toHaveCount(1);
  await expect(tutorialPanel.locator('img[src$="manual-02-01-five-file-upload-order.webp"]')).toHaveCount(0);
  await expect(tutorialPanel.locator('img[src$="manual-07-01-collinear-block-popup.webp"]')).toHaveCount(1);
  for (let idx = 0; idx < await mediaImages.count(); idx += 1) {
    const image = mediaImages.nth(idx);
    await image.scrollIntoViewIfNeeded();
    await expect.poll(() => image.evaluate((element) => element.complete && element.naturalWidth > 0)).toBe(true);
  }

  await page.getByRole('tab', { name: 'Files' }).click();
  const filesPanel = page.getByRole('tabpanel', { name: 'Files' });
  await expect(filesPanel.getByText('AP027078.gb')).toBeVisible();
  await expect(filesPanel.getByText('NZ_CP006932.gb')).toBeVisible();
  await expect(filesPanel.getByRole('link', { name: 'Session JSON' })).toBeVisible();

  expect(pageErrors).toEqual([]);
});

test('Gallery renders the Hepatoplasmataceae orthogroup tutorial and media', async ({ page }) => {
  const pageErrors = [];
  page.on('pageerror', (error) => pageErrors.push(error.message));

  await page.goto(`${baseUrl}/gallery/#hepatoplasmataceae_orthogroup`, { waitUntil: 'domcontentloaded' });
  await expect(page.getByRole('heading', { name: /Hepatoplasmataceae.*orthogroup matches/i })).toBeVisible();

  await page.getByRole('tab', { name: 'Tutorial' }).click();
  await expect(
    page.getByRole('heading', { name: 'Reproduce the Hepatoplasmataceae orthogroup comparison in the web app' })
  ).toBeVisible();
  await expect(page.getByText('blastp mode to Orthogroups')).toBeVisible();
  await expect(page.getByText('orthogroup ID, display name, member count')).toBeVisible();
  const tutorialPanel = page.getByRole('tabpanel', { name: 'Tutorial' });
  const mediaImages = tutorialPanel.getByRole('img');
  await expect(mediaImages).toHaveCount(13);
  await expect(tutorialPanel.locator('img[src$="manual-02-01-upload-row-context.webp"]')).toHaveCount(1);
  await expect(tutorialPanel.locator('img[src$="manual-02-01-five-file-upload-order.webp"]')).toHaveCount(0);
  await expect(tutorialPanel.locator('img[src$="manual-07-01-orthogroup-overview.webp"]')).toHaveCount(1);
  await expect(tutorialPanel.locator('img[src$="manual-04-01-orthogroups-mode.webp"]')).toHaveCount(1);
  for (let idx = 0; idx < await mediaImages.count(); idx += 1) {
    const image = mediaImages.nth(idx);
    await image.scrollIntoViewIfNeeded();
    await expect.poll(() => image.evaluate((element) => element.complete && element.naturalWidth > 0)).toBe(true);
  }

  await page.getByRole('tab', { name: 'Files' }).click();
  const filesPanel = page.getByRole('tabpanel', { name: 'Files' });
  await expect(filesPanel.getByText('AP027078.gb')).toBeVisible();
  await expect(filesPanel.getByRole('link', { name: 'Session JSON' })).toBeVisible();

  expect(pageErrors).toEqual([]);
});

test('Gallery renders the aminoglycoside BGC tutorial and media', async ({ page }) => {
  const pageErrors = [];
  page.on('pageerror', (error) => pageErrors.push(error.message));

  await page.goto(`${baseUrl}/gallery/#BGC0000708-BGC0000713`, { waitUntil: 'domcontentloaded' });
  await expect(page.getByRole('heading', { name: /Aminoglycoside biosynthetic gene clusters/i })).toBeVisible();

  await page.getByRole('tab', { name: 'Tutorial' }).click();
  await expect(
    page.getByRole('heading', { name: 'Reproduce the aminoglycoside BGC comparison in the web app' })
  ).toBeVisible();
  const tutorialPanel = page.getByRole('tabpanel', { name: 'Tutorial' });
  await expect(tutorialPanel.getByRole('heading', { name: 'Color rule basics' })).toBeVisible();
  await expect(
    tutorialPanel.getByRole('row', {
      name: 'CDS gene dnaA #d95f02 DnaA'
    })
  ).toBeVisible();
  await expect(page.getByText('antiSMASH gene-kind color rules')).toBeVisible();
  await expect(tutorialPanel.getByRole('columnheader', { name: 'Feature' }).first()).toBeVisible();
  await expect(tutorialPanel.getByRole('columnheader', { name: 'Qualifier' }).first()).toBeVisible();
  await expect(
    tutorialPanel.getByRole('row', {
      name: 'CDS gene_kind biosynthetic$ #d03535 Core biosynthetic genes'
    })
  ).toBeVisible();
  await expect(
    tutorialPanel.getByRole('row', {
      name: 'Reverse complement BGC0000713 only'
    })
  ).toBeVisible();
  const bgc0708MibigHref = 'https://mibig.secondarymetabolites.org/repository/BGC0000708.5/BGC0000708.gbk';
  await expect(tutorialPanel.locator(`a[href="${bgc0708MibigHref}"]`)).toHaveText('MIBiG repository');
  const mediaImages = tutorialPanel.getByRole('img');
  await expect(mediaImages).toHaveCount(20);
  await expect(tutorialPanel.locator('img[src$="manual-09-01-orthogroup-popup.webp"]')).toHaveCount(1);
  await expect(tutorialPanel.locator('img[src$="manual-10-01-feature-popup.webp"]')).toHaveCount(1);
  await expect(tutorialPanel.locator('img[src$="manual-04-03-track-layout-middle.webp"]')).toHaveCount(1);
  await expect(tutorialPanel.locator('img[src$="manual-04-04-pairwise-style-curve.webp"]')).toHaveCount(1);
  await expect(tutorialPanel.locator('img[src$="manual-07-01-specific-rules-all.webp"]')).toHaveCount(1);
  for (let idx = 0; idx < await mediaImages.count(); idx += 1) {
    const image = mediaImages.nth(idx);
    await image.scrollIntoViewIfNeeded();
    await expect.poll(() => image.evaluate((element) => element.complete && element.naturalWidth > 0)).toBe(true);
  }

  await page.getByRole('tab', { name: 'Files' }).click();
  const filesPanel = page.getByRole('tabpanel', { name: 'Files' });
  await expect(filesPanel.getByText('BGC0000708.gbk')).toBeVisible();
  await expect(filesPanel.locator(`a[href="${bgc0708MibigHref}"]`)).toHaveText('MIBiG repository');
  await expect(filesPanel.getByText('BGC0000713.gbk')).toBeVisible();
  await expect(filesPanel.getByRole('link', { name: 'Session JSON' })).toBeVisible();

  expect(pageErrors).toEqual([]);
});

test('Gallery renders the WSSV conservation tutorial and media', async ({ page }) => {
  const pageErrors = [];
  page.on('pageerror', (error) => pageErrors.push(error.message));

  await page.goto(`${baseUrl}/gallery/#WSSV_genome_comparison`, { waitUntil: 'domcontentloaded' });
  await expect(page.getByRole('heading', { name: /White spot syndrome virus genome comparison/i })).toBeVisible();

  await page.getByRole('tab', { name: 'Tutorial' }).click();
  await expect(
    page.getByRole('heading', { name: 'Reproduce the WSSV circular conservation ring comparison in the web app' })
  ).toBeVisible();
  const tutorialPanel = page.getByRole('tabpanel', { name: 'Tutorial' });
  await expect(tutorialPanel.getByText('browser LOSAT blastn results')).toBeVisible();
  await expect(
    tutorialPanel.getByRole('row', {
      name: 'Ring width 5'
    })
  ).toBeVisible();
  await expect(
    tutorialPanel.getByRole('row', {
      name: '20 Angostura2013 #c6bebb'
    })
  ).toBeVisible();
  await expect(tutorialPanel.getByRole('cell', { name: 'MG18PR-0187-N40S.fa' })).toBeVisible();
  const mediaImages = tutorialPanel.getByRole('img');
  await expect(mediaImages).toHaveCount(15);
  await expect(tutorialPanel.locator('img[src$="manual-04-02-comparison-fasta-series.webp"]')).toHaveCount(1);
  await expect(tutorialPanel.locator('img[src$="manual-04-02-upload-fasta-comparisons.webp"]')).toHaveCount(0);
  await expect(tutorialPanel.locator('img[src$="manual-09-01-conservation-rings.webp"]')).toHaveCount(1);
  await expect(tutorialPanel.locator('img[src$="manual-10-01-files-tab.webp"]')).toHaveCount(0);
  await expect(tutorialPanel.locator('img[src$="manual-11-01-feature-popup.webp"]')).toHaveCount(1);
  await expect(tutorialPanel.locator('img[src$="manual-08-01-browser-losat-run.webp"]')).toHaveCount(1);
  for (let idx = 0; idx < await mediaImages.count(); idx += 1) {
    const image = mediaImages.nth(idx);
    await image.scrollIntoViewIfNeeded();
    await expect.poll(() => image.evaluate((element) => element.complete && element.naturalWidth > 0)).toBe(true);
  }

  await page.getByRole('tab', { name: 'Files' }).click();
  const filesPanel = page.getByRole('tabpanel', { name: 'Files' });
  await expect(filesPanel.getByText('AP027280.gb')).toBeVisible();
  await expect(filesPanel.getByText('CN01.fasta')).toBeVisible();
  await expect(filesPanel.getByText('Angostura2013.fa')).toBeVisible();
  await expect(filesPanel.getByRole('link', { name: 'Session JSON' })).toBeVisible();

  expect(pageErrors).toEqual([]);
});

test('Gallery renders the human mitochondrial AT skew tutorial and media', async ({ page }) => {
  const pageErrors = [];
  page.on('pageerror', (error) => pageErrors.push(error.message));

  await page.goto(`${baseUrl}/gallery/#HmmtDNA_ATskew`, { waitUntil: 'domcontentloaded' });
  await expect(page.getByRole('heading', { name: /Homo sapiens.*AT skew/i })).toBeVisible();

  await page.getByRole('tab', { name: 'Tutorial' }).click();
  await expect(
    page.getByRole('heading', { name: 'Reproduce the human mitochondrial AT skew diagram in the web app' })
  ).toBeVisible();
  const tutorialPanel = page.getByRole('tabpanel', { name: 'Tutorial' });
  await expect(tutorialPanel.getByText('Open Custom Track Slots', { exact: true })).toBeVisible();
  await expect(
    tutorialPanel.getByRole('row', {
      name: 'gc_content Dinucleotide content inside Width: 0.1'
    })
  ).toBeVisible();
  await expect(
    tutorialPanel.getByRole('row', {
      name: 'Dinucleotide AT'
    })
  ).toBeVisible();
  await expect(
    tutorialPanel.getByRole('row', {
      name: 'Legend label AT skew'
    })
  ).toBeVisible();
  const mediaImages = tutorialPanel.getByRole('img');
  await expect(mediaImages).toHaveCount(16);
  await expect(tutorialPanel.locator('img[src$="manual-09-01-atskew-preview.webp"]')).toHaveCount(1);
  await expect(tutorialPanel.locator('img[src$="manual-10-01-feature-popup.webp"]')).toHaveCount(1);
  const atSkewSlotImage = tutorialPanel.locator('img[src$="manual-06-01-at-skew-slot-context.webp"]');
  const tickTrackImage = tutorialPanel.locator('img[src$="manual-07-01-tick-track-context.webp"]');
  await expect(atSkewSlotImage).toHaveCount(1);
  await expect(tickTrackImage).toHaveCount(1);
  await expect.poll(() => atSkewSlotImage.evaluate((element) => element.naturalHeight)).toBeGreaterThan(1000);
  await expect.poll(() => tickTrackImage.evaluate((element) => element.naturalHeight)).toBeGreaterThan(800);
  for (let idx = 0; idx < await mediaImages.count(); idx += 1) {
    const image = mediaImages.nth(idx);
    await image.scrollIntoViewIfNeeded();
    await expect.poll(() => image.evaluate((element) => element.complete && element.naturalWidth > 0)).toBe(true);
  }

  await page.getByRole('tab', { name: 'Files' }).click();
  const filesPanel = page.getByRole('tabpanel', { name: 'Files' });
  await expect(filesPanel.getByText('HmmtDNA.gbk')).toBeVisible();
  await expect(filesPanel.getByRole('link', { name: 'Session JSON' })).toBeVisible();

  expect(pageErrors).toEqual([]);
});

test('Gallery renders the majanivirus orthogroup tutorial and media', async ({ page }) => {
  const pageErrors = [];
  page.on('pageerror', (error) => pageErrors.push(error.message));

  await page.goto(`${baseUrl}/gallery/#majanivirus_orthogroup`, { waitUntil: 'domcontentloaded' });
  await expect(page.getByRole('heading', { name: /Large dsDNA viruses/i })).toBeVisible();

  await page.getByRole('tab', { name: 'Tutorial' }).click();
  await expect(
    page.getByRole('heading', { name: 'Reproduce the large majanivirus orthogroup comparison in the web app' })
  ).toBeVisible();
  const tutorialPanel = page.getByRole('tabpanel', { name: 'Tutorial' });
  await expect(
    tutorialPanel.getByRole('row', {
      name: 'protein blastp mode Orthogroups'
    })
  ).toBeVisible();
  await expect(
    tutorialPanel.getByRole('row', {
      name: 'Specific rule CDS product wsv.*-like protein #89d1fa WSSV-like proteins'
    })
  ).toBeVisible();
  await expect(tutorialPanel.getByRole('cell', { name: 'unmatched CDS' })).toBeVisible();
  await expect(tutorialPanel.getByText('Use 32 threads only deliberately')).toBeVisible();
  await expect(
    tutorialPanel.getByRole('row', {
      name: '1 MjeNMV.gb Marsupenaeus japonicus endogenous nimavirus'
    })
  ).toBeVisible();
  const mediaImages = tutorialPanel.getByRole('img');
  await expect(mediaImages).toHaveCount(15);
  await expect(tutorialPanel.locator('img[src$="manual-02-01-upload-row-label-context.webp"]')).toHaveCount(1);
  await expect(tutorialPanel.locator('img[src$="manual-02-01-nine-file-upload-order.webp"]')).toHaveCount(0);
  await expect(tutorialPanel.locator('img[src$="manual-05-01-record-labels.webp"]')).toHaveCount(0);
  await expect(tutorialPanel.locator('img[src$="manual-07-01-orthogroup-preview.webp"]')).toHaveCount(1);
  await expect(tutorialPanel.locator('img[src$="manual-08-01-orthogroup-popup.webp"]')).toHaveCount(1);
  await expect(tutorialPanel.locator('img[src$="manual-09-01-files-tab.webp"]')).toHaveCount(0);
  await expect(tutorialPanel.locator('img[src$="manual-03-03-thread-threshold-settings.webp"]')).toHaveCount(1);
  for (let idx = 0; idx < await mediaImages.count(); idx += 1) {
    const image = mediaImages.nth(idx);
    await image.scrollIntoViewIfNeeded();
    await expect.poll(() => image.evaluate((element) => element.complete && element.naturalWidth > 0)).toBe(true);
  }

  await page.getByRole('tab', { name: 'Files' }).click();
  const filesPanel = page.getByRole('tabpanel', { name: 'Files' });
  await expect(filesPanel.getByText('MjeNMV.gb')).toBeVisible();
  await expect(filesPanel.getByText('MejoMJNV.gb')).toBeVisible();
  await expect(filesPanel.getByRole('link', { name: 'Session JSON' })).toBeVisible();

  expect(pageErrors).toEqual([]);
});

test('Gallery renders the Vibrio multi-record tutorial and media', async ({ page }) => {
  const pageErrors = [];
  page.on('pageerror', (error) => pageErrors.push(error.message));

  await page.goto(`${baseUrl}/gallery/#Vnig_TUMSAT-TG-2018`, { waitUntil: 'domcontentloaded' });
  await expect(page.getByRole('heading', { name: /Vibrio nigripulchritudo/i })).toBeVisible();

  await page.getByRole('tab', { name: 'Tutorial' }).click();
  await expect(
    page.getByRole('heading', {
      name: 'Reproduce the Vibrio nigripulchritudo multi-record circular genome in the web app'
    })
  ).toBeVisible();
  const tutorialPanel = page.getByRole('tabpanel', { name: 'Tutorial' });
  await expect(tutorialPanel.getByRole('heading', { name: 'Requirements' })).toHaveCount(0);
  await expect(tutorialPanel.getByRole('heading', { name: 'Color rule basics' })).toHaveCount(0);
  await expect(tutorialPanel.getByRole('heading', { name: 'Color rules' })).toHaveCount(0);
  await expect(tutorialPanel.getByRole('heading', { name: 'Post-generation edits' })).toHaveCount(0);
  await expect(tutorialPanel.getByText('Enable Multi-record canvas')).toBeVisible();
  await expect(tutorialPanel.getByText('Upload the GenBank file')).toBeVisible();
  await expect(tutorialPanel.getByText('Set the circular layout')).toBeVisible();
  await expect(
    tutorialPanel.getByRole('row', {
      name: 'Track type tuckin'
    })
  ).toBeVisible();
  await expect(
    tutorialPanel.getByRole('row', {
      name: 'Record rows #1@1, #2@1, #3@2, #4@2, #5@2, #6@2'
    })
  ).toBeVisible();
  const mediaImages = tutorialPanel.getByRole('img');
  await expect(mediaImages).toHaveCount(11);
  await expect(tutorialPanel.locator('img[src$="manual-03-00-circular-layout-settings.webp"]')).toHaveCount(1);
  await expect(tutorialPanel.locator('img[src$="manual-03-03-selected-features.webp"]')).toHaveCount(1);
  await expect(tutorialPanel.locator('img[src$="manual-04-00-multirecord-toggle.webp"]')).toHaveCount(1);
  const multiRecordControlsImage = tutorialPanel.locator('img[src$="manual-04-01-multirecord-canvas.webp"]');
  await expect(tutorialPanel.locator('img[src$="manual-06-01-multirecord-preview.webp"]')).toHaveCount(1);
  await expect(tutorialPanel.locator('img[src$="manual-07-01-files-tab.webp"]')).toHaveCount(0);
  await expect(tutorialPanel.locator('img[src$="manual-08-01-feature-popup.webp"]')).toHaveCount(1);
  await expect(multiRecordControlsImage).toHaveCount(1);
  await expect(tutorialPanel.locator('img[src*="/post-"]')).toHaveCount(0);
  await expect(tutorialPanel.locator('img[src$="manual-04-02-record-positions.webp"]')).toHaveCount(0);
  for (let idx = 0; idx < await mediaImages.count(); idx += 1) {
    const image = mediaImages.nth(idx);
    await image.scrollIntoViewIfNeeded();
    await expect.poll(() => image.evaluate((element) => element.complete && element.naturalWidth > 0)).toBe(true);
  }

  await page.getByRole('tab', { name: 'Files' }).click();
  const filesPanel = page.getByRole('tabpanel', { name: 'Files' });
  await expect(filesPanel.getByText('GCF_015097735.1_ASM1509773v1_genomic.gbff')).toBeVisible();
  await expect(filesPanel.getByRole('link', { name: 'Session JSON' })).toBeVisible();

  expect(pageErrors).toEqual([]);
});

test('Gallery keeps command panel usable when a tutorial fetch fails', async ({ page }) => {
  await page.route('**/tutorials/Vnig_TUMSAT-TG-2018.json', (route) =>
    route.fulfill({
      status: 404,
      contentType: 'application/json; charset=utf-8',
      body: '{}'
    })
  );

  await page.goto(`${baseUrl}/gallery/#Vnig_TUMSAT-TG-2018`, { waitUntil: 'domcontentloaded' });
  await expect(page.getByRole('heading', { name: /Vibrio nigripulchritudo/i })).toBeVisible();

  await page.getByRole('tab', { name: 'Tutorial' }).click();
  await expect(page.getByText('Could not load tutorial (404).')).toBeVisible();

  await page.getByRole('tab', { name: 'Command' }).click();
  await expect(page.getByText('gbdraw circular -o Vnig_TUMSAT-TG-2018')).toBeVisible();
});

test('Gallery shows an inline fallback when tutorial media fails to load', async ({ page }) => {
  await page.route('**/media/hepatoplasmataceae_collinear/manual-01-01-linear-mode.webp', (route) =>
    route.fulfill({
      status: 404,
      contentType: 'image/webp',
      body: ''
    })
  );

  await page.goto(`${baseUrl}/gallery/#hepatoplasmataceae_collinear`, { waitUntil: 'domcontentloaded' });
  await page.getByRole('tab', { name: 'Tutorial' }).click();
  await expect(
    page.getByRole('heading', { name: 'Reproduce the Hepatoplasmataceae collinear analysis in the web app' })
  ).toBeVisible();
  await page.locator('#tutorial-panel .tutorial-media').first().scrollIntoViewIfNeeded();
  await expect(
    page.getByText('Media unavailable: ./media/hepatoplasmataceae_collinear/manual-01-01-linear-mode.webp')
  ).toBeVisible();
  await expect(page.getByText('Web app steps')).toBeVisible();
});

test('Gallery renders operation media and keeps it inside a mobile viewport', async ({ page }) => {
  await page.setViewportSize({ width: 390, height: 844 });
  await page.route('**/tutorials/hepatoplasmataceae_collinear.json', (route) => {
    route.fulfill({
      status: 200,
      contentType: 'application/json; charset=utf-8',
      body: readTutorialWithOperationMedia()
    });
  });

  await page.goto(`${baseUrl}/gallery/#hepatoplasmataceae_collinear`, { waitUntil: 'domcontentloaded' });
  await page.getByRole('tab', { name: 'Tutorial' }).click();

  const operation = page.locator('#tutorial-panel .tutorial-operation').first();
  await expect(operation.getByText('Select Linear mode')).toBeVisible();
  await expect(operation.getByText('Use Linear mode before adding the five GenBank inputs.')).toBeVisible();
  const image = operation.getByRole('img');
  await expect(image).toHaveAttribute('src', /manual-03-01-open-pairwise\.webp$/);
  await expect(
    operation.getByText('Operation-level media uses the same image path as the current tutorial.')
  ).toBeVisible();
  await image.scrollIntoViewIfNeeded();
  await expect.poll(() => image.evaluate((element) => element.complete && element.naturalWidth > 0)).toBe(true);

  const overflow = await image.evaluate((element) => {
    const rect = element.getBoundingClientRect();
    const parentRect = element.parentElement.getBoundingClientRect();
    return rect.width > parentRect.width + 1 || rect.right > document.documentElement.clientWidth + 1;
  });
  expect(overflow).toBe(false);
});

test('Gallery copy link button copies the selected sample URL', async ({ page }) => {
  await page.addInitScript(() => {
    let copiedText = '';
    Object.defineProperty(window, 'isSecureContext', { value: true, configurable: true });
    Object.defineProperty(navigator, 'clipboard', {
      value: {
        writeText: async (text) => {
          copiedText = String(text);
          window.__copiedGalleryLink = copiedText;
        },
        readText: async () => copiedText
      },
      configurable: true
    });
  });

  await page.goto(`${baseUrl}/gallery/#HmmtDNA_ATskew`, { waitUntil: 'domcontentloaded' });
  await page.getByRole('button', { name: 'Copy' }).click();

  await expect(page.getByRole('button', { name: 'Copied' })).toBeVisible();
  await expect.poll(() => page.evaluate(() => window.__copiedGalleryLink)).toBe(`${baseUrl}/gallery/#HmmtDNA_ATskew`);
});

test('Gallery tab controls support keyboard navigation', async ({ page }) => {
  await page.goto(`${baseUrl}/gallery/#Vnig_TUMSAT-TG-2018`, { waitUntil: 'domcontentloaded' });

  const previewTab = page.getByRole('tab', { name: 'Preview' });
  const tutorialTab = page.getByRole('tab', { name: 'Tutorial' });
  const commandTab = page.getByRole('tab', { name: 'Command' });
  const filesTab = page.getByRole('tab', { name: 'Files' });

  await previewTab.focus();
  await page.keyboard.press('ArrowRight');
  await expect(tutorialTab).toBeFocused();
  await expect(tutorialTab).toHaveAttribute('aria-selected', 'true');
  await expect(page.getByRole('tabpanel', { name: 'Tutorial' })).toBeVisible();

  await page.keyboard.press('ArrowRight');
  await expect(commandTab).toBeFocused();
  await expect(commandTab).toHaveAttribute('aria-selected', 'true');

  await page.keyboard.press('End');
  await expect(filesTab).toBeFocused();
  await expect(filesTab).toHaveAttribute('aria-selected', 'true');
  await expect(page.locator('#files-panel')).not.toHaveAttribute('hidden', '');
  await expect.poll(() => page.locator('#preview-panel').evaluate((panel) => getComputedStyle(panel).display)).toBe('none');

  await page.keyboard.press('Home');
  await expect(previewTab).toBeFocused();
  await expect(previewTab).toHaveAttribute('aria-selected', 'true');
  await expect.poll(() => page.locator('#files-panel').evaluate((panel) => getComputedStyle(panel).display)).toBe('none');

  await page.keyboard.press('ArrowLeft');
  await expect(filesTab).toBeFocused();
  await expect(filesTab).toHaveAttribute('aria-selected', 'true');
});

test('Gallery tutorial media fits a mobile viewport', async ({ page }) => {
  await page.setViewportSize({ width: 390, height: 844 });
  await page.goto(`${baseUrl}/gallery/#hepatoplasmataceae_collinear`, { waitUntil: 'domcontentloaded' });
  await page.getByRole('tab', { name: 'Tutorial' }).click();

  const mediaImages = page.locator('#tutorial-panel .tutorial-media img');
  await expect(mediaImages).toHaveCount(12);
  await mediaImages.last().scrollIntoViewIfNeeded();
  const overflowingImages = await mediaImages.evaluateAll((images) =>
    images
      .map((image) => {
        const rect = image.getBoundingClientRect();
        const parentRect = image.parentElement.getBoundingClientRect();
        return rect.width > parentRect.width + 1 || rect.right > document.documentElement.clientWidth + 1;
      })
      .filter(Boolean)
  );
  expect(overflowingImages).toEqual([]);
});

test('Gallery respects tutorialStatus draft and planned states', async ({ page }) => {
  let tutorialFetchCount = 0;
  await page.route('**/gallery/examples.json', (route) => {
    const examples = readGalleryExamples().map((entry) => {
      if (entry.id === 'hepatoplasmataceae_collinear') {
        return { ...entry, tutorialStatus: 'planned' };
      }
      if (entry.id === 'Vnig_TUMSAT-TG-2018') {
        return { ...entry, tutorialStatus: 'draft' };
      }
      return entry;
    });
    route.fulfill({
      status: 200,
      contentType: 'application/json; charset=utf-8',
      body: JSON.stringify(examples)
    });
  });
  await page.route('**/tutorials/hepatoplasmataceae_collinear.json', (route) => {
    tutorialFetchCount += 1;
    route.fulfill({
      status: 200,
      contentType: 'application/json; charset=utf-8',
      body: readTutorialJson('hepatoplasmataceae_collinear.json')
    });
  });

  await page.goto(`${baseUrl}/gallery/#hepatoplasmataceae_collinear`, { waitUntil: 'domcontentloaded' });
  await page.getByRole('tab', { name: 'Tutorial' }).click();
  await expect(page.getByRole('heading', { name: 'Tutorial planned' })).toBeVisible();
  await expect(page.getByText('A step-by-step web tutorial is not available')).toBeVisible();
  expect(tutorialFetchCount).toBe(0);

  await page.locator('[data-sample-id="Vnig_TUMSAT-TG-2018"]').click();
  await expect(page.getByRole('heading', { name: 'Tutorial draft' })).toBeVisible();
  await expect(page.getByText('not published in the public Gallery yet')).toBeVisible();
});

test('Gallery can preview draft tutorials when explicitly enabled', async ({ page }) => {
  await page.route('**/gallery/examples.json', (route) => {
    const examples = readGalleryExamples().map((entry) =>
      entry.id === 'Vnig_TUMSAT-TG-2018' ? { ...entry, tutorialStatus: 'draft' } : entry
    );
    route.fulfill({
      status: 200,
      contentType: 'application/json; charset=utf-8',
      body: JSON.stringify(examples)
    });
  });

  await page.goto(`${baseUrl}/gallery/?showDraftTutorials=1#Vnig_TUMSAT-TG-2018`, {
    waitUntil: 'domcontentloaded'
  });
  await page.getByRole('tab', { name: 'Tutorial' }).click();
  await expect(
    page.getByRole('heading', {
      name: 'Reproduce the Vibrio nigripulchritudo multi-record circular genome in the web app'
    })
  ).toBeVisible();
});

test('Gallery ignores stale tutorial fetch results after sample changes', async ({ page }) => {
  let releaseCollinearTutorial;
  let resolveCollinearRequestFinished;
  const collinearRequestStarted = new Promise((resolveStarted) => {
    page.route('**/tutorials/hepatoplasmataceae_collinear.json', async (route) => {
      resolveStarted();
      await new Promise((resolveRelease) => {
        releaseCollinearTutorial = resolveRelease;
      });
      await route.fulfill({
        status: 200,
        contentType: 'application/json; charset=utf-8',
        body: readTutorialJson('hepatoplasmataceae_collinear.json')
      });
      resolveCollinearRequestFinished();
    });
  });
  const collinearRequestFinished = new Promise((resolveFinished) => {
    resolveCollinearRequestFinished = resolveFinished;
  });

  await page.goto(`${baseUrl}/gallery/#hepatoplasmataceae_collinear`, { waitUntil: 'domcontentloaded' });
  await page.getByRole('tab', { name: 'Tutorial' }).click();
  await collinearRequestStarted;
  await expect(page.getByText('Loading tutorial for')).toBeVisible();

  await page.locator('[data-sample-id="Vnig_TUMSAT-TG-2018"]').click();
  const tutorialPanel = page.getByRole('tabpanel', { name: 'Tutorial' });
  await expect(
    tutorialPanel.getByRole('heading', {
      name: 'Reproduce the Vibrio nigripulchritudo multi-record circular genome in the web app'
    })
  ).toBeVisible();

  releaseCollinearTutorial();
  await collinearRequestFinished;
  await expect(
    tutorialPanel.getByRole('heading', {
      name: 'Reproduce the Vibrio nigripulchritudo multi-record circular genome in the web app'
    })
  ).toBeVisible();
  await expect(
    tutorialPanel.getByRole('heading', {
      name: 'Reproduce the Hepatoplasmataceae collinear analysis in the web app'
    })
  ).toHaveCount(0);
});

test('Gallery mobile controls and tutorial text do not overlap or overflow', async ({ page }) => {
  await page.setViewportSize({ width: 390, height: 844 });
  await page.goto(`${baseUrl}/gallery/#Vnig_TUMSAT-TG-2018`, { waitUntil: 'domcontentloaded' });
  await page.getByRole('tab', { name: 'Tutorial' }).click();
  await expect(
    page.getByRole('heading', {
      name: 'Reproduce the Vibrio nigripulchritudo multi-record circular genome in the web app'
    })
  ).toBeVisible();

  const layoutIssues = await page.evaluate(() => {
    const overlaps = [];
    const checkedParents = ['.viewer-actions', '.tab-bar'];
    checkedParents.forEach((selector) => {
      const parent = document.querySelector(selector);
      if (!parent) return;
      const children = Array.from(parent.children).filter((child) => !child.hidden);
      for (let leftIndex = 0; leftIndex < children.length; leftIndex += 1) {
        for (let rightIndex = leftIndex + 1; rightIndex < children.length; rightIndex += 1) {
          const leftRect = children[leftIndex].getBoundingClientRect();
          const rightRect = children[rightIndex].getBoundingClientRect();
          const hasOverlap =
            leftRect.left < rightRect.right - 1 &&
            leftRect.right > rightRect.left + 1 &&
            leftRect.top < rightRect.bottom - 1 &&
            leftRect.bottom > rightRect.top + 1;
          if (hasOverlap) {
            overlaps.push(`${selector}: ${children[leftIndex].textContent} / ${children[rightIndex].textContent}`);
          }
        }
      }
    });

    const viewportWidth = document.documentElement.clientWidth;
    const overflowing = Array.from(
      document.querySelectorAll(
        '.viewer-actions .action-button, .tab-button, #selected-title, .tutorial-step__title, .tutorial-step__body, .download-title, .file-title'
      )
    )
      .map((element) => {
        const rect = element.getBoundingClientRect();
        return rect.left < -1 || rect.right > viewportWidth + 1
          ? `${element.tagName.toLowerCase()}.${element.className || element.id}`
          : '';
      })
      .filter(Boolean);

    return {
      overlaps,
      overflowing,
      documentOverflow:
        document.documentElement.scrollWidth > document.documentElement.clientWidth + 1
    };
  });

  expect(layoutIssues).toEqual({
    overlaps: [],
    overflowing: [],
    documentOverflow: false
  });
});
