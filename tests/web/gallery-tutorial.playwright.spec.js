const { test, expect } = require('@playwright/test');
const { createReadStream, existsSync, statSync } = require('node:fs');
const { createServer } = require('node:http');
const { extname, join, normalize, resolve, sep } = require('node:path');

const repoRoot = resolve(process.env.GBDRAW_REPO || process.cwd());
const webRoot = join(repoRoot, 'gbdraw/web');

const contentTypes = {
  '.css': 'text/css; charset=utf-8',
  '.html': 'text/html; charset=utf-8',
  '.js': 'text/javascript; charset=utf-8',
  '.json': 'application/json; charset=utf-8',
  '.png': 'image/png',
  '.svg': 'image/svg+xml',
  '.webp': 'image/webp'
};

let server;
let baseUrl;

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
  await expect(page.getByText('Manual rebuild')).toBeVisible();
  await expect(page.getByText('Use browser LOSAT')).toBeVisible();
  await expect(page.getByText('Use Orientation + identity for Color mode.')).toBeVisible();
  const tutorialPanel = page.getByRole('tabpanel', { name: 'Tutorial' });
  const mediaImages = tutorialPanel.getByRole('img');
  await expect(mediaImages).toHaveCount(3);
  await expect(mediaImages.first()).toHaveAttribute('src', /01-gallery-session\.webp$/);
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
  await expect(mediaImages).toHaveCount(3);
  await expect(mediaImages.nth(1)).toHaveAttribute('src', /02-orthogroup-preview\.webp$/);
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
  await expect(tutorialPanel.getByText('CDS / gene / dnaA / #d95f02 / DnaA')).toBeVisible();
  await expect(page.getByText('antiSMASH gene-kind color rules')).toBeVisible();
  await expect(page.getByText('CDS / gene_kind / biosynthetic$')).toBeVisible();
  await expect(page.getByText('Reverse complement: BGC0000713 only')).toBeVisible();
  const mediaImages = tutorialPanel.getByRole('img');
  await expect(mediaImages).toHaveCount(4);
  await expect(mediaImages.nth(2)).toHaveAttribute('src', /03-orthogroup-popup\.webp$/);
  await expect(mediaImages.nth(3)).toHaveAttribute('src', /04-feature-popup\.webp$/);
  for (let idx = 0; idx < await mediaImages.count(); idx += 1) {
    const image = mediaImages.nth(idx);
    await image.scrollIntoViewIfNeeded();
    await expect.poll(() => image.evaluate((element) => element.complete && element.naturalWidth > 0)).toBe(true);
  }

  await page.getByRole('tab', { name: 'Files' }).click();
  const filesPanel = page.getByRole('tabpanel', { name: 'Files' });
  await expect(filesPanel.getByText('BGC0000708.gbk')).toBeVisible();
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
  await expect(tutorialPanel.getByText('Ring width: 5')).toBeVisible();
  await expect(tutorialPanel.getByText('MG18PR-0187-N40S.fa')).toBeVisible();
  const mediaImages = tutorialPanel.getByRole('img');
  await expect(mediaImages).toHaveCount(4);
  await expect(mediaImages.nth(1)).toHaveAttribute('src', /03-conservation-rings\.webp$/);
  await expect(mediaImages.nth(2)).toHaveAttribute('src', /02-input-files\.webp$/);
  await expect(mediaImages.nth(3)).toHaveAttribute('src', /04-feature-popup\.webp$/);
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
  await expect(tutorialPanel.getByText('Dinucleotide: AT')).toBeVisible();
  await expect(tutorialPanel.getByText('Legend label: AT skew')).toBeVisible();
  const mediaImages = tutorialPanel.getByRole('img');
  await expect(mediaImages).toHaveCount(3);
  await expect(mediaImages.nth(1)).toHaveAttribute('src', /02-atskew-preview\.webp$/);
  await expect(mediaImages.nth(2)).toHaveAttribute('src', /03-feature-popup\.webp$/);
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
  await expect(tutorialPanel.getByText('protein blastp mode: Orthogroups', { exact: true })).toBeVisible();
  await expect(tutorialPanel.getByText('WSSV-like proteins').first()).toBeVisible();
  await expect(tutorialPanel.getByText('Use 32 threads only deliberately')).toBeVisible();
  const mediaImages = tutorialPanel.getByRole('img');
  await expect(mediaImages).toHaveCount(4);
  await expect(mediaImages.nth(1)).toHaveAttribute('src', /02-orthogroup-preview\.webp$/);
  await expect(mediaImages.nth(2)).toHaveAttribute('src', /03-orthogroup-popup\.webp$/);
  await expect(mediaImages.nth(3)).toHaveAttribute('src', /04-files\.webp$/);
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
  await expect(tutorialPanel.getByText('Enable Multi-record canvas')).toBeVisible();
  await expect(tutorialPanel.getByText('NZ_AP024087.1: chromosome 1')).toBeVisible();
  await expect(tutorialPanel.getByText('Record positions: #1@1, #2@1, #3@2, #4@2, #5@2, #6@2')).toBeVisible();
  const mediaImages = tutorialPanel.getByRole('img');
  await expect(mediaImages).toHaveCount(4);
  await expect(mediaImages.nth(1)).toHaveAttribute('src', /02-multirecord-preview\.webp$/);
  await expect(mediaImages.nth(2)).toHaveAttribute('src', /03-files\.webp$/);
  await expect(mediaImages.nth(3)).toHaveAttribute('src', /04-feature-popup\.webp$/);
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

  await page.keyboard.press('Home');
  await expect(previewTab).toBeFocused();
  await expect(previewTab).toHaveAttribute('aria-selected', 'true');

  await page.keyboard.press('ArrowLeft');
  await expect(filesTab).toBeFocused();
  await expect(filesTab).toHaveAttribute('aria-selected', 'true');
});

test('Gallery tutorial media fits a mobile viewport', async ({ page }) => {
  await page.setViewportSize({ width: 390, height: 844 });
  await page.goto(`${baseUrl}/gallery/#hepatoplasmataceae_collinear`, { waitUntil: 'domcontentloaded' });
  await page.getByRole('tab', { name: 'Tutorial' }).click();

  const mediaImages = page.locator('#tutorial-panel .tutorial-media img');
  await expect(mediaImages).toHaveCount(3);
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
