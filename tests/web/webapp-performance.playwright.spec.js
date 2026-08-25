const { test, expect } = require('@playwright/test');
const { openApp } = require('./helpers/app-lifecycle.cjs');

test.describe.configure({ mode: 'serial' });

test('the browser SVG sanitizer retains its malicious-input security profile', async ({ page }) => {
  await openApp(page, { waitForPalette: false });
  await page.waitForFunction(() => window.DOMPurify?.sanitize);
  const sanitized = await page.evaluate(async () => {
    const { sanitizeSvgContent } = await import('/gbdraw/web/js/services/svg-sanitization.js');
    return sanitizeSvgContent([
      '<svg xmlns="http://www.w3.org/2000/svg" onload="window.__svg_xss = 1">',
      '<script>window.__svg_xss = 2</script>',
      '<foreignObject><div xmlns="http://www.w3.org/1999/xhtml">unsafe</div></foreignObject>',
      '<a href="javascript:window.__svg_xss = 3"><rect width="10" height="10" /></a>',
      '<rect id="safe-feature" data-gbdraw-feature-id="safe-feature" ',
      'width="10" height="10" fill="#54bcf8" onclick="window.__svg_xss = 4" />',
      '</svg>'
    ].join(''));
  });
  expect(sanitized).not.toMatch(/<script/i);
  expect(sanitized).not.toMatch(/foreignObject/i);
  expect(sanitized).not.toMatch(/\sonload=/i);
  expect(sanitized).not.toMatch(/\sonclick=/i);
  expect(sanitized).not.toMatch(/javascript:/i);
  expect(sanitized).toContain('data-gbdraw-feature-id="safe-feature"');
  expect(await page.evaluate(() => window.__svg_xss)).toBeUndefined();
});
