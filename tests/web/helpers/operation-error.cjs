const { expect } = require('@playwright/test');

// Change one real render input; the Worker, Python adapter and callers remain live.
const failNextNativeRender = (page) => page.evaluate(() => {
  const send = Worker.prototype.postMessage;
  Worker.prototype.postMessage = function (message, ...args) {
    if (message.type === 'run') {
      window.__S03_NATIVE_RENDER_SENT__ = true;
      Worker.prototype.postMessage = send;
      const bytes = new TextEncoder().encode('feature_type\tqualifier_key\tvalue\tcolor\tcaption\nCDS\tproduct\t😀[PRIVATE_PATTERN_SENTINEL\t#ff0000\t\n');
      const request = structuredClone(message.payload.request);
      request.diagramOptions.colors.colorTable = { resourceId: 's03-error-colors', representation: 'canonicalTsv' };
      message = { ...message, payload: { ...message.payload, request,
        resourceManifest: [...message.payload.resourceManifest, { resourceId: 's03-error-colors',
          cacheToken: 'render-resource-999999', name: 'PRIVATE_FILE_SENTINEL.tsv', kind: 'canonical-tsv',
          type: 'text/tab-separated-values', size: bytes.byteLength, lastModified: 0 }],
        stagedResources: [...message.payload.stagedResources, { resourceId: 's03-error-colors',
          cacheToken: 'render-resource-999999', bytes: bytes.buffer }] } };
    }
    return send.call(this, message, ...args);
  };
});

const inspectSafeDetails = async (page, alert, expected) => {
  const details = alert.locator('details');
  expect(await details.evaluate(element => element.open)).toBe(false);
  await expect(alert.getByRole('button', { name: 'Copy diagnostics', exact: true })).toBeHidden();
  await page.evaluate(() => {
    window.__S03_COPIED__ = [];
    Object.defineProperty(navigator, 'clipboard', { configurable: true, value: {
      writeText: async text => { window.__S03_COPIED__.push(text); }
    } });
  });
  expect(await page.evaluate(() => window.__S03_COPIED__)).toEqual([]);
  const summary = details.locator('summary');
  await summary.scrollIntoViewIfNeeded();
  await summary.focus();
  await summary.press('Enter');
  const diagnostics = alert.getByRole('textbox', { name: 'Safe diagnostics' });
  await expect(diagnostics).toBeVisible();
  await expect(diagnostics).toHaveAttribute('readonly', '');
  const text = await diagnostics.inputValue();
  expect(text.length).toBeLessThanOrEqual(4012);
  expect(text).toContain(`Code: ${expected.code}`);
  expect(text).toContain(`Operation: ${expected.operation}`);
  expect(text).toContain(`Stage: ${expected.stage}`);
  expect(text).not.toMatch(/PRIVATE_|Traceback|stdout|stderr|😀|<svg/);
  const copy = alert.getByRole('button', { name: 'Copy diagnostics', exact: true });
  await diagnostics.focus();
  await page.keyboard.press('Tab');
  await expect(copy).toBeFocused();
  await copy.press('Enter');
  await expect(alert.getByRole('status')).toHaveText('Diagnostics copied.');
  expect(await page.evaluate(() => window.__S03_COPIED__)).toEqual([text]);
  await page.evaluate(() => Object.defineProperty(navigator, 'clipboard', { configurable: true, value: undefined }));
  await copy.press('Enter');
  await expect(alert.getByRole('status')).toContainText('copy manually');
  await expect(diagnostics).toHaveValue(text);
  await page.evaluate(() => Object.defineProperty(navigator, 'clipboard', { configurable: true,
    value: { writeText: async () => { throw new Error('PRIVATE_CLIPBOARD_SENTINEL'); } } }));
  await copy.press('Enter');
  await expect(alert.getByRole('status')).toContainText('copy manually');
  const select = alert.getByRole('button', { name: 'Select diagnostics', exact: true });
  await select.scrollIntoViewIfNeeded();
  await select.focus();
  await select.press('Enter');
  await expect(diagnostics).toBeFocused();
  expect(await diagnostics.evaluate(element => element.selectionEnd - element.selectionStart)).toBe(text.length);
  const box = await copy.boundingBox();
  expect(box.x).toBeGreaterThanOrEqual(0);
  expect(box.x + box.width).toBeLessThanOrEqual(page.viewportSize().width);
  expect(await alert.textContent()).not.toMatch(/PRIVATE_|Traceback|stdout|stderr|😀|<svg/);
};
module.exports = { failNextNativeRender, inspectSafeDetails };
