// @ts-check
/** @param {unknown} text */
export const copyTextToClipboard = async (text) => {
  const value = String(text ?? '');
  if (globalThis.isSecureContext && navigator.clipboard?.writeText) {
    try {
      await navigator.clipboard.writeText(value);
      return;
    } catch {
      // Clipboard permission can still be denied in a secure context.
    }
  }

  const textarea = document.createElement('textarea');
  textarea.value = value;
  textarea.setAttribute('readonly', '');
  textarea.style.position = 'fixed';
  textarea.style.left = '-9999px';
  textarea.style.top = '0';
  document.body.appendChild(textarea);
  textarea.select();
  try {
    if (!document.execCommand('copy')) throw new Error('Clipboard copy failed.');
  } finally {
    textarea.remove();
  }
};

/**
 * Copies text that an asynchronous call makes. A click allows a clipboard write only while it is current, so the
 * write starts at once with a `ClipboardItem` that waits for the text; without one (or when that write is refused),
 * the text is copied once it is ready. A failure of the call itself is the caller's to report.
 * @param {Promise<string>} pending
 * @returns {Promise<string>} The copied text.
 */
export const copyPendingTextToClipboard = async (pending) => {
  const Item = globalThis.ClipboardItem;
  if (globalThis.isSecureContext && Item && navigator.clipboard?.write) {
    const blob = pending.then((text) => new Blob([text], { type: 'text/plain' }));
    try {
      await navigator.clipboard.write([new Item({ 'text/plain': blob })]);
      return await pending;
    } catch {
      // The call's own failure is thrown below; a refused write falls back to plain text.
    }
  }
  const text = await pending;
  await copyTextToClipboard(text);
  return text;
};
