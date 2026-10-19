// @ts-check
import { t } from '../utils/i18n.js';

/**
 * @param {Blob} blob
 * @param {string} filename
 */
export const downloadBlob = (blob, filename) => {
  const url = URL.createObjectURL(blob);
  const link = document.createElement('a');
  link.href = url;
  link.download = filename;
  link.addEventListener?.('click', (event) => event.stopPropagation(), { once: true });
  try {
    document.body?.appendChild?.(link);
    link.click();
  } finally {
    link.parentNode?.removeChild?.(link);
    URL.revokeObjectURL(url);
  }
};

/**
 * @param {string} filename
 * @param {BlobPart} text
 * @param {string} [type]
 */
export const downloadTextFile = (
  filename,
  text,
  type = 'text/tab-separated-values;charset=utf-8'
) => {
  downloadBlob(new Blob([text], { type }), filename);
};

/**
 * The app's choice dialog as a download port: resolves to the chosen action's id, or `cancel`.
 * @typedef {(request: {
 *   title: string, description?: string,
 *   actions: { id: string, label: string, tone?: 'primary' | 'secondary' | 'danger' }[]
 * }) => Promise<string>} DownloadChoicePort
 */

/**
 * A large download (a Session, LOSAT tables, CLI helper files) asks
 * first; it starts only when **Download** is chosen.
 * @param {DownloadChoicePort} ask @param {{ title: string, description: string }} question
 * @returns {Promise<boolean>} Whether to download.
 */
export const askBeforeLargeDownload = async (ask, { title, description }) => (await ask({
  title, description, actions: [{ id: 'download', label: t('Download'), tone: 'primary' }]
})) === 'download';
