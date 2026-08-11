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

export const downloadTextFile = (
  filename,
  text,
  type = 'text/tab-separated-values;charset=utf-8'
) => {
  downloadBlob(new Blob([text], { type }), filename);
};
