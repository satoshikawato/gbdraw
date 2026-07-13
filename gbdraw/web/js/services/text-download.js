export const downloadTextFile = (
  filename,
  text,
  type = 'text/tab-separated-values;charset=utf-8'
) => {
  const url = URL.createObjectURL(new Blob([text], { type }));
  const link = document.createElement('a');
  link.href = url;
  link.download = filename;
  link.addEventListener('click', (event) => event.stopPropagation(), { once: true });
  link.click();
  URL.revokeObjectURL(url);
};
