export const downloadTextFile = (
  filename,
  text,
  type = 'text/tab-separated-values;charset=utf-8'
) => {
  const url = URL.createObjectURL(new Blob([text], { type }));
  const link = document.createElement('a');
  link.href = url;
  link.download = filename;
  document.body.appendChild(link);
  link.addEventListener('click', (event) => event.stopPropagation(), { once: true });
  link.click();
  document.body.removeChild(link);
  URL.revokeObjectURL(url);
};
