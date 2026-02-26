(() => {
  const head = document.getElementsByTagName('head')[0];
  if (!head) return;

  const existing = document.querySelector('link[data-phosphor="regular"]');
  if (existing) return;

  const link = document.createElement('link');
  link.rel = 'stylesheet';
  link.type = 'text/css';
  link.href = './vendor/phosphor/regular/style.css';
  link.setAttribute('data-phosphor', 'regular');
  head.appendChild(link);
})();
