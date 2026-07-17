const paletteSelect = document.querySelector('#palette-select');
const previousButton = document.querySelector('#previous-palette');
const nextButton = document.querySelector('#next-palette');
const palettePosition = document.querySelector('#palette-position');
const paletteSwatches = document.querySelector('#palette-swatches');
const paletteCommand = document.querySelector('#palette-command');
const selectedPaletteName = document.querySelector('#selected-palette-name');
const paletteStatus = document.querySelector('#palette-status');
const palettePreview = document.querySelector('#palette-preview');

const ROLE_LABELS = {
  CDS: 'CDS',
  rRNA: 'rRNA',
  tRNA: 'tRNA',
  tmRNA: 'tmRNA',
  ncRNA: 'ncRNA',
  repeat_region: 'Repeat',
  skew_high: 'Skew high',
  skew_low: 'Skew low',
  gc_content: 'GC content'
};

let paletteData = null;
let paletteNames = [];
let diagramSvg = null;

const fetchText = async (url) => {
  const response = await fetch(url, { cache: 'no-cache' });
  if (!response.ok) throw new Error(`Could not load ${url} (${response.status})`);
  return response.text();
};

const fetchJson = async (url) => {
  const response = await fetch(url, { cache: 'no-cache' });
  if (!response.ok) throw new Error(`Could not load ${url} (${response.status})`);
  return response.json();
};

const selectedNameFromUrl = () => new URL(window.location.href).searchParams.get('palette') || '';

const colorForRole = (palette, role) => {
  const defaultPalette = paletteData.palettes[paletteData.defaultPalette] || {};
  return palette?.[role] || palette?.default || defaultPalette[role] || defaultPalette.default || '#d3d3d3';
};

const updateSwatches = (palette) => {
  paletteSwatches.replaceChildren();
  paletteData.roles.forEach((role) => {
    const color = colorForRole(palette, role);
    const item = document.createElement('div');
    item.className = 'palette-swatch';

    const chip = document.createElement('span');
    chip.className = 'palette-swatch__color';
    chip.style.backgroundColor = color;
    chip.setAttribute('aria-hidden', 'true');

    const text = document.createElement('span');
    text.className = 'palette-swatch__text';
    const label = document.createElement('strong');
    label.textContent = ROLE_LABELS[role] || role;
    const code = document.createElement('code');
    code.textContent = color;
    text.append(label, code);
    item.append(chip, text);
    paletteSwatches.appendChild(item);
  });
};

const applyPalette = (name, { updateUrl = true } = {}) => {
  if (!diagramSvg || !paletteData?.palettes?.[name]) return;
  const palette = paletteData.palettes[name];

  diagramSvg.querySelectorAll('[data-palette-key]').forEach((element) => {
    const role = element.dataset.paletteKey;
    const attribute = element.dataset.paletteAttribute || 'fill';
    element.setAttribute(attribute, colorForRole(palette, role));
  });

  const index = paletteNames.indexOf(name);
  paletteSelect.value = name;
  palettePosition.textContent = `${index + 1} / ${paletteNames.length}`;
  selectedPaletteName.textContent = name;
  paletteCommand.textContent = `gbdraw circular ... -p ${name}`;
  paletteStatus.textContent = `Showing ${name} on the same Circular SVG.`;
  updateSwatches(palette);

  if (updateUrl) {
    const url = new URL(window.location.href);
    url.searchParams.set('palette', name);
    window.history.replaceState(null, '', url);
  }
};

const stepPalette = (offset) => {
  const current = Math.max(0, paletteNames.indexOf(paletteSelect.value));
  const next = (current + offset + paletteNames.length) % paletteNames.length;
  applyPalette(paletteNames[next]);
  paletteSelect.focus();
};

const loadExplorer = async () => {
  try {
    const [data, svgText] = await Promise.all([
      fetchJson('./palettes.json'),
      fetchText('./circular.svg')
    ]);
    if (!data?.palettes || typeof data.palettes !== 'object') {
      throw new Error('Palette data is invalid.');
    }

    const parsed = new DOMParser().parseFromString(svgText, 'image/svg+xml');
    const parseError = parsed.querySelector('parsererror');
    const svg = parsed.documentElement;
    if (parseError || svg.localName !== 'svg') throw new Error('Circular palette SVG is invalid.');
    if (!svg.querySelector('[data-palette-key]')) {
      throw new Error('Circular palette SVG has no semantic color markers.');
    }

    paletteData = data;
    paletteNames = Object.keys(data.palettes);
    diagramSvg = document.importNode(svg, true);
    palettePreview.replaceChildren(diagramSvg);
    palettePreview.setAttribute('aria-busy', 'false');

    paletteSelect.replaceChildren();
    paletteNames.forEach((name) => {
      const option = document.createElement('option');
      option.value = name;
      option.textContent = name;
      paletteSelect.appendChild(option);
    });
    paletteSelect.disabled = false;
    previousButton.disabled = false;
    nextButton.disabled = false;

    const requested = selectedNameFromUrl();
    const initial = paletteNames.includes(requested) ? requested : data.defaultPalette;
    applyPalette(initial, { updateUrl: Boolean(requested) });
  } catch (error) {
    palettePreview.setAttribute('aria-busy', 'false');
    const message = document.createElement('p');
    message.className = 'palette-error';
    message.textContent = error?.message || 'Could not load the palette explorer.';
    palettePreview.replaceChildren(message);
    paletteStatus.textContent = 'Palette explorer unavailable.';
  }
};

paletteSelect.addEventListener('change', () => applyPalette(paletteSelect.value));
paletteSelect.addEventListener('keydown', (event) => {
  if (event.key === 'Home') applyPalette(paletteNames[0]);
  if (event.key === 'End') applyPalette(paletteNames[paletteNames.length - 1]);
});
previousButton.addEventListener('click', () => stepPalette(-1));
nextButton.addEventListener('click', () => stepPalette(1));
window.addEventListener('popstate', () => {
  const requested = selectedNameFromUrl();
  if (paletteNames.includes(requested)) applyPalette(requested, { updateUrl: false });
});

loadExplorer();
