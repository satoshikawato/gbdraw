const DEFAULT_SAMPLE_ID = 'lambda-phage-linear';

const sampleList = document.querySelector('#sample-list');
const sampleCount = document.querySelector('#sample-count');
const selectedTitle = document.querySelector('#selected-title');
const selectedDescription = document.querySelector('#selected-description');
const selectedTags = document.querySelector('#selected-tags');
const fileSize = document.querySelector('#file-size');
const frame = document.querySelector('#demo-frame');
const frameLoading = document.querySelector('#frame-loading');
const commandBlock = document.querySelector('#command-block code');
const interactiveStep = document.querySelector('#interactive-step');
const openLink = document.querySelector('#open-link');
const downloadLink = document.querySelector('#download-link');
const copyLink = document.querySelector('#copy-link');

let examples = [];
let selectedId = '';
let copyStatusTimer = 0;

const normalizeHash = () => decodeURIComponent(window.location.hash.replace(/^#/, '')).trim();

const getAbsoluteSampleUrl = (sample) => new URL(sample.svg, window.location.href).href;

const getCurrentSampleUrl = (sample) => {
  const url = new URL(window.location.href);
  url.hash = sample.id;
  return url.href;
};

const clearChildren = (element) => {
  while (element.firstChild) {
    element.removeChild(element.firstChild);
  }
};

const appendText = (parent, tagName, className, text) => {
  const element = document.createElement(tagName);
  if (className) element.className = className;
  element.textContent = text;
  parent.appendChild(element);
  return element;
};

const renderTags = (parent, tags) => {
  clearChildren(parent);
  tags.forEach((tag) => appendText(parent, 'span', 'tag', tag));
};

const setLoading = (sample) => {
  frameLoading.hidden = false;
  frame.title = `Interactive gbdraw SVG demo: ${sample.title}`;
};

const updateActionLinks = (sample) => {
  const href = getAbsoluteSampleUrl(sample);
  openLink.href = href;
  downloadLink.href = href;
  downloadLink.setAttribute('download', sample.svg.split('/').pop() || `${sample.id}.svg`);
};

const updateCopyButton = (text) => {
  window.clearTimeout(copyStatusTimer);
  copyLink.textContent = text;
  copyStatusTimer = window.setTimeout(() => {
    copyLink.textContent = 'Copy link';
  }, 1600);
};

const copyText = async (text) => {
  if (navigator.clipboard?.writeText && window.isSecureContext) {
    await navigator.clipboard.writeText(text);
    return;
  }
  const buffer = document.createElement('textarea');
  buffer.className = 'copy-buffer';
  buffer.value = text;
  document.body.appendChild(buffer);
  buffer.select();
  document.execCommand('copy');
  buffer.remove();
};

const renderSampleList = () => {
  clearChildren(sampleList);
  sampleCount.textContent = `${examples.length} examples`;

  examples.forEach((sample, index) => {
    const button = document.createElement('button');
    button.type = 'button';
    button.className = 'sample-card';
    button.setAttribute('role', 'listitem');
    button.setAttribute('aria-pressed', sample.id === selectedId ? 'true' : 'false');
    button.dataset.sampleId = sample.id;

    const thumbnail = document.createElement('img');
    thumbnail.className = 'sample-card__thumb';
    thumbnail.src = sample.thumbnail;
    thumbnail.alt = '';
    thumbnail.loading = index === 0 ? 'eager' : 'lazy';
    button.appendChild(thumbnail);

    const body = document.createElement('span');
    body.className = 'sample-card__body';
    appendText(body, 'span', 'sample-card__title', sample.title);
    appendText(body, 'span', 'sample-card__description', sample.description);
    const tagRow = document.createElement('span');
    tagRow.className = 'tag-row';
    renderTags(tagRow, sample.tags);
    body.appendChild(tagRow);
    appendText(body, 'span', 'sample-card__size', sample.fileSizeLabel);
    button.appendChild(body);

    button.addEventListener('click', () => selectSample(sample.id));
    sampleList.appendChild(button);
  });
};

const updatePressedState = () => {
  sampleList.querySelectorAll('.sample-card').forEach((button) => {
    button.setAttribute('aria-pressed', button.dataset.sampleId === selectedId ? 'true' : 'false');
  });
};

const selectSample = (id, { updateHash = true } = {}) => {
  const sample = examples.find((entry) => entry.id === id) || examples[0];
  if (!sample) return;

  selectedId = sample.id;
  selectedTitle.textContent = sample.title;
  selectedDescription.textContent = sample.description;
  renderTags(selectedTags, sample.tags);
  fileSize.textContent = sample.fileSizeLabel;
  commandBlock.textContent = sample.command;
  interactiveStep.textContent = sample.interactiveStep || '';
  updateActionLinks(sample);
  updatePressedState();
  setLoading(sample);

  if (updateHash && normalizeHash() !== sample.id) {
    window.history.replaceState(null, '', `#${encodeURIComponent(sample.id)}`);
  }

  window.requestAnimationFrame(() => {
    frame.src = sample.svg;
  });
};

const showError = (message) => {
  selectedTitle.textContent = 'Gallery unavailable';
  selectedDescription.textContent = message;
  frameLoading.textContent = message;
};

frame.addEventListener('load', () => {
  frameLoading.hidden = true;
});

copyLink.addEventListener('click', async () => {
  const sample = examples.find((entry) => entry.id === selectedId);
  if (!sample) return;
  try {
    await copyText(getCurrentSampleUrl(sample));
    updateCopyButton('Copied');
  } catch (error) {
    updateCopyButton('Copy failed');
  }
});

window.addEventListener('hashchange', () => {
  selectSample(normalizeHash(), { updateHash: false });
});

fetch('./examples.json', { cache: 'no-cache' })
  .then((response) => {
    if (!response.ok) throw new Error(`Could not load examples.json (${response.status})`);
    return response.json();
  })
  .then((data) => {
    examples = Array.isArray(data) ? data : [];
    if (!examples.length) throw new Error('No gallery examples are configured.');
    const initialId = normalizeHash() || DEFAULT_SAMPLE_ID;
    selectedId = initialId;
    renderSampleList();
    selectSample(initialId);
  })
  .catch((error) => {
    showError(error?.message || 'Could not load gallery examples.');
  });
