const DEFAULT_SAMPLE_ID = 'Vnig_TUMSAT-TG-2018';
const COMMON_COLOR_RULE_GUIDE_PATH = './tutorials/_common-color-rule-guide.json';

const sampleList = document.querySelector('#sample-list');
const sampleCount = document.querySelector('#sample-count');
const selectedTitle = document.querySelector('#selected-title');
const selectedDescription = document.querySelector('#selected-description');
const selectedMeta = document.querySelector('#selected-meta');
const selectedTags = document.querySelector('#selected-tags');
const fileSize = document.querySelector('#file-size');
const frame = document.querySelector('#demo-frame');
const frameLoading = document.querySelector('#frame-loading');
const previewNote = document.querySelector('#preview-note');
const commandBlock = document.querySelector('#command-block code');
const commandKind = document.querySelector('#command-kind');
const commandNote = document.querySelector('#command-note');
const copyCommandButton = document.querySelector('#copy-command');
const interactiveStep = document.querySelector('#interactive-step');
const openLink = document.querySelector('#open-link');
const downloadLink = document.querySelector('#download-link');
const sessionLink = document.querySelector('#session-link');
const copyLink = document.querySelector('#copy-link');
const tabButtons = Array.from(document.querySelectorAll('[data-gallery-tab]'));
const tabPanels = {
  preview: document.querySelector('#preview-panel'),
  tutorial: document.querySelector('#tutorial-panel'),
  command: document.querySelector('#command-panel'),
  files: document.querySelector('#files-panel')
};
const tutorialContent = document.querySelector('#tutorial-content');
const filesContent = document.querySelector('#files-content');
const filesSummary = document.querySelector('#files-summary');

const TAB_HASHES = {
  preview: '',
  tutorial: 'Tutorial',
  command: 'Command',
  files: 'Files'
};

const normalizeHash = () => decodeURIComponent(window.location.hash.replace(/^#/, '')).trim();

const getTabFromHash = () => {
  const value = normalizeHash().toLowerCase();
  return Object.keys(TAB_HASHES).find((tabName) => tabName === value) || '';
};

const getSampleIdFromPath = () => {
  const match = window.location.pathname.match(/\/gallery\/([^/]+)\/?$/i);
  if (!match) return '';
  const sampleId = decodeURIComponent(match[1]).trim();
  return sampleId.includes('.') ? '' : sampleId;
};

const getLegacySampleIdFromHash = () => (getTabFromHash() ? '' : normalizeHash());

const getGalleryRootPath = () => {
  const match = window.location.pathname.match(/^(.*\/gallery)(?:\/[^/]*)?\/?$/i);
  return match?.[1] || '/gallery';
};

let examples = [];
let selectedId = '';
let activeTab = getTabFromHash() || 'preview';
let copyStatusTimer = 0;
let tutorialRequestToken = 0;
let mediaLightbox = null;
let commonColorRuleGuide = null;
let commonColorRuleGuidePromise = null;
const tutorialCache = new Map();

const getAbsoluteSampleUrl = (sample) => new URL(sample.svg, window.location.href).href;

const getSampleRouteUrl = (sample, tabName = activeTab) => {
  const url = new URL(window.location.href);
  url.pathname = `${getGalleryRootPath()}/${encodeURIComponent(sample.id)}`;
  url.hash = TAB_HASHES[tabName] || '';
  return url;
};

const getCurrentSampleUrl = (sample) => getSampleRouteUrl(sample).href;

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

const appendHtml = (parent, tagName, className, html) => {
  const element = document.createElement(tagName);
  if (className) element.className = className;
  element.innerHTML = html;
  parent.appendChild(element);
  return element;
};

const appendLink = (parent, className, label, href, { download = '' } = {}) => {
  const link = document.createElement('a');
  if (className) link.className = className;
  link.textContent = label;
  link.href = href;
  if (download) link.setAttribute('download', download);
  if (/^https?:\/\//i.test(href)) {
    link.target = '_blank';
    link.rel = 'noopener';
  }
  parent.appendChild(link);
  return link;
};

const asArray = (value) => (Array.isArray(value) ? value : []);

const asText = (value) => (typeof value === 'string' ? value.trim() : '');

const getTutorialStatus = (sample) => {
  const status = asText(sample?.tutorialStatus).toLowerCase();
  if (['ready', 'draft', 'planned'].includes(status)) return status;
  return sample?.tutorial ? 'ready' : 'planned';
};

const isDraftTutorialPreviewEnabled = () => {
  const params = new URLSearchParams(window.location.search);
  return params.get('showDraftTutorials') === '1' || params.get('draftTutorials') === '1';
};

const canLoadTutorial = (sample) => {
  if (!sample?.tutorial) return false;
  const status = getTutorialStatus(sample);
  return status === 'ready' || (status === 'draft' && isDraftTutorialPreviewEnabled());
};

const selectedSample = () => examples.find((entry) => entry.id === selectedId);
const sampleIsInteractive = (sample) => asText(sample?.svgType).toLowerCase() !== 'static';

const pluralize = (count, singular, plural = `${singular}s`) =>
  `${count} ${count === 1 ? singular : plural}`;

const appendList = (parent, items, className = 'tutorial-list') => {
  const values = asArray(items).map(asText).filter(Boolean);
  if (!values.length) return null;
  const list = document.createElement('ul');
  list.className = className;
  values.forEach((item) => appendText(list, 'li', '', item));
  parent.appendChild(list);
  return list;
};

const appendTable = (parent, tableData) => {
  if (!tableData || typeof tableData !== 'object' || Array.isArray(tableData)) return null;
  const columns = asArray(tableData.columns).map(asText).filter(Boolean);
  const rows = asArray(tableData.rows)
    .map((row) => asArray(row).map(asText))
    .filter((row) => row.some(Boolean));
  if (!columns.length || !rows.length) return null;

  const wrap = document.createElement('div');
  wrap.className = 'tutorial-table-wrap';
  const table = document.createElement('table');
  table.className = 'tutorial-table';

  const caption = asText(tableData.caption);
  if (caption) appendText(table, 'caption', '', caption);

  const thead = document.createElement('thead');
  const headerRow = document.createElement('tr');
  columns.forEach((column) => {
    const cell = appendText(headerRow, 'th', '', column);
    cell.scope = 'col';
  });
  thead.appendChild(headerRow);
  table.appendChild(thead);

  const tbody = document.createElement('tbody');
  rows.forEach((row) => {
    const tr = document.createElement('tr');
    columns.forEach((_, columnIndex) => {
      appendText(tr, 'td', '', row[columnIndex] || '');
    });
    tbody.appendChild(tr);
  });
  table.appendChild(tbody);
  wrap.appendChild(table);
  parent.appendChild(wrap);
  return table;
};

const appendSection = (parent, title) => {
  const section = document.createElement('section');
  section.className = 'tutorial-section';
  appendText(section, 'h3', '', title);
  parent.appendChild(section);
  return section;
};

const tutorialSourceLabel = (download) => {
  const parts = [];
  if (download.source) parts.push(download.source);
  if (download.accession) parts.push(download.accession);
  return parts.join(': ');
};

const renderDownloads = (parent, downloads) => {
  const entries = asArray(downloads);
  if (!entries.length) return;

  const list = document.createElement('div');
  list.className = 'download-list';
  entries.forEach((download) => {
    if (!download || typeof download !== 'object') return;

    const item = document.createElement('div');
    item.className = 'download-item';
    const title = document.createElement('div');
    title.className = 'download-title';
    appendText(
      title,
      'span',
      '',
      asText(download.targetName) || asText(download.label) || asText(download.accession) || 'Input file'
    );

    const sourceLabel = tutorialSourceLabel(download);
    if (sourceLabel) {
      if (download.href) {
        appendLink(title, 'download-source', sourceLabel, download.href);
      } else if (download.accession) {
        appendLink(
          title,
          'download-source',
          sourceLabel,
          `https://www.ncbi.nlm.nih.gov/nuccore/${encodeURIComponent(download.accession)}`
        );
      } else {
        appendText(title, 'span', 'download-source', sourceLabel);
      }
    }

    item.appendChild(title);
    const note = asText(download.note);
    if (note) appendText(item, 'p', 'download-note', note);
    list.appendChild(item);
  });
  parent.appendChild(list);
};

const mediaEntryList = (media) => {
  if (!media) return [];
  return Array.isArray(media) ? media : [media];
};

const mediaIsVideo = (entry, source) => {
  const type = asText(entry.type).toLowerCase();
  return type === 'video' || /\.(mp4|webm|ogg)(?:[?#].*)?$/i.test(source);
};

const ensureMediaLightbox = () => {
  if (mediaLightbox) return mediaLightbox;

  const dialog = document.createElement('dialog');
  dialog.className = 'tutorial-lightbox';
  dialog.setAttribute('aria-label', 'Full-size tutorial image');
  dialog.tabIndex = -1;

  const image = document.createElement('img');
  image.className = 'tutorial-lightbox__image';
  dialog.appendChild(image);
  document.body.appendChild(dialog);

  image.addEventListener('click', () => dialog.close());
  dialog.addEventListener('click', (event) => {
    if (event.target === dialog) dialog.close();
  });
  dialog.addEventListener('close', () => {
    image.removeAttribute('src');
    image.alt = '';
  });

  mediaLightbox = { dialog, image };
  return mediaLightbox;
};

const openMediaLightbox = (source, label) => {
  const lightbox = ensureMediaLightbox();
  if (lightbox.dialog.open) lightbox.dialog.close();
  lightbox.image.src = source;
  lightbox.image.alt = label || '';
  if (typeof lightbox.dialog.showModal === 'function') {
    lightbox.dialog.showModal();
  } else {
    lightbox.dialog.setAttribute('open', '');
  }
  lightbox.dialog.focus();
};

const renderMedia = (parent, media) => {
  mediaEntryList(media).forEach((entry) => {
    const mediaObject = typeof entry === 'string' ? { src: entry } : entry;
    if (!mediaObject || typeof mediaObject !== 'object') return;

    const source = asText(mediaObject.src) || asText(mediaObject.href);
    if (!source) return;

    const caption = asText(mediaObject.caption);
    const altText = asText(mediaObject.alt) || caption || '';
    const figure = document.createElement('figure');
    figure.className = 'tutorial-media';
    const mediaElement = document.createElement(mediaIsVideo(mediaObject, source) ? 'video' : 'img');
    mediaElement.src = source;
    if (mediaElement.tagName === 'VIDEO') {
      mediaElement.controls = true;
      mediaElement.preload = 'metadata';
    } else {
      mediaElement.loading = 'lazy';
      mediaElement.alt = altText;
    }
    mediaElement.addEventListener('error', () => {
      clearChildren(figure);
      appendText(figure, 'p', 'panel-status panel-status--error', `Media unavailable: ${source}`);
    });

    if (mediaElement.tagName === 'IMG') {
      const openButton = document.createElement('button');
      openButton.type = 'button';
      openButton.className = 'tutorial-media__button';
      openButton.setAttribute('aria-label', `Open full-size image: ${caption || altText || source}`);
      openButton.appendChild(mediaElement);
      openButton.addEventListener('click', () => {
        openMediaLightbox(source, caption || altText);
      });
      figure.appendChild(openButton);
    } else {
      figure.appendChild(mediaElement);
    }

    if (caption) appendText(figure, 'figcaption', '', caption);
    parent.appendChild(figure);
  });
};

const renderOperationList = (parent, operations) => {
  const entries = asArray(operations);
  if (!entries.length) return;
  const numbered = entries.length > 1;

  const list = document.createElement('div');
  list.className = 'tutorial-operation-list';

  entries.forEach((operation, index) => {
    const operationObject = typeof operation === 'string' ? { body: operation } : operation;
    if (!operationObject || typeof operationObject !== 'object') return;

    const item = document.createElement('section');
    item.className = 'tutorial-operation';

    const title = asText(operationObject.title);
    if (title) {
      const titleRow = document.createElement('div');
      titleRow.className = 'tutorial-operation__title';
      if (numbered) appendText(titleRow, 'span', 'tutorial-operation__number', `${index + 1}.`);
      appendText(titleRow, 'span', '', title);
      item.appendChild(titleRow);
    }

    const body = asText(operationObject.body);
    if (body) appendText(item, 'p', 'tutorial-operation__body', body);
    appendList(item, operationObject.items, 'tutorial-operation__list');
    appendTable(item, operationObject.table);

    const note = asText(operationObject.note);
    if (note) appendText(item, 'p', 'tutorial-operation__note', note);
    renderMedia(item, operationObject.media);

    if (item.childElementCount) list.appendChild(item);
  });

  if (list.childElementCount) parent.appendChild(list);
};

const renderStepSection = (parent, title, steps, { numbered = false } = {}) => {
  const entries = asArray(steps);
  if (!entries.length) return;

  const section = appendSection(parent, title);
  const list = document.createElement('div');
  list.className = 'tutorial-step-list';

  entries.forEach((step, index) => {
    const stepObject = typeof step === 'string' ? { body: step } : step;
    if (!stepObject || typeof stepObject !== 'object') return;

    const item = document.createElement('section');
    item.className = 'tutorial-step';
    const stepTitle = asText(stepObject.title);
    if (stepTitle || numbered) {
      const titleRow = document.createElement('div');
      titleRow.className = 'tutorial-step__title';
      if (numbered) appendText(titleRow, 'span', 'tutorial-step__number', `${index + 1}.`);
      appendText(titleRow, 'span', '', stepTitle || `Step ${index + 1}`);
      item.appendChild(titleRow);
    }

    const body = asText(stepObject.body);
    if (body) appendText(item, 'p', 'tutorial-step__body', body);
    appendList(item, stepObject.items);
    appendTable(item, stepObject.table);

    const note = asText(stepObject.note);
    if (note) appendText(item, 'p', 'tutorial-step__note', note);
    renderMedia(item, stepObject.media);
    renderOperationList(item, stepObject.operations);

    list.appendChild(item);
  });

  section.appendChild(list);
};

const renderTextListSection = (parent, title, items) => {
  if (!asArray(items).length) return;
  const section = appendSection(parent, title);
  appendList(section, items);
};

const renderRelated = (parent, related) => {
  const entries = asArray(related);
  if (!entries.length) return;

  const section = appendSection(parent, 'Related');
  const list = document.createElement('ul');
  list.className = 'tutorial-list';
  entries.forEach((entry) => {
    if (!entry || typeof entry !== 'object') return;
    const label = asText(entry.label);
    const href = asText(entry.href);
    if (!label || !href) return;
    const item = document.createElement('li');
    appendLink(item, '', label, href);
    list.appendChild(item);
  });
  section.appendChild(list);
};

const renderTutorialUnavailable = (sample) => {
  clearChildren(tutorialContent);
  const status = getTutorialStatus(sample);
  const title = plainTitle(sample);

  if (status === 'draft' && sample?.tutorial) {
    appendText(tutorialContent, 'h2', '', 'Tutorial draft');
    appendText(
      tutorialContent,
      'p',
      'tutorial-summary',
      `A draft tutorial exists for ${title}, but it is not published in the public Gallery yet. Enable draft preview with ?showDraftTutorials=1 when checking local builds.`
    );
    return;
  }

  appendText(tutorialContent, 'h2', '', status === 'planned' ? 'Tutorial planned' : 'Tutorial coming soon');
  appendText(
    tutorialContent,
    'p',
    'tutorial-summary',
    `A step-by-step web tutorial is not available for ${title} yet. Use the Preview and Session controls to inspect the completed diagram, or open the Command tab for the CLI command.`
  );
};

const renderTutorialLoading = (sample) => {
  clearChildren(tutorialContent);
  appendText(tutorialContent, 'p', 'panel-status', `Loading tutorial for ${plainTitle(sample)}...`);
};

const renderTutorialError = (message) => {
  clearChildren(tutorialContent);
  appendText(tutorialContent, 'p', 'panel-status panel-status--error', message);
};

const validateTutorial = (data, sample) => {
  if (!data || typeof data !== 'object' || Array.isArray(data)) {
    throw new Error('Tutorial JSON must be an object.');
  }
  if (asText(data.id) !== sample.id) {
    throw new Error(`Tutorial ID mismatch for ${sample.id}.`);
  }
  if (typeof data.version !== 'number') {
    throw new Error(`Tutorial ${sample.id} is missing a numeric version.`);
  }
  if (!Array.isArray(data.manualSteps)) {
    throw new Error(`Tutorial ${sample.id} is missing manualSteps.`);
  }
  ['manualSteps', 'colorRules', 'postGenerationEdits', 'losatTips', 'troubleshooting'].forEach(
    (sectionName) => {
      asArray(data[sectionName]).forEach((step, index) => {
        if (!step || typeof step !== 'object' || Array.isArray(step)) return;
        if (step.operations !== undefined && !Array.isArray(step.operations)) {
          throw new Error(`Tutorial ${sample.id} ${sectionName}[${index}] operations must be an array.`);
        }
      });
    }
  );
  return data;
};

const validateCommonColorRuleGuide = (data) => {
  if (!data || typeof data !== 'object' || Array.isArray(data)) {
    throw new Error('Common color rule guide JSON must be an object.');
  }
  const colorRules = asArray(data.colorRules);
  if (!colorRules.length) {
    throw new Error('Common color rule guide is missing colorRules.');
  }
  return {
    title: asText(data.title) || 'Color rule basics',
    colorRules
  };
};

const loadCommonColorRuleGuide = async () => {
  if (commonColorRuleGuide) return commonColorRuleGuide;
  if (!commonColorRuleGuidePromise) {
    commonColorRuleGuidePromise = fetch(COMMON_COLOR_RULE_GUIDE_PATH, { cache: 'no-cache' })
      .then((response) => {
        if (!response.ok) {
          throw new Error(`Could not load common color rule guide (${response.status}).`);
        }
        return response.json();
      })
      .then(validateCommonColorRuleGuide)
      .then((guide) => {
        commonColorRuleGuide = guide;
        return guide;
      })
      .catch((error) => {
        commonColorRuleGuidePromise = null;
        throw error;
      });
  }
  return commonColorRuleGuidePromise;
};

const hydrateTutorial = async (tutorial) => {
  if (!asArray(tutorial.colorRules).length) return tutorial;

  try {
    return {
      ...tutorial,
      commonColorRuleGuide: await loadCommonColorRuleGuide()
    };
  } catch (error) {
    console.warn('Could not load common color rule guide.', error);
    return tutorial;
  }
};

const loadTutorial = async (sample) => {
  if (!sample?.tutorial) return null;
  if (tutorialCache.has(sample.id)) return tutorialCache.get(sample.id);

  const response = await fetch(sample.tutorial, { cache: 'no-cache' });
  if (!response.ok) {
    throw new Error(`Could not load tutorial (${response.status}).`);
  }
  const tutorial = await hydrateTutorial(validateTutorial(await response.json(), sample));
  tutorialCache.set(sample.id, tutorial);
  return tutorial;
};

const renderTutorial = (tutorial, sample) => {
  if (!tutorial) {
    renderTutorialUnavailable(sample);
    return;
  }

  clearChildren(tutorialContent);
  appendText(tutorialContent, 'h2', '', asText(tutorial.title) || plainTitle(sample));

  const summary = asText(tutorial.summary);
  if (summary) appendText(tutorialContent, 'p', 'tutorial-summary', summary);

  renderTextListSection(tutorialContent, 'Requirements', tutorial.requirements);
  if (asArray(tutorial.downloads).length) {
    const downloads = appendSection(tutorialContent, 'Input downloads');
    renderDownloads(downloads, tutorial.downloads);
  }
  renderStepSection(tutorialContent, 'Reproduce the figure', tutorial.manualSteps, { numbered: true });
  renderStepSection(
    tutorialContent,
    asText(tutorial.commonColorRuleGuide?.title) || 'Color rule basics',
    tutorial.commonColorRuleGuide?.colorRules
  );
  renderStepSection(tutorialContent, 'Color rules', tutorial.colorRules);
  renderStepSection(tutorialContent, 'Post-generation edits', tutorial.postGenerationEdits);
  renderStepSection(tutorialContent, 'Running LOSAT in the browser', tutorial.losatTips);
  renderStepSection(tutorialContent, 'Troubleshooting', tutorial.troubleshooting);
  renderRelated(tutorialContent, tutorial.related);
};

const renderArtifactList = (parent, sample) => {
  const artifacts = [
    {
      label: sampleIsInteractive(sample) ? 'Interactive SVG' : 'SVG',
      href: sample.svg,
      note: sample.sourceOutput || ''
    },
    {
      label: sample.session?.endsWith('.gz') ? 'Session JSON (gzip)' : 'Session JSON',
      href: sample.session,
      note: sample.sourceSession || ''
    },
    {
      label: 'Thumbnail',
      href: sample.thumbnail,
      note: 'Gallery thumbnail image.'
    }
  ].filter((item) => item.href);

  if (!artifacts.length) return;
  const section = document.createElement('section');
  appendText(section, 'h3', '', 'Gallery artifacts');
  const list = document.createElement('div');
  list.className = 'file-list';
  artifacts.forEach((artifact) => {
    const item = document.createElement('div');
    item.className = 'file-item';
    const title = document.createElement('div');
    title.className = 'file-title';
    appendLink(title, '', artifact.label, artifact.href, {
      download: artifact.href.split('/').pop() || ''
    });
    item.appendChild(title);
    if (artifact.note) appendText(item, 'p', 'file-note', artifact.note);
    list.appendChild(item);
  });
  section.appendChild(list);
  parent.appendChild(section);
};

const renderFeatureSourceList = (parent, featureSources) => {
  const sources = asArray(featureSources).map(asText).filter(Boolean);
  if (!sources.length) return;

  const section = document.createElement('section');
  appendText(section, 'h3', '', 'Input filenames');
  const list = document.createElement('div');
  list.className = 'file-list';
  sources.forEach((source) => {
    const item = document.createElement('div');
    item.className = 'file-item';
    appendText(item, 'div', 'file-title', source);
    appendText(item, 'p', 'file-note', 'Required input file listed in the gallery metadata.');
    list.appendChild(item);
  });
  section.appendChild(list);
  parent.appendChild(section);
};

const renderFilesPanel = (sample, tutorial = null) => {
  clearChildren(filesContent);
  const featureSources = asArray(sample.featureSources);
  const downloads = asArray(tutorial?.downloads);
  filesSummary.textContent = pluralize(featureSources.length, 'input file');

  if (downloads.length) {
    const section = document.createElement('section');
    appendText(section, 'h3', '', 'Input downloads');
    renderDownloads(section, downloads);
    filesContent.appendChild(section);
  } else {
    renderFeatureSourceList(filesContent, featureSources);
  }

  renderArtifactList(filesContent, sample);

  const note = asText(sample.sourceNote);
  if (note) appendText(filesContent, 'p', 'file-note', note);
};

const requestTutorialForSample = async (sample, { renderTutorialPanel = false, renderFiles = false } = {}) => {
  tutorialRequestToken += 1;
  const requestToken = tutorialRequestToken;

  if (!canLoadTutorial(sample)) {
    if (renderTutorialPanel) renderTutorialUnavailable(sample);
    if (renderFiles) renderFilesPanel(sample);
    return;
  }

  if (renderTutorialPanel) renderTutorialLoading(sample);

  try {
    const tutorial = await loadTutorial(sample);
    if (requestToken !== tutorialRequestToken || selectedId !== sample.id) return;
    if (renderTutorialPanel) renderTutorial(tutorial, sample);
    if (renderFiles) renderFilesPanel(sample, tutorial);
  } catch (error) {
    if (requestToken !== tutorialRequestToken || selectedId !== sample.id) return;
    if (renderTutorialPanel) {
      renderTutorialError(error?.message || 'Could not load tutorial.');
    }
    if (renderFiles) renderFilesPanel(sample);
  }
};

const plainTitle = (sample) => {
  const template = document.createElement('template');
  template.innerHTML = sample.title || '';
  return template.content.textContent || sample.title || sample.id;
};

const renderTags = (parent, tags) => {
  clearChildren(parent);
  tags.forEach((tag) => appendText(parent, 'span', 'tag', tag));
};

const setLoading = (sample) => {
  frameLoading.hidden = false;
  frame.title = `${sampleIsInteractive(sample) ? 'Interactive' : 'Static'} gbdraw SVG: ${plainTitle(sample)}`;
};

const updateActionLinks = (sample) => {
  const href = getAbsoluteSampleUrl(sample);
  openLink.href = href;
  downloadLink.href = href;
  downloadLink.setAttribute('download', sample.svg.split('/').pop() || `${sample.id}.svg`);

  if (sample.session) {
    const sessionHref = new URL(sample.session, window.location.href).href;
    sessionLink.hidden = false;
    sessionLink.href = sessionHref;
    sessionLink.setAttribute(
      'download',
      sample.session.split('/').pop() || `${sample.id}.gbdraw-session.json`
    );
  } else {
    sessionLink.hidden = true;
    sessionLink.removeAttribute('href');
  }
};

const updateCopyButton = (text) => {
  window.clearTimeout(copyStatusTimer);
  copyLink.textContent = text;
  copyStatusTimer = window.setTimeout(() => {
    copyLink.textContent = 'Copy';
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
    appendHtml(body, 'span', 'sample-card__title', sample.title);
    if (sample.description) {
      appendText(body, 'span', 'sample-card__description', sample.description);
    }
    appendText(
      body,
      'span',
      'sample-card__meta',
      [sample.difficulty, sample.workflow].filter(Boolean).join(' · ')
    );
    appendText(
      body,
      'span',
      'sample-card__meta sample-card__meta--secondary',
      [sample.inputSummary, sample.estimatedTime].filter(Boolean).join(' · ')
    );
    const tagRow = document.createElement('span');
    tagRow.className = 'tag-row';
    renderTags(tagRow, sample.tags);
    body.appendChild(tagRow);
    const tutorialStatus = getTutorialStatus(sample);
    if (tutorialStatus !== 'ready') {
      appendText(
        body,
        'span',
        `sample-card__status sample-card__status--${tutorialStatus}`,
        tutorialStatus === 'draft' ? 'Draft tutorial' : 'Tutorial planned'
      );
    }
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

const resetSamplePanels = (sample) => {
  tutorialRequestToken += 1;
  renderFilesPanel(sample);
  clearChildren(tutorialContent);
  if (canLoadTutorial(sample)) {
    appendText(
      tutorialContent,
      'p',
      'panel-status',
      `Open the Tutorial tab to load the step-by-step guide for ${plainTitle(sample)}.`
    );
  } else {
    renderTutorialUnavailable(sample);
  }
};

const setActiveTab = (tabName, { updateUrl = true } = {}) => {
  if (!tabPanels[tabName]) return;
  activeTab = tabName;

  tabButtons.forEach((button) => {
    const isActive = button.dataset.galleryTab === activeTab;
    button.setAttribute('aria-selected', isActive ? 'true' : 'false');
    button.tabIndex = isActive ? 0 : -1;
  });

  Object.entries(tabPanels).forEach(([name, panel]) => {
    if (!panel) return;
    panel.hidden = name !== activeTab;
  });

  const sample = selectedSample();
  if (!sample) return;
  if (updateUrl) {
    window.history.replaceState(null, '', getSampleRouteUrl(sample));
  }
  if (activeTab === 'tutorial') {
    requestTutorialForSample(sample, { renderTutorialPanel: true });
  } else if (activeTab === 'files') {
    requestTutorialForSample(sample, { renderFiles: true });
  }
};

const selectSample = (id, { updateUrl = true } = {}) => {
  const sample = examples.find((entry) => entry.id === id) || examples[0];
  if (!sample) return;

  selectedId = sample.id;
  selectedTitle.innerHTML = sample.title;
  selectedDescription.textContent = sample.description || '';
  selectedDescription.hidden = !sample.description;
  selectedMeta.textContent = [sample.difficulty, sample.workflow, sample.inputSummary, sample.estimatedTime]
    .filter(Boolean)
    .join(' · ');
  selectedMeta.hidden = !selectedMeta.textContent;
  renderTags(selectedTags, sample.tags);
  fileSize.textContent = sample.fileSizeLabel;
  commandBlock.textContent = sample.command;
  const isRunnable = sample.commandKind === 'runnable';
  commandKind.textContent = isRunnable ? 'Runnable' : 'Provenance';
  commandKind.className = `command-kind command-kind--${isRunnable ? 'runnable' : 'provenance'}`;
  commandNote.textContent = sample.commandNote || '';
  commandNote.hidden = !sample.commandNote;
  copyCommandButton.disabled = !isRunnable;
  copyCommandButton.textContent = isRunnable ? 'Copy command' : 'Provenance only';
  copyCommandButton.title = isRunnable
    ? 'Copy this runnable CLI command'
    : 'This command needs prepared inputs that are not fully public.';
  interactiveStep.textContent = sample.interactiveStep || '';
  interactiveStep.hidden = !sample.interactiveStep;
  previewNote.textContent = sampleIsInteractive(sample)
    ? 'This JavaScript-enabled SVG is embedded in a sandboxed iframe.'
    : 'This static SVG is embedded in a sandboxed iframe; it contains no JavaScript controls.';
  updateActionLinks(sample);
  resetSamplePanels(sample);
  updatePressedState();
  setLoading(sample);

  if (updateUrl) {
    window.history.replaceState(null, '', getSampleRouteUrl(sample));
  }

  if (activeTab === 'tutorial') {
    requestTutorialForSample(sample, { renderTutorialPanel: true });
  } else if (activeTab === 'files') {
    requestTutorialForSample(sample, { renderFiles: true });
  }

  window.requestAnimationFrame(() => {
    frame.src = sample.svg;
  });
};

const showError = (message) => {
  selectedTitle.textContent = 'Gallery unavailable';
  selectedDescription.textContent = message;
  selectedDescription.hidden = false;
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

copyCommandButton.addEventListener('click', async () => {
  const sample = examples.find((entry) => entry.id === selectedId);
  if (!sample || sample.commandKind !== 'runnable') return;
  try {
    await copyText(sample.command);
    copyCommandButton.textContent = 'Copied';
    window.setTimeout(() => {
      if (selectedId === sample.id) copyCommandButton.textContent = 'Copy command';
    }, 1600);
  } catch (error) {
    copyCommandButton.textContent = 'Copy failed';
  }
});

tabButtons.forEach((button, index) => {
  button.addEventListener('click', () => {
    setActiveTab(button.dataset.galleryTab);
  });
  button.addEventListener('keydown', (event) => {
    if (!['ArrowLeft', 'ArrowRight', 'Home', 'End'].includes(event.key)) return;
    event.preventDefault();
    let nextIndex = index;
    if (event.key === 'ArrowLeft') nextIndex = (index - 1 + tabButtons.length) % tabButtons.length;
    if (event.key === 'ArrowRight') nextIndex = (index + 1) % tabButtons.length;
    if (event.key === 'Home') nextIndex = 0;
    if (event.key === 'End') nextIndex = tabButtons.length - 1;
    tabButtons[nextIndex].focus();
    setActiveTab(tabButtons[nextIndex].dataset.galleryTab);
  });
});

const syncFromLocation = () => {
  const tabName = getTabFromHash() || 'preview';
  setActiveTab(tabName, { updateUrl: false });

  const legacySampleId = getLegacySampleIdFromHash();
  if (legacySampleId) {
    selectSample(legacySampleId, { updateUrl: false });
    return;
  }

  const pathSampleId = getSampleIdFromPath();
  if (pathSampleId && pathSampleId !== selectedId) {
    selectSample(pathSampleId, { updateUrl: false });
  }
};

window.addEventListener('hashchange', syncFromLocation);
window.addEventListener('popstate', syncFromLocation);

setActiveTab(activeTab, { updateUrl: false });

fetch('./examples.json', { cache: 'no-cache' })
  .then((response) => {
    if (!response.ok) throw new Error(`Could not load examples.json (${response.status})`);
    return response.json();
  })
  .then((data) => {
    examples = Array.isArray(data) ? data : [];
    if (!examples.length) throw new Error('No gallery examples are configured.');
    const pathSampleId = getSampleIdFromPath();
    const legacySampleId = getLegacySampleIdFromHash();
    const initialId = pathSampleId || legacySampleId || examples[0]?.id || DEFAULT_SAMPLE_ID;
    selectedId = initialId;
    renderSampleList();
    selectSample(initialId, { updateUrl: !legacySampleId });
  })
  .catch((error) => {
    showError(error?.message || 'Could not load gallery examples.');
  });
