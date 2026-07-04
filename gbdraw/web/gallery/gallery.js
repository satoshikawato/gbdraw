const DEFAULT_SAMPLE_ID = 'Vnig_TUMSAT-TG-2018';
const COMMON_COLOR_RULE_GUIDE_PATH = './tutorials/_common-color-rule-guide.json';

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

let examples = [];
let selectedId = '';
let activeTab = 'preview';
let copyStatusTimer = 0;
let tutorialRequestToken = 0;
let commonColorRuleGuide = null;
let commonColorRuleGuidePromise = null;
const tutorialCache = new Map();

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

const selectedSample = () => examples.find((entry) => entry.id === selectedId);

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

const renderMedia = (parent, media) => {
  mediaEntryList(media).forEach((entry) => {
    const mediaObject = typeof entry === 'string' ? { src: entry } : entry;
    if (!mediaObject || typeof mediaObject !== 'object') return;

    const source = asText(mediaObject.src) || asText(mediaObject.href);
    if (!source) return;

    const figure = document.createElement('figure');
    figure.className = 'tutorial-media';
    const mediaElement = document.createElement(mediaIsVideo(mediaObject, source) ? 'video' : 'img');
    mediaElement.src = source;
    if (mediaElement.tagName === 'VIDEO') {
      mediaElement.controls = true;
      mediaElement.preload = 'metadata';
    } else {
      mediaElement.loading = 'lazy';
      mediaElement.alt = asText(mediaObject.alt) || asText(mediaObject.caption) || '';
    }
    mediaElement.addEventListener('error', () => {
      clearChildren(figure);
      appendText(figure, 'p', 'panel-status panel-status--error', `Media unavailable: ${source}`);
    });
    figure.appendChild(mediaElement);

    const caption = asText(mediaObject.caption);
    if (caption) appendText(figure, 'figcaption', '', caption);
    parent.appendChild(figure);
  });
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

    const note = asText(stepObject.note);
    if (note) appendText(item, 'p', 'tutorial-step__note', note);
    renderMedia(item, stepObject.media);

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
  appendText(tutorialContent, 'h2', '', 'Tutorial coming soon');
  appendText(
    tutorialContent,
    'p',
    'tutorial-summary',
    `A step-by-step web tutorial is not available for ${plainTitle(sample)} yet. Use the Preview and Session controls to inspect the completed diagram, or open the Command tab for the CLI command.`
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

  const metaItems = [
    asText(tutorial.estimatedTime),
    asText(tutorial.difficulty),
    asText(tutorial.appMode),
    asText(tutorial.primaryWorkflow)
  ].filter(Boolean);
  if (metaItems.length) {
    const meta = document.createElement('div');
    meta.className = 'tutorial-meta';
    metaItems.forEach((item) => appendText(meta, 'span', 'tutorial-meta__item', item));
    tutorialContent.appendChild(meta);
  }

  renderTextListSection(tutorialContent, 'Requirements', tutorial.requirements);
  if (asArray(tutorial.downloads).length) {
    const downloads = appendSection(tutorialContent, 'Input downloads');
    renderDownloads(downloads, tutorial.downloads);
  }
  renderStepSection(tutorialContent, 'Quick reproduce', tutorial.quickReproduce);
  renderStepSection(tutorialContent, 'Manual rebuild', tutorial.manualSteps, { numbered: true });
  renderStepSection(
    tutorialContent,
    asText(tutorial.commonColorRuleGuide?.title) || 'Color rule basics',
    tutorial.commonColorRuleGuide?.colorRules
  );
  renderStepSection(tutorialContent, 'Color rules', tutorial.colorRules);
  renderStepSection(tutorialContent, 'Post-generation edits', tutorial.postGenerationEdits);
  renderStepSection(tutorialContent, 'LOSAT tips', tutorial.losatTips);
  renderStepSection(tutorialContent, 'Troubleshooting', tutorial.troubleshooting);
  renderRelated(tutorialContent, tutorial.related);
};

const renderArtifactList = (parent, sample) => {
  const artifacts = [
    {
      label: 'Interactive SVG',
      href: sample.svg,
      note: sample.sourceOutput || ''
    },
    {
      label: 'Session JSON',
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

  if (!sample?.tutorial) {
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
  frame.title = `Interactive gbdraw SVG demo: ${plainTitle(sample)}`;
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

const resetSamplePanels = (sample) => {
  tutorialRequestToken += 1;
  renderFilesPanel(sample);
  clearChildren(tutorialContent);
  if (sample.tutorial) {
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

const setActiveTab = (tabName) => {
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
  if (activeTab === 'tutorial') {
    requestTutorialForSample(sample, { renderTutorialPanel: true });
  } else if (activeTab === 'files') {
    requestTutorialForSample(sample, { renderFiles: true });
  }
};

const selectSample = (id, { updateHash = true } = {}) => {
  const sample = examples.find((entry) => entry.id === id) || examples[0];
  if (!sample) return;

  selectedId = sample.id;
  selectedTitle.innerHTML = sample.title;
  selectedDescription.textContent = sample.description || '';
  selectedDescription.hidden = !sample.description;
  renderTags(selectedTags, sample.tags);
  fileSize.textContent = sample.fileSizeLabel;
  commandBlock.textContent = sample.command;
  interactiveStep.textContent = sample.interactiveStep || '';
  interactiveStep.hidden = !sample.interactiveStep;
  updateActionLinks(sample);
  resetSamplePanels(sample);
  updatePressedState();
  setLoading(sample);

  if (updateHash && normalizeHash() !== sample.id) {
    window.history.replaceState(null, '', `#${encodeURIComponent(sample.id)}`);
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

window.addEventListener('hashchange', () => {
  selectSample(normalizeHash(), { updateHash: false });
});

setActiveTab(activeTab);

fetch('./examples.json', { cache: 'no-cache' })
  .then((response) => {
    if (!response.ok) throw new Error(`Could not load examples.json (${response.status})`);
    return response.json();
  })
  .then((data) => {
    examples = Array.isArray(data) ? data : [];
    if (!examples.length) throw new Error('No gallery examples are configured.');
    const initialId = normalizeHash() || examples[0]?.id || DEFAULT_SAMPLE_ID;
    selectedId = initialId;
    renderSampleList();
    selectSample(initialId);
  })
  .catch((error) => {
    showError(error?.message || 'Could not load gallery examples.');
  });
