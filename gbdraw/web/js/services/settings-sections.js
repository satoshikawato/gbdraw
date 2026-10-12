// @ts-check
// The sections of the settings pane, in card order. This list is the one owner
// of a section's key, title, icon, and modes. A card (or the Linear Comparison
// section) carries its key as `data-settings-section` and shows its title and
// icon through `SettingsSectionTitle`; `SettingsSectionIndex` lists the
// sections of the shown mode (js/components.js). Basic and About & Citation
// are cards without an entry.

/**
 * @template {string} Key
 * @typedef {object} SettingsSectionOf
 * @property {Key} key Identity of the section; the value of `data-settings-section`.
 * @property {string} title The heading text and the name of the index entry.
 * @property {string} icon Phosphor icon class of the heading and of the index entry.
 * @property {readonly ('circular' | 'linear')[]} modes The modes that show the section.
 */

/**
 * @template {string} Key
 * @param {readonly SettingsSectionOf<Key>[]} sections
 */
const catalog = (sections) => Object.freeze(sections.map((section) => Object.freeze(section)));

/** @type {readonly ('circular' | 'linear')[]} */
const BOTH_MODES = ['circular', 'linear'];

export const SETTINGS_SECTIONS = catalog([
  { key: 'input-genomes', title: 'Input Genomes', icon: 'ph-folder-open', modes: BOTH_MODES },
  { key: 'comparison', title: 'Comparison', icon: 'ph-arrows-left-right', modes: ['linear'] },
  { key: 'advanced-comparison', title: 'Advanced comparison and layout', icon: 'ph-gear-six', modes: ['linear'] },
  { key: 'region-annotations', title: 'Region Annotations', icon: 'ph-brackets-curly', modes: BOTH_MODES },
  { key: 'pairwise-comparisons', title: 'Pairwise Comparisons', icon: 'ph-arrows-left-right', modes: ['circular'] },
  { key: 'layout', title: 'Layout', icon: 'ph-layout', modes: BOTH_MODES },
  { key: 'titles', title: 'Titles & Record Labels', icon: 'ph-textbox', modes: BOTH_MODES },
  { key: 'legend', title: 'Legend', icon: 'ph-list-bullets', modes: BOTH_MODES },
  { key: 'colors', title: 'Colors', icon: 'ph-palette', modes: BOTH_MODES },
  { key: 'features', title: 'Features', icon: 'ph-shapes', modes: BOTH_MODES },
  { key: 'labels', title: 'Labels', icon: 'ph-text-t', modes: BOTH_MODES },
  // Axis & Scale is one card per mode, each with its own icon.
  { key: 'linear-axis', title: 'Axis & Scale', icon: 'ph-ruler', modes: ['linear'] },
  { key: 'circular-axis', title: 'Axis & Scale', icon: 'ph-paint-brush', modes: ['circular'] },
  { key: 'dinucleotide', title: 'Dinucleotide content/skew', icon: 'ph-dna', modes: BOTH_MODES }
]);

/** @typedef {(typeof SETTINGS_SECTIONS)[number]} SettingsSection */
/** @typedef {SettingsSection['key']} SettingsSectionKey */

/** @param {'circular' | 'linear'} mode */
export const settingsSectionsOf = (mode) => SETTINGS_SECTIONS.filter((section) => section.modes.includes(mode));

// What scrolls the settings: their scroll area, or the page in the phone
// layout, where the area does not scroll.
/** @param {Element} area The settings scroll area. */
export const settingsScrollOwner = (area) =>
  (getComputedStyle(area).overflowY !== 'visible' ? area : document.scrollingElement);
