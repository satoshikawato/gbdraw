// @ts-check
// The sections of the settings pane, in card order. This list is the one owner
// of a section's key and modes, and of what its index entry shows: the title
// and icon (the section heading shows them too, through
// `SettingsSectionTitle`) and the short name under the icon in the rail. An
// element carries its key through `v-settings-section`; `SettingsSectionIndex`
// lists the entries of the shown mode (js/components.js). A card without an
// entry of its own (Basic, the preserved-settings notice) names the entry that
// docks it in the rail. The list of drawn elements below names, for each
// element of the Result that has a section, the section a click on it opens.
import {
  DINUCLEOTIDE_CONTENT_GROUPS, DINUCLEOTIDE_SKEW_GROUPS, LENGTH_BAR_ID, PLOT_TITLE_ID, TICK_GROUPS, trackGroupSelector
} from './drawn-groups.js';

/** @import { SettingsPaneDisplay } from './settings-pane-state.js' */

/** @typedef {readonly ('circular' | 'linear')[]} SettingsSectionModes The modes that show a section. */

/**
 * A section with an index entry.
 * @template {string} Key
 * @typedef {object} SettingsEntryOf
 * @property {Key} key Identity of the section; the value of `data-settings-section`.
 * @property {string} title The heading text; the entry is named "Go to <title>".
 * @property {string} icon Phosphor icon class of the heading and of the entry.
 * @property {string} short The name under the icon of the rail entry.
 * @property {SettingsSectionModes} modes
 * @property {true} [input] An entry that chooses what is drawn; the rail sets these apart from the rest.
 * @property {true} [railOnly] An entry of the rail alone; the stacked row leaves it out.
 */

/**
 * A card that the rail docks with another section's entry.
 * @template {string} Key
 * @typedef {object} SettingsAttachedOf
 * @property {Key} key
 * @property {Key} dock The key of that entry.
 * @property {SettingsSectionModes} modes
 */

/**
 * @template {string} Key
 * @param {readonly (SettingsEntryOf<Key> | SettingsAttachedOf<Key>)[]} sections
 */
const catalog = (sections) => Object.freeze(sections.map((section) => Object.freeze(section)));

/** @type {SettingsSectionModes} */
const BOTH_MODES = ['circular', 'linear'];

export const SETTINGS_SECTIONS = catalog([
  { key: 'input-genomes', title: 'Input Genomes', icon: 'ph-folder-open', short: 'Input', modes: BOTH_MODES, input: true },
  { key: 'comparison', title: 'Comparison', icon: 'ph-arrows-left-right', short: 'Compare', modes: ['linear'], input: true },
  { key: 'basic', dock: 'input-genomes', modes: BOTH_MODES },
  { key: 'advanced-comparison', title: 'Advanced comparison and layout', icon: 'ph-gear-six', short: 'Advanced', modes: ['linear'] },
  { key: 'region-annotations', title: 'Region Annotations', icon: 'ph-brackets-curly', short: 'Regions', modes: BOTH_MODES },
  { key: 'pairwise-comparisons', title: 'Pairwise Comparisons', icon: 'ph-arrows-left-right', short: 'Compare', modes: ['circular'] },
  { key: 'preserved-settings', dock: 'input-genomes', modes: BOTH_MODES },
  { key: 'layout', title: 'Layout', icon: 'ph-layout', short: 'Layout', modes: BOTH_MODES },
  { key: 'font', title: 'Font', icon: 'ph-text-aa', short: 'Font', modes: BOTH_MODES },
  { key: 'titles', title: 'Titles & Record Labels', icon: 'ph-textbox', short: 'Titles', modes: BOTH_MODES },
  { key: 'legend', title: 'Legend', icon: 'ph-list-bullets', short: 'Legend', modes: BOTH_MODES },
  { key: 'colors', title: 'Colors', icon: 'ph-palette', short: 'Colors', modes: BOTH_MODES },
  { key: 'features', title: 'Features', icon: 'ph-shapes', short: 'Features', modes: BOTH_MODES },
  { key: 'labels', title: 'Labels', icon: 'ph-text-t', short: 'Labels', modes: BOTH_MODES },
  // Axis & Scale is one card per mode, each with its own icon.
  { key: 'linear-axis', title: 'Axis & Scale', icon: 'ph-ruler', short: 'Axis', modes: ['linear'] },
  { key: 'circular-axis', title: 'Axis & Scale', icon: 'ph-paint-brush', short: 'Axis', modes: ['circular'] },
  { key: 'dinucleotide', title: 'Dinucleotide content/skew', icon: 'ph-dna', short: 'GC / skew', modes: BOTH_MODES },
  { key: 'about', title: 'About & Citation', icon: 'ph-info', short: 'About', modes: BOTH_MODES, railOnly: true }
]);

/** @typedef {(typeof SETTINGS_SECTIONS)[number]} SettingsSection */
/** @typedef {SettingsSection['key']} SettingsSectionKey */
/** @typedef {Extract<SettingsSection, { short: string }>} SettingsEntry */

/** @param {SettingsSection} section @returns {section is SettingsEntry} */
const hasEntry = (section) => 'short' in section;

// The index entries of a mode: the stacked row, or the rail, which adds its own.
/** @param {'circular' | 'linear'} mode @param {SettingsPaneDisplay} display */
export const settingsEntriesOf = (mode, display) => SETTINGS_SECTIONS.filter(hasEntry)
  .filter((entry) => entry.modes.includes(mode) && (display === 'rail' || !entry.railOnly));

// The entry that docks a section in the rail: its own, or the one it names.
/** @param {SettingsSectionKey} key */
export const settingsEntryOf = (key) => {
  const section = SETTINGS_SECTIONS.find((candidate) => candidate.key === key);
  return section && 'dock' in section ? section.dock : key;
};

// The rail entry of a mode that shows a section: the section's entry, or the
// mode's entry with the same title (Axis & Scale is one card per mode). A mode
// without either has none.
/** @param {'circular' | 'linear'} mode @param {SettingsSectionKey} key @returns {SettingsSectionKey | undefined} */
export const settingsEntryIn = (mode, key) => {
  const own = SETTINGS_SECTIONS.filter(hasEntry).find((entry) => entry.key === settingsEntryOf(key));
  const entries = settingsEntriesOf(mode, 'rail');
  return (entries.find((entry) => entry === own) ?? entries.find((entry) => entry.title === own?.title))?.key;
};

// The element of a section in the settings scroll area.
/** @param {ParentNode} area The settings scroll area. @param {SettingsSectionKey} key */
export const settingsSectionElement = (area, key) =>
  /** @type {HTMLElement | null} */ (area.querySelector(`[data-settings-section="${key}"]`));

// The drawn elements of a Result whose settings have a section: the scale
// (the Circular ticks, the Linear bar or ruler), every dinucleotide track, and
// the plot title. The groups are those of services/drawn-groups.js; a section
// of one mode stands for the section of the same title in the other
// (`settingsEntryIn`).
/** @type {readonly { readonly selector: string, readonly section: SettingsSectionKey }[]} */
const DRAWN_ELEMENT_SECTIONS = Object.freeze([
  { selector: trackGroupSelector(TICK_GROUPS), section: 'circular-axis' },
  { selector: `#${LENGTH_BAR_ID}`, section: 'linear-axis' },
  { selector: `${trackGroupSelector(DINUCLEOTIDE_CONTENT_GROUPS)}, ${trackGroupSelector(DINUCLEOTIDE_SKEW_GROUPS)}`, section: 'dinucleotide' },
  { selector: `#${PLOT_TITLE_ID}`, section: 'titles' }
]);

export const DRAWN_ELEMENT_SELECTOR = DRAWN_ELEMENT_SECTIONS.map(({ selector }) => selector).join(', ');

/**
 * The section that holds the settings of the drawn element an event in the
 * Result reached: the element itself, or the group it is drawn in.
 * @param {EventTarget | null} target @param {Element} svg
 * @returns {SettingsSectionKey | null}
 */
export const drawnElementSection = (target, svg) => {
  const element = target instanceof Element ? target.closest(DRAWN_ELEMENT_SELECTOR) : null;
  if (!element || !svg.contains(element)) return null;
  return DRAWN_ELEMENT_SECTIONS.find(({ selector }) => element.matches(selector))?.section ?? null;
};

// What scrolls the settings: their scroll area, or the page in the phone
// layout, where the area does not scroll.
/** @param {Element} area The settings scroll area. */
export const settingsScrollOwner = (area) =>
  (getComputedStyle(area).overflowY !== 'visible' ? area : document.scrollingElement);

// What scrolls the settings, and the top and bottom edges of what it shows:
// the scroll area and its box, or the page between its sticky header and the
// Generate bar, which is fixed over the page's bottom edge in that layout.
/** @param {Element} area The settings scroll area. */
export const settingsView = (area) => {
  const owner = settingsScrollOwner(area);
  if (owner === area) {
    const { top, bottom } = area.getBoundingClientRect();
    return { owner, top, bottom };
  }
  return {
    owner,
    top: document.querySelector('.app-header')?.getBoundingClientRect().bottom ?? 0,
    bottom: document.querySelector('.generate-bar')?.getBoundingClientRect().top ?? window.innerHeight
  };
};
