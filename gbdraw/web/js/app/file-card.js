// @ts-check
// The File card shell that the Linear and the Circular input lists share: the
// card, its fold to one header line ("File N", the file chip, the moves and,
// in Circular, Remove), and what the header line says of a card's files. Each
// mode's list owner (app/linear-file-cards.js, app/circular-file-cards.js)
// feeds it a card's uid, place, line and state, and shows its own body in the
// default slot; the shell writes nothing and emits the user's presses.
import { linearCardFileFields } from '../services/linear-sources.js';
import { t } from '../utils/i18n.js';
/** @import { InputSlot } from '../services/input-types.js' */

const FILE_SLOT_NAMES = Object.freeze({ gb: 'GenBank', gff: 'GFF3', fasta: 'FASTA' });

/** @param {{ name?: string } | null} file */
export const chosenFileName = (file) => String(file?.name || t('Unnamed file'));

/**
 * What a File card's header line says of the slots its input type shows: each
 * file's name, "no FASTA" for a slot still empty, "No file" for none.
 * @param {string} inputType
 * @param {Readonly<Partial<Record<InputSlot, { name?: string } | null>>>} slots
 */
export const fileCardLine = (inputType, slots) => {
  const fields = linearCardFileFields(inputType);
  if (!fields.some((field) => slots[field])) return t('No file');
  return fields.map((field) => {
    const file = slots[field];
    return file ? chosenFileName(file) : t('no {file}', { file: FILE_SLOT_NAMES[field] });
  }).join(' + ');
};

/**
 * The props of one card.
 * @typedef {object} FileCardProps
 * @property {'linear' | 'circular'} mode Names the card's data attributes (`data-<mode>-source-card`).
 * @property {string} uid The card's identity.
 * @property {number} index Its 0-based place in the list.
 * @property {number} count The cards in the list; a card of a list of one has no moves.
 * @property {string} label The card's accessible name.
 * @property {string} line The header line (`fileCardLine`).
 * @property {boolean} ready Every slot the card shows holds a file: the line shows as the file chip.
 * @property {boolean} failed The card's file could not be read.
 * @property {boolean} folded
 * @property {boolean} disabled A session operation holds the inputs.
 * @property {boolean} canMoveUp
 * @property {boolean} canMoveDown
 * @property {string} moveBlockedReason Why the moves are off, '' when they are not.
 * @property {boolean} removable The header has Remove.
 * @property {string} describedby More ids that describe the header (the mode's own notes).
 */

export const FileCard = {
  template: '#file-card-template',
  props: {
    mode: { type: String, required: true },
    uid: { type: String, required: true },
    index: { type: Number, required: true },
    count: { type: Number, required: true },
    label: { type: String, required: true },
    line: { type: String, required: true },
    ready: Boolean,
    failed: Boolean,
    folded: Boolean,
    disabled: Boolean,
    canMoveUp: Boolean,
    canMoveDown: Boolean,
    moveBlockedReason: { type: String, default: '' },
    removable: Boolean,
    describedby: { type: String, default: '' }
  },
  // fold: the disclosure's state after the user's press or a path that opened
  // it (true: folded). move: -1 up, 1 down. remove: the pressed button, the
  // removal dialog's return focus.
  emits: ['fold', 'move', 'remove'],
  /** @param {FileCardProps} props */
  setup(props) {
    const failedId = () => `${props.mode}-file-failed-${props.uid}`;
    return {
      failedId,
      /** @param {string} part */
      dataOf: (part) => ({ [`data-${props.mode}-source-${part}`]: props.uid }),
      describedIds: () => [props.failed && failedId(), props.describedby].filter(Boolean).join(' ') || null
    };
  }
};
