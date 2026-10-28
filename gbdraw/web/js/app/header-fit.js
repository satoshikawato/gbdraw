// @ts-check
// #1383 (PD-OI-108): the one-row header gives way without narrowing the
// Session title (OV-368). Where the words of the Gallery & tutorials link do
// not fit, it shows its cap alone; where even the cap narrows the title, the
// file buttons (Reset Settings, Save Session, Load Session, Load example) show
// their icons alone too, as they do below 1280 px, and the link shows its
// words again if they fit then; with neither words, the toolbar's groups close
// up. Where even that narrows the title, the Session block shows the title
// without its "Session" caption, as the wrapped header does, and the groups
// close up further (TL-D06). OV-2177, OV-2178: where the title still narrows,
// the closed Language menu takes the width of the shown language's name
// instead of the longest one; where the title would still be narrower than
// `SESSION_TITLE_MIN_WIDTH` (the header overflows the window only past that),
// the header wraps to two rows, the logo above the toolbar. The room
// is the tagline room (`.app-tagline-room`), which takes the width
// the rest of the header leaves. The stylesheet applies these states only while
// the header is one row (above 768 px), and the wrapped header shows the cap
// alone; the buttons and the link keep their names in `aria-label`.
import { uiLocale } from '../utils/i18n.js';

/**
 * A state of the header: whether the file buttons and the Gallery link show
 * their words, whether the Session title has its caption, whether the closed
 * Language menu is as wide as its longest language name, and whether the
 * header keeps one row.
 * @typedef {{ fileWords: boolean, galleryWords: boolean, sessionCaption: boolean, languageWide: boolean, oneRow: boolean }} HeaderFit
 */

/**
 * What a state leaves: the tagline room's width, and the Session title
 * button's width (it narrows once the room is gone).
 * @typedef {{ room: number, title: number }} HeaderRoom
 */

/**
 * The states, in the order the header prefers them.
 * @type {readonly HeaderFit[]}
 */
export const HEADER_FITS = Object.freeze([
  { fileWords: true, galleryWords: true, sessionCaption: true, languageWide: true, oneRow: true },
  { fileWords: true, galleryWords: false, sessionCaption: true, languageWide: true, oneRow: true },
  { fileWords: false, galleryWords: true, sessionCaption: true, languageWide: true, oneRow: true },
  { fileWords: false, galleryWords: false, sessionCaption: true, languageWide: true, oneRow: true },
  { fileWords: false, galleryWords: false, sessionCaption: false, languageWide: true, oneRow: true },
  { fileWords: false, galleryWords: false, sessionCaption: false, languageWide: false, oneRow: true },
  { fileWords: false, galleryWords: false, sessionCaption: false, languageWide: false, oneRow: false }
].map((fit) => Object.freeze(fit)));

/**
 * The narrowest Session title (its pencil and about nine characters, 5rem)
 * the one-row header keeps; narrower, the header wraps (OV-2178).
 */
export const SESSION_TITLE_MIN_WIDTH = 80;

/**
 * The first one-row state that leaves the tagline room a pixel or more, so
 * that nothing else in the header narrows. The caption hides only where the
 * title narrows (TL-D06): where the icons alone leave the room any width,
 * nothing narrows, and the caption stays. Where every one-row state narrows
 * the title, the last of them while it keeps the title
 * `SESSION_TITLE_MIN_WIDTH` wide; else the two-row header.
 * @param {(fit: HeaderFit) => HeaderRoom} measure what the header leaves in a state
 * @returns {HeaderFit}
 */
export const chooseHeaderFit = (measure) => {
  const wrapped = HEADER_FITS.length - 1;
  const tightest = wrapped - 1;
  const captionKept = HEADER_FITS.findIndex((fit) => !fit.sessionCaption) - 1;
  for (let index = 0; index <= tightest; index += 1) {
    const { room, title } = measure(HEADER_FITS[index]);
    if (room >= 1 || (index === captionKept && room > 0)) return HEADER_FITS[index];
    if (index === tightest && title >= SESSION_TITLE_MIN_WIDTH) return HEADER_FITS[index];
  }
  return HEADER_FITS[wrapped];
};

/** @param {HTMLElement} header @param {HeaderFit} fit */
const showHeaderFit = (header, fit) => {
  header.toggleAttribute('data-file-icons', !fit.fileWords);
  header.toggleAttribute('data-gallery-icon', !fit.galleryWords);
  header.toggleAttribute('data-session-title-alone', !fit.sessionCaption);
  header.toggleAttribute('data-language-narrow', !fit.languageWide);
  header.toggleAttribute('data-header-wrapped', !fit.oneRow);
};

/** @type {WeakMap<HTMLElement, () => void>} */
const fitStops = new WeakMap();

// `v-header-fit` on the header: once per frame after the room changes size or
// the UI language changes, it lays out the states in order before the frame
// is painted and keeps the first that fits.
export const headerFit = {
  /** @param {HTMLElement} header */
  mounted: (header) => {
    const room = header.querySelector('.app-tagline-room');
    const title = header.querySelector('.app-session-button');
    if (!(room instanceof HTMLElement) || !(title instanceof HTMLElement) || typeof ResizeObserver !== 'function') return;
    let frame = 0;
    const fit = () => {
      frame = 0;
      showHeaderFit(header, chooseHeaderFit((state) => {
        showHeaderFit(header, state);
        return { room: room.getBoundingClientRect().width, title: title.getBoundingClientRect().width };
      }));
    };
    const schedule = () => { frame ||= requestAnimationFrame(fit); };
    const observer = new ResizeObserver(schedule);
    observer.observe(room);
    const stopLanguage = window.Vue.watch(() => uiLocale.value, schedule, { flush: 'post' });
    fitStops.set(header, () => {
      observer.disconnect();
      stopLanguage();
      if (frame) cancelAnimationFrame(frame);
    });
  },
  /** @param {HTMLElement} header */
  unmounted: (header) => {
    fitStops.get(header)?.();
    fitStops.delete(header);
  }
};
