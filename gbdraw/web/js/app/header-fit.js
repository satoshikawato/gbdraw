// @ts-check
// #1383 (PD-OI-108): the one-row header gives way without narrowing the
// Session title (OV-368). Where the words of the Gallery & tutorials link do
// not fit, it shows its cap alone; where even the cap narrows the title, the
// file buttons (Reset Settings, Save Session, Load Session, Load example) show
// their icons alone too, as they do below 1280 px, and the link shows its
// words again if they fit then; with neither words, the toolbar's groups close
// up. Where even that narrows the title, the Session block shows the title
// without its "Session" caption, as the wrapped header does, and the groups
// close up further (TL-D06). The room
// is the tagline room (`.app-tagline-room`), which takes the width
// the rest of the header leaves. The stylesheet applies these states only while
// the header is one row (above 768 px), and the wrapped header shows the cap
// alone; the buttons and the link keep their names in `aria-label`.
import { uiLocale } from '../utils/i18n.js';

/**
 * A state of the header: whether the file buttons and the Gallery link show
 * their words, and whether the Session title has its caption.
 * @typedef {{ fileWords: boolean, galleryWords: boolean, sessionCaption: boolean }} HeaderFit
 */

/**
 * The states, in the order the header prefers them.
 * @type {readonly HeaderFit[]}
 */
export const HEADER_FITS = Object.freeze([
  { fileWords: true, galleryWords: true, sessionCaption: true },
  { fileWords: true, galleryWords: false, sessionCaption: true },
  { fileWords: false, galleryWords: true, sessionCaption: true },
  { fileWords: false, galleryWords: false, sessionCaption: true },
  { fileWords: false, galleryWords: false, sessionCaption: false }
]);

/**
 * The first state that leaves the tagline room a pixel or more, so that nothing
 * else in the header narrows; the last state when none does. The caption hides
 * only where the title narrows (TL-D06): where the icons alone leave the room
 * any width, nothing narrows, and the caption stays.
 * @param {(fit: HeaderFit) => number} roomOf the tagline room's width in a state
 * @returns {HeaderFit}
 */
export const chooseHeaderFit = (roomOf) => {
  const last = HEADER_FITS.length - 1;
  for (let index = 0; index < last; index += 1) {
    const room = roomOf(HEADER_FITS[index]);
    if (room >= 1 || (index === last - 1 && room > 0)) return HEADER_FITS[index];
  }
  return HEADER_FITS[last];
};

/** @param {HTMLElement} header @param {HeaderFit} fit */
const showHeaderFit = (header, fit) => {
  header.toggleAttribute('data-file-icons', !fit.fileWords);
  header.toggleAttribute('data-gallery-icon', !fit.galleryWords);
  header.toggleAttribute('data-session-title-alone', !fit.sessionCaption);
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
    if (!(room instanceof HTMLElement) || typeof ResizeObserver !== 'function') return;
    let frame = 0;
    const fit = () => {
      frame = 0;
      showHeaderFit(header, chooseHeaderFit((state) => {
        showHeaderFit(header, state);
        return room.getBoundingClientRect().width;
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
