// A notice keeps how to word it (`Wording` of gbdraw/web/js/utils/i18n.js) and is worded where it shows (OV-2796,
// OV-3674). These word one as it shows in the shown UI language.

/** A notice request of the banner (`NoticeRequest` of app/app-messages.js) with its text and lines worded. */
export const wordedNotice = ({ text, items, ...rest }) => ({ text: text(), ...(items ? { items: items.map((item) => item()) } : {}), ...rest });

/** The text of a notice that may be cleared (`Wording | null`): empty when cleared. */
export const shownWording = (wording) => (wording ? wording() : '');
