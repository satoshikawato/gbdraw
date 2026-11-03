'use strict';

// A saved Session 49 as specs read it (DESIGN §1): each drawing is one
// `drawings[]` entry, found by its drawing ID, never by its mode or name.

/** The saved drawing `id`; throws when the Session holds none. */
const savedDrawing = (session, id) => {
  const drawing = session.drawings.find((entry) => entry.id === id);
  if (!drawing) throw new Error(`the saved Session holds no drawing "${id}"`);
  return drawing;
};

/** The drawing a saved Session opens: the one `ui.activeDrawingId` names. */
const shownSavedDrawing = (session) => savedDrawing(session, session.ui.activeDrawingId);

module.exports = { savedDrawing, shownSavedDrawing };
