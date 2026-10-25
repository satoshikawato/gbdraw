// A saved Session 49 read as the one-drawing views a test asserts on. The
// layout is the production reader's (`drawingViewOf`, `sessionDocumentOfViews`);
// a test finds a drawing by its ID.
import { drawingViewOf, sessionDocumentOfViews } from '../../../gbdraw/web/js/services/session-drawings-upgrade.js';

const drawingById = (document, id) => {
  const drawing = document.drawings.find((entry) => entry.id === id);
  if (!drawing) throw new Error(`the saved Session holds no drawing "${id}"`);
  return drawing;
};

/** The one-drawing view of the saved drawing `id` (`modes.<mode>`, `ui.mode`, its Result set at the top level). */
export const savedDrawingView = (document, id) => drawingViewOf(document, drawingById(document, id));

/**
 * A copy of the saved Session 49 whose drawing `id` is changed through its
 * one-drawing view: `mutate(view)` edits the view, and the document is built
 * again from every drawing's view in tab order. The project's own fields
 * (`ui`, `editorState`) come from the edited view.
 */
export const withSavedDrawingView = (document, id, mutate) => {
  const views = document.drawings.map((entry) => ({
    id: entry.id, name: entry.name, view: structuredClone(drawingViewOf(document, entry))
  }));
  const edited = views.find((entry) => entry.id === id);
  if (!edited) throw new Error(`the saved Session holds no drawing "${id}"`);
  mutate(edited.view);
  const rebuilt = sessionDocumentOfViews([edited, ...views.filter((entry) => entry !== edited)], document.ui.activeDrawingId);
  rebuilt.drawings = views.map((entry) => rebuilt.drawings.find((drawing) => drawing.id === entry.id));
  return rebuilt;
};
