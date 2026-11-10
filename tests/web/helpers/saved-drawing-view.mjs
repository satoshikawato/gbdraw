// A saved Session 49 read as the one-drawing views a test asserts on. The
// layout is the production reader's (`drawingViewOf`, `sessionDocumentOfViews`);
// a test finds a drawing by its ID.
import { drawingViewOf, sessionDocumentOfViews } from '../../../gbdraw/web/js/services/session-drawings-upgrade.js';
import { MODE_SCOPED_SETTINGS } from '../../../gbdraw/web/js/mode-scoped-settings.generated.js';
import { FEATURE_CATALOG_SCHEMA_BY_SESSION_VERSION } from '../../../gbdraw/web/js/contracts/schema.generated.js';

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

/**
 * Every saved drawing read as one 46-48-shaped view: the shown drawing's view
 * with `modes.<mode>` of each other drawing added (one drawing per mode in the
 * Sessions a Save writes from the two-mode state). For reads only; load the
 * saved document itself.
 */
export const savedModesView = (document) => {
  const view = savedDrawingView(document, document.ui.activeDrawingId);
  for (const entry of document.drawings) {
    const other = drawingViewOf(document, entry);
    view.modes = { ...other.modes, ...view.modes };
  }
  return view;
};

/**
 * `view` (a one-drawing view, or a 46-48-shaped view of every drawing) as a
 * Session `version` file holds it: `version`, no slice field a later Session
 * added (a registry row's `since`), which its writer did not write, and each
 * feature catalog (the view's own, the other mode's Result set's, each mode
 * slice's) labelled with the schema that version writes: a catalog with no
 * sequence (schema 6) is a valid schema 5 catalog
 * (tests/utils/session_drawing_view.py `view_at` does the same).
 */
export const drawingViewAt = (view, version) => {
  const at = structuredClone(view);
  const catalogSchema = FEATURE_CATALOG_SCHEMA_BY_SESSION_VERSION[version];
  const labelCatalog = (holder) => {
    const catalog = holder?.editorState?.featureCatalog;
    if (catalogSchema !== undefined && catalog && typeof catalog === 'object') catalog.schema = catalogSchema;
  };
  labelCatalog(at);
  labelCatalog(at.otherModeResult);
  for (const slice of Object.values(at.modes ?? {})) {
    labelCatalog(slice);
    for (const { domain, path, since } of MODE_SCOPED_SETTINGS.rows) {
      const container = domain.split('.').reduce((value, part) => value?.[part], slice);
      if (since > version && container && typeof container === 'object') delete container[path];
    }
  }
  return { ...at, version };
};
