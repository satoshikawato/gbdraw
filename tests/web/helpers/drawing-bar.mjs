// ESM twin of drawing-bar.cjs: a re-export, no second implementation.
import { createRequire } from 'node:module';

const { showDrawing, drawingTab, drawingIds, shownDrawing } = createRequire(import.meta.url)('./drawing-bar.cjs');

export { showDrawing, drawingTab, drawingIds, shownDrawing };
