// ESM twin of drawing-bar.cjs: a re-export, no second implementation.
import { createRequire } from 'node:module';

const { showDrawing, drawingIds, shownDrawing } = createRequire(import.meta.url)('./drawing-bar.cjs');

export { showDrawing, drawingIds, shownDrawing };
