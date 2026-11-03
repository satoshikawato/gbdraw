// The palette definitions the app loads at boot (app/palettes.js reads
// gbdraw/web/gallery/palettes/palettes.json). A new drawing takes its default
// palette's colors from them (`withNewDrawingDefaults`, services/reset.js), as
// tools/generate_mode_slice_defaults.mjs reads a mode's defaults, so a test
// whose Save meets an unused drawing runs with them.
import { readFileSync } from 'node:fs';
import { normalizePaletteDefinitions } from '../../../gbdraw/web/js/utils/color-utils.js';

/**
 * Runs `run` with the boot palette definitions in `state`, then puts back the ones it had.
 * @param {{ paletteDefinitions: { value: unknown } }} state
 * @param {() => Promise<void>} run
 */
export const withBootPalettes = async (state, run) => {
  const palettes = JSON.parse(readFileSync(new URL('../../../gbdraw/web/gallery/palettes/palettes.json', import.meta.url), 'utf8'));
  const previous = state.paletteDefinitions.value;
  state.paletteDefinitions.value = normalizePaletteDefinitions(palettes.palettes);
  try {
    await run();
  } finally {
    state.paletteDefinitions.value = previous;
  }
};
