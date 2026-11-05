// @ts-check
// The Rust drawing engine in the diagram Worker: Python's optional `gbdraw_core`
// (gbdraw/render/native.py). A module that is missing or fails to load leaves
// every group to Python; tools/build_gbdraw_core_wasm.py generates it.

const MODULE_URL = new URL('../../wasm/gbdraw-core/gbdraw_core.wasm', import.meta.url);

/**
 * The `gbdraw_core` module over one instantiated engine (rust/wasm/src/lib.rs).
 * @param {Uint8Array | ArrayBuffer} wasmBytes
 */
export const createGbdrawCore = async (wasmBytes) => {
  const { instance } = await WebAssembly.instantiate(wasmBytes);
  const engine = /** @type {Record<string, any>} */ (instance.exports);
  const encoder = new TextEncoder();
  const copyIn = (/** @type {Uint8Array} */ bytes) => {
    const pointer = engine.gbdraw_alloc(bytes.length);
    new Uint8Array(engine.memory.buffer, pointer, bytes.length).set(bytes);
    return pointer;
  };
  return {
    /**
     * @param {string} params
     * @param {{ toJs(): Uint8Array, destroy(): void }} columns Python bytes
     */
    pairwise_ribbons(params, columns) {
      const paramBytes = encoder.encode(params);
      const columnBytes = columns.toJs();
      columns.destroy();
      const paramPointer = copyIn(paramBytes);
      const columnPointer = copyIn(columnBytes);
      const status = engine.gbdraw_pairwise_ribbons(paramPointer, paramBytes.length, columnPointer, columnBytes.length);
      engine.gbdraw_free(paramPointer, paramBytes.length);
      engine.gbdraw_free(columnPointer, columnBytes.length);
      // A view of the engine's memory, which Python copies before its next call.
      const result = new Uint8Array(engine.memory.buffer, engine.gbdraw_result_pointer(), engine.gbdraw_result_length());
      if (status !== 0) throw new Error(new TextDecoder().decode(result));
      return result;
    }
  };
};

/**
 * Register the engine for Python's `import gbdraw_core`.
 * @param {{ registerJsModule(name: string, module: object): void }} pyodide
 */
export const registerGbdrawCore = async (pyodide) => {
  try {
    const response = await fetch(MODULE_URL);
    if (!response.ok) return;
    pyodide.registerJsModule('gbdraw_core', await createGbdrawCore(await response.arrayBuffer()));
  } catch (error) {
    console.warn('The Rust drawing engine did not load; Python draws every group.', error);
  }
};
