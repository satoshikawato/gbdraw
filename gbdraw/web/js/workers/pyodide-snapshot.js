// @ts-check
// The diagram Worker's Pyodide memory snapshot: the Wasm memory of a runtime
// that already holds the packages, the gbdraw wheel and the imported Worker
// helpers, so a Worker restores it instead of installing them (about 0.5 s
// instead of 5 s). Only a test server serves one (docs/internal/SELECTIVE_CI.md,
// "Pyodide snapshot"); config.js PYODIDE_SNAPSHOT_URL is null in the
// repository, and the Worker boots cold whenever a restore fails.
//
// Pyodide 0.29's makeMemorySnapshot() keeps only the Wasm memory. Two things
// live outside it and this module carries them in state.json and fs.bin:
//   1. the side modules (.so) loadPackage linked: their code, their function
//      table slots, the GOT globals and the dynamic-linker tables (LDSO);
//      without them the first numpy call traps ("table index is out of bounds");
//   2. MEMFS (site-packages files), which a later import or a package-data
//      read needs.
// The restore re-instantiates each side module at its recorded memory and
// table base, puts the snapshot bytes of its data region back, and rebuilds
// the table slots and GOT values. It uses Emscripten internals that Pyodide
// exports on its module (LDSO, GOT, GOTHandler, getDylinkMetadata, wasmTable,
// newDSO, mergeLibSymbols, setWasmTableEntry, updateTableMap) and the dlopen
// handle layout (HANDLE_OFFSETS): a Pyodide upgrade must rebuild and verify it.
//
// This module runs no Python and loads no runtime; the Worker and
// tools/build_pyodide_snapshot.mjs do that and pass the runtime in.

import { diagnosticError } from '../utils/error-normalization.js';

/** The Worker helper source, which runs before the snapshot is made. */
export const PYTHON_HELPER_FILE = new URL('../app/python-helpers.js', import.meta.url);

/** The files of one snapshot, by role. */
export const PYODIDE_SNAPSHOT_FILES = Object.freeze({
  memory: 'snapshot.bin',
  state: 'state.json',
  files: 'fs.bin'
});

// Emscripten's dlopen handle (`struct dso`): byte offsets of its fields.
const HANDLE_OFFSETS = Object.freeze({ allocated: 8, memAddr: 12, memSize: 16, tableAddr: 20, tableSize: 24 });

/** @typedef {(...args: unknown[]) => unknown} WasmFunction */
/**
 * Emscripten's record of one loaded library.
 * @typedef {object} DynamicLibrary
 * @property {Record<string, unknown>} exports
 * @property {boolean} [global]
 * @property {number} [refcount]
 */
/**
 * @typedef {object} EmscriptenFs
 * @property {(path: string) => string[]} readdir
 * @property {(path: string) => { mode: number }} stat
 * @property {(mode: number) => boolean} isDir
 * @property {(path: string) => Uint8Array} readFile
 * @property {(path: string) => void} mkdirTree
 * @property {(path: string, data: Uint8Array, options: { canOwn: boolean }) => void} writeFile
 */
/**
 * The Emscripten internals Pyodide 0.29 exports on `pyodide._module`.
 * @typedef {object} EmscriptenModule
 * @property {Int8Array} HEAP8
 * @property {Uint32Array} HEAPU32
 * @property {EmscriptenFs} FS
 * @property {WebAssembly.Table} wasmTable
 * @property {{ loadedLibsByName: Record<string, DynamicLibrary>, loadedLibsByHandle: Record<string, DynamicLibrary> }} LDSO
 * @property {Record<string, WebAssembly.Global>} GOT
 * @property {ProxyHandler<object>} GOTHandler
 * @property {(binary: Uint8Array) => { tableSize: number, memorySize: number }} getDylinkMetadata
 * @property {(symbol: string) => { sym: WasmFunction | undefined }} resolveGlobalSymbol
 * @property {(name: string, handle: undefined, state: string) => DynamicLibrary} newDSO
 * @property {(exports: Record<string, unknown>, name: string) => void} mergeLibSymbols
 * @property {(index: number, fn: unknown) => void} setWasmTableEntry
 * @property {(fn: WasmFunction) => number} getFunctionAddress
 * @property {(offset: number, count: number) => void} updateTableMap
 */
/**
 * The part of a Pyodide runtime the snapshot reads and writes.
 * @typedef {object} PyodideRuntime
 * @property {EmscriptenModule} _module
 * @property {{ sitePackages?: string, config: Record<string, unknown> } & Record<string, unknown>} _api
 * @property {Record<string, string>} loadedPackages
 * @property {(options: { serializer: (value: unknown) => HolderDescriptor }) => Uint8Array} makeMemorySnapshot
 */
/**
 * Where a JS object that Python still references lives: a named property of
 * globalThis, of the public Pyodide API or of its internal API, or the holder
 * itself. A null holder restores as null.
 * @typedef {{ holder: 'public' | 'api' | 'global' | null, name?: string }} HolderDescriptor
 */
/**
 * @typedef {object} SideModule
 * @property {string} name
 * @property {number[]} handles
 * @property {boolean} global
 * @property {number} memoryBase
 * @property {number} memorySize
 * @property {number} tableBase
 * @property {number} tableSize
 */
/**
 * state.json: the validity key and the dynamic-linker state.
 * @typedef {object} PyodideSnapshotState
 * @property {string} key See pyodideSnapshotKey.
 * @property {SideModule[]} libs
 * @property {Array<[number, string, string]>} extraSlots Table slots outside the side modules' ranges: [index, library, export].
 * @property {number[]} unknownSlots
 * @property {Record<string, number>} got
 * @property {number} tableLength
 * @property {number} firstLibTable
 * @property {Record<string, string>} loadedPackages
 * @property {string[]} unknownHiwire JS references that restore as null.
 */

/**
 * A snapshot that cannot be made or restored: a Worker start failure, after
 * which the page starts a Worker that boots cold.
 * @param {string} detail
 */
const snapshotFailure = (detail) => diagnosticError(
  'WORKER_INIT', { reason: 'pyodide-snapshot', detail }, { stage: 'initialization' }
);

const textEncoder = new TextEncoder();
const textDecoder = new TextDecoder();

/** @param {string} urlOrPath */
const fileName = (urlOrPath) => urlOrPath.split(/[?#]/, 1)[0].split('/').pop() || '';

/**
 * What a snapshot depends on besides the Pyodide build (whose BUILD_ID
 * makeMemorySnapshot writes into snapshot.bin and Pyodide checks on restore):
 * the dependency wheel files, the gbdraw wheel bytes and the bytes of
 * PYTHON_HELPER_FILE.
 * @param {{ localWheels: readonly string[], gbdrawWheel: Uint8Array, helperFile: Uint8Array }} inputs
 *   localWheels: the dependency wheel URLs or paths (their file names count).
 * @returns {Promise<string>} Hex SHA-256.
 */
export const pyodideSnapshotKey = async ({ localWheels, gbdrawWheel, helperFile }) => {
  const header = textEncoder.encode(`${JSON.stringify({ format: 1, localWheels: localWheels.map(fileName) })}\n`);
  const bytes = new Uint8Array(header.length + 4 + gbdrawWheel.length + helperFile.length);
  bytes.set(header);
  new DataView(bytes.buffer).setUint32(header.length, gbdrawWheel.length, true);
  bytes.set(gbdrawWheel, header.length + 4);
  bytes.set(helperFile, header.length + 4 + gbdrawWheel.length);
  const digest = new Uint8Array(await globalThis.crypto.subtle.digest('SHA-256', bytes));
  return Array.from(digest, (value) => value.toString(16).padStart(2, '0')).join('');
};

/**
 * @param {EmscriptenFs} FS
 * @param {string} root
 * @param {string[]} [out]
 * @returns {string[]}
 */
const walkFiles = (FS, root, out = []) => {
  FS.readdir(root).forEach((entry) => {
    if (entry === '.' || entry === '..') return;
    const path = `${root}/${entry}`;
    if (FS.isDir(FS.stat(path).mode)) walkFiles(FS, path, out);
    else out.push(path);
  });
  return out;
};

// fs.bin: a little-endian u32 header length, a JSON header of [path, size]
// pairs, then the file bytes in header order.
/** @param {EmscriptenFs} FS @param {string[]} paths */
const packFiles = (FS, paths) => {
  const chunks = paths.map((path) => FS.readFile(path));
  const header = textEncoder.encode(JSON.stringify(paths.map((path, i) => [path, chunks[i].length])));
  const out = new Uint8Array(4 + header.length + chunks.reduce((sum, chunk) => sum + chunk.length, 0));
  new DataView(out.buffer).setUint32(0, header.length, true);
  out.set(header, 4);
  let offset = 4 + header.length;
  chunks.forEach((chunk) => { out.set(chunk, offset); offset += chunk.length; });
  return out;
};

/**
 * @param {EmscriptenFs} FS
 * @param {Uint8Array} bundle
 * @returns {Map<string, Uint8Array>} The written files, as views of the bundle.
 */
const unpackFiles = (FS, bundle) => {
  const headerLength = new DataView(bundle.buffer, bundle.byteOffset, bundle.byteLength).getUint32(0, true);
  /** @type {Array<[string, number]>} */
  const entries = JSON.parse(textDecoder.decode(bundle.subarray(4, 4 + headerLength)));
  let offset = 4 + headerLength;
  /** @type {Map<string, Uint8Array>} */
  const files = new Map();
  entries.forEach(([path, size]) => {
    const data = bundle.subarray(offset, offset + size);
    offset += size;
    FS.mkdirTree(path.slice(0, path.lastIndexOf('/')));
    FS.writeFile(path, data, { canOwn: false });
    files.set(path, data);
  });
  return files;
};

/** @param {PyodideRuntime} pyodide @returns {Array<[Exclude<HolderDescriptor['holder'], null>, unknown]>} */
const namedHolders = (pyodide) => [['public', pyodide], ['api', pyodide._api], ['global', globalThis]];

/**
 * @param {PyodideRuntime} pyodide
 * @param {unknown} value
 * @param {string[]} unknown Collects the values no holder names.
 * @returns {HolderDescriptor}
 */
const describeHeldObject = (pyodide, value, unknown) => {
  const holders = namedHolders(pyodide);
  for (const [holder, object] of holders) {
    if (value === object) return { holder };
  }
  for (const [holder, object] of holders) {
    if (!object || typeof object !== 'object') continue;
    for (const name of Object.getOwnPropertyNames(object)) {
      let property;
      try { property = /** @type {Record<string, unknown>} */ (object)[name]; } catch { continue; }
      if (property === value) return { holder, name };
    }
  }
  unknown.push(`${typeof value} ${String(value).slice(0, 60)}`);
  return { holder: null };
};

/**
 * Records the dynamic-linker state, the site-packages files and the Wasm
 * memory of a runtime loaded with `_makeSnapshot: true`. The caller has
 * finished its setup and collected Python garbage.
 * @param {PyodideRuntime} pyodide
 * @param {string} key See pyodideSnapshotKey.
 * @returns {{ memory: Uint8Array, state: PyodideSnapshotState, files: Uint8Array }}
 */
export const capturePyodideSnapshot = (pyodide, key) => {
  const M = pyodide._module;
  const { LDSO } = M;
  /** @type {Map<DynamicLibrary, number[]>} */
  const handlesByLibrary = new Map();
  Object.entries(LDSO.loadedLibsByHandle).forEach(([handle, library]) => {
    handlesByLibrary.set(library, [...(handlesByLibrary.get(library) || []), Number(handle)]);
  });
  /** @type {SideModule[]} */
  const libs = [];
  Object.entries(LDSO.loadedLibsByName).forEach(([name, library]) => {
    if (name === '__main__') return;
    const handles = handlesByLibrary.get(library) || [];
    const handle = handles.find((h) => M.HEAP8[h + HANDLE_OFFSETS.allocated] === 1);
    if (handle === undefined) throw snapshotFailure(`No dlopen handle records the bases of ${name}`);
    /** @param {number} offset */
    const u32 = (offset) => M.HEAPU32[(handle + offset) >>> 2];
    libs.push({
      name, handles, global: Boolean(library.global),
      memoryBase: u32(HANDLE_OFFSETS.memAddr), memorySize: u32(HANDLE_OFFSETS.memSize),
      tableBase: u32(HANDLE_OFFSETS.tableAddr), tableSize: u32(HANDLE_OFFSETS.tableSize)
    });
  });
  // Every table slot outside the side modules' own ranges that holds a
  // side-module export (GOT.func and dlsym slots), by identity.
  /** @type {Map<unknown, [string, string]>} */
  const owners = new Map();
  libs.forEach((lib) => {
    Object.entries(LDSO.loadedLibsByName[lib.name].exports || {}).forEach(([symbol, value]) => {
      if (typeof value === 'function' && !owners.has(value)) owners.set(value, [lib.name, symbol]);
    });
  });
  /** @param {number} index */
  const inLibraryRange = (index) => libs.some((lib) => index >= lib.tableBase && index < lib.tableBase + lib.tableSize);
  const tabled = libs.filter((lib) => lib.tableSize > 0);
  const firstLibTable = tabled.length ? Math.min(...tabled.map((lib) => lib.tableBase)) : M.wasmTable.length;
  /** @type {Array<[number, string, string]>} */
  const extraSlots = [];
  /** @type {number[]} */
  const unknownSlots = [];
  for (let index = firstLibTable; index < M.wasmTable.length; index += 1) {
    if (inLibraryRange(index)) continue;
    const fn = M.wasmTable.get(index);
    if (!fn) continue;
    const owner = owners.get(fn);
    if (owner) extraSlots.push([index, ...owner]);
    else unknownSlots.push(index);
  }
  /** @type {Record<string, number>} */
  const got = {};
  Object.entries(M.GOT).forEach(([symbol, global]) => { got[symbol] = Number(global.value); });
  const files = packFiles(M.FS, walkFiles(M.FS, pyodide._api.sitePackages || '/lib/python3.13/site-packages'));

  // Pyodide 0.29 serializes hiwire slot 6 (pyodide.loadedPackages) only while
  // it is empty, and loadPackage fills it: empty it around the call.
  const loaded = pyodide.loadedPackages;
  const loadedPackages = { ...loaded };
  Object.keys(loaded).forEach((name) => { delete loaded[name]; });
  /** @type {string[]} */
  const unknownHiwire = [];
  let memory;
  try {
    memory = pyodide.makeMemorySnapshot({ serializer: (value) => describeHeldObject(pyodide, value, unknownHiwire) });
  } finally {
    Object.assign(loaded, loadedPackages);
  }
  return {
    memory,
    state: {
      key, libs, extraSlots, unknownSlots, got, tableLength: M.wasmTable.length, firstLibTable,
      loadedPackages, unknownHiwire
    },
    files
  };
};

/**
 * @param {EmscriptenModule} M
 * @param {Uint8Array} binary
 * @param {SideModule} lib
 */
const instantiateSideModule = async (M, binary, lib) => {
  const metadata = M.getDylinkMetadata(binary);
  const wasmImports = /** @type {Record<string, unknown>} */ (M.LDSO.loadedLibsByName.__main__.exports);
  const { memoryBase, tableBase } = lib;
  if (metadata.tableSize !== lib.tableSize || metadata.memorySize !== lib.memorySize) {
    throw snapshotFailure(`Dylink metadata of ${lib.name} differs from the snapshot`);
  }
  if (metadata.tableSize) {
    if (M.wasmTable.length > tableBase) {
      throw snapshotFailure(`Table already longer (${M.wasmTable.length}) than ${lib.name}'s base ${tableBase}`);
    }
    if (M.wasmTable.length < tableBase) M.wasmTable.grow(tableBase - M.wasmTable.length);
    M.wasmTable.grow(metadata.tableSize);
  }
  /** @type {Record<string, unknown> | null} */
  const localScope = lib.global ? null : {};
  /** @type {Record<string, unknown>} */
  let moduleExports = {};
  /** @param {string} symbol @returns {WasmFunction | undefined} */
  const resolveSymbol = (symbol) => {
    const resolved = M.resolveGlobalSymbol(symbol).sym || localScope?.[symbol] || moduleExports[symbol];
    return typeof resolved === 'function' ? /** @type {WasmFunction} */ (resolved) : undefined;
  };
  /** @type {Record<string, unknown>} */
  const stubs = {};
  const env = new Proxy(stubs, {
    get(target, property) {
      if (typeof property !== 'string') return undefined;
      if (property === '__memory_base') return memoryBase;
      if (property === '__table_base') return tableBase;
      const imported = /** @type {{ stub?: boolean } | undefined} */ (wasmImports[property]);
      if (property in wasmImports && !imported?.stub) return imported;
      if (!(property in target)) {
        /** @type {WasmFunction | undefined} */
        let resolved;
        target[property] = (/** @type {unknown[]} */ ...args) => {
          resolved ||= resolveSymbol(property);
          if (!resolved) throw snapshotFailure(`Dynamic linking error: cannot resolve symbol ${property}`);
          return resolved(...args);
        };
      }
      return target[property];
    }
  });
  // Proxies resolve each import by name when the module asks for it.
  const imports = /** @type {WebAssembly.Imports} */ (/** @type {unknown} */ ({
    'GOT.mem': new Proxy({}, M.GOTHandler),
    'GOT.func': new Proxy({}, M.GOTHandler),
    env,
    wasi_snapshot_preview1: env
  }));
  const { instance } = await WebAssembly.instantiate(/** @type {Uint8Array<ArrayBuffer>} */ (binary), imports);
  // Relocate the exports without updateGOT (its addFunction calls would add
  // table slots in another order); the GOT values come from the snapshot.
  moduleExports = {};
  Object.entries(instance.exports).forEach(([name, value]) => {
    const raw = value instanceof WebAssembly.Global ? value.value : value;
    moduleExports[name] = typeof raw === 'number' ? raw + memoryBase : raw;
  });
  if ('__start_em_asm' in moduleExports || Object.keys(moduleExports).some((name) => name.startsWith('__em_js__'))) {
    throw snapshotFailure(`${lib.name} uses EM_ASM/EM_JS, which the snapshot restore does not support`);
  }
  // __wasm_apply_data_relocs and __wasm_call_ctors stay uncalled: their
  // effects are in the snapshot memory. Instantiation rewrote the data
  // segments; the caller puts the snapshot bytes back.
  const library = M.newDSO(lib.name, undefined, 'loading');
  library.refcount = Infinity;
  library.global = lib.global;
  if (library.global) M.mergeLibSymbols(moduleExports, lib.name);
  library.exports = moduleExports;
  lib.handles.forEach((handle) => { M.LDSO.loadedLibsByHandle[handle] = library; });
};

/**
 * @param {string} url
 * @returns {Promise<Uint8Array>}
 */
const fetchBytes = async (url) => {
  const response = await fetch(url);
  if (!response.ok) throw snapshotFailure(`Pyodide snapshot input missing (${response.status}) at ${url}`);
  return new Uint8Array(await response.arrayBuffer());
};

/**
 * A fetched and checked snapshot, ready for Pyodide's loader.
 * @typedef {object} PendingPyodideSnapshot
 * @property {Uint8Array} memory the `_loadSnapshot` option of Pyodide's loader.
 * @property {(descriptor: HolderDescriptor) => unknown} deserialize the `_snapshotDeserializer` option of Pyodide's loader.
 * @property {(pyodide: PyodideRuntime) => Promise<void>} complete Restores the
 *   side modules and the files into the runtime the loader returned, then
 *   drops every download buffer.
 */

/**
 * Fetches the snapshot in `folderUrl` and checks it against the wheels and the
 * helper source this Worker would install. state.json is fetched first, so a
 * missing or blocked snapshot fails before any other request.
 * @param {{ folderUrl: string, gbdrawWheelUrl: string, localWheelUrls: readonly string[] }} inputs
 * @returns {Promise<PendingPyodideSnapshot>}
 */
export const fetchPyodideSnapshot = async ({ folderUrl, gbdrawWheelUrl, localWheelUrls }) => {
  /** @param {string} name */
  const file = (name) => new URL(name, folderUrl.endsWith('/') ? folderUrl : `${folderUrl}/`).toString();
  /** @type {PyodideSnapshotState} */
  const state = JSON.parse(textDecoder.decode(await fetchBytes(file(PYODIDE_SNAPSHOT_FILES.state))));
  const memoryPromise = fetchBytes(file(PYODIDE_SNAPSHOT_FILES.memory));
  const filesPromise = fetchBytes(file(PYODIDE_SNAPSHOT_FILES.files));
  // A rejected download must not surface as an unhandled rejection when the
  // key check below throws first.
  memoryPromise.catch(() => {});
  filesPromise.catch(() => {});
  const [gbdrawWheel, helperFile] = await Promise.all([fetchBytes(gbdrawWheelUrl), fetchBytes(PYTHON_HELPER_FILE.href)]);
  const key = await pyodideSnapshotKey({ localWheels: localWheelUrls, gbdrawWheel, helperFile });
  if (state.key !== key) throw snapshotFailure('Pyodide snapshot was built from another wheel or helper source');
  // Released by complete(): about 170 MB.
  /** @type {{ memory: Uint8Array, files: Uint8Array, state: PyodideSnapshotState } | null} */
  let downloads = { memory: await memoryPromise, files: await filesPromise, state };

  // The deserializer runs inside the loader, before it returns the API whose
  // objects some descriptors name: those resolve through a late proxy.
  /** @type {PyodideRuntime | null} */
  let runtime = null;
  /** @param {HolderDescriptor} descriptor */
  const resolveNow = (descriptor) => {
    const holder = new Map(namedHolders(/** @type {PyodideRuntime} */ (runtime))).get(/** @type {Exclude<HolderDescriptor['holder'], null>} */ (descriptor.holder));
    return descriptor.name ? /** @type {Record<string, unknown>} */ (holder)[descriptor.name] : holder;
  };
  /** @param {HolderDescriptor} descriptor */
  const deserialize = (descriptor) => {
    if (!descriptor?.holder) return null;
    if (runtime) return resolveNow(descriptor);
    if (descriptor.holder === 'global') {
      return descriptor.name ? /** @type {Record<string, unknown>} */ (globalThis)[descriptor.name] : globalThis;
    }
    const target = () => /** @type {object} */ (resolveNow(descriptor));
    return new Proxy(function lateBound() {}, {
      get: (_target, property) => Reflect.get(target(), property),
      has: (_target, property) => Reflect.has(target(), property),
      apply: (_target, thisArg, args) => Reflect.apply(/** @type {WasmFunction} */ (target()), thisArg, args)
    });
  };

  /** @param {PyodideRuntime} pyodide */
  const complete = async (pyodide) => {
    const { memory, files, state } = /** @type {NonNullable<typeof downloads>} */ (downloads);
    downloads = null;
    runtime = pyodide;
    const M = pyodide._module;
    const written = unpackFiles(M.FS, files);
    const heapOffset = new Uint32Array(memory.buffer, memory.byteOffset, 4)[1];
    for (const lib of state.libs) {
      await instantiateSideModule(M, written.get(lib.name) || M.FS.readFile(lib.name), lib);
    }
    // Put back the data regions the instantiations rewrote.
    state.libs.forEach((lib) => {
      if (!lib.memorySize) return;
      const start = heapOffset + lib.memoryBase;
      M.HEAP8.set(memory.subarray(start, start + lib.memorySize), lib.memoryBase);
    });
    if (M.wasmTable.length < state.tableLength) M.wasmTable.grow(state.tableLength - M.wasmTable.length);
    state.extraSlots.forEach(([index, library, symbol]) => {
      M.setWasmTableEntry(index, M.LDSO.loadedLibsByName[library].exports[symbol]);
    });
    Object.entries(state.got).forEach(([symbol, value]) => {
      M.GOT[symbol] ||= new WebAssembly.Global({ value: 'i32', mutable: true });
      M.GOT[symbol].value = value;
    });
    M.getFunctionAddress(() => {});
    M.updateTableMap(0, M.wasmTable.length);
    Object.assign(pyodide.loadedPackages, state.loadedPackages);
    // The loader keeps its options: drop its reference to the memory too.
    delete pyodide._api.config._loadSnapshot;
    delete pyodide._api.config._snapshotDeserializer;
  };

  return { memory: downloads.memory, deserialize, complete };
};
