// PROTOTYPE (proto/pyodide-snapshot, not for dev): capture and restore a
// Pyodide 0.29.0 memory snapshot that includes loaded C extension modules.
//
// Stock makeMemorySnapshot()/_loadSnapshot save only the Wasm linear memory.
// Two things live outside it and are lost on restore:
//   1. the side modules (.so) loaded by loadPackage: their code, their slots
//      in the function table, the GOT globals and the JS dynamic-linker
//      tables (LDSO); a call into numpy then traps ("table index is out of
//      bounds");
//   2. MEMFS (site-packages files): imported modules keep working, but a
//      later import or a package-data read (gbdraw/data/*.toml) fails.
// The capture records the dynamic-linker state and the site-packages files;
// the restore re-instantiates each side module at its recorded memory and
// table base, puts the snapshot bytes of its data region back, and rebuilds
// the table slots and GOT values exactly as they were.

const HANDLE_OFFSETS = Object.freeze({ allocated: 8, memAddr: 12, memSize: 16, tableAddr: 20, tableSize: 24 });

const walkFiles = (FS, root, out = []) => {
  FS.readdir(root).forEach((entry) => {
    if (entry === '.' || entry === '..') return;
    const path = `${root}/${entry}`;
    const stat = FS.stat(path);
    if (FS.isDir(stat.mode)) walkFiles(FS, path, out);
    else out.push(path);
  });
  return out;
};

const packFiles = (FS, paths) => {
  const chunks = paths.map((path) => FS.readFile(path));
  const header = new TextEncoder().encode(JSON.stringify(paths.map((path, i) => [path, chunks[i].length])));
  const total = 4 + header.length + chunks.reduce((sum, c) => sum + c.length, 0);
  const out = new Uint8Array(total);
  new DataView(out.buffer).setUint32(0, header.length, true);
  out.set(header, 4);
  let offset = 4 + header.length;
  chunks.forEach((c) => { out.set(c, offset); offset += c.length; });
  return out;
};

const unpackFiles = (FS, bundle) => {
  const view = new DataView(bundle.buffer, bundle.byteOffset, bundle.byteLength);
  const headerLength = view.getUint32(0, true);
  const entries = JSON.parse(new TextDecoder().decode(bundle.subarray(4, 4 + headerLength)));
  let offset = 4 + headerLength;
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

/**
 * Call after setup and before makeMemorySnapshot(). Returns the JSON-able
 * linker state and the site-packages bundle.
 */
export const captureSnapshotState = (pyodide) => {
  const M = pyodide._module;
  const { LDSO } = M;
  const handlesByDso = new Map();
  Object.entries(LDSO.loadedLibsByHandle).forEach(([handle, dso]) => {
    if (!handlesByDso.has(dso)) handlesByDso.set(dso, []);
    handlesByDso.get(dso).push(Number(handle));
  });
  const libs = [];
  Object.entries(LDSO.loadedLibsByName).forEach(([name, dso]) => {
    if (name === '__main__') return;
    const handles = handlesByDso.get(dso) || [];
    const handle = handles.find((h) => M.HEAP8[h + HANDLE_OFFSETS.allocated] === 1);
    if (handle === undefined) throw new Error(`No dlopen handle records the bases of ${name}`);
    const u32 = (offset) => M.HEAPU32[(handle + offset) >>> 2];
    libs.push({
      name, handles, global: Boolean(dso.global),
      memoryBase: u32(HANDLE_OFFSETS.memAddr), memorySize: u32(HANDLE_OFFSETS.memSize),
      tableBase: u32(HANDLE_OFFSETS.tableAddr), tableSize: u32(HANDLE_OFFSETS.tableSize)
    });
  });
  // Identify every table slot outside the side modules' own ranges that
  // holds a side-module export (GOT.func and dlsym slots).
  const owner = new Map();
  libs.forEach((lib) => {
    const exports = LDSO.loadedLibsByName[lib.name].exports || {};
    Object.entries(exports).forEach(([sym, value]) => {
      if (typeof value === 'function' && !owner.has(value)) owner.set(value, [lib.name, sym]);
    });
  });
  const inLibRange = (idx) => libs.some((l) => idx >= l.tableBase && idx < l.tableBase + l.tableSize);
  const tabled = libs.filter((l) => l.tableSize > 0);
  const firstLibTable = tabled.length ? Math.min(...tabled.map((l) => l.tableBase)) : M.wasmTable.length;
  const extraSlots = [];
  const unknownSlots = [];
  for (let idx = firstLibTable; idx < M.wasmTable.length; idx += 1) {
    if (inLibRange(idx)) continue;
    const fn = M.wasmTable.get(idx);
    if (!fn) continue;
    const id = owner.get(fn);
    if (id) extraSlots.push([idx, ...id]);
    else unknownSlots.push(idx);
  }
  const got = {};
  Object.entries(M.GOT).forEach(([sym, g]) => { got[sym] = g.value; });
  const sitePackages = pyodide._api?.sitePackages || '/lib/python3.13/site-packages';
  const fsBundle = packFiles(M.FS, walkFiles(M.FS, sitePackages));
  return {
    state: {
      libs, extraSlots, unknownSlots, got, tableLength: M.wasmTable.length, firstLibTable,
      loadedPackages: { ...pyodide.loadedPackages }
    },
    fsBundle
  };
};

/**
 * Pyodide 0.29.0 serializes hiwire slot 6 (pyodide.loadedPackages) only when
 * it is empty; loadPackage fills it. Empty it around makeMemorySnapshot().
 */
// JS objects that Python still references (micropip's pyfetch leftovers, the
// public API, bound API functions) need a serializer. Known objects are named
// by where they live; anything else restores as null and is counted.
const namedHolders = (pyodide) => [
  ['public', pyodide],
  ['api', pyodide._api],
  ['global', globalThis]
];
const describeObject = (pyodide, obj, unknown) => {
  for (const [holder, value] of namedHolders(pyodide)) {
    if (obj === value) return { holder };
  }
  for (const [holder, value] of namedHolders(pyodide)) {
    if (!value) continue;
    for (const name of Object.getOwnPropertyNames(value)) {
      let prop;
      try { prop = value[name]; } catch { continue; }
      if (prop === obj) return { holder, name };
    }
  }
  unknown.push(`${typeof obj} ${obj?.constructor?.name || ''} ${String(obj).slice(0, 60)}`);
  return { holder: null };
};
export const snapshotDeserializer = (pyodide) => (descriptor) => {
  if (!descriptor?.holder) return null;
  const holder = new Map(namedHolders(pyodide)).get(descriptor.holder);
  return descriptor.name ? holder[descriptor.name] : holder;
};

export const makeSnapshotWithPackages = (pyodide) => {
  pyodide.runPython('import gc; gc.collect()');
  const loaded = pyodide.loadedPackages;
  const saved = { ...loaded };
  Object.keys(loaded).forEach((k) => { delete loaded[k]; });
  const unknown = [];
  try {
    const bytes = pyodide.makeMemorySnapshot({ serializer: (obj) => describeObject(pyodide, obj, unknown) });
    return { bytes, unknown };
  } finally {
    Object.assign(loaded, saved);
  }
};

const instantiateAt = async (M, binary, lib, files) => {
  const metadata = M.getDylinkMetadata(binary);
  const wasmImports = M.LDSO.loadedLibsByName.__main__.exports;
  const { memoryBase, tableBase } = lib;
  if (metadata.tableSize !== lib.tableSize || metadata.memorySize !== lib.memorySize) {
    throw new Error(`Dylink metadata of ${lib.name} differs from the snapshot`);
  }
  if (metadata.tableSize) {
    if (M.wasmTable.length > tableBase) {
      throw new Error(`Table already longer (${M.wasmTable.length}) than ${lib.name}'s base ${tableBase}`);
    }
    if (M.wasmTable.length < tableBase) M.wasmTable.grow(tableBase - M.wasmTable.length);
    M.wasmTable.grow(metadata.tableSize);
  }
  const localScope = lib.global ? null : {};
  let moduleExports;
  const resolveSymbol = (sym) => {
    let resolved = M.resolveGlobalSymbol(sym).sym;
    if (!resolved && localScope) resolved = localScope[sym];
    if (!resolved) resolved = moduleExports[sym];
    return resolved;
  };
  const env = new Proxy({}, {
    get(stubs, prop) {
      if (prop === '__memory_base') return memoryBase;
      if (prop === '__table_base') return tableBase;
      if (prop in wasmImports && !wasmImports[prop].stub) return wasmImports[prop];
      if (!(prop in stubs)) {
        let resolved;
        stubs[prop] = (...args) => {
          resolved ||= resolveSymbol(prop);
          if (!resolved) throw new Error(`Dynamic linking error: cannot resolve symbol ${prop}`);
          return resolved(...args);
        };
      }
      return stubs[prop];
    }
  });
  const info = {
    'GOT.mem': new Proxy({}, M.GOTHandler),
    'GOT.func': new Proxy({}, M.GOTHandler),
    env,
    wasi_snapshot_preview1: env
  };
  const { instance } = await WebAssembly.instantiate(binary, info);
  // Relocate exports without updateGOT (its addFunction calls would add table
  // slots in a different order); GOT values come from the snapshot below.
  moduleExports = {};
  Object.entries(instance.exports).forEach(([name, value]) => {
    let v = value;
    if (typeof v === 'object') v = v.value;
    if (typeof v === 'number') v += memoryBase;
    moduleExports[name] = v;
  });
  if ('__start_em_asm' in moduleExports || Object.keys(moduleExports).some((n) => n.startsWith('__em_js__'))) {
    throw new Error(`${lib.name} uses EM_ASM/EM_JS; the prototype does not restore them`);
  }
  // Skip __wasm_apply_data_relocs and __wasm_call_ctors: their effects are in
  // the snapshot memory. Instantiation rewrote the data segments; the caller
  // puts the snapshot bytes back.
  const dso = M.newDSO(lib.name, undefined, 'loading');
  dso.refcount = Infinity;
  dso.global = lib.global;
  if (dso.global) M.mergeLibSymbols(moduleExports, lib.name);
  dso.exports = moduleExports;
  lib.handles.forEach((h) => { M.LDSO.loadedLibsByHandle[h] = dso; });
  return dso;
};

/**
 * Call right after loadPyodide({ _loadSnapshot }). `snapshot` is the same
 * byte array given to _loadSnapshot.
 */
export const restoreSnapshotState = async (pyodide, state, fsBundle, snapshot) => {
  const M = pyodide._module;
  const timings = {};
  let t0 = performance.now();
  const files = unpackFiles(M.FS, fsBundle);
  timings.fsUnpack = performance.now() - t0;
  const header = new Uint32Array(snapshot.buffer, snapshot.byteOffset, 4);
  const heapOffset = header[1];
  t0 = performance.now();
  if (M.wasmTable.length > state.firstLibTable) {
    throw new Error(`Restored table length ${M.wasmTable.length} exceeds the first side-module base ${state.firstLibTable}`);
  }
  for (const lib of state.libs) {
    const binary = files.get(lib.name) || M.FS.readFile(lib.name);
    await instantiateAt(M, binary, lib, files);
  }
  timings.instantiate = performance.now() - t0;
  t0 = performance.now();
  // Put back the data regions the instantiations rewrote.
  state.libs.forEach((lib) => {
    if (!lib.memorySize) return;
    const start = heapOffset + lib.memoryBase;
    M.HEAP8.set(snapshot.subarray(start, start + lib.memorySize), lib.memoryBase);
  });
  if (M.wasmTable.length < state.tableLength) M.wasmTable.grow(state.tableLength - M.wasmTable.length);
  state.extraSlots.forEach(([idx, libName, sym]) => {
    M.setWasmTableEntry(idx, M.LDSO.loadedLibsByName[libName].exports[sym]);
  });
  Object.entries(state.got).forEach(([sym, value]) => {
    M.GOT[sym] ||= new WebAssembly.Global({ value: 'i32', mutable: true });
    M.GOT[sym].value = value;
  });
  M.getFunctionAddress(() => {});
  M.updateTableMap(0, M.wasmTable.length);
  Object.assign(pyodide.loadedPackages, state.loadedPackages);
  timings.relink = performance.now() - t0;
  return timings;
};
