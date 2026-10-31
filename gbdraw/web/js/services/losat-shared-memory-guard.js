// @ts-check
// Guard every memory.copy and memory.fill of a threaded LOSAT module.
//
// A threaded module runs one instance per worker on one shared memory. In V8,
// a worker's view of that memory is updated some time after another worker
// grows it, and memory.copy and memory.fill are checked against that view, so
// a copy into memory that another thread has just grown traps "memory access
// out of bounds" (Chromium 156, the eight Vibrio chromosomes: 14 traps in 42
// jobs; with this guard, none in 42). The guard turns each bulk operation into
// a call of a helper that compares the range with `memory.size`, runs
// `memory.grow(0)` when the range lies outside it (which updates the calling
// worker's view and allocates nothing), and then runs the original operation
// once, so a range that is still invalid traps as before.
//
// Ported from LOSAT `LOSAT/tests/wasi_shared_memory.js` (`guardSharedMemory`,
// LOSAT dev 9e519214a), which LOSAT's Node host applies to the same modules.
// Binary encoding: WebAssembly Core spec, binary/modules and binary/instructions.
// The reader accepts the wasm32 instruction set of LOSAT's threaded builds and
// declines anything else instead of guessing instruction boundaries.

/** Thrown inside this module only: the module cannot be guarded. */
class GuardDeclined {}

/** @param {number} value */
const uleb = (value) => {
  /** @type {number[]} */
  const bytes = [];
  let rest = value;
  do {
    const byte = rest & 127;
    rest >>>= 7;
    bytes.push(byte | (rest ? 128 : 0));
  } while (rest);
  return Uint8Array.from(bytes);
};

/** @param {readonly Uint8Array<ArrayBuffer>[]} parts */
const concat = (parts) => {
  const out = new Uint8Array(parts.reduce((total, part) => total + part.length, 0));
  let offset = 0;
  for (const part of parts) {
    out.set(part, offset);
    offset += part.length;
  }
  return out;
};

/** @param {Uint8Array<ArrayBuffer>} bytes */
const sized = (bytes) => concat([uleb(bytes.length), bytes]);

const VALUE_TYPES = new Set([0x7f, 0x7e, 0x7d, 0x7c, 0x7b, 0x70, 0x6f]);
const BLOCK_OR_CONST = new Set([0x02, 0x03, 0x04, 0x41, 0x42, 0xd0]);
const ONE_INDEX = new Set([0x0c, 0x0d, 0x10, 0x12, 0x14, 0x15, 0xd2, 0xd5, 0xd6]);
const NO_IMMEDIATE = new Set([0x00, 0x01, 0x05, 0x0b, 0x0f, 0x1a, 0x1b, 0xd1, 0xd3, 0xd4]);
// Reserved SIMD holes (V8 wasm-opcodes.h): declined even if an engine assigns them.
const SIMD_HOLES = new Set([0x9a, 0xa2, 0xa5, 0xa6, 0xaf, 0xb0, 0xb2, 0xb3, 0xb4, 0xbb,
  0xc2, 0xc5, 0xc6, 0xcf, 0xd0, 0xd2, 0xd3, 0xd4, 0xe2, 0xee]);

class Reader {
  /** @param {Uint8Array<ArrayBuffer>} bytes */
  constructor(bytes) {
    this.bytes = bytes;
    this.pos = 0;
  }

  byte() {
    if (this.pos >= this.bytes.length) throw new GuardDeclined();
    return this.bytes[this.pos++];
  }

  /** @param {number} size */
  take(size) {
    if (size > this.bytes.length - this.pos) throw new GuardDeclined();
    const result = this.bytes.subarray(this.pos, this.pos + size);
    this.pos += size;
    return result;
  }

  uint() {
    let value = 0;
    for (let shift = 0; shift < 35; shift += 7) {
      const byte = this.byte();
      value += (byte & 127) * 2 ** shift;
      if (!(byte & 128) && value <= 0xffffffff) return value;
    }
    throw new GuardDeclined();
  }

  signed() {
    for (let i = 0; i < 10; i += 1) if (!(this.byte() & 128)) return;
    throw new GuardDeclined();
  }

  type() {
    if (!VALUE_TYPES.has(this.byte())) throw new GuardDeclined();
  }

  memory() {
    if (this.uint() !== 0) throw new GuardDeclined();
  }

  memarg() {
    if (this.uint() & ~0x3f) throw new GuardDeclined();
    this.uint();
  }

  end() {
    if (this.pos !== this.bytes.length) throw new GuardDeclined();
  }
}

/**
 * Read one instruction.
 * @param {Reader} r
 * @returns {'copy' | 'fill' | null} the bulk operation it is, if any
 */
const instruction = (r) => {
  const op = r.byte();
  if (BLOCK_OR_CONST.has(op)) r.signed();
  else if (ONE_INDEX.has(op) || (op >= 0x20 && op <= 0x26)) r.uint();
  else if (op === 0x0e) {
    const count = r.uint();
    for (let i = 0; i <= count; i += 1) r.uint();
  } else if (op === 0x11 || op === 0x13) { r.uint(); r.uint(); }
  else if (op === 0x1c) {
    const count = r.uint();
    for (let i = 0; i < count; i += 1) r.type();
  } else if (op >= 0x28 && op <= 0x3e) r.memarg();
  else if (op === 0x3f || op === 0x40) r.memory();
  else if (op === 0x43 || op === 0x44) r.take(op === 0x43 ? 4 : 8);
  else if (op === 0xfc) {
    const sub = r.uint();
    if (sub <= 7) return null;
    if (sub === 8) { r.uint(); r.memory(); }
    else if (sub === 9 || sub === 13 || (sub >= 15 && sub <= 17)) r.uint();
    else if (sub === 10) { r.memory(); r.memory(); return 'copy'; }
    else if (sub === 11) { r.memory(); return 'fill'; }
    else if (sub === 12 || sub === 14) { r.uint(); r.uint(); }
    else throw new GuardDeclined();
  } else if (op === 0xfd) {
    const sub = r.uint();
    if (SIMD_HOLES.has(sub)) throw new GuardDeclined();
    if (sub <= 11 || sub === 92 || sub === 93) r.memarg();
    else if (sub === 12 || sub === 13) r.take(16);
    else if (sub >= 21 && sub <= 34) r.take(1);
    else if (sub >= 84 && sub <= 91) { r.memarg(); r.take(1); }
    else if (sub > 0x113) throw new GuardDeclined();
  } else if (op === 0xfe) {
    const sub = r.uint();
    if (sub <= 2 || (sub >= 0x10 && sub <= 0x4e)) r.memarg();
    else if (sub === 3) {
      if (r.byte() !== 0) throw new GuardDeclined();
    } else throw new GuardDeclined();
  } else if (!NO_IMMEDIATE.has(op) && !(op >= 0x45 && op <= 0xc4)) {
    throw new GuardDeclined();
  }
  return null;
};

/**
 * The helper that replaces one bulk operation:
 * (param dst i32) (param src-or-value i32) (param len i32) (local bytes i64).
 * @param {boolean} copy
 */
const guardedBody = (copy) => {
  /** @param {number} index */
  const outside = (index) => [0x20, index, 0xad, 0x20, 2, 0xad, 0x7c, 0x20, 3, 0x56];
  return Uint8Array.from([
    1, 1, 0x7e,
    0x3f, 0, 0xad, 0x42, 16, 0x86, 0x21, 3, // bytes = u64(memory.size) << 16
    ...outside(0), ...(copy ? [...outside(1), 0x72] : []),
    0x04, 0x40, 0x41, 0, 0x40, 0, 0x1a, 0x0b, // if outside: grow(0); drop
    0x20, 0, 0x20, 1, 0x20, 2,
    0xfc, ...(copy ? [10, 0, 0] : [11, 0]), 0x0b
  ]);
};

/**
 * @typedef {object} SharedMemoryGuardResult
 * @property {Uint8Array<ArrayBuffer>} bytes the module to compile: the input itself when it
 *   imports no shared memory or has no bulk operation
 * @property {number} copy memory.copy operations guarded
 * @property {number} fill memory.fill operations guarded
 */

/**
 * @param {Uint8Array<ArrayBuffer>} bytes
 * @returns {SharedMemoryGuardResult}
 */
const guard = (bytes) => {
  const reader = new Reader(bytes);
  const header = reader.take(8);
  /** @type {{ id: number, bytes: Uint8Array<ArrayBuffer> }[]} */
  const sections = [];
  while (reader.pos < bytes.length) {
    const id = reader.byte();
    sections.push({ id, bytes: reader.take(reader.uint()) });
  }
  /** @param {number} id */
  const find = (id) => sections.find((s) => s.id === id);
  const importSection = find(2);
  if (!importSection) return { bytes, copy: 0, fill: 0 };
  let imports = 0;
  let sharedMemories = 0;
  const ir = new Reader(importSection.bytes);
  const importCount = ir.uint();
  for (let i = 0; i < importCount; i += 1) {
    ir.take(ir.uint()); ir.take(ir.uint());
    const kind = ir.byte();
    if (kind === 0) { imports += 1; ir.uint(); }
    else if (kind === 1 || kind === 2) {
      if (kind === 1) ir.type();
      const flags = ir.uint();
      if (kind === 2 && flags !== 3) return { bytes, copy: 0, fill: 0 }; // not a shared wasm32 memory
      if (kind === 1 && flags > 1) throw new GuardDeclined();
      ir.uint();
      if (flags & 1) ir.uint();
      if (kind === 2) sharedMemories += 1;
    } else if (kind === 3) { ir.type(); ir.byte(); }
    else throw new GuardDeclined();
  }
  ir.end();
  if (sharedMemories === 0) return { bytes, copy: 0, fill: 0 };
  const types = find(1);
  const functions = find(3);
  const code = find(10);
  if (sharedMemories !== 1 || find(5) || !types || !functions || !code) throw new GuardDeclined();
  const tr = new Reader(types.bytes);
  const typeCount = tr.uint();
  // One type index per entry (recursive and GC type groups are declined).
  for (let i = 0; i < typeCount; i += 1) {
    if (tr.byte() !== 0x60) throw new GuardDeclined();
    for (let j = 0; j < 2; j += 1) {
      const count = tr.uint();
      for (let k = 0; k < count; k += 1) tr.type();
    }
  }
  tr.end();
  const fr = new Reader(functions.bytes);
  const functionCount = fr.uint();
  const functionEntries = functions.bytes.subarray(fr.pos);
  const cr = new Reader(code.bytes);
  if (cr.uint() !== functionCount) throw new GuardDeclined();
  const counts = { copy: 0, fill: 0 };
  /** @type {Uint8Array<ArrayBuffer>[]} */
  const bodies = [];
  const helperCall = {
    fill: concat([Uint8Array.of(0x10), uleb(imports + functionCount)]),
    copy: concat([Uint8Array.of(0x10), uleb(imports + functionCount + 1)])
  };
  for (let i = 0; i < functionCount; i += 1) {
    const body = cr.take(cr.uint());
    const r = new Reader(body);
    const locals = r.uint();
    for (let j = 0; j < locals; j += 1) { r.uint(); r.type(); }
    /** @type {Uint8Array<ArrayBuffer>[]} */
    const chunks = [];
    let start = 0;
    while (r.pos < body.length) {
      const before = r.pos;
      const bulk = instruction(r);
      if (bulk) {
        chunks.push(body.subarray(start, before), helperCall[bulk]);
        counts[bulk] += 1;
        start = r.pos;
      }
    }
    chunks.push(body.subarray(start));
    bodies.push(sized(concat(chunks)));
  }
  cr.end();
  if (counts.copy + counts.fill === 0) return { bytes, ...counts };
  const originalTypes = new Reader(types.bytes);
  originalTypes.uint();
  types.bytes = concat([uleb(typeCount + 1), types.bytes.subarray(originalTypes.pos),
    Uint8Array.of(0x60, 3, 0x7f, 0x7f, 0x7f, 0)]);
  functions.bytes = concat([uleb(functionCount + 2), functionEntries, uleb(typeCount), uleb(typeCount)]);
  code.bytes = concat([uleb(functionCount + 2), ...bodies, sized(guardedBody(false)), sized(guardedBody(true))]);
  const guarded = concat([header, ...sections.map((s) => concat([Uint8Array.of(s.id), sized(s.bytes)]))]);
  if (!WebAssembly.validate(guarded)) throw new GuardDeclined();
  return { bytes: guarded, ...counts };
};

/**
 * Guard the bulk memory operations of a module that imports a shared memory.
 * A module without a shared memory or without bulk operations comes back as
 * it is.
 * @param {Uint8Array<ArrayBuffer>} bytes a valid Wasm module
 * @returns {SharedMemoryGuardResult | null} null when the reader declines the
 *   module (an encoding it does not know); the caller then compiles it as it is
 */
export const guardSharedMemoryBulkOps = (bytes) => {
  try {
    return guard(bytes);
  } catch (error) {
    if (error instanceof GuardDeclined) return null;
    throw error;
  }
};
