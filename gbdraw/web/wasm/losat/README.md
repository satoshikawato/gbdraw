The LOSAT WebAssembly files of the Web app:

- `losat.wasm`: serial WASI reactor with the direct `losat_web_*` API (LOSAT
  `docs/web/abi_v1.md`). The host calls `_initialize` once after instantiation.
- `losat-threaded.wasm`: threaded WASI command for cross-origin-isolated
  browsers. It exports `_start` and `wasi_thread_start`, imports shared memory
  from `env.memory`, and imports `wasi.thread-spawn`. LOSAT's shared-memory guard
  has rewritten its memory operations (its record is the custom section
  `losat.shared_memory_guard`); the Web app compiles the file as it is.

`losat-build.json` records their build: the LOSAT commit, the build commands,
the toolchain, each file's sha256 and size, and the page limits of the memory
the threaded file imports and its guard record. `tests/test_losat_wasm_build.py`
fails when a file does not match it or the threaded file has no current guard. The Web app records the build's `version` with each search,
and sizes each threaded job's shared memory from those limits
(`js/services/losat-build.generated.js`).

## Rebuild

From a clean checkout of `https://github.com/satoshikawato/LOSAT` at the chosen
commit (its `rust-toolchain.toml` selects the toolchain and both Wasm targets):

```bash
cd "$LOSAT_SRC/LOSAT"
export CARGO_TARGET_DIR="$TARGET_DIR"
cargo build --release --locked --lib --target wasm32-wasip1 --no-default-features
cargo build --release --locked --bin losat --target wasm32-wasip1-threads --features wasm-threads
cd "$GBDRAW_SRC"
python tools/record_losat_wasm_build.py --losat-src "$LOSAT_SRC" --target-dir "$TARGET_DIR"
```

The tool copies `wasm32-wasip1/release/LOSAT.wasm` to `losat.wasm`, writes
`wasm32-wasip1-threads/release/losat.wasm` through LOSAT's
`tests/wasi_shared_memory.js guard` to `losat-threaded.wasm` (cargo's output is
not guarded, and the tool stops when the result has no current guard record), and
rewrites `losat-build.json` and the generated module.

Panic messages in the files hold the absolute paths of the build (the Cargo
registry and LOSAT's `src/cli.rs`), so bytes depend on where it ran. For the same
bytes as LOSAT's `WASI threading integration` artifact
`losat-threaded-command.wasm`, build with `CARGO_HOME=/home/runner/.cargo` in a
checkout at `/home/runner/work/LOSAT/LOSAT`; `--remap-path-prefix` is not enough,
because it changes the file name in the module's name section.

Notes:

- The Web app keeps `losat.wasm` as the serial fallback.
- Threaded browser execution needs COOP/COEP response headers so that
  `SharedArrayBuffer` is available.
- For threaded LOSAT, `-num_threads N` is N compute threads. The browser
  wrapper runs the main job worker plus `N - 1` WASI thread workers.
- Both files read LOSAT's NCBI option spelling (`-task`, `-max_target_seqs`).
  The Web keeps the raw-cache spelling (`--task`) in its args and spells them
  with `LOSAT_SEARCH_FLAGS` (generated from `gbdraw/comparisons/losat_runtime.py`).
