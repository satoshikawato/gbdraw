Place the LOSAT WebAssembly binaries here:

- `losat.wasm`: serial browser artifact with the direct `losat_web_*` API.
- `losat-threaded.wasm`: threaded command artifact for cross-origin-isolated browsers.

Suggested serial build (from the LOSAT repo root):
```
cd /mnt/c/Users/genom/GitHub/LOSAT
rustup target add wasm32-wasip1
cargo build --release --target wasm32-wasip1 --no-default-features --lib
cp target/wasm32-wasip1/release/LOSAT.wasm /mnt/c/Users/genom/GitHub/gbdraw/gbdraw/web/wasm/losat/losat.wasm
```

Suggested threaded build:
```
cd /mnt/c/Users/genom/GitHub/LOSAT
rustup target add wasm32-wasip1-threads
cargo build --release --target wasm32-wasip1-threads --features wasm-threads
cp target/wasm32-wasip1-threads/release/LOSAT.wasm /mnt/c/Users/genom/GitHub/gbdraw/gbdraw/web/wasm/losat/losat-threaded.wasm
```

Notes:
- The web app keeps `losat.wasm` as the serial fallback.
- `losat-threaded.wasm` must export `_start` and `wasi_thread_start`, import
  shared memory from `env.memory`, and import `wasi.thread-spawn`.
- Threaded browser execution also requires COOP/COEP response headers so
  `SharedArrayBuffer` is available.
- For threaded LOSAT, `--num-threads N` is treated as N compute threads. The
  browser wrapper runs the main job worker plus `N - 1` WASI thread workers.
