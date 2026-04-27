Place the LOSAT WebAssembly binary here as `losat.wasm`.

Suggested build (from the LOSAT repo root):
```
cd /mnt/c/Users/kawato/Documents/GitHub/LOSAT/LOSAT
rustup target add wasm32-wasip1
cargo build --release --target wasm32-wasip1 --no-default-features --lib
cp target/wasm32-wasip1/release/LOSAT.wasm /mnt/c/Users/kawato/Documents/GitHub/gbdraw/gbdraw/web/wasm/losat/losat.wasm
```

Notes:
- The web app expects `losat.wasm` (lowercase filename).
- Newer builds export the direct `losat_web_*` API; the web app falls back to
  the WASI CLI path when those exports are absent.
- `tblastx` runs single-threaded in the browser build.
