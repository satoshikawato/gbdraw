Place the LOSAT WebAssembly binary here as `losat.wasm`.

Suggested build (from the LOSAT repo root):
```
cd /mnt/c/Users/kawato/Documents/GitHub/LOSAT/LOSAT
rustup target add wasm32-wasip1
cargo build --release --target wasm32-wasip1
cp target/wasm32-wasip1/release/LOSAT.wasm /mnt/c/Users/kawato/Documents/GitHub/gbdraw/gbdraw/web/wasm/losat/losat.wasm
```

Notes:
- The web app expects `losat.wasm` (lowercase filename).
- `tblastx` uses Rayon threads; ensure your WASI build supports threads and your hosting can set COOP/COEP if needed.
