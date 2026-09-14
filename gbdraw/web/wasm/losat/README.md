LOSAT WebAssembly assets, updated 2026-09-14:

- `losat.wasm`: serial command artifact with the direct `losat_web_*` API.
- `losat-threaded.wasm`: threaded command artifact for cross-origin-isolated browsers.

`-num_threads N` includes the calling worker. For N > 1, the parent participates
in computation alongside N-1 child workers. Thus n8 uses eight total compute
threads; n1 creates no children. The caller's Rayon registration is cleared
after each search.

The threaded module exports `_start` and `wasi_thread_start`, imports shared
`env.memory` and `wasi.thread-spawn`, and requires COOP/COEP response headers.
The serial module retains the ten-argument `losat_web_run_pair` API and remains
the fallback when the browser cannot use shared memory.

## Build identity

Both artifacts were built with Rust 1.92.0 from the LOSAT working-tree snapshot
identified below. The total-thread correction was subsequently committed as
`1de348ad` (`Count the caller in total LOSAT thread budgets`). The source
manifest identifies the exact build inputs, including changes present before
that commit; the pre-build HEAD alone does not identify these binaries.

- Pre-build LOSAT HEAD: `a4fddb4a279d6ce71ffd94002713ebf48c244276`.
- Source manifest SHA-256: `b403985aa02c5abadf6706e1777bf66e59bb7035c9e108f41add9e31d7fb4334`.
- `losat.wasm`: 2,092,934 bytes; SHA-256 `3bf733eacbef32733b6d5afd962bdccc291b59cd8d6addb24e1ec64965baa9d3`.
- `losat-threaded.wasm`: 2,348,248 bytes; SHA-256 `8240f41500994c2881297bb94b3c7d463e158a079f69a3b5e135b230e39d7149`.

Build from the `LOSAT/` crate directory within that repository, with the
recorded source snapshot and toolchain:

```bash
cargo build --release --locked --offline --bin LOSAT --no-default-features --target wasm32-wasip1 --target-dir target/serial-command
cargo build --release --locked --offline --bin LOSAT --features wasm-threads --target wasm32-wasip1-threads --target-dir target/threaded-command
```

Copy `target/serial-command/wasm32-wasip1/release/LOSAT.wasm` to `losat.wasm`
and `target/threaded-command/wasm32-wasip1-threads/release/LOSAT.wasm` to
`losat-threaded.wasm` together with the matching argument adapters.

The LOSAT checkout retains `artifacts/gbdraw-wasm-total-threads-20260914/`,
including the source snapshot, `build-identity.json`, binary/ABI inspection,
and verification records. Its `docs/wasm_total_threads_20260914.md` describes
the scheduler correction and its verification limits.

## Consumer integration

The runtime and threaded worker use canonical NCBI-style CLI options
(`-query`, `-subject`, `-outfmt`, `-num_threads`). Circular and linear
orchestration supplies `-task`, `-query_gencode`, `-db_gencode`, `-max_hsps`,
and `-max_target_seqs` as applicable. Option values and search policies remain
owned by their existing builders. Cache identities include these arguments,
so entries using the old spellings are not reused for new searches.

This implements the existing total-thread and search contracts. The same
orchestration, runtime adapters, N-1 child preparation, AUTO selection, and
serial fallback remain the owners and paths before and after this update.
The old argument spellings are replaced without adding a compatibility path.
Rollback must restore both binaries and their matching JavaScript adapters.

## Verification

Eight focused Node test files passed, covering command argument builders,
circular/linear orchestration, source batching, cache identity, settings and
session schemas. The command-builder test also passed native LOSAT argument
acceptance for BLASTN, BLASTP and TBLASTX. The Web change policy gate passed
with a clear review result.

Chromium 149.0.7827.55 passed 23 checks against these packaged files: all five
BLASTN, megablast, TBLASTX, genetic-code-4 TBLASTX and BLASTP fixtures matched
the recorded native bytes across serial fallback and threaded-command dispatch,
including repeated runs, cancellation/recovery and BLASTP n8. External network
requests and page errors were zero.

The existing dispatcher limits non-BLASTP jobs to one thread even when more
are requested. A separate three-run runtime-status probe confirmed BLASTN
requested n2/effective n1, and BLASTP effective n2 and n8. This update preserves
that dispatch policy; it makes no new speed or release-certification claim.
Results are retained as `browser-gbdraw-branch/results.json` and
`browser-gbdraw-budget-probe/results.json` in the LOSAT evidence directory above.
