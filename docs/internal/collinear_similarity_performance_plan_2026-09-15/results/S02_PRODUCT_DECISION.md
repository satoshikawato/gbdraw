# S02 follow-up — PATH-B Product decision receipt

**Product selection is complete. Authority serialization is complete; review
and integration into `origin/dev` remain pending.** Do not request the same
choice again. This handoff records the receipt and preparation work; it is not
another Product authority or an authorization to begin S06 runtime.

## 1. Exact human receipt

Received after the S02 investigation commit
`c46e55d14b6b8fae9e6778a7f04183a4d321185f` and explanation of the two options:

```text
PRODUCT_DECISION
Concern: protein-comparison.path-representation
Scenario revision: 1
Choice: PATH-B / lossless-graph
Rationale: 解析情報を維持しながら、通常の生成・保存での経路展開コストを減らしたい。
Must preserve: 図と解析情報、全経路の内容・順序・ID・shared情報、対応済み旧ファイルの読み込み、明示的な旧tuple形式での全経路取得。
May retire: 通常のAPI戻り値と保存形式が、常に全経路配列を含む仕様。
Accepted residual risk: 旧API依存コードの修正、新形式を旧バージョンで読めないこと、明示的な全量取得には大きな時間・メモリが必要になり得ること。
Owner: satoshikawato
Decision date: 2026-09-15
```

All required fields are present. The owner is in the base maintainer allowlist.
No rationale, preservation, retirement, risk, owner or date was filled in by
Codex. PATH-A is not selected. The original [Decision Pack](PATH_DECISION_PACK.md)
is retained as the pre-decision comparison; its historical intake wording and
the verification manifest's null receipt are not the current selection state.

## 2. Authority route and generated representation

Fetch confirmed `origin/dev` remains
`9a4f7e29ab1b99676bb783321f8a9f40f969d04c`. There is no new competing path
decision or occupied `PD-OI-023` on that base. This unmapped concern uses the
existing static OIPC route in `PRODUCT_IMPACT_RATCHET.md`.

| Item | Prepared state |
|---|---|
| Authority worktree | `/tmp/gbdraw-collinear-path-authority` |
| Branch / upstream | `perf/collinear-path-authority-20260915` / none |
| Authority-only commit | `77c8335925f872b83381b67ee04a29014205f776` |
| Sole changed file | `docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md` |
| Contract / record | OIPC revision **6**, **PD-OI-023**, concern revision **1** |
| Human outcome | PATH-B / lossless-graph, complete receipt accepted |
| Serialization | Exact supplied fields as JSON inside PD-OI-023; ready for review |
| Representation review | Pending; no review approval inferred from the receipt |
| Base integration | Not integrated; base remains OIPC revision 5 |
| Runtime / readers / writers | Unchanged |

The JSON is within the existing OIPC authority document. No new machine store,
evaluator, BD registry record, scenario revision, schema allocation or acceptance
waiver was created. Implementation names and version numbers in the
[S02 design](S02_PATH_CONTRACT.md) remain engineering proposals.

Inspect the exact generated representation independently of worktree paths:

```bash
git show 77c8335925f872b83381b67ee04a29014205f776:docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md
git diff 9a4f7e29ab1b99676bb783321f8a9f40f969d04c 77c8335925f872b83381b67ee04a29014205f776 --stat
```

The authority branch starts directly from fetched `origin/dev`. It does not
cherry-pick the S01/S02 evidence commits: the trusted checker requires the OIPC
file to be the only changed path. The existing S02 branch contains this
non-normative handoff update separately. The shared dirty checkout is preserved.

## 3. Focused verification and limits

- Parsed the generated JSON and compared every field with the supplied receipt.
- Compared the previous Product records, cross-surface clauses and acceptance
  catalog byte-for-byte; all are unchanged. Only revision metadata and the new
  PD-OI-023 record differ.
- Confirmed the authority diff contains exactly one file and no production,
  test, checker, schema or generated figure changes; `git diff --check` passes.
- Ran the unchanged trusted-base checker in the authority worktree:

```bash
node tools/check-web-change-budget.mjs --base origin/dev
```

Result: **Gate PASS / Review REQUIRED**, no blocking violations. Review is
required for the Product authority change; it is not an additional Product
choice or an architecture exception approval. The checker validates separation,
not the human meaning of the prose/JSON. The receipt comparison above checks
serialization fidelity; a reviewer still confirms it.

The first sandbox run failed with `spawnSync git EPERM`; the identical command
succeeded with the required sandbox permission. Both logs are retained:
[successful gate](data/s02-decision-authority-gate.log),
[sandbox failure](data/s02-decision-authority-sandbox.log).
The gate ran on the exact authority file later committed as `77c83359`.
There are no runtime or test edits, so the earlier 464 Python / 168 Node results
are retained as historical S02 evidence, not reported as newly rerun tests.

## 4. Remaining work and commit handoff

1. Review the generated PD-OI-023 representation against the receipt above.
2. Integrate the isolated authority change through the repository's normal
   authority-only route when publishing is authorized. No push, PR, merge, tag
   or deployment was requested or performed here.
3. S06 must find this outcome on its actual base and have completed S05 before
   starting dependent runtime. Retain the small equivalence oracles, published
   compatibility, large-integer transfer, explicit full-output Ω(L), and the
   separate S06 runtime/browser/architecture acceptance gates. This receipt does
   not waive those gates or authorize S03–S08 implementation in this session.

**Authority commit title:** `Record the lossless protein path Product decision`

**English summary:** Serialize the complete PATH-B receipt as PD-OI-023 in OIPC
revision 6, preserving the supplied terms and keeping runtime unchanged.

**S02 handoff commit title:** `Update S02 handoff after the PATH-B decision`

**English summary:** Record the received Product choice, isolated authority
commit, validation evidence and pending review and base integration.

The authority commit is local. The handoff update is a separate local commit on
the existing S02 branch. Revert each local commit independently if correcting
its transcription; changing the accepted outcome itself requires an explicit
new human decision under the existing lifecycle.
