# Session 01 Instruction Prompt — Baseline、Product preflight、architecture ledger

## Prompt

あなたはgbdraw Issue #563「feature popupからcircular recordをfeature基準で回転する」の
Session 01を担当する。過去の会話は存在しないものとして、この文書だけで開始する。

### 必須branch

`issue-563-feature-popup-record-rotation-20260922`を使用する。このbranchは2026-09-22に
最新`origin/dev`の`a9eaeadd105e0c26e46086626feaa695bdd33c94`から作成済みである。
別branchを作らず、`main`または`dev`へ直接commitしない。

最初に実行し、結果を記録する。

```bash
git branch --show-current
git status --short --branch
git rev-parse HEAD
git rev-parse --abbrev-ref --symbolic-full-name @{u}
```

upstreamは未設定が正しい。worktreeに無関係な差分があれば保持し、production、tests、docsに
分けて監査する。削除、reset、stash、rebaseをしない。

### 最初に全文を読む

1. `AGENTS.md`
2. `CLAUDE.md`
3. `gbdraw/web/CLAUDE.md`
4. `docs/internal/issue-563-feature-popup-record-rotation/IMPLEMENTATION_PLAN.md`
5. `docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md`
6. `docs/internal/PRODUCT_IMPACT_RATCHET.md`
7. `docs/internal/WEB_CHANGE_POLICY.md`
8. `docs/internal/PRODUCT_DECISION_PACKET_TEMPLATE.md`
9. GitHub Issue #563の最新本文

### Session goal

runtime実装前の事実、authority、owner/path/compatibility境界、acceptance evidenceを固定する。
このsessionでpopup機能を実装しない。推測に基づくschemaやUIを先行追加しない。

### 確認する現状

- `RecordPresentation.reverse_complement`、`RecordDisplayOptions.start_coordinate`、
  `RecordDisplayTransform`が既存表示変換を所有する。
- `record-display-options.js`がrecord表示draftを所有する。
- `session-request.js`がcanonical request schema 7の唯一ownerである。
- `run-analysis.js`と既存result admissionが生成経路を所有する。
- `services/history.js`がartifact replacement transactionを所有する。
- Session version 43、feature catalog schema 3である。
- 通常Generateでもdisplay start/reverse complement変更後のLOSAT executorはcache hit時に
  再実行されない。`LOSATP source jobs are reused after display start, reverse complement, and fresh Load`
  testを根拠にする。

### 作業

1. Issue #563の各要求を`AC-01`から`AC-20`へtraceするbaseline/preflight文書を、
   `docs/internal/issue-563-feature-popup-record-rotation/BASELINE_AND_PREFLIGHT.md`として作る。
2. Product Impact mandatory developer preflightを行う。
   - merged authorityがIssue記載の完全outcomeを既に選択しているなら
     `IMPLEMENT_EXISTING_AUTHORITY`と根拠を記す。
   - base-branch authority receiptが必要ならDecision/authority-only routeを記し、Session 02以降を
     `BLOCKED`とする。Product choice、rationale、retirement、riskを補完しない。
   - evidenceだけが不足するなら`EVIDENCE_REQUIRED`と不足証拠、取得方法を具体化する。
3. semantic owner、canonical path、privileged operatorをinventoryする。
   - canonical request
   - feature anchor calculation
   - record display effective state
   - candidate execution/admission
   - SVG admission
   - history
   - record coordinate transform
4. architecture ratchetの通常/例外判定を記録する。
   - owner/pathは既存owner内のprivate decompositionとして増加させない。
   - Session 44から見たreleased v42/catalog 3 readerはcompatibility pathになるため、namespace、
     positive fixture、removal condition、before/after CB setをpolicyどおり宣言する。計画時点で
     branch-onlyのv43はreader対象にせず、branch-owned artifactをcurrent writerへ書き換える。
   - request schema 7、Worker protocol、renderer pathは不変とする。
5. current source/testからfile-level implementation mapを作る。少なくとも次を含める。
   - `gbdraw/web_support/feature_metadata.py`
   - `gbdraw/web_support/feature_catalog.py`
   - `gbdraw/web/js/services/feature-catalog.js`
   - `gbdraw/web/js/app/record-display-options.js`
   - `gbdraw/web/js/services/session-request.js`
   - `gbdraw/web/js/app/run-analysis.js`
   - `gbdraw/web/js/services/history.js`
   - `gbdraw/web/index.html`
6. Issue acceptanceを満たすfixture inventoryを作る。
   plus/minus、unstranded、mixed、fuzzy、multipart、origin-spanning、duplicate record IDs、
   same-file multi-record、circular/linear mode、LOSAT comparison、depth/statistics trackを既存fixtureで
   再利用できるか確認し、足りない最小fixtureだけを列挙する。
7. 現状を守るcharacterization testを実行する。新規testを追加する場合、このsession終了時に
   passする既存挙動のtestだけにする。意図的なred testをbranchへ残さない。

### 最低限のverification

```bash
node --test tests/web/record-display-options.test.mjs
node --test tests/web/session-request.test.mjs
node --test tests/web/history.test.mjs tests/web/run-analysis-simple-path.test.mjs
npx playwright test tests/web/linear-multi-record.playwright.spec.js --grep "LOSATP source jobs are reused" --workers=1 --retries=0
node tests/web/architecture-contracts.test.mjs
git diff --check
```

browser依存が不足する場合はrepositoryのPlaywright guidanceに従い、Python Playwrightを含む代替で
同じ契約を確認する。sandbox failureは権限付き再実行を行う。

### 禁止事項

- popup、resolver、schema migrationのruntime実装
- IssueにないProduct outcomeの選択
- schema 3へ意味を黙って追加すること
- 「専用経路がないとLOSATが再実行される」という誤った前提
- 新しいrequest owner、render path、history engine

### 完了条件とhandoff

- `BASELINE_AND_PREFLIGHT.md`だけで新規参加者がauthority、owner、files、fixtures、risksを復元できる。
- Product classificationが一つに決まり、次sessionが進行可能か明記される。
- 実行commandと結果、未実行理由、working tree diffを報告する。
- 実装sessionでcommit/pushするかは、そのsessionを呼び出したuserの明示指示に従う。このprompt
  自体はpushまたはPR作成を許可しない。
