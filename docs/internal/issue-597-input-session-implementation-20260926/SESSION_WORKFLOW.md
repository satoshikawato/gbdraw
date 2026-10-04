# 全セッション共通作業規約

対象 repository は `satoshikawato/gbdraw`。対象 branch は `fix/issue-597-input-session-20260926`。
この branch は2026-09-26の最新 `origin/dev` から作成した。新sessionは最新 remote の同名branchを取得して続ける。
各promptはこの規約と [MASTER_PLAN.md](./MASTER_PLAN.md) を読むことを必須とする。

## 1. 独立 checkout に取得する

各sessionで次の手順を実行する。`SXX` は実行するsession番号に置き換える。
作業用変数は HOME/CODEX_HOME 等に流用しない。

```bash
ISSUE597_CHECKOUT=$(mktemp -d /tmp/gbdraw-issue597-SXX.XXXXXX)
git clone --single-branch --branch fix/issue-597-input-session-20260926 \
  https://github.com/satoshikawato/gbdraw.git "$ISSUE597_CHECKOUT"
cd "$ISSUE597_CHECKOUT"
git fetch origin
git switch fix/issue-597-input-session-20260926
git pull --ff-only origin fix/issue-597-input-session-20260926
git status --short --branch
git branch --show-current
git rev-parse --abbrev-ref --symbolic-full-name '@{upstream}'
```

期待branchは `fix/issue-597-input-session-20260926`、upstreamは `origin/fix/issue-597-input-session-20260926`。
shared checkout を switch/reset/clean しない。別sessionのworktree、node_modules、venv、サーバー、output filesを操作しない。
既存の自分専用checkoutを再開する場合はdirty diffを先に確認し、関係ない変更を残したまま pull/rebase しない。
originがsingle-branch cloneでも dev/main の明示取得は可能である。

```bash
git fetch origin dev:refs/remotes/origin/dev main:refs/remotes/origin/main --tags
git merge-base --is-ancestor d457b7189b137185a8dec800819a312c30b969fa HEAD
```

必要authorityがdevに入ったらruntime開始前にclean treeで取り込む。

```bash
git merge --no-edit origin/dev
```

このdev統合mergeは実装sessionのscope/resultsで報告し、同じimplementation branchへpushする。
無関係なdev変更と競合した場合はtargetedに解決し、devや別worktreeを変更しない。
最新devから別の実装branchを切り直して同名remoteへforce pushしない。
user指定のcontinuation branchを使うことが、このtaskの「fresh dev-derived work branch」方針の継続である。

## 2. 開始前の scope と readiness

1. `AGENTS.md`、`CLAUDE.md`、`gbdraw/web/CLAUDE.md`、本計画、担当prompt、前提sessionのresultsを読む。
2. `git status` と前提commitを確認。個人の作業開始HEAD、取得remote HEAD、dev HEADを記録。
3. scope、承認outcome、static/durable/privileged authority のbase merge SHAを確認。
4. 未完了の前提を隠さず、独立作業を完了した後で依存部分だけ止める。既存Product approvalを再要求しない。
5. repo外のtask専用output directoryと空きportを選ぶ。自分のserverだけ終了し、共通port/serverをkillしない。

セッションは順次実行する。一つのsessionを終えてpush・remote verificationを済ませてから次を開始する。
他担当者のcommitがremoteに増えていたら追加実装を始めず、latest stateと結果を再確認する。

## 3. 実装と verification

少ないowners/paths/compatibility branchesで契約を満たし、superseded codeを同じchangeで除く。
失敗をfallbackやgate緩和で隠さない。mapped contract、checker、authorityの変更は規範が指定した別deliveryにする。
unchanged code/input/environmentの証拠は再利用できる。変更/failure/未解決項目がなければ無目的にfull testsを繰返さない。

```bash
command -v playwright
playwright --version
python -c "from playwright.sync_api import sync_playwright; print('python playwright ok')"
node -e "console.log(require.resolve('@playwright/test'))"
```

Node specsには `@playwright/test` が必要。無ければPython Playwrightで同じfocused assertionを行う。
Chromiumのsandbox errorは同じcheckを適切なsandbox escalationで再実行する。
local dependency installationは自分のcheckout/environmentのみ。shared environmentの更新はしない。
wheelが必要なら `python tools/prepare_browser_wheel.py`、generated wheelはcommitしない。
portは例として `GBDRAW_WEB_TEST_PORT=...` をsession専用に指定する。
long testは30分以上の実行を許し、途中出力をmonitor。remote CI pollingは5分以上の間隔。

基本commands（担当promptの対象に絞る）:

```bash
node --test tests/web/session-file.test.mjs tests/web/session-request.test.mjs \
  tests/web/session-active-files.test.mjs tests/web/record-display-options.test.mjs \
  tests/web/record-metadata-inference.test.mjs
npx playwright test --config=playwright.functional.config.js --workers=1 \
  tests/web/circular-record-presentation.playwright.spec.js
npx playwright test --config=playwright.perf.config.js --workers=1 \
  tests/web/vibrio-session-save.performance.playwright.spec.js
pytest tests/ -v -m "not slow"
ruff check gbdraw/
```

production/tests/docs/generated diffsを別々にreviewする。
baseと同じchecker implementationによるlocal policy checkを実行し、trusted CIの結果と混同しない。
policy/checkerを変更しないimplementation sessionsは次を使える。

```bash
node tools/check-web-change-budget.mjs --base origin/dev
```

commit後は `--head HEAD` を指定して同じbaseとのcommitted rangeを確認できる。
required mapped contractの変更をruntimeと同じPRで自己証明に使わない。
`TestOutputComparison`はread-only。geometryを意図的変更した場合だけ別途reviewしてreferencesを更新する。

## 4. 結果文書と一commit

`results/SXX_RESULT.md` に以下を記入する。

- 完了/前提待ち/失敗の区別、scope、開始HEAD、source/fixtures/environmentのSHAやversion。
- 選んだengineering方式と理由、authority refs/base merge SHA、削除したold owner/path。
- 実行command、exit/result、artifact/log paths、未測定項目、残るdependency。
- Product checkpointsとacceptance IDの対応、production/tests/docs/generated review。
- 新compatibilityならnamespace別main/tag evidence、positive fixture、OE/PE/CBの完全なsets/実 arithmetic。
- 次sessionが読むべき具体的な結果と次のcommand。会話や外部の未保存メモを前提にしない。

巨大logsは既存CI artifact方式かrepo外のtask outputへ保存し、tracked resultsには再生成recipe/checksumsと必要なmetricsを残す。
private sequence/full rowsをlogsへ出さない。実測source SHAは実測時のHEADであり、結果文書を含むcommitの自己SHAを要求しない。
各promptのEnglish commit titleを使い、担当filesと結果だけを明示的にstageする。`git add .` は使わない。

```bash
git diff --check
git diff --cached --check
git diff --cached --stat
git branch --show-current
git rev-parse --abbrev-ref --symbolic-full-name '@{upstream}'
git commit -m "<担当promptのEnglish commit title>"
```

検証失敗を隠してsession完了にしない。前提待ちの候補/evidenceだけ完成したsessionも、その成果はcommit/pushし、runtime未着手を明記する。

## 5. 毎sessionのpushと確認

ユーザーはこの同名implementation branchへの各session commit/pushを指示している。再確認を求めない。
origin/dev/mainへ直接pushしない。`--force`/`--force-with-lease`は使わない。

```bash
git fetch origin
git merge-base --is-ancestor origin/fix/issue-597-input-session-20260926 HEAD
git push origin HEAD:refs/heads/fix/issue-597-input-session-20260926
git ls-remote --heads origin refs/heads/fix/issue-597-input-session-20260926
git rev-parse HEAD
git status --short --branch
```

remote SHAとlocal HEADを一致確認してからsession終了。remote statusを確かめず同じmutationを再試行しない。
non-fast-forwardでremoteだけ進んだ場合は実際のcommitとscopeを読み、相手の成果を保持して取り込む。
cleanな自分の未公開session commitだけならremoteの同名branchへrebase可能。実質競合はtargeted解決・関連checks再実行。
既にpush済みcommitをrewriteせず、必要ならmerge/追加修正で収束する。remote成功後のAPI errorを理由に新codeを足さない。
最終報告はbranch、commit、remote確認、checks、残るboundary、English summaryを含める。
PR作成/merge、別targetへのpush、deploy/tagはその対象の別authorizationを確認する。
