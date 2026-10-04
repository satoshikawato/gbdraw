# INSTRUCTION PROMPTS: 共通の開始・終了手順

対象ソフトはgbdraw（Python bioinformatics engine＋Vueのsingle-page Web UI）。Issue #619はCircularの保存Width/Radiusがobject textになる不具合。実装目標は数値＋明示unit controlsで、typed scientific meaning・保存draft・Result・Historyを保つ。

すべてのセッションは [総合計画](../MASTER_PLAN.md) と自分のpromptを読み、S00→S01→S02→S03の順に実行する。同じbranchの並行writerは作らない。別branchへruntimeを分岐して後でまとめる手順は採らない。

## 1. 対象checkoutの引継ぎ

実装branch: `fix/issue-619-circular-track-measure-inputs`。remoteから **このbranchを取得して使う**。S00〜S03は直列で、対象branchの既存checkoutを引き継ぐ。セッションごとのclone、worktree、lock directory作成は必須にしない。

まず現在のbranch、変更、worktreeを確認する。

```bash
git status --short --branch
git worktree list
git fetch origin
git branch --show-current
git rev-parse --abbrev-ref --symbolic-full-name '@{upstream}'
```

期待するbranchは対象名、upstreamは`origin/fix/issue-619-circular-track-measure-inputs`。そのcheckoutを他sessionが使用中でないことと、対象範囲の未完了編集がないことを確認してから更新する。無関係なdirty/untracked filesは保持し、pullと競合する場合にreset/cleanで消さない。

```bash
git pull --ff-only origin fix/issue-619-circular-track-measure-inputs
git merge-base --is-ancestor 88028fd242d263f0fe86aaf9da57b8dc9eb082f6 HEAD
```

別作業との同時進行により作業ディレクトリの隔離が必要な場合だけ、既存repositoryのGit履歴を共有するworktreeを使う。対象branchのcheckoutが既にあれば、その作業終了後に引き継ぐ。同じbranchを`--force`で二重checkoutしない。checkoutがなければ対象branch用worktreeを一つ作成し、後続sessionも使う。main/devからbranchを作り直したり、別sessionのcheckoutをswitch/reset/cleanしたりしない。

依存環境は利用可能な既存環境を使い、wheelは対象sourceに対応するものを必要時だけprepareする。browser検証では、serverが対象checkoutを配信していることとport/profile/outputの所有者を確認する。他sessionのserver/browser/outputを操作しない。新たなserverや出力先は必要な検証に限って用意する。

## 2. 読むものと前提

- `AGENTS.md` / `CLAUDE.md` / `gbdraw/web/CLAUDE.md`。
- `docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md`。
- `docs/internal/PRODUCT_IMPACT_RATCHET.md` / `WEB_CHANGE_POLICY.md`。
- 本計画、3 Decision Packs、前セッションの `SESSION_RESULTS/Sxx.md`。

S01以降は、署名済みoutcomeが実際に`origin/dev`の既存static Product Contractへ入っていることを確認する。task branchのcandidate文書だけを根拠にしない。

```bash
git fetch origin refs/heads/dev:refs/remotes/origin/dev
git show origin/dev:docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md
```

origin/devを確実に更新するため、上の明示refspecを使う。実装前にtrusted-base authorityを含む最新devが必要になった場合は、成果を保持して`git merge origin/dev`で取り込む。public branchをrebaseしてforce pushしない。conflictはin-scopeのowner/authorityを確認して解決し、Product choiceを推論で変えない。署名・evidence・必要なauthority-only統合が欠ける場合は依存runtimeを変更せず、独立成果と未成立の境界をcommit/pushする。

## 3. 実装と検証

自分のpromptのfile ownershipを守る。後続sessionが必要とするpublic interfacesは結果文書に残す。余分なstrategy/registry、parallel request path、schema migration、同期watcherを作らない。mapped contracts/checkerを変更して同じruntimeの唯一の安全根拠にしない。

テストは変更に対応するfocused checksとrequired gatesを行う。同じcode/input/environment/条件の証拠は再利用する。failure・material change・未解決の境界がない限り、通った全testを毎session繰り返さない。長いtestは少なくとも30分の余地を与えて監視し、既存assertion timeoutを緩めない。

browserが必要なsessionは両環境を調べる。

```bash
command -v playwright
playwright --version
python -c "from playwright.sync_api import sync_playwright; print('python playwright ok')"
node -e "console.log(require.resolve('@playwright/test'))"
```

Node spec runnerがなければPythonで等価checkを行う。sandboxのChromium起動失敗は同じcheckを必要なescalationで試す。wheelは`python tools/prepare_browser_wheel.py`で作り、commitしない。cache-bustはdeployable bundleの準備時のみ。tests/reference_outputsは通常read-only。

## 4. 各セッション終了時の必須 commit/push

**各セッションは、所有範囲の成果と結果文書を対象branchにcommitし、同名remote branchへpushして終了する。** この作業branchへのcommit/pushは本実装ワークフローのauthorized scope。別targetのpush、merge、deploymentへ広げない。

`SESSION_RESULTS/Sxx.md`に以下を記録する。

- branch、起点SHA、対象outcome/authority、変更owner/path。
- 実行commands/exits、input/environment、C619受入結果、artifact location。
- 未実行checkと理由、必要な残作業、次sessionの明示前提。
- proposed English commit titleと短いsummary。

stageは自分の変更pathだけを明示する。`git add .`／`git add -A`を使わない。production/test/docs/generated diffを別々にレビューし、generated wheel・private fixture・無関係な変更を混ぜない。失敗checkを直すか依存boundaryを明記した独立文書のみを残し、壊れたruntimeを完了と呼ばない。

```bash
git status --short
git diff --check
# git add <your explicit owned paths, including SESSION_RESULTS/Sxx.md>
git diff --cached --stat
git diff --cached
git branch --show-current
git rev-parse --abbrev-ref --symbolic-full-name '@{upstream}'
# git commit -m '<English title for the completed session>'
git fetch origin
git merge-base --is-ancestor origin/fix/issue-619-circular-track-measure-inputs HEAD
git push origin HEAD:refs/heads/fix/issue-619-circular-track-measure-inputs
git rev-parse HEAD
git ls-remote origin refs/heads/fix/issue-619-circular-track-measure-inputs
```

remoteが自分のHEADのancestorでなければ、pushせず新しい内容・ownerを確認して取り込む。retry前に実際のremote SHAを確認し、成功済みのpushをAPIエラーだけで繰り返さない。force pushはしない。commit後のlocal HEADとremote SHA一致を確認し、最終handoffでbranch、SHA、tests、残るboundaryを報告する。

終了時は、自分が起動したserver/browser processだけを停止する。対象checkoutは次sessionへ引き継ぎ、検証証拠を保持する。

## Prompt index

- [S00 — authorityと表現境界](S00_AUTHORITY_AND_BOUNDARIES.md)
- [S01 — scalar ownerとeditor model](S01_SCALAR_AND_EDITOR_MODEL.md)
- [S02 — controlsとlifecycle](S02_CONTROLS_AND_LIFECYCLE.md)
- [S03 — browserとhandoff](S03_BROWSER_AND_HANDOFF.md)
