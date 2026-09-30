# Web GUI 監査の修正 — 実装リファレンス

Status: 実装の入口。判断はすべて確定している（2026-09-30）。
読者: 監査で見つかった不具合をすべて直すセッションと、そのサブエージェント。
作成時の基準: `origin/dev` `4c89bab1`。行番号はこの commit のものなので、作業時は最新の dev で symbol から探し直す。

## 0. 目標と完了の定義

**目標:** 監査の 65 件と、設計段階で見つかった N-01〜N-20 を dev ですべて直す。あわせて、同じクラスの再発を防ぐガード（G-A〜G-J）を入れる。

**完了の定義:**

1. 第 5 節の PR（P00〜P20）が dev にマージされ、最後の dev SHA で `dev-staging-gate` と Gallery publication が緑になっている。
2. 各 ID について、次のどちらかが成り立っている。
   - assert を持つテストが入っていて、修正前の dev で失敗し、修正後に通る。known-defect の印は外してある。
   - 修正しないことにした理由が SESSION_LOG に記録されている（ER の結果、仮説が否定された、など）。
3. 第 2 節の receipt 30 件が Product Contract に登録され（P01）、各 PR の本文に、対応する PD-OI の番号が書かれている。
4. 挙動が変わった箇所の docs（`docs/REFERENCE/`、`docs/SESSION_COMPATIBILITY.md` など）と CHANGELOG が更新されている。
5. 出力が変わった参照出力と Gallery は、owner tool で作り直し、目視で確認したことが記録されている。
6. 最終報告（第 7 節）がある。

main（リリース版）への backport はしない。main には、次の dev → main の昇格で入る（W-1）。

## 1. 読む順序

1. リポジトリの `AGENTS.md`、`CLAUDE.md`、`gbdraw/web/CLAUDE.md`
2. この文書
3. [02_DECISION_PACK.md](02_DECISION_PACK.md): 承認済みの receipt の全文。この文書と矛盾する記述があれば 02 に従う
4. [01_REMEDIATION_PROPOSAL.md](01_REMEDIATION_PROPOSAL.md)
   - 第 3 節: 監査の訂正と N-01〜N-20
   - 第 4 節: 設計規則 R1〜R12
   - 第 7 節: バグ別の修正案
   - 第 8.2 節: ガード
5. 担当する PR の付録 [remediation/W*.md](remediation/): 行番号、試作、実験の記録
6. 監査の [README.md](README.md) と、baseline の `evidence/<area>/FINDINGS.full.md`: 再現の手順
7. 必要に応じて規約: `docs/internal/WEB_CHANGE_POLICY.md`、`docs/internal/PRODUCT_IMPACT_RATCHET.md`、`docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md`、`docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md`

## 2. 確定事項

**Product の決定（receipt の全文は 02）**

| ID | 対象 | 決定 |
|---|---|---|
| D-01 | FE-06 | All から配列の内容と `/translation` を外す |
| D-02 | TR-10(b) | help tip をキーボードとタップで開ける disclosure button にする |
| D-03 | TR-02 | 両モードで、論理 series が最初の source を得たときだけ managed 行を足す。series が source を失ったら除く |
| D-04 | IN-01 | Circular の定義系の設定を Applies on Generate にし、即時の定義再生成を削除する |
| D-05 | GE-02 | 全体の stroke 設定を Applies on Generate にする |
| D-06 | PV-04 | feature のある凡例項目の衝突にも、既存のダイアログを使う |
| D-07 | FE-02、PV-09 | batch の Result を表示するときに、正本の編集を投影する |
| D-08 | PV-03 | 凡例の順序を Generate 後も引き継ぐ |
| D-09 | PV-07 | canvas padding を Generate 後も引き継ぐ |
| D-10 | FE-01 | source を置き換えたときは、消えた target の label 編集だけを外す |
| D-11 | R1 | 即時編集 ≡ Generate を契約にする（OIC-027 を足す） |
| D-12 | IN-06 | record ごとに、その record 自身の推定値で定義を付ける |
| D-13 | IN-03 | GFF の行を持つ record だけを扱う（CLI と同じ） |
| D-14 | FE-09 | 常に stable hash を出す |
| D-15 | SE-10 | Reset で、Linear の record ごとの表示状態と alignment plan も初期化する |
| D-16 | FE-11 | 検索の Start/End の項目を削除する |
| D-17 | PV-05 | Web の PDF を px → pt で換算する（75%） |
| D-18 | Q-FRAME | すべての比較表を探索座標系（F）にし、planner が向きを投影する |
| D-19 | CO-04 | ファイル単位の database を維持し、要求していない自己検索をどのモードでも実行しない（PD-OI-018 revision 4） |
| D-20 | CO-06 | 矛盾はエラー、不明な ID は警告にする |
| D-21 | CO-08 | member 集合で付け替え、一致しないものは dormant として保存する |
| D-22 | CO-10 | source の座標を主にし、違うときだけ表の座標も出す |
| D-23 | CO-05 | 先頭 12 列が型どおりなら受け付ける |
| D-24 | PV-08 | 配置に失敗したときだけ、種名の行を折り返して配置し直す |
| D-25 | SE-05 | 旧形式の Session は、Generate してから Save する（案内を出す） |
| D-26 | Dinucleotide | ACGTU の 2 文字（大小を区別しない）。U は T として扱う |
| D-27 | 負の値 | フォントサイズ > 0 と stroke 幅 ≥ 0 だけを検証する |
| D-28 | GE-06 | Generate 中の Undo/Redo は busy として拒否する（draft の編集は PD-OI-051 のとおり許す） |
| D-29 | TR-03（Linear） | source のない series を有効な手動の行が参照していたら、row issue を出して Generate の前に止める |
| D-30 | PV-10 | 差を計測して JS を Python に合わせる。一致させられない場合だけ、Linear の即時の side 移動を Applies on Generate にする |

**現状を維持する項目（D-31〜D-40）:** PV-01 の delta、IN-04 の単一 → 単一の置き換え、BOM を Python で受け付けないこと、TR-09 の初回、SE-08(b)、SE-06 の読み取り専用、Generate 中の draft 編集、エラーの field へ移動するボタンを足さない、α 付きの hex を足さない、CLI の LOSAT database（Web との違いは docs に書く）。これらは変えない。

**作業の進め方（W）**

| ID | 決定 |
|---|---|
| W-1 | backport はしない |
| W-2 | 範囲は全件（N-01〜N-20 を含む） |
| W-3 | PR は workstream 単位にする。size の Review REQUIRED は受け入れる |
| W-4 | authority は最初に 1 本の PR で登録する |
| W-5 | dev へは、必須の check が通ったら自動でマージしてよい。main は対象外 |
| W-6 | 参照出力と Gallery は作り直してよい |
| W-7 | Workflow による並列化を使う |
| W-8 | ガードをすべて入れる |

## 3. 作業環境

- **文書の置き場所:** この文書群は main checkout（`/mnt/c/Users/genom/GitHub/gbdraw`）の `docs/internal/web-gui-audit-20260930/` に、未 commit のまま置いてある。P00 で dev に入れる。
  - 入れるもの: `README.md`、`01`〜`03`、`SESSION_LOG.md`、`remediation/`。
  - main checkout は別の作業ブランチで、未 commit の変更がある。その作業ツリーは変更しない。
- **証拠:** `/home/kawato/gbdraw-baselines/web-gui-audit-20260930/`（リポジトリの外）。
  - `evidence/<area>/`: FINDINGS、FINDINGS.full、fixture、出力。
  - `specs/audit-<area>/`: 監査の spec。
  - spec と script は scratchpad の絶対パス（`/tmp/claude-1000/...`）を参照しているので、使うときは baseline のパスに置き換える。
  - 取り込む fixture は小さな合成ファイルにして、`tests/fixtures/` に置く。
  - `design/`: 設計担当が作った probe、実験の script、試作。
    - `design/w4-history-prototype.patch` は、SE-01、SE-02〜04、GE-06、GE-09、SE-06 の試作を `4c89bab1` との差分にしたもの（7 files、316 行）。P09 と P10 の出発点に使える。
    - この patch は試作で、レビューを経ていない。規則（第 4 節）に照らして書き直す前提で使う。
- **ブランチ:** PR ごとに、最新の `origin/dev` から新しいブランチを作る（`git switch --no-track -c <branch> origin/dev`）。並列で進めるときは worktree を使い、`/mnt/c/Users/genom/GitHub/gbdraw/.worktrees/<name>` に置く。
- **worktree ごとの準備:**
  - `python tools/prepare_browser_wheel.py`
  - `ln -s /mnt/c/Users/genom/GitHub/gbdraw/node_modules node_modules`（`@playwright/test` 1.61）
  - Playwright はエージェントごとに別のポートで動かす（`GBDRAW_WEB_TEST_PORT=<port> npx playwright test <spec> --workers=1`）。
- **PR を出す前の手元の確認:**
  - `ruff check gbdraw/`
  - `pytest tests/ -m "not slow" -n auto`
  - `node --test tests/web/*.test.mjs`
  - 対象の Playwright spec
  - `node tools/check-web-change-budget.mjs --base origin/dev`
  - 参照出力を変えた場合は `pytest tests/test_output_comparison.py::TestOutputComparison -v`
- **PR の文面:**
  - `.agents/skills/write-clear-pull-request/SKILL.md` に従い、`.github/pull_request_template.md` を使う。
  - `node tools/check-pr-language.mjs --title "<title>" --body-file <path>` を 1 回実行する。
  - 本文の末尾には、指定された attribution の行を付ける。
- **CI:**
  - PR には `pr-gate` がかかり、dev へのマージ後には `dev-staging-gate`（functional Playwright を含む）と Gallery publication が走る。
  - 状態の確認は 5 分以上の間隔で行う。
- **自動マージ:** 必須の check がすべて通り、merge state が CLEAN になったら `gh pr merge <n> --merge --match-head-commit <sha>` でマージする。main への push とマージはしない。
- **参照出力:**
  - `pytest tests/test_output_comparison.py::TestGenerateReferences --update-reference-outputs -v` で作り直す。
  - `git diff -- tests/reference_outputs/` を確認し、比較テストを通す。
- **Gallery:**
  - owner tool（`tools/refresh_gallery_sessions.py`、`tools/prepare_interactive_gallery_assets.py`）で作り直す。
  - 手順は `.agents/skills/web-gallery-screenshot-maintenance/SKILL.md` に従う。
  - Gallery は generator が作るものなので、手で編集しない。
- **既知の flake:** 長い `page.evaluate` の「Execution context was destroyed」「Promise was collected」は、多くの場合 CDP による promise の回収で、画面遷移ではない。`evaluateWithRetainedPromise` を使うか、ボタンのクリックと待機で操作する。

## 4. すべての PR に適用する規則

### 4.1 手順

1. **失敗するテストを先に書く。** 監査の spec や script（または付録の試作）を assert を持つテストにし、最新の dev で失敗することを確かめる。P02 で known-defect の印を付けたテストがあれば、その印を外す。
2. **owner で直す。** 設計規則（01 の第 4 節）に従う。
   - R1: Result の書き手を限る。
   - R2: 表示の変化で編集を消さない。
   - R3: domain ごとに投影は 1 つ。
   - R4: fast path は完全一致か decline。
   - R5: 比較の座標系と再利用。
   - R6: 失敗の意味はエラーを出す側が持つ。
   - R7: 検証は Python の typed layer が持つ。
   - R8: Multi-Record Canvas は配置だけ。
   - R9: reader は 1 つ。
   - R10: watcher で修復しない。
   - R11: History の境界。
   - R12: Web↔CLI とアクセシビリティ。
3. **旧経路を同じ PR で消す。** 置き換えた parser、評価器、reader、helper、watcher、互換の分岐を残さない（DRY）。新しい module、framework、schema field は、2 つ以上の実経路を統合し、旧経路を消す場合だけ足す（YAGNI）。
4. **docs と CHANGELOG を更新する。**
5. 手元で確認し、PR を出す。マージした後に SESSION_LOG に記録する。

### 4.2 authority と分類

- P01 が入った後は、第 2 節の決定を受けた項目がすべて IMPLEMENT_EXISTING_AUTHORITY になる。PR の本文に、根拠の PD-OI の番号を書く。
- 次の ER の項目は、PR の中で証拠を取り、結果を SESSION_LOG に書く。
  - PV-10（D-30 の条件）
  - N-16（ラベルの reflow で draft が漏れるか）
  - FE-07（protein cache の昇格経路が古い結果を昇格しないか）
- **想定していない Product の判断が出たとき:** Owner の standing instruction（推奨案で進める）に従い、推奨案を Owner-delegated として採る。SESSION_LOG と最終報告に明記する。
  - ただし、承認済みの receipt と矛盾する場合や、receipt の範囲を超えて科学的出力を変える場合は、その項目を保留して報告する。
- **checker と authority の分離:** 次のものは別の PR にする。
  - checker（`tools/web-architecture-*.mjs`）の変更と、authority（`tools/web-change-policy.json`、規約文書）の変更。
  - guard（`tests/web/architecture-contracts.test.mjs` などの構造テスト）の変更と、runtime の変更。

### 4.3 書き換えが必要な既存テスト

現在の誤った挙動を正しいとしているテストがある。対応する修正と同じ PR で書き換える（OIPC-C08）。

| テスト | 関係する項目 |
|---|---|
| `tests/web/session-operation-consistency.test.mjs:28-33` | GE-06 |
| `tests/web/orthogroup-computation-cache.test.mjs:43-45` | CO-08 |
| `tests/web/linear-sources.test.mjs:190-238` | CO-04（新しい自己検索の規定） |
| `tests/web/annotations.test.mjs:392-402` | FE-05 |
| `tests/web/depth-track-session.playwright.spec.js:469-553` | TR-02 |
| `tests/web/track-slot-display.test.mjs:98-118` | TR-08 |
| `tests/web/history-generated-authority.playwright.spec.js:180-205` | GE-02 |
| `tests/web/match-sequences.test.mjs:393-432` | CO-10 |
| `tests/test_linear_multi_record_comparisons.py:35-37` の fixture | CO-06 |

`tests/web/composition-layout.test.mjs:796`（PV-01）は、D-31 で現状を維持するので変えない。

## 5. PR の計画

「主な対象」の詳しい内容は、01 の第 7 節と付録にある。ホットなファイルは次の 5 つ。

| 記号 | ファイル |
|---|---|
| RA | `gbdraw/web/js/app/run-analysis.js` |
| AS | `gbdraw/web/js/app/app-setup.js` |
| CF | `gbdraw/web/js/services/config.js` |
| SR | `gbdraw/web/js/services/session-request.js` |
| IX | `gbdraw/web/index.html` |

| PR | 内容 | 主な対象 | 前提 | ホットなファイル | 備考 |
|---|---|---|---|---|---|
| P00 | 文書を dev に入れる | README、01〜03、SESSION_LOG、remediation/ | — | — | docs だけ |
| P01 | authority の登録 | 第 2 節の receipt 30 件（PD-OI-056 以降、PD-OI-018 revision 4）、OIC-027、N-15（規約文書と実装のずれ） | P00 | — | authority だけ。Review REQUIRED。登録の方法は 02 の末尾 |
| P02 | 検査基盤（tests だけ） | 全項目の known-defect の印、共有ベクタ（G-D）、helper での自動検査（G-G(3)）、無断の状態変化の検出（G-C）、2 record の batch fixture、N-14、`architecture-contracts.test.mjs` の厳密な件数を上限に変える（P11 の削除に必要） | P01 | — | runtime は変えない |
| P03 | Python: 比較表 | CO-05（D-23）、N-03、CO-06（D-20） | P01 | — | 科学的出力の Review。CHANGELOG |
| P04 | Python: Circular の配置と凡例 | PV-08、N-01、D-24 の折り返し、TR-01、N-04、N-06 | P01 | — | Web 既定の Circular 出力がすべて変わる。Gallery（Vnig など）を作り直す |
| P05 | Python: 翻訳 | FE-07、N-05 | P01 | — | LOSATP と orthogroup の出力が変わる |
| P06 | Python: 検証とエラーの意味 | X-01（Python 側。`diagnostic=`、adapter、slotIndex、depth、config/modify）、X-02（PR-A。D-26、N-13）、D-27（legend の ValueError を含む）、SE-08(a)、GE-03（Python 側。interval の helper と codec の移行） | P01 | — | P07 の前提 |
| P07 | Web: エラーと数値 | X-01（JS の PR-1、PR-2 と ratchet）、X-02（PR-B、PR-C）、N-12、IN-08、SE-05（D-25） | P06 | RA、AS、CF、SR、IX | |
| P08 | Web: record の検出と構成 | IN-02、N-02、IN-03（D-13）、IN-04、IN-05、IN-06（D-12）、FE-05、SE-10（D-15）、FE-09（D-14） | P07 | RA、AS、CF、SR | 共有ベクタ（P02）を使う |
| P09 | Web: History | GE-09、SE-01、N-19、N-20、SE-02、SE-03、SE-04、N-18、GE-06（D-28） | P07 | CF | Undo の step が増える変化を docs に書く |
| P10 | Web: Session、CLI、Run Info | SE-06、N-17、SE-07、GE-03（JS 側。Linear の unavailable）、TR-07 | P08 | RA、CF、SR | main で作った v42 の CLI sidecar を fixture にする |
| P11 | Web: 即時反映の退役 | IN-01（D-04）、GE-02（D-05）、N-16 の検証 | P02 | AS、IX | `architecture-change`。net は約 −750 |
| P12 | Web: editor の状態と batch | FE-01（D-10）、FE-03、FE-04、FE-10、FE-02 と PV-09（D-07） | P11 | RA、AS、CF | `architecture-change` |
| P13 | Web: 凡例と配置の継承 | PV-02、PV-03（D-08）、PV-04（D-06）、PV-07（D-09）、PV-10（D-30）、GE-07、PV-11、PV-12、PV-01（D-31。挙動は変えない。D-09 の padding の継承で、はみ出しを直せることだけを確かめる） | P12 | AS、IX | |
| P14 | Web: 比較 | CO-01、CO-02、CO-03(a)、CO-04（D-19）、N-09、N-10、CO-08（D-21） | P08 | RA、SR | 科学的出力の Review |
| P15 | 比較の座標系（Q-FRAME） | CO-07（D-18）、CO-10（D-22）、N-07、N-08 | P03、P14 | RA、SR | Python の planner、JS の変換の削除、main の Session の reader での変換（fixture 付き）。CLI の release note。Gallery |
| P16 | Web: トラック | TR-02（D-03）、TR-03（D-29）、TR-06、TR-08、TR-09 | P07 | AS | `architecture-change`（watcher → 明示的な遷移） |
| P17 | 検索、出力、色 | FE-06（D-01）、FE-08、FE-11（D-16）と N-11、PV-05（D-17）、PV-06、FE-12 | P01 | — | Interactive SVG の runtime が変わるので Gallery を作り直す |
| P18 | アクセシビリティと画面 | TR-10(a)、TR-10(b)（D-02。パネル単位に分けてよい）、TR-11、TR-12 と link checker | P01 | IX | capture の撮り直し |
| P19 | ガードの仕上げと規約 | 残りのガード（G-A、G-B、G-E、G-F、G-H、G-I、G-J(2)〜(4)）、`gbdraw/web/CLAUDE.md` への R1〜R12 の不変条件、`tools/audit/` の sweep と README、昇格 PR のチェック項目 | P03〜P18 | — | 修正の PR に入れたガードは重複させない |
| P20 | Result への書き込みの allowlist | G-J(1): (a) authority の PR で `tools/web-change-policy.json` に新しい capability を足す → (b) checker の PR で detector を分ける → (c) allowlist を最終形まで狭める | P19 | — | 3 本に分ける。checker が key の先行追加を許すかを先に確かめる |

- **本数:** 21 本の目安（P20 を 3 本と数えれば 23 本）。W-3 の「workstream 単位」に従い、同じ wave の小さな PR はまとめてよい。
- **IN-01 の規模:** P11 の IN-01 は、01 の第 7 節（net 約 −600）と GE-02（net 約 −140）を合わせたもの。

## 6. 並列化と順序（W-7）

| Wave | PR | 進め方 |
|---|---|---|
| Wave 0 | P00 → P01 → P02 | 直列 |
| Wave 1 | Python の P03、P04、P05、P06。Web のうちホットなファイルを触らない P17、P18 | 並列 |
| Wave 2 | P07 → P08 → P10、P07 → P09、P07 → P16、P02 → P11 → P12 → P13、P08 → P14 → P15 | 依存に従って進める |
| Wave 3 | P19 → P20 | 直列 |

Wave 2 の規則:

- 同じホットなファイルを触る PR は、同時に開かない。
- 先にマージされた PR の上に rebase してから次を出す。

マージと同期の規則:

- dev の staging の 1 周期でマージする runtime の PR は 2〜3 本までにする。dev の先端で `dev-staging-gate` が緑になってから、同じクラスの次の PR をマージする。
- staging が赤になったら、マージを止める。その周期にマージした範囲の中で原因を切り分ける。
- Gallery を作り直す PR（P04、P15、P17、必要なら P14）は、互いに直列にする。
- index.html の衝突は、先にマージした側を正として rebase で解消する。

Workflow の使い方:

- worktree と port をエージェントごとに分ける。
- 各エージェントには、担当の PR の行と、01 と付録の該当項目だけを渡す。
- orchestrator は PR の順序、rebase、マージ、SESSION_LOG の更新を持つ。

## 7. 記録と最終報告

- **SESSION_LOG:** [SESSION_LOG.md](SESSION_LOG.md) に、PR をマージするたびに 1 行を追記する。書く内容は、日時、PR 番号、merge SHA、直した ID、テスト、生成物の更新、残った問題。
- **最終報告:** 次の内容を含める。
  - 全 ID の状態の表（修正済み／理由付きの見送り／ER の結果）と、各 ID の PR。
  - 最後の dev SHA と、その SHA の CI の結果。
  - Owner-delegated で決めた事項。
  - 作り直した参照出力と Gallery の一覧。
  - main への昇格の前に確かめる事項（PV-08 の出力の変化、CO-05、Q-FRAME の CLI の意味の変化、D-17 の PDF の大きさ）。
