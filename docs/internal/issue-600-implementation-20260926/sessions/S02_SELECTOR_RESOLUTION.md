# INSTRUCTION PROMPT — S02: 未一致注釈の隔離と warning 伝搬

あなたは gbdraw issue #600 の **S02** 担当です。この prompt の担当範囲を実装・検証し、
指定ブランチへ commit/push するところまで完了してください。他の session の作業は開始しない。
会話ログや未公開の提案を前提にせず、下記の committed files を入力に使ってください。

## Checkout と前提

他セッションの shared checkout を切り替えず、独立 clone に実装ブランチを取得する。
remote writer はこの session だけであることを確認し、前 session の pushed report があれば読む。

```bash
SESSION_DIR="$(mktemp -d /tmp/gbdraw-issue600-s02-XXXXXX)"
git clone --no-checkout https://github.com/satoshikawato/gbdraw.git "$SESSION_DIR"
cd "$SESSION_DIR"
git fetch origin
git switch --no-track -c fix/issue-600-annotations-styles-20260926 origin/fix/issue-600-annotations-styles-20260926
git status --short --branch
git rev-parse HEAD
```

この clone の docs を読み、次の内容を開始前に確認する。

- [総合計画](../MASTER_PLAN.md): scope、architecture、session 順、Git 運用、acceptance。
- [承認済み Product Decisions](../APPROVED_PRODUCT_DECISIONS.md): 完全な採用本文と署名者。
- `CLAUDE.md`、`gbdraw/web/CLAUDE.md`、`docs/internal/PRODUCT_IMPACT_RATCHET.md`、
  `ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md`、`WEB_CHANGE_POLICY.md`。
- 当該 clone にある nearest `AGENTS.md` があれば、それも読む。

ブランチは最新 dev から作成された `fix/issue-600-annotations-styles-20260926` を継続利用する。
新たに dev から作り直して前 session の成果を捨てない。承認済み４ outcome は選び直さない。
runtime に新しい Product outcome が必要になった場合は、その範囲だけ evidence/Decision
を用意し、選択を推測しない。無関係な編集・process・server・tmp files に触れない。

## 開始条件・所有範囲

S00 authority が `origin/dev` に存在し、S01 の commit/report が当該 remote branch にある。
担当は Python annotation resolver と、同じ解決結果を API/CLI/Worker/Result に届ける境界。
総合計画 §5 と **SEL-01～05** が仕様と acceptance owner。

主な files: `gbdraw/annotations/resolve.py`、`models.py`、`planning.py`、
`gbdraw/api/request_render.py`、`api/diagram.py` の既存 builder 入口、必要な assembler 引数、
`gbdraw/web_support/request_render.py`、既存 Result admission/通知と対応 tests。
Result の owner を増やさず、別 JS matching/exception catcher を作らない。

## 実装

1. record binding/syntax validation と feature matching の missing facts を区別する。
   binding 成功済み row に missing selector が１つでもあれば geometry を作らず row skip。
   partial union を作らない。完全一致時の envelope/segments/circular_path は維持。
2. `feature_selector_unmatched` の code と set/row/record 識別・欠落件数で１ row１ warning。
   不要な qualifier 値や genome content の console dump をしない。empty_span と区別する。
   record 欠落/曖昧/index 範囲外、省略した multi-record binding は fatal のまま。
3. materialized records/transforms に対する１回の解決 bundle を drawing/reporting が共有。
   native/typed、Circular single/grid/batch、Linear の既存 build/preparation を収束する。
   grid の旧事前解決、assembler の同義再解決を置換し、batch の global record index を保つ。
4. `PreparedDiagramRequest` → `RequestRenderResult` → Web metadata の explicit additive
   field で運ぶ。logger 捕捉、global collector、Drawing side-channel を使わない。
   API は warning を取得でき、CLI は当該 render の warning を１回表示する。
5. warning は admitted successful Result に帰属させ、failure/cancel/stale 候補から混入
   させない。Save/Load/History の既存 artifact lifecycle に従い、canonical row は削除しない。
   全 row skipped でも genome 図と skip notification を返す。request に含まれた slots/gaps
   は維持し、resolved mark から native が新規 auto slot を作るときだけ empty slot を作らない。
6. fixed-key metadata をこっそり緩めない。既存 namespace に必要な additive change を行い、
   persisted format/migration を追加する必要があるなら、その scope の evidence/Decision
   を用意する。未承認の reader を便宜的に増やさない。

## Focused verification

```bash
pytest tests/test_annotations.py tests/test_annotation_planning.py tests/test_circular_annotation_tracks.py tests/test_linear_annotation_tracks.py -q
node --test tests/web/annotations.test.mjs tests/web/session-request.test.mjs
```

matched/missed 混在、all missing、partial missing、wrong/duplicate record、crop/reverse/
rotation、single/grid/batch/Linear、empty/explicit slot、coordinate policy を検証。
Worker/client/admission の warning lifecycle、cancel/stale/fail、Undo/Redo、Session、
actual Generate と download を targeted browser check で確認する。最終 drawing のために
同じ生物入力を再 parse/resolve していないことを instrumentation または call-count で示す。

## 完了と publication

1. 担当 acceptance の実装・focused checks・必要な修正を完了し、scope の旧経路を除去する。
2. `docs/internal/issue-600-implementation-20260926/SESSION_RESULTS/S02.md` に start
   work SHA、authority base SHA、scope、owner/path before→after、実行 command/exit/result、
   artifacts、未実行範囲、remaining boundary を記す。未完了を completion と呼ばない。
3. production/test/docs/generated diff を別々に review し、`git diff --check` を実行する。
   command、事実、結論を区別し、過去 session の証拠は変更のない条件で再利用する。
4. `git fetch origin` と remote work SHA を確認。他 writer が進めていたら remote の実状態を
   調べてから統合方針を決める。force-push で他人の commits を消さない。
5. 現 branch と upstream を確認し、担当 path だけを `git add -- <明示 path>` で stage。
   session を１ commit とし、英語 title と短い summary を使う。**commit と push を実施して
   から session を終了する。**

```bash
git branch --show-current
git for-each-ref --format='%(upstream)' refs/heads/fix/issue-600-annotations-styles-20260926
git diff --cached --stat
git commit -m "fix: isolate unmatched annotation selectors"
git push --set-upstream origin HEAD:refs/heads/fix/issue-600-annotations-styles-20260926
git rev-parse HEAD
git ls-remote --heads origin refs/heads/fix/issue-600-annotations-styles-20260926
```

branch は指定名、upstream は未設定か同名 remote に限定。local HEAD と remote SHA の
一致を確認し、最終 handoff に pushed SHA、summary、checks、remaining boundary を記す。
push error が出たら先に実際の remote state を調べ、成功済み mutation を再実行しない。
main/dev への直接 push、force-push、他 branch への push、未承認の PR publication/merge/
deploy/tag は行わない。次 session の実装には着手しない。
