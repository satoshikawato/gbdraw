# INSTRUCTION PROMPT S00 — 承認済み製品結果の恒久記録

あなたは gbdraw の [Issue #599](https://github.com/satoshikawato/gbdraw/issues/599) の authority-only セッション担当者である。Web の配置リセット、Layout edit の発見困難、検索・toolbar 配置問題の3件について、承認済みの完全な製品結果を既存の恒久契約へ記録する。runtime は変更しない。

## ブランチ取得と必読

前セッションの公開完了と、同じブランチの writer が自分だけであることを確認する。共有作業ツリーを切り替えず、必ず実装ブランチを専用 clone へ取得する。

```bash
task_dir=$(mktemp -d /tmp/gbdraw-issue599-s00.XXXXXX)
git clone --branch fix/issue-599-preview-layout-20260926 \
  https://github.com/satoshikawato/gbdraw.git "$task_dir/repo"
cd "$task_dir/repo"
git fetch origin
git pull --ff-only origin fix/issue-599-preview-layout-20260926
test "$(git branch --show-current)" = "fix/issue-599-preview-layout-20260926"
test "$(git rev-parse --abbrev-ref '@{upstream}')" = "origin/fix/issue-599-preview-layout-20260926"
git status --short
```

入力 HEAD と origin/dev SHA を記録する。clean tree を確認し、repository `AGENTS.md`、`CLAUDE.md`、`gbdraw/web/CLAUDE.md`、[総合計画](MASTER_PLAN.md)、[承認一覧](00_APPROVED_PRODUCT_DECISIONS.md)、3つの個別 Pack、[進捗](SESSION_STATUS.md) を読む。
総合計画が参照する Product/Architecture/Web 規約を最新 dev で確認する。

## 選択済みの製品結果

- [Pack 01](DECISION_01_COMPOSITION_CONTINUITY.md): legend/title/Linear scale の差分を同じ図の新自動配置へ1回継承。非ゼロの対応不能は旧 Result/request/History を保持して公開前停止。
- [Pack 02](DECISION_02_LAYOUT_AFFORDANCE.md): Layout edit は明示切替。OFF は pan と有効化説明、ON は装飾 drag。説明を keyboard/touch からも読める。
- [Pack 03](DECISION_03_PREVIEW_CHROME.md): 検索上部行・toolbar下部行・canvas/editor workspace。検索自由 drag を退役。機能と到達性を保持。

Product Decision Owner `satoshikawato` は2026-09-26に各 Pack の完全な推奨 A を承認済み。rationale、must preserve、may retire、accepted residual risk、owner/date は Pack の承認本文そのものを使用する。再署名、未承認の補完、別 outcome の選択は行わない。

## 所有範囲と作業

所有する変更は、既存の適切な authority 記録、必要な receipt provenance、`results/S00.md` と進捗である。計画作成基準では unmapped concerns の正本は `docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md`、revision16、PD-OI-001–035。新 executable JSON registry、仮の BD 番号、checker、workflow、runtime、テスト期待値変更は追加しない。

1. fetched origin/dev の map、BD、OIPC、runtime contract、supported compatibility を調べ、各 concern の既存規定と差分を記録する。作業ツリーの他 issue の未追跡提案を base authority とみなさない。
2. 3件が unmapped のままなら OIPC に関心ごとの独立した accepted record を追加する。最新 base の番号と revision に従い採番し、他 concern を変更しない。同一 outcome が既に active なら再利用する。mapped へ変わっていれば既存 map/BD 記録先を使用し、第二 authority を作らない。
3. 各 record を完全な承認 outcome と receipt のフィールドへ照合し、総合計画 C/A/P/R の受入条件へ結び付ける。独立した preservation contributions がすべて成立することを AND-of-OR で確認する。
4. 既存 outcome と実質的に衝突する場合は該当収束だけ止め、対立する本文・影響・必要な人間判断を具体的に提示する。無言の上書き、risk の推定、承認範囲の拡張はしない。必要な supersession は規約どおり別 authority-only 変更に限定する。
5. 正式記録先以外に同じ active 規範を増やさない。計画の Pack は承認証跡で、runtime はそれを実行設定として読まない。

実装アーキテクチャやファイル配置を Product outcome として固定しない。単一の契約 owner、3つの concern の分離、最小の記録差分で SOLID/KISS/DRY/YAGNI を守る。

## 検証と dev merge 境界

各承認本文と JSON 表現の内容一致、owner/date、完全な preservation/retirement/risk、参照リンク、authority-only diff を確認する。既存の authority/schema validation と trusted-base policy の該当チェックを使う。runtime が不変なら不要な renderer/browser 再検査はしない。

依存 runtime は恒久 authority の dev merge 前に開始できない。S00 は concrete な差分を commit/push し、`results/S00.md` に authority commit、必要な PR 対象 `dev`、merge後確認項目を記録する。PR 作成・merge が別途明示承認されていなければ、レビュー可能な差分を提示した後その境界だけ承認待ちとして残す。PR タイトル・本文を作成する際は `write-clear-pull-request` skill と repository validator を適用する。

恒久記録が dev に存在することを確認できるまで、進捗を「authority push済み / dev merge待ち」とし、S01 の開始を許可しない。S01 は origin/dev の実状態から記録と receipt を確認して継続する。

## 完了・公開

`results/S00.md` と SESSION_STATUS を更新する。authority と docs の diff を別々に確認し、対象パスを列挙して stage する。他の担当者の変更を取り込んで消さない。branch/upstream を再確認し、fetch 後の remote branch が予期せず進んでいれば担当者と差分を確認する。必要な統合を行い、force push は使わない。

`git diff --cached --check` 後、必ず commit/push する。

```bash
git commit -m "docs: record approved preview layout decisions for issue 599"
git push origin HEAD:refs/heads/fix/issue-599-preview-layout-20260926
git rev-parse HEAD
git ls-remote --heads origin refs/heads/fix/issue-599-preview-layout-20260926
git status --short
```

local/remote SHA の一致を確認する。remote mutation のエラー時は実際の remote state を調べて再試行する。main/dev へ直接 push しない。終了報告に公開 SHA、検証結果、恒久 record ID、dev merge の状態、S01 の開始条件、英語 commit title と short summary を含める。
