# INSTRUCTION PROMPT S03 — 検索と toolbar を専用行へ配置

あなたは gbdraw の [Issue #599](https://github.com/satoshikawato/gbdraw/issues/599) BUG-18 の実装担当者である。現在の検索は自由 drag と固定360pxの drawer退避を持ち、CSS の transform override と競合して Preview 外へはみ出す。検索を上部専用行、toolbar を下部専用行へ置き、canvas/editor workspace と構造で分離する。

## 取得・開始条件

S02 の commit/push 完了と、自分が同じ remote branch の唯一の writer であることを確認する。他セッションの共有作業ツリーを変更せず、必ず指定実装ブランチを専用 clone へ取得する。

```bash
task_dir=$(mktemp -d /tmp/gbdraw-issue599-s03.XXXXXX)
git clone --branch fix/issue-599-preview-layout-20260926 \
  https://github.com/satoshikawato/gbdraw.git "$task_dir/repo"
cd "$task_dir/repo"
git fetch origin
git pull --ff-only origin fix/issue-599-preview-layout-20260926
test "$(git branch --show-current)" = "fix/issue-599-preview-layout-20260926"
test "$(git rev-parse --abbrev-ref '@{upstream}')" = "origin/fix/issue-599-preview-layout-20260926"
git status --short
```

repository AGENTS/CLAUDE、Web CLAUDE、[総合計画](MASTER_PLAN.md)、[Pack 03](DECISION_03_PREVIEW_CHROME.md)、[承認一覧](00_APPROVED_PRODUCT_DECISIONS.md)、[進捗](SESSION_STATUS.md)、S00–S02 の結果を読む。最新 origin/dev の authority を確認し、入力 SHA・authority base/record を記録する。検索自由 drag の退役は承認済みだが、candidate authority で runtime を自己承認しない。

## 承認結果と所有範囲

選択 `A / DOCKED-SEARCH-AND-CONTROLS` は2026-09-26に `satoshikawato` が完全な本文で承認済み。検索 drag は全幅で退役。query/field/regex/active match/focus、Prev/Next/Open/Enter、toolbar、drawer tab/Close/Escape、pan/zoom、同じ canvas/editor、Session/Export/History は維持する。

所有範囲は `gbdraw/web/index.html` の Preview DOM/CSS、`js/app/feature-search/preview-actions.js`、`app-setup.js` の撤去する検索配置 exports、既存 feature-search/navigation/responsive tests。drawer の visibility/tab transitions は `app/right-drawer.js` の単一 owner を維持する。S01 の delta provider と S02 の hint/mode を消さず、新 DOM でも再bindを確認する。

## 実装タスク

1. 既存 search DOM を上部 auto row、toolbar DOM を下部 auto row、同じ canvas と editor drawer を中央 workspace に置く。drawer の containing block を workspace にする。editor の第二コピーや別 mobile editor を作らない。
2. geometry は index.html の単一 CSS owner へ収束させる。container 幅で wrap、min-width:0、必要な scroll を制御する。操作行を drawer/canvas の下に潜り込ませない。
3. 固定 `RIGHT_DRAWER_WIDTH_PX`、dragOffset/activeDrag、search position computed、start/move/stopDrag、search global drag listener、関連 Vue export、template style/drag binding を同じ変更で削除する。検索の意味を扱う機能は維持する。
4. 旧 search/toolbar absolute/translate と container-query の競合 override を削除する。新 grid/flex に旧 `!important` solver を残さない。無関係な feature popup 等の360pxは変更しない。座標 ref、ResizeObserver、collision solver は追加しない。
5. 通常高さ740px以上の受入条件で workspace>=200px。短い viewport、soft keyboard、200% zoom では Preview/card または文書の scroll で全操作へ到達できるようにする。fixed height + overflow clip で隠さない。
6. drawer開閉/animation、settings幅変更、resize で検索 state や focus を reset しない。Close/Escape は visibility-only の既存意味を維持する。狭幅 editor 上下 dock と alignment review、PD-OI-035 の mobile palette coverage は変更しない。

SRP は search意味/CSS geometry/drawer transitions の分離、DRY は同一DOMと単一geometry owner、KISS/YAGNI は専用行と旧JS位置経路の撤去で守る。互換性のために検索 drag を残す branch や新設定 toggle は不要。

## 検証

総合計画 P01–P03 と S02 の説明到達性を検査する。viewport 幅1440/1024/900/768/390、高さ844/740/480、container境界前後、settings resize、drawer open/close系列、200% zoom、soft keyboard の条件を整理して実ブラウザーで検査する。すべてを無意味に直積化せず、境界・系列・短画面をカバーする matrix と結果を明示する。

矩形非交差だけでは合格にしない。Previewを可視位置へscrollし、elementFromPoint と実 click/keyboard、clipping/到達性、740px以上条件のworkspace高さを検査する。keyboard/zoom はCSS縮小だけで代用せず、実表示倍率やvisual viewportに相当する操作の手順・制限を証拠に残す。

```bash
node --test tests/web/preview-feature-search.test.mjs
```

既存 `preview-navigation.playwright.spec.js` 等へ browser regression を追加し、全search機能、toolbar、drawer tab/Close/Escape、pan/zoom、同じDOM identity、S01配置継承・S02hintを確認する。Node/Python Playwright を調べ、必要な source一致 wheel を準備する。旧listener/computed/exports/CSSが残っていないことを `rg` と diff で確認する。新規 dependencies、科学的geometry/reference SVG変更は不要。

## 完了・commit/push

`results/S03.md` に入力/authority/検査 SHA、matrix、hit target/scroll/DOM/検索stateの結果、コマンド・artifact、制限、旧geometry paths撤去を記録し、進捗を更新する。production/tests/docs/generated を別々に review する。公開操作文書は S04 の既存 docs owner へ変更点を引き継ぐ。

branch/upstream を再確認して fetch。予期しない remote 更新を調べ、他者変更を上書きしない。対象パスのみ stage して `git diff --cached --check`。**終了時は必ず commit と push。**

```bash
git commit -m "fix(web): dock preview search and toolbar in dedicated rows"
git push origin HEAD:refs/heads/fix/issue-599-preview-layout-20260926
git rev-parse HEAD
git ls-remote --heads origin refs/heads/fix/issue-599-preview-layout-20260926
git status --short
```

local/remote SHA一致を確認。エラー時はremote実状態を調べてから再試行する。main/dev直接 push・force push・無承認のPR作成/merge/deployは禁止。公開 SHA、検証結果、未達条件、S04 開始条件、英語 commit title と short summary を報告する。
