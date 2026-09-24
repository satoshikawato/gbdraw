# Issue #581 — S02 Non-modal floating palette INSTRUCTION PROMPT

以下を一つの独立した作業指示として使用する。過去の会話は前提にしない。

## 目的と開始条件

gbdraw Issue #581 の Similarity Alignment 選択 UI を、比較図を見られる non-modal な
floating palette に変更する。設計提案は
https://github.com/satoshikawato/gbdraw/issues/581#issuecomment-5811037046、
完全な実装境界は
`docs/internal/issue-581-feature-popup-alignment-ux/IMPLEMENTATION_PLAN.md` にある。
S01 のローカル選択・Apply batch 検証が完了したこと、2026-09-24 に承認された
`A / FLOATING_GUIDE_CANVAS_PICK`（PD-OI-035）が PR #582 で
`origin/dev` にマージされ、固定実装ブランチへ取り込まれたことを確認する。

## ブランチと資料

- runtime 作業は **`issue-581-feature-popup-alignment-ux-20260924`** に積む。
  branch、HEAD、upstream、worktree、総合計画第10節、S01 の diff と test を確認する。
- `AGENTS.md`、`CLAUDE.md`、`gbdraw/web/CLAUDE.md`、
  `docs/internal/PRODUCT_IMPACT_RATCHET.md`、
  `docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md` と、
  `docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md` の PD-OI-026、033〜035 を読む。
- `gbdraw/web/index.html` の現行 dialog、`gbdraw/web/js/app/app-setup.js` の
  focus trap／return、`gbdraw/web/js/app/similarity-alignment.js` の state、
  既存の popup drag、responsive layout、関連 browser tests を読む。

## 実装

1. 全画面の暗い backdrop、`aria-modal="true"`、Tab trap を外す。
   名前と説明を持つ `role="dialog"` は維持し、初回 open では palette の選択 control へ
   focus を移す。閉じたときは、元の invoker がまだある場合にだけ focus を戻す。
2. palette の header を drag handle にする。位置は画面内に clamp し、browser resize
   でも操作可能にする。位置は Session、render request、History に保存しない。
   drag は選択 radio、Apply、Close、text selection、図の pan／zoom と競合させない。
3. desktop では比較図の有効幅を減らさず、palette 本体だけに適切な最大幅・最大高と
   scroll を付ける。390 px の narrow viewport では画面外に control が出ず、選択、
   Skip、Apply が到達可能な配置にする。最初の配置と reset で mouse による drag を
   必須にしない。
4. background の diagram pan／zoom を可能にする。palette 外の操作による
   source／crop／group／committed Result の変更には S01 の stale 判定を適用し、
   古い draft を黙って Apply しない。別 dialog との focus と Escape の優先順位を確認する。
5. candidate の主行、Details disclosure、resolution count、Apply reason を読みやすく
   整える。HTML template に Resolver 規則や候補順位の判定を入れない。

## 検証と終了

- 実ブラウザで desktop（例: 1366 px）と 390 px を目視し、図の width と vertical
  record context、palette の drag と clamp、スクロール、Apply の可視性を確認する。
- keyboard だけで candidate／Skip／Apply／Cancel を操作でき、Tab が palette に
  閉じ込められず、Escape と focus return が既存 popup と衝突しないことを確認する。
- palette を動かしても diagram Result／Session／History が変わらず、invalidating
  mutation で古い選択を Apply できないことを確認する。
- Node の `@playwright/test` と Python Playwright の両方の可用性を確認する。
  Node spec が動かない環境でも Python Playwright で実ブラウザを確認する。
- production、test、visual diff を別々に監査し、総合計画第10節へ結果と S03 開始可否を
  記す。S03 の guide／badge が未実装なら、Issue 全体の完成とは報告しない。

SOLID: palette は表示・focus のみ所有し、alignment controller が選択を所有する。
KISS/DRY: 既存 Vue template と drag の小さい再利用範囲で済ませる。
YAGNI: layout framework、永続化、第二の modal 管理器は作らない。

## コミット・push と次セッションへの引き継ぎ

本依頼は担当分のコミットと同名 remote work branch への push を明示的に許可する。
このセッションの作業は固定ブランチ `issue-581-feature-popup-alignment-ux-20260924` 内で完了する。
検証と差分監査の後、担当分を英語の題名で１コミットにまとめ、同名の remote work
branch へ `git push origin HEAD:refs/heads/issue-581-feature-popup-alignment-ux-20260924` で push する。
開始前と push 前に branch、upstream、作業ツリー、remote の状態を確認する。

回答の最後に、次の S03 `SESSION_03_CANVAS_INTEGRATION.md` を新規参加者が単独で実行できる
完全な INSTRUCTION PROMPT として提示する。そのプロンプトにも、同じ固定ブランチで
作業・検証・コミット・push まで行い、さらに次セッション用の完全なプロンプトを
回答末尾に提示する指示を含める。S04 が完了した場合は追加セッションを作らない。
