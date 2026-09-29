# S06 INSTRUCTION PROMPT: Preview配置と必要時ヘルプ

検索バーを小さく戻し、通常幅のEditorをPreview上端へ配置して検索を退避させてください。Align/Reviewを横並びにし、指定説明を必要時のhelpへ移します。

## 前提と場所

**fix/gui-feedback-remediation-20260928**、/mnt/c/Users/genom/GitHub/gbdraw/.worktrees/gui-feedback-remediation-20260928を使用。
[総合計画](../00_MASTER_PLAN.md)、[SESSION_LOG](../SESSION_LOG.md)、S02/S03/S05の結果、AGENTS/CLAUDE/Web CLAUDEを読む。
S02の同時更新例外がbaseへ反映済みであること。S05のcontroller/表示名を利用し、別actionを作らない。

## 所有と実装

所有はgbdraw/web/index.htmlの検索/Editor/Align/Review/説明markupとCSS、既存components.jsのHelpTip、必要なpreview/search配線、対応tests/docs。
S03が担当したPendingとoperation/errorの意味は維持する。

1. PD-OI-054の検索専用行、PD-OI-024の常時Lock説明を今回の結果に改訂し、runtime/testsと同じ差分でレビュー可能にする。既に選ばれた結果を再質問しない。
2. 検索はmainの最大39.5remを基準に利用可能幅内へ収める。通常幅ではEditorがPreview上端から開き、検索を残り幅へ退避する。
3. 同じdrawer幅CSS変数を利用する。固定360pxのJS移動、新しい座標ref/observer、検索の自由drag、旧commitの全revertを追加しない。
4. 狭幅の既存Editor/review契約を保持する。検索の折返し/スクロールとcanvas/Close/toolbar到達を成立させる。
5. popupとdrawerの両方でAlignの右隣にReview alignment optionsを置く。busy/disabled/exact-reference actionは既存controllerを使い、狭幅で切断しない。
6. Alignの長い常時説明とON aligns definitions…を既存HelpTipへ集約する。hover/focus/tap、Escape/blur、accessible descriptionを扱い、一つのtext sourceを使う。
7. Current: Run LOSAT…とprogram/threshold要約を削除する。選択済みボタンにaria-pressed等で状態を示し、Customと実エラーは必要時に保持する。
8. capture/testが削除した文字列を待つ場合、実際の選択状態とoperation settlementへ変更する。単なるsleepや内部refの直書きで置き換えない。

## 受入

viewport 1440×900、1024×740、768×740、390×844、390×740と200% zoomを確認する。
通常幅では検索とdrawerが重ならず、検索が場外へ出ず、Editorが検索より下の専用行から始まらない。
狭幅はcanvasの利用可能幅と通常条件で200px以上の高さ、独立scroll、短い高さ/soft keyboard時の操作到達を保持する。

検索の全field、regex、query、Prev/Next/Open/Enter、active match、focus、camera、Editor tab/Close/Escapeを確認する。配置変更でSVGをcloneしない。query/Result/History/Exportを変更しない。

helpのpointer/keyboard/touch操作を確認し、hover専用へ退行しない。説明の削除を理由に実エラーやrecoveryを消さない。CW-01の「表示変更だけで全feature計算しない」も検証する。

## 検証とhandoff

preview-navigation、right-drawer、similarity-alignment-ui、comparison-uiの対応Playwrightとunit testsを実行。Gallery/tutorialのcapture変更時は該当skillを適用し、生成ownerから再現する。完成画面を目視確認する。

results/S06.mdにviewport結果、tooltip、契約と差分の対応、スクリーンショット根拠を記録する。SESSION_LOGを更新し、production/tests/docs/generatedを分けてレビューしてcommit、英語title/summaryを残す。
