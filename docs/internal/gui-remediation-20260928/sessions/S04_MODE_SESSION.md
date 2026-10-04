# S04 INSTRUCTION PROMPT: モード別設定とSession復元

Circular SessionからLinearへ切り替えたときのtitle・行配置・表示設定の混入を修正してください。両モードを編集して保存した設定の互換性も検証します。

## 前提と場所

**fix/gui-feedback-remediation-20260928**、/mnt/c/Users/genom/GitHub/gbdraw/.worktrees/gui-feedback-remediation-20260928を使用。clone/別worktree/別branchを作らない。
[総合計画](../00_MASTER_PLAN.md)、[SESSION_LOG](../SESSION_LOG.md)、results/S00.mdとAGENTS/CLAUDE/Web CLAUDE、docs/SESSION_COMPATIBILITY.mdを読む。
同時更新ルールはS02がbaseへ反映済みであること。旧flat Sessionの非active値の扱いとper-mode対象は、S00で記録された具体的なProduct判断を使う。

## 実装

1. 既存session-request/config/active-config-contract境界で、欠落と明示false/Show/Autoを区別する。Circular requestが持たないLinear設定をfalse/Showとして作らない。
2. 未訪問Linearはtitle空、Arrange in rows=true、Replicon=false、Accession/Length=autoを維持する。Autoは共有行topologyに従って解決する。
3. titleと決定済み表示設定を既存mode-profilesまたはlayout-preferencesの一つの責任箇所へ統合する。第三のmode snapshot managerや全設定resetは作らない。
4. writer/reader/History/reset/mode往復を同じownerの遷移へ集約する。保存形式の変更が必要なら公開済みschemaの証拠とfixtureを示し、branch-only migrationを積み重ねない。
5. 明示inactive draftの扱いは採用された契約の範囲に従う。filename/default一致/source不存在をユーザー意図の判定に使わない。
6. gallery-session-publicationの未使用mode既定値を修正する。固定した旧Vnigの入力は保存し、修正版fixtureで上書きしない。
7. Product Contract/Session説明とruntime/testsを同じ差分で整合させる。per-mode化の範囲を勝手に全formへ拡張しない。

## 判断が不足する場合

旧flat Sessionのinactive Showと手動保存Showを区別できない場合に、片方を都合よく消さない。S00のDecision Packへ具体的な欠落判断を返す。独立した欠落既定値修正・テストは進めてよいが、未承認の互換変更を実装せず、G03全体を完了としない。

## 受入

- main v41/dev v44のVnigをLoad→Linear→指定Vibrioの2入力→Generate。選択した旧形式方針に基づく期待値を明記する。
- 新規/更新Galleryでtitle空、行group維持、Accession/Length auto、Replicon false。
- Linear→CircularとCircular→Linearの往復、両modeを明示編集したsettings-only Session、Save→fresh Load。
- 明示OFF/Show、title/fonts、既存layoutPreferences、比較設定を保持。
- Load時は保存Resultを勝手に再生成せず、失敗importは元のsource/config/Resultへrollback。
- Replicon ONの報告はcheckbox/model/描画を別々に観測する。未再現を修正済みと書かない。

既存mode-profiles、session-request、settings-only-session、linear-record-layoutのunit/browser testsへ必要な回帰ケースを追加する。

## 終了

results/S04.mdにfield-owner表、旧/新fixtureの結果、compatibility evidence、Product判断と実装の対応を記録。Gallery生成物が必要なら所定skillとowner toolを使い、S07へ視覚確認の残件を明記する。SESSION_LOGを更新して担当差分をcommitし、英語title/summaryを残す。
