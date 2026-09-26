# Decision Pack 02 — 原因のある復旧案内と安全な診断情報

承認者・author: `satoshikawato`。承認日: `2026-09-26`。
Concern: `web.errors.user-facing-diagnostic-disclosure`。Scenario revision: `1`。
選択: **A / ACTIONABLE_SUMMARY_WITH_SAFE_DETAILS**。

## 承認した製品動作

summaryには失敗した操作、既知の原因、利用可能な次の操作を表示する。Detailsは初期collapseの任意表示とし、code、stage、row、field、position、codepointなど許可済みで有界の診断情報だけを公開する。summaryは常に可視で、実操作に合う見出し、keyboardと390 pxの表示を提供する。

既存の上限であるsummary 1,000字、detail text 4,000字、section 8個を維持する。paramsの型と長さをallowlistで検査する。raw内部例外、stack/traceback、自由形式stdout/stderr、入力sequence、ファイル全文、SVG、任意exception messageは画面や通常consoleへ公開しない。未知原因は不明と表示し、既知の修正情報をunknownへ潰さない。

失敗時にはdraft、最後の成功Result、canonical request、record方向、Historyを保持する。cancel/stale/supersededをerrorへ変換しない。error情報はtransientであり、保存Sessionや成功artifactへ追加しない。native CLI/APIの例外捕捉とPython rule validationは維持する。

## 実装前のauthorityと検証

製品動作は承認済み。S01でこのoutcomeをPDF対応とは独立したstatic Product Contract recordへ記録し、baseへ統合した後にS03で既存error経路を修正する。source→Worker→caller→UIの情報保持、raw非公開、rollback/retryと操作可能な案内が実装受入条件となる。

## 承認本文

```text
PRODUCT_DECISION
Concern: web.errors.user-facing-diagnostic-disclosure
Scenario revision: 1
Choice: A / ACTIONABLE_SUMMARY_WITH_SAFE_DETAILS
Rationale: 利用者が失敗した操作と既知の原因から修正・再試行でき、解決しない場合には入力内容を公開せず問題報告に必要な診断を確認できるようにする。
Must preserve: 既知のunderlying failureと利用可能な次の操作、Alignのdraftとretry、最後の成功Resultとcanonical request、record方向とHistory、Pythonによるrule validation、cancel/stale/supersededとerrorの区別、native CLI/APIの例外捕捉、local-only処理、actual stageの表示、unknown原因の非捏造、keyboard/390 px対応。Detailsは初期collapseの任意表示とし、code/stage/row/field/position/codepointなど許可済みで有界の診断だけを公開する。
May retire: 利用者向け画面と通常consoleへのraw内部例外文、stack/traceback、自由形式stdout/stderr、raw cleanup noteの直接公開、Error型prefixの重複、一律のGeneration Error表示、原因を消すgeneric fallbackだけの案内。rawの除去は既知原因の消失を許可しない。
Accepted residual risk: 未知の例外では実際の操作stage、原因不明という説明、安全な次の操作と有界diagnostic codeまでの案内になり、raw traceを利用者から回収できないため、開発者側で別途再現調査が必要になる。private data漏出、原因の捏造、validation無視、失敗時のdraft/Result/History損失、required acceptance evidenceの不合格は認めない。
Owner: satoshikawato
Decision date: 2026-09-26
```

PDFのfont fallbackと出力品質は[Decision Pack 01](./DECISION_PACK_01_PDF_JA_ZH.md)が別に扱う。
