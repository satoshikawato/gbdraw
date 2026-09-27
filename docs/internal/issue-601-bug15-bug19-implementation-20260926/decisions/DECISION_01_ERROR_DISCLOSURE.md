# 承認済みProduct Decision — エラー診断情報の公開範囲

Status: Product Decision Ownerが2026-09-26に選択Aを承認済み。
Product Decision Owner: satoshikawato。
対象: Issue #601 の BUG-15 / BUG-19。B案は採択していない。

## 選択された製品動作

短い修正案内と、利用者が任意に開ける有界のDetails、表示中の安全な内容だけを手動コピーするCopy diagnosticsを提供する。

本文の各条件が実装の境界である。別の関心への承認を推測したり、残余リスクを追加したりしない。
このファイルは承認内容のレビュー・引継ぎ記録であり、新しいCI decision storeや正式Product authorityではない。
正式な未登録concernのauthorityは、既存OPTION_INTEGRITY_PRODUCT_CONTRACT.mdへの記録がdevに統合された後に有効となる。
[総合計画書](../MASTER_PLAN.md)のS01→S02の順序を守る。

## 承認された本文

以下の本文は選択Aの記入済みreceiptをそのまま保存したもの。
receipt SHA-256（UTF-8、末尾改行を除く）: 8617888cb1838521f78640db4653e2dcf88dd427cd466bbbc33d59efacaceae8。

~~~text
PRODUCT_DECISION
Concern: web.errors.diagnostic-disclosure
Scenario revision: 1
Choice: A / GUIDANCE_WITH_BOUNDED_DIAGNOSTICS
Rationale: 利用者が短い修正案内から作業を続けられ、必要な場合は入力内容を公開せずに安全な失敗種別と段階を調査へ渡せるようにする。
Must preserve: すべての移行対象で既知 validation の修正情報を保持する。Generate/Align の以前の Result、canonical request、draft、orientation、History、retry、Save/Export、cancel/stale/superseded を保つ。Details は keyboard で開け、Copy diagnostics は表示中の bounded code/operation/stage/許可 context/副因だけを手動コピーする。unknown は stage と stable code を示し、元 pattern、sequence、file/record 名、path、SVG、自由な exception/traceback/stdout/stderr を画面・Copy・console に自動公開しない。初回と旧 Result 保持を区別する。
May retire: user-facing raw exception/traceback と自由な stdout/stderr、個別の例外型 prefix の直接表示。安全な修正情報、Details 入口、既存 recovery は退役しない。
Accepted residual risk: bounded 診断だけでは稀な未知例外を特定できず、利用者の明示的な Session 保存・別途再現情報が必要になる場合がある。Clipboard 不可時も表示情報の手動選択コピーと通常 recovery を維持する。
Owner: satoshikawato
Decision date: 2026-09-26
~~~

## 機械表現

本文の九項目だけを、同じ意味のJSONへ機械的に変換した。
scenarioRevisionだけを整数にし、rationale・mustPreserve・mayRetire・acceptedResidualRiskは原文のまま保存した。
これはレビュー用の表示であり、独立したJSON registryではない。

~~~json
{
  "concern": "web.errors.diagnostic-disclosure",
  "scenarioRevision": 1,
  "choice": "A / GUIDANCE_WITH_BOUNDED_DIAGNOSTICS",
  "rationale": "利用者が短い修正案内から作業を続けられ、必要な場合は入力内容を公開せずに安全な失敗種別と段階を調査へ渡せるようにする。",
  "mustPreserve": "すべての移行対象で既知 validation の修正情報を保持する。Generate/Align の以前の Result、canonical request、draft、orientation、History、retry、Save/Export、cancel/stale/superseded を保つ。Details は keyboard で開け、Copy diagnostics は表示中の bounded code/operation/stage/許可 context/副因だけを手動コピーする。unknown は stage と stable code を示し、元 pattern、sequence、file/record 名、path、SVG、自由な exception/traceback/stdout/stderr を画面・Copy・console に自動公開しない。初回と旧 Result 保持を区別する。",
  "mayRetire": "user-facing raw exception/traceback と自由な stdout/stderr、個別の例外型 prefix の直接表示。安全な修正情報、Details 入口、既存 recovery は退役しない。",
  "acceptedResidualRisk": "bounded 診断だけでは稀な未知例外を特定できず、利用者の明示的な Session 保存・別途再現情報が必要になる場合がある。Clipboard 不可時も表示情報の手動選択コピーと通常 recovery を維持する。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-26"
}
~~~
