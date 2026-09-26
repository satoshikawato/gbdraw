# 承認内容の機械表現

承認者・author: `satoshikawato`。承認日: `2026-09-26`。
以下は承認済みの各`PRODUCT_DECISION`本文から機械的に生成したレビュー用JSONである。各recordは選択された一つのoutcomeだけを含む。
これは承認内容の提示であり、CIが読む新しいdecision storeではない。両concernはdeveloper preflightのLane Bとして扱い、S01が既存の`docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md`へ独立recordとして記録する。正式なProduct authorityはbaseへの統合後に有効となる。

## PDF対応

原文: [Decision Pack 01](./DECISION_PACK_01_PDF_JA_ZH.md)。PDFの実現性証拠はS00で取得する。

```json
{
  "concern": "web.export.pdf-unsupported-glyph-output",
  "scenarioRevision": 4,
  "choice": "A",
  "optionId": "SELECTIVE_JA_ZH_VECTOR_TEXT_FALLBACK",
  "rationale": "日本語、中国語（簡体・繁体）を含む図も文字として保持したPDFで提出できるよう、既存fontで表現できるアルファベットなどを維持し、不足部分だけ対応fontへ自動切替する。",
  "mustPreserve": "既存fontで正しく表現できる部分のfamily/style、原文textと全feature/label/legend/比較関係、page寸法とanchor、クリック時snapshotとfilename、PNG DPI、成功Resultとcanonical request、Session/History、local-only処理、必要fontの遅延取得と失敗後retry、fallback後のglyph validation。検証対象は日本語・簡体/繁体中国語と既存文字の混在で、surrogate pairや結合単位を分断せず、文字の位置と原文のUnicode対応を保持する。文字として抽出・検索・選択できるfont埋込みPDFを通常操作で取得し、文字削除・画像化・path化・style低下を黙った代用にしない。language別の手動font選択を通常操作に要求しない。",
  "mayRetire": "同梱fallbackと検証済み配置で表現できる不足文字があるだけでPDF全体を拒否する動作。一次fontに不足がある単位について、対応fontへの部分的な代替を認める。",
  "acceptedResidualRisk": "対応fontによって配布容量と初回load/parse時間が増え、不足部分のtypefaceとmetricが一次fontとは異なる場合がある。明示text languageがない漢字は固定manifestの既定字形になり地域字形の最適化は保証しない。検証済みrepertoire/style外は原因と同snapshotのSVG/PNG回復を提示する。誤配置・文字化け・無告知の画像化・geometry破損、license/packaging/資源/独立renderer/原文text等のrequired evidence不合格、architecture/vendor/privilege Gate違反は認めない。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-26"
}
```

## 診断情報公開

原文: [Decision Pack 02](./DECISION_PACK_02_ERROR_DISCLOSURE.md)。実装の情報保持・非公開・復旧検証はS03/S04で実施する。

```json
{
  "concern": "web.errors.user-facing-diagnostic-disclosure",
  "scenarioRevision": 1,
  "choice": "A",
  "optionId": "ACTIONABLE_SUMMARY_WITH_SAFE_DETAILS",
  "rationale": "利用者が失敗した操作と既知の原因から修正・再試行でき、解決しない場合には入力内容を公開せず問題報告に必要な診断を確認できるようにする。",
  "mustPreserve": "既知のunderlying failureと利用可能な次の操作、Alignのdraftとretry、最後の成功Resultとcanonical request、record方向とHistory、Pythonによるrule validation、cancel/stale/supersededとerrorの区別、native CLI/APIの例外捕捉、local-only処理、actual stageの表示、unknown原因の非捏造、keyboard/390 px対応。Detailsは初期collapseの任意表示とし、code/stage/row/field/position/codepointなど許可済みで有界の診断だけを公開する。",
  "mayRetire": "利用者向け画面と通常consoleへのraw内部例外文、stack/traceback、自由形式stdout/stderr、raw cleanup noteの直接公開、Error型prefixの重複、一律のGeneration Error表示、原因を消すgeneric fallbackだけの案内。rawの除去は既知原因の消失を許可しない。",
  "acceptedResidualRisk": "未知の例外では実際の操作stage、原因不明という説明、安全な次の操作と有界diagnostic codeまでの案内になり、raw traceを利用者から回収できないため、開発者側で別途再現調査が必要になる。private data漏出、原因の捏造、validation無視、失敗時のdraft/Result/History損失、required acceptance evidenceの不合格は認めない。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-26"
}
```
