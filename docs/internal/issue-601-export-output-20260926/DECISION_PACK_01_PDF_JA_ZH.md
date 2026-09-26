# Decision Pack 01 — 日本語・中国語のPDF出力

承認者・author: `satoshikawato`。承認日: `2026-09-26`。
Concern: `web.export.pdf-unsupported-glyph-output`。Scenario revision: `4`。
選択: **A / SELECTIVE_JA_ZH_VECTOR_TEXT_FALLBACK**。

## 承認した製品動作

日本語と中国語（簡体字・繁体字）を含む図のPDF出力に対応する。既存fontで表せるアルファベット・数字・記号とそのstyleを維持し、不足する部分だけ同梱の対応fontへ自動切替する。PDF内に文字を保持し、検索・選択できる出力を提供する。

font選択は通常のPDFボタンの処理で自動実行する。クリック時の図とfilenameを同期captureし、成功ResultやSession、Historyを変更せず、一時的なexport cloneを使う。fallback後にもglyph coverageを検査する。font取得失敗・未対応文字・未対応style・配置失敗を区別し、原因と同snapshotのSVG/PNG取得または再試行を案内する。

明示されたtext languageは地域字形選択へ利用する。情報がない漢字は固定manifestのcoverage優先順でfontを選ぶ。UIの言語やOS fontの偶然に依存させない。地域字形の言語別最適化は、text languageがない場合には保証しない。

必要fontだけsame-originから遅延取得し、失敗した取得を再試行できるようにする。font原本、版、checksum、license、styleと生成assetの対応を一つの準備経路で管理する。既存のベクターPDF経路でrunを切り替え、原文Unicodeとglyph位置を保つ。

## 実装前の証拠とauthority

製品動作は承認済み。日本語の小規模なtext PDF probeは成功している。中国語のfont選定、製品の自動fallback、独立PDF renderer、混在runの配置、資源・license・source/wheel配信の受入は必要な確認として残る。

[総合計画書](./MASTER_PLAN.md)のS00で実現性証拠を取得し、S01で承認内容を独立したstatic Product Contractへ記録する。承認は欠字、文字化け、原文抽出不一致、required evidenceの失敗やarchitecture Gate違反を許可するものではない。

## 承認本文

```text
PRODUCT_DECISION
Concern: web.export.pdf-unsupported-glyph-output
Scenario revision: 4
Choice: A / SELECTIVE_JA_ZH_VECTOR_TEXT_FALLBACK
Rationale: 日本語、中国語（簡体・繁体）を含む図も文字として保持したPDFで提出できるよう、既存fontで表現できるアルファベットなどを維持し、不足部分だけ対応fontへ自動切替する。
Must preserve: 既存fontで正しく表現できる部分のfamily/style、原文textと全feature/label/legend/比較関係、page寸法とanchor、クリック時snapshotとfilename、PNG DPI、成功Resultとcanonical request、Session/History、local-only処理、必要fontの遅延取得と失敗後retry、fallback後のglyph validation。検証対象は日本語・簡体/繁体中国語と既存文字の混在で、surrogate pairや結合単位を分断せず、文字の位置と原文のUnicode対応を保持する。文字として抽出・検索・選択できるfont埋込みPDFを通常操作で取得し、文字削除・画像化・path化・style低下を黙った代用にしない。language別の手動font選択を通常操作に要求しない。
May retire: 同梱fallbackと検証済み配置で表現できる不足文字があるだけでPDF全体を拒否する動作。一次fontに不足がある単位について、対応fontへの部分的な代替を認める。
Accepted residual risk: 対応fontによって配布容量と初回load/parse時間が増え、不足部分のtypefaceとmetricが一次fontとは異なる場合がある。明示text languageがない漢字は固定manifestの既定字形になり地域字形の最適化は保証しない。検証済みrepertoire/style外は原因と同snapshotのSVG/PNG回復を提示する。誤配置・文字化け・無告知の画像化・geometry破損、license/packaging/資源/独立renderer/原文text等のrequired evidence不合格、architecture/vendor/privilege Gate違反は認めない。
Owner: satoshikawato
Decision date: 2026-09-26
```

診断情報の公開範囲は[Decision Pack 02](./DECISION_PACK_02_ERROR_DISCLOSURE.md)が別に扱う。
