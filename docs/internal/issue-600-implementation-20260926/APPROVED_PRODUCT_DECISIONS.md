# Issue #600 — 承認済み Product Decisions

Decision date: **2026-09-26**。下記４つの outcome はすべて Choice A として承認済み。
この文書は受領内容を保持する実装 handoff であり、既存静的 Product Contract の代替
authority や新しい decision store ではない。S00 がこの complete outcome を既存
`OPTION_INTEGRITY_PRODUCT_CONTRACT.md` の４つの独立 record として登録する。

Approval receipt: repository task の Product approval は「**すべて推奨案で承認します。**」。
この approval が選ぶ内容は、下記の complete outcome / rationale / must preserve /
may retire / accepted residual risk と stable choice ID で固定する。別の新 outcome、
public-contract extension、追加 retirement は含めない。
Product Decision Owner: **satoshikawato**（明示回答により確認済み）。
４ outcome の再承認や署名者の再確認は不要。S00 はこの受領内容を既存 authority へ
転記し、生成表現と対応を review 可能にする。

| Concern | Approved stable choice |
| --- | --- |
| `annotations.table-auxiliary-columns` | `A / AUX_COLUMNS_WARN_IGNORE` |
| `annotations.feature-selector-miss` | `A / MISSING_SELECTOR_SKIP_ROW_WARN` |
| `styles.specific-color-caption-multiplicity` | `A / CAPTION_AUTO_DISAMBIGUATE_SOLID_ROWS` |
| `tracks.pixel-text-input-domain` | `A / PIXEL_TEXT_OPTIONAL_PX` |

## annotations.table-auxiliary-columns

- Outcome status: APPROVED / A / AUX_COLUMNS_WARN_IGNORE
- Durable registration status: S00 で base authority を照合・登録する
- Complete outcome: Web/CLI/Python の Annotation TSV は任意の未知 header を受理し、その列を捨てて既知列のみ import する。１表につき列名を集約して「無視され、Session/TSV 再出力に保存されない」と通知。fill_colour 等の typo も未知列として通知する。header にない余剰 cell、必須欠落・正規化後重複・不正な既知値は全 import 拒否。
- User access / feedback: Web import 操作直後に読み上げ可能な status と列名一覧。CLI logger の集約 warning。通知には cell contents を含めない。
- Session / regeneration: annotation 値のみ保存。Load で未知列復元や自動 Generate をしない。再生成は既知列だけの import と同じ。
- Export / artifact: TSV writer は既存 column inventory のみ。付加列の lossless export はしない。SVG に未知 metadata を入れない。
- Failure / recovery: import は成功、利用者は通知を確認して編集・Generate へ進める。誤字だった場合は原 TSV を修正して再 import。known-invalid 時は直前 state のまま。

```text
PRODUCT_DECISION
Concern: annotations.table-auxiliary-columns
Scenario revision: 1
Choice: A / AUX_COLUMNS_WARN_IGNORE
Rationale: 生物学的な annotation に使う列の意味を検証しつつ、解析 TSV の付加 metadata だけで作図を止めない。取り込まれない列を明示し、利用者が誤字や非保存を判断できるようにする。
Must preserve: 既知 annotation の値、行/集合順序、strict typed schema、valid-input export、失敗時の既存 draft/Result。 unknown field を typed annotation に通さない。必須列、duplicate、target/known enum/数値/style を検証する。malformed row を付加列と誤認しない。未知列の cell contents を console に出さない。
May retire: Annotation TSV に対する「unknown header はすべて fatal」の契約のみ。records/track 等の他の表の unknown policy は退役しない。
Accepted residual risk: optional typo が無視され、デフォルト style になる可能性。列名と非保存の通知、known-required/known-value strict 検証で範囲を制限。
Owner: satoshikawato
Decision date: 2026-09-26
```

Inert machine representation（既存静的契約への転記用。独立した JSON registry ではない）:

```json
{
  "concern": "annotations.table-auxiliary-columns",
  "scenarioRevision": 1,
  "choice": "A / AUX_COLUMNS_WARN_IGNORE",
  "rationale": "生物学的な annotation に使う列の意味を検証しつつ、解析 TSV の付加 metadata だけで作図を止めない。取り込まれない列を明示し、利用者が誤字や非保存を判断できるようにする。",
  "mustPreserve": "既知 annotation の値、行/集合順序、strict typed schema、valid-input export、失敗時の既存 draft/Result。 unknown field を typed annotation に通さない。必須列、duplicate、target/known enum/数値/style を検証する。malformed row を付加列と誤認しない。未知列の cell contents を console に出さない。",
  "mayRetire": "Annotation TSV に対する「unknown header はすべて fatal」の契約のみ。records/track 等の他の表の unknown policy は退役しない。",
  "acceptedResidualRisk": "optional typo が無視され、デフォルト style になる可能性。列名と非保存の通知、known-required/known-value strict 検証で範囲を制限。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-26"
}
```

## annotations.feature-selector-miss

- Outcome status: APPROVED / A / MISSING_SELECTOR_SKIP_ROW_WARN
- Durable registration status: S00 で base authority を照合・登録する
- Complete outcome: binding 成功済みの annotation に１件でも未一致 feature selector があれば、その record に対する annotation 行全体を skip し、code/set/annotation/record 識別と欠落件数を持つ warning を返す。他の行/record は継続。全注釈 missing でも genome 図は正常に返し、skip 件数を表示する。
- User access / feedback: Generate 成功後の status に skip 件数と row/record 識別を表示。CLI warning、API の structured warning を提供。未一致の qualifier 値を console に dump しない。
- Session / regeneration: selector と row を保存し、次回はそのときの record に再解決。Session Load で自動 Generate をしない。保存 preview は保持する。
- Export / artifact: SVG/PNG/PDF に skipped mark を出さない。annotation TSV は元 row を含み、次の入力で再利用できる。empty mark の legend は作らない。request に含まれた explicit slot と Web の既存自動 projection の配置/gap は維持し、skip を理由に縮小しない。Python が resolved marks から新規 auto slot を作る場合は empty set の slot を作らない。
- Failure / recovery: 成功図を確認して selector を修正・削除・別 record を明示して再 Generate できる。構造エラー時は直前 Result/draft を保ち修正へ。

```text
PRODUCT_DECISION
Concern: annotations.feature-selector-miss
Scenario revision: 1
Choice: A / MISSING_SELECTOR_SKIP_ROW_WARN
Rationale: gene の欠落で比較図全体を失敗させず、複数 anchor で指定した annotation の意味も保つ。部分的な範囲の図示を自動で選ばず、欠落行のスキップを利用者に明示する。
Must preserve: 完全一致行の geometry、既存 record 意味、coordinate policy、crop/reverse/rotation、他の注釈、failure/cancel/stale 隔離。 record 欠落/曖昧/index 範囲外、multi-record の record 省略、malformed selector は fatal。coordinate clip/skip/error、transform、selector matching の意味を維持。任意 exception を skip にしない。
May retire: feature selector miss の blanket fatal だけ。record/syntax/coordinate error の fatal は維持。
Accepted residual risk: gene typo でも図が成功する。skip を表示することで隠れた欠落を防ぐ。一部 anchor が正しくてもその行の有用な mark は表示されない。request に含まれた注釈 slot は空き領域として残り得る。
Owner: satoshikawato
Decision date: 2026-09-26
```

Inert machine representation（既存静的契約への転記用。独立した JSON registry ではない）:

```json
{
  "concern": "annotations.feature-selector-miss",
  "scenarioRevision": 1,
  "choice": "A / MISSING_SELECTOR_SKIP_ROW_WARN",
  "rationale": "gene の欠落で比較図全体を失敗させず、複数 anchor で指定した annotation の意味も保つ。部分的な範囲の図示を自動で選ばず、欠落行のスキップを利用者に明示する。",
  "mustPreserve": "完全一致行の geometry、既存 record 意味、coordinate policy、crop/reverse/rotation、他の注釈、failure/cancel/stale 隔離。 record 欠落/曖昧/index 範囲外、multi-record の record 省略、malformed selector は fatal。coordinate clip/skip/error、transform、selector matching の意味を維持。任意 exception を skip にしない。",
  "mayRetire": "feature selector miss の blanket fatal だけ。record/syntax/coordinate error の fatal は維持。",
  "acceptedResidualRisk": "gene typo でも図が成功する。skip を表示することで隠れた欠落を防ぐ。一部 anchor が正しくてもその行の有用な mark は表示されない。request に含まれた注釈 slot は空き領域として残り得る。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-26"
}
```

## styles.specific-color-caption-multiplicity

- Outcome status: APPROVED / A / CAPTION_AUTO_DISAMBIGUATE_SOLID_ROWS
- Durable registration status: S00 で base authority を照合・登録する
- Complete outcome: 同 caption・異色 rule を受け付け、全色に lowercase normalized hex を付けた caption を canonical rule として採用する（例 Transporter [#112233] / Transporter [#445566]）。同名同色は共有、空 caption は凡例なし。既存 literal caption/legend key に衝突する場合は予約後に決定的な追加 suffix で区別。first/last-wins は廃止。各実際に使用された色を別の solid 凡例行で示す。
- User access / feedback: import/manual edit の正常完了時に caption 変更を通知。利用者は生成した solid 行を既存 editor から編集できる。
- Session / regeneration: admitted caption は普通の文字列として保存。Load は preview/draft を保ち自動 Generate しない。過去の同名異色 draft は次の rule edit/Generate の通常 preparation で通知付き正規化。
- Export / artifact: 新 TSV は区別した canonical caption。SVG/PNG/PDF は各色の solid 行。元ファイルの同名 caption のままの lossless 復元は約束しない。
- Failure / recovery: 自動区別後すぐ図を使える。必要なら caption を編集して再生成。stale/preparation failure は直前 rules/Result に戻す。

```text
PRODUCT_DECISION
Concern: styles.specific-color-caption-multiplicity
Scenario revision: 1
Choice: A / CAPTION_AUTO_DISAMBIGUATE_SOLID_ROWS
Rationale: 近い色を使う rule を拒否せず、実際の各色を凡例に表示する。今回は既存 solid 行を再利用する自動 caption 区別を採用し、複数 swatch 用の renderer・editor・保存形式を追加せずに fresh/live/native の意味を揃える。
Must preserve: feature 色、rule 順序/precedence、同名同色共有、single-color caption、unused rule の凡例除外、既存 solid editor・保存 preview・failure/History 契約。 rule order/regex/precedence/visibility を変えない。使用色の忠実な図示、stable identity、file/manual provenance、Result rollback、既存 SVG sanitizer を維持。
May retire: caption 衝突による Web 拒否、過去の last-wins/上書き。退役は specific-color rule の同名異色 scope に限定。
Accepted residual risk: 凡例が長くなり元の同名文字列は変わる。multi-swatch grouping を望む利用者には複数行になる。hex suffix、既存 edit、layout 再計測で扱いを明確にする。
Owner: satoshikawato
Decision date: 2026-09-26
```

Inert machine representation（既存静的契約への転記用。独立した JSON registry ではない）:

```json
{
  "concern": "styles.specific-color-caption-multiplicity",
  "scenarioRevision": 1,
  "choice": "A / CAPTION_AUTO_DISAMBIGUATE_SOLID_ROWS",
  "rationale": "近い色を使う rule を拒否せず、実際の各色を凡例に表示する。今回は既存 solid 行を再利用する自動 caption 区別を採用し、複数 swatch 用の renderer・editor・保存形式を追加せずに fresh/live/native の意味を揃える。",
  "mustPreserve": "feature 色、rule 順序/precedence、同名同色共有、single-color caption、unused rule の凡例除外、既存 solid editor・保存 preview・failure/History 契約。 rule order/regex/precedence/visibility を変えない。使用色の忠実な図示、stable identity、file/manual provenance、Result rollback、既存 SVG sanitizer を維持。",
  "mayRetire": "caption 衝突による Web 拒否、過去の last-wins/上書き。退役は specific-color rule の同名異色 scope に限定。",
  "acceptedResidualRisk": "凡例が長くなり元の同名文字列は変わる。multi-swatch grouping を望む利用者には複数行になる。hex suffix、既存 edit、layout 再計測で扱いを明確にする。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-26"
}
```

## tracks.pixel-text-input-domain

- Outcome status: APPROVED / A / PIXEL_TEXT_OPTIONAL_PX
- Durable registration status: S00 で base authority を照合・登録する
- Complete outcome: 純 pixel track geometry（Linear height/spacing、Circular inner_gap_px/outer_gap_px）の文字列入口は finite decimal/exponent と optional px（大文字小文字・前後空白可）を受理。trim 後空欄/null は auto。10、10px、10PX、10 px は同値。height は正、gap/spacing は非負。不正文字列/単位/非有限を拒否し、draft で保持して row error を出す。
- User access / feedback: 対象 field の help/placeholder を「px optional」に統一。field 名と正/非負条件を row error と CLI error で示す。
- Session / regeneration: 現行 canonical 型だけ保存。同値 input は同じ geometry。Load は既存 preview を保つ。text acceptance は新 migration ではない。
- Export / artifact: 同値 input の SVG/download は同じ。TSV/CLI の書き出しは既存 canonical 数値形式でよい。
- Failure / recovery: row error の値を編集して再 submission。失敗時に直前 Result を保つ。

```text
PRODUCT_DECISION
Concern: tracks.pixel-text-input-domain
Scenario revision: 1
Choice: A / PIXEL_TEXT_OPTIONAL_PX
Rationale: 利用者が pixel 値を単位付きで paste できるようにし、検証・正規化・request の値を揃える。物理 pixel と factor scalar は分けたまま、保存形式を増やさずに入力の一貫性を改善する。
Must preserve: 既存 valid 数値、Linear px acceptance、auto、physical pixel 意味、現行 typed request/Session、Circular radius/width factor/%、retired key 拒否。 typed JSON の gaps は数値、Linear は既存 ScalarSpec。Circular ratio/% semantics を保つ。不正値を null/0 化しない。不要な arbitrary CSS unit conversion を作らない。
May retire: pure pixel 対象の without-a-unit restriction と、invalid→null/zero の黙示的変換。一般 dimension input の制限は退役しない。
Accepted residual risk: trim空欄はauto。decimal/exponent以外のJS Number形式を使っていた入力は拒否され得るが、Pythonと一致しない隠れた入力経路を支持しない。scope は listed slot fields に限る。
Owner: satoshikawato
Decision date: 2026-09-26
```

Inert machine representation（既存静的契約への転記用。独立した JSON registry ではない）:

```json
{
  "concern": "tracks.pixel-text-input-domain",
  "scenarioRevision": 1,
  "choice": "A / PIXEL_TEXT_OPTIONAL_PX",
  "rationale": "利用者が pixel 値を単位付きで paste できるようにし、検証・正規化・request の値を揃える。物理 pixel と factor scalar は分けたまま、保存形式を増やさずに入力の一貫性を改善する。",
  "mustPreserve": "既存 valid 数値、Linear px acceptance、auto、physical pixel 意味、現行 typed request/Session、Circular radius/width factor/%、retired key 拒否。 typed JSON の gaps は数値、Linear は既存 ScalarSpec。Circular ratio/% semantics を保つ。不正値を null/0 化しない。不要な arbitrary CSS unit conversion を作らない。",
  "mayRetire": "pure pixel 対象の without-a-unit restriction と、invalid→null/zero の黙示的変換。一般 dimension input の制限は退役しない。",
  "acceptedResidualRisk": "trim空欄はauto。decimal/exponent以外のJS Number形式を使っていた入力は拒否され得るが、Pythonと一致しない隠れた入力経路を支持しない。scope は listed slot fields に限る。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-26"
}
```
