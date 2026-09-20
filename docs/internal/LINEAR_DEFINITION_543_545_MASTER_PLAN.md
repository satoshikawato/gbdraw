# Linear Definition・Replicon表示修正 — 総合計画書

状態: S0・S1・S2完了（2026-09-20 JST）。D1-A・D2-P・D3-Aは契約PR [#546](https://github.com/satoshikawato/gbdraw/pull/546)でdevへ反映済み。修正、実ブラウザのGenerate／Save／Load、文書・図、必要なローカルgateを検証済み。runtime差分は`fix/linear-definition-543-545-20260920`で未commit・未push。対象版、受入結果、制約と再現コマンドは第10節のS1/S2記録にある。

対象Issue:

- [#543 Replicon labels persist despite being turned off](https://github.com/satoshikawato/gbdraw/issues/543)
- [#545 Definition is aligned to the right in linear mode](https://github.com/satoshikawato/gbdraw/issues/545)

実装者向けの開始指示は、別紙の[INSTRUCTION PROMPTS](LINEAR_DEFINITION_543_545_INSTRUCTION_PROMPTS.md)にある。本書は問題、根拠、表示契約、設計、受入条件の管理先とする。チャット履歴や一時ディレクトリ内の監査ファイルを前提にしない。

## 1. 目的と対象

gbdrawは、GenBankまたはGFF3＋FASTAからゲノム図を生成するPythonソフトウェアである。Web版はVueの単一ページアプリで、ブラウザ内のWorkerがPyodideを使って同じPython描画処理を実行する。

Linearモードでは、一つの行に一つ、または複数の生物学的レコードを配置できる。本計画は、レコードに付く文字の表示条件と横配置を修正し、Web・CLI・Pythonで共有する描画処理の責務を明確にする。

| 用語 | 意味 |
| --- | --- |
| レコード | 一本の染色体、プラスミドなどの配列と注釈。入力ファイル一つに複数含まれる場合がある。 |
| 行 | Linear図でレコードを横方向に並べる配置単位。入力ファイルとは別の概念。 |
| Definition | Name、Subtitle、Replicon、Accession、Lengthなどの文字を積み重ねるブロック。 |
| 行Definition | 行全体に共通するNameなどを、行の左側へ一度だけ描く部分。 |
| レコード固有Definition | レコードに固有の文字。複数レコード行では各配列の上に描く。 |
| Subtitle | レコード固有の明示入力、または入力ファイルの既定値から得る副題。 |
| Replicon行 | 生物学的メタデータから得た名称を、Show Repliconの設定に従って描く行。 |
| Lock Definition Column | `form.keep_definition_left_aligned`。Definitionを共通左列に固定する設定。 |
| canonical request | Web状態をPythonへ渡し、Sessionにも使用する共通の型付き描画要求。 |
| Session | 入力、設定、保存済みResultなどを再利用するための保存ファイル。 |

完成条件は、設定から最終SVGまで表示の意味が一致し、文字の配置と衝突判定が同じ幾何情報を使うことである。中央揃えの基準、Organelle、保存済みデータの扱いは第4節の明示決定で確定している。本書は計画と実施記録の管理先であり、正式なdurable authorityはbaseの`PD-OI-024`、実装・検証の対象は第10節の版と差分で識別する。

## 2. 適用する規則

着手時に以下を読む。パスは本書からの相対パスである。

- [AGENTS.md](../../AGENTS.md)、[CLAUDE.md](../../CLAUDE.md)、[Web CLAUDE.md](../../gbdraw/web/CLAUDE.md)
- [Product Impact Ratchet](PRODUCT_IMPACT_RATCHET.md)と[Product Decision Packet Template](PRODUCT_DECISION_PACKET_TEMPLATE.md)
- [Architecture Fitness Function Ratchet](ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md)
- [Session互換性](../SESSION_COMPATIBILITY.md)、[Web仕様](../REFERENCE/web-app.md)、[CLI仕様](../CLI_Reference.md)

初回の実装着手では作業ツリーを確認し、`git fetch origin`の後、最新`origin/dev`から追跡先なしの作業ブランチを作る。既存の無関係な変更とworktreeを保持する。継続セッションは実施記録に記載された同じ作業を引き継ぎ、ブランチの由来と変更内容を確認する。ブランチを作り直す必要がある場合も最新`origin/dev`を起点とし、必要な既存変更を重複なく引き継ぐ。

本計画はpush、PR作成、merge、tag、公開、deployを許可する文書ではない。実装セッションではローカル変更と必要な検証を進め、外部操作はその時点の明示的な許可に従う。

## 3. 監査結果と再現可能な根拠

### 3.1 監査した版と検証の限界

監査対象は `100ac32ad5a73cae22e7b43a2dbbc0bb2126e699`。両Issueの画像の表示版も`v0.14.0+100ac32`である。実装開始時にはこのSHAを最新と仮定せず、最新の作業対象との差分を確認する。

元Issueには入力GenBankとSessionが添付されていないため、元データそのものの再実行ではなく、下記の独立した入力で再現した。Webのアップロード状態、現在のPythonソースで生成したSVG、そのSVGのChromiumでの文字座標を確認している。

監査時のローカルbrowser wheelでは`builders.py`と`assemble.py`が監査対象ソースと一致しなかった。したがって、そのwheelによるPyodideのGenerate全経路を合格の根拠にしていない。このローカルの差異をホスト済みWebの不具合原因とは断定しない。実装完了には、実装対象ソースから用意したwheelでの実ブラウザ生成が必要である。

### 3.2 #543の原因

複数レコードを展開するWeb処理が、推論したReplicon名を`record_subtitle`へ書き込む。その後は通常のSubtitleとしてrequestへ渡る。Python側は非空のSubtitleを表示し、`show_replicon`は別のReplicon行にだけ適用する。

| 境界 | 監査対象コードと確認事項 |
| --- | --- |
| Web入力 | [app-setup.js](../../gbdraw/web/js/app/app-setup.js)の`expandDiscoveredLinearRecords()`が`inferredSubtitleFromReplicon`を条件に`record_subtitle`を設定する。 |
| 明示入力の解決 | [linear-sources.js](../../gbdraw/web/js/app/linear-sources.js)の`resolveLinearRecordEffectiveSubtitle()`がレコード値、ファイル既定値の順に解決する。 |
| request | [session-request.js](../../gbdraw/web/js/services/session-request.js)の`buildRecords()`が`presentation.subtitle`へ投影する。Show Repliconのboolも正しく投影されている。 |
| 描画 | [definition.py](../../gbdraw/render/groups/linear/definition.py)の`DefinitionGroup._build_definition_lines()`がSubtitleとRepliconを別の種類として扱う。 |

これはフラグの受け渡し漏れではなく、自動推論値を明示入力のフィールドへ移したことで表示制御を迂回する不具合である。

再現入力F1を次のように構成する。全レコードを1,000 bpの`ATGC`反復配列、異なるID、`molecule_type=DNA`とし、全長のsource featureを付ける。sourceには共通の`/organism="Aeromonas hydrophila"`、`/strain="A1"`を付ける。手入力のlabel/subtitleは設定しない。

| レコード | 追加source qualifier |
| --- | --- |
| C | chromosome、plasmid、organelle qualifierなし |
| P1 | `/plasmid="p1"` |
| P2 | `/plasmid="p2"` |

必要なら50–950 bpのCDSを各レコードに加える。Biopythonの`SeqIO.write(..., "genbank")`で一つのファイルへ出力する。

監査で確認した結果:

- C、P1、P2を一つのファイルとしてアップロードすると、Show Replicon=falseでもSubtitleは空、`p1`、`p2`となる。
- P1だけを単独ファイルでアップロードすると、Subtitleは空のまま。
- P1、C、P2に並べ替えてアップロードしても、先頭のP1に`p1`が入る。「2番目以降」は必要条件ではない。
- Webが補ったSubtitleを持つレコードをPythonで描くと、falseでも`p1`と`p2`が`data-definition-line-kind="subtitle"`として残る。
- 同じ入力でtrueにすると、それぞれがSubtitleとRepliconの2行になる。

### 3.3 #545の原因

[builders.py](../../gbdraw/diagrams/linear/builders.py)の`add_record_definition_group()`は、複数レコード行の行Definitionについて次を指定する。

```python
text_anchor="start" if keep_definition_left_aligned else "end"
```

falseの横位置も`horizontal_offset + placement.x - definition_gap`という右端基準になる。通常の一レコード／行経路は既定値`middle`を使うため、経路によって文字揃えが異なる。

再現入力F2は2行×2レコードとし、1行目の両レコードにName=`Aeromonas hydrophila`、Subtitle=`A1`、2行目にはName=`Aeromonas sp.`、Subtitle=`B`を明示する。Replicon、Accession、Length、GC、Skew、feature labelsをオフにする。Lock=falseではNameとSubtitleが右寄せ、trueでは左寄せになる。

通常の一レコード／行経路でも中心座標は各Definition自身の幅から決まる。このため、`middle`にするだけでは、異なる長さのNameを共通の中心軸へそろえられない。[assemble.py](../../gbdraw/diagrams/linear/assemble.py)の`_record_collision_bands()`も横位置を独自に計算しており、描画側だけの変更は境界の不一致を生む。

### 3.4 関連性と既存テストの不足

両不具合はPR #541に含まれる変更群に関連するが、直接原因は別々である。

| commit | 関連する変更 |
| --- | --- |
| `643aa450` | メタデータ推論とSubtitleへのコピーを導入。当初は呼出元へ渡らないplasmidフィールドを条件にしていた。 |
| `24df536f` | 行Definitionの分割をLock設定から分離し、falseで`end`を指定した。 |
| `57459b85` | 分割をassemble側へ集約。Replicon判定情報の受け渡しを直し、Subtitleへのコピーが有効になった。 |
| `35558856` | 行全体に共通するName/Subtitleだけを行Definitionへ振り分けるよう変更した。 |

監査では関連Pythonテスト8件とJSテスト3ファイルが通った。これらは正しい修正の証明ではなく、検出範囲の不足を示す。入力推論、文字内容の振り分け、固定列の配置を個別に確認しても、アップロード後の表示スイッチと最終SVG、複数行の共通中心までは確認できていなかった。

## 4. 確定した表示契約とauthorityへの反映

### 4.1 修正で守る要件

- 自動Replicon表示はShow Repliconに従い、オンでは対象レコードに一度だけ描く。
- ユーザーのSubtitleをShow Repliconで一括非表示にせず、文字列一致だけで削除しない。
- 行共通のName/Subtitleは一度だけ表示する。異なる明示値、先頭だけに設定された行識別、後続だけに設定された固有値が欠落しない。
- レコード固有の情報と行全体の情報を、レコード数や先頭位置だけで混同しない。
- 表示行の種類とLine Stylesの対象を一致させる。Repliconのスタイルで自動Replicon行を編集できる。
- 測定、描画、衝突判定、最終キャンバスの境界を一致させる。
- 左側の設定変更は、既存のGenerateによる適用手順で検証する。本修正のために新しい自動再生成やwatcherを追加しない。

### 4.2 採用した結果

以下の選択は第4.3節のProduct Decisionによって確定している。D1–D3は本計画内の管理IDであり、登録済みの`BD-###`ではない。実装者は同じ選択を再質問しない。

| 採用ID | 完成後の表示と適用範囲 |
| --- | --- |
| **D1-A: 共通幅の中央** | Lock=falseでは、同じ開始位置の行のDefinitionが共通幅の中心にそろう。行を移動すると対応するDefinitionも追従する。単一／複数／混在行に適用する。Lock=trueの共通左列は維持する。 |
| **D2-P: 保存値を保持** | 保存済みSubtitleは自動／手入力を推測して削除しない。読み込みだけでは保存Resultを変えず、Generateで新しい表示契約を適用する。不要なSubtitleは利用者が明示的にクリアする。ファイル既定値へ戻る既存の継承規則を維持する。 |
| **D3-A: OrganelleもReplicon行で制御** | chromosome、plasmid、organelle由来の自動名をShow Repliconで制御する。対象名はオンで一つ、オフでゼロ。複数qualifierの候補はchromosome→plasmid→organelleの順で一つを選ぶ。organelleの表記は既存の自動Subtitle表記を引き継ぐ。Web・CLI・Pythonの共通描画に適用する。 |

D1-Aは既定配置の選択である。既にサポートする明示的な`text_anchor`設定は上書きせず、その受入範囲を調査してテストする。明示設定を新しい経路へ拡張することや既存対応を廃止することは、この決定に含めない。

D2-Pでは現行Sessionの保存済みResultと設定ドラフトを区別する。互換readerやmigrationが必要と判断する前に、`main`のfirst-parent履歴またはrelease tagの実在形式とpositive fixtureを確認する。dev上だけの中間形式を新しい互換契約として増やさない。

D3-Aでは[record_metadata.py](../../gbdraw/core/record_metadata.py)の`infer_record_source_metadata()`と`format_inferred_subtitle()`を確認する。現行Pythonの`metadata.replicon`はchromosome/plasmidのみなので、選択されたorganelle表示を既存の生物学的名称解決とReplicon行で実現する。Show Repliconの既定値falseは維持する。自動名のオン／オフは手入力Subtitleの表示を変更しない。

保存済みSessionに格納された自動SubtitleがShow Replicon=falseでも残りうることと、再Generate後の配置が保存済みプレビューと変わりうることは、明示的に受け入れられた制約である。

### 4.3 Product Decision原文と機械可読な転記

Product Decision Ownerが明示した回答を、そのまま記録する。

```text
PRODUCT_DECISION
Concern: linear.definition-display
Scenario revision: 1
Choice: D1-A, D2-P, D3-A
Rationale: 名前の比較をしやすくし、自動の生物学的名称を一つの表示スイッチで制御する。
Must preserve: 手入力Subtitle、保存済みSessionの値とプレビュー、Lock=trueの共通左列、行共通・レコード固有ラベルの区別。
May retire: Replicon/Organelle名のSubtitleへの自動コピー、Lock=falseで各Definition自身の幅に基づく既定横配置。
Accepted residual risk: 保存済みの自動Subtitleはオフでも残りうる。再Generate後の配置は保存済みプレビューと変わりうる。
Owner: satoshikawato
Decision date: 2026-09-19
```

次のJSONは回答のフィールドを転記した確認用の表現である。選択肢一覧だけを配列へ変換し、理由・維持範囲・廃止範囲・リスクの文言は変更していない。CIが読み込むdecision schemaや新しいauthority storeを定義するものではない。

```json
{
  "concern": "linear.definition-display",
  "scenarioRevision": 1,
  "choices": ["D1-A", "D2-P", "D3-A"],
  "rationale": "名前の比較をしやすくし、自動の生物学的名称を一つの表示スイッチで制御する。",
  "mustPreserve": "手入力Subtitle、保存済みSessionの値とプレビュー、Lock=trueの共通左列、行共通・レコード固有ラベルの区別。",
  "mayRetire": "Replicon/Organelle名のSubtitleへの自動コピー、Lock=falseで各Definition自身の幅に基づく既定横配置。",
  "acceptedResidualRisk": "保存済みの自動Subtitleはオフでも残りうる。再Generate後の配置は保存済みプレビューと変わりうる。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-19"
}
```

### 4.4 実装前に残る手続き

ユーザー可視の結果は確定している。S0で残るのは最新baseへの適用確認、必要な再現、互換性の実在根拠、durable authorityへの反映である。人の選択が未解決であるかのように扱わない。

`linear.definition-display`、scenario revision `1`は手動preflightの識別子であり、Product Impact mapへ登録されたことを意味しない。[Product Impact Ratchet](PRODUCT_IMPACT_RATCHET.md)に従い、未登録のmaterial Product変更の正式な反映先として、既存の[Option Integrity Product Contract](OPTION_INTEGRITY_PRODUCT_CONTRACT.md)とそのbase側の状態を確認する。mapped concernへ該当することが判明した場合は既存のauthority経路を使い、並行した権威を作らない。

選択された結果、理由、維持・廃止範囲、受け入れられたリスクだけを反映するauthority-only変更を準備し、生成した正式表現をレビュー可能にする。runtimeは別に扱い、必要なauthorityがそのbaseへ反映されてから実装する。候補のauthorityは同じ候補のruntimeを許可しない。本書の受領記録だけをdurable base authorityと呼ばない。

同じ選択の再承認は不要である。外部へのpushやmergeが必要な段階では、その対象についての明示的な許可を確認する。通常のローカル編集・調査・検証を停止する理由にしない。

最新baseとの衝突や、決定に含まれない新しいユーザー可視差が発見された場合に限り、追加の根拠調査と必要なProduct手続きを行う。その場合も確定済みの結果は維持し、依存する変更だけを保留して独立作業を進める。

## 5. 修正アーキテクチャ

### 5.1 所有者と依存方向

| 責務 | 既存所有者 | 修正後の役割 |
| --- | --- | --- |
| ファイル選択とレコード展開 | `web/js/app/app-setup.js`、`record-discovery.js`、`python-helpers.js` | 入力と選択肢を扱う。自動Replicon名を手入力Subtitleへ昇格させない。 |
| 明示値・ファイル既定値 | `web/js/app/linear-sources.js` | SubtitleとNameの明示値・継承を解決する。 |
| request／Session投影 | `web/js/services/session-request.js` | 同じ型付き表現へ一度だけ投影する。表示内容を重複管理しない。 |
| 生物学的名称 | `core/record_metadata.py` | 選ばれた契約に従って元データから名称を解決する。 |
| 行とレコードの振り分け | `diagrams/linear/assemble.py::_split_definition_line_kinds()` | 行共通／レコード固有の種類を決定する。 |
| Definitionの測定と配置 | 既存のLinear geometry／layout境界 | 表示行、測定幅、anchor、位置から共有する境界を得る。 |
| SVG描画 | `render/groups/linear/definition.py`、`diagrams/linear/builders.py` | 決定された種類・配置を描画する。別の表示方針を作らない。 |

Web状態→canonical request→Pythonの意味解決・配置→SVGという既存の依存方向を維持する。描画処理からVue、入力ファイルカード、Session読込処理へ依存させない。Pyodideは既存Workerのままとし、メインスレッドへPythonを追加しない。

### 5.2 RepliconとSubtitle

`expandDiscoveredLinearRecords()`によるReplicon→`record_subtitle`の自動コピーを廃止し、自動Replicon名を既存Replicon行で描く。D3のorganelle方針を満たしてから、関連経路を収束させる。

コピー廃止後に利用者がなくなる`inferredSubtitleFromReplicon`、対応するdiscovery payload、Web専用のSubtitle推論処理は、参照検索で確認して同じ変更内で削除する。GFF/FASTAのhelper経路も確認する。Organism/strain推論、他の実利用者があるコア関数、必要な互換readerを巻き込まない。

手入力と自動入力を後から区別するためだけの新しいフラグ、Sessionスキーマ、文字列推測処理は追加しない。DOMから重複文字を消す方法も、生成物とcanonical stateが食い違うため採用しない。

### 5.3 横配置と境界

既存の`_LinearRecordDefinitionGeometry`または同じLinear layout境界にある小さなprivate処理を使い、横位置と境界の計算を共有する。`builders.py`と`_record_collision_bands()`が独立して同じ式を持つ状態を解消する。新しい公開APIや別のレイアウト処理系を作る必要はない。

D1-Aの共通幅による中央揃えを、次の座標で実現する。

- `H`: レコード配置領域のキャンバス上の左端。
- `X`: 行先頭レコードの、配置領域内の横offset。
- `G`: Definition領域と行の間に確保するgap。
- `W`: 表示される左側Definition群の実測最大幅。
- `w`: 対象Definitionの実測幅。

| 状態 | anchorと基準座標 | 対象Definitionの横境界 |
| --- | --- | --- |
| Lock=false、D1-A | `middle`、`H + X - G - W/2` | 基準座標の左右に`w/2` |
| Lock=true、現行の共通左列 | `start`、`H - G - W` | 基準座標から右へ`w` |
| 複数レコード行のレコード固有文字 | 配列区間中央の`middle` | その中央の左右に実測幅の半分 |

D1-Aの`G`は予約領域と行のgapであり、短い文字自体と行の距離は`G`より大きくなる。テストとhelp文でこの違いを明記する。既定配置はこの一つの契約で実装する。

行Definitionの垂直位置、レコード固有文字の配列上配置、行振り分けは保つ。各layout計算でoffsetまたは幅が変わったら配置結果を更新し、同じ結果を描画・衝突判定・キャンバス境界へ渡す。最終位置を用いた収まりを検証する。

### 5.4 原則を設計と作業へ適用する

| 原則 | 実装での適用 | ワークフローでの適用 |
| --- | --- | --- |
| SRP | 入力、名称解決、行振り分け、配置、描画に所有者を定める。 | 表示契約の判断と実装案の判断を分ける。 |
| OCP | 既存の行種類と配置結果で必要な挙動を表す。 | 将来機能を想定した拡張基盤を作らない。 |
| LSP | 単一／複数レコードが同じ配置契約を満たす。 | 経路を切り替える再現例で意味の変化を確認する。 |
| ISP | builderへ必要な行種類・位置・境界を渡す。 | セッションごとに必要な入力、成果、未決事項を明記する。 |
| DIP | 描画は解決済み入力とgeometryに依存する。 | 証拠はチャット記憶でなく版と再現手順に結び付ける。 |
| DRY | Repliconの二重経路と横位置の重複式を削る。 | 本書を契約の管理先とし、指示書は本書を参照する。 |
| KISS | 既存境界に小さな変更を加える。 | 同じ修正を三つの連続工程で完了し、独立した管理システムを増やさない。 |
| YAGNI | 未使用経路を同時に削り、汎用policy engineや不要な互換層を追加しない。 | 関係しないGallery刷新、全Web再設計、包括的ベンチマークを必須依存にしない。 |

## 6. 受入条件とテスト配置

全組合せを無条件に増やさず、各境界を識別できる代表例を使う。コードの分岐数ではなく、利用者が得る結果をassertする。

| ID | 入力・操作 | 受入条件 | 主な検証先 |
| --- | --- | --- | --- |
| R1 | F1をWebへupload、off→Generate→on→Generate→off→Generate | 自動名はoffで0行、onでレコードごとにReplicon行1行。自動Subtitleを作らない。 | `tests/web/linear-multi-record.playwright.spec.js` |
| R2 | 単独P1、P1を先頭とした複数入力、並べ替え、ファイル差替え | 位置とレコード数に依存したSubtitle生成がない。差替え後に旧自動名が残らない。 | 同上、`linear-source-defaults.test.mjs` |
| R3 | Replicon名と同一文字列の手入力Subtitle、異なる手入力、ファイル既定値 | Show Repliconで明示値を削らない。手入力とRepliconが同名の2行はユーザー指定として区別する。 | JS、request、Python Definition tests |
| R4 | chromosome、plasmid、organelle、qualifierなし、競合qualifier、GFF＋FASTA | D3の対象・表記・優先順を満たす。対象外の自動名を作らない。 | `test_record_metadata.py`、discovery/helper、browser |
| R5 | Line StylesでRepliconとSubtitleを別々に指定 | 正しい種類にだけfont、weight、colorが適用される。 | `test_definition_line_styles.py`、browser |
| A1 | F2、Nameだけの場合、単一／複数／混在行、Lock両値 | D1の基準で中央または左列をそろえる。文字属性だけでなく表示座標を測る。 | `test_linear_definition_alignment.py`、browser |
| A2 | 長いName、短いSubtitle、行offset、中央／類似性に基づく配列整列 | 非固定は行に追従、固定は共通左列。予約gap、衝突境界、キャンバス内の収まりが一致する。 | alignment／multi-record layout tests、browser |
| A3 | 共通Name、異なるSubtitle、先頭だけの行識別、後続だけの固有Name、空行見出し | 現在の行振り分けで必要な各文字が一度表示され、欠落しない。 | `test_linear_multi_record_layout.py` |
| P1 | 新規SessionをSave→新しいページでLoad→Generate | フラグ、明示Subtitle、配置契約が維持される。保存Resultと再生成Resultを区別する。 | session-request tests、browser |
| P2 | 実在する対応済み旧Session、明示Subtitleクリア→Save→Load | D2を満たし、読み込みだけで保存ResultやSubtitleを破壊しない。クリア後の継承規則も維持する。 | Session互換tests、browser |
| X1 | 同じ表示契約のWeb request、CLI、typed Python入力 | 共通Python描画で行種類と配置が一致する。入力adapterの差を意図せず広げない。 | Python API／CLIの既存test所有者 |
| X2 | 共通state／metadata変更が到達するCircular操作 | Name等の既存結果とrequest契約を維持する。 | 影響する既存Circular test |

ブラウザの文字中心は固定フォントのDOM座標または`getBBox()`＋CTMで比較する。字形のはみ出しを考慮した小さな許容差を事前に定め、失敗を通すために緩めない。解析モデルの境界とSVGの実測値も比較する。

## 7. 実装セッション

| セッション | 入力 | 作業と成果 | 完了条件 |
| --- | --- | --- | --- |
| S0: 基準・authority反映 | 第4節の明示決定、最新base、現行仕様・履歴 | F1/F2の基準確認、互換性と到達範囲の調査、確定したD1-A・D2-P・D3-Aのauthority-only変更準備。実施記録へ追記する。 | 必要なauthority経路が満たされ、runtimeを実装できるbaseが特定されている。外部操作待ちは選択未決と区別する。 |
| S1: 修正とfocused検証 | S0の結果と対象ブランチ | 失敗する回帰例を既存test所有者へ追加し、二重表示経路と重複配置式を解消する。関連focused testsを通す。 | 対象runtimeとtestsが整合し、不要経路を削除済み。S2未実施を明示して引き継ぐ。 |
| S2: 実ブラウザ・文書・図・最終確認 | S1の変更、版と検証記録 | 現在のwheelでR/A/P受入操作、Xの到達範囲、必要な図・文書・必須gateを確認。見つかったin-scope不具合を直す。 | 第9節の条件がすべて成立するか、未完了の境界を具体的に報告する。 |

各セッションは別プロンプトで開始できる。同じ担当者が連続実行してもよい。セッション分割は外部承認の追加を意味しない。通常の編集・build・検証のたびに許可を取り直さない。

## 8. 検証、文書、生成物

### 8.1 focused checksと必須gate

コマンドはリポジトリルートで実行する。変更範囲に合わせて既存test所有者を利用し、消した内部helperを試すだけのtestは整理する。

```bash
pytest tests/test_linear_multi_record_layout.py tests/test_linear_definition_alignment.py tests/test_definition_line_styles.py tests/test_record_metadata.py -q
node --test tests/web/linear-source-defaults.test.mjs tests/web/record-metadata-inference.test.mjs tests/web/session-request.test.mjs
ruff check gbdraw/
```

関連するrequest、Session、CLI、Circular testsを第6節の到達範囲に従って追加する。基準と変更後を比較するWeb policy／architecture checksと、CIのimpact判定が要求するgateは省略しない。既存checkerの失敗を新しい許容設定で回避しない。

```bash
command -v playwright
playwright --version
python -c "from playwright.sync_api import sync_playwright; print('python playwright ok')"
node -e "console.log(require.resolve('@playwright/test'))"
python tools/prepare_browser_wheel.py
npm run test:web:functional-full -- tests/web/linear-multi-record.playwright.spec.js --workers=1 --retries=0
pytest tests/test_output_comparison.py::TestOutputComparison -v
```

browser specは必要な追加ケースと影響する既存ケースを対象に絞ってよい。その場合は選んだテスト名と対象外の理由を記録する。NodeのPlaywrightが使えない場合もPython Playwrightで等価の実操作を確認する。Chromiumが`sandbox_host_linux.cc ... Operation not permitted`で失敗した場合は、同じローカル確認を必要なsandbox権限で再実行する。

テスト全体を30分未満の外部timeoutで打ち切らず、長い実行は途中経過を確認する。テスト自身の短いtimeoutは変更しない。同じソース、入力、環境、受入条件の有効な証拠は再利用する。失敗、変更、未解決事項がある範囲だけを再検証する。

### 8.2 文書と図

設定の意味を変える場合は、既存のWeb仕様、CLI help／reference、Session互換説明の該当箇所を更新する。public仕様文書を分散させず、説明の管理先を維持する。Galleryの画面説明やcaptureを変更する必要が生じたときは、その作業に対応する既存skillを読む。

実際に影響する既存のmulti-record例を選び、記録されたrecipeまたはSessionから描く。候補はVibrioのmulti-record／collinear例である。必要なラベル、legend、色、track、比較情報を保ち、読める倍率で目視する。F1/F2の最小入力を公開の完成例へ転用しない。

`tests/reference_outputs/`は通常検証では読み取り専用。意図した幾何変更をレビューした後に限り、`--update-reference-outputs`で必要な参照を更新し、SVG差分と比較テストを確認する。browser wheelは生成されたgitignored assetとして扱い、手編集・commitしない。配布用bundleを作らない通常検証ではcache-bust更新を目的にしない。`examples/gbdraw_social_preview.png`は変更しない。

## 9. 完了条件とレビュー

- D1–D3について、選択内容、根拠、適用範囲、必要なProduct手続きが記録されている。
- #543と#545の受入条件が実際の生成SVGで成立し、残る既存Subtitleの扱いが説明されている。
- 単一／複数／混在行、明示Subtitle、保存・再読込・再生成が選択契約に一致する。
- 自動Repliconの表示経路が一つになり、消費者のない推論payloadが残っていない。
- 横位置と境界の重複した判断を削除し、測定・描画・衝突判定・キャンバス境界が一致する。
- 変更後のPythonソースとbrowser wheelの対応を確認し、そのwheelでブラウザ生成を検証している。
- focused checks、影響する既存回帰、必要なpolicy／CI gateを通している。未実施のものを合格扱いしない。
- production、tests、docs、生成図の差分をそれぞれ確認し、公開例を再現・目視している。
- 変更ファイル、検証、残る制約、ロールバック方法、英語のproposed commit titleと短いsummaryを記録している。

通常のアーキテクチャ変更では、責務／経路の前後、削除した旧経路、挙動検証、適用gate、rollbackを簡潔に記録する。本計画の想定は「Replicon表示の二重経路を削除」「行振り分け所有者は維持」「横配置結果を共有」である。所有者・経路・互換分岐が増えるなど[例外条件](ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md)に該当するときだけ、完全なOE／PE／CBの集合と必要な例外判断を準備する。

ロールバックは本修正に属するruntime・tests・docs・生成物を一貫した単位で戻せるようにする。保存済みSessionの一括書換えをrollbackの前提にしない。実際のrevert、破壊的な操作、公開は許可された範囲で行う。

## 10. 実施記録

Product Decision Owner `satoshikawato`による`2026-09-19`のD1-A・D2-P・D3-A選択を受領し、第4節へ原文とJSON転記を保存した。正式なauthorityはPR #546で反映済み。S1/S2の修正後の証拠は本節末尾に記録する。第3節とS0の再現は修正前の観測として保持する。

後続セッションはこの節に必要な行を追記し、同じ内容の報告書を量産しない。ログやfixtureは担当する既存test／artifactの保存先を使い、リンクを記載する。チャットや消失する一時パスだけを根拠にしない。

| セッション | base／対象版・ブランチ | 判断・実施範囲 | 変更／再現証拠 | 検証結果・制約 | 次の作業 |
| --- | --- | --- | --- | --- | --- |
| 決定記録 | `100ac32ad5a73cae22e7b43a2dbbc0bb2126e699` / `docs/linear-definition-543-545-plan` | D1-A・D2-P・D3-A受領済み | 第4.3節に原文とJSON転記 | runtime・authority storeは未変更 | 確定した決定をS0へ引き継ぐ |
| S0 | `100ac32ad5a73cae22e7b43a2dbbc0bb2126e699` / `docs/linear-definition-543-545-authority-20260919` | F1/F2再現、main互換履歴、明示anchor範囲を確認。契約1ファイルを許可後にPR #546で反映 | 下記S0記録と再現コード、PRのCI結果 | merge SHA `9789a2cc6a9df131caacee2caa8e88c1c2747ef8`。完了 | なし |
| S1 | `9789a2cc6a9df131caacee2caa8e88c1c2747ef8` / `fix/linear-definition-543-545-20260920` | 自動名の所有者、配置計算を集約。明示値と既存anchor範囲を保持 | 下記S1/S2記録の変更一覧・差分SHA | focused、到達範囲、core、JS回帰を通過。完了 | なし |
| S2 | S1と同じ未commit候補 | 現行wheel、実Generate／Save／Load、文書・図、policy、CI相当のローカル検証 | 下記受入表、実行表、再生成コード | 受入R/A/P/Xを充足。残存失敗なし。remote runtime CIは未実施 | runtimeの外部操作は別途の対象別許可に従う |

未commitの変更を引き継ぐ場合は、HEADだけでなく変更ファイル一覧と差分を識別できる情報を添える。各検証に対象ソース、入力、環境、コマンド、結果を結び付ける。推薦、明示決定、実装済み、検証済みを別々に記録する。


### S0記録（2026-09-19～20 JST）

以下は契約候補の公開許可を受ける前の調査記録。未merge・未実装という記述はこの時点の状態であり、最新状態は末尾のS1/S2記録を参照する。

#### 対象・正式反映先・実装開始条件

- `git fetch origin`後の`origin/dev`は監査版と同じ
  `100ac32ad5a73cae22e7b43a2dbbc0bb2126e699`。追跡先なしの
  `docs/linear-definition-543-545-authority-20260919`をこのSHAから作成した。
  `.git`の書込制約でbranch作成を一度失敗した後、権限拡張で同じ操作を完了した。
- 元Issueは`gh issue view 543/545 --repo satoshikawato/gbdraw --json title,body,url`
  で確認した。本文の添付は画像のみで、GenBank／Sessionはない。F1/F2は独立した合成入力。
- baseのProduct Impact mapにあるconcernは
  `product.canonical-render-request-boundary`と
  `product.saved-session-regeneration-continuity`。いずれも維持対象であり、
  今回の表示選択を所有していない。baseの`tools/web-product-decisions.json`は
  `maintainerLogins=["satoshikawato"]`、`decisions=[]`。
- `linear.definition-display`は未登録のLane B。正式反映先は既存の
  [Option Integrity Product Contract](OPTION_INTEGRITY_PRODUCT_CONTRACT.md#pd-oi-024-linear-definition-alignment-and-automatic-replicon-visibility)。
  `PD-OI-018/019`はレコード範囲・配置・初期値の契約で、中央揃えの基準や
  Organelle表示を選択しない。既存仕様・main履歴にもD1/D3の完全な選択はなく、
  新しい衝突は見つからなかった。
- 手続き分類は`PRODUCT_DECISION_REQUIRED`の受領後段階：人の選択は第4.3節で
  完了し、routeは`DURABLE_AUTHORITY_REQUIRED`。未mergeなので
  `IMPLEMENT_EXISTING_AUTHORITY`ではない。決定に必要な追加証拠や未記入項目はなく、
  `EVIDENCE_REQUIRED`／`NOT_ALLOWED`でもない。決定の再質問は不要。
- 契約revision 7、`PD-OI-024`、scenario revision 1を候補として追加した。
  JSON全フィールドは第4.3節と機械的に同一比較済み。契約ファイルSHA-256は
  `1820649d884e768aab41f7a2ced8eb1f15c335794435674e4f2ffb769e8108ef`。
  これは未commit・未push・未mergeの候補であり、base authorityではない。
- checkerは契約変更を**その1ファイルだけ**に分離する。既存の未追跡`.worktrees/`と
  本計画／指示書を含む作業ツリー全体はGate FAILになるため、それらを消さずに
  最新baseの隔離コピーへ契約1ファイルだけを反映して検証した。
  計画書と指示書をauthority-only commit／PRへ一緒に追加してはいけない。
- S1開始条件は、この契約候補のauthority-only反映が`origin/dev`へmerge済みであること。
  同じ作業ブランチ上のローカルcommitだけでは代用できない。外部操作は未許可のため、
  push・PR・mergeは行っていない。許可後はremoteの実状態とcandidate差分を再確認し、
  契約だけを反映する。S1はそのbaseから開始し、本書の未追跡記録を保持する。

#### F1/F2の観測と維持範囲

再現時のPythonソースとWeb JSはbaseと同一。F1は1,000 bp、source featureのみ、
CDSなし。手入力Name／Subtitleなし、ファイルSubtitle既定値も空。
Webでは実際のfile inputからアップロードし、discovery完了を待った。
`adv.linear_show_replicon=false`を観測した。

| F1の入力順 | 展開後のSubtitle | Python描画off | Python描画on |
| --- | --- | --- | --- |
| C,P1,P2 | 空,p1,p2 | p1/p2が各subtitle 1行 | 各subtitle 1行＋replicon 1行 |
| P1,C,P2 | p1,空,p2 | 先頭p1もsubtitleとして残る | 先頭p1も2行に重複 |
| P1単独 | 空 | 自動名0行 | p1がreplicon 1行 |

Python描画はWebで観測したNameとSubtitleをannotationsへ渡し、現行assemblerを
実行した。これはuploadと共通描画を結んだ基準再現であり、Pyodide Generate全経路の
合格証拠ではない。browser wheelは作成・使用していない。

F2は4レコード、`#1@1,#2@1,#3@2,#4@2`。1行目は
Name=`Aeromonas hydrophila`／Subtitle=`A1`、2行目は
Name=`Aeromonas sp.`／Subtitle=`B`。Replicon／Accession／Length／GC／Skew／Depth／
feature labelsをオフにした。Chromiumで`document.fonts.ready`を待ち、
`getBBox()`の両端を`getCTM()`で変換して測った。以下はSVG内px、3桁丸め。

| F2経路／Lock | 種類・文字 | anchor | 左端 | 中心 | 右端 |
| --- | --- | --- | ---: | ---: | ---: |
| 複数／false | name: Aeromonas hydrophila | end | 17.039 | 137.362 | 257.685 |
| 複数／false | subtitle: A1 | end | 227.696 | 242.367 | 257.039 |
| 複数／false | name: Aeromonas sp. | end | 95.727 | 176.383 | 257.039 |
| 複数／false | subtitle: B | end | 241.039 | 249.039 | 257.039 |
| 複数／true | name: Aeromonas hydrophila | start | 16.050 | 136.373 | 256.695 |
| 複数／true | name: Aeromonas sp. | start | 16.050 | 96.706 | 177.363 |
| 単一／false | name: Aeromonas hydrophila | middle | 16.859 | 137.256 | 257.652 |
| 単一／false | name: Aeromonas sp. | middle | 95.582 | 176.293 | 257.004 |

F2の各textの`x`属性は`0.0`。複数／falseの行groupのtranslate x合計は257、
複数／trueは16。単一／falseは136.9296875と176.29296875。
属性だけでなく、短いSubtitleが右へ寄ることと、長短Nameの中心の不一致を確認した。
raw観測JSONのSHA-256は
`04c5d4657f7e5a11193394c21a3a209887033a944a7224eef50006e329e74b1f`。

明示的な`objects.definition.linear.text_anchor=start/middle/end`も全経路で試した。
通常の単一レコード／行かつLock=falseだけが明示値を消費する。Lock=trueはstart、
複数レコード行の行DefinitionはLock=falseでend／trueでstart、上部のレコード固有行は
middleを指定する。D1は既存の明示値消費範囲を保ち、未対応経路へ拡張しない。

#### Session・到達範囲・S1の準備

- `origin/main`のfirst-parent先頭は
  `f41a12d2c4323d70d6f5edc96d46d18c0084e2c9`。
  その`services/session-request.js`は`presentation.subtitle`の保存・復元を持ち、
  `data/config.toml`はLinearの`text_anchor=middle`、`show_replicon=false`を持つ。
- 正の旧形式fixtureは
  [`BGC0000708-BGC0000713.v40-schema5.json`](../../tests/fixtures/sessions/BGC0000708-BGC0000713.v40-schema5.json)。
  Session 40／request 5、明示Subtitle 5件と保存Result 1件がある。
  main上のbytesと現行fixtureが一致し、first-parentの
  `4e8c93804186f9c4b163b584bd81d759b1b3522d`に存在する。
  `gallery-session-migration.test.mjs`の既存正例も通過した。
  新規のreader／migration／schemaは不要。mainにないファイル既定値の中間形式を
  新しい互換契約にしない。
- `linear-sources.js::resolveLinearRecordEffectiveSubtitle()`は非空のレコード値、
  ファイル既定値の順。空文字は継承復帰であり、既定値がある場合の「完全非表示」ではない。
  同名の明示値を保存する既存JS正例も通過した。D2ではこの区別を維持する。
- 自動コピーの唯一の消費箇所は`app-setup.js::expandDiscoveredLinearRecords()`。
  `record-discovery.js`の`inferredSubtitle`／`inferredSubtitleFromReplicon`と
  `python-helpers.js::list_sequence_records()`の対応payloadをS1で同時に整理する。
  GFF/FASTA helperの`list_gff_fasta_records()`はselector・ID・長さ・topologyのみを返す。
  Organism／strain推論は保持する。
- Pythonの`infer_record_source_metadata()`はrepliconとorganelleを別々に返し、
  Linear Definitionは現状repliconだけを参照する。Circularも両方を使うので、
  `metadata.replicon`へ無条件にorganelleを追加するとCircularで重複するおそれがある。
  D3は既存metadata境界を用い、Circularの既存2フィールドの意味を維持して実現する。
  現行organelleの自動Subtitle表記は`strip().capitalize()`。
- CLIの`--record_subtitle`／`--show_replicon`、typed requestのpresentationと
  definition config、Webのcanonical requestは同じPython描画へ到達する。
  adapter側に別の自動名処理を設けない。
- S1で追加する失敗例は第6節のR1/R2（upload状態とSVG行数）、R4（organelleと
  優先順）、A1/A2（単一／複数／混在の実測中心・gap・offsetと衝突帯）。
  R3/R5/A3は既存の正しい明示値・style・行振り分けを維持する回帰例として扱う。
  P1/P2は保存ResultとGenerate後を区別し、X1/X2は共通Pythonへの到達を確認する。
  既存の通過テストをD1/D3の合格と取り違えない。

#### 実行コマンドと結果

環境はPython 3.13.3、Python Playwright 1.61.0、Chromium headless shell 1228、
Node 26.8.2。Nodeの`@playwright/test`も解決可能だった。

| コマンド | 結果・意味 |
| --- | --- |
| `pytest tests/test_linear_multi_record_layout.py tests/test_linear_definition_alignment.py tests/test_definition_line_styles.py tests/test_record_metadata.py -q` | 82 passed／browser 3 failed（sandbox起動制約）。305.68秒。 |
| `pytest tests/test_linear_definition_alignment.py -q -m browser`（同じテストを権限拡張） | 3 passed、15 deselected、1.60秒。上記の失敗3件を解消。timeoutや受入式は未変更。 |
| `node --test tests/web/linear-source-defaults.test.mjs tests/web/record-metadata-inference.test.mjs tests/web/session-request.test.mjs` | 3ファイル成功。 |
| `node --test tests/web/gallery-session-migration.test.mjs` | 1ファイル成功。 |
| 下記の再現コード（Python Playwrightを権限拡張） | 実upload 3ケース、F1 SVG 6点、F2 SVG 12点を記録。修正前の観測。 |
| `node tools/check-web-change-budget.mjs --base origin/dev`（作業ツリー全体） | Gate FAIL。既存未追跡worktreeと2計画ファイルがauthority単独条件に抵触。 |
| `node tools/check-web-change-budget.mjs --base HEAD`（同じbaseの隔離コピー、契約1ファイルだけ変更） | Gate PASS／Review REQUIRED。baseとcandidateのauthority validationとseparationはVALID。 |
| `git diff --check` | 成功。 |
| 第4.3節JSONと`PD-OI-024`JSONの`json.loads()`後の同値比較 | 全フィールド一致。 |

隔離検証は`git clone --shared --no-checkout . /tmp/gbdraw-linear-definition-authority-s0`、
`git -C /tmp/gbdraw-linear-definition-authority-s0 checkout --detach 100ac32ad5a73cae22e7b43a2dbbc0bb2126e699`
の後、契約文書だけをコピーして実行した。checker自体は変更していない。
checkerの`spawnSync git EPERM`も同じコマンドを権限拡張して解消した。

production差分とtests差分はそれぞれ空であることを確認済み。docsは契約1ファイルの
候補と本書の実施記録を別々に確認した。生成SVGはテスト用のみで、参照SVG・Gallery・
public図・browser wheelを変更していない。runtime owner/path/compatibilityは全て不変で、
完全OE/PE/CB集合が必要な例外条件は発生していない。S1/S2のfocused後検証、wheel生成、
実Generate、保存再読込、文書・公開図更新は未実施。受入R/A/P/Xの最終合格はない。

未commitの差分は契約文書だけがtracked変更、本計画と指示書は開始時からuntracked。
他のworktreeは保全した。rollbackは未公開の`PD-OI-024`追加とrevision metadataだけを
取り消す範囲で足り、Session変換やruntime rollbackを要しない。

Proposed commit title: `docs: record Linear definition display decisions`

Summary: `Record D1-A, D2-P, and D3-A as PD-OI-024 before dependent runtime changes.`

#### S0再現コード

以下を`/tmp/linear_definition_s0.py`へ保存し、リポジトリrootで別端末から
`python -m http.server 8765 --bind 127.0.0.1`を起動した上で、
`PYTHONPATH=. python /tmp/linear_definition_s0.py`を実行する。
出力先`test-results/linear-definition-s0/`はgitignoredの使い捨てテスト生成物。
このコードは修正前の観測器であり、期待する修正結果を合格判定するtestではない。
観測コードSHA-256:
`c0da810f90dc4e6b47c185f5cff014468a37b7fd6aa4428b4b71f09cc4f7e2a5`。

```python
import json
from pathlib import Path
from copy import deepcopy

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from playwright.sync_api import sync_playwright
from gbdraw.api import LinearMultiRecordOptions
from gbdraw.api.diagram import assemble_linear_diagram_from_records
from gbdraw.config.models import GbdrawConfig
from gbdraw.config.toml import load_config_toml

out = Path('test-results/linear-definition-s0')
out.mkdir(parents=True, exist_ok=True)
records = []
for rid, extra in [('C', {}), ('P1', {'plasmid': ['p1']}), ('P2', {'plasmid': ['p2']})]:
    record = SeqRecord(Seq('ATGC' * 250), id=rid, name=rid, description='S0 synthetic record')
    record.annotations['molecule_type'] = 'DNA'
    record.features = [SeqFeature(FeatureLocation(0, 1000), type='source', qualifiers={
        'organism': ['Aeromonas hydrophila'], 'strain': ['A1'], **extra})]
    records.append(record)
for name, indices in [('F1', [0, 1, 2]), ('F1-first-plasmid', [1, 0, 2]), ('F1-single', [1])]:
    SeqIO.write([records[i] for i in indices], out / f'{name}.gb', 'genbank')

def render(name, items, lock=False, show=False, positions=None, anchor='middle'):
    config = load_config_toml('gbdraw.data', 'config.toml')
    config['canvas'].update(show_gc=False, show_skew=False, show_depth=False)
    config['labels']['linear']['scope'] = 'none'
    config['canvas']['linear']['keep_definition_left_aligned'] = lock
    config['objects']['definition']['linear'].update(
        show_replicon=show, show_accession=False, show_length=False, text_anchor=anchor)
    kwargs = {} if positions is None else {'layout': LinearMultiRecordOptions(multi_record_positions=positions)}
    svg = assemble_linear_diagram_from_records(items, cfg=GbdrawConfig.from_dict(config),
        selected_features_set=[], legend='none', **kwargs).tostring()
    (out / f'{name}.svg').write_text(svg)
    return svg

measure = '''async () => {
  await document.fonts.ready;
  return [...document.querySelectorAll('text[data-definition-line-kind]')].map(e => {
    const b = e.getBBox(), m = e.getCTM();
    const p = new DOMPoint(b.x, b.y).matrixTransform(m);
    const q = new DOMPoint(b.x + b.width, b.y + b.height).matrixTransform(m);
    return {kind:e.dataset.definitionLineKind, text:e.textContent,
      anchor:e.getAttribute('text-anchor'), x:e.getAttribute('x'),
      group:e.parentElement.id, transform:e.parentElement.getAttribute('transform'),
      left:p.x, right:q.x, center:(p.x+q.x)/2, top:p.y, bottom:q.y};
  });
}'''
report = {'uploads': {}, 'svg': {}}
with sync_playwright() as pw:
    browser = pw.chromium.launch()
    for name, indices in [('F1', [0, 1, 2]), ('F1-first-plasmid', [1, 0, 2]), ('F1-single', [1])]:
        page = browser.new_page()
        page.goto('http://127.0.0.1:8765/gbdraw/web/index.html')
        page.wait_for_function('Boolean(window.__GBDRAW_APP__)')
        page.evaluate("window.__GBDRAW_APP__.mode = 'linear'")
        page.locator('[data-linear-source-card] input[type=file]').first.set_input_files(out / f'{name}.gb')
        page.wait_for_function('n => window.__GBDRAW_APP__.linearSeqs.length === n && Boolean(window.__GBDRAW_APP__.linearSeqs[0].file_definition)', arg=len(indices))
        state = page.evaluate('''() => ({show:window.__GBDRAW_APP__.adv.linear_show_replicon,
          records:window.__GBDRAW_APP__.linearSeqs.map(s=>({selector:s.region_record_id,
          subtitle:s.record_subtitle,fileSubtitle:s.file_subtitle,label:s.file_definition}))})''')
        report['uploads'][name] = state
        projected = [deepcopy(records[i]) for i in indices]
        for record, value in zip(projected, state['records']):
            record.annotations.update(gbdraw_record_label=value['label'], gbdraw_record_subtitle=value['subtitle'])
        for show in [False, True]:
            key = f'{name}-replicon-{show}'
            svg = render(key, projected, show=show, positions=tuple(f'#{i+1}@1' for i in range(len(projected))))
            page.set_content(svg)
            report['svg'][key] = page.evaluate(measure)
        page.close()
    page = browser.new_page()
    f2 = [deepcopy(records[0]) for _ in range(4)]
    for i, record in enumerate(f2):
        record.id = f'F2_{i+1}'
        record.annotations.update(gbdraw_record_label='Aeromonas hydrophila' if i < 2 else 'Aeromonas sp.',
            gbdraw_record_subtitle='A1' if i < 2 else 'B')
    for mode, positions in [('multi', ('#1@1', '#2@1', '#3@2', '#4@2')), ('single', None)]:
        items = f2 if mode == 'multi' else [f2[0], f2[2]]
        for lock in [False, True]:
            for anchor in ['middle', 'start', 'end']:
                key = f'F2-{mode}-lock-{lock}-anchor-{anchor}'
                page.set_content(render(key, items, lock=lock, positions=positions, anchor=anchor))
                report['svg'][key] = page.evaluate(measure)
    browser.close()
(out / 'observations.json').write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n')
print(json.dumps(report, ensure_ascii=False, indent=2))
```

### S1/S2完了記録（2026-09-20 JST）

#### Authority、対象版と公開範囲

ユーザーは契約文書1ファイルのpush・PR作成・必要なレビューとCI確認後のdevへのmergeを明示許可した。契約commit `823ddfe56d02285ffd1536fc4b620ca31748c8dc`を同名の作業ブランチへpushし、[PR #546](https://github.com/satoshikawato/gbdraw/pull/546)を作成した。差分と決定の転記を確認し、必須CIの成功後、2026-09-19 15:24:49 UTCにmergeした。追加の必須reviewer人数設定はなかった。契約revision 7／`PD-OI-024`の内容とSHA-256は上記S0候補と同一。

runtimeはmerge済みbase `9789a2cc6a9df131caacee2caa8e88c1c2747ef8`から、追跡先なしの`fix/linear-definition-543-545-20260920`を作成して実装した。手続きは`IMPLEMENT_EXISTING_AUTHORITY`。D1/D2/D3を再選択せず、既存のcanonical request、Worker、保存形式を使う。runtimeのcommit・push・PR作成・merge・deployは行っていない。上記の許可は契約候補に限定される。

最終tracked差分は40ファイル、662 additions／512 deletions。`git diff --binary HEAD`のSHA-256は`dfb6ef403ad9d543b198bfc808ed89e00dd046c1f8536a225c739188d9d8b460`。本計画と指示書は開始時から未追跡であり、このSHAには含まない。既存の未追跡`.worktrees/`も保持した。変更一覧は[paths JSON](../../test-results/linear-definition-543-545/linear-definition-paths.json)、再現可能な所有者と変更内容は以下のとおり。

| 対象 | 変更と責務 |
| --- | --- |
| `gbdraw/web/js/app/app-setup.js`、`record-discovery.js`、`python-helpers.js` | Repliconから明示Subtitleへの自動コピーと、そのためだけのpayload・フラグ・推論関数を削除。Organism／strain、GenBankとGFF/FASTAの既存discoveryは維持。 |
| `gbdraw/core/record_metadata.py`、`render/groups/linear/definition.py` | 自動名を`format_replicon_label()`へ集約し、Replicon行だけで表示。chromosome→plasmid→organelleの優先順、organelleの既存capitalize表記を維持。Circularが使うreplicon／organelleの2フィールドは統合しない。 |
| `gbdraw/layout/linear.py` | `place_linear_definition()`の小さな配置結果にx・left・rightを集約。Lock ONは共通左列、OFFの既定値は共通幅の中央＋行offset、レコード固有行は配列中央。 |
| `gbdraw/diagrams/linear/{precalc,assemble,builders,positioning}.py` | 一度の事前計測から正確な個別幅と共通予約幅を得る。描画と衝突帯が同じ配置結果を使う。重複計測・未使用positioning関数を削除。行振り分けは既存`_split_definition_line_kinds()`に保持。 |
| `gbdraw/web/js/services/session-resources.js` | Save時のbindingに既存readerが受理する`file_definition`・`file_subtitle`を保存。明示クリア→ファイル既定値へ継承→Save→Loadの失敗で、writerの欠落2フィールドを特定した。schema・migration・新しい継承規則は追加しない。 |
| `gbdraw/linear.py`、`gbdraw/web/index.html` | Lockとgap、自動Repliconの対象・Subtitleとの独立性についてhelpを更新。 |
| `tests/test_{linear_definition_alignment,definition_line_styles,record_metadata}.py`、`tests/web/{linear-multi-record.playwright.spec,record-metadata-inference.test,session-resources.test}.js/mjs` | 最終SVG、実Generate、Save／Load、旧Session、明示値、style、配置境界の回帰を既存test所有者へ追加。metadata推論fixtureから未使用Subtitleケースを削除。 |
| `tests/fixtures/public_contract.json` | CLI helpの2説明文に対応するLinear actionsのhashだけを更新。API、default、代表描画hashは不変。 |
| `docs/CLI_Reference.md`、`docs/REFERENCE/{web-app,session-and-request-compatibility}.md` | 上記の表示と保存互換性を既存の説明へ反映。 |
| `docs/images/`10 SVG、`tests/reference_outputs/`6 Linear SVG | 実測幅を共有する意図した幾何変更に合わせ、既存recipeと正式な参照更新コマンドで再生成。下記の構造比較・目視・再実行で検証。 |

通常の一レコード／行・Lock OFFにおける明示`objects.definition.linear.text_anchor=start/end`の消費範囲は維持する。複数／混在行や配列上のレコード固有Definitionへ、この設定の意味を拡張していない。Show Repliconの既定値はfalseのまま。保存済みSubtitleが自動名と同じでも削除しない。

#### 失敗例から修正後の受入まで

最初に回帰例を追加し、修正前のPythonでは9件失敗（共通中央6、organelle3）、36件中27件成功だった。実ブラウザでもupload直後のSubtitleが`["p1", "", "p2", "Plastid:chloroplast"]`となる失敗を確認した。修正後は自動Subtitleが空で、Show Repliconに従う独立の行として生成される。

S2ではファイル既定値を設定して明示Subtitleを消去した後のSave／Loadが失敗し、上記writerの欠落を修正した。旧Sessionの比較は既存Result ingestionがXML宣言やeditor属性を正規化するため、全図形・文字・幾何を比較し、LoadだけではGenerate実行数が0であることも確認した。受入の「保存プレビューを変えない」をバイト列の同一性と混同しない。browser wheelのHEAD応答は空なので、実GET応答だけをhash対象にした。テスト所有のtimeoutと1 pxの文字実測許容差は緩めていない。

| ID | 修正後の証拠と結果 |
| --- | --- |
| R1 | `Linear automatic replicon names follow Generate and preserve saved subtitles`。P1,C,P2,Oを実uploadし、off/on/offそれぞれの実Pyodide Generateで行種類・文字・個数を確認。成功。 |
| R2 | 同testの並べ替え、先頭P1、単独P1へのファイル差替えで自動Subtitleがないことを確認。既存Automatic Linearとdiscovery JSも成功。 |
| R3 | 同名p1と異なる手入力Subtitle、ファイル既定値を維持。Show OFFでも明示値を表示し、ONでは別種類の2行を許す。Python／JS／実ブラウザで成功。 |
| R4 | Pythonのqualifier優先順・organelle・qualifierなし、browserのorganelle、GFF/FASTAのCLIとtyped renderで成功。GFFは既存adapterのfeature選択を保ち、sourceを明示選択した入力で検証。 |
| R5 | SubtitleとRepliconの異なるfont／weight／colorをPythonと実ブラウザで確認。成功。 |
| A1 | 単一／複数／混在行、Nameのみ／Name＋Subtitle、Lock両値をSVGとChromium実測で確認。UIの`Linear Lock Definition Column applies common centers and left edges after Generate`も成功。ONでNameとSubtitleの共通左端、OFFで行offsetを引いた共通中央。 |
| A2 | 不等長配列、align_center両値、既存類似性整列、gap・canvas境界・最終collision帯を確認。新規browser 12ケースと既存gap 3ケースが成功。`getBBox()`＋CTMをSVG座標へ戻し、1 px許容差で比較。 |
| A3 | 既存`test_linear_multi_record_layout.py`の共通／固有／空見出し／先頭・後続の振り分けを維持。focusedとcoreで成功。 |
| P1 | Show値、同名の明示Subtitle、Line Styles、Lock両値をSave→新しいページ→Load→Generateで確認。保存時のpreviewと未生成のdraftフラグを別々に保持。成功。 |
| P2 | 実在v40/schema5 fixtureのSubtitle5件と保存図形をLoad時に保持し、Generate後も明示値を維持。明示クリア後のファイル既定値もSave／Loadを通過。成功。 |
| X1 | GenBank／GFF、CLI／typed Python／Web requestの共通描画、API・request・Sessionの574件、公開契約、文書recipeを確認。成功。 |
| X2 | Circular metadata／bounds／request、共有Canvasの実ブラウザ、Circular参照10件と関連recipeを確認。既存結果を維持。 |

#### 実行環境、wheelと検証コマンド

Python 3.13.3、pytest 9.0.2、Node 26.8.2、Python Playwright 1.61.0、Node `@playwright/test`とChromiumを使用。両方のPlaywright導入経路を確認した。sandboxによるChromium起動拒否／Node `spawnSync git EPERM`は同じローカルcheckを必要な権限で実行して解消した。以下はremote CIの代用結果を偽装するものではなく、同じ必要範囲のローカル検証記録である。

`python tools/prepare_browser_wheel.py`で生成した`gbdraw/web/gbdraw-0.14.0-py3-none-any.whl`のSHA-256は`5ae81e16570192dd0ec024c6637b6c7433dcc016462d6003d9d960a2f713211f`。変更したPython8モジュールについてwheel内bytesとソースbytesを比較し、全一致した。実Generateのブラウザが要求したwheelのGETを`route.fetch()`で取得し、同じ応答をそのまま返した上でSHA一致をassertした。描画のmock置換はない。cache-bustは配布準備ではないため更新していない。wheelはgitignoredで、commit対象外。

[wheel identity](../../test-results/linear-definition-543-545/linear-definition-wheel-identity.json)に各モジュールのSHAを保存した。以下でも同じ対応を再確認できる。

```python
from pathlib import Path
from zipfile import ZipFile
import hashlib
import subprocess
wheel = Path("gbdraw/web/gbdraw-0.14.0-py3-none-any.whl")
assert hashlib.sha256(wheel.read_bytes()).hexdigest() == "5ae81e16570192dd0ec024c6637b6c7433dcc016462d6003d9d960a2f713211f"
paths = subprocess.check_output(["git", "diff", "--name-only", "HEAD"], text=True).splitlines()
with ZipFile(wheel) as archive:
    for path in paths:
        if path.startswith("gbdraw/") and path.endswith(".py"):
            assert archive.read(path) == Path(path).read_bytes(), path
```

ログは[test-results/linear-definition-543-545](../../test-results/linear-definition-543-545/)へ保存した。ログは補助artifactであり、本書のコマンドと結果、既存のtest／recipeから再現できる。

| コマンド（repo root、特記以外） | 結果／ログ名（`linear-definition-`接頭辞を省略） |
| --- | --- |
| `pytest tests/test_linear_multi_record_layout.py tests/test_linear_definition_alignment.py tests/test_definition_line_styles.py tests/test_record_metadata.py -q -m 'not browser'` | focused 110 passed、3 deselected（追加browser matrix前）。`focused.log`等。最終追加分は下記core／browserで確認。 |
| `pytest tests/test_api_request_render.py tests/test_api_requests.py tests/test_web_request_render.py tests/test_session_request_codec.py tests/test_api_session.py tests/test_session_io.py tests/test_session_compat.py tests/test_circular_definition_bounds.py tests/test_linear_vertical_layout.py -q` | 574 passed。`integration.log`。 |
| `pytest tests/ -m 'not slow and not (recipe or gallery or browser)' -n 4 --dist loadfile --durations=20 -q` | 5733 passed、1 failed、17 skipped。失敗は意図したCLI help変更によるsnapshotだけ。`core.log`。 |
| `pytest tests/test_public_contract.py -q` | help以外のhashが不変であることを確認してfixture更新後、2 passed。上記core失敗を解消。`public-contract.log`。未変更の5733件を無条件に再実行していない。 |
| `node --test --test-concurrency=4 tests/web/*.test.mjs` | 606 passed、0 failed／skipped。`web-contracts.log`。 |
| `pytest tests/ -m 'browser and not slow' --durations=20 -q` | 37 passed、6052 deselected。Node comparison-contractsを起動する既存wrapperも含む。`python-browser.log`。 |
| `npm run test:web:pr-smoke -- --workers=1 --retries=0` | 12 passed。`pr-smoke.log`。 |
| `npx playwright test --config=playwright.functional.config.js tests/web/linear-multi-record.playwright.spec.js --workers=1 --retries=0 --grep 'Linear automatic replicon names'` | 最終wheel GET一致を含め1 passed。`browser-identity-final.log`。 |
| 同コマンドの`--grep 'Released Linear session'`相当（新規2件の選択実行内） | 旧Session 1 passed。`browser-new-final.log`。同実行のもう1件のHEAD hash失敗は上記GET対象化で解消。 |
| 同コマンドの`--grep 'Linear Lock Definition Column'` | Lock ON／OFFと再生成1 passed。`browser-columns.log`。 |
| 同コマンドの`--grep 'Automatic Linear renders every record\|Automatic Linear per-record rows\|GFF annotation targets follow FASTA\|multi-record defaults render a shared Circular\|File-level default organism and subtitle\|GenBank file upload infers'`（正規表現のOR） | 既存6件passed。`browser-final.log`。同ログの新規2件の初回失敗は上記で解消。 |
| `pytest tests/test_output_comparison.py::TestGenerateReferences --update-reference-outputs -k linear -v` | 目視した意図的な幾何変更のLinear 6件を正式更新。`reference-update.log`。 |
| `pytest tests/test_output_comparison.py::TestOutputComparison -v` | Linear 6＋Circular 10＝16 passed。`reference-final.log`。 |
| `pytest tests/ -m 'not slow and (recipe or gallery)' --durations=20 -q` | 282 passed、8 failed。失敗は今回の幅変更で古くなった文書SVG。`recipes-gallery.log`。 |
| `pytest tests/test_cli_comparison_how_to_recipe_contracts.py tests/test_cli_how_to_recipe_contracts.py tests/test_cli_tables_tracks_sessions_exports_recipe_contracts.py tests/test_onboarding_recipe_contracts.py tests/test_python_howto_recipe_contracts.py tests/test_python_tutorial_recipe_contracts.py -q --durations=10` | 図更新後56 passed、1 failed。順次検証するH-PY-04で追加の古いSVGを検出。`recipes-final.log`。 |
| `pytest tests/test_python_howto_recipe_contracts.py::test_python_evidence_recipes_regenerate_from_a_clean_external_context -q` | H-PY-04更新後1 passed。最初の8件をすべて解消し、recipe／galleryの未解決失敗は0。`recipes-hpy-final.log`。 |
| `ruff check gbdraw/`、`git diff --check` | 成功。`ruff-final.log`。 |
| `node tools/check-web-change-budget.mjs --base 9789a2cc6a9df131caacee2caa8e88c1c2747ef8`（下記隔離候補） | Gate PASS、Review REQUIRED。参照SVGとexport削除のレビューは下記で記録。`policy-final.log`。 |

上記のgrep ORはシェルへ渡すとき`|`を使う（Markdown表内のエスケープを外す）。最終候補の残存test失敗はない。coreの17 skippedは成功件数に含めない。slow、全functional／perf／release用ブラウザ、remote runtime CIはこのPR相当scopeでは実行していない。契約PRのremote CI成功をruntime候補へ流用していない。

#### Policyと差分レビュー

`tools/ci-impact-policy.mjs`をbaseのものから使用し、最終40 tracked pathを分類した。`documentation, tests-only, python-core, renderer, web-runtime, session-persistence`、valid=true。architecture変更としてPR full範囲の`web-change-budget, core-pr, recipes-standard, gallery, lint, web-contracts-pr, web-pr-smoke`が必要であり、上記のローカル実行で確認した。[判定JSON](../../test-results/linear-definition-543-545/linear-definition-impact.json)に保存。

既存の未追跡worktreeを候補へ混入させないため、同じHEADのローカルshared clone `/tmp/gbdraw-linear-definition-runtime-review`の作業ブランチへ40 tracked pathだけをコピーした。checkerは変更していない。本書・指示書は未追跡の作業記録として判定scope外。再現は`git clone --shared . <隔離先>`後に同じbaseからbranchを作り、`git diff --name-only HEAD`にある変更ファイルを同じ相対パスへコピーして上記checkerを実行する。

production、tests、docs、生成SVGを別々にレビューした。Web inventoryはexports 882→881（未使用`formatInferredSubtitle`の削除）、canonical owner 2、authority location 2、watcher 52、import graph 151 modules／472 edges／0 cyclesを維持。privilege entries 48／54、依存関係、Worker、schema、互換branchの増加なし。Pythonは自動名と横配置の所有者へ既存責務を集約し、旧コピー・payload・計測・positioningを同時削除した。通常の非増加変更であり、完全OE/PE/CB集合を要する例外はない。

registered Product Impactの2 concernはUNCHANGED。Lane Bの表示差はbaseの`PD-OI-024`に従う。writer2フィールドの復元は既存readerと明示クリア時の継承契約を満たす修正であり、新しい保存形式／別のProduct optionではない。異なる明示Subtitle、行識別、レコード固有文字、saved previewという独立の寄与は受入表でそれぞれ確認した。

#### 文書図・実用図と再生成

`love-me-love-my-docs`を文書図の既存recipe再実行へ適用した。ページ判断はS2の明示範囲どおり、CLI Reference、Web Reference、Session互換説明と既存CLI／Python Tutorial・How-toをkeep。新規ページ、統合、削除は0。新しいGUI screenshot／captureは0で、Gallery文書・captureは変更していない。

再生成コマンドは以下。各runnerは既存manifestから入力をそろえたclean temporary directoryで実際のCLIコマンドまたは公開Python literal snippetを実行する。PNGはSVGの目視用artifactだけで、公開SVGを手編集していない。

```bash
for scenario in H-CLI-01 H-CLI-02 H-CLI-05 H-CLI-06 T-CLI-02; do
  python docs/recipes/run_cli_scenarios.py --scenario "$scenario"
done
for scenario in H-PY-02 H-PY-04 T-PY-03 T-PY-04; do
  python docs/recipes/run_python_scenarios.py --scenario "$scenario"
done
```

更新した公開SVGは`h-cli-01/{lambda_genbank,lambda_gff3}.svg`、`h-cli-02/comparison_table.svg`、`h-cli-05/linear_precomputed_comparison.svg`、`h-cli-06/cli_losatp_pairwise.svg`、`h-py-02/python_linear_comparison.svg`、`h-py-04/python_gff3.svg`、`t-cli-02/lambda_linear.svg`、`t-py-03/python_lambda_linear.svg`、`t-py-04/python_lambda_de3_losatn.svg`。別のCircular生成物は変化しなかった。6参照SVGを含む16点のXMLについて要素数、tag、text、tailが不変で、変化した属性は`transform, data-gbdraw-composition, width, viewBox`だけ。path形状、feature、font、色、label内容は不変。共通幅のceilと実幅のずれを除いたことで、既存例では主に0.254 px（Hepatoplasma例0.476 px）の水平補正が生じた。

[構造比較JSON](../../test-results/linear-definition-543-545/linear-definition-generated-review.json)、[図確認1](../../test-results/linear-definition-543-545/linear-definition-doc-review-1.png)、[図確認2](../../test-results/linear-definition-543-545/linear-definition-doc-review-2.png)、[H-PY-04](../../test-results/linear-definition-543-545/linear-definition-doc-hpy04-review.png)、[BLAST参照図](../../test-results/linear-definition-543-545/linear-definition-reference-review.png)で結果を保持。PNGを読み込んで文字、凡例、feature、比較リンク、canvas端を目視した。参照SVGの通常比較も16件成功。

実用的な複数レコード図は既存Vibrio Gallery Session（11 records／5 rows）からLockだけを変えて生成した。Name、Subtitle、色規則、legend、定量track、比較情報を保持し、2400 pxの白背景PNGでONの共通左端とOFFの共通中央を確認した。[ON図](../../test-results/linear-definition-543-545/showcase/vibrio-locked-white.png)、[OFF図](../../test-results/linear-definition-543-545/showcase/vibrio-centered-white.png)。元のGallery assetやsocial previewは変更していない。これは内部受入の再現用Sessionであり、公開Tutorialの読者入力をSessionへ置き換えるものではない。

再現コード（repo rootで実行）:

```python
from pathlib import Path
from dataclasses import replace
import cairosvg
from gbdraw.session import materialize_session, session_to_request, with_request_output
from gbdraw.api.session_compat import render_session_compatible_request

out = Path("test-results/linear-definition-543-545/showcase")
out.mkdir(parents=True, exist_ok=True)
source = Path("gbdraw/web/gallery/sessions/vibrio-harveyi-group-collinear.gbdraw-session.json.gz")
with materialize_session(source, output_directory=out) as materialized:
    request = session_to_request(materialized)
    for locked in (True, False):
        name = "vibrio-locked" if locked else "vibrio-centered"
        options = replace(request.options, config_overrides={
            **dict(request.options.config_overrides or {}),
            "canvas.linear.keep_definition_left_aligned": locked,
        })
        candidate = with_request_output(
            replace(request, options=options), output_prefix=name,
            output_directory=out, formats=("svg",), overwrite=True,
        )
        result = render_session_compatible_request(candidate, materialized.document.to_dict())
        svg = result.drawing.tostring()
        (out / f"{name}.svg").write_text(svg)
        cairosvg.svg2png(bytestring=svg.encode(), output_width=2400,
                        background_color="white", write_to=str(out / f"{name}-white.png"))
```

文書recipeの入力、literal snippet、artifact参照は既存contract testsで確認した。新しいcapture selectorや認証状態はない。既存recipe／Gallery CIがこの再現を所有しているため、別のharnessやCIを追加していない。

Docs Progress:
- [x] Step 1: Frame — 既存英語CLI／Python文書と3つのReferenceを対象化
- [x] Step 2: Flow census — 既存manifest／ページ所有者をkeep、追加ページなし
- [x] Step 3: Demo data — 下記の公式source比較、既存宣言入力
- [x] Step 4: Smoke proof — 既存recipeの実行と生成図の表示を確認
- [x] Step 5: Execution harness — 既存CLI／Python runnerを使用
- [x] Step 6: Evidence run — 必要な生成物を更新、失敗を修正して再実行
- [x] Step 7: Manual written — 既存help／reference／互換説明を更新
- [x] Step 8: Verify + report — SVG構造・目視・recipe contracts、再生成コマンドを記録

#### 公開入力mirrorのsource比較

再生成対象のLambda、DE3、MIBiG5件は公式sourceを再取得して比較した。UTC 2026-09-19に実行。MIBiGは`.5`のannotation release、NCBIはversioned accessionの`gbwithparts`、retmode=textを指定した。6件は直接byte一致。Lambdaだけは現行taxonomy行の改訂があり、元source hashを残して次の一行だけを既存mirrorの表記へ正規化した結果、全bytesが一致した。配列と全featureのbytesは不変である。旧taxonomyを現行分類とする主張ではない。読者の入力元は引き続き各ページの公式databaseで、リポジトリmirrorはoffline automation用に限定する。

```python
source_line = b"            Caudoviricetes; Zimmerviridae; Jacobvirinae; Lambdavirus;\n"
mirror_line = b"            Caudoviricetes; Caudoviricetes incertae sedis; Lambdavirus;\n"
assert source.count(source_line) == 1
normalized = source.replace(source_line, mirror_line)
assert normalized == mirror
```

以下の表と再検証コードは補助JSONがなくてもsourceとmirrorを照合できる。将来sourceが変わってhashが違う場合、過去の結果を流用せず変更を調査する。

| Fixture／format | UTC取得時刻 | Source URL | Source bytes／SHA-256 | Mirror bytes／SHA-256 | 比較 |
| --- | --- | --- | --- | --- | --- |
| `lambda-genbank`／GenBank | 2026-09-19T15:56:44.075839+00:00 | [公式 NC_001416.1](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_001416.1&rettype=gbwithparts&retmode=text) | 176721／`3c624302adeeb3c00649f549903ab781b9e75bab16069ae655833d536407367f` | 176723／`4b76b8bacc8026aac3f19a4a915f4ac772ad61e7ec18f0e2cc859229f95a66e7` | 直接不一致、上記正規化後一致 |
| `de3-genbank`／GenBank | 2026-09-19T15:56:46.047324+00:00 | [公式 NC_042057.1](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_042057.1&rettype=gbwithparts&retmode=text) | 111686／`288eb87480f8fe6eab6246fe1fc4af78c85a6cb7e591c47fd7d7f0170c932e09` | 111686／`288eb87480f8fe6eab6246fe1fc4af78c85a6cb7e591c47fd7d7f0170c932e09` | 直接一致 |
| `bgc-0000708-genbank`／GenBank | 2026-09-19T15:56:47.446322+00:00 | [公式 BGC0000708](https://mibig.secondarymetabolites.org/repository/BGC0000708.5/BGC0000708.gbk) | 105185／`9a5f971c5ed8c406b20574fb50aac567609deb787eb1e8d4635050aa264a04b0` | 105185／`9a5f971c5ed8c406b20574fb50aac567609deb787eb1e8d4635050aa264a04b0` | 直接一致 |
| `bgc-0000709-genbank`／GenBank | 2026-09-19T15:56:50.164867+00:00 | [公式 BGC0000709](https://mibig.secondarymetabolites.org/repository/BGC0000709.5/BGC0000709.gbk) | 121511／`4b66b7e4b78d429d12176e1e36d0e48178c562a9d128d4308b38753af9995255` | 121511／`4b66b7e4b78d429d12176e1e36d0e48178c562a9d128d4308b38753af9995255` | 直接一致 |
| `bgc-0000711-genbank`／GenBank | 2026-09-19T15:56:53.064853+00:00 | [公式 BGC0000711](https://mibig.secondarymetabolites.org/repository/BGC0000711.5/BGC0000711.gbk) | 72291／`32393648f6a91166444331b83687f1b9b7b24c60553a7ddcb677dfe207736789` | 72291／`32393648f6a91166444331b83687f1b9b7b24c60553a7ddcb677dfe207736789` | 直接一致 |
| `bgc-0000712-genbank`／GenBank | 2026-09-19T15:56:55.540633+00:00 | [公式 BGC0000712](https://mibig.secondarymetabolites.org/repository/BGC0000712.5/BGC0000712.gbk) | 134734／`705104a0daa5c44981b0a1e5352d3e56f012dd2e3ae94c98c85cd0ee9198bf94` | 134734／`705104a0daa5c44981b0a1e5352d3e56f012dd2e3ae94c98c85cd0ee9198bf94` | 直接一致 |
| `bgc-0000713-genbank`／GenBank | 2026-09-19T15:56:58.036847+00:00 | [公式 BGC0000713](https://mibig.secondarymetabolites.org/repository/BGC0000713.5/BGC0000713.gbk) | 79197／`bf182663de453f4a3fc30ed0aa8f040a164eeab1c98e604983844994996e58fb` | 79197／`bf182663de453f4a3fc30ed0aa8f040a164eeab1c98e604983844994996e58fb` | 直接一致 |

取得・比較の再実行（repo root、ネットワーク使用）:

```python
from pathlib import Path
import hashlib
import json
import time
import urllib.request

manifest = json.loads(Path("gbdraw/web/tutorial-data/manifest.json").read_text())
keys = ["lambda-genbank", "de3-genbank"] + [
    key for key in manifest["files"] if key.startswith("bgc-") and key.endswith("genbank")
]
for key in keys:
    item = manifest["files"][key]
    accession = item["records"][0]["id"]
    url = item["provenance"]["sourceUrl"]
    if key in ("lambda-genbank", "de3-genbank"):
        url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
               f"?db=nuccore&id={accession}&rettype=gbwithparts&retmode=text")
    with urllib.request.urlopen(url, timeout=60) as response:
        source = response.read()
    mirror = Path("gbdraw/web/tutorial-data", item["relativePath"]).read_bytes()
    normalized = source
    if key == "lambda-genbank":
        source_line = b"            Caudoviricetes; Zimmerviridae; Jacobvirinae; Lambdavirus;\n"
        mirror_line = b"            Caudoviricetes; Caudoviricetes incertae sedis; Lambdavirus;\n"
        assert source.count(source_line) == 1
        normalized = source.replace(source_line, mirror_line)
    assert normalized == mirror, key
    print(key, url, len(source), hashlib.sha256(source).hexdigest(),
          len(mirror), hashlib.sha256(mirror).hexdigest(), source == mirror)
    time.sleep(0.4)
```

GFF3／FASTAは既存`tools/build_lambda_gff3_fixture.py`によるLambda派生物で、recipeの入力hash検証を維持した。source比較の詳細JSONも[artifact](../../test-results/linear-definition-543-545/linear-definition-source-verification.json)へ保存した。公開入力やannotationを別の版に置き換えていない。

#### 完了判定・制約・rollback

第9節のローカル完了条件はすべて成立。未解決のProduct判断、in-scopeのtest失敗、未確認のwheel差異はない。Lock ONは行共通Name／Subtitleの左揃えで、レコード固有行は従来どおり配列中央。保存済みの自動由来Subtitleは明示値と識別できないため保持する。消す場合はユーザーがSubtitleを明示クリアするが、ファイル既定値がある場合はその値へ継承する。

rollbackはこの候補の40 tracked pathをbaseへ戻す（将来commit後ならそのruntime commitをrevertする）範囲。作業ツリー全体のreset、無関係worktreeの削除、保存Sessionの変換は不要。元の描画へ戻す場合は対応する文書図／参照SVGも同じ変更単位で戻す。merge済みの契約PRは独立のauthorityであり、runtime rollbackに伴って自動撤回しない。

Proposed commit title: `Fix Linear definition alignment and automatic Replicon visibility`

Summary: `Keep explicit and saved subtitles, share Linear definition geometry, preserve file defaults through sessions, and verify rendering with browser, CLI, Python, and documentation regressions.`
