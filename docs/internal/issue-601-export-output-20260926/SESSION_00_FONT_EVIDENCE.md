# INSTRUCTION PROMPT — S00 日中fontとPDF実現性の検証

あなたはgbdraw Issue #601のPDF実現性検証担当です。日本語と中国語（簡体字・繁体字）のPDF出力に対応し、既存fontで表せるアルファベット・数字・記号とstyleを維持して不足部分だけ対応fontへ切り替える方針が、`satoshikawato`に承認されています。文字を保持したPDFを提供するため、font asset、配置、原文抽出と配布コストの証拠を取得してください。

## branch取得と必読資料

SESSION_IDは`s00`です。[総合計画書](./MASTER_PLAN.md)の「セッション共通の取得・終了手順」に従って、remoteから**`fix/issue-601-export-output-20260926`**を取得し、未使用pathの専用cloneで作業してください。他sessionのcheckoutを使わず、このbranchの最新push済みcommitから始めます。

取得後、AGENTS.md、CLAUDE.md、Web CLAUDE、総合計画書、[Decision Pack 01](./DECISION_PACK_01_PDF_JA_ZH.md)、[修正前の証拠](./CURRENT_BEHAVIOR_EVIDENCE.md)、[承認の機械表現](./APPROVED_DECISIONS.md)を読みます。

## 作業所有範囲

この計画directoryの`evidence/`に再現script・入力・結果を保存し、`SESSION_00_RESULT.md`を作成してください。font原本やTTF変換output、PDF、renderer出力は専用temporary directoryへ置き、取得元・版・checksum・生成commandsで再生成できるようにします。必要な証拠artifactは適切な保存先とdigestを記録します。production source、正式Product authority、expected reference SVGを変更しません。

## 実施内容

1. source基点、既存jsPDF/svg2pdf版、font loader、source/wheel配布経路を確認します。現在のexport pipeline、`preparePdfFonts`、`read_pdf_font`、helper operation境界を図示し、同じowner/pathを使う方針を具体化します。
2. 日本語・簡体中国語・繁体中国語のstatic PDF-compatible font候補を比較します。版・出所・license・checksum、normal/bold、既存italic混在とfamily対応、cmap、地域字形、parse結果と容量を表にします。原本と必要representationを一つの有限manifestから生成できる案を確定します。現行JP subsetだけで簡体字を保証しません。
3. disposable browser harnessで現在のexport librariesへfontを登録し、文字列の不足部分だけfontを切り替えます。製品runtimeを編集せず、自動run選択を検証可能な内部probeで試します。probe source、commands、environmentを保存し、一時scriptの存在だけに依存させません。
4. `プラスミドA`、`日本語`、`简体中文`、`繁體中文`、`ABC プラスミドA 简体中文 繁體中文 β ≥ 95%`を、normal/boldと既存Latin italic/bolditalicで検査します。tspan、anchor、baseline、spacing、rotation、textPathとstyle継承も含めます。必要fallbackにstyleがない場合は明示非対応として記録し、無告知normal化をしません。
5. 明示`lang`ありの地域字形と、language情報がない場合の固定順・coverage選択を比較します。code unitで文字を割らず、結合文字、異体字、補助平面`𠮷`は表示・原文抽出の成功または明示的な非対応として記録します。一般的な日中対応を確認せず非対応宣言だけで合格にしません。
6. 実PDF bytesの原文text、font埋込み、image代用なし、page寸法を検査します。既存`pdf-text.cjs`だけで全glyph配置が正しいとは判断せず、別のPDF rendererをdisposable環境へ用意してbrowser referenceと字形・anchor・run配置を比較します。rendererと版、render commands、比較図と目視所見を記録します。
7. 実用的なGallery図でもfont-matched DOMとPDFを比較します。小さいprobeは内部fixtureです。asset payload、初回load/parse、反復export、PDF glyph subset、多数labelと最大図の時間・メモリを測り、有限cache/cleanup設計を決めます。自動判断用の資源境界が未確定なら必要evidenceとして明記し、暫定値で縮小しません。
8. source配信と生成wheelが同じmanifest assetを取得できる案、必要fontだけsame-origin lazy loadする案、取得失敗時retry案を具体化します。新しいruntime libraryやshaperは先回りして導入しません。既存経路で不可能なcheckpointは具体的な再現と必要boundaryを示します。

## 完了条件と停止条件

日本語・簡体/繁体中国語と既存文字の混在でfont保持、coverage、字形・配置、原文抽出が検証でき、font/配布案とコストが再現可能になったらS00合格です。必要asset・版・styleと対応範囲を明示し、実装が選ぶべき有限データを引き渡してください。

required evidenceが不足する場合は不足checkpointだけを保留し、未検証を合格としません。製品outcomeの変更や追加dependencyが必要なら具体的な境界を提示します。既知のerror情報修正は独立したconcernとして扱えます。

## 終了時のcommit・push

`SESSION_00_RESULT.md`に再現commands、結果、limits、S01へ渡す証拠とS02のfont案を保存します。総合計画書の終了手順で差分reviewと検証を行い、**commitして同名remote branchへpushしてください**。完了/保留はpush済み結果fileで区別します。

Commit title例: `Record Japanese and Chinese PDF font feasibility evidence`。
remote SHA一致を確認し、次sessionにbranch、commit SHA、結果fileを引き渡してください。
