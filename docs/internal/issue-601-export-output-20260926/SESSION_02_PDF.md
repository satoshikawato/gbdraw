# INSTRUCTION PROMPT — S02 日本語・中国語の部分font切替

あなたはgbdraw Issue #601のPDF実装担当です。日本語と中国語（簡体字・繁体字）のPDF出力に対応し、既存fontで表せるアルファベット・数字・記号とstyleを維持して、不足部分だけ対応fontへ自動切替してください。文字検索・選択できるfont埋込みPDFを既存のPDFボタンで取得できるようにします。

## branch取得と必読資料

SESSION_IDは`s02`です。[総合計画書](./MASTER_PLAN.md)の共通取得手順で**`fix/issue-601-export-output-20260926`**をremoteから取得し、専用cloneで作業してください。S00のpush済みevidence、S01のauthority commitとmerge状態、両Decision Packs、AGENTS.md、CLAUDE.md、Web CLAUDE、両ratchetを読みます。

S01のPDF契約が`origin/dev`へ統合されたことをsourceとSHAで確認し、共通手順の通常mergeで実装branchへ取り込んでください。未mergeのcandidateや本計画の承認記録だけで依存runtimeを開始しません。

## 所有範囲

省略したJS pathは`gbdraw/web/js/`からの相対pathです。

- `gbdraw/web/js/services/pdf-fonts.js`: primary維持、coverage、fallback/run選択、必要font登録。
- `services/export.js`: 同snapshotのPDF処理、staging計測とcleanup。
- `app/python-helpers.js::read_pdf_font`と既存READ_PDF_FONT glue: 検証済みassetの取得allowlist。
- `gbdraw/_web_assets.py`、既存package準備owner、必要なfont原本/manifest、packaging/include/license: 一つの再現可能準備経路。
- `app-setup.js`: 限定loader wiringとexport recoveryの必要変更だけ。
- 関連PDF/export/packaging testsと`SESSION_02_RESULT.md`。

新しい正式authority、guard/CIをruntimeと同時に変更しません。shared filesの先行変更を保ち、このsessionの範囲に合わせて統合します。一般errorのsource/transport/文言はS03へ引き渡し、並行error基盤を実装しません。

## 実施内容

1. S00で選定したasset/style/coverageとsource→生成representationの経路を実装します。原本、版、checksum、license、region、coverageと生成先を有限manifestで管理し、sourceとwheelの両配信へ収束させます。font upload/汎用plugin基盤を追加しません。
2. `preparePdfFonts`で既存computed styleをcaptureし、従来のprimary family/styleで表現できる文字を維持します。不足部分だけfontを選び、surrogate pairや結合単位を壊さず、連続する同font/styleをrunへまとめます。各tspanの文字とstyleを一つの親family設定で上書きしないでください。
3. 明示text languageを地域字形へ利用し、情報がない場合はS00で確定したmanifestの固定順とcoverageを使います。UI language、OS font、file列挙順、漢字codepointからの言語推測で変えません。
4. font registration前後のcoverageとstyleを検査し、fallback後も欠字がないことを確かめます。取得失敗、未対応glyph/style、layout failureの識別可能な情報を既存error境界へ渡します。raw文言を公開する専用UIを作らず、S03がnormalizerへ統合できるようcode/stage/safe paramsを定義します。
5. 必要fontだけsame-origin lazy loadし、既存single diagram Workerのhelperを使います。cacheを固定asset数で有界にし、取得失敗Promiseを再利用せずretryできるようにします。DOM計測でも同じfontをload完了させ、文字符号とPDF glyphの対応を保ちます。
6. `captureSvgExport`の同期snapshot、filename、page寸法、anchor、baseline、spacing、rotation、textPath flatteningの既存順序を維持します。文字run分割後も単一labelの全体anchorを保持してください。save前のエラーは部分PDFを保存せず、同snapshotのSVG/PNG取得へ回復可能にします。
7. 必要資産以外を先読みせず、PDF glyph subset、初回/反復時間、最大図のparse/encode/cache/cleanupを測ります。font table/派生glyphの削除で字形を壊さないでください。
8. 置き換えたLiberationだけの強制選択を削除し、font-selection ownerとPDF pipelineを一つに保ちます。new dependency、schema、compat readerを先回りして作りません。

## 必須検証

S00 corpusの日中・混在、normal/bold、既存italic/bolditalic、Sans/Serif/Mono、tspan/textPath/anchor等を実PDFで検査します。原文text抽出一致、font埋込み、image代用なし、primary部分のfamily/style維持、同page寸法を自動assertしてください。supplementary/結合・異体字のsuccessまたは安全な明示非対応も確認します。

同snapshotの内容と名前、出力途中のResult変更、font/library取得失敗後retry、staging/object URLのcleanup、PNG DPIとSVGの不変を検査します。Gallery品質の実図を別PDF rendererで表示し、文字の可読性と配置、legendを目視確認します。

source/wheel配布、license/include、same-origin/offline/no-data-uploadも確認します。実際のoffline auditではbrowser-offline-qa skillを読み、wheelを専用cloneで生成し、sourceとwheelのevidenceを別に記録します。generated browser wheelやdist/egg-infoはcommitしません。

新しいblocking assertionsは総合計画書の既存PR inventoryへ登録します。Playwright local specだけをPR automatic safetyの証拠にしません。trusted-base checkerとarchitecture contractsを通し、production/tests/docs/generated差分を別々にreviewします。新しいdependencyやauthorityが必要なGate失敗は原因boundaryを提示し、実装でGateを弱めません。

## 終了時のcommit・push

`SESSION_02_RESULT.md`にasset/版/容量/coverage、manifest/owner/path、commandsとresults、artifact digest、required acceptance、S03へ渡すerror契約、authority base SHAとmerge状態を記録します。共通終了手順で**commitし、同名remote branchへpushしてください**。
Commit title例: `Support Japanese and Chinese text in Web PDF exports`。
remote SHA一致を確認してS03へhandoffします。
