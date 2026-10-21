# Web UI translation style: ja, zh-Hans, zh-Hant, ko, es, fr, pt

Revision 4, 2026-10-21. Companion of `I18N_GLOSSARY.tsv`. Translators and LLM drafts
follow this file and the glossary together. The glossary decides words; this
file decides how sentences are written around them. Where the two disagree,
the glossary wins and this file is corrected.

## 1. All locales

### 1.1 What is never translated

- **The figure.** The rendered SVG and every export (SVG, Interactive SVG,
  PNG, PDF) keep their English legend words, axis labels, titles, and popup
  script. So do the Source recipe command, the Exact replay file, Session
  JSON, the diagnostics text inside the Details block, console output, and
  `data-*` hooks.
- **Identifiers.** Setting names and paths (`large_tick_interval`,
  `record_selector`, `losatp_max_hits`, `config.form.keep_definition_left_aligned`),
  CLI options (`--feature_override_table`), option values shown as code
  (`above_feature`, `one_to_one`, `exclude_matching`, `overflow=error`),
  renderer ids (`features`, `ticks`, `dinucleotide_content`,
  `dinucleotide_skew`, `depth`, `annotations`), error codes
  (`LOSAT_THREADING_UNAVAILABLE`), Python exception names (`ValueError`), and
  the record selector syntax (`#1`, `#<number>`, `record_id:start-end[:rc]`).
  The sentence around an identifier is translated; the identifier is
  interpolated as a parameter.
- **Biological and file names.** Feature types (`CDS`, `tRNA`, `rRNA`,
  `misc_feature`), qualifier keys with their slash (`/gene`, `/product`,
  `/organism`), record IDs, accessions, file names and extensions, species
  names in placeholders (*Escherichia coli* O157:H7 Sakai), and user-typed
  names (legend captions, slot ids, Session titles).
- **Product, format and program names.** gbdraw (always lower case), LOSAT,
  LOSATN, LOSATP, TLOSATX, BLAST, NCBI BLAST+, blastn, blastp, megablast,
  dc-megablast, GenBank, DDBJ, GFF3, FASTA, TSV, JSON, SVG, PNG, PDF,
  BLAST outfmt 6/7, DPI, CLI, Python, JavaScript, WebAssembly, Pyodide,
  SharedArrayBuffer, Web Workers, IGV, Bakta, Pharokka, Liberation Sans /
  Serif / Mono, and the Circular presets Tuckin, Middle, Spreadout.
- **The Language menu entries** (English, Español, Français, Português,
  日本語, 简体中文, 繁體中文, 한국어) are endonyms and are not translated; the
  menu's own label is (言語 / 语言 / 語言 / 언어 / Idioma / Langue / Idioma
  for ja / zh-Hans / zh-Hant / ko / es / fr / pt). Gallery example titles from `tutorial-data/manifest.json` stay
  English too.
- **Names fixed in English by PD-OI-093 and PD-OI-105** stay English in the
  English UI only. Other locales translate them with the glossary (D-11).

### 1.2 Numbers, units, symbols

- ASCII digits only: no full-width digits and no kanji or Chinese numerals for
  quantities. Decimal point and digit grouping stay as the English source
  shows them, because the app formats numbers with `en-US`.
- Units stay Latin and keep the source's spacing between number and unit:
  `100 bp`, `1kb`, `20px`, `0.65 ×R`, `8 pt`, `10x`. Never convert a unit.
  Time words are words, not units: 30 seconds is 30秒 in `ja`, `zh-Hans` and `zh-Hant`.
- Keep `→ ← · × + − %` and the en dash of ranges (`Rows 3–4`) as written.
- Ellipsis is one `…` (U+2026) whether the source writes `...` or `…`, both
  for "more input follows" (Choose records…) and for progress (Saving
  session…).

### 1.3 Placeholders and keys

- The English source text is the key (D-02). There is no English catalog:
  `t('Save Session')`, `tn(count, one, other)` (the `other` text is the key),
  `tx(context, text)` (key `context + '\u0004' + text`, for homonyms), `msg()`
  to keep a text until `tm()` translates it where it is shown. A catalog
  `<locale>.js` maps those keys to plain strings.
- `{name}` placeholders are copied verbatim: the same set as the English key,
  never translated, never added or dropped. They may move inside the
  sentence.
- A `tn` key has one string in `ja`, `zh-Hans`, `zh-Hant` and `ko`, which have
  no plural forms (`Intl.PluralRules` returns `other`), and a `one` and an
  `other` string in `es`, `fr` and `pt` (§5.5 to §8.5).
- One key is one sentence or label. Never build a sentence from translated
  fragments; ask for a key with placeholders instead. A sentence with an
  inline element (`<strong>`, a swatch, a bound name) is one key; the
  `<i18n-text>` component places the element at its `{name}` slot (D-04).
- Two sentences shown together are two keys. The code joins them with a
  language-dependent separator: a space in English, none after 。 in `ja` and
  `zh-*` (D-18). Do not add a trailing space to a value.
- Quotes follow the English: a `{name}` that the English key quotes is quoted
  in the target language's quotes (`ja` 「」, `zh-Hans` “ ”, `zh-Hant` 「」,
  `ko` ‘ ’, `es` and `pt` “ ”, `fr` « »);
  one the English leaves bare stays bare (D-25). English double quotes never
  appear in `ja` or `zh-*` prose.

### 1.4 Words and capitals

- Use the glossary rendering for every English term it lists, in every key,
  including compounds: Record label is レコードラベル because Record is
  レコード and Label is ラベル. A row tagged `[review]` marks a judgment call;
  it is still the value to use until the glossary changes.
- Title Case, ALL CAPS headings (PRESET SCHEMES, DEFAULT COLORS, SPECIFIC
  RULES, LABEL FILTERING, QUALIFIER PRIORITY, LABEL GEOMETRY, ADD KEY) and
  letter spacing carry no meaning in CJK; translate them as ordinary
  headings. Outside the figure, the CJK rules in `gbdraw/web/index.html` turn off
  `uppercase` and `tracking-*`.
- A word that the source capitalizes as a name of a control (Generate,
  Apply, Save Session) is written as that control's translated label, so a
  reader can find it: 「生成」をクリック / 点击“生成” / 點選「產生」.

### 1.5 Length

- CJK glyphs are about twice as wide as Latin letters at the same size, and
  the CJK CSS raises 10–11 px text to 12 px. Count width, not characters.
- Rail short names: at most 5 full-width glyphs, 4 preferred (the 60 px icon
  column has 52 px buttons at 10 px). Header and bar buttons: no wider than
  the English label. Option labels inside 10 px rows: as short as the
  glossary form.
- Help texts and notices may be as long as needed; they wrap.
- `es`, `fr` and `pt` labels run 15 to 35% longer than English; their length
  rules are §6.6, §7.6 and §8.6.

## 2. Japanese (`ja`)

### 2.1 Vocabulary

- **Katakana** for established loanwords, with the JIS/Microsoft long-vowel
  style: words ending in -er, -or, -ar take ー (エディター,
  フィルター, レンダラー, レイヤー, ユーザー, サーバー, パラメーター,
  フォルダー, ヘッダー, マーカー). Words that are established without it keep
  their form (データ, プロトコル, インターフェイス), and -y endings follow
  Microsoft (メモリ, ライブラリ, クエリ, エントリ; but ギャラリー, ファミリー
  are established with ー). Exceptions: the seed クオリファイア is the
  maintainer's form; ブラウザ is written without ー, as in docs/internal (R1).
- **Kanji** for native concepts: 凡例, 軸, 目盛り, 鎖, 逆相補, 配列, 遺伝子,
  結果, 設定, 履歴, 生成, 比較, 深度, 環状, 直鎖状, 余白, 切り出し.
- **既定**, not デフォルト, for Default (D-13). **直鎖状**, not 線状, for
  Linear (D-20); the Record Size Mode option Linear is 線形. As a modifier of
  a diagram: 環状図の, 直鎖状図の.
- **Depth**: the Depth input file is リード深度TSV; every other Depth (track,
  axis, series, window) is 深度 (D-27). 一致率 for Identity (D-31).
- **DDBJ terms** for INSDC fields: 位置情報 (Location), 遺伝子産物
  (/product), アクセッション番号, 備考 (/note). Bitscore stays Latin (D-31).
- Never 記録 or 配列 for Record, never 特徴 or 機能 for Feature (seeds).
- Okurigana as in Office and the glossary: 目盛り, 読み込み, 貼り付け,
  並べ替え, 切り出し, 元に戻す, 折りたたむ.
- Prefer Office command words: 元に戻す, やり直し, 貼り付け, 名前を変更, 削除,
  追加, 適用, キャンセル, 閉じる, 参照 (Browse).
- No ・ or spaces inside the katakana compounds the glossary lists
  (カスタムトラックスロット, フィーチャータイプ); use の where Japanese wants
  it (フィーチャーの配置).

### 2.2 Sentence endings

- Labels, buttons, headings, options, column headers, tabs, badges: 体言止め
  or the Office command form (セッションを保存, サンプルを読み込む, 元に戻す).
  No です・ます and no final 。 in a label.
- Help texts, notices, dialogs and error messages: です・ます in complete
  sentences ending with 。. An English imperative becomes 〜してください
  (hiragana ください). No 〜しましょう, no 〜して下さい.
- Status and progress: 〜中… for work in progress (保存中…, 読み込み中…,
  生成中…, 適用中…; with an object: セッションを保存中…, 図を生成中…).
  Completion: 〜しました。 Busy reasons: 図を生成中です。生成が終わってから再試行してください。

### 2.3 Punctuation and spacing

- Prose uses full-width 、。（）「」：. No space follows 。 between two
  sentences (D-18). Parentheses are full-width when they
  hold Japanese (領域（任意）, 自動（既定）) and half-width, with no
  space, when they hold Latin, a number or a unit (許容値(bp), 最低一致率(%)).
  A label colon is full-width (状態：) unless Latin follows it directly. The spaced slash of
  the source stays: 生物名 / 株名, as in the seed 二塩基含量 / スキュー.
- **No half-width space between Latin or digits and Japanese**, matching the
  maintainer's seeds: GC含量, LOSATを実行, 10件, BLAST TSVを選択 (the space
  inside the English term BLAST TSV stays). The number-unit spacing of the
  source stays: 100 bp.

### 2.4 Placeholders and identifiers

- `{n}` attaches without spaces and takes a counter: 件 for records,
  features, rows, items and edits ({n}件のレコード, {drawn}/{total}件を描画);
  個 for other objects; 回 for repetitions.
- Interpolated identifiers attach without spaces and without quotes:
  tick_intervalではなくlarge_tick_intervalを使用してください。 After a locator
  label use the full-width colon with no space: 設定：{path}。 行：{n}。
- User-typed names take 「」: 「{name}」に統合.

### 2.5 Length

Japanese is usually 0.5–0.8× the English character count, which is about
the same pixel width. Buttons stay within 8 full-width glyphs where
possible; rail short names within 5 (Features is the known exception).

## 3. Simplified Chinese (`zh-Hans`, mainland usage)

### 3.1 Vocabulary

- 设置 for Settings, 默认 for Default, 预设 for Preset, 会话 for Session,
  生成 for Generate, 加载 for Load, 保存 for Save, 导入 / 导出 for Import /
  Export, 上传 / 下载 for Upload / Download, 文件 for File, 数据 for data.
- Mainland computing terms: 撤销 / 重做, 添加, 选择, 应用, 重置, 自定义,
  高级, 布局, 选项卡, 字段, 缓存, 线程, 串行 / 并行.
- Biology: 基因组, 链 (strand), 氨基酸, 细胞器, 登录号, 遗传密码表.
- Terms settled in rounds 1 and 2 (see the glossary): 子轨道 (Lane), 指定参考
  (Exact reference), 共线性 (Collinear), 软件包, 几何参数, 片段 (feature
  part), 单元 (Collinear unit), 测序深度 / 深度, 比例尺 for scale bar and
  length bar. Round 3 (D-33): GC偏斜 (GC skew), 运行时 (Runtime), 部分
  (Section), 成比例 (the Record Size Mode option Linear). The glossary holds
  the current value.

### 3.2 Buttons and sentences

- Buttons and menu items are a verb or verb-object without 请: 保存会话,
  加载示例, 应用, 取消, 重试.
- Help texts and messages are complete sentences ending with 。. 请 may open
  an instruction in a message (请选择一条记录。) but not more than once per
  message; avoid 您 unless the sentence needs a pronoun.
- Progress: 正在保存…, 正在加载会话…, 正在生成图谱…; completion: 已保存。;
  busy reasons: 正在生成图谱。生成完成后请重试。

### 3.3 Punctuation and spacing

- Full-width ，。：；！？（）、 everywhere, including labels with Latin; no
  space follows 。 between two sentences (D-18):
  最小一致性（%）. Quotes “ ”, nested ‘ ’.
- No space between Chinese and Latin or digits: GC含量, 运行LOSAT,
  上传BLAST TSV, 10条记录. The exception is an interpolated identifier, which
  takes a half-width space on each side inside a sentence:
  请使用 large_tick_interval 代替 tick_interval。 After a full-width colon no
  space: 设置：{path}。

### 3.4 Placeholders and measure words

- {n} attaches without a space and takes a measure word: 条 for records,
  rules, tracks and edits ({n}条记录), 个 for files, features and items,
  行 for table lines (第{n}行).

### 3.5 Length

Usually 0.4–0.6× the English character count. Rail short names 2–4 glyphs.

## 4. Traditional Chinese (`zh-Hant`, Taiwan usage)

### 4.1 Locale label and vocabulary

- The Language menu shows plain **繁體中文**, with no 台灣 suffix (Owner
  decision, 2026-10-19). The vocabulary is nevertheless Taiwan usage
  throughout; a Hong Kong catalog would be a separate `zh-Hant-HK` later.
- Required words: 設定 (Settings), 預設 (Default), and therefore **預設組**
  (Preset), 工作階段 (Session), 產生 (Generate), 檔案 (File), 資料 (data),
  介面 (interface), 執行 (run), 載入 (Load), 儲存 (Save), 匯入 / 匯出
  (Import / Export), 下載 / 上傳 (Download / Upload).
- Microsoft Taiwan UI words: 新增, 選取 (Select; Choose stays 選擇), 套用,
  重設, 復原 / 重做, 進階, 版面配置, 分頁 (tab), 說明 (Help), 自訂, 欄位,
  程式, 指令, 執行緒, 快取, 預覽, 搜尋, 字型, 對話方塊, 點選, 拖曳, 略過,
  取代, 內嵌, 覆寫, 位移.
- Table structure: a row is **列** and a column is **欄** (the reverse of the
  mainland 行 / 列).
- Biology: 基因體 (genome), 股 (strand: 雙股、正股), 胺基酸, 胞器, 登錄號,
  遺傳密碼表.

### 4.2 Punctuation and spacing

- Full-width ，。：；！？（）、 as in `zh-Hans`; quotes **「」**, nested 『』.
- Spacing rules are those of §3.3: no space between Chinese and Latin or
  digits, a half-width space around interpolated identifiers.

### 4.3 Where Taiwan and Hong Kong differ

The catalog follows the Taiwan column. Hong Kong readers will notice:

| Concept | Taiwan (this catalog) | Hong Kong |
| --- | --- | --- |
| Session | 工作階段 | 會話 |
| data | 資料 | 數據 |
| interface | 介面 | 界面 |
| software | 軟體 | 軟件 |
| network | 網路 | 網絡 |
| genome | 基因體 | 基因組 |
| strand | 股 | 鏈 |
| row / column | 列 / 欄 | 行 / 列 |
| thread | 執行緒 | 線程 |
| cache | 快取 | 緩存 |
| print | 列印 | 打印 |
| amino acid | 胺基酸 | 氨基酸 |
| Annotation | 註解 | 註釋 |

### 4.4 How the `zh-Hant` catalog is produced

1. Start from the `zh-Hans` catalog and convert it with OpenCC
   **`s2twp`** (script plus Taiwan phrases). Never `s2t` (characters only),
   and never the reverse direction.
2. Correct the output against this file: the table of §4.4.1 lists every
   glossary term whose `zh-Hant` value differs from the `s2twp` output of its
   `zh-Hans` value. Replace each by hand, in context, and check any term the
   table does not cover against `I18N_GLOSSARY.tsv`.
3. A native Taiwan reader checks the whole catalog.

#### 4.4.1 Words `s2twp` gets wrong for this app

Measured on 2026-10-19 by converting the `zh-Hans` column of the glossary:
87 of 372 rows come out different from the `zh-Hant` decision. The run used
a Python re-implementation of the `s2twp.json` chain (STPhrases with the
regional phrases, STCharacters, then TWPhrases, TWVariantsPhrases,
TWVariants) over OpenCC's text dictionaries at commit `2a09f3a7`
(2026-10-18), because the OpenCC binary could not be installed in this
environment. Re-run the real `opencc -c s2twp.json` over the `zh-Hans` strings
before relying on the list; segmentation may differ in a few places.

| `s2twp` output | Correct for this app | Why |
| --- | --- | --- |
| 會話 | 工作階段 | Session (product decision) |
| 生成 | 產生 | Generate, Regenerate, Generate Diagram |
| 預設 (from 预设) | 預設組 | Preset; 預設 is Default |
| 註釋 | 註解 | Annotation, Region Annotations, Annotation set |
| 示例, 示例庫 | 範例, 範例庫 | Example, Gallery |
| 特徵型別 | 特徵類型 | 型別 is a programming data type |
| 應用 | 套用 | Apply, Applied, Applies on Generate, Pending |
| 重置 | 重設 | Reset, Reset Settings |
| 撤銷 | 復原 | Undo |
| 選擇 | 選取 | Select only; Choose keeps 選擇 |
| 高級 | 進階 | Advanced |
| 佈局 | 版面配置 | Layout, Layout edit, Track layout |
| 選項卡 | 分頁 | Tab |
| 幫助 | 說明 | Help |
| 欄目 | 區段 | Section (and 分區 becomes 分割槽, a disk partition) |
| 自定義 | 自訂 | Custom, Custom Track Slots |
| 固定, 停靠 | 釘選, 停駐 | Pin, Dock |
| 拖動, 跳過, 替換 | 拖曳, 略過, 取代 | Drag, Skip, Replace |
| 重用 | 重複使用 | Reuse |
| 修復 | 修正 | Fix |
| 嵌入, 覆蓋, 偏移 | 內嵌, 覆寫, 位移 | Embed, Override, Offset |
| 歸一化 | 正規化 | Normalize |
| 命令, 請求 | 指令, 要求 | Command, Request |
| 執行時, 配置 | 執行階段, 組態 | Runtime, Configuration |
| 後設資料, 架構 | 中繼資料, 結構描述 | Metadata, Schema |
| 字首 | 前綴 | Prefix (R2); Suffix is 後綴 |
| 許可證, 開源 | 授權條款, 開放原始碼 | License, Open source |
| 不可用, 可選 | 無法使用, 選填 | Unavailable, Optional |
| 必需, 缺失 | 必要, 缺少 | Required, Missing |
| 推薦 | 建議 | Recommended |
| 序列 (from 串行) | 循序 | Serial execution |
| 並行 | 平行 | Parallel runs |
| 跨源 | 跨來源 | Cross-origin isolated |
| 高亮 | 醒目提示 | Highlight |
| 標尺 | 尺規 | Ruler |
| 描邊, 填充 | 筆畫, 填色 | Stroke, Fill |
| 渲染, 渲染器 | 繪製方式 (noun) or 繪製 (verb), 繪製器 | Rendering, Renderer (R2) |
| 泳道 | 子軌道 | Lane (R1) |
| 字型系列 | 字型家族 | Font family |
| 圖例項 | 圖例項目 | Legend item |
| 行, 列 | 列, 欄 | Row, Column: swap both in one pass, never one replace after the other |
| 基因組 | 基因體 | Genome, Input Genomes |
| 鏈 | 股 | Strand, Separate strands, Forward strand |
| 細胞器 | 胞器 | Organelle |
| 遺傳密碼錶 | 遺傳密碼表 | Genetic code (錶 is a watch) |
| 登入號 | 登錄號 | Accession (登入 is log in) |
| 標識 | 識別資訊 | Identity of a record or feature |
| 匹配 | 符合 | Search matches only; a comparison match is 比對 (連線 in the figure, 結果 in counts and lists) (R2) |
| 溢位 | 溢出 | Layout overflow; 溢位 reads as arithmetic overflow |
| 原始檔 | 來源檔案 | Source recipe (源文件 means input files here) |

`s2twp` already gets these right, so they need no correction: 檔案, 資料,
介面, 執行, 載入, 儲存, 匯入, 匯出, 新增, 搜尋, 字型, 快取, 執行緒, 欄位,
對話方塊, 點選, 預覽, 視窗, 即時, 互動式, 剪貼簿, 貼上, 重新命名,
正規表示式, 調色盤, 支援, 摺疊, 停用, 胺基酸, 擷取.

### 4.5 Linear and Collinear

The diagram mode Linear is 線狀 (D-33). The Record Size Mode option Linear
has its own key, `tx('size scaling', 'Linear')`, and is 線性 (proportional
to record length). Collinear is 共線性 (R2).

### 4.6 Length

As `zh-Hans`. 工作階段 is four glyphs where 会话 is two; check the header
buttons 儲存工作階段 and 載入工作階段 at 1,024 px in PR-1.

## 5. Korean (`ko`)

Standard Korean (표준어) as South Korean software writes it: Microsoft and
Google Korean UIs for command words, Korean molecular-biology textbooks and
papers for the science. The Language menu entry is **한국어** (D-03).

### 5.1 Vocabulary

- **Field terms**: 유전체 (genome), 서열 (sequence), 염기서열 (the Nucleotide
  tab), 아미노산 서열 (the Amino acid tab), 가닥 (strand; 정방향 가닥), 역상보
  (reverse complement), 주석 (annotation), 유전자, 단백질, 유전자좌 (locus),
  균주 (strain), 생물명 (organism), 종명 (species), 세포소기관, 등록 번호
  (accession), 유전 암호표 (genetic code), 산물 (/product), 비고 (/note), 정의
  (DEFINITION), 위치 (Location), GC 함량, 정렬 (sequence alignment), 동일성
  (identity), 유사성, 보존성, 공선성 (collinear), 쌍별 (pairwise), 깊이 / 리드
  깊이 (depth), 원형 / 선형 (Circular / Linear), 주 눈금 / 보조 눈금, 범례,
  축척 (scale), 축척 막대.
- **Kept in Latin**, as Korean papers write them: GC skew, AT skew, Bitscore,
  E-value, orthogroup, RBH, and every §1.1 name. Skew alone is 스큐.
- **Loanwords** in the standard spelling (외래어 표기법) Microsoft uses: 피처,
  레코드, 트랙, 레인, 슬롯, 앵커, 히트, 쿼리, 블록, 리본, 링, 밴드, 마크,
  레이블, 라이선스, 메타데이터, 스키마, 런타임, 렌더러, 렌더링, 오버레이,
  오프셋, 팔레트, 프리셋, 캔버스, 레이어, 레이아웃, 세션, 섹션, 필드, 필터,
  패턴, 스택, 스텝, 유닛, 다이어그램, 다이뉴클레오타이드. Exception: 윈도우
  (sliding window), as bioinformatics writes it, not the standard 윈도. Add no
  loanword the glossary does not list. The glossary's 0.15.0 rows add
  these, which are allowed: 템플릿, 시드, 패널, 라이브러리, 인코딩, 파이프라인,
  워크플로, 해시, 플랫 파일, 포커스, 그라데이션, 픽셀, 인트론, 레플리콘,
  프리프린트, 캐시, 사이드바, 호스팅, 이스케이프 문자, 웹 앱, and the
  loanword parts of 파라로그, 오솔로그 and 갭 열림 수 (rows of
  `I18N_GLOSSARY.tsv`). The compound 소스 코드 (Source code) is allowed.
- **Command words** (Microsoft Korean): 저장, 열기, 닫기, 취소, 추가, 삭제
  (Delete), 제거 (Remove), 지우기 (Clear), 복사 (Copy), 복제 (Duplicate),
  붙여넣기, 편집, 적용, 선택, 모두 선택 해제, 표시 / 숨기기, 켬 / 끔, 활성화 /
  비활성화, 실행 취소 / 다시 실행, 되돌리기 (Revert), 복원 (Restore), 초기화
  (Reset), 다시 시도, 건너뛰기, 바꾸기, 이름 바꾸기, 병합, 가져오기 /
  내보내기, 불러오기 (Load), 찾아보기 (Browse), 고정, 도킹, 접기 / 펼치기,
  확대/축소, 맞춤 (Align); nouns 미리 보기, 대화 상자, 편집기, 도움말, 세부
  정보, 사용자 지정, 기본값, 재정의, 임계값, 정규식, 글꼴, 오류, 기록 (History).
- **Never**: 게놈, 시퀀스, 스트랜드, 어노테이션, 특징 or 기능 for Feature,
  기록 for Record, 라벨, 라이센스, 디폴트, 컬러 / 색상 (use 색), 폰트, 에러,
  리셋, 쓰레드, 끌기 (use 드래그), 소스 (also "sauce"; use 원본, except in the compound 소스 코드), 당신 / 귀하
  (leave the subject out), 유저.

### 5.2 Buttons and sentences

- **Labels, buttons, headings, menu commands, options, tabs, badges**: noun
  forms, no final period. A Sino-Korean action noun alone or after its object,
  without a particle: 저장, 세션 저장, LOSAT 실행, 파일 제거, 다이어그램 생성.
  A native verb in its -기 form: 지우기, 숨기기, 이름 바꾸기, 건너뛰기,
  되돌리기, 예제 불러오기. A command that opens a further step keeps `…`
  (맞춤…, 레코드 선택…).
- **Sentences** (help, notices, dialogs, errors): 합니다체 (~합니다,
  ~했습니다, ~할 수 없습니다). Never 해요체 (~해요, ~하세요), 해라체 or
  ~하시오. Instructions: ~하십시오. Dialog questions: ~하시겠습니까?
- **Progress**: a Sino-Korean noun takes 중… (저장 중…, 적용 중…, 준비 중…;
  with its object: 세션 저장 중…, 다이어그램 생성 중…); a native verb takes
  -는 중… with its particle (세션을 불러오는 중…, 원본 레코드를 읽는 중…).
  Busy reasons are two keys: 다이어그램을 생성하는 중입니다. 생성이 끝난 후
  다시 시도하십시오.
- **Completion**: ~했습니다 in a sentence (저장했습니다.); -됨 in a state
  label (적용됨, 준비됨, 활성화됨), negative ~ 안 됨 (적용 안 됨).
- **Errors** say what failed (~에 실패했습니다. / ~할 수 없습니다.), then
  what to do in ~하십시오. A control name in a sentence is its translated
  label in ‘ ’: ‘생성’을 클릭하십시오.

### 5.3 Punctuation and spacing

- Half-width `. , : ; ? ! ( )`. A sentence ends with `.`; two sentences are
  separated by one space (D-06). A sentence whose English has no final period
  keeps none (help hints, tooltips, titles, and the `{effect}` fragments that
  an outer key ends).
- **Quotes**: ‘ ’ (U+2018, U+2019) around a placeholder the English quotes
  and around a control name: ‘{name}’(으)로 병합. Windows and Google Korean
  messages quote names in single quotes; this is their typographic form.
  English `"…"` never appears; “ ” is not used.
- **Label colon**: no space before, one after (상태: 준비됨). **Parentheses**
  attach to what they gloss: 영역(선택 사항), 레코드 간격(px), 자동(기본값).
- **Word spacing** follows 한글 맞춤법; compound UI terms are spaced as
  Microsoft spaces them (글꼴 크기, 실행 취소, 미리 보기, 대화 상자), and the
  glossary's spacing is binding. These stay joined: 염기서열, 세포소기관,
  유전자좌, 눈금자, 기본값, 임계값, 정규식, 우선순위, 축상.
- **Latin and digits**: a particle attaches with no space (LOSAT을, TSV를,
  SVG에); a Korean noun after a Latin word takes a space (GC 함량, BLAST TSV
  업로드). A counter attaches to its number ({count}개, 3행); a unit keeps the
  source spacing (100 bp, 1kb, 20px).
- Ranges keep the source's en dash (3–4행), never `~`; a slash keeps the
  source spacing (생물명 / 균주, GC/스큐); one `…`, attached.

### 5.4 Words and capitals

- Hangul has no case. Title Case and ALL CAPS headings become ordinary
  headings: PRESET SCHEMES 프리셋 배색, DEFAULT COLORS 기본 색, SPECIFIC
  RULES 개별 규칙, LABEL FILTERING 레이블 필터링, QUALIFIER PRIORITY 한정자
  우선순위, LABEL GEOMETRY 레이블 치수, ADD KEY 키 추가. The CSS that turns
  off `uppercase` and `tracking-*` for CJK must match `:lang(ko)` too.
- Latin keeps its own case (gbdraw, LOSAT, orthogroup). Compounds reuse the
  glossary words: 레코드 레이블 (Record + Label), 트랙 프리셋 (Track + Preset).
- One English word keeps one Korean word: Remove is always 제거, Delete 삭제,
  Clear 지우기, Reset 초기화. No Hanja, not even in parentheses.

### 5.5 Placeholders and plurals

- `{name}` is copied verbatim and moves to Korean order (object before verb).
- **Particles**: the value's last syllable is unknown, so an alternating
  particle is written in both forms, post-consonant first: 을(를), 이(가),
  은(는), 과(와), (으)로. Non-alternating particles attach plainly (의, 에,
  에서, 도, 만, 까지); the copula 입니다 fits both.
  `Remove {name}?` → {name}을(를) 제거하시겠습니까?;
  `Rename "{from}" to "{to}"?` → ‘{from}’의 이름을 ‘{to}’(으)로
  바꾸시겠습니까?; `{first} and {second}` → {first}과(와) {second}.
- Prefer a known noun after the placeholder, so the particle attaches to
  Korean: {file} 파일을, {track} 트랙을, {names} 글꼴은. Use the paired form
  only when no such noun fits.
- **Plurals**: `Intl.PluralRules('ko')` always returns `other`, so there is
  one form. Translate the `other` key; the build script gives the `one` key
  the same string (D-05). No plural marker 들 on a counted noun; English
  agreement (its / their, was / were) disappears.
- **Counters**: noun first, counter attached to the number: 레코드 {count}개,
  not {count}개의 레코드. 개 for records, features, files, tracks, colors,
  entries; 건 for issues, results, matches, edits; 쌍 for pairs; 행 for table
  lines; 번째 for an ordinal. The particle follows the counter ({count}개를).
- Worked examples (`plural-pairs.tsv`):
  1. `{count} record` / `{count} records` → 레코드 {count}개
  2. `{count} gene pair` / `… pairs` → 유전자 {count}쌍
  3. `{count} issue to fix before Generate` / `… issues …` → 생성 전에
     수정할 문제 {count}건
  4. `{current} of {count} match` / `… matches` → {count}건 중 {current}번째
  5. `Keep your {count} changed color from the "{palette}" palette ({keys})?`
     → ‘{palette}’ 팔레트에서 변경한 색 {count}개({keys})를 유지하시겠습니까?

### 5.6 Length

- A Hangul syllable is full width; a Korean label has 0.5–0.7× the English
  character count and about the same pixel width. Count width.
- Rail short names: at most 4 syllables, 2 preferred: 입력, 피처, 레이블,
  범례, 색, 축, 비교, 영역, 제목, 고급, 정보, 글꼴, GC/스큐; 레이아웃 (4) is
  the longest.
- Header and bar buttons no wider than the English label: 세션 저장, 설정
  초기화, 다이어그램 생성. 세션 불러오기 and 예제 불러오기 (7 glyphs) are the
  widest; check them at 1,024 px. A shorter form must keep the glossary word
  (never 열기 for Load).
- Labels take the short forms: 적용 안 됨 (not 적용되지 않음), 자동(기본값),
  필수, 누락, 없음, 켬 / 끔; full sentences belong in help and notices, which
  may be as long as needed.
- Korean wraps at word spaces: give `:lang(ko)` `word-break: keep-all` and
  the CJK 12 px floor, so a word is never split between syllables.

## 6. Spanish (`es`)

International Spanish (D-03): words a reader understands in Spain and in the
Americas; no voseo, no tuteo, no vosotros. Where Spain and Latin America
differ, take the word Microsoft uses in both markets or the RAE form
(archivo, not fichero; agregar; orgánulo). The Language menu entry is
**Español**; the menu label is Idioma.

### 6.1 Vocabulary

- **Field terms**: genoma, secuencia, registro (Record), característica
  (Feature; tabla de características), tipo de característica, calificador
  (the keys /gene, /product stay), anotación, número de acceso, locus (plural
  loci), cadena (Strand: cadena +, cadena −, cadena directa), complemento
  inverso, contenido GC, sesgo GC, profundidad (Depth), profundidad de
  lectura (read depth in sentences), cobertura, identidad, homología,
  ortólogo, parálogo, ortogrupo, bloque colineal, código genético, orgánulo,
  organismo, cepa, especie, producto, gen; the tabs Nucleótidos / Aminoácidos.
- **Alignment has two words**: alineamiento for sequences (BLAST, LOSAT,
  longitud del alineamiento), alineación for the layout action of Similarity
  groups (Revisar alineación, Restablecer alineación).
- **Diagram terms**: pista (Track), carril (Lane), ranura (track slot),
  posición (File card slot), marca de graduación (Tick; marca principal /
  secundaria), marca (annotation mark), regla graduada (Ruler; regla is a
  color rule), margen (Padding), separación (Gap), espaciado (Spacing), trazo,
  relleno, muestra de color, fondo / primer plano, preajuste / predefinido
  (Preset), predeterminado (Default).
- **Command words** (Microsoft Spanish): Restablecer (Reset), Eliminar
  (Delete), Quitar (Remove), Borrar (Clear), Agregar, Cambiar nombre,
  Examinar (Browse), Subir (Upload), Cargar (Load), Importar, Seleccionar
  (Select), Elegir (Choose), Contraer, Omitir (Skip), Combinar (Merge),
  Acoplar (Dock), Girar, Mantener (Keep), Volver a generar; Habilitar /
  Deshabilitar (Enable / Disable); Activado / Desactivado (On / Off).
- **Loanwords**: only hit (plural hits), query, subject, E-value, Bitscore,
  zoom, renderizador and the names of §1.1; never track, feature, slot,
  preset, lane, popup, default or reset.
- **Never use**: fichero, ordenador, computadora, salvar, remover, resetear,
  setear, chequear, clicar, pinchar, pulsar (write hacer clic), similaridad,
  accesión, cualificador, hebra (cadena only), organelo, rasgo, coger, vale.

### 6.2 Buttons and sentences

- Register: usted, implicit; su only where a possessive is needed. Never
  tú, vos or vosotros forms (selecciona, seleccioná, seleccionad).
- Buttons, menu commands, links: the infinitive, no article where the English
  has none: Guardar sesión, Cargar ejemplo, Cambiar nombre, Aplicar y generar.
- Labels, headings, tabs, column headers, badges: a noun phrase without a
  final period. Options are nouns or adjectives agreeing with the setting
  (Automático, Activado, En serie). A dialog title that asks is ¿ + infinitive:
  ¿Quitar todo el archivo?, ¿Restablecer toda la configuración?
- Instructions (help texts, notices, dialogs, errors): formal imperative.
  Seleccione un registro. Haga clic en Generar. Vuelva a generar el diagrama.
  Leave blank for auto: Déjelo vacío para usar el valor automático.
- Descriptions: third person present or impersonal se (Muestra la
  profundidad de lectura como una pista. Se aplica al generar.).
- Progress (…ing): gerund + … (Guardando sesión…, Generando diagrama…).
  Completion: a participle agreeing with its noun in a status line (Sesión
  guardada.), se + preterite in a sentence (Se aplicó la alineación.).
- Errors: what failed, then what to do: No se pudo leer el archivo. Compruebe
  que sea un archivo FASTA. Headings: Error de generación, Error de
  exportación SVG. Busy reasons: Generando el diagrama. Vuelva a intentarlo
  cuando termine la generación. Retry is Reintentar on a button and Vuelva a
  intentarlo in a sentence.

### 6.3 Punctuation and spacing

- ¿ and ¡ open every question and exclamation, also inside a sentence (Si no
  coincide, ¿desea quitarlo?).
- Quotes “ ” (U+201C, U+201D), nested ‘ ’. A {name} the English quotes is
  “{name}”; a UI word the English quotes ('Bar') is “Barra”. No straight " or
  ' in prose. A period or comma follows the closing quote.
- Colon: no space before, one after (Campo: {path}). A label ending with a
  colon keeps it.
- Parentheses, slashes and ranges as the source writes them: Región
  (opcional), Organismo / cepa, Filas 3–4.
- ASCII digits with the source's decimal point and grouping (0.65, 1,000),
  because the app formats numbers with en-US; never 0,65. Units and % keep
  the source spacing (100 bp, 50%), although the RAE writes 50 %.
- One … (U+2026), no space before it. e.g. is p. ej.; i.e. is es decir. No
  comma before y / o in a list. Sentences are separated by one space (§1.3).

### 6.4 Words and capitals

- Sentence case for every label, button, heading, tab, option and column
  header: the first word and proper names only (Ranuras de pista
  personalizadas, Comparaciones por pares, Eje y escala). Never copy English
  Title Case.
- A key the English writes in capitals stays in capitals, with accents:
  ESQUEMAS PREDEFINIDOS, COLORES PREDETERMINADOS, REGLAS ESPECÍFICAS,
  FILTRADO DE ETIQUETAS, PRIORIDAD DE CALIFICADORES, GEOMETRÍA DE ETIQUETAS,
  AGREGAR CLAVE, COLOR; the badge OFF is DESACTIVADO.
- The action Generate in running text is the verb in lower case: generar, al
  generar, antes de generar, pendiente de generar. Only the button itself is
  cited with its label Generar (capital, no quotes, as the next bullet says;
  “Generar” where the English quotes it).
- Nouns the English capitalizes as product terms are lower case in a
  sentence: la sesión, el resultado, el registro. A control, section or
  option cited in a sentence keeps its translated label and capital without
  quotes, as Microsoft writes it: haga clic en Generar; en Ranuras de pista
  personalizadas; elija Automático. Where the English quotes it, use “ ”.
- Compounds put the head noun first, joined with de: tipo de característica,
  TSV de profundidad, TSV de BLAST, caché de LOSAT. Genders: el ancla
  (feminine: el ancla habilitada), la CDS, el ID, el TSV, el hit, la query,
  el subject, el locus.

### 6.5 Placeholders and plurals

- A placeholder's gender is unknown. Put the class noun before it and let
  articles and adjectives agree with that noun (el registro {record}, la
  pista {track}, el archivo “{file}”); never make a word agree with the
  placeholder itself ({name} is hidden → El elemento {name} está oculto). A
  placeholder as subject takes a third-person singular verb.
- {count} precedes the noun; the noun and every word agreeing with it follow
  the number. tn takes the `one` string only when Intl.PluralRules('es')
  returns one, which is the count 1 alone; 0 and 1.5 take the plural (0
  registros). Intl returns many for 1,000,000 and other exact millions; tn
  must use `other` there. A `one` string drops {count} only where the English
  one does.
- Worked pairs, `one` / `other` (English `other` key in parentheses):
  1. ({count} records): ({count} registro) / ({count} registros)
  2. ({count} features ready): {count} característica lista / {count}
     características listas
  3. (From these {count} color rules:): De esta regla de color: / De estas
     {count} reglas de color:
  4. (Reset {count} placements to Auto): Restablecer {count} colocación a
     Automático / Restablecer {count} colocaciones a Automático
  5. ({current} of {count} matches): {current} de {count} coincidencia /
     {current} de {count} coincidencias

### 6.6 Length

- Spanish runs 15–30% longer than English; help texts and notices wrap. Rail
  short names: at most 9 letters, about 50 px at 10 px semibold in the
  52 px button, which does not wrap: Entrada, Comparar, Avanzado, Regiones,
  Diseño, Fuente, Títulos, Leyenda, Colores, Etiquetas, Eje, GC / sesgo,
  Acerca de; Generar beside the rail. Features shares one key with its
  section title, and Características (15 letters) overflows; see the notes.
- Header buttons stay within about 1.3× the English label; check Restablecer
  configuración at 1,024 px. Options in 10 px rows use the glossary's short forms: Automático, En serie,
  Multihilo, Arriba / Abajo, Dentro / Fuera, Sobre el eje, Fondo, Cinta /
  Curva. To shorten a label, drop articles first (Ancho de pista, not Ancho de
  la pista); never switch terms, and abbreviate with a period only where a
  screenshot shows an overflow.

## 7. French (`fr`)

Standard French of France, vouvoiement, French typography (D-03, D-06). The Language menu
entry is the endonym **Français**; the menu label itself is Langue.

### 7.1 Vocabulary

- **Biology and bioinformatics**, as French scientific writing has them: génome, séquence, gène,
  protéine, brin (brin +, brin −), complément inverse, annotation, élément (Feature), type
  d’élément, enregistrement (Record), numéro d’accession, locus (invariable), qualificateur,
  produit (/product), localisation (Location), organisme, souche, espèce, organite, teneur en
  GC, dinucléotide, profondeur, profondeur de séquençage (Read depth), couverture, code
  génétique, identité, seuil, requête / sujet (BLAST query / subject), hit, homologie,
  orthologue, paralogue, co-orthologue, inparalogue, colinéaire, meilleur hit réciproque.
- **Figure words**: piste, sous-piste (lane), emplacement (track slot), empilement (stack),
  anneau, règle graduée (ruler), graduation principale / secondaire (Excel French ticks), entrée
  de légende, étiquette, libellé (caption), trait (line), contour (stroke), remplissage (fill),
  pastille (swatch), calque, marge (padding), écart (gap), décalage.
- **Command words** of Microsoft and Apple French: Enregistrer, Charger, Ouvrir, Fermer, Copier,
  Coller, Dupliquer, Supprimer (Delete), Retirer (Remove), Renommer, Effacer, Vider le cache,
  Réinitialiser, Restaurer, Appliquer, Réessayer, Ignorer (Skip), Conserver, Afficher / Masquer,
  Activer / Désactiver, Réduire / Développer, Monter / Descendre, Parcourir, Importer /
  Exporter, Télécharger, Sélectionner (Select) / Choisir (Choose), Faire glisser, Faire pivoter;
  Settings are paramètres (Windows French), never réglages.
- **Undo / Redo is Annuler / Rétablir**, the pair of Office, macOS, LibreOffice and Google Docs
  French. Annuler is also Cancel, as in every French program: the header Undo is an icon button
  whose tooltip names the step (Annuler : {action}) and never sits beside a dialog Annuler.
  Rétablir means Redo only (Reset to default is Revenir à la valeur par défaut, Revert is
  Revenir à la valeur d’origine).
- **Loanwords kept**: hit, skew (GC skew, le skew positif), orthogroup, open source, thread /
  multithread, cross-origin, Bitscore, E-value, the names of 1.1, tRNA and rRNA (prose ARNt, ARNr)
  and the unit bp (prose pb).
- **One word, one concept**: élément is Feature only (a list item is entrée, a thing on the
  canvas is objet); entrée is Input and Entry, saisie a typed value; enregistrement is Record
  only, never the act of saving (write enregistrer, enregistrée, en cours d’enregistrement);
  règle is Rule, règle graduée Ruler; requête is Query, demande Request; cible is a target,
  sujet the BLAST subject; tableau is a UI table, table de couleurs the color-table file;
  recadrage (Crop), rognage (Clip) and écrêtage (value clipping) differ; importer serves Import
  and Upload, charger Load.
- **Never**: caractéristique or fonctionnalité (Feature), sauvegarder (Save), téléverser,
  record, track or feature in French text, optionnel (facultatif), invalide (non valide),
  supporter (prendre en charge), éditer for modify (modifier; éditeur and édition stay for the
  Editor and the edit mode), éligible (admissible), pop-up, Veuillez, S’il vous plaît.

### 7.2 Buttons and sentences

- Register: vous everywhere, never tu. Buttons, menu commands, links and button tooltips:
  infinitive, no final period (Enregistrer la session, Charger un exemple, Copier le
  diagnostic). Checkboxes follow the English form (Afficher la profondeur, Échelle
  logarithmique). Headings, field labels, column headers, tabs and options: noun phrases without
  an initial article (Couleur de remplissage, Mode d’exécution LOSAT).
- Instructions (notices, remedies, help addressed to the user): imperative with vous
  (Sélectionnez un fichier.); a help text that describes its control is in the third person
  present (Affiche la ligne de l’axe de profondeur.).
- Work in progress: action noun and …, no verb (Chargement de la session…, Génération du
  diagramme…, Réinitialisation…); Saving session… is Session en cours d’enregistrement… (7.1);
  in a sentence, … en cours. (Annulation ou rétablissement en cours.)
- Completion: the noun and its past participle (Session enregistrée. Alignement appliqué : …).
  Errors: what failed, then what to do, as two sentences; Could not X is Impossible de X
  (Impossible d’importer les annotations : {reason}), X failed is X a échoué or Échec de X,
  Retry after it finishes is Réessayez une fois l’opération terminée.
- Dialog questions: infinitive (Retirer tout le fichier ?). Negation is always ne … pas. A
  control or section named in a sentence takes its French label in guillemets: Cliquez sur
  « Générer », dans « Annotations de régions ».

### 7.3 Punctuation and spacing

- U+00A0 before a colon, in a label too (Paramètre : {path}, Libellé :); U+202F before a
  semicolon, question mark or exclamation mark (Retirer ?); U+00A0 inside guillemets
  (« {name} »). No space before a period, comma, ellipsis or closing parenthesis.
- Quotes are guillemets with “ ” nested. An English quoted "{name}", “{name}” or '{name}'
  becomes « {name} »; a bare placeholder stays bare. A quoted option name takes the French label
  ('Bar' draws a simple line → « Barre » trace une ligne simple).
- Apostrophe: the typographic ’ (U+2019) in every string (l’axe, S’applique), as Apple French
  and French typography write it, and it needs no escaping in JavaScript; ASCII ' only in code.
- Parentheses as in English, no inner space: Région (facultatif). (facultatif)
  is invariable and never agrees with its noun. Ranges keep the source en dash
  ({inner}–{outer} px); in prose from … to … is de … à …. The slash keeps the source spacing:
  Organisme / souche, GenBank/DDBJ.
- Numbers and units keep the en-US form the app prints (0.65 ×R, 1,024 px, 100 bp, 20px, 65%):
  never 0,65, 1 024 or 65 %, because the formatter writes them and users type them back. One …
  (U+2026), never three dots. The em dash keeps the source spacing.
- The code joins two sentences with an ordinary space. Abbreviations take a period: min., max.,
  séq. (Identité min. (%)); ID is invariable (des ID).

### 7.4 Words and capitals

- Sentence case for every label, as Microsoft and Apple French UIs write it: first word
  capitalized, the rest lower case except names and identifiers (Emplacements de piste
  personnalisés, Paramètres de session conservés).
- Capitals keep their accents: Éléments, Étiquettes, À propos, Édition de la disposition. An
  English key in capitals stays in capitals, with accents: SCHÉMAS PRÉDÉFINIS, COULEURS PAR
  DÉFAUT, RÈGLES SPÉCIFIQUES, PRIORITÉ DES QUALIFICATEURS, AJOUTER UNE CLÉ, DÉSACTIVÉ. Latin
  script keeps the CSS `uppercase` and `tracking-*` rules, so a sentence-case key that CSS
  capitalizes stays in sentence case.
- Concepts the English capitalizes (Session, Result, Legend) are lower case in prose (la
  session, le résultat) unless they name a control (7.2); gbdraw stays lower case. Compounds
  follow the glossary words, head noun first: Record label is étiquette d’enregistrement, Depth
  TSV is TSV de profondeur, Feature type is type d’élément.

### 7.5 Placeholders and plurals

- A placeholder value has no known gender (record IDs, file names, setting paths, track names,
  user-typed names). Put a French noun before it and make articles, adjectives and participles
  agree with that noun: le fichier {file}, l’enregistrement {record}, la piste {track}, le
  paramètre {setting}, la palette « {palette} ». Never elide or inflect against the placeholder
  itself: de {name}, que {name}, never d’{name}. A placeholder holding a translated phrase
  ({action}, {reason}) follows a colon: Annuler : {action}.
- None is one shared key and is always Aucun (legend position, label mode,
  color); a context key is needed before a feminine Aucune can exist.
- Masculine: enregistrement, élément, fichier, résultat, hit, TSV, diagramme; feminine: piste,
  étiquette, entrée, séquence, correspondance, paire, règle, session; no article on CDS.
- `tn(count, one, other)` shows `one` when `Intl.PluralRules('fr')` returns one, for 0 and for
  1, so the singular must read right for 0 (0 enregistrement); 1000000 returns many, which `tn`
  treats as other. Keep {count} in the one form whenever the English one form has it; where it
  has none (From this color rule:), the calls pass 1 or more and French omits it too.
- A key outside `tn` that writes "(s)" in English writes none in French: the count follows a
  colon (Correspondances sans clé d’élément, non appliquées : {count}.).
- Worked examples (one / other):
  1. {count} record → {count} enregistrement / {count} enregistrements.
  2. {count} LOSAT result is cached and reused when the inputs and settings match. → {count}
     résultat LOSAT est en cache et sera réutilisé si les entrées et les paramètres
     correspondent. / {count} résultats LOSAT sont en cache et seront réutilisés si …
  3. Input files: {count} entry of {path} could not be read and was dropped. → Fichiers
     d’entrée : {count} entrée de {path} n’a pas pu être lue et a été ignorée. / Fichiers
     d’entrée : {count} entrées de {path} n’ont pas pu être lues et ont été ignorées.
  4. {current} of {count} match → {current} sur {count} résultat / … résultats.
  5. From this color rule: / From these {count} color rules: → De cette règle de couleur : / De
     ces {count} règles de couleur :
- A list placeholder agrees the same way: La paire {pairs} est… / Les paires {pairs} sont….

### 7.6 Length

- French labels run 15–30% longer than English, a button with an object about twice as long
  (Enregistrer la session); take the short form the glossary gives.
- Rail short names (52 px buttons, 10 px semibold, no wrapping: about nine Latin letters):
  Entrée, Comparer, Avancé, Régions, Disposition (about 52 px, the limit), Police, Titres,
  Légende, Couleurs, Éléments, Étiquettes (about 48 px), Axe, GC / skew, À propos, Générer.
- The four header buttons (Réinitialiser les paramètres, Enregistrer la session, Charger une
  session, Charger un exemple) are about 1.7× the English width: check them at 1,024 px in the
  PR-B screenshots and fix an overflow in CSS, since the same keys name them in prose.
- Short forms: Monter / Descendre, Tout effacer, min. / max., Auto, ID, séq., no initial article
  in a label. Help texts and notices may be as long as needed.

## 8. Portuguese (`pt`)

### 8.1 Vocabulary

- **Norm.** Brazilian Portuguese spelled to the 1990 Orthographic Agreement
  (sequência, ideia, sub-região, pré-visualização), one catalog under `pt`
  (D-03). No slang or Brazilian-only idioms, so a reader in Portugal is not put
  off. Where the norms use different words, take the Brazilian one: arquivo
  (not ficheiro), tela (ecrã), salvar (guardar), baixar (transferir), excluir
  (eliminar), registro (registo), seção (secção), aba (separador), usuário
  (utilizador), aplicativo (aplicação), configurações (definições), mesclar
  (unir), organela, nucleotídeo (organelo, nucleótido), cepa (estirpe), fita
  (cadeia), and the gerund for progress (Salvando…, not A guardar…). Where one
  word serves both norms, use it: ignorar (not pular), pré-visualização (not
  prévia).
- **Field terms**: genoma, gene, sequência, registro, número de acesso, locus
  (pl. loci), qualificador, anotação, fita (fita direta, fita + / fita −),
  complemento reverso, homologia, ortólogo, parálogo, similaridade, identidade,
  profundidade de leitura, cobertura, conteúdo de GC, código genético,
  organismo, espécie, cepa, produto, localização, definição.
- **Diagram terms**: diagrama, trilha (Track), subtrilha (Lane), banda (Band),
  faixa (Ribbon), anel (Ring), eixo, escala, marca de escala (Tick; marca
  principal / secundária), régua, barra de escala, legenda, rótulo (Label),
  camada, sobreposição (Overlay), fundo / primeiro plano (Underlay /
  Foreground), realce, tela (Canvas), margem (Padding), espaço (Gap),
  espaçamento (Spacing), intervalo, posição / posicionamento / localização
  (Position / Placement / Location).
- **Command words** (Microsoft and Google pt-BR): Carregar (Load), Enviar
  (Upload), Baixar, Redefinir (Reset), Restaurar, Reverter, Remover, Excluir
  (Delete), Limpar (Clear), Mesclar, Ignorar (Skip), Tentar novamente,
  Procurar (Browse), Escolher, Selecionar, Buscar, Recolher, Habilitar /
  Desabilitar, Ativado / Desativado (On / Off), Gerar.
- **English kept** where the field says it in Portuguese: feature (f.), slot,
  layout, hit, subject, skew (with GC skew, AT skew), cache, pop-up, zoom,
  orthogroup (m.), query, thread (f.), E-value, Bitscore. Write them as
  Portuguese nouns: lower case in a sentence, plural -s, an article of the
  stated gender, no italics or quotes. Never mix one with a translation of it
  (no característica beside feature, no consulta beside query).
- **Never use**: característica or elemento for Feature; padrão for Pattern
  (padrão is Default, Pattern is expressão); faixa for Track or Lane; cadeia
  for strand; pb or kpb for bp, kb; deletar, resetar, setar, printar, upar;
  ficheiro, ecrã, guardar, registo, secção, utilizador; por favor.

### 8.2 Buttons and sentences

- **Register**: você, always implicit; never tu or o senhor.
- **Buttons, menu commands, action labels, tooltips of actions**: infinitive
  (Salvar sessão, Gerar diagrama, Desfazer, Refazer, Renomear), never an
  imperative (Salve) or a noun (Salvamento). Labels, headings, tabs and
  options are noun phrases with no final period (Cor da linha).
- **Instructions** in help, notices, dialogs, errors: imperative of você
  (Selecione um arquivo. Clique em Gerar. Deixe em branco para usar o valor
  automático.). **Descriptions** ("Shows…") take the present: Mostra…,
  Desenha….
- **Progress**: gerund, no article in a short status (Salvando sessão…,
  Carregando sessão…, Gerando diagrama…, Aplicando…). **Completion**: a
  participle agreeing with its noun (Sessão salva. Alinhamento aplicado.).
- **Errors**: "Could not" is Não foi possível + infinitive; "Failed to" in a
  heading is Falha ao + infinitive; "X failed." is X falhou; cause and remedy
  follow as sentences. Busy reasons: Gerando o diagrama. Tente novamente
  quando a geração terminar. A dialog title asks with the infinitive (Remover
  o arquivo inteiro?); Deseja …? only in body text.
- **A control named in a sentence** takes its translated label with that
  label's capitals, no quotes, after em: clique em Gerar diagrama, use Salvar
  sessão, em Configurações.

### 8.3 Punctuation and spacing

- Quotes “ ” (U+201C, U+201D), nested ‘ ’; never straight quotes. A word the
  English puts in single quotes ('Bar', '{selector}') also takes “ ”; a quoted
  option label is translated inside them (“Barra”).
- A label colon has no space before and one after: Campo: {label}. Inside a
  sentence, continue after a colon in lower case unless a name, a control
  label or a placeholder follows.
- No space before ; ? ! (unlike fr); no ¿ ¡. Sentences are separated by one
  space (D-06).
- Parentheses, slashes and ranges as in English: Região (opcional), Espaço
  entre registros (px), Organismo / cepa, Linhas 3–4 (en dash, no spaces).
- Numbers keep the source's decimal point and grouping (0.65, 1,000; never
  0,65 or 1.000); units and percent keep the source spacing (100 bp, 20px,
  65%, 0.65 ×R). No ordinals (1º): primeira linha. One ellipsis character ….

### 8.4 Words and capitals

- **Sentence case** for every label, button, heading, tab and option, as in
  Microsoft and Google pt-BR: only the first word and names are capitalized.
  Save Session → Salvar sessão; Custom Track Slots → Slots de trilha
  personalizados; Pairwise Comparisons → Comparações par a par.
- ALL CAPS headings (PRESET SCHEMES) are written in sentence case (Esquemas
  predefinidos); the CSS draws the capitals as for English.
- Names keep their form: gbdraw (lower case), LOSAT, GenBank, BLAST, Tuckin,
  Middle, Spreadout. Language names in prose are lower case (português); the
  menu endonym Português is not translated.
- Compounds are built from glossary words, head noun first: Record label →
  rótulo do registro; Feature type → tipo de feature; Linear track slot Lane
  gap → Espaço entre subtrilhas do slot de trilha linear. de + article for a
  specific thing (da legenda), bare de for a kind (tabela de cores).
- Adjectives agree (posição personalizada); a standalone option label is
  invariable masculine: Automático, Nenhum, Ativado, Desativado, Padrão.

### 8.5 Placeholders and plurals

- Placeholders are copied verbatim and may move. A value's gender is unknown,
  so put a noun before it and make articles and adjectives agree with that
  noun: o arquivo “{file}”, o registro {record}, a fonte {names}, a
  configuração {setting}, o valor {value}, a trilha {track}. Never put o, a,
  do, da, no, na, ao or à straight before such a placeholder: Mesclar em
  “{name}”, not no {name}. A user-typed name the English quotes is “{name}”.
- `Intl.PluralRules('pt')` selects `one` for 0 and 1 (D-05): write `one` in
  the singular and make it read right with 0 (0 registro is the CLDR pt form).
  A `one` text without {count} (From this color rule:) is shown for a single
  item only. The count is digits; noun, article, adjective and verb agree with
  the form, and so do lists ({items}, {names}): o item {items} / os itens
  {items}.
- Worked examples (one / other):
  1. `{count} record` → {count} registro / {count} registros
  2. `{count} needs a choice` → {count} precisa de uma escolha / {count}
     precisam de uma escolha
  3. `{count} chosen feature keeps its current direction.` → {count} feature
     escolhida mantém a direção atual. / {count} features escolhidas mantêm a
     direção atual.
  4. `Keep your {count} changed color from the "{palette}" palette ({keys})?`
     → Manter {count} cor alterada da paleta “{palette}” ({keys})? / Manter
     {count} cores alteradas da paleta “{palette}” ({keys})?
  5. `Legend: entry {items} … was dropped.` → Legenda: o item {items} tinha
     uma cor que o aplicativo não aceita e foi descartado. / Legenda: os itens
     {items} tinham uma cor que o aplicativo não aceita e foram descartados.

### 8.6 Length

- Portuguese runs 20–35% longer than English. Rail short names fit the
  52 px button at 10 px in about 9 letters: Entrada, Features, Cores,
  Rótulos, Legenda, Layout, Eixo, GC / skew, Regiões, Títulos, Comparar (the
  verb; Comparação does not fit), Avançado, Fonte, Sobre.
- Header and bar buttons stay within about 1.3× the English (Salvar sessão,
  Carregar exemplo, Gerar diagrama). Redefinir configurações and Tentar
  novamente have no shorter form; the PR screenshots at 1,024 px check them.
- Inside 10 px rows use the short glossary form: marca principal (not marca
  de escala principal), fundo (not plano de fundo), Automático (never Auto).
  Abbreviate only where English does (Mín. / Máx.). Help texts wrap.
