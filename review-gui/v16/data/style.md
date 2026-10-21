# Web UI translation style: ja, zh-Hans, zh-Hant

Revision 3, 2026-10-21. Companion of `I18N_GLOSSARY.tsv`. Translators and LLM drafts
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
- **The Language menu entries** (English, 日本語, 简体中文, 繁體中文) are
  endonyms and are not translated; the menu's own label is (言語 / 语言 /
  語言). Gallery example titles from `tutorial-data/manifest.json` stay
  English too.
- **Names fixed in English by PD-OI-093 and PD-OI-105** stay English in the
  English UI only. Other locales translate them with the glossary (D-11).

### 1.2 Numbers, units, symbols

- ASCII digits only: no full-width digits and no kanji or Chinese numerals for
  quantities. Decimal point and digit grouping stay as the English source
  shows them, because the app formats numbers with `en-US`.
- Units stay Latin and keep the source's spacing between number and unit:
  `100 bp`, `1kb`, `20px`, `0.65 ×R`, `8 pt`, `10x`. Never convert a unit.
  Time words are words, not units: 30 seconds is 30秒 in all three locales.
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
- A `tn` key has one string per locale: `ja`, `zh-Hans` and `zh-Hant` have no
  plural forms.
- One key is one sentence or label. Never build a sentence from translated
  fragments; ask for a key with placeholders instead. A sentence with an
  inline element (`<strong>`, a swatch, a bound name) is one key; the
  `<i18n-text>` component places the element at its `{name}` slot (D-04).
- Two sentences shown together are two keys. The code joins them with a
  language-dependent separator: a space in English, none after 。 in `ja` and
  `zh-*` (D-18). Do not add a trailing space to a value.
- Quotes follow the English: a `{name}` that the English key quotes is quoted
  in the target language's quotes (`ja` 「」, `zh-Hans` “ ”, `zh-Hant` 「」);
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
