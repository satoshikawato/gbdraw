# Web UI translation style: ja, zh-Hans, zh-Hant, ko, es, fr, pt, id, vi, th, de, it, pl, tr, ru, uk

Revision 11, 2026-10-31 (German, Italian, Polish, Turkish, Russian and Ukrainian: §12–§17; plural objects for few/many locales and the de, it, pl, tr plural and separator lists: §1.3; the cited Auto option in de: §12.2; tr dialog titles and singular `one` forms: §15.2, §15.5; ru and uk time abbreviations: §1.2, §16.3, §17.3; ru tab names: §16.2; uk Load / Download / Upload sources and plural worked pair 5: §17.1, §17.5; ru and uk in the separator and plural lists: §1.3). Revision 10 was the 2026-10-22 revision (cross-locale fixes: the drawing mode versus topology and dimension words with their track in every
locale, rail names wrap only between words in every locale, ja/zh comparative-genomics terms with no established use stay English,
ja Attach, zh Either and axis line, es/fr/pt counts and exact matches, fr stroke and Linear scale, id/vi/th record selector
and feature wording: §1.5, §2, §3, §4, §5.1, §6, §7, §8, §9, §10, §11. Revision 9: counted sentences are `tn` pairs, and
track, mode, control and setting names are no longer placeholders inside sentences (§1.3, §7.5); fr Feature and MB confirmed
by the Owner (P-07): §7.1, §7.3; this file and the glossary, not the Product Contract, hold the term rules: §2.1, §3, §9.1.
Revision 8 was the ko, es, fr and pt round-2 revision; Revision 7 the round-1 term audit of Indonesian, Vietnamese and Thai:
§9, §10, §11). Companion of `I18N_GLOSSARY.tsv`. Translators and LLM
drafts follow this file and the glossary together. The glossary decides words;
this file decides how sentences are written around them. Where the two
disagree, the glossary wins and this file is corrected.

When keys are translated: a pull request that adds or changes UI text adds the
English key only; a catalog keeps the old translation of an edited or removed
text as a stale key, and the translation batch before each release fills the
new keys and drops the stale ones. The batch follows
`.claude/skills/translate-web-ui/SKILL.md` with `node tools/i18n/sync.mjs`
(`status`, `export`, `apply`), and only the release tier of `Tests` fails while a
catalog lacks a key of the sources or holds a stale one
(`node tools/i18n/sync.mjs status --require-complete`). Pull-request and
ordinary `dev` runs fail on a non-string or empty value or a placeholder
mismatch, and print missing and stale keys without failing. Each locale has a lint
list, `tests/web/i18n-lint/<locale>.tsv`, run by
`tests/web/i18n-catalogs.test.mjs` as §9.7 describes. Translations are drafts;
nothing states that they are reviewed.

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
  (`above_feature`, `exclude_matching`, `overflow=error`),
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
- **The Language menu entries** (English, Bahasa Indonesia, Deutsch, Español,
  Français, Italiano, Polski, Português, Tiếng Việt, Türkçe, Русский,
  Українська, ไทย, 日本語, 简体中文, 繁體中文, 한국어; the menu order) are
  endonyms and are not translated; the menu's own label is (言語 / 语言 / 語言 /
  언어 / Idioma / Langue / Idioma / Bahasa / Ngôn ngữ / ภาษา / Sprache /
  Lingua / Język / Dil / Язык / Мова for ja / zh-Hans / zh-Hant / ko / es /
  fr / pt / id / vi / th / de / it / pl / tr / ru / uk). Gallery
  example titles from `tutorial-data/manifest.json` stay English too.
- **Names fixed in English by PD-OI-093 and PD-OI-105** stay English in the
  English UI only. Other locales translate them with the glossary (D-11).

### 1.2 Numbers, units, symbols

- ASCII digits only: no full-width digits and no kanji or Chinese numerals for
  quantities. Decimal point and digit grouping stay as the English source
  shows them, because the app formats numbers with `en-US`.
- Units stay Latin and keep the source's spacing between number and unit:
  `100 bp`, `1kb`, `20px`, `0.65 ×R`, `8 pt`, `10x`. Never convert a unit (bp, kb and Mb
  stay in every locale; `fr` writes byte units in French, §7.3).
  Time words are words, not units: 30 seconds is 30秒 in `ja`, `zh-Hans` and `zh-Hant`,
  and `ru` and `uk` write мс, с and мин / хв (§16.3, §17.3).
- Keep `→ ← · × + − %` and the en dash of ranges (`Rows 3–4`) as written.
- Ellipsis is one `…` (U+2026) whether the source writes `...` or `…`, both
  for "more input follows" (Choose records…) and for progress (Saving
  session…).

### 1.3 Placeholders and keys

- The English source text is the key (D-02). There is no English catalog:
  `t('Save Session')`, `tn(count, one, other)` (the English text of each form is a key),
  `tx(context, text)` (key `context + '\u0004' + text`, for homonyms), `msg()`
  to keep a text until `tm()` translates it where it is shown. A catalog
  `<locale>.js` maps those keys to plain strings; only the `other` key of a
  `tn` pair may map to plural forms (below).
- `{name}` placeholders are copied verbatim: the same set as the English key,
  never translated, never added or dropped. They may move inside the
  sentence.
- A `tn` key has one string in `ja`, `zh-Hans`, `zh-Hant` and `ko`, which have
  no plural forms (`Intl.PluralRules` returns `other`), and a `one` and an
  `other` string in `es`, `fr`, `pt`, `de`, `it` and `tr` (§5.5 to §11.5,
  §12.5, §13.5, §15.5). `th`, `id` and `vi`
  have no plural forms either: `Intl.PluralRules` selects `other` for every
  count, so their `one` rows are never shown (they carry the `other` text or a
  natural singular; §9.5, §10.5, §11.5).
- A `tn` key's translation is one string, except that in a language whose
  plural rule has `few` or `many` categories (Polish, Russian, Ukrainian)
  the `other` key holds an object of the CLDR categories the language uses
  (`pl`: few, many, other; `ru` and `uk`: one, few, many, other),
  with `other` required, for example
  `{"one": "…", "few": "…", "many": "…", "other": "…"}`. The `one` key is
  then the text for exactly 1 and stays a string; the object's `one` is the
  text of the other counts the rule calls `one` (21, 31, … in Russian and
  Ukrainian), so it must not read as "this one thing". A category the object
  leaves out shows its `other` text.
- A sentence that counts things is a `tn` pair; no key writes a plural as "(s)".
  A track, mode, control or setting name is never a placeholder inside a
  sentence: each name has a whole-sentence key (Hide GC content?, No Circular
  Result yet. …), so the sentence can take the article, case or particle the
  name needs.
- One key is one sentence or label. Never build a sentence from translated
  fragments; ask for a key with placeholders instead. A sentence with an
  inline element (`<strong>`, a swatch, a bound name) is one key; the
  `<i18n-text>` component places the element at its `{name}` slot (D-04).
- Two sentences shown together are two keys. The code joins them with a
  language-dependent separator: a space in English and in `ko`, `es`,
  `fr`, `pt`, `id`, `vi`, `th`, `de`, `it`, `pl`, `tr`, `ru` and `uk`, none after 。 in
  `ja` and `zh-*` (D-18). Do
  not add a trailing space to a value.
- Quotes follow the English: a `{name}` that the English key quotes is quoted
  in the target language's quotes (`ja` 「」, `zh-Hans` “ ”, `zh-Hant` 「」,
  `ko` ‘ ’, `es` and `pt` “ ”, `fr` « », `id` straight double
  quotes "…", `vi` “ ”, `th` “ ” with nested ‘ ’);
  one the English leaves bare stays bare (D-25). English double quotes never
  appear in `ja` or `zh-*` prose.

### 1.4 Words and capitals

- Use the glossary rendering for every English term it lists, in every key,
  including compounds: Record label is レコードラベル because Record is
  レコード and Label is ラベル. A row tagged `[review]` marks a judgment call;
  it is still the value to use until the glossary changes.
- Do not coin a target-language term. A domain term with no established use
  in that language's molecular biology or bioinformatics writing stays in
  English (lower case in running text unless it is an acronym or a name); an
  established native term or loanword stays. Example: reciprocal best hit.
- Title Case, ALL CAPS headings (PRESET SCHEMES, DEFAULT COLORS, SPECIFIC
  RULES, LABEL FILTERING, QUALIFIER PRIORITY, LABEL GEOMETRY, ADD KEY) and
  letter spacing carry no meaning in CJK; translate them as ordinary
  headings. The same rules cover `th` (Thai has no letter case) and `vi`
  (capitals with stacked tone marks collide). Outside the figure, the CJK,
  `th` and `vi` rules in `gbdraw/web/index.html` turn off `uppercase` and
  `tracking-*`.
- A word that the source capitalizes as a name of a control (Generate,
  Apply, Save Session) is written as that control's translated label, so a
  reader can find it: 「生成」をクリック / 点击“生成” / 點選「產生」.

### 1.5 Length

- CJK glyphs are about twice as wide as Latin letters at the same size, and
  the CJK CSS raises 10–11 px text to 12 px. Count width, not characters.
- Rail short names (52 px entries; 10 px, CJK 12 px, Thai 12 px): prefer a
  name that fits on one line (CJK at most 4 full-width glyphs, 5 at most;
  Latin about nine letters). A name never breaks inside a word in any
  locale: CJK words keep their glyphs together (`keep-all`) and Latin words
  are not hyphenated (no soft hyphens in rail values). A name that does not
  fit takes two lines, breaking only between words (GC / sesgo, Acerca de),
  and its entry grows. A single word wider than the entry stays on one line
  and may run into the rail column's padding (58 px inside the column at most;
  ja フィーチャー uses proportional kana to fit); a word wider than that needs a
  shorter rail key (`tx('rail', '…')`). Rail keys: `rail\u0004Features`,
  `rail\u0004Layout`, `rail\u0004Labels`; a catalog value equal to the
  section title is correct when it fits. Header and bar buttons: no wider than
  the English label. Option labels inside 10 px rows: as short as the
  glossary form.
- Help texts and notices may be as long as needed; they wrap.
- `es`, `fr` and `pt` labels run 15 to 35% longer than English; their length
  rules are §6.6, §7.6 and §8.6.
- `id` and `vi` short verbs and phrases can run up to twice the English width;
  see §9.6 and §10.6.
- Thai UI text has a 12 px floor, also for the rail names, the hints and the
  header buttons, and line-height 1.5 on the rail (§11.6).

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
- **Diagram mode versus topology**: a record that belongs to the Circular or
  Linear diagram is written 環状図の… / 直鎖状図の… (環状図のレコード,
  直鎖状図の配列{n}). 環状レコード is only the topology (Circular record
  checkbox: the record is circular); 直鎖状配列 is never used for a Linear
  diagram's sequence.
- **Comparative-genomics terms (§1.4)**: オーソログ, パラログ, 相互ベストヒット
  (reciprocal best hit; U-Tokyo bioinformatics textbook) and グラフのエッジ
  (エッジ) are established and stay. Terms with no established Japanese use
  stay English in lower case: near-reciprocal hit, coortholog, inparalog. Do
  not write 準相互ヒット, 共オーソログ or インパラログ. Compounds of an
  established word keep the Japanese pattern: レコード内パラログ, RBHシード.
- **Attach / Attached** (a Depth TSV to a series or record) is 割り当て /
  割り当て済み in every key, not 関連付け.
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

- フィーチャー on the rail keeps one line at every window size (never
  フィー / チャー, Owner 2026-10-22); no shorter rail form.

## 3. Simplified Chinese (`zh-Hans`, mainland usage)

Open zh-Hans and zh-Hant terms (glossary rows marked [review]) are decided by
zh reviewers; this file and the glossary record each decision.

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

### 3.6 Circular and Linear as the name of a drawing

- When Circular or Linear names the diagram mode and is followed by a record,
  sequence, file, region, scale or row, say 图谱: 线性图谱的记录间隙,
  环状图谱的记录标签, 序列{n}在线性图谱中的记录行, 线性图谱的比例尺,
  清除线性图谱的文件. 线性记录, 线性序列, 线性文件, 线性比例尺, 环状区域,
  环状序列 read as molecule topology (or a mathematical linear scale) and are
  not used.
- A bare 环状 stays for the topology of a record (环状记录 for the key
  `Circular record`), for the mode name itself, and for things that are
  ring-shaped (环状轨道, 环状比较环).

### 3.7 Axis and line, not an outline

- The axis line and the intron line are lines, not outlines: 轴线颜色 /
  轴线宽度, 线条颜色 / 线条宽度 (as zh-Hant 軸線, 線條). 描边 stays for the
  outline of blocks, features and legend entries.

### 3.8 Established Chinese terms for comparative genomics (§1.4)

- Kept: 直系同源 and 旁系同源 (ortholog, paralog; zh.wikipedia, Chinese
  textbooks), 命中 for a BLAST hit, 共线性, 锚点, 区块.
- English, lower case in running text (no established Chinese term found):
  reciprocal best hit, near-reciprocal hit, coortholog, inparalog. 双向最佳命中,
  近似双向最佳命中, 共直系同源 and 内旁系同源 are not used. A Chinese gloss may
  follow once in a help text if reviewers want one (双向最优匹配 for BBH appears
  in one Chinese course page; zh.wikipedia writes 后重复 for inparalog).
- Sentence or label start keeps the English capital (Reciprocal best hit：…);
  Chinese has no plural, so the option Reciprocal best hits is Reciprocal best
  hit.

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

### 4.7 Circular and Linear as the name of a drawing

- The rule of §3.6 holds with 圖譜 (the glossary's Diagram): when Circular or
  Linear names the diagram mode and is followed by a record, sequence, file,
  region, scale or row, say 線狀圖譜的… / 環狀圖譜的… (線狀圖譜的比例尺, not
  線狀圖 or 線狀比例尺). A bare 環狀 / 線狀 record, sequence,
  file or scale reads as molecule topology and is not used for the mode. A
  bare 環狀 stays for the topology of a record, the mode name itself, and
  ring-shaped things.

### 4.8 Dimension words of a track

- Width, Height and Radius of a GC, skew or depth track carry 軌道:
  深度軌道寬度, GC含量軌道半徑, GC skew軌道寬度, GC軌道高度 (zh-Hans 轨道
  already). A bare 深度寬度 / GC含量寬度 reads as the width of a value.
  Sentences that quote these labels quote the new wording.

### 4.9 Established Chinese terms for comparative genomics (§1.4)

- The rule of §3.8 holds: 直系同源 and 旁系同源, 命中, 共線性, 錨點, 區塊 stay;
  reciprocal best hit, near-reciprocal hit, coortholog and inparalog stay
  English in lower case, with the English capital at a label or sentence start.

## 5. Korean (`ko`)

Standard Korean (표준어) as South Korean software writes it: Microsoft and
Google Korean UIs for command words, Korean molecular-biology textbooks and
papers for the science. The Language menu entry is **한국어** (PD-OI-106).

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
- **Kept in Latin**, as Korean papers write them: GC skew, AT skew, skew,
  Bitscore, E-value, orthogroup, RBH, reciprocal best hit, near-reciprocal
  hit, coortholog, inparalog, locus tag, and every §1.1 name.
- **Feature** is written in English, lower case, also at the start of a label
  (feature 유형), and its particle follows the Korean reading: feature를,
  feature의. Never 피처 or 특징.
- **Loanwords** in the standard spelling (외래어 표기법) Microsoft uses:
  레코드, 트랙, 레인, 슬롯, 앵커, 히트, 쿼리, 블록, 리본, 링, 밴드, 마크,
  레이블, 라이선스, 메타데이터, 스키마, 런타임, 렌더러, 렌더링, 오버레이,
  오프셋, 팔레트, 프리셋, 캔버스, 레이어, 레이아웃, 세션, 섹션, 필드, 필터,
  패턴, 스택, 스텝, 유닛, 다이어그램, 다이뉴클레오타이드. Exception: 윈도우
  (sliding window), as bioinformatics writes it, not the standard 윈도. Add no
  loanword the glossary does not list. The glossary's 0.15.0 rows add
  these, which are allowed: 템플릿, 시드, 패널, 라이브러리, 인코딩, 파이프라인,
  워크플로, 해시, 플랫 파일, 포커스, 그라데이션, 픽셀, 인트론, 레플리콘,
  프리프린트, 캐시, 사이드바, 호스팅, 웹 앱, 리소스, 컨트롤, and the
  loanword parts of 파라로그, 오솔로그 and 갭 개시 수 (rows of
  `I18N_GLOSSARY.tsv`). The compound 소스 코드 (Source code) is allowed.
- **Track dimensions** carry 트랙: the Width, Height and Radius of a GC, skew
  or depth track are GC 함량 트랙 너비, GC skew 트랙 반지름, GC 트랙 높이,
  깊이 트랙 너비. Without 트랙 the label reads as the width of the GC content
  itself. A control name for a slot is 트랙 슬롯 {id} …, never 트랙 {id} …
  (선형 다이어그램의 트랙 슬롯 {id} 렌더러).
- **Command words** (Microsoft Korean): 저장, 열기, 닫기, 취소, 추가, 삭제
  (Delete), 제거 (Remove), 지우기 (Clear), 복사 (Copy), 복제 (Duplicate),
  붙여넣기, 편집, 적용, 선택, 모두 선택 해제, 표시 / 숨기기, 켬 / 끔, 활성화 /
  비활성화, 실행 취소 / 다시 실행, 되돌리기 (Revert), 복원 (Restore), 초기화
  (Reset), 다시 시도, 건너뛰기, 바꾸기, 이름 바꾸기, 병합, 가져오기 /
  내보내기, 불러오기 (Load), 찾아보기 (Browse), 고정, 도킹, 접기 / 펼치기,
  확대/축소, 맞춤 (Align); nouns 미리 보기, 대화 상자, 편집기, 도움말, 세부
  정보, 사용자 지정, 기본값, 재정의, 임계값, 정규식, 글꼴, 오류, 기록 (History).
- **Never**: 게놈, 시퀀스 (except 이스케이프 시퀀스, regular-expression escape
  sequences), 스트랜드, 어노테이션, 피처, 특징 or 기능 for Feature,
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
- **Recurring phrasings**: Leave blank for Auto is 비워 두면 자동으로
  설정됩니다; a label that does not fit is 배치할 공간이 부족합니다 (heading:
  레이블 배치 공간 부족), never 들어가지 않습니다.

### 5.3 Punctuation and spacing

- Half-width `. , : ; ? ! ( )`. A sentence ends with `.`; two sentences are
  separated by one space. A sentence whose English has no final period
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
- A Latin word takes the particle of its Korean reading (feature를, LOSAT을,
  px이).
- Prefer a known noun after the placeholder, so the particle attaches to
  Korean: {file} 파일을, {track} 트랙을, {names} 글꼴은. Use the paired form
  only when no such noun fits.
- **Plurals**: `Intl.PluralRules('ko')` always returns `other`, so there is
  one form. Translate the `other` key; the build script gives the `one` key
  the same string. No plural marker 들 on a counted noun; English
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
- Rail short names: at most 4 syllables, 2 preferred: 입력, feature, 레이블,
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

International Spanish (PD-OI-106): words a reader understands in Spain and in the
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
  groups (Revisar alineación, Restablecer alineación; Similarity alignment is
  alineación por similitud).
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
- **A slot id is a slot, not a track**: a label whose `{id}` is a track slot
  id names the ranura (Renderizador de la ranura de pista lineal {id}), so it
  cannot read as a track id.
- **A count of members is "Número de …"**: Número de genes del grupo de la
  query; a bare plural (Genes del grupo) reads as the list of genes.
- **Loanwords**: only hit (plural hits), query, subject, E-value, Bitscore,
  zoom, renderizador and the names of §1.1; also the English terms reciprocal
  best hit, near-reciprocal hit, coortholog, inparalog and orthogroup (el
  coortholog, masculine; plural -s); never track, feature, slot,
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
  Leave blank for auto: Déjelo vacío para usar el valor automático. Name a
  button as an action: haga clic en Restablecer. Do not mix the third-person
  indicative and the usted imperative in one sentence.
- Descriptions: third person present or impersonal se (Muestra la
  profundidad de lectura como una pista. Se aplica al generar.).
- An option, button or setting name is not the grammatical subject of a
  sentence: Con Automático, …; Al guardar, se almacenan …, not Automático
  elige … or Guardar almacena ….
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
- **An exact-match condition is parenthesized after the field**: Producto
  (coincidencia exacta): {product}, never Producto exacto, which reads as an
  adjective of the value.
- **Avoid su for "its"** where the sentence has no clear owner (hover hints,
  per-item settings): name the thing (la configuración de este elemento),
  because su is also "your" under usted.

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
  one does. A string without plural forms that counts things is a noun phrase
  and a colon: Etiqueta: {count}.
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
  52 px entry; a name that does not fit takes two lines (GC / sesgo takes two
  lines at 1,440 px): Entrada, Comparar, Avanzado, Regiones,
  Diseño, Fuente, Títulos, Leyenda, Colores, Etiquetas, Eje, GC / sesgo,
  Acerca de; Generar beside the rail. The rail short name of Features has its
  own key, `tx('rail', 'Features')`: Caract. stays its rail form (the section
  title Características, 15 letters, does not fit).
- Header buttons stay within about 1.3× the English label; check Restablecer
  configuración at 1,024 px. Options in 10 px rows use the glossary's short forms: Automático, En serie,
  Multihilo, Arriba / Abajo, Dentro / Fuera, En el eje, Fondo, Cinta /
  Curva. To shorten a label, drop articles first (Ancho de pista, not Ancho de
  la pista); never switch terms, and abbreviate with a period only where a
  screenshot shows an overflow.

## 7. French (`fr`)

Standard French of France, vouvoiement, French typography (PD-OI-106). The Language menu
entry is the endonym **Français**; the menu label itself is Langue.

### 7.1 Vocabulary

- **Biology and bioinformatics**, as French scientific writing has them: génome, séquence, gène,
  protéine, brin (brin +, brin −), complément inverse, annotation, élément (Feature; Owner P-07,
  2026-10-22, after the French GFF and annotation sources of the glossary row), type
  d’élément, enregistrement (Record), numéro d’accession, locus (invariable), qualificateur,
  produit (/product), localisation (Location), organisme, souche, espèce, organite, teneur en
  GC, dinucléotide, profondeur, profondeur de séquençage (Read depth), couverture, code
  génétique, identité, seuil, requête / sujet (BLAST query / subject), hit, homologie,
  orthologue, paralogue, co-orthologue, inparalogue, colinéaire, meilleur hit réciproque.
- **Figure words**: piste, sous-piste (lane), emplacement (track slot), empilement (stack),
  anneau, règle graduée (ruler), graduation principale / secondaire (Excel French ticks), entrée
  de légende, étiquette, libellé (caption), trait (line), contour (stroke), remplissage (fill),
  pastille (swatch), calque, marge (padding), écart (gap), décalage.
- **Stroke has two words**: trait for the axis and for lines (trait de l’axe,
  Épaisseur du trait: the line itself); contour for the outline of a shape
  (blocs, éléments, pastilles de légende). The axis stroke is never « contour
  de l’axe », which reads as an outline around the axis.
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
- **Loanwords kept**: hit, skew (GC skew, le skew positif), orthogroupe (prose; option values such as
  orthogroup stay as written), open source, thread /
  multithread, cross-origin, Bitscore, E-value, the names of 1.1, tRNA and rRNA (prose ARNt, ARNr)
  and the unit bp (prose pb). Byte units are French (§7.3).
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
  present (Affiche le trait de l’axe de profondeur.).
- Work in progress: action noun and …, no verb (Chargement de la session…, Génération du
  diagramme…, Réinitialisation…); Saving session… is Session en cours d’enregistrement… (7.1);
  in a sentence, … en cours. (Annulation ou rétablissement en cours.)
- Completion: the noun and its past participle (Session enregistrée. Alignement appliqué : …).
  Errors: what failed, then what to do, as two sentences; Could not X is Impossible de X
  (Impossible d’importer les annotations : {reason}), X failed is X a échoué or Échec de X,
  Retry after it finishes is Réessayez une fois l’opération terminée.
- A number that must be above 0 is strictement positif / positive; one that may be 0 is
  positif ou nul (positif alone includes 0 in France).
- An option or button name as the subject is quoted: « Avec « Enregistrer », … », or
  « Auto » masque …, so it is not read as an adverb. A setting that is Enabled is réglé(e) sur
  « Activé »; an instruction is Définissez … or Choisissez ….
- Linear / Circular qualifying a noun is du diagramme linéaire / circulaire (étiquettes du
  diagramme circulaire); « en mode circulaire » only for the app mode.
- The Linear scale (scale bar and ruler of the drawing) is « l’échelle du
  diagramme linéaire », never « l’échelle linéaire » (linear versus
  logarithmic scale, as the Échelle logarithmique checkbox shows); the same
  for the ruler: la règle graduée du diagramme linéaire.
- A count of members is « Nombre de … » : Nombre de gènes du groupe (requête).
- An exact-match condition: Produit (correspondance exacte) : {product}.
- The « mode » of a legend position is the diagram’s mode: des positions du
  mode du diagramme, not de son mode.
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
  never 0,65, 1 024 or 65 %, because the formatter writes them and users type them back.
  Byte units are the exception to §1.2: Mo, Go, ko (as octets), never MB, GB, kB; bp, kb
  and Mb stay. A file size in megabytes is Mo; a sequence length in megabases is Mb (Owner
  P-07, 2026-10-22). One …
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
  paramètre {setting}, la palette « {palette} ». At a sentence start or as subject: Le paramètre {setting} …, la piste « {track} ».
  Never elide or inflect against the placeholder
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
- A count outside `tn` (a short label such as {count} aligned, whose English has one form)
  follows a colon, so the participle need not agree with it (Alignés : {count}).
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
- Rail short names (52 px entries, 10 px semibold; about nine Latin letters on one line):
  Entrée, Comparer, Avancé, Régions, Police, Titres, Légende, Couleurs, Éléments, Axe,
  GC / skew, À propos, Générer. Disposition and Étiquettes have rail keys with a soft
  hyphen, Dispo­sition and Éti­quettes; they take two lines (Dispo- / sition,
  Éti- / quettes) where the one-line word does not fit (1,440 px). The section titles
  keep the plain words.
- The four header buttons (Réinitialiser les paramètres, Enregistrer la session, Charger une
  session, Charger un exemple) are about 1.7× the English width: check them at 1,024 px and
  fix an overflow in CSS, not with a shorter key, since the same keys name them in prose.
- Short forms: Monter / Descendre, Tout effacer, min. / max., Auto, ID, séq., no initial article
  in a label. Help texts and notices may be as long as needed.

## 8. Portuguese (`pt`)

### 8.1 Vocabulary

- **Norm.** Brazilian Portuguese spelled to the 1990 Orthographic Agreement
  (sequência, ideia, sub-região, pré-visualização), one catalog under `pt`
  (PD-OI-106). No slang or Brazilian-only idioms, so a reader in Portugal is not put
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
- **Connective and examples**: so is por isso, not então. e.g. is
  (por exemplo, …) inside a help sentence and ex.: in a placeholder or a short
  field hint. Nonnegative is não negativo / não negativa, agreeing with its
  noun, not maior ou igual a 0.
- **Diagram mode**: Circular / Linear qualifying a noun is diagrama circular /
  diagrama linear (registros do diagrama linear); modo stays for modo simples.
  The rule covers every noun, not only registros: a session or a comparison of
  a mode is sessão do diagrama circular, comparação do diagrama linear; sessão
  circular reads as a circular topology. The trilha / slot compounds of the
  glossary (slot de trilha linear) stay.
- **An option or button name is not the sentence subject**: Ao salvar a
  sessão, … são armazenados; Com Nenhum, …; never Salvar armazena … or
  Desativado as centraliza.
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
  space.
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
- **Track dimensions**: the Width, Height and Radius of a GC, skew or depth
  track carry trilha: Largura da trilha de conteúdo de GC, Raio da trilha de
  GC skew, Altura das trilhas de GC. A slot control name is Renderizador do
  slot de trilha linear {id}, not Renderizador da trilha {id}.
- **Counts and exact matches**: a value that counts members is Número de genes
  do grupo … (a bare plural noun reads as a list of genes); an exact-match
  condition is written Produto (correspondência exata): {product}, because
  exato beside a noun describes the noun.
- Adjectives agree (posição personalizada); a standalone option label is
  invariable masculine: Automático, Nenhum, Ativado, Desativado, Padrão.

### 8.5 Placeholders and plurals

- Placeholders are copied verbatim and may move. A value's gender is unknown,
  so put a noun before it and make articles and adjectives agree with that
  noun: o arquivo “{file}”, o registro {record}, a fonte {names}, a
  configuração {setting} (when {setting} is a setting's label; a {setting} that
  can be a button name or Esta alteração opens the sentence bare; at a sentence
  start: A configuração {setting} …), o valor
  {value}, a trilha {track}. Never put o, a,
  do, da, no, na, ao or à straight before such a placeholder: Mesclar em
  “{name}”, not no {name}. A user-typed name the English quotes is “{name}”.
- `Intl.PluralRules('pt')` selects `one` for 0 and 1: write `one` in
  the singular and make it read right with 0 (0 registro is the CLDR pt form).
  A `one` text without {count} (From this color rule:) is shown for a single
  item only. The count is digits; noun, article, adjective and verb agree with
  the form, and so do lists ({items}, {names}): o item {items} / os itens
  {items}.
- Worked examples (one / other):
  1. `{count} record` → {count} registro / {count} registros
  2. `{count} record awaiting a choice` → {count} registro aguardando uma
     escolha / {count} registros aguardando uma escolha
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
  52 px button at 10 px in about 9 letters; a rail name of two words takes two
  lines when it does not fit, and a single word must fit or get a rail key: Entrada, Features, Cores,
  Rótulos, Legenda, Layout, Eixo, GC / skew, Regiões, Títulos, Comparar (the
  verb; Comparação does not fit), Avançado, Fonte, Sobre.
- Header and bar buttons stay within about 1.3× the English (Salvar sessão,
  Carregar exemplo, Gerar diagrama). Redefinir configurações and Tentar
  novamente have no shorter form; the PR screenshots at 1,024 px check them.
- Inside 10 px rows use the short glossary form: marca principal (not marca
  de escala principal), fundo (not plano de fundo), Automático (never Auto).
  Abbreviate only where English does (Mín. / Máx.). Help texts wrap.

## 9. Indonesian (`id`)

Standard Indonesian in the register of Indonesian software and lecture notes
(WordPress, Proton and Google), not maximal Badan Bahasa coinage, never Malay.
Menu entry **Bahasa Indonesia**, menu label Bahasa (PD-OI-106).

### 9.1 Vocabulary

- **Provisional** (Owner decision D-05 of the Indonesian, Vietnamese and Thai
  translation work, 2026-10-21, written "id D-05" in the glossary; the vendors
  differ and the choice is to be reconsidered): Pengaturan (Settings), default, file, Unduh / Unggah,
  Urungkan (Undo), Atur ulang (Reset), Lanjutan (Advanced), penjajaran
  (alignment), strain, sekuens. The glossary marks each row.
- **Field terms** (spelling fixed): genom, gen, sekuens, nukleotida, asam
  amino, organel, lokus, anotasi, fitur (Feature; jenis fitur), untai (Strand;
  untai forward, untai +/−), kode genetik, nomor aksesi, kandungan GC,
  dinukleotida, kedalaman (Depth, depth coverage; kedalaman baca only for the
  read-depth table), identitas, kemiripan, ortolog, penjajaran (sequence
  alignment), kecocokan (comparison match), sirkular, linear (topology; the
  drawing mode is diagram linear / diagram sirkular), Kolinear.
- **Kept in English** : E-value, Bitscore, hit, query, subject, reverse
  complement, locus tag, upstream, downstream, RBH, reciprocal best hit,
  orthogroup, synteny, GC skew, AT skew, skew, qualifier, record, track, tick,
  ribbon, preset, offset, padding, thread, runtime and the names of §1.1:
  lower case in prose, no English plural ({count} record, Hit maksimum), no
  Indonesian affix. A help text may gloss once: qualifier (atribut fitur),
  Strain (galur), default (bawaan). Adapted loanwords, no italics: sesi,
  galeri, pratinjau, fitur, kustom, draf, filter, font, file, browser, pop-up,
  templat.
- **Diagram terms**: lajur (Lane), posisi (File card slot), garis skala (scale
  bar; Bar = Garis), penggaris (Ruler), celah (Gap), jarak (Spacing), pita
  (Band), tanda (annotation mark), latar belakang (Underlay), garis tepi
  (Stroke), isian (Fill), legenda, tumpukan, jangkar (anchor), perataan /
  Ratakan (Similarity group alignment), grafik (plot), gambar (Session
  drawing).
- **Command words**: Simpan, Muat (Load), Impor, Ekspor, Unduh, Unggah, Buat
  (Generate), Terapkan, Atur ulang, Urungkan, Ulangi, Batal, Hapus (Delete,
  Remove, Clear), Tambahkan, Pilih (Select, Choose), Telusuri (Browse), Cari
  (Search), Lewati (Skip), Kecualikan (Omit), Kembalikan (Revert), Periksa
  (Inspect); states Aktif / Nonaktif, Otomatis.
- **Never use** (checked in 9.7): Malay (muat turun, muat naik, fail, tetapan,
  pelayar, peranti, pangkalan data, perisian, boleh, bila, kerana); galat,
  peramban, prasetel, tetikus, utas; Setelan, Setel ulang, Tingkat lanjut,
  Download; sekuen, sirkuler, linier, unting, eror (variants or misspellings;
  one form is kept); komplemen balik, tag lokus, bilah skala, kueri, subjek,
  hulu / hilir, tumpang-tindih; khusus for Custom, rekaman for record, nilai
  E.

### 9.2 Buttons and sentences

- Register: formal, impersonal, subject omitted where possible; **Anda**
  (capital A) where a person is needed; never kamu, kau, gak, udah; no kami or
  saya.
- Buttons, menu commands, links: the bare base verb or its -kan imperative,
  sentence case, no period, no Silakan (Simpan sesi, Muat contoh, Terapkan);
  never nominalized (Penyimpanan sesi) and never a meN- form.
- Labels, headings, tabs, column headers: noun phrases, head noun first (Ukuran
  font label), impersonal (Pemberitahuan, not Beri tahu saya); a heading may be
  a pe-…-an noun (Pemfilteran label). A dialog title that asks: verb + object +
  ? (Hapus seluruh file?).
- Instructions: imperative without a subject (Pilih record. Klik Buat.
  Kosongkan untuk otomatis.). Descriptions: an active verb, subject implied
  (Menampilkan kedalaman.), or a di- passive. Progress: meN- verb + object +
  one … (Menyimpan sesi…). Completion: di- passive with a period (Sesi
  disimpan.).
- Failure: Gagal + verb (Gagal memuat sesi.) or Tidak dapat + verb (Tidak dapat
  membaca file.), never Gagal untuk or Tidak bisa; also … tidak tersedia /
  tidak valid / tidak didukung. Chips Kesalahan:, Peringatan:, Catatan:. Busy:
  Membuat diagram. Coba lagi setelah pembuatan selesai.
- Required, optional, recommended: Wajib ("… wajib diisi" in a sentence),
  Opsional and "(opsional)" after a field name, Direkomendasikan. Generate as
  an event: saat diagram dibuat, sebelum diagram dibuat, pembuatan diagram
  berikutnya; a cited button is Buat (klik Buat). An option that selects what a
  button does is cited with its label (klik Buat); "manually" for a generation
  started by that button is written with the button (dibuat dengan Buat), not
  secara manual.

### 9.3 Punctuation and spacing

- A period ends each sentence of help texts, notices, dialogs and errors;
  labels, buttons, headings, options and links take none; no "!". Straight
  double quotes wherever the English quotes a name or a UI word, also for its
  'Bar' style (Pilih "Garis" atau "Penggaris"); straight single quotes only
  for a gloss (EYD); never “ ”, ‘ ’ or « ».
- Colon: no space before, one after (Bidang: {path}); keep a label colon; no
  NBSP. Spacing around Latin terms, digits, placeholders and units is as in
  English; a placeholder takes no affix. Lists: no comma before dan / atau
  joining two items, a comma before the last dan of three or more (A, B, dan C;
  PUEBI); "&" is dan; e.g. is misalnya.
- ASCII digits, English point and grouping (0.65, 1,000), not 0,65 /
  1.000; units as the source (100 bp, 20px); bp, never pb.
- Parentheses, slashes and ranges as the source (Daerah (opsional), Baris 3–4,
  Kanan →); one … (U+2026), no space before it. di and ke as prepositions
  stand apart (di atas, ke bawah); the prefix di- attaches (dihapus). Hyphens:
  huruf besar-kecil; tumpang tindih is two words (KBBI); prefixes attach
  (nonaktif) except before English (multi-record); never di-reset.

### 9.4 Words and capitals

- Sentence case for every label, button, heading, tab, option and column header
  (Slot track kustom, Atur ulang pengaturan), never Title Case (Simpan Sesi).
  Slash alternatives follow the English capitals (Atas / Bawah, Organisme /
  strain).
- A key the English writes in capitals stays in capitals (TAMBAHKAN KUNCI,
  WARNA DEFAULT, SKEMA PRESET, ATURAN SPESIFIK; OFF is NONAKTIF).
- Lower case in prose: sesi, hasil, record, track, file, gbdraw. A cited
  control or section keeps its label and capital without quotes (klik Buat; di
  Slot track kustom); where the English quotes it, straight double quotes.
  One exception: the LOSATP mode name is quoted ("Kecocokan berpasangan") when
  help or message prose cites it (di mode "Kecocokan berpasangan"), because the
  bare phrase is also the ordinary noun phrase "pairwise matches". On option
  labels, headings and buttons it stays bare (Gunakan Kecocokan berpasangan).
  Other cited controls stay unquoted.
- Head noun first: TSV kedalaman, TSV BLAST, cache LOSAT, ID record, lebar
  slot track linear {id}, warna default. Never split or abbreviate a fixed
  compound (tata letak, atur ulang, ambang batas, asam amino, kandungan GC,
  nomor aksesi, kedalaman baca).
- A cited control that is a select or radio option and the subject of the
  sentence is written Opsi + label (Opsi Otomatis …); a button, a pressed
  toggle button included, is Tombol + label (Tombol Simpan sesi …).
- A mode word after a noun means the topology (record sirkular); for the
  Linear or Circular drawing, write the noun phrase + diagram linear / diagram
  sirkular (Label record diagram sirkular), or … di diagram linear in prose.
  The drawing mode is attached to the diagram also after file, reference and
  task nouns: file diagram linear, Kode genetik referensi di diagram sirkular,
  Tugas LOSATN diagram sirkular. The names of records that the record selector
  lists use the same form (Record diagram sirkular); Record sirkular is the
  topology.
- A dimension word (Lebar, Tinggi, Radius) for a depth or GC value takes track,
  because the dimension belongs to the drawn track, not to the value: Lebar
  track kedalaman, Tinggi track GC skew, Radius track kandungan GC.
- A common noun kept for a feature is lower case inside a sentence or phrase
  and capital only at the start of a label: Aturan warna untuk fitur apa pun;
  Fitur apa pun dengan {qualifier} yang cocok dengan "{pattern}".

### 9.5 Placeholders and counts

- No inflection or agreement: a placeholder sits like a noun, after its class
  noun (record {record}, file "{file}").
- The count precedes an unchanged noun ({count} record, {count} pasangan gen):
  no reduplication after a number ({count} fitur-fitur), no English plural, no
  classifier (buah); {current} dari {count}, {drawn}/{total}.
- tn: no plural forms; the `one` string repeats the `other` and must
  read well for 1. Locators: Baris {n}, File {n}; ke-{n} only for "the nth".
- Worked pairs (`one` = `other`):
  1. ({count} records) → ({count} record)
  2. From these {count} color rules: → Dari {count} aturan warna ini:
  3. Reset {count} placements to Auto → Atur ulang {count} penempatan ke
     Otomatis
  4. {current} of {count} matches → {current} dari {count} hasil
  5. Input files: {count} entries of {path} could not be read and were dropped.
     → File input: {count} entri dari {path} tidak dapat dibaca dan dibuang.

### 9.6 Length

- About 1.05–1.15× English overall, but short verbs grow
  1.4–2.25× (Terapkan 1.65×, Urungkan 2.0×) and meN- progress forms run 8–12
  letters: allow +50% for labels of up to 10 characters, +15–20% for sentences.
- Rail short names (52 px buttons, 10 px semibold) hold about nine letters on
  one line; Pengaturan (54 px) overflows. Use Input, Banding, Lanjutan,
  Daerah, Tata letak, Font, Judul, Legenda, Warna, Label, Sumbu, GC / skew,
  Fitur, Tentang; Buat beside the rail. A rail name of two words takes two
  lines when it does not fit; a single word must fit (about nine letters) or
  get a rail key. Pengaturan still does not fit as one word.
- Header buttons within about 1.5× the English: check Atur ulang pengaturan (21
  letters for 14) and Urungkan at 1,024 px. To shorten, drop yang, untuk and
  the subject; never switch terms; abbreviate only after a screenshot shows
  overflow.

### 9.7 Machine checks

The rule lists live in `tests/web/i18n-lint/id.tsv` and run in
`tests/web/i18n-catalogs.test.mjs` on every catalog value, never on the
key, after each `{name}` placeholder is masked (as `{X}`); `error`
fails the test and `warn` is printed. Each rule flags:

- Malay and vendor forms: `id-malay-download` (error) muat turun / muat naik;
  `id-malay-words` (error) the other Malay words (fail, tetapan, pelayar,
  boleh, kerana …) and -iti forms (identiti, kualiti);
  `id-official-coinage` (warn) galat, peramban, bawaan, prasetel, utas;
  `id-d05-terms` (error) the rejected forms of the provisional terms (Setelan,
  Setel ulang, Tingkat lanjut, Download, Upload); `id-d05-file` (warn) berkas.
- Spelling and register: `id-spelling` (error) sekuen, sirkuler, linier,
  unting, helai, eror, histori, trek, penyejajaran, persejajaran, pensejajaran
  (not the glossary form; some are attested variants, not misspellings);
  `id-anda-case` (error) lower-case anda; `id-informal` (error) kamu, gak,
  udah; `id-register-phrases` (warn) tidak bisa, gagal untuk, silakan;
  `id-tumpang-tindih` (error) tumpang-tindih; `id-r1-terms` (error) the terms
  retired in round 1 (komplemen balik, tag lokus, bilah skala, kueri, subjek,
  ke hulu / ke hilir).
- Word forms: `id-di-ke-preposition` (error) diatas, kebawah;
  `id-di-prefix-split` (error) di hapus; `id-men-assimilation` (warn)
  mensimpan; `id-affixed-loanword` (error) di-reset; `id-nominalized-command`
  (warn) a short value opening with Penyimpanan and the like;
  `id-progress-ellipsis` (warn) a bare meN- phrase without … or period;
  `id-title-case` (warn) Title Case; `id-english-plural` (warn) records,
  tracks; `id-count-reduplication` (warn) {count} fitur-fitur.
- Punctuation, numbers and script: `id-ellipsis` (error) ... and a space before
  …; `id-decimal-comma` (error) 0,65; `id-thousands-dot` (warn) 1.000;
  `id-unit-pb` (error) pb; `id-curly-quotes` (error) “ ” ‘ ’ « »;
  `id-foreign-letter` (error) accented letters and other scripts;
  `id-unexpected-char` (warn) other non-ASCII, NBSP included.

## 10. Vietnamese (`vi`)

Standard northern written Vietnamese: the UI chrome in the Web-app register
of Firefox, Google and Android, help texts in the register of Vietnamese
journals (English in parentheses at a term's first use). The Language menu
entry is **Tiếng Việt**. **Every value is NFC** (`normalize('NFC')`, also
for text pasted from PDFs): decomposed tone marks look the same but break
search and the lint (`vi-nfc`). The UI shows Vietnamese in Inter: the five vendored
Inter `vietnamese` subset files (the `unicode-range` blocks in
`gbdraw/web/index.html`) carry the Vietnamese letters, so a word never mixes
typefaces.

### 10.1 Vocabulary

- **New orthography** (rule QĐ 1989/QĐ-BGDĐT of 2018, 2024 Sinh học 12
  textbooks): gene, DNA, RNA, mRNA, nucleotide, dinucleotide, amino acid,
  protein, enzyme, codon, allele, plasmid; never gen, ADN, ARN, nuclêôtit,
  axit amin, prôtêin, enzim, côđon. No plural `-s`.
- **Field terms** (the glossary decides; never swap them): **hệ gene**
  (genome; as the 2018-program textbooks and papers that spell gene write it),
  trình tự (sequence; never chuỗi), bản ghi (record), **mạch** (strand: mạch
  xuôi (+), mạch ngược), bổ sung ngược, **chú thích** (annotation) and **chú
  giải** (legend), **sắp gióng cột (alignment)** for sequences ("(alignment)"
  once per string) and **căn chỉnh** for aligning records, **độ đồng nhất
  (identity)** and **độ tương đồng** (similarity), hàm lượng GC, độ sâu / độ
  sâu đọc / độ bao phủ, vạch chia, thang đo, sơ đồ dạng vòng / dạng thẳng,
  bảng mã di truyền, chủng, mã số (accession), **mô tả** (a record definition,
  the GenBank DEFINITION line; never định nghĩa), đoạn khớp (a drawn match),
  kết quả phù hợp (a search result).
- **Kept in English**, lower case inside a phrase, after the head noun, no
  plural: feature, qualifier, track, query, subject, hit, E-value, Bitscore,
  RBH, reciprocal best hit, near-reciprocal hit, orthogroup, ortholog,
  paralog, synteny, Collinear, GC skew, AT skew, skew, locus, locus tag,
  runtime, seed (loại feature, bố cục track, số hit tối đa). A help text may
  gloss one once, never a label: feature (như gene, CDS, tRNA), qualifier
  (thuộc tính), track (dải / vòng), query (trình tự truy vấn), subject (trình
  tự đích), GC skew (độ lệch GC).
- **Diagram mode after a noun**: a record, sequence, file, region, row or
  Session of the Circular or Linear diagram takes sơ đồ: bản ghi trong sơ đồ
  dạng thẳng, Nhãn bản ghi (sơ đồ dạng vòng). bản ghi / trình tự / hệ gene /
  tệp / vùng / tham chiếu + dạng vòng or dạng thẳng states the molecule
  topology (Display start requires a circular record = … bản ghi dạng vòng).
  Nouns without a topology keep the short form: track dạng vòng, khe track
  dạng thẳng {id}. The record selector of the Circular diagram is Bản ghi
  trong sơ đồ dạng vòng (the selector's aria name), and bản ghi dạng vòng is
  the topology only. A task or reference of a diagram takes the bracket form:
  Tác vụ LOSATN (sơ đồ dạng vòng), Bảng mã di truyền tham chiếu (sơ đồ dạng
  vòng).
- A stroke-width placeholder under a thickness label uses the label's word:
  Độ dày (not Chiều rộng).
- **Command words** (Firefox, Google, Android): Cài đặt, Mặc định, Đặt sẵn
  (Preset; after a head noun: Chế độ track đặt sẵn, Bộ màu đặt sẵn), Phiên,
  Tạo, Lưu, Tải (Load; Mở phiên), Nhập / Xuất, Tải lên / Tải xuống, Hoàn tác /
  Làm lại, Áp dụng, Đặt lại, **Huỷ**, Đóng, Xoá (Delete, Clear), Gỡ bỏ
  (Remove), Nhân bản, Bỏ qua, **tệp**, Thẻ, Bật / Tắt, Nhấp.
- **Loanwords**: Vietnamese where Firefox, Google or the journals have a word
  (tệp, bộ nhớ đệm); ID, web, pixel, pipeline, intron as written; never a
  phonetic respelling (pích-xen, lô-cút).
- **Never** (beyond the lint of §10.7): biểu đồ for the map, chuỗi for a
  sequence, đặc điểm for a feature, tuỳ chọn for "optional".

### 10.2 Buttons and sentences

- **Register**: bạn when a pronoun is needed, usually none; của tôi only
  where the English speaks as the user (Giữ {count} màu của tôi).
- **Buttons, menu commands, tabs, chips**: bare verb or verb + object, no
  period, no Vui lòng or hãy: Lưu, Huỷ, Áp dụng, Lưu phiên, Thử tạo lại.
- **Labels, headings, column headers, options**: noun phrases, head first:
  Cỡ chữ nhãn, Bố cục track, Độ đồng nhất tối thiểu (%).
- **Help and instructions**: full sentences with a period; an English
  imperative is a bare verb (Chọn một bản ghi.); Hãy may open an instruction
  after a reason; Vui lòng only opens a request in a blocking or error
  message, once at most.
- **Progress** Đang + verb + … (Đang lưu phiên…); busy reasons Đang tạo sơ
  đồ. Hãy thử lại sau khi tạo xong. **Completion** Đã + verb (Đã sao chép).
- **Errors, "cannot"**: active Không thể + verb (Không thể tải tệp.) or
  Không + verb + được (Không đọc được tệp.); Không tìm thấy {name}.; generic
  Đã xảy ra lỗi khi …. Never "… không được tìm thấy", "… thất bại", "được …
  bởi". A short error label has no period (Lỗi không xác định); no "!".
- **Active voice**: Đã lưu tệp., not Tệp đã được lưu.; được without an agent
  stays where it is the natural adjective (được hỗ trợ, được giữ lại).
- **Required / optional / recommended**: Bắt buộc, Cần … (Cần áp dụng),
  (không bắt buộc) after a field name, Khuyến nghị. Quote a control: “Tạo”.

### 10.3 Punctuation and spacing

- Syllables are space-separated, so Latin terms, digits, placeholders and
  units take one ordinary space (Chạy LOSAT, Trình tự {n}, 100 bp); no extra
  rule at script boundaries.
- No space before `, . ; : ! ?`, one after. A label colon follows the word:
  Nhãn:, Trường: {label}. Sentences end with a period; labels, buttons,
  headings and chips do not; a dialog title keeps the English question mark.
- Quotes “ ” for user-typed names and quoted labels; a `{name}` the
  English quotes is quoted, a bare one stays bare (STYLE §1.3). Straight
  quotes only in literal code (`/key="value"`).
- Parentheses: no inner spaces, one space before (Vùng (không bắt buộc)).
  Ranges and symbols as the source (Hàng 3–4; → ← × −); in prose an em dash
  becomes a comma or colon.
- One `…`. Numbers keep the English format: 1,000 and 0.65, never
  1.000 or 0,65. `&` in a title becomes và (Trục và thang đo).

### 10.4 Words and capitals

- **Sentence case** for labels, buttons, headings, menu items, tabs and
  options: Cài đặt nâng cao, Tải xuống TSV. No Title Case (Cài Đặt Nâng Cao),
  no ALL CAPS (stacked tone marks collide). The English capital headings
  become sentence case: Bộ màu đặt sẵn (PRESET SCHEMES), Màu mặc định, Quy tắc
  riêng, Lọc nhãn, Thứ tự ưu tiên qualifier, Vị trí nhãn (LABEL GEOMETRY),
  Thêm khoá (ADD KEY). `:lang(vi)` turns off `uppercase` and `tracking-*`.
- Names keep their capitals (GenBank, LOSAT); gbdraw stays lower case. A
  common noun kept in English is capitalized only at the start of a label
  (Feature, Track) and is lower case inside a sentence or phrase: Quy tắc màu
  cho feature bất kỳ.
- **Tone marks in the new position** in every string (QĐ 1989/QĐ-BGDĐT Điều 8,
  2018 textbooks, Google; many papers still use the old position): hoà, hoá,
  xoá, khoá, toạ, khoẻ, tuỳ, huỷ, thuỷ, so Xoá, Huỷ, Tuỳ chỉnh, Tuỳ chọn, Toạ
  độ, Chuẩn hoá, Mã hoá. Only open syllables in oa, oe, uy change; closed ones
  are the same in both systems (hoàn, khoảng, thuyết). Cancel is Huỷ: not Hủy
  bỏ, and not Hủy (old position).
- **i / y** as Google, Firefox and most Vietnamese journal papers: tỷ lệ, lý
  do, kỹ thuật, ký tự, ký hiệu. This deliberately differs from the textbook
  rule QĐ 1989/QĐ-BGDĐT Điều 9, which writes i (tỉ lệ, lí luận).

### 10.5 Placeholders and counts

- Count before the noun, no plural marker, and these nouns need no
  classifier: {count} bản ghi, {count} feature, {count} luồng. các / những
  only without a number; "all {count}" is tất cả {count} ….
- No plural forms: the `one` row repeats the `other` text.
- "Records: {n}" becomes {n} bản ghi or Số bản ghi: {n}; locators put the
  number last (Trình tự {n}, Dòng {n}); an ordinal is bản ghi thứ {n}.
- Placeholders stay verbatim at the end of their noun phrase; identifiers take
  plain spaces, no quotes: Dùng large_tick_interval thay cho tick_interval.
- Worked examples (both rows get the same text):
  1. `{count} record` / `{count} records` → `{count} bản ghi`
  2. `{count} gene pair` / `… pairs` → `{count} cặp gene`
  3. `{current} of {count} match` / `… matches` → `{current} trên {count} kết quả phù hợp`
  4. `{count} feature updating...` / `… features …` → `Đang cập nhật {count} feature…`
  5. `Reset {count} placement to Auto` / `… placements …` → `Đặt lại {count} bố trí về Tự động`

### 10.6 Length

- Letters are Latin-width, so character counts are a fair guide: the core
  words measured 0.98× the English width in Inter 10 px; a word rendered as a
  phrase reaches 1.5–1.7× (Áp dụng 1.47×; Cài đặt sẵn, the earlier Preset
  form, 1.74×).
- Rail short names (52 px buttons at 9–10 px), about 9 characters: Đầu vào,
  Feature, Nhãn, Màu, Trục, Vùng, So sánh, Nâng cao, Tiêu đề, Bố cục, GC /
  skew. Phông chữ and Giới thiệu are the widest: Phông chữ takes two lines at
  1,440 px; both fit. The short form Phông is not used. The rail names have
  line-height 1.3 in vi so a two-line name keeps its diacritics clear.
- Header and bar buttons: within about 1.3× the English (Lưu phiên, Tải ví
  dụ, Đặt lại cài đặt); Huỷ, not Huỷ bỏ. No abbreviations with periods
  (Trình đ.khiển). Help texts wrap. Headings need line-height ≥ 1.4 and no
  negative letter-spacing, or Ấ, Ổ, Ữ clip.

### 10.7 Machine checks

The rule lists live in `tests/web/i18n-lint/vi.tsv` and run in
`tests/web/i18n-catalogs.test.mjs` on every catalog value, never on the
key, after each `{name}` placeholder is masked (as `{X}`); `error`
fails the test and `warn` is printed. Each rule flags:

- `vi-old-gen` (error): the word gen anywhere, also in hệ gen (genome is hệ gene).
- `vi-old-adn-arn` (error): ADN, ARN, mARN, tARN or rARN instead of DNA / RNA.
- `vi-old-acid` (error): axit amin, acid amin or axit nuclêic.
- `vi-old-loanword` (error): nuclêôtit, prôtêin, côđon, enzim, pôlipeptit, ribôxôm, plasmit, alen.
- `vi-genome-term` (error): bộ gen or bộ gene instead of hệ gene.
- `vi-tap-tin` (error): tập tin instead of tệp.
- `vi-huy-bo` (error): Hủy bỏ or Huỷ bỏ instead of Huỷ.
- `vi-tai-ve` (error): Tải về instead of Tải xuống.
- `vi-thiet-dat` (error): Thiết đặt instead of Cài đặt.
- `vi-that-bai` (error): thất bại instead of Không thể … or không thành công.
- `vi-passive-not-found` (error): … không được tìm thấy instead of Không tìm thấy ….
- `vi-passive-by` (warn): được … bởi with up to six words between (bởi vì excluded).
- `vi-tone-old-position` (error): a tone mark on the o of a final oa / oe or the u of a final uy (Hòa, Xóa, Tùy, Hủy, Tọa).
- `vi-y-spelling` (warn): tỉ lệ, tỉ số, lí do, kĩ thuật, kí tự, kí hiệu.
- `vi-nfc` (error): a combining mark U+0300–U+036F, so the value is not NFC.
- `vi-ascii-word` (error): one of 28 common UI words without diacritics (khong, cai dat, luu, tep, …).
- `vi-all-caps` (error): an upper-case word of two or more letters with a Vietnamese capital (CÀI ĐẶT, TẮT); ASCII acronyms such as LOSAT pass.
- `vi-title-case` (warn): two capitalized words in a row, the second with a Vietnamese letter (Tải Xuống); Tiếng Việt is the known exception.
- `vi-register` (error): xin vui lòng, quý vị, quý khách or chúng tôi.
- `vi-vui-long-label` (warn): a value that starts with Vui lòng and has no sentence punctuation.
- `vi-straight-quote` (error): a straight " that opens a quoted word or placeholder.
- `vi-space-before-punct` (warn): a space between a word and `, . ; : ! ?`.
- `vi-decimal-comma` (warn): a decimal comma (0,65) where the English format keeps 0.65.
- `vi-english-plural` (warn): features, tracks, hits, queries, subjects, qualifiers, orthogroups, genes, reads or loci.
- `vi-ellipsis` (error): three periods instead of one ….
- `vi-preset-term` (error): cài đặt sẵn instead of đặt sẵn after a head noun (chế độ track đặt sẵn, bộ màu đặt sẵn).

## 11. Thai (`th`)

Neutral written Thai as Google software writes it (Chromium, AOSP, Flutter)
for command words, IPST / NSTDA usage for biology, and English for the BLAST
and comparative-genomics words Thai papers leave in English. Royal Institute
coinages appear only as alternatives in glossary notes. The Language menu shows
**ไทย**. No Thai font file is bundled: Thai uses the OS font stack (a bare
Linux without a Thai font shows missing-glyph boxes) at 12 px or more, with
`line-height` 1.5 or more.

### 11.1 Vocabulary

- **Field terms**: จีโนม, ยีน, ลำดับ, นิวคลีโอไทด์, กรดอะมิโน, สาย (strand),
  สายพันธุ์ (strain), สิ่งมีชีวิต (organism), ออร์แกเนลล์, คำอธิบายประกอบ
  (annotation), ฟีเจอร์, ระเบียน (record), บริเวณ (region), ปริมาณ GC,
  การเทียบลำดับ (sequence alignment), การจัดแนว (Similarity-group alignment),
  ความคล้ายคลึง (similarity), ความเหมือน (sequence identity),
  ความอนุรักษ์ของลำดับ (conservation), ผลิตภัณฑ์ของยีน (/product),
  การเปรียบเทียบรายคู่ (pairwise comparison), แผนภาพ, แทร็ก. The stroke of a
  line element is the line itself (สีเส้นแกน, ความหนาเส้น); เส้นขอบ is the
  outline of a block or swatch. A dimension word for GC content or GC skew
  takes แทร็ก, as it already does for depth: ความกว้างแทร็กปริมาณ GC. A
  stroke-width placeholder under a thickness label uses the label's word:
  ความหนา (not ความกว้าง).
- **Kept in English**, no plural s: E-value, Bitscore, Query, Subject, hit,
  RBH, orthogroup, ortholog, paralog, coortholog, inparalog, upstream,
  downstream, Collinear, GC skew, AT skew, coverage, qualifier, accession
  (หมายเลข accession), locus tag, LOCUS, the §1.1 names, units. Help may gloss
  one once: ค่า E; ความเหมือน may take “(identity)” once.
- **Command words** (Google): บันทึก, โหลด, ยกเลิก, เพิ่ม, ลบ (Delete), นำออก
  (Remove), ล้าง (Clear), ทำสำเนา, แก้ไข, นำไปใช้, เลิกทำ / ทำซ้ำ, เปลี่ยนกลับ,
  กู้คืน, รีเซ็ต, นำเข้า / ส่งออก, อัปโหลด, สร้าง (Generate), เรียกใช้ (Run).
- **UI nouns** (Google): การตั้งค่า, ค่าเริ่มต้น, ค่าที่ตั้งไว้ล่วงหน้า /
  พรีเซ็ต (not ค่าที่กำหนดล่วงหน้า: กำหนดเอง is Custom), เซสชัน, ไฟล์,
  ตัวอย่าง (Example), ภาพตัวอย่าง (Preview), แกลเลอรี, ขั้นสูง, เลย์เอาต์,
  แบบอักษร, ชุดสี, ป้ายกำกับ, ประกาศ (Notice; หมายเหตุ is /note; a notice that
  refers to itself says ข้อความนี้), การลบล้าง (override, as a noun or label;
  a help sentence says ใช้ … แทน).
- **Legend**: คำอธิบายสัญลักษณ์ in sentences, titles and panel labels; the
  short form คำอธิบาย only on the rail, the Legend tab and buttons, so it
  never stands next to คำอธิบายประกอบ (annotation) or คำจำกัดความ
  (Definition).
- **Loanwords**: only those the glossary lists, in the Google spelling
  (§11.4); otherwise a native word. Prefer แบบอักษร to ฟอนต์.
- **Never** (official coinages and known mistranslations): ค่าโดยปริยาย, บรรณนิทัศน์, เอกลักษณ์,
  การปรับแนว, จีน for gene, อินทรีย์, ภูมิภาค, ข้อคำถาม, บรรจุลง / บรรจุขึ้น,
  แฟ้ม, สายโยงใย, ช่วงเวลาสื่อสาร, ทำกลับ, จุดตั้งระยะ, คะแนนบิต, ความเบ้.
- **Retired** (r1 term audit, lint `th-retired-term`): แอนโนเทชัน, ความคล้าย
  without คลึง, การเปรียบเทียบเป็นคู่, การจับคู่แบบเป็นคู่, ออร์โธล็อก /
  พาราล็อก and their compounds, การอนุรักษ์ for sequence conservation, ต้นน้ำ
  / ปลายน้ำ.

### 11.2 Buttons and sentences

- **Register**: neutral, as in Chromium; no particles (ครับ, ค่ะ, นะ), no ท่าน
  or เรา. Leave the subject out, or write คุณ (user) or ระบบ (app) if needed.
- **Buttons, menu commands, toggles**: a bare verb or verb + object (บันทึก,
  บันทึกเซสชัน, นำไฟล์ออก); a command that opens a further step keeps `…`.
- **Headings, tabs, column headers**: nouns; an action becomes a การ- noun
  (การส่งออก, การแก้ไขเลย์เอาต์). The key's kind decides ส่งออก or การส่งออก.
- **Instructions and help**: complete sentences with no final period; an
  imperative is a bare verb (เลือกระเบียนที่จะครอบตัด). "Please" is โปรด, at
  most once per message, never กรุณา.
- **Status "…ing"**: กำลัง + verb + `…`: กำลังโหลด…, กำลังสร้างแผนภาพ…
- **Completion**: …แล้ว (บันทึกแล้ว); "not yet" is ยังไม่ได้… (ยังไม่ได้นำไปใช้).
- **Cannot, failed, errors**: verb phrase + ไม่ได้ (อ่านไฟล์ไม่ได้); X +
  ไม่สำเร็จ for a failed operation (ส่งออก PNG ไม่สำเร็จ); เกิดข้อผิดพลาด,
  ไม่ถูกต้อง (invalid), ไม่พร้อมใช้งาน (unavailable), ไม่รองรับ (unsupported).
- **Required, optional, recommended**: the badge จำเป็น, or ต้องระบุ X in a
  sentence; the tags (ไม่บังคับ), (แนะนำ), (ค่าเริ่มต้น).
- **Questions** end with หรือไม่ and no `?`; no `!` anywhere.
- **Active voice**: avoid the adverse passive ถูก…; leave the agent out or
  name ระบบ (ทิ้งบรรทัด {lines} แล้ว, not บรรทัด {lines} ถูกทิ้ง).
- A control named in a sentence is its Thai label in “ ”: คลิก “สร้าง”.
- **Relative clauses**: never stack ที่…ที่; the second clause takes ซึ่ง, or
  a redundant modifier is dropped (ฟีเจอร์ที่เลือก). At most two การ- nouns in
  a row.
- **Generate again**: สร้างแผนภาพอีกครั้ง, without quotes; คลิก “สร้าง” only
  when the button is clicked.
- **Any feature**: a whole-sentence phrase writes ฟีเจอร์ใดก็ได้ with no space
  before its ที่ clause (ฟีเจอร์ใดก็ได้ที่ {qualifier} ตรงกับ “{pattern}”);
  กฎสีของฟีเจอร์ใดก็ได้ is one phrase, not two fragments.
- **Drawing mode**: Linear / Circular after a record, file, sequence, axis or
  scale attaches to the diagram: ระเบียนในแผนภาพแบบเส้นตรง,
  ลำดับที่แสดงในแผนภาพแบบวงกลม. ระเบียนแบบวงกลม is the topology of the
  molecule only; tracks and slots keep แทร็กแบบวงกลม. The record selector of
  the Circular diagram is ระเบียนในแผนภาพแบบวงกลม (the selector's aria name).

### 11.3 Punctuation and spacing

- **Sentence ends**: no period, comma or semicolon after Thai text. One space
  separates clauses and sentences; keys shown together join with `' '`. No
  space inside a clause: never between a noun and its ที่-clause, a subject
  and its verb, or a verb and its object; this rule wins over th-long-run, so
  a sentence whose closed clause would run 60 or more Thai characters is
  restructured (a fronted เมื่อ / หาก / จาก phrase, ระบบ as subject, a break
  before เพื่อ(ให้), โดย, ใน…, และ or หรือ, or a shorter synonym) instead of
  keeping the space. Periods stay inside Latin names, numbers and ranges
  (0.65, 1..200).
- **Lists**: Thai items are separated by spaces, with และ / หรือ before the last
  item (ฟีเจอร์ ป้ายกำกับ และคำอธิบายสัญลักษณ์). Latin-only lists (formats,
  identifiers, code values) keep commas: A, B หรือ C. Run-time
  `Intl.ListFormat('th')` output stays as it is ("A B และC", no comma).
- **Latin, digits, placeholders, units**: one space between a Thai letter and
  a Latin letter, digit or `{placeholder}`, on either side: ปริมาณ GC,
  {count} ระเบียน, ไฟล์ {name} ไม่ถูกต้อง. Units keep the source spacing
  (100 bp, 20px); % attaches to its number. No NBSP or ZWSP.
- **Colon** in a label: attached, one space after: สถานะ: {state}.
- **Parentheses**: one space outside, none inside: บริเวณ (ไม่บังคับ).
- **Quotes**: “ ”, nested ‘ ’; never ASCII quotes.
- **Ellipsis**: one `…`, attached. **ๆ** attaches to its word and takes one
  space after it (ต่างๆ ที่ใช้).
- **Ranges and symbols** keep the source form (3–4, GC / skew, → × −).
- **Digits**: ASCII only, never ๐–๙; English decimal point and grouping.

### 11.4 Spelling of transliterations

- Thai has no case: Title Case and ALL CAPS headings become ordinary text
  (PRESET SCHEMES ชุดสีที่ตั้งไว้ล่วงหน้า); `:lang(th)` turns off `uppercase`
  and `tracking-*`. Latin words keep their own case.
- **No tone marks** in transliterations, as Google spells them: เซสชัน,
  เวอร์ชัน, อัปโหลด, อัปเดต, แกลเลอรี, เลย์เอาต์, ไฮไลต์, รีเซ็ต, คลิก, ลิงก์,
  ไอคอน, อีเมล, อินเทอร์เฟซ. Never เซสชั่น, เวอร์ชั่น, อัพโหลด, อัพเดท,
  แกลเลอรี่, ไฮไลท์, ลิ้งค์. The one exception is Chromium's ป๊อปอัป.
- Mai taikhu ็ is not a tone mark and stays: บล็อก, แท็บ, ล็อก. Biology
  follows IPST / NSTDA: จีโนม, ยีน, นิวคลีโอไทด์, โครโมโซม, ออร์แกเนลล์.
- **Storage order**: consonant, upper or lower vowel, tone mark, then ำ
  (น้ำ). แ and ำ are one character each, never เเ or ํา.
- **Compounds** put the head noun first and reuse glossary words:
  ป้ายกำกับระเบียน (Record label), พรีเซ็ตแทร็ก (Track Preset). One English
  word keeps one Thai word (Remove นำออก, Delete ลบ, Clear ล้าง).

### 11.5 Placeholders and counts

- `{name}` is copied verbatim, placed in Thai order, with one space on each
  side except next to `( ) : “ ” /` or the edge of the value.
- **No plural forms**: `Intl.PluralRules('th')` selects `other`, and the
  `one` row repeats it.
- **Counts**: number + classifier. In labels a noun that is its own classifier
  follows the number ({count} ระเบียน, {count} ไฟล์); in sentences the noun
  comes first (ฟีเจอร์ {count} รายการ, สี {count} สี). รายการ counts features,
  entries, matches, rules and issues; status chips that count records put the
  status first: จัดแนวแล้ว {count} ระเบียน. An ordinal is ที่ + number.
- Examples:
  1. `{count} record` → {count} ระเบียน
  2. `{count} feature` → ฟีเจอร์ {count} รายการ
  3. `{count} gene pair` → ยีน {count} คู่
  4. `{current} of {count} match` → รายการที่ {current} จาก {count}
  5. `Keep your {count} changed color from the "{palette}" palette ({keys})?`
     → เก็บสีที่เปลี่ยนแปลง {count} สีจากชุดสี “{palette}” ({keys}) ไว้หรือไม่

### 11.6 Length

- Thai rarely widens a label: in Chromium EN→TH the median rendered width is
  0.85–0.92× the English. The risk is a short label whose English word
  becomes a phrase: Legend → คำอธิบายสัญลักษณ์, Preset → ค่าที่ตั้งไว้ล่วงหน้า
  (about 2.4×).
- **Code points overstate width**: about a fifth of Thai code points are
  zero-width marks (ค่าที่ตั้งไว้ล่วงหน้า: 21 code points, 13 base glyphs).
  Measure, or count base glyphs; never apply CJK character limits. Thai
  needs about 150% of the Latin line height.
- **Rail** (52 px buttons): rail names are 12 px with line-height 1.5; two
  lines are the rule of every locale (§1.5). แบบอักษร, คำอธิบาย, ป้ายกำกับ and GC / skew take two lines
  and fit the entry at 1,024 and 1,280 px (Chromium on Linux, Loma), so ฟอนต์
  and ป้าย are not used. Compare is เทียบ, not เปรียบเทียบ. The inline button
  รีเซ็ตตำแหน่งเริ่มแสดงผล (Reset start) keeps the glossary form.
- **Header buttons** come out narrower than English: สร้างแผนภาพ is 72 px
  against 105 px for Generate Diagram.
- **Short forms** (the glossary's "long / short" cells) go on buttons, the
  rail and option rows: คำอธิบาย, พรีเซ็ต. Titles and help use the long form.

### 11.7 Machine checks

The rule lists live in `tests/web/i18n-lint/th.tsv` and run in
`tests/web/i18n-catalogs.test.mjs` on every catalog value, never on the
key, after each `{name}` placeholder is masked (as `{X}`); `error`
fails the test and `warn` is printed. Each rule flags:

- `th-final-period` (error): a value ends with a period right after a Thai letter.
- `th-period` (warn): a period ends a Thai sentence inside a value, or ends a Thai value after a Latin word.
- `th-thai-digit` (error): a Thai digit ๐–๙.
- `th-comma-thai` (error): a comma between two Thai words.
- `th-comma-after-thai` (warn): any other comma right after a Thai letter.
- `th-space-thai-latin` (error): a Thai letter directly before a Latin letter, a digit or `{`.
- `th-space-latin-thai` (error): a Latin letter, a digit, `%` or `}` directly before a Thai letter.
- `th-translit-tone` (error): a tone-marked transliteration (…ชั่น, แกลเลอรี่, ลิ้งค์, คลิ้ก, ไอค่อน).
- `th-translit-spelling` (error): a transliteration spelled against Google (อัพ…, อัปเดท, เลย์เอาท์, ไฮไลท์).
- `th-coinage` (error): a banned coinage or known mistranslation (ค่าโดยปริยาย, บรรณนิทัศน์, เอกลักษณ์, แฟ้ม …).
- `th-gene-china` (error): จีน, which means China, written for gene.
- `th-kept-english` (warn): a Thai rendering of a kept-English term, or ฟอนต์ for แบบอักษร.
- `th-polite` (error): ครับ, ค่ะ or กรุณา.
- `th-pronoun` (warn): ท่าน, ข้าพเจ้า, กระผม, ดิฉัน or พวกเรา.
- `th-tone-before-vowel` (error): a tone mark U+0E48–0E4B stored before U+0E31, U+0E34–0E37 or U+0E47.
- `th-mark-order` (error): any other misordered, doubled or orphaned mark, เเ for แ, or ํา for ำ.
- `th-maiyamok` (warn): a space before ๆ, or a Thai letter right after it.
- `th-question-exclaim` (warn): `?` or `!` after Thai text.
- `th-punct-space` (warn): no space around `( )` between Thai words, a space before a colon or `…`, or none after a colon.
- `th-foreign-script` (error): Arabic, Devanagari, Lao, Myanmar, Khmer, CJK, kana, Hangul or full-width letters.
- `th-ellipsis` (error): three periods instead of `…`.
- `th-quotes` (warn): an ASCII quote instead of “ ” or ‘ ’.
- `th-long-run` (warn): 60 or more Thai characters without a space.
- `th-whitespace` (warn): a leading, trailing or doubled space, or NBSP, ZWSP, ZWNJ, ZWJ, WJ or BOM.
- `th-retired-term` (error): a term the r1 term audit retired (แอนโนเทชัน, ความคล้าย without คลึง, …เป็นคู่, ออร์โธล็อก / พาราล็อก, การอนุรักษ์, ต้นน้ำ / ปลายน้ำ).
- `th-identity-thai` (warn): the English word identity outside a “(identity)” gloss or a code value.
- `th-line-stroke` (warn): เส้นขอบของแกน / เส้นขอบของเส้น, or ไม่ยอมรับ.

## 12. German (`de`)

Standard German of Germany, formal Sie, spelling of the 1996 reform with ß (Duden). No Swiss or
Austrian variants. The Language menu entry is the endonym **Deutsch** (page `lang` de); the menu
label itself is Sprache. This section is a draft; no text states that it is reviewed.

### 12.1 Vocabulary

- **Field terms**: Genom, Sequenz, Datensatz (Record; Microsoft Datensatz), Feature (see below),
  Feature-Typ, Qualifier (the keys /gene, /product stay), Annotation, Zugangsnummer (Accession),
  Locus (plural Loci), Locus-Tag, Strang (Strand: Vorwärtsstrang; de.wikipedia GC-Gehalt writes
  DNA-Strang), reverses Komplement, GC-Gehalt, GC-Skew, AT-Skew, Skew, Dinukleotid, Tiefe (Depth),
  Lesetiefe (read depth in sentences), Abdeckung (Coverage), Identität, Homologie, Ortholog (plural
  Orthologe), Paralog (Paraloge), Treffer (BLAST hit; de.wikipedia BLAST-Algorithmus), E-Wert,
  Alignment (plural Alignments), Lücke (alignment gap), Fehlpaarung (mismatch), genetischer Code,
  Organell, Organismus, Stamm, Art (species), Produkt, Gen, Protein.
- **Kept in English** (no established German term, or the field keeps the English noun; D-12):
  Feature (GenBank's own table is FEATURES; Merkmal means a phenotypic trait), Track (UCSC and IGV
  usage), Slot, Qualifier, Query, Subject, Skew, Locus-Tag, Orthogroup, Reciprocal Best Hit,
  Near-reciprocal Hit, Coortholog, Inparalog, Ribbon, Flatfile, Popup, Seed, Pipeline, Workflow,
  Cache, Preprint, Runtime, Task, Job, Markup, Lookbehind, Issues (GitHub page name), and the rail
  name Labels. Nouns are capitalized and take German compounds with a hyphen (Feature-Typ).
- **Diagram terms**: Diagramm (the generated figure; never Karte), Track (a ring or band),
  Spur (Lane; a gel lane is a Spur, so it never means Track), Track-Slot, Platz (File card
  slot), Teilstrich (Tick; Hauptteilstrich, Nebenteilstrich), Markierung (annotation mark),
  Lineal (Ruler), Skala (Scale), Maßstabsbalken (Length bar), Innenabstand (Padding),
  Zwischenraum (Gap), Abstand (Spacing), Versatz (Offset), Kontur (Stroke of shapes), Linie
  (lines, the axis line), Füllung, Farbfeld (Swatch), Vordergrund / Hintergrund, Voreinstellung
  (Preset), Standard (Default), Beschriftung (Label), Legende, Legendeneintrag, Legendentext.
- **Two-sense words**: Alignment is the sequence alignment, Ausrichtung the layout action of
  Similarity groups (Ausrichtung prüfen), Orientierung the record orientation; Treffer is the BLAST
  hit and the search match, Übereinstimmung the drawn pairwise match; Lücke is the alignment gap,
  Zwischenraum the layout gap; Lage is the INSDC location, Position the legend position; Kontur
  is the outline, Linie the line itself (Achsenlinie, Linienfarbe); Wählen is Choose, Auswählen is
  Select; Zuklappen / Aufklappen are Collapse / Expand (Erweitern would echo Erweitert).
- **Command words** (Microsoft German): Speichern, Laden, Öffnen, Schließen, Kopieren, Einfügen,
  Duplizieren, Löschen (Delete), Entfernen (Remove), Leeren (Clear), Zurücksetzen (Reset),
  Wiederherstellen (Restore), Rückgängig (Undo), Wiederholen (Redo), Erneut versuchen (Retry),
  Anwenden, Abbrechen, Hinzufügen, Umbenennen, Durchsuchen, Hochladen, Herunterladen, Importieren,
  Exportieren, Überspringen, Beibehalten, Zusammenführen, Andocken, Drehen, Ausblenden, Anzeigen,
  Aktivieren / Deaktivieren; Ein / Aus (On / Off), Ohne (None).
- **Generate is Erzeugen** (Diagramm erzeugen, beim Erzeugen, Neu erzeugen); never generieren. It
  is 8 letters beside the rail where Generieren has 10.
- **Never use**: Merkmal or Funktion for Feature, Spur for Track, generieren, sichern for Save,
  Satz for Record, Karte for the figure, kollinear with two l
  (write kolinear, as de.wikipedia Syntenie writes Kolinearität), Deppenleerzeichen (Locus Tag),
  du forms, Bitte in instructions, English leftovers (Record, Legend, Session, Settings).

### 12.2 Buttons and sentences

- Register: Sie, written with a capital (Sie, Ihr, Ihnen). Never du, dein, euch. Microsoft and
  Apple German use Sie for software text; no source was fetched this session (from memory).
- Buttons, menu commands, links, checkboxes: the infinitive, no article where the English has
  none: Sitzung speichern, Beispiel laden, Datei hochladen, Tiefe anzeigen. A reflexive object
  stays short: Zurücksetzen, Neu erzeugen.
- Labels, headings, tabs, column headers, badges, options: noun phrases, no final period
  (Benutzerdefinierte Track-Slots, Paarweise Vergleiche). Options are nouns or adjectives that
  agree with the setting (Automatisch, Seriell, Multithread). A dialog title that asks is a
  question with the infinitive: Datei vollständig entfernen?
- Instructions (help, notices, dialogs, errors): the Sie imperative (Wählen Sie einen
  Datensatz. Klicken Sie auf Erzeugen.). Leave blank for auto: Zur automatischen Bestimmung leer
  lassen. Separable verbs put the particle last: Wählen Sie die Datei aus.
- Progress ("…ing"): passive with werden and … (Sitzung wird gespeichert…, Diagramm wird
  erzeugt…). Completion: the noun and its participle (Sitzung gespeichert.), or a passive sentence
  (Die Ausrichtung wurde angewendet.). Errors: what failed, then what to do (Die Datei konnte nicht
  gelesen werden. Prüfen Sie, ob es sich um eine FASTA-Datei handelt.).
- "Cannot": a state is nicht möglich / kann nicht (Rückgängig ist derzeit nicht möglich.); a
  failure is konnte nicht (… konnte nicht gespeichert werden.). Retry after it finishes: Versuchen
  Sie es erneut, wenn der Vorgang abgeschlossen ist.
- Required / optional / recommended: erforderlich (invariable after a field: Region (optional)),
  optional, empfohlen. Not Pflichtfeld, not optionell.
- An option or button name is never the grammatical subject: Bei „Auto“ wird …; Beim
  Speichern wird …, not Auto wählt … . A text cites the label the control shows: the options
  read Auto, so help and notices cite „Auto“; „Automatisch“ only where the control itself shows
  automatisch.

### 12.3 Punctuation and spacing

- Quotes „…“ (U+201E, U+201C), nested ‚…‘ (U+201A, U+2018). A {name} the English quotes is
  „{name}“; a UI word the English quotes ('Bar') is „Balken“. Never straight " or ', never “…”,
  ”, « » in German prose. A period or comma follows the closing quote.
- Colon: no space before, one after (Feld: {path}). A label that ends with a colon keeps it.
  No space before ? ! ; either.
- Ellipsis: one … (U+2026) directly after the word (Laden…, Sitzung wird geladen…).
- Dashes: the en dash with spaces for a break (Wert – Maximum), the en dash without spaces for
  ranges and between numbers (Zeilen 3–4), the hyphen only inside compounds. Slash as in the source.
- Digits are ASCII with the English decimal point and grouping (0.65, 1,000.5), never 0,65 or
  1.000 or 1 000 (D-05); units and % keep the source spacing (100 bp, 50%).
- Abbreviations: z. B., d. h., u. a. with a no-break space (U+00A0); min., max.,
  Abb. take a period; ID is invariable (die IDs). Sentences are joined with one space (D-09).
- Parentheses as in English, no inner space: Region (optional).

### 12.4 Words and capitals

- Every noun is capitalized; labels are sentence case, never English Title Case: Benutzerdefinierte
  Track-Slots, Paarweise Vergleiche, Achse und Skala.
  Adjectives, verbs and the words of an option written in lower case in the source stay lower case
  (hoch, mittel, niedrig; the Python-owned member values).
- A key the English writes in capitals stays in capitals; Latin script keeps the CSS `uppercase`
  rule, which writes ß as SS and keeps Ä, Ö, Ü (STANDARDFARBEN, SPEZIFISCHE REGELN,
  BESCHRIFTUNGSFILTER, QUALIFIER-PRIORITÄT, SCHLÜSSEL HINZUFÜGEN, AUS).
- The action Generate in running text is the verb or noun in German case rules (beim Erzeugen,
  vor dem Erzeugen, nach dem Erzeugen). Only the button itself is cited with its label Erzeugen
  and without quotes (Klicken Sie auf Erzeugen); where the English quotes it, use „Erzeugen“.
- Concepts the English capitalizes (Session, Result, Legend) follow normal German noun capitals
  (die Sitzung, das Ergebnis, die Legende); gbdraw stays lower case.
- Compounds put the head last: Datensatzbeschriftung, Legendeneintrag, Tiefen-Track; the
  Fugen-s follows Duden (Sitzungseinstellungen, Regionsannotationen, Zugangsnummer). Hyphenate
  a compound with an acronym, a number or a kept English noun (Feature-Typ, GFF3-Datei, Track-Slot);
  never write the parts as two words.
- Genders: der Datensatz, das Feature, der Track, der Slot, die Spur, die Legende, der Strang,
  die Sitzung, die CDS, die Region, das Set, der Anker, die Beschriftung, das Alignment, der
  Treffer, die Zugangsnummer, die Einstellung, das Ergebnis, der Qualifier, der Locus, die Query.
- A count of members is "Anzahl der …" (Anzahl der Gene in der Gruppe); exact match is parenthesized
  after the field: Produkt (genaue Übereinstimmung): {product}.

### 12.5 Placeholders and counts

- A placeholder's gender, number and case are unknown and it is never inflected, never carries
  a suffix, and never takes an article that depends on its value. Put a German noun before it
  and decline the noun: die Datei {file}, der Datensatz {record}, des Datensatzes {record}, dem
  Track {track}, die Option {setting}, die Palette „{palette}“. Not des {record}, der {track}.
- A placeholder that holds a translated phrase ({action}, {reason}) follows a colon:
  Rückgängig: {action}. Lists ({pairs}, {names}) follow a number-bearing noun that agrees:
  Paar {pairs} ist … / Die Paare {pairs} sind ….
- {count} stands before the noun and is not declined; the noun and the words agreeing with it
  follow. `Intl.PluralRules('de')` returns `one` for the integer 1 alone; 0, 1.5, 21, 22, 25 and
  every other count take `other`: 0 Datensätze, 1.5 Datensätze, 21 Datensätze. A `one` string
  drops {count} only where the English `one` does.
- Worked pairs, `one` / `other`:
  1. ({count} records): ({count} Datensatz) / ({count} Datensätze)
  2. {count} features ready: {count} Feature bereit / {count} Features bereit
  3. From this color rule: / From these {count} color rules: Aus dieser Farbregel: / Aus diesen
     {count} Farbregeln:
  4. Reset {count} placements to Auto: {count} Platzierung auf Automatisch zurücksetzen /
     {count} Platzierungen auf Automatisch zurücksetzen
  5. {count} chosen feature keeps its current direction.: {count} gewähltes Feature behält seine
     aktuelle Richtung. / {count} gewählte Features behalten ihre aktuelle Richtung.

### 12.6 Length

- German runs 20–30% longer than English (compounds, longer verbs); a button with an object about
  1.5 times. Take the short form of the glossary first; abbreviate only with an accepted form.
- Rail short names (52 px entries, 10 px semibold, about nine Latin letters; a word is never
  split, there are no soft hyphens): Eingabe (Input), Vergleich (Compare), Erweitert (Advanced,
  Linear only), Regionen, Layout, Schrift, Titel, Legende, Farben, Features, Labels (own rail key;
  Beschriftungen has 14 letters), Achse, GC / Skew (two words, wraps between them), Info. The
  Generate button beside the rail reads Erzeugen. The section titles keep the plain words
  (Beschriftungen, Eingabe-Genome).
- Header buttons: prefer the short forms Zurücksetzen (Reset Settings), Sitzung speichern,
  Sitzung laden, Beispiel laden. Einstellungen zurücksetzen is the long form and is used only if
  the 1,024 px fit check passes (D-11). Options in 10 px rows use Auto, Seriell, Multithread,
  Oben / Unten, Innen / Außen, Auf der Achse, Hintergrund, Ribbon / Kurve.
- To shorten a label drop the article first (Track-Breite, not Breite des Tracks); never switch terms.

### 12.7 Machine checks

The rules are in `tests/web/i18n-lint/de.tsv` (draft: `lint-de.tsv`); each runs on one catalog value.

- `de-sie-lowercase` (error): lower-case sie after a sentence-initial imperative infinitive.
- `de-du-pronoun`, `de-du-imperative` (error): du / dein / euch forms and du-imperatives.
- `de-bitte` (warn): Bitte at the start of an instruction.
- `de-spelling` (error): Swiss or transliterated spelling (ss for ß, ae / oe / ue) and old spellings.
- `de-quote-straight`, `de-quote-closing`, `de-quote-other` (error / error / warn): wrong quotes.
- `de-space-before-punct` (error): space before : ; ? !. `de-hyphen-dash` (warn): " - " as a dash.
- `de-three-dots` (error): "..." or a space before the final ellipsis.
- `de-decimal-comma` (error), `de-thousand-space` (warn), `de-nonascii-digits` (error): numbers.
- `de-placeholder-glue` (error): a letter glued to a placeholder. `de-apostrophe-genitive` (error).
- `de-compound-space` (warn): GC Gehalt, Locus Tag written as two words.
- `de-feature-merkmal` (warn), `de-generieren` (error), `de-kollinear` (warn): glossary terms.
- `de-english-leftover` (warn): Record, Legend, Settings and similar English words.
- `de-abbrev-spacing` (warn): z. B. without the no-break space. `de-soft-hyphen` (error).

## 13. Italian (`it`)

Standard Italian, informal second person (tu) as in Microsoft, Apple and Google
software (a draft). The Language menu entry is the endonym **Italiano**; the
menu label is Lingua; the page `lang` is `it`. Sources: `glossary-it.notes.md`.

### 13.1 Vocabulary

- **Field terms**: genoma, sequenza, gene, annotazione, numero di accesso
  (Accession; it.wikipedia GenBank), locus (plural loci), filamento (Strand:
  filamento +, filamento −, filamento diretto / inverso), complementare inverso,
  contenuto in GC, profondità (Depth), profondità di lettura (read depth in
  sentences), copertura, identità, similarità, omologia, ortologo, paralogo,
  allineamento, blocco colineare, codice genetico, organello, organismo, ceppo,
  specie, prodotto, localizzazione (Location), amminoacido, nucleotide,
  a monte / a valle (upstream / downstream), gruppo di similarità.
- **Diagram terms**: traccia (Track), corsia (Lane), slot (track slot), posizione
  (File card slot), tacca (Tick; principale / secondaria), righello (Ruler),
  barra di scala, legenda, voce della legenda, anello, banda, nastro, margine
  (Padding), spazio (Gap), spaziatura (Spacing), contorno (Stroke), riempimento,
  campione di colore, sfondo / primo piano, area di disegno (Canvas), livello
  (Layer), preimpostazione (Preset), predefinito (Default).
- **Command words** (Microsoft Italian): Annulla (Cancel and Undo), Ripeti (Redo),
  Elimina (Delete), Rimuovi (Remove), Cancella (Clear), Reimposta (Reset),
  Ripristina (Restore), Sfoglia (Browse), Carica (Load and Upload), Scarica,
  Importa, Esporta, Seleziona, Scegli, Mostra / Nascondi, Abilita / Disabilita
  (verbs), Attivato / Disattivato (On / Off), Riduci / Espandi (Collapse /
  Expand), Comprimi (Compress), Salta (Skip), Unisci (Merge), Mantieni (Keep),
  Genera, Rigenera, Fai clic.
- **Loanwords kept in English** because Italian bioinformatics and software keep
  them (lower case, invariable plural, an article of the stated gender): feature
  (f.), record (m.), hit (m.), query (f.), subject (m.), qualifier (m.), mismatch,
  skew, locus tag, slot (m.), card (f.), layout, popup, file, input / output,
  browser, runtime, thread, cache, hash, pipeline, tutorial, zoom, reciprocal
  best hit, near-reciprocal hit, coortholog, inparalog, orthogroup (§13.5). Never
  mix one with a translation of it (no caratteristica beside feature).
- **Never use**: cliccare / clicca (write Fai clic), settare / setta,
  resettare, uploadare, downloadare, editare, linkare, loggare, customizzare,
  accessione, etichetta di locus, caratteristica or funzionalità for Feature,
  registro or registrazione for Record, elica or trefolo for strand, mappa for
  Diagram, inclinazione or asimmetria for skew, opzionale, pb or kpb for bp,
  Si prega di, Per favore.

### 13.2 Buttons and sentences

- Register: tu, implicit. Never Lei (Faccia clic, Scelga), never the impersonal
  infinitive of technical manuals (Fare clic su, Selezionare un file), never
  Si prega di. Source: Microsoft Italian Style Guide examples (Vuoi continuare?,
  fai clic su Annulla, Per confermare … guarda la riproduzione).
- Buttons, menu commands, links, tooltips of actions: imperative of tu, no final
  period, no article where the English has none: Salva sessione, Carica esempio,
  Reimposta impostazioni, Genera diagramma. Checkboxes follow the English
  (Mostra la profondità, Scala logaritmica).
- Labels, headings, tabs, column headers, badges: noun phrases without a final
  period (Colore del riempimento). Options are nouns or adjectives (Automatico,
  Attivato, Seriale). Descriptions: a noun phrase or the third person present
  (Prefisso per i file di output. Mostra la profondità di lettura come traccia.).
- Instructions (help, notices, dialogs, errors): imperative of tu. Seleziona un
  record. Fai clic su Genera. Lascia vuoto per usare il valore automatico.
  A dialog title that asks takes the infinitive (Rimuovere l’intero file?); a
  question in the body: Vuoi …?
- Progress: noun and … (Generazione del diagramma…, Salvataggio della
  sessione…, Caricamento della sessione…); in a sentence …in corso.
  Completion: noun and participle (Sessione salvata. Allineamento applicato.).
- Errors: what failed, then what to do. Cannot: Impossibile + infinitive
  (Impossibile leggere il file. Controlla che sia un file FASTA.); X failed is
  Errore di X. Busy: Generazione del diagramma in corso. Riprova al termine.
- Required is Obbligatorio, optional (facoltativo) (invariable, after the field
  name), recommended Consigliato; non valido, non supportato, non disponibile.
- An option or button name is not the subject: Con Automatico, …; Al salvataggio, ….

### 13.3 Punctuation and spacing

- Quotes “ ” (U+201C, U+201D), nested ‘ ’. A {name} the English quotes is
  “{name}”; a UI word the English quotes ('Bar') is “Barra”. Never « », „ or
  straight quotes in prose. Punctuation follows the closing quote.
- Apostrophe: the typographic ’ (U+2019) in every elision (l’asse, dell’asse,
  un’annotazione); the straight ' only in code. An accent is a real letter: è,
  perché, più, così, già, può, attività, visibilità, È at a sentence start;
  never e', perche', piu', E'.
- Colon: no space before, one after (Campo: {path}); after a colon a lower-case
  word unless a name, control label or placeholder follows. No space before
  ; ? ! (unlike fr): Vuoi continuare?
- Parentheses, slashes and ranges as the source writes them: Regione
  (facoltativo), Organismo / ceppo, Righe 3–4. A spaced hyphen is not a dash.
  Sentences are joined by an ordinary space (D-09); no non-breaking space.
- Digits and separators stay en-US (D-05): 0.65, 1,000.5, never 0,65 or
  1.000,5; units and % keep the source spacing (100 bp, 20px, 65%). One …
  (U+2026), no space before it. ad es. for e.g., cioè for i.e. No comma before
  e / o; ed / ad only before the same vowel.

### 13.4 Words and capitals

- Sentence case for every label, button, heading, tab and option, as in Microsoft
  Italian: the first word and names only (Slot di traccia personalizzati).
- A key the English writes in capitals stays in capitals with accents:
  SCHEMI PREIMPOSTATI, COLORI PREDEFINITI, REGOLE SPECIFICHE, FILTRO DELLE
  ETICHETTE, PRIORITÀ DEI QUALIFIER, AGGIUNGI CHIAVE; the badge OFF is
  DISATTIVATO.
- The action Generate in running text is lower case (genera, alla generazione,
  prima di generare); only the button is cited as Genera, without quotes
  (fai clic su Genera), or “Genera” where the English quotes it. A cited control
  keeps its translated label and capital: in Slot di traccia personalizzati.
- Compounds put the head noun first, joined with di: tipo di feature, TSV di
  profondità. Genders: la feature, la CDS, la query, la card; il record, il file,
  l’hit, il locus, lo slot, il subject, il popup.
- An exact-match condition is parenthesized after the field: Prodotto
  (corrispondenza esatta): {product}. A count of members is Numero di …
  (Numero di geni del gruppo), not a bare plural.

### 13.5 Placeholders and counts

- A placeholder value has no known gender. Put a noun before it and make the
  article agree with the noun: il file {file}, il record {record}, la traccia
  {track}, lo slot {id}, l’impostazione {setting}. Never glue an elision,
  suffix or prefix to a placeholder (d’{name}, {name}i); never put an article or
  an articulated preposition (del, nella) straight before one. A placeholder that
  holds a translated phrase ({action}, {reason}) follows a colon.
- `Intl.PluralRules('it')` gives `one` for 1 only. 0, 1.5, 21, 22, 25 and 1000
  take `other`; 1000000 gives `many`, which falls back to `other`. The `one`
  string drops {count} only where the English `one` does. Kept English nouns do
  not change in the plural (1 record, 0 record, 2 record).
- Worked pairs, `one` / `other`:
  1. ({count} record) / ({count} records): ({count} record) / ({count} record)
  2. {count} feature ready / {count} features ready: {count} feature pronta /
     {count} feature pronte
  3. From this color rule: / From these {count} color rules: → Da questa regola
     di colore: / Da queste {count} regole di colore:
  4. Remove {count} unmatched feature edit / edits → Rimuovi {count} modifica
     di feature senza corrispondenza / Rimuovi {count} modifiche di feature …
  5. {count} LOSAT result is cached and reused … → {count} risultato LOSAT è in
     cache e viene riutilizzato quando input e impostazioni corrispondono. /
     {count} risultati LOSAT sono in cache e vengono riutilizzati quando input e
     impostazioni corrispondono.

### 13.6 Length

- Italian runs 15–25% longer than English; help texts and notices wrap. Rail
  short names fit about nine letters in the 52 px entry and break only between
  words: Input, Confronto, Avanzate, Regioni, Layout, Carattere, Titoli,
  Legenda, Colori, Feature, Etichette, Asse, GC / skew, Info (About), with Genera
  beside the rail. Carattere and Etichette (9 letters) are the widest. The rail
  keys of Features, Layout and Labels are Feature, Layout and Etichette.
- Header buttons stay within about 1.4× the English label: Salva sessione,
  Carica sessione, Carica esempio, Reimposta impostazioni (1.6×: check at 1,024
  px). Options in 10 px rows use the short forms: Automatico, Seriale,
  Multithread, Sopra / Sotto, Dentro / Fuori, Sull’asse. Drop an article before
  switching terms.

### 13.7 Machine checks

`tests/web/i18n-lint/it.tsv` (draft `lint-it.tsv`, 25 rules). Errors: `it-register`
(Lei forms, Fare clic, Si prega di), `it-never-words` (clicca, settare,
resettare, uploadare, editare, accessione, etichetta di locus), `it-skew-translated`
(inclinazione / asimmetria with GC or AT), `it-optional` ((facoltativa),
opzionale), `it-quotes` (« » „ ‚), `it-apostrophe` (straight apostrophe, e' or e’
as an accent), `it-missing-accent` (più, così, perché, già, può, -ità nouns,
future forms), `it-decomposed` (not NFC), `it-digits` (non-ASCII digits),
`it-number-format` (decimal comma, dot or space grouping), `it-space-before-punct`
(a space before : ; ? !), `it-ellipsis` (three dots, or a space before …),
`it-placeholder-suffix` (suffix, prefix or elision glued to {name}),
`it-loanword-plural` (English -s on feature, record, hit, query, file),
`it-auto-label` (the label Auto), `it-bp-unit` (pb, kpb, Mpb). Warnings:
`it-english-ui-words` (upload, download, reset, settings, click, track, lane,
legend, preset, default, tick, ribbon, strand), `it-feature-words`
(caratteristica, funzionalità), `it-record-words` (registro, registrazione),
`it-strand-words` (elica, trefolo, strand), `it-straight-quotes`,
`it-placeholder-article` (an article straight before a placeholder),
`it-hyphen-dash` (a spaced hyphen as a dash), `it-title-case` (three or more
capitalized words in a label), `it-diagram-map` (mappa for diagram).

## 14. Polish (`pl`)

Standard Polish in the register of Polish software (Microsoft, Google, Apple): the bare imperative
without a pronoun, no Pan/Pani. The Language menu entry is the endonym **Polski**; the menu label is
Język; the page `lang` is `pl`. Numbers keep the English format (D-05). Translations are drafts; nothing
states that they are reviewed.

### 14.1 Vocabulary

- **Field terms**: genom, sekwencja, rekord (Record; numer akcesyjny for Accession), cecha (Feature; typ
  cechy; INSDC "tabela cech"), adnotacja, kwalifikator (the keys /gene, /product stay), locus (plural
  loci), nić (Strand: nić +, nić −, nić dodatnia), odwrotna komplementarna, zawartość GC, głębokość
  (Depth), głębokość odczytów (read depth), pokrycie (coverage), identyczność, podobieństwo, homologia,
  ortolog, paralog, blok kolinearny, kod genetyczny, organellum (plural organella), organizm, szczep,
  gatunek, produkt, gen.
- **Alignment has two words**: dopasowanie for sequences (BLAST, LOSAT, długość dopasowania; a drawn
  match is also dopasowanie), wyrównanie for the layout action of Similarity groups (Wyrównaj, Przejrzyj
  wyrównanie, wyrównanie według podobieństwa). **Hit** is trafienie.
- **Diagram terms**: ścieżka (Track; plural ścieżki, ścieżek), tor (Lane), slot (track slot), miejsce
  (File card slot), karta (File card), zakładka (Tab), znacznik podziałki (Tick; główny / pomocniczy
  znacznik), oznaczenie (annotation mark), linijka (Ruler; reguła is a color rule), pasek skali (Length
  bar), margines (Padding), odstęp (Gap), rozstaw (Spacing), interwał, kontur (Stroke), wypełnienie
  (Fill), próbka koloru, podkład / pierwszy plan (Underlay / Foreground), ustawienie wstępne (Preset),
  domyślny (Default).
- **Command words** (imperative): Zapisz, Wczytaj (Load), Importuj, Eksportuj, Prześlij (Upload),
  Pobierz (Download), Wygeneruj ponownie, Zastosuj, Resetuj, Cofnij, Ponów, Anuluj, Usuń (Delete and
  Remove), Wyczyść (Clear), Dodaj, Zmień nazwę, Przeglądaj (Browse), Zaznacz (Select), Wybierz (Choose),
  Zwiń, Rozwiń, Pomiń (Skip), Wyklucz (Omit), Scal (Merge), Zadokuj, Przypnij, Zachowaj (Keep), Obróć,
  Spróbuj ponownie; states Włączone / Wyłączone (Wł. / Wył. on badges), Automatycznie.
- **Kept in English** (D-12, evidence in the glossary notes): GC skew, AT skew, skew, query, subject,
  locus tag, orthogroup, coortholog, inparalog, reciprocal best hit, near-reciprocal hit, upstream,
  downstream, open source, cross-origin isolated and the names of §1.1. They are lower case in prose,
  take no Polish case ending and no plural (sekwencja query; rephrase a count). Loanwords in use: slot,
  diagram, interfejs, piksel, intron, replikon, preprint, hosting.
- **Never use**: numer dostępu (numer akcesyjny), pasmo, łańcuch, włókno (nić), mapa (diagram), trasa,
  utwór (ścieżka), kwerenda (query), wgraj, ściągnij, załaduj, ładowanie (Prześlij, Pobierz, Wczytaj),
  wartość E, wynik bitowy, skos GC, nachylenie GC, asymetria GC, znacznik locus, Pan, Pani, proszę,
  kołowy (kolisty), and the English feature, track, preset, lane, popup, default, reset, record, legend,
  label.

### 14.2 Buttons and sentences

- Register (Microsoft, Google, Apple Polish): the second person singular imperative with no pronoun; no
  Pan/Pani/Państwo, no proszę, no infinitive on a button. Where a possessive is unavoidable write twój
  in lower case; usually drop it (Zastąp zmienione kolory, not Zastąp twoje zmienione kolory).
- Buttons, menu commands, links: perfective imperative: Zapisz sesję, Wczytaj przykład, Zmień nazwę,
  Zastosuj i generuj. Generate is the imperfective Generuj (Generuj diagram); Regenerate is Wygeneruj
  ponownie.
- Labels, headings, tabs, column headers, badges: a noun phrase without a final period (Zawartość GC,
  Typ cechy). Options are nouns, adverbs or adjectives that agree with the setting (Automatycznie,
  Włączone, Szeregowo). A dialog title that asks is the perfective infinitive + ?: Usunąć cały plik?,
  Zresetować wszystkie ustawienia?
- Standalone state chips use the neuter singular (Dostępne, Włączone, Gotowe, Zalecane, Wymagane);
  inside a phrase the adjective agrees with its noun (slot włączony, ścieżka włączona). Glossary
  adjectives are listed in the masculine.
- Required / optional / recommended: Wymagane, (opcjonalnie) after a field name (invariable), Zalecane.
- Instructions (help, notices, dialogs, errors): the same imperative. Wybierz rekord. Kliknij Generuj.
  Wygeneruj diagram ponownie. Leave blank for auto: Pozostaw puste, aby użyć wartości automatycznej.
  Name a button as an action: kliknij Resetuj.
- Descriptions: third person or the impersonal reflexive (Pokazuje głębokość odczytów. Stosuje się przy
  generowaniu.).
- Progress (…ing): the verbal noun + … (Zapisywanie sesji…, Generowanie diagramu…). Completion: the
  impersonal -no / -to form in a sentence (Zapisano sesję. Zastosowano wyrównanie.), the participle as a
  status label (Zapisane).
- Errors: what failed, then what to do. Nie można odczytać pliku. Sprawdź, czy to plik FASTA. Nie można
  … says the action is impossible; Nie udało się … says an attempt failed (Nie udało się zapisać
  sesji.). Headings: Błąd generowania, Błąd eksportu SVG. Busy reason: Trwa generowanie diagramu.
  Spróbuj ponownie po zakończeniu. Retry is Spróbuj ponownie.

### 14.3 Punctuation and spacing

- Quotes „…” (U+201E opens, U+201D closes), nested «…». A `{name}` the English quotes is „{name}”; a UI
  word the English quotes ('Bar') is „Pasek”. No straight " and no “ ‘ ’ in prose. A period or comma
  follows the closing ”.
- Colon: no space before, one after (Pole: {path}). A label that ends with a colon keeps it. No space
  before ? ! , . ;
- Dashes: the en dash with spaces for a pause (Plik – brak), without spaces for a range (Wiersze 3–4);
  never a spaced hyphen. The minus − and → × stay.
- ASCII digits with the source's decimal point and grouping (0.65, 1,000), never 0,65 or 1 000. Units
  and % keep the source spacing (100 bp, 50%). One plain space everywhere: the catalog uses no
  non-breaking space, although print Polish puts one after single-letter words (a, i, o, u, w, z) and
  between a number and its unit; the tests compare plain text.
- One … (U+2026) directly after the word. e.g. is np.; i.e. is tzn.; etc. is itd. Sentences are
  separated by one space (§1.3). Diacritics are NFC: ą ć ę ł ń ó ś ź ż are precomposed.

### 14.4 Words and capitals

- Sentence case for every label, button, heading, tab, option and column header: the first word and
  proper names only (Niestandardowe sloty ścieżek, Porównania parami). Never copy English Title Case.
- A key written in capitals stays in capitals, with diacritics: KOLORY DOMYŚLNE, REGUŁY SZCZEGÓŁOWE,
  FILTROWANIE ETYKIET, PRIORYTET KWALIFIKATORÓW, DODAJ KLUCZ; the badge OFF is WYŁ.
- The action Generate in running text is the verb in lower case (przy generowaniu, wygeneruj ponownie);
  only the button itself is cited with its label Generuj (capital, no quotes, „Generuj” where the
  English quotes it).
- Nouns the English capitalizes as product terms are lower case in a sentence: sesja, wynik, rekord. A
  control, section or option cited in a sentence keeps its translated label and capital without quotes
  (kliknij Generuj; w sekcji Etykiety; wybierz Automatycznie).
- Compounds put the head noun first and use the genitive: typ cechy, TSV głębokości, plik BLAST TSV,
  pamięć podręczna LOSAT. Genders: ta CDS (ta sekwencja), ten TSV, ten ID, ten locus; query and subject
  have no ending.

### 14.5 Placeholders and counts

- A placeholder is never inflected and carries no suffix. Put a class noun before it, in the case the
  sentence needs, and let adjectives agree with that noun (rekord {record}, ścieżka {track}, plik
  „{file}”); never write {name}u or {n}-ego. {name} is hidden becomes Element {name} jest ukryty.
- Categories (Intl.PluralRules('pl')): one is 1 only; few is 2–4, 22–24, 32–34 … (not 12–14); many is 0,
  5–21, 25–31 … and 112; other is a fraction (1.5). The `one` key holds the text for exactly one (count 1); the `other` key holds an object with one entry per plural category of the other counts: `few`, `many` and `other` (fractions); the runtime shows the object's category text for every count but 1. 21 takes many (21 rekordów), 22 few (22 rekordy), 25
  many, 0 many (0 rekordów), 1.5 other (1.5 rekordu: the genitive singular). A `one` string drops
  {count} only where the English one does.
- After a numeral the noun is nominative plural for few and genitive plural for many. In a frame that
  already asks for the genitive (z, od, do, dla, po liczbie) the noun is genitive after one as well,
  singular, and genitive plural after few and many. A verb after few is plural; after many it is the
  third person singular neuter (5 edycji wskazuje); in the other case it is singular. Prefer an
  impersonal verb (Pominięto 5 wierszy) to avoid the agreement.
- Worked sets, `one` / few / many / other (English key in parentheses):
  1. ({count} records): ({count} rekord) / ({count} rekordy) / ({count} rekordów) / ({count} rekordu)
  2. ({count} features ready): {count} cecha gotowa / {count} cechy gotowe / {count} cech gotowych /
     {count} cechy gotowej
  3. ({current} of {count} matches): {current} z {count} dopasowania / {current} z {count} dopasowań /
     {current} z {count} dopasowań / {current} z {count} dopasowania (z takes the genitive, so few and
     many are equal)
  4. (From these {count} color rules:): Z tej reguły koloru: / Z tych {count} reguł koloru: / Z tych
     {count} reguł koloru: / Z tych {count} reguł koloru: (tych needs the plural, so other equals many)
  5. ({count} feature edits name features the source does not have): {count} edycja cechy wskazuje
     cechę, której nie ma w źródle. / {count} edycje cech wskazują cechy, których nie ma w źródle. /
     {count} edycji cech wskazuje cechy, których nie ma w źródle. / {count} edycji cech wskazuje cechy,
     których nie ma w źródle.
- A string without plural forms that counts things is a noun phrase and a colon: Etykiety: {count}.
  Never build a plural with (y) or (s).

### 14.6 Length

- Polish runs 10–30% longer than English; help texts and notices wrap. Rail short names, at most 9
  letters in the 52 px entry (10 px semibold): Wejście, Porównaj (Compare), Zaawans. (Advanced;
  Zaawansowane has 12 letters and does not fit; the period is the accepted abbreviation), Regiony,
  Układ, Czcionka, Tytuły, Legenda, Kolory, Cechy, Etykiety, Oś, GC / skew, O programie (two lines,
  between words); Generuj beside the rail. A word is never hyphenated.
- Header buttons stay within about 1.35× the English label: Wczytaj przykład (16 for 12), Wczytaj sesję,
  Zapisz sesję, Cofnij, Ponów, Resetuj ustawienia (18 for 14); check the first and the last at 1,024 px.
  Generate bar: Generuj diagram, Anuluj, Generowanie….
- Options in 10 px rows use the short forms: Automatycznie / Auto, Wł. / Wył., Szeregowo, Wielowątkowo,
  Nad / Pod, Na osi, Na zewnątrz / Wewnątrz, Wstążka / Krzywa. To shorten a label drop a genitive
  attribute first (Szerokość ścieżki, not Szerokość tej ścieżki); never switch terms.

### 14.7 Machine checks

The lint list is `tests/web/i18n-lint/pl.tsv`; each rule reads one catalog value
with the `{name}` placeholders masked as `{X}`. `error` fails, `warn` is printed.

- Typography: `pl-nfc` (error) a combining mark, so not NFC; `pl-ellipsis-dots` (error) three periods;
  `pl-ellipsis-space` (error) a space before …; `pl-straight-quotes` (error) a straight double quote;
  `pl-wrong-quotes` (error) “ ‘ ’ ‚ or a ” before any „; `pl-space-before-punct` (warn) a space before ,
  . ; : ! ?; `pl-hyphen-dash` (warn) a spaced hyphen; `pl-nbsp` (warn) a non-breaking or thin space.
- Numbers and script: `pl-decimal-comma` (warn) 0,65; `pl-thousands-space` (warn) 1 000;
  `pl-foreign-script` (error) Cyrillic, CJK, Thai, Arabic or full-width characters.
- Register: `pl-formal-address` (error) Pan, Pani, Państwo, proszę; `pl-pronoun` (warn) ty, ciebie,
  twój; `pl-infinitive-button` (warn) a value that opens with an infinitive and is not a question;
  `pl-title-case` (warn) two capitalized words.
- Terms: `pl-english-words` (warn) feature, track, preset, lane, popup, default, reset, record, legend,
  label, strand, session, settings, canvas; `pl-strand-words` (warn) pasmo, łańcuch, włókno;
  `pl-accession` (warn) numer dostępu; `pl-kept-latin` (error) the Polish renderings of GC skew, locus
  tag, E-value, Bitscore, query, subject and reciprocal best hit; `pl-map-track-words` (warn) mapa,
  trasa, utwór; `pl-vendor-verbs` (warn) wgraj, ściągnij, załaduj, ładowanie.
- Words and placeholders: `pl-placeholder-suffix` (error) a letter glued to a placeholder;
  `pl-plural-hack` (error) rekord(y) or file(s); `pl-ascii-diacritics` (error) błąd written as blad,
  ścieżka as sciezka, and similar.

## 15. Turkish (`tr`)

Standard Turkish of Turkey, as Microsoft tr-TR and the TDK spelling guide write it. Turkish is agglutinative: suffixes
carry case and possession, vowel harmony picks their form, and a suffix can never be attached to a `{placeholder}`,
because the value is unknown. §15.5 shows how to rephrase around that. The Language menu entry is **Türkçe**; the
page `lang` is `tr`; the menu label is Dil. This section is a draft and nothing in it states that it is reviewed.

### 15.1 Vocabulary

- **Field terms**: genom, dizi (sequence), kayıt (Record), özellik (Feature; özellik tipi, never tür, which is the
  Species), anotasyon, erişim numarası (Accession), gen, ürün, protein, nükleotit, amino asit, organizma, tür (Species),
  suş, izolat, organel, replikon, iplik (Strand: ileri iplik; tr.wikipedia DNA and Genetik kod write çift iplikli DNA),
  ters tamamlayıcı (Reverse complement), GC içeriği, okuma derinliği (Read depth), kapsama, özdeşlik (Identity),
  homoloji, ortolog, paralog, hizalama (Alignment), benzerlik grubu (Similarity group), kolinear blok, genetik kod,
  topoloji, sorgu (Query), yuva (track slot), iz (Track), şerit (Lane), halka (Ring), gösterge (Legend, Microsoft
  Excel Turkish), işaret (Tick; ana işaret, ara işaret), cetvel (Ruler), kontur (Stroke), dolgu (Fill), kenar
  boşluğu (Padding), boşluk (Gap), aralık (Spacing and Interval), ofset (Offset), kaydırma (Pan), hazır ayar (Preset).
- **Command words** (Microsoft tr-TR): Kaydet, Yükle (Load), Karşıya yükle (Upload), İndir, İçeri aktar, Dışarı aktar,
  Sil (Delete), Kaldır (Remove), Temizle (Clear), Sıfırla (Reset), Geri yükle (Restore), Seç (Select and Choose),
  Gözat, Ekle, Uygula, İptal, Kapat, Kopyala, Çoğalt (Duplicate), Yeniden adlandır, Birleştir, Atla (Skip), Hariç tut
  (Omit), Göster, Gizle, Etkinleştir, Devre dışı bırak, Genişlet, Daralt, Geri al, Yinele, Yeniden dene, Oluştur,
  Yeniden oluştur; Açık / Kapalı (On / Off); Etkin / Devre dışı (Enabled / Disabled).
- **Kept in English** (no settled Turkish form; lower case in running text): anchor, qualifier, locus tag, hit,
  subject, skew (GC skew, AT skew), seed, reciprocal best hit, near-reciprocal hit, coortholog, inparalog,
  orthogroup, upstream / downstream, flat file, hash, lookbehind, Halo, and the names of §1.1. The Latin word
  never takes a suffix: write a compound with the head noun (§15.4).
- **Two words for one English word**: Hizalama is both the sequence alignment and the layout action; the layout
  method is "benzerliğe göre hizalama". Kaynak is both Source and Resource. Konum is both Position and Location.
  Çöz is the button Resolve; belirle is "determine automatically".
- **Never use**: zincir for a DNA strand, öznitelik / nitelik for feature, çapa for anchor, lokus etiketi, opsiyonel,
  tavsiye edilen, giriş for input (it is login; write girdi), font (write yazı tipi), lütfen, the informal
  imperative tıkla in a sentence, the -mak / -mek infinitive as a label, English plural -s (hits, anchors).

### 15.2 Buttons and sentences

- Register: the polite plural imperative in instructions, help and errors (Dosyayı seçin. Oluştur düğmesine
  tıklayın. Tanılama bilgilerini kopyalayın.), as Microsoft tr-TR writes it. Never the singular sen form
  (seçebilirsin), and no "lütfen". Descriptions are third person present or impersonal (Etiketleri gösterir.
  Oluşturulurken uygulanır.).
- Buttons, menu commands, links: the bare imperative, which is also Microsoft's form (Kaydet, Uygula, İptal,
  Oturum kaydet, Örnek yükle). A single generic object takes no accusative (Oturum yükle); a specific whole takes it
  (Ayarları sıfırla, Tüm kayıtları seç). Never the infinitive (Kaydetmek).
- Labels, headings, tabs, column headers, badges: a noun phrase without a final period, with the possessive
  compound Turkish needs (Yazı tipi boyutu, Etiket görünürlüğü). Options are nouns or adjectives (Otomatik, Etkin).
  A dialog title that asks ends in mI: Bu dosya kaldırılsın mı? A dialog title that names its command stays the
  command, as the button it opens from (Varsayılan renkleri sıfırla, Kayıtları seçin): the -mA noun (sıfırlama)
  would read as the negative imperative ("do not reset").
- Required, optional, recommended: Gerekli, İsteğe bağlı, Önerilen. A field marker is (isteğe bağlı).
- Progress (…ing): -iliyor / -ıyor + … (Oturum kaydediliyor…, Diyagram oluşturuluyor…). Completion: -ildi in a
  status line (Oturum kaydedildi.; Alignment applied: Hizalama uygulandı.).
- "Cannot": a failure that happened is -amadı / -ilemedi (Dosya okunamadı.); a rule is -ilemez / -amaz (Etiket boş
  olamaz.); a feature that is off is Kullanılamıyor.
- Errors: what failed, then what to do: Dosya okunamadı. Dosyanın FASTA olduğunu denetleyin. Headings: Oluşturma
  hatası, SVG dışa aktarma hatası. Retry is Yeniden dene on a button and yeniden deneyin in a sentence.
- A control is cited by its label plus the noun düğme, sekme or bölüm, never by a suffix on the label: Oluştur
  düğmesine tıklayın; Gösterge sekmesi. A setting is not the subject: write "Otomatik seçeneğinde, …", not
  "Otomatik …seçer". Where the English quotes a label, use “ ”.

### 15.3 Punctuation and spacing

- Sentences end with a period; labels and buttons do not. Colon: no space before, one after (Etiket: {path}). A
  question mark or exclamation mark has no space before it.
- Quotes “ ” (U+201C, U+201D), nested ‘ ’. A {name} the English quotes is “{name}”. No straight " in prose, and no
  « » or „ “.
- The apostrophe is ’ (U+2019) and is used only where a suffix must follow a literal Latin name in static text
  (GenBank’ten); prefer a compound (GenBank dosyası). Never an apostrophe after a placeholder or a number.
- ASCII digits with the English point and grouping (0.65, 1,000); never 0,65. Units and % keep the source spacing
  (100 bp, 50%, 8 pt) although Turkish writes %50 in prose.
- One … (U+2026) with no space before it. Ranges keep the en dash (Satır 3–4); a spaced hyphen is not a dash.
  Slashes as the source writes them (Üstte / önde). No non-breaking spaces. Sentences are joined by one space (§1.3).

### 15.4 Words and capitals

- Sentence case for every label, button, heading and tab: the first word and proper names only (Özel iz yuvaları,
  İkili karşılaştırmalar). Never copy English Title Case.
- İ and I: Turkish has dotted i / İ and dotless ı / I. The capital of i is İ (İptal, İşlem, İndir, İz) and the
  lower case of I is ı (Kayıt, Sıfırla). Only a few words start with a plain I (Isı, Işık, Irmak). A capital
  inside a lower-case word is an error (kayIt, kayİt).
- Keys the English writes in capitals stay in capitals, with İ where the word has i: HAZIR ŞEMALAR, VARSAYILAN
  RENKLER, BELİRLİ KURALLAR, ETİKET FİLTRELEME, QUALIFIER ÖNCELİĞİ, ETİKET GEOMETRİSİ, ANAHTAR EKLE, RENK; the
  badge OFF is KAPALI. CSS uppercase on lang tr maps i to İ, but a typed value must do it itself.
- The action Generate in running text is the lower-case noun or verb form (oluşturma, oluşturulurken, oluşturmadan
  önce); only the button itself is cited by its label Oluştur, followed by düğmesi.
- Nouns the English capitalizes as product terms are lower case in a sentence: oturum, sonuç, kayıt.
- Compounds put the modifier first and mark the head with the possessive: özellik tipi, kayıt etiketi, yazı tipi
  boyutu, qualifier önceliği, locus tag değeri. A kept English word is always the modifier, so it needs no suffix.
- Diacritics are never dropped: ı İ ç Ç ğ Ğ ö Ö ş Ş ü Ü. A word typed without them (sifirla) is an error.

### 15.5 Placeholders and counts

- A `{placeholder}` is never inflected: no case, possessive, plural, apostrophe or ordinal dot after it (`}'yi`,
  `}'nı`, `}ı`, `}’ı`, `{n}.`). Vowel harmony depends on a value the translator cannot see. Rephrase:
  1. Put a class noun after the placeholder: {file} dosyası, {record} kaydı, {track} izi, {id} kimliği, {section} bölümü.
     The suffix sits on the class noun, which is known.
  2. Use a colon: Şuraya git: {section}. Yeniden sırala: {name}. Seçili kayıtlar: {records}.
  3. Name the thing with adlı: “{name}” adlı yardımcı dosya oluşturulamadı.
  4. Use a postposition that takes no suffix (için, ile, kadar, önce, sonra): {count} kayıt için; {name} ile.
  5. A pair is an arrow: #{from} → #{to}. An ordinal is a class noun: Dosya {n} (not {n}. dosya).
  6. A quoted placeholder is followed by a space and a noun, never by a suffix: “{name}” dosyası.
- {count} precedes a singular noun (a noun after a numeral is singular): {count} kayıt, never {count} kayıtlar. Turkish
  has `one` and `other` only; Intl.PluralRules('tr') returns one for 1 alone, and 0 or 1.5 take other. The two forms
  are usually the same string (D-09). Because `one` is exactly 1, the `one` value may be a singular where the English
  `one` form names a single thing (Seçilen özellik zaten sola bakıyor. / Seçilen {count} özelliğin tümü zaten sola
  bakıyor.). A string that counts things without
  `tn` is a noun phrase and a colon: Gerekli sütun sayısı: {count}.
- Worked pairs, `one` / `other` (English `other` key in parentheses):
  1. ({count} source records inspected): {count} kaynak kayıt incelendi / {count} kaynak kayıt incelendi
  2. ({drawn} of {count} records drawn): Çizilen kayıt: {drawn} / {count} (both forms)
  3. ({count} issues to fix before Generate): Oluşturmadan önce düzeltilecek {count} sorun (both forms)
  4. (Auto ({count} parallel runs)): Otomatik ({count} paralel çalıştırma) (both forms)
  5. ({current} of {count} matches): Eşleşme {current} / {count} (both forms)

### 15.6 Length

- Turkish runs 10 to 25% longer than English; help texts and notices wrap. Compounds are long words and are never
  hyphenated. Prefer the short native word (Düzen, Eksen, Girdi) and drop the possessive on a button where the
  noun stays clear.
- Rail names (52 px entries; one word may run into the padding, a name wraps only between words): Girdi, Karşılaştır
  (11 letters; if it overflows at 1,024 px use the rail key Kıyasla), Gelişmiş, Bölgeler, Düzen, Yazı tipi, Başlıklar,
  Gösterge, Renkler, Özellikler (10 letters; fallback Özellik), Etiketler, Eksen, GC / skew, Hakkında; Oluştur
  beside the rail. Rail keys that need their own value: `rail\u0004Features`, `rail\u0004Layout`, `rail\u0004Labels`.
- Header buttons stay within about 1.3 times the English label: Oturum kaydet, Oturum yükle, Örnek yükle,
  Ayarları sıfırla, Geri al, Yinele, Dil. Options in 10 px rows use the short forms: Otomatik, Seri, Yukarı / Aşağı,
  Dışta / İçte, Eksen üzerinde, Arka plan, Kurdele / Eğri.

### 15.7 Machine checks

The rule list is `tests/web/i18n-lint/tr.tsv`; each rule runs on one catalog value and `error` fails the test.

- Suffixes on placeholders: `tr-placeholder-apostrophe-suffix` (error) `}'yi`, `}'nı`; `tr-placeholder-glued-suffix` (error)
  `}ı`, `}yi`; `tr-placeholder-quote-suffix` (error) `”ı` after a quoted placeholder; `tr-count-ordinal-dot` (warn) `{n}.`;
  `tr-count-plural-after-numeral` (warn) `{count} kayıtlar`; `tr-apostrophe-on-kept-word` (warn) anchor’ı.
- Case and diacritics: `tr-capital-i-inside-word` (error) kayIt; `tr-capital-i-initial` (warn) Iptal;
  `tr-caps-heading-dotless` (error) ETIKET; `tr-ascii-folded-words` (error) sifirla, yukle.
- Register: `tr-informal-sen` (error) sen and -sin forms; `tr-please-word` (warn) lütfen; `tr-infinitive-label` (warn)
  a value ending in -mak / -mek.
- Vocabulary: `tr-english-word-left` (warn) track, feature, layout …; `tr-retired-terms` (warn) zincir, öznitelik, çapa,
  opsiyonel, giriş, font; `tr-english-plural-s` (warn) hits, anchors.
- Punctuation: `tr-straight-double-quote` (error); `tr-foreign-quotes` (error) « » „; `tr-straight-apostrophe` (warn);
  `tr-space-before-punct` (error) before : ; ? !; `tr-ellipsis` (error) ... or a space before …; `tr-spaced-hyphen-dash`
  (warn); `tr-edge-space` (error) leading, trailing or doubled whitespace.
- Numbers: `tr-nonascii-digits` (error); `tr-decimal-comma` (warn) 0,65.

## 16. Russian (`ru`)

Standard Russian as Microsoft, Apple and Google Russian UIs write command words
and as Russian molecular-biology writing has the science. The Language menu
entry is **Русский**; the menu label is Язык; the page `lang` is `ru`.

### 16.1 Vocabulary

- **Field terms**: геном, последовательность, запись (Record; ru.wikipedia
  GenBank), элемент (Feature; тип элемента), квалификатор (keys /gene,
  /product stay), аннотация, номер доступа (accession), локус (plural локусы),
  цепь (Strand; ru.wikipedia GFF), обратный комплемент, GC-состав, глубина
  (Depth), глубина покрытия (read depth in sentences), покрытие, идентичность,
  сходство, гомология, ортолог, паралог, ортогруппа, коллинеарный блок,
  генетический код, органелла, организм, штамм, вид, продукт, ген.
- **Feature is элемент**, an ordinary lower-case noun. A list item is пункт
  (пункт легенды), never элемент; элементы управления is Microsoft's Control.
- **Diagram terms**: трек (Track), дорожка (Lane), слот (track slot), позиция
  (File card slot), деление (Tick; основное / промежуточное), метка (annotation
  mark), линейка (Ruler; правило is a color rule), масштабный отрезок (Length
  bar), отступ (Padding), промежуток (Gap), расстояние (Spacing), заливка,
  обводка (outline of a shape; the axis line is линия оси), образец цвета,
  подложка, предустановка (Preset), по умолчанию, макет (Layout), подпись
  (Label), легенда.
- **Command words** (Microsoft Russian): Применить, Добавить, Удалить (Delete),
  Убрать (Remove), Очистить (Clear), Сбросить (Reset), Восстановить, Обзор,
  Загрузить (Load and Upload), Скачать (Download), Выбрать (Select, Choose;
  Выделить for text), Развернуть / Свернуть, Оставить (Keep), Отменить (Undo),
  Вернуть (Redo), Повторить (Retry), Откатить (Revert), Включить / Отключить,
  Вкл. / Выкл. (On / Off). Generate is Построить (построение, построить
  заново), never Сгенерировать. Session is сеанс, Settings is Настройки, one
  setting is параметр.
- **Kept in English**, lower case in running text and never declined: query,
  subject, coortholog, inparalog, near-reciprocal hit, reciprocal best hit,
  seed, locus tag, GC skew, AT skew, skew, upstream / downstream, flat file,
  lookbehind, and every name of STYLE 1.1. Declined loanwords: хит (hit),
  рендерер, хеш, кэш, репликон, конвейер.
- **Never use**: кликнуть, клик, фича, признак or особенность (Feature), юзер,
  конфиг, ресет, сессия, нить, аксессия, учетный / регистрационный /
  инвентарный номер, метка локуса, локус-тег, GC-асимметрия, ГЦ, рекорд,
  линейная / круговая запись (write запись линейной диаграммы; кольцевая
  запись is only the topology), capital Вы / Ваш, пожалуйста, any ты form.

### 16.2 Buttons and sentences

- Register: вы, lower case, implicit (Microsoft, Apple); never ты or capital Вы.
- Buttons, menu commands, links, button tooltips: the infinitive, no final
  period (Сохранить сеанс, Загрузить пример, Применить и построить).
- Labels, headings, tabs, column headers, badges: a nominative noun phrase
  without a final period (Пользовательские слоты треков). Options are nouns or
  adjectives (Авто, Последовательно). A dialog title that asks is the
  infinitive and a question mark (Убрать весь файл?).
- Instructions (help, notices, dialogs, errors): the plural imperative: Выберите
  запись. Нажмите кнопку Построить. Leave blank for auto: Оставьте пустым для
  автоматического значения.
- Descriptions: third person present or impersonal (Показывает глубину покрытия
  как трек. Применяется при построении.). An option or button name is not a
  subject: Если выбрано Авто, …, not Авто выбирает ….
- Progress (…ing): verbal noun + … (Сохранение сеанса…, Построение диаграммы…).
  Completion: a neuter short participle as a status (Применено. Сохранено.) or
  noun and participle (Сеанс сохранен.). A lone state label is neuter: Включено,
  Отключено, Недоступно, Не применено.
- Errors: what failed, then what to do: Не удалось прочитать файл. Убедитесь,
  что это файл FASTA. Cannot is не удалось (failed) or невозможно (not
  allowed). Retry is Повторить on a button, Повторите попытку in a sentence.
- Required / optional / recommended: Обязательно, необязательно (invariable:
  Область (необязательно)), Рекомендуется. Above 0 is положительное (больше 0);
  possibly 0 is неотрицательное.
- A control named in a sentence takes a class noun and its label: нажмите кнопку
  Построить. A name the English quotes is quoted in «…»; one it leaves bare
  stays bare, except a tab name, which is always quoted after вкладка: на
  вкладке «Легенда», вкладки «Подробности» и «Квалификаторы».

### 16.3 Punctuation and spacing

- Sentences end with a period; labels and button texts have none. A label
  ending with a colon keeps it. No space before a colon, semicolon, comma,
  question mark or exclamation mark (Метка: {name}).
- Quotes are guillemets «…» with no space inside. Never straight " or English
  “ ”. Do not nest quotes: rephrase.
- Hyphen only inside words (GC-состав). The em dash — has a space on each side;
  the en dash – marks a range (Строки 3–4). One … (U+2026), no space before it.
  e.g. is например, i.e. is то есть; abbreviations: мин., макс. (minimum,
  maximum).
- Digits are ASCII with the source's point and grouping (0.65, 1,000), never
  0,65 or 1 000, because the app formats numbers with en-US. Units and % keep
  the source spacing (100 bp, 50%, 20px). Sentences are joined by one space.
- Time abbreviations are Russian words, as in the ja, zh, ko, th, id and vi
  catalogs: {value} мс, {value} с, {minutes} мин {seconds} с (no period). The
  measurement units bp, kb, Mb, px, pt, ×R and DPI stay Latin (STYLE 1.2).
- No non-breaking or soft hyphen characters: a word is never split. The letter
  ё is never written, only е (все, учет, отсчет), as Microsoft Russian does.
  Separate Latin and Cyrillic words by a space or a hyphen (файл LOSAT).

### 16.4 Words and capitals

- Sentence case for every label, button, heading, tab, option and column
  header: the first word, proper names and acronyms only (Заголовки и подписи
  записей); a concept the English capitalizes (Session) is lower case in prose.
- A key the English writes in capitals stays in capitals: ЦВЕТА ПО УМОЛЧАНИЮ,
  ИНДИВИДУАЛЬНЫЕ ПРАВИЛА, ФИЛЬТРАЦИЯ ПОДПИСЕЙ, ПРИОРИТЕТ КВАЛИФИКАТОРОВ,
  ГЕОМЕТРИЯ ПОДПИСЕЙ, ДОБАВИТЬ КЛЮЧ; the badge OFF is ВЫКЛ. Cyrillic keeps the
  CSS `uppercase` and `tracking-*` rules.
- Compounds put the head noun first, genitive after it: тип элемента, TSV
  глубины. Gender: запись, цепь, подпись, глубина feminine; элемент, трек, слот,
  ген, хит, сеанс masculine; покрытие, сходство, выравнивание neuter; CDS
  feminine, TSV masculine. Exact match: Продукт (точное совпадение): {product}.

### 16.5 Placeholders and counts

- A placeholder's gender, number and case are unknown and never inflected. Put a
  class noun before it and let agreeing words follow that noun: запись {record},
  файл {file}, трек {track}. No ending or hyphenated suffix ({count}-й). A
  placeholder as subject takes a neuter or impersonal verb.
- `Intl.PluralRules('ru')`: one for 1, 21, 31, 101; few for 2–4, 22–24, 32–34;
  many for 0, 5–20 (11–14 included), 25–30, 100; other only for fractions (1.5).
  After the number: one nominative singular, few genitive singular, many
  genitive plural, other genitive singular: 0 записей, 1 запись, 1.5 записи,
  21 запись, 22 записи, 25 записей.
- The `one` key holds the text for exactly one (count 1); the `other` key holds an object with one entry per plural category of the other counts: `few`, `many` and `other` (fractions), and also `one` (the form of 21, 31, 101 …, which Intl classes as one); the runtime shows the object's category text for every count but 1.
  Prefer a construction with the noun in the genitive plural after a colon
  (Пропущено строк аннотации: {count}): all four forms are then the same text.
  Where agreement is needed: a verb is singular in one, plural in few, neuter
  singular in many; an adjective before a feminine noun is nominative plural in
  few, before a masculine or neuter noun genitive plural.
- When the English `one` drops {count} (From this color rule:), Russian `one`
  still covers 21, 31 and 101: use a count-neutral wording if one exists,
  else the text is right for 1 only (listed for PR-B in the notes).
- Worked pairs (ru one / few / many / other):
  1. ({count} record) / ({count} records): ({count} запись) / ({count} записи)
     / ({count} записей) / ({count} записи)
  2. {count} gene / {count} genes: {count} ген / {count} гена / {count} генов
     / {count} гена
  3. Skipped {count} annotation row. / rows.: Пропущена {count} строка
     аннотации. / Пропущены {count} строки аннотации. / Пропущено {count} строк
     аннотации. / Пропущено {count} строки аннотации.
  4. {current} of {count} match / matches: {current} из {count} совпадения /
     {current} из {count} совпадений / {current} из {count} совпадений /
     {current} из {count} совпадения
  5. Remove {count} unmatched feature edit / edits: Убрать {count}
     несопоставленную правку / … несопоставленные правки / … несопоставленных
     правок / … несопоставленной правки

### 16.6 Length

- Russian runs 20–30% longer than English. Help texts and notices wrap; a label
  may be about 1.3× the English label, a header button about 1.4×.
- Rail short names (52 px entries, 10 px semibold, about nine Cyrillic letters
  on one line): Ввод, Сравнение, Доп., Области, Макет, Шрифт, Заголовки, Легенда,
  Цвета, Элементы, Подписи, Ось, GC / skew, О программе (two words, wraps
  between them); Построить beside the rail. Fallbacks if a screenshot shows an
  overflow: Сравн. for Сравнение, Названия for Заголовки.
- Header buttons: Сбросить настройки, Сохранить сеанс, Загрузить сеанс, Загрузить
  пример, Построить диаграмму; Undo and Redo are icons with the tooltips
  Отменить: {action} and Вернуть: {action}. Check them at 1,024 px and fix an
  overflow in CSS, not with a shorter key, since the keys name them in prose.
- Option labels in 10 px rows use the glossary short forms: Авто, Вкл. / Выкл.,
  Снаружи / Внутри, Сверху / Снизу, Попарно, Последовательно, Многопоточно,
  Лента / Кривая. Never switch terms to shorten.

### 16.7 Machine checks

`tests/web/i18n-lint/ru.tsv` holds these rules; each matches one catalog value.

- Letters: `ru-ukrainian-letters` (і ї є ґ), `ru-yo` (ё), `ru-combining-marks`
  (stress accents), `ru-mixed-script` (warn), `ru-non-ascii-digits`.
- Punctuation: `ru-quotes-not-guillemets`, `ru-space-before-punct`,
  `ru-hyphen-as-dash`, `ru-three-dots`, `ru-decimal-comma` (warn).
- Placeholders and plurals: `ru-placeholder-suffix`, `ru-placeholder-prefix`,
  `ru-plural-parens` (no "(и)" plural; use a tn pair).
- Register: `ru-informal-ty`, `ru-vy-capital`, `ru-singular-imperative`.
- Words: `ru-never-slang`, `ru-session-word` (warn), `ru-feature-synonyms`
  (warn), `ru-accession-alt`, `ru-locus-tag-translated`,
  `ru-gc-skew-translated`, `ru-strand-nit` (warn), `ru-mode-record`,
  `ru-english-loanwords` (warn).

## 17. Ukrainian (`uk`)

Standard Ukrainian of Ukraine, written from the English source and the glossary
only (D-07: no Russian draft is read or reused). The Language menu entry is the
endonym **Українська**; the menu label itself is Мова. Where Ukrainian software
and science have a word, it is used even when the Russian cognate is shorter
(§17.1 lists the calques); where neither has one, the English word stays (D-12).

### 17.1 Vocabulary

- **Biology and bioinformatics**: геном, послідовність, ген, білок, ланцюг (ланцюг +, ланцюг −; прямий ланцюг;
  "ланцюжок" is not used), зворотний комплемент, анотація, запис (Record), ознака (Feature; тип ознаки, таблиця
  ознак), реєстраційний номер (Accession), локус, кваліфікатор (the keys /gene, /product stay), розташування
  (Location), продукт, організм, штам, вид, органела, вміст GC, динуклеотид, глибина (Depth), глибина секвенування
  (read depth in sentences), покриття (Coverage), ідентичність, поріг, гомологія, ортолог, паралог, колінеарний,
  блок, вирівнювання, генетичний код.
- **Figure words**: доріжка (Track), рівень (Lane), слот (track slot), позиція (File card slot), кільце (Ring),
  легенда, підпис (Label), заголовок, поділка (Tick; основна / проміжна поділка), позначка (annotation mark),
  лінійка (Ruler; "правило" is a color rule), масштабний відрізок (scale bar, Length bar), шкала (scale of
  coordinates), заливка, контур (Stroke of shapes), лінія (axis and intron lines), зразок кольору (Swatch),
  підкладка / передній план (Underlay / Foreground), відступ (Padding), проміжок (Gap), відстань (Spacing), інтервал
  (Interval), зсув (Offset), пресет (Preset), типовий (Default).
- **Command words**: Зберегти, Завантажити (Load), Звантажити (Download), Передати (Upload), Відкрити, Закрити, Копіювати, Вставити, Дублювати, Видалити (Delete), Вилучити
  (Remove), Очистити (Clear), Скинути (Reset), Відновити (Restore), Застосувати, Скасувати (Cancel and Undo),
  Повторити (Redo), Спробувати ще раз (Retry), Пропустити (Skip), Залишити (Keep), Перейменувати, Додати, Вибрати
  (Select), Обрати (Choose), Показати / Приховати, Увімкнути / Вимкнути (Enable / Disable), Увімкнено / Вимкнено (On
  / Off), Розгорнути / Згорнути, Імпортувати / Експортувати, Створити (Generate), Створити знову (Regenerate),
  Натиснути (click, tap), Перетягнути, Огляд (Browse).
- **Load, Download and Upload** need three verbs, because Load Session, Download and Upload appear on the same
  screens. The sources disagree: Microsoft writes Завантажити for Download and Передати for Upload; Google writes
  Завантажити for Upload; Mozilla writes Завантажити for Download; KDE and Debian write Звантажити for Download.
  The split Завантажити (Load) / Звантажити (Download) / Передати (Upload) is a decision of this catalog on weak
  evidence, listed for the Owner; no vendor uses all three this way.
- **One word, one concept**: ознака is Feature only (a list item is елемент, a thing on the canvas is об’єкт); запис
  is Record only, never a log line; запит is Request, and Query stays query; збіг is Match; hit stays hit; таблиця
  is a UI table, таблиця кольорів the color-table file; доріжка is Track, смуга is Band, рівень is Lane, слот is a
  track slot, позиція is a File card slot or Position; вираз is a regex Pattern, шаблон a Template.
- **Kept in English** (no established Ukrainian term found, D-12; lower case, indeclinable): GC skew and skew, hit,
  query, subject, locus tag, orthogroup, coortholog, inparalog, reciprocal best hit, near-reciprocal hit, seed,
  upstream / downstream, cross-origin isolated, lookbehind, params, hash; plus the names of §1.1. The Ukrainian
  loanwords пресет, слот, escape-послідовність and хеш (kept as hash in the hash= syntax) are accepted.
- **Never use** (Russian calques where a Ukrainian term exists; also in the lint list): сесія / сессия (Session is
  сеанс), відмінити (скасувати), по замовчуванню and дефолт (типовий, типово), включити / відключити / виключити
  (увімкнути / вимкнути), слідуючий (наступний), приймати участь (брати участь), на протязі (протягом), настройки /
  настроюваний (налаштування, власний), колонка (стовпець), трек (доріжка), фіча (ознака), закладка (вкладка),
  кликнути (натиснути), вірний (правильний), любий (будь-який), додаток for an app (застосунок), "при" for "when"
  (під час, після, у разі), the letters ы э ъ ё.

### 17.2 Buttons and sentences

- Register: ви, written in lower case; never ти. Buttons, menu commands, links and button tooltips: the infinitive
  (Зберегти сеанс, Завантажити приклад, Скасувати), as Microsoft, Google and Apple Ukrainian do; no article exists.
  Checkboxes follow the English form (Показати глибину).
- Labels, headings, tabs, column headers, options, badges: a noun phrase in the nominative without a final period.
  An option agrees with its setting (Авто, Увімкнено, Послідовно, Вище / Нижче). A dialog title that asks is the
  infinitive plus ?: Вилучити весь файл? Скинути всі налаштування?
- Instructions (help, notices, dialogs, errors): the polite plural imperative (Виберіть запис. Натисніть «Створити».
  Створіть діаграму знову.). A help text that describes its control is in the third person present (Показує вісь
  глибини.).
- Work in progress: a verbal noun and … (Збереження сеансу…, Створення діаграми…). Completion: an impersonal -но /
  -то form (Сеанс збережено. Вирівнювання застосовано: …). Prefer this form for counted sentences: it does not agree
  with the number (Прочитано {count} вихідний запис).
- Errors: what failed, then what to do. "Could not X" is Не вдалося X (Не вдалося імпортувати анотації: {reason});
  "cannot" in a rule is не можна or неможливо. "X failed" is Помилка: X or X не вдалося. Busy reason: Створення
  діаграми. Повторіть після завершення.
- Required, optional, recommended: Обов’язково, (необов’язково), Рекомендовано. Auto is Авто on a button and in 10
  px rows, автоматично in a sentence; Default is типовий (Авто (типово)).
- An option, button or setting name is quoted when it is cited as a word: «Авто» приховує…, Натисніть «Створити»,
  not Авто приховує….
- Circular and Linear name the diagram: for a mode noun use the feminine (Кільцева, Лінійна, to agree with
  діаграма); as a modifier of a record or scale say кільцевої / лінійної діаграми (запис кільцевої діаграми,
  масштабний відрізок лінійної діаграми); кільцевий запис is only the topology.
- An exact-match condition follows the field in parentheses: Продукт (точний збіг): {product}. A count of members is
  "Кількість …": Кількість генів у групі.

### 17.3 Punctuation and spacing

- Quotes are « » with „ “ nested; a {name} the English quotes becomes «{name}»; one the English leaves bare stays
  bare. No straight " and no “ ” in prose. A period or comma follows the closing quote.
- Apostrophe: the typographic ’ (U+2019) in every word (об’єкт, необов’язково, ім’я, пов’язаний, Об’єднати, п’ять),
  never ASCII ' or ʼ. Letters ґ є і ї are used as spelled; і / й, у / в and з / із / зі follow euphony in running
  text, but before a placeholder the sentence is rephrased (see 17.5).
- Colon: no space before, one after (Поле: {path}); a label ending with a colon keeps it. No space before ; ? !
  either.
- Dash: the em dash — with U+00A0 before it and an ordinary space after; a hyphen is not a dash. The en dash of
  ranges stays (Рядки 3–4). The minus sign − and + in ланцюг + / ланцюг − stay.
- Digits are ASCII with the English decimal point and grouping (0.65, 1,000), never 0,65 or 1 000 (D-05). Units and
  % keep the source spacing (100 bp, 50%, 8 pt, 0.65 ×R); one … (U+2026), no space before it.
- Time abbreviations are Ukrainian words, as in the ja, zh, ko, th, id and vi catalogs: {value} мс, {value} с,
  {minutes} хв {seconds} с, without a period. The measurement units bp, kb, Mb, px, pt, ×R and DPI stay Latin
  (STYLE 1.2).
- Parentheses and slashes as in the source: Ділянка (необов’язково), Організм / штам. "from … to" is від … до (від
  #{from} до #{to}).
- The code joins sentences with one ordinary space (D-09). e.g. is наприклад, i.e. is тобто; abbreviations take a
  period (макс., мін., Порівн.), except the time units above. ID is invariable.

### 17.4 Words and capitals

- Sentence case for every label, button, heading, tab, option and column header: the first word and proper names
  only (Власні слоти доріжок, Попарні порівняння, Вісь і шкала). Never copy English Title Case. Cyrillic keeps the
  CSS `uppercase` rule, so a sentence-case key that CSS capitalizes stays in sentence case; a key the English writes
  in capitals stays in capitals: ТИПОВІ КОЛЬОРИ, ОКРЕМІ ПРАВИЛА, ФІЛЬТРАЦІЯ ПІДПИСІВ,
  ПРІОРИТЕТ КВАЛІФІКАТОРІВ, ДОДАТИ КЛЮЧ; the badge OFF is ВИМКНЕНО.
- Concepts the English capitalizes (Session, Result, Legend) are lower case in a sentence (сеанс, результат,
  легенда); gbdraw stays lower case. A control or section cited in a sentence keeps its label with a capital in « »:
  натисніть «Створити», у розділі «Власні слоти доріжок». The action Generate in running text is створення, створити
  (під час створення, перед створенням).
- Compounds put the head noun first: тип ознаки, TSV глибини, кеш LOSAT, підпис запису, ID запису. Genders:
  діаграма, доріжка, ознака (f.), запис, файл, слот, результат (m.), розміщення (n.), CDS (f., like послідовність),
  TSV and BLAST (m.). The adjective after a Latin name agrees with the Ukrainian head noun: таблиця BLAST, формат
  GenBank.
- An English word kept in the text (hit, query, GC skew) is lower case, never declined, and has no hyphenated
  suffix: кількість hit, не hit-ів.

### 17.5 Placeholders and counts

- A placeholder's gender, number and first sound are unknown. Put a class noun before it and let prepositions, adjectives
  and case endings agree with that noun: запис {record}, файл «{file}», доріжка {track}, налаштування
  {setting}. Never decline or suffix a placeholder ({name}-ом, {n}-й), never put у / в, з / із / зі or і / й
  directly against one: use від … до, "із" before numerals ({current} із {count}), "та" instead of і / й before a
  placeholder.
- {count} precedes the noun; the noun and its adjectives follow the number by `Intl.PluralRules('uk')`: `one` for 1,
  21, 31, 101 (not 11); `few` for 2–4, 22–24, 32–34 (not 12–14); `many` for 0, 5–20 (so 11–14 too), 25–30; `other`
  for fractions such as 1.5. The `one` key holds the text for exactly one (count 1); the `other` key holds an object with one entry per plural category of the other counts: `few`, `many` and `other` (fractions), and also `one` (the form of 21, 31, 101 …, which Intl classes as one); the runtime shows the object's category text for every count but 1. After a numeral the noun is nominative singular (one), nominative plural (few),
  genitive plural (many) and genitive singular (other); a neuter noun such as розміщення has the same form for one
  and few. Keep {count} in a form wherever the English one form has it; "From this color rule:" has none.
- Prefer the impersonal -но form or a noun phrase over a verb that must agree with the number. Where a verb is
  needed use a form that does not change (Прочитано, Вилучено), or put the count after a colon (Обов’язкові стовпці:
  {count}).
- Worked pairs, `one` / `few` / `many` / `other` (English `other` key in parentheses):
  1. ({count} source records inspected): Прочитано {count} вихідний запис / Прочитано {count} вихідні записи /
     Прочитано {count} вихідних записів / Прочитано {count} вихідного запису
  2. ({count} LOSAT results are cached and reused …): У кеші {count} результат LOSAT; його використовують повторно,
     якщо вхідні дані й налаштування збігаються. / … {count} результати LOSAT; їх використовують … / … {count}
     результатів LOSAT; їх використовують … / … {count} результату LOSAT; його використовують …
  3. (Reset {count} placements to Auto): Скинути {count} розміщення до значення Авто / Скинути {count} розміщення до
     значення Авто / Скинути {count} розміщень до значення Авто / Скинути {count} розміщення до значення Авто
  4. ({current} of {count} matches): {current} із {count} збігу / {current} із {count} збігів / {current} із {count}
     збігів / {current} із {count} збігу
  5. (From these {count} color rules:): the `one` key From this color rule: is З цього правила кольору: (count 1
     only); the object is З {count} правила кольору: / З цих {count} правил кольору: / З цих {count} правил
     кольору: / З {count} правила кольору: (the catalog's forms)
- Few / many of the counted nouns (also in the glossary notes): запис записи / записів, ознака ознаки / ознак, ген гени /
  генів, пара пари / пар, рядок рядки / рядків, збіг збіги / збігів, доріжка доріжки / доріжок, результат результати /
  результатів, редагування редагування / редагувань.
- Rule check: 0, 5, 12 and 25 take `many` (0 записів, 12 записів, 25 записів); 21 takes the `one` entry of the `other` object (21 запис); 22 takes
  `few` (22 записи); 1.5 takes `other` (1.5 запису). A `one` string never says "один", since 21 and 101 use it.

### 17.6 Length

- Ukrainian runs 20–30 % longer than English, and Cyrillic letters at 10 px semibold are about 10 % wider than Latin
  ones (about eight letters on one line in the 52 px rail entry). Shorter wording first (D-11): drop a noun, then
  use the short form in the glossary; never hyphenate or split a word.
- The 14 rail names (rail keys where the section title is longer): Вхід, Порівн., Додатково, Ділянки, Макет, Шрифт,
  Заголовки, Легенда, Кольори, Ознаки, Підписи, Вісь, GC / skew, Про програму (two lines, between words); Створити
  beside the rail. Додатково and Заголовки (9 letters) are the widest words; a word wider than 58 px needs a shorter
  rail key chosen at the fit check.
- Header buttons (Скинути налаштування, Зберегти сеанс, Завантажити сеанс, Завантажити приклад) are 1.2–1.6 times
  the English width: check them at 1,024 px. If one overflows, the fit check decides a shorter glossary form
  (candidates: Скинути, Відкрити сеанс, Приклад), not a syllable break.
- Options in 10 px rows use the glossary's short forms: Авто, Послідовно, Багатопотоково, Вище / Нижче, Зовні /
  Всередині, На осі, Підкладка, Стрічка / Крива. Help texts and notices wrap and may be as long as needed.

### 17.7 Machine checks (`tests/web/i18n-lint/uk.tsv`)

- `uk-russian-letters`: a string with ы э ъ ё. `uk-ascii-apostrophe`, `uk-modifier-apostrophe`: an apostrophe that
  is not ’. `uk-straight-quote`, `uk-wrong-quotes`, `uk-guillemet-space`: quotes other than « » or a space inside
  them.
- `uk-space-before-punct`, `uk-three-dots`, `uk-spaced-hyphen`: a space before : ; ? !, three dots, a hyphen used as
  a dash. `uk-nonascii-digits`: digits other than 0-9. `uk-placeholder-suffix`: Cyrillic letters glued to a
  placeholder. `uk-mixed-script`: Latin and Cyrillic letters touching.
- `uk-register-vy-capital`, `uk-register-ty`: capitalised Ви / Ваш mid-sentence, and the ти form or a singular
  imperative (натисни, вибери).
- Russian calques (17.1): `uk-calque-session`, `uk-calque-undo`, `uk-calque-default` (warn variant
  `uk-default-variant`), `uk-calque-enable`, `uk-calque-participle`, `uk-calque-settings`, `uk-calque-pry` (warn),
  `uk-calque-column`, `uk-track-trek-fich`, `uk-calque-words` (warn).
