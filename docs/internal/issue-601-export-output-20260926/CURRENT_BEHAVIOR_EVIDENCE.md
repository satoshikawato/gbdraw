# Issue 601 — 調査結果と検証の範囲

調査基準: `origin/dev` / `d457b7189b137185a8dec800819a312c30b969fa`。
確認日: `2026-09-26`。以下は修正前sourceの証拠であり、修正実装の合格証拠ではない。

## Issueの内容と現行実装

[Issue #601](https://github.com/satoshikawato/gbdraw/issues/601)はPDFの文字対応、エラーの使いやすさ、Python regexの拒否を報告している。

| 報告 | 確認結果 | 実装への含意 |
| --- | --- | --- |
| BUG-07 | `pdf-fonts.js`はLiberationにないglyphでPDF全体を拒否する。非ASCII全般の拒否ではない | glyph検査を維持して対応fontで不足部分を補う |
| BUG-15 | normalizer単体で内部ValueErrorをsummaryへ直接公開し、stderrのtraceback除去では末尾causeまで消すことを再現。Alignにはerrorなしreturnによるgeneric fallbackもある | producerの識別情報をtransportで保ち、一つのnormalizerで表示する |
| BUG-19 | Color/Labelは[PR #538](https://github.com/satoshikawato/gbdraw/pull/538)でPython評価へ統一済み。Feature SearchはJavaScript regex | 既存評価を保護し、入力欄ごとにsyntaxを案内する |

Issueに記載された`GUI_AUDIT_DEV_20260926.md`は基準treeに存在しない。報告者が使った入力欄・pattern・buildは未特定であり、BUG-19をすべてのregex画面の回帰と断定しない。
比較identity不一致の例外発生源は`gbdraw/render/groups/linear/pairwise_match.py`。そのraw例外文によるnormalizer再現は行ったが、実Genome入力から不一致を誘発する再現は未実施。

## 実行済みの既存検査

```bash
node --test tests/web/error-normalization.test.mjs tests/web/rule-matching.test.mjs
pytest tests/test_web_rule_matching.py -q
GBDRAW_WEB_TEST_PORT=4189 npx playwright test tests/web/python-rule-parity.playwright.spec.js tests/web/gui-audit-regressions.playwright.spec.js --grep 'Python-only color|Python label TSV|PDF preserves Greek' --workers=1 --retries=0
```

| 検査 | 結果 | 範囲 |
| --- | --- | --- |
| Node | 8 passed | 既存normalizer、Python-only pattern、catalog/stale/priority |
| Python | 9 passed | Python ruleとnativeの一致、invalid syntax、empty catalog |
| Chromium | 3 passed / 24.0 s | Greek/math PDF text、Python-only ColorとLabelのatomic commit/History/Generate |

上記のsourceとinputsが変われば関連検査を再実行する。全CI、Firefox/WebKit、deployの検証結果ではない。

## font coverageとPDF probe

fontToolsの`TTFont(...).getBestCmap()`によるRegularの確認結果:

| font / text | 結果 |
| --- | --- |
| Liberation Sans/Serif/Mono / `β-lactamase α ≥ 95%` | 不足なし |
| Liberation / `プラスミドA` | 日本語glyph不足 |
| 現行Noto Sans JP Japanese / `日本語` | 不足なし |
| 現行Noto Sans JP Japanese / `简体中文` | U+7B80不足 |
| 現行Noto Sans JP Japanese / `繁體中文` | coverage不足なし、地域字形は未検証 |
| 上記font / `𠮷` U+20BB7 | 不足 |

現行Result/export経路へ内部fixture `プラスミドA`を入れると、PDFはU+30D7のerror、download 0件。同snapshotのSVG（172 bytes）とPNG（16,393 bytes）は取得できた。PNGの字形を全端末で保証する証拠ではない。

既存Noto Japanese Regular/Bold WOFF2を一時TTFに変換し、同じjsPDF/svg2pdfへ手動登録して500×160 ptの内部fixtureを出力した。Regular 2,363,384 bytes、Bold 2,363,108 bytes、計4,726,492 bytes。PDFは284,415 bytes、FontFile2あり、image XObjects 0件。既存`tests/web/helpers/pdf-text.cjs`による抽出は`プラスミドA 日本語日本語 太字プラスミドA β ≥ 95%`。

これは手動run指定による日本語text PDFの実現性証拠。自動fallback、Latinの自動維持、製品経路、font-matched metrics、独立renderer、italic、textPath、補助平面、最大図は未検証。中国語や日中font全体の容量を保証しない。一時probeは永続的な再現recipeになっていないため、S00はsourceとcommandsを保存して再現可能な証拠を取得する。

vendored jsPDFは3.0.3。UTF-8 pathに`characterToGlyph(t.charCodeAt(s))`と`toUnicode[n] = t.charCodeAt(s)`を使う箇所がある。supplementary/結合・異体字をcmap coverageだけで保証しない。

## authorityの基準

- `OPTION_INTEGRITY_PRODUCT_CONTRACT.md`: revision 16。C03/C04/C06/C07、PD-OI-016、PD-OI-031 revision 4、PD-OI-034 revision 4の既存failure isolationとcause/retry要件を確認。
- `tools/web-product-impact-map.json`: canonical-render-request-boundary / saved-session-regeneration-continuity。PDF fallbackとdiagnostic disclosure concernは未登録。
- `tools/web-product-decisions.json`: active decisionsは空。存在しないBD番号を引用しない。
- [承認内容](./APPROVED_DECISIONS.md)は製品判断を確定する。実現性証拠と独立authorityのbase統合は[総合計画書](./MASTER_PLAN.md)の手順で扱う。

## 一次資料

- [jsPDFのUTF-8 / custom TTF](https://github.com/parallax/jsPDF#use-of-unicode-characters--utf-8)
- [svg2pdfのcustom fonts](https://github.com/yWorks/svg2pdf.js#concerning-custom-fonts-and-non-us-ascii-characters)
- [NotoのfallbackとCJK地域font](https://notofonts.github.io/noto-docs/website/use/)
- [文字clusterとglyphの対応](https://harfbuzz.github.io/clusters.html)
