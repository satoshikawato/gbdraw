# Issue #600 — 変更前の再現根拠

調査・実行: 2026-09-26。
Source: `origin/dev@d457b7189b137185a8dec800819a312c30b969fa`。
`git fetch origin` 後も同 SHA。`origin/main` は
`4556e04e929a4a85ad28d1833ce7304bd764881c`。
この baseline 記録の作成時点では runtime code / existing tests / authority / references は未変更。

## Source availability and authority

[Issue #600](https://github.com/satoshikawato/gbdraw/issues/600) の本文とコメントを
GitHub connector の読取で取得した（取得時 comments 0）。Issue が挙げる
`docs/internal/GUI_AUDIT_DEV_20260926.md` はこの checkout にない。未読の監査を
根拠として補完しない。source/runtime の既知挙動は以下の probe で独立に再現した。

`tools/web-product-decisions.json` の decision 配列は空。
`tools/web-product-impact-map.json` の concern は canonical request と Result admission。
`OPTION_INTEGRITY_PRODUCT_CONTRACT.md` には本４ outcome の decision がない。
`PD-OI-016` の last-successful-Result isolation は横断して維持する。
基準と `origin/main` の CLI Reference は unknown columns の拒否、Circular gap の
without-a-unit を明記。これらは S00 の authority-only 改訂で対象 scope に限定して更新する。
main 第一親の関連 path history も確認した。新 compatibility reader は本案にはない。

## Regenerate observations

repo root から実行する。production の読み取りと、使い捨て record/table/ESM copy のみ。
JSON はこの evidence ディレクトリへ書く。

```bash
PYTHONPATH=. python docs/internal/issue-600-implementation-20260926/evidence/reproduce.py > docs/internal/issue-600-implementation-20260926/evidence/python-observations.json
node docs/internal/issue-600-implementation-20260926/evidence/reproduce.mjs > docs/internal/issue-600-implementation-20260926/evidence/javascript-observations.json
```

- [Python probe](evidence/reproduce.py) / [raw output](evidence/python-observations.json)
- [JavaScript probe](evidence/reproduce.mjs) / [raw output](evidence/javascript-observations.json)

| Case | Observed result | Interpretation / limit |
| --- | --- | --- |
| BUG-08 notes column | JS/Python とも unknown column error | Auxiliary header 単独で拒否。browser file upload は未実行 |
| BUG-09 absent / present+absent | どちらも ValidationError | 未一致１件で row geometry を返せない。完全一致 control は１ annotation、warning なし |
| BUG-10 two used colors | Web import error、Python 凡例は最後の #445566 のみ | Web の拒否を外すだけでは色を説明できない。Python は legend function probe で、full SVG generation は未実施 |
| BUG-14 10 | raw validation 成功、payload 10 | control |
| BUG-14 10px | raw validation 拒否、payload 10 | validator と projection が不同意 |
| BUG-14 10PX / oops | raw validation 拒否、payload null | payload 経路で不正値が auto に見える |
| BUG-14 px | raw validation 拒否、payload 0 | 単位だけの文字列が zero になる |
| BUG-14 linear 10px | raw validation 成功 | Circular と異なる |
| BUG-14 Circular CLI 10px | parse error: without a unit | 正しい slot 構文 `gc:dinucleotide_content@inner_gap_px=10px` で確認 |

payload-only の probe は validator を意図的に呼ばない直接関数テスト。
通常の Generate が不正値を黙って通すと主張する証拠ではなく、normalizer/
serializer を単独利用した場合の差を示す。browser journey の未検証と分けて扱う。

## Existing focused checks

```bash
pytest tests/test_annotations.py tests/test_annotation_planning.py -q
node --test tests/web/annotations.test.mjs tests/web/file-imports.test.mjs tests/web/track-slot-validation.test.mjs
```

- Python: **15 passed**, 2.55s。Python 3.13.3、pytest 9.0.2。
- Node: **30 passed**, failures/skips 0、869.546049ms。
- これらは strict な現状を含む既存 tests の基準確認。承認済み新 outcome が
  既に実装されて通ったという意味ではない。
- Playwright CLI と Python package、Node `@playwright/test` を確認済み。
  browser unavailable ではない。baseline は code review + function probe の範囲で、
  browser/画面確認は実施していない。runtime 実装後の acceptance に残す。
- full suite、Ruff production、output comparison、browser wheel build は今回未実施。
  この baseline から新 outcome の acceptance 成功を推論しない。

## Evidence checksums

観察結果は生成時点の値。将来コードが変われば同じ出力は保証しない。

| File | SHA-256 |
| --- | --- |
| `reproduce.py` | `5dd869df1b5e9a9f9b90e16039f8bfffbc9e0779783252fcec6cee7889371b29` |
| `reproduce.mjs` | `175c16e9d2eeb840b1b87bf61cf35ab6bcecbaee0beb9976a2d2f6ee390532b7` |
| `python-observations.json` | `001f3fee44ab7885d0673d412b18b0254a5393f770177fbfceb5a8aff68d328b` |
| `javascript-observations.json` | `a557a87906552c1e1030872fbc1a56de3ad81f4b4f052f0b4782cbdbc0f8fa95` |
