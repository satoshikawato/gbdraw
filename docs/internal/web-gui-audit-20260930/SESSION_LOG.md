# 実行状況

計画: [03_IMPLEMENTATION_REFERENCE.md](03_IMPLEMENTATION_REFERENCE.md)。判断: [02_DECISION_PACK.md](02_DECISION_PACK.md)。

| 日時 | PR | merge SHA | 直した ID | テスト | 生成物 | 残った問題 |
|---|---|---|---|---|---|---|
| 2026-09-30 | — | — | — | — | — | 監査（README）、修正提案（01）、Decision Pack（02。Owner が確定）、実装リファレンス（03）を作成した。runtime、Product Contract、CI は未変更。基準は origin/dev 4c89bab1 |
| 2026-09-30 17:28 JST | #649 (P00) | e97d90fe | — | docs tests 21 passed | CHANGELOG `[Unreleased]` に PR ごとの slot（P03〜P20）を追加 | なし |
| 2026-09-30 17:47 JST | #650 (P01) | a0cb65eb | — (authority) | 02 の 30 receipt が verbatim、SHA-256 と JSON が一致することを script で照合。architecture-contracts 0 fail | Product Contract revision 29: PD-OI-056〜084、PD-OI-018 revision 4、OIC-027 | N-15 は Product Contract と同じ PR に入れられない（isolation 規則）ため、別の docs PR にした |
| 2026-09-30 18:14 JST | #651 (N-15) | bfa18b8a | N-15 | ci 229/229、docs contracts 21 passed。コード側（smoke 上限 19、ratchet の MAXIMUM_RULE_COUNT=4）が正しく、文書を合わせた | なし | 監査に無かったずれ（gallery parity command は web-pr-smoke ではなく gallery job で走る）も直した |
| 2026-09-30 21:45 JST | #653 (P05) | 7545b074 | FE-07、N-05（ER FE-07: 旧い翻訳の結果は protein cache から昇格しない。cache key が protein ごとの aaSha256 を含むため。修正不要） | 修正前に 13 fail + ImportError。GTG/TTG→M、GFF3 phase、MG1655 の /translation と一致。dev staging（Tests、Gallery publication）緑 | なし（Gallery と参照出力は変わらない） | Owner が GitHub から直接マージした（レビュー担当は上限で途中終了。orchestrator が差分を事後に確認した）。直後に Owner の操作で dev→main の PR #654 が開かれたが、main は対象外なので触れていない |
| 2026-10-01 01:28 JST | #656 (P17) | 95228d4f | FE-06（D-01）、FE-11 と N-11（D-16）、FE-08、FE-12、PV-05（D-17）、PV-06 | 全 ID で修正前に失敗（dev の production 14 files で wheel を作って確認）。node 1126、pytest、Playwright 10 spec 92 passed | Gallery の interactive SVG 10 件、examples.json、artifact-manifest.json、docs/images の h-cli-13・t-py-08・t-cli-11 の interactive SVG を owner tool で作り直した。tobacco-chloroplast の検索と popup を Chromium で目視確認 | Owner-delegated 2 件（read_color_table は None/NA/null を文字として読む、match popup も共通の location formatter を使う）。残: app-setup.js の到達しない envelope fallback、t-cli-11 の session が dev の時点で schema 7 のまま（テストなし）。read_color_table の diagnostic= 化は P07 へ |
| 2026-10-01 02:46 JST | #655 (P02) | abf25c8f | N-14（検査基盤: G-G(3) の helper 検査、G-C、G-D の record 検出ベクタ、2 record の batch fixture、architecture-contracts の上限化） | Wave 2 の 50 件余りに known-defect の印（test.fail、assertKnownDefect、strict xfail、vector）。新 spec 48 件が dev 7545b074 と 95228d4f で想定どおり。既存 33 spec 202 tests が helper 変更後も通る | なし | tests だけの PR。staging が P17 で赤の間だが runtime を変えず原因と無関係なので先にマージした。印のない ID: N-09、N-10、N-18、N-19、N-20（再現なし）、N-17（P10 が fixture を足す） |
| 2026-10-01 12:31 JST | #657 (P06) | ceeb38d3 | X-01（Python 側）、X-02 PR-A（D-26 ACGTU、N-13）、D-27、SE-08(a)、GE-03（Python 側） | 修正前に失敗（共有ベクタ option_domain_vectors.json、test_dinucleotide_domain、test_web_error_diagnostics、test_label_filtering_derived_keys、test_circular_definition_interval）。別のレビュー担当が 191 の config 葉に -1/0/-0.5/NaN を与えて、拒否に変わったのが font size・stroke width・tick_width の 33 件だけと確認。Gallery 4 Session と main の v42 sidecar 7 件が dev と同じ SVG に replay | なし（参照出力と Gallery は不変）。mode-profiles.generated.js は owner tool で再生成 | Owner-delegated: GE-03 の間隔の導出を codec ではなく override の適用時に置いた（0.13.0 の Python API の挙動に戻る。W3 の注記とは異なる）、tick_width を stroke 幅として扱う。レビューで Linear 側に残っていた slot の dinucleotide 解決関数と ratchet の抜け道（diagnostic=None）を直した。P07 へ: Web 側の field・reason・context key の定義と、option_domain_vectors.json の Web 側 |

## 進め方と Owner-delegated の記録

表の行は PR をマージするたびに追記する。orchestrator が SESSION_LOG だけを変える docs PR でまとめて dev に入れる（並列の PR が同じ行を触って衝突しないようにするため）。

- 2026-09-30（進め方）: CHANGELOG の `[Unreleased]` に PR ごとの slot を置き、各 PR は自分の slot だけを置き換える。dev の branch protection が `strict`（最新への追従が必須）なので、マージは直列になる。
- 2026-09-30（進め方）: P02 と Wave 1（P03〜P06、P17、P18）を同時に実装する。P02 の known-defect の印は Wave 2 以降（P07〜P16）の ID に付け、Wave 1 の PR は自分の ID の失敗するテストを自分で書く（03 §4.1 の「失敗するテストを先に書く」は各 PR で満たす）。
- 2026-09-30（進め方）: N-15 は P01 から外し、別の docs PR にした。Product Contract の変更は他のパスと同じ PR に入れられない（`tools/check-web-change-budget.mjs` の isolation 規則）。
- 2026-10-01（進め方）: 使用量の上限で Workflow のエージェントが 2 回まとめて止まった（文脈を失い、再開のたびに読み直しが要る）。以後は Agent を最大 4 本並列に動かし、止まったら同じエージェントに続きを送って再開する。実装担当が自己レビューと CI まで持ち、科学的出力を変える PR（P03、P04、P06）だけ別のレビュー担当を付ける。
- 2026-10-01（進め方）: P09 を P07 より先に進める。03 の前提 P07 はホットなファイル（services/config.js）の共有によるもので、P07 はまだ始まっていないため。P09 が先にマージされれば、P07 がその上に rebase する。GE-06 の拒否の理由は既存の notice の経路で出し、P07 が diagnosticError の code に移す。
- 2026-10-01 02:45（staging）: #656（P17）のマージ後の dev staging（Tests run 36744424557）が赤。feature-fill-scope と active-result-edit-transaction（dev staging だけで走る functional Playwright）が feature popup の周りで失敗し、performance は cancelled。規則どおりマージを止め、P17 の担当が dev から新しいブランチで修正する。
- 2026-10-01（Owner-delegated、P04 / D-24）: 「definition_font_size を明示した場合は折り返さない」を「既定値 18 以外のとき明示とみなす」と近似した。明示的に 18 を指定した Web・Python の要求では折り返しが起きる（receipt の Must preserve との差）。正確に区別するには CLI の Session が保存する完全な config の形式を変える必要があり、D-35（CLI Session は現状維持）と R4 の Web↔CLI 一致に反する。D-24 の目的（Web 既定で長い学名の図を作れる）を保つため近似を採る。最終報告で Owner の確認事項にする。
- 2026-10-01 12:31（staging）: P17 の 3 件で staging が赤のまま P06 をマージした。原因は P17 の feature popup に特定済みで修正中、P06 は Python だけの変更で無関係、かつ P07 の前提のため。以後のマージでも、staging の失敗がこの 3 件だけであることを確かめる。
- 2026-10-01 12:50（staging）: P17 後の staging の失敗の原因は、FE-11 で Features drawer の行の title が 1 始まりの location になったのに、dev staging だけで走る 3 つの spec が古い 0 始まりの start..end で行を探していたこと（runtime の不具合ではない）。テストだけを直す PR #660 を出した。performance の cancelled は依存の install に 12 分かかった runner の遅れで、無関係。
- 2026-10-01 12:50（P03）: #658 を P06 の上に追従させると、P06 の producer coverage の ratchet が gbdraw/io/comparisons.py の新しい ValidationError 8 か所（diagnostic= なし）を拒否した。次のセッションで diagnostic= を付ける（Web で表の読み込みエラーが汎用表示になる P03 の残りも同時に直る）。
