# S08 follow-up — 重複計算の修正候補

## 2026-09-17 RW-04完了 — aliasの二重正規化

commit `6d361253` の後続を局所調査し、`validateProteinIdentityManifest` 内で
同じprimitive stringの`displayAlias.normalize('NFC').trim()`を、空値判定と
直後のalias集約用変数のために2回実行していると確認した。間にawait・内容変更はない。
既存の変数で一度計算し、その値の空値判定を残した。正常なP featureを訪れる
manifest検証1回につき2P→P、fixtureでは4→2。manifest検証回数は変えない。
Unicode・ordinal・不正入力・失敗順、RW-01〜03とSessionの既存回帰、実Workerでの
保存raw・失敗隔離・retryを確認した。[調査・証拠・handoff](S08_FOLLOWUP_ALIAS_VALIDATION.md)。
秒数改善は未測定。別境界の検証は維持し、次候補は未登録。この1件で今回の調査を終了する。

## 2026-09-17 RW-03完了

[merge検証の修正・証拠・handoff](S08_FOLLOWUP_MANIFEST_MERGE.md)。Generateの入力検証を
merge owner内へ集約し、正常時2R→R、必要な統合後検証1回を維持した。
不正入力・conflictの優先順、reload案内、公開helperの既定動作を回帰確認した。
source/installed packageの実Worker・desktop/mobileで保存raw検索0回、結果/provenance一致、
不正manifest時のResult/History/cache保持と再試行を確認。RW-01/RW-02は維持する。
秒数改善は未測定。近傍に同じ入力・寿命・結果と証明できる新しい重複は見つからなかった。

## 2026-09-17後続更新

RW-01→RW-02を実装し、getterと処理区間の回帰検証を行った。
[修正・検証・handoff](S08_FOLLOWUP_RAW_VALIDATION.md)を参照。
RW-01は後段の重複検証のみ削除、RW-02はGenerateでdeep copyした非共有manifestを
prepared-pair loop内だけ索引化し、成功・例外・Cancelでreleaseする。
Session reader/writerの検証、entry固有TSV検証と公開helperを維持した。
H件の正常hit・R recordsの区間はmanifest検証2H→1、runtime ID Set 4H→R、
TSV走査2H→H。秒数の改善を示すものではない。

引き継ぎsourceはcommit `8bf4b98ae02ef14f72caaebe3114f7a6fc3e0c75`。
以下の「02ca8f95＋未コミット差分」はcommit前の観測記録として残す。
S08は完了済み、S07承認・S06/S08 writerの将来merge前Review・S05却下を維持する。

## 初版の調査記録

2026-09-17。ユーザーの「明らかに無駄な計算とかあったらそれをドキュメントして。それを次なおそうや片っ端から」という追加依頼に対する記録。
今回は読み取りと文書化のみ。次回は下記の順で扱い、実装前に最終sourceとの差を確認する。
S07.5〜S07.8の追加削減は終了のまま、S07の承認・S06のmerge前Review・S05却下を維持する。
この一覧をS08の未完了条件にせず、cache/runtime再設計を承認したものとも扱わない。

## Sourceと証拠

専用branch `verify/collinear-s08-integration-20260916` のHEAD
`02ca8f950a7561ae9569c66f3f91d30702b53a17` **と未コミット差分**が対象。
全体は[S08 handoff](S08.md)・[最終hash一覧](data/s08-review.json)を参照。
今回読み取った二ownerはS08では変更していない。

| Owner | SHA-256 |
| --- | --- |
| `gbdraw/web/js/app/losat-cache.js` | `24b489d989f5fb396e2695529f6b7b6c7540c8c3bb2cbac3e2e0b23d6b6a69cf` |
| `gbdraw/web/js/app/run-analysis.js` | `cdd420cf73b9512bc5a28450f69ee39c89bc18710d083d7d566aa0b99c33957e` |

以下の回数はsourceからの演算数であり、新しいprofileの実測値ではない。
時間は[保存済みbrowser観測](data/s08-browser-summary.json)を参照する。
各区間には別の処理も含まれるため、重複処理の所要秒数・削減可能秒数とはしない。

## RW-01 — 同じraw entryを1回取得する間に検証が2回走る【後続で修正済み】

到達経路は `run-analysis.js:3732` のprepared-job loop →
`getReusableLosatCacheEntry`（1279）→ `getRawLosatCacheEntry`（260）→
`losat-cache.js` の `getCurrentRawLosatCacheEntry`（581）。wrapperに検証済み結果のmemoはない。

protein-currentの正常hitは、同じ同期呼び出し内で次を実行する。

1. 599行の `validateProteinRawEntryReferences(entry, manifest)` が、manifest全体を検証し、
   query/subjectのruntime ID Setを作り、`rawProteinTextMatchesBindings` でTSV全行を検査する
   （291–318行）。
2. 600行の `proteinRuntimeIdSets` が同じmanifest全体を再検証し、同じ2つのSetを再生成する
   （321–330行）。
3. 605行が、同じ `entry.text` と同じ内容のSetで `rawProteinTextMatchesBindings` を再実行する。

この間にawaitや状態変更はない。通常のSession/plain-object入力では2回目は同じ結果になる。
正常hit **1件につきmanifest検証2回、query/subject Set計4個、TSV全走査2回**。
VibrioのS08保存raw再解析は47件hitなので、その47件のgetter成功経路だけで
**94回のmanifest検証、188個のSet、94回のTSV走査**となる。
他の境界での検証はこの数に含めない。derived-warmについて同じ47件を通るとは主張しない。

**次回の最小変更:** 既存getter内の後段 `proteinRuntimeIdSets` と直接のTSV検査を除き、
最初の `validateProteinRawEntryReferences` に集約する。公開helperは他の利用者を調べ、
この修正のために新owner、永続cache、Worker protocol、compatibility分岐を作らない。
このgetter内ではmanifest検証1回・Set2個・TSV走査1回になるという演算数の削減を確認できる。

**守る契約:** schema・program・outfmt・args・search context、方向別binding hash、
protein-set参照、runtime ID、厳密な12列/numeric/finite検証、空rawの正常hitを維持する。
不正entryをhitに変えない。S08のSave時の検証とreaderの検証は別の境界であり、削除対象ではない。

**検証:** 既存 `tests/web/losat-cache.test.mjs` の正常hit、逆方向、未知ID、manifest不整合、
strict numeric fixtureをgetter経由でも確認する。raw/derived identityとSession exportの既存gateを実行し、
保存rawからの実Generateで検索0回・結果/provenance一致を確認する。
壁時計の閾値をテストにせず、重複呼び出し削除と既存oracleで正確性を示す。

## RW-02 — 同じmanifestを各pairで検証し直す【後続で寿命確認・修正済み】

RW-01を直しても、prepared-job loopは各entryに同じ
`workingProteinIdentityManifest` を渡し、validatorがmanifest全体を1回ずつ検証する。
また、同じrecordがself/forward/reverseの複数pairに現れるたびにruntime ID Setを作り直す。
これは全manifestの構造確認を必要なentry固有のbinding/TSV検査と一緒に繰り返している。

既存ownerには `buildValidatedProteinIdentityIndex`（`losat-cache.js:250`）、
`releaseValidatedProteinIdentityIndex`（266）、validatorの `{ identityIndex }` 引数がある。
Session readerとS08 writerもこの既存経路を使う。次回は、**manifestの内容が変わらない処理区間だけ**
同じindexを渡せるかを確認する。新しい共有cacheへ置き換える案ではない。

**未確認点:** loopにはawait・旧Session昇格・cancel/retryがある。object identityが同じだけでは
内容不変の証拠にならない。manifest更新点、昇格、source/History/Session置換を確認し、
検証済みindexが古い内容を承認しない寿命を既存owner内で表せる場合のみ採用する。
try/finallyで成功・例外・Cancelのreleaseを保証する。公開getterの直接利用者も維持する。

**検証:** RW-01のgateに加え、内容変更後に古いindexを使わないこと、旧schema昇格、
source差し替え、Cancel→member変更→retry、Clear Cache、Session/History置換を確認する。
TSVはentryごとに少なくとも1回検証する。RW-01と削減回数・秒数を二重計上しない。
全47件の同一manifestが不変なら、その範囲のmanifest検証を47回から1回にできるという
条件付きの演算数であり、現時点の採用・実時間短縮の証明ではない。

## 別枠 — 既知の再parseと、大きいが無駄と断定できない区間

共有64-entry LRUの49/64/81-table既存観測はwarm hit/missが49/0、0/64、0/81。
64/81では同じ入力のparse/filterを繰り返すことが確認済み
（[source一致と再利用根拠](data/s08-contract-reuse.json)）。
ただし現行の有限cacheとownerの問題を含むため、RW-01/RW-02の局所修正とは分ける。
**S05で却下された共有cache/共通解析class/runtime再設計を復活させない。**
今回の47-table Vibrio全体待ち時間をこの64/81-table現象のせいにしない。

Vibrioのrender往復と、render返却後のResult admission/History/DOM区間はそれぞれ大きいが、
今回のhookでは内部を分離していない。Historyだけの時間や、全部削除可能な処理とはしない。
`feature-catalog.js:499` のadmissionには既存cacheがあり、
`preview-runtime.js:798` のfeature indexも保持される。単に複数のcallerがあることを根拠に、
catalog全走査やSVG全走査が毎回重複すると決めつけない。

raw新規検索、異なる設定の再解析、必要なself/reverse evidence、Session reader/writerの
別境界の検証、明示的全path取得を「無駄」として削らない。
追加の明確な重複が見つかった場合は、同じ入力・同じ寿命・同じ結果という根拠を付けて
この一覧へ追加する。新しい測定は必要な境界のみ各3回まで、warmup最大1回、noise追加なし。


## RW-03 — merge直前のper-record manifest検証【後続で修正済み】

以下は候補登録時の調査記録。実装・保持契約・検証結果は冒頭のRW-03結果文書を参照。

`run-analysis.js`で`proteinEntries.map(entry => entry.identityManifest)`の直後に
`manifests.some(manifest => !validateProteinIdentityManifest(manifest))`を実行し、
そのまま`mergeProteinIdentityManifests(manifests)`を呼ぶ。merge ownerは各入力manifestに
同じvalidatorを再実行する。正常入力では同じ配列・同じ内容であり、この二呼び出し間に
await・mutationはない。R入力なら同じper-record manifest検証がR回余分にある。
merge後の**統合manifest**検証は別の必要な検証であり、この数には含めない。

ただし不正入力時の前段はユーザー向けreload案内を含む固有error、merge側は一般的な
invalid-manifest errorを返す。単に前段を削除するとerror wordingが変わるため、既存の
失敗契約・他のmerge callerを調べてから既存owner内で整理する。今回は実装・時間測定せず、
RW-01/RW-02の削減回数に加算しない。全体改善率や削減可能秒数は未測定。
