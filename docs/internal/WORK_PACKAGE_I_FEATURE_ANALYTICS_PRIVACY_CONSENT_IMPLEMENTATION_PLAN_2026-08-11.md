# Work package I: Feature analytics and privacy consent implementation plan

- Date: 2026-08-11
- Status: ready for implementation
- Baseline: `docs_renovation` / `1a1780d3a9a7` と文書作成時点の inventoried working tree
- Source: [`gbdraw_v0.14.0_codex_roadmap.md`](./gbdraw_v0.14.0_codex_roadmap.md) の Work package I
- Target release: `v0.14.0`
- Prerequisites: A0 baseline、H1 hosted-artifact topology、I2については Work package E の accepted-Result owner と export owner の凍結

Related contracts:

- [Web application maintenance contract](../../gbdraw/web/CLAUDE.md)
- [Browser/offline verifier](../../tools/verify_gui_offline.py)
- [Cloudflare bundle builder](../../tools/prepare_cloudflare_pages.py)
- [GitHub Pages workflow](../../.github/workflows/deploy_web.yml)
- [Web packaging tests](../../tests/test_web_packaging.py)

External implementation references (recheck at Phase 0 because provider
behaviour and settings can change):

- [Google Consent Mode concepts](https://developers.google.com/tag-platform/security/concepts/consent-mode)
- [Google privacy controls for tags](https://developers.google.com/tag-platform/security/guides/privacy)
- [GA4 page-view measurement](https://developers.google.com/analytics/devguides/collection/ga4/views)
- [GA4 automatically collected events](https://support.google.com/analytics/answer/9234069)
- [GA4 data collection overview](https://support.google.com/analytics/answer/11593727?hl=en)
- [GA4 granular location and device collection](https://support.google.com/analytics/answer/12002752?hl=en)

この文書は Work package I を実行可能な phase と証拠 gate に分解する。
文書作成時点では production code、test、public documentation、deployment、
GA property、generated artifact を変更していない。working tree には本計画と
無関係な変更があるため、実装開始時に差分を再取得し、関係しない変更を
巻き戻さない。

## 1. 結論

Work package I は一つの analytics feature ではなく、次の二つに分ける。

1. **I1: mandatory privacy remediation**
   - 現在の eager `gtag.js` と `gtag('config')` を撤去する。
   - source/local/wheel/analytics-disabled hosted artifact を Google に接続不能にする。
   - hosted bundle と production deployment の owner を一つにする。
2. **I2/I3: conditional analytics**
   - 明示的な許可後だけ GA4 を動的 load する。
   - closed registry の outcome event だけを送る。
   - notice、property 設定、実 transport、withdrawal、browser matrix の全証拠がそろった場合だけ有効化する。

I1 は release blocker である。Analytics 自体は release blocker ではない。
Enabled branch のどれか一つでも証明できなければ、`v0.14.0` は analytics-disabled
profile で出す。Privacy を妥協して analytics を残す分岐は作らない。

## 2. 現状と根拠

| Area | Current owner/evidence | 計画への影響 |
| --- | --- | --- |
| Source CSP と marker | `gbdraw/web/index.html` の CSP、`GOOGLE_ANALYTICS_SCRIPT`、`GOOGLE_ANALYTICS_NOTICE` | source は same-origin のまま維持する。hosted builder は script ではなく passive config だけを注入する。 |
| Cloudflare bundle | `tools/prepare_cloudflare_pages.py` の `_render_google_analytics_script()` と `build_cloudflare_pages_bundle()` | 現在は consent 前に外部 script と `config` を注入している。最初に反転する contract である。 |
| Cloudflare runtime | `wrangler.toml` と `gbdraw/web/cloudflare-worker.js` | `dist/cloudflare-pages` と remote Gallery proxy は一組であり、canonical hosted path の第一候補である。 |
| GitHub Pages | `.github/workflows/deploy_web.yml` の `Prepare Web Directory` | `gbdraw/web` を直接 copy する第二の production recipe で、Cloudflare builder を迂回する。 |
| Packaging test | `tests/test_web_packaging.py::test_cloudflare_bundle_includes_google_analytics_and_hosted_notice` | eager load を正として固定している。parser-time GA absence と passive config を正とする test に置換する。 |
| Local GUI | `gbdraw/cli.py` と `gbdraw/_build_support.py` | packaged `gbdraw/web` を配る。hosted config がない限り、forged storage でも hard no-op にする。 |
| Consent persistence | production Web に first-party `localStorage` / `sessionStorage` owner はない | versioned consent store は新規 owner になる。diagram/session state に混ぜない。 |
| Session and Reset | `gbdraw/web/js/services/config.js`、`history-snapshot.js`、`app/app-setup.js::resetSettings()` | consent は列挙 field に追加しない。negative contract test で分離を証明する。 |
| Accepted Result | `gbdraw/web/js/app/run-analysis.js` の dispatch、candidate staging、stale checks、atomic adoption | `render_result=success` と rendered feature snapshot は complete current commit 後にだけ出す。 |
| Export | `gbdraw/web/js/services/export.js` | SVG/interactive SVG、PNG callback、PDF save の完了点が異なる。各完了後にだけ event を出す。 |
| Error text | `gbdraw/web/js/services/error-normalization.js` | UI 用 text/stdout/stderr を analytics に再利用しない。control-flow 由来の enum を別に作る。 |
| Gallery | `gbdraw/web/gallery/index.html` と `gallery.js` | 別 document/CSP で、tutorial completion owner もない。v0.14 analytics scope から外す。 |
| Offline evidence | `tools/verify_gui_offline.py` | external request の block/record が既にある。source と wheel で forged consent matrix を追加する。 |

## 3. 固定する製品判断

### D1. Canonical hosted path は一つだけ

Phase 0 で live DNS、Cloudflare Worker、GitHub Pages custom-domain の実際の
owner を read-only に確認する。本計画の選択は次である。

```text
tools/prepare_cloudflare_pages.py
  -> dist/cloudflare-pages
  -> gbdraw/web/cloudflare-worker.js
  -> production hosted origin
```

この選択が実運用と一致すれば、GitHub Pages workflow は production publisher
から外すか、明示的な non-production hostname 用にする。実運用が逆なら作業を
止め、H1 で一度だけ canonical path を変更する。二つの production recipe を
同期し続ける設計は採用しない。

### D2. Analytics profile は disabled が default

Builder API の default は `analytics_enabled=False`、measurement ID は `None`
とする。CLI は `--enable-analytics` と valid な `--google-analytics-id` の両方を
要求する。ID の環境変数を許しても、明示 enable flag の代わりにはしない。
repository default の measurement ID は削除する。

Enabled build は first-party consent markup と、enabled flag、validated
measurement ID、notice version を表す passive same-origin metadata を marker に
注入できる。Parser-time の external `<script>`、inline `gtag()`、`dataLayer`
初期化は禁止する。

### D3. Artifact profile ごとに capability を閉じる

- source/local/wheel: hosted metadata なし、Google CSP origin なし、consent UI なし、loader は hard no-op。
- hosted disabled: measurement ID なし、Google CSP origin なし、consent UI なし、event path なし。
- hosted enabled: passive metadata と最小 CSP のみ。valid `allowed` 後に dynamic load できる。

Shared JavaScript が存在することだけで local artifact を hosted-capable と判定しない。
ただし storage、query、hash、session import だけで hosted mode を偽装できる実装は
禁止する。Capability root は builder が注入した validated metadata だけである。

### D4. Consent は三状態、versioned、fail-closed

Runtime state は `unknown | allowed | declined` だけとする。missing、expired、
old notice version、malformed value、invalid timestamp、storage exception は
`unknown` へ正規化し、外部 load を行わない。Consent の最大有効期間は180日とし、
notice version 変更時には再同意を求める。

Planning baseline の eager deployment が残した既知の `_ga` cookies は、hosted
capability の `unknown` / `declined` boot で tag を loadせず削除する。Legacy cookie
を新しい opt-in 後の client identity に持ち越さない。Prior collection data の
retention/deletion disposition は property operator と privacy reviewer が記録する。

### D5. Pre-consent queue と backfill は作らない

`unknown` / `declined` では `dataLayer` 自体を作らない。Event は破棄し、
後で許可されても再送しない。App の起動時刻、同意前の render、export、
session load を buffer しない。

### D6. Consent は project state ではない

Consent record と in-memory dedupe は canonical request、`.gbdraw-session`、
history、Undo/Redo、saved navigation、Reset Settings の対象外である。Session import
が同意を変更する code path は存在させない。

### D7. Basic tag gating を使う

Unknown/declined のまま tag を load して consent mode の cookieless ping を使う
方式は採用しない。Allowed decision を保存した後に初めて `dataLayer` と `gtag`
を初期化し、script を一回だけ append する。

Initialization order は次で固定する。

1. current capability と versioned `allowed` を再検証する。
2. loader generation を取得する。
3. `dataLayer` と local `gtag` shim を作る。
4. analytics storage を granted、ad storage/ad user data/ad personalization を denied にする。
5. 許可後の時点で `gtag('js', new Date())` command を一回だけ初期化する。
6. `gtag.js` を一回だけ append する。
7. current generation と consent を再検証する。
8. `config` を `send_page_view:false`、Google Signals/personalization false、固定 page context で呼ぶ。
9. その後の approved custom event だけを dispatch する。

### D8. Withdrawal は live tag を残さない

Allow から decline へ変える action は、最初に wrapper を disabled にし、
`window['ga-disable-MEASUREMENT_ID']` を設定し、loader generation を invalidate
する。次に configured GA cookie を clear し、stored `allowed` を remove または
`declined` で overwrite して read-back する。Decoded state が `allowed` でないことを
確認できた場合だけ reload する。Declined write が失敗しても、remove と read-back が
成功して missing/`unknown` なら安全に reload できる。

Read-back がなお `allowed`、または removal/overwrite/verification がすべて失敗する
場合は reload しない。Current document を hard-disabled のまま保ち、cookie を消し、
privacy-state error と site-data cleanup guidance を表示して enabled-profile gate を
fail にする。Old allowed record を残したまま「declined」と表示して reload しない。

別 tab の decline、policy expiry、pending loader 中の decline も同じ coordinator
を通す。Observation は二つに分ける。Withdrawal click から unload まで新しい
request が開始されないこと、新 document で script/`dataLayer`/GA cookie/request が
ゼロであることを別々に証明する。どちらかが real-tag capture で fail した場合は
enabled profile を release しない。

### D9. GA wrapper の schema と transport payload は別々に証明する

Closed wrapper は custom values だけを制御する。GA が付加する page/device/network
metadata、cookie、automatic event も capture して承認対象にする。Enhanced
Measurement、Google Signals、advertising、cross-domain linking、user-provided
data、不要な data sharing は property 側で無効化する。Page context は固定し、
hosted page は strict referrer policy を使う。

### D10. v0.14 は main app の非生物学的な bounded outcome だけを測る

Gallery、Palette Explorer、tutorial、session save/load、label edit、legend edit は
初回 registry に入れない。`multi_record_layout` も「record count >= 2」という
biological-input-derived bucket になるため外す。Gallery に event を追加する将来変更は、
直接訪問者を含む同じ consent UI/state/CSP/transport gate を別途満たす。

### D11. Event registry は closed かつ result-owned

Analytics transport API は enum primitive だけを受ける。`Error`、request、
resources、reactive state、session object、File、SVG、任意 text を受け取らない。
Unknown event/key/value、required key 不足、extra key は fail-closed で送らない。

Feature classification は二段階にする。Genuine renderer dispatch 時に、validated
current candidate から allowlisted capability booleans だけへ eligibility を射影する。
Atomic acceptance 後に、accepted Result owner が active/rendered booleans を射影する。
Transport から raw request/resource を到達不能にし、export は live draft ではなく
committed active snapshot を使う。

### D12. Analytics failure は application failure ではない

Tag block、timeout、collector error、storage failure、sanitizer rejection は、
upload、ready、Generate、preview、save、load、export の state と error UI を変更しない。
Analytics 用 Promise を user workflow から await しない。

### D13. Enabled profile には三者 sign-off が要る

- maintainer: product purpose、feature/KPI registry、release decision
- privacy/legal reviewer: controller/contact、notice、retention、recipient/transfer、withdrawal/deletion wording
- property operator: account ownership、access、retention、Enhanced Measurement/Signals/ads/data-sharing、deletion、kill switch

Property operator はさらに event-scoped custom dimensions と schema-filtered report/query
を所有する。Event を送信できることだけで KPI を運用可能とみなさない。

いずれかが不在または未承認なら disabled profile を選ぶ。実装者は法的適合性を
推定して checkmark を付けない。

## 4. Artifact and network matrix

| Profile | Build input | Consent surface | CSP | Forged stored `allowed` | External-blocked expectation |
| --- | --- | --- | --- | --- | --- |
| Source checkout | `gbdraw/web` | hidden | same-origin | ignored | zero external attempts; core workflow passes |
| Installed local GUI/wheel | packaged `gbdraw/web` | hidden | same-origin | ignored | zero external attempts; core workflow passes |
| Hosted disabled | canonical hosted bundle, flag off | hidden | no Google origin | ignored | zero external attempts; core workflow passes |
| Hosted enabled, unknown | flag on + valid ID | visible | captured Google allowlist | still unknown if record invalid/expired | zero Google attempts/cookies |
| Hosted enabled, declined | same | persistent settings available | same | remains declined | zero Google attempts/cookies |
| Hosted enabled, allowed | same | persistent withdrawal available | same | valid record activates | consented GA attempt may occur; blocking it cannot affect app |
| Hosted enabled, revoked | reload into declined | persistent settings available | same | declined wins | no request initiated after withdrawal; GA cookies cleared |

Hosted allowed + collector blocked は zero-attempt offline case ではない。許可済みの
analytics request attempt は認めるが、application workflow は完全に成功させる。

## 5. Consent storage and transition contract

### 5.1 Wire shape

一つの namespaced key を使う。

```text
key: gbdraw.analyticsConsent.v1
```

JSON shape:

```json
{
  "schema": 1,
  "noticeVersion": 1,
  "status": "allowed",
  "decidedAtMs": 1786377600000,
  "expiresAtMs": 1801929600000
}
```

`status` は persisted record では `allowed | declined` のみで、missing record が
`unknown` を表す。整数 timestamp、順序、最大180日、current notice version を
厳密に検証する。未知 key は将来の誤読を避けるため invalid とする。Invalid record
は可能なら削除するが、削除失敗時も runtime は `unknown` のままにする。

### 5.2 Transition table

| Current | Input | Next | Side effects |
| --- | --- | --- | --- |
| boot | missing/invalid/expired/storage error | `unknown` | known legacy GA cookie を削除し、notice を表示。Google/dataLayer/event はゼロ |
| `unknown` | allow | `allowed` | valid record 保存後に new in-memory visit と loader 開始 |
| `unknown` | decline | `declined` | declined record を保存、notice を閉じる、外部 side effect なし |
| `declined` | allow from Privacy settings | `allowed` | new visit、dynamic load。過去 event は送らない |
| `allowed` | decline; stored state verified non-allowed | `declined` または `unknown` + reload | dispatch block、GA disable、loader invalidate、cookie cleanup、remove/overwrite/read-back、reload |
| `allowed` | decline; stored state cannot be verified non-allowed | in-memory declined; no reload | hard-disable と cookie cleanup を維持し、privacy error。Enabled gate fails |
| `allowed`, loader pending | decline | 上の persistence result に従う | generation mismatch で onload/config/event を無効化 |
| any | cross-tab valid allow | `allowed` | policy 再検証後にその tab の new visit を開始 |
| any | cross-tab valid decline | `declined` + reload if needed | 同じ withdrawal coordinator を通す |
| `allowed` | policy expiry timer | `unknown` + reload | dispatch block、cleanup、stored record を invalidation |
| any | duplicate action | unchanged | idempotent。script、event、reload を重複させない |
| `allowed` | script/collector failure | `allowed` | app は継続、analytics は no-op。Privacy settings は使用可能 |

Withdrawal 後も declined preference は残す。Cookie cleanup と consent record
cleanup を同じ操作として扱わない。

## 6. Event, feature, and KPI registry

### 6.1 Analytics visit

Analytics visit は main app の一つの top-level document における uninterrupted
`allowed` epoch とする。Reload または decline 後の再 allow は新しい visit である。
Visit ID、user ID、session file ID は送らない。Dedupe set は memory-only とする。

### 6.2 Event registry v1

| Event | Required params | Values | Dedupe | Truth owner |
| --- | --- | --- | --- | --- |
| `workflow_start` | `event_schema` | `1` | once/visit | first genuine renderer dispatch after preflight |
| `feature_progress` | `event_schema`, `feature_id`, `stage` | closed feature enum; `eligible\|rendered\|exported` | once/feature/stage/visit | genuine-dispatch eligibility; accepted active snapshot; export completion |
| `render_result` | `event_schema`, `status`, `trigger`, `error_code` | `success\|failure`; trigger/error enums below | none; one/genuine attempt | current accepted commit or current genuine failure |
| `export_figure` | `event_schema`, `format` | `svg\|interactive_svg\|png\|pdf` | none; one/completed app handoff | format-specific committed-Result Blob/save handoff; OS disk completion is not observable |

Trigger enum:

```text
analysis_generate
manual_render
automatic_render
session_replay
```

Initial error enum:

```text
none
worker_unavailable
request_rejected
render_failed
unsafe_result
commit_failed
```

Phase 4 で実 control flow と対応しない code は削除し、不足する stage は自由 text
ではなく enum review で追加する。Pre-dispatch validation、debounce replacement、
pending supersession、stale response、cancellation、direct patch は
`render_result` を出さない。

### 6.3 Initial feature registry

| `feature_id` | Eligible predicate at genuine dispatch | Rendered/active predicate | Safe owner |
| --- | --- | --- | --- |
| `record_display_start` | validated candidate に retained Work C の display-start control が利用可能な record がある | 少なくとも一つの explicit start override が committed | Work C safe dispatch/accepted flags |
| `manual_feature_placement` | validated candidate に supported placement target/slot がある | non-Auto placement override が一つ以上 committed | Work D safe dispatch/accepted flags |
| `region_annotations` | current mode の annotation surface が candidate に適用可能 | effective region annotation が一つ以上 committed | annotation capability/summary booleans |
| `custom_tracks` | current mode の custom-track surface が candidate に適用可能 | explicit custom track-slot layout が committed | track capability/summary booleans |
| `depth_tracks` | current mode の depth-track surface が candidate に適用可能 | depth track が一つ以上 committed | track capability/summary booleans |

Phase 0/4 で各 predicate を retained feature contract と照合する。Genuine dispatch
時の safe eligibility と accepted Result の active state を truthful に作れない ID は
registry から外す。Record count、sequence size、filename、record identity、mode、
input type を parameter として追加しない。

LOSAT program/presentation ID は initial registry に含めない。LOSAT の truthful な
attempt boundary は renderer より前の validated analysis-job dispatch にあり、その失敗を
renderer-owned denominator へ同じ意味で入れられないためである。Post-analysis survivor
だけを数える近似は採らない。将来追加する場合は、LOSAT job dispatch、bounded outcome、
failure semantics、KPI consumer を一組で設計する。

### 6.4 KPI definitions

| KPI | Formula | Decision use | Caveat |
| --- | --- | --- | --- |
| Successful adoption | `feature_progress(rendered) / feature_progress(eligible)` per feature | genuine dispatch に到達した eligible visit から成功利用への比率 | consenting visits のみ。preflight-only failure は範囲外 |
| Export conversion | `feature_progress(exported) / feature_progress(rendered)` per feature | committed use が output に到達したか | one-per-visit dedupe の比率 |
| Render failure rate | failed `render_result` / all `render_result` | renderer reliability guardrail | scheduler churn と cancellation を除外 |
| Privacy violations | prohibited request/cookie/payload count | release guardrail | target は常に exactly zero |
| Collector independence | analytics-blocked core workflow failures | release guardrail | target は常に exactly zero |

Baseline が得られるまで adoption target は置かない。Active feature を attempt outcome
へ安全に相関する別契約ができるまで per-feature failure rate は報告しない。

## 7. Ownership and change surface

| Responsibility | File/module owner | Planned change |
| --- | --- | --- |
| Hosted profile/config/CSP | `tools/prepare_cloudflare_pages.py` | eager injection と hard-coded default ID を削除。disabled default、passive config、profile-specific CSP を生成 |
| Production path | `.github/workflows/deploy_web.yml`, `wrangler.toml`, `gbdraw/web/cloudflare-worker.js` | Phase 0 decisionに従い一つだけを production owner にする |
| Markup and source CSP | `gbdraw/web/index.html` | inert hosted-build markers、strict source CSP/referrer policy。source/local/disabled artifact に consent markup を置かない |
| Consent UI adapter | new `gbdraw/web/js/app/analytics-consent.js` | Vue refs、notice/details、allow/decline/withdraw actions |
| Consent/registry/loader | new `gbdraw/web/js/services/analytics.js` | capability、storage reducer、dedupe、schema validation、dynamic loader、disable/cookie cleanup |
| Composition | `gbdraw/web/js/app/app-setup.js` | service と UI adapter を一回だけ wiring。diagram state owner にしない |
| Render outcomes | `gbdraw/web/js/app/run-analysis.js` と Work E の final accepted-Result owner | genuine dispatch/result hook と safe feature facts を注入 |
| Candidate staging | `gbdraw/web/js/app/candidate-render.js` | event owner にしない。staging success は committed success ではない |
| Export outcomes | `gbdraw/web/js/services/export.js` | format別 completion 後に reporter callback。committed snapshot を使用 |
| Negative persistence contract | `config.js`, `history-snapshot.js`, Reset owner | consent field を追加しない。test だけで absence を固定 |
| Pure unit tests | new `tests/web/analytics.test.mjs` | reducer、registry、dedupe、loader race、no-backfill、cookie plan |
| Packaging tests | `tests/test_web_packaging.py` | three profile output、deployment owner、source CSP、parser-time script absence |
| Browser privacy acceptance | new focused browser test under `tests/` or `tests/web/` | exact built artifact、request/cookie/sentinel/state matrix。CIで確実に実行される runner を選ぶ |
| Offline matrix | `tools/verify_gui_offline.py` とその contract test | required subcommands を使い、source に加えて exact outer wheel を temp に展開して checkout 外の `gbdraw/web` を runtime smoke する |
| Web maintenance contract | `gbdraw/web/CLAUDE.md` | same-origin source/local rule と enabled hosted-build-only analytics exception を明示 |
| Reporting/property contract | authorised GA property と internal evidence/query specification | event-scoped dimensions、`event_schema=1` filter、KPI formula/timezone、synthetic reconciliation、access owner |
| Public privacy docs | `README.md`, `docs/INSTALL.md`, `docs/REFERENCE/web-app.md`, `docs/FAQ.md` | enabled/disabled/local boundaries、sent/unsent、withdrawal、provider metadata |
| Generated notice | `tools/generate_open_source_notices.py` と generated notice | generator owner を更新して再生成。generated HTML を手編集しない |

Gallery file、canonical request schema、session schema、Python API、CLI analytics flag は
この Work package の change surface にしない。

## 8. Execution phases

### Phase 0: Baseline, deployment owner, and retain/disable entry gate

Status: pending

#### Depends on

- A0 integration baseline
- H1 owner が参加できること

#### Owners

- `tools/prepare_cloudflare_pages.py`
- `.github/workflows/deploy_web.yml`
- `wrangler.toml`
- Cloudflare Worker/Pages と GitHub Pages の read-only deployment metadata
- Work package I decision/evidence ledger

#### Work

1. 開始 commit、branch、dirty-tree inventory、package version、source CSP、current builder output を記録する。
2. `gbdraw.app` の DNS、実 response headers/build stamp、Cloudflare deployment、GitHub Pages custom-domain、現在の GA/cookie behaviour を read-only に確認する。
3. Canonical production path を H1 と一緒に承認する。本計画の Cloudflare 選択と異なる場合はここで計画を改訂する。
4. Current eager script/config と二重 recipe を reproducible test/output で記録する。現行 eager test の pass を privacy evidence と扱わない。
5. Enabled branch の maintainer、privacy/legal reviewer、property operator を記録し、prior eager-deployment data/cookie の disposition を decision item にする。担当者または staging property がない場合は disabled decision を記録する。
6. Initial feature IDs、eligibility/active predicates、error enum、notice version、180-day consent term を retained contracts と照合する。
7. Enabled/disabled の選択期限を Gate 0 より前に設定する。未決定のまま release candidate を作らない。

#### Tests and evidence

```bash
git status --short --untracked-files=all
git rev-parse HEAD
python -m pytest tests/test_web_packaging.py -q -k "cloudflare or deploy"
python tools/verify_gui_offline.py check-assets
python tools/verify_gui_offline.py smoke-test
```

Current baseline で privacy target が fail することは defect evidence であり、test を
green に見せるため期待値を維持しない。External ownership check は command、日時、
observed target/build stamp を ledger に残す。

#### Completion gate

- Canonical production owner と bundle recipe が一つに決まっている。
- Enabled branch に必要な三者 owner と staging capture route があるか、disabled decision が記録されている。
- Event/feature registry の各 candidate に truthful predicate と owner がある。
- Baseline evidence が exact commit に結び付く。

#### Stop condition

Live production owner、DNS authority、または Gallery remote-asset dependency が不明なまま
deployment workflow を変更しない。Read-only evidence で解消できなければ maintainer に
判断を求める。

### Phase 1: Canonical hosted artifact and mandatory analytics-disabled profile

Status: pending

#### Depends on

- Phase 0 canonical path decision

#### Owners

- `tools/prepare_cloudflare_pages.py`
- selected production workflow/config
- `gbdraw/web/index.html` markers/CSP
- `gbdraw/web/CLAUDE.md`
- `README.md`, `docs/INSTALL.md`, `docs/REFERENCE/web-app.md`, `docs/FAQ.md`
- `tests/test_web_packaging.py`
- documentation contract tests

#### Work

1. `_render_google_analytics_script()` と eager `gtag('config')` injection を削除する。
2. Hard-coded default measurement ID を削除し、builder API を disabled default にする。
3. `--enable-analytics` と valid `--google-analytics-id` の AND gate を追加する。ID は strict `G-...` pattern と HTML escaping を通す。
4. Disabled output では Google origin、measurement ID、consent notice/config を一切生成しない。
5. Enabled output では passive hosted config と notice-version metadata だけを注入し、external script/inline gtag/dataLayer を生成しない。Consent markup は Phase 2 で first-party marker injection として追加する。
6. Source CSP を変更せず、enabled built artifact だけに最小 Google host allowlist を与える。実 capture 前の broad wildcard を最終値にしない。
7. Canonical path が Cloudflare の場合、GitHub Pages の automatic production publication と `gbdraw.app` CNAME ownership を外す。Non-production use を残すなら hostname、trigger、artifact purpose を明記する。
8. Canonical builder が browser wheel、build stamp、Gallery remote manifest、headers を一度だけ生成することを固定する。
9. Maintainer approval の下で `gbdraw/web/CLAUDE.md` を更新し、source/local の same-origin rule を維持しつつ、retained enabled canonical hosted build だけに許す consent-gated exception と必須 test を記載する。承認前に Phase 3 の external loader を実装しない。
10. 現行 public docs の unconditional「hosted app は GA を使う」という主張を削除する。Phase 1 時点では disabled-safe behavior を記載し、enabled evidence が後で揃った場合だけ Phase 5 で opt-in 説明へ更新する。

#### Tests

- disabled builder output の rendered HTML/CSP に measurement ID、Google origin、hosted analytics config、consent markup がなく、shared service が存在しても capability detector が常に no-op になる。
- enabled builder output に passive config と approved CSP があるが、external `<script>`、`gtag(`、`dataLayer` がない。
- flag だけ、ID だけ、invalid ID は fail closed または disabled で、暗黙 enable しない。
- source `index.html` と installed Web tree は same-origin CSP のまま。
- selected production workflow は canonical bundle だけを publish し、別 recipe で rebuild しない。
- Cloudflare Worker を選ぶ場合、remote Gallery assets が Worker proxy なしの Pages artifact に誤って出ない。
- Web maintenance contract が source/local same-origin invariant と conditional hosted exception を矛盾なく記載する。
- Public docs に実 artifact より先行した GA 利用 claim がない。Disabled branch はこの copy で release可能である。

#### Completion gate

I1 の最小安全状態として、ordinary build と canonical hosted disabled build のどちらも
Google に接続できず、eager analytics contract が production/test/public docs から消えて
いる。Enabled runtime exception は authoritative guidance で承認済みか、enabled branch
が停止されている。

#### Stop condition

Canonical bundle を二つの deploy target で同じと仮定しない。Remote Gallery behavior、
headers、CNAME の差が一つの artifact で表現できなければ H1 の topology decision に戻る。

### Phase 2: Consent owner and UI with a no-op transport

Status: pending; skip when analytics is disabled at scope freeze

#### Depends on

- Phase 1
- Enabled branch の owner assignment

#### Owners

- new `gbdraw/web/js/services/analytics.js`
- new `gbdraw/web/js/app/analytics-consent.js`
- `gbdraw/web/js/app/app-setup.js`
- `gbdraw/web/index.html`
- `tools/prepare_cloudflare_pages.py` の enabled-only first-party markup injection
- new `tests/web/analytics.test.mjs`

#### Work

1. Hosted capability detector、strict storage decoder、three-state reducer、expiry timer、legacy GA-cookie cleanup、in-memory visit/dedupe owner を pure service に実装する。
2. Transport はこの phase では injectable no-op/test spy のままとし、Google URL を呼ばない。
3. Source HTML の inert marker に対し、builder が enabled profile だけへ first-party consent/details/footer markup を注入する。
4. Enabled hosted build の `unknown` だけに non-modal first-choice notice を表示する。App の操作/focus を trap せず、無視した状態は `unknown`/analytics-off のままにする。Allow と decline を同等に操作可能にし、silent-allow の dismiss を置かず、keyboard order、screen-reader labels、narrow layout を実装する。
5. Footer の persistent `Privacy settings` から details、allow/decline/withdraw を再表示できるようにする。
6. Storage failure は UI を壊さず `unknown`/no analytics に戻す。Allow を persistence できない場合は外部 load しない。
7. `storage` event、duplicate action、notice-version change、clock expiry を同じ reducer に通す。
8. Consent refs を Vue diagram state に追加せず、app composition から明示 dependency として渡す。

#### Tests

- valid allowed/declined、missing、extra key、wrong type、NaN/fractional timestamp、future/expired/overlong duration、old notice version。
- `getItem`/`setItem`/`removeItem` exception を個別および組合せで試す。特に current `allowed` に対する (a) remove成功+declined write失敗、(b) remove失敗+declined overwrite/read-back成功、(c) remove/overwrite/read-backで non-allowed を確認不能、を固定する。
- allow/decline/re-allow、reload hydration、expiry、cross-tab、double click。
- Unknown/declined boot に pre-seed した既知の legacy GA cookies が tag/request なしで消える。
- unknown/declined の report call が spy transport、queue、dedupe を変更しない。
- Allow 前に起きた event が Allow 後に出ない。
- session save/load、history snapshot、Undo/Redo、Reset Settings に storage key/value が入らず、session import で consent が変わらない。
- Source/local capability では forged valid record が UI/transport を起動しない。

#### Completion gate

全 consent transition が pure unit test と UI state test で再現でき、external transport を
一切追加せずに persistence/isolation/accessibility contract が成立する。

#### Stop condition

Consent を `state.js`、session writer、history snapshot に入れないと UI を組めない設計に
なった場合は composition boundary を直す。Schema/session version は上げない。

### Phase 3: Consent-gated GA loader, property-safe defaults, and withdrawal

Status: pending; retained analytics only

#### Depends on

- Phase 2
- staging measurement property

#### Owners

- `gbdraw/web/js/services/analytics.js`
- hosted CSP/referrer policy in `tools/prepare_cloudflare_pages.py`
- `tests/web/analytics.test.mjs`
- focused browser test transport fixture

#### Work

1. Dynamic script loader を generation token と one-load Promise で実装する。
2. D7 の initialization order と consent/ad defaults を実装する。
3. `config` に `send_page_view:false`、signals/personalization false、fixed page context を設定する。
4. Loader timeout/error を reporter no-op に変え、application state/error banner に伝播させない。
5. Withdrawal coordinator に dispatch block、GA disable flag、pending generation invalidation、cookie cleanup、allowed-record remove/overwrite/read-back、verified-non-allowed の場合だけの reload を実装する。Verification不能時は reload せず enabled gate を fail する。
6. `_ga` と measurement-specific GA cookie の host/path/domain cleanup policy を明示し、未知 cookie を broad delete しない。
7. Cross-tab decline と expiry を withdrawal coordinator に接続する。Pending script の `onload` は consent/generation を再検証してから config する。
8. Final CSP host set は test stub ではなく staging real-tag capture の endpoint inventory に基づき最小化する。

#### Tests

- Unknown/declined では `document.createElement('script')`、`dataLayer`、cookie、transport call がゼロ。
- Allowed の同時/重複 call でも script append と config は一回。
- Loader success、error、timeout、CSP block が app behavior を変えない。
- Allow 直後の decline、script fetch 中の decline、onload と decline の競合で stale config/event が出ない。
- Withdrawal の action order、partial storage failure matrix、verified reload/no-reload branch、cookie cleanup。
- Click から unload までの zero-new-request と、reload 後 new document の zero script/dataLayer/cookie/request を別 oracle にする。
- Cross-tab decline と expiry で同じ停止 contract が成立する。
- Automatic page view、enhanced-measurement surrogate event、unapproved endpoint を test fixture が検出する。

#### Completion gate

Stubbed browser transport で全 race が deterministic に通り、staging real-tag capture で
withdrawal 後の新規 request がゼロ、automatic pageview がゼロ、cookie/endpoint inventory
が approved contract と一致する。

#### Stop condition

Real tag が withdrawal 後に beacon/request を開始する、unapproved automatic event/field を
抑止できない、必要 endpoint が過度に broad な CSP を要求する、または current
COOP/COEP cross-origin isolation と両立しない場合は Phase 1 の disabled profile に戻る。
Analytics のために worker/isolation headers を弱めない。

### Phase 4: Closed event registry and committed-outcome hooks

Status: pending; retained analytics only

#### Depends on

- Phase 3
- Work package E の current accepted-Result/renderer-attempt owner
- retained export owner
- feature predicates frozen against final C/D/G contracts

#### Owners

- `gbdraw/web/js/services/analytics.js`
- a narrow safe-feature projection at the accepted-Result boundary
- `gbdraw/web/js/app/run-analysis.js` / final Work E coordinator
- `gbdraw/web/js/services/export.js`
- focused Node tests for run/export hooks

#### Work

1. Event registry v1 を data として一か所に定義し、required/exact keys と enum values を検証する。
2. Transport API を event-specific methods または discriminated primitive payload に限定する。Generic object spread API は作らない。
3. First genuine renderer dispatch で `workflow_start` を once/visit 送り、同じ validated candidate から safe feature eligibility snapshot を送る。
4. Current candidate の full atomic adoption 後に `render_result(success)` と safe active/rendered snapshot を送る。
5. Dispatched current attempt の genuine failure だけに fixed error enum を付ける。UI error normalizer や exception message を渡さない。
6. Debounce replacement、pending replacement、stale/superseded、cancel、direct patch、preflight validation は event を出さない。
7. Accepted Result と一緒に memory-only safe feature snapshot を更新する。Raw request/resources、feature names/IDs、record IDs は保持しない。
8. SVG/interactive SVG、PNG callback、PDF save の application-level Blob/save handoff 後に `export_figure` と active feature の exported stage を送る。失敗/キャンセルでは送らず、browser/OS の disk completion を主張しない。

#### Tests

- Unknown event、extra/missing key、unknown enum、non-primitive、free text、Error/object input を reject する。
- Dispatch eligibility と accepted active snapshot に feature IDs と booleans 以外が存在しない。
- Success は commit 前/staging 後には出ず、complete current adoption 後に一回だけ出る。
- Failure/cancel/stale/supersede/debounce/direct-patch matrix と exact counts。
- Genuine dispatch 後に失敗した eligible visit は denominator に入り、preflight-only failure は入らない。
- Feature eligible/rendered/exported は once/feature/stage/visit。
- Export formatごとの async completion ordering と committed Result ownership。
- Raw error、filename、session title、record/accession/gene/coordinate sentinel を reporter に渡せない。

#### Completion gate

KPI の各 numerator/denominator が emitted registry だけから計算でき、failed genuine
dispatch を含む synthetic truth table と scheduler/transaction tests の期待値が一致する。
Truth owner のない event/feature は registry から削除されている。

#### Stop condition

Metric を作るために raw canonical object、biological count bucket、free text、Gallery event、
または click-level instrumentation が必要になった場合は metric/ID を defer する。

### Phase 5: Enabled privacy details, reporting, and property operations

Status: pending; retained analytics only. Disabled-copy cleanup is mandatory in Phase 1

#### Depends on

- Phases 3 and 4 payload inventory

#### Owners

- `gbdraw/web/index.html` user-facing copy
- `README.md`
- `docs/INSTALL.md`
- `docs/REFERENCE/web-app.md`
- `docs/FAQ.md`
- `gbdraw/web/CLAUDE.md`
- analytics-property operations record
- generated open-source notice owner if affected

#### Work

1. Main notice と Privacy details を actual custom/automatic payload inventory に合わせる。
2. Controller/contact、purpose、provider、cookie/client identifier、browser/device/network/approximate-location metadata、retention、recipient/transfer、withdrawal、deletion/contact を記載する。
3. Sequence、file/file name、record/accession、organism、gene/product/locus、annotation label、coordinate、LOSAT raw data、figure、raw error、session title が送信禁止であることを明記する。
4. Local/source/wheel は telemetry-incapable、hosted disabled は analytics UI/path なし、hosted enabled は opt-in のみという artifact distinction を public docs に反映する。
5. Phase 1 で承認済みの Web maintenance contract を actual endpoint/payload/test inventory と照合し、必要なら狭める。External exception をこの phase で初めて導入しない。
6. Property operator が Enhanced Measurement、Signals、ads、cross-domain、user-provided data、data sharing、retention、access、deletion、kill switch と prior eager-deployment data disposition を screen/export evidence で確認する。
7. `event_schema`、`feature_id`、`stage`、`status`、`trigger`、`error_code`、`format` を retained event-scoped reporting dimensions として登録し、`event_schema=1` で filter する saved Exploration または Data API query specification を一つ owner-owned にする。
8. Reporting specification は property timezone、date range、numerator、denominator、dedupe assumption、zero-denominator handling、consenting-sample caveat を固定する。Synthetic event truth table の期待値と report output を照合する。
9. Privacy/legal review で applicable jurisdictions、DNT/GPC、minor users、international transfer、consent record、prior unconsented collection の扱いを決める。
10. Privacy/legal reviewer と maintainer の approval、date、notice version を ledger に残す。

#### Tests and review

- Documentation contract tests が enabled/disabled/local の主張と build behavior を照合する。
- UI copy と reference docs の sent/unsent/event/version が一致する。
- Synthetic registry fixture から successful adoption、export conversion、render failure rate を計算し、saved report/query と同じ値になる。
- Generated notice を変更する場合は generator command から再生成し、generated diff を別 review する。
- External legal/property approval は automated test の代替にせず、evidence record として扱う。

#### Completion gate

Observed transport、property settings、UI、public docs が同じ notice version と data inventory
を記述し、三者 approval がある。Approval がない場合は disabled branch を選ぶ。

#### Stop condition

「bio data を送らない」という狭い説明だけで GA の automatic metadata/cookie を省略しない。
Controller/contact、retention、transfer、withdrawal 文面を実装者が推測しない。

### Phase 6: Exact-artifact browser privacy and offline evidence

Status: pending

#### Depends on

- Disabled branch: Phase 1
- Enabled branch: Phases 2–5
- H1 の provisional artifact recipe と clean-staging owner。Phase 6 の開始時に、
  tested hosted bundle と同一 source commit から fresh outer wheel を一度 build し、
  exact path と SHA-256 を固定できること

#### Owners

- canonical bundle builder
- H1 provisional wheel builder/manifest owner
- `tests/web/analytics.test.mjs`
- focused Playwright/Python-Playwright privacy acceptance
- `tools/verify_gui_offline.py`
- `tests/test_web_packaging.py`
- CI workflow that actually runs the selected browser test

#### Work

1. H1 の declared recipe を使い、新規の空 staging/output directory へ同一 source commit の provisional outer wheel を build する。最低限の wheel-only invocation は `python -m build --wheel --outdir "${GBDRAW_I_ARTIFACT_DIR}"` とし、実際の H1 wrapper/command、commit、絶対 path、SHA-256 を ledger に残す。既存 `dist/` や glob で選んだ wheel を入力にしない。
2. Fresh temporary directories に source/local、hosted disabled、hosted enabled bundle を build し、provisional wheel と同じ source commit であること、および各 path/hash を記録する。
3. Browser context を各 scenario ごとに新規作成し、storage/cookie/request logs を分離する。
4. Google script/collector は local stub/interception で応答し、live property へ automated data を送らない。
5. `fetch`/XHR/script/image/beacon/pagehide を含む request の URL、method、body、headers、referrer と前後 cookie を記録する。
6. Filename、accession、organism、gene、session title、coordinate、raw error、sequence fragment、SVG text に一意 sentinel を入れる。
7. Source/wheel は forged allowed と corrupt storage を seed し、Circular/Linear full workflow で zero external attempts を証明する。
8. Enabled allowed + collector blocked は analytics attempt を許すが、ready/generate/save/reopen/export が成功することを別 oracle で確認する。
9. Desktop と narrow viewport、mouse と keyboard の privacy controls を確認する。
10. `verify_gui_offline.py` に exact outer wheel path を受け、hash を記録し、temporary directory へ展開した `gbdraw/web` を checkout 外から serve する `smoke-wheel` subcommand を追加する。既存 source `smoke-test` の pass を wheel runtime evidence に流用しない。

#### Required browser scenarios

| Scenario | Script/request/cookie oracle | App oracle |
| --- | --- | --- |
| unknown first visit | zero/zero/zero | notice usable; app initializes locally |
| decline | zero/zero/zero | all features remain available |
| decline + reload | zero/zero/zero | choice retained; settings available |
| allow | exactly approved script/config/events/cookies | app unchanged |
| allow + reload | no pre-consent backfill; approved traffic only | app unchanged |
| allow -> revoke, persistence verified | click-to-unload new request zero; GA cookies cleared | reload returns usable declined/unknown app with no tag/dataLayer |
| allow -> revoke, persistence unverified | current document starts no new request; cookies cleared | no reload; visible privacy-state error; enabled gate fails |
| revoke + reload | zero/zero/zero | choice retained |
| decline -> allow | only post-allow events | app unchanged |
| corrupt/old/expired storage | zero/zero/zero | resolves unknown |
| unknown/declined + legacy GA cookies | no request; cookies cleared | app/choice state unchanged |
| storage exception | zero/zero/zero | app usable; consent unavailable/fail-closed |
| double action | one transition/load at most | no duplicate UI/reload |
| revoke during load | stale loader cannot config/send | usable declined reload |
| cross-tab revoke | all allowed tabs stop/reload | no project-state loss |
| collector blocked/timeout | consented attempt only | ready/generate/save/reopen/export pass |
| disabled hosted profile | zero/zero/zero; no UI/config/CSP path | core workflow passes |
| local/source/wheel + forged allow | zero attempts and cookies | Circular/Linear workflow passes |

#### Completion gate

Source smoke、exact outer-wheel inspection、checkout 外 extracted-wheel runtime smoke を
別 evidence として持つ。Exact artifact hash に結び付いた request/cookie/sentinel
evidence が全 scenario で pass し、
CI が同じ focused browser gate を実際に実行する。Code inspection や static string test
だけで network gate を完了扱いにしない。この時点の wheel evidence は H1 provisional
qualification と明記し、release-final evidence には昇格させない。H-final/J は final
candidate artifact に対して同じ inspection、checkout 外 smoke、privacy scenario を
再実行する。

#### Stop condition

Playwright runner が CI に含まれていない場合は test file の存在を evidence にしない。
Node runner を追加するか、既存 pytest browser lane で実行される Python Playwright test
へ移す。

### Phase 7: Broad verification, release handoff, and diff audit

Status: pending

#### Depends on

- Phase 6
- Gate 0 の final enabled/disabled decision

#### Focused gate

最終 test 名に合わせて command を ledger で固定する。最低限:

```bash
node tests/web/analytics.test.mjs
python -m pytest tests/test_web_packaging.py -q -k "analytics or cloudflare or deploy or offline"
python tools/verify_gui_offline.py check-assets
python tools/verify_gui_offline.py smoke-test
GBDRAW_I_EXACT_WHEEL_PATH=/absolute/path/from-the-candidate-manifest/gbdraw.whl
python tools/verify_gui_offline.py inspect-wheel "${GBDRAW_I_EXACT_WHEEL_PATH}"
python tools/verify_gui_offline.py smoke-wheel "${GBDRAW_I_EXACT_WHEEL_PATH}"
```

Browser privacy test の exact command と artifact path/hash を同じ record に残す。
ここで参照する wheel は Phase 6 で固定した H1 provisional wheel であり、下記
`python -m build` の出力へ暗黙に差し替えない。

#### Broad gate

```bash
python -m pytest tests/ -v -m "not slow"
python -m build
```

加えて、repository が宣言する complete Node unit suite、supported Chromium project、
additional supported browser smoke、built wheel/local GUI offline smoke を実行する。
30分級の test は repository timeout contract に従い、途中出力を監視する。

#### Diff audit

次を別々に review する。

1. Production: source CSP、analytics service/UI、render/export hooks、builder/deploy path
2. Tests: static、Node、browser、offline、CI invocation
3. Documentation: UI copy、README/INSTALL/reference/FAQ、roadmap/plan
4. Generated artifacts: hosted bundles、browser wheel、notices。tracked/generated owner を確認し、`dist/` や browser wheel を commit しない

最後に `git diff --check`、unexpected file inventory、secret/measurement-ID scan、
production bundle hash comparison を行う。

#### Completion gate

- Enabled branch: I1/I2/I3、all browser/payload/property/copy gates が pass。
- Disabled branch: I1、complete-absence/offline/docs gates が passし、enabled claim/path がない。
- Work package J traceability table に exact commands、environment、artifact hashes、results が登録される。
- Production/test/docs/generated diff に unexplained change がない。

#### Stop condition

Focused gate だけで release-ready としない。Broad regression、exact artifact、CI invocation、
diff category のどれかが未確認なら status は pending のままにする。

## 9. Dependencies and parallelism

```text
A0
 └─ H1 canonical hosted topology
     └─ I Phase 0
         └─ Phase 1 mandatory disabled-safe artifact
             ├─ disabled release branch ───────────────┐
             └─ retained branch: Phase 2 -> Phase 3   │
                                      E/C/D/G freeze -> Phase 4
                                      Phase 3 + 4 -> Phase 5
                                                     └─ Phase 6
disabled branch absence evidence ──────────────────────┘
                                                        -> Phase 7
                                                        -> Gate 0 / A1 / H-final / J
```

- Phase 0/1 は B–G と並行できる。Current eager GA remediation を E の完成まで待たない。
- Phase 2/3 の consent infrastructure は feature registry と独立して進められるが、enabled branch owner/sign-off がない場合は作らない。
- Phase 4 は Work E の final current-attempt/accepted-commit semantics と export completion owner を待つ。
- Phase 5 copy は observed Phase 3/4 payload より先に final 承認しない。
- A1 は enabled/disabled decision と matching public contract が freeze した後に同期する。
- Phase 6 は同一 commit の H1 provisional wheel/hosted bundle を qualification input にする。
- H-final は停止した final candidate から exact artifacts を一度 build し、J はそれを
  rebuild せずに受け取る。J は Phase 6 と同じ wheel/privacy commands を final hashes
  に対して再実行する。

## 10. Risks and mandatory stop conditions

| Risk | Detection | Required response |
| --- | --- | --- |
| Live deployment owner が repository assumptions と違う | DNS/headers/build stamp/workflow audit | Workflow を触らず H1 decision を改訂 |
| GitHub Pages と Cloudflare が同じ production hostname を publish | deployment metadata と post-deploy hash mismatch | 一方を production から外す。同期 workaround は作らない |
| Consent Mode denied でも ping が出る | unknown/declined real/browser capture | Tag を load しない basic gating に戻す。解消不能なら disable |
| Withdrawal 後に beacon/request が出る | click timestamp 以降の pagehide/beacon capture | disable flag/reload boundary を修正。証明不能なら disable |
| Withdrawal storage が部分失敗し stale `allowed` が残る | operation-by-operation exception + read-back tests | verified non-allowed まで reload 禁止。Supported pathで解消不能なら disable |
| GA automatic fields/events が registry を越える | real transport URL/body/header/cookie inventory | property/config を縮小。承認・開示不能なら disable |
| User data が URL/title/referrer/error から漏れる | multi-surface sentinel test | P0。Transport を停止し root owner で修正 |
| Consent が session/history/reset に混入 | snapshot/session diff と import tests | P0。Persistence field を削除し owner を分離 |
| Event metric の denominator が観測不能 | Registry-to-KPI traceability audit | KPI/feature ID を削除。click proxy で代用しない |
| Scheduler churn を failure と数える | deterministic attempt-count tests | accepted attempt owner へ hook を移す |
| Analytics failure が app Promise/error に伝播 | blocked/timeout workflow matrix | fire-and-forget/no-op boundary を修正 |
| CSP allowlist が broad/unstable | captured endpoint と bundle CSP diff | endpoint を縮小。安定しなければ disable |
| GA script が COOP/COEP isolation と両立しない | exact hosted headers での loader capture | isolation/worker security を維持し analytics を disable |
| Legal/property approval が間に合わない | Phase 5 ledger blank | disabled profile を選ぶ |
| Browser test が CI で実行されない | workflow job evidence | Runner を追加するまで gate 未完了 |

Privacy P0 が見つかった場合、enabled branch の追加実装を続けて coverage を増やす前に
transport を default-disabled に戻す。

## 11. Scope cuts and non-goals

### 11.1 Cut order

Enabled analytics の工数が膨らむ場合は次の順で削る。

1. Feature ID 数を減らす。
2. `feature_progress(exported)` と format breakdown の一方を、KPI decision に不要なら減らす。
3. Per-feature analysis を止め、overall `render_result` reliability だけにする。
4. Custom event をすべて削り、analytics-disabled profile にする。

Consent correctness、withdrawal、local telemetry incapability、canonical artifact、payload
inspection、absence/offline tests は削らない。Analytics を切った場合は consent UI、
measurement ID、Google CSP origin、event/metric claims を shipped artifact/docs から外す。

### 11.2 Explicit non-goals

- Gallery/tutorial/Palette Explorer analytics
- session save/load event
- click、panel-open、configured-control analytics
- LOSAT program/presentation analytics。Validated analysis-job dispatch と bounded outcome
  の専用 owner がない v0.14 では post-analysis survivor を代理母数にしない
- record-count/sequence-size bucket
- raw mode/input type、record/feature identity parameter
- user ID、account、cross-device identity、custom session ID
- heatmap、replay recording、performance tracing、crash report upload
- server-side proxy、first-party analytics backend、Cloudflare Analytics への置換
- A/B test、feature gating、consent-dependent product behavior
- GA property の実変更、production deploy、external publicationを本計画作成だけで行うこと

Provider を将来変更する場合も consent/event contract を迂回せず、transport adapter と
observed payload/property gate を置き換える別計画にする。

## 12. Acceptance matrix

| Gate | Enabled profile oracle | Disabled/local oracle | Evidence owner |
| --- | --- | --- | --- |
| I-ART-01 canonical artifact | selected production path が exact bundle hash を配る | second production recipe がない | H1 provisional manifest + H-final/J final manifest + packaging/deployment record |
| I-ART-02 parser-time absence | passive config のみ、external GA script/inline config なし | measurement ID/Google CSP/config/consent markup なし。shared no-op module だけでは activate 不可 | packaging test |
| I-CONSENT-01 fail-closed | invalid/old/expired/storage-error は unknown、zero contact | forged allowed は hard no-op | Node + browser + offline verifier |
| I-CONSENT-02 lifecycle | allow/reload/re-allow/revoke/cross-tab/race/partial storage failure が transition table と一致し、unverified allowed は reload しない | consent UI/path なし | browser state matrix |
| I-CONSENT-03 isolation | request/session/history/Undo/Reset に consent がない | 同左 | session/history unit tests |
| I-NET-01 zero-contact | unknown/declined/revoked で script/request/cookie/dataLayer/backfill がゼロ。legacy GA cookie は無通信で削除 | 全状態でゼロ | browser interception + cookies |
| I-NET-02 allowed transport | approved endpoint/event/fields/cookies だけ | not applicable | staging capture + browser stub |
| I-DATA-01 sentinel | URL/body/header/referrer/cookie/dataLayer に sentinel なし | external attempts 自体ゼロ | browser capture |
| I-EVENT-01 schema | exact event/param/enum registry、unknown は reject | no effective event path | Node unit tests |
| I-EVENT-02 truth | success は atomic commit後、failure は genuine current attemptのみ | not applicable | run-analysis deterministic tests |
| I-EVENT-03 export | format completion後、committed snapshotのみ | not applicable | export unit/browser tests |
| I-KPI-01 derivability | event registry と `event_schema=1` reporting specification から3 KPI が一意に計算でき、synthetic truth table と一致 | analytics claim なし | property/query owner + traceability review |
| I-FAIL-01 independence | tag/collector block/timeoutでも全 core workflow pass | external-blocked full workflow pass | browser matrix |
| I-DOC-01 disclosure | actual custom/automatic data、cookie、withdrawal、ownerを記載 | analytics disabled/local absenceを記載 | docs tests + approvals |
| I-OPS-01 property | property config、access、retention、deletion、kill switch evidence | not applicable | operator record |
| I-RC-01 release choice | all enabled gates pass | all absence gates pass | Gate 0 + J traceability |

Enabled/disabled の双方を一つの曖昧な checkbox で pass させない。Gate 0 で選んだ branch
の行だけが release oracle になるが、I1 の local/disabled safety evidence は常に残す。

## 13. Evidence and implementation ledger

### 13.1 Phase ledger

| Phase | Status | Behaviour implemented | Evidence | Deviations | Remaining risk |
| --- | --- | --- | --- | --- | --- |
| 0 | pending |  |  |  |  |
| 1 | pending |  |  |  |  |
| 2 | pending/conditional |  |  |  |  |
| 3 | pending/conditional |  |  |  |  |
| 4 | pending/conditional |  |  |  |  |
| 5 | pending/conditional |  |  |  |  |
| 6 | pending |  |  |  |  |
| 7 | pending |  |  |  |  |

Skipped conditional phase は空欄にせず `not retained` と Gate 0 decision/evidence を記載する。

### 13.2 Release-choice record

```text
Decision: enabled | disabled
Decision date:
Candidate commit:
Canonical hosted builder:
Production deployment owner:
Measurement property owner (enabled only):
Notice version (enabled only):
Approvers:
Reason:
Invalidated evidence after later changes:
```

### 13.3 Artifact evidence record

各 profile について次を残す。

```text
Profile:
Qualification class: H1 provisional | H-final/J final
Source commit:
Build command:
Artifact path:
SHA-256:
Index/CSP/config summary:
Browser/runtime version:
Fixture and sentinels:
Consent scenario:
Observed script requests:
Observed collection requests:
Observed cookies before/after:
Observed pagehide/beacon:
App workflow oracle:
Result:
Log/capture path:
Reviewer/date:
```

### 13.4 Property and copy record

```text
Property/staging identifier:
Enhanced Measurement:
Google Signals:
Advertising/ad personalization:
Cross-domain linking:
User-provided data:
Data sharing:
Retention:
Access owners:
Registered event-scoped dimensions:
Saved report/Data API query owner and reference:
Property timezone and zero-denominator rule:
Synthetic KPI reconciliation evidence:
Deletion procedure:
Prior eager-deployment data disposition:
DNT/GPC and minor-user decision:
Kill-switch procedure:
Observed unavoidable metadata:
Privacy-details version:
Maintainer approval:
Privacy/legal approval:
Property-owner approval:
```

Measurement ID や account detail を repository secret として扱う必要がある場合、ledger
には secret 本文ではなく secret owner/reference と review evidence を記録する。

### 13.5 Test evidence rule

各 test record は command、commit、artifact hash、OS/browser/Python/Node version、
fixture、expected oracle、exit status、duration、unexpected skip/retry を含む。Static string
assertionだけを runtime network/cookie/payload evidence にしない。Manual staging capture
だけを deterministic consent reducer evidence にしない。両方が必要である。

## 14. Definition of done

Work package I は次の共通条件を満たしたときだけ完了する。

1. Eager GA script/config が source と canonical hosted builder から消えている。
2. Production hosted artifact と deployment owner が一つで、exact bundle hash を追跡できる。
3. Source/local/wheel/hosted-disabled artifact は forged storage を含む全状態で Google attempt/cookie がゼロである。
4. Consent は request/session/history/Undo/Reset から分離され、invalid storage は fail-closed である。
5. External analytics failure が core application workflow を一つも壊さない。
6. Focused、browser、offline、broad、artifact、diff-audit evidence が exact commit/hash に結び付いている。Work package I の handoff では H1 provisional hash を明示し、release qualification は H-final/J が final hash に対して同じ artifact/privacy gates を再実行するまで未完了とする。

さらに、Gate 0 の選択に応じて次のどちらか一方を満たす。

### Enabled completion

- allow 前/decline 後は script/request/cookie/dataLayer/backfill がゼロ。Revoke click から unload までは新規 request がゼロで cookie が消え、verified reload 後の new document は script/request/cookie/dataLayer がゼロ。
- withdrawal、partial persistence failure、pending load、cross-tab、expiry、reload の全 state matrix が passし、stored non-allowed を確認できない branch は reload しない。
- closed event/feature registry が accepted Result/export truth owner と一致し、KPI が導出可能。
- event-scoped dimensions と schema-filtered report/query が登録・所有され、synthetic KPI reconciliation が pass。
- real GA transport と property settings が approved/disclosed inventory に一致し、sentinel leak がゼロ。
- maintainer、privacy/legal reviewer、property operator の approval がある。

### Disabled completion

- Hosted artifact に measurement ID、Google CSP origin、consent notice、effective loader/event path がない。
- Public docs と release claims は analytics を使用すると記載しない。
- Enabled-only conditional phases は `not retained` decision と absence evidence を ledger に持つ。

この document の作成だけでは上記条件を一つも完了扱いにしない。実装時は phase ledger
を observed evidence に基づいて更新し、最後に Work package J の traceability table へ
handoff する。
