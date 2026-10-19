# Option Integrity Product Contract

Status: active Product authority

## Authority metadata

- Contract ID: `OIPC`
- Contract revision: `68`
- Product Decision Owner: `satoshikawato`
- Decision date: `2026-08-28`
- Decision source: explicit Product Decision Owner selection of one (`1`) after
  review of the merged deterministic evidence
- Reviewed candidate:
  `03_INITIAL_OPTION_INTEGRITY_PRODUCT_CONTRACT_CANDIDATE.md`
- Candidate SHA-256:
  `26b41219ca04ff26b56a29e11aa4be74c6030b0a22e386332afc28ac7a80623f`
- Reviewed evidence:
  [`COLLINEAR_MAX_CONFLICTS_COMPARISON_EVIDENCE_2026-08-28.md`](./COLLINEAR_MAX_CONFLICTS_COMPARISON_EVIDENCE_2026-08-28.md),
  merged by PR `#422` at `878c62ba17c61c45cd0adbd05cbd9fb36306db9d`
- Approved decision IDs: `PD-OI-001`, `PD-OI-002`, `PD-OI-003`,
  `PD-OI-004`, `PD-OI-005`, `PD-OI-006`, `PD-OI-007`, `PD-OI-008`,
  `PD-OI-009`, `PD-OI-010`, `PD-OI-011`, `PD-OI-012`, `PD-OI-013`,
  `PD-OI-014`, `PD-OI-015`, `PD-OI-016`, and `PD-OI-017`
- Initial candidate modification: `PD-OI-014`, as recorded below
- Revision 2 change: `PD-OI-007`, as recorded below
- Additional approved decision IDs: `PD-OI-018`, `PD-OI-019`, `PD-OI-020`,
  `PD-OI-021`, `PD-OI-022`, `PD-OI-023`, `PD-OI-024`, `PD-OI-025`,
  `PD-OI-026`, `PD-OI-027`, `PD-OI-028`, `PD-OI-029`, `PD-OI-030`,
  `PD-OI-031`, `PD-OI-032`, `PD-OI-033`, `PD-OI-034`, `PD-OI-035`,
  `PD-OI-036`, `PD-OI-037`, `PD-OI-038`, `PD-OI-039`, `PD-OI-040`,
  `PD-OI-041`, `PD-OI-042`, `PD-OI-043`, `PD-OI-044`, `PD-OI-045`,
  `PD-OI-046`, `PD-OI-047`, `PD-OI-048`, `PD-OI-049`, `PD-OI-050`,
  `PD-OI-051`, `PD-OI-052`, `PD-OI-053`, `PD-OI-054`, `PD-OI-055`,
  `PD-OI-056`, `PD-OI-057`, `PD-OI-058`, `PD-OI-059`, `PD-OI-060`,
  `PD-OI-061`, `PD-OI-062`, `PD-OI-063`, `PD-OI-064`, `PD-OI-065`,
  `PD-OI-066`, `PD-OI-067`, `PD-OI-068`, `PD-OI-069`, `PD-OI-070`,
  `PD-OI-071`, `PD-OI-072`, `PD-OI-073`, `PD-OI-074`, `PD-OI-075`,
  `PD-OI-076`, `PD-OI-077`, `PD-OI-078`, `PD-OI-079`, `PD-OI-080`,
  `PD-OI-081`, `PD-OI-082`, `PD-OI-083`, `PD-OI-084`, `PD-OI-085`,
  `PD-OI-086`, `PD-OI-087`, `PD-OI-088`, `PD-OI-089`, `PD-OI-090`, `PD-OI-091`,
  `PD-OI-092`, `PD-OI-093`, `PD-OI-094`, `PD-OI-095`, `PD-OI-096`,
  `PD-OI-097`, `PD-OI-098`, `PD-OI-099`, `PD-OI-100`, `PD-OI-101`,
  `PD-OI-102`, `PD-OI-103`, `PD-OI-104`, and `PD-OI-105`
- Revision 3 addition: `PD-OI-018`, accepted by `satoshikawato` on
  `2026-09-13` after confirming the complete record/search outcome, no feature
  retirement, and the runtime/memory cost of complete comparisons. The initial
  approval and its date above continue to describe `PD-OI-001`–`PD-OI-017`.
  The maintainer subsequently specified the fresh defaults in `PD-OI-019`.
  On the same date, the maintainer replaced the default-five member cap in
  `PD-OI-004` with an unbounded default and clarified `PD-OI-018` to require
  source-file execution, reusable raw evidence after display transforms,
  compact record disclosures, and execution controls in LOSAT Settings.
- Revision 4 addition: `PD-OI-020`, requested by `satoshikawato` on
  `2026-09-14` to restore automatic thread allocation when Total threads changes
  and protect that behavior against regression.
- Revision 5 changes: `PD-OI-001`, `PD-OI-002`, `PD-OI-004`, and `PD-OI-018`
  are replaced for the maintainer's `2026-09-14` LOSATP follow-up. New
  `PD-OI-021` records optional Collinear self-search and inference;
  `PD-OI-022` records reuse after cancellation. The final clarification requires
  each mode's edited limits to return when that mode is selected again, not a
  reset on every switch. The maintainer explicitly requested these decisions
  and contracts be recorded to prevent regression. Earlier approvals above
  retain their original scope; this amendment records the new instructions.
- Revision 6 addition: `PD-OI-023`, selected by `satoshikawato` on
  `2026-09-15` through the complete `PRODUCT_DECISION` response for
  `protein-comparison.path-representation`, scenario revision `1`.
  Only PATH-B and the supplied preservation, retirement and risk terms are
  recorded. The earlier decisions retain their scope. Dependent runtime still
  requires this authority on its base; this amendment contains no runtime.
- Revision 7 addition: `PD-OI-024`, selected by `satoshikawato` on
  `2026-09-19` through the complete `PRODUCT_DECISION` response for
  `linear.definition-display`, scenario revision `1`. Only D1-A, D2-P,
  D3-A and the supplied preservation, retirement and risk terms are recorded.
  Earlier decisions retain their scope. Dependent runtime requires this
  authority merged into its base; this amendment contains no runtime.
- Revision 8 change: `PD-OI-018` is replaced for scenario revision `3`,
  selected by `satoshikawato` on `2026-09-20` through the complete
  `PRODUCT_DECISION` response for
  `diagram-generation.linear-record-universe-and-search-scope`. The selected
  `LINEAR-FILE-ROW-BLOCK` outcome makes normal-layout File-card order and
  visual row order one operation, blocks File-card moves for custom layouts,
  and records the supplied preservation, retirement, and risk terms. Earlier
  decisions retain their scope. Dependent runtime requires this authority
  merged into its base; this amendment contains no runtime.
- Revision 9 addition: `PD-OI-025`, selected by `satoshikawato` on
  `2026-09-21` through the complete `PRODUCT_DECISION` response for
  `diagram-generation.linear-depth-source-scope-and-discoverability`, scenario
  revision `1`. The selected `FILE-BULK-WITH-RECORD-OVERRIDES` outcome exposes
  common Depth TSV assignment on each Linear File card while preserving sparse
  per-record bindings. Earlier decisions retain their scope. Dependent runtime
  requires this authority merged into its base; this amendment contains no
  runtime.
- Revision 10 additions: `PD-OI-026` through `PD-OI-031`, selected by
  `satoshikawato` on `2026-09-22` through six complete `PRODUCT_DECISION`
  responses for issue `#561`. These additions record only the supplied choices,
  preservation requirements, retirement permissions, and accepted residual
  risks for deterministic Similarity Group alignment. Earlier decisions retain
  their scope. Dependent runtime requires this authority merged into its base;
  this amendment contains no runtime.
- Revision 11 addition: `PD-OI-032`, selected by `satoshikawato` on
  `2026-09-22` through the complete `PRODUCT_DECISION` response for issue
  `#563`. This addition records the selected feature-popup record-rotation
  outcome, preservation requirements, lack of retirement permission, and
  accepted residual risk. Earlier decisions retain their scope. Dependent
  runtime requires this authority merged into its base; this amendment contains
  no runtime.
- Revision 12 additions: `PD-OI-033` through `PD-OI-035`, approved by
  `satoshikawato` on `2026-09-24` through explicit approval of the exact three
  `PRODUCT_DECISION` texts presented for issue `#581`. The receipts below are
  the full scope of these additions. Dependent runtime requires this authority
  merged into its base; this amendment contains no runtime.
- Revision 13 changes: `PD-OI-026`, `PD-OI-027`, `PD-OI-031`, and `PD-OI-034`
  are replaced for scenario revision `2`, selected by `satoshikawato` on
  `2026-09-24` through four complete `PRODUCT_DECISION` responses for issue
  `#586`. The receipts below are the full scope of these replacements. Earlier
  unaffected decisions retain their scope. Dependent runtime requires this
  authority merged into its base; this amendment contains no runtime.
- Revision 14 changes: `PD-OI-026`, `PD-OI-027`, `PD-OI-031`, and `PD-OI-034`
  are replaced for scenario revision `3`. `satoshikawato` explicitly
  approved the four complete Choice A `PRODUCT_DECISION` texts on
  `2026-09-25` for the issue `#586` follow-up. The receipts below are the
  full scope of these replacements. Other decisions retain their scope.
  Dependent runtime requires this authority merged into its base; this
  amendment contains no runtime.
- Revision 15 change: `PD-OI-035` is replaced for scenario revision `2`.
  `satoshikawato` supplied the complete `A / RETAIN_MOBILE_PALETTE_COVERAGE`
  `PRODUCT_DECISION` response on `2026-09-25`. The receipt below is the full
  scope of this replacement. Other decisions retain their scope. This
  amendment contains no runtime.
- Revision 16 changes: `PD-OI-027`, `PD-OI-028`, `PD-OI-029`, `PD-OI-031`,
  and `PD-OI-034` are replaced for scenario revisions `4`, `2`, `2`, `4`,
  and `4`. `satoshikawato` explicitly approved the five complete
  `PRODUCT_DECISION` texts as written on `2026-09-26` for the record-owned
  orientation follow-up to issue `#586`. The receipts below are the full
  scope of these replacements. `PD-OI-026`, `PD-OI-030`, `PD-OI-035`, and
  other decisions retain their scope. Dependent runtime requires this
  authority merged into its base; this amendment contains no runtime.
- Revision 17 changes: `satoshikawato` signed five complete
  `PRODUCT_DECISION` responses for issue `#602` on `2026-09-26`:
  01=B and 02–05=A. `PD-OI-024` scenario revision `2` replaces only the
  D1 Web fresh/reset default and retains D2-P/D3-A. New `PD-OI-036` through
  `PD-OI-039` serialize the other four receipts. `PD-OI-035` scenario
  revision `3` retains its independent canvas-interaction guarantees and
  supersedes only the scenario-2 mobile coverage exception and close-review
  continuation, using the signed review-presentation receipt in `PD-OI-039`.
  The five receipts below preserve all nine supplied fields exactly. Earlier
  unaffected outcomes retain their scope. This amendment contains no runtime;
  dependent runtime requires all five outcomes and the limited supersession
  merged into its base.
- Revision 18 additions: `PD-OI-040` through `PD-OI-043`, approved by
  `satoshikawato` on `2026-09-26` as the four complete Choice A outcomes for
  issue `#600`. The approval receipt is 「すべて推奨案で承認します。」 and
  its complete scope is retained in
  [`APPROVED_PRODUCT_DECISIONS.md`](./issue-600-implementation-20260926/APPROVED_PRODUCT_DECISIONS.md).
  The four independent records below preserve every supplied receipt field
  and outcome without additional retirement or risk terms. Earlier decisions
  retain their scope. This amendment contains no runtime; dependent runtime
  requires these records merged into its base.
- Revision 19 changes: `PD-OI-027`, `PD-OI-029`, `PD-OI-031`, and
  `PD-OI-034` are replaced for scenario revisions `5`, `3`, `5`, and `5`.
  `satoshikawato` explicitly approved the four complete issue `#598` receipts
  on `2026-09-26`. They define exclusive displayed-direction modes and
  selectable Reset scope. `PD-OI-026`, `PD-OI-028`, `PD-OI-030`,
  `PD-OI-032`, `PD-OI-033`, `PD-OI-035`, and other decisions retain their
  scope. This amendment contains no runtime; dependent implementation must
  have this authority merged into its base. Issue `#598` BUG-17 is excluded
  from implementation at the owner's instruction, with no circular-rotation
  authority change.
- Revision 20 change: `PD-OI-039` is replaced for scenario revision `2`.
  `satoshikawato` explicitly approved the complete nine-field
  `A / EXCLUSIVE_DIRECTIONS_WITHOUT_MATCH` receipt on `2026-09-26`
  (UTF-8 receipt SHA-256
  `06a9d2fe9b1d1406f6f8e04c23a9ca031133b9fae24b0683403ec2c6cae55270`).
  Only the old Match affordance is additionally retired to reconcile issue
  `#598` with issue `#602`. All independent `PD-OI-035`/`PD-OI-039`
  requirements and all four accepted issue `#598` decisions retain their
  scope. This authority-only amendment contains no runtime.
- Revision 21 additions: `PD-OI-044` and `PD-OI-045`, approved by
  `satoshikawato` on `2026-09-26` as the two complete Choice A outcomes for
  issue `#597` BUG-02 and BUG-20. The approval receipt is
  「推奨案で承認します」 for each outcome. The complete supplied receipts
  are serialized below without inferred rationale, preservation, retirement,
  or risk terms. BUG-01 is outside this delivery scope at the owner's
  instruction and receives no authority record. Earlier outcomes retain
  their scope. This amendment contains no runtime; dependent implementation
  requires this authority merged into its base.
- Revision 22 additions: `PD-OI-046` and `PD-OI-047`, selected by
  `satoshikawato` on `2026-09-26` through the two complete Choice A
  `PRODUCT_DECISION` receipts for issue `#601` BUG-15 and BUG-19.
  The two independent records below preserve all nine supplied fields and
  their original concern keys exactly. Earlier decisions and acceptance
  conditions retain their scope. This authority-only amendment contains no
  runtime; dependent implementation requires these records merged into its
  base. It neither registers nor supersedes the separate export-plan
  diagnostic-disclosure candidate.
- Revision 23 additions: `PD-OI-048` through `PD-OI-050`, approved by
  `satoshikawato` on `2026-09-27` as the three complete Choice A receipts for
  issue `#619`. The owner explicitly answered 「署名します。」 to confirmation
  of all three Choice A texts, Owner `satoshikawato`, and Decision date
  `2026-09-27`. The independent records below preserve all nine supplied
  fields. Earlier decisions and acceptance conditions retain their scope.
  This authority-only amendment contains no runtime or runtime acceptance
  evidence; dependent implementation requires these records merged into its
  base.
- Revision 24 addition: `PD-OI-051`, the complete Choice A receipt for
  `diagram-generation.inflight-comparison-draft`, approved and signed by
  `satoshikawato` on `2026-09-27`. The record and source receipt preserve all
  nine approved fields. Earlier decisions and acceptance conditions retain
  their scope. This authority-only amendment contains no runtime or runtime
  acceptance evidence; dependent implementation requires the record merged
  into its base.
- Revision 25 additions: `PD-OI-052` through `PD-OI-054`, approved by
  `satoshikawato` on `2026-09-26` as the three complete independent Choice A
  outcomes for issue `#599`. The approval receipt is
  「すべて推奨案で承認します。」; each published Decision Pack and its exact
  receipt are identified below. These additions preserve all supplied fields
  without extending rationale, retirement, or risk. Earlier outcomes and
  their independent preservation requirements retain their scope. This
  authority-only amendment contains no runtime or runtime acceptance evidence;
  dependent implementation requires these records merged into its base.
- Revision 26 addition: `PD-OI-055`, selected as
  `A / RETAIN_VALIDATED_BINDING_ENRICHMENT` by `satoshikawato` on
  `2026-09-28` through explicit approval of the complete nine-field
  `PRODUCT_DECISION` text below for Issue `#619` finding 5. The receipt and
  JSON preserve only that outcome and its supplied preservation, retirement,
  and risk terms. Earlier decisions retain their scope. This authority-only
  amendment contains no runtime or runtime acceptance evidence; dependent
  runtime requires this authority merged into its base.
- Revision 27 changes: `PD-OI-037` and `PD-OI-049` are replaced for scenario
  revision `2`.
  - `PD-OI-037`: the Product Decision Owner `satoshikawato` approved the
    receipt text on `2026-09-29` (GUI remediation S00 decision 4), covering
    Rationale and Accepted residual risk. Choice, Must preserve, and May
    retire come from the Owner's confirmed requirement to remove the always-on
    Pending display and the named always-on explanations.
  - `PD-OI-049`: on the same date, the Owner approved a help-only mitigation
    in place of help and Pending. The other fields are unchanged.
  - Both revisions use the reviewed static Product Contract co-change route
    and merge with their implementation. The lifecycle text below records
    that route.
- Revision 28 changes: `PD-OI-024` is replaced for scenario revision `3` and
  `PD-OI-054` for scenario revision `2`, through the same co-change route.
  - Both use the receipt text that `satoshikawato` approved on `2026-09-29`
    (GUI remediation S00 decision 4) for Rationale and Accepted residual risk.
  - `PD-OI-024`: the always-on Lock explanation moves into the existing
    help-tip and the checkbox's accessible description. D1-B, D2-P, and D3-A
    are otherwise unchanged.
  - `PD-OI-054`: the approved additional scope places search and toolbar in
    the width left by the Editor, which opens from the Preview top edge and
    covers neither.
- Revision 29 changes: `PD-OI-056` through `PD-OI-084` are added and
  `PD-OI-018` is replaced for scenario revision `4`, from the 30 receipts of the
  2026-09-30 Web GUI audit Decision Pack
  ([`02_DECISION_PACK.md`](./web-gui-audit-20260930/02_DECISION_PACK.md)).
  `OIC-027` is added for `PD-OI-066`.
  - The Product Decision Owner `satoshikawato` answered on `2026-09-30`, first
    to `D-01`–`D-39` and `W-1`–`W-8`: 「だいたいそのままで承認。けど、D-12: 1ファイルに複数生物あるときは、Definitionはそれぞれにつけてくれるとうれしいかも。Circularでしょ？ D-19: 今まではBだったんじゃないの？たとえばmultiple repliconのゲノムだったら、染色体ごとじゃなくてゲノム=ファイル単位のE-valueが欲しいんじゃないかな？ D-26: UはOKにして。」
    and then to the rewritten `D-12`, `D-19`, and `D-26` receipts and the `D-40`
    CLI scope: 「OK,これで承認します。」
  - Receipt mapping: `D-01` → `PD-OI-056`; `D-02` → `PD-OI-057`; `D-03` →
    `PD-OI-058`; `D-04` → `PD-OI-059`; `D-05` → `PD-OI-060`; `D-06` →
    `PD-OI-061`; `D-07` → `PD-OI-062`; `D-08` → `PD-OI-063`; `D-09` →
    `PD-OI-064`; `D-10` → `PD-OI-065`; `D-11` → `PD-OI-066`; `D-12` →
    `PD-OI-067`; `D-13` → `PD-OI-068`; `D-14` → `PD-OI-069`; `D-15` →
    `PD-OI-070`; `D-16` → `PD-OI-071`; `D-17` → `PD-OI-072`; `D-18` →
    `PD-OI-073`; `D-20` → `PD-OI-074`; `D-21` → `PD-OI-075`; `D-22` →
    `PD-OI-076`; `D-23` → `PD-OI-077`; `D-24` → `PD-OI-078`; `D-25` →
    `PD-OI-079`; `D-26` → `PD-OI-080`; `D-27` → `PD-OI-081`; `D-28` →
    `PD-OI-082`; `D-29` → `PD-OI-083`; `D-30` → `PD-OI-084`; `D-19` →
    `PD-OI-018` scenario revision `4`.
  - Each record reproduces its receipt without translation or additional terms.
    `D-31`–`D-40` keep current behavior and receive no record. Earlier
    decisions retain their scope. This authority-only amendment contains no
    runtime; dependent implementation requires it merged into its base.
- Revision 30 changes: `PD-OI-033` is replaced for scenario revision `2` and
  `PD-OI-085` is added, from the Product Decision Owner `satoshikawato`'s
  answers of `2026-10-04` recorded in section 1 of
  [`POPUP_RECORD_ACTIONS_AND_VIBRIO_SESSION_PLAN_2026-10-04.md`](./POPUP_RECORD_ACTIONS_AND_VIBRIO_SESSION_PLAN_2026-10-04.md)
  and its Appendix A.
  - The Owner answered the plan's open questions: `OD-1` (move Record actions
    from the top of the Edit tab into a Layout group): 「移しましょう。」
    `OD-2` (name the staging button **Apply on Generate**) and `OD-3` (keep
    **Apply and regenerate** target-only so staged changes on other records
    stay pending): 「推奨通りでお願いします。」 `OD-4` (proceed):
    「計画書をまとめた後、実装に移ってください。」
  - `PD-OI-033` scenario revision `2` (`B / LAYOUT_GROUP_DISCLOSURE`) follows
    from `OD-1`. `PD-OI-085` (`A / APPLY_ON_GENERATE`) follows from `OD-2` and
    `OD-3`.
  - The receipt fields restate those answers and the plan that was presented
    in the same session. The Owner did not separately review the receipt
    wording; the Owner's approval covers the answers quoted above and the
    plan. Earlier decisions, including `PD-OI-032`, retain their scope. This
    authority-only amendment contains no runtime; dependent runtime requires
    it merged into its base.
- Revision 31 changes: `PD-OI-018` is replaced for scenario revision `5`, from
  the Product Decision Owner `satoshikawato`'s answer of `2026-10-04`, given in
  that session through a multiple-choice question.
  - The question: 「9/20 に承認した PD-OI-018 には「1つの File が複数の行にまたがるときは File を移動できない（custom layout 扱い）」と書かれています。CLI で書いたセッション（1レコード＝1行）を Load すると、このルールのため警告が出て並べ替えもできません。R9（#780）はルールを「各 File のレコードが連続した行にまとまり、別の File と行を共有していなければ、その行のまとまりごと移動できる」に広げる変更です。どうしますか？」
  - The Owner selected 「ルールを広げる (Recommended)」, whose option text was: 「PD-OI-018 を revision 5 に改訂する PR を先に出し、そのあと #780 をマージします。CLI セッションでも警告が出ず、File を並べ替えられるようになります。複数の File が1行を共有するレイアウトや、File の行が交互に入り組んだレイアウトは、今までどおり custom 扱いです。」
  - `PD-OI-018` scenario revision `5` (`C / CONSECUTIVE-FILE-ROW-BLOCKS`)
    follows from that answer. Pull request `#780`
    (`fix/web-linear-file-block-moves`) implements it.
  - The receipt fields restate that answer and the option text. The Owner did
    not separately review the receipt wording; the Owner's approval covers the
    answer quoted above. Earlier decisions retain their scope. This
    authority-only amendment contains no runtime; dependent runtime (`#780`)
    requires it merged into its base.
- Revision 32 changes: `PD-OI-061` is replaced for scenario revision `2`
  (`B / MERGE-SAME-FEATURE-TYPE-ONLY`), from the Product Decision Owner
  `satoshikawato`'s replies of `2026-10-06` while reviewing OV-62 (a Legend row
  without features renamed onto a drawn row's caption).
  - Reply 1, on the screenshot of the Merge, Suffix, and Cancel dialog for
    GC content renamed onto `CDS`: 「それはマズいね。CDSとGC contentは全く別物じゃん。トラックとか、同じフィーチャートラックでも別タイプのフィーチャーの場合は混ぜちゃいけないよね基本的に。」
  - Reply 2, selecting an option of the multiple-choice question
    「Merge を出してよい条件」: 「同じ track の同じ feature type 同士だけ (Recommended)」
  - Reply 3, on the full receipt text recorded below: 「OKです」
  - The Owner approved the receipt text as written. Earlier decisions retain
    their scope. This change is a static Product Contract co-change: the
    runtime, tests, and documentation that implement it are in the same pull
    request, and the Review is `REQUIRED`.
- Revision 33 changes: `PD-OI-086` is added
  (`A / ALL-DIAGRAM-SETTINGS-PER-MODE`), from the Product Decision Owner
  `satoshikawato`'s reply of `2026-10-07` for OV-80, OV-82, and related
  findings (settings and edits made in one diagram mode changed or failed the
  other mode's diagram and Generate).
  - The question, asked with the full receipt text recorded below: 「上の PD-OI-086 の文面を、このまま Product Contract に記録してよいですか？（Contract だけの PR として先に入れます。Review が必要で、自動 merge はしません）」
  - The Owner's reply: 「OK、この文面で」
  - The Owner approved the receipt text as written. The receipt text
    incorporates the Owner's earlier answers of the same date.
  - No other record changes. `PD-OI-061`, `PD-OI-062`, `PD-OI-063`, and
    `PD-OI-084` retain their scope within each mode. `PD-OI-063` names no
    mode and is read per mode: each mode carries its own Legend order, as
    `PD-OI-052` already matches decoration deltas by mode. `PD-OI-002`
    retains its scope within Linear. `PD-OI-044`, `PD-OI-066`, and
    `PD-OI-070` are cited by the receipt and retain their scope.
  - Outside this Contract, the receipt retires the scope of GUI remediation
    S00 decision 2 of `2026-09-29`
    ([`S00.md`](./gui-remediation-20260928/results/S00.md) section 11.4),
    which made only `plot_title`, `plot_title_font_size`, and `def_font_size`
    per mode. The receipt keeps the rule of S00 decisions 1 and 2 for old
    Session values that were per mode.
  - Earlier decisions retain their scope. This authority-only amendment
    contains no runtime; dependent runtime (the OV-80 and OV-82
    implementation) requires it merged into its base, and the Review is
    `REQUIRED`.
- Revision 34 changes: `PD-OI-087` is added
  (`A / POPUP_ONLY_FEATURE_ROTATION`), from the Product Decision Owner
  `satoshikawato`'s replies of `2026-10-08` for OV-212 (two routes set a
  record's display start from a feature: the sidebar buttons need a Ctrl/⌘/Shift
  selection, and their 5′ end splits a − strand feature across both ends of the
  record).
  - Reply 1, selecting option A of the OV-212 question: 「A（推奨） 2ボタンを廃止し、遺伝子を起点にした回転はポップアップに一本化する。数値欄と Reset start は残す。にします」
  - Reply 2, on the proposed Rationale and Accepted residual risk wording
    recorded below: 「案のとおり」
  - Must preserve and May retire restate option A. The Owner adopted the
    Rationale and Accepted residual risk as written.
  - The receipt retires two controls that earlier receipts preserve as sidebar
    operations: 「既存 sidebar 操作」 of `PD-OI-033` scenario revision `2`,
    「sidebar の record display 操作とその意味」 of `PD-OI-085`, and
    「既存sidebar workflow」 of `PD-OI-032`. It narrows those items by the two
    buttons only; their other items and the other fields of those records
    retain their scope.
  - This change is a static Product Contract co-change: the runtime, tests, and
    documentation that implement it are in the same pull request, and the
    Review is `REQUIRED`.
- Revision 35 changes: `PD-OI-088` is added
  (`B / DIALOG_FIRST_BUSY_CHOICE`), with the acceptance contract `OIC-028`,
  from the Product Decision Owner `satoshikawato`'s replies of `2026-10-09`
  for OV-225 (on a Session just loaded, a popup color pick waited about 6 s
  for the Python runtime before its scope dialog opened, with no visible
  status).
  - Reply 1, selecting the recommended outcome of the OV-225 question (keep
    the current order in the 0.14.0 performance pull request, and option (b),
    dialog first with the dialog busy until its choice commits, as its own
    change): 「OV-225 のダイアログの順序: 推奨:どおり」
  - Reply 2, on the complete proposed receipt recorded below:
    「OV-225: これでいいです」
  - Reply 3 (`2026-10-09`), to the question whether the receipt covers the
    same dialogs opened from the Features drawer color input, popup stroke
    edits, and the Legend panel rename: 「F4: A」. The receipt and its SHA-256
    are unchanged; `OIC-028` names these entrances.
  - The receipt narrows no earlier record. It realizes, for these three
    dialogs, the 「live edit適用中/失敗の通知」 that `PD-OI-037` scenario
    revision `2` preserves.
  - This change is a static Product Contract co-change: the runtime, tests, and
    documentation that implement it are in the same pull request, and the
    Review is `REQUIRED`.
- Revision 36 changes: `PD-OI-089` is added
  (`A / USER_DEFAULT_COLORS_OVER_PALETTE`), with the acceptance contract
  `OIC-029`, from the Product Decision Owner `satoshikawato`'s replies of
  `2026-10-09` on palette changes and edited Default colors (a palette switch
  or the Default colors Reset discarded the colors a user had edited without
  asking, and popup Apply to all on a Legend row that no Specific color rule
  draws wrote one hash rule per feature, so a feature that joined the row
  later kept the palette color and formed a second row).
  - Reply 1, the proposal (a palette change adds `-d` entries, and Apply to
    all or a palette color that would overwrite an existing `-d` entry asks in
    a dialog): 「-pを変えるというのは、-dにエントリを追加するという挙動として再定義するのがいいんじゃないかな？そのうえで、Apply to all CDSでもパレットから直接色を変えるのでも、それが既存の-dのエントリを上書きするものであれば、「設定上書き/消えちゃうけどどうしますか？」っていう選択肢のダイアログのポップアップを出すのがいいんじゃないかな？」
  - Reply 2, accepting the refined direction (a user default color is a
    Default colors value that differs from the selected palette's color;
    Apply to all on such a row sets the default color; a palette switch and
    the Default colors Reset ask before they discard user colors), in v0.14.0
    as a separate pull request: 「OK、その方向で行きましょう。v0.14.0に入れる。別のPRにする。」
  - Reply 3, on the open outcomes Q1-Q6, the complete proposed receipt
    recorded below, and F4 (recorded in the Revision 35 entry):
    「Q1=B Q2=B Q3=B Q4=A Q5=A Q6=A」「受領書: これでいい」「F4: A」
  - Reply 4, withdrawing Q3 because its premise did not hold, and on the
    Q2 line wording: 「Q3: A'」「Q2 文言: 推奨どおり」
  - Q1=B: Apply to all on a palette row while a palette is queued sets the
    type's color in the queued palette's colors and in the shown Result's
    colors, so the Result shows it now.
  - Q2=B: the scope dialog always shows one line under Apply to all on a
    palette row: "Sets the <type> default color (<type> features without
    their own color or rule, also hidden ones)." or, when the type has a user
    default color, "Replaces your <type> default color ■ #hex." (Reply 4
    replaced the approved "all <type> features", which overstated).
  - Q3 is withdrawn: it assumed a Legend-tab fill color control for a palette
    row, and the GUI has none, so nothing is built and `OIC-029` has no Q3
    term; the receipt is recorded exactly as approved.
  - Q4=A: Specific color rule rows that older Sessions saved for Apply to all
    stay rule rows; nothing migrates them.
  - Q5=A: switching back to the applied palette while another is queued asks
    like any switch when user colors exist; the queued palette is dropped and
    the kept colors apply now.
  - Q6=A: dialog bodies name up to three keys.
  - The receipt retires no item of an earlier record. The records its Must
    preserve cites (`PD-OI-037`, `PD-OI-066`, `PD-OI-069`, `PD-OI-070`,
    `PD-OI-086`, and `PD-OI-088`) retain their scope.
  - This change is a static Product Contract co-change: the runtime, tests, and
    documentation that implement it are in the same pull request, and the
    Review is `REQUIRED`.
- Revision 37 changes: `PD-OI-090` is added
  (`A / PROJECT-WITH-DRAWINGS`, concern `session.project-drawings`), with the
  acceptance contract `OIC-030`, and replaces `PD-OI-086` in full;
  `PD-OI-044` is replaced for scenario revision `2`
  (`A / REVEAL_APPLICABLE_SINGLE_RECORD_CONTROLS`). Both come from the Product
  Decision Owner `satoshikawato`'s replies of `2026-10-07` and `2026-10-08`
  for the v0.15.0 multi-drawing project (one Session holds shared input files
  and several drawings, each with its own settings, edits, History, and
  Result).
  - Replies (question text abbreviated, option labels as selected): R1-1
    「Drawing tabは基本的にまっさらなものを新規作成がいいのでは。スタイルとか入力をコピーする場合はそういうexplicitな選択肢があるといいと思う。新規作成、複製、スタイルをコピー、みたいな複数の選択肢。そのうえで、たとえば既存のdrawingから特定のエントリーとか特定の領域を切り出して新規drawingにする場合は、スタイルを引き継ぐ（数mbのレコードから数十kbを切り出す場合は、フィーチャーの縦幅とかラベルのフォントサイズとかストロークの太さとかGC content/skeｗの幅とかを適宜気を利かせて変える必要があると思う）。」;
    R1-2 「使った mode だけ (推奨)」; R1-3 「Result の drawing だけ」;
    R1-4 「起動を認める (推奨)」; R2-1 「中身は1回保存、選択はdrawingごと (推奨)」;
    R2-2 「その都度聞く (推奨)」; R2-3 「gbdraw render を追加 (推奨)」;
    R2-4 「Contract だけの PR を先に (推奨)」; R3-1 「上部に drawing バー (推奨)」;
    R3-2 「この 4 種類 (推奨)」; R3-3 「確認＋Undo で戻せる (推奨)」;
    R4-1 「それは完全にregressionじゃん。ダメだよこれは。Session fileはすべてのResultを保存して。」;
    R4-3 「切り出しも最初に」; R6-1 「v0.14.0 に入れる (推奨)」;
    R7-1 「その record の見た目すべて (推奨)」; R7-2 「Auto が変わる値だけ戻す (推奨)」;
    R7-3 「最初の版では扱わない (推奨)」; R11-1 「ボタン 3 つ (推奨)」 (Reset this
    drawing / Reset all drawings / Cancel; replaces the checkbox wording of
    R3-4 「表示中の drawing＋他も選べる (推奨)」); and the Owner's replies of
    `2026-10-07` on the Reset dialog 「このモードだけか全部か選択可能にする」 and
    「「現在の描画消えるけどいいですか？」ていうダイアログをつけて、すぐにリセット後の見た目にする（即時反映する）か、キャンセルする選択肢を与える」.
  - The question asked with the full receipt texts recorded below:
    「AUTHORITY-15 の受領書 2 通を承認しますか？ 対象は PD-OI-087 と PD-OI-044 rev 2。本文 owner/receipt-draft.r2.txt、owner/receipt-draft-pd-oi-044-rev2.r2.txt」
    The Owner's reply (`2026-10-08`) selected 「承認 (推奨)」, whose option
    text reads 「決定日は 2026-10-08。正式な記録（SHA と JSON）を作って見せる。Contract の PR は 0.14.0 のタグの後。」
    The Owner approved the receipt texts as written.
  - The question named the drawings receipt by its draft ID `PD-OI-087`, and
    its draft acceptance contract was `OIC-028`. Revisions 34 and 35 took
    those IDs first, so they are recorded as `PD-OI-090` and `OIC-030`. The
    `PD-OI-090` receipt names no ID of its own and is recorded byte for byte
    as approved. The `PD-OI-044` receipt names the drawings record twice; both
    read `PD-OI-090` in place of `PD-OI-087`, no other character differs from
    the approved text, and its SHA-256 below is that of the renumbered text.
  - `PD-OI-086` is replaced in full: its record is removed from the active
    Contract (Git history retains it) because `PD-OI-090` restates each of
    its terms per drawing and changes four of them: one History in which a
    mode switch is a step, Reset Settings resetting both modes, Session 27-44
    Load without Python, and the limit of two diagrams. GUI remediation S00
    decision 2
    ([`S00.md`](./gui-remediation-20260928/results/S00.md) section 11.4) stays
    retired; the rule of S00 decisions 1 and 2 for old Session values stays as
    the migration rule in `PD-OI-090`.
  - `PD-OI-044` changes only its Python clause: the preview-only Load without
    Python holds for Sessions of version 46 and later; Sessions of version 27
    to 44 are converted by Python once on Load.
  - `PD-OI-089` cites `PD-OI-086` in its Must preserve for separate palettes
    and colors in Circular and Linear. That item is read per drawing:
    `PD-OI-090` keeps the palette and color rules of each drawing separately.
    Its 「Session の形式」 is the Session format that `PD-OI-089` itself does
    not change; the format change of `PD-OI-090` is outside that item.
  - Persisted-format evidence: version 46 is the first-parent `dev` writer
    since `acd35a17` (PR `#945`), which release `0.14.0` ships, and this
    amendment merges only after the `0.14.0` tag; version 44 is the
    first-parent `main` writer `fe6861f0`; release `0.13.0` writes 30.
    Version 45 was development-only and is not read.
  - Cited records retain their scope and are read per drawing: `PD-OI-002`,
    `PD-OI-004` (LOSATP limits within a Linear drawing), `PD-OI-009`,
    `PD-OI-016`, `PD-OI-018`, `PD-OI-019`, `PD-OI-020`, `PD-OI-027`,
    `PD-OI-030`, `PD-OI-032`, `PD-OI-037`, `PD-OI-038`, `PD-OI-045`,
    `PD-OI-046`, `PD-OI-052`, `PD-OI-061` to `PD-OI-066`, `PD-OI-070`,
    `PD-OI-073`, `PD-OI-079`, and `PD-OI-089`. Earlier decisions retain their
    scope. This authority-only amendment contains no runtime; dependent
    runtime requires it merged into its base, and the Review is `REQUIRED`.
- Revision 38 changes: `PD-OI-091` is added
  (`A / RECORD-DRAW-ON-OFF`, concern `diagram-generation.record-draw-selection`).
  It comes from the Product Decision Owner `satoshikawato`'s replies of
  `2026-10-10` for drawing a chosen subset of a file's records in the Web app.
  - Replies (question text abbreviated, option labels as selected): D-01
    「特定のレコードだけ描画する機能をv0.14.0に絶対に入れなければならない。例えば100本以上コンティグがあるドラフトアセンブリとかだったら絶対に必要。」;
    D-02 「個々のレコードのプルタブを開くとチェックを外して/入れて描画OFF/ONできるようにするほかに、レコードのリストをポップアップしてそこからチェックを外して/入れて描画OFF/ONできるようにする2ルートが必要」;
    D-03 「範囲：Circular の Multi-Record Canvas の部分選択も 0.14.0 に入れる。」;
    D-04 「アップロード直後の既定： レコード数が閾値（例：20）を超えたらピッカーを自動で開く」;
    D-05 「外したレコードのカード設定: 外しても保持しよう。オプションで削除できるようにしよう。」 and
    「Session には基本保存するけど、さっきの回答のように、完全削除するオプションもつける」;
    D-06 「これはポップアップで質問して「ファイルごと消しますか？」って出すのがいいと思う」;
    D-07 「P-1：推奨どおり」 (the LOSAT E-value database stays the whole subject source file);
    D-08 「P-4: 推奨どおり」 (Delete settings also removes the record's feature edits, Feature placements, annotations that target it, and comparison pairs that name it);
    D-09 「残りも推奨通り。」;
    D-10 「一覧ではレコードのソートもできるようにしたい。IDとか、長さ順とか、オリジナルの順に戻すとか」.
  - The question asked with the full receipt text recorded below was answered
    by the Owner on `2026-10-10` with 「文面承認」. The Owner approved the
    receipt text as written; it is recorded byte for byte.
  - `PD-OI-010` accepted residual risk "Arbitrary multi-record subset editing
    remains unavailable" and the exclusion of record-subset editing in
    `PD-OI-044` scenario revision `2` `May retire` are retired by the
    `PD-OI-091` `May retire`. Their texts are unchanged; the later record
    controls where they differ, and Git history retains the approved receipts.
  - `PD-OI-018` is not amended: the E-value database stays the whole subject
    source file. `PD-OI-019`, `PD-OI-044`, `PD-OI-070`, and `PD-OI-090` keep
    their scope and are read together with `PD-OI-091`. Session version `46`
    is unchanged; the new field is optional.
  - This change is a static Product Contract co-change: the runtime, tests, and
    documentation that implement it are in the same pull request, and the
    Review is `REQUIRED`.
- Revision 39 changes: `PD-OI-028` gains scenario revision `3`
  (`A / KEEP_ANCHOR_PLAN`, regrouping comparison settings) from the Product
  Decision Owner `satoshikawato`'s replies of `2026-10-11` for OV-341 (after a
  Result filter such as **Minimum identity** regrouped the plan's Similarity
  Group, every later Generate was blocked with a stale-reference repair).
  - Reply 1, selecting the recommended outcome: 「A（推奨）整列を保つ」.
  - Reply 2, on the complete receipt recorded below: 「下書きの文面を承認する」.
    The receipt is recorded as approved.
  - The scenario-revision-`2` receipt keeps its full scope; scenario revision
    `3` adds the regrouping outcome and retires nothing.
  - This change is a static Product Contract co-change: the runtime, tests, and
    documentation that implement it are in the same pull request, and the
    Review is `REQUIRED`.
- Revision 40 changes: `PD-OI-085` gains scenario revision `2`
  (`RENAME_ROTATION_APPLY_BUTTONS`), `PD-OI-089` gains scenario revision `2`,
  and `PD-OI-061` gains scenario revision `3`. Each renames visible labels
  and changes no behavior. They come from the Product Decision Owner
  `satoshikawato`'s replies of `2026-10-11` on the redesigned feature popup
  and choice dialogs of the Web app.
  - Reply 1, to the review questions on the Before/After pairs of the
    redesign, one of which recommended the two new rotation button names:
    「全部推奨通りでいいです。」
  - Reply 2, to the question 「回転ボタンと合わせて、この内容でよいか確認をお願いします。」
    asked with the receipt fields recorded below for the three records:
    「変えていいです。」
  - Reply 3, to the question 「受領の文面に、私の書き間違いが 1 か所あります。確認をお願いします。」
    asked with one corrected sentence of the `PD-OI-061` receipt: 「承認します。」
    The sentence Reply 2 was asked with said that **Merge** kept its name.
    Before this revision that button read
    `Merge into existing "<name>" and adopt its color`, so the `Choice` field
    of the receipt below carries the corrected sentence in its place. No other
    field of any receipt differs from the text Reply 2 confirmed.
  - `PD-OI-085`: the rotation buttons **Apply on Generate** and **Apply and
    regenerate** are named **Apply when I Generate** and **Apply now**. The
    record title and the concern key keep their wording.
  - `PD-OI-089` and `OIC-029`: the Color/Stroke Change Scope option
    `Apply to all "<legend name>" (N)` is named
    `Legend item "<legend name>" (N)`. The default-color line of `OIC-029` is
    unchanged.
  - `PD-OI-061` and `OIC-028`: the Legend Name Conflict choices **Suffix** and
    **Merge** are named `Keep both: name it "<name 2>"` and
    `Merge into "<name>"`. `PD-OI-088`, whose acceptance contract is
    `OIC-028`, is unchanged; the row changes only those two names.
  - The confirmed fields of each receipt are its choice, rationale,
    preservation, retirement scope, and accepted residual risk; the three
    receipts share the owner and the decision date of that confirmation. The
    concern key is the record's and the scenario revision is its next one.
    The confirmed fields for `PD-OI-089` and `PD-OI-061` name no choice code;
    their `Choice` field carries the confirmed statement of the label change,
    and no code is added. Wrapped lines of the confirmed text are joined; no
    other character differs.
  - Each earlier receipt keeps its text and its full scope. Where it or an
    earlier revision entry names an old label, it means the control that now
    carries the new label. The new receipts retire the old visible labels
    only.
  - This change is a static Product Contract co-change: the runtime, tests, and
    documentation that show the new labels are in the same pull request, and
    the Review is `REQUIRED`.
- Revision 41 changes: `PD-OI-067` gains scenario revision `2`
  (`B / PER-RECORD-INFERRED-DEFINITION, one-record File shows record names on the card`)
  from the Product Decision Owner `satoshikawato`'s replies of `2026-10-11`
  and `2026-10-12` on the names area of a Linear File card of the Web app.
  Each quotation below is the part of a reply that carries the answer; what
  it leaves out is stated.
  - Reply 1 (`2026-10-11`), to the review questions on the Before/After
    images of the redesigned Input Genomes sidebar, one of which recommended
    that a File with one record hide the File-level names block and show the
    record's name fields on its card. The reply first says that the Owner had
    looked at the images, then: 「推奨通りでいいです。」
  - Reply 2 (`2026-10-11`), to the request to confirm the draft receipt for a
    File with one record. The message before that request had named the
    draft file that holds the nine fields recorded below and had summarized
    what changes, what stays and the accepted risk; it did not show the nine
    fields themselves. The reply answered two requests, each with
    「承認します。」; the other one confirmed a correction of the `PD-OI-061`
    receipt (Revision 40).
  - Reply 3 (`2026-10-12`), to three questions asked with Before/After images
    of a prototype of this outcome. The reply also covers other matters; its
    line for these questions begins with the address of the page that asked
    them and continues:
    「A か B か:B.空のカードで先に入力した名前: 消さない; BGC の Gallery Session の保存し直し:保存し直す」
    - Its first answer selects option B of the question 「A か B か」. The
      question gave the rule of option B as:
      「カードに出す名前のブロックは常に 1 つ。File の名前があるあいだは dev と同じカード（File のブロックがカード、レコードの欄は Record options の中）で、File の名前が空ならレコードの欄をカードに出す。」
      In it `dev` is the Web app before this change. The option's text went
      on with one sentence that compared card heights; that sentence is a
      statement by the question's author, not part of the rule, and is not
      reproduced.
    - Its second answer keeps a name that was typed on a File card before
      its file is chosen. It is implemented in the same pull request and adds
      no term to this record.
    - Its third answer is about a Gallery example Session and adds no term to
      this record.
  - The scenario-2 receipt holds the nine fields of the draft that Reply 2
    approved, with no character changed. The first selection of Reply 3 is recorded in
    `PD-OI-067` below its receipts, in the words of its question. It is not a
    receipt field, and no rationale, preservation, retirement scope, or
    residual risk is added for it.
  - The scenario-1 receipt keeps its text and its SHA-256.
  - This change is a static Product Contract co-change: the runtime, tests, and
    documentation that implement it are in the same pull request, and the
    Review is `REQUIRED`.
- Revision 42 change: `PD-OI-039` is replaced for scenario revision `3`
  (`A / DOCKED-REVIEW-AT-EVERY-WIDTH`) from the Product Decision Owner
  `satoshikawato`'s replies of `2026-10-12` on the Review alignment panel of
  the Web app. Each quotation below is the part of a reply that carries the
  answer; what it leaves out is stated.
  - Reply 1 (`2026-10-12`), to four questions asked with same-state
    Before/After images of a mock of the redesigned panel: strips in the
    panel that show each record after Apply, rings and arrows on the figure
    in place of the candidate numbers, the panel docked at the right of the
    figure, and the records that need a choice first. The reply, complete,
    with the questions it quotes (its line breaks are shown as ` / `):
    「>どこまで変えるか。 / 全部。 / >右に固定する間、遺伝子のポップアップを隠してよいか。 / 隠していいです。 / >英語の文言。 / 任せます。 / >時期。 / 0.14.0に入れる。」
    Its first answer selects the whole proposal, which docks the panel. Its
    second answer lets the feature popup be hidden while the panel is docked.
    Its third and fourth answers add no term to this record.
  - Reply 2 (`2026-10-12`), to the request to approve the receipt recorded
    below. The request showed the nine fields in full, twice: a first
    wording, and a second wording that was named as the current one and
    differed in two passages about the feature popup (in Must preserve, that
    a click on the figure neither opens, switches nor closes the popup during
    the review and that the search **Open** switches the popup's target and
    shows it when the review ends; in May retire, the matching split). The
    reply: 「文面を承認します」
    The shown text had a placeholder for the decision date; the date recorded
    is the day of this reply.
  - The receipt below holds the nine fields of the second wording, with no
    character changed apart from that date.
  - The runtime that realizes this outcome merged before this amendment, as
    pull request `#1020` (merge commit `d8af5a7f`), which the Owner merged on
    `2026-10-12` before Reply 2. This amendment contains no runtime.
  - Scope: scenario revision `2` (`A / EXCLUSIVE_DIRECTIONS_WITHOUT_MATCH`,
    receipt SHA-256
    `06a9d2fe9b1d1406f6f8e04c23a9ca031133b9fae24b0683403ec2c6cae55270`) is
    retained in Git history, and what it retired stays retired. `PD-OI-035`
    keeps its text; only its reference to the governing `PD-OI-039` scenario
    revision now reads `3`. The `OIC-026` row is reworded to the receipt's
    dock, Editor, popup and view terms. `PD-OI-031`, `PD-OI-037`,
    `PD-OI-038`, `PD-OI-054` and all four accepted issue `#598` decisions
    are unchanged.
- Revision 43 changes: `PD-OI-092` is added
  (`A / PINNED-ICON-SECTION-INDEX`, concern `web.settings-pane.section-index`).
  It comes from the Product Decision Owner `satoshikawato`'s replies of
  `2026-10-12` on reaching a settings section in the left settings pane of
  the Web app. Each quotation below is the part of a reply that carries the
  answer.
  - Reply 1, to a measurement shown with images of the pane: in Linear with
    five records the Labels section was 4,009 px and the Comparison section
    2,070 px down the pane.
    「これも重要だね。もしかすると「今」もこの部分は改善が必要かも。」
  - Reply 2, to images of the pane with one pinned row above it: a section
    index whose icons jump to a section and open it, with an entry of its own
    for the Linear Comparison section. 「この形で進めていいです。」
    The same images showed a button that switches the pane to a second
    display; that part is not in this record.
  - Reply 3, to the question whether the index shows icons only or icons with
    names.
    「目次はアイコンでいいです。マウスオーバーすると説明(あるいはタイトル）が浮かび上がるようにしてほしい。」
  - Reply 4, to the request to approve the receipt text as written or to name
    what to change, asked with Before/After images of the built index.
    「そのまま承認します。」
    The page that asked showed the Rationale and the Accepted residual risk as
    recorded below, the Must preserve text as a list whose `〔…〕` headings were
    set in bold without the brackets, the Choice as its identifier with a
    description, May retire as 「なし（none）」, and the Owner. It gave the
    decision date as the day of the confirmation and did not show the concern
    key or the scenario revision. The receipt is recorded byte for byte.
  - The distances in the Rationale (3,909 px and 1,972 px) are those of the
    Web app on the day of the confirmation; a change to the Linear File cards
    had shortened the pane since Reply 1.
  - No other record is amended or retired. The record names what it leaves
    unchanged in its Must preserve.
  - This change is a static Product Contract co-change: the runtime, tests, and
    documentation that implement it are in the same pull request, and the
    Review is `REQUIRED`.
- Revision 44 changes: `PD-OI-093` is added
  (`A / STACKED-DEFAULT-WITH-ICON-RAIL-SWITCH`, concern
  `web.settings-pane.display`). It comes from the Product Decision Owner
  `satoshikawato`'s replies of `2026-10-12` on a second display of the left
  settings pane of the Web app. Each quotation below is the part of a reply
  that carries the answer.
  - Reply 1, to images of today's pane beside three other places for the
    settings: `A` a column of icons with one section beside it, `B` a grid of
    options with a floating popup, `C` a toolbar above the Result.
    「今とAを行き来できるようにするのがいい気がする。デフォは今のままで、ボタンひとつで切り替えられる形。」
  - Reply 2, to images of the pane with a button at the right end of the
    section index row that switches it to the column of icons, a button on
    top of the column that switches back, and the choice remembered by the
    browser and not stored in the Session. 「この形で進めていいです。」
    This is Reply 2 of the Revision 43 entry; `PD-OI-092` recorded its section
    index part.
  - Reply 3, given on a page of Before/After images of the built display,
    after the Owner asked where the button that folds the pane is and was
    told that it is at the right end of the section index row.
    「OKOK,これでいきましょう。全部承認します。」
  - Reply 4, to the three requests of the same page: approve the receipt text
    as written or name what to change, choose one of three outcomes for the
    figure while Review alignment is open, and confirm the four English
    button names. 「推奨通りでいきます。Product Decisionも承認します」
    The page showed the nine fields of the receipt as recorded below, the
    Must preserve text as a list whose `〔…〕` headings were set in bold
    without the brackets, and the decision date as the day of the
    confirmation. The receipt is recorded byte for byte.
  - The recommended outcome that Reply 4 selects for the figure while Review
    alignment is open (fit the figure again whenever the room of the Result
    changes) is not part of this record. It is a later change; until it
    lands, the last sentences of the Accepted residual risk below describe
    the Web app.
  - `PD-OI-092` keeps its text. The receipt below places the section
    index row of `PD-OI-092` in the stacked display (「今の表示」),
    which stays the default, and keeps it unchanged there (`〔変えないもの〕`).
    The Accepted residual risk of `PD-OI-092` gives the index buttons
    as 22 × 26 px; this record's Accepted residual risk states that the
    Linear ones become 21.2 px wide (20.5 px in a browser that always shows
    its scrollbar) to make room for the switch button.
  - No other record is amended or retired.
  - This change is a static Product Contract co-change: the runtime, tests, and
    documentation that implement it are in the same pull request, and the
    Review is `REQUIRED`.
- Revision 45 change: `PD-OI-061` is replaced for scenario revision `4`
  (`B / MERGE-SAME-FEATURE-TYPE-ONLY`, extended to the names of color rules
  and to the list under the feature popup's **Legend name**). It comes from
  the Product Decision Owner `satoshikawato`'s replies of `2026-10-13` on a
  user report (OV-580, OV-581, OV-582): in the feature popup of the Web app,
  choosing another color rule's name from the list under **Legend name** kept
  the feature's own color, so one Legend name had two colors and both rows
  were drawn as `<name> [#<hex>]`; renaming a color rule's row onto another
  color rule's row ended the same way without asking. Each quotation below is
  the part of a reply that carries the answer.
  - Reply 1, on the reported Legend:
    「凡例に色のIDがつくのは許されないので直して。Ownership couplingとか境界型付けの問題でもあると思う。」
  - Reply 2, on the list under **Legend name**:
    「Select legend listがまずバグってる気がする」
  - Reply 3, on merging a whole Legend item into another one:
    「前は普通にマージできたんだよ全部Transcription regulationに。」
    Release `0.13.0` asked the Legend Name Conflict dialog for every rename
    onto an existing name of another color. The scenario-1 receipt
    (`2026-09-30`) preserved the `PD-OI-042` naming for a target that is a
    color rule's caption, the scenario-2 receipt kept that term, and the
    runtime has carried it since commit `4b7b8825` (`2026-10-02`).
  - Reply 4, to the request to confirm the ready-to-copy `PRODUCT_DECISION`
    response for this concern and scenario revision, prepared with the
    decision pack for the report:
    「Web の操作(PD-OI-061 改訂 4): それでいいです」
  - The receipt below is that response, byte for byte.
  - The decision pack described the same outcome in English, with one phrase
    the receipt does not repeat: a typed rename onto a color rule's name asks
    the dialog "also when the target is a color rule's name, listed in the
    legend or not". That phrase is the pack's description, not a receipt
    field. The receipt's 「衝突先が色ルールの名前でも」 names a color rule's name
    without a condition on a drawn row, and its Accepted residual risk names
    a color rule that the figure does not draw.
  - Reply 5, after the review of the implementation, to the question
    「全 feature type 向けのルール(`*`)の扱い」 asked with the recommended fix
    「一覧に出る項目は必ず入れる。`*` のルールはどの type の一覧にも出す」:
    「全 feature type 向けのルール(*)の扱い: 推奨通りなおしましょう」
    No field of the receipt differs from the text Reply 4 confirmed.
  - Scope: scenario revision `4` keeps the outcome
    `B / MERGE-SAME-FEATURE-TYPE-ONLY` of scenario revision `2` and the
    choice names of scenario revision `3`. It retires one preservation term
    of the scenario-2 receipt, 「衝突先が色ルールの caption のときの PD-OI-042 の区別」,
    and adds the outcome of a pick from the popup's list. The scenario-2 and
    scenario-3 receipts keep their text and their SHA-256.
  - `PD-OI-042` is unchanged. One Legend name written with several colors in
    a color table, by **Add rule**, or on the command line is still drawn as
    `<name> [#<hex>]` rows; the Owner's decision on that concern is not part
    of this revision. `PD-OI-088` and `OIC-028` are unchanged.
  - This change is a static Product Contract co-change: the runtime and tests
    that implement it are in the same pull request, and the Review is
    `REQUIRED`.
- Revision 46 changes: `PD-OI-094` is added
  (`A / FOLDABLE-FILE-CARD, FOLDED-AFTER-LOAD`, concern
  `web.linear.file-card-fold`). It comes from the Product Decision Owner
  `satoshikawato`'s replies of `2026-10-12` and `2026-10-13` on folding the
  File cards of Linear mode in the Web app to one line each. Each quotation
  below is the part of a reply that carries the answer.
  - Reply 1 (`2026-10-12`), to mock images of the settings pane in which the
    five File cards of a Linear Session are one line each.
    「ファイルのカードを畳む案を入れる。」
  - Reply 2 (`2026-10-13`), to a page of Before/After images of the built
    cards that asked for one choice: `A`, a card that holds its files is
    folded after **Load Session** or **Load example** (recommended, built), or
    `B`, no card folds by itself. 「PR-B: 推奨通り。」
  - Reply 3 (`2026-10-13`), to the request to approve the receipt text of
    that page as written or to name what to change.
    「PR-B の文面の承認: 承認します。」
    The page showed the nine fields of the receipt as recorded below, the
    Must preserve text as a list whose `〔…〕` headings were set in bold
    without the brackets, and the decision date as the day of the approval.
    The receipt is recorded byte for byte.
  - `PD-OI-018`, `PD-OI-025`, and `PD-OI-067` keep their text. They describe
    the contents of a File card, which the receipt below leaves as they are
    in an open card (`〔畳んでも変わらないもの〕`). A folded card shows its
    header line only; the Accepted residual risk below states what is then
    out of sight until the card is opened, the common Depth TSV assignment
    of `PD-OI-025` included.
  - No other record is amended or retired.
  - This change is a static Product Contract co-change: the runtime, tests, and
    documentation that implement it are in the same pull request, and the
    Review is `REQUIRED`.
- Revision 47 changes: `PD-OI-042` is replaced for scenario revision `2`
  (`B / NUMBERED-NAMES`, concern
  `styles.specific-color-caption-multiplicity`), and `PD-OI-095` is added
  (concern `legend.row-name-multiplicity`; its receipt names no choice code).
  They come from the Product Decision Owner `satoshikawato`'s replies of
  `2026-10-13` on the user report of the Revision 45 entry and on three
  findings made with it (OV-583, OV-584, OV-585): one Legend name written
  with several colors was drawn as `<name> [#<hex>]` rows by the CLI, the
  Python API, and the Web app; the rule list, the popup, and the downloaded
  color table held the rewritten captions; and the second track row of one
  label was named `<label> (2)` by another mechanism. Each quotation below is
  the part of a reply that carries the answer.
  - Reply 1, on the reported Legend. This is Reply 1 of the Revision 45
    entry.
    「凡例に色のIDがつくのは許されないので直して。Ownership couplingとか境界型付けの問題でもあると思う。」
  - Reply 2, to the question on one Legend name written with several colors,
    asked with three outcomes: `A`, one row per color, each drawn with the
    name as written; `B`, numbered names (`name (1)`, `name (2)`), which the
    question described as still rewriting the names and the rule captions;
    `C`, one row with a swatch per color.
    「1 つの名前に複数の色があるとき(PD-OI-042 改訂 2):色が違う場合はB.色が一緒の場合はマージ。てかそうだったでしょ今まで？前はPharokkaでもできたんだから。」
  - Reply 3, on the report that two mechanisms make a Legend name unique
    (numbers and color IDs): 「これはだめだろ。」
  - Reply 4, to the question whether the first color keeps the name without
    a number while later colors get `(2)`, `(3)`:
    「2色目が出るなら全部番号アリがいいのでは？今までそうだっただろ？」
  - Reply 5, to the question whether the Legend Name Conflict choice
    **Keep both** starts at `(2)` instead of `(1)`:
    「何言ってるのかわからないけど、(1), (2), ってなるべきでは？」
    That choice belongs to `PD-OI-061`, which this revision does not change.
  - Reply 6, to the question whether track rows with one label, drawn
    `GC content` and `GC content (2)` since release `0.13.0`, follow the same
    rule:
    「色が同じなら凡例は1個のままでいいよね。色が違うならGC content (1) / GC content (2) に揃える」
  - Reply 7, to the question on the feature popup's **Legend name** field,
    whose suggestions and the list under it offered the same names with two
    results (OV-587), asked with two options: remove the field's suggestions
    (recommended), or make a chosen suggestion join the group.
    「入力欄の候補と一覧の食い違い(OV-587)：候補を消す(推奨)でいきましょう」
    The removal is implemented in the same pull request and adds no term to
    a record.
  - Reply 8, to the request to confirm three receipt texts (`PD-OI-042`
    scenario revision `2`; the new concern `legend.row-name-multiplicity`,
    scenario revision `1`; `PD-OI-061` scenario revision `5`) and two
    choices: 「選択 1: 保存済みの色 ID の扱い」 with
    「A（推奨）読み込むときに外す / B そのまま読む」, and 「選択 2: OV-587」 with
    「A（推奨）入力候補を外す / B 入力候補を選んでも一覧と同じ動作にする」. The request offered
    the reply form 「文面 1〜3 は承認。選択 1 は A。選択 2 は A。」
    「すべて推奨通り。」
  - The two receipts recorded below are the first and the second of those
    three texts, byte for byte. Option `A` of 選択 1 is in the Accepted
    residual risk of the `PD-OI-042` receipt (「読み込むときにその ID を外す」).
    Unlike the description of `B` in the question of Reply 2, the
    `PD-OI-042` receipt keeps a rule's caption as written
    (「rule の caption は書かれたとおりに保つ」). The third text is not recorded by this revision: `PD-OI-061` keeps its
    text and its scenario revision `4`.
  - Scope: scenario revision `2` of `PD-OI-042` replaces scenario revision
    `1` (`A / CAPTION_AUTO_DISAMBIGUATE_SOLID_ROWS`) as a whole. The
    scenario-1 outcome fields, its receipt, and its Acceptance IDs (`CLR-01`
    to `CLR-05`) are not carried into the record; Git history retains them.
  - No other record is amended or retired.
  - This change is a static Product Contract co-change: the runtime, tests, and
    documentation that implement it are in the same pull request, and the
    Review is `REQUIRED`.
- Revision 48 changes: `PD-OI-096` is added
  (`B / MODE-DEFAULT-LISTED-AT-LOAD`, concern `web.legend.position-of-mode`). It comes from the
  Product Decision Owner `satoshikawato`'s replies of `2026-10-13` on a legend
  position that is not one of its diagram mode's positions in the Web app:
  the Legend position setting, and the side a loaded Result's Legend was laid
  out with. Each quotation below is the part of a reply that carries the
  answer.
  - Reply 1, to a report line that named three outcomes for such a setting:
    「選択肢は、モードの集合に正規化、Load で拒否、現状維持の 3 つです。」
    「推奨どおりでよろしくです」
  - Reply 2, to the question 「OV-652(BA-D07)の確認: そのモードの集合に無い凡例位置(例: Linear で "upper_left" や "banana")は、どう扱いますか?」
    with the options 「モードの既定値にする(推奨)」 (「プロットタイトル位置と同じ扱いです。集合に無い値はそのモードの既定の凡例位置になります。」)
    and 「近い位置へ対応付ける」. The Owner selected 「モードの既定値にする(推奨)」.
  - Reply 3, a message on loading older files:
    「てか古いファイルを拒否するよりもさ、「これこれこういう不具合というか、バグがあるから一部完全には復元できないけどいい？次回保存するときは細心のフォーマットになるから」っていうひとことnoticeを出して読み込めるようにするのがいいんじゃないかな。」
  - Reply 4, to four questions that confirmed Reply 3. The three answers this
    record uses:
    - 「一言 notice の出し方はどちらにしますか?」: the Owner selected
      「読み込む前に確認(推奨)」 (a dialog before the load, with Load and
      Cancel; Cancel leaves the current screen unchanged).
    - 「拒否をやめて「読めるところだけ読んで notice」にする範囲はどこまでですか?」:
      the Owner selected 「gbdraw のファイルなら全部(推奨)」.
    - 「意味が変わらない自動変換(古いスキーマの昇格、"RBH" を "rbh" にするような綴りの正規化)でも notice を出しますか?」:
      the Owner selected 「出さない(推奨)」.
  - Reply 5, to three outcomes for the side a loaded Result's Legend was laid
    out with when that side is not one of the Result's mode's positions (`A`,
    `none`; `B`, the mode's default; `C`, keep the stored corner). The Owner
    quoted the recommended line: "B (recommended): use the mode's default, so
    the legend moves below the figure, and list it in the dialog. This is the
    same handling as OV-652, and Sessions 31 and later already read it this
    way."
  - Reply 6, to the request to approve the receipt text below as written or to
    name what to change. The request showed the nine fields of the receipt as
    recorded below and said that one meaning differs from an earlier draft: a
    Legend position setting outside its mode reads as an unset value does
    (usually the mode's default).
    「この文面のまま承認します」
    The receipt is recorded byte for byte.
  - `OIPC-C01` keeps its text. The values this record decides are stored in a
    file that is being loaded, not supplied in a control or a request; the
    receipt lists them before the load, and its Accepted residual risk states
    the interim in which they are replaced without a notice.
  - No other record is amended or retired.
  - This change is a static Product Contract co-change: the runtime, tests, and
    documentation that implement it are in the same pull request, and the
    Review is `REQUIRED`.
- Revision 49 changes: `PD-OI-061` is replaced for scenario revision `5`
  (concern `web.legend.rename-collision`; its receipt names no choice code).
  It comes from the Product Decision Owner `satoshikawato`'s replies of
  `2026-10-13` on the Legend Name Conflict choice **Keep both** after
  Revision 47: `PD-OI-042` scenario revision `2` numbers every drawn row of a
  name shared by rows of different swatches, while **Keep both** gave the
  renamed row alone a numbered name and left the other row's name without a
  number. Each quotation below is the part of a reply that carries the
  answer.
  - Reply 1, to the request to confirm three receipt texts (`PD-OI-042`
    scenario revision `2`; the new concern `legend.row-name-multiplicity`,
    scenario revision `1`; `PD-OI-061` scenario revision `5`) and two
    choices. This is Reply 8 of the Revision 47 entry, which records the
    request in full.
    「すべて推奨通り。」
    The receipt recorded below is the third of those texts, byte for byte.
  - Reply 2, to three options for this change:
    「A（推奨）PR 3 を 0.14.0 から外して止める / B 小さい PR 3（色ルールと feature の改名だけ rev 5）/ C Session の形式変更まで含めて全部やる（一番大きい）」
    「C: Session の形式変更まで含めて全部やります（一番大きい）。もう始めたんでしょ？あとどれくらいかかるの？」
  - Reply 3, to the question
    「P-8（Session の版）: 推奨は、この PR で version 47 にすることです。ただし、複数 drawing の決定（PD-OI-090）の文面が「47」を予定しているので、そちらを 48 に読み替える必要があります。」
    「P-8（Session の版）: 推奨通り」
    Session version `47` is written by this change.
  - Reply 4, to the question
    「P-9（同じ色の行で Keep both を選んだとき）: 推奨は、「同名同色は 1 行」の規則どおり 1 行にまとめることです。保留中の P-5 も同じ扱いになります。」
    「P-9：推奨通り」
    P-5 is the item of that name in the change's decision record
    (`pharokka-legend-merge-20261013/DECISIONS.md`); this entry does not
    restate it.
  - Reply 5, to the question
    「P-10（feature を 1 つだけ切り出して Keep both にしたときの新しい行の位置）: 推奨は、同じ名前の行のすぐ後ろに置くことです。こうすると、ダイアログに出る (2) が凡例の表示と一致します。」
    「P-10:　推奨通り」
  - Reply 6 (`2026-10-14`), to the question
    「残る制限（D-33）: 同じ名前・同じ色の 2 グループは 1 行にまとめて描きます。その行の色を変えると、変わるのは見えている方のグループだけです。もう一方は元の色のまま別の行として戻り、凡例は「名前 (1)」「名前 (2)」の 2 行になります。また、Generate の後にこれが起きると、戻った行が凡例タブの一覧に出るのは次の描画のときです。どうしますか？」
    the selected option and its description:
    「PR 3 の中で直す (Recommended)」
    「まとめた行の色を変えると、その行にまとめたすべてのグループの色が変わり、1 行のまま残る。凡例タブもすぐ更新される。PR 3 のマージは半日ほど遅れますが、0.14.0 は boundary-audit を待つので、リリース日は変わりません。」
  - Scope: scenario revision `5` keeps the outcome
    `B / MERGE-SAME-FEATURE-TYPE-ONLY`, the choices of the dialog, and the
    conditions of **Merge**. The scenario-2, scenario-3, and scenario-4
    receipts keep their text and their SHA-256.
  - `PD-OI-042`, `PD-OI-095`, and `PD-OI-090` keep their text. No other
    record is amended or retired.
  - This change is a static Product Contract co-change: the runtime, tests, and
    documentation that implement it are in the same pull request, and the
    Review is `REQUIRED`.
- Revision 50 changes: `PD-OI-097` is added
  (concern `product.circular-outer-label-fit`, scenario revision `7`). It comes
  from the Product Decision Owner `satoshikawato`'s replies of `2026-10-13` and
  `2026-10-14` on Circular and Linear labels that do not fit: labels that
  overlap, leader lines that cross, and labels over the definition text in the
  middle of a Circular map. Each quotation below is the part of a reply that
  carries the answer.
  - Reply 1 (scenario revision 1): 「決めていただくこと: DEF-CROSS、DO-ASK、TIGHT-FIX-FIRST。」
  - Reply 2 (scenario revision 2): 「SHOW-ALL,、REL-0140だと思うけど、 >0.13 では既定設定で描けていた WSSV、153 kb
    プラスミドの 8 pt、Epifagus が、0.14.0 では止まります なんで？WSSVはかけてたじゃん。」
    「今まで描けてたものは今まで通り描けるべきだろ。」
  - Reply 3 (scenario revision 3): 「CROSS-WARN、っていうか引き出し線が交差しているってどういうこと？」
    「owner-seen-fit-def_crossing_n4-20261013.png これはダメだね。WSSVは問題ないけどこれはダメ。これは警告入りしていい」
  - Reply 4 (scenario revision 4): the selected options 「開いてよい (推奨)」, 「通知の帯 (推奨)」,
    「時間を区切る (推奨)」 and 「全部この PR に」, and the reply 「まあええよ」.
  - Reply 5 (scenario revision 5): 「全部推奨」
  - Reply 6 (scenario revision 6, `2026-10-13`): 「E-ROW（推奨）」
  - Reply 7 (`2026-10-14`), to four questions on the receipt fields, each
    answered by option: Rationale 「案のままでよい (推奨)」, Must preserve
    「案のままでよい (推奨)」, May retire 「案のままでよい (推奨)」, and Accepted
    residual risk 「交差したりラベルが定義文にかかるのも、やっぱりラベルの配置を微調整することで直すことを提案してほしい」.
  - Reply 8 (scenario revision 7, `2026-10-14`), to the automatic repair of
    slight defects and to a fix row that leaves crossings: 「入れる、上限 100 (推奨)」
    (`REPAIR-100`) and 「交差が残る行も出す (推奨)」 (`CROSS-LEFT-ROW`).
  - Reply 9 (`2026-10-14`), to the rewritten Accepted residual risk
    「微調整で全部直る重なり・交差・定義文へのかぶり（欠陥 100 組まで）は、自動で直して描く。直しきれないもの（極の内側ラベルの積み重なり、長い線のもつれ、MjeNMV など）は、重なりなら描かずに直し方を選ばせ、交差やかぶりだけなら描いて帯で知らせる。直し方の候補は、帯のボタンを押したときに計算する。この結果は Session に保存しない。CLI と Python は、描けた図なら 5 秒、描けない図なら 30 秒まで探し、間に合わない候補は『未計算』と示す。」:
    「これでよい (推奨)」.
  - Reply 10 (`2026-10-14`), to two questions after the review of this
    change. To the question whether to record Reply 9's sentence on search
    times as written or to correct it to the option text of Reply 4's
    「時間を区切る (推奨)」, which says 「重なりで止まる場合は値を必ず出す」:
    「Q9 に合わせて直す (推奨)」, whose option text replaces that sentence with
    「CLI と Python は、描けた図なら 5 秒まで探す。描けない図は、外側ラベルだけの重なりなら値が出るまで探し、内側・埋め込み・放射ラベルの重なりなら 30 秒まで探す。間に合わない候補は『未計算』と示す。」.
    To the question what Apply and Cancel do while a fix row is being
    computed in the Web: 「1 候補 30 秒まで (推奨)」, whose option text adds
    「Web の候補は 1 つ 30 秒まで探し、見つからなければそう示す。」.
  - The receipt below combines these replies: Choice lists the options of
    scenario revisions 1 to 7, Rationale, Must preserve, and May retire are the
    texts of Reply 7, and Accepted residual risk is the text of Reply 9 as
    corrected by Reply 10.
  - `PD-OI-088` (a dialog choice is one History step, Cancel records none, the
    dialog is busy while applying) applies to the new dialog unchanged. No
    other record is amended or retired.
  - This change is a static Product Contract co-change: the runtime, tests, and
    documentation that implement it are in the same pull request, and the
    Review is `REQUIRED`.
- Revision 51 changes: `PD-OI-098` is added
  (`A+B / CHOSEN-FONT-MEASURED-AND-EMBEDDED`, concern
  `diagram-generation.diagram-font`), and `PD-OI-092` (concern
  `web.settings-pane.section-index`) and `PD-OI-093` (concern
  `web.settings-pane.display`) are replaced for scenario revision `2`. They
  come from the Product Decision Owner `satoshikawato`'s reply of
  `2026-10-13` on choosing the font of the diagram text in the Web app, the
  CLI, and the Python API, asked with the three receipt texts below: the
  new concern, and the two settings pane records that gain a Font section and
  a Font item.
  - Reply, to the request to approve the three receipt texts as written or to
    name what to change:
    「お願い（1 つ）: 三つとも承認します」
  - The receipts recorded below are those three texts, byte for byte. The
    SHA-256 of each is over its UTF-8 bytes without the final newline:
    `PD-OI-098` `14dc5e981573dce50e9815f4457d34de6786c819ae31a5b4672e178359af3675`,
    `PD-OI-092` `0fba5433eeb10e8c606b967f142c4a33b938d23db7fc3d6dc7fde5396e636c92`,
    `PD-OI-093` `81bf15fcf9788503b438545f62256e6eb72f8e21f9e0f1213245ea3f858316c0`.
  - Scope: scenario revision `2` of `PD-OI-092` keeps the outcome
    `A / PINNED-ICON-SECTION-INDEX` and adds the **Font** section (icon
    `Aa`, between Layout and Titles & Record Labels) to the index, so the
    Circular index has 12 buttons and the Linear index 13. Scenario revision
    `2` of `PD-OI-093` keeps the outcome
    `A / STACKED-DEFAULT-WITH-ICON-RAIL-SWITCH` and adds the Font item to the
    icon rail and the Linear index button width of its Accepted residual
    risk. Each replaces scenario revision `1` as a whole; Git history retains
    the former receipt.
  - `PD-OI-098` does not supersede another decision. Sessions already kept
    the font name as the preserved override `objects.text.font_family`; the
    Web setting is saved as `modes.<mode>.config.adv.font_family` in Session
    47 (see `docs/internal/SESSION_FORMAT_HISTORY.md`), and no Session stores
    a font file.
  - No other record is amended or retired.
  - This change is a static Product Contract co-change: the runtime, tests, and
    documentation that implement it are in the same pull request, and the
    Review is `REQUIRED`.
- Revision 52 changes: `PD-OI-099` is added
  (`A / CONFIRM-REPAIRS-BEFORE-LOAD`, concern
  `web.session.load-repair-confirmation`), with acceptance contract `OIC-031`.
  It comes from the Product Decision Owner `satoshikawato`'s replies of
  `2026-10-13` and `2026-10-14` on how the Web app loads a gbdraw file that is
  old or holds values it cannot use. Each quotation below is the part of a
  reply that carries the answer.
  - Reply 1 (`2026-10-13`), a message on loading older files. It is Reply 3 of
    the Revision 48 entry.
    「てか古いファイルを拒否するよりもさ、「これこれこういう不具合というか、バグがあるから一部完全には復元できないけどいい？次回保存するときは細心のフォーマットになるから」っていうひとことnoticeを出して読み込めるようにするのがいいんじゃないかな。」
  - Reply 2 (`2026-10-13`), to four questions that confirmed Reply 1, each with
    the recommended option first. Reply 4 of the Revision 48 entry quotes the
    questions of (a), (b), and (d).
    - (a) how the notice is shown: the Owner selected
      「読み込む前に確認(推奨)」.
    - (b) the scope of loading as far as possible instead of refusing: the
      Owner selected 「gbdraw のファイルなら全部(推奨)」.
    - (c) files the Web refused and Python reads: the Owner selected
      「Web も notice 付きで読む(推奨)」.
    - (d) a notice for conversions that do not change a meaning: the Owner
      selected 「出さない(推奨)」.
  - Reply 3 (`2026-10-14`), to the status line
    「PR-N（Load の修復ダイアログ）の receipt は Owner の承認待ちです」 of another
    Claude Code session, which relayed the reply; the reply is in that
    session's transcript.
    「承認します」
    The status line names this pull request by its working name. The receipt
    approved is draft v3.6 of the receipt text, recorded below byte for byte.
    That session did not show the receipt text. Earlier, the session that
    drafted the receipt had shown the full text of draft v3.1 on `2026-10-13`.
    For drafts v3.2 to v3.6 it gave a summary of each draft's changes, not the
    full text, and with drafts v3.2, v3.3, and v3.5 the path of the draft file.
    It had not yet asked for approval of v3.6. After the reply, the relaying
    session told the Owner that a later change to the text needs a new
    approval.
  - `OIPC-C01` and `OIPC-C03` keep their text. The values this record decides
    are stored in a file that is being loaded, not supplied in a control or a
    request; the receipt lists them before the load, except the values its
    〔出さないもの〕 names, and its Accepted residual risk states the interim.
  - `PD-OI-096` keeps its text; the confirmation dialog its Accepted residual
    risk names is the dialog of this record. No other record is amended or
    retired.
  - This change is a static Product Contract co-change: the runtime, tests, and
    documentation that implement it are in the same pull request, and the
    Review is `REQUIRED`.
- Revision 53 changes: `PD-OI-097` is replaced for scenario revision `8`
  (concern `product.circular-outer-label-fit`). It comes from the Product
  Decision Owner `satoshikawato`'s reply of `2026-10-14` on where the labels
  around a repaired defect yield room after the automatic repair of scenario
  revision `7`. Each quotation below is the part of a reply that carries the
  answer; Replies 1 to 10 are those of the Revision 50 entry.
  - Reply 11 (scenario revision `8`, `2026-10-14`), to the question
    「【譲り合いを効かせる範囲】今収まっている Epifagus 14 pt に、譲り合いを効かせた場合を 4 つ並べました（logs/owner-seen-P1_epifagus14_scope-20261014.png。tobacco は P2）。1 は今の #1053、2 は接点の規則だけ（ラベルは動かず、見た目はほぼ同じ。追加 0.01 s）、3 は譲り合い（50 個が動き、追加 2.0 s）、4 は譲り合い＋上へ譲る（51 個が動き、追加 2.8 s）です。小さな図でも、譲り合いを効かせると全ラベルが動き、約 1 s 増えます。4 では maturase の線が tRNA-His の文字を貫いています（拡大：owner-seen-P1_panel4_maturase_through_tRNA-His-20261014.png）。この種類の欠陥は、どの選択肢でも起きないようにします。譲り合いをどこに効かせますか？」
    the selected option and its description:
    「微調整した図だけ (推奨)」
    「修正が入った図（N4 など）だけ、修正の後で周りを整える。収まっている図は 0.13.0 と同じまま（②を変えない）。速さも変わらない」
    (`YIELD-REPAIRED-ONLY`)
  - Reply 12 (`2026-10-14`), approving the receipt text below:
    「Contract 受領文の案をOKします」
  - The receipt of this change (SHA-256, UTF-8, excluding the final newline:
    `38b212f4fa546885f1edf3767718cb8b15811b6dc30b195b1f7f4312c57eba3b`):

    ```text
    PRODUCT_DECISION
    Concern: product.circular-outer-label-fit
    Scenario revision: 8
    Choice: YIELD-REPAIRED-ONLY
    Rationale: not stated
    Must preserve: not stated
    May retire: not stated
    Accepted residual risk: not stated
    Owner: Satoshi Kawato
    Decision date: 2026-10-14
    ```

  - Scope: the record's receipt for scenario revision `8` is recorded below;
    the scenario revision `7` receipt keeps its text and its SHA-256.
    `PD-OI-088` applies unchanged. No other record is amended or retired.
  - This change is a static Product Contract co-change: the runtime, tests, and
    documentation that implement it are in the same pull request, and the
    Review is `REQUIRED`.
- Revision 54 changes: `PD-OI-100` is added (`A / NO-UNDO-STEP`, concern
  `web.history.rich-feature-popup`). It comes from the Product Decision Owner
  `satoshikawato`'s replies of `2026-10-16` on whether the **Rich Feature
  Popup** checkbox records an Undo step (it did, one step per toggle, and Undo
  or Redo of other steps changed it back). The questions and options are
  recorded in campaign `band-hover-perf-20261015`, `DECISIONS.md` `D-03`
  (`Q-R1`; A: it records no Undo step, B: keep recording a step) and `D-07`.
  - Reply 1, selecting option A:
    「D-03: Undoにすべきではない。」
  - Reply (a), to the question asked again with examples after the Owner's
    「具体的にどういうシチュエーション？」:
    「D-03 補足(a) の具体例です。
    例1: Popup ON → 色を変更（Undo 1段）→ Popup を OFF（Undo にならない）→ Undo を押す。色は戻る。Popup は OFF のまま？それとも色を変える前の ON に戻る？
    例2: Popup OFF → Reset Settings（Popup を含む全設定が既定に戻り Popup ON）→ Undo を押す。ほかの設定は Reset 前に戻る。Popup は ON のまま？OFF に戻る？」
    the selected option and its description:
    「どちらも今の値のまま」
    「例1 は OFF のまま、例2 は ON のまま。Popup は Undo/Redo では一切変わらない」
  - Reply (b), to the question
    「D-03 補足(b): Popup だけを切り替えた後に Load Session / Load example を押したとき、「Replace the current work?」の確認を出しますか？（Popup の値は Session に保存されるので、出さないと黙って上書きされます。今の PR-C は出しません）」
    the selected option and its description:
    「確認を出す」
    「Undo の段は作らず「未保存の変更あり」だけ記録して、確認ダイアログを出す」
  - Reply (c), to the question
    「D-03 補足(c): Auto Reflow、PNG DPI、palette の Instant Preview、LOSAT の実行方法とスレッド数は、今も切り替えるたびに Undo の段になります。どうしますか？」
    the selected option and its description:
    「今回は対象外」
    「今回は Rich Feature Popup だけ。ほかは別途決める」
  - Reply (d), to the question
    「D-03 補足(d): Product Contract の記録に残す理由はありますか？（なければ「not stated」と記録します。理由を書く場合は「その他」に入力してください）」
    the selected option:
    「not stated で記録」
  - The receipt of `PD-OI-100` below takes effect with the Owner's approval
    on the pull request that adds this revision; that pull request does not
    merge without it.
  - Scope: only the **Rich Feature Popup** preference. Per Reply (c), the
    other preferences are outside this decision and this revision does not
    change them. No other record is amended or retired.
- Revision 55 changes: `PD-OI-096` is replaced for scenario revision `2`
  (concern `web.legend.position-of-mode`). It comes from the Product Decision
  Owner `satoshikawato`'s replies of `2026-10-15` and `2026-10-16` on a
  Session that release `0.13.0` saved after **Generate** and a switch to the
  other diagram mode (OV-677): the Web app opened its Result in the shown
  mode. Each quotation below is the part of a reply that carries the answer;
  Replies 1 to 6 are those of the Revision 48 entry.
  - Reply 7 (`2026-10-15`), to the question
    「OV-677: 0.13.0 で Circular を描いた後に Linear に切り替えて保存した Session は、0.14.0 で開くと、描いたモード(Circular)ではなく表示中のモード(Linear)の結果として読まれます。0.13.0 自身も同じ読み方をしていたので、退行ではありません。どう扱いますか?」
    with the options 「A. 0.14.0 後に回す (Recommended)」,
    「B. 0.14.0 前に直す」, and 「C. 0.14.0 前に古い結果を捨てる」:
    「B. 0.14.0 前に直す」
    The Result of such a Session opens in the mode it was drawn in, which
    `0.13.0` recorded only as the mode of the last **Generate**
    (`cliInvocation.mode`).
  - Reply 8 (`2026-10-16`), to the question
    「OV-677 の修正で、PD-OI-096 の受領文(Product Contract)を scenario revision 2 に改めます。変わるのは「Accepted residual risk」だけで、他の欄はすべて今のままです。どちらの文で承認しますか？」
    with the options 「案1 終わった暫定の文も消す (Recommended)」 and
    「案2 OV-677 の文だけ直す」:
    「案1 終わった暫定の文も消す (Recommended)」
    The preview of each option showed the new and the old Accepted residual
    risk in full and said that every other field stays, with Scenario
    revision `2` and the approval date as the decision date. The other fields
    were not shown again. The receipt recorded under `PD-OI-096` below is the
    scenario revision `1` receipt with that Accepted residual risk, Scenario
    revision `2`, and Decision date `2026-10-16` (SHA-256, UTF-8, excluding
    the final newline:
    `abb6fe7056b5b15c82a883c3b5e3d9f13892a472a499e9fbca3986cdda2f6a87`).
  - Scope: only the Accepted residual risk changes. It drops the two
    sentences on the interim before the Load confirmation dialog, which
    Revision 52 brought (`PD-OI-099`; its entry says "`PD-OI-096` keeps its
    text; the confirmation dialog its Accepted residual risk names is the
    dialog of this record"). Its sentence on a Legend that `0.13.0` drew in a
    Linear corner no longer applies only to Sessions saved while Linear was
    shown, and it names the Sessions whose Result's mode is not recorded. The
    scenario revision `1` receipt is not kept in the record; Git history
    retains it. `PD-OI-090` keeps its text (its 〔以前の Session の移行〕
    states the outcome of Reply 7). No other record is amended or retired.
  - This change is a static Product Contract co-change: the runtime, tests, and
    documentation that implement it are in the same pull request, and the
    Review is `REQUIRED`.
- Revision 56 changes: `PD-OI-101` is added (`O-B3 / TASK-YIELD`, concern
  `web.generate.later-status-paint-yields`). It comes from the Product
  Decision Owner `satoshikawato`'s reply of `2026-10-16` to question
  `S3C-Q1` of the decision packet in campaign `incremental-render-20261013`,
  `design/STEP3C-DESIGN.md` revision 3, section 9 (outcomes and their
  consequences in sections 2 B and 5; `DECISIONS.md` `D-53`).
  - The question: what a Generate does at its two later paint waits, after
    "Preparing render inputs and session..." and "Preparing preview...". The
    outcomes: **O-B3** replace both frame waits by a task yield (recommended),
    **O-B4** one frame instead of two, **O-B1** drop the yield after
    "Preparing render inputs and session...", **O-B2** drop both, **O-B0**
    keep both.
  - Reply:
    「S3C-Q1: B3」
  - The reply names the outcome. The field texts of the receipt below are the
    O-B3 parts of the receipt fields that the packet proposed (its section 5),
    with nothing added. The receipt takes effect with the Owner's approval on
    the pull request that adds this revision; that pull request does not merge
    without it.
  - Scope: only the two later status paint waits of a Generate. The paint
    wait after the first status and the post-bind frame are unchanged.
    `PD-OI-016` and every other record keep their text; no record is amended
    or retired.
  - This change is a static Product Contract co-change: the runtime and tests
    that implement it are in the same pull request, and the Review is
    `REQUIRED`.
- Revision 57 changes: `PD-OI-102`, `PD-OI-103`, and `PD-OI-104` are added,
  and `PD-OI-066` and `PD-OI-059` are replaced for scenario revision `2`. They
  come from the Product Decision Owner `satoshikawato`'s replies of
  `2026-10-17` to the questions of the decision packet in campaign
  `incremental-render-20261013`, `design/STEP4-DECISION-PACKET.md` revision
  2.2 (design `design/STEP4-DESIGN.md` revision 2.2; `DECISIONS.md` `D-67`,
  `D-68`, and `D-70`). The questions were asked after the Owner invited questions.
  - Reply 1, one selection per question; each question and the selected label:
    - `S4-Q1` 「いつ反映しますか？」: 「A 確定時に裏で再描画 (Recommended)」
    - `S4-Q2` 「どのセクションを即時にしますか？」: 「A G1 だけ (Recommended)」
    - `S4-Q3` 「GC・skew・depth の表示切替をどうしますか？（OV-1041 と OV-1043 もこの選択に従う）」:
      「B すぐ隠す＋正確な注記 (Recommended)」
    - `S4-Q5` 「既存の即時編集が設定反映の終わりを待つことを認めますか？」:
      「A G1 なら約 0.9 s まで (Recommended)」
    - the `PD-OI-059` revision `2` question 「定義のどの設定を即時にしますか？」:
      「A 題名だけ即時 (Recommended)」
  - Reply 2, to the question
    「receipt 下書き 7 件（Stage D DL、S4-Q1 A、S4-Q2 A、S4-Q3 B、S4-Q4 A、S4-Q5 A、PD-OI-059 rev 2 A、Owner: Satoshi Kawato、2026-10-17）をこのまま receipt として確定してよいですか？」:
    「このまま確定 (Recommended)」
    The drafts filled Rationale, Must preserve, May retire, and Accepted
    residual risk from the packet text of each selected option; Reply 2
    confirmed them as written.
  - Reply 3, two corrections of the confirmed text after a review found that
    it repeated errors of the packet; each question and the selected label:
    - 「S4-Q3 B の receipt（PD-OI-104）に事実の誤りがありました。「表示に戻すのは Generate で」と書きましたが、今の dev では描画済みのトラックは隠すのも戻すのも即時で、注記もそう書いてあります（packet の記述が原因）。receipt をどうしますか？」:
      「今の動作どおりに直す (Recommended)」
      It corrects the Rationale and Must preserve of the `S4-Q3` receipt: a
      drawn track is hidden or shown at once, and a track the Result lacks
      and the layout apply on Generate, with the two notes the Web app shows
      today.
    - 「S4-Q4 の効果文に、選ばなかった S4-Q1 C の仕組み「a section's Apply now」が入っています。文を直しますか？」:
      「選んだ仕組みだけにする (Recommended)」
      It corrects the `S4-Q4` receipt, which is not recorded here.
  - The receipts below hold the confirmed text with the Reply 3 correction.
    Wrapped lines of the confirmed text are joined; no other character
    differs.
  - Receipt mapping: `S4-Q1` → `PD-OI-102`; `S4-Q2` → `PD-OI-103`; `S4-Q3` →
    `PD-OI-104`; `S4-Q5` → `PD-OI-066` scenario revision `2`; the `PD-OI-059`
    question → `PD-OI-059` scenario revision `2`. The receipts name the
    questions by these IDs and the setting groups by the packet's names
    (G1, G2, G3, GT), which `PD-OI-103` quotes. `S4-Q4` (concern
    `product.saved-session-regeneration-continuity`) belongs to the Product
    Impact map, and the stage D receipt of the same
    reply (concern `architecture.r1.generate-preview-mount`) keeps today's
    behaviour; neither is recorded here. `PD-OI-060` was not asked and keeps
    scenario revision `1`.
  - Scope: the `PD-OI-066` scenario revision `1` receipt is kept as it was
    approved; scenario revision `2` changes only what its receipt states. The
    `PD-OI-059` scenario revision `1` receipt is not kept in the record; Git
    history retains it. The `PD-OI-059` record title changes to state
    scenario revision `2` (no link targets its anchor); its concern key keeps
    its wording. `OIC-027` is reworded to the parity of `PD-OI-066`
    scenario revision `2`. `PD-OI-037`, `PD-OI-052`, `PD-OI-060`,
    `PD-OI-084`, `PD-OI-085`, and every other record keep their text; no
    other record is amended or retired.
  - This authority-only amendment contains no runtime; dependent runtime
    requires it merged into its base.
- Revision 58 change: `PD-OI-039` is replaced for scenario revision `4`
  (`A / DOCKED-REVIEW-AT-EVERY-WIDTH`, the same outcome) from the Product
  Decision Owner `satoshikawato`'s replies of `2026-10-12` and `2026-10-13` on
  fitting the figure again during Review alignment in the Web app. Each
  quotation below is the part of a reply that carries the answer; what it
  leaves out is stated.
  - Reply 1 (`2026-10-12`), to a page that asked, among other things, what
    the figure does when its room changes while Review alignment is open:
    `A`, recommended, fit the figure again whenever the room of the Result
    changes, with the stated cost that a zoom or position changed by hand
    during the review is replaced at such a change; `B`, on a display switch
    only; or `C`, leave it. 「推奨通りでいきます。Product Decisionも承認します」
    The reply selects `A`. It does not word the receipt.
  - Reply 2 (`2026-10-13`), to a page of Before/After images of the built
    behaviour (nine causes of a room change, each from the fit of the
    opening and after a zoom made by hand) and the request to approve the
    receipt recorded below as written or to name what to change. Section 6
    of that page showed the fields in full with the passages added to
    scenario revision `3` marked.
    「OV-536（上のページの 6）: そのまま承認」
    The shown text had a placeholder for the decision date; the date recorded
    is the day of this reply. The receipt is recorded byte for byte.
  - Scope: scenario revision `3` (same outcome, receipt SHA-256
    `3d1ec79f071f7f6f5ea1c407853d8550df2491085954f253f1b00aa9cb62eb44`) is
    retained in Git history, and what it retired stays retired. Every
    sentence of its receipt is kept; the receipt below adds one sentence to
    the Rationale, the passage `〔review 中の Fit のし直し〕` to Must preserve,
    one sentence to May retire, and the sentences after 「…区別する必要がある。」
    to Accepted residual risk. `PD-OI-035` keeps its text; only its reference
    to the governing `PD-OI-039` scenario revision now reads `4`.
  - `PD-OI-093` keeps its receipt byte for byte. The last sentences of its
    Accepted residual risk
    (「Review alignment を開いたまま表示を切り替えるか section を開閉すると、図は自動では合わせ直されない（窓の幅を変えたときの今の動きと同じ）。」)
    described the Web app until this change, as its
    Revision 44 entry announced; from this revision the figure is fitted
    again in those cases and that risk no longer exists.
  - No other record is amended or retired.
  - This change is a static Product Contract co-change: the runtime and tests
    that implement it are in the same pull request, and the Review is
    `REQUIRED`.
- Revision 59 change: `PD-OI-053` is replaced for scenario revision `2`
  (`A / EXPLICIT-MODE-WITH-DISCOVERABLE-TARGETS`, the same outcome) from the
  Product Decision Owner `satoshikawato`'s replies of `2026-10-13` on the
  pointer cursor of a Legend row in the Web app. Each quotation below is the
  part of a reply that carries the answer; what it leaves out is stated.
  - The finding: a Legend row (its colour swatch and its name) showed the
    pointing-hand cursor in every state, where scenario revision `1` names
    the help cursor with Layout edit off and grab / grabbing with it on; and
    with Layout edit off a drag that starts on a row did not move the canvas.
    A row is a button that opens its popup.
  - Reply 1 (`2026-10-13`), to two outcomes: `(a)`, recommended, correct only
    the state with Layout edit on (grab, and grabbing during the drag) and
    record that with Layout edit off a row stays a button; or `(b)`, make the
    product follow scenario revision `1` in both states.
    「OV-603：a)（推奨）：ON のときだけ直して「つかめる手／つかんだ手」にします。」
    The reply's other line approved the receipt text of `PD-OI-094`
    (Revision 46).
  - Reply 2 (`2026-10-13`), to the request to approve the receipt recorded
    below as written or to name what to change. The request showed the
    receipt in full and named the parts added to scenario revision `1`.
    「OV-603（前のメッセージの文面）: そのまま承認」
    The reply's two other lines approved the receipts of the other two
    records of the same campaign, a click on a drawn element that opens its
    settings section and the refit during Review alignment, which are
    recorded in their own revisions.
    The shown text had a placeholder for the decision date; the date recorded
    is the day of this reply. The receipt is recorded byte for byte.
  - The same request offered one more sentence, which would have recorded
    that a record definition of a Linear figure is, with Layout edit off, a
    button like a Legend row. The reply does not take it up, so the receipt
    does not contain it and this record says nothing about record
    definitions.
  - Scope: scenario revision `1` (same outcome, receipt SHA-256
    `0365e0057b61e6b606c6714dc17f0e5c89d755ce651c1a8eec66d0139e0c65ec`) is retained in Git history. Every sentence of its receipt is
    kept; the receipt below adds one sentence to the Rationale, the two
    `〔凡例の項目…〕` passages to Must preserve, and one sentence each to May
    retire and Accepted residual risk. For a Legend row with Layout edit off,
    the help cursor and the canvas pan from a drag that starts on the row are
    retired; for every other supported target scenario revision `1` holds
    unchanged.
  - No other record is amended or retired.
  - This change is a static Product Contract co-change: the runtime and tests
    that implement it are in the same pull request, and the Review is
    `REQUIRED`.
- Revision 60 changes: `PD-OI-105` is added
  (`B / CLICK-DRAWN-ELEMENT-SHOWS-ITS-SECTION-IN-POPUP`, concern
  `web.result.element-settings-entry`). It comes from the Product Decision
  Owner `satoshikawato`'s replies of `2026-10-12`, `2026-10-13`, and
  `2026-10-18` on reaching the settings of a drawn element from the figure in
  the Web app. Each quotation below is the part of a reply that carries the
  answer.
  - Reply 1 (`2026-10-12`), to mock images in which a click on the Circular
    scale ticks, the GC skew track, the Linear scale, or the plot title shows
    that element's settings section beside the figure.
    「図の要素をクリックして設定を開く案(D、1 回目の最後)を足す。」
  - Reply 2 (`2026-10-13`), to a page of Before/After images of a first built
    form, in which the section opened in the settings pane. The page asked
    three choices: where the settings open (`A` in the settings pane, built;
    `B` beside the figure, as in the mock images), what the plot title and
    the Linear scale say on hover (both sentences; the present sentence only;
    or not clickable), and the English wording **Click to open its settings**.
    「PR-C（図の要素から設定を開く）: 図の横に浮かべる（最初の案の画像の形）。設定はポップアップする方がいいと思う。選んでいただく点 2：タイトルと Linear の scale の説明: 2文とも出す。選んでいただく点 3：英語の文言：このままでよい。」
  - Reply 3 (`2026-10-13`), to a page of Before/After images of the rebuilt
    form `B` and the request to approve its receipt text as written or to
    name what to change.
    「PR-C（https://claude.ai/artifact/RDBmSKD27nNQ6Xd21ivmrb の 12）: そのまま承認」
    Section 12 of that page showed the nine fields of the receipt as recorded
    below, the Must preserve text as a list whose `〔…〕` headings were set in
    bold without the brackets, and the decision date as the day of the
    approval. The receipt is recorded byte for byte.
  - Reply 4 (`2026-10-18`), to the question what the receipt's
    「値は今までどおり Generate で反映される。」 means now that `PD-OI-102`,
    `PD-OI-103`, and `PD-OI-104` make some rows of these sections apply at
    once: `B`, recommended, keep the receipt as approved and record the
    reading in the item on those records below; or `A`, amend that sentence
    to say that a value applies when the same change in the settings pane
    applies. 「D-41: B、D-44: B」
    「D-41: B」 selects `B`. 「D-44: B」 answers a question on the record of
    Review alignment, recorded in its own revision.
  - `PD-OI-053` keeps its text. With Layout edit off, the plot title and the
    Linear scale keep its help cursor, dashed frame, and "Turn on Layout edit
    to move this item"; the receipt below adds a second line to that hint
    (`〔押せることの表示〕`). With Layout edit on, nothing of this record applies
    (`〔出ないとき〕`).
  - `PD-OI-039` keeps its text: while Review alignment is open a click on the
    figure opens, switches, and closes nothing, this popup included
    (`〔出ないとき〕`, `〔閉じるとき〕`).
  - `PD-OI-092` and `PD-OI-093` keep their text: the section index and the
    icon column work as recorded, and using either on the section that the
    popup holds returns the section to the settings pane first
    (`〔閉じるとき〕`, `〔アイコンの列〕`).
  - `PD-OI-102`, `PD-OI-103`, and `PD-OI-104` keep their text. The receipt
    below was approved on `2026-10-13`, before they were added (Revision
    57). As Reply 4 selects, the popup holds the section itself, so a row
    of that section applies when those records say it applies, and the
    receipt's 「値は今までどおり Generate で反映される」 does not change when
    a row applies.
  - No other record is amended or retired.
  - This change is a static Product Contract co-change: the runtime, tests, and
    documentation that implement it are in the same pull request, and the
    Review is `REQUIRED`.
- Revision 61 change: `PD-OI-039` keeps scenario revision `4` and its
  receipt byte for byte (SHA-256
  `8f755440a6655699b1341600bd2aba0d17ef80172490fe600a7cf25d5e16b899`). The
  Product Decision Owner `satoshikawato`'s reply of `2026-10-18` reads one
  sentence of its Accepted residual risk. The quotation below is the part of
  the reply that carries the answer; what it leaves out is stated.
  - The reply (`2026-10-18`), to the question whether an Undo or Redo of a
    change that applies at once (`PD-OI-102`, such as a feature colour)
    replaces the Result in the Accepted residual risk sentence
    「reviewを開いたままSessionの読み込み、Undo、CircularとLinearの切り替えでResultが置き換わったあとは、panelが開いていても図はFitし直されない（今と同じ）。」:
    `B`, recommended, it does not: the review still shows the figure it
    fitted, keeps its rings and the view from before the review, and fits
    the figure again when the room changes; or `A`, after any Undo or Redo
    the figure is not fitted again. 「D-41: B、D-44: B」
    「D-44: B」 selects `B`. A Session load, a Generate, and an Undo or Redo
    across a Generate still replace the Result as that sentence says.
    「D-41: B」 answers a question on the record of a click on a drawn
    element that opens its settings section, recorded in its own revision.
  - The Decision source of `PD-OI-039` now names this reply too. No other
    record is amended or retired.
  - This change is a static Product Contract co-change: the runtime and tests
    that keep the rings and the held view across a change that applies at
    once are in the same pull request, and the Review is `REQUIRED`.
- Revision 62 changes: `PD-OI-090` is replaced for scenario revision `2`
  (concern `session.project-drawings`; the outcome
  `A / PROJECT-WITH-DRAWINGS` is kept), and `PD-OI-044` is replaced for
  scenario revision `3` (`A / REVEAL_APPLICABLE_SINGLE_RECORD_CONTROLS` is
  kept). They come from two replies of the Product Decision Owner
  `satoshikawato`. Each quotation below is the part of a reply that carries
  the answer.
  - Session version: Reply 3 of the Revision 49 entry (`2026-10-13`)
    accepted the reading 「そちらを 48 に読み替える」 of the drawings record
    (「P-8（Session の版）: 推奨通り」). Session version `47` is now the Legend
    row keys format that `main` and `dev` write, so the project and drawing
    format is Session version `48`, and version `47` Sessions are an earlier
    format that moves to it without Python, as version `46` Sessions do.
    Every Session version in the two receipts that named the drawings format
    reads accordingly; every other number is unchanged. In `PD-OI-090`:
    「version 46 と 47 の Session の preview だけの Load は Python を起動しない（PD-OI-044 revision 2）」
    reads
    「version 46、47 と 48 の Session の preview だけの Load は Python を起動しない（PD-OI-044 revision 3）」;
    「version 46（0.14.0）は Python を起動せずに移し」 reads
    「version 46 と 47（0.14.0 は 47）は Python を起動せずに移し」;
    「version 46 の modes と otherModeResult。version 46 は読むだけ」 reads
    「version 46 と 47 の modes と otherModeResult。読むだけ」; and
    「新しい形式（version 47）」 reads 「新しい形式（version 48）」. In
    `PD-OI-044`: 「（version 47、PD-OI-090）」 reads
    「（version 48、PD-OI-090）」.
  - Drawing switch: the reply (`2026-10-18`) to the question
    「OV-1100（Product の判断）: Rules 表の編集をルールと照合している途中に図を切り替えたとき、どうするか。」
    with the choices 「A: 編集を捨てる。」, 「B: 照合が終わるまで切り替えを止める（推奨）。」, and
    「C: 編集を裏で元の図に入れる。」:
    「OV-1100（Product の判断）:B: 照合が終わるまで切り替えを止める」
    It selects `B / SWITCH-WAITS-FOR-RULE-COMMIT`. `PD-OI-090` adds one
    sentence after 「Generate と Session の保存・読み込みの間は切り替えられない（PD-OI-045）。」:
    「色ルールの表（Colors の SPECIFIC RULES）の編集が照合を待っている間も切り替えられず、照合が終わるとその編集はその drawing の History の 1 つの step になる。」
    The reply states no rationale, preserved or retired effect, or accepted
    risk for this choice, and the receipt adds none.
  - Receipt approval: the reply (`2026-10-18`) to the two receipt texts
    above, shown as written with the `OIC-030` wording: 「C2 承認」. It
    approves both receipts as recorded.
  - New receipt SHA-256 values: `PD-OI-090`
    `a69007c6ea265135e50df20e3780ffa78b5bd3cab035603bb05dba96c685484e`;
    `PD-OI-044`
    `5a26a1e530c80e9cfaa44c7f39f6971bcb95b90e94adea1b0ffd0ccaa6e0a813`.
    Both receipts carry the decision date `2026-10-18`.
  - `OIC-030` reads 「each Session 46 and 47 fixture」 and 「a Session 46,
    47, or 48 preview-only Load」, and adds after 「a drawing switch writes no
    setting and is not an Undo step;」 the obligation 「a drawing switch is
    refused while a color rule table edit is being matched, and that edit
    becomes one History step of its drawing;」.
  - `PD-OI-091` cites 「PD-OI-044 revision 2 の May retire」; scenario
    revision `3` keeps that May retire byte for byte, so the citation keeps
    its meaning. No other record is amended or retired.
  - This change is Contract-only: the Contract is the only changed path, and
    the Review is `REQUIRED`. The runtime that implements it follows in
    separate pull requests, which cite this authority on their base: the
    drawing switch in DRAWINGS-A (step A3b), and the Session `48` writer in
    DRAWINGS-B (step X).
- Revision 63 change: `PD-OI-090` is replaced for scenario revision `3`
  (concern `session.project-drawings`; the outcome
  `A / PROJECT-WITH-DRAWINGS` is kept), and `OIC-030` is amended to match.
  They come from two replies of the Product Decision Owner `satoshikawato`.
  Each quotation below is the part of a reply that carries the answer.
  - Project steps: the reply (`2026-10-18`) to the questions
    「Q-DR-5 OV-1160 Instant Preview applies only the shown drawing's queued palette (dev too): A turning it on applies every drawing's queue in the toggle's one History step (Recommended) / B drops the hidden drawings' queues / C keep.」
    and
    「Q-DR-8 (A3 review M2, intent form) Which History holds edits of Session-level values (Session title, Auto Reflow, DPI, Instant Preview, LOSAT execution, drafts)? PD-OI-090 says Undo never undoes another drawing's operation but does not place these. A the project stream, undone newest-first together with the shown drawing's steps (Recommended; Circular rename A, Linear rename B, Undo anywhere undoes B first) / B the drawing shown when the edit was made (A3 today: Undo in Circular can drop Linear's later rename; conflicts with the Contract's intent) / C not undoable.」:
    「Q-DR-5（OV-1160）: A」 and 「Q-DR-8: A」.
    Q-DR-8 `A` makes an edit of a value that does not depend on a drawing a
    project step, undone newest first together with the shown drawing's
    steps; Q-DR-5 `A` applies every drawing's queued palette in the one step
    that turns Instant Preview on. Its answers to Q-DR-4, Q-DR-6, and
    Q-DR-7 fall within the existing records and change no text.
  - Owner-delegated parts: two choices taken under the Owner's standing
    instruction (`2026-09-29`, restated `2026-10-07`) are written into the
    receipt, and the approval below covers them. D-15-51: Reset all drawings
    is a project operation, as 〔Reset〕「実行は 1 回の Undo で戻る」 requires.
    D-15-53: an operation that also changes a drawing value is one step of
    that drawing, as 「1 操作 1 step」 and 「他の drawing の操作は戻さない」
    require.
  - Receipt change: in 〔Undo と削除〕 of `PD-OI-090`,
    「Undo/Redo は表示中の drawing の操作と、drawing の追加、複製、名前の変更、並べ替え、削除を、新しい順に 1 操作 1 step で戻し、他の drawing の操作は戻さない。」
    reads
    「Undo/Redo は表示中の drawing の操作と project の操作を、新しい順に 1 操作 1 step で戻し、他の drawing の操作は戻さない。project の操作は、drawing の追加、複製、名前の変更、並べ替え、削除、Reset all drawings と、drawing によらず 1 つの値（Session の title、Auto Reflow、PNG DPI、palette の Instant Preview、LOSAT の実行方法とスレッド数、編集の入力欄に書きかけた値）だけを変える操作で、どの drawing を表示していても戻せる。drawing の値も一緒に変える操作（色ルールを追加して入力欄の書きかけが消えるなど）は、その drawing の操作として 1 step で戻る。Rich Feature Popup は今までどおり Undo/Redo では変わらない（PD-OI-100）。palette の Instant Preview を On にすると、すべての drawing で反映を待っている palette をその 1 step で反映し、その Undo で待っている状態に戻る。」.
    Every other term of scenario revision `2` is unchanged. `PD-OI-100`
    keeps its text; the receipt cites it.
  - Receipt approval: the reply (`2026-10-18`) to the receipt text above,
    shown as written with the `OIC-030` wording: 「C3（Contract の文面）の承認: 承認」.
    It approves the receipt as recorded.
  - The replies state no rationale, preserved or retired effect, or accepted
    risk, and the receipt adds none.
  - New receipt SHA-256: `PD-OI-090`
    `e3c534983538a18884e9e5e3400657ccdb14e6c0b75679b4e5eeaa4bc040e795`.
    The receipt keeps the decision date `2026-10-18`.
  - `OIC-030` reads 「Undo/Redo of a step in one drawing」 in place of
    「Undo/Redo in one drawing」, so that the Undo of Reset all drawings,
    which restores every drawing, keeps the obligation, and adds after
    「Reset all drawings is the only action that resets more than the shown drawing;」
    the obligation
    「Undo/Redo in any drawing walks the shown drawing's steps and the project steps (drawing add, duplicate, rename, reorder, delete, Reset all drawings, and edits that change only the Session title, Auto Reflow, PNG DPI, palette Instant Preview, LOSAT execution and threads, or input drafts) newest first, one step per operation; an operation that also changes a drawing value is one step of that drawing; turning Instant Preview on applies every drawing's queued palette in that one step;」.
  - No other record is amended or retired.
  - This change is Contract-only: the Contract is the only changed path, and
    the Review is `REQUIRED`. The runtime that implements it follows in
    DRAWINGS-A (the project steps of History), which cites this authority on
    its base.
- Revision 64 changes: `PD-OI-044` is replaced for scenario revision `4`
  (`A / REVEAL_APPLICABLE_SINGLE_RECORD_CONTROLS` is kept) and `PD-OI-090`
  for scenario revision `4` (`A / PROJECT-WITH-DRAWINGS` is kept), and
  `OIC-030` is amended to match; `PD-OI-103` and `PD-OI-104` are replaced for
  scenario revision `2`, `PD-OI-066` and `PD-OI-059` for scenario revision
  `3`, and `PD-OI-060` for scenario revision `2`. They come from replies of
  the Product Decision Owner `satoshikawato` of `2026-10-18` in campaign
  `incremental-render-20261013` (`DECISIONS.md` `D-128` to `D-134`), asked
  after the Owner invited questions (「質問して」). Each quotation below is the
  part of a reply that carries the answer; each question is given as
  `DECISIONS.md` records it.
  - Session version (`D-129`, OV-1178): the reply to the question "re-read
    Contract revision 62's Session numbers now (0.14.0 = Session 48,
    drawings 49, the Python-free preview Load covers 46-49; Contract-only PR
    before the tag) or leave it to the drawings campaign":
    「A 読み替える（推奨）」. It accepts that reading of the drawings record.
    Session version `48` is the format the 0.14.0 runtime writes, which
    stores protein comparison results in the search frame (OV-1145, D-124).
    So the project and drawing format is Session version `49`, and version
    `48` Sessions are an earlier format that moves to it without Python, as
    version `46` and `47` Sessions do. Every Session version in the two
    receipts that named the drawings format reads accordingly; every other
    number is unchanged.
  - Load conversion (`D-128`, OV-1179): D-126 (`2026-10-18`) replied
    「A. Python Worker 経由」 to a question that said "no authority PR". The
    reply to Q-O2, "B-LOAD (keep D-126; a Session 46/47 with protein results
    stored on a reversed record shows its stored figure at once, the Python
    Worker converts once at Load; PD-OI-044 rev 4 / PD-OI-090 rev 3 rows) or
    A-SAVE (Load stays Python-free; conversion at the first Generate or
    Save)", is: 「B-LOAD（推奨）」.
    `PD-OI-044` and `PD-OI-090` therefore let a Session up to 47 start
    Python once at Load, but only when it holds protein comparison results
    stored on a reversed record. The figure shows first, as revision 3
    requires for version 27-44. The replies state no rationale, and the
    receipts add none.
  - Step 5 (`D-130` to `D-133`): the questions of the decision packet
    `design/STEP5-G2GT-DECISION-PACKET.md` revision 3 of that campaign (base
    `332e3ad5f`); each question and the reply:
    - `S5-Q1` (`web.settings.live-sections` revision 2), how G2 sections go
      live (A evidence rule / B wave receipts / C all after SW1 / D stay on
      Generate), with the reading (MARGIN / SENSITIVITY) and the evidence
      unit (U1 / U2): 「A・MARGIN・U2（推奨）」
    - `S5-Q2` (`web.layout.track-toggle-application` revision 2), GC / skew
      / Depth toggles (A redraw when within target / B redraw now), reading
      (MARGIN / SENSITIVITY / AS S5-Q1), set (TOGETHER / PER-TOGGLE):
      「A・S5-Q1 と同じ読み方・まとめて（推奨）」
    - `S5-Q3` (`web.live-edit.regeneration-parity` revision 3), how long an
      editor edit may wait for live applies (A about two applies within the
      rule / B the measured cost): 「A 目標内の 2 回分（推奨）」
    - the packet's conditional receipts `web.circular.definition-application`
      revision 3 (A live with S5-Q1 / B keep revision 2) and
      `web.stroke.application` revision 2 (A live with S5-Q1 / B keep
      revision 1): 「A / A（推奨）」
  - Receipt text (`D-134`): the reply to the question of the Rationale, Must
    preserve, May retire, and Accepted residual risk fields of `D-128` to
    `D-133`: 「下書きをそのまま採用」. It adopts the drafts as written:
    `OWNER-RECEIPTS-S18-DRAFT.md` revision 2.1 for `PD-OI-044` and
    `PD-OI-090` (its rows for both variants and for B-LOAD; no A-SAVE row),
    and for the five step 5 receipts the packet's text of the chosen option:
    Must preserve, May retire, and Accepted residual risk are the option's
    Preserved, Retired, and Residual risk text, the `S5-Q1` and `S5-Q2`
    receipts add to Must preserve the packet's "Evidence rule" text, and each
    Rationale quotes the packet's stated reason for the option. Where the
    packet states a field as "as S5-Q1 A for the unit", the receipt says so.
    The serialized receipts are shown to the Owner for review before this
    change merges.
  - Numbering: `PD-OI-090` scenario revision `3` (Revision 63) reached the
    base after the draft was written against scenario revision `2`. The
    drafted scenario revision `3` is recorded as scenario revision `4`,
    superseding scenario revision `3`; its text changes are applied to
    scenario revision `3`'s text, whose 〔Undo と削除〕 is kept.
  - Receipt mapping: `S5-Q1` → `PD-OI-103` scenario revision `2`; `S5-Q2` →
    `PD-OI-104` scenario revision `2`, which supersedes the `S4-Q3` outcome;
    `S5-Q3` → `PD-OI-066` scenario revision `3`;
    `web.circular.definition-application` revision 3 → `PD-OI-059` scenario
    revision `3`; `web.stroke.application` revision 2 → `PD-OI-060`
    scenario revision `2`. Each Choice names the packet's option of that ID.
    The receipts name the setting groups (G1, G2, G3, GT) as `PD-OI-103`
    lists them and the stage work as the packet does (SW1, SW2, SW3).
  - New receipt SHA-256 values; every receipt carries the decision date
    `2026-10-18`:
    - `PD-OI-044` scenario revision `4`:
      `896078be8ce11203eedaddae8d048ae312bceb46f4f89f63f19ff14baadcb398`
    - `PD-OI-090` scenario revision `4`:
      `dbc29cc437fd62eb2aad80cbbe7872f7ab701896b9a7ccf54609ca5b9c168f81`
    - `PD-OI-103` scenario revision `2`:
      `95f13d32a53bccd0b4997663e78ff5acb731971220cc3fa36e8ba426267ab6ca`
    - `PD-OI-104` scenario revision `2`:
      `55b1a9d81c4eada55480cc440f892fb6d4058e6ead85dcbbf681cf67ccfcbaaa`
    - `PD-OI-066` scenario revision `3`:
      `d051f8b4a669d2be36d0f91c06867d8c4079fdd5605451583d3029ccd3b150a7`
    - `PD-OI-059` scenario revision `3`:
      `9d9f749fff0f9a8f24a74d912b707a21712649d674c0a06f183cb0ef90ffc4cf`
    - `PD-OI-060` scenario revision `2`:
      `c11a661c2922c3bda0a37c77a94dc0015dc4a9323bad86ef78f4c9aefb20ae0a`
  - `OIC-030` reads 「each Session 46, 47, and 48 fixture」 and 「a Session 46,
    47, 48, or 49 preview-only Load starts no Python Worker, except that a
    Session 46 or 47 holding protein comparison results stored on a reversed
    record shows its stored figure first and then starts it once to convert
    them」.
  - `PD-OI-091` cites 「PD-OI-044 revision 2 の May retire」. Scenario
    revision `4` keeps that May retire byte for byte, as revision `3` did, so
    the citation keeps its meaning.
  - Scope: `PD-OI-044` and `PD-OI-090` keep only their new receipts; Git
    history retains the former ones. `PD-OI-103`, `PD-OI-104`, `PD-OI-059`,
    and `PD-OI-060` keep their former receipt, which holds for each unit or
    toggle until it goes live, and `PD-OI-066` keeps its scenario revision
    `1` and `2` receipts. The titles of `PD-OI-060` and `PD-OI-104` change
    to state when the settings apply (no link targets their anchors).
    `PD-OI-102`, `OIC-024`, `OIC-027`, and every other record keep their
    text, and no other record is amended or retired. The R1 exception
    sentence of `gbdraw/web/CLAUDE.md` that cites `PD-OI-104` is unchanged
    here; the packet drops it in a later authority-only pull request after
    the GT runtime change.
  - This change is Contract-only: the Contract is the only changed path, and
    the Review is `REQUIRED`. The runtime that implements it follows in
    separate pull requests, which cite this authority on their base: the
    Web Load conversion of OV-1145 (its step 2), and the step 5 work after
    the 0.14.0 tag and DRAWINGS A1.
- Revision 65 change: `PD-OI-104` is replaced for scenario revision `3`
  (`A / TOGGLE-LIVE-REDRAW-WHEN-WITHIN-TARGET`, reading `AS S5-Q1`, set
  `TOGETHER`, is kept). It comes from the Product Decision Owner
  `satoshikawato`'s reply of `2026-10-18` (campaign
  `incremental-render-20261013`, `DECISIONS.md` `D-135`) after Revision 64
  merged in pull request #1193.
  - The reply, to the question whether the May retire of `PD-OI-104`
    scenario revision `2`, serialized from the packet's "Retired" text only,
    also holds the "Lost" text of `S5-Q2` A in
    `design/STEP5-G2GT-DECISION-PACKET.md` revision 3 (the instant hide and
    show, replaced by a redraw within the rule; Show Depth turned on after a
    Depth source change no longer shows the old drawn track at once):
    「足す（推奨）」.
  - Receipt change: May retire
    「`applyTrackVisibility`, its watcher, its History restore call (`app/app-setup.js` 3904), the R1 exception sentence.」
    reads
    「`applyTrackVisibility`, its watcher, its History restore call (`app/app-setup.js` 3904), the R1 exception sentence; the instant hide and show (a redraw within the rule instead); Show Depth turned on after a Depth source change no longer shows the old drawn track at once.」.
    Every other field of the scenario revision `2` receipt is unchanged; the
    scenario revision `1` receipt stays in the record, as before.
  - The reply states no rationale or accepted risk for this change, and the
    receipt adds none.
  - New receipt SHA-256: `PD-OI-104` scenario revision `3`
    `8d2096e0ea8236387807688b956388ce96e13f005ada0a8e08f130364137e3cd`.
    The receipt keeps the decision date `2026-10-18`.
  - No other record is amended or retired.
  - This change is Contract-only: the Contract is the only changed path, and
    the Review is `REQUIRED`.
- Revision 66 changes: `PD-OI-097` is replaced for scenario revision `9`
  (concern `product.circular-outer-label-fit`). It comes from the Product
  Decision Owner `satoshikawato`'s replies of `2026-10-18` and `2026-10-19`
  on Linear region (annotation) labels that the placer cannot move clear of
  feature labels and leader lines (OV-1155). Each quotation below is the part
  of a reply that carries the answer; Replies 1 to 12 are those of the
  Revision 50 and 53 entries.
  - Reply 13 (scenario revision `9`, `2026-10-18`), to the question
    「OV-1155: Linear の領域(annotation)ラベルを、自分の範囲の中でどこへ動かしても引き出し線や feature ラベルを避けられないとき、今は重ねたまま、警告なしで描いています。どうしますか？」
    the selected option and its description:
    「警告に載せる (推奨)」
    「図は今のまま描き、既存の Linear の収まり警告(CLI の警告、Web の通知帯)に「領域ラベルの重なり N 件」を加え、直し方の選択肢も出します。利点: LP-D42・LP-D48 の「収まらなければ警告して選んでもらう」と同じ扱いになり、図が勝手に変わらない。欠点: 通知文と選択肢を足す必要がある(Contract の変更を伴うかもしれない)。重なり自体は残る。」
    (`REGION-WARN`)
  - Reply 14 (`2026-10-18`), the message
    「cli-hmmtdna-20.png てかこれ"rRNA genes"のラベル下に置けばよくね？紫で塗られた空白の領域が空いてるじゃん。ここにうまく配置すれば、unambiguousだしきれいだしフォントサイズも変えなくていい。」
  - Reply 15 (`2026-10-18`), to the question
    「crowded のように、2行にしても重なりが残る（別のラベルにかぶる）とき、どうしますか？」
    the selected option and its description:
    「全部解消するときだけ (推奨)」
    「そのラベルの重なりが 0 になるときだけ枠の中や 2 行を使う。解消しなければ今の 1 行のままにして警告と「領域ラベルを小さく」に任せる。予測しやすく、重なりを新しく作らない。」
  - Reply 16 (`2026-10-18`), the message
    「二行とフォントサイズ小さくするのやめたい。ほかのラベルに動いてもらうのも考えたら？」
  - Reply 17 (`2026-10-18`), the message
    「ラベルの背景を半透明の塗りつぶしにして、レイヤーの前後関係を明確にすれば，フィーチャーラベルのリーダー線の前にregion labelが来れるんじゃないかな」
  - Reply 18 (scenario revision `9`, `2026-10-19`), to the question
    「OV-1155：領域ラベルとほかのラベルが重なるとき、どう直しますか？（上の画像。どの案でも領域ラベル同士の重なりは直らず、警告が残ります）」
    the selected option and its description:
    「(e')+(b) (推奨)」
    「リーダー線とだけ重なる領域ラベルは、白背景つきで線の前に出す。文字同士が重なるときは、そのフィーチャーラベルを同じ側の別の段へ移す（空き段がなければ1段足す）。ラベルは今の側のまま、線はまっすぐ。left_end__e / hmmtdna20__b / combo20__eb」
    (`RGN-FRONT`, `RGN-YIELD-ROW`)
  - Reply 19 (`2026-10-19`), to the question
    「背景の裏を通るリーダー線を「重なり」として数えますか？」
    the selected option and its description:
    「数えない (推奨)」
    「背景で線が隠れるなら重なりとしない。left_end は警告が消える。文字と文字、文字とフィーチャーの図形は、これまでどおり重なりとして数える。」
  - Reply 20 (`2026-10-19`), to the question
    「領域ラベルの背景の見た目はどれにしますか？（backing__left_end / backing__crowded / backing__combo14）」
    the selected option and its description:
    「白の角丸ボックス 70%」
    「背景が薄く、裏の線が透けて見える。」
  - Reply 21 (`2026-10-19`), to the question
    「ほかのラベルを動かすのは自動にしますか？」
    the selected option and its description:
    「自動 (推奨)」
    「重なりが全部消えるときだけ自動で動かす（LP-D77）。消えないときは今の位置のまま警告。2行と「小さく」をやめるので、警告に直し方の候補は出ない（領域ラベル同士の重なりは位置を手で直す）。」
  - Replies the Owner gave on `2026-10-18` that chose region labels on two
    lines and a row that makes region labels smaller were withdrawn by
    Reply 16 before they reached `dev`; they are not part of this record.
  - Replies 13 to 21 state no rationale, preserved or retired effect, or
    accepted risk for these choices, and the change receipt below adds none.
    Replies 14, 15, 17, and 19 to 21 say how the choices are realized (the
    label inside its own filled box; a place used only when it clears every
    conflict of the label; a leader line behind the backing not counted; a
    white rounded box at 70 % opacity; automatic); they add no choice token.
  - Reply 22 (`2026-10-19`), approving the receipt texts below, to the
    question
    「PD-OI-097 rev 9 の受領文（作り直し）：Choice に REGION-WARN ; RGN-FRONT ; RGN-YIELD-ROW を加え、Rationale の末尾に「ただし Linear の領域ラベルの重なりは、図を描いたうえで警告する（LP-D68）。」、May retire の末尾に「（Linear の領域ラベルの重なりを除く）」を足し、Decision date を 2026-10-19 にします。この受領文を承認しますか？」
    the selected option and its description:
    「承認 (推奨)」
    「このとおり Contract に記録し、承認の返答を Reply 22 として引用して SHA-256 を計算する。そのあと PR を開く（マージは Owner）。」
    The record's scenario revision `9` receipt and the change receipt below
    are the approved texts.
    So the record's scenario revision `9` receipt appends
    「ただし Linear の領域ラベルの重なりは、図を描いたうえで警告する（LP-D68）。」 to the Rationale and 「（Linear の領域ラベルの重なりを除く）」 to the May retire of scenario
    revision `8`; its other fields are those of scenario revision `8`.
  - The receipt of this change (SHA-256, UTF-8, excluding the final newline:
    `133d121c2a05096bfa6ed961d40991dba56bac8f41b412e1039d0921539ca27b`):

    ```text
    PRODUCT_DECISION
    Concern: product.circular-outer-label-fit
    Scenario revision: 9
    Choice: REGION-WARN ; RGN-FRONT ; RGN-YIELD-ROW
    Rationale: not stated
    Must preserve: not stated
    May retire: not stated
    Accepted residual risk: not stated
    Owner: Satoshi Kawato
    Decision date: 2026-10-19
    ```

  - Scope: the record's receipt for scenario revision `9` is recorded below;
    the scenario revision `7` and `8` receipts keep their text and their
    SHA-256. `PD-OI-088` applies unchanged. No other record is amended or
    retired.
  - This change is a static Product Contract co-change: the runtime, tests, and
    documentation that implement it are in the same pull request, and the
    Review is `REQUIRED`.
- Revision 67 changes: `PD-OI-097` is replaced for scenario revision `10`
  (concern `product.circular-outer-label-fit`). It comes from the Product
  Decision Owner `satoshikawato`'s replies of `2026-10-19` on Circular region
  (annotation) labels drawn along their circle that are longer than the
  circle allows (OV-1210). Each quotation below is the part of a reply that
  carries the answer; Replies 1 to 22 are those of the Revision 50, 53, and 66
  entries.
  - Reply 23 (scenario revision `10`, `2026-10-19`), to the question
    「Q2（OV-1210）：円のかなりの部分を回るほど長いラベルはどう描きますか？（ページの f2, f3）」
    the selected option and its description:
    「常に曲線＋欠けたら警告」
    「常に曲線で描き、円周に収まらない分は欠ける。欠けたときは警告を出す。」
    and the answer:
    「常に曲線＋欠けたら警告」
    (`RGN-CUT-WARN`)
  - Reply 23 states no rationale, preserved or retired effect, or accepted
    risk for this choice, and the change receipt below adds none.
  - Reply 24 (`2026-10-19`), approving the receipt texts below, to the
    question
    「Q5（OV-1210）：PD-OI-097 rev 10 の受領文です。Choice に RGN-CUT-WARN を加え、Rationale の末尾に「Circular の領域ラベルが円に収まらないときは、収まる分を曲線で描いたうえで警告する（LP-D87）。」を足し、Decision date を 2026-10-19 にします。ほかの欄（Must preserve を含む）は rev 9 のままです。この受領文を承認しますか？（ページの「実装と受領文の確認」）」
    the selected option and its description:
    「承認 (推奨)」
    「このとおり Contract rev 67 に記録し、この返答を Reply 24 として引用して SHA-256 を計算する。そのあと review を経て PR を開く（merge は Owner）。末尾の「…」もこのまま。」
    and the answer:
    「承認 (推奨)」
    The record's scenario revision `10` receipt and the change receipt below
    are the approved texts.
    So the record's scenario revision `10` receipt appends
    「Circular の領域ラベルが円に収まらないときは、収まる分を曲線で描いたうえで警告する（LP-D87）。」 to the Rationale of scenario revision `9`; its other
    fields are those of scenario revision `9`, with the decision date of the
    last reply.
  - The receipt of this change (SHA-256, UTF-8, excluding the final newline:
    `6192643855d405eb55379ad7f92cae6f12745c54938f4ae9355ccc495219c81e`):

    ```text
    PRODUCT_DECISION
    Concern: product.circular-outer-label-fit
    Scenario revision: 10
    Choice: RGN-CUT-WARN
    Rationale: not stated
    Must preserve: not stated
    May retire: not stated
    Accepted residual risk: not stated
    Owner: Satoshi Kawato
    Decision date: 2026-10-19
    ```

  - Scope: the record's receipt for scenario revision `10` is recorded below;
    the scenario revision `7`, `8`, and `9` receipts keep their text and their
    SHA-256. `PD-OI-088` applies unchanged. No other record is amended or
    retired.
  - This change is a static Product Contract co-change: the runtime, tests, and
    documentation that implement it are in the same pull request, and the
    Review is `REQUIRED`.
- Revision 68 changes: `PD-OI-097` is replaced for scenario revision `11`
  (concern `product.circular-outer-label-fit`), and `PD-OI-098` for scenario
  revision `2` (concern `diagram-generation.diagram-font`; the outcome
  `A+B / CHOSEN-FONT-MEASURED-AND-EMBEDDED` is kept). They come from the
  Product Decision Owner `satoshikawato`'s replies of `2026-10-19` on the
  style of feature labels: what keeps the leader lines that pass behind a
  label off its letters, a Linear underline joined to the leader line, the
  width of leader lines, and a minimum gap between labels in one Linear row
  (OV-1355). Each quotation below is the part of a reply that carries the
  answer; Replies 1 to 24 are those of the Revision 50, 53, 66, and 67
  entries.
  - Reply 25 (scenario revision `11`, `2026-10-19`), to the question
    「Q9: フィーチャーラベルを、後ろを通るリーダー線からどう守りますか？（ページの図1・2・6・8）」
    the selected option and its description:
    「ハロー 0.08-0.12em (推奨)」
    「文字の形に沿った白い縁取り。リーダー線は文字の間で元の濃さのまま残る。領域ラベルの背景 0.7 は今のまま」
    and the answer:
    「ハロー 0.08-0.12em (推奨)」
    (`LBL-HALO`)
  - Reply 26 (scenario revision `11`, `2026-10-19`), to the question
    「Q10: Linear のデフォルトはどれにしますか？（図1・3・4）」
    the selected option and its description:
    「下線+ハロー (推奨)」
    「LP-D93 の「下線+半透明背景」と同じ目的を、リーダー線を隠さないハローで果たす」
    and the answer:
    「下線+ハロー (推奨)」
    (`LIN-UNDER`)
  - Reply 27 (scenario revision `11`, `2026-10-19`), to the question
    「Q11: Circular のデフォルトはどれにしますか？（図6・7）」
    the selected option and its description:
    「ハロー、下線なし (推奨)」
    「Circular でリーダー線が文字を通る所は少なく（lambda で3か所）、違いは小さい」
    and the answer:
    「ハロー、下線なし (推奨)」
    (`CIR-NO-UNDER`)
  - Reply 28 (scenario revision `11`, `2026-10-19`), to the question
    「Q12: 上から来るリーダー線（ラベルがトラックの下にあるときなど）はどう終えますか？（図5）」
    the selected option and its description:
    「0.20em 上で止める (推奨)」
    「文字の 0.20 em 上で止める。線の交差も折れ点も増えない（LP-D90 と両立）」
    and the answer:
    「0.20em 上で止める (推奨)」
    (`GAP-STOP`)
  - Reply 29 (scenario revision `11`, `2026-10-19`), to the question
    「Q13: リーダー線（と下線）の太さをどうしますか？（図4・7・9・10）」
    the selected option and its description:
    「両モード 0.07倍・≧1px (推奨)」
    「文字サイズの 0.07 倍、ただし 1.0 px 以上。Linear と長い配列の Circular は 0.5→1.0 px 以上に太く、短い配列の Circular（14 pt）は 1.5→1.0 px に細くなる。下線と文字の間は線の縁から測る」
    and the answer:
    「両モード 0.07倍・≧1px (推奨)」
    (`LEADER-007`)
  - Reply 30 (scenario revision `11`, `2026-10-19`), to the question
    「Q14: 選べるようにするオプションの形は？」
    the selected option and its description:
    「2つに分ける (推奨)」
    「label_line（none / under）と label_protection（none / halo / box）。下線と守りを独立に選べ、「下線+背景」も選べる」
    and the answer:
    「2つに分ける (推奨)」
    (`LBL-OPTS-2`)
  - Reply 31 (scenario revision `11`, `2026-10-19`), to the question
    「Q15: この変更（Q9〜Q14）をいつ入れますか？」
    the selected option and its description:
    「0.14.0 に入れる」
    「タグの前に入れる（P-03 の凍結の例外になる）」
    and the answer:
    「0.14.0 に入れる」
    (`LBL-0140`)
  - Reply 32 (scenario revision `11`, `2026-10-19`), to the question
    「Q-A: フィーチャーラベルのハローの幅は？（ページ https://claude.ai/artifact/NN3ZELDPbp3u6pFqP8TsTz Version 6 の図 Q-A 1・2）」
    the selected option and its description:
    「0.08 em (推奨)」
    「読む大きさでは差が見えず、リーダー線を一番多く残す（ラベル内で見える割合 89 %、dev 83 %）」
    and the answer:
    「0.08 em (推奨)」
    (`HALO-008`)
  - Reply 33 (scenario revision `11`, `2026-10-19`), to the question
    「Q-B: 名前と置き場所はこれでよいですか？ CLI --label_line {none,under} / --label_protection {none,halo,box}（両モード）、config [labels.linear]/[labels.circular] の line / protection、Python は config_overrides のみ、GUI は Labels の Label Rendering の後に Label Line（None/Underline）と Label Protection（None/Halo/Box）。ハローの幅と箱の不透明度は固定。」
    the selected option and its description:
    「このまま (推奨)」
    「上の名前と置き場所。GUI の画面は PR でスクリーンショットをお見せします」
    and the answer:
    「このまま (推奨)」
    (`OPT-NAMES`)
  - Reply 34 (scenario revision `11`, `2026-10-19`), to the question
    「Q-D: この変更より前に保存した Session の見た目は？（図 Q-D 1・2）」
    the selected option and its description:
    「A 次の Generate で新しく (推奨)」
    「保存された図は Generate までそのまま。Python で labels.stroke_width を指定するとエラー（CHANGELOG に書く）」
    and the answer:
    「A 次の Generate で新しく (推奨)」
    (`SESS-NEXT-GEN`)
  - Reply 35 (`2026-10-19`), to the question
    「Q-F: 全体の図（図 Q-D 2）を見て、Linear のリーダー線の太さは Q13 のままでよいですか？ 24 pt で 0.5 → 1.68 px になり、下線とあわせて灰色の格子が目立ちます。」
    the selected option and its description:
    「Q13 のまま」
    「文字サイズの 0.07 倍、1.0 px 以上。24 pt で 1.68 px、50 % に縮めてもコントラスト 2.49:1」
    and the answer:
    「Q13 のまま」
    Reply 35 confirms the width of Reply 29 and adds no choice token.
  - Reply 36 (scenario revision `11`, `2026-10-19`), to the question
    「S1 Web PDF のハロー：Web の PDF 書き出し（svg2pdf）は、ハローを黒い文字の重複として描き、ラベルを 2 回取り出してしまいます（ページの図 S1）。どうしますか？」
    the selected option and its description:
    「B ハローを描かない (推奨)」
    「Web PDF だけハローを外す。SVG・PNG・CLI の PDF にはハローがある。Web PDF の見た目は今と同じで、ラベルは 1 回ずつ取り出せる。」
    and the answer:
    「B ハローを描かない (推奨)」
    (`WEBPDF-NO-HALO`)
  - Reply 37 (scenario revision `11`, `2026-10-19`), to the question
    「S2 Linear の段と下線：計画 (a) では全ラベルが動いて図が高くなり、BGC で段が変わり、HmmtDNA 20 pt で警告なしの重なりが出ました。(b') なら実データは dev と同じ位置・大きさです（ページの図 S2 1〜6）。どうしますか？」
    the selected option and its description:
    「(b') 段に数えない (推奨)」
    「下線は段の高さに数えず、領域名だけが下線をよける。実データの図は dev と同じ。領域名が下線に当たるときだけ、ラベルが新しい段に移るか警告。背景で下線を隠さない。」
    and the answer:
    「(b') 段に数えない (推奨)」
    (`UNDER-NOT-IN-ROW`)
  - Reply 38 (scenario revision `11`, `2026-10-19`), to the question
    「Q-G 下線の白い縁：上の段へ伸びるリーダー線が、ほかのラベルの下線を横切るとつなぎ目に見えます（lambda で 293 か所）。下線に白い縁を付けると、横切る線は下線の所で途切れます（新しいページの図 G 1〜3）。どうしますか？」
    the selected option and its description:
    「付ける (推奨)」
    「ほかのラベルのリーダー線は下線の所で途切れ、後ろを通るように見える。自分のリーダー線はつながったまま。文字や形は隠れない。PDF +0.5 %。」
    and the answer:
    「付ける (推奨)」
    (`UNDER-CASING`)
  - Reply 39 (scenario revision `11`, `2026-10-19`), to the question
    「Q-H 下線の端：同じ段で隣り合うラベルの下線が、密な所では 1 本につながって見えます（tobacco IR 5 pt で 0.43 px）。下線の両端を 0.12 em ずつ短くしますか（図 H 1〜3）？」
    the selected option and its description:
    「両端を 0.12 em 短く (推奨)」
    「隣どうしの下線が離れる（0.43 → 1.63 px、24 pt で 1.07 → 6.83 px）。どの下線も文字より少し短くなる（24 pt で片側 2.9 px）。ラベルの位置は変わらない。」
    and the answer:
    「両端を 0.12 em 短く (推奨)」
    (`UNDER-TRIM`)
  - Reply 40 (scenario revision `11`, `2026-10-19`), to the question
    「Q-I 同じ段の最小のすき間：Linear では、同じ段の隣どうしのラベルがすき間なしで並び、1 語に見えることがあります（HmmtDNA 24 pt の "tRNA-SertRNA-Lys"。dev にもある不具合 OV-1355）。すき間を入れると、詰まった段のラベルが別の段に移ります（新しいページの図 I 1〜3）。どうしますか？」
    the selected option and its description:
    「0.25 em、LS の PR で (推奨)」
    「ほぼ単語 1 つ分の空白。近すぎる組はどの図でもなくなる。lambda 21、tobacco IR 22、HmmtDNA 4 ラベルが別の段へ。図の高さは変わらない。出力の作り直しは 1 回で済む。」
    and the answer:
    「0.25 em、LS の PR で (推奨)」
    (`ROW-GAP-025`)
  - Reply 41 (`2026-10-19`), the message on figure E1 of the page that showed
    the built figures:
    「実装した図の確認: やっぱ塗りつぶし前の方がいいかも。Lambda,ラベルの背景のリーダー線目立ちすぎ。下線＋塗りつぶしがいいかな。Thoughts?」
  - Reply 42 (scenario revision `11`, `2026-10-19`), to the question
    「Q-K1：Linear のフィーチャーラベルの既定をどれにしますか？白い背景は文字の外 0.18 em までで、ラベルの位置・段・図の大きさはどの案でも変わりません。白い背景なら PDF はほとんど大きくならず（lambda 24 pt でハロー 355 KB → 27 KB）、Web の PDF も SVG と同じ見た目になります。背景の濃さは Linear と Circular の Box で共通の 1 つの値です。（ページの図 K 1〜K 5）」
    the selected option and its description:
    「下線＋白い背景 0.8」
    「文字の上でも線が 2 割の濃さで続いて見える。」
    and the answer:
    「下線＋白い背景 0.8」
    (`LIN-BOX-080`)
  - Reply 43 (`2026-10-19`), to the question
    「Q-K2：Circular の既定はどうしますか？白い背景にすると、タバコ葉緑体のように色の付いた帯の上で、ラベルごとに白い四角が見えます（ページの図 K 6）。」
    the selected option and its description:
    「ハローのまま (推奨)」
    「Circular でリーダー線が文字を通る所は少ない（lambda で 3 か所）。帯の上に白い四角が出ない。」
    and the answer:
    「ハローのまま (推奨)」
  - Reply 42 did not select the recommended option. For Linear it replaces
    the halo of Replies 25 and 26 with a white box behind the letters at
    opacity `0.8`, the one fixed Box opacity of Reply 33 for both modes. Reply 43
    keeps the halo of Reply 27 for Circular and adds no choice token.
  - Approvals the Owner gave on `2026-10-19` of an earlier wording of both
    receipts, with a halo as the Linear default, were made void by Reply 42
    before they reached `dev`; they are not part of this record.
  - Replies 25 to 43 state no rationale, preserved or retired effect, or
    accepted risk as a receipt field; the descriptions of their options say
    how each choice is realized. The approved receipts of Replies 44 to 46
    word those fields, and the change receipt below adds none.
  - Reply 44 (`2026-10-19`), approving the `PD-OI-097` receipt texts below,
    to the question
    「Q-E1（書き直し）：PD-OI-097 rev 11 の受領文です。前の文面から、Linear の守りを白い背景（不透明度 0.8）、Circular をハローにし、Choice に LIN-BOX-080 を足し、Must preserve に「領域ラベルは下線か白い背景にかかるとき別の段に移すか警告する」「Linear の白い背景はほかのラベルの文字と下線を覆わない」を入れ、Accepted residual risk の PDF/EPS の数字を測り直した値にしました。Decision date は 2026-10-19 です。足す文と全文はページの一番上の節にあります。この受領文を承認しますか？」
    the selected option and its description:
    「承認 (推奨)」
    「このとおり Contract rev 68 に記録し、この返答を Reply として引用して SHA-256 を計算する。そのあと出力の作り直しとレビューを経て PR を開く（merge は Owner）。」
    and the answer:
    「承認 (推奨)」
    The change receipt below is the approved text. Reply 46 corrected one
    Must preserve sentence of the record's receipt that Reply 44 approved,
    before that wording reached `dev`.
    So the record's scenario revision `11` receipt appends the choices of
    Replies 25 to 34, 36 to 40, and 42 to the Choice; the sentences from
    「フィーチャーラベルは、後ろを通るリーダー線から守る。」 to the end to the Rationale;
    the sentences from 「0.13.0 のフィーチャーラベルの見た目（細いリーダー線、守りと下線なし）を既定で出すこと。」
    to the end to the May retire; and the sentences from
    「Circular の PDF と EPS はハローの分だけ大きくなる（タバコ葉緑体の PDF は 56 KB から 235 KB、EPS は 1.32 MB から 1.59 MB）。」
    to the end to the Accepted residual risk of scenario revision `10`. In
    its Must preserve, 「重なりも交差もない図は 0.13.0 と同じに描く。」 reads the
    sentences from 「重なりも交差もない図は 0.13.0 と同じ位置にラベルを描く。」 to
    「Linear の白い背景は、ほかのラベルの文字と下線を覆わない。」. Its other fields
    are those of scenario revision `10`, with the decision date of the last
    reply.
  - The receipt of this change (SHA-256, UTF-8, excluding the final newline:
    `09af42ad585ac8add25662150d394b6e45b4ea647fa325a261d773c650fcf5a5`):

    ```text
    PRODUCT_DECISION
    Concern: product.circular-outer-label-fit
    Scenario revision: 11
    Choice: LBL-HALO ; LIN-UNDER ; CIR-NO-UNDER ; GAP-STOP ; LEADER-007 ; LBL-OPTS-2 ; LBL-0140 ; HALO-008 ; OPT-NAMES ; SESS-NEXT-GEN ; WEBPDF-NO-HALO ; UNDER-NOT-IN-ROW ; UNDER-CASING ; UNDER-TRIM ; ROW-GAP-025 ; LIN-BOX-080
    Rationale: not stated
    Must preserve: not stated
    May retire: not stated
    Accepted residual risk: not stated
    Owner: Satoshi Kawato
    Decision date: 2026-10-19
    ```

  - Reply 45 (`2026-10-19`), approving the `PD-OI-098` receipt text recorded
    under that record below, to the question
    「Q-E2（書き直し）：PD-OI-098 rev 2 の受領文です。〔変えないもの〕の「font を選ばない図の文字の配置と出力」の後に「（ただし、フィーチャーラベルのハローと白い背景、下線とその白い縁、リーダー線の太さ、Linear の同じ段のすき間（LP-D98、LP-D100、LP-D102〜LP-D104、LP-D106）と、Web の PDF で中央揃えと右揃えの文字の位置を直したこと（OV-1303）は除く）」を足し、Scenario revision 2、Decision date 2026-10-19 にします。ほかの欄は rev 1 のままです。この受領文を承認しますか？」
    the selected option and its description:
    「承認 (推奨)」
    「このとおり Contract rev 68 に記録する。」
    and the answer:
    「承認 (推奨)」
  - Receipt change (`PD-OI-098`): in 〔変えないもの〕 of the Must preserve,
    「font を選ばない図の文字の配置と出力、」
    reads
    「font を選ばない図の文字の配置と出力（ただし、フィーチャーラベルのハローと白い背景、下線とその白い縁、リーダー線の太さ、Linear の同じ段のすき間（LP-D98、LP-D100、LP-D102〜LP-D104、LP-D106）と、Web の PDF で中央揃えと右揃えの文字の位置を直したこと（OV-1303）は除く）、」.
    Every other field of the scenario revision `1` receipt is unchanged,
    with Scenario revision `2` and the decision date of Reply 45; the
    scenario revision `1` receipt stays in the record.
  - Reply 45 states no rationale, preserved or retired effect, or accepted
    risk for this change, and the receipt adds none.
  - New receipt SHA-256: `PD-OI-098` scenario revision `2`
    `a305e5d6eefb69914e20e1397cca98298c3050315dff6d488638c2dff5160cc8`.
  - Reply 46 (`2026-10-19`), approving the corrected `PD-OI-097` receipt text
    recorded under that record below, to the question
    「Q-E1（訂正）：PR 前のレビューで、承認いただいた PD-OI-097 rev 11 の Must preserve の 1 文が動きと合っていませんでした（「領域ラベルを別の段に移す」と書いたが、実際は文字にかかるときと同じくフィーチャーラベルの方を別の段へ移し、解消しなければ描いて警告する）。動きは変えず、その 1 文を「Linear の領域ラベルは、フィーチャーラベルの下線と白い背景を文字と同じによける。よけられないときは、文字にかかるときと同じく、そのフィーチャーラベルを同じ側の別の段へ移す（RGN-YIELD-ROW）。移しても解消しないときは、そのまま描いて警告する。領域ラベルの背景は、リーダー線の上には出るが、下線と白い背景は覆わない（RGN-FRONT）。」に直します。ほかの文と PD-OI-098 rev 2 はそのままです（ページの一番上の節）。この訂正した受領文を承認しますか？」
    the selected option and its description:
    「承認 (推奨)」
    「訂正した文面を Contract rev 68 に記録し、この返答を Reply として引用して SHA-256 を計算し直す。そのあと PR を開く（merge は Owner）。」
    and the answer:
    「承認 (推奨)」
    The record's scenario revision `11` receipt is the approved text. In its
    Must preserve, the sentence that Reply 44 approved as
    「ただし Linear で同じ段の隣のラベルとの間が 0.25 em に満たないときはそのラベルを別の段に移し、Linear の領域ラベルがフィーチャーラベルの下線か白い背景にかかるときは領域ラベルを別の段に移すか警告する。」
    reads
    「ただし Linear で同じ段の隣のラベルとの間が 0.25 em に満たないときは、そのラベルを別の段に移す。Linear の領域ラベルは、フィーチャーラベルの下線と白い背景を文字と同じによける。よけられないときは、文字にかかるときと同じく、そのフィーチャーラベルを同じ側の別の段へ移す（RGN-YIELD-ROW）。移しても解消しないときは、そのまま描いて警告する。領域ラベルの背景は、リーダー線の上には出るが、下線と白い背景は覆わない（RGN-FRONT）。」.
    The behavior is unchanged, and every other sentence is as Reply 44
    approved it. The wording that Reply 44 approved was corrected by Reply 46
    before it reached `dev`; it is not part of this record.
  - Reply 46 states no rationale, preserved or retired effect, or accepted
    risk beyond the corrected sentence, and the change receipt above is
    unchanged.
  - Scope: the record receipts for `PD-OI-097` scenario revision `11` and
    `PD-OI-098` scenario revision `2` are recorded below; the `PD-OI-097`
    scenario revision `7`, `8`, `9`, and `10` receipts and the `PD-OI-098`
    scenario revision `1` receipt keep their text and their SHA-256.
    `PD-OI-088` applies unchanged. No other record is amended or retired.
  - This change is a static Product Contract co-change: the runtime, tests, and
    documentation that implement it are in the same pull request, and the
    Review is `REQUIRED`.
- Records remaining `EVIDENCE_REQUIRED`: none
- Excluded records: none

This contract owns the user-observable outcomes recorded below. It does not
select source files, classes, module paths, canonical call edges, delivery
order, cache implementation, or test design. Code, tests, fixtures,
screenshots, and historical behavior are evidence, not Product authority.

The mapped concern `product.canonical-render-request-boundary`, scenario
revision `1`, remains governed by its existing selected option
`canonical-typed-request-boundary`. This contract references that authority
where applicable and does not replace or duplicate it. No active `BD-###`
decision governs the concerns recorded here at this revision.

## Interpretation and lifecycle

Active records use only these statuses:

- `ACCEPTED`: the complete outcome is normative.
- `EVIDENCE_REQUIRED`: evidence-only work may proceed, but implementation must
  not select the pending outcome.
- `DEFERRED`: the named capability is outside the current delivery scope.
- `UNSUPPORTED`: the named input or journey must be rejected explicitly.

Authority precedes dependent runtime implementation. A runtime change cites
authority already present on its base and does not change the decision it
implements. Evidence precedes a decision when a record is
`EVIDENCE_REQUIRED`; evidence does not select its own outcome.

To correct an active outcome, increment the scenario revision, identify the
prior decision and revision in `Supersedes`, and record the complete
replacement outcome with its explicit receipt. The replacement either merges
as an authority-only change before dependent runtime, or merges together with
its implementation through the reviewed static Product Contract co-change
route in [`WEB_CHANGE_POLICY.md`](./WEB_CHANGE_POLICY.md#static-product-contract-co-change).
Where a record says that dependent runtime requires its authority merged into
its base, a reviewed co-change that merges the record and its implementation
together satisfies that requirement. That route changes only this procedure.
It does not change any recorded outcome, preservation condition, or
compatibility commitment. Git history retains the former text; the active
contract does not accumulate superseded records.

## Cross-surface clauses

### OIPC-C01: Omitted and explicit values

- Omission uses the public typed default.
- An explicitly supplied default is execution-equivalent to omission.
- An explicit valid non-default changes the documented execution or
  presentation behavior.
- Invalid explicit values are rejected, not silently coerced.

### OIPC-C02: Requested and effective values

When automatic or context-dependent resolution exists, requested intent and
effective execution are both retained. Requested intent is not overwritten by
the effective subtype.

### OIPC-C03: Consume or reject

Every accepted public value reaches its real consumer or is rejected before
execution. A surface must not accept and silently ignore a value.

### OIPC-C04: Request, execution, cache, and artifact agreement

For each applicable field, canonical request intent, resolved execution
values, actual helper invocation, correct stage-specific cache identity,
Session data, and artifact metadata agree. Requested/effective differences are
explicit. The mapped `product.canonical-render-request-boundary` concern
continues to own canonical Web request continuity.

### OIPC-C05: Preservation of valid intent

A valid value is not deleted because a surface lacks an editor. The generic
surface disposition is `EDITABLE`, `READ_ONLY`, `PASS_THROUGH`, or
`UNSUPPORTED`. Imported comparison reconstruction uses the more specific
states `EDITABLE`, `PRESERVED_READ_ONLY`, and `DECISION_REQUIRED` defined in
`PD-OI-008`.

### OIPC-C06: Explicit replacement and clearing

Replacement and removal are explicit actions. Empty controls, missing
properties, inactive modes, failed reconstruction, and failed generation do
not imply deletion.

### OIPC-C07: Failure isolation

Failed, canceled, superseded, or stale generation does not replace the last
successful Result or committed request.

### OIPC-C08: Evidence is not authority

A test or historical implementation that contradicts an active decision is
corrected. Passing evidence does not make incorrect behavior normative.

## Product Decision records

Record requirement from 2026-10-05 (`gbdraw/web/CLAUDE.md` R13): a decision
in which one editor domain follows another (a hidden feature hides its label,
a label action shows its feature, a slot input routes through placement)
names, in its Choice and in its receipt's `Must preserve` or `Normative
outcome`, the **reaction owner** (the module that owns the reaction), the
**channel** (`port`, `root-projection`, or `event`), and the **projection
function**. Records accepted before this date keep their text; their reaction
owners are recorded by the runtime pull requests of the owner-coupling plan
(`docs/internal/WEB_OWNER_COUPLING_PREVENTION_IMPLEMENTATION_PLAN_2026-10-05.md`,
Phase E): F-3 and `PD-OI-066` label-follows-feature → `app/feature-editor/visibility-actions.js`,
`port` (`applyFeatureVisibilityToLabels`), `projectFeatureVisibility`;
Owner Q2 show-feature-and-label → `app/feature-editor/visibility-actions.js`,
`port` (`setFeatureVisibility`), `projectFeatureVisibility`; Owner Q3
slot inputs → `app/feature-editor/placement-actions.js`, `port`
(`changeTrackLayout`), the transition itself; E7 (#822) kept #805's direction,
so placement remains the writer of the slot inputs.


### PD-OI-001: LOSATP raw-search limit fresh defaults

- Concern key: `diagram-generation.losatp-candidate-limit-default`
- Scenario revision: `2`
- Supersedes: `PD-OI-001`, scenario revision `1`.
- Status: `ACCEPTED`
- Normative outcome: Web exposes the raw LOSATP `max_target_seqs` limit as
  **Max target seqs**. Fresh and reset Collinear starts at `5`; Similarity
  groups starts unbounded (`None`). Pairwise, CLI, and Python omission defaults
  retain their existing meanings. A blank Web value explicitly means unbounded,
  including in Collinear. No hidden cap substitutes for an unbounded request.
- Rationale: The maintainer requested exposure of the actual raw-search limit
  and distinct Collinear/Similarity defaults, with regression protection.
- Must preserve: Explicit finite and unbounded values; truthful requested and
  effective metadata; raw-cache identity; cancellation and errors; saved values.
- May retire: The unbounded fresh Web Collinear default. Unbounded search itself
  remains available.
- Accepted residual risk: The existing unbounded-work cost remains. This
  amendment adds no performance guarantee or additional risk waiver.
- Acceptance contracts: `OIC-001`, `OIC-005`, `OIC-013`, `OIC-017`.
- Decision source: Explicit maintainer instruction to expose Max target seqs,
  default Collinear to `5`, and leave Similarity groups unbounded.
- Owner and decision date: `satoshikawato`, `2026-09-14`.

### PD-OI-002: LOSATP mode-specific limit retention and GUI placement

- Concern key: `diagram-generation.losatp-candidate-limit-scope`
- Scenario revision: `2`
- Supersedes: `PD-OI-002`, scenario revision `1`.
- Status: `ACCEPTED`
- Normative outcome: **Max target seqs** is directly editable in LOSATP
  Settings, alongside **Member hits per protein** when member selection applies.
  Each mode remembers its own raw and member limits. The first visit to
  Collinear uses `5`/`5`; the first visit to Similarity groups uses unbounded/
  unbounded. After editing either mode, switching away and back restores its
  edited values. This works repeatedly in both directions and preserves blanks.
  Switching modes neither resets the returning mode to defaults nor copies the
  departing mode's limits into it. Save and fresh Load preserve active and
  inactive limits; Reset Settings clears the drafts to their declared defaults.
- Rationale: The maintainer clarified the required sequence as edit Similarity,
  edit Collinear, return to Similarity's values, then return to Collinear's values.
- Must preserve: Independently edited limits; existing Session values;
  Pairwise display-limit independence; discoverability and keyboard operation.
- May retire: A single raw-limit value shared across presentations and
  Advanced-only placement. There must not be duplicate controls for one active
  raw-search setting.
- Accepted residual risk: No additional risk waiver was supplied. Search
  changes caused by selecting a different saved raw limit remain visible;
  matching raw evidence remains reusable.
- Acceptance contracts: `OIC-001`, `OIC-002`, `OIC-006`, `OIC-017`.
- Decision source: Explicit maintainer correction to remember each mode's
  settings, superseding the earlier same-session request to reset on every switch.
- Owner and decision date: `satoshikawato`, `2026-09-14`.

### PD-OI-003: Pairwise display max hits

- Concern key: `diagram-generation.pairwise-display-max-hits`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Normative outcome: Pairwise display max hits has a fresh default of `5` and
  is applied after threshold filtering to Pairwise result selection only.
- Rationale: Presentation density must not change raw search evidence.
- Must preserve: Raw-search cache reuse when only this value changes; explicit
  Session values.
- May retire: Aliasing Pairwise display max hits to Candidate limit or member
  hits.
- Accepted residual risk: Additional qualified hits remain absent from the
  Pairwise view while retained in raw evidence.
- Acceptance contracts: `OIC-002`, `OIC-005`.
- Owner and decision date: `satoshikawato`, `2026-08-28`.

### PD-OI-004: Similarity/Collinear member hits per protein

- Concern key: `diagram-generation.member-hits-per-protein`
- Scenario revision: `3`
- Supersedes: `PD-OI-004`, scenario revision `2`.
- Status: `ACCEPTED`
- Normative outcome: Fresh and reset Web **Member hits per protein** starts at
  `5` in Collinear and unbounded (`None`) in Similarity groups. CLI/Python
  omission semantics remain unchanged. A blank Web control always means
  unbounded, including in Collinear; an explicit positive integer limits the
  distinct directional subject candidates retained after result filtering.
  Member selection remains independent of raw Max target seqs and Pairwise
  display matches. Collinear consumes it with inference both OFF and ON.
  Per-mode retention follows `PD-OI-002`. Request, helper execution, provenance,
  and Session replay distinguish unbounded and finite choices.
- Rationale: The maintainer requested the same mode defaults and retention for
  both exposed limits, while preserving their distinct scientific roles.
- Must preserve: Raw-search reuse when only member hits changes; correct
  derived invalidation; explicit saved values; threshold filtering.
- May retire: The unbounded fresh Web Collinear member default. Blank-to-five
  coercion and aliasing this field to the raw or Pairwise limits remain prohibited.
- Accepted residual risk: Member selection can change blocks or groups and
  remains visible in provenance. The existing unbounded computation/memory
  allowance is unchanged; no additional waiver is recorded.
- Acceptance contracts: `OIC-002`, `OIC-004`, `OIC-005`, `OIC-006`, `OIC-017`.
- Decision source: Explicit maintainer instruction that Member hits per protein
  receive the same mode defaults and remembered-value behavior as Max target seqs.
- Owner and decision date: `satoshikawato`, `2026-09-14`.

### PD-OI-005: Supported Collinear value domains

- Concern key: `diagram-generation.collinear-value-domains`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Normative outcome: Web execution consumes every currently public Collinear
  enum value, including the current equivalents of unit mode `auto`, `cds`,
  `locus`; anchor mode `all`, `one_to_one`, `rbh`; and merge orientation
  `strand`, `order`, `either`. Exact names are verified against the current
  typed API before implementation.
- Rationale: A valid public value must reach execution or be rejected
  explicitly.
- Must preserve: Typed Python validation; requested `auto`; separately
  reported effective resolution.
- May retire: Browser or Worker branches that coerce one valid enum to another.
- Accepted residual risk: Combinations are covered primarily by unit/contract
  tests rather than browser cases.
- Acceptance contracts: `OIC-003`, `OIC-005`.
- Owner and decision date: `satoshikawato`, `2026-08-28`.

### PD-OI-006: Collinear search-scope fresh default

- Concern key: `diagram-generation.collinear-search-scope-default`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Normative outcome: The fresh search-scope default is `adjacent`; `all`
  remains an explicit supported value.
- Rationale: Adjacent comparison is the primary Web journey while explicit
  all-pairs analysis remains available.
- Must preserve: Explicit `all` across CLI, Python API, Web, and imported
  Session.
- May retire: Conflicting fresh defaults across surfaces.
- Accepted residual risk: Fresh output can differ from a historical surface
  that used implicit `all`; released Sessions must be preserved or migrated
  explicitly, not reinterpreted.
- Acceptance contracts: `OIC-003`, `OIC-004`, `OIC-005`.
- Owner and decision date: `satoshikawato`, `2026-08-28`.

### PD-OI-007: Collinear merge-conflict fresh default

- Concern key: `diagram-generation.collinear-max-conflicts-default`
- Scenario revision: `2`
- Status: `ACCEPTED`
- Normative outcome: The fresh `max_conflicts` default is one (`1`). Explicit
  zero (`0`) remains supported. Omission and explicit one are
  execution-equivalent.
- Rationale: Merged deterministic evidence established the consumer behavior
  of both values, and the Product Decision Owner selected one after reviewing
  that evidence.
- Must preserve: Explicit zero and one; retained singleton anchors; agreement
  among request intent, execution, round trip, and provenance; reproducible
  evidence for the merge-threshold effect.
- May retire: Conflicting fresh omission defaults across surfaces.
- Accepted residual risk: At one, compatible clusters may merge across one
  retained interior singleton where zero keeps them separate. The singleton
  remains in the result, and the selected behavior remains visible in
  provenance.
- Acceptance contracts: `OIC-004`, `OIC-005`.
- Evidence: [`COLLINEAR_MAX_CONFLICTS_COMPARISON_EVIDENCE_2026-08-28.md`](./COLLINEAR_MAX_CONFLICTS_COMPARISON_EVIDENCE_2026-08-28.md),
  merged by PR `#422` at `878c62ba17c61c45cd0adbd05cbd9fb36306db9d`.
- Supersedes: `PD-OI-007`, scenario revision `1`.
- Owner and decision date: `satoshikawato`, `2026-08-28`.

### PD-OI-008: Imported comparison reconstruction and resolution

- Concern key: `diagram-generation.imported-comparison-state`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Normative outcome: An imported committed comparison is `EDITABLE`,
  `PRESERVED_READ_ONLY`, or `DECISION_REQUIRED`. Exact reconstruction becomes
  editable without marking a user edit. A valid executable but non-projectable
  comparison remains usable through explicit inheritance. An ambiguous,
  incomplete, or non-executable comparison blocks Generate until explicit
  replacement or clearing. No comparison is silently converted to no
  comparison.
- Rationale: GUI limitations must not destroy valid work or make regeneration
  untrustworthy.
- Must preserve: Last successful Result; committed request; Save/export; exact
  executable data; explicit continuation; failure isolation.
- May retire: Silent empty-comparison fallback, mode-change resets, and
  reconstruction that invents defaults.
- Accepted residual risk: Some imported Sessions require read-only disclosure
  or explicit resolution before regeneration.
- Acceptance contracts: `OIC-007`, `OIC-013`.
- Owner and decision date: `satoshikawato`, `2026-08-28`.

### PD-OI-009: GUI-unmanaged configuration overrides

- Concern key: `diagram-generation.gui-unmanaged-config-overrides`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Normative outcome: Schema-known, valid, safe, active-mode-compatible
  configuration leaves without a GUI editor are preserved losslessly and
  disclosed. Unknown paths, non-leaf paths, unsafe keys, invalid literals, and
  active-mode-incompatible values are rejected explicitly.
- Rationale: Lack of an editor is not authorization to discard a valid typed
  value; generic unknown pass-through is unsafe.
- Must preserve: Managed siblings; valid imported leaves; typed validation;
  explicit reset.
- May retire: Silent dropping of valid GUI-unmanaged values and permissive
  unknown-path pass-through.
- Accepted residual risk: Users can see read-only settings they cannot edit in
  the current GUI.
- Acceptance contracts: `OIC-008`.
- Owner and decision date: `satoshikawato`, `2026-08-28`.

### PD-OI-010: Circular record selection, crop, and reverse display

- Concern key: `diagram-generation.circular-record-transforms`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Normative outcome: Circular Web workflows support selection of one record,
  valid region crop, and reverse-complement display. Arbitrary multi-record
  subset editing is outside this program.
- Rationale: These operations form one coherent single-record preparation
  journey.
- Must preserve: Source record order; stable record identity; deterministic
  selection; coordinate validation; Session round trip; no double reverse.
- May retire: Hidden accepted fields that never reach rendering and duplicate
  display-direction state.
- Accepted residual risk: Arbitrary multi-record subset editing remains
  unavailable.
- Acceptance contracts: `OIC-009`.
- Owner and decision date: `satoshikawato`, `2026-08-28`.

### PD-OI-011: Linear scale and ruler-label font interaction

- Concern key: `diagram-generation.linear-scale-ruler-fonts`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Normative outcome: Scale font size and ruler-label font size are linked by
  default and can be unlinked in Advanced settings. Imported Sessions with
  unequal values open unlinked and retain both. Explicit relink states the
  effect and copies the current scale font size to ruler-label font size once.
- Rationale: The public values are independent while linked fresh state keeps
  the common UI simple.
- Must preserve: Both explicit values; Session round trip; deterministic
  relink behavior.
- May retire: One-field aliasing that overwrites the other value.
- Accepted residual risk: One small linked/unlinked UI state transition is
  added.
- Acceptance contracts: `OIC-011`.
- Owner and decision date: `satoshikawato`, `2026-08-28`.

### PD-OI-012: Circular record label and subtitle

- Concern key: `diagram-generation.circular-label-subtitle`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Normative outcome: Explicit record label and subtitle values affect Circular
  rendering and are editable in the applicable Web journey. Explicit values
  override corresponding inferred title lines; empty values preserve inferred
  output.
- Rationale: Accepted public presentation fields must work or be rejected.
- Must preserve: Existing output when both values are empty; explicit values
  through Session round trip.
- May retire: Parsing and serialization paths that accept fields without
  consuming them.
- Accepted residual risk: Explicit text can change layout and requires targeted
  geometry review.
- Acceptance contracts: `OIC-010`.
- Owner and decision date: `satoshikawato`, `2026-08-28`.

### PD-OI-013: Web compatibility meaning of `grid_column`

- Concern key: `diagram-generation.grid-column-compatibility`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Normative outcome: The Web normalizes `grid_column` to row-internal render
  ordering. Render-equivalent placement is guaranteed; exact numeric identity
  is not a Web round-trip promise.
- Rationale: Preserving a number with no distinct supported visual effect
  would add compatibility complexity without Product value.
- Must preserve: Relative order and rendered placement.
- May retire: Exact-number assertions without a supported visual difference.
- Accepted residual risk: A Web-exported Session can use different numeric
  columns while rendering equivalently.
- Acceptance contracts: `OIC-012`.
- Owner and decision date: `satoshikawato`, `2026-08-28`.

### PD-OI-014: Direct non-adjacent Similarity and Collinear display links

- Concern key: `diagram-generation.non-adjacent-grouping-topology`
- Scenario revision: `1`
- Status: `DEFERRED`
- Normative outcome: Similarity groups continue to use all-vs-all
  protein-search evidence across every loaded record. Collinear continues to
  support both **Adjacent pairs** and **All records** evidence scope. An
  uploaded BLAST TSV remains assignable to every supported comparison edge.
  The only deferred capability is drawing a direct link or ribbon between
  arbitrary user-selected non-adjacent display rows for a Similarity-group or
  Collinear result. This Web display deferral does not make the typed API as a
  whole unsupported.
- Rationale: Direct links between arbitrary non-adjacent result rows require
  separate interaction, persistence, and recovery design. Existing evidence
  scopes and supported comparison-edge inputs are independent capabilities.
- Must preserve: Similarity-group all-vs-all evidence; both Collinear evidence
  scopes; uploaded BLAST TSV assignment to supported comparison edges;
  existing adjacent/all display workflows; accurate typed API documentation.
- May retire: Claims that the typed API lacks a capability solely because the
  Web does not expose it.
- Accepted residual risk: The Web cannot directly draw a Similarity-group or
  Collinear link/ribbon between arbitrary user-selected non-adjacent display
  rows until that interaction is designed.
- Acceptance contracts: `OIC-012`.
- Owner and decision date: `satoshikawato`, `2026-08-28`.

### PD-OI-015: Advanced GUI exposure policy

- Concern key: `diagram-generation.gui-exposure-policy`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Normative outcome: API field existence alone does not require a GUI control.
  A control is added for a defined journey or credible usage evidence. Valid
  GUI-unmanaged values use `READ_ONLY` or `PASS_THROUGH` behavior.
- Rationale: GUI completeness is measured by supported journeys, not field
  count.
- Must preserve: Explicit disposition; lossless valid values; accurate
  unsupported/deferred scope.
- May retire: Field-count parity as an acceptance criterion and a
  repository-wide registry used only to enforce it.
- Accepted residual risk: Some advanced values remain non-editable in the Web.
- Acceptance contracts: `OIC-006`, `OIC-008`, `OIC-012`.
- Owner and decision date: `satoshikawato`, `2026-08-28`.

### PD-OI-016: Generation failure, cancellation, and stale-result isolation

- Concern key: `diagram-generation.failure-isolation`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Normative outcome: A failed, canceled, superseded, or stale Generate attempt
  does not replace the last successful Result or committed canonical request.
  The draft remains available, and the UI identifies that the displayed Result
  is the last successful one.
- Rationale: Failed work must not destroy valid output or present candidate
  state as committed.
- Must preserve: User inputs; last successful artifact; actionable error;
  ability to correct and retry.
- May retire: Pre-validation committed mutation, clearing output on failure,
  and partial-state admission.
- Accepted residual risk: Displayed Result can differ from the current draft;
  the UI discloses this.
- Acceptance contracts: `OIC-013`.
- Owner and decision date: `satoshikawato`, `2026-08-28`.

### PD-OI-017: Active comparison appearance controls remain reachable

- Concern key: `diagram-generation.comparison-appearance-reachability`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Normative outcome: When a Linear comparison presentation is active,
  applicable Comparison Appearance controls remain reachable. Supported match
  styles include `ribbon` and `curve`, and match-height adjustment remains
  available where visibly effective. Switching analysis modes does not reset
  independent appearance values. Appearance settings affect rendering, not
  raw protein search or grouping semantics.
- Rationale: A supported rendering choice must remain discoverable and
  operable.
- Must preserve: Control reachability; accessibility; Session round trip;
  independence from analysis defaults and cache stages.
- May retire: Hidden active controls and unauthorized mode-change resets.
- Accepted residual risk: Advanced placement can require one additional
  disclosure action if the control remains discoverable.
- Acceptance contracts: `OIC-006`.
- Owner and decision date: `satoshikawato`, `2026-08-28`.

### PD-OI-018: Complete Linear records, placement, and comparison scope

- Concern key: `diagram-generation.linear-record-universe-and-search-scope`
- Scenario revision: `5`
- Supersedes: `PD-OI-018`, scenario revision `4` (`B / FILE-DATABASE-WITHOUT-UNREQUESTED-SELF`).
- Status: `ACCEPTED`
- Selected outcome: `C / CONSECUTIVE-FILE-ROW-BLOCKS`
- Normative outcome: exactly the complete approved scenario revision `5`
  `PRODUCT_DECISION` receipt and its JSON representation below, together with
  the scenario revision `4` outcome and the scenario revision `3` outcome
  retained below. The revision `5` receipt's Must preserve retains the revision
  `4` receipt and its Must preserve (including items 1–7 of revision `3`). Its
  May retire narrows revision `3` item 6's custom-layout paragraph and the
  revision `3` JSON `customLayoutRule` only for a File whose records occupy
  consecutive rows that no other File shares: such a File moves as one block of
  rows. A File that shares a row with another File, or whose rows interleave
  with another File's rows, stays custom and its move stays unavailable. The
  other revision `3` and `4` terms remain in force.
- Decision source: `satoshikawato`'s answer of `2026-10-04`, quoted verbatim in
  the Revision 31 entry above. The receipt and JSON restate that answer and its
  option text; the Owner did not separately review the receipt wording.
  Dependent runtime (`#780`) requires this authority merged into its base; this
  amendment supplies no runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `ae2b8a91a6c99c0eb2cc4d320bd4d461cbd05c2eb38e5d11c148f52e2702aab4`.
- Acceptance contracts: `OIC-005`, `OIC-006`, `OIC-007`, `OIC-013`, `OIC-015`,
  `OIC-018`.

```text
PRODUCT_DECISION
Concern: diagram-generation.linear-record-universe-and-search-scope
Scenario revision: 5
Supersedes: PD-OI-018, scenario revision 4
Choice: C / CONSECUTIVE-FILE-ROW-BLOCKS
Rationale: CLI で書いた Session は 1 record = 1 行で並ぶため、複数 record のファイルが複数の行にまたがり、Web で Load すると File の並べ替えができず custom layout の警告が出る。各 File の record が連続した行にまとまり、別の File と行を共有していなければ、その行のまとまりごと File を移動できるようにする。
Must preserve: revision 4 の receipt とその Must preserve（revision 3 の item 1〜7 を含む）。File の移動はその File のすべての record を 1 つのまとまりとして動かし、File 内の record の相対的な行・列・順序と record に付いた設定を保つ。移動は 1 回で Undo できる draft 操作である。1 File = 1 行のレイアウトでの動作は変えない。複数の File が 1 行を共有するレイアウトと、File の行が別の File の行と入り組んだレイアウトでは、従来どおり移動できず、Record Layout へ案内する。止めた移動は File の順、行の割り当て、比較、cache、現在の Result を何も変えない。
May retire: 1 つの File が複数の行にまたがるときは File を移動できないという制限（revision 3 の item 6 と customLayoutRule のうち、その File の行が連続していて別の File と共有されていない場合）。
Accepted residual risk: 行のまとまりの並びが File カードの順と違うとき、最初の移動で File カードの順に並べ直される。Arrange in rows が OFF のときも、移動は表示されていない行の割り当てを書き換える（既存の動作）。
Owner: satoshikawato
Decision date: 2026-10-04
```

```json
{
  "concern": "diagram-generation.linear-record-universe-and-search-scope",
  "scenarioRevision": 5,
  "supersedes": "PD-OI-018, scenario revision 4",
  "choice": "C / CONSECUTIVE-FILE-ROW-BLOCKS",
  "rationale": "CLI で書いた Session は 1 record = 1 行で並ぶため、複数 record のファイルが複数の行にまたがり、Web で Load すると File の並べ替えができず custom layout の警告が出る。各 File の record が連続した行にまとまり、別の File と行を共有していなければ、その行のまとまりごと File を移動できるようにする。",
  "mustPreserve": "revision 4 の receipt とその Must preserve（revision 3 の item 1〜7 を含む）。File の移動はその File のすべての record を 1 つのまとまりとして動かし、File 内の record の相対的な行・列・順序と record に付いた設定を保つ。移動は 1 回で Undo できる draft 操作である。1 File = 1 行のレイアウトでの動作は変えない。複数の File が 1 行を共有するレイアウトと、File の行が別の File の行と入り組んだレイアウトでは、従来どおり移動できず、Record Layout へ案内する。止めた移動は File の順、行の割り当て、比較、cache、現在の Result を何も変えない。",
  "mayRetire": "1 つの File が複数の行にまたがるときは File を移動できないという制限（revision 3 の item 6 と customLayoutRule のうち、その File の行が連続していて別の File と共有されていない場合）。",
  "acceptedResidualRisk": "行のまとまりの並びが File カードの順と違うとき、最初の移動で File カードの順に並べ直される。Arrange in rows が OFF のときも、移動は表示されていない行の割り当てを書き換える（既存の動作）。",
  "owner": "satoshikawato",
  "decisionDate": "2026-10-04"
}
```

#### Retained scenario revision `4` outcome

- Scenario revision `4` (`B / FILE-DATABASE-WITHOUT-UNREQUESTED-SELF`)
  normative outcome, retained: exactly the scenario revision `4`
  `PRODUCT_DECISION` receipt and its JSON representation below, together with
  the scenario revision `3` outcome retained below. The revision `4` receipt's
  Must preserve retains items 1–7 of revision `3`; its May retire removes only
  the limitation of self-search exclusion to Collinear inference OFF and the
  loss of links when unrequested self-hits fill Max target seqs.
- Decision source (revision `4`): the complete `D-19` receipt in [`02_DECISION_PACK.md`](./web-gui-audit-20260930/02_DECISION_PACK.md) at P00 merge
  commit `e97d90fecfb327135eb50e85cfa9a87be3145823`, approved by `satoshikawato` on `2026-09-30` through
  the two Owner replies quoted verbatim in the Revision 29 entry above. The
  receipt and JSON reproduce all supplied fields without translation or
  additional terms. The CLI LOSAT database scope is unchanged (`D-40` in the
  same Decision Pack).
- Receipt SHA-256 (revision `4`, UTF-8, excluding the final newline):
  `39df7d14e08db9cf2872da02436c49a11bb7e2d2105cf256dfe0a86fd0786dbb`.

```text
PRODUCT_DECISION
Concern: diagram-generation.linear-record-universe-and-search-scope
Scenario revision: 4
Supersedes: PD-OI-018, scenario revision 3
Choice: B / FILE-DATABASE-WITHOUT-UNREQUESTED-SELF
Rationale: 1 つのファイルを 1 つのゲノムとして扱い、複数 replicon のゲノムでもゲノム単位の E-value で比較する。そのうえで、要求していない自己一致が Max target seqs を埋めてリンクが消えることを防ぐ。
Must preserve: revision 3 の item 1〜7（source ファイル単位の batch、database の範囲を raw cache の identity に含めることを含む）。表示に使わない自己検索（record 自身への検索）は、どのモードでも実行しない。同じ source の中の record どうしの比較は、query の record を除いた database で検索する。E-value の database は subject 側の source ファイル（query と同じ source なら query の record を除いたもの）であることを、docs と Run Info に明記する。job 数の見積もりは、実行の計画と同じ関数から出す。
May retire: 自己検索を除くのが Collinear の inference が OFF のときだけという限定と、要求していない自己一致で Max target seqs が埋まりリンクが消える動作。
Accepted residual risk: 同じ record でも、ファイルの分け方を変えると E-value が変わる（1 ファイル = 1 ゲノムという前提）。record の対ごとに検索する CLI とは、複数 record のファイルで E-value が一致しない。同じ source の中の比較は別の job になり、job 数が増える。
Owner: satoshikawato
Decision date: 2026-09-30
```

```json
{
  "concern": "diagram-generation.linear-record-universe-and-search-scope",
  "scenarioRevision": 4,
  "supersedes": "PD-OI-018, scenario revision 3",
  "choice": "B / FILE-DATABASE-WITHOUT-UNREQUESTED-SELF",
  "rationale": "1 つのファイルを 1 つのゲノムとして扱い、複数 replicon のゲノムでもゲノム単位の E-value で比較する。そのうえで、要求していない自己一致が Max target seqs を埋めてリンクが消えることを防ぐ。",
  "mustPreserve": "revision 3 の item 1〜7（source ファイル単位の batch、database の範囲を raw cache の identity に含めることを含む）。表示に使わない自己検索（record 自身への検索）は、どのモードでも実行しない。同じ source の中の record どうしの比較は、query の record を除いた database で検索する。E-value の database は subject 側の source ファイル（query と同じ source なら query の record を除いたもの）であることを、docs と Run Info に明記する。job 数の見積もりは、実行の計画と同じ関数から出す。",
  "mayRetire": "自己検索を除くのが Collinear の inference が OFF のときだけという限定と、要求していない自己一致で Max target seqs が埋まりリンクが消える動作。",
  "acceptedResidualRisk": "同じ record でも、ファイルの分け方を変えると E-value が変わる（1 ファイル = 1 ゲノムという前提）。record の対ごとに検索する CLI とは、複数 record のファイルで E-value が一致しない。同じ source の中の比較は別の job になり、job 数が増える。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-30"
}
```

#### Retained scenario revision `3` outcome

- Scenario revision `3` normative outcome, retained:
  1. Without an explicit record selector or crop, Linear includes every record
     from each GenBank or paired GFF3/FASTA source. Enabling comparisons does
     not shrink that set or select only the first record.
  2. Record Layout exposes one independently editable placement per biological
     record, with an identifiable record label. File count is not record count.
     Exposing the controls preserves the existing row placement. Each uploaded
     source has one file-input card regardless of its record count; GFF3 and its
     paired FASTA remain one source. Per-record controls belong under that
     source or in Record Layout and must not appear as repeated file uploads.
     Each multi-record source's record list starts collapsed behind a compact, single-line
     `Number of records: N` disclosure. Expanding it exposes every record's
     controls; collapsing it changes no record selection or placement.
     Removing or replacing a source updates all of its records without leaving
     hidden records from the former source. Save and fresh Load preserve this
     distinction between source files and biological records.
  3. Adjacent comparison selects the Cartesian product between neighboring
     occupied display rows. Two records above three records means six pairs.
     This applies to nucleotide, translated-nucleotide, and Pairwise protein
     comparisons. Explicit selected subsets and uploaded TSV associations
     retain their meanings.
  4. Similarity groups and Collinear with scope `all` retain every ordered
     between-record comparison, including same-row and non-adjacent pairs.
     Similarity groups always includes within-record evidence. Collinear
     includes it only when **Infer orthogroups with self-comparisons** is ON,
     as specified by `PD-OI-021`; OFF excludes it. Display placement, including
     a single occupied row, does not restrict the selected search scope.
  5. Explicit comparison endpoints stay explicit through decoding and
     rendering, including endpoints whose numeric indices are consecutive.
  6. A File card number and its up/down controls represent both input-source
     order and visual row order in a normal Linear layout. Moving a File moves
     every record belonging to that source as one block. All records in the
     moved File remain in the same File-owned row, and File-owned rows follow
     File-card order. The move is one atomic, undoable draft operation.
     Record identity, source association, selector, crop, reverse complement,
     definition, subtitle, depth binding, feature state, and within-File
     record order remain attached to the same record.

     Scenario revision `5` narrows the custom-layout paragraph below: a File
     whose records occupy consecutive rows that no other File shares is
     movable as one block of rows.

     A normal layout has exactly one occupied row per File and no row shared
     by records from another File. When one File spans multiple rows or one row
     is shared by multiple Files, the File move is unavailable. The application
     explains that custom Record Layout controls visual placement and directs
     the user to Record Layout. A blocked move makes no partial change to File
     order, row assignments, comparisons, cache metadata, or the current
     Result.

     Explicit comparison endpoints remain attached by stable record UID and
     are reindexed without changing their biological endpoints. Adjacent
     comparison is resolved from the new occupied-row adjacency after the File
     move. Compatible raw LOSAT evidence remains reusable; incompatible
     derived comparison artifacts are recalculated.

     Save, fresh Load, regeneration, keyboard operation, and Session replay
     preserve the resulting File and row order. The last successful Result
     remains unchanged until Generate succeeds. Failed, canceled, superseded,
     or stale Generate does not replace it. Explicit selection/cropping,
     comparison omission, and imported read-only intent remain supported.
     `OIPC-C07` governs failed, canceled, and stale work. No new File-order
     state, Session field, compatibility migration, request schema, Worker
     protocol, or rendering path is introduced.
  7. Shared source bytes remain shared. Complete comparisons may require more
     jobs, but no hidden record or pair cap is permitted.
     LOSAT batches compatible records by input source file. A source job searches
     multi-sequence FASTA inputs containing the selected records, then routes
     hits to their record endpoints. Two sources
     in all-vs-all with within-record evidence enabled require four directed
     source jobs, including self and reverse searches; eight total records must
     not become 64 LOSAT invocations merely because records are expanded.
     TLOSATX records with different explicit genetic codes use compatible
     subsets within each source because each invocation accepts one query and
     one subject genetic code; records sharing those settings remain batched.
     Adjacent and explicit selections restrict retained record-pair evidence.
     When Collinear inference is OFF, batching must not search any record
     against itself, even inside a multi-record source. Comparisons within one
     source may therefore require separate jobs; compatible between-source
     comparisons remain batched. Search arguments, source contents,
     and the actual searched database scope are part of raw-cache identity.
     Progress reports actual source jobs separately from biological record
     pairs. File batching must preserve cancellation and exact hit routing.
     Drawing-start changes and reverse-complement display reuse raw LOSATP
     results whenever source contents, selected biological regions, and search
     settings are unchanged. Only display coordinates and derived presentation
     are updated; display transforms do not become raw-search identity.
- Rationale: Prevent recurrence of incomplete-record search and per-record
  placement regressions reported by the maintainer. With Arrange in rows
  enabled by default, changing File order while retaining the File's previous
  absolute row produces a visible result that contradicts the explicit reorder
  action. Source order and visual row order therefore change as one operation.
- Must preserve: All seven outcomes above; `PD-OI-006` keeps the fresh Collinear
  scope `adjacent`, and `PD-OI-014` continues to govern evidence versus displayed
  links. A fresh no-comparison document does not gain comparison intent.
  Comparison places Run LOSAT across the top, with No comparison on the left
  and Upload BLAST TSV on the right of the row below. DOM and keyboard Tab
  order follow that visual order on desktop and narrow screens.
  Existing LOSAT Execution, Total threads, Parallel runs, and Threads per run
  controls remain editable in Comparison Settings when LOSAT is active.
  Selecting Run LOSAT opens Settings immediately, exposing the mode and its
  settings without another disclosure click. Restoring active LOSAT intent
  also starts with Settings open. Users can collapse it manually; selecting
  Run LOSAT again reopens it without changing the chosen mode or thread values.
  Fresh and reset Web state defaults Execution to `threaded`. Explicit saved
  `auto`, `serial`, or `threaded` choices remain authoritative on Load.
  The maintainer specified these command order, default, and disclosure
  requirements on `2026-09-14`; the existing execution modes and their support
  checks remain.
- May retire: Mandatory Collinear within-record evidence when inference is OFF;
  preservation of a moved File's absolute numeric row in a normal
  one-File-per-row layout; File-card movement that changes input source order
  while leaving visible row order unchanged; and any acceptance test that
  disables Arrange in rows before proving the primary default-layout
  File-reorder outcome. Other record coverage, placement, and comparison
  capabilities remain supported.
- Accepted residual risk: Increased computation time and memory from complete
  all-record comparisons. The maintainer explicitly accepted this cost as the
  original behavior. This does not permit silent truncation or hidden caps.
  Moving a File can change the derived Adjacent comparison pair set and can
  require comparison recomputation on the next Generate. Custom layouts require
  the user to edit or normalize Record Layout before File-card reordering is
  available. No silent loss of custom placement or comparison intent is
  accepted.
- Acceptance contracts: `OIC-005`, `OIC-006`, `OIC-007`, `OIC-013`, `OIC-015`.
- Original decision source: The maintainer explicitly specified Cartesian Adjacent and
  complete all-record scopes, requested durable regression protection, accepted
  the computation/memory cost, and confirmed `May retire: none`, the Owner, and
  the decision date.
- Amendment source: The maintainer requested Collinear self-comparison and
  orthogroup inference be optional and default OFF on `2026-09-14`.
- File-row amendment source: The maintainer supplied the complete
  `PRODUCT_DECISION` response reproduced below on `2026-09-20`.
- Owner and decision date: `satoshikawato`, `2026-09-20`.

Scenario revision `5` narrows the `customLayoutRule` below in the same way: a File whose records occupy consecutive rows that no other File shares is movable as one block of rows.

```json
{
  "concern": "diagram-generation.linear-record-universe-and-search-scope",
  "scenarioRevision": 3,
  "choice": "LINEAR-FILE-ROW-BLOCK",
  "rationale": "A File card number and its up/down controls communicate the source's visual order in a Linear diagram. With Arrange in rows enabled by default, changing File order while retaining the File's previous absolute row number produces a visible result that contradicts the user's explicit reorder action. File reordering must therefore update source order and visual row order as one operation.",
  "mustPreserve": [
    "Each uploaded GenBank file or paired GFF3/FASTA source remains one File card.",
    "Moving a File moves every record belonging to that source as one block.",
    "In the normal layout, every record from the moved File remains in the same File-owned row, and File rows follow File-card order.",
    "Record identity, source association, selector, crop, reverse complement, definition, subtitle, depth binding, feature state, and within-File record order remain attached to the same record.",
    "Explicit comparison endpoints remain attached by stable record UID and are reindexed without changing their biological endpoints.",
    "Compatible raw LOSAT evidence remains reusable; incompatible derived comparison artifacts are recalculated.",
    "Adjacent comparison is resolved from the new occupied-row adjacency after the File move.",
    "The move is one atomic, undoable draft operation.",
    "The last successful Result remains unchanged until Generate succeeds. Failed, canceled, superseded, or stale Generate does not replace it.",
    "Save, fresh Load, regeneration, keyboard operation, and Session replay preserve the resulting File and row order.",
    "No new File-order state, Session field, compatibility migration, request schema, Worker protocol, or rendering path is introduced."
  ],
  "customLayoutRule": "A normal layout means that every File occupies exactly one row and no row is shared by records from another File. If one File spans multiple rows or one row is shared by multiple Files, the File move is unavailable. The application must explain that the custom Record Layout controls visual placement and direct the user to Record Layout. A blocked move makes no partial change to File order, row assignments, comparisons, cache metadata, or the current Result.",
  "mayRetire": [
    "Preservation of a moved File's absolute numeric row in a normal one-File-per-row layout.",
    "File-card movement that changes input source order while leaving the visible row order unchanged.",
    "Any acceptance test that disables Arrange in rows before proving the primary default-layout File-reorder outcome."
  ],
  "acceptedResidualRisk": "Moving a File can change the derived Adjacent comparison pair set and can therefore require comparison recomputation on the next Generate. Custom layouts require the user to edit or normalize Record Layout before File-card reordering becomes available. No silent loss of custom placement or comparison intent is accepted.",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-20"
}
```

### PD-OI-019: Fresh Web multi-record layout defaults

- Concern key: `diagram-generation.web-multi-record-layout-defaults`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Normative outcome: In a fresh Web document and after Reset Settings, Linear
  defaults to Arrange in rows, grouping records from each source file into the
  same row while retaining independent per-record placement controls. Circular
  defaults to Multiple records in a single canvas (Multi-Record Canvas).
  Repeated record accessions retain distinct feature identities on that canvas.
  Explicitly saved layout choices take precedence on Session import; disabling
  either setting remains supported.
- Rationale: The maintainer explicitly requested these two defaults as part of
  restoring complete multi-record workflows.
- Must preserve: All loaded records; independent placement; explicit saved
  settings; the ability to choose separate rows or disable Circular shared
  canvas; comparison scope as specified in `PD-OI-018`.
- May retire: none. Both existing layout choices remain available.
- Accepted residual risk: Complete all-record comparisons retain the computation
  and memory cost accepted under `PD-OI-018`.
- Acceptance contracts: `OIC-004`, `OIC-006`, `OIC-015`.
- Decision source: Explicit maintainer follow-up requesting same-file Linear
  rows and Circular shared canvas as defaults in the same decision session.
- Owner and decision date: `satoshikawato`, `2026-09-13`.


### PD-OI-020: LOSAT total-thread allocation and control agreement

- Concern key: `diagram-generation.losat-total-thread-allocation`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Normative outcome:
  1. Total threads sets the effective total budget. Safe uses half the available
     hardware threads, with a minimum of one; Available uses the hardware
     count. Explicit numeric budgets cannot exceed the available hardware.
  2. Auto Parallel runs and Auto Threads per run recalculate when the budget
     or pending source-job count changes. Auto Threads per run distributes the
     budget across the selected simultaneous runs; multiple jobs do not impose
     a hidden fixed two-thread limit. With four pending LOSATP source jobs and
     both controls on Auto, numeric totals 32, 16, and 2 produce respectively
     4 runs × 8 threads, 4 × 4, and 2 × 1, when hardware permits those totals.
  3. The displayed plan and actual execution use the same allocation rules.
     Before Generate, the controls may use the estimated source-job count;
     cached jobs are omitted from actual execution. Actual run information
     reports the effective execution allocation. Simultaneous runs multiplied
     by threads per run never exceeds the effective total budget.
  4. Explicit Parallel runs and Threads per run choices remain independent
     editable intent. Auto adjusts around the manual choice. A temporary
     effective clamp does not replace the saved manual value or Auto with its
     computed value; any requested/effective difference is visible. Increasing
     the budget makes a preserved manual choice effective again when feasible.
     Save and fresh Load preserve these choices.
- Rationale: Restore the maintainer-requested linkage between Total threads,
  simultaneous runs, and threads per run, and prevent the controls from
  promising an allocation that execution does not use.
- Must preserve: Existing execution modes, browser support checks, fixed
  single-thread-per-run constraints for LOSATN and TLOSATX, source-file job
  batching, cancellation, and explicit saved choices. Scheduling changes do
  not alter search parameters, member-hit limits, or raw-search identity.
- May retire: none.
- Accepted residual risk: The existing computation and memory allowance in
  `PD-OI-018` remains; this decision grants no exception to the selected total
  thread budget and does not promise a fixed speedup for every workload.
- Acceptance contracts: `OIC-005`, `OIC-013`, `OIC-016`.
- Decision source: Explicit maintainer request to restore or implement linked
  automatic allocation and record this regression as a Contract.
- Owner and decision date: `satoshikawato`, `2026-09-14`.

### PD-OI-021: Optional Collinear self-search and orthogroup inference

- Concern key: `diagram-generation.collinear-optional-orthogroup-inference`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Normative outcome: Web Collinear exposes **Infer orthogroups with
  self-comparisons** as a checkbox, OFF in fresh and reset state. OFF submits
  only between-record searches and builds blocks without orthogroup or paralog
  inference. This exclusion applies even when query and subject records share
  one source file. The chosen anchor, unit, member-hit, threshold, scope, and
  block settings still apply. ON enables within-record evidence and the existing
  orthogroup/paralog inference. Similarity groups retains its inference behavior.
  Session save/load preserves explicit ON/OFF. Older Collinear Sessions that
  omit the choice preserve their historical ON behavior; their saved preview
  and explicit limits are not reinterpreted as fresh defaults. Typed request,
  actual helper call, derived identity, and provenance agree on the choice.
  Matching raw evidence can be reused; inference changes cannot reuse an
  incompatible derived result.
- Rationale: The maintainer identified mandatory self-search and paralog-aware
  grouping during Collinear work and requested a default-OFF checkbox.
- Must preserve: Direct block construction when OFF; existing inference when
  ON; complete cross-record scope; explicit saved choices; raw/derived cache
  correctness. An inactive paralog-link setting does not claim an effect in OFF.
- May retire: Mandatory self-search and orthogroup inference in Web Collinear.
  Both remain available through the checkbox.
- Accepted residual risk: No additional risk waiver was supplied. ON retains
  the existing complete-evidence computation and memory allowance in `PD-OI-018`.
- Acceptance contracts: `OIC-003`, `OIC-005`, `OIC-006`, `OIC-015`, `OIC-018`.
- Decision source: Explicit maintainer request to make Collinear self-comparison
  and orthogroup inference optional with default `false`/OFF.
- Owner and decision date: `satoshikawato`, `2026-09-14`.

### PD-OI-022: Reuse completed LOSAT searches after cancellation

- Concern key: `diagram-generation.completed-losat-search-retry`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Normative outcome: After LOSAT raw search completes, canceling downstream
  Collinear work does not require the same raw search to run again on Generate.
  A member-only change recomputes the dependent analysis and reuses matching
  completed raw evidence. Cache reuse still validates the actual inputs and
  raw search settings; changing Max target seqs can require a new search.
  Clear Cache invalidates retained evidence. Session/History replacement must
  not revive unrelated retry data. Failure or cancellation preserves the last
  successful Result and committed request under `OIPC-C07`.
- Rationale: The maintainer reported a completed LOSATP search restarting after
  canceling Collinear work and changing only member hits, and requested a fix
  with durable regression protection.
- Must preserve: Correct raw and derived identities; deterministic cancellation;
  explicit clearing; last successful output. Partial raw batches are not falsely
  admitted as completed searches.
- May retire: Unnecessary raw reruns caused solely by downstream cancellation
  or a member-limit edit.
- Accepted residual risk: No new persistence or cross-reload guarantee is
  introduced; no additional risk waiver was supplied.
- Acceptance contracts: `OIC-002`, `OIC-005`, `OIC-013`, `OIC-019`.
- Decision source: Explicit maintainer request to fix the reported cancellation
  cache behavior and add Product Decisions/Contracts against regression.
- Owner and decision date: `satoshikawato`, `2026-09-14`.

### PD-OI-023: Lossless protein paths during normal generation and saving

- Concern key: `protein-comparison.path-representation`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `PATH-B / lossless-graph`
- Normative outcome: 通常の生成・保存では、全経路を再現できるgraphを保持する。
  図と解析情報、全経路の内容・順序・ID・shared情報を維持し、対応済み旧ファイルを
  読み込めることと、明示的な旧tuple形式での全経路取得を維持する。
  通常のAPI戻り値と保存形式が常に全経路配列を含む仕様は廃止してよい。
- Decision source: The complete maintainer response to S02 Decision Pack
  revision `1`, prepared at `c46e55d14b6b8fae9e6778a7f04183a4d321185f`.
  The receipt below preserves every supplied field without translating or
  extending its rationale, retirement, risk, owner or date. It is a reviewable
  serialization within this existing authority document, not a new decision
  store or evaluator. Implementation names and schema allocations in the S02
  design remain engineering proposals; this receipt does not freeze them.

```json
{
  "concern": "protein-comparison.path-representation",
  "scenarioRevision": 1,
  "choice": "PATH-B / lossless-graph",
  "rationale": "解析情報を維持しながら、通常の生成・保存での経路展開コストを減らしたい。",
  "mustPreserve": "図と解析情報、全経路の内容・順序・ID・shared情報、対応済み旧ファイルの読み込み、明示的な旧tuple形式での全経路取得。",
  "mayRetire": "通常のAPI戻り値と保存形式が、常に全経路配列を含む仕様。",
  "acceptedResidualRisk": "旧API依存コードの修正、新形式を旧バージョンで読めないこと、明示的な全量取得には大きな時間・メモリが必要になり得ること。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-15"
}
```

### PD-OI-024: Linear Definition alignment and automatic Replicon visibility

- Concern key: `linear.definition-display`
- Scenario revision: `3`
- Supersedes: `PD-OI-024`, scenario revision `2`, only the requirement that
  Linear Layout always explains ON/OFF and application on Generate. D1-B,
  D2-P, and D3-A otherwise remain required in full.
- Status: `ACCEPTED`
- Selected outcomes: `D1-B-WEB-LOCKED-FRESH`, `D2-P`, `D3-A`
- Normative outcome:
  1. **D1-B-WEB-LOCKED-FRESH:** Web fresh/reset selects Lock Definition
     Column=true. Lock=true retains the common left edge and configured gap;
     explicit Lock=false retains common-width centering and row following.
     Single, shared, and mixed rows, the existing accepted `text_anchor`
     domain, explicit saved values, supported old omission meanings, the
     saved Result on Load, and CLI/Python omitted defaults remain supported.
     The Lock Definition Column help-tip (hover, keyboard focus, tap) and the
     checkbox's accessible description explain ON/OFF and application on
     Generate from one text source.
  2. **D2-P — 保存値を保持:** 保存済みSubtitleは、自動／手入力を推測したり、
     Replicon名と文字列が一致したりすることを理由に削除しない。
     読み込みだけでは保存Resultを変えず、Generateで新しい表示契約を適用する。
     不要なSubtitleは利用者が明示的にクリアし、ファイル既定値へ戻る既存の
     継承規則を維持する。手入力Subtitleと行共通・レコード固有ラベルの区別を保つ。
  3. **D3-A — OrganelleもReplicon行で制御:** chromosome、plasmid、organelle
     由来の自動名をShow Repliconで制御し、対象名はオンで一つ、オフでゼロとする。
     候補が競合するときはchromosome→plasmid→organelleの順で一つを選ぶ。
     organelleの表記は既存の自動Subtitle表記を引き継ぐ。
     Web・CLI・Pythonの共通描画に適用し、Show Repliconの既定値falseを維持する。
     自動名のオン／オフは手入力Subtitleの表示を変更しない。
- Decision source: the receipt text approved by `satoshikawato` on
  `2026-09-29` (GUI remediation S00 decision 4) supplies the Rationale and
  the Accepted residual risk. The Owner's confirmed requirement supplies the
  change: move the always-on Lock explanation into the existing help-tip.
  Every other clause repeats scenario revision `2` (signed on `2026-09-26`)
  without change. This serialization adds no other terms.
- Acceptance contracts: `OIC-004`, `OIC-005`, `OIC-006`, `OIC-023`.

```json
{
  "concern": "linear.definition-display",
  "scenarioRevision": 3,
  "choice": "A / D1-B-WEB-LOCKED-FRESH; retain D2-P and D3-A; explain Lock in its help-tip",
  "rationale": "派生 status と常時説明が操作応答を損ない（Result 後の比較切替 約 1.3 s）、画面を圧迫するため削除・help-tip 化する。",
  "mustPreserve": "Lock=trueの共通左端とconfigured gap、明示Lock=falseの共通幅中央とrow追従、単一/共有/混在行、既存text_anchorの受入範囲、保存Sessionの明示値と対応済み旧省略意味、読込時の保存Result、CLI/Python省略default。PD-OI-024のD2-Pの保存/手入力Subtitleと継承・ラベル区別、およびD3-AのReplicon/Organelle選択順・独立制御・既定falseをすべて維持する。Lock Definition ColumnのON/OFFの違いとGenerate適用は、既存help-tip（hover・keyboard focus・tap）とcheckboxのaccessible descriptionで一つのtext sourceから説明する。",
  "mayRetire": "Linear Layoutの常時段落によるON/OFFとGenerate適用の説明だけ。説明はhelp-tipとaccessible descriptionに残す。",
  "acceptedResidualRisk": "help-tip を開かない利用者は Generate/Save/Lock の事前説明を見ない。実 error・Processing/Canceling・recovery は保持する。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-29"
}
```

### PD-OI-025: Linear Depth source scope and discoverability

- Concern key: `diagram-generation.linear-depth-source-scope-and-discoverability`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `FILE-BULK-WITH-RECORD-OVERRIDES`
- Normative outcome:
  1. Each Linear File card exposes the common Depth TSV assignment without
     requiring the user to expand its record list or any record options.
     Applying or clearing a File-level value updates every record binding for
     that File and logical series as one undoable operation.
  2. Sparse per-record Depth bindings remain editable. A File card distinguishes
     empty, common, and mixed record bindings. Applying a File-level value to a
     mixed series replaces every record binding in that File and series; the UI
     discloses this effect before the action.
  3. Logical-series settings shared across records are presented once rather
     than repeated inside every record card. Per-record controls expose only the
     record-specific source assignment. A single-record File does not receive a
     duplicate record-level uploader for the same binding.
  4. The canonical state remains the record-major Depth matrix. Null cells,
     logical series indexes, per-record overrides, source isolation, one-step
     Undo/Redo, Session round trips, regeneration, and canonical render-request
     semantics remain supported. The outcome introduces no new render path or
     requirement for a new persisted Depth-default field.
- Decision source: The complete `PRODUCT_DECISION` response from
  `satoshikawato` dated `2026-09-21` for issue `#554`, reproduced below. The
  receipt preserves the supplied fields without extending its rationale,
  preservation, retirement, risk, owner, or date. This is a reviewable
  serialization in the existing static authority document, not a new decision
  store or a `BD-###` record. It cannot authorize dependent runtime until
  merged into that runtime's base.
- Acceptance contracts: `OIC-004`, `OIC-005`, `OIC-006`, `OIC-020`.

```json
{
  "concern": "diagram-generation.linear-depth-source-scope-and-discoverability",
  "scenarioRevision": 1,
  "choice": "FILE-BULK-WITH-RECORD-OVERRIDES",
  "rationale": "Common Depth TSV input should be available once on the File card, while supported sparse per-record bindings remain editable.",
  "mustPreserve": [
    "Record-major Depth matrices",
    "Null cells and logical series indexes",
    "Per-record overrides",
    "Source isolation",
    "One-step Undo/Redo",
    "Session round trips",
    "Regeneration",
    "Existing canonical request semantics"
  ],
  "mayRetire": [
    "Duplicated global Depth settings inside every record card",
    "The need to expand every record before finding Depth input"
  ],
  "acceptedResidualRisk": "Applying or clearing a File-level value replaces or clears every record binding for that File and series; the UI must disclose this and Undo must restore the previous matrix.",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-21"
}
```

### PD-OI-026: Deterministic Similarity Group anchor resolution with optional review

- Concern key: `diagram-generation.similarity-alignment.anchor-resolution`
- Scenario revision: `3`
- Supersedes: `PD-OI-026`, scenario revision `2` (`A / AUTO_PRESELECTED_SUGGESTIONS`).
- Status: `ACCEPTED`
- Selected outcome: `A / DETERMINISTIC_RESOLUTION_WITH_OPTIONAL_REVIEW`
- Normative outcome: The exact selected reference and stable biological
  identities govern independent target resolution. Python retains
  only-usable-candidate and unique-direct-RBH automatic resolution and
  disclosed transient recommendations for ambiguity. A fully resolved default
  alignment may commit without opening review. An explicit review still
  exposes candidate reasons, replacement, and Skip before commitment;
  ambiguous suggestions still require review. Missing, unusable, and skipped
  targets remain unchanged, and resolution is independent of viewport,
  geometry, confidence, supporting-edge count, and multi-hop evidence.
- Decision source: The complete four Choice A `PRODUCT_DECISION` texts
  explicitly approved by `satoshikawato` on `2026-09-25` for the issue
  `#586` follow-up, with this concern's receipt reproduced below. This
  serialization adds no terms to the approved receipt and becomes runtime
  authority only after merge into the runtime base.

```json
{
  "concern": "diagram-generation.similarity-alignment.anchor-resolution",
  "scenarioRevision": 3,
  "choice": "A / DETERMINISTIC_RESOLUTION_WITH_OPTIONAL_REVIEW",
  "rationale": "Python should keep its deterministic automatic resolutions and disclosed ambiguity recommendations, while a fully resolved default alignment can commit without opening candidate review unless the user requests it.",
  "mustPreserve": "The exact selected reference; stable record and biological-feature identity; only-usable-candidate and unique-direct-RBH automatic resolution; independent displayed-record treatment; unchanged missing, unusable, and skipped records; visible recommendation reasons and local replacement or Skip in every opened review; final shared Python validation; and independence from viewport, scroll, ribbon geometry, confidence score, supporting-edge count, and multi-hop evidence. Unique representative and deterministic candidate 1 remain transient recommendations, not automatic resolutions.",
  "mayRetire": "The requirement that replacement or Skip be presented before every Python-resolved default alignment. The explicit review action must still present those choices before commitment.",
  "acceptedResidualRisk": "A resolved default plan may commit without per-target inspection. The named explicit review action and reliable Undo/Reset must remain available; ambiguous suggestions still require the full disclosed review.",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-25"
}
```

### PD-OI-027: Explicit display directions with record-owned orientation

- Concern key: `diagram-generation.similarity-alignment.transform-semantics`
- Scenario revision: `5`
- Supersedes: `PD-OI-027`, scenario revision `4` (`A / RECORD_OWNED_ORIENTATION_WITH_REVIEW_MATCH`).
- Status: `ACCEPTED`
- Selected outcome: `A / EXPLICIT_DISPLAY_DIRECTION_MODES`
- Normative outcome: exactly the approved `PRODUCT_DECISION` receipt below.
- Decision source: `satoshikawato` explicitly approved the four complete
  issue `#598` direction and Reset receipts on `2026-09-26`. The complete
  approval receipt for this concern is reproduced below.
  This serialization adds no terms to the accepted receipt. It becomes
  dependent runtime authority only after merge into the runtime base.

```json
{
  "concern": "diagram-generation.similarity-alignment.transform-semantics",
  "scenarioRevision": 5,
  "choice": "A / EXPLICIT_DISPLAY_DIRECTION_MODES",
  "rationale": "Users should choose the final displayed direction of selected alignment features instead of making all targets follow a potentially minority-direction reference. Reference identity defines positioning, while record direction remains independent record state.",
  "mustPreserve": "Default Keep current directions; one exclusive Keep, All selected features right-facing, All selected features left-facing or Custom direction mode; Custom per-record Keep/right/left choices; exact reference identity and selected target anchors; bulk direction scope including the reference and only known-strand selected anchors; explicit unchanged reasons for unknown directions and unchanged missing, unusable and skipped records; record-wide absolute orientation updates without editing biological source strands; a fixed pre-Align canvas x of the reference feature center, adjusted record placement as necessary, unchanged vertical placement and exact idempotent selected-anchor center alignment; per-record before/after arrow previews and truthful scope coverage; readable text and feature/label/annotation/ribbon geometry in the same record transform; one atomic validated orientations/placement/plan/Result commit; orientation-independent plans, ordinary Reverse, accepted Reset scope and complete artifact Undo/Redo. Changed final validated directions or reference-center placement require a refreshed preview and another Apply before commitment. Capture actual direction deltas including the reference if a separately authorized reset receipt is enabled.",
  "mayRetire": "The rule that alignment always preserves the reference record's direction and left-edge position; the single reference-relative Match checkbox; the restriction that target exceptions require post-Apply sidebar Reverse. Do not retire exact reference selection, its fixed anchor-center position, the default Keep mode or ordinary Reverse.",
  "acceptedResidualRisk": "Users can confuse display arrows with biological strand annotations or interpret all as every feature in a source. Show selected-anchor scope, reference participation, record names, before/after arrows and unknown exclusions; preserve sources and preview reference left-edge movement while its feature center stays fixed. Custom increases review density and requires keyboard/390 px acceptance. Majority inference, guessed unknown directions, source annotation changes, hidden flips and separate orientation/render owners are not accepted.",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-26"
}
```

### PD-OI-028: Orientation-independent active alignment plan

- Concern key: `diagram-generation.similarity-alignment.plan-lifecycle`
- Scenario revision: `3`
- Supersedes: `PD-OI-028`, scenario revision `2`
  (`A / ORIENTATION_INDEPENDENT_ACTIVE_PLAN`), only by adding the scenario-3
  outcome for regrouping comparison settings; every scenario-2 term stays in
  force. Scenario revision `2` superseded scenario revision `1`
  (`A / PERSISTED_ACTIVE_ALIGNMENT_PLAN`).
- Status: `ACCEPTED`
- Selected outcomes: `A / ORIENTATION_INDEPENDENT_ACTIVE_PLAN` (scenario
  revision `2`) and `A / KEEP_ANCHOR_PLAN` (scenario revision `3`)
- Normative outcome: exactly the two approved `PRODUCT_DECISION` receipts
  below, read together.
- Scenario of revision `3`: an active plan exists; the user changes a
  Generate-time comparison setting (Bitscore, E-value, Minimum identity,
  Minimum length, or a LOSATP grouping setting) so that a Generate regroups
  and the plan's reference or a target is no longer a unique member of the
  plan's Similarity Group; every anchor feature still exists in its record.
- Decision source (scenario revision `2`): The complete five `PRODUCT_DECISION` texts for the
  record-owned orientation follow-up to issue `#586`, explicitly approved
  as written by `satoshikawato` on `2026-09-26`, with this concern's
  receipt reproduced below. This serialization adds no terms to the
  approved receipt and becomes runtime authority only after merge into
  the runtime base.

```json
{
  "concern": "diagram-generation.similarity-alignment.plan-lifecycle",
  "scenarioRevision": 2,
  "choice": "A / ORIENTATION_INDEPENDENT_ACTIVE_PLAN",
  "rationale": "An active plan identifies anchors, not directions. A manual orientation change leaves every anchor valid, so keeping the plan spares the user a second alignment after reversing a record.",
  "mustPreserve": "Survival of the active plan across regeneration after style, label, and canvas-size changes, and across record reorder when stable record identities remain; survival across manual record orientation changes, with the same anchors aligned in the new orientation when the diagram is next generated; explicit clearing and notification after manual record movement, source replacement, crop changes, or record-selector changes; validation before regeneration; the last successful Result while a stale plan is repaired; and explicit reselect, Skip, or Clear actions without automatic anchor substitution.",
  "mayRetire": "Clearing the active plan when the user manually changes a record's orientation.",
  "acceptedResidualRisk": "After a manual orientation change, the reversed record moves horizontally on the next generation so that its anchor stays aligned. A user who wanted the previous position uses Undo or Reset Align.",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-26"
}
```

Scenario revision `3` receipt, approved by `satoshikawato` on `2026-10-11`
(OV-341):

```json
{
  "concern": "diagram-generation.similarity-alignment.plan-lifecycle",
  "scenarioRevision": 3,
  "choice": "A / KEEP_ANCHOR_PLAN",
  "rationale": "An alignment plan identifies anchor features, not groups (PD-OI-028); regrouping alone does not make it stale.",
  "mustPreserve": "no anchor substitution; repair still blocks when an anchor feature is gone; a no-edit Generate after a successful Generate succeeds.",
  "mayRetire": "none",
  "acceptedResidualRisk": "the drawer's Active plan may name a group that no longer lists the reference until the user aligns again.",
  "owner": "satoshikawato",
  "decisionDate": "2026-10-11"
}
```

### PD-OI-029: Selectable position and alignment-direction Reset with artifact history

- Concern key: `diagram-generation.similarity-alignment.reset-and-history`
- Scenario revision: `3`
- Supersedes: `PD-OI-029`, scenario revision `2` (`A / IMMEDIATE_PREALIGN_POSITION_BASELINE`).
- Status: `ACCEPTED`
- Selected outcome: `A / SELECTABLE_RESET_WITH_ALIGNMENT_DIRECTION_RESTORE`
- Normative outcome: exactly the approved `PRODUCT_DECISION` receipt below.
- Decision source: `satoshikawato` explicitly approved the four complete
  issue `#598` direction and Reset receipts on `2026-09-26`. The complete
  approval receipt for this concern is reproduced below.
  This serialization adds no terms to the accepted receipt. It becomes
  dependent runtime authority only after merge into the runtime base.

```json
{
  "concern": "diagram-generation.similarity-alignment.reset-and-history",
  "scenarioRevision": 3,
  "choice": "A / SELECTABLE_RESET_WITH_ALIGNMENT_DIRECTION_RESTORE",
  "rationale": "Alignment can change both placement and record direction, so users must be able to choose whether Reset removes only positioning or also restores the directions actually changed by the latest Align. The restoration scope must be explicit without turning alignment plans into direction owners or replaying unrelated edits.",
  "mustPreserve": "An explicit default Reset positions command that restores positions immediately before the latest successful Align, clears its active plan and keeps all current directions; an additional Reset positions and alignment direction changes command that uses the same position baseline and restores absolute before-Align direction only for records whose direction actually changed in that Align; unchanged direction for every record not reversed by that Align, including the reference when unchanged and all later manual edits on unaffected records; include a reference record in restoration only if an independently authorized alignment direction choice actually reversed it; visible target names, count, current and restored directions, and disclosure that later manual direction edits on restoration targets are replaced by the combined command; replacement of the reset receipt by each new successful Align, with no first-Align or source-orientation fallback; an orientation-independent plan and ordinary record-owned current orientation; source-bound validated restoration information captured from actual successful before/after states, retained across style regeneration and stable reorder, saved and freshly loaded with new Sessions, and cleared atomically with plan invalidation or successful Reset. Missing old restoration information leaves positions Reset available and direction restoration unavailable with an explicit reason, never guessed; an empty modern delta means no Align direction changes. Both Reset scopes consume the active plan and receipt, preserve unrelated settings and pending form edits, and use one canonical artifact transaction. Undo/Redo restores or reapplies complete artifacts including directions, plan and receipt; failed, canceled, stale, superseded, preview-readiness or History-finalization work creates no committed history and preserves the prior artifact and receipt. Original sources, biological identities, target-external settings, readable text, record-consistent feature/label/ribbon geometry, existing canonical request, rendering, Worker, sanitizer/admission and History owners, compatible comparison reuse and zero additional Reset LOSAT jobs remain intact. Restoration receipts never select rendering orientation.",
  "mayRetire": "The rule that Reset always retains directions changed by alignment and that their restoration is available only through ordinary Undo or manual Reverse. Do not retire the position-only choice or either existing recovery workflow.",
  "acceptedResidualRisk": "The combined scope intentionally replaces later manual direction edits on records actually reversed by the latest Align. Bound this to a visible target and before-direction preview, an explicit position-only alternative and one-operation Undo. A positions-only Reset consumes the same active plan and receipt; switching to combined afterward requires Undo of that Reset first. Sessions without historical direction evidence cannot restore it and must show that limit while retaining positions Reset. A compact Session/History receipt adds bounded state-maintenance cost; require source/plan binding, atomic lifecycle, old-information/no-op/re-Align/session/failure/geometry/job regression coverage, and keyboard/390 px acceptance. Guessing missing history, reversing unrelated records, partial restoration and new direction/render/History owners are not accepted.",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-26"
}
```

### PD-OI-030: Reader-only legacy Similarity Group alignment compatibility

- Concern key: `diagram-generation.similarity-alignment.session-compatibility`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `A / LEGACY_READER_ONLY`
- Normative outcome:
  1. New Sessions persist the exact resolved alignment plan and do not write the
     legacy `alignOrthogroupFeature` group-ID string. Saving and loading a new
     Session preserves the exact reference, per-record resolved anchors, skips,
     effective transform intent, and reset baseline required by the other
     accepted decisions.
  2. Existing Sessions remain loadable through a bounded reader-only adapter
     that reproduces their historical implicit group-resolution behavior. The
     legacy resolver is unavailable to new alignment requests and to the normal
     writer path.
  3. A new alignment or supported edit converts the loaded state to the new
     resolved representation. The writer never downgrades a resolved plan to the
     legacy representation.
  4. Malformed or unsupported legacy values produce an explicit error without
     automatic substitution. The last successful Result remains visible when
     migration or validation fails.
- Decision source: The complete `PRODUCT_DECISION` response from
  `satoshikawato` dated `2026-09-22` for issue `#561`, reproduced below.
  The receipt preserves the supplied fields without extending its rationale,
  preservation, retirement, risk, owner, or date. This is a reviewable
  serialization in the existing static authority document, not a new decision
  store or a `BD-###` record. It cannot authorize dependent runtime until
  merged into that runtime's base.

```json
{
  "concern": "diagram-generation.similarity-alignment.session-compatibility",
  "scenarioRevision": 1,
  "choice": "A / LEGACY_READER_ONLY",
  "rationale": "Existing Sessions must remain loadable, but new Sessions must not perpetuate the ambiguous group-ID representation. Compatibility therefore belongs in a bounded reader-only adapter rather than the normal writer and runtime path.",
  "mustPreserve": "Reader-only reproduction of existing Sessions; the exact resolved plan in new Session round trips; explicit errors for malformed or unsupported legacy values; the last successful Result on migration or validation failure; and conversion to the new representation after a new Align or supported edit.",
  "mayRetire": "Writing the legacy alignOrthogroupFeature string in new Sessions; normal-runtime use of the legacy group resolver; and downgrade writing from a resolved plan to the ambiguous legacy representation.",
  "acceptedResidualRisk": "Replaying an old Session remains dependent on an isolated legacy resolver and can retain its historical implicit selection until the user creates a new alignment.",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-22"
}
```

### PD-OI-031: Automatic resolved alignment with explicit display-direction review

- Concern key: `diagram-generation.similarity-alignment.surface-scope`
- Scenario revision: `5`
- Supersedes: `PD-OI-031`, scenario revision `4` (`A / AUTO_APPLY_RESOLVED_WITH_EXPLICIT_REVIEW_SINGLE_MATCH`).
- Status: `ACCEPTED`
- Selected outcome: `A / AUTO_APPLY_WITH_EXPLICIT_DIRECTION_REVIEW`
- Normative outcome: exactly the approved `PRODUCT_DECISION` receipt below.
- Decision source: `satoshikawato` explicitly approved the four complete
  issue `#598` direction and Reset receipts on `2026-09-26`. The complete
  approval receipt for this concern is reproduced below.
  This serialization adds no terms to the accepted receipt. It becomes
  dependent runtime authority only after merge into the runtime base.

```json
{
  "concern": "diagram-generation.similarity-alignment.surface-scope",
  "scenarioRevision": 5,
  "choice": "A / AUTO_APPLY_WITH_EXPLICIT_DIRECTION_REVIEW",
  "rationale": "Keep automatic resolved alignment uncomplicated while opened reviews expose mutually exclusive final display direction outcomes, including reversing only a minority reference through an all-right or all-left choice.",
  "mustPreserve": "Exact popup/drawer reference selection; automatic Python-resolved default alignment in Keep mode without forced review; an accessible ambiguity-required or explicit review with candidate facts, Select/Skip, resolution summary, one exclusive Keep/right/left/Custom direction selection and per-record resulting arrows; clear known-strand selected-anchor scope including the reference and unknown/skipped exclusion; shared typed validation, actionable underlying errors, strict CLI ambiguity rejection, existing CLI/API defaults, ordinary record controls, keyboard operation and the existing accepted narrow-screen palette coverage limitation. Apply persists absolute record transforms and placements, never review policies in the alignment plan.",
  "mayRetire": "An opened review exposing only a single reference-relative Match checkbox, permanent preservation of reference direction during an explicit direction operation, and mandatory post-Apply correction for per-record exceptions. Keep automatic/default and explicit review entry points.",
  "acceptedResidualRisk": "Explicit direction choices have more outcomes than the default path and can move the reference record's left edge while its selected feature center stays fixed. Use arrow-based labels, a single radio group, truthful target previews and Custom disclosure. This decision does not redesign narrow-screen palette coverage or introduce new CLI direction flags.",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-26"
}
```

### PD-OI-032: Feature-popup rotation for one circular record

- Concern key: `diagram-generation.feature-popup-record-rotation`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `A / POPUP-RECORD-ROTATION`
- Normative outcome:
  1. The open feature popup targets exactly one source-bound feature through its
     explicit stable record and biological-feature identity. It never falls
     back to a global selection. The operation is available for a complete,
     effectively circular record in either Circular or Linear diagram mode.
  2. Record actions expose the selected feature's 5-prime base, covered
     midpoint, and 3-prime base; a signed strand-relative offset; optional
     absolute forward orientation; and a distinct feature-end placement.
     Preview and resolution use original source coordinates, exact multipart
     traversal, and non-negative circular wrapping. Feature-end placement does
     not collapse into the 3-prime-base anchor.
  3. Apply updates only the target record's absolute display start and, when
     requested, absolute reverse-complement state. Leaving orientation off
     preserves its current value. Repeating the same operation is idempotent;
     every target-external record and layout value remains unchanged.
  4. Apply derives a target-only candidate from the last committed request, so
     unrelated pending form edits remain pending and are neither applied nor
     discarded. The existing sidebar workflow remains available and resolves
     the same display-transform meaning.
  5. Successful Apply admits the fresh Result, absolute transform, and
     non-authoritative provenance as one artifact-history transaction. One Undo
     or Redo restores or reapplies them together. Cancel and failed, stale, or
     superseded work leave the prior Result, transform, provenance, and History
     unchanged.
  6. New Sessions persist the absolute transform and provenance in Session 44,
     catalog 4. A bounded reader conservatively accepts released catalog 3;
     provenance never becomes rendering authority. Manual display-start or
     orientation changes clear stale anchor provenance without changing the
     effective transform.
  7. Invalid offsets and unsafe, ambiguous, fuzzy, cropped, linear, stale, or
     otherwise unsupported operations expose operation-specific reasons instead
     of truncating, guessing, or substituting another target. Duplicate record
     identifiers and split feature fragments retain stable source-bound
     identity.
  8. Feature, label, tick, depth, statistics, and comparison geometry follow the
     same record display transform. Source sequence, annotation, qualifiers,
     biological identity, and source-file export remain unchanged. Compatible
     LOSAT evidence is reused with zero additional executor jobs for a
     transform-only operation.
  9. Feature search and post-generation continuation remain available. Record
     actions are keyboard-operable in both rich and simple popup surfaces, show
     visible reason text, preserve the search query and stable target across
     Result replacement, and remain usable at a 390 px viewport.
  10. The request remains schema 7, and the existing Worker protocol and
      rendering path remain unchanged. The implementation adds no second
      request owner, Worker path, SVG admission path, History engine, or record
      rotation engine.
- Decision source: The complete `PRODUCT_DECISION` response from
  `satoshikawato` dated `2026-09-22` for issue `#563`, reproduced below.
  The receipt preserves the supplied fields without extending its rationale,
  preservation, retirement, risk, owner, or date. This is a reviewable
  serialization in the existing static authority document, not a new decision
  store or a `BD-###` record. It cannot authorize dependent runtime until
  merged into that runtime's base.

```json
{
  "concern": "diagram-generation.feature-popup-record-rotation",
  "scenarioRevision": 1,
  "choice": "A / POPUP-RECORD-ROTATION",
  "rationale": "Feature popupから対象featureを基準にrecordを直接回転できるようにし、sidebarとの往復や手動座標計算を減らす。source-coordinate preview、target-only適用、atomic Undo/Redoによって、操作結果を予測可能かつ安全にする。",
  "mustPreserve": "Source sequence、annotation、qualifiers、biological identity、対象外recordのtransformとlayout、未適用のform edits、既存sidebar workflow、canonical request owner、Worker経路、SVG sanitizer/admission経路、ResultとHistoryのowner、RecordDisplayTransform、LOSAT evidence reuse、searchおよびpost-generation workflow、failure/cancel/stale/superseded時の直前Resultとrecord transform。",
  "mayRetire": "none",
  "acceptedResidualRisk": "Popup UIおよびSession catalog compatibility pathの追加に伴う限定的なUI・保守負担を受容する。この負担は既存ownerの再利用、catalog 3からcatalog 4への単一のbounded reader、390 px・keyboard acceptance、AC-01～AC-20、およびfull regression gatesで制限する。科学的意味の変更、source dataの変更、global-selection fallback、追加LOSAT executor jobは受容しない。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-22"
}
```

### PD-OI-033: Feature-popup record-actions presentation

- Concern key: `web.feature-popup.record-actions-presentation`
- Scenario revision: `2`
- Supersedes: `PD-OI-033`, scenario revision `1` (`A / EDIT_DISCLOSURE`).
- Status: `ACCEPTED`
- Selected outcome: `B / LAYOUT_GROUP_DISCLOSURE`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its nine-field JSON representation below.
- Decision source: the Owner answers of `2026-10-04` quoted verbatim in the
  Revision 30 entry above (`OD-1`, `OD-4`), and Appendix A of
  [`POPUP_RECORD_ACTIONS_AND_VIBRIO_SESSION_PLAN_2026-10-04.md`](./POPUP_RECORD_ACTIONS_AND_VIBRIO_SESSION_PLAN_2026-10-04.md).
  The receipt fields restate those answers and the plan presented in that
  session; the Owner did not separately review the receipt wording. Dependent
  runtime requires this authority merged into its base; this amendment supplies
  no runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `2abe2f6f57152e61eb0034d77f6f0a9ac734022fb352814e51512db9ecc61a13`.
- Acceptance contracts: `OIC-021`, `OIC-024`. These obligations and the
  complete selected outcome are jointly required; their citation does not claim
  completed dependent-runtime checks.

```text
PRODUCT_DECISION
Concern: web.feature-popup.record-actions-presentation
Scenario revision: 2
Choice: B / LAYOUT_GROUP_DISCLOSURE
Rationale: popup の Edit を、現在の Result に反映する Appearance と、Generate で適用する Layout の 2 つに分ける。record 回転は Feature placement と同じ Layout グループの開閉セクションに置く。重複した表示と、同じ理由文の繰り返しをなくし、popup を読みやすくする。
Must preserve: 開いた feature だけを対象とする回転、既存の anchor・offset・orientation・feature-end 操作、適用前 preview、操作できない理由の表示（1 か所）、rich と simple の両 popup での同じ操作、keyboard と 390 px での到達性、既存 sidebar 操作、成功時の一体的な Result と Undo/Redo、Cancel・失敗時の直前 Result と record transform、Feature placement の次回 Generate 適用、fill color・stroke・label・feature visibility・legend name・similarity group の操作。
May retire: Edit タブ先頭の Record actions 配置、セクション内の重複見出し・Record と Feature の行・状態バッジ、同じ理由文の複数表示、ヘッダの fill color 入力（Fill Color と重複）、ヘッダの similarity group 行（Similarity group セクションと重複）、タブより上の Feature placement 配置とその個別注記、Edit 上部の「Live edit: …」説明文（グループ見出しで置き換える）。
Accepted residual risk: record 回転は Edit の上部から下へ移るため、見つけにくくなる。Layout 見出しと開閉ボタンの名前で補い、keyboard と 390 px で到達できることを確認する。
Owner: satoshikawato
Decision date: 2026-10-04
```

```json
{
  "concern": "web.feature-popup.record-actions-presentation",
  "scenarioRevision": 2,
  "choice": "B / LAYOUT_GROUP_DISCLOSURE",
  "rationale": "popup の Edit を、現在の Result に反映する Appearance と、Generate で適用する Layout の 2 つに分ける。record 回転は Feature placement と同じ Layout グループの開閉セクションに置く。重複した表示と、同じ理由文の繰り返しをなくし、popup を読みやすくする。",
  "mustPreserve": "開いた feature だけを対象とする回転、既存の anchor・offset・orientation・feature-end 操作、適用前 preview、操作できない理由の表示（1 か所）、rich と simple の両 popup での同じ操作、keyboard と 390 px での到達性、既存 sidebar 操作、成功時の一体的な Result と Undo/Redo、Cancel・失敗時の直前 Result と record transform、Feature placement の次回 Generate 適用、fill color・stroke・label・feature visibility・legend name・similarity group の操作。",
  "mayRetire": "Edit タブ先頭の Record actions 配置、セクション内の重複見出し・Record と Feature の行・状態バッジ、同じ理由文の複数表示、ヘッダの fill color 入力（Fill Color と重複）、ヘッダの similarity group 行（Similarity group セクションと重複）、タブより上の Feature placement 配置とその個別注記、Edit 上部の「Live edit: …」説明文（グループ見出しで置き換える）。",
  "acceptedResidualRisk": "record 回転は Edit の上部から下へ移るため、見つけにくくなる。Layout 見出しと開閉ボタンの名前で補い、keyboard と 390 px で到達できることを確認する。",
  "owner": "satoshikawato",
  "decisionDate": "2026-10-04"
}
```

### PD-OI-034: Local exclusive display-direction review with retry

- Concern key: `web.similarity-alignment.choice-and-retry`
- Scenario revision: `5`
- Supersedes: `PD-OI-034`, scenario revision `4` (`A / AUTO_APPLY_RESOLVED_WITH_REVIEW_RETRY_SINGLE_MATCH`).
- Status: `ACCEPTED`
- Selected outcome: `A / LOCAL_EXCLUSIVE_DIRECTION_REVIEW_WITH_RETRY`
- Normative outcome: exactly the approved `PRODUCT_DECISION` receipt below.
- Decision source: `satoshikawato` explicitly approved the four complete
  issue `#598` direction and Reset receipts on `2026-09-26`. The complete
  approval receipt for this concern is reproduced below.
  This serialization adds no terms to the accepted receipt. It becomes
  dependent runtime authority only after merge into the runtime base.

```json
{
  "concern": "web.similarity-alignment.choice-and-retry",
  "scenarioRevision": 5,
  "choice": "A / LOCAL_EXCLUSIVE_DIRECTION_REVIEW_WITH_RETRY",
  "rationale": "One exclusive direction selection should govern local candidate and direction previews, and retry must retain user intent without combining a global policy with per-row overrides or committing an unseen reference reversal.",
  "mustPreserve": "Local no-Worker mode/custom/candidate/Skip editing; one direction tagged-union state with custom row values only in Custom; Python ownership of candidate eligibility and final source/display facts; one final batch validation per Apply attempt; one resolver for visible preview and final absolute orientations/reference-center placement; unchanged unknown/skipped/missing/unusable targets with reasons; an updated review and another Apply if final validated output differs; editable mode/custom choices after validation or render failure with the underlying message; prior Result, directions, placement and History after failed, canceled, stale or superseded work; canvas interaction, retained initially-Keep review after automatic render failure, atomic artifact Undo/Redo, Session/regeneration and the independently accepted Reset contract. Policies remain transient, and successful committed direction changes are actual record-state deltas, including the reference when changed.",
  "mayRetire": "The single draft-level Match flag, reference-relative direction as the only bulk operation, and post-Apply-only per-record exceptions. Do not add a second validation, rendering or History path or persisted direction policy.",
  "acceptedResidualRisk": "A final changed direction/placement receipt may require another Apply. Keep prior artifacts and local intent, show the new preview and limit re-review to actual output differences. Unknown directions are never guessed and stale reference/source bindings reject explicitly.",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-26"
}
```

### PD-OI-035: Similarity-alignment canvas interaction

- Concern key: `web.similarity-alignment.canvas-interaction`
- Scenario revision: `3`
- Supersedes: `PD-OI-035`, scenario revision `2`
  (`A / RETAIN_MOBILE_PALETTE_COVERAGE`), only its permission for review
  coverage of the 390 px Preview, retirement of simultaneous canvas picking
  and pan/zoom, and continuation requiring review closure to inspect the
  diagram. Those permissions are no longer active.
- Status: `ACCEPTED`
- Normative outcome: the independent scenario-2 preservation requirements
  remain unchanged:

  正確な feature identity に基づく選択、候補一覧からの keyboard radio 選択と Skip、390 px での操作と適切な focus、描画されない候補の一覧からの選択、図上位置だけによる自動選択の禁止、デスクトップでの canvas 選択と手動 pan／zoom、レビューを閉じた後のプレビュー、ガイド・番号・draft が Result・download・Session に混入しないこと。

  The compact review presentation, visible canvas and simultaneous canvas
  interaction are governed by `PD-OI-039`, scenario revision `4`, in full.
  This record and `PD-OI-039` are jointly required; presentation does not
  replace identity, keyboard/Skip, non-rendered-candidate, desktop canvas,
  focus, or overlay-exclusion guarantees. No candidate may be selected from
  canvas position alone. The previous mobile coverage rationale and risk
  are retained in Git history, not as an active coverage exception.
- Decision source: the signed issue `#602` response for
  `web.similarity-alignment.review-presentation`, scenario revision `1`,
  retained in Git history and preserved by `PD-OI-039` scenario revision `2`,
  supplies exactly the rationale, preservation,
  limited retirement, accepted residual risk, owner, and date for this
  limited supersession. No separate human choice or rationale is inferred
  for the canvas-interaction concern. Other transform, plan, Reset, History,
  and alignment outcomes, including `PD-OI-031`/`PD-OI-034`, retain their scope.
  Dependent runtime requires this supersession merged into its base.
- Acceptance contracts: `OIC-006`, `OIC-013`, `OIC-014`, `OIC-026`.

### PD-OI-036: Linear record-label Auto visibility and disclosure

- Concern key: `linear.record-label-auto-visibility`
- Scenario revision: `2`
- Status: `ACCEPTED`
- Selected outcome: `B / AUTO-FRESH-RESET-WITH-DISCLOSURE`
- Scenario revision `2` is the signed scenario; no prior record for this
  concern exists in the base Contract. Fresh/reset remains independently
  Auto for Accession and Length; default Show (01-A) is not selected.
- Normative outcome: exactly the signed `PRODUCT_DECISION` receipt below.
- Decision source: the complete issue `#602` response signed in full by
  `satoshikawato` on `2026-09-26`. This serialization preserves all nine
  supplied fields without translation or extension. It is not a new decision
  store or a `BD-###` record and cannot authorize dependent runtime until
  merged into its base.
- Acceptance contracts: `OIC-004`, `OIC-005`, `OIC-006`, `OIC-022`.

```json
{
  "concern": "linear.record-label-auto-visibility",
  "scenarioRevision": 2,
  "choice": "B / AUTO-FRESH-RESET-WITH-DISCLOSURE",
  "rationale": "共有行の図では簡潔な既定表示を維持し、情報が非表示になる理由とShowへの変更先を配置操作の場所で明示する。",
  "mustPreserve": "fresh/resetの独立Auto、Show/Hideの明示値、diagram-wide Auto解決、休眠行除外、既存Sessionと保存Result、Undo/Redo、GenerateとExportの区別。Auto非表示時はLayoutにも理由・対象field・図全体の範囲・次回Generateの効果・Record Labelsへの変更先を表示する。",
  "mayRetire": "なし。",
  "acceptedResidualRisk": "共有行でAccession/Lengthが非表示になる結果自体は残る。説明を見落とす可能性があるため、layout操作場所とLabelsの両方で実効値と変更先を示す。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-26"
}
```

### PD-OI-037: Operation feedback without derived application status

- Concern key: `web.edit-application-feedback`
- Scenario revision: `2`
- Supersedes: `PD-OI-037`, scenario revision `1` (`A / DERIVED-APPLICATION-STATUS`).
- Status: `ACCEPTED`
- Selected outcome: `A / OPERATION-FEEDBACK-WITHOUT-DERIVED-STATUS`
- Normative outcome: exactly the receipt below.
- Decision source: the receipt text approved by `satoshikawato` on
  `2026-09-29` (GUI remediation S00 decision 4) supplies the Rationale and
  the Accepted residual risk. The Owner's confirmed requirement supplies the
  Choice, Must preserve, and May retire: remove the always-on Pending display
  and the named always-on explanations, and keep Processing/Canceling, real
  errors, and recovery. This serialization adds no other terms. It is not a
  new decision store or a `BD-###` record.
- Acceptance contracts: `OIC-005`, `OIC-013`, `OIC-014`, `OIC-024`.

```json
{
  "concern": "web.edit-application-feedback",
  "scenarioRevision": 2,
  "choice": "A / OPERATION-FEEDBACK-WITHOUT-DERIVED-STATUS",
  "rationale": "派生 status と常時説明が操作応答を損ない（Result 後の比較切替 約 1.3 s）、画面を圧迫するため削除・help-tip 化する。",
  "mustPreserve": "操作単位のLive edit、Applies on Generate、Apply requiredの分類表示、live edit適用中/失敗の通知、Generateボタン内のProcessing/Canceling、実際のerrorとrecovery、canonical即時commitと必要時自動rerender、reviewのlocal draft、atomic Generate、失敗/Cancel/stale時の旧ResultとHistory、Undo/Redo、SessionのResult/draft分離、Exportの現在Result出力を維持する。",
  "mayRetire": "Result上とGenerate上の常時Pending/Applied/Unknown/Invalid/Not generated表示とその派生計算、生成intentの比較基準と更新処理、Generate recalculates placement…、Supported color…、Save stores…の常時説明を退役する。製品の適用タイミング・保存・復旧・編集機能は退役しない。",
  "acceptedResidualRisk": "help-tip を開かない利用者は Generate/Save/Lock の事前説明を見ない。実 error・Processing/Canceling・recovery は保持する。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-29"
}
```

### PD-OI-038: Compact Editor presentation

- Concern key: `web.editor.compact-presentation`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `A / DOCKED-COMPACT-EDITOR`
- Normative outcome: exactly the signed `PRODUCT_DECISION` receipt below.
- Decision source: the complete issue `#602` response signed in full by
  `satoshikawato` on `2026-09-26`. This serialization preserves all nine
  supplied fields without translation or extension. It is not a new decision
  store or a `BD-###` record and cannot authorize dependent runtime until
  merged into its base.
- Acceptance contracts: `OIC-006`, `OIC-013`, `OIC-014`, `OIC-025`.

```json
{
  "concern": "web.editor.compact-presentation",
  "scenarioRevision": 1,
  "choice": "A / DOCKED-COMPACT-EDITOR",
  "rationale": "狭いPreviewでも即時編集の変化を図で確認できるよう、図とEditorを上下の領域へ配置する。",
  "mustPreserve": "同じSVGとEditor、全tabと同期可用性、canonical live commitと必要時rerender、既存History/Session/Export、camera操作、keyboard、Close/Escapeのvisibility-only意味、選択tab、Result置換/失敗復旧。390×844/740では利用可能幅全体かつ高さ200px以上のcanvasを確保し、Editor内容を独立scrollさせ、Close/headerとtoolbarを操作可能にする。短いviewport/soft keyboardでは全操作へscrollで到達できる。wideのside drawerを維持する。",
  "mayRetire": "狭いPreviewでEditorが横から全面高さを覆う表示配置だけ。編集機能や保存意味は退役しない。",
  "acceptedResidualRisk": "上下分割で図とEditor listの縦領域が短くなり、list scrollが増える。実操作のpointer/keyboard/browser検証を必須とし、複製Preview・SVG clone・第二editorによる回避は受け入れない。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-26"
}
```

### PD-OI-039: Compact review with exclusive alignment directions

- Concern key: `web.similarity-alignment.review-presentation`
- Scenario revision: `4`
- Status: `ACCEPTED`
- Selected outcome: `A / DOCKED-REVIEW-AT-EVERY-WIDTH`
- Supersedes: `PD-OI-039`, scenario revision `3`
  (`A / DOCKED-REVIEW-AT-EVERY-WIDTH`, the same outcome), adding one rule: while
  the review shows the figure it fitted, the figure is fitted again whenever
  the width or height of its room in the Result changes, so a zoom or
  position changed by hand during the review lasts only until such a change.
  Every sentence of the scenario revision `3` receipt is kept in the receipt
  below: the review docked at every width, its Editor, feature-popup and view
  terms, the exclusive direction controls and all independent identity,
  keyboard/Skip, non-rendered-candidate, canvas, focus, validation/retry,
  artifact, Session and History requirements remain. What scenario revisions
  `2` and `3` retired stays retired. The scenario-2 mobile coverage exception
  and close-review continuation in `PD-OI-035` remain superseded; its
  scenario-3 independent guarantees remain jointly required. Keeping a zoom
  or position changed by hand during a review after the room changes may
  retire only as stated in this receipt. Non-modal canvas interaction remains
  required.
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  below. The four issue `#598` decisions are unchanged.
- Decision source: the Owner replies of `2026-10-12` and `2026-10-13` quoted
  in the Revision 58 entry above, which states what each quotation leaves
  out. The Owner approved the receipt text below as written. The reply of
  `2026-10-18` quoted in the Revision 61 entry above reads one sentence of
  its Accepted residual risk on an Undo without changing it. The receipt and
  JSON below reproduce its fields without translation or additional terms.
  Dependent runtime requires this authority merged into its base, or merged
  together with it through the reviewed static Product Contract co-change
  route; this record supplies no runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `8f755440a6655699b1341600bd2aba0d17ef80172490fe600a7cf25d5e16b899`.
- Acceptance contracts: `OIC-006`, `OIC-013`, `OIC-014`, `OIC-026`.

```text
PRODUCT_DECISION
Concern: web.similarity-alignment.review-presentation
Scenario revision: 4
Supersedes: PD-OI-039, scenario revision 3 (web.similarity-alignment.review-presentation)
Choice: A / DOCKED-REVIEW-AT-EVERY-WIDTH
Rationale: 候補の番号はrecordごとに振り直されて図上の位置と対応せず、浮動panelは候補を隠すため、reviewをすべての幅でResult内に固定する（広い幅では図の右、狭い幅では図の下）。Apply後の並びをrecordごとの帯で示し、図上では番号の代わりにリング・矢印・reference線で候補を示し、選択が必要なrecordを先頭に置く。review中に図を表示する領域（残りの領域）の大きさが変わると、開始時の大きさのままの図は切れるか余白が残るので、そのたびに図をFitし直す。
Must preserve: PD-OI-031/034と現行transform/plan/reset/historyのすべての結果。resolvedの通常自動Apply、ambiguousと明示reviewのlocal draft、独立Select/Skip、候補根拠とreference identity、Issue #598で承認済みの排他的Keep/right/left/Custom、referenceを含むselected known-strand anchorsの方向選択と各recordのbefore/after矢印、unknown/skipped/missing/unusableの理由付き不変、canvas操作（図上の候補選択と手動pan/zoom）、local編集でWorkerを呼ばないこと、Applyの共有Python batch validationとatomic Result/History。失敗時draft/error/retry、Cancel/stale/superseded時の以前のResult/orientation/History、Session/regeneration/Export、focus復帰を維持する。390×844/740では利用可能幅全体かつ高さ200px以上のcanvasを確保し、候補listをscroll、Apply/Cancelを到達可能にする。すべての幅で、review開始時はEditorをownerで閉じ、tabを保持し、review中は理由付きでopenをdisableし、開始時に開いていたEditorは終了時に同じtabで開き直す。review中はfeature popupを隠し、終了時（Cancel・閉じる・Escape・Apply成功）に戻す。review中、図のクリックはpopupを開かず、切り替えず、閉じない。検索のOpenは表示位置を動かしてpopupの対象を切り替え、そのpopupは終了時に表示する。review開始時に図を残りの領域へFitし、終了時に開始前のzoomと位置へ戻す。review中に保存したSessionとApplyのHistoryには開始前のzoomと位置を残す。リング・矢印・reference線・帯・draftはResult・download・Sessionに混入しない。非モーダルcanvasを維持する。〔review 中の Fit のし直し〕review が Fit した図を表示している間、残りの領域の幅か高さが変わるたびに、Fit ボタンと同じ合わせ方で図をその領域へ Fit し直す。領域が変わるのは、窓の大きさを変えたとき、設定ペインの幅をドラッグしたとき（ドラッグの間、動きに合わせて Fit し続ける）、設定ペインの表示（今の表示とアイコンの列）を切り替えたとき、アイコンの列の横の section を開閉したとき、review の panel が図の右と下のあいだで移ったとき、Result の上の通知やバナーが出たとき・消えたとき。領域が変わらない間は、手で変えた zoom と位置を残す（候補を選ぶ、行にマウスを載せる、キーを押す、Run info のタブに切り替えて戻る）。Fit ボタンは review 中も今までどおり使える。Fit のし直しは History に入らず、panel の内容とリング・矢印・reference 線を変えない。review 中に保存した Session と Apply の History に開始前の zoom と位置を残すこと、終了時に開始前の zoom と位置へ戻すことは変えない。review を開いていないときは、領域が変わっても zoom と位置を変えない。
May retire: review panelの自由drag（すべての幅）、review中にEditorを同時openする継続（すべての幅）、review中にfeature popupを表示する継続、review中に図のクリックでpopupを開く・切り替える継続、候補の番号バッジ、review中に変えたzoomと位置を終了後に残すこと。これ以外のPD-OI-035/039の独立要求とIssue #598の承認済み4決定は退役しない。加えて、review中に変えたzoomと位置を、残りの領域の大きさが変わったあとまで残すことは退役する。
Accepted residual risk: 広い幅でもreview中はEditorとfeature popupを使えず、panelを動かせない。review中に変えたzoomと位置は終了時に失われる。候補は番号ではなく名前・座標・図上のリングで見分けるため、同名の候補は座標と図上の位置で区別する必要がある。review中に変えたzoomと位置は、終了時のほかに、残りの領域の大きさが変わるたびに失われる。窓や設定ペインの操作のほか、Resultの上の通知やバナーが出る・消えるとき、Applyを押したとき（「Applying similarity alignment…」が出ている間、領域が低くなる。Applyが成功すれば開始前のzoomと位置に戻る）にも起きる。スマートフォンでは、ブラウザのツールバーが出入りして領域の高さが変わるときにも起きうる（コードから。実機では確かめていない）。図をドラッグして動かしている途中に領域が変わると、その移動はそこで終わる。reviewを開いたままSessionの読み込み、Undo、CircularとLinearの切り替えでResultが置き換わったあとは、panelが開いていても図はFitし直されない（今と同じ）。
Owner: satoshikawato
Decision date: 2026-10-13
```

```json
{
  "concern": "web.similarity-alignment.review-presentation",
  "scenarioRevision": 4,
  "supersedes": "PD-OI-039, scenario revision 3 (web.similarity-alignment.review-presentation)",
  "choice": "A / DOCKED-REVIEW-AT-EVERY-WIDTH",
  "rationale": "候補の番号はrecordごとに振り直されて図上の位置と対応せず、浮動panelは候補を隠すため、reviewをすべての幅でResult内に固定する（広い幅では図の右、狭い幅では図の下）。Apply後の並びをrecordごとの帯で示し、図上では番号の代わりにリング・矢印・reference線で候補を示し、選択が必要なrecordを先頭に置く。review中に図を表示する領域（残りの領域）の大きさが変わると、開始時の大きさのままの図は切れるか余白が残るので、そのたびに図をFitし直す。",
  "mustPreserve": "PD-OI-031/034と現行transform/plan/reset/historyのすべての結果。resolvedの通常自動Apply、ambiguousと明示reviewのlocal draft、独立Select/Skip、候補根拠とreference identity、Issue #598で承認済みの排他的Keep/right/left/Custom、referenceを含むselected known-strand anchorsの方向選択と各recordのbefore/after矢印、unknown/skipped/missing/unusableの理由付き不変、canvas操作（図上の候補選択と手動pan/zoom）、local編集でWorkerを呼ばないこと、Applyの共有Python batch validationとatomic Result/History。失敗時draft/error/retry、Cancel/stale/superseded時の以前のResult/orientation/History、Session/regeneration/Export、focus復帰を維持する。390×844/740では利用可能幅全体かつ高さ200px以上のcanvasを確保し、候補listをscroll、Apply/Cancelを到達可能にする。すべての幅で、review開始時はEditorをownerで閉じ、tabを保持し、review中は理由付きでopenをdisableし、開始時に開いていたEditorは終了時に同じtabで開き直す。review中はfeature popupを隠し、終了時（Cancel・閉じる・Escape・Apply成功）に戻す。review中、図のクリックはpopupを開かず、切り替えず、閉じない。検索のOpenは表示位置を動かしてpopupの対象を切り替え、そのpopupは終了時に表示する。review開始時に図を残りの領域へFitし、終了時に開始前のzoomと位置へ戻す。review中に保存したSessionとApplyのHistoryには開始前のzoomと位置を残す。リング・矢印・reference線・帯・draftはResult・download・Sessionに混入しない。非モーダルcanvasを維持する。〔review 中の Fit のし直し〕review が Fit した図を表示している間、残りの領域の幅か高さが変わるたびに、Fit ボタンと同じ合わせ方で図をその領域へ Fit し直す。領域が変わるのは、窓の大きさを変えたとき、設定ペインの幅をドラッグしたとき（ドラッグの間、動きに合わせて Fit し続ける）、設定ペインの表示（今の表示とアイコンの列）を切り替えたとき、アイコンの列の横の section を開閉したとき、review の panel が図の右と下のあいだで移ったとき、Result の上の通知やバナーが出たとき・消えたとき。領域が変わらない間は、手で変えた zoom と位置を残す（候補を選ぶ、行にマウスを載せる、キーを押す、Run info のタブに切り替えて戻る）。Fit ボタンは review 中も今までどおり使える。Fit のし直しは History に入らず、panel の内容とリング・矢印・reference 線を変えない。review 中に保存した Session と Apply の History に開始前の zoom と位置を残すこと、終了時に開始前の zoom と位置へ戻すことは変えない。review を開いていないときは、領域が変わっても zoom と位置を変えない。",
  "mayRetire": "review panelの自由drag（すべての幅）、review中にEditorを同時openする継続（すべての幅）、review中にfeature popupを表示する継続、review中に図のクリックでpopupを開く・切り替える継続、候補の番号バッジ、review中に変えたzoomと位置を終了後に残すこと。これ以外のPD-OI-035/039の独立要求とIssue #598の承認済み4決定は退役しない。加えて、review中に変えたzoomと位置を、残りの領域の大きさが変わったあとまで残すことは退役する。",
  "acceptedResidualRisk": "広い幅でもreview中はEditorとfeature popupを使えず、panelを動かせない。review中に変えたzoomと位置は終了時に失われる。候補は番号ではなく名前・座標・図上のリングで見分けるため、同名の候補は座標と図上の位置で区別する必要がある。review中に変えたzoomと位置は、終了時のほかに、残りの領域の大きさが変わるたびに失われる。窓や設定ペインの操作のほか、Resultの上の通知やバナーが出る・消えるとき、Applyを押したとき（「Applying similarity alignment…」が出ている間、領域が低くなる。Applyが成功すれば開始前のzoomと位置に戻る）にも起きる。スマートフォンでは、ブラウザのツールバーが出入りして領域の高さが変わるときにも起きうる（コードから。実機では確かめていない）。図をドラッグして動かしている途中に領域が変わると、その移動はそこで終わる。reviewを開いたままSessionの読み込み、Undo、CircularとLinearの切り替えでResultが置き換わったあとは、panelが開いていても図はFitし直されない（今と同じ）。",
  "owner": "satoshikawato",
  "decisionDate": "2026-10-13"
}
```

### PD-OI-040: Annotation TSV auxiliary columns

- Concern key: `annotations.table-auxiliary-columns`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `A / AUX_COLUMNS_WARN_IGNORE`
- Normative outcome: Web/CLI/Python の Annotation TSV は任意の未知 header を受理し、その列を捨てて既知列のみ import する。１表につき列名を集約して「無視され、Session/TSV 再出力に保存されない」と通知。fill_colour 等の typo も未知列として通知する。header にない余剰 cell、必須欠落・正規化後重複・不正な既知値は全 import 拒否。
- User access / feedback: Web import 操作直後に読み上げ可能な status と列名一覧。CLI logger の集約 warning。通知には cell contents を含めない。
- Session / regeneration: annotation 値のみ保存。Load で未知列復元や自動 Generate をしない。再生成は既知列だけの import と同じ。
- Export / artifact: TSV writer は既存 column inventory のみ。付加列の lossless export はしない。SVG に未知 metadata を入れない。
- Failure / recovery: import は成功、利用者は通知を確認して編集・Generate へ進める。誤字だった場合は原 TSV を修正して再 import。known-invalid 時は直前 state のまま。
- Approval receipt: 「すべて推奨案で承認します。」, selecting exactly this
  complete Choice A outcome, approved by `satoshikawato` on `2026-09-26`.
- Decision source: The complete outcome and `PRODUCT_DECISION` receipt in
  [`APPROVED_PRODUCT_DECISIONS.md`](./issue-600-implementation-20260926/APPROVED_PRODUCT_DECISIONS.md),
  reproduced without inferred rationale, preservation, retirement, or risk
  terms. This is an inert serialization in the existing static authority;
  it cannot authorize dependent runtime until merged into its base.
- Acceptance IDs: `TSV-01`, `TSV-02`, `TSV-03`; their definitions remain owned by
  [`MASTER_PLAN.md`](./issue-600-implementation-20260926/MASTER_PLAN.md).

```json
{
  "concern": "annotations.table-auxiliary-columns",
  "scenarioRevision": 1,
  "choice": "A / AUX_COLUMNS_WARN_IGNORE",
  "rationale": "生物学的な annotation に使う列の意味を検証しつつ、解析 TSV の付加 metadata だけで作図を止めない。取り込まれない列を明示し、利用者が誤字や非保存を判断できるようにする。",
  "mustPreserve": "既知 annotation の値、行/集合順序、strict typed schema、valid-input export、失敗時の既存 draft/Result。 unknown field を typed annotation に通さない。必須列、duplicate、target/known enum/数値/style を検証する。malformed row を付加列と誤認しない。未知列の cell contents を console に出さない。",
  "mayRetire": "Annotation TSV に対する「unknown header はすべて fatal」の契約のみ。records/track 等の他の表の unknown policy は退役しない。",
  "acceptedResidualRisk": "optional typo が無視され、デフォルト style になる可能性。列名と非保存の通知、known-required/known-value strict 検証で範囲を制限。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-26"
}
```

### PD-OI-041: Unmatched annotation feature selectors

- Concern key: `annotations.feature-selector-miss`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `A / MISSING_SELECTOR_SKIP_ROW_WARN`
- Normative outcome: binding 成功済みの annotation に１件でも未一致 feature selector があれば、その record に対する annotation 行全体を skip し、code/set/annotation/record 識別と欠落件数を持つ warning を返す。他の行/record は継続。全注釈 missing でも genome 図は正常に返し、skip 件数を表示する。
- User access / feedback: Generate 成功後の status に skip 件数と row/record 識別を表示。CLI warning、API の structured warning を提供。未一致の qualifier 値を console に dump しない。
- Session / regeneration: selector と row を保存し、次回はそのときの record に再解決。Session Load で自動 Generate をしない。保存 preview は保持する。
- Export / artifact: SVG/PNG/PDF に skipped mark を出さない。annotation TSV は元 row を含み、次の入力で再利用できる。empty mark の legend は作らない。request に含まれた explicit slot と Web の既存自動 projection の配置/gap は維持し、skip を理由に縮小しない。Python が resolved marks から新規 auto slot を作る場合は empty set の slot を作らない。
- Failure / recovery: 成功図を確認して selector を修正・削除・別 record を明示して再 Generate できる。構造エラー時は直前 Result/draft を保ち修正へ。
- Approval receipt: 「すべて推奨案で承認します。」, selecting exactly this
  complete Choice A outcome, approved by `satoshikawato` on `2026-09-26`.
- Decision source: The complete outcome and `PRODUCT_DECISION` receipt in
  [`APPROVED_PRODUCT_DECISIONS.md`](./issue-600-implementation-20260926/APPROVED_PRODUCT_DECISIONS.md),
  reproduced without inferred rationale, preservation, retirement, or risk
  terms. This is an inert serialization in the existing static authority;
  it cannot authorize dependent runtime until merged into its base.
- Acceptance IDs: `SEL-01`, `SEL-02`, `SEL-03`, `SEL-04`, `SEL-05`; their definitions remain owned by
  [`MASTER_PLAN.md`](./issue-600-implementation-20260926/MASTER_PLAN.md).

```json
{
  "concern": "annotations.feature-selector-miss",
  "scenarioRevision": 1,
  "choice": "A / MISSING_SELECTOR_SKIP_ROW_WARN",
  "rationale": "gene の欠落で比較図全体を失敗させず、複数 anchor で指定した annotation の意味も保つ。部分的な範囲の図示を自動で選ばず、欠落行のスキップを利用者に明示する。",
  "mustPreserve": "完全一致行の geometry、既存 record 意味、coordinate policy、crop/reverse/rotation、他の注釈、failure/cancel/stale 隔離。 record 欠落/曖昧/index 範囲外、multi-record の record 省略、malformed selector は fatal。coordinate clip/skip/error、transform、selector matching の意味を維持。任意 exception を skip にしない。",
  "mayRetire": "feature selector miss の blanket fatal だけ。record/syntax/coordinate error の fatal は維持。",
  "acceptedResidualRisk": "gene typo でも図が成功する。skip を表示することで隠れた欠落を防ぐ。一部 anchor が正しくてもその行の有用な mark は表示されない。request に含まれた注釈 slot は空き領域として残り得る。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-26"
}
```

### PD-OI-042: Specific-color caption disambiguation

- Concern key: `styles.specific-color-caption-multiplicity`
- Scenario revision: `2`
- Supersedes: `PD-OI-042`, scenario revision `1`
  (`A / CAPTION_AUTO_DISAMBIGUATE_SOLID_ROWS`).
- Status: `ACCEPTED`
- Selected outcome: `B / NUMBERED-NAMES`, as the Choice field of the receipt
  below states it
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its nine-field JSON representation below.
- Decision source: the Owner replies of `2026-10-13` quoted verbatim in the
  Revision 47 entry above. The receipt and JSON below reproduce all nine
  supplied fields without translation or additional terms. The runtime,
  tests, and documentation that implement it merge with it (static Product
  Contract co-change); this record supplies no runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `506a9f6df4ebd2536a70714ff666ccf1e1acaffa0eb17ac78b115d2b0766350b`.

```text
PRODUCT_DECISION
Concern: styles.specific-color-caption-multiplicity
Scenario revision: 2
Choice: B / NUMBERED-NAMES。同じ名前で見本（塗りと線）が同じ行は 1 つの凡例行にまとめる。同じ名前で見本が違う行が 2 つ以上描かれるときは、その全部の名前に、凡例に並ぶ順で (1), (2), … を付ける（例 Packaging (1) / Packaging (2)）。描画が名付けた行（feature type、other …、track）と同じ名前の rule の行にも同じ規則を使う（例 tRNA (1) / tRNA (2)）。番号は凡例に描く名前だけに付け、rule の caption は書かれたとおりに保つ。
Rationale: 凡例の名前に色の ID を付けない。名前を区別する仕組みは番号付けの 1 つだけにする。利用者が書いた名前を gbdraw が書き換えない。
Must preserve: feature の色、rule の順序と precedence、同名同色の共有、空 caption は凡例なし、使われていない rule を凡例に出さないこと、使った色をすべて凡例に示すこと、Undo/Redo、Session と Generate での保持、CLI・Python API・Web で同じ規則。
May retire: 名前に hex suffix を付ける区別（同名異色の caption、描画が名付けた行と同じ名前の rule）。rule の caption の書き換えとその通知（Updated N specific-color caption(s) to distinguish their colors.）。区別済みの caption を色表の TSV に書き出すこと（書かれたとおりの caption を書き出す）。
Accepted residual risk: 番号は凡例に並ぶ順なので、行が増える・減る・並び替わると番号が変わる。1 行だけになった名前は番号なしに戻る。保存済みの Session と色表で、caption の末尾がその rule 自身の色の ID のとき（例: 色が #ffe700 の rule の Packaging [#ffe700]）は、読み込むときにその ID を外す。利用者が自分でその形の名前を書いていた場合も外れる。
Owner: satoshikawato
Decision date: 2026-10-13
```

```json
{
  "concern": "styles.specific-color-caption-multiplicity",
  "scenarioRevision": 2,
  "choice": "B / NUMBERED-NAMES。同じ名前で見本（塗りと線）が同じ行は 1 つの凡例行にまとめる。同じ名前で見本が違う行が 2 つ以上描かれるときは、その全部の名前に、凡例に並ぶ順で (1), (2), … を付ける（例 Packaging (1) / Packaging (2)）。描画が名付けた行（feature type、other …、track）と同じ名前の rule の行にも同じ規則を使う（例 tRNA (1) / tRNA (2)）。番号は凡例に描く名前だけに付け、rule の caption は書かれたとおりに保つ。",
  "rationale": "凡例の名前に色の ID を付けない。名前を区別する仕組みは番号付けの 1 つだけにする。利用者が書いた名前を gbdraw が書き換えない。",
  "mustPreserve": "feature の色、rule の順序と precedence、同名同色の共有、空 caption は凡例なし、使われていない rule を凡例に出さないこと、使った色をすべて凡例に示すこと、Undo/Redo、Session と Generate での保持、CLI・Python API・Web で同じ規則。",
  "mayRetire": "名前に hex suffix を付ける区別（同名異色の caption、描画が名付けた行と同じ名前の rule）。rule の caption の書き換えとその通知（Updated N specific-color caption(s) to distinguish their colors.）。区別済みの caption を色表の TSV に書き出すこと（書かれたとおりの caption を書き出す）。",
  "acceptedResidualRisk": "番号は凡例に並ぶ順なので、行が増える・減る・並び替わると番号が変わる。1 行だけになった名前は番号なしに戻る。保存済みの Session と色表で、caption の末尾がその rule 自身の色の ID のとき（例: 色が #ffe700 の rule の Packaging [#ffe700]）は、読み込むときにその ID を外す。利用者が自分でその形の名前を書いていた場合も外れる。",
  "owner": "satoshikawato",
  "decisionDate": "2026-10-13"
}
```

### PD-OI-043: Pure pixel track text inputs

- Concern key: `tracks.pixel-text-input-domain`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `A / PIXEL_TEXT_OPTIONAL_PX`
- Normative outcome: 純 pixel track geometry（Linear height/spacing、Circular inner_gap_px/outer_gap_px）の文字列入口は finite decimal/exponent と optional px（大文字小文字・前後空白可）を受理。trim 後空欄/null は auto。10、10px、10PX、10 px は同値。height は正、gap/spacing は非負。不正文字列/単位/非有限を拒否し、draft で保持して row error を出す。
- User access / feedback: 対象 field の help/placeholder を「px optional」に統一。field 名と正/非負条件を row error と CLI error で示す。
- Session / regeneration: 現行 canonical 型だけ保存。同値 input は同じ geometry。Load は既存 preview を保つ。text acceptance は新 migration ではない。
- Export / artifact: 同値 input の SVG/download は同じ。TSV/CLI の書き出しは既存 canonical 数値形式でよい。
- Failure / recovery: row error の値を編集して再 submission。失敗時に直前 Result を保つ。
- Approval receipt: 「すべて推奨案で承認します。」, selecting exactly this
  complete Choice A outcome, approved by `satoshikawato` on `2026-09-26`.
- Decision source: The complete outcome and `PRODUCT_DECISION` receipt in
  [`APPROVED_PRODUCT_DECISIONS.md`](./issue-600-implementation-20260926/APPROVED_PRODUCT_DECISIONS.md),
  reproduced without inferred rationale, preservation, retirement, or risk
  terms. This is an inert serialization in the existing static authority;
  it cannot authorize dependent runtime until merged into its base.
- Acceptance IDs: `PX-01`, `PX-02`, `PX-03`; their definitions remain owned by
  [`MASTER_PLAN.md`](./issue-600-implementation-20260926/MASTER_PLAN.md).

```json
{
  "concern": "tracks.pixel-text-input-domain",
  "scenarioRevision": 1,
  "choice": "A / PIXEL_TEXT_OPTIONAL_PX",
  "rationale": "利用者が pixel 値を単位付きで paste できるようにし、検証・正規化・request の値を揃える。物理 pixel と factor scalar は分けたまま、保存形式を増やさずに入力の一貫性を改善する。",
  "mustPreserve": "既存 valid 数値、Linear px acceptance、auto、physical pixel 意味、現行 typed request/Session、Circular radius/width factor/%、retired key 拒否。 typed JSON の gaps は数値、Linear は既存 ScalarSpec。Circular ratio/% semantics を保つ。不正値を null/0 化しない。不要な arbitrary CSS unit conversion を作らない。",
  "mayRetire": "pure pixel 対象の without-a-unit restriction と、invalid→null/zero の黙示的変換。一般 dimension input の制限は退役しない。",
  "acceptedResidualRisk": "trim空欄はauto。decimal/exponent以外のJS Number形式を使っていた入力は拒否され得るが、Pythonと一致しない隠れた入力経路を支持しない。scope は listed slot fields に限る。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-26"
}
```

### PD-OI-044: Applicable Circular controls and truthful record discovery

- Concern key: `diagram-generation.circular-transform-discoverability`
- Scenario revision: `4`
- Supersedes: `PD-OI-044`, scenario revision `3` (`A / REVEAL_APPLICABLE_SINGLE_RECORD_CONTROLS`).
  Scenario revision `4` changes two things only: the Session version of the
  project and drawing format in its Rationale, from `48` to `49` (`PD-OI-090`
  scenario revision `4`), and a Session `46` or `47` that holds protein
  comparison results stored on a reversed record, which shows its stored
  figure at once and starts Python once at Load to move those results to the
  search coordinates (Rationale, Must preserve, Accepted residual risk).
  Scenario revision `3` superseded scenario revision `2`, which superseded
  scenario revision `1`.
- Status: `ACCEPTED`
- Selected outcome: `A / REVEAL_APPLICABLE_SINGLE_RECORD_CONTROLS`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its JSON representation below.
- Decision source: the Owner replies quoted verbatim in the Revision 37 entry
  above (R1-4 and the receipt approval), and for scenario revision `3` those
  quoted in the Revision 62 entry above, and for scenario revision `4` those
  quoted in the Revision 64 entry above. The Owner approved the receipt text
  below as written. Its two references to the drawings record carry the
  recorded ID `PD-OI-090` in place of the draft ID `PD-OI-087` (Revision 37).
  Dependent runtime requires this authority merged into its base; this record
  supplies no runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `896078be8ce11203eedaddae8d048ae312bceb46f4f89f63f19ff14baadcb398`.

```text
PRODUCT_DECISION
Concern: diagram-generation.circular-transform-discoverability
Scenario revision: 4
Supersedes: PD-OI-044, scenario revision 3
Choice: A / REVEAL_APPLICABLE_SINGLE_RECORD_CONTROLS
Rationale: source の探索状況と一件用 crop の適用条件を区別して示し、編集可能になった一件用 controls は selection の直後に見えるようにする。通常 upload の自動探索と保存済みプレビューの軽い閲覧を両立する。保存済みプレビューの軽い閲覧は、現在の形式（version 46 以降）の Session で保つ。version 46 と 47 の Session で、逆向きにした record を含むタンパク質比較の保存された結果を持つものも、保存された図はすぐに示し、その結果は Python が読み込み時に一度だけ探索の座標に移す（D-124、D-126）。version 27〜44 の Session は、Session を project と drawing の形式（version 49、PD-OI-090）へ移すために、Python が一度だけ読み込み時に変換する。
Must preserve: valid native upload の自動 record discovery と Generate 前の rotation controls、既存 parser/helper 境界、exact source-bound identity、explicit single/grid/batch、fresh shared-canvas default、saved explicit choices、一件 crop と topology/start/reverse の適用条件、手動 close/expand、元の focus、keyboard/390 px、現在の形式（version 46 以降）の Session の preview-only Load の Python Worker 0（version 46 と 47 の Session で、逆向きにした record を含むタンパク質比較の保存された結果を持つものを除く）、active draft と saved artifact の分離、失敗時の旧 Result、Retry/Replace/Remove/Inspect/Generate の継続。version 27〜44 の Session と、その version 46 と 47 の Session の Load でも、保存された図はすぐに示し、変換に失敗したときは何も変えない（PD-OI-045、PD-OI-046）。
May retire: 適用可能になった一件用 section が常に collapsed で始まる挙動、valid fresh upload に manual Load が必須であるかのような prompt、実行中でない deferred discovery を loading と表す UI。全 record subset editing や複数 source support の選択は含まない。以前の形式（version 27〜44）の Session の Load で Python を起動しないこと。
Accepted residual risk: applicable になった時に一件用 section が展開されて pane 高さが変わる。操作元の focus と scroll anchor を維持し、無関係な更新で再展開しない。科学的意味の変更、先頭 record の自動選択、grouping の自動切替、現在の形式（version 46 以降）の Session の preview-only Load による Python 初期化は受容しない（version 46 と 47 の Session で、逆向きにした record を含むタンパク質比較の保存された結果を持つものを除く）。その Session を開くときは、保存された図をすぐに示したうえで Python を一度起動して結果を探索の座標に移し、その間（約 8 秒）は編集を待つ。どちらの座標にも合わない結果を持つ Session は、読み込みを取り消して理由を示す（PD-OI-046）。新しい形式で保存し直すと起こらない。以前の形式（version 27〜44）の Session を開くときは、保存された図をすぐに示したうえで Python を一度起動して project の形式に移し、その間（約 8 秒）は編集を待つ（PD-OI-090）。新しい形式で保存し直すと起こらない。version 45 は開発版だけで、読まない。
Owner: satoshikawato
Decision date: 2026-10-18
```

```json
{
  "concern": "diagram-generation.circular-transform-discoverability",
  "scenarioRevision": 4,
  "supersedes": "PD-OI-044, scenario revision 3",
  "choice": "A / REVEAL_APPLICABLE_SINGLE_RECORD_CONTROLS",
  "rationale": "source の探索状況と一件用 crop の適用条件を区別して示し、編集可能になった一件用 controls は selection の直後に見えるようにする。通常 upload の自動探索と保存済みプレビューの軽い閲覧を両立する。保存済みプレビューの軽い閲覧は、現在の形式（version 46 以降）の Session で保つ。version 46 と 47 の Session で、逆向きにした record を含むタンパク質比較の保存された結果を持つものも、保存された図はすぐに示し、その結果は Python が読み込み時に一度だけ探索の座標に移す（D-124、D-126）。version 27〜44 の Session は、Session を project と drawing の形式（version 49、PD-OI-090）へ移すために、Python が一度だけ読み込み時に変換する。",
  "mustPreserve": "valid native upload の自動 record discovery と Generate 前の rotation controls、既存 parser/helper 境界、exact source-bound identity、explicit single/grid/batch、fresh shared-canvas default、saved explicit choices、一件 crop と topology/start/reverse の適用条件、手動 close/expand、元の focus、keyboard/390 px、現在の形式（version 46 以降）の Session の preview-only Load の Python Worker 0（version 46 と 47 の Session で、逆向きにした record を含むタンパク質比較の保存された結果を持つものを除く）、active draft と saved artifact の分離、失敗時の旧 Result、Retry/Replace/Remove/Inspect/Generate の継続。version 27〜44 の Session と、その version 46 と 47 の Session の Load でも、保存された図はすぐに示し、変換に失敗したときは何も変えない（PD-OI-045、PD-OI-046）。",
  "mayRetire": "適用可能になった一件用 section が常に collapsed で始まる挙動、valid fresh upload に manual Load が必須であるかのような prompt、実行中でない deferred discovery を loading と表す UI。全 record subset editing や複数 source support の選択は含まない。以前の形式（version 27〜44）の Session の Load で Python を起動しないこと。",
  "acceptedResidualRisk": "applicable になった時に一件用 section が展開されて pane 高さが変わる。操作元の focus と scroll anchor を維持し、無関係な更新で再展開しない。科学的意味の変更、先頭 record の自動選択、grouping の自動切替、現在の形式（version 46 以降）の Session の preview-only Load による Python 初期化は受容しない（version 46 と 47 の Session で、逆向きにした record を含むタンパク質比較の保存された結果を持つものを除く）。その Session を開くときは、保存された図をすぐに示したうえで Python を一度起動して結果を探索の座標に移し、その間（約 8 秒）は編集を待つ。どちらの座標にも合わない結果を持つ Session は、読み込みを取り消して理由を示す（PD-OI-046）。新しい形式で保存し直すと起こらない。以前の形式（version 27〜44）の Session を開くときは、保存された図をすぐに示したうえで Python を一度起動して project の形式に移し、その間（約 8 秒）は編集を待つ（PD-OI-090）。新しい形式で保存し直すと起こらない。version 45 は開発版だけで、読まない。",
  "owner": "satoshikawato",
  "decisionDate": "2026-10-18"
}
```

### PD-OI-045: Exclusive semantic Session operations with responsive browsing

- Concern key: `web.session-operation-consistency`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `A / EXCLUSIVE_SEMANTIC_SESSION_OPERATION`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt below.
- Approval receipt: 「推奨案で承認します」, selecting this complete outcome,
  approved by `satoshikawato` on `2026-09-26`.
- Reviewed outcome SHA-256: `7cef6306aefc5ff935d949efef7f76964eae2d24df8281c6ca371df5c2e87633`.
- Decision source: the complete receipt in
  `docs/internal/issue-597-input-session-implementation-20260926/decisions/03_SESSION_OPERATIONS.md`
  at commit `51a786086dc2777e5fa375e5f94d8b7ac7deeedc` of
  `fix/issue-597-input-session-20260926`. Its nine supplied fields are
  reproduced below without translation or additional terms. The receipt
  document is evidence of the human choice; this record is its static authority.
  Dependent runtime requires this authority merged into its base.

```json
{
  "concern": "web.session-operation-consistency",
  "scenarioRevision": 1,
  "choice": "A / EXCLUSIVE_SEMANTIC_SESSION_OPERATION",
  "rationale": "Save/Load は一つの整合した document に対する操作として完了させ、異なる時点の Result、draft、source、cache を混合しない。処理中は閲覧を維持し、semantic edits を終了後に再開する明確な workflow を優先する。",
  "mustPreserve": "主スレッドの応答と閲覧・scroll・pan/zoom・検索、visible pending/busy reasons、同時 Save の join と一度の download、title/size/repeat-download 取消、committed Result と editable draft の分離、supported Sessions と settings-only、atomic Load、failed/canceled/stale/teardown recovery、旧 request/resources/Result/History、source bytesと全 cache/evidence/provenance、JSON/gzip と CLI/Python replay、privacyとsize/sanitization constraints、現行 performance gates。",
  "mayRetire": "Save/Load pending 中の source/editor/History/Reset/Generate 等の semantic mutation と、mutation entry point によって偶然編集可能または silent no-op になる振る舞い。Generate/automatic reflow 中の Save/Load 開始も busy reason 付きで停止し、完了後の再試行を提供する。通常編集・閲覧・成功後の操作は廃止しない。",
  "acceptedResidualRisk": "長い Save/Load の間、document 編集は一時停止する。閲覧、status、bounded completion、error/retry を維持する。無期限 lock、main-thread freeze、データの省略、checkpoint混合、追加memoryの未計測、既存gateの弱化は受容しない。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-26"
}
```

### PD-OI-046: Guidance with bounded diagnostics

- Concern key: `web.errors.diagnostic-disclosure`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `A / GUIDANCE_WITH_BOUNDED_DIAGNOSTICS`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its nine-field JSON representation below.
- Decision source: the complete Choice A receipt approved by `satoshikawato`
  on `2026-09-26`, retained in
  [`DECISION_01_ERROR_DISCLOSURE.md`](./issue-601-bug15-bug19-implementation-20260926/decisions/DECISION_01_ERROR_DISCLOSURE.md)
  at S00 commit `202fe9de554aaa70dc731deb80bf032f26d80061`.
  All nine supplied fields are reproduced without translation or additional
  terms. This static record does not change the source receipt's concern key
  or retire another receipt. Dependent runtime requires this authority merged
  into its base; this amendment supplies no runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `8617888cb1838521f78640db4653e2dcf88dd427cd466bbbc33d59efacaceae8`.

```text
PRODUCT_DECISION
Concern: web.errors.diagnostic-disclosure
Scenario revision: 1
Choice: A / GUIDANCE_WITH_BOUNDED_DIAGNOSTICS
Rationale: 利用者が短い修正案内から作業を続けられ、必要な場合は入力内容を公開せずに安全な失敗種別と段階を調査へ渡せるようにする。
Must preserve: すべての移行対象で既知 validation の修正情報を保持する。Generate/Align の以前の Result、canonical request、draft、orientation、History、retry、Save/Export、cancel/stale/superseded を保つ。Details は keyboard で開け、Copy diagnostics は表示中の bounded code/operation/stage/許可 context/副因だけを手動コピーする。unknown は stage と stable code を示し、元 pattern、sequence、file/record 名、path、SVG、自由な exception/traceback/stdout/stderr を画面・Copy・console に自動公開しない。初回と旧 Result 保持を区別する。
May retire: user-facing raw exception/traceback と自由な stdout/stderr、個別の例外型 prefix の直接表示。安全な修正情報、Details 入口、既存 recovery は退役しない。
Accepted residual risk: bounded 診断だけでは稀な未知例外を特定できず、利用者の明示的な Session 保存・別途再現情報が必要になる場合がある。Clipboard 不可時も表示情報の手動選択コピーと通常 recovery を維持する。
Owner: satoshikawato
Decision date: 2026-09-26
```

```json
{
  "concern": "web.errors.diagnostic-disclosure",
  "scenarioRevision": 1,
  "choice": "A / GUIDANCE_WITH_BOUNDED_DIAGNOSTICS",
  "rationale": "利用者が短い修正案内から作業を続けられ、必要な場合は入力内容を公開せずに安全な失敗種別と段階を調査へ渡せるようにする。",
  "mustPreserve": "すべての移行対象で既知 validation の修正情報を保持する。Generate/Align の以前の Result、canonical request、draft、orientation、History、retry、Save/Export、cancel/stale/superseded を保つ。Details は keyboard で開け、Copy diagnostics は表示中の bounded code/operation/stage/許可 context/副因だけを手動コピーする。unknown は stage と stable code を示し、元 pattern、sequence、file/record 名、path、SVG、自由な exception/traceback/stdout/stderr を画面・Copy・console に自動公開しない。初回と旧 Result 保持を区別する。",
  "mayRetire": "user-facing raw exception/traceback と自由な stdout/stderr、個別の例外型 prefix の直接表示。安全な修正情報、Details 入口、既存 recovery は退役しない。",
  "acceptedResidualRisk": "bounded 診断だけでは稀な未知例外を特定できず、利用者の明示的な Session 保存・別途再現情報が必要になる場合がある。Clipboard 不可時も表示情報の手動選択コピーと通常 recovery を維持する。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-26"
}
```

### PD-OI-047: Rejected Color rule pattern edit recovery

- Concern key: `web.rules.rejected-pattern-edit-recovery`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `A / KEEP_REJECTED_PATTERN_DRAFT`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its nine-field JSON representation below.
- Decision source: the complete Choice A receipt approved by `satoshikawato`
  on `2026-09-26`, retained in
  [`DECISION_02_REGEX_EDIT_RECOVERY.md`](./issue-601-bug15-bug19-implementation-20260926/decisions/DECISION_02_REGEX_EDIT_RECOVERY.md)
  at S00 commit `202fe9de554aaa70dc731deb80bf032f26d80061`.
  All nine supplied fields are reproduced without translation or additional
  terms. This static record does not change the source receipt's concern key
  or retire another receipt. Dependent runtime requires this authority merged
  into its base; this amendment supplies no runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `9cc66bc30f97cc995b9b0002b094ba74cfa00a8e77b163c33074689f0e73b496`.

```text
PRODUCT_DECISION
Concern: web.rules.rejected-pattern-edit-recovery
Scenario revision: 1
Choice: A / KEEP_REJECTED_PATTERN_DRAFT
Rationale: 有効なルールと図を保護しながら、入力ミスや一時的な検証失敗から同じ文字列を修正・再試行できるようにする。
Must preserve: Color/Label の Python regex semantics、既存一 Worker・一 preparation・atomic live commit、優先順位、prepared reuse、valid target と Generate の同値、canonical rule/Result/History、stale/cancel 隔離を維持する。対象は既存 Color rule の pattern field。拒否された text を field に保持して原因と Not applied、Save/Generate は last accepted rule、Export は現在 Result を使うことを示す。keyboard の編集/Retry/Revert、正しい syntax/runtime 分類、同 document の drawer close/reopen・一時 mode 切替での draft 保持、row/revision の現在性を保つ。対象 rule の Undo/Redo 置換、row 削除、document/session 成功置換、reset で draft を解放する。成功 edit だけ History へ記録し、未確定 text は Session/diagnostics/console に自動保存・公開しない。TSV、新規 rule、preset、Search の意味は変えない。
May retire: 対象 field の failure 後に未確定 pattern text を無条件で accepted 値へ戻す表示。不正 rule の拒否、last accepted rule の保護、正常な live edit は退役しない。
Accepted residual risk: 表示 text と accepted rule が一時的に異なり、Save/Generate は accepted rule、Export は現在 Result を使う。Not applied と対象説明、Revert を提供する。Session/document や対象 rule の History 置換後に未確定 draft は保持しない。
Owner: satoshikawato
Decision date: 2026-09-26
```

```json
{
  "concern": "web.rules.rejected-pattern-edit-recovery",
  "scenarioRevision": 1,
  "choice": "A / KEEP_REJECTED_PATTERN_DRAFT",
  "rationale": "有効なルールと図を保護しながら、入力ミスや一時的な検証失敗から同じ文字列を修正・再試行できるようにする。",
  "mustPreserve": "Color/Label の Python regex semantics、既存一 Worker・一 preparation・atomic live commit、優先順位、prepared reuse、valid target と Generate の同値、canonical rule/Result/History、stale/cancel 隔離を維持する。対象は既存 Color rule の pattern field。拒否された text を field に保持して原因と Not applied、Save/Generate は last accepted rule、Export は現在 Result を使うことを示す。keyboard の編集/Retry/Revert、正しい syntax/runtime 分類、同 document の drawer close/reopen・一時 mode 切替での draft 保持、row/revision の現在性を保つ。対象 rule の Undo/Redo 置換、row 削除、document/session 成功置換、reset で draft を解放する。成功 edit だけ History へ記録し、未確定 text は Session/diagnostics/console に自動保存・公開しない。TSV、新規 rule、preset、Search の意味は変えない。",
  "mayRetire": "対象 field の failure 後に未確定 pattern text を無条件で accepted 値へ戻す表示。不正 rule の拒否、last accepted rule の保護、正常な live edit は退役しない。",
  "acceptedResidualRisk": "表示 text と accepted rule が一時的に異なり、Save/Generate は accepted rule、Export は現在 Result を使う。Not applied と対象説明、Revert を提供する。Session/document や対象 rule の History 置換後に未確定 draft は保持しない。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-26"
}
```

### PD-OI-048: Circular Width/Radius input representation

- Concern key: `tracks.circular-measure-input-representation`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `A / NUMERIC_PX_FACTOR_WITH_LEGACY_INPUT`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its nine-field JSON representation below.
- Decision source: the complete Choice A text reviewed in
  [01_INPUT_REPRESENTATION.md](https://github.com/satoshikawato/gbdraw/blob/727b876214b58d4233790ce1feb4913eaf48bf4d/docs/internal/issue-619-implementation-plan-20260927/DECISION_PACKS/01_INPUT_REPRESENTATION.md)
  at S00 commit `727b876214b58d4233790ce1feb4913eaf48bf4d`. On `2026-09-27`,
  `satoshikawato` explicitly confirmed signing all three Choice A texts with
  that Owner and Decision date. All nine supplied fields are reproduced
  without translation or additional terms. This record does not supersede
  another decision. Dependent runtime requires this authority merged into
  its base; this amendment supplies no runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `b405f184de5362132bb62451eb0e8ab3d28969809e8a3c9b12cf854d2c1ae55d`.

```text
PRODUCT_DECISION
Concern: tracks.circular-measure-input-representation
Scenario revision: 1
Choice: A / NUMERIC_PX_FACTOR_WITH_LEGACY_INPUT
Rationale: 数値の意味を明示しつつ、通常操作の選択肢をpxと倍率の二つに絞る。percentによる既存入力と保存値の意味は維持する。
Must preserve: 既存px/factor/%の値とunit、precision、typed request/Session、Auto、invaliddraft、draft/Result分離、適用時点、History、privacy、失敗復旧。
May retire: 数値欄内に単位を恒常表示する旧UIと、percentをliteral spellingのまま通常表示することのみ。percent入力や既存Sessionの受理は退役しない。
Accepted residual risk: percent入力を倍率表示へまとめるため、65%が0.65と読めることをhelpで説明する必要がある。suffix入力の途中と確定を区別する。
Owner: satoshikawato
Decision date: 2026-09-27
```

```json
{
  "concern": "tracks.circular-measure-input-representation",
  "scenarioRevision": 1,
  "choice": "A / NUMERIC_PX_FACTOR_WITH_LEGACY_INPUT",
  "rationale": "数値の意味を明示しつつ、通常操作の選択肢をpxと倍率の二つに絞る。percentによる既存入力と保存値の意味は維持する。",
  "mustPreserve": "既存px/factor/%の値とunit、precision、typed request/Session、Auto、invaliddraft、draft/Result分離、適用時点、History、privacy、失敗復旧。",
  "mayRetire": "数値欄内に単位を恒常表示する旧UIと、percentをliteral spellingのまま通常表示することのみ。percent入力や既存Sessionの受理は退役しない。",
  "acceptedResidualRisk": "percent入力を倍率表示へまとめるため、65%が0.65と読めることをhelpで説明する必要がある。suffix入力の途中と確定を区別する。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-27"
}
```

### PD-OI-049: Circular Width/Radius unit changes

- Concern key: `tracks.circular-measure-unit-change`
- Scenario revision: `2`
- Supersedes: `PD-OI-049`, scenario revision `1`. Only the Accepted residual
  risk and the decision date change.
- Status: `ACCEPTED`
- Selected outcome: `A / KEEP_NUMBER_CHANGE_UNIT`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its nine-field JSON representation below.
- Decision source: `PD-OI-037` revision 2 retires the always-on Pending
  display, which removes the Pending half of this record's mitigation. On
  `2026-09-29`, `satoshikawato` explicitly approved the help-only mitigation
  below and kept every other field unchanged. This revision merges with its
  implementation through the reviewed co-change route.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `b879318bb23587d4695706d0102b7e2d0a211229467fcef07658ba407f77c596`.

```text
PRODUCT_DECISION
Concern: tracks.circular-measure-unit-change
Scenario revision: 2
Choice: A / KEEP_NUMBER_CHANGE_UNIT
Rationale: 単位選択を、入力した数値の意味を明示的に変更する編集として統一する。現在の円半径や古いResultに依存せず、生成前や複数recordでも同じ操作を使える。
Must preserve: numericdraftとunitの明示、Auto、Generateまで旧Resultを保つこと、History/Session、失敗復旧、既存scalarの科学的意味、local-only。
May retire: なし。
Accepted residual risk: 倍率からpxへ切り替えると図上の大きさが変わる。helpに、数値維持・次回Generate反映を明記する。
Owner: satoshikawato
Decision date: 2026-09-29
```

```json
{
  "concern": "tracks.circular-measure-unit-change",
  "scenarioRevision": 2,
  "choice": "A / KEEP_NUMBER_CHANGE_UNIT",
  "rationale": "単位選択を、入力した数値の意味を明示的に変更する編集として統一する。現在の円半径や古いResultに依存せず、生成前や複数recordでも同じ操作を使える。",
  "mustPreserve": "numericdraftとunitの明示、Auto、Generateまで旧Resultを保つこと、History/Session、失敗復旧、既存scalarの科学的意味、local-only。",
  "mayRetire": "なし。",
  "acceptedResidualRisk": "倍率からpxへ切り替えると図上の大きさが変わる。helpに、数値維持・次回Generate反映を明記する。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-29"
}
```

### PD-OI-050: Circular Auto unit preference lifecycle

- Concern key: `tracks.circular-measure-auto-unit-lifecycle`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `A / TRANSIENT_AUTO_UNIT`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its nine-field JSON representation below.
- Decision source: the complete Choice A text reviewed in
  [03_AUTO_UNIT_LIFECYCLE.md](https://github.com/satoshikawato/gbdraw/blob/727b876214b58d4233790ce1feb4913eaf48bf4d/docs/internal/issue-619-implementation-plan-20260927/DECISION_PACKS/03_AUTO_UNIT_LIFECYCLE.md)
  at S00 commit `727b876214b58d4233790ce1feb4913eaf48bf4d`. On `2026-09-27`,
  `satoshikawato` explicitly confirmed signing all three Choice A texts with
  that Owner and Decision date. All nine supplied fields are reproduced
  without translation or additional terms. This record does not supersede
  another decision. Dependent runtime requires this authority merged into
  its base; this amendment supplies no runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `352ce20ffa4f81c9a7bc147d0275b8edad8b2a281311abf0dba01415217b40d6`.

```text
PRODUCT_DECISION
Concern: tracks.circular-measure-auto-unit-lifecycle
Scenario revision: 1
Choice: A / TRANSIENT_AUTO_UNIT
Rationale: Autoにgeometry上のunitはないため、その選択を次回入力用の小さなtransient preferenceとして扱う。manual値の意味と保存は保ち、追加の永続schemaやunit mirrorを避ける。
Must preserve: 空欄/Autoのnull意味、unitを先に選ぶ操作、manual値のunitとHistory/Session、既存preview/request、invaliddraftと失敗復旧。
May retire: なし。Autoのunit preferenceのHistory/Session保証は新設しない。
Accepted residual risk: 空欄時だけのunit選択はpanel再マウントやLoadで忘れられる。manual scalarのunitは必ず残り、Auto geometryは変わらない。
Owner: satoshikawato
Decision date: 2026-09-27
```

```json
{
  "concern": "tracks.circular-measure-auto-unit-lifecycle",
  "scenarioRevision": 1,
  "choice": "A / TRANSIENT_AUTO_UNIT",
  "rationale": "Autoにgeometry上のunitはないため、その選択を次回入力用の小さなtransient preferenceとして扱う。manual値の意味と保存は保ち、追加の永続schemaやunit mirrorを避ける。",
  "mustPreserve": "空欄/Autoのnull意味、unitを先に選ぶ操作、manual値のunitとHistory/Session、既存preview/request、invaliddraftと失敗復旧。",
  "mayRetire": "なし。Autoのunit preferenceのHistory/Session保証は新設しない。",
  "acceptedResidualRisk": "空欄時だけのunit選択はpanel再マウントやLoadで忘れられる。manual scalarのunitは必ず残り、Auto geometryは変わらない。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-27"
}
```

### PD-OI-051: In-flight comparison draft changes

- Concern key: `diagram-generation.inflight-comparison-draft`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `A / RUN_SNAPSHOT_COMMIT_DRAFT_RETAINED`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its nine-field JSON representation below.
- Decision source: the complete nine-field receipt and JSON below, reviewed,
  approved, and signed by `satoshikawato` on `2026-09-27` with
  「承認、署名します。devに統合してください。」 The owner selected Choice A
  and requested preparation of the remaining wording before that approval.
  All nine approved fields are reproduced without translation or additional
  terms. This record does not supersede another decision. Dependent runtime
  requires this authority merged into its base; this amendment supplies no
  runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding final newline):
  `b23439b7d9a744f1683735cc97ae2e86091febdf4eb3e11cfe967952f1c78776`.

```text
PRODUCT_DECISION
Concern: diagram-generation.inflight-comparison-draft
Scenario revision: 1
Choice: A / RUN_SNAPSHOT_COMMIT_DRAFT_RETAINED
Rationale: 開始時に固定したrequestの生成を完了させ、利用者が後から変更した比較draftは次のGenerateに残すことで、進行中の作業と次の操作を両立する。
Must preserve: 開始時のmixed upload/LOSAT比較snapshotとstable source-bound query/subject endpoints・ordinalsを使い、biological inputs・rules・current artifactが不変でCancel/new Generateがなければ生成を完了する。後から選んだnone draftと保持された3 edge draftsを上書きせず、次のGenerateではnoneを適用する。source bytes/exact identity、mixed upload/LOSATの区別、content-addressed raw cacheと互換cache再利用、Last successful Result/request、draft、Cancel・true stale・new-run隔離、atomic artifact/History/rollback、SessionのResult/draft分離、current Result export、既存keyboard/focusを維持する。rule/catalog/Result/editor/source変更による旧candidate拒否を維持する。#598の4決定、Decision Packs 01–05、PD-OI-035 revision 3、EXCLUSIVE_DIRECTIONS_WITHOUT_MATCH、logical pre-Align reference center、canvas/keyboard/focus/Editor/exact identity/Select/Skipの独立要求を変更しない。
May retire: none
Accepted residual risk: none
Owner: satoshikawato
Decision date: 2026-09-27
```

```json
{
  "concern": "diagram-generation.inflight-comparison-draft",
  "scenarioRevision": 1,
  "choice": "A / RUN_SNAPSHOT_COMMIT_DRAFT_RETAINED",
  "rationale": "開始時に固定したrequestの生成を完了させ、利用者が後から変更した比較draftは次のGenerateに残すことで、進行中の作業と次の操作を両立する。",
  "mustPreserve": "開始時のmixed upload/LOSAT比較snapshotとstable source-bound query/subject endpoints・ordinalsを使い、biological inputs・rules・current artifactが不変でCancel/new Generateがなければ生成を完了する。後から選んだnone draftと保持された3 edge draftsを上書きせず、次のGenerateではnoneを適用する。source bytes/exact identity、mixed upload/LOSATの区別、content-addressed raw cacheと互換cache再利用、Last successful Result/request、draft、Cancel・true stale・new-run隔離、atomic artifact/History/rollback、SessionのResult/draft分離、current Result export、既存keyboard/focusを維持する。rule/catalog/Result/editor/source変更による旧candidate拒否を維持する。#598の4決定、Decision Packs 01–05、PD-OI-035 revision 3、EXCLUSIVE_DIRECTIONS_WITHOUT_MATCH、logical pre-Align reference center、canvas/keyboard/focus/Editor/exact identity/Select/Skipの独立要求を変更しない。",
  "mayRetire": "none",
  "acceptedResidualRisk": "none",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-27"
}
```

### PD-OI-052: Matched decoration deltas across regeneration

- Concern key: `web.composition-decoration-continuity`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `A / CARRY-MATCHED-DECORATION-DELTAS`
- Normative outcome: 同じ図の legend/title/Linear scale の delta を新 automatic 配置へ1回加算。mode/grouping、validated source/region、record identity で照合。prefix/配列順/DOM順を使わず、batchの各出力も別々に対応。非ゼロ target の消失・未知対応は候補を公開せず、旧 Result を保持して対象 Reset または設定修正を案内。zero/fresh は自動 Generate
- Discoverability/accessibility / immediate feedback: 通常は自動継承。対応不能の理由・対象・Reset/修正を keyboard/touch から到達可能に表示
- Canonical state update: Result SVG を差分の正本とし、transaction-local snapshot を candidate に適用。UI refs は同期値。latent/global map なし
- Undo/Redo: Generate は継承を含む1 replacement。保存 Result を復元し delta を再加算しない
- Session / regeneration: 保存 Result から次の Generate の差分を取得。通常 load は保存 Result を維持。未知対応は無言破棄しない
- Export/artifact: 適用済みの current Result を各形式へ出力。raw Python recipeだけで手動位置を再現する保証なし
- Validation/error: finite delta、一意 target、図の同一性を検査。source/region/mode/grouping/record集合変更や欠落/重複/未知が転用不能なら候補公開前のエラー
- Failure/recovery / next available action: render/transform/bind失敗、Cancel/staleは旧Result/request/History保持。対応不能は対象 Reset または設定修正→Generate
- Scientific-output: 装飾位置のみ。recordTranslations/active alignment の record delta を二重加算しない
- Cache/provenance: 既存 validated identity/digest 使用。raw search/cache key に delta を追加しない。candidate と保存Resultを一致
- Performance: 非ゼロ対象のみ。候補の既存parse/serialize共用、batch旧SVGは必要分のみ、zero fast path維持
- Compatibility: writer/readerを維持。照合できない保存図は明示回復
- Decision source: the complete Choice A outcome and nine-field receipt in
  [DECISION_01_COMPOSITION_CONTINUITY.md](https://github.com/satoshikawato/gbdraw/blob/bcc4e0aa5ffcf4bdf9a12952caa8bf5de2cfee08/docs/internal/issue-599-preview-layout-20260926/DECISION_01_COMPOSITION_CONTINUITY.md)
  at published planning commit `bcc4e0aa5ffcf4bdf9a12952caa8bf5de2cfee08`.
  `satoshikawato` explicitly approved all three independent recommended A
  outcomes with 「すべて推奨案で承認します。」 on `2026-09-26`.
  The receipt and JSON below reproduce this concern's supplied fields without
  translation or additional rationale, retirement, or risk terms. This record
  does not supersede another decision. Dependent runtime requires this
  authority merged into its base; this amendment contains no runtime or
  runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `c7668d708ac6cd90b4373df0a685e418f21be7e653e99cbe899b82fba25e5e7d`.
- Acceptance contracts: `OIC-005`, `OIC-006`, `OIC-013`, `OIC-014`. These existing obligations and the
  complete selected outcome are jointly required; their citation does not
  claim completed dependent-runtime checks.

```text
PRODUCT_DECISION
Concern: web.composition-decoration-continuity
Scenario revision: 1
Choice: A / CARRY-MATCHED-DECORATION-DELTAS
Rationale: 色やfontを直すたびに装飾の配置をやり直す負担をなくし、別の図に位置を誤転用しない。
Must preserve: legend/title/Linear scaleのdrag・適用・Reset・History・Session・current Result export、source/region/record同一性、既存record/alignmentの意味、zero fast path、失敗/Cancel/stale時の旧Resultとcommitted requestを保持する。未知/欠落からの無言削除、別sourceへの誤転用、record deltaの二重加算、新schema/Worker/全History cloneを認めない。通常Generate/committed-candidate/automatic reflowに同じ候補境界を使い、batch全出力を対応identityへだけ適用する。
May retire: 通常Generateがlegend/title/Linear scaleの非ゼロdeltaを無言で捨てる動作だけ。diagram全体、個別record、padding、legend順の新しい継承保証は含めない。
Accepted residual risk: 新automatic配置に同じdeltaを加えるので絶対位置は変わり、clipping/overlapが残りうる。自動clampせずpadding/Resetで調整する。対応不能は候補公開前に止まり、明示Reset/設定修正が必要。
Owner: satoshikawato
Decision date: 2026-09-26
```

```json
{
  "concern": "web.composition-decoration-continuity",
  "scenarioRevision": 1,
  "choice": "A / CARRY-MATCHED-DECORATION-DELTAS",
  "rationale": "色やfontを直すたびに装飾の配置をやり直す負担をなくし、別の図に位置を誤転用しない。",
  "mustPreserve": "legend/title/Linear scaleのdrag・適用・Reset・History・Session・current Result export、source/region/record同一性、既存record/alignmentの意味、zero fast path、失敗/Cancel/stale時の旧Resultとcommitted requestを保持する。未知/欠落からの無言削除、別sourceへの誤転用、record deltaの二重加算、新schema/Worker/全History cloneを認めない。通常Generate/committed-candidate/automatic reflowに同じ候補境界を使い、batch全出力を対応identityへだけ適用する。",
  "mayRetire": "通常Generateがlegend/title/Linear scaleの非ゼロdeltaを無言で捨てる動作だけ。diagram全体、個別record、padding、legend順の新しい継承保証は含めない。",
  "acceptedResidualRisk": "新automatic配置に同じdeltaを加えるので絶対位置は変わり、clipping/overlapが残りうる。自動clampせずpadding/Resetで調整する。対応不能は候補公開前に止まり、明示Reset/設定修正が必要。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-26"
}
```

### PD-OI-053: Explicit Layout edit with discoverable targets

- Concern key: `web.layout-edit-affordance`
- Scenario revision: `2`
- Status: `ACCEPTED`
- Selected outcome: `A / EXPLICIT-MODE-WITH-DISCOVERABLE-TARGETS`
- Supersedes: `PD-OI-053`, scenario revision `1`
  (`A / EXPLICIT-MODE-WITH-DISCOVERABLE-TARGETS`, the same outcome),
  keeping every sentence of its receipt and adding the pointer cursors of
  a Legend row: with Layout edit off the pointing hand, a button that opens
  its popup and is not a canvas pan start; with Layout edit on grab and
  grabbing. For every other supported target scenario revision `1` holds
  unchanged.
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its JSON representation below.
- Discoverability/accessibility / feedback: 常設説明、toggle aria-pressed/説明、focus/touch。大量の SVG tab stop は追加しない
- Canonical state update: 既存mode ref。hintは派生表示でSVG/History/canonicalを変更しない
- Undo/Redo: hover/hint は履歴なし。実 drag のみ既存1操作
- Session / regeneration: mode/Result復元後に表示をrebind。hintは保存しない。Generate継承はPack01が決める
- Export/artifact: cursor/outline/hint は Preview 専用。plain/interactive SVG、PNG/PDF、保存Resultへ入れない
- Validation/error: 既存 composition eligibility で対象限定。未対応targetを動かせると説明しない
- Failure/recovery / next action: Result/load/Historyの既存bind。hintが出せない場合も常設説明とtoggleを使用可能
- Scientific-output: 発見方法だけ。生物学的意味/comparison/alignment不変
- Cache/provenance: hintをrequest/cache keyに含めず、Preview transientをclean serializationで除去
- Performance: 既存bindで対象限定、hoverで全走査/Worker/History cloneなし
- Compatibility: Session/modeの意味、writer/readerを維持
- Decision source: the Owner replies of `2026-10-13` quoted verbatim in the
  Revision 59 entry above. The Owner approved the receipt text below as
  written. The receipt reproduces its fields without translation or
  additional terms. Dependent runtime requires this authority merged into its
  base, or merged together with it through the reviewed static Product
  Contract co-change route; this record supplies no runtime acceptance
  evidence. Scenario revision `1` (receipt SHA-256
  `0365e0057b61e6b606c6714dc17f0e5c89d755ce651c1a8eec66d0139e0c65ec`) is retained in
  Git history.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `8cf2674991d1054424e478f7e5b1d5c31c1eb19467bbe68eab88d6827e2b23a7`.
- Acceptance contracts: `OIC-006`, `OIC-013`, `OIC-014`. These existing obligations and the
  complete selected outcome are jointly required; their citation does not
  claim completed dependent-runtime checks.

```text
PRODUCT_DECISION
Concern: web.layout-edit-affordance
Scenario revision: 2
Supersedes: PD-OI-053, scenario revision 1 (web.layout-edit-affordance)
Choice: A / EXPLICIT-MODE-WITH-DISCOVERABLE-TARGETS
Rationale: canvas panと配置編集の区別を保ちながら、対象と有効化方法を初めての利用者にも示す。凡例の項目（色の四角と名前）は、Layout edit が OFF のとき押すと popup が開くので指差しの手のポインターが合い、ON のときは押しても popup は開かずドラッグで凡例が動くので「つかめる手／つかんだ手」のポインターが合う。
Must preserve: OFFの従来panと明示toggle、ONのtarget drag、supported targetへのhelp/hover説明とtoolbarの常設説明、keyboard focus/touchで同じ説明へ到達すること、feature/label/legend個別編集とShift/Ctrlの優先順位、record/alignmentの既存動作、実dragの1 History操作、Session/Exportを維持する。hintだけでcanonical値を変えず、成果物へhintを保存しない。modeを自動ONにしない。〔凡例の項目、Layout edit が OFF〕項目はボタンである。載せると指差しの手のポインターになり、点線の枠と「Turn on Layout edit to move this item」は今までどおり出る。クリック、または Enter か Space で項目の popup が開く。項目から始めたドラッグでは図全体は動かず、popup も開かない。〔凡例の項目、Layout edit が ON〕載せると「つかめる手」（grab）、押している間とドラッグ中は「つかんだ手」（grabbing）、離すと「つかめる手」に戻る。ドラッグで凡例が動く。クリック、Enter、Space では popup は開かない。
May retire: OFF/ONの意味を区別できないcursorと説明不足だけ。gesture/編集機能/保存意味は退役しない。ただし凡例の項目の上に限り、Layout edit が OFF のときの「?」のポインターと、そこから始めたドラッグで図全体が動くことは退役する。
Accepted residual risk: 移動前に有効化の1手順が残る。常設説明とkeyboard/touch検証で発見可能性を補う。Layout edit が OFF のとき、凡例の項目の上から始めたドラッグでは図全体を動かせない（項目のない場所からは動かせる）。
Owner: satoshikawato
Decision date: 2026-10-13
```

```json
{
  "concern": "web.layout-edit-affordance",
  "scenarioRevision": 2,
  "supersedes": "PD-OI-053, scenario revision 1 (web.layout-edit-affordance)",
  "choice": "A / EXPLICIT-MODE-WITH-DISCOVERABLE-TARGETS",
  "rationale": "canvas panと配置編集の区別を保ちながら、対象と有効化方法を初めての利用者にも示す。凡例の項目（色の四角と名前）は、Layout edit が OFF のとき押すと popup が開くので指差しの手のポインターが合い、ON のときは押しても popup は開かずドラッグで凡例が動くので「つかめる手／つかんだ手」のポインターが合う。",
  "mustPreserve": "OFFの従来panと明示toggle、ONのtarget drag、supported targetへのhelp/hover説明とtoolbarの常設説明、keyboard focus/touchで同じ説明へ到達すること、feature/label/legend個別編集とShift/Ctrlの優先順位、record/alignmentの既存動作、実dragの1 History操作、Session/Exportを維持する。hintだけでcanonical値を変えず、成果物へhintを保存しない。modeを自動ONにしない。〔凡例の項目、Layout edit が OFF〕項目はボタンである。載せると指差しの手のポインターになり、点線の枠と「Turn on Layout edit to move this item」は今までどおり出る。クリック、または Enter か Space で項目の popup が開く。項目から始めたドラッグでは図全体は動かず、popup も開かない。〔凡例の項目、Layout edit が ON〕載せると「つかめる手」（grab）、押している間とドラッグ中は「つかんだ手」（grabbing）、離すと「つかめる手」に戻る。ドラッグで凡例が動く。クリック、Enter、Space では popup は開かない。",
  "mayRetire": "OFF/ONの意味を区別できないcursorと説明不足だけ。gesture/編集機能/保存意味は退役しない。ただし凡例の項目の上に限り、Layout edit が OFF のときの「?」のポインターと、そこから始めたドラッグで図全体が動くことは退役する。",
  "acceptedResidualRisk": "移動前に有効化の1手順が残る。常設説明とkeyboard/touch検証で発見可能性を補う。Layout edit が OFF のとき、凡例の項目の上から始めたドラッグでは図全体を動かせない（項目のない場所からは動かせる）。",
  "owner": "satoshikawato",
  "decisionDate": "2026-10-13"
}
```

### PD-OI-054: Docked Preview search and controls

- Concern key: `web.preview-search-placement`
- Scenario revision: `2`
- Supersedes: `PD-OI-054`, scenario revision `1` (`A / DOCKED-SEARCH-AND-CONTROLS`).
- Status: `ACCEPTED`
- Selected outcome: `A / DOCKED-SEARCH-WITH-TOP-EDITOR`
- Normative outcome: exactly the receipt below.
- Discoverability/accessibility / feedback: 安定した順序、keyboard、wrap/scroll、focus維持。short高さでoverflow clipによる隠れなし
- Canonical state update: query等は既存search owner。geometryはCSSのみ、新座標ref/observerなし
- Undo/Redo: chromeはartifact Historyに入れず、図の履歴を維持
- Session / regeneration: chrome位置は新たに保存しない。既存Session/Result/draftを維持
- Export/artifact: HTML chromeは図/SVG/PNG/PDF/保存Resultへ入らない
- Validation/error: 通常幅でsearch/toolbarとEditorの重なりなし、EditorはPreview上端から開始、検索は利用可能幅内。短い高さで到達性検査
- Failure/recovery / next action: drawer/resizeでquery/active/focus保持。Result消失は既存visibilityで閉じる
- Scientific-output: chromeのみ、生物学的値/比較/scale不変
- Cache/provenance: chromeをWorker/cache/requestに入れない
- Performance: CSSのみ。位置computed/global drag listenerなし
- Compatibility: 保存形式維持、自由drag退役を維持
- Decision source: the receipt text approved by `satoshikawato` on
  `2026-09-29` (GUI remediation S00 decision 4) supplies the Rationale, the
  Accepted residual risk, and the additional scope (search and toolbar in the
  width left by the Editor, which covers neither). The Owner's confirmed
  requirement supplies the rest: search at most 39.5 rem within the available
  width and the Editor from the Preview top edge. The remaining fields repeat
  scenario revision `1`. This serialization adds no other terms.
- Acceptance contracts: `OIC-006`, `OIC-013`, `OIC-014`, `OIC-025`, `OIC-026`.

```json
{
  "concern": "web.preview-search-placement",
  "scenarioRevision": 2,
  "choice": "A / DOCKED-SEARCH-WITH-TOP-EDITOR",
  "rationale": "派生 status と常時説明が操作応答を損ない（Result 後の比較切替 約 1.3 s）、画面を圧迫するため削除・help-tip 化する。",
  "mustPreserve": "searchの全field/query/regex/Prev/Next/Open/Enter、active match/focus、全toolbar操作、drawer tab/Close/Escape、同じsearch/canvas/SVG/editor DOM、既存Session/Export/Historyを維持する。通常幅では検索（最大39.5rem、利用可能幅内）とtoolbarをEditor幅を除いた残り幅に置き、EditorはPreview上端から開いて両方を覆わない。drawer幅は同一CSS変数を共有する。通常高さ740px以上の受入条件でworkspace200px以上、short viewport/keyboard/200%zoomで全操作へscroll到達可能にする。狭幅のEditor/reviewのcanvas確保、scroll、Close、keyboard操作とalignment reviewの仕様は変更しない。",
  "mayRetire": "通常幅での全幅の検索専用rowと、EditorがPreviewの検索rowより下から始まる配置だけ。検索の自由drag、固定360pxのJS退避、新座標ref/observerは導入しない。",
  "acceptedResidualRisk": "help-tip を開かない利用者は Generate/Save/Lock の事前説明を見ない。実 error・Processing/Canceling・recovery は保持する。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-29"
}
```

### PD-OI-055: Valid bindings after failed Generate

- Concern key: `web.generate.failed-source-binding-continuation`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `A / RETAIN_VALIDATED_BINDING_ENRICHMENT`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its nine-field JSON representation below.
- Decision source: the complete Choice A text presented to the Product Decision
  Owner for Issue `#619` finding 5, explicitly approved in full, including
  Owner and Decision date, by `satoshikawato` on `2026-09-28`. The receipt
  and JSON below reproduce the approved fields without translation or
  additional terms. This record does not supersede another decision.
  Dependent runtime requires this authority merged into its base; this
  amendment supplies no runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline): `472b82ef4731c8656827a1c0903e80cc3b46d60e0b424ea0568de9f95fdb7f7b`.
- Acceptance contracts: `OIC-005`, `OIC-013`, `OIC-014`, and `OIC-015` remain
  jointly required with the selected outcome; citation does not claim the
  remaining multi-record, cancel/stale, Export, or retry evidence is complete.

```text
PRODUCT_DECISION
Concern: web.generate.failed-source-binding-continuation
Scenario revision: 1
Choice: A / RETAIN_VALIDATED_BINDING_ENRICHMENT
Rationale: Preserve the continuation observed in S01. Valid source bindings resolved during a failed Generate may remain in the editable document and a subsequently saved Session, although no new Result was admitted.
Must preserve: The previous Result and canonical request, editable drafts, History, exact source bytes, record and annotation identity, coherent Save/Load, Export of the previous Result, retry, and cancel, stale, and superseded recovery. Only complete, validated bindings to the same source may persist. Independent source discovery completed before Generate remains valid.
May retire: None of the existing supported behavior. A strict guarantee that every saved binding field remains unchanged after a failed Generate is not adopted.
Accepted residual risk: A failed Generate may change binding metadata in a later saved Session while the displayed Result remains unchanged; this may surprise someone comparing Session files. Incorrect, incomplete, stale, dangling, or wrong-source bindings are not accepted.
Owner: satoshikawato
Decision date: 2026-09-28
```

```json
{
  "concern": "web.generate.failed-source-binding-continuation",
  "scenarioRevision": 1,
  "choice": "A / RETAIN_VALIDATED_BINDING_ENRICHMENT",
  "rationale": "Preserve the continuation observed in S01. Valid source bindings resolved during a failed Generate may remain in the editable document and a subsequently saved Session, although no new Result was admitted.",
  "mustPreserve": "The previous Result and canonical request, editable drafts, History, exact source bytes, record and annotation identity, coherent Save/Load, Export of the previous Result, retry, and cancel, stale, and superseded recovery. Only complete, validated bindings to the same source may persist. Independent source discovery completed before Generate remains valid.",
  "mayRetire": "None of the existing supported behavior. A strict guarantee that every saved binding field remains unchanged after a failed Generate is not adopted.",
  "acceptedResidualRisk": "A failed Generate may change binding metadata in a later saved Session while the displayed Result remains unchanged; this may surprise someone comparing Session files. Incorrect, incomplete, stale, dangling, or wrong-source bindings are not accepted.",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-28"
}
```

### PD-OI-056: Feature search All-field scope

- Concern key: `web.feature-search.all-field-scope`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `A / ALL-EXCLUDES-SEQUENCE-CONTENT`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its nine-field JSON representation below.
- Decision source: the complete `D-01` receipt in [`02_DECISION_PACK.md`](./web-gui-audit-20260930/02_DECISION_PACK.md) at P00 merge commit
  `e97d90fecfb327135eb50e85cfa9a87be3145823`, approved by `satoshikawato` on `2026-09-30` through the two
  Owner replies quoted verbatim in the Revision 29 entry above. The receipt and
  JSON below reproduce all nine supplied fields without translation or
  additional terms. This record does not supersede another decision.
  Dependent runtime requires this authority merged into its base; this
  amendment supplies no runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `6d866b12aea2e0e5de149f1fdad7802d42e6ac68e91983c1e7bb98b530094ec4`.
- Acceptance contracts: `OIC-006`. These obligations and the complete
  selected outcome are jointly required; their citation does not claim
  completed dependent-runtime checks.

```text
PRODUCT_DECISION
Concern: web.feature-search.all-field-scope
Scenario revision: 1
Choice: A / ALL-EXCLUDES-SEQUENCE-CONTENT
Rationale: 既定の All で遺伝子名を検索したとき、名前の一致だけが返るようにする。配列への偶然の一致で結果が埋まらないようにする。
Must preserve: 専用の Nucleotide sequence と Amino acid sequence の field による配列検索（IUPAC の展開を含む）、Label・qualifier・Location など他の field、編集後のラベルの検索、Preview と Interactive SVG の検索結果の一致、件数の表示。
May retire: All が配列の内容と /translation の値に一致する動作。
Accepted residual risk: 配列の断片を All で検索していた利用者は、専用の field を選ぶ必要がある。作り直す前の Interactive SVG は旧挙動のまま。
Owner: satoshikawato
Decision date: 2026-09-30
```

```json
{
  "concern": "web.feature-search.all-field-scope",
  "scenarioRevision": 1,
  "choice": "A / ALL-EXCLUDES-SEQUENCE-CONTENT",
  "rationale": "既定の All で遺伝子名を検索したとき、名前の一致だけが返るようにする。配列への偶然の一致で結果が埋まらないようにする。",
  "mustPreserve": "専用の Nucleotide sequence と Amino acid sequence の field による配列検索（IUPAC の展開を含む）、Label・qualifier・Location など他の field、編集後のラベルの検索、Preview と Interactive SVG の検索結果の一致、件数の表示。",
  "mayRetire": "All が配列の内容と /translation の値に一致する動作。",
  "acceptedResidualRisk": "配列の断片を All で検索していた利用者は、専用の field を選ぶ必要がある。作り直す前の Interactive SVG は旧挙動のまま。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-30"
}
```

### PD-OI-057: Keyboard- and touch-reachable help tips

- Concern key: `web.help-tip.reachability`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `A / FOCUSABLE-DISCLOSURE-TIPS`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its nine-field JSON representation below.
- Decision source: the complete `D-02` receipt in [`02_DECISION_PACK.md`](./web-gui-audit-20260930/02_DECISION_PACK.md) at P00 merge commit
  `e97d90fecfb327135eb50e85cfa9a87be3145823`, approved by `satoshikawato` on `2026-09-30` through the two
  Owner replies quoted verbatim in the Revision 29 entry above. The receipt and
  JSON below reproduce all nine supplied fields without translation or
  additional terms. This record does not supersede another decision.
  Dependent runtime requires this authority merged into its base; this
  amendment supplies no runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `d61be689354b873b898680fae92e59e64672d02fbe9d1540cb70c5d8faa7b5c9`.
- Acceptance contracts: `OIC-006`. These obligations and the complete
  selected outcome are jointly required; their citation does not claim
  completed dependent-runtime checks.

```text
PRODUCT_DECISION
Concern: web.help-tip.reachability
Scenario revision: 1
Choice: A / FOCUSABLE-DISCLOSURE-TIPS
Rationale: キーボードとタッチの利用者が、hover と同じ説明に届くようにする。常時表示の説明を help tip に移した（PD-OI-054）ので、tip に届くことが必要になった。
Must preserve: 各 tip の文言、hover での表示、周囲の label と control の accessible name、Escape で閉じられること、390 px での操作、既存の id 付き tip の挙動。
May retire: id のない tip を hover 専用の aria-hidden の icon にする設計（js/components.js の意図のコメント）。
Accepted residual risk: tab stop が最大 175 個増え、キーボードでの移動が長くなる。Gallery と docs の capture を撮り直すことがある。
Owner: satoshikawato
Decision date: 2026-09-30
```

```json
{
  "concern": "web.help-tip.reachability",
  "scenarioRevision": 1,
  "choice": "A / FOCUSABLE-DISCLOSURE-TIPS",
  "rationale": "キーボードとタッチの利用者が、hover と同じ説明に届くようにする。常時表示の説明を help tip に移した（PD-OI-054）ので、tip に届くことが必要になった。",
  "mustPreserve": "各 tip の文言、hover での表示、周囲の label と control の accessible name、Escape で閉じられること、390 px での操作、既存の id 付き tip の挙動。",
  "mayRetire": "id のない tip を hover 専用の aria-hidden の icon にする設計（js/components.js の意図のコメント）。",
  "acceptedResidualRisk": "tab stop が最大 175 個増え、キーボードでの移動が長くなる。Gallery と docs の capture を撮り直すことがある。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-30"
}
```

### PD-OI-058: Managed Depth track-row lifecycle

- Concern key: `web.depth.managed-slot-lifecycle`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `A / ADD-ON-FIRST-SOURCE-BOTH-MODES`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its nine-field JSON representation below.
- Decision source: the complete `D-03` receipt in [`02_DECISION_PACK.md`](./web-gui-audit-20260930/02_DECISION_PACK.md) at P00 merge commit
  `e97d90fecfb327135eb50e85cfa9a87be3145823`, approved by `satoshikawato` on `2026-09-30` through the two
  Owner replies quoted verbatim in the Revision 29 entry above. The receipt and
  JSON below reproduce all nine supplied fields without translation or
  additional terms. This record does not supersede another decision.
  Dependent runtime requires this authority merged into its base; this
  amendment supplies no runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `3dbd45fd6b6f047f6c5cdda98f1723330bf22be0802f49d8051726260f03609b`.
- Acceptance contracts: `OIC-006`, `OIC-020`. These obligations and the complete
  selected outcome are jointly required; their citation does not claim
  completed dependent-runtime checks.

```text
PRODUCT_DECISION
Concern: web.depth.managed-slot-lifecycle
Scenario revision: 1
Choice: A / ADD-ON-FIRST-SOURCE-BOTH-MODES
Rationale: Circular と Linear で、Depth ファイルを設定したときの track 行の振る舞いを揃える。利用者が消した行や無効にした行を勝手に戻さない。
Must preserve: 明示の track slot が有効なときの authority、利用者が削除・無効化・移動した行、行の params と legend_label、Reset による再生成、Undo/Redo、Session の往復、PD-OI-025 の論理 series の範囲。論理 series が最初の source を得たとき、その index を参照する行（有効・無効を問わない）がなければ managed 行を 1 つ足し、series が source を失ったら managed 行を除く。
May retire: 無関係な切り替えのたびに Circular の watcher が depth 行を作り直す・再び有効にする・付け替える動作。Linear の "Add Depth TSV series" が無効な行を再び有効にする動作。
Accepted residual risk: 無効な stack に depth 行を持たない旧 Session は、読み込んで stack を有効にしても行が自動では足されない（Reset で作れる）。
Owner: satoshikawato
Decision date: 2026-09-30
```

```json
{
  "concern": "web.depth.managed-slot-lifecycle",
  "scenarioRevision": 1,
  "choice": "A / ADD-ON-FIRST-SOURCE-BOTH-MODES",
  "rationale": "Circular と Linear で、Depth ファイルを設定したときの track 行の振る舞いを揃える。利用者が消した行や無効にした行を勝手に戻さない。",
  "mustPreserve": "明示の track slot が有効なときの authority、利用者が削除・無効化・移動した行、行の params と legend_label、Reset による再生成、Undo/Redo、Session の往復、PD-OI-025 の論理 series の範囲。論理 series が最初の source を得たとき、その index を参照する行（有効・無効を問わない）がなければ managed 行を 1 つ足し、series が source を失ったら managed 行を除く。",
  "mayRetire": "無関係な切り替えのたびに Circular の watcher が depth 行を作り直す・再び有効にする・付け替える動作。Linear の \"Add Depth TSV series\" が無効な行を再び有効にする動作。",
  "acceptedResidualRisk": "無効な stack に depth 行を持たない旧 Session は、読み込んで stack を有効にしても行が自動では足されない（Reset で作れる）。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-30"
}
```

### PD-OI-059: When the definition settings apply

- Concern key: `web.circular.definition-application`
- Scenario revision: `3`
- Supersedes: `PD-OI-059`, scenario revision `2` (`A / TITLES-LIVE-BOTH-MODES`).
  Scenario revision `3` changes only what its receipt states; in the decision
  packet's words, Species, Strain, Default font size, the definition rows,
  Keep Full Definition and the Circular record label go live "when the S5-Q1
  outcome makes their unit live; until then on Generate with the note".
  Scenario revision `2` superseded `PD-OI-059`, scenario revision `1`
  (`A / CIRCULAR-DEFINITION-APPLIES-ON-GENERATE`), in both diagram modes:
  Plot title, Title position, and Title font apply live (`PD-OI-102`,
  `PD-OI-103`), and Species, Strain, Default font size, and Keep Full
  Definition stay on Generate with an Applies on Generate note, as the
  receipt below states.
- Status: `ACCEPTED`
- Selected outcomes: `A / TITLES-LIVE-BOTH-MODES` (scenario revision `2`)
  and `A / DEFINITION-SETTINGS-LIVE-WITH-S5-Q1` (scenario revision `3`)
- Normative outcome: exactly the complete approved scenario revision `3`
  `PRODUCT_DECISION` receipt and its nine-field JSON representation below,
  with the effects of the scenario revision `2` receipt that it does not
  change. The scenario revision `2` receipt is kept as it was approved.
- Decision source (scenario revision `2`): the Owner replies quoted verbatim
  in the Revision 57 entry above. The receipt and JSON below reproduce all
  nine confirmed fields without translation or additional terms.
- Receipt SHA-256 (scenario revision `2`; UTF-8, excluding the final newline):
  `eaf258733eaaa93868f28d06b55f8fa7bf8d4ff7cd6afda8f38ac5b0c96c83a3`.
- Decision source (scenario revision `3`): the Owner replies quoted
  verbatim in the Revision 64 entry above. The scenario revision `3`
  receipt and JSON below reproduce all nine fields without translation or
  additional terms. Dependent runtime requires this authority merged into its
  base; this record supplies no runtime acceptance evidence.
- Receipt SHA-256 (scenario revision `3`; UTF-8, excluding the final newline):
  `9d9f749fff0f9a8f24a74d912b707a21712649d674c0a06f183cb0ef90ffc4cf`.
- Acceptance contracts: `OIC-009`, `OIC-010`, `OIC-024`, `OIC-027`. These obligations and the complete
  selected outcome are jointly required; their citation does not claim
  completed dependent-runtime checks.

Scenario revision `2` receipt, confirmed on `2026-10-17`:

```text
PRODUCT_DECISION
Concern: web.circular.definition-application
Scenario revision: 2
Choice: A / TITLES-LIVE-BOTH-MODES
Rationale: Plot Title is in both modes and in G1; a full rerender keeps every definition line consistent with crop length, GC% and record label without the retired regenerate_definition_svgs helper.
Must preserve: editing, saving and History of Species, Strain, Plot title, Title position, Title font, Default font size, Keep Full Definition; the correct definition after Generate (region length and GC%, record label and subtitle, reverse complement, grid order); the other live edits; Species, Strain, Default font size and Keep Full Definition stay on Generate with an Applies on Generate note.
May retire: the Generate-only timing of Plot Title, Title position and Title font in both modes, and revision 1's reason "align with Linear's on-Generate timing".
Accepted residual risk: as S4-Q2 A for Plot Title (+54-85 ms over 250 ms on large Sessions).
Owner: Satoshi Kawato
Decision date: 2026-10-17
```

```json
{
  "concern": "web.circular.definition-application",
  "scenarioRevision": 2,
  "choice": "A / TITLES-LIVE-BOTH-MODES",
  "rationale": "Plot Title is in both modes and in G1; a full rerender keeps every definition line consistent with crop length, GC% and record label without the retired regenerate_definition_svgs helper.",
  "mustPreserve": "editing, saving and History of Species, Strain, Plot title, Title position, Title font, Default font size, Keep Full Definition; the correct definition after Generate (region length and GC%, record label and subtitle, reverse complement, grid order); the other live edits; Species, Strain, Default font size and Keep Full Definition stay on Generate with an Applies on Generate note.",
  "mayRetire": "the Generate-only timing of Plot Title, Title position and Title font in both modes, and revision 1's reason \"align with Linear's on-Generate timing\".",
  "acceptedResidualRisk": "as S4-Q2 A for Plot Title (+54-85 ms over 250 ms on large Sessions).",
  "owner": "Satoshi Kawato",
  "decisionDate": "2026-10-17"
}
```

Scenario revision `3` receipt, confirmed on `2026-10-18`:

```text
PRODUCT_DECISION
Concern: web.circular.definition-application
Scenario revision: 3
Choice: A / DEFINITION-SETTINGS-LIVE-WITH-S5-Q1
Rationale: they follow S5-Q1 without another round.
Must preserve: the correct definition after each apply and Generate (region length and GC%, record label and subtitle, reverse complement, grid order); their editing, saving, History; rev 2's G1 titles.
May retire: rev 2's "stay on Generate with an Applies on Generate note", per live unit.
Accepted residual risk: as S5-Q1 A for the unit.
Owner: satoshikawato
Decision date: 2026-10-18
```

```json
{
  "concern": "web.circular.definition-application",
  "scenarioRevision": 3,
  "choice": "A / DEFINITION-SETTINGS-LIVE-WITH-S5-Q1",
  "rationale": "they follow S5-Q1 without another round.",
  "mustPreserve": "the correct definition after each apply and Generate (region length and GC%, record label and subtitle, reverse complement, grid order); their editing, saving, History; rev 2's G1 titles.",
  "mayRetire": "rev 2's \"stay on Generate with an Applies on Generate note\", per live unit.",
  "acceptedResidualRisk": "as S5-Q1 A for the unit.",
  "owner": "satoshikawato",
  "decisionDate": "2026-10-18"
}
```

### PD-OI-060: When the diagram-wide stroke settings apply

- Concern key: `web.stroke.application`
- Scenario revision: `2`
- Supersedes: `PD-OI-060`, scenario revision `1`
  (`C / STROKE-APPLIES-ON-GENERATE`). Scenario revision `2` changes only what
  its receipt states; in the decision packet's words, the global block and
  line strokes (Feature stroke) and the axis and scale strokes (Axis) go live
  "when the S5-Q1 outcome makes their unit live; until then on Generate with
  the note".
- Status: `ACCEPTED`
- Selected outcomes: `C / STROKE-APPLIES-ON-GENERATE` (scenario revision `1`)
  and `A / GLOBAL-STROKES-LIVE-WITH-S5-Q1` (scenario revision `2`)
- Normative outcome: exactly the complete approved scenario revision `2`
  `PRODUCT_DECISION` receipt and its nine-field JSON representation below,
  with the effects of the scenario revision `1` receipt that it does not
  change. The scenario revision `1` receipt is kept as it was approved.
- Decision source (scenario revision `1`): the complete `D-05` receipt in [`02_DECISION_PACK.md`](./web-gui-audit-20260930/02_DECISION_PACK.md) at P00 merge commit
  `e97d90fecfb327135eb50e85cfa9a87be3145823`, approved by `satoshikawato` on `2026-09-30` through the two
  Owner replies quoted verbatim in the Revision 29 entry above. The receipt and
  JSON below reproduce all nine supplied fields without translation or
  additional terms.
- Receipt SHA-256 (scenario revision `1`; UTF-8, excluding the final newline):
  `f95d7c6e1b7354f805d186d2709637b246ef2f15d416bc1143e1c46816ef14af`.
- Decision source (scenario revision `2`): the Owner replies quoted
  verbatim in the Revision 64 entry above. The scenario revision `2`
  receipt and JSON below reproduce all nine fields without translation or
  additional terms. Dependent runtime requires this authority merged into its
  base; this record supplies no runtime acceptance evidence.
- Receipt SHA-256 (scenario revision `2`; UTF-8, excluding the final newline):
  `c11a661c2922c3bda0a37c77a94dc0015dc4a9323bad86ef78f4c9aefb20ae0a`.
- Acceptance contracts: `OIC-024`, `OIC-027`. These obligations and the complete
  selected outcome are jointly required; their citation does not claim
  completed dependent-runtime checks.

Scenario revision `1` receipt, approved by `satoshikawato` on `2026-09-30`:

```text
PRODUCT_DECISION
Concern: web.stroke.application
Scenario revision: 1
Choice: C / STROKE-APPLIES-ON-GENERATE
Rationale: 空欄・不正な値・Auto への戻しが、Generate と違う stroke として Result に残らないようにする。個々の feature の stroke 指定を全体の設定で上書きしないようにする。
Must preserve: 全体の stroke 設定（block、line、axis、scale の幅と色）の編集・保存・Generate での適用、個々の feature の stroke 編集の即時反映と Auto への復元、不正な値の検証。各設定には Applies on Generate の表示を付ける。
May retire: 全体の stroke 設定の即時反映。
Accepted residual risk: 全体の stroke を変えたとき、Generate するまで preview が変わらない。
Owner: satoshikawato
Decision date: 2026-09-30
```

```json
{
  "concern": "web.stroke.application",
  "scenarioRevision": 1,
  "choice": "C / STROKE-APPLIES-ON-GENERATE",
  "rationale": "空欄・不正な値・Auto への戻しが、Generate と違う stroke として Result に残らないようにする。個々の feature の stroke 指定を全体の設定で上書きしないようにする。",
  "mustPreserve": "全体の stroke 設定（block、line、axis、scale の幅と色）の編集・保存・Generate での適用、個々の feature の stroke 編集の即時反映と Auto への復元、不正な値の検証。各設定には Applies on Generate の表示を付ける。",
  "mayRetire": "全体の stroke 設定の即時反映。",
  "acceptedResidualRisk": "全体の stroke を変えたとき、Generate するまで preview が変わらない。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-30"
}
```

Scenario revision `2` receipt, confirmed on `2026-10-18`:

```text
PRODUCT_DECISION
Concern: web.stroke.application
Scenario revision: 2
Choice: A / GLOBAL-STROKES-LIVE-WITH-S5-Q1
Rationale: they follow S5-Q1 without another round.
Must preserve: rev 1's reasons (a blank, invalid or Auto value fails the apply or draws as Generate does, through one projection); per-feature stroke edits stay instant with their Auto restore and are not overridden (editor intent applies after the rerender).
May retire: rev 1's Generate-only timing per live unit.
Accepted residual risk: as S5-Q1 A for the unit.
Owner: satoshikawato
Decision date: 2026-10-18
```

```json
{
  "concern": "web.stroke.application",
  "scenarioRevision": 2,
  "choice": "A / GLOBAL-STROKES-LIVE-WITH-S5-Q1",
  "rationale": "they follow S5-Q1 without another round.",
  "mustPreserve": "rev 1's reasons (a blank, invalid or Auto value fails the apply or draws as Generate does, through one projection); per-feature stroke edits stay instant with their Auto restore and are not overridden (editor intent applies after the rerender).",
  "mayRetire": "rev 1's Generate-only timing per live unit.",
  "acceptedResidualRisk": "as S5-Q1 A for the unit.",
  "owner": "satoshikawato",
  "decisionDate": "2026-10-18"
}
```

### PD-OI-061: Legend rename collisions for every entry

- Concern key: `web.legend.rename-collision`
- Scenario revision: `5`
- Supersedes: `PD-OI-061`, scenario revision `4`. Scenario revision `5`
  keeps the outcome `B / MERGE-SAME-FEATURE-TYPE-ONLY`, the choices, and the
  conditions of **Merge**, and retires the naming in which the original row
  keeps its name without a number and only the renamed row gets a numbered
  name. Scenario revision `4`
  keeps the outcome `B / MERGE-SAME-FEATURE-TYPE-ONLY` and the choice names,
  retires the scenario-2 preservation term
  「衝突先が色ルールの caption のときの PD-OI-042 の区別」, and adds the outcome of a
  pick from the list under the feature popup's **Legend name**. Scenario
  revision `3` superseded scenario revision `2` only by renaming the two
  choices of the Legend Name Conflict dialog. Scenario revision `2`
  (`B / MERGE-SAME-FEATURE-TYPE-ONLY`) superseded scenario revision `1`
  (`A / DIALOG-FOR-ALL-ENTRIES`).
- Status: `ACCEPTED`
- Selected outcomes: `B / MERGE-SAME-FEATURE-TYPE-ONLY` (scenario revision
  `2`, and scenario revision `4` with the names of color rules and the
  popup's list), the label change of scenario revision `3`, and the naming of
  **Keep both** of scenario revision `5`; the receipts of scenario revisions
  `3` and `5` name no choice code
- Normative outcome: exactly the four approved `PRODUCT_DECISION` receipts
  and their nine-field JSON representations below, read together, and the
  Owner's Replies 4, 5, and 6 of the Revision 49 entry with their questions.
  Where the scenario-5 receipt differs from an earlier receipt, the
  scenario-5 receipt governs: **Keep both** keeps the typed name, and every
  row of that name is numbered in legend order. Where the scenario-4 receipt
  differs from an earlier receipt, the scenario-4 receipt
  governs: a rename onto the name of a color rule asks the same dialog as a
  rename onto any other name, and the scenario-2 preservation of the
  `PD-OI-042` naming for such a target is retired. Where the scenario-2
  receipt names **Suffix** and **Merge**, it means the choices now labelled
  `Keep both: name it "<name 2>"` and `Merge into "<name>"`.
- Decision source (scenario revision `2`): the Owner replies of `2026-10-06`
  quoted verbatim in the Revision 32 entry above. The Owner approved the
  receipt text below as written ("OKです"). Dependent runtime (OV-62) is
  changed in the same pull request, by the static Product Contract co-change
  route; this record supplies no runtime acceptance evidence.
- Receipt SHA-256 (scenario revision `2`; UTF-8, excluding the final newline):
  `783436b02c70f043f89dc07b4bf8d56e0869e996507f72fedbbecd0301c8f322`.
- Decision source (scenario revision `3`): the Owner replies of `2026-10-11`
  quoted verbatim in the Revision 40 entry above. The second receipt and JSON
  below reproduce the confirmed fields without translation or additional
  terms. The runtime that shows the new names merges with it (static Product
  Contract co-change).
- Receipt SHA-256 (scenario revision `3`; UTF-8, excluding the final newline):
  `eb19efb8436cfee30d6955f19fc9b50433c24b5112ef20a7088d4a988a3f5da2`.
- Decision source (scenario revision `4`): the Owner replies of `2026-10-13`
  quoted verbatim in the Revision 45 entry above. The third receipt and JSON
  below reproduce all nine supplied fields without translation or additional
  terms. The runtime and tests that implement it merge with it (static
  Product Contract co-change); this record supplies no runtime acceptance
  evidence.
- Receipt SHA-256 (scenario revision `4`; UTF-8, excluding the final newline):
  `3f96ff5e0e94f2689233e2bde57b17ff945130e3ad92f04070eced277d80bf74`.
- Decision source (scenario revision `5`): the Owner replies of `2026-10-13`
  quoted verbatim in the Revision 49 entry above. The fourth receipt and JSON
  below reproduce all nine supplied fields without translation or additional
  terms. The runtime, tests, and documentation that implement it merge with
  it (static Product Contract co-change); this record supplies no runtime
  acceptance evidence.
- Receipt SHA-256 (scenario revision `5`; UTF-8, excluding the final newline):
  `fe82de9744104885fea0d9f44209565f879373c6b9462ae6f44df99f06247f63`.
- Acceptance contracts: `OIC-006`. These obligations and the complete
  selected outcome are jointly required; their citation does not claim
  completed dependent-runtime checks.

```text
PRODUCT_DECISION
Concern: web.legend.rename-collision
Scenario revision: 2
Choice: B / MERGE-SAME-FEATURE-TYPE-ONLY
Rationale: 凡例の名前を既存の名前に変えたとき、Merge は両方の項目が同じ track の同じ feature type を表すときだけ示す。それ以外は Suffix と Cancel を示す。CDS と GC content のように別のものを 1 つの凡例項目にまとめない。
Must preserve: 衝突しない rename の即時反映、すべての項目で衝突時に選択のダイアログを示すこと（原因の分からないエラーで止めない）、衝突先が色ルールの caption のときの PD-OI-042 の区別、Undo/Redo、Generate と Session での保持。
May retire: feature のない項目や別の feature type の項目を、Merge で別の項目にまとめる動作。
Accepted residual risk: 別の feature type の項目に同じ名前を付けたいときは、Suffix の付いた名前（例: "CDS (1)"）になる。描画が名付けた項目（例: "other proteins"）は、その feature type が分かるまで Suffix と Cancel だけを示す。
Owner: satoshikawato
Decision date: 2026-10-06
```

```json
{
  "concern": "web.legend.rename-collision",
  "scenarioRevision": 2,
  "choice": "B / MERGE-SAME-FEATURE-TYPE-ONLY",
  "rationale": "凡例の名前を既存の名前に変えたとき、Merge は両方の項目が同じ track の同じ feature type を表すときだけ示す。それ以外は Suffix と Cancel を示す。CDS と GC content のように別のものを 1 つの凡例項目にまとめない。",
  "mustPreserve": "衝突しない rename の即時反映、すべての項目で衝突時に選択のダイアログを示すこと（原因の分からないエラーで止めない）、衝突先が色ルールの caption のときの PD-OI-042 の区別、Undo/Redo、Generate と Session での保持。",
  "mayRetire": "feature のない項目や別の feature type の項目を、Merge で別の項目にまとめる動作。",
  "acceptedResidualRisk": "別の feature type の項目に同じ名前を付けたいときは、Suffix の付いた名前（例: \"CDS (1)\"）になる。描画が名付けた項目（例: \"other proteins\"）は、その feature type が分かるまで Suffix と Cancel だけを示す。",
  "owner": "satoshikawato",
  "decisionDate": "2026-10-06"
}
```

Scenario revision `3` receipt, confirmed by `satoshikawato` on `2026-10-11`:

```text
PRODUCT_DECISION
Concern: web.legend.rename-collision
Scenario revision: 3
Choice: Legend Name Conflict choice "Suffix" -> "Keep both: name it "<name 2>"" (shows the resulting name); "Merge" は `Merge into existing "<name>" and adopt its color` から `Merge into "<name>"` に変える（引き継ぐ色はボタンの色見本で示す）。選択肢と効果は変えない。
Rationale: UX-A10 (show what the suffix is).
Must preserve: Options and effects unchanged.
May retire: the old visible labels only.
Accepted residual risk: older docs/videos show the old labels.
Owner: satoshikawato
Decision date: 2026-10-11
```

```json
{
  "concern": "web.legend.rename-collision",
  "scenarioRevision": 3,
  "choice": "Legend Name Conflict choice \"Suffix\" -> \"Keep both: name it \"<name 2>\"\" (shows the resulting name); \"Merge\" は `Merge into existing \"<name>\" and adopt its color` から `Merge into \"<name>\"` に変える（引き継ぐ色はボタンの色見本で示す）。選択肢と効果は変えない。",
  "rationale": "UX-A10 (show what the suffix is).",
  "mustPreserve": "Options and effects unchanged.",
  "mayRetire": "the old visible labels only.",
  "acceptedResidualRisk": "older docs/videos show the old labels.",
  "owner": "satoshikawato",
  "decisionDate": "2026-10-11"
}
```

Scenario revision `4` receipt, confirmed by `satoshikawato` on `2026-10-13`:

```text
PRODUCT_DECISION
Concern: web.legend.rename-collision
Scenario revision: 4
Choice: B / MERGE-SAME-FEATURE-TYPE-ONLY。色ルールの名前も同じ扱いにする。Feature のポップアップの一覧（Choose from color rules…）から選んだ項目には、名前と色の両方を引き継いで入る。
Rationale: 凡例の名前に色の ID を付けない。一覧から既存の項目を選ぶのは「その項目に入れる」操作なので、色も引き継いで 1 つの項目にまとめる。名前を入力して既存の名前に変えたときは、衝突先が色ルールの名前でも、ほかの項目と同じ選択のダイアログ（Merge into "<name>" / Keep both: name it "<name> (1)" / Cancel）を示す。
Must preserve: 衝突しない rename の即時反映、すべての項目で衝突時に選択のダイアログを示すこと、Merge は同じ track の同じ feature type の項目どうしだけ、This feature only / All N の範囲の選択、Undo/Redo、Generate と Session での保持。
May retire: 衝突先が色ルールの caption のときに確認なしで名前だけを変え、PD-OI-042 の色 ID で区別する動作。一覧に別の feature type の色ルールを示すこと。
Accepted residual risk: 一覧から選ぶと feature の色が変わる（Undo で戻せる）。図に描かれていない色ルール（例: Connector）を選んだときも、そのルールの色になる。
Owner: satoshikawato
Decision date: 2026-10-13
```

```json
{
  "concern": "web.legend.rename-collision",
  "scenarioRevision": 4,
  "choice": "B / MERGE-SAME-FEATURE-TYPE-ONLY。色ルールの名前も同じ扱いにする。Feature のポップアップの一覧（Choose from color rules…）から選んだ項目には、名前と色の両方を引き継いで入る。",
  "rationale": "凡例の名前に色の ID を付けない。一覧から既存の項目を選ぶのは「その項目に入れる」操作なので、色も引き継いで 1 つの項目にまとめる。名前を入力して既存の名前に変えたときは、衝突先が色ルールの名前でも、ほかの項目と同じ選択のダイアログ（Merge into \"<name>\" / Keep both: name it \"<name> (1)\" / Cancel）を示す。",
  "mustPreserve": "衝突しない rename の即時反映、すべての項目で衝突時に選択のダイアログを示すこと、Merge は同じ track の同じ feature type の項目どうしだけ、This feature only / All N の範囲の選択、Undo/Redo、Generate と Session での保持。",
  "mayRetire": "衝突先が色ルールの caption のときに確認なしで名前だけを変え、PD-OI-042 の色 ID で区別する動作。一覧に別の feature type の色ルールを示すこと。",
  "acceptedResidualRisk": "一覧から選ぶと feature の色が変わる（Undo で戻せる）。図に描かれていない色ルール（例: Connector）を選んだときも、そのルールの色になる。",
  "owner": "satoshikawato",
  "decisionDate": "2026-10-13"
}
```

Scenario revision `5` receipt, confirmed by `satoshikawato` on `2026-10-13`:

```text
PRODUCT_DECISION
Concern: web.legend.rename-collision
Scenario revision: 5
Choice: Legend Name Conflict の "Keep both" を選ぶと、名前を変えた行は入力した名前をそのまま持ち、同じ名前の行はすべて凡例に並ぶ順で番号付きになる（"<name> (1)" と "<name> (2)"）。ボタンは今の形のまま、名前を変えた行が受け取る名前を示す（例 Keep both: name it "<name> (2)"）。ダイアログの説明に、もう一方の行の新しい名前を 1 文で示す（The existing "<name>" becomes "<name> (1)".）。番号付きで表示されている名前（例 "<name> (2)"）を入力したときは、その行を衝突先として同じダイアログを示す。選択肢（Merge / Keep both / Cancel）と Merge の条件は変えない。
Rationale: 凡例の名前の番号付けを 1 つの仕組みに揃える。
Must preserve: scenario revision 4 の内容すべて。
May retire: 元の行を番号なしのまま残し、名前を変えた行だけを "<name> (1)" にする名前の付け方。
Accepted residual risk: 以前の付け方で保存した名前（例 "CDS (1)"）は、その文字列のまま読む。
Owner: satoshikawato
Decision date: 2026-10-13
```

```json
{
  "concern": "web.legend.rename-collision",
  "scenarioRevision": 5,
  "choice": "Legend Name Conflict の \"Keep both\" を選ぶと、名前を変えた行は入力した名前をそのまま持ち、同じ名前の行はすべて凡例に並ぶ順で番号付きになる（\"<name> (1)\" と \"<name> (2)\"）。ボタンは今の形のまま、名前を変えた行が受け取る名前を示す（例 Keep both: name it \"<name> (2)\"）。ダイアログの説明に、もう一方の行の新しい名前を 1 文で示す（The existing \"<name>\" becomes \"<name> (1)\".）。番号付きで表示されている名前（例 \"<name> (2)\"）を入力したときは、その行を衝突先として同じダイアログを示す。選択肢（Merge / Keep both / Cancel）と Merge の条件は変えない。",
  "rationale": "凡例の名前の番号付けを 1 つの仕組みに揃える。",
  "mustPreserve": "scenario revision 4 の内容すべて。",
  "mayRetire": "元の行を番号なしのまま残し、名前を変えた行だけを \"<name> (1)\" にする名前の付け方。",
  "acceptedResidualRisk": "以前の付け方で保存した名前（例 \"CDS (1)\"）は、その文字列のまま読む。",
  "owner": "satoshikawato",
  "decisionDate": "2026-10-13"
}
```

### PD-OI-062: Batch live edits projected on Result mount

- Concern key: `web.batch.live-edit-projection`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `B / PROJECT-ON-MOUNT`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its nine-field JSON representation below.
- Decision source: the complete `D-07` receipt in [`02_DECISION_PACK.md`](./web-gui-audit-20260930/02_DECISION_PACK.md) at P00 merge commit
  `e97d90fecfb327135eb50e85cfa9a87be3145823`, approved by `satoshikawato` on `2026-09-30` through the two
  Owner replies quoted verbatim in the Revision 29 entry above. The receipt and
  JSON below reproduce all nine supplied fields without translation or
  additional terms. This record does not supersede another decision.
  Dependent runtime requires this authority merged into its base; this
  amendment supplies no runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `2125389ed12b1b9fd2a874b726ebfc4f03d0a6e518342bec2b953fe4d41bd5b1`.
- Acceptance contracts: `OIC-013`, `OIC-027`. These obligations and the complete
  selected outcome are jointly required; their citation does not claim
  completed dependent-runtime checks.

```text
PRODUCT_DECISION
Concern: web.batch.live-edit-projection
Scenario revision: 1
Choice: B / PROJECT-ON-MOUNT
Rationale: batch の各 Result を選んだとき、その Result の preview と出力に、すでに行った色・非表示・凡例・ラベルの編集が反映されているようにする。編集のたびに全 Result を処理し直すことは避ける。
Must preserve: 表示中の Result への即時反映、Undo/Redo で全 Result の見え方が戻ること、Save → Load → Result 選択、各 Result の export、編集がないときの zero fast path、stale・cancel のときの旧 Result の保持、label の DOM identity。
May retire: 表示していない Result に、Generate まで古い色・非表示・凡例が残る動作。
Accepted residual risk: 一度も表示していない Result は、Session に保存される SVG の bytes が、次に表示するか Generate するまで古い（Load 後に選べば投影される）。大きな batch では表示のたびに投影のコストがかかる。
Owner: satoshikawato
Decision date: 2026-09-30
```

```json
{
  "concern": "web.batch.live-edit-projection",
  "scenarioRevision": 1,
  "choice": "B / PROJECT-ON-MOUNT",
  "rationale": "batch の各 Result を選んだとき、その Result の preview と出力に、すでに行った色・非表示・凡例・ラベルの編集が反映されているようにする。編集のたびに全 Result を処理し直すことは避ける。",
  "mustPreserve": "表示中の Result への即時反映、Undo/Redo で全 Result の見え方が戻ること、Save → Load → Result 選択、各 Result の export、編集がないときの zero fast path、stale・cancel のときの旧 Result の保持、label の DOM identity。",
  "mayRetire": "表示していない Result に、Generate まで古い色・非表示・凡例が残る動作。",
  "acceptedResidualRisk": "一度も表示していない Result は、Session に保存される SVG の bytes が、次に表示するか Generate するまで古い（Load 後に選べば投影される）。大きな batch では表示のたびに投影のコストがかかる。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-30"
}
```

### PD-OI-063: Legend order continuity across Generate

- Concern key: `web.legend.order-continuity`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `B / CARRY-LEGEND-ORDER`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its nine-field JSON representation below.
- Decision source: the complete `D-08` receipt in [`02_DECISION_PACK.md`](./web-gui-audit-20260930/02_DECISION_PACK.md) at P00 merge commit
  `e97d90fecfb327135eb50e85cfa9a87be3145823`, approved by `satoshikawato` on `2026-09-30` through the two
  Owner replies quoted verbatim in the Revision 29 entry above. The receipt and
  JSON below reproduce all nine supplied fields without translation or
  additional terms. This record does not supersede another decision.
  Dependent runtime requires this authority merged into its base; this
  amendment supplies no runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `8c1695739454244618f079ea84586600d5e2f165d7a0d0aa2cfb8de4b167691f`.
- Acceptance contracts: `OIC-027`. These obligations and the complete
  selected outcome are jointly required; their citation does not claim
  completed dependent-runtime checks.

```text
PRODUCT_DECISION
Concern: web.legend.order-continuity
Scenario revision: 1
Choice: B / CARRY-LEGEND-ORDER
Rationale: Sort や Move で整えた凡例の順序を、色や font を直して Generate するたびにやり直さなくてよいようにする（PD-OI-052 と同じ負担をなくす）。
Must preserve: 編集がないときの既定の順序、Sort と Move の即時反映、Undo/Redo、Session の往復、batch の全出力への適用、新しく現れた項目の表示、PD-OI-052 の装飾 delta。
May retire: Generate が凡例の順序を既定に戻す動作。
Accepted residual risk: 並べ替えの後に新しく現れた項目は末尾に置かれる。消えた項目の順序の情報は捨てる。
Owner: satoshikawato
Decision date: 2026-09-30
```

```json
{
  "concern": "web.legend.order-continuity",
  "scenarioRevision": 1,
  "choice": "B / CARRY-LEGEND-ORDER",
  "rationale": "Sort や Move で整えた凡例の順序を、色や font を直して Generate するたびにやり直さなくてよいようにする（PD-OI-052 と同じ負担をなくす）。",
  "mustPreserve": "編集がないときの既定の順序、Sort と Move の即時反映、Undo/Redo、Session の往復、batch の全出力への適用、新しく現れた項目の表示、PD-OI-052 の装飾 delta。",
  "mayRetire": "Generate が凡例の順序を既定に戻す動作。",
  "acceptedResidualRisk": "並べ替えの後に新しく現れた項目は末尾に置かれる。消えた項目の順序の情報は捨てる。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-30"
}
```

### PD-OI-064: Canvas padding continuity across Generate

- Concern key: `web.canvas.padding-continuity`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `B / CARRY-CANVAS-PADDING`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its nine-field JSON representation below.
- Decision source: the complete `D-09` receipt in [`02_DECISION_PACK.md`](./web-gui-audit-20260930/02_DECISION_PACK.md) at P00 merge commit
  `e97d90fecfb327135eb50e85cfa9a87be3145823`, approved by `satoshikawato` on `2026-09-30` through the two
  Owner replies quoted verbatim in the Revision 29 entry above. The receipt and
  JSON below reproduce all nine supplied fields without translation or
  additional terms. This record does not supersede another decision.
  Dependent runtime requires this authority merged into its base; this
  amendment supplies no runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `680081a403e4faca4b75c6824345322645103e3b8e5f44b2d24993040136843d`.
- Acceptance contracts: `OIC-027`. These obligations and the complete
  selected outcome are jointly required; their citation does not claim
  completed dependent-runtime checks.

```text
PRODUCT_DECISION
Concern: web.canvas.padding-continuity
Scenario revision: 1
Choice: B / CARRY-CANVAS-PADDING
Rationale: PD-OI-052 が clipping の緩和策として示す padding を、Generate のたびに入れ直さなくてよいようにする。
Must preserve: padding の編集と即時反映、Reset、Undo/Redo、Session の往復、batch の全出力、export。padding を二重に適用しないこと。
May retire: Generate が canvas padding を 0 に戻す動作。
Accepted residual risk: 図の大きさが大きく変わる設定変更の後も同じ padding が残るので、余白が合わないことがある。
Owner: satoshikawato
Decision date: 2026-09-30
```

```json
{
  "concern": "web.canvas.padding-continuity",
  "scenarioRevision": 1,
  "choice": "B / CARRY-CANVAS-PADDING",
  "rationale": "PD-OI-052 が clipping の緩和策として示す padding を、Generate のたびに入れ直さなくてよいようにする。",
  "mustPreserve": "padding の編集と即時反映、Reset、Undo/Redo、Session の往復、batch の全出力、export。padding を二重に適用しないこと。",
  "mayRetire": "Generate が canvas padding を 0 に戻す動作。",
  "acceptedResidualRisk": "図の大きさが大きく変わる設定変更の後も同じ padding が残るので、余白が合わないことがある。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-30"
}
```

### PD-OI-065: Label edits after source replacement

- Concern key: `web.labels.source-replacement-reconciliation`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `A / PRUNE-UNMATCHED-TARGETS`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its nine-field JSON representation below.
- Decision source: the complete `D-10` receipt in [`02_DECISION_PACK.md`](./web-gui-audit-20260930/02_DECISION_PACK.md) at P00 merge commit
  `e97d90fecfb327135eb50e85cfa9a87be3145823`, approved by `satoshikawato` on `2026-09-30` through the two
  Owner replies quoted verbatim in the Revision 29 entry above. The receipt and
  JSON below reproduce all nine supplied fields without translation or
  additional terms. This record does not supersede another decision.
  Dependent runtime requires this authority merged into its base; this
  amendment supplies no runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `81d6b93f60092310c78f8bcfe0e196329f9b6426c1824783ba761723e6e6cfbe`.

```text
PRODUCT_DECISION
Concern: web.labels.source-replacement-reconciliation
Scenario revision: 1
Choice: A / PRUNE-UNMATCHED-TARGETS
Rationale: 別のゲノムに置き換えて Generate したとき、もう存在しない feature への label の編集だけを外し、残る feature への編集は保つ。
Must preserve: 表示の変化（Result 選択、mount、record 選択、mode、非表示、reflow）では label の override を作成・削除しないこと、bulk の label override（matcher として残す）、Undo による復元、Session の往復。
May retire: source の置き換えのとき、target が 1 つでも消えると label の override をすべて消す動作。
Accepted residual risk: 置き換えた後のゲノムに同じ identity の feature があれば、その override はそのまま適用される。
Owner: satoshikawato
Decision date: 2026-09-30
```

```json
{
  "concern": "web.labels.source-replacement-reconciliation",
  "scenarioRevision": 1,
  "choice": "A / PRUNE-UNMATCHED-TARGETS",
  "rationale": "別のゲノムに置き換えて Generate したとき、もう存在しない feature への label の編集だけを外し、残る feature への編集は保つ。",
  "mustPreserve": "表示の変化（Result 選択、mount、record 選択、mode、非表示、reflow）では label の override を作成・削除しないこと、bulk の label override（matcher として残す）、Undo による復元、Session の往復。",
  "mayRetire": "source の置き換えのとき、target が 1 つでも消えると label の override をすべて消す動作。",
  "acceptedResidualRisk": "置き換えた後のゲノムに同じ identity の feature があれば、その override はそのまま適用される。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-30"
}
```

### PD-OI-066: Live-edit and regeneration parity

- Concern key: `web.live-edit.regeneration-parity`
- Scenario revision: `3`
- Supersedes: `PD-OI-066`, scenario revision `2` (`A / WAIT-FOR-G1-APPLIES`).
  Scenario revision `3` keeps every effect of scenario revision `2` except
  what its receipt states: the wait of an existing live edit is bounded by
  the live applies within the rule of `PD-OI-103` scenario revision `2`, not
  by the G1 applies only. Scenario revision `2` superseded `PD-OI-066`,
  scenario revision `1` (`A / LIVE-EDIT-EQUALS-REGENERATION`), and kept every
  effect of scenario revision `1` except what its receipt states: "the
  response speed of each live edit" is no longer an unconditional rule, only
  for the wait for G1 applies (`PD-OI-103`), and the parity is stated with its
  condition on pending Applies on Generate values.
- Status: `ACCEPTED`
- Selected outcome: `A / WAIT-FOR-LIVE-APPLIES-WITHIN-TARGET`
- Normative outcome: exactly the complete approved scenario revision `3`
  `PRODUCT_DECISION` receipt and its nine-field JSON representation below,
  with the effects of the scenario revision `2` receipt that it does not
  change, and those of the scenario revision `1` receipt that scenario
  revision `2` does not change. The scenario revision `1` and `2` receipts
  are kept as they were approved.
- Decision source (scenario revision `1`): the complete `D-11` receipt in [`02_DECISION_PACK.md`](./web-gui-audit-20260930/02_DECISION_PACK.md) at P00 merge commit
  `e97d90fecfb327135eb50e85cfa9a87be3145823`, approved by `satoshikawato` on `2026-09-30` through the two
  Owner replies quoted verbatim in the Revision 29 entry above. The receipt and
  JSON below reproduce all nine supplied fields without translation or
  additional terms.
- Receipt SHA-256 (scenario revision `1`; UTF-8, excluding the final newline):
  `9298e2bb07d617f9f2562e8168d412a867082f4130b60ec709ef2d7ff005e9e8`.
- Decision source (scenario revision `2`): the Owner replies quoted verbatim
  in the Revision 57 entry above. The second receipt and JSON below reproduce
  all nine confirmed fields without translation or additional terms.
- Receipt SHA-256 (scenario revision `2`; UTF-8, excluding the final newline):
  `8027018ef3be9188c563c3c8a122a9c01b6512a04cb71347a45542bd93e915c0`.
- Decision source (scenario revision `3`): the Owner replies quoted
  verbatim in the Revision 64 entry above. The scenario revision `3`
  receipt and JSON below reproduce all nine fields without translation or
  additional terms. Dependent runtime requires this authority merged into its
  base; this record supplies no runtime acceptance evidence.
- Receipt SHA-256 (scenario revision `3`; UTF-8, excluding the final newline):
  `d051f8b4a669d2be36d0f91c06867d8c4079fdd5605451583d3029ccd3b150a7`.
- Acceptance contracts: `OIC-024`, `OIC-027`. These obligations and the complete
  selected outcome are jointly required; their citation does not claim
  completed dependent-runtime checks.

Scenario revision `1` receipt, approved by `satoshikawato` on `2026-09-30`:

```text
PRODUCT_DECISION
Concern: web.live-edit.regeneration-parity
Scenario revision: 1
Choice: A / LIVE-EDIT-EQUALS-REGENERATION
Rationale: 即時の編集で見えている図が、次の Generate、Session の読み込み、export でも同じになることを保証する。
Must preserve: 各即時編集の応答の速さ、Live edit と Applies on Generate の表示の正確さ（OIC-024）、PD-OI-052 の装飾 delta。即時に編集した Result は、同じ draft から新しく Generate した Result と、対象要素の意味（位置、色、文字、表示）で一致する。Applies on Generate の設定は、Generate の前に Result を変えない。
May retire: 即時の編集と Generate で結果が違ってよいという暗黙の扱い。Generate の compiler で再現できない即時編集は、Applies on Generate に切り替える。
Accepted residual risk: 即時に反映される設定が減ることがある（D-04 と D-05 と同じ方向）。parity を取れない即時編集は退役しうる（D-30）。
Owner: satoshikawato
Decision date: 2026-09-30
```

```json
{
  "concern": "web.live-edit.regeneration-parity",
  "scenarioRevision": 1,
  "choice": "A / LIVE-EDIT-EQUALS-REGENERATION",
  "rationale": "即時の編集で見えている図が、次の Generate、Session の読み込み、export でも同じになることを保証する。",
  "mustPreserve": "各即時編集の応答の速さ、Live edit と Applies on Generate の表示の正確さ（OIC-024）、PD-OI-052 の装飾 delta。即時に編集した Result は、同じ draft から新しく Generate した Result と、対象要素の意味（位置、色、文字、表示）で一致する。Applies on Generate の設定は、Generate の前に Result を変えない。",
  "mayRetire": "即時の編集と Generate で結果が違ってよいという暗黙の扱い。Generate の compiler で再現できない即時編集は、Applies on Generate に切り替える。",
  "acceptedResidualRisk": "即時に反映される設定が減ることがある（D-04 と D-05 と同じ方向）。parity を取れない即時編集は退役しうる（D-30）。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-30"
}
```

Scenario revision `2` receipt, confirmed on `2026-10-17`:

```text
PRODUCT_DECISION
Concern: web.live-edit.regeneration-parity
Scenario revision: 2
Choice: A / WAIT-FOR-G1-APPLIES
Rationale: an apply is not cancelled mid-flight, so an existing live edit made during a G1 apply waits for it; the wait is short (about two G1 applies) with S4-Q2 A and S4-Q3 B.
Must preserve: accurate Live edit and Applies on Generate notes (OIC-024); the PD-OI-052 decoration delta; Applies on Generate settings unchanged before Generate; parity stated as: a live Result equals a fresh Generate from the same draft in the target elements' meaning when no pending Applies on Generate value changes those elements; otherwise it equals a Generate of the committed Session plus the applied live values.
May retire: "the response speed of each live edit" as an unconditional rule, only for the G1 wait.
Accepted residual risk: a wait of up to about 0.9 s after a committed G1 value (about 0.5-0.6 s on tobacco), longer while values keep being committed.
Owner: Satoshi Kawato
Decision date: 2026-10-17
```

```json
{
  "concern": "web.live-edit.regeneration-parity",
  "scenarioRevision": 2,
  "choice": "A / WAIT-FOR-G1-APPLIES",
  "rationale": "an apply is not cancelled mid-flight, so an existing live edit made during a G1 apply waits for it; the wait is short (about two G1 applies) with S4-Q2 A and S4-Q3 B.",
  "mustPreserve": "accurate Live edit and Applies on Generate notes (OIC-024); the PD-OI-052 decoration delta; Applies on Generate settings unchanged before Generate; parity stated as: a live Result equals a fresh Generate from the same draft in the target elements' meaning when no pending Applies on Generate value changes those elements; otherwise it equals a Generate of the committed Session plus the applied live values.",
  "mayRetire": "\"the response speed of each live edit\" as an unconditional rule, only for the G1 wait.",
  "acceptedResidualRisk": "a wait of up to about 0.9 s after a committed G1 value (about 0.5-0.6 s on tobacco), longer while values keep being committed.",
  "owner": "Satoshi Kawato",
  "decisionDate": "2026-10-17"
}
```

Scenario revision `3` receipt, confirmed on `2026-10-18`:

```text
PRODUCT_DECISION
Concern: web.live-edit.regeneration-parity
Scenario revision: 3
Choice: A / WAIT-FOR-LIVE-APPLIES-WITHIN-TARGET
Rationale: the choice matching S5-Q1 / S5-Q2 (A with A / A).
Must preserve: rev 2 (notes, PD-OI-052 delta, parity text).
May retire: the G1-only bound.
Accepted residual risk: an editor edit may wait about two live applies, each within the rule: up to 0.6 s (P2) or 1.5 s (P3) with the margin (0.8 s or 2 s without), longer while values keep coming.
Owner: satoshikawato
Decision date: 2026-10-18
```

```json
{
  "concern": "web.live-edit.regeneration-parity",
  "scenarioRevision": 3,
  "choice": "A / WAIT-FOR-LIVE-APPLIES-WITHIN-TARGET",
  "rationale": "the choice matching S5-Q1 / S5-Q2 (A with A / A).",
  "mustPreserve": "rev 2 (notes, PD-OI-052 delta, parity text).",
  "mayRetire": "the G1-only bound.",
  "acceptedResidualRisk": "an editor edit may wait about two live applies, each within the rule: up to 0.6 s (P2) or 1.5 s (P3) with the margin (0.8 s or 2 s without), longer while values keep coming.",
  "owner": "satoshikawato",
  "decisionDate": "2026-10-18"
}
```

### PD-OI-067: Per-record inferred Linear definitions

- Concern key: `web.linear.file-default-definition`
- Scenario revision: `2`
- Supersedes: `PD-OI-067`, scenario revision `1`
  (`B / PER-RECORD-INFERRED-DEFINITION`), only for the names area of the
  card of a File with one record; every scenario-1 term stays in force.
- Status: `ACCEPTED`
- Selected outcomes: `B / PER-RECORD-INFERRED-DEFINITION` (scenario revision
  `1`) and
  `B / PER-RECORD-INFERRED-DEFINITION, one-record File shows record names on the card`
  (scenario revision `2`)
- Normative outcome: exactly the two approved `PRODUCT_DECISION` receipts and
  their nine-field JSON representations below, read together, and the Owner
  selection of `2026-10-12` recorded below them.
- Decision source (scenario revision `1`): the complete `D-12` receipt in [`02_DECISION_PACK.md`](./web-gui-audit-20260930/02_DECISION_PACK.md) at P00 merge commit
  `e97d90fecfb327135eb50e85cfa9a87be3145823`, approved by `satoshikawato` on `2026-09-30` through the two
  Owner replies quoted verbatim in the Revision 29 entry above. The receipt and
  JSON below reproduce all nine supplied fields without translation or
  additional terms. This record does not supersede another decision.
  Dependent runtime requires this authority merged into its base; this
  amendment supplies no runtime acceptance evidence.
- Receipt SHA-256 (scenario revision `1`; UTF-8, excluding the final newline):
  `daf932cb4892019717ea3ae082d1c799b58a436f909738fbc6239df5d43d265b`.
- Decision source (scenario revision `2`): the Owner replies of `2026-10-11`
  and `2026-10-12` quoted in the Revision 41 entry above, which states what
  each quotation leaves out. The second
  receipt and JSON below reproduce all nine fields of the approved draft
  without translation or additional terms. The runtime that implements it merges with
  it (static Product Contract co-change).
- Receipt SHA-256 (scenario revision `2`; UTF-8, excluding the final newline):
  `93b14a1df028f81ff573dc9c1a60f312b8160e6a7416fe14d0bce8c85be67ee2`.
- Acceptance contracts: `OIC-015`. These obligations and the complete
  selected outcome are jointly required; their citation does not claim
  completed dependent-runtime checks.

```text
PRODUCT_DECISION
Concern: web.linear.file-default-definition
Scenario revision: 1
Choice: B / PER-RECORD-INFERRED-DEFINITION
Rationale: 1 つのファイルに別の生物の record が入っていても、各 record に自分の生物名の定義を付ける。
Must preserve: record ごとの Definition の編集、利用者が file に入力した Definition を全 record に適用すること、全 record が同じ生物のときの現在の表示、Reset で推定値に戻ること、Session の往復、Circular の現在の挙動。定義の優先順位は、record に入力した値 → file に入力した値 → その record 自身の推定値。
May retire: 1 番目の record の推定値を file の既定値として全 record に使う動作。
Accepted residual risk: file 欄の「Using file default」は、利用者が file に入力した値だけを指すようになる。推定値を保存するために Session の形式が変わる場合がある（旧 Session は読み込み時に推定し直す）。
Owner: satoshikawato
Decision date: 2026-09-30
```

```json
{
  "concern": "web.linear.file-default-definition",
  "scenarioRevision": 1,
  "choice": "B / PER-RECORD-INFERRED-DEFINITION",
  "rationale": "1 つのファイルに別の生物の record が入っていても、各 record に自分の生物名の定義を付ける。",
  "mustPreserve": "record ごとの Definition の編集、利用者が file に入力した Definition を全 record に適用すること、全 record が同じ生物のときの現在の表示、Reset で推定値に戻ること、Session の往復、Circular の現在の挙動。定義の優先順位は、record に入力した値 → file に入力した値 → その record 自身の推定値。",
  "mayRetire": "1 番目の record の推定値を file の既定値として全 record に使う動作。",
  "acceptedResidualRisk": "file 欄の「Using file default」は、利用者が file に入力した値だけを指すようになる。推定値を保存するために Session の形式が変わる場合がある（旧 Session は読み込み時に推定し直す）。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-30"
}
```

Scenario revision `2` receipt (the draft approved by `satoshikawato` on
`2026-10-11`; Reply 2 of the Revision 41 entry above):

```text
PRODUCT_DECISION
Concern: web.linear.file-default-definition
Scenario revision: 2
Choice: B / PER-RECORD-INFERRED-DEFINITION, one-record File shows record names on the card
Rationale: record が 1 つの File では「この File の全 record の名前」と record の名前が同じものを指すので、2 か所に分けず File カードに record の名前を直接出す（UX 監査 UX-C06）。
Must preserve: revision 1 の Must preserve すべて（record ごとの Definition の編集、File に入力した Definition を全 record に適用すること、全 record が同じ生物のときの表示、Reset で推定値に戻ること、Session の往復、Circular の挙動、優先順位 record → file → 推定値）。record が 2 つ以上の File では今の File 全体の名前欄をそのまま出す。
May retire: record が 1 つの File で、値の入っていない File 全体の名前欄（Organism / strain, Subtitle）を表示すること。
Accepted residual risk: record が 1 つの File で File 全体の名前を入力するには、File の値の入った Session を読み込む以外の入口がなくなる（record の欄で同じ結果になる）。古いドキュメントと動画は今の配置を示す。
Owner: satoshikawato
Decision date: 2026-10-11
```

```json
{
  "concern": "web.linear.file-default-definition",
  "scenarioRevision": 2,
  "choice": "B / PER-RECORD-INFERRED-DEFINITION, one-record File shows record names on the card",
  "rationale": "record が 1 つの File では「この File の全 record の名前」と record の名前が同じものを指すので、2 か所に分けず File カードに record の名前を直接出す（UX 監査 UX-C06）。",
  "mustPreserve": "revision 1 の Must preserve すべて（record ごとの Definition の編集、File に入力した Definition を全 record に適用すること、全 record が同じ生物のときの表示、Reset で推定値に戻ること、Session の往復、Circular の挙動、優先順位 record → file → 推定値）。record が 2 つ以上の File では今の File 全体の名前欄をそのまま出す。",
  "mayRetire": "record が 1 つの File で、値の入っていない File 全体の名前欄（Organism / strain, Subtitle）を表示すること。",
  "acceptedResidualRisk": "record が 1 つの File で File 全体の名前を入力するには、File の値の入った Session を読み込む以外の入口がなくなる（record の欄で同じ結果になる）。古いドキュメントと動画は今の配置を示す。",
  "owner": "satoshikawato",
  "decisionDate": "2026-10-11"
}
```

Owner selection of `2026-10-12` within scenario revision `2` (Reply 3 of the
Revision 41 entry above; not a receipt field):

- A File with one record that holds a File-level name: option B of the
  question 「A か B か」, whose rule the question gave as
  「カードに出す名前のブロックは常に 1 つ。File の名前があるあいだは dev と同じカード（File のブロックがカード、レコードの欄は Record options の中）で、File の名前が空ならレコードの欄をカードに出す。」
  In it `dev` is the Web app before this change.

### PD-OI-068: GFF3+FASTA record universe

- Concern key: `diagram-generation.gff-fasta-record-universe`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `A / GFF-ANNOTATED-RECORDS-ONLY`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its nine-field JSON representation below.
- Decision source: the complete `D-13` receipt in [`02_DECISION_PACK.md`](./web-gui-audit-20260930/02_DECISION_PACK.md) at P00 merge commit
  `e97d90fecfb327135eb50e85cfa9a87be3145823`, approved by `satoshikawato` on `2026-09-30` through the two
  Owner replies quoted verbatim in the Revision 29 entry above. The receipt and
  JSON below reproduce all nine supplied fields without translation or
  additional terms. This record does not supersede another decision.
  Dependent runtime requires this authority merged into its base; this
  amendment supplies no runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `a4a04f680c03b1e1a13ccf82fd35d62b10c01700546c147ab567d28c47ebd14c`.
- Acceptance contracts: `OIC-015`. These obligations and the complete
  selected outcome are jointly required; their citation does not claim
  completed dependent-runtime checks.

```text
PRODUCT_DECISION
Concern: diagram-generation.gff-fasta-record-universe
Scenario revision: 1
Choice: A / GFF-ANNOTATED-RECORDS-ONLY
Rationale: GFF3+FASTA の record の集合を CLI と同じにし、科学的な出力を変えない。
Must preserve: GFF の行を持つ record の表示（feature が 0 でも region 行や埋め込み ##FASTA を持つものを含む）、FASTA の順序、CLI の出力、PD-OI-018 の他の項目。
May retire: GFF の行を 1 つも持たない FASTA の配列を、record の候補として一覧に出す動作（Generate できない候補）。
Accepted residual risk: FASTA だけにある配列は描けない。描くには GFF に region 行を足す必要がある。
Owner: satoshikawato
Decision date: 2026-09-30
```

```json
{
  "concern": "diagram-generation.gff-fasta-record-universe",
  "scenarioRevision": 1,
  "choice": "A / GFF-ANNOTATED-RECORDS-ONLY",
  "rationale": "GFF3+FASTA の record の集合を CLI と同じにし、科学的な出力を変えない。",
  "mustPreserve": "GFF の行を持つ record の表示（feature が 0 でも region 行や埋め込み ##FASTA を持つものを含む）、FASTA の順序、CLI の出力、PD-OI-018 の他の項目。",
  "mayRetire": "GFF の行を 1 つも持たない FASTA の配列を、record の候補として一覧に出す動作（Generate できない候補）。",
  "acceptedResidualRisk": "FASTA だけにある配列は描けない。描くには GFF に region 行を足す必要がある。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-30"
}
```

### PD-OI-069: Stable-hash single-feature colors with duplicate record IDs

- Concern key: `web.feature-color.duplicate-record-instance`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `A / STABLE-HASH-ONLY`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its nine-field JSON representation below.
- Decision source: the complete `D-14` receipt in [`02_DECISION_PACK.md`](./web-gui-audit-20260930/02_DECISION_PACK.md) at P00 merge commit
  `e97d90fecfb327135eb50e85cfa9a87be3145823`, approved by `satoshikawato` on `2026-09-30` through the two
  Owner replies quoted verbatim in the Revision 29 entry above. The receipt and
  JSON below reproduce all nine supplied fields without translation or
  additional terms. This record does not supersede another decision.
  Dependent runtime requires this authority merged into its base; this
  amendment supplies no runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `8d0e04078648f77ebe017d95a929eb66616e78aea8d736d42575b418977e6f9d`.
- Acceptance contracts: `OIC-005`. These obligations and the complete
  selected outcome are jointly required; their citation does not claim
  completed dependent-runtime checks.

```text
PRODUCT_DECISION
Concern: web.feature-color.duplicate-record-instance
Scenario revision: 1
Choice: A / STABLE-HASH-ONLY
Rationale: Web が作る色ルールを、Python が必ず照合できる値にする（OIPC-C03）。
Must preserve: 重複しない feature の「This feature only」、label の instance 単位の編集、Undo/Redo、Session。
May retire: record ID が重複するとき、rendered instance id を色ルールに書く動作。
Accepted residual risk: 同じ record ID を持つ同一の複製がある場合、「This feature only」の色は両方の複製に適用される。
Owner: satoshikawato
Decision date: 2026-09-30
```

```json
{
  "concern": "web.feature-color.duplicate-record-instance",
  "scenarioRevision": 1,
  "choice": "A / STABLE-HASH-ONLY",
  "rationale": "Web が作る色ルールを、Python が必ず照合できる値にする（OIPC-C03）。",
  "mustPreserve": "重複しない feature の「This feature only」、label の instance 単位の編集、Undo/Redo、Session。",
  "mayRetire": "record ID が重複するとき、rendered instance id を色ルールに書く動作。",
  "acceptedResidualRisk": "同じ record ID を持つ同一の複製がある場合、「This feature only」の色は両方の複製に適用される。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-30"
}
```

### PD-OI-070: Reset Settings scope for Linear record display

- Concern key: `web.reset.linear-record-display`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `A / RESET-LINEAR-RECORD-DISPLAY`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its nine-field JSON representation below.
- Decision source: the complete `D-15` receipt in [`02_DECISION_PACK.md`](./web-gui-audit-20260930/02_DECISION_PACK.md) at P00 merge commit
  `e97d90fecfb327135eb50e85cfa9a87be3145823`, approved by `satoshikawato` on `2026-09-30` through the two
  Owner replies quoted verbatim in the Revision 29 entry above. The receipt and
  JSON below reproduce all nine supplied fields without translation or
  additional terms. This record does not supersede another decision.
  Dependent runtime requires this authority merged into its base; this
  amendment supplies no runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `1000473579a0d63c4d1321bac369987dfc8124e15568ed0ce3215560020942ce`.

```text
PRODUCT_DECISION
Concern: web.reset.linear-record-display
Scenario revision: 1
Choice: A / RESET-LINEAR-RECORD-DISPLAY
Rationale: Reset Settings の範囲を Circular と Linear で揃える。
Must preserve: ファイル、展開された行の record 選択（region_record_id）、file の既定値、depth の割り当て、Undo による復元。
May retire: Reset Settings の後も、Linear の record ごとの Definition・Subtitle・region・逆相補と alignment plan が残る動作。
Accepted residual risk: Reset で Linear の record ごとの表示設定が消える（Undo で戻せる）。
Owner: satoshikawato
Decision date: 2026-09-30
```

```json
{
  "concern": "web.reset.linear-record-display",
  "scenarioRevision": 1,
  "choice": "A / RESET-LINEAR-RECORD-DISPLAY",
  "rationale": "Reset Settings の範囲を Circular と Linear で揃える。",
  "mustPreserve": "ファイル、展開された行の record 選択（region_record_id）、file の既定値、depth の割り当て、Undo による復元。",
  "mayRetire": "Reset Settings の後も、Linear の record ごとの Definition・Subtitle・region・逆相補と alignment plan が残る動作。",
  "acceptedResidualRisk": "Reset で Linear の record ごとの表示設定が消える（Undo で戻せる）。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-30"
}
```

### PD-OI-071: Location-only feature-search positions

- Concern key: `web.feature-search.location-fields`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `A / LOCATION-FIELD-ONLY`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its nine-field JSON representation below.
- Decision source: the complete `D-16` receipt in [`02_DECISION_PACK.md`](./web-gui-audit-20260930/02_DECISION_PACK.md) at P00 merge commit
  `e97d90fecfb327135eb50e85cfa9a87be3145823`, approved by `satoshikawato` on `2026-09-30` through the two
  Owner replies quoted verbatim in the Revision 29 entry above. The receipt and
  JSON below reproduce all nine supplied fields without translation or
  additional terms. This record does not supersede another decision.
  Dependent runtime requires this authority merged into its base; this
  amendment supplies no runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `7fd68a5b264ba1a9b941d592dcfd75a8bf371df9687a431d6983ecc617bcc78b`.

```text
PRODUCT_DECISION
Concern: web.feature-search.location-fields
Scenario revision: 1
Choice: A / LOCATION-FIELD-ONLY
Rationale: 位置の検索と表示を 1 始まりの INSDC 形式に揃え、誤った一致をなくす。
Must preserve: Location での検索（原点をまたぐ feature と分割された feature を含む）、drawer と popup の位置と長さの表示。
May retire: 検索の Start と End の項目（0 始まりの生の値）。
Accepted residual risk: 開始位置の数値だけで検索していた場合は、Location の値で検索し直す必要がある。
Owner: satoshikawato
Decision date: 2026-09-30
```

```json
{
  "concern": "web.feature-search.location-fields",
  "scenarioRevision": 1,
  "choice": "A / LOCATION-FIELD-ONLY",
  "rationale": "位置の検索と表示を 1 始まりの INSDC 形式に揃え、誤った一致をなくす。",
  "mustPreserve": "Location での検索（原点をまたぐ feature と分割された feature を含む）、drawer と popup の位置と長さの表示。",
  "mayRetire": "検索の Start と End の項目（0 始まりの生の値）。",
  "acceptedResidualRisk": "開始位置の数値だけで検索していた場合は、Location の値で検索し直す必要がある。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-30"
}
```

### PD-OI-072: Web PDF physical size

- Concern key: `web.export.pdf-physical-size`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `A / CSS-PX-TO-PT`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its nine-field JSON representation below.
- Decision source: the complete `D-17` receipt in [`02_DECISION_PACK.md`](./web-gui-audit-20260930/02_DECISION_PACK.md) at P00 merge commit
  `e97d90fecfb327135eb50e85cfa9a87be3145823`, approved by `satoshikawato` on `2026-09-30` through the two
  Owner replies quoted verbatim in the Revision 29 entry above. The receipt and
  JSON below reproduce all nine supplied fields without translation or
  additional terms. This record does not supersede another decision.
  Dependent runtime requires this authority merged into its base; this
  amendment supplies no runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `bb555cc2d77b10dffe8ade436876858b54612dcbb17b6ae41df6aac7172bd6ad`.

```text
PRODUCT_DECISION
Concern: web.export.pdf-physical-size
Scenario revision: 1
Choice: A / CSS-PX-TO-PT
Rationale: Web の PDF の物理的な大きさを、CLI（CairoSVG）の PDF と、PNG の DPI に揃える。
Must preserve: PDF の見た目、文字の抽出、ページが 1 枚であること、Web の PNG と SVG の大きさ。
May retire: Web の PDF を 1 px = 1 pt で作る動作。
Accepted residual risk: Web で作る PDF の物理サイズは、これまでの 75% になる。
Owner: satoshikawato
Decision date: 2026-09-30
```

```json
{
  "concern": "web.export.pdf-physical-size",
  "scenarioRevision": 1,
  "choice": "A / CSS-PX-TO-PT",
  "rationale": "Web の PDF の物理的な大きさを、CLI（CairoSVG）の PDF と、PNG の DPI に揃える。",
  "mustPreserve": "PDF の見た目、文字の抽出、ページが 1 枚であること、Web の PNG と SVG の大きさ。",
  "mayRetire": "Web の PDF を 1 px = 1 pt で作る動作。",
  "acceptedResidualRisk": "Web で作る PDF の物理サイズは、これまでの 75% になる。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-30"
}
```

### PD-OI-073: Comparison table coordinate frame

- Concern key: `comparison.table-coordinate-frame`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `C / SEARCH-FRAME-EVERYWHERE`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its nine-field JSON representation below.
- Decision source: the complete `D-18` receipt in [`02_DECISION_PACK.md`](./web-gui-audit-20260930/02_DECISION_PACK.md) at P00 merge commit
  `e97d90fecfb327135eb50e85cfa9a87be3145823`, approved by `satoshikawato` on `2026-09-30` through the two
  Owner replies quoted verbatim in the Revision 29 entry above. The receipt and
  JSON below reproduce all nine supplied fields without translation or
  additional terms. This record does not supersede another decision.
  Dependent runtime requires this authority merged into its base; this
  amendment supplies no runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `586edec020dd3da7e9084881ec83a6f42b254b3c353abb70a076a638d83bc730`.
- Acceptance contracts: `OIC-005`, `OIC-021`. These obligations and the complete
  selected outcome are jointly required; their citation does not claim
  completed dependent-runtime checks.

```text
PRODUCT_DECISION
Concern: comparison.table-coordinate-frame
Scenario revision: 1
Choice: C / SEARCH-FRAME-EVERYWHERE
Rationale: 比較表の座標を表示の向きに関係なく同じ意味にし、逆相補・回転・再アップロード・CLI で、同じ表が同じ相同領域を指すようにする。
Must preserve: LOSAT の raw cache と Save Raw の内容、crop の意味（表は crop 後の record 内の座標）、feature binding を持つ protein 比較の投影、向きを変えない場合の既存の BLAST 表の結果、main で保存された Session の読み込み（読み込み時に一度だけ変換する）。表の座標が record の範囲外なら、検証で止めるか警告する。
May retire: アップロードした表と CLI の -b を表示座標（逆相補の後）として読む動作、JS 側の探索座標から表示座標への変換、範囲外の行を record の外に描く動作。
Accepted residual risk: 逆相補と -b を組み合わせていた CLI の利用者にとって、表の座標の意味が変わる（release note に書く）。crop より前の全長の座標で作った表は、crop すると範囲外として止まる。
Owner: satoshikawato
Decision date: 2026-09-30
```

```json
{
  "concern": "comparison.table-coordinate-frame",
  "scenarioRevision": 1,
  "choice": "C / SEARCH-FRAME-EVERYWHERE",
  "rationale": "比較表の座標を表示の向きに関係なく同じ意味にし、逆相補・回転・再アップロード・CLI で、同じ表が同じ相同領域を指すようにする。",
  "mustPreserve": "LOSAT の raw cache と Save Raw の内容、crop の意味（表は crop 後の record 内の座標）、feature binding を持つ protein 比較の投影、向きを変えない場合の既存の BLAST 表の結果、main で保存された Session の読み込み（読み込み時に一度だけ変換する）。表の座標が record の範囲外なら、検証で止めるか警告する。",
  "mayRetire": "アップロードした表と CLI の -b を表示座標（逆相補の後）として読む動作、JS 側の探索座標から表示座標への変換、範囲外の行を record の外に描く動作。",
  "acceptedResidualRisk": "逆相補と -b を組み合わせていた CLI の利用者にとって、表の座標の意味が変わる（release note に書く）。crop より前の全長の座標で作った表は、crop すると範囲外として止まる。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-30"
}
```

### PD-OI-074: Uploaded comparison table record binding

- Concern key: `comparison.uploaded-table-record-binding`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `B / CONTRADICTION-ERROR-UNKNOWN-WARN`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its nine-field JSON representation below.
- Decision source: the complete `D-20` receipt in [`02_DECISION_PACK.md`](./web-gui-audit-20260930/02_DECISION_PACK.md) at P00 merge commit
  `e97d90fecfb327135eb50e85cfa9a87be3145823`, approved by `satoshikawato` on `2026-09-30` through the two
  Owner replies quoted verbatim in the Revision 29 entry above. The receipt and
  JSON below reproduce all nine supplied fields without translation or
  additional terms. This record does not supersede another decision.
  Dependent runtime requires this authority merged into its base; this
  amendment supplies no runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `a9fbfd8905669bfe9e9d2c511f517e5190a62c42b742074b52893b380333a825`.
- Acceptance contracts: `OIC-005`. These obligations and the complete
  selected outcome are jointly required; their citation does not claim
  completed dependent-runtime checks.

```text
PRODUCT_DECISION
Concern: comparison.uploaded-table-record-binding
Scenario revision: 1
Choice: B / CONTRADICTION-ERROR-UNKNOWN-WARN
Rationale: query と subject を取り違えた表で、誤った相同領域を黙って描かないようにする。
Must preserve: record ID と一致する表の結果、version 接尾辞の違い（.1 など）を許すこと、ID が record と無関係な表の位置による割り当て（警告付き）、CLI と Web で同じ結果。
May retire: 端点と逆の record を指す行や、ほかの record を指す行を、位置のまま描く動作と、metadata の ID と index の矛盾。
Accepted residual risk: ID が record と無関係な表は、今と同じく位置で割り当てられる（警告は出る）。
Owner: satoshikawato
Decision date: 2026-09-30
```

```json
{
  "concern": "comparison.uploaded-table-record-binding",
  "scenarioRevision": 1,
  "choice": "B / CONTRADICTION-ERROR-UNKNOWN-WARN",
  "rationale": "query と subject を取り違えた表で、誤った相同領域を黙って描かないようにする。",
  "mustPreserve": "record ID と一致する表の結果、version 接尾辞の違い（.1 など）を許すこと、ID が record と無関係な表の位置による割り当て（警告付き）、CLI と Web で同じ結果。",
  "mayRetire": "端点と逆の record を指す行や、ほかの record を指す行を、位置のまま描く動作と、metadata の ID と index の矛盾。",
  "acceptedResidualRisk": "ID が record と無関係な表は、今と同じく位置で割り当てられる（警告は出る）。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-30"
}
```

### PD-OI-075: Similarity group name and description identity

- Concern key: `web.similarity-group.override-identity`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `A / REKEY-BY-MEMBERSET-KEEP-DORMANT`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its nine-field JSON representation below.
- Decision source: the complete `D-21` receipt in [`02_DECISION_PACK.md`](./web-gui-audit-20260930/02_DECISION_PACK.md) at P00 merge commit
  `e97d90fecfb327135eb50e85cfa9a87be3145823`, approved by `satoshikawato` on `2026-09-30` through the two
  Owner replies quoted verbatim in the Revision 29 entry above. The receipt and
  JSON below reproduce all nine supplied fields without translation or
  additional terms. This record does not supersede another decision.
  Dependent runtime requires this authority merged into its base; this
  amendment supplies no runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `a4d825e0c630f8a4b147e5a3ae3f179a466a973d7b1b1030c75a022c2e05c62d`.

```text
PRODUCT_DECISION
Concern: web.similarity-group.override-identity
Scenario revision: 1
Choice: A / REKEY-BY-MEMBERSET-KEEP-DORMANT
Rationale: 利用者が付けた group の名前と説明を同じ member の group に付け続け、別の group へ移したり黙って消したりしない。
Must preserve: group の名前と説明の編集、Session の往復、Undo/Redo、Interactive SVG への出力、og_* の ID の表示。
May retire: ID の文字列だけを頼りに名前を残す動作と、残らない名前を黙って消す動作。
Accepted residual risk: member が 1 つでも変わった group には名前が付かない（dormant として保存し、一覧と Clear から扱える）。Session に項目が 1 つ増える。
Owner: satoshikawato
Decision date: 2026-09-30
```

```json
{
  "concern": "web.similarity-group.override-identity",
  "scenarioRevision": 1,
  "choice": "A / REKEY-BY-MEMBERSET-KEEP-DORMANT",
  "rationale": "利用者が付けた group の名前と説明を同じ member の group に付け続け、別の group へ移したり黙って消したりしない。",
  "mustPreserve": "group の名前と説明の編集、Session の往復、Undo/Redo、Interactive SVG への出力、og_* の ID の表示。",
  "mayRetire": "ID の文字列だけを頼りに名前を残す動作と、残らない名前を黙って消す動作。",
  "acceptedResidualRisk": "member が 1 つでも変わった group には名前が付かない（dormant として保存し、一覧と Clear から扱える）。Session に項目が 1 つ増える。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-30"
}
```

### PD-OI-076: Match popup coordinates

- Concern key: `web.match-popup.coordinates`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `B / SOURCE-PRIMARY-WITH-TABLE`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its nine-field JSON representation below.
- Decision source: the complete `D-22` receipt in [`02_DECISION_PACK.md`](./web-gui-audit-20260930/02_DECISION_PACK.md) at P00 merge commit
  `e97d90fecfb327135eb50e85cfa9a87be3145823`, approved by `satoshikawato` on `2026-09-30` through the two
  Owner replies quoted verbatim in the Revision 29 entry above. The receipt and
  JSON below reproduce all nine supplied fields without translation or
  additional terms. This record does not supersede another decision.
  Dependent runtime requires this authority merged into its base; this
  amendment supplies no runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `8986c071c749e19bbd68cd51d6a913f7fb902e4cbaf967d4229334027c113bc9`.

```text
PRODUCT_DECISION
Concern: web.match-popup.coordinates
Scenario revision: 1
Choice: B / SOURCE-PRIMARY-WITH-TABLE
Rationale: match popup と FASTA ヘッダの座標を、feature popup と同じ入力ファイルの座標にする。
Must preserve: 取り出す配列そのもの、逆鎖の扱い、表の座標の参照（違うときだけ併記）、Interactive SVG の popup。
May retire: match popup と FASTA ヘッダが crop 後の表示座標だけを出す動作。
Accepted residual risk: popup の行が 1 行増えることがある。
Owner: satoshikawato
Decision date: 2026-09-30
```

```json
{
  "concern": "web.match-popup.coordinates",
  "scenarioRevision": 1,
  "choice": "B / SOURCE-PRIMARY-WITH-TABLE",
  "rationale": "match popup と FASTA ヘッダの座標を、feature popup と同じ入力ファイルの座標にする。",
  "mustPreserve": "取り出す配列そのもの、逆鎖の扱い、表の座標の参照（違うときだけ併記）、Interactive SVG の popup。",
  "mayRetire": "match popup と FASTA ヘッダが crop 後の表示座標だけを出す動作。",
  "acceptedResidualRisk": "popup の行が 1 行増えることがある。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-30"
}
```

### PD-OI-077: outfmt 6 tables with extra columns

- Concern key: `comparison.outfmt6-extra-columns`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `B / FIRST-12-COLUMNS`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its nine-field JSON representation below.
- Decision source: the complete `D-23` receipt in [`02_DECISION_PACK.md`](./web-gui-audit-20260930/02_DECISION_PACK.md) at P00 merge commit
  `e97d90fecfb327135eb50e85cfa9a87be3145823`, approved by `satoshikawato` on `2026-09-30` through the two
  Owner replies quoted verbatim in the Revision 29 entry above. The receipt and
  JSON below reproduce all nine supplied fields without translation or
  additional terms. This record does not supersede another decision.
  Dependent runtime requires this authority merged into its base; this
  amendment supplies no runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `6c49ba56b74a05b6e19ae38422befad218884821281c48e71eed26b6d49217e6`.
- Acceptance contracts: `OIC-005`. These obligations and the complete
  selected outcome are jointly required; their citation does not claim
  completed dependent-runtime checks.

```text
PRODUCT_DECISION
Concern: comparison.outfmt6-extra-columns
Scenario revision: 1
Choice: B / FIRST-12-COLUMNS
Rationale: -outfmt "6 std qlen slen" のように列を足した表を、CLI と Web でそのまま使えるようにする。
Must preserve: 12 列の表の結果、outfmt 7 のコメント行、空のファイル、先頭 12 列の型の検証、CLI・Web・conservation で同じ規則。
May retire: 列を足した表を黙って誤読する動作。存在しないファイルや読めないファイルを飛ばして、後ろの比較をずらす動作。
Accepted residual risk: 13 列目以降は使わずに捨てる（INFO ログを出す）。存在しないファイルを渡していた CLI の実行は失敗に変わる。
Owner: satoshikawato
Decision date: 2026-09-30
```

```json
{
  "concern": "comparison.outfmt6-extra-columns",
  "scenarioRevision": 1,
  "choice": "B / FIRST-12-COLUMNS",
  "rationale": "-outfmt \"6 std qlen slen\" のように列を足した表を、CLI と Web でそのまま使えるようにする。",
  "mustPreserve": "12 列の表の結果、outfmt 7 のコメント行、空のファイル、先頭 12 列の型の検証、CLI・Web・conservation で同じ規則。",
  "mayRetire": "列を足した表を黙って誤読する動作。存在しないファイルや読めないファイルを飛ばして、後ろの比較をずらす動作。",
  "acceptedResidualRisk": "13 列目以降は使わずに捨てる（INFO ログを出す）。存在しないファイルを渡していた CLI の実行は失敗に変わる。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-30"
}
```

### PD-OI-078: Circular definition wrap on fit failure

- Concern key: `diagram-generation.circular-definition-fit`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `B / WRAP-DEFINITION-ON-FIT-FAILURE`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its nine-field JSON representation below.
- Decision source: the complete `D-24` receipt in [`02_DECISION_PACK.md`](./web-gui-audit-20260930/02_DECISION_PACK.md) at P00 merge commit
  `e97d90fecfb327135eb50e85cfa9a87be3145823`, approved by `satoshikawato` on `2026-09-30` through the two
  Owner replies quoted verbatim in the Revision 29 entry above. The receipt and
  JSON below reproduce all nine supplied fields without translation or
  additional terms. This record does not supersede another decision.
  Dependent runtime requires this authority merged into its base; this
  amendment supplies no runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `700a70feb21766e8d634b765ff9c5519cebe7871110a2037f7ec3f41cdfe1e80`.

```text
PRODUCT_DECISION
Concern: diagram-generation.circular-definition-fit
Scenario revision: 1
Choice: B / WRAP-DEFINITION-ON-FIT-FAILURE
Rationale: よくある細菌の長い学名（subsp.、serovar、str. などを含むもの）でも、Web の既定の設定で図を作れるようにする。
Must preserve: これまで成功していた出力（折り返さない）、定義の文字そのもの、center_reserved_radius と definition_font_size を明示した場合の扱い（折り返しを適用しない）、失敗したときの案内。
May retire: 定義の円が入りきらないとき、配置し直さずに失敗する動作。
Accepted residual risk: 折り返しても入らない場合は今までどおり失敗し、改善した案内を出す。
Owner: satoshikawato
Decision date: 2026-09-30
```

```json
{
  "concern": "diagram-generation.circular-definition-fit",
  "scenarioRevision": 1,
  "choice": "B / WRAP-DEFINITION-ON-FIT-FAILURE",
  "rationale": "よくある細菌の長い学名（subsp.、serovar、str. などを含むもの）でも、Web の既定の設定で図を作れるようにする。",
  "mustPreserve": "これまで成功していた出力（折り返さない）、定義の文字そのもの、center_reserved_radius と definition_font_size を明示した場合の扱い（折り返しを適用しない）、失敗したときの案内。",
  "mayRetire": "定義の円が入りきらないとき、配置し直さずに失敗する動作。",
  "acceptedResidualRisk": "折り返しても入らない場合は今までどおり失敗し、改善した案内を出す。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-30"
}
```

### PD-OI-079: Generate before saving a legacy Session

- Concern key: `web.session.legacy-save`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `A / REQUIRE-GENERATE-BEFORE-SAVE`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its nine-field JSON representation below.
- Decision source: the complete `D-25` receipt in [`02_DECISION_PACK.md`](./web-gui-audit-20260930/02_DECISION_PACK.md) at P00 merge commit
  `e97d90fecfb327135eb50e85cfa9a87be3145823`, approved by `satoshikawato` on `2026-09-30` through the two
  Owner replies quoted verbatim in the Revision 29 entry above. The receipt and
  JSON below reproduce all nine supplied fields without translation or
  additional terms. This record does not supersede another decision.
  Dependent runtime requires this authority merged into its base; this
  amendment supplies no runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `46e30d472ee232ae1228fb2fc17f85f2e9509fc19348236f4e53b9307e7d9150`.

```text
PRODUCT_DECISION
Concern: web.session.legacy-save
Scenario revision: 1
Choice: A / REQUIRE-GENERATE-BEFORE-SAVE
Rationale: 旧形式の Session から、feature の identity が確かでない状態のまま現行の形式を書き出さない。
Must preserve: 旧形式の Session の読み込みと preview、Generate 後の Save、v40 以降の Session の Save、PD-OI-045 の Session 操作。Save が必要とする Generate を案内し、エラーパネルから Generate を実行できるようにする。
May retire: 0.13.0 で可能だった「旧形式の Session を読み込んで、そのまま Save する」操作。
Accepted residual risk: 旧形式の Session を保存し直すには 1 回 Generate が必要。
Owner: satoshikawato
Decision date: 2026-09-30
```

```json
{
  "concern": "web.session.legacy-save",
  "scenarioRevision": 1,
  "choice": "A / REQUIRE-GENERATE-BEFORE-SAVE",
  "rationale": "旧形式の Session から、feature の identity が確かでない状態のまま現行の形式を書き出さない。",
  "mustPreserve": "旧形式の Session の読み込みと preview、Generate 後の Save、v40 以降の Session の Save、PD-OI-045 の Session 操作。Save が必要とする Generate を案内し、エラーパネルから Generate を実行できるようにする。",
  "mayRetire": "0.13.0 で可能だった「旧形式の Session を読み込んで、そのまま Save する」操作。",
  "acceptedResidualRisk": "旧形式の Session を保存し直すには 1 回 Generate が必要。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-30"
}
```

### PD-OI-080: Dinucleotide alphabet

- Concern key: `options.dinucleotide-alphabet`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `B / ACGTU-PAIRS`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its nine-field JSON representation below.
- Decision source: the complete `D-26` receipt in [`02_DECISION_PACK.md`](./web-gui-audit-20260930/02_DECISION_PACK.md) at P00 merge commit
  `e97d90fecfb327135eb50e85cfa9a87be3145823`, approved by `satoshikawato` on `2026-09-30` through the two
  Owner replies quoted verbatim in the Revision 29 entry above. The receipt and
  JSON below reproduce all nine supplied fields without translation or
  additional terms. This record does not supersede another decision.
  Dependent runtime requires this authority merged into its base; this
  amendment supplies no runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `0461b20190328105e23fb0557a22ae7f0643adb679eecb92e09e553bf78f6d50`.

```text
PRODUCT_DECISION
Concern: options.dinucleotide-alphabet
Scenario revision: 1
Choice: B / ACGTU-PAIRS
Rationale: 無効な指定で空や平坦な track を黙って描かないようにし、CLI の traceback もなくす。RNA の表記（U）でも指定できるようにする。
Must preserve: ACGT の 2 文字の指定（大小を区別しない）、slot の nt、CLI と Web で同じ検証。U は T と同じ塩基として扱う（AU は AT と同じ結果になり、配列中の U も T として数える）。
May retire: 2 文字でない指定を黙って既定値に戻す動作、XY のような塩基でない文字の受理、G での IndexError。
Accepted residual risk: N などの曖昧な塩基の記号は指定できない。凡例などの表示名は入力した文字（AU）のまま出す。
Owner: satoshikawato
Decision date: 2026-09-30
```

```json
{
  "concern": "options.dinucleotide-alphabet",
  "scenarioRevision": 1,
  "choice": "B / ACGTU-PAIRS",
  "rationale": "無効な指定で空や平坦な track を黙って描かないようにし、CLI の traceback もなくす。RNA の表記（U）でも指定できるようにする。",
  "mustPreserve": "ACGT の 2 文字の指定（大小を区別しない）、slot の nt、CLI と Web で同じ検証。U は T と同じ塩基として扱う（AU は AT と同じ結果になり、配列中の U も T として数える）。",
  "mayRetire": "2 文字でない指定を黙って既定値に戻す動作、XY のような塩基でない文字の受理、G での IndexError。",
  "acceptedResidualRisk": "N などの曖昧な塩基の記号は指定できない。凡例などの表示名は入力した文字（AU）のまま出す。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-30"
}
```

### PD-OI-081: Validated font sizes and stroke widths

- Concern key: `options.nonnegative-style-values`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `A / SPEC-DETERMINED-ONLY`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its nine-field JSON representation below.
- Decision source: the complete `D-27` receipt in [`02_DECISION_PACK.md`](./web-gui-audit-20260930/02_DECISION_PACK.md) at P00 merge commit
  `e97d90fecfb327135eb50e85cfa9a87be3145823`, approved by `satoshikawato` on `2026-09-30` through the two
  Owner replies quoted verbatim in the Revision 29 entry above. The receipt and
  JSON below reproduce all nine supplied fields without translation or
  additional terms. This record does not supersede another decision.
  Dependent runtime requires this authority merged into its base; this
  amendment supplies no runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `17c67282ae607cc5a4c36af8945337e6835ffb3758c3b5cd2350ac896100fb8b`.

```text
PRODUCT_DECISION
Concern: options.nonnegative-style-values
Scenario revision: 1
Choice: A / SPEC-DETERMINED-ONLY
Rationale: SVG と CSS の仕様で意味が決まる値だけを検証し、意味を確かめていない値は変えない。
Must preserve: offset、spacing、track_axis_gap、label_rotation の現在の受理範囲。CLI、Python API、Web、Session で同じ検証とエラー。
May retire: 0 以下のフォントサイズと負の stroke 幅の受理、描画の途中の ValueError による traceback。
Accepted residual risk: 負の offset などの意味は、今回は確かめない。
Owner: satoshikawato
Decision date: 2026-09-30
```

```json
{
  "concern": "options.nonnegative-style-values",
  "scenarioRevision": 1,
  "choice": "A / SPEC-DETERMINED-ONLY",
  "rationale": "SVG と CSS の仕様で意味が決まる値だけを検証し、意味を確かめていない値は変えない。",
  "mustPreserve": "offset、spacing、track_axis_gap、label_rotation の現在の受理範囲。CLI、Python API、Web、Session で同じ検証とエラー。",
  "mayRetire": "0 以下のフォントサイズと負の stroke 幅の受理、描画の途中の ValueError による traceback。",
  "acceptedResidualRisk": "負の offset などの意味は、今回は確かめない。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-30"
}
```

### PD-OI-082: Undo and Redo during Generate

- Concern key: `web.generation.in-flight-history`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `A / BUSY-UNDO-DURING-GENERATE`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its nine-field JSON representation below.
- Decision source: the complete `D-28` receipt in [`02_DECISION_PACK.md`](./web-gui-audit-20260930/02_DECISION_PACK.md) at P00 merge commit
  `e97d90fecfb327135eb50e85cfa9a87be3145823`, approved by `satoshikawato` on `2026-09-30` through the two
  Owner replies quoted verbatim in the Revision 29 entry above. The receipt and
  JSON below reproduce all nine supplied fields without translation or
  additional terms. This record does not supersede another decision.
  Dependent runtime requires this authority merged into its base; this
  amendment supplies no runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `e8605fbcc01935848811be174be8adc383ce46cd9abad58f3162ce5372749758`.
- Acceptance contracts: `OIC-013`. These obligations and the complete
  selected outcome are jointly required; their citation does not claim
  completed dependent-runtime checks.

```text
PRODUCT_DECISION
Concern: web.generation.in-flight-history
Scenario revision: 1
Choice: A / BUSY-UNDO-DURING-GENERATE
Rationale: Generate 中に Undo や Redo を押しても、確定済みの request と Result が古いものに戻ったり、実行中の Generate が黙って捨てられたりしないようにする。
Must preserve: Generate の Cancel、処理中の表示、Generate が終わった後の Undo/Redo、Save と Load の拒否（既存）、PD-OI-051 による Generate 中の draft の編集。Undo と Redo のボタンとショートカットは同じ判定を使い、拒否した理由を示す。
May retire: Generate 中の Undo と Redo。
Accepted residual risk: 長い LOSAT の実行中は、Cancel するか終わるまで Undo できない。
Owner: satoshikawato
Decision date: 2026-09-30
```

```json
{
  "concern": "web.generation.in-flight-history",
  "scenarioRevision": 1,
  "choice": "A / BUSY-UNDO-DURING-GENERATE",
  "rationale": "Generate 中に Undo や Redo を押しても、確定済みの request と Result が古いものに戻ったり、実行中の Generate が黙って捨てられたりしないようにする。",
  "mustPreserve": "Generate の Cancel、処理中の表示、Generate が終わった後の Undo/Redo、Save と Load の拒否（既存）、PD-OI-051 による Generate 中の draft の編集。Undo と Redo のボタンとショートカットは同じ判定を使い、拒否した理由を示す。",
  "mayRetire": "Generate 中の Undo と Redo。",
  "acceptedResidualRisk": "長い LOSAT の実行中は、Cancel するか終わるまで Undo できない。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-30"
}
```

### PD-OI-083: Linear Depth series without a source

- Concern key: `web.depth.linear-sourceless-series`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `A / ROW-ISSUE-BEFORE-GENERATE`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its nine-field JSON representation below.
- Decision source: the complete `D-29` receipt in [`02_DECISION_PACK.md`](./web-gui-audit-20260930/02_DECISION_PACK.md) at P00 merge commit
  `e97d90fecfb327135eb50e85cfa9a87be3145823`, approved by `satoshikawato` on `2026-09-30` through the two
  Owner replies quoted verbatim in the Revision 29 entry above. The receipt and
  JSON below reproduce all nine supplied fields without translation or
  additional terms. This record does not supersede another decision.
  Dependent runtime requires this authority merged into its base; this
  amendment supplies no runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `a79ed36bc2c8edc21df4a5bc4e39e33f96a53adcf4a9b4b9018d1cd25c6f0da6`.
- Acceptance contracts: `OIC-020`. These obligations and the complete
  selected outcome are jointly required; their citation does not claim
  completed dependent-runtime checks.

```text
PRODUCT_DECISION
Concern: web.depth.linear-sourceless-series
Scenario revision: 1
Choice: A / ROW-ISSUE-BEFORE-GENERATE
Rationale: Linear で File の Depth を消した後、source を持たない series を有効な手動の行が参照していても、原因の分からない失敗にしない。
Must preserve: PD-OI-025（論理 series の保持、File 単位の apply と clear が 1 つの undoable 操作）、OIC-020、Circular の row issue と同じ文言。
May retire: この状態の Generate が汎用のエラーで失敗する動作。
Accepted residual risk: 利用者が行を無効にするか削除する必要がある。
Owner: satoshikawato
Decision date: 2026-09-30
```

```json
{
  "concern": "web.depth.linear-sourceless-series",
  "scenarioRevision": 1,
  "choice": "A / ROW-ISSUE-BEFORE-GENERATE",
  "rationale": "Linear で File の Depth を消した後、source を持たない series を有効な手動の行が参照していても、原因の分からない失敗にしない。",
  "mustPreserve": "PD-OI-025（論理 series の保持、File 単位の apply と clear が 1 つの undoable 操作）、OIC-020、Circular の row issue と同じ文言。",
  "mayRetire": "この状態の Generate が汎用のエラーで失敗する動作。",
  "acceptedResidualRisk": "利用者が行を無効にするか削除する必要がある。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-30"
}
```

### PD-OI-084: Linear live legend side moves

- Concern key: `web.legend.live-side-move`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `A / PARITY-OR-APPLY-ON-GENERATE`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its nine-field JSON representation below.
- Decision source: the complete `D-30` receipt in [`02_DECISION_PACK.md`](./web-gui-audit-20260930/02_DECISION_PACK.md) at P00 merge commit
  `e97d90fecfb327135eb50e85cfa9a87be3145823`, approved by `satoshikawato` on `2026-09-30` through the two
  Owner replies quoted verbatim in the Revision 29 entry above. The receipt and
  JSON below reproduce all nine supplied fields without translation or
  additional terms. This record does not supersede another decision.
  Dependent runtime requires this authority merged into its base; this
  amendment supplies no runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `e534794a0cb4ae542c9f313f1e3dce666e7f121ddebe360b03a59420112fe6bb`.
- Acceptance contracts: `OIC-024`, `OIC-027`. These obligations and the complete
  selected outcome are jointly required; their citation does not claim
  completed dependent-runtime checks.

```text
PRODUCT_DECISION
Concern: web.legend.live-side-move
Scenario revision: 1
Choice: A / PARITY-OR-APPLY-ON-GENERATE
Rationale: 即時に見えた凡例の配置と、Generate 後の配置が食い違わないようにする（D-11 の契約）。
Must preserve: Circular の凡例の side の即時移動、凡例の drag、PD-OI-052 の装飾 delta、Undo/Redo。Linear で parity を取れる場合は、即時移動も残す。
May retire: parity を取れない場合に限り、Linear の凡例の side の即時移動（凡例を持たない図での side 変更の例外を含む）。
Accepted residual risk: 退役した場合、Linear では side を変えても Generate まで preview が変わらない。
Owner: satoshikawato
Decision date: 2026-09-30
```

```json
{
  "concern": "web.legend.live-side-move",
  "scenarioRevision": 1,
  "choice": "A / PARITY-OR-APPLY-ON-GENERATE",
  "rationale": "即時に見えた凡例の配置と、Generate 後の配置が食い違わないようにする（D-11 の契約）。",
  "mustPreserve": "Circular の凡例の side の即時移動、凡例の drag、PD-OI-052 の装飾 delta、Undo/Redo。Linear で parity を取れる場合は、即時移動も残す。",
  "mayRetire": "parity を取れない場合に限り、Linear の凡例の side の即時移動（凡例を持たない図での side 変更の例外を含む）。",
  "acceptedResidualRisk": "退役した場合、Linear では side を変えても Generate まで preview が変わらない。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-30"
}
```

### PD-OI-085: Feature-popup record rotation Apply on Generate

- Concern key: `web.feature-popup.record-rotation-apply-on-generate`
- Scenario revision: `2`
- Supersedes: `PD-OI-085`, scenario revision `1` (`A / APPLY_ON_GENERATE`),
  only by renaming the two rotation buttons; every scenario-1 term stays in
  force.
- Status: `ACCEPTED`
- Selected outcomes: `A / APPLY_ON_GENERATE` (scenario revision `1`) and
  `RENAME_ROTATION_APPLY_BUTTONS` (scenario revision `2`)
- Normative outcome: exactly the two approved `PRODUCT_DECISION` receipts and
  their nine-field JSON representations below, read together. Where the
  scenario-1 receipt names **Apply on Generate** and **Apply and regenerate**,
  it means the buttons now labelled **Apply when I Generate** and
  **Apply now**.
- Decision source (scenario revision `1`): the Owner answers of `2026-10-04`
  quoted verbatim in the Revision 30 entry above (`OD-2`, `OD-3`, `OD-4`), and
  Appendix A of
  [`POPUP_RECORD_ACTIONS_AND_VIBRIO_SESSION_PLAN_2026-10-04.md`](./POPUP_RECORD_ACTIONS_AND_VIBRIO_SESSION_PLAN_2026-10-04.md).
  The receipt fields restate those answers and the plan presented in that
  session; the Owner did not separately review the receipt wording. This record
  does not supersede another decision; `PD-OI-032` retains its scope. Dependent
  runtime requires this authority merged into its base; this amendment supplies
  no runtime acceptance evidence.
- Receipt SHA-256 (scenario revision `1`; UTF-8, excluding the final newline):
  `38f89876744c098aa5695b570a0250c6f672272594aafc6c1e1e739b69712f40`.
- Decision source (scenario revision `2`): the Owner replies of `2026-10-11`
  quoted verbatim in the Revision 40 entry above. The second receipt and JSON
  below reproduce the confirmed fields without translation or additional
  terms. The runtime that shows the new names merges with it (static Product
  Contract co-change).
- Receipt SHA-256 (scenario revision `2`; UTF-8, excluding the final newline):
  `2942593046d32fe65bc3ed1c858fc57073518154adaa961be45e2f0b9b23c66d`.
- Acceptance contracts: `OIC-021`, `OIC-024`. These obligations and the
  complete selected outcome are jointly required; their citation does not claim
  completed dependent-runtime checks.

```text
PRODUCT_DECISION
Concern: web.feature-popup.record-rotation-apply-on-generate
Scenario revision: 1
Choice: A / APPLY_ON_GENERATE
Rationale: Record actions に Apply on Generate を加える。開いた feature の record の表示開始位置と、指定したときの向きを、sidebar と同じ未適用の record display 設定に書き込むだけで、Result は再生成しない。Generate Diagram を 1 回押すと、予約したすべての record がまとめて描画される。複数 record の回転を予約してから一度に描画したいという要望（画像編集ソフトの「適用」と「OK」の区別）に応える。
Must preserve: PD-OI-032 の対象特定・anchor・offset・orientation・feature-end・適用前 preview・理由表示。Apply and regenerate は最後の committed request から対象 record だけの candidate を作り、他の record の予約は予約のまま残す（PD-OI-032 item 4）。sidebar の record display 操作とその意味。Generate Diagram による未適用設定の一括適用。同じ record への再予約は後の値が有効。予約は 1 回の Undo/Redo で戻せる。popup を開き直すと予約済みの値が分かる。既存の pending Generate 通知のほかに、Result と Generate への常時 Pending/Applied 表示を加えない（PD-OI-037 revision 2）。
May retire: none
Accepted residual risk: 予約したまま Generate しないと、表示中の Result と設定がずれたままになる。既存の pending Generate 通知と popup での予約値表示で補う。
Owner: satoshikawato
Decision date: 2026-10-04
```

```json
{
  "concern": "web.feature-popup.record-rotation-apply-on-generate",
  "scenarioRevision": 1,
  "choice": "A / APPLY_ON_GENERATE",
  "rationale": "Record actions に Apply on Generate を加える。開いた feature の record の表示開始位置と、指定したときの向きを、sidebar と同じ未適用の record display 設定に書き込むだけで、Result は再生成しない。Generate Diagram を 1 回押すと、予約したすべての record がまとめて描画される。複数 record の回転を予約してから一度に描画したいという要望（画像編集ソフトの「適用」と「OK」の区別）に応える。",
  "mustPreserve": "PD-OI-032 の対象特定・anchor・offset・orientation・feature-end・適用前 preview・理由表示。Apply and regenerate は最後の committed request から対象 record だけの candidate を作り、他の record の予約は予約のまま残す（PD-OI-032 item 4）。sidebar の record display 操作とその意味。Generate Diagram による未適用設定の一括適用。同じ record への再予約は後の値が有効。予約は 1 回の Undo/Redo で戻せる。popup を開き直すと予約済みの値が分かる。既存の pending Generate 通知のほかに、Result と Generate への常時 Pending/Applied 表示を加えない（PD-OI-037 revision 2）。",
  "mayRetire": "none",
  "acceptedResidualRisk": "予約したまま Generate しないと、表示中の Result と設定がずれたままになる。既存の pending Generate 通知と popup での予約値表示で補う。",
  "owner": "satoshikawato",
  "decisionDate": "2026-10-04"
}
```

Scenario revision `2` receipt, confirmed by `satoshikawato` on `2026-10-11`:

```text
PRODUCT_DECISION
Concern: web.feature-popup.record-rotation-apply-on-generate
Scenario revision: 2
Choice: RENAME_ROTATION_APPLY_BUTTONS — "Apply on Generate" -> "Apply when I Generate"; "Apply and regenerate" -> "Apply now"
Rationale: 「regenerate」が Generate ボタンと紛らわしい（UX 監査 UX-B09）。ボタン名でそれぞれの結果を示す。
Must preserve: 両ボタンの動作は変えない（Apply when I Generate は sidebar と同じ未適用の record display 設定に予約するだけで Result を再生成しない。Apply now は最後の committed request から対象 record だけの candidate を作り、他の record の予約は予約のまま残す）。PD-OI-032 の対象特定・anchor・offset・orientation・feature-end・適用前 preview・理由表示。各ボタンのアクセシブル名は表示名と同じ。
May retire: 旧ボタン名「Apply on Generate」「Apply and regenerate」の表示のみ（動作は退役しない）。
Accepted residual risk: 旧ボタン名で書かれた既存の資料（v0.14.0 リリースノートの記述、収録済みのチュートリアル動画）と画面の表記が一致しない。
Owner: satoshikawato
Decision date: 2026-10-11
```

```json
{
  "concern": "web.feature-popup.record-rotation-apply-on-generate",
  "scenarioRevision": 2,
  "choice": "RENAME_ROTATION_APPLY_BUTTONS — \"Apply on Generate\" -> \"Apply when I Generate\"; \"Apply and regenerate\" -> \"Apply now\"",
  "rationale": "「regenerate」が Generate ボタンと紛らわしい（UX 監査 UX-B09）。ボタン名でそれぞれの結果を示す。",
  "mustPreserve": "両ボタンの動作は変えない（Apply when I Generate は sidebar と同じ未適用の record display 設定に予約するだけで Result を再生成しない。Apply now は最後の committed request から対象 record だけの candidate を作り、他の record の予約は予約のまま残す）。PD-OI-032 の対象特定・anchor・offset・orientation・feature-end・適用前 preview・理由表示。各ボタンのアクセシブル名は表示名と同じ。",
  "mayRetire": "旧ボタン名「Apply on Generate」「Apply and regenerate」の表示のみ（動作は退役しない）。",
  "acceptedResidualRisk": "旧ボタン名で書かれた既存の資料（v0.14.0 リリースノートの記述、収録済みのチュートリアル動画）と画面の表記が一致しない。",
  "owner": "satoshikawato",
  "decisionDate": "2026-10-11"
}
```

### PD-OI-087: Feature-based display start only in the feature popup

- Concern key: `web.record-display.feature-start-shortcuts`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `A / POPUP_ONLY_FEATURE_ROTATION`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its nine-field JSON representation below.
- Decision source: the Owner replies of `2026-10-08` quoted verbatim in the
  Revision 34 entry above. The receipt and JSON below reproduce all nine
  fields without translation or additional terms. This record does not
  supersede another decision. It narrows the preserved sidebar operations of
  `PD-OI-032`, `PD-OI-033` scenario revision `2`, and `PD-OI-085` by the two
  retired buttons only, as the Revision 34 entry states. The runtime that
  implements it merges with it (static Product Contract co-change).
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `91529d4e7a6d0729497d13bf8b552be59e4bac39b641e9aabfdcc13f3c06b6c1`.
- Acceptance contracts: `OIC-021`. These obligations and the complete selected
  outcome are jointly required.

```text
PRODUCT_DECISION
Concern: web.record-display.feature-start-shortcuts
Scenario revision: 1
Choice: A / POPUP_ONLY_FEATURE_ROTATION
Rationale: 遺伝子を起点にした回転の入口を、開いた遺伝子に効き、鎖の向きを正しく扱うポップアップの 1 か所にまとめる。
Must preserve: sidebar の Display start 数値欄と Reset start。
May retire: sidebar の「Use selected feature 5′ end」と「Use selected feature midpoint」の 2 ボタン。遺伝子を起点にした回転はポップアップ（Rotate record using this feature）に一本化する。
Accepted residual risk: Ctrl/Shift で選択した遺伝子から、サイドバーで開始点を決める手段がなくなる。今後は遺伝子を開いて、ポップアップで操作する。
Owner: satoshikawato
Decision date: 2026-10-08
```

```json
{
  "concern": "web.record-display.feature-start-shortcuts",
  "scenarioRevision": 1,
  "choice": "A / POPUP_ONLY_FEATURE_ROTATION",
  "rationale": "遺伝子を起点にした回転の入口を、開いた遺伝子に効き、鎖の向きを正しく扱うポップアップの 1 か所にまとめる。",
  "mustPreserve": "sidebar の Display start 数値欄と Reset start。",
  "mayRetire": "sidebar の「Use selected feature 5′ end」と「Use selected feature midpoint」の 2 ボタン。遺伝子を起点にした回転はポップアップ（Rotate record using this feature）に一本化する。",
  "acceptedResidualRisk": "Ctrl/Shift で選択した遺伝子から、サイドバーで開始点を決める手段がなくなる。今後は遺伝子を開いて、ポップアップで操作する。",
  "owner": "satoshikawato",
  "decisionDate": "2026-10-08"
}
```

### PD-OI-088: Popup dialog busy until its choice commits

- Concern key: `web.popup-dialog-choice-application`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `B / DIALOG_FIRST_BUSY_CHOICE`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its nine-field JSON representation below.
- Decision source: the Owner replies of `2026-10-09` quoted verbatim in the
  Revision 35 entry above. The receipt and JSON below reproduce all nine
  fields without translation or additional terms. This record does not
  supersede or narrow another decision. The runtime that implements it merges
  with it (static Product Contract co-change).
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `b856d6da3f8f3c9b079284e6f542fde09b39d4ce19740aada5912053092264d3`.
- Acceptance contracts: `OIC-028`. These obligations and the complete selected
  outcome are jointly required.

```text
PRODUCT_DECISION
Concern: web.popup-dialog-choice-application
Scenario revision: 1
Choice: B / DIALOG_FIRST_BUSY_CHOICE
Rationale: 待ち時間を利用者が操作した場所（ダイアログ）に示し、適用中のもう一度のクリックが別の選択を始める競合をなくす。
Must preserve: ポップアップの色・Legend 名・Reset fill color の各ダイアログの選択肢と結果、1 回の選択が History の 1 ステップであること、Cancel は History に記録しないこと、live edit 適用中/失敗の通知（PD-OI-037）。
May retire: 色を選んでからダイアログが開くまでの待ち（追加する規則の準備は選択の後に行う）。選択の適用中に、そのダイアログの選択肢・Cancel・Escape・背景クリックで操作できること。
Accepted residual risk: Session 読み込み直後の最初の選択では、ダイアログが数秒（手元で約 5〜6 s）「Applying an edit…」のまま操作できず、その間は取り消せない。
Owner: satoshikawato
Decision date: 2026-10-09
```

```json
{
  "concern": "web.popup-dialog-choice-application",
  "scenarioRevision": 1,
  "choice": "B / DIALOG_FIRST_BUSY_CHOICE",
  "rationale": "待ち時間を利用者が操作した場所（ダイアログ）に示し、適用中のもう一度のクリックが別の選択を始める競合をなくす。",
  "mustPreserve": "ポップアップの色・Legend 名・Reset fill color の各ダイアログの選択肢と結果、1 回の選択が History の 1 ステップであること、Cancel は History に記録しないこと、live edit 適用中/失敗の通知（PD-OI-037）。",
  "mayRetire": "色を選んでからダイアログが開くまでの待ち（追加する規則の準備は選択の後に行う）。選択の適用中に、そのダイアログの選択肢・Cancel・Escape・背景クリックで操作できること。",
  "acceptedResidualRisk": "Session 読み込み直後の最初の選択では、ダイアログが数秒（手元で約 5〜6 s）「Applying an edit…」のまま操作できず、その間は取り消せない。",
  "owner": "satoshikawato",
  "decisionDate": "2026-10-09"
}
```

### PD-OI-089: Default colors a user edits win over the palette

- Concern key: `web.default-colors-over-palette`
- Scenario revision: `2`
- Supersedes: `PD-OI-089`, scenario revision `1`
  (`A / USER_DEFAULT_COLORS_OVER_PALETTE`), only by renaming one option of
  the Color/Stroke Change Scope dialog; every scenario-1 term stays in force.
- Status: `ACCEPTED`
- Selected outcomes: `A / USER_DEFAULT_COLORS_OVER_PALETTE` (scenario revision
  `1`) and the label change of scenario revision `2`, whose receipt names no
  choice code
- Normative outcome: exactly the two approved `PRODUCT_DECISION` receipts and
  their nine-field JSON representations below, read together. Where the
  scenario-1 receipt names **Apply to all** on a Legend row, it means the
  option now labelled `Legend item "<legend name>" (N)`.
- Decision source (scenario revision `1`): the Owner replies of `2026-10-09`
  quoted verbatim in the Revision 36 entry above. The receipt and JSON below
  reproduce all nine fields without translation or additional terms. This
  record does not supersede or narrow another decision. The runtime that
  implements it merges with it (static Product Contract co-change).
- Receipt SHA-256 (scenario revision `1`; UTF-8, excluding the final newline):
  `131d1bb0cacdcdc031387e6bff6c35672c4cb468a64220dba04dc70460e0ef78`.
- Decision source (scenario revision `2`): the Owner replies of `2026-10-11`
  quoted verbatim in the Revision 40 entry above. The second receipt and JSON
  below reproduce the confirmed fields without translation or additional
  terms. The runtime that shows the new name merges with it (static Product
  Contract co-change).
- Receipt SHA-256 (scenario revision `2`; UTF-8, excluding the final newline):
  `68eb7c40ac71b1363789062e85b82f7171603a9a026a331fc4019d862346c07d`.
- Acceptance contracts: `OIC-029`. These obligations and the complete selected
  outcome are jointly required.

```text
PRODUCT_DECISION
Concern: web.default-colors-over-palette
Scenario revision: 1
Choice: A / USER_DEFAULT_COLORS_OVER_PALETTE
Rationale: パレットは土台で、利用者が変えた既定の色（パレットの色と違う色）がその上に勝つ。CLI の -p と -d と同じ関係にし、パレットの切り替えや Reset で利用者の色が黙って消えたり、Apply to all で一部の feature だけが新しい色になって行が二つに分かれたりしないようにする。
Must preserve: Default colors の一覧での直接の編集は確認なしに反映されること、Palette instant preview とキューの表示の正確さ（PD-OI-037、OIC-024）、即時の編集と次の Generate の一致（PD-OI-066）、Specific color rule が描く行の編集は rule を変えること、「This feature only」と label・source label の範囲は rule を書くこと（PD-OI-069）、ダイアログの選択は History の 1 ステップで Cancel は記録せず適用中は busy（PD-OI-088）、Circular と Linear で別々の palette と色（PD-OI-086）、Session の形式、Generate と Run Info の CLI recipe が -p とパレットと違う色だけの -d を書くこと、Reset Settings（PD-OI-070）。
May retire: パレットを切り替えると利用者が変えた既定の色が確認なしに消える動作（キューから適用中のパレットに戻すと適用中の色に戻る動作を含む）。Default colors の Reset が確認なしに利用者の色を消す動作。Specific color rule が描かない Legend 行で Apply to all を選ぶと feature ごとの hash rule を書き、後から行に加わる feature がパレットの色のまま別の行になる動作。
Accepted residual risk: 利用者の色は保存せず導出するので、新しいパレットの色と同じになった色はそれ以降パレットに従う。パレットの定義が版で変わると、以前の既定値と同じ色も利用者の色として数える。Specific color rule と既定の色の feature が同じ色で一つの行にあるときと、その型の feature が一つだけのとき（Apply to all が出ない）は、今までどおり rule を書く。Apply to all は画面にない（非表示や他の record・Result の）その型の feature の色も変える。
Owner: satoshikawato
Decision date: 2026-10-09
```

```json
{
  "concern": "web.default-colors-over-palette",
  "scenarioRevision": 1,
  "choice": "A / USER_DEFAULT_COLORS_OVER_PALETTE",
  "rationale": "パレットは土台で、利用者が変えた既定の色（パレットの色と違う色）がその上に勝つ。CLI の -p と -d と同じ関係にし、パレットの切り替えや Reset で利用者の色が黙って消えたり、Apply to all で一部の feature だけが新しい色になって行が二つに分かれたりしないようにする。",
  "mustPreserve": "Default colors の一覧での直接の編集は確認なしに反映されること、Palette instant preview とキューの表示の正確さ（PD-OI-037、OIC-024）、即時の編集と次の Generate の一致（PD-OI-066）、Specific color rule が描く行の編集は rule を変えること、「This feature only」と label・source label の範囲は rule を書くこと（PD-OI-069）、ダイアログの選択は History の 1 ステップで Cancel は記録せず適用中は busy（PD-OI-088）、Circular と Linear で別々の palette と色（PD-OI-086）、Session の形式、Generate と Run Info の CLI recipe が -p とパレットと違う色だけの -d を書くこと、Reset Settings（PD-OI-070）。",
  "mayRetire": "パレットを切り替えると利用者が変えた既定の色が確認なしに消える動作（キューから適用中のパレットに戻すと適用中の色に戻る動作を含む）。Default colors の Reset が確認なしに利用者の色を消す動作。Specific color rule が描かない Legend 行で Apply to all を選ぶと feature ごとの hash rule を書き、後から行に加わる feature がパレットの色のまま別の行になる動作。",
  "acceptedResidualRisk": "利用者の色は保存せず導出するので、新しいパレットの色と同じになった色はそれ以降パレットに従う。パレットの定義が版で変わると、以前の既定値と同じ色も利用者の色として数える。Specific color rule と既定の色の feature が同じ色で一つの行にあるときと、その型の feature が一つだけのとき（Apply to all が出ない）は、今までどおり rule を書く。Apply to all は画面にない（非表示や他の record・Result の）その型の feature の色も変える。",
  "owner": "satoshikawato",
  "decisionDate": "2026-10-09"
}
```

Scenario revision `2` receipt, confirmed by `satoshikawato` on `2026-10-11`:

```text
PRODUCT_DECISION
Concern: web.default-colors-over-palette
Scenario revision: 2
Choice: Color/Stroke Change Scope option "Apply to all "<legend name>" (N)" on a Legend row -> "Legend item "<legend name>" (N)" under the group heading for features that share it.
Rationale: UX-A05 (one label per group, no repeated "Apply to all").
Must preserve: Options and effects unchanged; OIC-029's default-colour line ("Sets the <type> default color …" / "Replaces your <type> default color ■ #hex.") kept verbatim.
May retire: the old visible labels only.
Accepted residual risk: older docs/videos show the old labels.
Owner: satoshikawato
Decision date: 2026-10-11
```

```json
{
  "concern": "web.default-colors-over-palette",
  "scenarioRevision": 2,
  "choice": "Color/Stroke Change Scope option \"Apply to all \"<legend name>\" (N)\" on a Legend row -> \"Legend item \"<legend name>\" (N)\" under the group heading for features that share it.",
  "rationale": "UX-A05 (one label per group, no repeated \"Apply to all\").",
  "mustPreserve": "Options and effects unchanged; OIC-029's default-colour line (\"Sets the <type> default color …\" / \"Replaces your <type> default color ■ #hex.\") kept verbatim.",
  "mayRetire": "the old visible labels only.",
  "acceptedResidualRisk": "older docs/videos show the old labels.",
  "owner": "satoshikawato",
  "decisionDate": "2026-10-11"
}
```

### PD-OI-090: Project with drawings

- Concern key: `session.project-drawings`
- Scenario revision: `4`
- Supersedes: `PD-OI-090`, scenario revision `3` (`A / PROJECT-WITH-DRAWINGS`).
  Scenario revision `4` keeps the outcome and every term of scenario
  revision `3` except the Session versions of 〔Save と Load〕, 〔以前の Session
  の移行〕, May retire, and Accepted residual risk (the drawings format is
  Session version `49`; `48` is the 0.14.0 format, which moves to it without
  Python), and one added sentence in 〔Save と Load〕: a Session `46` or `47`
  that holds protein comparison results stored on a reversed record shows
  its stored figure at once, and Python moves those results to the search
  coordinates once at Load.
  Scenario revision `3` kept every term of scenario revision `2` except the
  first sentence of 〔Undo と削除〕: Undo/Redo also
  walks the project operations, which now include Reset all drawings and the
  edits that change only values that do not depend on a drawing, and turning
  Instant Preview on applies every drawing's queued palette in its one step.
  Scenario revision `2` kept every term of scenario revision `1` except
  two: the drawings format is Session version `48` (scenario revision `1`
  named `47`, which the Legend row keys of Revision 49 took), and a drawing switch also waits while a color rule table edit is
  being matched. Scenario revision `1` superseded `PD-OI-086`, scenario
  revision `1` (`A / ALL-DIAGRAM-SETTINGS-PER-MODE`, concern
  `web.mode.scoped-settings`), and removed that record from the active
  Contract.
- Status: `ACCEPTED`
- Selected outcome: `A / PROJECT-WITH-DRAWINGS`; at a drawing switch while a
  color rule table edit is being matched, `B / SWITCH-WAITS-FOR-RULE-COMMIT`
  of the OV-1100 question (Revision 62 entry); for the Undo/Redo of project
  operations, the answers `A` to Q-DR-5 and Q-DR-8 (Revision 63 entry); for
  the Load of a Session `46` or `47` with protein comparison results stored
  on a reversed record, `B-LOAD` of Q-O2 (Revision 64 entry). The receipt
  names the record's outcome only.
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its JSON representation below.
- Decision source: the Owner replies quoted verbatim in the Revision 37 entry
  above, for scenario revision `2` those quoted in the Revision 62 entry
  above, for scenario revision `3` those quoted in the Revision 63 entry
  above, and for scenario revision `4` those quoted in the Revision 64 entry
  above. The Owner approved the receipt text below as written. The receipt
  reproduces its fields without translation or additional terms. Dependent
  runtime requires this authority merged into its base; this record supplies
  no runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `dbc29cc437fd62eb2aad80cbbe7872f7ab701896b9a7ccf54609ca5b9c168f81`.
- Acceptance contracts: `OIC-004`, `OIC-005`, `OIC-006`, `OIC-013`,
  `OIC-027`, `OIC-030`. These obligations and the complete selected outcome
  are jointly required; their citation does not claim completed
  dependent-runtime checks.

```text
PRODUCT_DECISION
Concern: session.project-drawings
Scenario revision: 4
Supersedes: PD-OI-090, scenario revision 3
Choice: A / PROJECT-WITH-DRAWINGS
Rationale: 1 つの Session（project）に、共有の入力ファイルと複数の drawing（図）を持たせる。PD-OI-086 が Circular と Linear に分けた設定、凡例の編集、History、Result を、drawing ごとに持つ形へ広げる。ある図の編集や Generate が別の図の設定や Result を変えたり失敗させたりしないようにし（OV-80、OV-82、OV-104）、同じ mode の drawing も複数持てるようにして、同じゲノムの全体図と遺伝子クラスターの拡大図のように、同じ入力から複数の図を作れるようにする。Session の形式は Python が定め（読み、検証し、移行し、書く）、Web、CLI、Python API が同じ形式を使う。
Must preserve: 〔構造〕Session は共有の入力ファイルと 1 つ以上の drawing を持つ。drawing は作るときに mode（Circular または Linear）を決め、あとから変えない。同じ mode の drawing を複数持てる。drawing の ID は作るときの種類から決まり（circular、linear、linear-2 …）、変わらない。Session の title、アプリの設定（Auto Reflow、PNG DPI、palette の Instant Preview など）、LOSAT の実行方法とスレッド数は drawing によらず 1 つ（PD-OI-020、OIC-015、OIC-016）。〔drawing ごとに持つもの〕図の設定（Depth、GC content と GC skew、ラベル、軸と目盛り、フォント、線、トラック、凡例と title の位置、比較の閾値、LOSAT の検索の設定を含むすべての設定）、palette と色ルール（「この feature だけ」の色を含む）、qualifier priority、label の表と filter、注釈セット、凡例の編集（色、線、名前の変更、項目の削除と追加、順序）、feature と label の編集、Feature placement、描く record と region と表示（crop、開始位置、向き）、canvas の余白、画面にない設定の上書き（PD-OI-009）、Depth の割り当て、History、Result。drawing どうしで共有するのは、入力ファイルの中身と、結果を変えない計算の cache だけ。Show Depth は、その drawing に最初の Depth ファイルが入ると On、最後のファイルがなくなると Off になり、他の drawing は変わらない。〔drawing bar と切り替え〕画面上部の drawing bar が Circular/Linear の切り替えのボタンに替わる。drawing の切り替えは移動であり、どの drawing の値も書き換えず、他の drawing へ写さず、既定に戻さず、Undo の step にもならない（PD-OI-002、PD-OI-065）。Generate と Session の保存・読み込みの間は切り替えられない（PD-OI-045）。色ルールの表（Colors の SPECIFIC RULES）の編集が照合を待っている間も切り替えられず、照合が終わるとその編集はその drawing の History の 1 つの step になる。新しい Session は Circular の drawing 1 つで始まり、Linear は空の Linear の drawing を作る案内の tab から作る。ファイル、設定の変更、編集、Result のどれもない drawing は保存しない。tab は Result があるかどうか以外の導出された状態を常時表示しない（PD-OI-037）。bar は 390 px の幅でも操作でき、Preview と Editor の条件を保つ（PD-OI-038、OIC-025）。〔Result と失敗の分離〕Generate は表示中の drawing の Result だけを置き換え、他の drawing の Result はそのまま残る。戻ると、その drawing の Result がドラッグや編集を保って表示される。ある drawing の編集、Generate、失敗、Cancel、Undo/Redo は、他の drawing の設定、凡例、History、Result を変えない（OIPC-C07、PD-OI-016）。Generate、即時の編集、Session の保存と読み込み、export は、その drawing の値だけを使い（PD-OI-066、OIC-027）、export は表示中の drawing の現在の Result を出す。〔作成〕「+ New drawing」から次の 4 種類を選べる。（1）新規 Circular または新規 Linear：まっさら（ファイルなし、その mode の既定値、凡例と feature の編集なし。PD-OI-019、OIC-004）。（2）複製：ファイル、設定、凡例と feature の編集、Result を写し、History は空から。（3）スタイルをコピーして新規：設定とスタイルを写し、ファイル、凡例の編集、feature ごとの編集は写さない。別の mode へは両方の mode にある設定だけを写す。（4）選んだ領域から新規：領域は feature、record 全体、座標のいずれかで選ぶ。同じファイルの同じ record の領域を、元の座標で 1 始まり、両端を含む範囲として描く。複数の record を選んだときは、record ごとの領域を 1 つの Linear の drawing に並べる。元の drawing のその record の見た目、つまり設定、palette と色ルール、label の規則、凡例の編集、領域の中の feature の編集と注釈、その record の Depth を引き継ぐ。Auto の値は新しい drawing の長さで決まる。自分で入れた値は、新しい長さと mode での Auto が元と同じときだけ残し、違うときは Auto に戻して作成の画面に示し、戻さないことも選べる。アップロードした比較表は写さず、作成の画面に示す。元の drawing は変わらない。〔入力ファイル〕入力ファイルの中身は project に 1 回だけ保存し、同じ内容のファイルは 1 つにまとめる。どのファイルのどの record と region を描くかは drawing ごとに選び、別の drawing は project にあるファイルを再 upload なしで選べる。Linear の File の順序は行の順序と一体で drawing ごとに持つ（PD-OI-018）。ある drawing でファイルを外したり置き換えたりしても、他の drawing の選択は変わらない（次の「すべての drawing」を選んだときを除く）。別の drawing も使う入力ファイルを置き換えるときは、「この drawing だけ / このファイルを使うすべての drawing」をその都度選ばせ、選んだ drawing の record を更新して古い source の record を残さない（反応の owner は app/drawing-list.js、channel は port の applyInputReplacementToDrawings、投影の関数は projectInputReplacement）。どの drawing も使わなくなったファイルは保存しない。〔Undo と削除〕Undo/Redo は表示中の drawing の操作と project の操作を、新しい順に 1 操作 1 step で戻し、他の drawing の操作は戻さない。project の操作は、drawing の追加、複製、名前の変更、並べ替え、削除、Reset all drawings と、drawing によらず 1 つの値（Session の title、Auto Reflow、PNG DPI、palette の Instant Preview、LOSAT の実行方法とスレッド数、編集の入力欄に書きかけた値）だけを変える操作で、どの drawing を表示していても戻せる。drawing の値も一緒に変える操作（色ルールを追加して入力欄の書きかけが消えるなど）は、その drawing の操作として 1 step で戻る。Rich Feature Popup は今までどおり Undo/Redo では変わらない（PD-OI-100）。palette の Instant Preview を On にすると、すべての drawing で反映を待っている palette をその 1 step で反映し、その Undo で待っている状態に戻る。Result か編集のある drawing の削除は確認し、Undo で戻せる。最後の 1 つの drawing は削除できない。〔Reset〕Reset Settings は確認の dialog を出し、Reset this drawing / Reset all drawings / Cancel の 3 つのボタンから選ぶ。Reset this drawing は表示中の drawing の、Reset all drawings はすべての drawing の設定と編集を既定に戻す。編集（色、線、非表示、凡例、label）は表示中の図にすぐ反映し、Generate で反映する設定は次の Generate で変わる。表示していない drawing の図には、その drawing を次に表示したときに反映する。どちらでも入力ファイル、record の選択、ファイルの既定値、Depth の割り当て、Result は残る（PD-OI-070）。Cancel は何も変えず、履歴にも残さない。実行は 1 回の Undo で戻る。〔Save と Load〕Save Session は保存する全 drawing と、表示していた drawing を 1 つの整合した document として保存し、Load は全 drawing を一度に戻すか、何も変えない（変換に失敗したときも。PD-OI-045、PD-OI-046）。Load は各 drawing の保存された Result を Generate せずに示し、保存したときに表示していた drawing を開く。version 46、47、48 と 49 の Session の preview だけの Load は Python を起動しない（PD-OI-044 revision 4）。ただし、version 46 と 47 の Session で、逆向きにした record を含むタンパク質比較の保存された結果を持つものは、保存された図をすぐに示したうえで、読み込み時に Python が一度だけその結果を探索の座標に移す（PD-OI-044 revision 4）。〔既存の決定〕PD-OI-002 と PD-OI-004 の LOSATP の上限（Similarity と Collinear の比較 mode ごと）は、Linear の drawing の中で今までどおり。PD-OI-052 の装飾の delta の照合は drawing をまたがない（同じ mode と record の drawing の間も）。PD-OI-061 の衝突の選択、PD-OI-062 の batch の投影、PD-OI-063 の凡例の順序、PD-OI-064 の余白、PD-OI-065 の label の override（「表示の変化」に drawing の切り替えを含む）、PD-OI-027 と PD-OI-032 の record の向きと回転と provenance は、それぞれの drawing の中で今までどおり。〔以前の Session の移行〕以前の Session は、読み込むときに一度だけ project に移し、保存された値は表示していなかった mode のものも消さない（OIPC-C05、OIPC-C06、S00 の判断 1）。version 46、47 と 48（0.14.0 は 48）は Python を起動せずに移し、mode ごとの設定、凡例の編集、feature の編集、Result をそれぞれの mode の drawing に入れる。Result のない mode の drawing は、使った形跡（その mode の入力ファイル、その mode を表示したままの保存、既定値とも他方の mode の値とも違うその mode の設定）があるときだけ作る。version 27〜44 は Python が移し、まず次の規則で mode ごとに分けてから、同じ規則で drawing にする。両方の mode で共有されていた設定は作るすべての drawing に写す。mode ごとだった値（mode profile、片方の mode 専用の設定、Circular と Linear で分かれていた設定、mode 付きの行）はその mode の drawing に入れ、mode profile にない他方の mode の値は既定値にする。mode を持たずに保存された title と font は、保存された Result の drawing だけに入れる。Show Depth はその drawing に Depth ファイルがあるときだけ On にする。凡例の編集と、保存された Result に結びついた feature ごとの編集は、保存された Result の drawing だけに入れる。record を選んだ注釈は、その record を選んだ drawing だけに入れる。version 45 は開発版だけで、読まない。以前の Session の再生は今までどおりで、drawing が 1 つの project は以前の Session と同じ図を描く。PD-OI-030 の旧 alignment、PD-OI-073 の main の Session の比較表の変換、PD-OI-079 の version 40 より前の Session は Result のある drawing を Generate してから Save する規則は、今までどおり。〔CLI と Python API〕gbdraw render は project の 1 つまたはすべての drawing を描き（--drawing ID、--list_drawings）、複数を描くときは出力ファイルの名前に drawing の ID を付けて互いに上書きしない。drawing が 1 つの mode は今までどおり gbdraw circular / gbdraw linear の --session で再生でき、同じ mode の drawing が複数あるときはエラーで候補を示して --drawing で選ばせる。Python API は drawing を選んで、またはすべてを描く。
May retire: 図が Circular と Linear の 2 つに限られ、mode ごとに設定、凡例の編集、Result、入力ファイルを 1 組ずつ持つ形（PD-OI-086 と、version 46、47 と 48 の modes と otherModeResult。読むだけ）。画面上部の Circular/Linear の切り替えのボタン（drawing の bar に置き換える）。両方の mode が 1 つの Undo の履歴を共有し、mode の切り替えが Undo の step になる動作。Reset Settings が確認の window.confirm 1 つで、表示していない mode の設定も既定に戻す動作。同じ内容の入力ファイルを mode ごとに別々に保存する動作。version 27〜44 の Session を Python を起動せずに読み込む動作（PD-OI-044 revision 1 と PD-OI-086 の範囲）。
Accepted residual risk: 環状の record の原点をまたぐ領域は、最初の版では切り出せない（作成の画面で理由を示す）。Auto の大きさは record の長さの区分（50 kb など）で切り替わるため、区分の境をまたぐ領域では見た目が急に変わる。2 つの図を同じ設定にしたいときは、複製かスタイルのコピーを使うか、両方の drawing で設定する。以前の Session から作った他方の drawing に Result がなければ Generate が必要で、version 27〜44 の Session の凡例の編集は保存された Result の drawing だけに入る。version 27〜44 の Session で両方の mode に写した共有の値（Depth の最小・最大、window など）が、他方の mode のデータに合わないことがある。0.14.0 で他方の mode を表示も Generate もせず、ファイルも置かず、設定を一方の mode と同じにしていた Session では、その mode の drawing は作らない（スタイルをコピーして新規で作り直せる）。version 27〜44 の Session を開くときは、保存された図をすぐ示したうえで Python を一度起動して移し、その間（約 8 秒）は編集を待つ（新しい形式で保存し直すと起こらない）。新しい形式（version 49）で保存した Session は以前の版の gbdraw では読めない。drawing が増えると Session とメモリが大きくなる（Result と History は drawing ごと）。
Owner: satoshikawato
Decision date: 2026-10-18
```

```json
{
  "concern": "session.project-drawings",
  "scenarioRevision": 4,
  "supersedes": "PD-OI-090, scenario revision 3",
  "choice": "A / PROJECT-WITH-DRAWINGS",
  "rationale": "1 つの Session（project）に、共有の入力ファイルと複数の drawing（図）を持たせる。PD-OI-086 が Circular と Linear に分けた設定、凡例の編集、History、Result を、drawing ごとに持つ形へ広げる。ある図の編集や Generate が別の図の設定や Result を変えたり失敗させたりしないようにし（OV-80、OV-82、OV-104）、同じ mode の drawing も複数持てるようにして、同じゲノムの全体図と遺伝子クラスターの拡大図のように、同じ入力から複数の図を作れるようにする。Session の形式は Python が定め（読み、検証し、移行し、書く）、Web、CLI、Python API が同じ形式を使う。",
  "mustPreserve": "〔構造〕Session は共有の入力ファイルと 1 つ以上の drawing を持つ。drawing は作るときに mode（Circular または Linear）を決め、あとから変えない。同じ mode の drawing を複数持てる。drawing の ID は作るときの種類から決まり（circular、linear、linear-2 …）、変わらない。Session の title、アプリの設定（Auto Reflow、PNG DPI、palette の Instant Preview など）、LOSAT の実行方法とスレッド数は drawing によらず 1 つ（PD-OI-020、OIC-015、OIC-016）。〔drawing ごとに持つもの〕図の設定（Depth、GC content と GC skew、ラベル、軸と目盛り、フォント、線、トラック、凡例と title の位置、比較の閾値、LOSAT の検索の設定を含むすべての設定）、palette と色ルール（「この feature だけ」の色を含む）、qualifier priority、label の表と filter、注釈セット、凡例の編集（色、線、名前の変更、項目の削除と追加、順序）、feature と label の編集、Feature placement、描く record と region と表示（crop、開始位置、向き）、canvas の余白、画面にない設定の上書き（PD-OI-009）、Depth の割り当て、History、Result。drawing どうしで共有するのは、入力ファイルの中身と、結果を変えない計算の cache だけ。Show Depth は、その drawing に最初の Depth ファイルが入ると On、最後のファイルがなくなると Off になり、他の drawing は変わらない。〔drawing bar と切り替え〕画面上部の drawing bar が Circular/Linear の切り替えのボタンに替わる。drawing の切り替えは移動であり、どの drawing の値も書き換えず、他の drawing へ写さず、既定に戻さず、Undo の step にもならない（PD-OI-002、PD-OI-065）。Generate と Session の保存・読み込みの間は切り替えられない（PD-OI-045）。色ルールの表（Colors の SPECIFIC RULES）の編集が照合を待っている間も切り替えられず、照合が終わるとその編集はその drawing の History の 1 つの step になる。新しい Session は Circular の drawing 1 つで始まり、Linear は空の Linear の drawing を作る案内の tab から作る。ファイル、設定の変更、編集、Result のどれもない drawing は保存しない。tab は Result があるかどうか以外の導出された状態を常時表示しない（PD-OI-037）。bar は 390 px の幅でも操作でき、Preview と Editor の条件を保つ（PD-OI-038、OIC-025）。〔Result と失敗の分離〕Generate は表示中の drawing の Result だけを置き換え、他の drawing の Result はそのまま残る。戻ると、その drawing の Result がドラッグや編集を保って表示される。ある drawing の編集、Generate、失敗、Cancel、Undo/Redo は、他の drawing の設定、凡例、History、Result を変えない（OIPC-C07、PD-OI-016）。Generate、即時の編集、Session の保存と読み込み、export は、その drawing の値だけを使い（PD-OI-066、OIC-027）、export は表示中の drawing の現在の Result を出す。〔作成〕「+ New drawing」から次の 4 種類を選べる。（1）新規 Circular または新規 Linear：まっさら（ファイルなし、その mode の既定値、凡例と feature の編集なし。PD-OI-019、OIC-004）。（2）複製：ファイル、設定、凡例と feature の編集、Result を写し、History は空から。（3）スタイルをコピーして新規：設定とスタイルを写し、ファイル、凡例の編集、feature ごとの編集は写さない。別の mode へは両方の mode にある設定だけを写す。（4）選んだ領域から新規：領域は feature、record 全体、座標のいずれかで選ぶ。同じファイルの同じ record の領域を、元の座標で 1 始まり、両端を含む範囲として描く。複数の record を選んだときは、record ごとの領域を 1 つの Linear の drawing に並べる。元の drawing のその record の見た目、つまり設定、palette と色ルール、label の規則、凡例の編集、領域の中の feature の編集と注釈、その record の Depth を引き継ぐ。Auto の値は新しい drawing の長さで決まる。自分で入れた値は、新しい長さと mode での Auto が元と同じときだけ残し、違うときは Auto に戻して作成の画面に示し、戻さないことも選べる。アップロードした比較表は写さず、作成の画面に示す。元の drawing は変わらない。〔入力ファイル〕入力ファイルの中身は project に 1 回だけ保存し、同じ内容のファイルは 1 つにまとめる。どのファイルのどの record と region を描くかは drawing ごとに選び、別の drawing は project にあるファイルを再 upload なしで選べる。Linear の File の順序は行の順序と一体で drawing ごとに持つ（PD-OI-018）。ある drawing でファイルを外したり置き換えたりしても、他の drawing の選択は変わらない（次の「すべての drawing」を選んだときを除く）。別の drawing も使う入力ファイルを置き換えるときは、「この drawing だけ / このファイルを使うすべての drawing」をその都度選ばせ、選んだ drawing の record を更新して古い source の record を残さない（反応の owner は app/drawing-list.js、channel は port の applyInputReplacementToDrawings、投影の関数は projectInputReplacement）。どの drawing も使わなくなったファイルは保存しない。〔Undo と削除〕Undo/Redo は表示中の drawing の操作と project の操作を、新しい順に 1 操作 1 step で戻し、他の drawing の操作は戻さない。project の操作は、drawing の追加、複製、名前の変更、並べ替え、削除、Reset all drawings と、drawing によらず 1 つの値（Session の title、Auto Reflow、PNG DPI、palette の Instant Preview、LOSAT の実行方法とスレッド数、編集の入力欄に書きかけた値）だけを変える操作で、どの drawing を表示していても戻せる。drawing の値も一緒に変える操作（色ルールを追加して入力欄の書きかけが消えるなど）は、その drawing の操作として 1 step で戻る。Rich Feature Popup は今までどおり Undo/Redo では変わらない（PD-OI-100）。palette の Instant Preview を On にすると、すべての drawing で反映を待っている palette をその 1 step で反映し、その Undo で待っている状態に戻る。Result か編集のある drawing の削除は確認し、Undo で戻せる。最後の 1 つの drawing は削除できない。〔Reset〕Reset Settings は確認の dialog を出し、Reset this drawing / Reset all drawings / Cancel の 3 つのボタンから選ぶ。Reset this drawing は表示中の drawing の、Reset all drawings はすべての drawing の設定と編集を既定に戻す。編集（色、線、非表示、凡例、label）は表示中の図にすぐ反映し、Generate で反映する設定は次の Generate で変わる。表示していない drawing の図には、その drawing を次に表示したときに反映する。どちらでも入力ファイル、record の選択、ファイルの既定値、Depth の割り当て、Result は残る（PD-OI-070）。Cancel は何も変えず、履歴にも残さない。実行は 1 回の Undo で戻る。〔Save と Load〕Save Session は保存する全 drawing と、表示していた drawing を 1 つの整合した document として保存し、Load は全 drawing を一度に戻すか、何も変えない（変換に失敗したときも。PD-OI-045、PD-OI-046）。Load は各 drawing の保存された Result を Generate せずに示し、保存したときに表示していた drawing を開く。version 46、47、48 と 49 の Session の preview だけの Load は Python を起動しない（PD-OI-044 revision 4）。ただし、version 46 と 47 の Session で、逆向きにした record を含むタンパク質比較の保存された結果を持つものは、保存された図をすぐに示したうえで、読み込み時に Python が一度だけその結果を探索の座標に移す（PD-OI-044 revision 4）。〔既存の決定〕PD-OI-002 と PD-OI-004 の LOSATP の上限（Similarity と Collinear の比較 mode ごと）は、Linear の drawing の中で今までどおり。PD-OI-052 の装飾の delta の照合は drawing をまたがない（同じ mode と record の drawing の間も）。PD-OI-061 の衝突の選択、PD-OI-062 の batch の投影、PD-OI-063 の凡例の順序、PD-OI-064 の余白、PD-OI-065 の label の override（「表示の変化」に drawing の切り替えを含む）、PD-OI-027 と PD-OI-032 の record の向きと回転と provenance は、それぞれの drawing の中で今までどおり。〔以前の Session の移行〕以前の Session は、読み込むときに一度だけ project に移し、保存された値は表示していなかった mode のものも消さない（OIPC-C05、OIPC-C06、S00 の判断 1）。version 46、47 と 48（0.14.0 は 48）は Python を起動せずに移し、mode ごとの設定、凡例の編集、feature の編集、Result をそれぞれの mode の drawing に入れる。Result のない mode の drawing は、使った形跡（その mode の入力ファイル、その mode を表示したままの保存、既定値とも他方の mode の値とも違うその mode の設定）があるときだけ作る。version 27〜44 は Python が移し、まず次の規則で mode ごとに分けてから、同じ規則で drawing にする。両方の mode で共有されていた設定は作るすべての drawing に写す。mode ごとだった値（mode profile、片方の mode 専用の設定、Circular と Linear で分かれていた設定、mode 付きの行）はその mode の drawing に入れ、mode profile にない他方の mode の値は既定値にする。mode を持たずに保存された title と font は、保存された Result の drawing だけに入れる。Show Depth はその drawing に Depth ファイルがあるときだけ On にする。凡例の編集と、保存された Result に結びついた feature ごとの編集は、保存された Result の drawing だけに入れる。record を選んだ注釈は、その record を選んだ drawing だけに入れる。version 45 は開発版だけで、読まない。以前の Session の再生は今までどおりで、drawing が 1 つの project は以前の Session と同じ図を描く。PD-OI-030 の旧 alignment、PD-OI-073 の main の Session の比較表の変換、PD-OI-079 の version 40 より前の Session は Result のある drawing を Generate してから Save する規則は、今までどおり。〔CLI と Python API〕gbdraw render は project の 1 つまたはすべての drawing を描き（--drawing ID、--list_drawings）、複数を描くときは出力ファイルの名前に drawing の ID を付けて互いに上書きしない。drawing が 1 つの mode は今までどおり gbdraw circular / gbdraw linear の --session で再生でき、同じ mode の drawing が複数あるときはエラーで候補を示して --drawing で選ばせる。Python API は drawing を選んで、またはすべてを描く。",
  "mayRetire": "図が Circular と Linear の 2 つに限られ、mode ごとに設定、凡例の編集、Result、入力ファイルを 1 組ずつ持つ形（PD-OI-086 と、version 46、47 と 48 の modes と otherModeResult。読むだけ）。画面上部の Circular/Linear の切り替えのボタン（drawing の bar に置き換える）。両方の mode が 1 つの Undo の履歴を共有し、mode の切り替えが Undo の step になる動作。Reset Settings が確認の window.confirm 1 つで、表示していない mode の設定も既定に戻す動作。同じ内容の入力ファイルを mode ごとに別々に保存する動作。version 27〜44 の Session を Python を起動せずに読み込む動作（PD-OI-044 revision 1 と PD-OI-086 の範囲）。",
  "acceptedResidualRisk": "環状の record の原点をまたぐ領域は、最初の版では切り出せない（作成の画面で理由を示す）。Auto の大きさは record の長さの区分（50 kb など）で切り替わるため、区分の境をまたぐ領域では見た目が急に変わる。2 つの図を同じ設定にしたいときは、複製かスタイルのコピーを使うか、両方の drawing で設定する。以前の Session から作った他方の drawing に Result がなければ Generate が必要で、version 27〜44 の Session の凡例の編集は保存された Result の drawing だけに入る。version 27〜44 の Session で両方の mode に写した共有の値（Depth の最小・最大、window など）が、他方の mode のデータに合わないことがある。0.14.0 で他方の mode を表示も Generate もせず、ファイルも置かず、設定を一方の mode と同じにしていた Session では、その mode の drawing は作らない（スタイルをコピーして新規で作り直せる）。version 27〜44 の Session を開くときは、保存された図をすぐ示したうえで Python を一度起動して移し、その間（約 8 秒）は編集を待つ（新しい形式で保存し直すと起こらない）。新しい形式（version 49）で保存した Session は以前の版の gbdraw では読めない。drawing が増えると Session とメモリが大きくなる（Result と History は drawing ごと）。",
  "owner": "satoshikawato",
  "decisionDate": "2026-10-18"
}
```

### PD-OI-091: Per-record draw ON/OFF

- Concern key: `diagram-generation.record-draw-selection`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `A / RECORD-DRAW-ON-OFF`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its JSON representation below.
- Decision source: the Owner replies quoted verbatim in the Revision 38 entry
  above. The Owner approved the receipt text below as written. The receipt
  reproduces its fields without translation or additional terms. Dependent
  runtime requires this authority merged into its base, or merged together
  with it through the reviewed static Product Contract co-change route; this
  record supplies no runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `3eba5a6876af4b5ab575ab586291516e42bc298e53ef61b80a29093e2d5f49d6`.

```text
PRODUCT_DECISION
Concern: diagram-generation.record-draw-selection
Scenario revision: 1
Choice: A / RECORD-DRAW-ON-OFF
Rationale: 100 本以上のコンティグがあるドラフトアセンブリのように record の多いファイルから、見たい record だけを描けるようにする。0.13.0 の Web の Linear では record を 1 つ選んで描けたが、複数 record のファイルを record ごとの card に展開してから（#526）、record を外す方法がなくなった。CLI（--record_id、--records_table）と Python API はすでに一部の record だけを描ける。
Must preserve: 〔状態〕record ごとの「描く／描かない（ON/OFF）」を drawing ごとに 1 つ持つ（PD-OI-090）。既定はすべて ON（PD-OI-019）。〔2 つの入口〕Linear の record の card の中の Draw this record のチェックと、ファイルの record の一覧（Choose records…）のチェックの 2 か所から切り替え、どちらも同じ状態を読み書きする。Circular は Record の select の横の一覧から切り替える。一覧は record ID で検索でき、長さを示し、ファイルの順（元の順。既定）、record ID の順（昇順と降順。数字は数として比べ、contig_2 は contig_10 より前）、長さの順（長い順と短い順）に並べ替えられる。並べ替えは一覧の表示だけで、描く順は変えない。表示中の行をまとめて ON または OFF にできる。変更はすぐ反映し、1 つの操作は 1 回の Undo で戻る。〔描画〕OFF の record は Generate に送らず、図、比較のペア（All adjacent pairs を含む）、行の配置、凡例、注釈から外れる。比較の画面には ON の record だけを示し、Record Layout と Circular の Record Order には OFF の record を灰色で示す。OFF の card は閉じた状態でも灰色と OFF の印で分かる。ON に戻すと、ファイルの中の元の位置（Linear の行、Circular の grid の位置）に戻る。Circular の Multi-Record Canvas と record ごとの別々の図は ON の record だけを描く。Record の select で 1 つの record を選んだときは、ON/OFF によらずその record を描く。ON の record が 1 つでも、All records は別々の図のままで、自動で 1 つの図に切り替えない（PD-OI-044）。ON/OFF を変えると、次の Generate まで比較の結果と Similarity alignment の計画は無効になる。〔LOSAT〕E-value の database は今までどおり subject 側の source ファイル全体で（PD-OI-018）、record を ON/OFF しても他の record の E-value は変わらない。ペアは ON の record どうしだけを作る。〔設定の保持と削除〕OFF の record の card の設定（Organism、Subtitle、Region、Reverse complement、LOSAT Gencode、開始位置）、feature の編集、注釈、比較のペア、行は残り、Session の保存と読み込みでも残る。Delete settings（OFF の card と一覧にある）は、その record の card の設定を既定に戻し、その record の feature の編集と Feature placement、その record を選んだ注釈、その record を含む比較のペアも消す。1 回の Undo で戻る。〔最後の record〕ファイルの最後の ON の record を OFF にしようとすると、Remove the whole file? と聞く。Remove File はそのファイルを外し（ほかにファイルがなければ空にする。Circular は GenBank、または GFF3 と FASTA の組を外す）、Cancel は何も変えない（まとめて OFF にする操作なら、その操作全体を取り消す）。〔upload〕record が 20 を超えるファイルを upload すると、すべて ON のまま一覧を自動で開く。複数のファイルなら 1 つずつ開く。Session の読み込み、取り消し、Undo では開かない。〔Reset と Session〕Reset Settings は ON/OFF を残す（PD-OI-070）。Session は drawing ごとに OFF の record を保存し、その項目のない Session（以前の Web、CLI、Gallery の Session）はすべて ON として読む。Session の version は 46 のまま。CLI の --session の再生は ON の record だけを描き、保存し直しても OFF の record とその設定を残す。あるファイルの record がすべて OFF の Session はそのまま読み込み、そのファイルは描かない。すべての record が OFF なら、Generate は Turn on at least one record to draw. と示す。
May retire: 複数 record のファイルの全 record を必ず描く Web の Linear の動作（#526 以降）。Multi-Record Canvas と record ごとの別々の図が全 record を必ず描く動作。PD-OI-010 の受け入れた残余リスク（Arbitrary multi-record subset editing remains unavailable）と、PD-OI-044 revision 2 の May retire が record の一部を選ぶ編集を除外していること。
Accepted residual risk: record が数百あると Linear の card が数百並び、画面が重くなることがある（#526 から変わらない）。OFF の record の設定は Session に残るので、Session が大きくなる（Delete settings で消せる）。CLI には ON/OFF の option はなく、--record_id と --records_table で同じ図を描く。
Owner: satoshikawato
Decision date: 2026-10-10
```

```json
{
  "concern": "diagram-generation.record-draw-selection",
  "scenarioRevision": 1,
  "choice": "A / RECORD-DRAW-ON-OFF",
  "rationale": "100 本以上のコンティグがあるドラフトアセンブリのように record の多いファイルから、見たい record だけを描けるようにする。0.13.0 の Web の Linear では record を 1 つ選んで描けたが、複数 record のファイルを record ごとの card に展開してから（#526）、record を外す方法がなくなった。CLI（--record_id、--records_table）と Python API はすでに一部の record だけを描ける。",
  "mustPreserve": "〔状態〕record ごとの「描く／描かない（ON/OFF）」を drawing ごとに 1 つ持つ（PD-OI-090）。既定はすべて ON（PD-OI-019）。〔2 つの入口〕Linear の record の card の中の Draw this record のチェックと、ファイルの record の一覧（Choose records…）のチェックの 2 か所から切り替え、どちらも同じ状態を読み書きする。Circular は Record の select の横の一覧から切り替える。一覧は record ID で検索でき、長さを示し、ファイルの順（元の順。既定）、record ID の順（昇順と降順。数字は数として比べ、contig_2 は contig_10 より前）、長さの順（長い順と短い順）に並べ替えられる。並べ替えは一覧の表示だけで、描く順は変えない。表示中の行をまとめて ON または OFF にできる。変更はすぐ反映し、1 つの操作は 1 回の Undo で戻る。〔描画〕OFF の record は Generate に送らず、図、比較のペア（All adjacent pairs を含む）、行の配置、凡例、注釈から外れる。比較の画面には ON の record だけを示し、Record Layout と Circular の Record Order には OFF の record を灰色で示す。OFF の card は閉じた状態でも灰色と OFF の印で分かる。ON に戻すと、ファイルの中の元の位置（Linear の行、Circular の grid の位置）に戻る。Circular の Multi-Record Canvas と record ごとの別々の図は ON の record だけを描く。Record の select で 1 つの record を選んだときは、ON/OFF によらずその record を描く。ON の record が 1 つでも、All records は別々の図のままで、自動で 1 つの図に切り替えない（PD-OI-044）。ON/OFF を変えると、次の Generate まで比較の結果と Similarity alignment の計画は無効になる。〔LOSAT〕E-value の database は今までどおり subject 側の source ファイル全体で（PD-OI-018）、record を ON/OFF しても他の record の E-value は変わらない。ペアは ON の record どうしだけを作る。〔設定の保持と削除〕OFF の record の card の設定（Organism、Subtitle、Region、Reverse complement、LOSAT Gencode、開始位置）、feature の編集、注釈、比較のペア、行は残り、Session の保存と読み込みでも残る。Delete settings（OFF の card と一覧にある）は、その record の card の設定を既定に戻し、その record の feature の編集と Feature placement、その record を選んだ注釈、その record を含む比較のペアも消す。1 回の Undo で戻る。〔最後の record〕ファイルの最後の ON の record を OFF にしようとすると、Remove the whole file? と聞く。Remove File はそのファイルを外し（ほかにファイルがなければ空にする。Circular は GenBank、または GFF3 と FASTA の組を外す）、Cancel は何も変えない（まとめて OFF にする操作なら、その操作全体を取り消す）。〔upload〕record が 20 を超えるファイルを upload すると、すべて ON のまま一覧を自動で開く。複数のファイルなら 1 つずつ開く。Session の読み込み、取り消し、Undo では開かない。〔Reset と Session〕Reset Settings は ON/OFF を残す（PD-OI-070）。Session は drawing ごとに OFF の record を保存し、その項目のない Session（以前の Web、CLI、Gallery の Session）はすべて ON として読む。Session の version は 46 のまま。CLI の --session の再生は ON の record だけを描き、保存し直しても OFF の record とその設定を残す。あるファイルの record がすべて OFF の Session はそのまま読み込み、そのファイルは描かない。すべての record が OFF なら、Generate は Turn on at least one record to draw. と示す。",
  "mayRetire": "複数 record のファイルの全 record を必ず描く Web の Linear の動作（#526 以降）。Multi-Record Canvas と record ごとの別々の図が全 record を必ず描く動作。PD-OI-010 の受け入れた残余リスク（Arbitrary multi-record subset editing remains unavailable）と、PD-OI-044 revision 2 の May retire が record の一部を選ぶ編集を除外していること。",
  "acceptedResidualRisk": "record が数百あると Linear の card が数百並び、画面が重くなることがある（#526 から変わらない）。OFF の record の設定は Session に残るので、Session が大きくなる（Delete settings で消せる）。CLI には ON/OFF の option はなく、--record_id と --records_table で同じ図を描く。",
  "owner": "satoshikawato",
  "decisionDate": "2026-10-10"
}
```

### PD-OI-092: Settings pane section index

- Concern key: `web.settings-pane.section-index`
- Scenario revision: `2`
- Supersedes: `PD-OI-092`, scenario revision `1` (`A / PINNED-ICON-SECTION-INDEX`).
- Status: `ACCEPTED`
- Selected outcome: `A / PINNED-ICON-SECTION-INDEX`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its JSON representation below.
- Decision source: the Owner reply of `2026-10-13` quoted verbatim in the
  Revision 51 entry above. The Owner approved the receipt text below as
  written. The receipt reproduces its fields without translation or
  additional terms. Scenario revision `2` keeps the outcome and adds the Font section of the Revision 51 entry to the index. The runtime, tests, and documentation that
  implement it merge with it (static Product Contract co-change); this record
  supplies no runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `0fba5433eeb10e8c606b967f142c4a33b938d23db7fc3d6dc7fde5396e636c92`.

```text
PRODUCT_DECISION
Concern: web.settings-pane.section-index
Scenario revision: 2
Choice: A / PINNED-ICON-SECTION-INDEX
Rationale: 設定ペインは長く、下の section までのスクロールが長い（Linear で record が 5 つのとき、Labels は 3,909 px、Comparison は 1,972 px 下）。今の縦並びのまま、どの section にも 1 回のクリックで届くようにする。図の文字の font の設定は図全体にかかるので、section を 1 つ加えて置く。Titles & Record Labels の Default font size は record label の大きさだけなので、その横に置くと図全体の設定に見えない。
Must preserve: 〔目次〕設定ペインの一番上に、設定をスクロールしても動かない目次の行を置く。目次は section ごとのアイコンのボタンで、card と同じ順に並ぶ。Circular は Input Genomes、Region Annotations、Pairwise Comparisons、Layout、Font、Titles & Record Labels、Legend、Colors、Features、Labels、Axis & Scale、Dinucleotide content/skew の 12 個。Linear は Input Genomes、Comparison、Advanced comparison and layout、Region Annotations、Layout、Font、Titles & Record Labels、Legend、Colors、Features、Labels、Axis & Scale、Dinucleotide content/skew の 13 個。Basic（Input Genomes のすぐ下）と About & Citation は目次に入れない。〔押したとき〕その section が閉じていれば開き、設定の表示範囲の上端までスクロールし（設定の終わり近くの section は、スクロールできるところまで）、キーボードのフォーカスをその section の見出しに移す。目次は section を閉じない。〔名前と印〕ボタンにマウスを載せたとき、またはキーボードでフォーカスしたときに、section の名前がボタンの下に浮かび上がる。ボタンのアクセシブル名は「Go to ＜section の名前＞」。表示範囲の上端にある section のボタンに印が付く（押した直後は、次にスクロールするまで押したボタンに印が付く）。〔Font section〕図の文字の font の設定（Font）は、Layout と Titles & Record Labels のあいだの新しい section「Font」に置く。目次のアイコンは Aa。〔変えないもの〕Font section を加えること以外の card の順と中身、既定で開いている section、Generate Diagram が設定の最後の Tab stop で下に固定されていること、設定ペインの幅のドラッグ変更。目次の操作は設定の値、Session、Undo の履歴、図を変えない。
May retire: none
Accepted residual risk: 設定の表示範囲が 44 px 低くなる。Tab stop が Circular で 12、Linear で 13 増える（設定の最初の入力まで Tab が長くなる）。ボタンは 22 × 26 px と小さい。設定ペインを狭めると目次は 2 行になる。スマートフォンの幅では目次は固定されず、設定の先頭に置かれる。設定ペインの上端が写っている Gallery と docs の capture は、目次のアイコンと section が 1 つ増えるので撮り直すまで古い。設定ペインが 1 section ぶん長くなる。
Owner: satoshikawato
Decision date: 2026-10-13
```

```json
{
  "concern": "web.settings-pane.section-index",
  "scenarioRevision": 2,
  "choice": "A / PINNED-ICON-SECTION-INDEX",
  "rationale": "設定ペインは長く、下の section までのスクロールが長い（Linear で record が 5 つのとき、Labels は 3,909 px、Comparison は 1,972 px 下）。今の縦並びのまま、どの section にも 1 回のクリックで届くようにする。図の文字の font の設定は図全体にかかるので、section を 1 つ加えて置く。Titles & Record Labels の Default font size は record label の大きさだけなので、その横に置くと図全体の設定に見えない。",
  "mustPreserve": "〔目次〕設定ペインの一番上に、設定をスクロールしても動かない目次の行を置く。目次は section ごとのアイコンのボタンで、card と同じ順に並ぶ。Circular は Input Genomes、Region Annotations、Pairwise Comparisons、Layout、Font、Titles & Record Labels、Legend、Colors、Features、Labels、Axis & Scale、Dinucleotide content/skew の 12 個。Linear は Input Genomes、Comparison、Advanced comparison and layout、Region Annotations、Layout、Font、Titles & Record Labels、Legend、Colors、Features、Labels、Axis & Scale、Dinucleotide content/skew の 13 個。Basic（Input Genomes のすぐ下）と About & Citation は目次に入れない。〔押したとき〕その section が閉じていれば開き、設定の表示範囲の上端までスクロールし（設定の終わり近くの section は、スクロールできるところまで）、キーボードのフォーカスをその section の見出しに移す。目次は section を閉じない。〔名前と印〕ボタンにマウスを載せたとき、またはキーボードでフォーカスしたときに、section の名前がボタンの下に浮かび上がる。ボタンのアクセシブル名は「Go to ＜section の名前＞」。表示範囲の上端にある section のボタンに印が付く（押した直後は、次にスクロールするまで押したボタンに印が付く）。〔Font section〕図の文字の font の設定（Font）は、Layout と Titles & Record Labels のあいだの新しい section「Font」に置く。目次のアイコンは Aa。〔変えないもの〕Font section を加えること以外の card の順と中身、既定で開いている section、Generate Diagram が設定の最後の Tab stop で下に固定されていること、設定ペインの幅のドラッグ変更。目次の操作は設定の値、Session、Undo の履歴、図を変えない。",
  "mayRetire": "none",
  "acceptedResidualRisk": "設定の表示範囲が 44 px 低くなる。Tab stop が Circular で 12、Linear で 13 増える（設定の最初の入力まで Tab が長くなる）。ボタンは 22 × 26 px と小さい。設定ペインを狭めると目次は 2 行になる。スマートフォンの幅では目次は固定されず、設定の先頭に置かれる。設定ペインの上端が写っている Gallery と docs の capture は、目次のアイコンと section が 1 つ増えるので撮り直すまで古い。設定ペインが 1 section ぶん長くなる。",
  "owner": "satoshikawato",
  "decisionDate": "2026-10-13"
}
```

### PD-OI-093: Settings pane display switch and icon rail

- Concern key: `web.settings-pane.display`
- Scenario revision: `2`
- Supersedes: `PD-OI-093`, scenario revision `1` (`A / STACKED-DEFAULT-WITH-ICON-RAIL-SWITCH`).
- Status: `ACCEPTED`
- Selected outcome: `A / STACKED-DEFAULT-WITH-ICON-RAIL-SWITCH`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its JSON representation below.
- Decision source: the Owner reply of `2026-10-13` quoted verbatim in the
  Revision 51 entry above. The Owner approved the receipt text below as
  written. The receipt reproduces its fields without translation or
  additional terms. Scenario revision `2` keeps the outcome and adds the Font item of the Revision 51 entry to the icon rail. The runtime, tests, and documentation that
  implement it merge with it (static Product Contract co-change); this record
  supplies no runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `81bf15fcf9788503b438545f62256e6eb72f8e21f9e0f1213245ea3f858316c0`.

```text
PRODUCT_DECISION
Concern: web.settings-pane.display
Scenario revision: 2
Choice: A / STACKED-DEFAULT-WITH-ICON-RAIL-SWITCH
Rationale: 設定ペインが図の幅を取る（1600 px の窓で 320 px）。図を広く見たいときに、設定をアイコンの列だけ、または 1 つの section だけに絞れるようにする。今の表示に慣れた人のために、既定は今の表示のままにする。図の文字の font の設定のために加えた Font section も、目次と同じ順で列に置く。
Must preserve: 〔2 つの表示〕設定ペインは「今の表示（縦に並んだ card と目次の行）」か「アイコンの列の表示（アイコンと短い名前の縦の列と、その横に開いた 1 つの section）」のどちらか。既定は今の表示。目次の行の右端のボタン（Collapse settings to icons）でアイコンの列の表示に、列の一番上のボタン（Expand settings）で今の表示に切り替える。設定の値と section の中身はどちらの表示でも同じ 1 つのもので、切り替えても入力途中の値と section の開閉は変わらない。〔アイコンの列〕列は 60 px、横に開く section は 320 px。項目は上から Input（Input Genomes と Basic。Linear では Comparison を除く）、Linear では Compare（Comparison だけ）と Advanced、Regions、Circular では Compare（Pairwise Comparisons）、Layout、Font、Titles、Legend、Colors、Features、Labels、Axis、GC / skew、一番下に About。項目を押すと、その section を横に開き（閉じていれば開く）、キーボードのフォーカスを見出しに移し、その項目に印を付ける。開いている項目をもう一度押しても開いたまま。Close section のボタンで列だけにし、フォーカスは開いていた項目に戻る。項目にマウスを載せたとき、またはキーボードでフォーカスしたときに、section の名前が項目の右に出る。項目のアクセシブル名は目次と同じ「Go to ＜section の名前＞」。section 自身の開閉は、横に開いた状態でも使える。窓が低くて項目が収まらないときは列がスクロールし、細いスクロールバーが出る。〔切り替えたとき〕アイコンの列の表示にすると、目次で印の付いていた section が横に開く。今の表示に戻すと、開いていた section が設定の一番上に来る（列だけだったときは、スクロール位置はそのまま）。フォーカスは新しい表示の切り替えボタンに移る。Circular と Linear を切り替えたとき、開いていた section と同じ名前の section が新しいモードにあればそれを開いたままにし（Axis & Scale）、無ければ（Compare、Advanced）Input を開く。列だけのときは列だけのまま。〔設定へ案内する操作〕比較のエラーの原因を示す操作、比較の問題で Generate が止まったときの案内、レコードのポップアップから設定の欄へ移る操作、Circular の canvas の設定、Record Labels の表示の設定、色ルールのエラーの Edit pattern は、アイコンの列の表示ではその section を横に開いてから欄にフォーカスする。〔Generate〕Generate Diagram は 1 つのボタンで、どの表示でも設定の最後の Tab stop で、下に固定する。section を開いているときは列と section を合わせた幅、列だけのときは 60 px の正方形（アクセシブル名は Generate Diagram のまま）。処理中はローダーだけを表示する。Circular の layout の警告は、列だけのときは Generate の上の印（Layout warning という名前のボタン）になり、その文は Generate Diagram の説明として読み上げられ、印にマウスを載せる、フォーカスする、押す、タップすると表示される（印を押しても Generate は始まらない）。〔記憶〕選んだ表示（今の表示か、アイコンの列の表示か）だけをブラウザの保存領域に書き、次に開いたときにその表示にする（アイコンの列の表示で開いたときは Input を開く）。開いていた section、設定ペインの幅、設定、ファイルは書かない。Session の保存と読み込み、example の読み込み、Reset Settings は表示を変えない。表示の切り替え、項目、Close section は Undo の履歴に入らない。保存領域が使えないときは既定の表示で動き、メッセージは出さない。〔スマートフォンの幅〕幅 768 px 以下では常に今の表示で、切り替えボタンは出ない。覚えている選択は変えない。〔変えないもの〕今の表示での目次の項目、順、名前、押したときの動き、印、浮かび上がる名前。card の順と中身、既定で開いている section、設定の値と既定値、Session、図。今の表示での設定ペインの幅のドラッグ変更。
May retire: none
Accepted residual risk: Web アプリが初めてブラウザの保存領域に書く（表示の選択 1 つだけ）。選択はブラウザごと、アドレスごとで、gbdraw gui は起動のたびにポートが変わるため、ローカルではサーバーを止めると選択は消える。アイコンの列の表示で section を開くと図の幅は今より 62 px 狭い（1600 px の窓で 1,228 → 1,166 px。section を閉じると 1,494 px）。アイコンの列の表示では設定ペインの幅をドラッグで変えられない。今の表示の Tab stop が 1 つ増える（切り替えボタン）。Linear の目次のボタンは、切り替えボタンの場所を空けるために幅が 22 px から 21.2 px（スクロールバーが常に出るブラウザでは 20.5 px）になる。列の短い名前のうち Regions、Compare、GC / skew はアクセシブル名に含まれない。列だけのとき、layout の警告は印だけになり、文は印を開くまで見えない。高さ 800 px 以下の窓では列の項目が小さくなる（1280 × 720 で高さ 34〜39 px、名前は 9 px）。さらに低い窓では列がスクロールする（1366 × 650 の Linear で section を開いたとき）。設定へ案内する操作は、今の表示でも 1 つの規則になる（欄がすでに見えていれば動かさず、見えていなければ設定の中央に出す。これまでは操作ごとに違った）。横に開いた section をスクロールすると、右端の部品が Close section のボタンの下を通る。Review alignment を開いたまま表示を切り替えるか section を開閉すると、図は自動では合わせ直されない（窓の幅を変えたときの今の動きと同じ）。幅が約 1280 px の窓で section を開くと、Review alignment のパネルは図の下に移る。設定ペインが写っている Gallery と docs の capture は、撮り直すまで古い。列の項目が 1 つ（Font）増えるので、低い窓では列がスクロールしやすくなる。
Owner: satoshikawato
Decision date: 2026-10-13
```

```json
{
  "concern": "web.settings-pane.display",
  "scenarioRevision": 2,
  "choice": "A / STACKED-DEFAULT-WITH-ICON-RAIL-SWITCH",
  "rationale": "設定ペインが図の幅を取る（1600 px の窓で 320 px）。図を広く見たいときに、設定をアイコンの列だけ、または 1 つの section だけに絞れるようにする。今の表示に慣れた人のために、既定は今の表示のままにする。図の文字の font の設定のために加えた Font section も、目次と同じ順で列に置く。",
  "mustPreserve": "〔2 つの表示〕設定ペインは「今の表示（縦に並んだ card と目次の行）」か「アイコンの列の表示（アイコンと短い名前の縦の列と、その横に開いた 1 つの section）」のどちらか。既定は今の表示。目次の行の右端のボタン（Collapse settings to icons）でアイコンの列の表示に、列の一番上のボタン（Expand settings）で今の表示に切り替える。設定の値と section の中身はどちらの表示でも同じ 1 つのもので、切り替えても入力途中の値と section の開閉は変わらない。〔アイコンの列〕列は 60 px、横に開く section は 320 px。項目は上から Input（Input Genomes と Basic。Linear では Comparison を除く）、Linear では Compare（Comparison だけ）と Advanced、Regions、Circular では Compare（Pairwise Comparisons）、Layout、Font、Titles、Legend、Colors、Features、Labels、Axis、GC / skew、一番下に About。項目を押すと、その section を横に開き（閉じていれば開く）、キーボードのフォーカスを見出しに移し、その項目に印を付ける。開いている項目をもう一度押しても開いたまま。Close section のボタンで列だけにし、フォーカスは開いていた項目に戻る。項目にマウスを載せたとき、またはキーボードでフォーカスしたときに、section の名前が項目の右に出る。項目のアクセシブル名は目次と同じ「Go to ＜section の名前＞」。section 自身の開閉は、横に開いた状態でも使える。窓が低くて項目が収まらないときは列がスクロールし、細いスクロールバーが出る。〔切り替えたとき〕アイコンの列の表示にすると、目次で印の付いていた section が横に開く。今の表示に戻すと、開いていた section が設定の一番上に来る（列だけだったときは、スクロール位置はそのまま）。フォーカスは新しい表示の切り替えボタンに移る。Circular と Linear を切り替えたとき、開いていた section と同じ名前の section が新しいモードにあればそれを開いたままにし（Axis & Scale）、無ければ（Compare、Advanced）Input を開く。列だけのときは列だけのまま。〔設定へ案内する操作〕比較のエラーの原因を示す操作、比較の問題で Generate が止まったときの案内、レコードのポップアップから設定の欄へ移る操作、Circular の canvas の設定、Record Labels の表示の設定、色ルールのエラーの Edit pattern は、アイコンの列の表示ではその section を横に開いてから欄にフォーカスする。〔Generate〕Generate Diagram は 1 つのボタンで、どの表示でも設定の最後の Tab stop で、下に固定する。section を開いているときは列と section を合わせた幅、列だけのときは 60 px の正方形（アクセシブル名は Generate Diagram のまま）。処理中はローダーだけを表示する。Circular の layout の警告は、列だけのときは Generate の上の印（Layout warning という名前のボタン）になり、その文は Generate Diagram の説明として読み上げられ、印にマウスを載せる、フォーカスする、押す、タップすると表示される（印を押しても Generate は始まらない）。〔記憶〕選んだ表示（今の表示か、アイコンの列の表示か）だけをブラウザの保存領域に書き、次に開いたときにその表示にする（アイコンの列の表示で開いたときは Input を開く）。開いていた section、設定ペインの幅、設定、ファイルは書かない。Session の保存と読み込み、example の読み込み、Reset Settings は表示を変えない。表示の切り替え、項目、Close section は Undo の履歴に入らない。保存領域が使えないときは既定の表示で動き、メッセージは出さない。〔スマートフォンの幅〕幅 768 px 以下では常に今の表示で、切り替えボタンは出ない。覚えている選択は変えない。〔変えないもの〕今の表示での目次の項目、順、名前、押したときの動き、印、浮かび上がる名前。card の順と中身、既定で開いている section、設定の値と既定値、Session、図。今の表示での設定ペインの幅のドラッグ変更。",
  "mayRetire": "none",
  "acceptedResidualRisk": "Web アプリが初めてブラウザの保存領域に書く（表示の選択 1 つだけ）。選択はブラウザごと、アドレスごとで、gbdraw gui は起動のたびにポートが変わるため、ローカルではサーバーを止めると選択は消える。アイコンの列の表示で section を開くと図の幅は今より 62 px 狭い（1600 px の窓で 1,228 → 1,166 px。section を閉じると 1,494 px）。アイコンの列の表示では設定ペインの幅をドラッグで変えられない。今の表示の Tab stop が 1 つ増える（切り替えボタン）。Linear の目次のボタンは、切り替えボタンの場所を空けるために幅が 22 px から 21.2 px（スクロールバーが常に出るブラウザでは 20.5 px）になる。列の短い名前のうち Regions、Compare、GC / skew はアクセシブル名に含まれない。列だけのとき、layout の警告は印だけになり、文は印を開くまで見えない。高さ 800 px 以下の窓では列の項目が小さくなる（1280 × 720 で高さ 34〜39 px、名前は 9 px）。さらに低い窓では列がスクロールする（1366 × 650 の Linear で section を開いたとき）。設定へ案内する操作は、今の表示でも 1 つの規則になる（欄がすでに見えていれば動かさず、見えていなければ設定の中央に出す。これまでは操作ごとに違った）。横に開いた section をスクロールすると、右端の部品が Close section のボタンの下を通る。Review alignment を開いたまま表示を切り替えるか section を開閉すると、図は自動では合わせ直されない（窓の幅を変えたときの今の動きと同じ）。幅が約 1280 px の窓で section を開くと、Review alignment のパネルは図の下に移る。設定ペインが写っている Gallery と docs の capture は、撮り直すまで古い。列の項目が 1 つ（Font）増えるので、低い窓では列がスクロールしやすくなる。",
  "owner": "satoshikawato",
  "decisionDate": "2026-10-13"
}
```

### PD-OI-094: Linear File cards fold to one line

- Concern key: `web.linear.file-card-fold`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `A / FOLDABLE-FILE-CARD, FOLDED-AFTER-LOAD`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its JSON representation below.
- Decision source: the Owner replies quoted verbatim in the Revision 46 entry
  above. The Owner approved the receipt text below as written. The receipt
  reproduces its fields without translation or additional terms. Dependent
  runtime requires this authority merged into its base, or merged together
  with it through the reviewed static Product Contract co-change route; this
  record supplies no runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `633f30bac2795b4df63eca6f2fc82fedf3a199af14f9eaf4bd7f94bde8f93026`.

```text
PRODUCT_DECISION
Concern: web.linear.file-card-fold
Scenario revision: 1
Choice: A / FOLDABLE-FILE-CARD, FOLDED-AFTER-LOAD
Rationale: Linear の File カードは 1 枚が長く、ファイルが増えると Comparison や下の section までのスクロールが長くなる（ファイル 5 つのとき、ファイルの一覧は 1,778 px、Comparison は 1,972 px 下）。カードを 1 行に畳めるようにして、ファイルの一覧を 1 画面で見渡せるようにする（畳むと一覧は 218 px、Comparison は 412 px 下）。
Must preserve: 〔畳む・開く〕各 File カードの見出しの行を押すとカードを 1 行に畳み、もう一度押すと開く。上下の移動ボタンを押したときは畳まれない。畳んだ行には File の番号、ファイル名、上下の移動ボタンが出る。〔最初の状態〕Session または example を読み込んだとき、ファイルがそろっているカードは畳んだ状態で出る。ファイルのないカードと、GFF3 + FASTA の片方だけのカードは開いている。自分で追加したカードは開いている。ファイルを選んでもカードは畳まれない。〔畳んでも変わらないもの〕カードの中身（ファイル、名前、Depth TSV の割り当て、Record options、レコードの一覧）と設定の値。開いたカードの中身は今と同じ（入力欄からはみ出す長いファイル名だけは、両端で切れる代わりに末尾を「…」で省略する）。畳んだ状態は Session に保存せず、Undo の履歴と図にも入らない。カードを上下に移動すると、畳んだ状態も一緒に動く。〔問題の表示〕畳んだカードに直すべき問題（読み取れなかったファイル、足りないファイル、Gencode required）があるときは、畳んだ行にも出る。設定へ案内する操作（レコードのポップアップの Show in settings）は、行き先のカードを開く。
May retire: none
Accepted residual risk: 畳んだカードでは、名前、Depth TSV の割り当て（付いているかどうかを含む）、Record options、レコードの一覧は、開くまで見えない（1 回押す手間が増える）。Session や example を読み込んだ直後はカードが畳まれている。一番狭い設定ペイン（240 px）では、畳んだ行のファイル名は数文字まで省略される。File カードが写っている Gallery と docs の capture は、撮り直すまで見出しの行に開閉の印がない。
Owner: satoshikawato
Decision date: 2026-10-13
```

```json
{
  "concern": "web.linear.file-card-fold",
  "scenarioRevision": 1,
  "choice": "A / FOLDABLE-FILE-CARD, FOLDED-AFTER-LOAD",
  "rationale": "Linear の File カードは 1 枚が長く、ファイルが増えると Comparison や下の section までのスクロールが長くなる（ファイル 5 つのとき、ファイルの一覧は 1,778 px、Comparison は 1,972 px 下）。カードを 1 行に畳めるようにして、ファイルの一覧を 1 画面で見渡せるようにする（畳むと一覧は 218 px、Comparison は 412 px 下）。",
  "mustPreserve": "〔畳む・開く〕各 File カードの見出しの行を押すとカードを 1 行に畳み、もう一度押すと開く。上下の移動ボタンを押したときは畳まれない。畳んだ行には File の番号、ファイル名、上下の移動ボタンが出る。〔最初の状態〕Session または example を読み込んだとき、ファイルがそろっているカードは畳んだ状態で出る。ファイルのないカードと、GFF3 + FASTA の片方だけのカードは開いている。自分で追加したカードは開いている。ファイルを選んでもカードは畳まれない。〔畳んでも変わらないもの〕カードの中身（ファイル、名前、Depth TSV の割り当て、Record options、レコードの一覧）と設定の値。開いたカードの中身は今と同じ（入力欄からはみ出す長いファイル名だけは、両端で切れる代わりに末尾を「…」で省略する）。畳んだ状態は Session に保存せず、Undo の履歴と図にも入らない。カードを上下に移動すると、畳んだ状態も一緒に動く。〔問題の表示〕畳んだカードに直すべき問題（読み取れなかったファイル、足りないファイル、Gencode required）があるときは、畳んだ行にも出る。設定へ案内する操作（レコードのポップアップの Show in settings）は、行き先のカードを開く。",
  "mayRetire": "none",
  "acceptedResidualRisk": "畳んだカードでは、名前、Depth TSV の割り当て（付いているかどうかを含む）、Record options、レコードの一覧は、開くまで見えない（1 回押す手間が増える）。Session や example を読み込んだ直後はカードが畳まれている。一番狭い設定ペイン（240 px）では、畳んだ行のファイル名は数文字まで省略される。File カードが写っている Gallery と docs の capture は、撮り直すまで見出しの行に開閉の印がない。",
  "owner": "satoshikawato",
  "decisionDate": "2026-10-13"
}
```

### PD-OI-095: One naming rule for Legend rows that share a name

- Concern key: `legend.row-name-multiplicity`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: the Choice field of the receipt below, which names no
  choice code
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its nine-field JSON representation below.
- Decision source: the Owner replies of `2026-10-13` quoted verbatim in the
  Revision 47 entry above. The receipt and JSON below reproduce all nine
  supplied fields without translation or additional terms. This record does
  not supersede another decision. The runtime, tests, and documentation that
  implement it merge with it (static Product Contract co-change); this record
  supplies no runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `52dc43606bf87dfd8b1dd53395522fb9d52f43e7b7ac36cd781a3efb8e0db9bc`.

```text
PRODUCT_DECISION
Concern: legend.row-name-multiplicity
Scenario revision: 1
Choice: 凡例のすべての行に 1 つの規則を使う。同じ名前で見本（塗りと線）が同じ行は 1 つにまとめる。同じ名前で見本が違う行が 2 つ以上描かれるときは、その全部の名前に、凡例に並ぶ順で (1), (2), … を付ける。1 回だけ使う名前には番号を付けない。track の行（GC content、GC skew、Depth、Conservation）と region annotation の行にも使う（例: 同じ label の GC content が 2 つあり、色が同じなら GC content の 1 行、色が違えば GC content (1) / GC content (2)）。
Rationale: 凡例の名前を区別する仕組みを番号付けの 1 つに揃える。同じ見本の行を 2 つに分けない。
Must preserve: track と annotation の legend_label、slot の順、見本の色と線、保存済みの Session（0.13.0 を含む）が読めて、凡例の行に付けた編集が同じ track の行に残ること。
May retire: 同じ label の 2 つ目以降の行だけに (2), (3) を付け、1 つ目を番号なしにする名前の付け方。色が同じでも行を分ける動作。
Accepted residual risk: 0.13.0 が GC content / GC content (2) と描いた図は、色が違えば GC content (1) / GC content (2) に、色が同じなら GC content の 1 行に変わる。1 行にまとめられた 2 つ目の行に付けていた編集は使われなくなる。
Owner: satoshikawato
Decision date: 2026-10-13
```

```json
{
  "concern": "legend.row-name-multiplicity",
  "scenarioRevision": 1,
  "choice": "凡例のすべての行に 1 つの規則を使う。同じ名前で見本（塗りと線）が同じ行は 1 つにまとめる。同じ名前で見本が違う行が 2 つ以上描かれるときは、その全部の名前に、凡例に並ぶ順で (1), (2), … を付ける。1 回だけ使う名前には番号を付けない。track の行（GC content、GC skew、Depth、Conservation）と region annotation の行にも使う（例: 同じ label の GC content が 2 つあり、色が同じなら GC content の 1 行、色が違えば GC content (1) / GC content (2)）。",
  "rationale": "凡例の名前を区別する仕組みを番号付けの 1 つに揃える。同じ見本の行を 2 つに分けない。",
  "mustPreserve": "track と annotation の legend_label、slot の順、見本の色と線、保存済みの Session（0.13.0 を含む）が読めて、凡例の行に付けた編集が同じ track の行に残ること。",
  "mayRetire": "同じ label の 2 つ目以降の行だけに (2), (3) を付け、1 つ目を番号なしにする名前の付け方。色が同じでも行を分ける動作。",
  "acceptedResidualRisk": "0.13.0 が GC content / GC content (2) と描いた図は、色が違えば GC content (1) / GC content (2) に、色が同じなら GC content の 1 行に変わる。1 行にまとめられた 2 つ目の行に付けていた編集は使われなくなる。",
  "owner": "satoshikawato",
  "decisionDate": "2026-10-13"
}
```

### PD-OI-096: Legend position of a diagram mode at Load

- Concern key: `web.legend.position-of-mode`
- Scenario revision: `2`
- Supersedes: `PD-OI-096`, scenario revision `1`. Scenario revision `2`
  keeps every field of scenario revision `1` except its Accepted residual
  risk, its scenario revision, and its decision date.
- Status: `ACCEPTED`
- Selected outcome: `B / MODE-DEFAULT-LISTED-AT-LOAD`
- Normative outcome: exactly the complete approved scenario revision `2`
  `PRODUCT_DECISION` receipt and its JSON representation below.
- Decision source: the Owner replies quoted verbatim in the Revision 48 entry
  (scenario revision `1`) and the Revision 55 entry (scenario revision `2`)
  above. The Owner approved the scenario revision `1` receipt as written, and
  the Accepted residual risk of scenario revision `2` as written; the other
  fields are those of scenario revision `1`. The receipt reproduces its
  fields without translation or additional terms. Dependent runtime requires
  this authority merged into its base, or merged together with it through the
  reviewed static Product Contract co-change route; this record supplies no
  runtime acceptance evidence.
- Receipt SHA-256 (scenario revision `2`; UTF-8, excluding the final newline):
  `abb6fe7056b5b15c82a883c3b5e3d9f13892a472a499e9fbca3986cdda2f6a87`.

```text
PRODUCT_DECISION
Concern: web.legend.position-of-mode
Scenario revision: 2
Choice: B / MODE-DEFAULT-LISTED-AT-LOAD
Rationale: 0.13.0 までの版が保存したファイルにも、その mode で使えない凡例の位置（Linear の upper_left など）がありうる。拒否せずに読み込み、Generate が失敗しないようにする。直せなかった値は読み込む前に知らせる。Session 31 以降はすでに描画済みの側をこう読んでいるので、古い Session もそろえる。
Must preserve: 〔mode の位置〕mode の集合にある位置（Circular: left、right、top、bottom、upper_left、upper_right、lower_left、lower_right、none。Linear: left、right、top、bottom、none）の読み込み、保存、Generate、History。大文字や前後の空白を含む値は小文字にして位置として読み、知らせない。〔集合の外の値〕その mode で使えない値は、設定の凡例の位置なら値がないときと同じ位置（ふつうはその mode の既定の位置。Circular は left、Linear は bottom）として、読み込んだ Result に保存された描画済みの凡例の側ならその mode の既定の位置として読み、読み込む前の確認 dialog に載せる。構成 metadata に保存された側が空のときも、その mode の既定の位置として読む。〔値がないとき〕値がないときは現在の既定の位置のまま（構成 metadata のない古い Result の描画済みの側は、空白だけの値も含めて none）。〔図〕読み込んだ直後の図は保存されたときのまま。〔CLI と Python〕CLI と Python API の検証とエラー文は変えない。
May retire: 集合の外の凡例の位置をそのまま保持し、Generate が "Input validation failed" で失敗する動作。集合の外の描画済みの側で、最初の凡例の編集が Linear の角へ凡例を置く動作。大文字や空白を含む、または集合の外や空の描画済みの側で、読み込みや表示が失敗する動作。
Accepted residual risk: 0.13.0 が Linear の角に描いた凡例は、最初の凡例の編集で図の下へ移る。0.13.0 までの版が mode を切り替えてから保存した Session のうち、最後の Generate の mode を記録していないものは、保存したときに表示していた mode の Result として読む。
Owner: satoshikawato
Decision date: 2026-10-16
```

```json
{
  "concern": "web.legend.position-of-mode",
  "scenarioRevision": 2,
  "choice": "B / MODE-DEFAULT-LISTED-AT-LOAD",
  "rationale": "0.13.0 までの版が保存したファイルにも、その mode で使えない凡例の位置（Linear の upper_left など）がありうる。拒否せずに読み込み、Generate が失敗しないようにする。直せなかった値は読み込む前に知らせる。Session 31 以降はすでに描画済みの側をこう読んでいるので、古い Session もそろえる。",
  "mustPreserve": "〔mode の位置〕mode の集合にある位置（Circular: left、right、top、bottom、upper_left、upper_right、lower_left、lower_right、none。Linear: left、right、top、bottom、none）の読み込み、保存、Generate、History。大文字や前後の空白を含む値は小文字にして位置として読み、知らせない。〔集合の外の値〕その mode で使えない値は、設定の凡例の位置なら値がないときと同じ位置（ふつうはその mode の既定の位置。Circular は left、Linear は bottom）として、読み込んだ Result に保存された描画済みの凡例の側ならその mode の既定の位置として読み、読み込む前の確認 dialog に載せる。構成 metadata に保存された側が空のときも、その mode の既定の位置として読む。〔値がないとき〕値がないときは現在の既定の位置のまま（構成 metadata のない古い Result の描画済みの側は、空白だけの値も含めて none）。〔図〕読み込んだ直後の図は保存されたときのまま。〔CLI と Python〕CLI と Python API の検証とエラー文は変えない。",
  "mayRetire": "集合の外の凡例の位置をそのまま保持し、Generate が \"Input validation failed\" で失敗する動作。集合の外の描画済みの側で、最初の凡例の編集が Linear の角へ凡例を置く動作。大文字や空白を含む、または集合の外や空の描画済みの側で、読み込みや表示が失敗する動作。",
  "acceptedResidualRisk": "0.13.0 が Linear の角に描いた凡例は、最初の凡例の編集で図の下へ移る。0.13.0 までの版が mode を切り替えてから保存した Session のうち、最後の Generate の mode を記録していないものは、保存したときに表示していた mode の Result として読む。",
  "owner": "satoshikawato",
  "decisionDate": "2026-10-16"
}
```

### PD-OI-097: Labels that do not fit

- Concern key: `product.circular-outer-label-fit`
- Scenario revision: `11`
- Supersedes: `PD-OI-097`, scenario revision `10`. Scenario revision `11`
  keeps every outcome of scenario revision `10` and adds `LBL-HALO`,
  `LIN-UNDER`, `CIR-NO-UNDER`, `GAP-STOP`, `LEADER-007`, `LBL-OPTS-2`,
  `LBL-0140`, `HALO-008`, `OPT-NAMES`, `SESS-NEXT-GEN`, `WEBPDF-NO-HALO`,
  `UNDER-NOT-IN-ROW`, `UNDER-CASING`, `UNDER-TRIM`, `ROW-GAP-025`, and
  `LIN-BOX-080`; it keeps every field of scenario revision `10` except the
  Rationale, the May retire, and the Accepted residual risk, to each of which
  Replies 44 and 46 of the Revision 68 entry append the feature label
  style, and the Must preserve, in which they replace the sentence on
  figures without overlaps or crossings. Scenario revision `10`
  keeps every outcome of scenario revision `9` and adds `RGN-CUT-WARN`; it
  keeps every field of scenario revision `9` except the Rationale, to which
  Reply 24 of the Revision 67 entry appends the Circular region labels cut to
  fit their circle, and the decision date. Scenario revision `9`
  keeps every outcome of scenario revision `8` and adds `REGION-WARN`,
  `RGN-FRONT`, and `RGN-YIELD-ROW`; it keeps every field of scenario revision
  `8` except the Rationale and the May retire, to each of which Reply 22 of
  the Revision 66 entry appends the exception for Linear region labels. Scenario revision `8` kept every outcome
  and field of scenario revision `7` and added `YIELD-REPAIRED-ONLY`.
- Status: `ACCEPTED`
- Selected outcome: `DEF-CROSS ; DO-ASK ; TIGHT-FIX-FIRST ; SHOW-ALL ; REL-0140 ; CROSS-WARN ; ROWS-OPEN ; WARN-NOTICE ; SEARCH-BOUNDED-5S ; SCOPE-ALL-LABELS ; INNER-3 ; SAME-RULES ; SWITCHES-FIRST ; RADIAL-JOIN ; DEF-WARN ; E-ROW ; REPAIR-100 ; CROSS-LEFT-ROW ; YIELD-REPAIRED-ONLY ; REGION-WARN ; RGN-FRONT ; RGN-YIELD-ROW ; RGN-CUT-WARN ; LBL-HALO ; LIN-UNDER ; CIR-NO-UNDER ; GAP-STOP ; LEADER-007 ; LBL-OPTS-2 ; LBL-0140 ; HALO-008 ; OPT-NAMES ; SESS-NEXT-GEN ; WEBPDF-NO-HALO ; UNDER-NOT-IN-ROW ; UNDER-CASING ; UNDER-TRIM ; ROW-GAP-025 ; LIN-BOX-080`
- Normative outcome: exactly the complete approved scenario revision `11`
  `PRODUCT_DECISION` receipt and its JSON representation below. The
  scenario revision `7`, `8`, `9`, and `10` receipts are kept as they were
  approved.
- Decision source (scenario revision `7`): the Owner replies quoted verbatim in
  the Revision 50 entry above. The receipt reproduces the approved field texts
  without translation or additional terms. Dependent runtime requires this
  authority merged into its base, or merged together with it through the
  reviewed static Product Contract co-change route; this record supplies no
  runtime acceptance evidence.
- Receipt SHA-256 (scenario revision `7`; UTF-8, excluding the final newline):
  `e77d81603000fc47273f5bc86fc1331bf27b666be988cef11cc388301a7324e8`.
- Decision source (scenario revision `8`): the Owner replies of `2026-10-14`
  quoted verbatim in the Revision 53 entry above. The second receipt and JSON
  below reproduce all nine approved fields without translation or additional
  terms. The runtime, tests, and documentation that implement it merge with it
  (static Product Contract co-change); this record supplies no runtime
  acceptance evidence.
- Receipt SHA-256 (scenario revision `8`; UTF-8, excluding the final newline):
  `405597f512e865299b743d4307924a44955721a391917fa7360230ba4004417b`.
- Decision source (scenario revision `9`): the Owner replies of `2026-10-18`
  and `2026-10-19` quoted verbatim in the Revision 66 entry above; Reply 22
  approves the third receipt and JSON below. They keep the fields of scenario
  revision `8`, add the three choices and the two exceptions of Reply 22, and
  carry the decision date of the last reply, without translation or additional terms. The runtime, tests, and documentation that implement it
  merge with it (static Product Contract co-change); this record supplies no
  runtime acceptance evidence.
- Receipt SHA-256 (scenario revision `9`; UTF-8, excluding the final newline):
  `65994020bd2ce347daa2d39f5e7fb470b7fe1d3efdf5600fcd43a4201ad3065d`.
- Decision source (scenario revision `10`): the Owner replies of `2026-10-19`
  quoted verbatim in the Revision 67 entry above; Reply 24 approves the fourth
  receipt and JSON below. They keep the fields of scenario revision `9`, add
  the choice and the Rationale sentence of Reply 24, and carry the decision
  date of the last reply, without translation or additional terms. The
  runtime, tests, and documentation that implement it merge with it (static
  Product Contract co-change); this record supplies no runtime acceptance
  evidence.
- Receipt SHA-256 (scenario revision `10`; UTF-8, excluding the final newline):
  `ae0bca6a9a6e4778ebc0a1bc9b12976a34eb389d874eab7c8e8857d502d5ee16`.
- Decision source (scenario revision `11`): the Owner replies of `2026-10-19`
  quoted verbatim in the Revision 68 entry above; Reply 46 approves the fifth
  receipt and JSON below, which correct one Must preserve sentence of the text
  that Reply 44 approved. They keep the fields of scenario revision `10`, add
  the choices and the sentences of Replies 44 and 46, replace its Must preserve
  sentence on figures without overlaps or crossings, and carry the decision
  date of the last reply, without translation or additional terms. The
  runtime, tests, and documentation that implement it merge with it (static
  Product Contract co-change); this record supplies no runtime acceptance
  evidence.
- Receipt SHA-256 (scenario revision `11`; UTF-8, excluding the final newline):
  `21bb6008d55e26f39220929f1f0047a4a2dcdd5ee9b10610015c4b512bb5a6ef`.

```text
PRODUCT_DECISION
Concern: product.circular-outer-label-fit
Scenario revision: 7
Choice: DEF-CROSS ; DO-ASK ; TIGHT-FIX-FIRST ; SHOW-ALL ; REL-0140 ; CROSS-WARN ; ROWS-OPEN ; WARN-NOTICE ; SEARCH-BOUNDED-5S ; SCOPE-ALL-LABELS ; INNER-3 ; SAME-RULES ; SWITCHES-FIRST ; RADIAL-JOIN ; DEF-WARN ; E-ROW ; REPAIR-100 ; CROSS-LEFT-ROW
Rationale: ラベルが重なった図を黙って出さない。重なるときは描かずに直し方を選ばせ、引き出し線の交差や定義文へのかぶりだけなら描いたうえで知らせる。
Must preserve: 今まで重ならずに描けていた図は今まで通り描く（rev 2 の『今まで描けてたものは今まで通り描けるべき』）。重なりも交差もない図は 0.13.0 と同じに描く。ダイアログでの選択は History の 1 ステップ、Cancel は記録しない、適用中は操作を受け付けない（PD-OI-088）。
May retire: ラベルが重なったまま図を出す動作。
Accepted residual risk: 微調整で全部直る重なり・交差・定義文へのかぶり（欠陥 100 組まで）は、自動で直して描く。直しきれないもの（極の内側ラベルの積み重なり、長い線のもつれ、MjeNMV など）は、重なりなら描かずに直し方を選ばせ、交差やかぶりだけなら描いて帯で知らせる。直し方の候補は、帯のボタンを押したときに計算する。この結果は Session に保存しない。CLI と Python は、描けた図なら 5 秒まで探す。描けない図は、外側ラベルだけの重なりなら値が出るまで探し、内側・埋め込み・放射ラベルの重なりなら 30 秒まで探す。間に合わない候補は『未計算』と示す。Web の候補は 1 つ 30 秒まで探し、見つからなければそう示す。
Owner: satoshikawato
Decision date: 2026-10-14
```

```json
{
  "concern": "product.circular-outer-label-fit",
  "scenarioRevision": 7,
  "choice": "DEF-CROSS ; DO-ASK ; TIGHT-FIX-FIRST ; SHOW-ALL ; REL-0140 ; CROSS-WARN ; ROWS-OPEN ; WARN-NOTICE ; SEARCH-BOUNDED-5S ; SCOPE-ALL-LABELS ; INNER-3 ; SAME-RULES ; SWITCHES-FIRST ; RADIAL-JOIN ; DEF-WARN ; E-ROW ; REPAIR-100 ; CROSS-LEFT-ROW",
  "rationale": "ラベルが重なった図を黙って出さない。重なるときは描かずに直し方を選ばせ、引き出し線の交差や定義文へのかぶりだけなら描いたうえで知らせる。",
  "mustPreserve": "今まで重ならずに描けていた図は今まで通り描く（rev 2 の『今まで描けてたものは今まで通り描けるべき』）。重なりも交差もない図は 0.13.0 と同じに描く。ダイアログでの選択は History の 1 ステップ、Cancel は記録しない、適用中は操作を受け付けない（PD-OI-088）。",
  "mayRetire": "ラベルが重なったまま図を出す動作。",
  "acceptedResidualRisk": "微調整で全部直る重なり・交差・定義文へのかぶり（欠陥 100 組まで）は、自動で直して描く。直しきれないもの（極の内側ラベルの積み重なり、長い線のもつれ、MjeNMV など）は、重なりなら描かずに直し方を選ばせ、交差やかぶりだけなら描いて帯で知らせる。直し方の候補は、帯のボタンを押したときに計算する。この結果は Session に保存しない。CLI と Python は、描けた図なら 5 秒まで探す。描けない図は、外側ラベルだけの重なりなら値が出るまで探し、内側・埋め込み・放射ラベルの重なりなら 30 秒まで探す。間に合わない候補は『未計算』と示す。Web の候補は 1 つ 30 秒まで探し、見つからなければそう示す。",
  "owner": "satoshikawato",
  "decisionDate": "2026-10-14"
}
```

Scenario revision `8` receipt, confirmed by `satoshikawato` on `2026-10-14`:

```text
PRODUCT_DECISION
Concern: product.circular-outer-label-fit
Scenario revision: 8
Choice: DEF-CROSS ; DO-ASK ; TIGHT-FIX-FIRST ; SHOW-ALL ; REL-0140 ; CROSS-WARN ; ROWS-OPEN ; WARN-NOTICE ; SEARCH-BOUNDED-5S ; SCOPE-ALL-LABELS ; INNER-3 ; SAME-RULES ; SWITCHES-FIRST ; RADIAL-JOIN ; DEF-WARN ; E-ROW ; REPAIR-100 ; CROSS-LEFT-ROW ; YIELD-REPAIRED-ONLY
Rationale: ラベルが重なった図を黙って出さない。重なるときは描かずに直し方を選ばせ、引き出し線の交差や定義文へのかぶりだけなら描いたうえで知らせる。
Must preserve: 今まで重ならずに描けていた図は今まで通り描く（rev 2 の『今まで描けてたものは今まで通り描けるべき』）。重なりも交差もない図は 0.13.0 と同じに描く。ダイアログでの選択は History の 1 ステップ、Cancel は記録しない、適用中は操作を受け付けない（PD-OI-088）。
May retire: ラベルが重なったまま図を出す動作。
Accepted residual risk: 微調整で全部直る重なり・交差・定義文へのかぶり（欠陥 100 組まで）は、自動で直して描く。直しきれないもの（極の内側ラベルの積み重なり、長い線のもつれ、MjeNMV など）は、重なりなら描かずに直し方を選ばせ、交差やかぶりだけなら描いて帯で知らせる。直し方の候補は、帯のボタンを押したときに計算する。この結果は Session に保存しない。CLI と Python は、描けた図なら 5 秒まで探す。描けない図は、外側ラベルだけの重なりなら値が出るまで探し、内側・埋め込み・放射ラベルの重なりなら 30 秒まで探す。間に合わない候補は『未計算』と示す。Web の候補は 1 つ 30 秒まで探し、見つからなければそう示す。
Owner: satoshikawato
Decision date: 2026-10-14
```

```json
{
  "concern": "product.circular-outer-label-fit",
  "scenarioRevision": 8,
  "choice": "DEF-CROSS ; DO-ASK ; TIGHT-FIX-FIRST ; SHOW-ALL ; REL-0140 ; CROSS-WARN ; ROWS-OPEN ; WARN-NOTICE ; SEARCH-BOUNDED-5S ; SCOPE-ALL-LABELS ; INNER-3 ; SAME-RULES ; SWITCHES-FIRST ; RADIAL-JOIN ; DEF-WARN ; E-ROW ; REPAIR-100 ; CROSS-LEFT-ROW ; YIELD-REPAIRED-ONLY",
  "rationale": "ラベルが重なった図を黙って出さない。重なるときは描かずに直し方を選ばせ、引き出し線の交差や定義文へのかぶりだけなら描いたうえで知らせる。",
  "mustPreserve": "今まで重ならずに描けていた図は今まで通り描く（rev 2 の『今まで描けてたものは今まで通り描けるべき』）。重なりも交差もない図は 0.13.0 と同じに描く。ダイアログでの選択は History の 1 ステップ、Cancel は記録しない、適用中は操作を受け付けない（PD-OI-088）。",
  "mayRetire": "ラベルが重なったまま図を出す動作。",
  "acceptedResidualRisk": "微調整で全部直る重なり・交差・定義文へのかぶり（欠陥 100 組まで）は、自動で直して描く。直しきれないもの（極の内側ラベルの積み重なり、長い線のもつれ、MjeNMV など）は、重なりなら描かずに直し方を選ばせ、交差やかぶりだけなら描いて帯で知らせる。直し方の候補は、帯のボタンを押したときに計算する。この結果は Session に保存しない。CLI と Python は、描けた図なら 5 秒まで探す。描けない図は、外側ラベルだけの重なりなら値が出るまで探し、内側・埋め込み・放射ラベルの重なりなら 30 秒まで探す。間に合わない候補は『未計算』と示す。Web の候補は 1 つ 30 秒まで探し、見つからなければそう示す。",
  "owner": "satoshikawato",
  "decisionDate": "2026-10-14"
}
```

Scenario revision `9` receipt, approved by `satoshikawato` on `2026-10-19` (Reply 22 of the Revision 66 entry):

```text
PRODUCT_DECISION
Concern: product.circular-outer-label-fit
Scenario revision: 9
Choice: DEF-CROSS ; DO-ASK ; TIGHT-FIX-FIRST ; SHOW-ALL ; REL-0140 ; CROSS-WARN ; ROWS-OPEN ; WARN-NOTICE ; SEARCH-BOUNDED-5S ; SCOPE-ALL-LABELS ; INNER-3 ; SAME-RULES ; SWITCHES-FIRST ; RADIAL-JOIN ; DEF-WARN ; E-ROW ; REPAIR-100 ; CROSS-LEFT-ROW ; YIELD-REPAIRED-ONLY ; REGION-WARN ; RGN-FRONT ; RGN-YIELD-ROW
Rationale: ラベルが重なった図を黙って出さない。重なるときは描かずに直し方を選ばせ、引き出し線の交差や定義文へのかぶりだけなら描いたうえで知らせる。ただし Linear の領域ラベルの重なりは、図を描いたうえで警告する（LP-D68）。
Must preserve: 今まで重ならずに描けていた図は今まで通り描く（rev 2 の『今まで描けてたものは今まで通り描けるべき』）。重なりも交差もない図は 0.13.0 と同じに描く。ダイアログでの選択は History の 1 ステップ、Cancel は記録しない、適用中は操作を受け付けない（PD-OI-088）。
May retire: ラベルが重なったまま図を出す動作。（Linear の領域ラベルの重なりを除く）
Accepted residual risk: 微調整で全部直る重なり・交差・定義文へのかぶり（欠陥 100 組まで）は、自動で直して描く。直しきれないもの（極の内側ラベルの積み重なり、長い線のもつれ、MjeNMV など）は、重なりなら描かずに直し方を選ばせ、交差やかぶりだけなら描いて帯で知らせる。直し方の候補は、帯のボタンを押したときに計算する。この結果は Session に保存しない。CLI と Python は、描けた図なら 5 秒まで探す。描けない図は、外側ラベルだけの重なりなら値が出るまで探し、内側・埋め込み・放射ラベルの重なりなら 30 秒まで探す。間に合わない候補は『未計算』と示す。Web の候補は 1 つ 30 秒まで探し、見つからなければそう示す。
Owner: satoshikawato
Decision date: 2026-10-19
```

```json
{
  "concern": "product.circular-outer-label-fit",
  "scenarioRevision": 9,
  "choice": "DEF-CROSS ; DO-ASK ; TIGHT-FIX-FIRST ; SHOW-ALL ; REL-0140 ; CROSS-WARN ; ROWS-OPEN ; WARN-NOTICE ; SEARCH-BOUNDED-5S ; SCOPE-ALL-LABELS ; INNER-3 ; SAME-RULES ; SWITCHES-FIRST ; RADIAL-JOIN ; DEF-WARN ; E-ROW ; REPAIR-100 ; CROSS-LEFT-ROW ; YIELD-REPAIRED-ONLY ; REGION-WARN ; RGN-FRONT ; RGN-YIELD-ROW",
  "rationale": "ラベルが重なった図を黙って出さない。重なるときは描かずに直し方を選ばせ、引き出し線の交差や定義文へのかぶりだけなら描いたうえで知らせる。ただし Linear の領域ラベルの重なりは、図を描いたうえで警告する（LP-D68）。",
  "mustPreserve": "今まで重ならずに描けていた図は今まで通り描く（rev 2 の『今まで描けてたものは今まで通り描けるべき』）。重なりも交差もない図は 0.13.0 と同じに描く。ダイアログでの選択は History の 1 ステップ、Cancel は記録しない、適用中は操作を受け付けない（PD-OI-088）。",
  "mayRetire": "ラベルが重なったまま図を出す動作。（Linear の領域ラベルの重なりを除く）",
  "acceptedResidualRisk": "微調整で全部直る重なり・交差・定義文へのかぶり（欠陥 100 組まで）は、自動で直して描く。直しきれないもの（極の内側ラベルの積み重なり、長い線のもつれ、MjeNMV など）は、重なりなら描かずに直し方を選ばせ、交差やかぶりだけなら描いて帯で知らせる。直し方の候補は、帯のボタンを押したときに計算する。この結果は Session に保存しない。CLI と Python は、描けた図なら 5 秒まで探す。描けない図は、外側ラベルだけの重なりなら値が出るまで探し、内側・埋め込み・放射ラベルの重なりなら 30 秒まで探す。間に合わない候補は『未計算』と示す。Web の候補は 1 つ 30 秒まで探し、見つからなければそう示す。",
  "owner": "satoshikawato",
  "decisionDate": "2026-10-19"
}
```

Scenario revision `10` receipt, approved by `satoshikawato` on `2026-10-19` (Reply 24 of the Revision 67 entry):

```text
PRODUCT_DECISION
Concern: product.circular-outer-label-fit
Scenario revision: 10
Choice: DEF-CROSS ; DO-ASK ; TIGHT-FIX-FIRST ; SHOW-ALL ; REL-0140 ; CROSS-WARN ; ROWS-OPEN ; WARN-NOTICE ; SEARCH-BOUNDED-5S ; SCOPE-ALL-LABELS ; INNER-3 ; SAME-RULES ; SWITCHES-FIRST ; RADIAL-JOIN ; DEF-WARN ; E-ROW ; REPAIR-100 ; CROSS-LEFT-ROW ; YIELD-REPAIRED-ONLY ; REGION-WARN ; RGN-FRONT ; RGN-YIELD-ROW ; RGN-CUT-WARN
Rationale: ラベルが重なった図を黙って出さない。重なるときは描かずに直し方を選ばせ、引き出し線の交差や定義文へのかぶりだけなら描いたうえで知らせる。ただし Linear の領域ラベルの重なりは、図を描いたうえで警告する（LP-D68）。Circular の領域ラベルが円に収まらないときは、収まる分を曲線で描いたうえで警告する（LP-D87）。
Must preserve: 今まで重ならずに描けていた図は今まで通り描く（rev 2 の『今まで描けてたものは今まで通り描けるべき』）。重なりも交差もない図は 0.13.0 と同じに描く。ダイアログでの選択は History の 1 ステップ、Cancel は記録しない、適用中は操作を受け付けない（PD-OI-088）。
May retire: ラベルが重なったまま図を出す動作。（Linear の領域ラベルの重なりを除く）
Accepted residual risk: 微調整で全部直る重なり・交差・定義文へのかぶり（欠陥 100 組まで）は、自動で直して描く。直しきれないもの（極の内側ラベルの積み重なり、長い線のもつれ、MjeNMV など）は、重なりなら描かずに直し方を選ばせ、交差やかぶりだけなら描いて帯で知らせる。直し方の候補は、帯のボタンを押したときに計算する。この結果は Session に保存しない。CLI と Python は、描けた図なら 5 秒まで探す。描けない図は、外側ラベルだけの重なりなら値が出るまで探し、内側・埋め込み・放射ラベルの重なりなら 30 秒まで探す。間に合わない候補は『未計算』と示す。Web の候補は 1 つ 30 秒まで探し、見つからなければそう示す。
Owner: satoshikawato
Decision date: 2026-10-19
```

```json
{
  "concern": "product.circular-outer-label-fit",
  "scenarioRevision": 10,
  "choice": "DEF-CROSS ; DO-ASK ; TIGHT-FIX-FIRST ; SHOW-ALL ; REL-0140 ; CROSS-WARN ; ROWS-OPEN ; WARN-NOTICE ; SEARCH-BOUNDED-5S ; SCOPE-ALL-LABELS ; INNER-3 ; SAME-RULES ; SWITCHES-FIRST ; RADIAL-JOIN ; DEF-WARN ; E-ROW ; REPAIR-100 ; CROSS-LEFT-ROW ; YIELD-REPAIRED-ONLY ; REGION-WARN ; RGN-FRONT ; RGN-YIELD-ROW ; RGN-CUT-WARN",
  "rationale": "ラベルが重なった図を黙って出さない。重なるときは描かずに直し方を選ばせ、引き出し線の交差や定義文へのかぶりだけなら描いたうえで知らせる。ただし Linear の領域ラベルの重なりは、図を描いたうえで警告する（LP-D68）。Circular の領域ラベルが円に収まらないときは、収まる分を曲線で描いたうえで警告する（LP-D87）。",
  "mustPreserve": "今まで重ならずに描けていた図は今まで通り描く（rev 2 の『今まで描けてたものは今まで通り描けるべき』）。重なりも交差もない図は 0.13.0 と同じに描く。ダイアログでの選択は History の 1 ステップ、Cancel は記録しない、適用中は操作を受け付けない（PD-OI-088）。",
  "mayRetire": "ラベルが重なったまま図を出す動作。（Linear の領域ラベルの重なりを除く）",
  "acceptedResidualRisk": "微調整で全部直る重なり・交差・定義文へのかぶり（欠陥 100 組まで）は、自動で直して描く。直しきれないもの（極の内側ラベルの積み重なり、長い線のもつれ、MjeNMV など）は、重なりなら描かずに直し方を選ばせ、交差やかぶりだけなら描いて帯で知らせる。直し方の候補は、帯のボタンを押したときに計算する。この結果は Session に保存しない。CLI と Python は、描けた図なら 5 秒まで探す。描けない図は、外側ラベルだけの重なりなら値が出るまで探し、内側・埋め込み・放射ラベルの重なりなら 30 秒まで探す。間に合わない候補は『未計算』と示す。Web の候補は 1 つ 30 秒まで探し、見つからなければそう示す。",
  "owner": "satoshikawato",
  "decisionDate": "2026-10-19"
}
```

Scenario revision `11` receipt, approved by `satoshikawato` on `2026-10-19` (Reply 46 of the Revision 68 entry):

```text
PRODUCT_DECISION
Concern: product.circular-outer-label-fit
Scenario revision: 11
Choice: DEF-CROSS ; DO-ASK ; TIGHT-FIX-FIRST ; SHOW-ALL ; REL-0140 ; CROSS-WARN ; ROWS-OPEN ; WARN-NOTICE ; SEARCH-BOUNDED-5S ; SCOPE-ALL-LABELS ; INNER-3 ; SAME-RULES ; SWITCHES-FIRST ; RADIAL-JOIN ; DEF-WARN ; E-ROW ; REPAIR-100 ; CROSS-LEFT-ROW ; YIELD-REPAIRED-ONLY ; REGION-WARN ; RGN-FRONT ; RGN-YIELD-ROW ; RGN-CUT-WARN ; LBL-HALO ; LIN-UNDER ; CIR-NO-UNDER ; GAP-STOP ; LEADER-007 ; LBL-OPTS-2 ; LBL-0140 ; HALO-008 ; OPT-NAMES ; SESS-NEXT-GEN ; WEBPDF-NO-HALO ; UNDER-NOT-IN-ROW ; UNDER-CASING ; UNDER-TRIM ; ROW-GAP-025 ; LIN-BOX-080
Rationale: ラベルが重なった図を黙って出さない。重なるときは描かずに直し方を選ばせ、引き出し線の交差や定義文へのかぶりだけなら描いたうえで知らせる。ただし Linear の領域ラベルの重なりは、図を描いたうえで警告する（LP-D68）。Circular の領域ラベルが円に収まらないときは、収まる分を曲線で描いたうえで警告する（LP-D87）。フィーチャーラベルは、後ろを通るリーダー線から守る。Linear は文字の後ろに白い背景（不透明度 0.8、文字の外 0.18 em）を敷き、Circular は文字の形の白いハロー（0.08 em）で守る（LP-D98、LP-D106）。Linear は下線をリーダー線につなぎ、下線に白い縁を付けて、ほかのラベルのリーダー線を下線の所で途切れさせる。下線の両端は 0.12 em ずつ短くする。Circular は下線を引かない。文字の上から来て下線につながるリーダー線は、文字の 0.20 em 上で止める。リーダー線と下線は同じ太さ max(1.0 px, 0.07 × 文字の大きさ) で描く。下線（label_line）と守り（label_protection：なし、ハロー、白い背景）はモードごとに選べる（LP-D98、LP-D100、LP-D103）。Linear の下線と白い背景はラベルの段の高さに数えず、領域ラベルだけがそれをよける（LP-D102）。Linear の同じ段のラベルは 0.25 em 以上離す（LP-D104）。
Must preserve: 今まで重ならずに描けていた図は今まで通り描く（rev 2 の『今まで描けてたものは今まで通り描けるべき』）。重なりも交差もない図は 0.13.0 と同じ位置にラベルを描く。ただし Linear で同じ段の隣のラベルとの間が 0.25 em に満たないときは、そのラベルを別の段に移す。Linear の領域ラベルは、フィーチャーラベルの下線と白い背景を文字と同じによける。よけられないときは、文字にかかるときと同じく、そのフィーチャーラベルを同じ側の別の段へ移す（RGN-YIELD-ROW）。移しても解消しないときは、そのまま描いて警告する。領域ラベルの背景は、リーダー線の上には出るが、下線と白い背景は覆わない（RGN-FRONT）。フィーチャーラベルの守り（Linear は白い背景、Circular はハロー）、下線とその白い縁、リーダー線の太さと止め方は LP-D98 と LP-D106 の既定で描く。Linear の白い背景は、ほかのラベルの文字と下線を覆わない。ダイアログでの選択は History の 1 ステップ、Cancel は記録しない、適用中は操作を受け付けない（PD-OI-088）。
May retire: ラベルが重なったまま図を出す動作。（Linear の領域ラベルの重なりを除く）0.13.0 のフィーチャーラベルの見た目（細いリーダー線、守りと下線なし）を既定で出すこと。labels.stroke_width の設定（Python の config_overrides で指定するとエラーになる）。Linear で同じ段のラベルをすき間なしで並べること。
Accepted residual risk: 微調整で全部直る重なり・交差・定義文へのかぶり（欠陥 100 組まで）は、自動で直して描く。直しきれないもの（極の内側ラベルの積み重なり、長い線のもつれ、MjeNMV など）は、重なりなら描かずに直し方を選ばせ、交差やかぶりだけなら描いて帯で知らせる。直し方の候補は、帯のボタンを押したときに計算する。この結果は Session に保存しない。CLI と Python は、描けた図なら 5 秒まで探す。描けない図は、外側ラベルだけの重なりなら値が出るまで探し、内側・埋め込み・放射ラベルの重なりなら 30 秒まで探す。間に合わない候補は『未計算』と示す。Web の候補は 1 つ 30 秒まで探し、見つからなければそう示す。Circular の PDF と EPS はハローの分だけ大きくなる（タバコ葉緑体の PDF は 56 KB から 235 KB、EPS は 1.32 MB から 1.59 MB）。Linear の白い背景では PDF はほとんど変わらず（lambda Linear 24 pt で 17.6 KB から 27 KB）、EPS は半透明の背景を画像で描くので大きくなる（65 KB から 664 KB）。Web の PDF には Circular のハローを描かない（Linear の白い背景、下線とその白い縁は描く。LP-D102）。Web の PDF は文字をカーニングせずに描くので、ラベルが SVG より最大 0.1 em ほど広くなることがある（中心の位置は同じ）。保存済みの Session は次の Generate で新しい見た目になる（それまでは保存した図のまま）。Linear の詰まった段では、ラベルが別の段に移る（lambda 24 pt で 73 のうち 21、タバコ IR 5 pt で 143 のうち 22）。
Owner: satoshikawato
Decision date: 2026-10-19
```

```json
{
  "concern": "product.circular-outer-label-fit",
  "scenarioRevision": 11,
  "choice": "DEF-CROSS ; DO-ASK ; TIGHT-FIX-FIRST ; SHOW-ALL ; REL-0140 ; CROSS-WARN ; ROWS-OPEN ; WARN-NOTICE ; SEARCH-BOUNDED-5S ; SCOPE-ALL-LABELS ; INNER-3 ; SAME-RULES ; SWITCHES-FIRST ; RADIAL-JOIN ; DEF-WARN ; E-ROW ; REPAIR-100 ; CROSS-LEFT-ROW ; YIELD-REPAIRED-ONLY ; REGION-WARN ; RGN-FRONT ; RGN-YIELD-ROW ; RGN-CUT-WARN ; LBL-HALO ; LIN-UNDER ; CIR-NO-UNDER ; GAP-STOP ; LEADER-007 ; LBL-OPTS-2 ; LBL-0140 ; HALO-008 ; OPT-NAMES ; SESS-NEXT-GEN ; WEBPDF-NO-HALO ; UNDER-NOT-IN-ROW ; UNDER-CASING ; UNDER-TRIM ; ROW-GAP-025 ; LIN-BOX-080",
  "rationale": "ラベルが重なった図を黙って出さない。重なるときは描かずに直し方を選ばせ、引き出し線の交差や定義文へのかぶりだけなら描いたうえで知らせる。ただし Linear の領域ラベルの重なりは、図を描いたうえで警告する（LP-D68）。Circular の領域ラベルが円に収まらないときは、収まる分を曲線で描いたうえで警告する（LP-D87）。フィーチャーラベルは、後ろを通るリーダー線から守る。Linear は文字の後ろに白い背景（不透明度 0.8、文字の外 0.18 em）を敷き、Circular は文字の形の白いハロー（0.08 em）で守る（LP-D98、LP-D106）。Linear は下線をリーダー線につなぎ、下線に白い縁を付けて、ほかのラベルのリーダー線を下線の所で途切れさせる。下線の両端は 0.12 em ずつ短くする。Circular は下線を引かない。文字の上から来て下線につながるリーダー線は、文字の 0.20 em 上で止める。リーダー線と下線は同じ太さ max(1.0 px, 0.07 × 文字の大きさ) で描く。下線（label_line）と守り（label_protection：なし、ハロー、白い背景）はモードごとに選べる（LP-D98、LP-D100、LP-D103）。Linear の下線と白い背景はラベルの段の高さに数えず、領域ラベルだけがそれをよける（LP-D102）。Linear の同じ段のラベルは 0.25 em 以上離す（LP-D104）。",
  "mustPreserve": "今まで重ならずに描けていた図は今まで通り描く（rev 2 の『今まで描けてたものは今まで通り描けるべき』）。重なりも交差もない図は 0.13.0 と同じ位置にラベルを描く。ただし Linear で同じ段の隣のラベルとの間が 0.25 em に満たないときは、そのラベルを別の段に移す。Linear の領域ラベルは、フィーチャーラベルの下線と白い背景を文字と同じによける。よけられないときは、文字にかかるときと同じく、そのフィーチャーラベルを同じ側の別の段へ移す（RGN-YIELD-ROW）。移しても解消しないときは、そのまま描いて警告する。領域ラベルの背景は、リーダー線の上には出るが、下線と白い背景は覆わない（RGN-FRONT）。フィーチャーラベルの守り（Linear は白い背景、Circular はハロー）、下線とその白い縁、リーダー線の太さと止め方は LP-D98 と LP-D106 の既定で描く。Linear の白い背景は、ほかのラベルの文字と下線を覆わない。ダイアログでの選択は History の 1 ステップ、Cancel は記録しない、適用中は操作を受け付けない（PD-OI-088）。",
  "mayRetire": "ラベルが重なったまま図を出す動作。（Linear の領域ラベルの重なりを除く）0.13.0 のフィーチャーラベルの見た目（細いリーダー線、守りと下線なし）を既定で出すこと。labels.stroke_width の設定（Python の config_overrides で指定するとエラーになる）。Linear で同じ段のラベルをすき間なしで並べること。",
  "acceptedResidualRisk": "微調整で全部直る重なり・交差・定義文へのかぶり（欠陥 100 組まで）は、自動で直して描く。直しきれないもの（極の内側ラベルの積み重なり、長い線のもつれ、MjeNMV など）は、重なりなら描かずに直し方を選ばせ、交差やかぶりだけなら描いて帯で知らせる。直し方の候補は、帯のボタンを押したときに計算する。この結果は Session に保存しない。CLI と Python は、描けた図なら 5 秒まで探す。描けない図は、外側ラベルだけの重なりなら値が出るまで探し、内側・埋め込み・放射ラベルの重なりなら 30 秒まで探す。間に合わない候補は『未計算』と示す。Web の候補は 1 つ 30 秒まで探し、見つからなければそう示す。Circular の PDF と EPS はハローの分だけ大きくなる（タバコ葉緑体の PDF は 56 KB から 235 KB、EPS は 1.32 MB から 1.59 MB）。Linear の白い背景では PDF はほとんど変わらず（lambda Linear 24 pt で 17.6 KB から 27 KB）、EPS は半透明の背景を画像で描くので大きくなる（65 KB から 664 KB）。Web の PDF には Circular のハローを描かない（Linear の白い背景、下線とその白い縁は描く。LP-D102）。Web の PDF は文字をカーニングせずに描くので、ラベルが SVG より最大 0.1 em ほど広くなることがある（中心の位置は同じ）。保存済みの Session は次の Generate で新しい見た目になる（それまでは保存した図のまま）。Linear の詰まった段では、ラベルが別の段に移る（lambda 24 pt で 73 のうち 21、タバコ IR 5 pt で 143 のうち 22）。",
  "owner": "satoshikawato",
  "decisionDate": "2026-10-19"
}
```

### PD-OI-098: Font of the diagram text

- Concern key: `diagram-generation.diagram-font`
- Scenario revision: `2`
- Supersedes: `PD-OI-098`, scenario revision `1`. Scenario revision `2`
  keeps the outcome and every field of scenario revision `1` except the Must
  preserve, to whose 〔変えないもの〕 Reply 45 of the Revision 68 entry adds the
  exception for the feature label style and the Web PDF text anchors
  (OV-1303), and the decision date.
- Status: `ACCEPTED`
- Selected outcome: `A+B / CHOSEN-FONT-MEASURED-AND-EMBEDDED`
- Normative outcome: exactly the complete approved scenario revision `2`
  `PRODUCT_DECISION` receipt and its JSON representation below. The
  scenario revision `1` receipt is kept as it was approved.
- Decision source (scenario revision `1`): the Owner reply of `2026-10-13` quoted verbatim in the
  Revision 51 entry above. The Owner approved the receipt text below as
  written. The receipt reproduces its fields without translation or
  additional terms. This record does not supersede another decision. The
  runtime, tests, and documentation that implement it merge with it (static
  Product Contract co-change); this record supplies no runtime acceptance
  evidence.
- Receipt SHA-256 (scenario revision `1`; UTF-8, excluding the final newline):
  `14dc5e981573dce50e9815f4457d34de6786c819ae31a5b4672e178359af3675`.
- Decision source (scenario revision `2`): the Owner replies of `2026-10-19`
  quoted verbatim in the Revision 68 entry above; Reply 45 approves the
  second receipt and JSON below. They keep the fields of scenario revision
  `1`, add the exception of Reply 45 to its Must preserve, and carry the
  decision date of the reply, without translation or additional terms. The
  runtime, tests, and documentation that implement it merge with it (static
  Product Contract co-change); this record supplies no runtime acceptance
  evidence.
- Receipt SHA-256 (scenario revision `2`; UTF-8, excluding the final newline):
  `a305e5d6eefb69914e20e1397cca98298c3050315dff6d488638c2dff5160cc8`.

```text
PRODUCT_DECISION
Concern: diagram-generation.diagram-font
Scenario revision: 1
Choice: A+B / CHOSEN-FONT-MEASURED-AND-EMBEDDED
Rationale: 図の文字の font を、システムにインストールされた font も含めて自由に選べるようにする。Web UI では font file の読み込み（A）とインストール済み font の選択（B）の両方を出す。読み込んだ font は Session に保存せず、gbdraw に含まれない font は警告する。システムにない font を使う Session は読み込むときに警告し、読み込む場合は近い font に置き換える。SVG には PowerPoint などと同じように font を埋め込む。
Must preserve: 〔設定〕図の文字の font は 1 つの設定にする（Python の font_family、CLI の --font_family、Web の Font section の Font）。Web では Circular と Linear の drawing がそれぞれ持つ。既定は Liberation Sans で、font を選ばない図は今と同じ。〔Web での選び方〕Font の欄に Liberation Sans（既定）、Liberation Serif、Liberation Mono、この page で読み込んだ font、Load font file…（.ttf と .otf を選ぶ）を出す。Chrome と Edge では Choose installed font… も出し、インストール済みの font を名前で検索して選べる（ブラウザが許可を求める）。Firefox と Safari では Choose installed font… は出さない。Font は Applies on Generate（Generate で図に入る）。Undo、Redo、Reset Settings は他の設定と同じ。〔測り方と表示〕label の配置、legend、title の文字の幅と高さは、選んだ font そのもので測る。preview も選んだ font で表示する。〔Session〕Session には font の名前（CSS の font-family の並び）だけを保存し、font file は保存しない。Session に入る preview にも font を埋め込まない。〔警告〕gbdraw に含まれない font を選ぶと、Font の欄の下に注意を出す（Session には font の名前だけが残る。別の環境で開くには同じ font が必要で、ないときは近い font に置き換わる）。Web で Session を保存したときも、保存の通知に同じ注意を出す。CLI は --save_session のときに同じ注意を log に出す。〔読み込むとき〕Session の font がその環境にないとき、Web は Session を開く前に dialog で知らせる。ボタンは Load font file…、Use installed font…（Chrome と Edge だけ）、Load with ＜近い font の名前＞、Cancel。font を読み込めばその font で開く。Load with … では Session を開き、文字を近い font で描く。Font の設定には Session の font の名前が残り、あとでその font を読み込めばその font で描く。保存すると元の名前のまま保存する。Cancel では何も変えない。この確認では描画用の Worker を起動しない。CLI と Python は警告を log に出し、近い font で描く。〔近い font〕名前の並びの中で、その環境で使える最初の font。どれも使えないときは、その font の種類（serif、sans-serif、monospace）に合う Liberation（Serif、Sans、Mono）。〔書き出し〕SVG と Interactive SVG は、gbdraw に含まれない font を使っていれば、その font の使った文字だけ（subset）を埋め込む。既定でオンで、Web と CLI のどちらでも切れる。埋め込みを禁止している font（fsType）は埋め込まず、そう知らせる。Liberation は埋め込まない。PNG も同じ埋め込みを使い、選んだ font で描く。PDF も選んだ font で描く。〔変えないもの〕font を選ばない図の文字の配置と出力、Session の保存と読み込みの手順、Gallery の Session。
May retire: Web で font_family を GUI で編集できない設定（読み取り専用の表示）として扱うこと。
Accepted residual risk: 別の環境で同じ font がないと、近い font で描くので文字の幅が変わり、label の配置や legend の大きさが変わる。保存した preview は、font がない環境では Generate するまでブラウザの代わりの font で表示される。Web で読み込んだ font は page を閉じると消え、毎回読み込み直す必要がある。インストール済み font の選択は Chrome と Edge だけ。埋め込んだ SVG と PNG は font の subset のぶん大きくなる。font の license は利用者の責任になる（埋め込みを禁止している font は埋め込まない）。
Owner: satoshikawato
Decision date: 2026-10-13
```

```json
{
  "concern": "diagram-generation.diagram-font",
  "scenarioRevision": 1,
  "choice": "A+B / CHOSEN-FONT-MEASURED-AND-EMBEDDED",
  "rationale": "図の文字の font を、システムにインストールされた font も含めて自由に選べるようにする。Web UI では font file の読み込み（A）とインストール済み font の選択（B）の両方を出す。読み込んだ font は Session に保存せず、gbdraw に含まれない font は警告する。システムにない font を使う Session は読み込むときに警告し、読み込む場合は近い font に置き換える。SVG には PowerPoint などと同じように font を埋め込む。",
  "mustPreserve": "〔設定〕図の文字の font は 1 つの設定にする（Python の font_family、CLI の --font_family、Web の Font section の Font）。Web では Circular と Linear の drawing がそれぞれ持つ。既定は Liberation Sans で、font を選ばない図は今と同じ。〔Web での選び方〕Font の欄に Liberation Sans（既定）、Liberation Serif、Liberation Mono、この page で読み込んだ font、Load font file…（.ttf と .otf を選ぶ）を出す。Chrome と Edge では Choose installed font… も出し、インストール済みの font を名前で検索して選べる（ブラウザが許可を求める）。Firefox と Safari では Choose installed font… は出さない。Font は Applies on Generate（Generate で図に入る）。Undo、Redo、Reset Settings は他の設定と同じ。〔測り方と表示〕label の配置、legend、title の文字の幅と高さは、選んだ font そのもので測る。preview も選んだ font で表示する。〔Session〕Session には font の名前（CSS の font-family の並び）だけを保存し、font file は保存しない。Session に入る preview にも font を埋め込まない。〔警告〕gbdraw に含まれない font を選ぶと、Font の欄の下に注意を出す（Session には font の名前だけが残る。別の環境で開くには同じ font が必要で、ないときは近い font に置き換わる）。Web で Session を保存したときも、保存の通知に同じ注意を出す。CLI は --save_session のときに同じ注意を log に出す。〔読み込むとき〕Session の font がその環境にないとき、Web は Session を開く前に dialog で知らせる。ボタンは Load font file…、Use installed font…（Chrome と Edge だけ）、Load with ＜近い font の名前＞、Cancel。font を読み込めばその font で開く。Load with … では Session を開き、文字を近い font で描く。Font の設定には Session の font の名前が残り、あとでその font を読み込めばその font で描く。保存すると元の名前のまま保存する。Cancel では何も変えない。この確認では描画用の Worker を起動しない。CLI と Python は警告を log に出し、近い font で描く。〔近い font〕名前の並びの中で、その環境で使える最初の font。どれも使えないときは、その font の種類（serif、sans-serif、monospace）に合う Liberation（Serif、Sans、Mono）。〔書き出し〕SVG と Interactive SVG は、gbdraw に含まれない font を使っていれば、その font の使った文字だけ（subset）を埋め込む。既定でオンで、Web と CLI のどちらでも切れる。埋め込みを禁止している font（fsType）は埋め込まず、そう知らせる。Liberation は埋め込まない。PNG も同じ埋め込みを使い、選んだ font で描く。PDF も選んだ font で描く。〔変えないもの〕font を選ばない図の文字の配置と出力、Session の保存と読み込みの手順、Gallery の Session。",
  "mayRetire": "Web で font_family を GUI で編集できない設定（読み取り専用の表示）として扱うこと。",
  "acceptedResidualRisk": "別の環境で同じ font がないと、近い font で描くので文字の幅が変わり、label の配置や legend の大きさが変わる。保存した preview は、font がない環境では Generate するまでブラウザの代わりの font で表示される。Web で読み込んだ font は page を閉じると消え、毎回読み込み直す必要がある。インストール済み font の選択は Chrome と Edge だけ。埋め込んだ SVG と PNG は font の subset のぶん大きくなる。font の license は利用者の責任になる（埋め込みを禁止している font は埋め込まない）。",
  "owner": "satoshikawato",
  "decisionDate": "2026-10-13"
}
```

Scenario revision `2` receipt, approved by `satoshikawato` on `2026-10-19` (Reply 45 of the Revision 68 entry):

```text
PRODUCT_DECISION
Concern: diagram-generation.diagram-font
Scenario revision: 2
Choice: A+B / CHOSEN-FONT-MEASURED-AND-EMBEDDED
Rationale: 図の文字の font を、システムにインストールされた font も含めて自由に選べるようにする。Web UI では font file の読み込み（A）とインストール済み font の選択（B）の両方を出す。読み込んだ font は Session に保存せず、gbdraw に含まれない font は警告する。システムにない font を使う Session は読み込むときに警告し、読み込む場合は近い font に置き換える。SVG には PowerPoint などと同じように font を埋め込む。
Must preserve: 〔設定〕図の文字の font は 1 つの設定にする（Python の font_family、CLI の --font_family、Web の Font section の Font）。Web では Circular と Linear の drawing がそれぞれ持つ。既定は Liberation Sans で、font を選ばない図は今と同じ。〔Web での選び方〕Font の欄に Liberation Sans（既定）、Liberation Serif、Liberation Mono、この page で読み込んだ font、Load font file…（.ttf と .otf を選ぶ）を出す。Chrome と Edge では Choose installed font… も出し、インストール済みの font を名前で検索して選べる（ブラウザが許可を求める）。Firefox と Safari では Choose installed font… は出さない。Font は Applies on Generate（Generate で図に入る）。Undo、Redo、Reset Settings は他の設定と同じ。〔測り方と表示〕label の配置、legend、title の文字の幅と高さは、選んだ font そのもので測る。preview も選んだ font で表示する。〔Session〕Session には font の名前（CSS の font-family の並び）だけを保存し、font file は保存しない。Session に入る preview にも font を埋め込まない。〔警告〕gbdraw に含まれない font を選ぶと、Font の欄の下に注意を出す（Session には font の名前だけが残る。別の環境で開くには同じ font が必要で、ないときは近い font に置き換わる）。Web で Session を保存したときも、保存の通知に同じ注意を出す。CLI は --save_session のときに同じ注意を log に出す。〔読み込むとき〕Session の font がその環境にないとき、Web は Session を開く前に dialog で知らせる。ボタンは Load font file…、Use installed font…（Chrome と Edge だけ）、Load with ＜近い font の名前＞、Cancel。font を読み込めばその font で開く。Load with … では Session を開き、文字を近い font で描く。Font の設定には Session の font の名前が残り、あとでその font を読み込めばその font で描く。保存すると元の名前のまま保存する。Cancel では何も変えない。この確認では描画用の Worker を起動しない。CLI と Python は警告を log に出し、近い font で描く。〔近い font〕名前の並びの中で、その環境で使える最初の font。どれも使えないときは、その font の種類（serif、sans-serif、monospace）に合う Liberation（Serif、Sans、Mono）。〔書き出し〕SVG と Interactive SVG は、gbdraw に含まれない font を使っていれば、その font の使った文字だけ（subset）を埋め込む。既定でオンで、Web と CLI のどちらでも切れる。埋め込みを禁止している font（fsType）は埋め込まず、そう知らせる。Liberation は埋め込まない。PNG も同じ埋め込みを使い、選んだ font で描く。PDF も選んだ font で描く。〔変えないもの〕font を選ばない図の文字の配置と出力（ただし、フィーチャーラベルのハローと白い背景、下線とその白い縁、リーダー線の太さ、Linear の同じ段のすき間（LP-D98、LP-D100、LP-D102〜LP-D104、LP-D106）と、Web の PDF で中央揃えと右揃えの文字の位置を直したこと（OV-1303）は除く）、Session の保存と読み込みの手順、Gallery の Session。
May retire: Web で font_family を GUI で編集できない設定（読み取り専用の表示）として扱うこと。
Accepted residual risk: 別の環境で同じ font がないと、近い font で描くので文字の幅が変わり、label の配置や legend の大きさが変わる。保存した preview は、font がない環境では Generate するまでブラウザの代わりの font で表示される。Web で読み込んだ font は page を閉じると消え、毎回読み込み直す必要がある。インストール済み font の選択は Chrome と Edge だけ。埋め込んだ SVG と PNG は font の subset のぶん大きくなる。font の license は利用者の責任になる（埋め込みを禁止している font は埋め込まない）。
Owner: satoshikawato
Decision date: 2026-10-19
```

```json
{
  "concern": "diagram-generation.diagram-font",
  "scenarioRevision": 2,
  "choice": "A+B / CHOSEN-FONT-MEASURED-AND-EMBEDDED",
  "rationale": "図の文字の font を、システムにインストールされた font も含めて自由に選べるようにする。Web UI では font file の読み込み（A）とインストール済み font の選択（B）の両方を出す。読み込んだ font は Session に保存せず、gbdraw に含まれない font は警告する。システムにない font を使う Session は読み込むときに警告し、読み込む場合は近い font に置き換える。SVG には PowerPoint などと同じように font を埋め込む。",
  "mustPreserve": "〔設定〕図の文字の font は 1 つの設定にする（Python の font_family、CLI の --font_family、Web の Font section の Font）。Web では Circular と Linear の drawing がそれぞれ持つ。既定は Liberation Sans で、font を選ばない図は今と同じ。〔Web での選び方〕Font の欄に Liberation Sans（既定）、Liberation Serif、Liberation Mono、この page で読み込んだ font、Load font file…（.ttf と .otf を選ぶ）を出す。Chrome と Edge では Choose installed font… も出し、インストール済みの font を名前で検索して選べる（ブラウザが許可を求める）。Firefox と Safari では Choose installed font… は出さない。Font は Applies on Generate（Generate で図に入る）。Undo、Redo、Reset Settings は他の設定と同じ。〔測り方と表示〕label の配置、legend、title の文字の幅と高さは、選んだ font そのもので測る。preview も選んだ font で表示する。〔Session〕Session には font の名前（CSS の font-family の並び）だけを保存し、font file は保存しない。Session に入る preview にも font を埋め込まない。〔警告〕gbdraw に含まれない font を選ぶと、Font の欄の下に注意を出す（Session には font の名前だけが残る。別の環境で開くには同じ font が必要で、ないときは近い font に置き換わる）。Web で Session を保存したときも、保存の通知に同じ注意を出す。CLI は --save_session のときに同じ注意を log に出す。〔読み込むとき〕Session の font がその環境にないとき、Web は Session を開く前に dialog で知らせる。ボタンは Load font file…、Use installed font…（Chrome と Edge だけ）、Load with ＜近い font の名前＞、Cancel。font を読み込めばその font で開く。Load with … では Session を開き、文字を近い font で描く。Font の設定には Session の font の名前が残り、あとでその font を読み込めばその font で描く。保存すると元の名前のまま保存する。Cancel では何も変えない。この確認では描画用の Worker を起動しない。CLI と Python は警告を log に出し、近い font で描く。〔近い font〕名前の並びの中で、その環境で使える最初の font。どれも使えないときは、その font の種類（serif、sans-serif、monospace）に合う Liberation（Serif、Sans、Mono）。〔書き出し〕SVG と Interactive SVG は、gbdraw に含まれない font を使っていれば、その font の使った文字だけ（subset）を埋め込む。既定でオンで、Web と CLI のどちらでも切れる。埋め込みを禁止している font（fsType）は埋め込まず、そう知らせる。Liberation は埋め込まない。PNG も同じ埋め込みを使い、選んだ font で描く。PDF も選んだ font で描く。〔変えないもの〕font を選ばない図の文字の配置と出力（ただし、フィーチャーラベルのハローと白い背景、下線とその白い縁、リーダー線の太さ、Linear の同じ段のすき間（LP-D98、LP-D100、LP-D102〜LP-D104、LP-D106）と、Web の PDF で中央揃えと右揃えの文字の位置を直したこと（OV-1303）は除く）、Session の保存と読み込みの手順、Gallery の Session。",
  "mayRetire": "Web で font_family を GUI で編集できない設定（読み取り専用の表示）として扱うこと。",
  "acceptedResidualRisk": "別の環境で同じ font がないと、近い font で描くので文字の幅が変わり、label の配置や legend の大きさが変わる。保存した preview は、font がない環境では Generate するまでブラウザの代わりの font で表示される。Web で読み込んだ font は page を閉じると消え、毎回読み込み直す必要がある。インストール済み font の選択は Chrome と Edge だけ。埋め込んだ SVG と PNG は font の subset のぶん大きくなる。font の license は利用者の責任になる（埋め込みを禁止している font は埋め込まない）。",
  "owner": "satoshikawato",
  "decisionDate": "2026-10-19"
}
```

### PD-OI-099: Confirm what Load cannot restore before it loads

- Concern key: `web.session.load-repair-confirmation`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `A / CONFIRM-REPAIRS-BEFORE-LOAD`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its JSON representation below.
- Decision source: the Owner replies quoted verbatim in the Revision 52 entry
  above. The Owner approved the receipt text below in Reply 3 of that entry.
  The receipt reproduces its fields without translation or additional terms.
  Dependent runtime requires this authority merged into its base, or merged
  together with it through the reviewed static Product Contract co-change
  route.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `d646d6704e71e3066315abdba14cf762e847653a7041982359bdb98711452245`.
- Acceptance contracts: `OIC-031`. These obligations and the complete selected
  outcome are jointly required; their citation does not claim completed
  dependent-runtime checks.

```text
PRODUCT_DECISION
Concern: web.session.load-repair-confirmation
Scenario revision: 1
Choice: A / CONFIRM-REPAIRS-BEFORE-LOAD
Rationale: 古いファイルや、バグなどで不正な値を含むファイルを拒否すると、ユーザーは何も取り戻せない。拒否するよりも、読み込めるところまで読み込み、一部を完全には復元できないことを読み込む前に一言知らせ、次に保存すると最新の形式になると伝える。
Must preserve: 〔拒否〕次のファイルだけは今までどおり読み込まずにエラーを示す。JSON として読めないファイル（壊れた圧縮ファイル、UTF-8 でない文字を含むファイルを含む）。gbdraw の Session でも設定ファイルでもないファイル（開発版だけの Session 34〜38、43、45 と、開発版だけの request の形式 3、4 を含む）。版 26 以下の Session（gbdraw 0.11.0 以前と、0.12.0 より前の gbdraw.app が保存したもの）。このときはエラーで古い版だと示す。新しい版の Session や request の形式。大きすぎるファイル（200 MB、展開して 512 MB を超えるもの）。読み込みに要る機能（gzip の展開、module Worker）がない browser で開いたファイル。安全でない名前（__proto__、constructor、prototype）を含むファイル。何も復元できないファイル（ファイルにある設定、request、Result、入力ファイル、編集、cache のすべてを捨てることになるもの）。読み込みの途中で browser の機能や Worker が失敗したときも、読み込まずにエラーを示す。〔確認〕それ以外の Session と古い設定ファイルは、読み込めるところまで読む。値を捨てる、既定値にする、意味を変えて読むところが 1 つでもあれば、読み込む前に 1 つの dialog でその一覧（何をどう読んだか。設定の値の行はファイルの中の場所も示す）を示し、次に Save Session すると最新の形式で保存されることを伝える。[Load] で読み込み、[Cancel] では設定、Result、入力ファイル、History、画面を何も変えない。ファイル自体は変えない。〔出さないもの〕古い版の形式からの変換や、大文字、前後の空白、別名の違いのように、意味を変えずに読めるものでは dialog を出さない。計算し直せる cache の項目、同じ Legend の行に後の編集が上書きしていた Session 46 の色や線の編集、Result にない group の名前や説明の編集、Instant Preview が on のときの当てていない palette の選択、画面の表示の状態、もう一方の図の型（Circular か Linear）だけで使う設定、ファイルの中の別の有効な値で置き換わる値のように、図に効いていなかった値を捨てるだけのときも dialog を出さない。何も捨てずに読めるファイル（Gallery の Session、現在の Web が保存した Session を含む）では今までどおり dialog を出さない。〔今までの知らせ〕今は読み込んだ後の知らせに出ている、合わせられずに捨てた feature の編集、範囲を狭めた Feature visibility の編集、hash= を付け替えられなかった値、読み替えたり捨てたりした表の行、色を受け付けずに捨てた Legend の項目は、この dialog に移す。意味を変えない変換の知らせ（annotation の hash= を source で指す形に変えた数）と、古い設定ファイルを読んだという知らせは、今までどおり読み込んだ後に出す。Linear の Comparisons 欄の、保存された比較を読み取り専用で残した、または決定が要るという枠（選ぶ操作つき）は、今までどおり読み込んだ後に出す。feature の情報を復元できないときは、dialog の一覧に示し、Features 欄の表示（SVG の preview は残る）も今までどおり出す。〔ほかの確認〕読み込む前の「Replace the current work?」の確認は今までどおり先に出す。〔読み込んだ後〕Session を読み込んだ状態が Undo の起点になり、変更のない状態として始まる。保存すると直した値が最新の形式で書かれる。〔Python と CLI〕CLI と Python API は、request の形式が定めない値（文字や真偽値や名前の欄の null や別の種類の値、形式 9 では定義にない key や名前の別の綴り）、比較の表のファイルの path（以前の Python API が一部の request に書いたもの。今は表を resource として保存する）、ASCII でない数字を使った record の selector を拒否する。それ以外に Session を読む範囲は変えない（古い形式の request の大文字や別名などは今までどおり読む）。
May retire: Web が、上の〔拒否〕に当たらない Session（古い Session、不正な値を含む Session、Python は読める Session）と古い設定ファイルを、読み込みのときに拒否する動作。完全には復元できなかったことを、読み込んだ後の知らせで伝える動作。知らせずに捨てたり既定値にしたりする動作（凡例の位置を既定値にする今の暫定の動作を含む）。
Accepted residual risk: 直した部分は保存された Result（図）と合わないことがあり、次の Generate で直した設定の図になる。dialog の一覧には、ファイルの中の場所（modes.circular.config.adv.depth_tracks[0].height など）をそのまま示す行がある。一覧が長いときは 12 行までを示し、残りの数を示す。Web では dialog を出して読めるが、CLI では拒否される Session がある（手で書き換えた Session の不正な値、request の形式が定めない値）。request を直せない Session は、設定だけの Session として読み、保存された Result、GenBank、GFF、FASTA の入力ファイル、feature の一覧、cache を捨てて一覧に示すことがある。この動作はいくつかの変更に分けて入る。最初の変更（この dialog）より後の変更が入るまでは、次のことが残る。古い版の Session や古い設定ファイルの一部の値や不正な値（track の並びの設定の値、受け付けない Legend の色の編集、設定ファイルの図の型に合わない Legend の位置や知らない key など）は、今までどおり dialog に示さずに直されるか、そのまま読み込まれる（知らない track の型など）。不正な値を含む Session や設定ファイルの一部（定義にない key や欄、知らない値、図の型にない値、範囲外の値、種類の違う値、足りない resource や欄、読めない track の並びや cache、確認を通らない保存された Result を含むもの）は、今までどおり拒否される。何も復元できないファイルを拒否する確認はまだなく、後の変更で入る。保存された request に形式が定めない値を含む Session のうち、今までの Web が読めたものは、その値を dialog に示さずに読み込まれ、Generate や、保存された request で描き直す操作（label の配置し直し、alignment、record の変換）で「An input value is invalid」と出る。それ以外（feature の編集の表や grouping が null のもの、layout に定義にない key があるものなど）は今までどおり拒否される。
Owner: satoshikawato
Decision date: 2026-10-14
```

```json
{
  "concern": "web.session.load-repair-confirmation",
  "scenarioRevision": 1,
  "choice": "A / CONFIRM-REPAIRS-BEFORE-LOAD",
  "rationale": "古いファイルや、バグなどで不正な値を含むファイルを拒否すると、ユーザーは何も取り戻せない。拒否するよりも、読み込めるところまで読み込み、一部を完全には復元できないことを読み込む前に一言知らせ、次に保存すると最新の形式になると伝える。",
  "mustPreserve": "〔拒否〕次のファイルだけは今までどおり読み込まずにエラーを示す。JSON として読めないファイル（壊れた圧縮ファイル、UTF-8 でない文字を含むファイルを含む）。gbdraw の Session でも設定ファイルでもないファイル（開発版だけの Session 34〜38、43、45 と、開発版だけの request の形式 3、4 を含む）。版 26 以下の Session（gbdraw 0.11.0 以前と、0.12.0 より前の gbdraw.app が保存したもの）。このときはエラーで古い版だと示す。新しい版の Session や request の形式。大きすぎるファイル（200 MB、展開して 512 MB を超えるもの）。読み込みに要る機能（gzip の展開、module Worker）がない browser で開いたファイル。安全でない名前（__proto__、constructor、prototype）を含むファイル。何も復元できないファイル（ファイルにある設定、request、Result、入力ファイル、編集、cache のすべてを捨てることになるもの）。読み込みの途中で browser の機能や Worker が失敗したときも、読み込まずにエラーを示す。〔確認〕それ以外の Session と古い設定ファイルは、読み込めるところまで読む。値を捨てる、既定値にする、意味を変えて読むところが 1 つでもあれば、読み込む前に 1 つの dialog でその一覧（何をどう読んだか。設定の値の行はファイルの中の場所も示す）を示し、次に Save Session すると最新の形式で保存されることを伝える。[Load] で読み込み、[Cancel] では設定、Result、入力ファイル、History、画面を何も変えない。ファイル自体は変えない。〔出さないもの〕古い版の形式からの変換や、大文字、前後の空白、別名の違いのように、意味を変えずに読めるものでは dialog を出さない。計算し直せる cache の項目、同じ Legend の行に後の編集が上書きしていた Session 46 の色や線の編集、Result にない group の名前や説明の編集、Instant Preview が on のときの当てていない palette の選択、画面の表示の状態、もう一方の図の型（Circular か Linear）だけで使う設定、ファイルの中の別の有効な値で置き換わる値のように、図に効いていなかった値を捨てるだけのときも dialog を出さない。何も捨てずに読めるファイル（Gallery の Session、現在の Web が保存した Session を含む）では今までどおり dialog を出さない。〔今までの知らせ〕今は読み込んだ後の知らせに出ている、合わせられずに捨てた feature の編集、範囲を狭めた Feature visibility の編集、hash= を付け替えられなかった値、読み替えたり捨てたりした表の行、色を受け付けずに捨てた Legend の項目は、この dialog に移す。意味を変えない変換の知らせ（annotation の hash= を source で指す形に変えた数）と、古い設定ファイルを読んだという知らせは、今までどおり読み込んだ後に出す。Linear の Comparisons 欄の、保存された比較を読み取り専用で残した、または決定が要るという枠（選ぶ操作つき）は、今までどおり読み込んだ後に出す。feature の情報を復元できないときは、dialog の一覧に示し、Features 欄の表示（SVG の preview は残る）も今までどおり出す。〔ほかの確認〕読み込む前の「Replace the current work?」の確認は今までどおり先に出す。〔読み込んだ後〕Session を読み込んだ状態が Undo の起点になり、変更のない状態として始まる。保存すると直した値が最新の形式で書かれる。〔Python と CLI〕CLI と Python API は、request の形式が定めない値（文字や真偽値や名前の欄の null や別の種類の値、形式 9 では定義にない key や名前の別の綴り）、比較の表のファイルの path（以前の Python API が一部の request に書いたもの。今は表を resource として保存する）、ASCII でない数字を使った record の selector を拒否する。それ以外に Session を読む範囲は変えない（古い形式の request の大文字や別名などは今までどおり読む）。",
  "mayRetire": "Web が、上の〔拒否〕に当たらない Session（古い Session、不正な値を含む Session、Python は読める Session）と古い設定ファイルを、読み込みのときに拒否する動作。完全には復元できなかったことを、読み込んだ後の知らせで伝える動作。知らせずに捨てたり既定値にしたりする動作（凡例の位置を既定値にする今の暫定の動作を含む）。",
  "acceptedResidualRisk": "直した部分は保存された Result（図）と合わないことがあり、次の Generate で直した設定の図になる。dialog の一覧には、ファイルの中の場所（modes.circular.config.adv.depth_tracks[0].height など）をそのまま示す行がある。一覧が長いときは 12 行までを示し、残りの数を示す。Web では dialog を出して読めるが、CLI では拒否される Session がある（手で書き換えた Session の不正な値、request の形式が定めない値）。request を直せない Session は、設定だけの Session として読み、保存された Result、GenBank、GFF、FASTA の入力ファイル、feature の一覧、cache を捨てて一覧に示すことがある。この動作はいくつかの変更に分けて入る。最初の変更（この dialog）より後の変更が入るまでは、次のことが残る。古い版の Session や古い設定ファイルの一部の値や不正な値（track の並びの設定の値、受け付けない Legend の色の編集、設定ファイルの図の型に合わない Legend の位置や知らない key など）は、今までどおり dialog に示さずに直されるか、そのまま読み込まれる（知らない track の型など）。不正な値を含む Session や設定ファイルの一部（定義にない key や欄、知らない値、図の型にない値、範囲外の値、種類の違う値、足りない resource や欄、読めない track の並びや cache、確認を通らない保存された Result を含むもの）は、今までどおり拒否される。何も復元できないファイルを拒否する確認はまだなく、後の変更で入る。保存された request に形式が定めない値を含む Session のうち、今までの Web が読めたものは、その値を dialog に示さずに読み込まれ、Generate や、保存された request で描き直す操作（label の配置し直し、alignment、record の変換）で「An input value is invalid」と出る。それ以外（feature の編集の表や grouping が null のもの、layout に定義にない key があるものなど）は今までどおり拒否される。",
  "owner": "satoshikawato",
  "decisionDate": "2026-10-14"
}
```

### PD-OI-100: Rich Feature Popup outside Undo

- Concern key: `web.history.rich-feature-popup`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `A / NO-UNDO-STEP`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its JSON representation below.
- Decision source: the Owner replies quoted verbatim in the Revision 54 entry
  above (Reply 1 and Replies (a) to (d)); `(a)` and `(b)` in the receipt name
  those replies. Fields for which the Owner gave nothing read `not stated`.
  The runtime, tests, and documentation that implement it merge with it
  (static Product Contract co-change); this record supplies no runtime
  acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `35f45f5ae5bfbf8ddef85d1649f4e518f90c4610a0ac04a8f0bde2a01a1f49d4`.

```text
PRODUCT_DECISION
Concern: web.history.rich-feature-popup
Scenario revision: 1
Choice: A / NO-UNDO-STEP
Rationale: not stated
Must preserve: (a) 「どちらも今の値のまま」 (b) 「確認を出す」
May retire: not stated
Accepted residual risk: not stated
Owner: satoshikawato
Decision date: 2026-10-16
```

```json
{
  "concern": "web.history.rich-feature-popup",
  "scenarioRevision": 1,
  "choice": "A / NO-UNDO-STEP",
  "rationale": "not stated",
  "mustPreserve": "(a) 「どちらも今の値のまま」 (b) 「確認を出す」",
  "mayRetire": "not stated",
  "acceptedResidualRisk": "not stated",
  "owner": "satoshikawato",
  "decisionDate": "2026-10-16"
}
```

### PD-OI-101: Task yield after the later Generate statuses

- Concern key: `web.generate.later-status-paint-yields`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `O-B3 / TASK-YIELD`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its JSON representation below.
- Decision source: the Owner reply quoted verbatim in the Revision 56 entry
  above. The runtime and tests that implement it merge with it (static
  Product Contract co-change); this record supplies no runtime acceptance
  evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `7f295b16535ec23ede9fdafca1391754fce7f7db250617ce3679b1e563cb658c`.

```text
PRODUCT_DECISION
Concern: web.generate.later-status-paint-yields
Scenario revision: 1
Choice: O-B3 / TASK-YIELD
Rationale: warm Generates wait up to 21-50 ms of frame waits for text that lives about one frame
Must preserve: overlay and first status paint; Worker progress texts; every status value in order; `generate.paint-opportunity-completed` before capture; GE-09 and the pre-commit Cancel in Chromium
May retire: the paint opportunity of "Preparing render inputs and session..." and "Preparing preview..."
Accepted residual risk: texts missed on fast Generates; in Firefox or WebKit a queued Cancel may be dispatched after the continuation, giving B1's consequence (cold start after a Cancel in a 42-69 ms window) or B2's (a candidate mounted then restored) in that window
Owner: satoshikawato
Decision date: 2026-10-16
```

```json
{
  "concern": "web.generate.later-status-paint-yields",
  "scenarioRevision": 1,
  "choice": "O-B3 / TASK-YIELD",
  "rationale": "warm Generates wait up to 21-50 ms of frame waits for text that lives about one frame",
  "mustPreserve": "overlay and first status paint; Worker progress texts; every status value in order; `generate.paint-opportunity-completed` before capture; GE-09 and the pre-commit Cancel in Chromium",
  "mayRetire": "the paint opportunity of \"Preparing render inputs and session...\" and \"Preparing preview...\"",
  "acceptedResidualRisk": "texts missed on fast Generates; in Firefox or WebKit a queued Cancel may be dispatched after the continuation, giving B1's consequence (cold start after a Cancel in a 42-69 ms window) or B2's (a candidate mounted then restored) in that window",
  "owner": "satoshikawato",
  "decisionDate": "2026-10-16"
}
```

### PD-OI-102: Live settings redraw the Result when a value is committed

- Concern key: `web.settings.live-application`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `A / LIVE-ON-COMMIT`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its nine-field JSON representation below.
- Decision source: the Owner replies quoted verbatim in the Revision 57 entry
  above. The receipt and JSON below reproduce all nine confirmed fields
  without translation or additional terms. This record does not supersede
  another decision. Dependent runtime requires this authority merged into its
  base; this record supplies no runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `3593074885d0092e44834791a2cfaa3de11d88001352203542149c594e81a599`.

```text
PRODUCT_DECISION
Concern: web.settings.live-application
Scenario revision: 1
Choice: A / LIVE-ON-COMMIT
Rationale: the campaign goal (every setting the Result shows applies in real time) with one mechanism: committing a value redraws the current Result in the background through the same renderer and stages as Generate.
Must preserve: Generate (atomic, applies every other setting, resets zoom); Processing / Canceling only in the Generate button; the live edit notices; the previous Result on failure; Undo / Redo (one step per change); the Session Result / draft separation for settings that apply on Generate; export of the current Result; every existing live edit; a live Result identical to a Generate.
May retire: the Generate-only timing of the sections chosen in S4-Q2 (G1); the "Applies on Generate" note on those controls; staging several of those changes to see them only at Generate.
Accepted residual risk: Save, Load, a mode switch and the Undo/Redo of a mode switch are refused during an apply; editor edits wait until the apply loop is idle (about two applies, longer if values keep being committed); Reset Settings or Import settings show a half-applied state until Generate; a valid value coalesced with an invalid one is not applied; an Undo whose re-apply fails leaves the control and the figure on different values.
Owner: Satoshi Kawato
Decision date: 2026-10-17
```

```json
{
  "concern": "web.settings.live-application",
  "scenarioRevision": 1,
  "choice": "A / LIVE-ON-COMMIT",
  "rationale": "the campaign goal (every setting the Result shows applies in real time) with one mechanism: committing a value redraws the current Result in the background through the same renderer and stages as Generate.",
  "mustPreserve": "Generate (atomic, applies every other setting, resets zoom); Processing / Canceling only in the Generate button; the live edit notices; the previous Result on failure; Undo / Redo (one step per change); the Session Result / draft separation for settings that apply on Generate; export of the current Result; every existing live edit; a live Result identical to a Generate.",
  "mayRetire": "the Generate-only timing of the sections chosen in S4-Q2 (G1); the \"Applies on Generate\" note on those controls; staging several of those changes to see them only at Generate.",
  "acceptedResidualRisk": "Save, Load, a mode switch and the Undo/Redo of a mode switch are refused during an apply; editor edits wait until the apply loop is idle (about two applies, longer if values keep being committed); Reset Settings or Import settings show a half-applied state until Generate; a valid value coalesced with an invalid one is not applied; an Undo whose re-apply fails leaves the control and the figure on different values.",
  "owner": "Satoshi Kawato",
  "decisionDate": "2026-10-17"
}
```

### PD-OI-103: Settings sections that apply live

- Concern key: `web.settings.live-sections`
- Scenario revision: `2`
- Supersedes: `PD-OI-103`, scenario revision `1` (`A / G1-LIVE`). Scenario
  revision `2` changes only what its receipt states, for the G2 units and GT;
  in the decision packet's words: "the rule gates G2 units and GT only; G1
  stays live under PD-OI-103 rev 1, whose Accepted residual risk covers Plot
  Title over P1 (+54-85 ms). Under A it is a gate: PD-OI-103 rev 1's question
  trigger becomes a go-live."
- Status: `ACCEPTED`
- Selected outcomes: `A / G1-LIVE` (scenario revision `1`) and
  `A / LIVE-BY-EVIDENCE-RULE`, reading `MARGIN`, unit `U2 CONTROL GROUPS`
  (scenario revision `2`)
- Normative outcome: exactly the complete approved scenario revision `2`
  `PRODUCT_DECISION` receipt and its nine-field JSON representation below,
  with the effects of the scenario revision `1` receipt that it does not
  change. The scenario revision `1` receipt is kept as it was approved.
- Decision source (scenario revision `1`): the Owner replies quoted verbatim
  in the Revision 57 entry above. The receipt and JSON below reproduce all
  nine confirmed fields without translation or additional terms.
- Receipt SHA-256 (scenario revision `1`; UTF-8, excluding the final newline):
  `d30daa401717a284efa14fd6e9ba211617c3e51af15499371127b2bc71e87fe2`.
- Decision source (scenario revision `2`): the Owner replies quoted
  verbatim in the Revision 64 entry above. The scenario revision `2`
  receipt and JSON below reproduce all nine fields without translation or
  additional terms. Dependent runtime requires this authority merged into its
  base; this record supplies no runtime acceptance evidence.
- Receipt SHA-256 (scenario revision `2`; UTF-8, excluding the final newline):
  `95f13d32a53bccd0b4997663e78ff5acb731971220cc3fa36e8ba426267ab6ca`.
- Setting groups the receipts name (G1, G2, G3, GT): the membership quoted
  from the packet's "Groups used below". It names the groups; it adds no
  outcome.
  - G1: Legend (position, swatch size, font size); Plot Title, its Position
    and Font Size (both modes); Circular multi-record Record Size Mode, Min
    Radius, Column and Row Gap Ratio, Row for record; Linear Lock Definition
    Column.
  - G2: Species, Strain, Default font size, definition rows, Keep Full
    Definition, record names; Axis; Dinucleotide; Region Annotations;
    comparison match height and style; Linear record gap and rows; feature
    block and line stroke; custom track stacks.
  - GT: Show / Hide GC Content, GC Skew, Depth.
  - G3 (stays on Generate in every choice): searches, record set, comparison
    filters, Labels, Font, feature types, rendering, geometry and placement,
    layout presets, reverse, region, display start, depth and Circular
    conservation display settings.

Scenario revision `1` receipt, confirmed on `2026-10-17`:

```text
PRODUCT_DECISION
Concern: web.settings.live-sections
Scenario revision: 1
Choice: A / G1-LIVE
Rationale: G1 applies in 0.25-0.45 s on every measured Session; G2 costs 0.8-6.3 s on large Sessions, so each G2 section is asked later on its own evidence (the packet's later-question rule).
Must preserve: every G2 setting's Generate-only timing and the Session separation for it; G3 unchanged; GT per S4-Q3.
May retire: the Generate-only timing and the note of Legend, Plot Title and its Position and Font Size, the Circular multi-record layout and the Linear Lock Definition Column, in both modes.
Accepted residual risk: Plot Title +55 / +85 / +54 ms over 250 ms on Vnig / hepa / vibrio, Lock Definition Column about +30 ms, tobacco +138 / +157 ms over 100 ms; the S4-Q1 risks for G1.
Owner: Satoshi Kawato
Decision date: 2026-10-17
```

```json
{
  "concern": "web.settings.live-sections",
  "scenarioRevision": 1,
  "choice": "A / G1-LIVE",
  "rationale": "G1 applies in 0.25-0.45 s on every measured Session; G2 costs 0.8-6.3 s on large Sessions, so each G2 section is asked later on its own evidence (the packet's later-question rule).",
  "mustPreserve": "every G2 setting's Generate-only timing and the Session separation for it; G3 unchanged; GT per S4-Q3.",
  "mayRetire": "the Generate-only timing and the note of Legend, Plot Title and its Position and Font Size, the Circular multi-record layout and the Linear Lock Definition Column, in both modes.",
  "acceptedResidualRisk": "Plot Title +55 / +85 / +54 ms over 250 ms on Vnig / hepa / vibrio, Lock Definition Column about +30 ms, tobacco +138 / +157 ms over 100 ms; the S4-Q1 risks for G1.",
  "owner": "Satoshi Kawato",
  "decisionDate": "2026-10-17"
}
```

Scenario revision `2` receipt, confirmed on `2026-10-18`:

```text
PRODUCT_DECISION
Concern: web.settings.live-sections
Scenario revision: 2
Choice: A / LIVE-BY-EVIDENCE-RULE; reading MARGIN, unit U2 CONTROL GROUPS
Rationale: D-116's order ("first those meeting the targets") in one receipt; MARGIN because the gate has no Owner round, so a unit near its target comes back to the Owner; cost: Axis stays on Generate unless stage work brings it under 192 ms. U2: U1 keeps a fast control (window, hepa 535 ms) behind the slowest; U2 adds no owner.
Must preserve: every G2 setting's Generate-only timing and Session separation until its unit is live; G1, G3; GT per S5-Q2. Evidence rule: "A unit (a G2 unit, or the GT set or toggle as chosen in S5-Q2) is within the rule when its acceptance report shows, with the live mode of the step 4 spec (controls selected by CSS), the unit's representative change (U1: the section's slowest control; U2: each control group's representative change; GT: hide and show of each toggle), warm, at least 3 interleaved rounds with alternating Session order, each run started at a 1-minute load at or below 8 and a 5-minute load at or below 6, that in every round the median live apply (MARGIN: multiplied by 1.3; SENSITIVITY: itself, the x 1.3 result reported) is at or below its D-03 class target (P1 250, P2 400, P3 1000 ms) on each deciding Session. The deciding Sessions are the Gallery Vnig, hepatoplasmataceae orthogroup and vibrio collinear Sessions where the control is enabled and its change draws something different, with these cases fixed: hepa's record gap (one record) is left out; record names are judged on hepa and vibrio (Vnig's Multi-Record Canvas disables them); Dinucleotide controls on vibrio with Show GC Content turned on first; a toggle a Session has off is turned on first for its hide; Show Depth on all three Sessions, each with an added per-base Depth TSV (`samtools depth -a` text: one row per position of every record). The tobacco result, the stage D state and the stage work are reported. After SW1 is measured (and SW2 / SW3 where designed), every unit still over returns to the Owner as a 'stay on Generate (GT: stay instant) / accept the measured cost' question with its numbers."
May retire: per unit at go-live, its Generate-only timing, its note, and staging its changes to see them only at Generate.
Accepted residual risk: units may wait for the return question; MARGIN may hold back a unit under target (Axis).
Owner: satoshikawato
Decision date: 2026-10-18
```

```json
{
  "concern": "web.settings.live-sections",
  "scenarioRevision": 2,
  "choice": "A / LIVE-BY-EVIDENCE-RULE; reading MARGIN, unit U2 CONTROL GROUPS",
  "rationale": "D-116's order (\"first those meeting the targets\") in one receipt; MARGIN because the gate has no Owner round, so a unit near its target comes back to the Owner; cost: Axis stays on Generate unless stage work brings it under 192 ms. U2: U1 keeps a fast control (window, hepa 535 ms) behind the slowest; U2 adds no owner.",
  "mustPreserve": "every G2 setting's Generate-only timing and Session separation until its unit is live; G1, G3; GT per S5-Q2. Evidence rule: \"A unit (a G2 unit, or the GT set or toggle as chosen in S5-Q2) is within the rule when its acceptance report shows, with the live mode of the step 4 spec (controls selected by CSS), the unit's representative change (U1: the section's slowest control; U2: each control group's representative change; GT: hide and show of each toggle), warm, at least 3 interleaved rounds with alternating Session order, each run started at a 1-minute load at or below 8 and a 5-minute load at or below 6, that in every round the median live apply (MARGIN: multiplied by 1.3; SENSITIVITY: itself, the x 1.3 result reported) is at or below its D-03 class target (P1 250, P2 400, P3 1000 ms) on each deciding Session. The deciding Sessions are the Gallery Vnig, hepatoplasmataceae orthogroup and vibrio collinear Sessions where the control is enabled and its change draws something different, with these cases fixed: hepa's record gap (one record) is left out; record names are judged on hepa and vibrio (Vnig's Multi-Record Canvas disables them); Dinucleotide controls on vibrio with Show GC Content turned on first; a toggle a Session has off is turned on first for its hide; Show Depth on all three Sessions, each with an added per-base Depth TSV (`samtools depth -a` text: one row per position of every record). The tobacco result, the stage D state and the stage work are reported. After SW1 is measured (and SW2 / SW3 where designed), every unit still over returns to the Owner as a 'stay on Generate (GT: stay instant) / accept the measured cost' question with its numbers.\"",
  "mayRetire": "per unit at go-live, its Generate-only timing, its note, and staging its changes to see them only at Generate.",
  "acceptedResidualRisk": "units may wait for the return question; MARGIN may hold back a unit under target (Axis).",
  "owner": "satoshikawato",
  "decisionDate": "2026-10-18"
}
```

### PD-OI-104: When the GC, skew and depth toggles apply

- Concern key: `web.layout.track-toggle-application`
- Scenario revision: `3`
- Supersedes: `PD-OI-104`, scenario revision `2`
  (`A / TOGGLE-LIVE-REDRAW-WHEN-WITHIN-TARGET`). Scenario revision `3` keeps
  every field of the scenario revision `2` receipt and adds to its May retire
  the packet's "Lost" text of that option (Revision 65 entry). Scenario
  revision `2` superseded scenario revision `1`
  (`B / TOGGLE-HIDE-AT-ONCE-NOTED`, the outcome of `S4-Q3`) and changed only
  what its receipt states; in the decision
  packet's words, "PD-OI-104 rev 1 holds for a toggle until its redraw is
  within the evidence rule (P3; reading and set as chosen) both ways".
- Status: `ACCEPTED`
- Selected outcomes: `B / TOGGLE-HIDE-AT-ONCE-NOTED` (scenario revision `1`)
  and `A / TOGGLE-LIVE-REDRAW-WHEN-WITHIN-TARGET`, reading `AS S5-Q1`, set
  `TOGETHER` (scenario revision `3`)
- Normative outcome: exactly the complete approved scenario revision `3`
  `PRODUCT_DECISION` receipt and its nine-field JSON representation below,
  with the effects of the scenario revision `1` receipt that it does not
  change. The scenario revision `1` receipt is kept as it was approved.
- Decision source (scenario revision `1`): the Owner replies quoted verbatim
  in the Revision 57 entry above. The receipt and JSON below reproduce all
  nine confirmed fields without translation or additional terms.
- Receipt SHA-256 (scenario revision `1`; UTF-8, excluding the final newline):
  `23aa2302bfdbf219d08593e7e37f484ab039efb7b030bc53a98d916228883aec`.
- Decision source (scenario revision `3`): the Owner replies quoted
  verbatim in the Revision 64 and Revision 65 entries above. The scenario
  revision `3` receipt and JSON below reproduce all nine fields without translation or
  additional terms. Dependent runtime requires this authority merged into its
  base; this record supplies no runtime acceptance evidence.
- Receipt SHA-256 (scenario revision `3`; UTF-8, excluding the final newline):
  `8d2096e0ea8236387807688b956388ce96e13f005ada0a8e08f130364137e3cd`.

Scenario revision `1` receipt, confirmed on `2026-10-17`:

```text
PRODUCT_DECISION
Concern: web.layout.track-toggle-application
Scenario revision: 1
Choice: B / TOGGLE-HIDE-AT-ONCE-NOTED
Rationale: hiding or showing a drawn GC, skew or depth track is instant today; redrawing on every toggle costs seconds on large Sessions. Keep the behaviour and state it truthfully.
Must preserve: all of today's behaviour (hiding or showing a drawn track applies at once in the shown Result; a track the Result lacks, and the layout, apply on Generate), with today's notes: "Hiding or showing a drawn track applies at once; a new track and the layout apply on Generate" (track visibility) and "Adding a Depth file shows and removing the last one hides the drawn Depth track at once; the rest applies on Generate" (Depth files). OV-1041 and OV-1043 are fixed to match this outcome (truthful notes or behaviour).
May retire: none.
Accepted residual risk: a gap where a hidden track was, until Generate.
Owner: Satoshi Kawato
Decision date: 2026-10-17
```

```json
{
  "concern": "web.layout.track-toggle-application",
  "scenarioRevision": 1,
  "choice": "B / TOGGLE-HIDE-AT-ONCE-NOTED",
  "rationale": "hiding or showing a drawn GC, skew or depth track is instant today; redrawing on every toggle costs seconds on large Sessions. Keep the behaviour and state it truthfully.",
  "mustPreserve": "all of today's behaviour (hiding or showing a drawn track applies at once in the shown Result; a track the Result lacks, and the layout, apply on Generate), with today's notes: \"Hiding or showing a drawn track applies at once; a new track and the layout apply on Generate\" (track visibility) and \"Adding a Depth file shows and removing the last one hides the drawn Depth track at once; the rest applies on Generate\" (Depth files). OV-1041 and OV-1043 are fixed to match this outcome (truthful notes or behaviour).",
  "mayRetire": "none.",
  "acceptedResidualRisk": "a gap where a hidden track was, until Generate.",
  "owner": "Satoshi Kawato",
  "decisionDate": "2026-10-17"
}
```

Scenario revision `3` receipt, confirmed on `2026-10-18`:

```text
PRODUCT_DECISION
Concern: web.layout.track-toggle-application
Scenario revision: 3
Choice: A / TOGGLE-LIVE-REDRAW-WHEN-WITHIN-TARGET; reading AS S5-Q1, set TOGETHER
Rationale: the G2 rule, and the instant response until the redraw is within it. TOGETHER: one PR retires the second writer and its exception whole; no state where one toggle redraws and its neighbour is patched.
Must preserve: R1 (c) parity; hiding stays live (Show Depth off after a removal too); OV-1043; Depth display settings (G3). Evidence rule: "A unit (a G2 unit, or the GT set or toggle as chosen in S5-Q2) is within the rule when its acceptance report shows, with the live mode of the step 4 spec (controls selected by CSS), the unit's representative change (U1: the section's slowest control; U2: each control group's representative change; GT: hide and show of each toggle), warm, at least 3 interleaved rounds with alternating Session order, each run started at a 1-minute load at or below 8 and a 5-minute load at or below 6, that in every round the median live apply (MARGIN: multiplied by 1.3; SENSITIVITY: itself, the x 1.3 result reported) is at or below its D-03 class target (P1 250, P2 400, P3 1000 ms) on each deciding Session. The deciding Sessions are the Gallery Vnig, hepatoplasmataceae orthogroup and vibrio collinear Sessions where the control is enabled and its change draws something different, with these cases fixed: hepa's record gap (one record) is left out; record names are judged on hepa and vibrio (Vnig's Multi-Record Canvas disables them); Dinucleotide controls on vibrio with Show GC Content turned on first; a toggle a Session has off is turned on first for its hide; Show Depth on all three Sessions, each with an added per-base Depth TSV (`samtools depth -a` text: one row per position of every record). The tobacco result, the stage D state and the stage work are reported. After SW1 is measured (and SW2 / SW3 where designed), every unit still over returns to the Owner as a 'stay on Generate (GT: stay instant) / accept the measured cost' question with its numbers."
May retire: `applyTrackVisibility`, its watcher, its History restore call (`app/app-setup.js` 3904), the R1 exception sentence; the instant hide and show (a redraw within the rule instead); Show Depth turned on after a Depth source change no longer shows the old drawn track at once.
Accepted residual risk: the gap after hiding until go-live (PD-OI-104's accepted risk); TOGETHER: GC waits.
Owner: satoshikawato
Decision date: 2026-10-18
```

```json
{
  "concern": "web.layout.track-toggle-application",
  "scenarioRevision": 3,
  "choice": "A / TOGGLE-LIVE-REDRAW-WHEN-WITHIN-TARGET; reading AS S5-Q1, set TOGETHER",
  "rationale": "the G2 rule, and the instant response until the redraw is within it. TOGETHER: one PR retires the second writer and its exception whole; no state where one toggle redraws and its neighbour is patched.",
  "mustPreserve": "R1 (c) parity; hiding stays live (Show Depth off after a removal too); OV-1043; Depth display settings (G3). Evidence rule: \"A unit (a G2 unit, or the GT set or toggle as chosen in S5-Q2) is within the rule when its acceptance report shows, with the live mode of the step 4 spec (controls selected by CSS), the unit's representative change (U1: the section's slowest control; U2: each control group's representative change; GT: hide and show of each toggle), warm, at least 3 interleaved rounds with alternating Session order, each run started at a 1-minute load at or below 8 and a 5-minute load at or below 6, that in every round the median live apply (MARGIN: multiplied by 1.3; SENSITIVITY: itself, the x 1.3 result reported) is at or below its D-03 class target (P1 250, P2 400, P3 1000 ms) on each deciding Session. The deciding Sessions are the Gallery Vnig, hepatoplasmataceae orthogroup and vibrio collinear Sessions where the control is enabled and its change draws something different, with these cases fixed: hepa's record gap (one record) is left out; record names are judged on hepa and vibrio (Vnig's Multi-Record Canvas disables them); Dinucleotide controls on vibrio with Show GC Content turned on first; a toggle a Session has off is turned on first for its hide; Show Depth on all three Sessions, each with an added per-base Depth TSV (`samtools depth -a` text: one row per position of every record). The tobacco result, the stage D state and the stage work are reported. After SW1 is measured (and SW2 / SW3 where designed), every unit still over returns to the Owner as a 'stay on Generate (GT: stay instant) / accept the measured cost' question with its numbers.\"",
  "mayRetire": "`applyTrackVisibility`, its watcher, its History restore call (`app/app-setup.js` 3904), the R1 exception sentence; the instant hide and show (a redraw within the rule instead); Show Depth turned on after a Depth source change no longer shows the old drawn track at once.",
  "acceptedResidualRisk": "the gap after hiding until go-live (PD-OI-104's accepted risk); TOGETHER: GC waits.",
  "owner": "satoshikawato",
  "decisionDate": "2026-10-18"
}
```

### PD-OI-105: A click on a drawn element shows its settings in a popup

- Concern key: `web.result.element-settings-entry`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `B / CLICK-DRAWN-ELEMENT-SHOWS-ITS-SECTION-IN-POPUP`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its JSON representation below.
- Decision source: the Owner replies quoted verbatim in the Revision 60 entry
  above. The Owner approved the receipt text below as written. The receipt
  reproduces its fields without translation or additional terms. Dependent
  runtime requires this authority merged into its base, or merged together
  with it through the reviewed static Product Contract co-change route; this
  record supplies no runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `0ed87f6520aeac57b438bed6bd00b609083e36fd6d03af4c28a94cb957caf82d`.

```text
PRODUCT_DECISION
Concern: web.result.element-settings-entry
Scenario revision: 1
Choice: B / CLICK-DRAWN-ELEMENT-SHOWS-ITS-SECTION-IN-POPUP
Rationale: 図の目盛りやタイトルを変えたいとき、設定ペインのどの section にあるかを探さなくても、図のその要素を押せば設定が出るようにする。設定は図の横に popup で出す。
Must preserve: 〔押せる要素と出る section〕Circular の目盛り → Axis & Scale。Linear の scale → Axis & Scale。GC content、GC skew などの dinucleotide の track → Dinucleotide content/skew。図のタイトル → Titles & Record Labels。〔押したとき〕Result の中に popup が出て、その section が入る。popup は押した要素に重ならない場所（右、左、下、上の順で、入るところ）に出る。どこにも入らないとき（Result いっぱいの輪の track、上下に余白のない横長の要素、スマートフォンの幅）は、押した場所の横に出て要素の一部を覆う。入るのは section そのもので、写しではない。popup に出ている間、その section は設定ペインにはなく、目次のボタンは残る。フォーカスは popup の見出しに移る。値は今までどおり Generate で反映される。〔popup が出ている間〕popup の中の編集は設定ペインでの編集と同じで、Undo の履歴に 1 つずつ入る。Undo、Redo、Generate、図の移動と拡大縮小、Fit、Result の外のクリックでは閉じない。見出しをドラッグすると popup が動き、見出しを押しても section は畳まれない。幅は 320 px で、長い section は popup の中でスクロールする。別の押せる要素を押すと、同じ popup の中身がその section に替わる。窓の大きさが変わって popup が画面からはみ出すときは、閉じずに、画面の中の一番近い位置へ動く。〔Result の popup は 1 つ〕feature、比較の帯、凡例の行、レコードの定義の popup が開いているときに押せる要素を押すと、その 1 回のクリックで前の popup が閉じて設定の popup が出る。それらの popup が開くと、設定の popup は閉じる。feature の選択は残る。〔閉じるとき〕popup の Close section のボタン。Escape（feature などの popup が開いていないとき。右の Editor の引き出しが開いていても、先に popup だけが閉じる）。Result の何もない場所、または popup も動作も持たない要素のクリック。ほかの popup が開いたとき。Layout edit を ON にしたとき。Review alignment を開いたとき。Circular と Linear を切り替えたとき（Undo で切り替わったときも）。Session を読み込んだとき（読み込みに失敗したときは、popup は読み込む前のまま残る）。Reset Settings。Run info のタブに切り替えたとき。Result がなくなったとき（最初の Generate より前まで Undo したとき）。目次のボタンまたはアイコンの列でその section を開いたとき。閉じると section は設定ペインに戻り、開いていたか畳んでいたかは popup に出す前のまま。フォーカスは popup を出す前にあった場所に戻る。〔アイコンの列〕押した要素の section が、横に開いている section と同じときは、popup は出ず、フォーカスがその見出しに移る。ほかの section のときは popup が出て、横に開いている section はそのまま。〔押せることの表示〕Layout edit が OFF のとき、これらの要素にマウスを載せると「Click to open its settings」が出る。目盛りと track ではポインターが手の形になる。タイトルと Linear の scale では、「Click to open its settings」と今の「Turn on Layout edit to move this item」が 2 行で出て、「?」のポインターと点線の枠は今のまま。〔出ないとき〕Layout edit が ON のあいだ。Review alignment を開いているあいだ（説明も手の形のポインターも出ない）。図をドラッグして動かしたとき。〔変えないもの〕feature、比較の帯、凡例の行、レコードの定義を押したときの動きと、それらの popup（窓の大きさが変わって画面からはみ出すときに、画面の中の一番近い位置へ動くようになる点だけが変わる）。popup を出す・閉じるだけでは、設定の値、図、Session、Undo の履歴、書き出した SVG・PNG・PDF は変わらない。popup とその位置は保存しない。
May retire: none
Accepted residual risk: 目盛りは線と小さい数字なので、押せる範囲が狭い。キーボードでは図の要素からは出せない（目次のボタンで同じ section を開ける）。Result いっぱいの輪の track や、上下に popup の入る余白がない横長の要素（図の下の端にある Linear の scale など）を押したとき、およびスマートフォンの幅では、popup が押した要素の一部を覆う。長い section の popup は Result の高さいっぱいになり、図の一部を覆う。popup はドラッグすると Result の外へも動かせる。popup に出ている間、設定ペインはその section のぶん短くなり、ペインのスクロールの位置がずれることがある。アイコンの列で、押した要素の section がすでに横に開いているときは、押しても見た目は変わらない。popup の見出しは、読み上げではボタンと読まれるが、popup の中では畳めない。Custom Track Slots が有効な図では、section の一部の欄が無効になっている。
Owner: satoshikawato
Decision date: 2026-10-13
```

```json
{
  "concern": "web.result.element-settings-entry",
  "scenarioRevision": 1,
  "choice": "B / CLICK-DRAWN-ELEMENT-SHOWS-ITS-SECTION-IN-POPUP",
  "rationale": "図の目盛りやタイトルを変えたいとき、設定ペインのどの section にあるかを探さなくても、図のその要素を押せば設定が出るようにする。設定は図の横に popup で出す。",
  "mustPreserve": "〔押せる要素と出る section〕Circular の目盛り → Axis & Scale。Linear の scale → Axis & Scale。GC content、GC skew などの dinucleotide の track → Dinucleotide content/skew。図のタイトル → Titles & Record Labels。〔押したとき〕Result の中に popup が出て、その section が入る。popup は押した要素に重ならない場所（右、左、下、上の順で、入るところ）に出る。どこにも入らないとき（Result いっぱいの輪の track、上下に余白のない横長の要素、スマートフォンの幅）は、押した場所の横に出て要素の一部を覆う。入るのは section そのもので、写しではない。popup に出ている間、その section は設定ペインにはなく、目次のボタンは残る。フォーカスは popup の見出しに移る。値は今までどおり Generate で反映される。〔popup が出ている間〕popup の中の編集は設定ペインでの編集と同じで、Undo の履歴に 1 つずつ入る。Undo、Redo、Generate、図の移動と拡大縮小、Fit、Result の外のクリックでは閉じない。見出しをドラッグすると popup が動き、見出しを押しても section は畳まれない。幅は 320 px で、長い section は popup の中でスクロールする。別の押せる要素を押すと、同じ popup の中身がその section に替わる。窓の大きさが変わって popup が画面からはみ出すときは、閉じずに、画面の中の一番近い位置へ動く。〔Result の popup は 1 つ〕feature、比較の帯、凡例の行、レコードの定義の popup が開いているときに押せる要素を押すと、その 1 回のクリックで前の popup が閉じて設定の popup が出る。それらの popup が開くと、設定の popup は閉じる。feature の選択は残る。〔閉じるとき〕popup の Close section のボタン。Escape（feature などの popup が開いていないとき。右の Editor の引き出しが開いていても、先に popup だけが閉じる）。Result の何もない場所、または popup も動作も持たない要素のクリック。ほかの popup が開いたとき。Layout edit を ON にしたとき。Review alignment を開いたとき。Circular と Linear を切り替えたとき（Undo で切り替わったときも）。Session を読み込んだとき（読み込みに失敗したときは、popup は読み込む前のまま残る）。Reset Settings。Run info のタブに切り替えたとき。Result がなくなったとき（最初の Generate より前まで Undo したとき）。目次のボタンまたはアイコンの列でその section を開いたとき。閉じると section は設定ペインに戻り、開いていたか畳んでいたかは popup に出す前のまま。フォーカスは popup を出す前にあった場所に戻る。〔アイコンの列〕押した要素の section が、横に開いている section と同じときは、popup は出ず、フォーカスがその見出しに移る。ほかの section のときは popup が出て、横に開いている section はそのまま。〔押せることの表示〕Layout edit が OFF のとき、これらの要素にマウスを載せると「Click to open its settings」が出る。目盛りと track ではポインターが手の形になる。タイトルと Linear の scale では、「Click to open its settings」と今の「Turn on Layout edit to move this item」が 2 行で出て、「?」のポインターと点線の枠は今のまま。〔出ないとき〕Layout edit が ON のあいだ。Review alignment を開いているあいだ（説明も手の形のポインターも出ない）。図をドラッグして動かしたとき。〔変えないもの〕feature、比較の帯、凡例の行、レコードの定義を押したときの動きと、それらの popup（窓の大きさが変わって画面からはみ出すときに、画面の中の一番近い位置へ動くようになる点だけが変わる）。popup を出す・閉じるだけでは、設定の値、図、Session、Undo の履歴、書き出した SVG・PNG・PDF は変わらない。popup とその位置は保存しない。",
  "mayRetire": "none",
  "acceptedResidualRisk": "目盛りは線と小さい数字なので、押せる範囲が狭い。キーボードでは図の要素からは出せない（目次のボタンで同じ section を開ける）。Result いっぱいの輪の track や、上下に popup の入る余白がない横長の要素（図の下の端にある Linear の scale など）を押したとき、およびスマートフォンの幅では、popup が押した要素の一部を覆う。長い section の popup は Result の高さいっぱいになり、図の一部を覆う。popup はドラッグすると Result の外へも動かせる。popup に出ている間、設定ペインはその section のぶん短くなり、ペインのスクロールの位置がずれることがある。アイコンの列で、押した要素の section がすでに横に開いているときは、押しても見た目は変わらない。popup の見出しは、読み上げではボタンと読まれるが、popup の中では畳めない。Custom Track Slots が有効な図では、section の一部の欄が無効になっている。",
  "owner": "satoshikawato",
  "decisionDate": "2026-10-13"
}
```

## Acceptance contract catalog

| Contract | Required meaning |
| --- | --- |
| `OIC-001` | Candidate limit is truthful; `None` remains unbounded and no hidden cap is applied. |
| `OIC-002` | Candidate, Pairwise display, and member-hit limits are independent and invalidate only the correct stages. |
| `OIC-003` | Every supported Collinear enum reaches the real typed Python analysis path. |
| `OIC-004` | Fresh defaults follow the declared surface and mode rules; explicit imported values win. |
| `OIC-005` | Canonical request, resolved values, actual helper invocation, stage cache identities, and artifact provenance agree. |
| `OIC-006` | Required controls are discoverable, operable, mode-safe, persistent, and accessible. |
| `OIC-007` | Imported comparison intent is never silently cleared; unresolved state has explicit actions. |
| `OIC-008` | Valid GUI-unmanaged config survives and is disclosed; invalid or unknown input is rejected. |
| `OIC-009` | Circular selection, crop, and reverse use one effective path and round-trip without double reversal. |
| `OIC-010` | Circular label/subtitle are effective; empty values preserve inferred output. |
| `OIC-011` | Scale and ruler-label font sizes remain independently public with explicit linked-default behavior. |
| `OIC-012` | `grid_column`, deferred direct-link topology, and GUI surface scope are represented accurately. |
| `OIC-013` | Failed, canceled, or stale Generate preserves the committed request and last successful Result. |
| `OIC-014` | Product authority, Product Impact mapping, Architecture Ratchet, runtime owners, and evidence remain separate. |
| `OIC-015` | Linear discovery, per-record placement, actual comparison jobs, explicit request endpoints, Session replay, and SVG endpoints retain the complete selected record universe. Adjacent uses neighboring-row Cartesian products; Similarity and Collinear `all` retain complete directed between-record evidence regardless of display rows; Collinear within-record evidence follows the inference checkbox. |
| `OIC-016` | LOSAT Auto allocation follows Total threads; displayed and executed budgets agree, and manual intent survives temporary clamps and Session replay. |
| `OIC-017` | Web raw/member defaults are 5/5 in Collinear and unbounded/unbounded in Similarity. Each mode restores its own edits repeatedly, including blanks; Session round trips retain both modes; Reset Settings restores defaults. |
| `OIC-018` | Collinear inference defaults OFF; actual raw jobs exclude every self-comparison, including within multi-record source batches, and the real Python path skips orthogroup inference. ON retains the existing inference; request, cache, provenance, and legacy Session interpretation agree. |
| `OIC-019` | Completed raw searches survive downstream cancellation for matching retries; member-only edits do not rerun LOSAT. Raw-setting/input changes, Clear Cache, and Session/History replacement prevent incompatible reuse; the committed Result remains intact. |
| `OIC-020` | Linear File cards expose common Depth TSV assignment without expanding records. File-level apply and clear update only that File and logical series as one undoable operation; empty, common, and mixed states remain truthful. Per-record sparse overrides, logical indexes, canonical requests, Session replay, and regeneration remain unchanged. |
| `OIC-021` | Feature-popup record rotation uses the explicit popup target and source coordinates, changes only one effectively circular record through a target-only atomic transaction, preserves pending edits and the prior artifact on every no-op path, round-trips the absolute transform, and reuses compatible LOSAT evidence without additional executor jobs. |
| `OIC-022` | `PD-OI-036`: independent fresh/reset Auto, explicit Show/Hide, diagram-wide resolution excluding dormant rows, saved Session/Result and Undo/Redo remain supported; Layout and Labels disclose the effective Auto result, affected fields, reason, scope, next Generate effect, and route to Record Labels. |
| `OIC-023` | `PD-OI-024`: only Web fresh/reset selects Lock ON; explicit OFF, saved values and supported omission meanings, saved Result on Load, CLI/Python defaults, accepted anchors, D2-P and D3-A remain supported. ON/OFF and Generate application are always explained in Linear Layout. |
| `OIC-024` | `PD-OI-037` revision 2: operation-level Live edit, Applies on Generate, and Apply required remain truthful, including Palette Instant Preview and Alignment review. No always-on derived Pending/Applied status or its computation remains. Live applying/error, Processing/Canceling, real errors, and recovery are shown when they occur, with no status-only Worker, genome-byte read/hash, canonical request projection, generated-table build, or SVG/checkpoint clone. |
| `OIC-025` | `PD-OI-038`: one SVG and Editor retain live commit/rerender, all tabs, availability, History/Session/Export, camera, keyboard, visibility-only Close/Escape, selected tab and Result recovery. At 390×844/740 the canvas uses the available full width and at least 200 px height; content scroll, reachable header/Close/toolbar, short-viewport/soft-keyboard access, and wide side drawer remain required. Pointer/keyboard/browser verification is required; duplicated Preview/SVG/editor is not accepted. |
| `OIC-026` | `PD-OI-035` and `PD-OI-039`: identity, keyboard Select/Skip, non-rendered candidates, no position-only selection, desktop canvas, focus and transient overlay exclusion remain required with all PD-OI-031/034 outcomes. Compact review retains visible, operable canvas at full available width and at least 200 px height at 390×844/740, scrollable candidates and reachable Apply/Cancel, local no-Worker draft edits, atomic batch validation, failure/error/retry and artifact/orientation/History recovery. At every width the review is docked in the Result, closes Editor through its owner while retaining tab, disables reopening with a reason until review ends, then reopens an Editor that was open on the same tab; the review cannot be dragged. During the review the feature popup is hidden, a click on the figure neither opens, switches nor closes it, search Open switches its target, and it returns when the review ends. The figure is fitted when the review opens and the earlier zoom and pan return when it ends; a Session saved during the review and the History of its Apply keep the earlier zoom and pan. While the review shows the figure it fitted, the figure is fitted again, as by the Fit button and without a History step, whenever the width or height of its room changes; a zoom and pan changed by hand stay while the room does not change, and outside a review a change of the room keeps the zoom and pan. Browser verification must show presentation changes leave draft and Result unchanged. |
| `OIC-027` | `PD-OI-066` scenario revision `2`: a live-edited Result, including after Session load and export, agrees with a Result freshly generated from the same draft in the meaning (position, color, text, and visibility) of every edited element when no pending Applies on Generate value changes those elements; otherwise it agrees with a Generate of the committed Session plus the applied live values. Settings shown as Applies on Generate do not change the Result before Generate. A live edit that the Generate compiler cannot reproduce is shown as Applies on Generate instead. |
| `OIC-028` | `PD-OI-088`: the color scope, Legend name, and Reset fill color dialogs, whether opened from the feature popup's fill color, the Features drawer color input, a popup stroke edit, or a Legend panel rename (`Merge into "<name>"` or `Keep both: name it "<name 2>"`), open from the saved rules before the rules a choice may add are prepared. A choice is one History step and Cancel records none. From a choice until its step ends, the dialog stays open, its choices and Cancel are disabled, it states that the edit is applying, and no second choice, Cancel, Escape, or backdrop click starts or closes anything; the dialog closes when the step ends. |
| `OIC-029` | `PD-OI-089`: a user default color is a Default colors value that differs from the selected palette's color (case, named colors, and `#rgb` compared as one color; Auto is none). A palette switch with user colors opens one dialog: Keep my N colors (the new palette plus the user colors), Use the palette's colors, or Cancel (no change, no History step, the select shows the old palette); without user colors it switches as one step. Switching back to the applied palette while another is queued asks the same when user colors exist; the queued palette is dropped and the kept colors apply now. Default colors Reset with user colors asks first. Dialog bodies name up to three keys. `Legend item "<legend name>" (N)` on a Legend row that no Specific color rule draws sets that feature type's default color as one History step and writes no Specific color rule; the scope dialog shows under it "Sets the <type> default color (<type> features without their own color or rule, also hidden ones)." or, when the type has a user color, "Replaces your <type> default color ■ #hex." While a palette is queued, that choice sets the color in the queued palette's colors and in the shown Result now. Generate and the Source recipe write `-p` and only the user colors as `-d`. |
| `OIC-030` | `PD-OI-090`: an edit, Generate, failure, cancel, or Undo/Redo of a step in one drawing leaves every other drawing's settings, Legend edits, History, and Result unchanged, and Reset all drawings is the only action that resets more than the shown drawing; Undo/Redo in any drawing walks the shown drawing's steps and the project steps (drawing add, duplicate, rename, reorder, delete, Reset all drawings, and edits that change only the Session title, Auto Reflow, PNG DPI, palette Instant Preview, LOSAT execution and threads, or input drafts) newest first, one step per operation; an operation that also changes a drawing value is one step of that drawing; turning Instant Preview on applies every drawing's queued palette in that one step; a drawing switch writes no setting and is not an Undo step; a drawing switch is refused while a color rule table edit is being matched, and that edit becomes one History step of its drawing; Save and a fresh Load restore every saved drawing, its Result, and the drawing that was shown; each Session 27-44 fixture and each Session 46, 47, and 48 fixture converts with no saved value lost and its first drawing renders as the old Session replays; a Session 46, 47, 48, or 49 preview-only Load starts no Python Worker, except that a Session 46 or 47 holding protein comparison results stored on a reversed record shows its stored figure first and then starts it once to convert them; the Python and CLI replay and render one or all drawings with distinct output paths; a region drawing keeps the source coordinates and carries over exactly what `PD-OI-090` lists. |
| `OIC-031` | `PD-OI-099`: a Load of a file with something it cannot fully restore lists it in one dialog before any setting, Result, input file, History, or screen changes, and Cancel changes none of them. Every Gallery Session, and every Session fixture that loses nothing, loads with no row; a fixture that is lossy on purpose lists exactly its loss. A refusal of the Web Load comes before the dialog, never after [Load]. No file the Web loaded before the dialog existed is refused. The receipt's 〔拒否〕 classes are refused; once the interim of its Accepted residual risk ends, they are the only refusals of the Web Load. |

These new acceptance entries are obligations for dependent runtime work, not
claims of completed runtime or browser verification by this authority amendment.
The signed receipts remain the complete outcome; existing acceptance contracts
and independent preserved guarantees remain jointly required.

### OIC-020 required regression coverage

The normal automated PR gate must observe all of the following:

- A multi-record Linear GenBank File exposes its Depth TSV assignment while its
  record list and record options remain closed. Applying one file binds the same
  logical series to every record in that File and does not affect another File.
- One per-record replacement produces a truthful mixed File state. Applying or
  clearing the File-level value then replaces or clears every record cell in
  that File and series, and one Undo restores the complete prior matrix.
- Empty cells and later logical columns do not shift when a source is cleared.
  Same-named independent files remain distinct, including after Session
  restoration.
- Save, fresh Load, canonical request construction, generation, and subsequent
  regeneration preserve common and mixed bindings without a new Session schema,
  request schema, Worker protocol, or rendering path.

### OIC-021 required regression coverage

| ID | Required observation |
| --- | --- |
| `AC-01` | The popup alone rotates the target feature's record and never consults another or global selection. |
| `AC-02` | An effectively circular record resolves the same source request in Circular and Linear diagram modes. |
| `AC-03` | Every non-target record and layout value remains unchanged, including same-file multi-record inputs. |
| `AC-04` | Positive and negative offsets resolve relative to feature direction and wrap correctly. |
| `AC-05` | Orientation intent is absolute and idempotent; leaving it off preserves the current value. |
| `AC-06` | The 3-prime base anchor and feature-end placement remain distinct operations. |
| `AC-07` | Multipart, origin-spanning, and odd/even covered midpoints follow exact covered traversal. |
| `AC-08` | Unstranded, ambiguous, fuzzy, cropped, linear-topology, and stale cases expose operation-specific reasons. |
| `AC-09` | Duplicate record IDs and split fragments retain stable source-bound identity. |
| `AC-10` | One Undo or Redo restores or reapplies origin, orientation, provenance, and Result together. |
| `AC-11` | Save and fresh Load restore the absolute transform and provenance. |
| `AC-12` | Failed, canceled, stale, and superseded rendering preserves the previous Result and transform. |
| `AC-13` | Unrelated pending edits remain pending and are neither applied nor discarded. |
| `AC-14` | A transform-only operation adds zero LOSAT executor jobs and reuses compatible raw evidence. |
| `AC-15` | Feature, label, tick, depth, statistics, and comparison geometry use the same transform. |
| `AC-16` | Manual display-start or orientation changes clear stale anchor provenance. |
| `AC-17` | Cancel changes no draft, Result, transform, provenance, or History state. |
| `AC-18` | Search query and stable target re-identification survive replacement; popup actions remain keyboard- and 390 px-accessible. |
| `AC-19` | Request schema 7, the Worker protocol, and the renderer path do not expand. |
| `AC-20` | Product Impact and Architecture Ratchet evidence remain reviewable and all required gates pass. |

### OIC-015 required regression coverage

The normal automated PR gate must observe all of the following:

- Multi-record GenBank and GFF3/FASTA discovery retains every record with
  comparisons enabled and disabled; explicit selectors remain exact.
- A two-record source and a three-record source expose two file-input cards
  and five independent placements. A single two-record upload has one file
  card and two record controls. Assert both counts after Save/fresh Load,
  source replacement, source removal, and per-record row changes; filenames
  alone must not collapse distinct input sources. The neighboring 2-by-3 rows
  retain all six Pairwise record combinations through one source-file search
  job. The number of record pairs must not be reported as LOSAT job count.
  The record list initially shows only its single-line count; pointer and
  keyboard expansion expose every record without changing its state.
  Fresh and reset documents use same-file Linear rows and Circular shared
  canvas; explicit saved opt-outs survive Load.
- Similarity and Collinear `all` with inference ON retain the complete directed
  record-pair matrix through source-file jobs, including self, same-row, and
  non-adjacent evidence. Two sources containing five records require four
  source jobs and cover 25 directed record pairs. Collinear `all` with inference
  OFF retains all 20 between-record directions and searches no record against
  itself, including within multi-record sources. A single-row layout retains
  the selected analysis even though it has no between-row links.
- Explicit endpoints survive the typed decoder and reach the intended SVG
  records, including a numerically consecutive pair beside a multi-record row.
- Save, fresh Load, regeneration, and reordering preserve record identity,
  placement, pair mapping, and shared source resources. Repeated execution
  reuses only semantically equivalent cached searches.
- Every LOSAT mode exposes the existing execution and thread controls in
  Comparison Settings; mode switches and Session replay preserve their values.
  Verify Run LOSAT occupies the full top row, with No comparison and Upload
  BLAST TSV side by side below it. At desktop and narrow viewport widths,
  DOM order and actual keyboard Tab traversal match that visual order.
  Pointer and keyboard activation of Run LOSAT open Settings without moving
  focus away from the command; reopening it requires no mode change. Fresh
  and reset Execution is `threaded`, and saved explicit execution modes survive
  Load with active LOSAT settings immediately exposed.
- Changing only a record's drawing start or reverse-complement display issues
  no additional LOSATP search jobs, updates the rendered coordinates/orientation,
  and preserves that reuse through Save and fresh Load.
- Unbounded member selection retains more than five threshold-qualified hits.
  Explicit finite member limits still apply, and changing only that limit
  invalidates derived output while preserving raw search cache reuse.

These are observations of jobs, controls, requests, Sessions, and SVG results;
checking an `All` label or counting uploaded files is insufficient. Restoring
first-record truncation, zipped Adjacent pairing, or positional endpoint
coercion must fail the corresponding regression test. Source changes to these
boundaries require this coverage in the normal PR gate.

### OIC-031 required regression coverage

These tests protect the outcome; they do not supply Product authority.

| Observation | Executable regression owners |
| --- | --- |
| The dialog lists what Load cannot restore before anything changes; Cancel changes nothing. | [`load-repair-dialog.test.mjs`](../../tests/web/load-repair-dialog.test.mjs), `a lossy Load asks once before it changes anything, and Cancel leaves everything as it was (PLAN 12.7 (1), (2))`; [`session-loading-feedback.playwright.spec.js`](../../tests/web/session-loading-feedback.playwright.spec.js), `a lossy Session asks before it loads: Cancel changes nothing, Load loads it`. |
| Every Gallery Session, and every Session fixture that loses nothing, loads with no row; a fixture that is lossy on purpose lists exactly its loss. | [`load-repair-dialog.test.mjs`](../../tests/web/load-repair-dialog.test.mjs), `every Gallery Session and Session fixture loads without a repair (G3)`. |
| A refusal comes before the dialog, never after [Load]. | [`load-repair-dialog.test.mjs`](../../tests/web/load-repair-dialog.test.mjs), `a configuration field Load refuses is refused before the dialog (BA-OD96)`. |
| No file the Web loaded before the dialog existed is refused. | [`load-repair-dialog.test.mjs`](../../tests/web/load-repair-dialog.test.mjs), `no file the Web loaded before the dialog existed is refused (BA-OD100)`. |
| The 〔拒否〕 classes are refused. | [`session-import-codec.test.mjs`](../../tests/web/session-import-codec.test.mjs), `malformed JSON/gzip and fatal UTF-8 stay structured failures`, `unsafe own keys remain untrusted after the production whole-object reply`, and `codec rejects exact file and expanded caps before reading/decoding oversize data`; [`session-import-client.test.mjs`](../../tests/web/session-import-client.test.mjs), `Worker unavailable, constructor and postMessage failures have no main parser fallback`; [`session-authority.test.mjs`](../../tests/web/session-authority.test.mjs), the `readSessionInput` refusals; [`mode-scoped-session.test.mjs`](../../tests/web/mode-scoped-session.test.mjs), `a file that is not a gbdraw Session is refused as such before its contents are read`; [`session-losat-cache-validation.test.mjs`](../../tests/web/session-losat-cache-validation.test.mjs), `session import rejects branch-internal version 34` to `38`; [`load-repair-dialog.test.mjs`](../../tests/web/load-repair-dialog.test.mjs), `a Session older than 27 is refused as an old version (B51)`; [`session-import-worker.playwright.spec.js`](../../tests/web/session-import-worker.playwright.spec.js), `unsafe keys, malformed JSON/gzip, fatal UTF-8 and Worker crash preserve the old document and History`. Request formats 3, 4 and newer: no test yet. |

## Residual-risk boundary

Accepted residual risk never authorizes a security vulnerability, silent
scientific-output corruption, loss of a must-preserve effect, deterministic
Architecture Ratchet failure, undocumented unbounded performance regression,
cache reuse across different execution semantics, artifact provenance that
disagrees with actual execution, or failure of a required acceptance contract.


## OIC-016 acceptance evidence

- With four source jobs, change Total threads between 32, 16, and 2 through the
  visible control and observe both Auto labels recalculate as specified above.
- Exercise one explicit concurrency value with automatic per-run threads and
  one explicit per-run thread value with automatic concurrency. Verify the
  effective budget bound after decreasing and increasing Total threads.
- Preserve explicit choices and Auto through Save and fresh Load, including a
  temporarily clamped value. Show the effective clamp without rewriting intent.
- Observe real threaded LOSATP dispatch and completion on a small workload;
  compare the worker allocation and runtime report with the displayed plan.
- Keep fixed one-thread-per-run programs fixed. Changing only scheduling
  settings preserves the existing raw-search cache identity.

## OIC-017–OIC-019 regression evidence owners

These tests protect the outcomes; they do not supply Product authority.
The mode-restoration, inference-toggle, and cancellation/retry browser
observations below must run in the normal automated PR gate. Full dev staging
provides additional coverage and does not replace this pre-merge requirement.

| Contract | Executable regression owners |
| --- | --- |
| `OIC-017` | [`comparison-ui.playwright.spec.js`](../../tests/web/comparison-ui.playwright.spec.js), `comparison controls drive appearance and current Session round trips`: initial values, repeated mode restoration, saved inactive values, blanks, and Reset Settings. |
| `OIC-018` | [`linear-multi-record.playwright.spec.js`](../../tests/web/linear-multi-record.playwright.spec.js), `Collinear inference checkbox skips self searches and reuses matching evidence`; [`linear-sources.test.mjs`](../../tests/web/linear-sources.test.mjs); [`losat-settings.test.mjs`](../../tests/web/losat-settings.test.mjs); [`test_collinearity.py`](../../tests/test_collinearity.py); [`test_session_request_codec.py`](../../tests/test_session_request_codec.py); [`test_api_request_render.py`](../../tests/test_api_request_render.py). Protect actual job endpoints, no inference call when OFF, member selection, both scopes, all anchor modes, and request/derived identities. |
| `OIC-019` | [`linear-multi-record.playwright.spec.js`](../../tests/web/linear-multi-record.playwright.spec.js), `protein raw cache survives cancellation and derived options preserve search identity`; [`run-analysis-simple-path.test.mjs`](../../tests/web/run-analysis-simple-path.test.mjs). Protect cancellation, repeated retry, member edits, raw-setting/input changes, clearing, and late artifact failure. |
