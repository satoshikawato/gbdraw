# Option Integrity Product Contract

Status: active Product authority

## Authority metadata

- Contract ID: `OIPC`
- Contract revision: `35`
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
  `PD-OI-086`, `PD-OI-087`, and `PD-OI-088`
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
- Scenario revision: `2`
- Supersedes: `PD-OI-028`, scenario revision `1` (`A / PERSISTED_ACTIVE_ALIGNMENT_PLAN`).
- Status: `ACCEPTED`
- Selected outcome: `A / ORIENTATION_INDEPENDENT_ACTIVE_PLAN`
- Normative outcome: exactly the approved `PRODUCT_DECISION` receipt below.
- Decision source: The complete five `PRODUCT_DECISION` texts for the
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
  interaction are governed by `PD-OI-039`, scenario revision `2`, in full.
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
- Scenario revision: `2`
- Status: `ACCEPTED`
- Selected outcome: `A / EXCLUSIVE_DIRECTIONS_WITHOUT_MATCH`
- Supersedes: `PD-OI-039`, scenario revision `1`
  (`A / DOCKED-COMPACT-ALIGNMENT-REVIEW`), replacing its Match retention
  with the explicitly approved exclusive direction controls. All independent
  identity, keyboard/Skip, non-rendered-candidate, canvas, focus, Editor,
  validation/retry, artifact, Session and History requirements remain.
  The scenario-2 mobile coverage exception and close-review continuation in
  `PD-OI-035` remain superseded; its scenario-3 independent guarantees remain
  jointly required. Narrow free drag and concurrent Editor opening during
  review may retire only as stated in this receipt. Wide drag and non-modal
  canvas interaction remain required.
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  below. The four issue `#598` decisions are unchanged.
- Decision source: `satoshikawato` selected A and explicitly approved all nine
  fields of the displayed receipt on `2026-09-26`, then authorized this
  authority-only update through dev. Receipt UTF-8 SHA-256:
  `06a9d2fe9b1d1406f6f8e04c23a9ca031133b9fae24b0683403ec2c6cae55270`.
  This serialization adds no retirement or risk terms to that approval.
  Dependent runtime requires the supersession merged into its base.
- Acceptance contracts: `OIC-006`, `OIC-013`, `OIC-014`, `OIC-026`.

```json
{
  "concern": "web.similarity-alignment.review-presentation",
  "scenarioRevision": 2,
  "choice": "A / EXCLUSIVE_DIRECTIONS_WITHOUT_MATCH",
  "rationale": "狭いPreviewでもalignment候補をcanvasで確認できるよう、reviewを図の下段に固定し、候補比較へ操作を集中させる。表示方向はIssue #598のKeep/right/left/Customへ統一し、reference相対のMatch操作による結果との混同を避ける。",
  "mustPreserve": "PD-OI-031/034と現行transform/plan/reset/historyのすべての結果。resolvedの通常自動Apply、ambiguousと明示reviewのlocal draft、独立Select/Skip、候補根拠とreference identity、Issue #598で承認済みの排他的Keep/right/left/Custom、referenceを含むselected known-strand anchorsの方向選択と各recordのbefore/after矢印、unknown/skipped/missing/unusableの理由付き不変、canvas操作、local編集でWorkerを呼ばないこと、Applyの共有Python batch validationとatomic Result/History。失敗時draft/error/retry、Cancel/stale/superseded時の以前のResult/orientation/History、Session/regeneration/Export、focus復帰を維持する。390×844/740では利用可能幅全体かつ高さ200px以上のcanvasを確保し、候補listをscroll、Apply/Cancelを到達可能にする。狭いreview開始時はEditorをownerで閉じ、tabを保持し、review中は理由付きでopenをdisable、終了後は明示reopen可能。wideのdragと非モーダルcanvasを維持する。",
  "mayRetire": "旧Match reference direction checkbox・flag・操作affordance。狭いPreviewでreviewを自由にdragする操作、およびreview中にEditorを同時openする継続。これ以外のPD-OI-035/039の独立要求とIssue #598の承認済み4決定は退役しない。",
  "acceptedResidualRisk": "狭いreviewではlist scrollが増え、自由に位置を動かせなくなる。開始時Editorは閉じるがtabは保持し、終了後再openできる。位置変更で候補draftやResultを変えないことをbrowserで確認する。旧Match利用者はright/leftまたはCustomで表示方向を明示的に選ぶ必要がある。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-26"
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
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `A / CAPTION_AUTO_DISAMBIGUATE_SOLID_ROWS`
- Normative outcome: 同 caption・異色 rule を受け付け、全色に lowercase normalized hex を付けた caption を canonical rule として採用する（例 Transporter [#112233] / Transporter [#445566]）。同名同色は共有、空 caption は凡例なし。既存 literal caption/legend key に衝突する場合は予約後に決定的な追加 suffix で区別。first/last-wins は廃止。各実際に使用された色を別の solid 凡例行で示す。
- User access / feedback: import/manual edit の正常完了時に caption 変更を通知。利用者は生成した solid 行を既存 editor から編集できる。
- Session / regeneration: admitted caption は普通の文字列として保存。Load は preview/draft を保ち自動 Generate しない。過去の同名異色 draft は次の rule edit/Generate の通常 preparation で通知付き正規化。
- Export / artifact: 新 TSV は区別した canonical caption。SVG/PNG/PDF は各色の solid 行。元ファイルの同名 caption のままの lossless 復元は約束しない。
- Failure / recovery: 自動区別後すぐ図を使える。必要なら caption を編集して再生成。stale/preparation failure は直前 rules/Result に戻す。
- Approval receipt: 「すべて推奨案で承認します。」, selecting exactly this
  complete Choice A outcome, approved by `satoshikawato` on `2026-09-26`.
- Decision source: The complete outcome and `PRODUCT_DECISION` receipt in
  [`APPROVED_PRODUCT_DECISIONS.md`](./issue-600-implementation-20260926/APPROVED_PRODUCT_DECISIONS.md),
  reproduced without inferred rationale, preservation, retirement, or risk
  terms. This is an inert serialization in the existing static authority;
  it cannot authorize dependent runtime until merged into its base.
- Acceptance IDs: `CLR-01`, `CLR-02`, `CLR-03`, `CLR-04`, `CLR-05`; their definitions remain owned by
  [`MASTER_PLAN.md`](./issue-600-implementation-20260926/MASTER_PLAN.md).

```json
{
  "concern": "styles.specific-color-caption-multiplicity",
  "scenarioRevision": 1,
  "choice": "A / CAPTION_AUTO_DISAMBIGUATE_SOLID_ROWS",
  "rationale": "近い色を使う rule を拒否せず、実際の各色を凡例に表示する。今回は既存 solid 行を再利用する自動 caption 区別を採用し、複数 swatch 用の renderer・editor・保存形式を追加せずに fresh/live/native の意味を揃える。",
  "mustPreserve": "feature 色、rule 順序/precedence、同名同色共有、single-color caption、unused rule の凡例除外、既存 solid editor・保存 preview・failure/History 契約。 rule order/regex/precedence/visibility を変えない。使用色の忠実な図示、stable identity、file/manual provenance、Result rollback、既存 SVG sanitizer を維持。",
  "mayRetire": "caption 衝突による Web 拒否、過去の last-wins/上書き。退役は specific-color rule の同名異色 scope に限定。",
  "acceptedResidualRisk": "凡例が長くなり元の同名文字列は変わる。multi-swatch grouping を望む利用者には複数行になる。hex suffix、既存 edit、layout 再計測で扱いを明確にする。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-26"
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
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `A / REVEAL_APPLICABLE_SINGLE_RECORD_CONTROLS`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt below.
- Approval receipt: 「推奨案で承認します」, selecting this complete outcome,
  approved by `satoshikawato` on `2026-09-26`.
- Reviewed outcome SHA-256: `57debd7625f007c99849b8de5f992508b369b7a2277c01ed1ff51e7c95354d93`.
- Decision source: the complete receipt in
  `docs/internal/issue-597-input-session-implementation-20260926/decisions/02_RECORD_DISCOVERY.md`
  at commit `51a786086dc2777e5fa375e5f94d8b7ac7deeedc` of
  `fix/issue-597-input-session-20260926`. Its nine supplied fields are
  reproduced below without translation or additional terms. The receipt
  document is evidence of the human choice; this record is its static authority.
  Dependent runtime requires this authority merged into its base.

```json
{
  "concern": "diagram-generation.circular-transform-discoverability",
  "scenarioRevision": 1,
  "choice": "A / REVEAL_APPLICABLE_SINGLE_RECORD_CONTROLS",
  "rationale": "source の探索状況と一件用 crop の適用条件を区別して示し、編集可能になった一件用 controls は selection の直後に見えるようにする。通常 upload の自動探索と保存済みプレビューの軽い閲覧を両立する。",
  "mustPreserve": "valid native upload の自動 record discovery と Generate 前の rotation controls、既存 parser/helper 境界、exact source-bound identity、explicit single/grid/batch、fresh shared-canvas default、saved explicit choices、一件 crop と topology/start/reverse の適用条件、手動 close/expand、元の focus、keyboard/390 px、preview-only Load の Python Worker 0、active draft と saved artifact の分離、失敗時の旧 Result、Retry/Replace/Remove/Inspect/Generate の継続。",
  "mayRetire": "適用可能になった一件用 section が常に collapsed で始まる挙動、valid fresh upload に manual Load が必須であるかのような prompt、実行中でない deferred discovery を loading と表す UI。全 record subset editing や複数 source support の選択は含まない。",
  "acceptedResidualRisk": "applicable になった時に一件用 section が展開されて pane 高さが変わる。操作元の focus と scroll anchor を維持し、無関係な更新で再展開しない。科学的意味の変更、先頭 record の自動選択、grouping の自動切替、preview-only Load による Python 初期化は受容しない。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-26"
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
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `A / EXPLICIT-MODE-WITH-DISCOVERABLE-TARGETS`
- Normative outcome: OFF の drag は従来の canvas pan。supported target は help cursor/hover枠と「Turn on Layout edit to move this item」を表示。toolbar の常設説明、keyboard focus、touch でも同じ情報へ到達。ON は grab、drag 中は grabbing。target 上で mode を自動ONにしない
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
- Decision source: the complete Choice A outcome and nine-field receipt in
  [DECISION_02_LAYOUT_AFFORDANCE.md](https://github.com/satoshikawato/gbdraw/blob/bcc4e0aa5ffcf4bdf9a12952caa8bf5de2cfee08/docs/internal/issue-599-preview-layout-20260926/DECISION_02_LAYOUT_AFFORDANCE.md)
  at published planning commit `bcc4e0aa5ffcf4bdf9a12952caa8bf5de2cfee08`.
  `satoshikawato` explicitly approved all three independent recommended A
  outcomes with 「すべて推奨案で承認します。」 on `2026-09-26`.
  The receipt and JSON below reproduce this concern's supplied fields without
  translation or additional rationale, retirement, or risk terms. This record
  does not supersede another decision. Dependent runtime requires this
  authority merged into its base; this amendment contains no runtime or
  runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `0365e0057b61e6b606c6714dc17f0e5c89d755ce651c1a8eec66d0139e0c65ec`.
- Acceptance contracts: `OIC-006`, `OIC-013`, `OIC-014`. These existing obligations and the
  complete selected outcome are jointly required; their citation does not
  claim completed dependent-runtime checks.

```text
PRODUCT_DECISION
Concern: web.layout-edit-affordance
Scenario revision: 1
Choice: A / EXPLICIT-MODE-WITH-DISCOVERABLE-TARGETS
Rationale: canvas panと配置編集の区別を保ちながら、対象と有効化方法を初めての利用者にも示す。
Must preserve: OFFの従来panと明示toggle、ONのtarget drag、supported targetへのhelp/hover説明とtoolbarの常設説明、keyboard focus/touchで同じ説明へ到達すること、feature/label/legend個別編集とShift/Ctrlの優先順位、record/alignmentの既存動作、実dragの1 History操作、Session/Exportを維持する。hintだけでcanonical値を変えず、成果物へhintを保存しない。modeを自動ONにしない。
May retire: OFF/ONの意味を区別できないcursorと説明不足だけ。gesture/編集機能/保存意味は退役しない。
Accepted residual risk: 移動前に有効化の1手順が残る。常設説明とkeyboard/touch検証で発見可能性を補う。
Owner: satoshikawato
Decision date: 2026-09-26
```

```json
{
  "concern": "web.layout-edit-affordance",
  "scenarioRevision": 1,
  "choice": "A / EXPLICIT-MODE-WITH-DISCOVERABLE-TARGETS",
  "rationale": "canvas panと配置編集の区別を保ちながら、対象と有効化方法を初めての利用者にも示す。",
  "mustPreserve": "OFFの従来panと明示toggle、ONのtarget drag、supported targetへのhelp/hover説明とtoolbarの常設説明、keyboard focus/touchで同じ説明へ到達すること、feature/label/legend個別編集とShift/Ctrlの優先順位、record/alignmentの既存動作、実dragの1 History操作、Session/Exportを維持する。hintだけでcanonical値を変えず、成果物へhintを保存しない。modeを自動ONにしない。",
  "mayRetire": "OFF/ONの意味を区別できないcursorと説明不足だけ。gesture/編集機能/保存意味は退役しない。",
  "acceptedResidualRisk": "移動前に有効化の1手順が残る。常設説明とkeyboard/touch検証で発見可能性を補う。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-26"
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

### PD-OI-059: Circular definition settings apply on Generate

- Concern key: `web.circular.definition-application`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `A / CIRCULAR-DEFINITION-APPLIES-ON-GENERATE`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its nine-field JSON representation below.
- Decision source: the complete `D-04` receipt in [`02_DECISION_PACK.md`](./web-gui-audit-20260930/02_DECISION_PACK.md) at P00 merge commit
  `e97d90fecfb327135eb50e85cfa9a87be3145823`, approved by `satoshikawato` on `2026-09-30` through the two
  Owner replies quoted verbatim in the Revision 29 entry above. The receipt and
  JSON below reproduce all nine supplied fields without translation or
  additional terms. This record does not supersede another decision.
  Dependent runtime requires this authority merged into its base; this
  amendment supplies no runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `6380086684b70c8be665093da9549d884082413cbec22de6c26bd1238dc0ce8c`.
- Acceptance contracts: `OIC-009`, `OIC-010`, `OIC-024`, `OIC-027`. These obligations and the complete
  selected outcome are jointly required; their citation does not claim
  completed dependent-runtime checks.

```text
PRODUCT_DECISION
Concern: web.circular.definition-application
Scenario revision: 1
Choice: A / CIRCULAR-DEFINITION-APPLIES-ON-GENERATE
Rationale: Result の定義行に、crop の長さ・GC%・record label と食い違う値が書き込まれないようにする。Linear と同じ「Applies on Generate」に揃える。
Must preserve: Species、Strain、Plot title、Title position、Title font、Default font size、Keep Full Definition の編集・保存・History、Generate 後の正しい定義（region の長さと GC%、record label と subtitle、逆相補、grid の順序）、Linear の現在の挙動、他の即時編集（色、ラベル、凡例など）。各設定には Applies on Generate の表示を付ける。
May retire: 上の Circular の設定の即時反映と、そのための helper（regenerate_definition_svgs）。
Accepted residual risk: これらを変えたとき、Generate するまで preview が変わらない。
Owner: satoshikawato
Decision date: 2026-09-30
```

```json
{
  "concern": "web.circular.definition-application",
  "scenarioRevision": 1,
  "choice": "A / CIRCULAR-DEFINITION-APPLIES-ON-GENERATE",
  "rationale": "Result の定義行に、crop の長さ・GC%・record label と食い違う値が書き込まれないようにする。Linear と同じ「Applies on Generate」に揃える。",
  "mustPreserve": "Species、Strain、Plot title、Title position、Title font、Default font size、Keep Full Definition の編集・保存・History、Generate 後の正しい定義（region の長さと GC%、record label と subtitle、逆相補、grid の順序）、Linear の現在の挙動、他の即時編集（色、ラベル、凡例など）。各設定には Applies on Generate の表示を付ける。",
  "mayRetire": "上の Circular の設定の即時反映と、そのための helper（regenerate_definition_svgs）。",
  "acceptedResidualRisk": "これらを変えたとき、Generate するまで preview が変わらない。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-30"
}
```

### PD-OI-060: Diagram-wide stroke settings apply on Generate

- Concern key: `web.stroke.application`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `C / STROKE-APPLIES-ON-GENERATE`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its nine-field JSON representation below.
- Decision source: the complete `D-05` receipt in [`02_DECISION_PACK.md`](./web-gui-audit-20260930/02_DECISION_PACK.md) at P00 merge commit
  `e97d90fecfb327135eb50e85cfa9a87be3145823`, approved by `satoshikawato` on `2026-09-30` through the two
  Owner replies quoted verbatim in the Revision 29 entry above. The receipt and
  JSON below reproduce all nine supplied fields without translation or
  additional terms. This record does not supersede another decision.
  Dependent runtime requires this authority merged into its base; this
  amendment supplies no runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `f95d7c6e1b7354f805d186d2709637b246ef2f15d416bc1143e1c46816ef14af`.
- Acceptance contracts: `OIC-024`, `OIC-027`. These obligations and the complete
  selected outcome are jointly required; their citation does not claim
  completed dependent-runtime checks.

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

### PD-OI-061: Legend rename collisions for every entry

- Concern key: `web.legend.rename-collision`
- Scenario revision: `2`
- Supersedes: `PD-OI-061`, scenario revision `1` (`A / DIALOG-FOR-ALL-ENTRIES`).
- Status: `ACCEPTED`
- Selected outcome: `B / MERGE-SAME-FEATURE-TYPE-ONLY`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its nine-field JSON representation below.
- Decision source: the Owner replies of `2026-10-06` quoted verbatim in the
  Revision 32 entry above. The Owner approved the receipt text below as
  written ("OKです"). Dependent runtime (OV-62) is changed in the same pull
  request, by the static Product Contract co-change route; this record supplies
  no runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `783436b02c70f043f89dc07b4bf8d56e0869e996507f72fedbbecd0301c8f322`.
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
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `A / LIVE-EDIT-EQUALS-REGENERATION`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its nine-field JSON representation below.
- Decision source: the complete `D-11` receipt in [`02_DECISION_PACK.md`](./web-gui-audit-20260930/02_DECISION_PACK.md) at P00 merge commit
  `e97d90fecfb327135eb50e85cfa9a87be3145823`, approved by `satoshikawato` on `2026-09-30` through the two
  Owner replies quoted verbatim in the Revision 29 entry above. The receipt and
  JSON below reproduce all nine supplied fields without translation or
  additional terms. This record does not supersede another decision.
  Dependent runtime requires this authority merged into its base; this
  amendment supplies no runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `9298e2bb07d617f9f2562e8168d412a867082f4130b60ec709ef2d7ff005e9e8`.
- Acceptance contracts: `OIC-024`, `OIC-027`. These obligations and the complete
  selected outcome are jointly required; their citation does not claim
  completed dependent-runtime checks.

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

### PD-OI-067: Per-record inferred Linear definitions

- Concern key: `web.linear.file-default-definition`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `B / PER-RECORD-INFERRED-DEFINITION`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its nine-field JSON representation below.
- Decision source: the complete `D-12` receipt in [`02_DECISION_PACK.md`](./web-gui-audit-20260930/02_DECISION_PACK.md) at P00 merge commit
  `e97d90fecfb327135eb50e85cfa9a87be3145823`, approved by `satoshikawato` on `2026-09-30` through the two
  Owner replies quoted verbatim in the Revision 29 entry above. The receipt and
  JSON below reproduce all nine supplied fields without translation or
  additional terms. This record does not supersede another decision.
  Dependent runtime requires this authority merged into its base; this
  amendment supplies no runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `daf932cb4892019717ea3ae082d1c799b58a436f909738fbc6239df5d43d265b`.
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
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `A / APPLY_ON_GENERATE`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its nine-field JSON representation below.
- Decision source: the Owner answers of `2026-10-04` quoted verbatim in the
  Revision 30 entry above (`OD-2`, `OD-3`, `OD-4`), and Appendix A of
  [`POPUP_RECORD_ACTIONS_AND_VIBRIO_SESSION_PLAN_2026-10-04.md`](./POPUP_RECORD_ACTIONS_AND_VIBRIO_SESSION_PLAN_2026-10-04.md).
  The receipt fields restate those answers and the plan presented in that
  session; the Owner did not separately review the receipt wording. This record
  does not supersede another decision; `PD-OI-032` retains its scope. Dependent
  runtime requires this authority merged into its base; this amendment supplies
  no runtime acceptance evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `38f89876744c098aa5695b570a0250c6f672272594aafc6c1e1e739b69712f40`.
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

### PD-OI-086: Mode-scoped diagram settings and edits

- Concern key: `web.mode.scoped-settings`
- Scenario revision: `1`
- Status: `ACCEPTED`
- Selected outcome: `A / ALL-DIAGRAM-SETTINGS-PER-MODE`
- Normative outcome: exactly the complete approved `PRODUCT_DECISION` receipt
  and its nine-field JSON representation below.
- Decision source: the Owner reply of `2026-10-07` quoted verbatim in the
  Revision 33 entry above. The Owner approved the receipt text below as
  written ("OK、この文面で"). The receipt and JSON below reproduce all nine
  fields without translation or additional terms. This record does not
  supersede another decision. `PD-OI-061`, `PD-OI-062`, `PD-OI-063`, and
  `PD-OI-084` retain their scope within each mode, with `PD-OI-063` read per
  mode, and `PD-OI-002` within Linear. The Revision 33 entry states these
  readings and the scope of GUI remediation S00 decision 2 that the receipt
  retires outside this Contract. Dependent runtime requires this authority
  merged into its base; this amendment supplies no runtime acceptance
  evidence.
- Receipt SHA-256 (UTF-8, excluding the final newline):
  `ee7ef6a59d18e149670340a3b2aac94c76f036e3f72029b74767a3dc30d9dd92`.
- Acceptance contracts: `OIC-004`, `OIC-027`. These obligations and the
  complete selected outcome are jointly required; their citation does not
  claim completed dependent-runtime checks.

```text
PRODUCT_DECISION
Concern: web.mode.scoped-settings
Scenario revision: 1
Choice: A / ALL-DIAGRAM-SETTINGS-PER-MODE
Rationale: Circular と Linear で、図の設定と編集をすべて別々に持つ。片方のモードで行った設定や編集が、もう片方のモードの図や Generate を変えたり失敗させたりしないようにする（OV-80、OV-82 ほか）。v0.15.0 の drawing はこの形をそのまま引き継ぐ。
Must preserve: 次のものは Circular と Linear で別々に持ち、モードを何度切り替えても残る。切り替えは他方のモードの値を消さず、写さない。図の設定（Depth、GC content と GC skew、ラベル、軸と目盛り、フォント、線、トラック、比較の閾値、LOSAT の検索の設定を含むすべての設定）、palette と色、色ルール（「この feature だけ」の色を含む）、qualifier priority、ラベルの表と filter、注釈セット、凡例の編集（色・線・名前・削除・追加・順序）、feature ごとの編集（表示、ラベル、塗り、線、配置）、record の表示、canvas の余白、画面にない設定の上書き。両方のモードで共通のままのものは、入力ファイル（今までどおりモードごと）、LOSAT の結果の cache、LOSAT の実行方法とスレッド数、Auto Reflow・PNG DPI・palette の Instant Preview などのアプリの設定、Session の title。初めて使うモードは既定値で始まり、凡例と feature の編集はない。凡例の編集は表示中の Result のモードの値を変える。Show Depth はそのモードに最初の Depth ファイルが入ると On、最後のファイルがなくなると Off になり、他方のモードは変わらない。Generate、即時の編集、Session の保存と読み込み、export は、そのモードの値だけを使う（PD-OI-066）。PD-OI-061、062、063、084 は各モードの中で今までどおり。PD-OI-002 の LOSATP の上限は Linear の中で今までどおり。Undo/Redo は今までどおり 1 つの履歴で、モードの切り替えも 1 step。Reset Settings は今までどおり両方のモードを既定に戻し、ファイル、Depth の割り当て、Result は残す（PD-OI-070）。Session は両方のモードの値を保存する（Session 46）。以前の Session（27〜44）は Python を起動せずに読み込み（PD-OI-044）、次の規則で両方のモードに分ける。両モード共通だった設定は両方のモードに写す。モードごとだった値（mode profile、片方のモード専用の設定、Circular と Linear で分かれていた設定、モード付きの行）はそのモードに入れ、mode profile にない他方のモードの値は既定値にする（S00 の判断 1・2 と同じ）。Show Depth はそのモードに Depth ファイルがあるときだけ On にする。凡例の編集と、保存された Result に結びついた feature ごとの編集は、その Result のモードだけに入れる。record を選んだ注釈は、その record を選んだモードだけに入れる。保存された値は失わない（OIPC-C05、C06、OIC-004）。CLI と Python API の描画と再現は今までどおり。
May retire: 図の設定と編集が両方のモードで共有される動作と、S00 の判断 2 の「モードごとにするのは title と font だけ」という範囲。モードの切り替えで設定を入れ替える仕組み（mode profile）と、他方のモードの Show Depth を Off にする動作。あるモードの凡例・feature の編集・色ルール・注釈が、別のモードの Generate に持ち込まれる動作。開発版だけの Session 45 の読み込み。
Accepted residual risk: 両方のモードで同じ設定や色にしたいときは、両方のモードで設定する（コピーする操作は v0.15.0 の drawing で入る）。以前の Session から両方のモードに写した値（Depth の最小・最大、window など）が、他方のモードのデータに合わないことがある。開発版で保存した Session 45 は読み込めない。変更が大きく、0.14.0 のリリース前の検証の期間が短くなる。
Owner: satoshikawato
Decision date: 2026-10-07
```

```json
{
  "concern": "web.mode.scoped-settings",
  "scenarioRevision": 1,
  "choice": "A / ALL-DIAGRAM-SETTINGS-PER-MODE",
  "rationale": "Circular と Linear で、図の設定と編集をすべて別々に持つ。片方のモードで行った設定や編集が、もう片方のモードの図や Generate を変えたり失敗させたりしないようにする（OV-80、OV-82 ほか）。v0.15.0 の drawing はこの形をそのまま引き継ぐ。",
  "mustPreserve": "次のものは Circular と Linear で別々に持ち、モードを何度切り替えても残る。切り替えは他方のモードの値を消さず、写さない。図の設定（Depth、GC content と GC skew、ラベル、軸と目盛り、フォント、線、トラック、比較の閾値、LOSAT の検索の設定を含むすべての設定）、palette と色、色ルール（「この feature だけ」の色を含む）、qualifier priority、ラベルの表と filter、注釈セット、凡例の編集（色・線・名前・削除・追加・順序）、feature ごとの編集（表示、ラベル、塗り、線、配置）、record の表示、canvas の余白、画面にない設定の上書き。両方のモードで共通のままのものは、入力ファイル（今までどおりモードごと）、LOSAT の結果の cache、LOSAT の実行方法とスレッド数、Auto Reflow・PNG DPI・palette の Instant Preview などのアプリの設定、Session の title。初めて使うモードは既定値で始まり、凡例と feature の編集はない。凡例の編集は表示中の Result のモードの値を変える。Show Depth はそのモードに最初の Depth ファイルが入ると On、最後のファイルがなくなると Off になり、他方のモードは変わらない。Generate、即時の編集、Session の保存と読み込み、export は、そのモードの値だけを使う（PD-OI-066）。PD-OI-061、062、063、084 は各モードの中で今までどおり。PD-OI-002 の LOSATP の上限は Linear の中で今までどおり。Undo/Redo は今までどおり 1 つの履歴で、モードの切り替えも 1 step。Reset Settings は今までどおり両方のモードを既定に戻し、ファイル、Depth の割り当て、Result は残す（PD-OI-070）。Session は両方のモードの値を保存する（Session 46）。以前の Session（27〜44）は Python を起動せずに読み込み（PD-OI-044）、次の規則で両方のモードに分ける。両モード共通だった設定は両方のモードに写す。モードごとだった値（mode profile、片方のモード専用の設定、Circular と Linear で分かれていた設定、モード付きの行）はそのモードに入れ、mode profile にない他方のモードの値は既定値にする（S00 の判断 1・2 と同じ）。Show Depth はそのモードに Depth ファイルがあるときだけ On にする。凡例の編集と、保存された Result に結びついた feature ごとの編集は、その Result のモードだけに入れる。record を選んだ注釈は、その record を選んだモードだけに入れる。保存された値は失わない（OIPC-C05、C06、OIC-004）。CLI と Python API の描画と再現は今までどおり。",
  "mayRetire": "図の設定と編集が両方のモードで共有される動作と、S00 の判断 2 の「モードごとにするのは title と font だけ」という範囲。モードの切り替えで設定を入れ替える仕組み（mode profile）と、他方のモードの Show Depth を Off にする動作。あるモードの凡例・feature の編集・色ルール・注釈が、別のモードの Generate に持ち込まれる動作。開発版だけの Session 45 の読み込み。",
  "acceptedResidualRisk": "両方のモードで同じ設定や色にしたいときは、両方のモードで設定する（コピーする操作は v0.15.0 の drawing で入る）。以前の Session から両方のモードに写した値（Depth の最小・最大、window など）が、他方のモードのデータに合わないことがある。開発版で保存した Session 45 は読み込めない。変更が大きく、0.14.0 のリリース前の検証の期間が短くなる。",
  "owner": "satoshikawato",
  "decisionDate": "2026-10-07"
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
| `OIC-026` | `PD-OI-035` and `PD-OI-039`: identity, keyboard Select/Skip, non-rendered candidates, no position-only selection, desktop canvas, focus and transient overlay exclusion remain required with all PD-OI-031/034 outcomes. Compact review retains visible, operable canvas at full available width and at least 200 px height at 390×844/740, scrollable candidates and reachable Apply/Cancel, local no-Worker draft edits, atomic batch validation, failure/error/retry and artifact/orientation/History recovery. Narrow review closes Editor through its owner while retaining tab, disables reopening with a reason until review ends, then permits explicit reopen; wide drag remains. Browser verification must show presentation changes leave draft and Result unchanged. |
| `OIC-027` | `PD-OI-066`: a live-edited Result agrees with a Result freshly generated from the same draft in the meaning (position, color, text, and visibility) of every edited element, including after Session load and export. Settings shown as Applies on Generate do not change the Result before Generate. A live edit that the Generate compiler cannot reproduce is shown as Applies on Generate instead. |
| `OIC-028` | `PD-OI-088`: the color scope, Legend name, and Reset fill color dialogs, whether opened from the feature popup's fill color, the Features drawer color input, a popup stroke edit, or a Legend panel rename (Merge/Suffix), open from the saved rules before the rules a choice may add are prepared. A choice is one History step and Cancel records none. From a choice until its step ends, the dialog stays open, its choices and Cancel are disabled, it states that the edit is applying, and no second choice, Cancel, Escape, or backdrop click starts or closes anything; the dialog closes when the step ends. |

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
