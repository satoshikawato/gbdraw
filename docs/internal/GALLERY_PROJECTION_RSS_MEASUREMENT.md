# Gallery projection RSS measurement

S04 release-gate remediation, 2026-09-08. The projection benchmark uses
`--no-incremental-marking` for both its fixed baseline and candidate processes.
The application, fixture, 1.25 ratio, three-trial count, wall-time median, and
whole-process maximum-RSS assertion are unchanged. Projection always measures
the fixed baseline locally; historical Node 18 measurements cannot stand in for
that baseline under a new Node/GC profile. Each completed trial is printed, and
a failed gate writes its report before returning a nonzero exit code.

## Source and environment

| Revision | SHA | Production difference from preceding row |
| --- | --- | --- |
| Fixed benchmark base | `574b33b83962949397839e2aaa862a8b96667625` | Historical benchmark owner |
| A, S02 merged | `d7345f5e280366835f16f2c4295cba763acaf14e` | Not part of this remediation's production comparison |
| B, S03 merged | `67bab304ef6d990a1e69f9713937def30893eb7b` | None; only replay tests and their fixture changed |
| C, S04 merged | `d5c5e1b0296d6b402fbc41cee700aa5458b01250` | `gbdraw/web/js/app/ui.js`, five additions and one deletion |

The existing tool's `_operation_command('projection')` runs
`tests/web/session-request.test.mjs --project-session` against the existing
Vibrio Gallery gzip session. It copies JS files into a temporary ES-module
package. An ESM load hook observed 33 production modules; `ui.js` was copied but
never loaded. All loaded production modules, the test entry, input, measurement
tool, workflow and fixed baseline have identical bytes across A/B/C. The
[module inventory and SHA-256 values](gallery_projection_rss_evidence.json)
record that closure. The fixed benchmark revision has its own historical
modules and input; it is measured as committed, not replaced with C's sources.

Comparisons were serial, isolated exports of each revision with all production
JS, the original test entry and input. Initial and confirmation comparisons
used three prespecified blocks of three trials, rotating revision order. No
trial was excluded. Local Python was 3.11.14, Node 20.20.2 / V8
11.3.244.8-node.38, WSL2 Linux x86_64, glibc 2.39, with no NODE_OPTIONS or
NODE_PATH. CI uses Ubuntu 24.04.4 image 20260831.293.1, Python 3.11.16 and the
same Node 20.20.2. Python runs the standard-library measurement controller;
the child uses Node built-ins and the recorded repository modules, with no
Python packages, npm dependencies, browser, wheel or Worker. The workflow's
editable Python export dependencies are not imported by this child. Node 24
used to execute GitHub Actions is distinct from the selected Node 20 workload.

## Observations and intervention

All numbers below are KiB. The historical description of approximately
1.294 to 1.016 refers to millions of KiB, not GiB.

The retained run 34088103527 on the same SHA
`2efb21e57a5db81344930c54928fac5654e73b2b` failed at 1,293,696 KiB on attempt 1
and passed at 1,016,492 KiB on its authorized attempt 2. Failed historical runs
have no individual trial values: the former reporter raised before emitting
them. The original S04 attempts likewise remain FAIL at 1,294,136 and
1,294,532 KiB. These historical observations motivate investigation; they are
not a diagnosis by themselves.

| Original workload, initial nine trials | RSS minimum–maximum | Wall median (seconds) |
| --- | --- | --- |
| Fixed base | 1,006,244–1,016,328 | 1.260 |
| A | 1,005,172–1,016,772 | 1.295 |
| B | 999,616–1,016,368 | 1.278 |
| C | 1,006,508–1,020,748 | 1.322 |

Adding only `--trace-gc` to the original workload reproduced both modes on B
(1,294,976; 1,286,120; 1,016,764) and C (1,010,916; 1,016,496; 1,294,484).
Separately, checkpoints after read, decode, parse and projection reproduced
high RSS on A/B/C: 1,287,048–1,287,800 / 1,271,632–1,287,568 /
1,286,328–1,295,096 in three trials each. These probes affect timing and are
reported separately from uninstrumented measurements; they do not estimate
natural high-mode frequency. A one-CPU control instead gave approximately
1,083,000 KiB for A/B/C. Thus C's UI change is not needed to produce the modes.

The high peak occurs during JSON parsing, before the production projection
call. For example, C's checkpoint run ended decoding at 976,472 KiB maximum,
parsing at 1,294,884 KiB maximum, and projection at that same maximum. Current
RSS had already fallen to about 1,015,000 KiB when parsing returned. The decoded
session is approximately 285 MB; GC traces and Buffer counters show a varying
overlap of its temporary inflated Buffer with parser allocations. Measuring
only process-exit RSS would miss this peak.

The modes differ by about 27%, enough for the current three-trial maximum
ratio to reject unchanged executed production when the baseline samples the
low mode and a head trial samples the high mode. Nine low samples cannot
establish that the high mode is absent. We therefore control the measured GC
condition instead of changing the estimator or removing observations.

Disabling concurrent ArrayBuffer sweeping alone did not resolve the modes.
`--predictable-gc-schedule` also left high checkpoint measurements. Disabling
incremental marking resolved the measured mode split for both original and
checkpoint workloads. This is an intervention in the measurement process,
not an application memory fix. V8 documents incremental marking as an engine
flag in [Node 20.20.2's flag definitions](https://github.com/nodejs/node/blob/v20.20.2/deps/v8/src/flags/flag-definitions.h);
its [GC description](https://v8.dev/blog/trash-talk) explains the overlap of
incremental/concurrent GC with application execution.

| Confirmation, nine uninstrumented trials with the selected profile | RSS minimum–maximum | Wall median (seconds) |
| --- | --- | --- |
| Fixed base | 1,081,148–1,082,140 | 1.215 |
| A | 1,069,080–1,086,160 | 1.187 |
| B | 1,069,428–1,086,216 | 1.237 |
| C | 1,069,116–1,086,372 | 1.220 |

Three additional checkpoint trials per revision also remained between
1,068,828 and 1,085,868 KiB. All diagnostic trials, including unsuccessful
alternative controls, are retained in the linked JSON. The profile tests
allocation/time regressions under a shared GC condition; it does not claim to
benchmark default browser GC scheduling or prove all hosts deterministic.

## Regression detection and scope

A temporary 512 MiB allocation overlapping C's real projection produced
1,530,412 / 1,530,580 / 1,530,624 KiB and failed the unchanged maximum-RSS gate
against the contemporaneous 1,357,795 KiB limit. Combining just one of those
spikes with two control trials also failed. A one-second delay at the same
boundary failed the wall-time gate (median 2.78 > 1.79 seconds). The injected
allocation was released after projection; it was not retained until process
exit. These are positive controls in isolated diagnostic copies, not changes
to production or the committed workload.

Focused tests cover lossless reporting of RSS and time failures, mandatory
same-profile fixed-base measurement, a single RSS spike, and rejection of
trial counts other than three. The existing ratio assertions still apply.
Application owner/path counts do not change; the measurement controller
remains the sole gate owner. No CI authority, product behavior, public figure,
session fixture or fixed baseline was changed. Existing policy requires human
review because registered performance evidence paths change.
