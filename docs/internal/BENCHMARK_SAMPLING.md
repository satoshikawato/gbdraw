# Benchmark sampling

Why the timing comparisons sample the way they do, and the measurements behind
each number. The rules live in code. Change a number there and here in the same
pull request.

| Tool | Where the rule is |
| --- | --- |
| `tools/benchmark_protein_comparison.py` | `POLICY`, `PROCESS_POLICY`, `spread()`, `compare_reports()`, subcommand `compare-sources` |
| `tools/benchmark_diagram_layout.py` | `DEFAULT_THRESHOLD_PCT`, `compare_reports()`, subcommand `compare-sources` |
| `tools/measure_gallery_publication_performance.py` | three trials per side; `baseSha` and `maximumRegressionRatio` in `tests/performance/gallery-publication-baseline.json`; run by `.github/workflows/gallery-publication.yml` |

## Terms

- Noise: 100 × MAD / median, where MAD is the median absolute deviation from
  the median. One outlier among three samples does not move it.
- CV: standard deviation / mean, in percent.
- False regression rate: how often a comparison of the same code with itself
  reads "regression".
- Slowdown caught: the slowdown that a comparison reads as a regression 80% of
  the time.
- Model: the logarithm of a time is a within-process part (σw) plus a
  between-process part (σb). The spread of a ratio of medians follows from
  both.

## Three samples per side, each from a fresh process, taking turns

`compare-sources` in both local tools runs base, current, base, current, base,
current, one sample per process.

Samples taken in one process share its state, so most of the spread is
between processes. More samples in one process barely narrow the median; the
same three samples from separate processes do. Taking turns spreads a slow
period of the machine over both sides.

Measured from 104 saved protein reports (September 2026, the maintainer's
workstation under ordinary load; σw 0.035, σb 0.034 to 0.064). Per stage, for
σb 0.064 (0.034 in brackets):

| Sampling | Threshold for 5% false regressions | Slowdown caught | Observed false regressions (difference over 10%) |
| --- | --- | --- | --- |
| 3 samples in one process (the earlier rule) | 17.4% (10.4%) | 27.4% (16.1%) | 13.4% |
| 7 samples in one process | 16.6% (9.2%) | 26.2% (14.2%) | not measured |
| 21 samples in one process | 16.2% (8.6%) | 25.6% (13.2%) | not measured |
| **3 processes, taking turns (the rule)** | 10.3% (6.8%) | 16.0% (10.4%) | 6.5% |
| 5 processes, taking turns | 7.9% (5.2%) | 12.2% (8.0%) | not measured |

The observed column draws three reports per side from the 16 cases that had
six or more reports of the same code. A comparison reads "regression" when any
of its stages does, so the rate per comparison is higher than the rate per
stage.

The count stays three, the project's rule for timing comparisons. Five
processes would lower the threshold for 5% false regressions to about 8%, at
two thirds more run time per comparison.

## Noise limit 10% in protein `compare-sources`

A stage whose noise exceeds the limit on either side is "inconclusive" and not
judged on time. Samples from separate processes carry the spread between
processes, so the 5% limit of one-process comparisons left many same-code
comparisons inconclusive. Per comparison, from the same saved reports (the four
cases with enough reports):

| Sampling | Noise limit | Pass | False regression | Inconclusive |
| --- | --- | --- | --- | --- |
| 3 samples in one process | 5% | 67% | 24% | 8% |
| 3 processes, taking turns | 5% | 49% | 10% | 42% |
| 3 processes, taking turns | 7.5% | 70% | 11% | 19% |
| **3 processes, taking turns** | **10%** | 80% | 12% | 8% |

10% equals the regression threshold: a stage goes unjudged only when its
spread is larger than the difference the comparison looks for. One-process
comparisons keep 5%. Each report saves its policy, and `compare` judges by the
policy saved in the report, so reports saved earlier keep their limit.

## Protein stages under 2 ms not judged on time

In `compare-sources`, a stage whose base and current medians are both under
2 ms reads `pass_below_floor`. Its time is recorded and its result is still
compared, but its time difference decides nothing.

Measured on the quiet workstation on 2026-11-06 (dev `d1ee38719`; four blocks
with a one-minute load average of at most 1.7 and no other heavy job), the
spread between processes (robust σ) by stage length:

| Stage median | Robust σ between processes |
| --- | --- |
| 100 ms or more | 1.5% |
| 10 to 100 ms | 2.3% |
| 1 to 10 ms | 3.4% |
| under 1 ms | 4.4% |

Per stage, 1 of 116 same-code pairs differed by more than 10%. A comparison
has 29 stages and any one decides, so the model gave pass 88%, regression 10%
and inconclusive 2% per comparison. Every miss came from three `path-16` stages
with medians of 0.18, 0.51 and 0.94 ms; at 0.18 ms a delay of 0.03 ms is
already over 15%. With the floor the model passes 97% of same-code
comparisons, and all five measured blocks pass.

Stage medians jump from 0.97 ms to 3.7 ms, so a boundary at 2 ms does not move
a stage from one side to the other between runs. A stage that reaches 2 ms on
either side is judged again, so a tenfold slowdown of a 0.2 ms stage is still
judged.

## Layout `compare-sources`: 5% threshold

A case is a regression when the current median exceeds the base median by more
than 5% and by at least the larger noise of the two sides. A difference over 5%
but smaller than that noise reads `noise_exempt`.

Measured on the quiet workstation, same code: spread between processes 2.2% to
3.2% over the three cases; differences −7.1% to +1.5%; no regression in 12
pairs (none in 15 with a block under load). The model gives 1% to 5% false
regressions per case and about 8% per three-case comparison, so the threshold
stays 5%.

## Gallery projection gate: fixed base, three trials each, ratio 1.25

The `Gallery publication` workflow measures the fixed base commit (`baseSha`)
and the head in the same job, three trials each, base first. It fails when the
head's median wall time exceeds 1.25 times the base median, or the head's
largest peak RSS exceeds 1.25 times the base's largest.
Each revision projects its own file of the Vibrio Gallery Session
(`vibrio-harveyi-group-collinear.gbdraw*`): a base before #1760 stores it as
`.gbdraw-session.json.gz`, later commits as a `.gbdraw` file. The report lists
the command of each trial.

Measured in CI over 73 runs (2026-10-31 to 2026-11-05): trials within one job
had noise 0.5% (σ 0.8%), but base medians across runners fell into three groups
from 1.30 to 2.09 s (CV 11.6%). Measuring the base in the same job cancels the
runner difference; a fixed time limit would not.

The base before `8b8afec97` was `574b33b8`, which was slower than current code
(head / base 0.595, range 0.56 to 0.64), so 1.25 let current code become about
2.1 times slower before the gate failed. With a recent base, 1.25 means a 25%
slowdown again. The model gives no false failures at 1.25 with the CI spread.
On the quiet workstation (σ 2.7% within a job) it gives none at 1.10 and 2.5%
at 1.05. Head peak RSS in the 73 CI runs was 571,604 to 683,108 KiB (largest
over smallest within one run up to 1.19); against the new base, the RSS limit
went from about 2.4 times to 1.25 times current code.

## Fixed limits in the Playwright performance specs

The Vibrio save time, the longest heartbeat gap during that save, the longest
task during mouse-wheel input on a large preview, and the search times have fixed limits. They stay as
guards against stalls, not as slowdown checks: CI runners fall into three speed
groups up to 1.6 times apart, so tighter fixed limits would fail on the slow
ones. A check meant to catch small slowdowns measures its base in the same job,
as the Gallery projection gate does.

## Measured spread

| Measurement | Where | Count | Noise (median) | CV (median) |
| --- | --- | --- | --- | --- |
| Protein stage under 10 ms, within one process | workstation, saved reports, September 2026 | 168 stages | 3.6% | 6.3% |
| Protein stage 10 to 100 ms, within one process | same | 113 stages | 2.0% | 3.4% |
| Protein stage 0.1 to 1 s, within one process | same | 249 stages | 1.4% | 2.6% |
| Protein stage 1 s or more, within one process | same | 76 stages | 1.3% | 3.0% |
| Gallery projection, base trials within one job | CI, 2026-10-31 to 2026-11-05 | 73 runs | 0.5% (largest 2.2%) | not computed |
| Gallery projection, base medians across runs | same | 73 runs | not computed | 11.6% |
| Vibrio save wall time (median 5,816 ms) | CI performance job | 70 runs | 9.9% | 15.9% |
| Vibrio save wall time (median 3,100 ms) | quiet workstation, 2026-11-06 | 4 runs | not computed | 2.4% |
| Longest heartbeat gap during the Vibrio save (median 250 ms) | CI performance job | 70 runs | 9.3% | 13.4% |
| Longest task during mouse-wheel input on a large preview | CI performance job | 70 runs, 210 measurements | every measurement 0 ms | every measurement 0 ms |

## Re-measuring

Compare a tree with itself on a quiet machine, several times, and count the
decisions:

```bash
python tools/benchmark_protein_comparison.py compare-sources \
  --baseline-root <tree> --current-root <tree> --output-dir <folder>
python tools/benchmark_diagram_layout.py compare-sources \
  --baseline-root <tree> --current-root <tree> --output-dir <folder>
```

Use two separate checkouts of the same commit, so that each side imports its
own `gbdraw`. A rule change needs measurements of the kind above for the
sampling it replaces and the one it proposes.
