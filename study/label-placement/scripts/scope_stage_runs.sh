#!/bin/bash
# Stage check runs (scope study), 4 processes; output logs/scope-stage-20261013.jsonl
cd /home/kawato/gbdraw-baselines/label-placement-20261011
export TMPDIR=/home/kawato/gbdraw-baselines/tmp PYTHONPATH=/home/kawato/gbdraw-work/.worktrees/lp-scope
printf '%s\n' c_t5_both_f14 c_t5_both_f14:f10 c_tobacco_both_f14 c_tobacco_both_f14:f10 c_tobacco_both_def c_p1_both_f14 c_p1_both_f14:f12 c_p1_both_f14:f11 \
  c_epifagus_both_def c_epifagus_both_def:f9 c_epifagus_both_def:f12 c_plasmid37_both_def c_plasmid37_both_def:f10 c_mjenmv_both_def c_mjenmv_both_def:f6.5 \
  c_lambda_both_def c_wssv_both_def c_mjenmv_ex_both c_wssv_ex_both c_p1_both_def |
  xargs -P 4 -I{} sh -c 'mkdir -p work/SCOPE/stage/$(echo {} | tr ":" "_") && cd work/SCOPE/stage/$(echo {} | tr ":" "_") && /home/kawato/micromamba/bin/python3.13 -P ../../../../scripts/scope_stage.py {} 2>/dev/null | tail -1' >> logs/scope-stage-20261013.jsonl
