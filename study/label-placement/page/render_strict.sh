#!/bin/bash
# LP-D07: C and G with the crossing fallback off (fallback=0) as the main panels; the defaults kept as *_cross on three cases.
P=/home/kawato/gbdraw-baselines/label-placement-20261011/page
PY=/home/kawato/micromamba/bin/python3.13
cd "$P" || exit 1
ALL=$($PY -c "import json; print(' '.join(c['id'] for c in json.load(open('page.json'))['cases']))")
OPEN="p1_open_f14 lambda_open_f14 pbr322_open_f14"
runs() { for c in $1; do printf '%s:%s ' "$c" "$2"; done; }
$PY run_extra.py --out renders --jobs 4 --suffix _cross $(runs "p1_f14 $OPEN" C) $(runs "p1_f14 $OPEN" G) > renders/batch_cross.log 2>&1
$PY run_extra.py --out renders --jobs 4 --opts fallback=0 $(runs "$ALL" C) > renders/batch_C_strict.log 2>&1
$PY run_extra.py --out renders --jobs 4 --opts fallback=0 $(runs "$ALL" G) > renders/batch_G_strict.log 2>&1
$PY run_extra.py --out renders --jobs 3 --opts fallback=0,inside=0 --suffix _in0 $(runs "$OPEN" G) > renders/batch_G_var.log 2>&1
$PY run_extra.py --out renders --jobs 3 --opts fallback=0,clear=6 --suffix _clear6 $(runs "$OPEN" G) >> renders/batch_G_var.log 2>&1
$PY run_extra.py --out renders --jobs 3 --opts fallback=0,merge=1 --suffix _merge $(runs "$OPEN" G) >> renders/batch_G_var.log 2>&1
$PY run_extra.py --out renders --jobs 3 --opts fallback=0,color=1 --suffix _color $(runs "$OPEN" G) >> renders/batch_G_var.log 2>&1
PYTHONPATH=/home/kawato/gbdraw-work/.worktrees/label-study TMPDIR=/home/kawato/gbdraw-baselines/tmp $PY backfill_gap.py renders > renders/batch_gap.log 2>&1
echo "strict renders done: $(grep -c 'seconds=' renders/batch_C_strict.log) C, $(grep -c 'seconds=' renders/batch_G_strict.log) G, $(grep -c 'seconds=' renders/batch_G_var.log) G variants, $(grep -c 'seconds=' renders/batch_cross.log) cross; failed: $(cat renders/batch_cross.log renders/batch_C_strict.log renders/batch_G_strict.log renders/batch_G_var.log | grep -c FAILED)"
