"""Compose the study figures from the page config: figures/<case>.png and figures/<case>--crop<i>.png.

Usage: python compose_figures.py <page.json> <figures_dir> [case id ...]
<case>.png: every design of the case on one shared canvas (map centres aligned, one scale), captioned with the
design, font size, labels placed and dropped, and canvas size. <case>--crop<i>.png: the same local box of the map
at 1:1 under every design (the box of the page config, or the densest label region of the first design).
Prints one README table row per figure.
"""
import json, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import figs

cfg = json.loads(Path(sys.argv[1]).read_text())
out = Path(sys.argv[2]); out.mkdir(parents=True, exist_ok=True)
only = set(sys.argv[3:])
designs = {d['key']: d for d in cfg['designs']}
for case in cfg['cases']:
    cid = case['id']
    if only and cid not in only:
        continue
    panels, same = figs.load_case(cfg['render_dirs'], cid, designs)
    keys = list(panels)
    caption = lambda k: f"{designs[k]['label']} {designs[k]['en']} | {figs.summary(panels[k])}"
    cols = 3 if len(keys) > 4 else len(keys)
    w, _, _ = figs.union(list(panels.values()))
    scale = min(1.0, cfg.get('figure_cell_px', 1150) / w)
    figs.grid([(caption(k), im) for k, im in zip(keys, figs.aligned(list(panels.values()), scale))], cols).save(out / f'{cid}.png')
    print(f"| {cid}.png | {case.get('en', cid)} | {', '.join(designs[k]['label'] for k in keys)} | whole map, scale {scale:.2f} | |")
    for i, spec in enumerate(case.get('crops', []), 1):
        box = spec.get('box') or figs.densest_box(panels[spec.get('from', keys[0])]['labels'], *spec.get('size', (520, 380)))
        figs.grid([(f"{designs[k]['label']} {designs[k]['en']}", figs.crop_local(panels[k], box)) for k in keys], cols).save(out / f'{cid}--crop{i}.png')
        print(f"| {cid}--crop{i}.png | {case.get('en', cid)} | {', '.join(designs[k]['label'] for k in keys)} | 1:1, local box {list(box)} | |")
