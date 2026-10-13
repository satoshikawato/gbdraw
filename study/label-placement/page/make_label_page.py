"""Build the Owner comparison page for the Circular label placement study.

Usage: python make_label_page.py <page.json> <out_dir>
Writes <out_dir>/index.html and the figures to <out_dir>/img/ (lossless WebP on a white ground).

The page is image-first (Owner 2026-10-12, LP-D03):
- per case, one canvas whose placement is switched in place: the panels are aligned on the map centre at one
  scale, so only the labels move;
- per case, 1:1 crops of the same part of the map for every design, side by side;
- per design, one contact sheet of that design on every case.
Renders are found as <render_dir>/<case id>__<design key>.{png,svg,metrics.json,labels.json} (scripts/run_case.py).
Text fields in the config are HTML (as in gallery-video-20261010/s4/pages/make_page.py).
"""
import html, json, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import figs

cfg = json.loads(Path(sys.argv[1]).read_text())
out = Path(sys.argv[2]); (out / 'img').mkdir(parents=True, exist_ok=True)
esc = html.escape
designs = {d['key']: d for d in cfg['designs']}
columns = cfg.get('metric_columns', [])
dirs = cfg['render_dirs']


def save(im, name):
    im.save(out / 'img' / name, format='WEBP', lossless=True, method=4)
    return name


def metric(panel, key):
    m = panel['metrics']
    if key == 'dropped':
        v = (m.get('proto') or {}).get('dropped')
        return len(v) if isinstance(v, list) else (v or 0)
    if key == 'canvas_px':
        return ' × '.join(str(v) for v in panel['canvas'])
    v = m.get(key)
    if key == 'leader_len_sum' and v is not None:
        return f'{round(v):,}'
    if isinstance(v, list):
        return ' / '.join(f'{x:g}' if isinstance(x, (int, float)) else str(x) for x in v)
    return '' if v is None else (f'{v:g}' if isinstance(v, float) else v)


def cell(key, value):
    cls = ''
    if key in cfg.get('zero_is_good', []) and value != '':
        cls = ' class="ok"' if value in (0, '0') else ' class="bad"'
    return f'<td{cls}>{esc(str(value))}</td>'


def crop_grid(cid, index, box, panels, caption):
    cells = ''
    for k, p in panels.items():
        im = figs.crop_local(p, box)
        name = save(im, f'{cid}--crop{index}--{k}.webp')
        cells += (f'<figure><div class="frame"><img src="img/{esc(name)}" width="{im.width}" height="{im.height}" loading="lazy" '
                  f'alt="{esc(designs[k]["label"])} {esc(designs[k]["name"])}: 同じ部分の原寸"></div>'
                  f'<figcaption><b>{esc(designs[k]["label"])}</b>{esc(designs[k]["name"])}</figcaption></figure>')
    return f'<div class="cropset">{f"<p>{caption}</p>" if caption else ""}<div class="crops">{cells}</div></div>'


loaded = {}
case_html = []
hero = ''
for case in cfg['cases']:
    cid = case['id']
    panels, same = figs.load_case(dirs, cid, designs)
    if not panels:
        raise SystemExit(f'no render for case {cid}')
    loaded[cid] = panels
    keys = list(panels)
    W, H, offsets = figs.union(list(panels.values()))
    first = case.get('default', keys[0])
    tabs = ''.join(
        f'<button type="button" role="radio" aria-checked="{"true" if k == first else "false"}" '
        f'tabindex="{"0" if k == first else "-1"}" data-k="{esc(k)}" id="{esc(cid)}-tab-{esc(k)}">'
        f'<b>{esc(designs[k]["label"])}</b><span>{esc(designs[k]["name"])}</span></button>' for k in keys)
    imgs = ''
    for (k, p), (ox, oy) in zip(panels.items(), offsets):
        name = save(p['img'], f'{cid}--{k}.webp')
        w, h = p['img'].size
        imgs += (f'<img data-k="{esc(k)}" src="img/{esc(name)}" width="{w}" height="{h}" loading="lazy" '
                 f'style="left:{ox / W * 100:.3f}%;top:{oy / H * 100:.3f}%;width:{w / W * 100:.3f}%" '
                 f'alt="{esc(case["title"])}: {esc(designs[k]["label"])} {esc(designs[k]["name"])}"{"" if k == first else " hidden"}>')
    notes = ''.join(
        f'<p data-k="{esc(k)}"{"" if k == first else " hidden"}><b>{esc(designs[k]["label"])}</b>{note}</p>'
        for k, note in case.get('notes', {}).items() if k in panels)
    table = ''
    if columns:
        head = ''.join(f'<th scope="col">{esc(label)}</th>' for _, label in columns)
        body = ''.join(
            f'<tr data-k="{esc(k)}"{" class=on" if k == first else ""}><th scope="row"><button type="button" data-k="{esc(k)}">'
            f'{esc(designs[k]["label"])}</button></th>' + ''.join(cell(key, metric(p, key)) for key, _ in columns) + '</tr>'
            for k, p in panels.items())
        table = (f'<div class="scroll"><table class="metrics"><thead><tr><th scope="col">配置</th>{head}</tr></thead>'
                 f'<tbody>{body}</tbody></table></div>')
    crops = ''
    for i, spec in enumerate(case.get('crops', []), 1):
        box = spec.get('box') or figs.densest_box(panels[spec.get('from', keys[0])]['labels'], *spec.get('size', (520, 380)))
        block = crop_grid(cid, i, box, panels, spec.get('caption', ''))
        crops += block
        if cfg.get('hero', {}).get('case') == cid and cfg['hero'].get('crop', 1) == i:
            hero = block.replace('<div class="cropset">', '<div class="cropset hero">', 1)
    case_html.append(f'''<section class="case" id="{esc(cid)}" data-first="{esc(first)}">
  <header><h3>{esc(case["title"])}</h3><p class="meta">{case.get("meta", "")}</p>{f'<p>{case["note"]}</p>' if case.get("note") else ''}</header>
  <div class="bar"><div class="tabs" role="radiogroup" aria-label="{esc(case["title"])} の配置を選ぶ">{tabs}</div>
    <button type="button" class="hold" aria-label="押している間 {esc(designs[keys[0]]["label"])} を表示">押している間 {esc(designs[keys[0]]["label"])}</button></div>
  <div class="frame"><div class="stack" style="--w:{W:.0f}px;aspect-ratio:{W:.0f}/{H:.0f}">{imgs}</div></div>
  <div class="notes" aria-live="polite">{notes}</div>
  {f'<p class="same">このケースでは今の配置と同じ結果になるため図を省略: {esc("、".join(designs[k]["label"] + " " + designs[k]["name"] for k in same))}</p>' if same else ''}
  {table}
  {crops}
</section>''')

sheets = ''
if cfg.get('sheets', True):
    scale = cfg.get('sheet_scale', 0.3)
    for k, d in designs.items():
        cells = []
        for case in cfg['cases']:
            panels = loaded[case['id']]
            if k in panels:
                index = list(panels).index(k)
                cells.append((case.get('en', case['id']), figs.aligned(list(panels.values()), scale)[index]))
        if not cells:
            continue
        sheet = figs.grid(cells, cfg.get('sheet_cols', 4), pad=10, caption_px=17)
        name = save(sheet, f'sheet--{k}.webp')
        sheets += (f'<figure class="sheet" id="sheet-{esc(k)}"><figcaption><b>{esc(d["label"])}</b>{esc(d["name"])}'
                   f'<span>{d["desc"]}</span></figcaption><div class="frame"><img class="zoom" style="--w:{sheet.width}px" '
                   f'src="img/{esc(name)}" width="{sheet.width}" height="{sheet.height}" loading="lazy" '
                   f'alt="{esc(d["label"])} {esc(d["name"])} を全ケースに適用した一覧"></div></figure>')

def picture_groups(groups, prefix):
    """One case under several placements or settings: whole maps on one shared canvas, or one map region at 1:1."""
    html = ''
    for vi, group in enumerate(groups, 1):
        loaded_v = [(v, figs.load(figs.find_render(dirs, group['case'], v['render']))) for v in group['panels']]
        if group.get('box'):
            images = [figs.crop_local(pv, group['box']) for _, pv in loaded_v]
        else:
            images = figs.aligned([pv for _, pv in loaded_v], 1.0)
        cells = ''
        for (v, pv), im in zip(loaded_v, images):
            name = save(im, f'{prefix}{vi}--{v["render"]}.webp')
            cells += (f'<figure><div class="frame"><img class="zoom" style="--w:{im.width}px" src="img/{esc(name)}" width="{im.width}" '
                      f'height="{im.height}" loading="lazy" alt="{esc(v["label"])} {esc(v["name"])}"></div>'
                      f'<figcaption><b>{esc(v["label"])}</b>{esc(v["name"])}<span>{esc(figs.summary_ja(pv))}</span></figcaption></figure>')
        html += (f'<section class="variant" id="{prefix}{vi}"><h3>{esc(group["title"])}</h3><p>{group.get("note", "")}</p>'
                 f'<div class="crops{"" if group.get("box") else " whole"}">{cells}</div></section>')
    return html


variants = picture_groups(cfg.get('variants', []), 'variant')
docs = picture_groups(cfg.get('docs', []), 'doc')

# LP-D11 (Owner 2026-10-12 on H2): its figures live on a second page (one version holds at most 511 files)
d11 = ''
if cfg.get('d11_url'):
    d11 = (f'<section class="wide" id="d11"><h2>{esc(cfg["d11_title"])}</h2><p class="hint">{cfg["d11_intro"]} '
           f'<a href="{esc(cfg["d11_url"])}">図と結果は別のページにあります。</a></p></section>')

design_list = ''.join(
    f'<div><dt><b>{esc(d["label"])}</b>{esc(d["name"])}</dt><dd>{d["desc"]}</dd></div>' for d in cfg['designs'])
compare = ''
if cfg.get('compare'):
    head = ''.join(f'<th scope="col">{esc(c)}</th>' for c in cfg['compare']['columns'])
    body = ''.join(
        f'<tr{" class=pick" if r.get("pick") else ""}><th scope="row"><b>{esc(designs[r["design"]]["label"])}</b>'
        f'{esc(designs[r["design"]]["name"])}</th>' + ''.join(f'<td>{c}</td>' for c in r['cells']) + '</tr>'
        for r in cfg['compare']['rows'])
    compare = (f'<section class="wide"><h2>候補の比較</h2><div class="scroll"><table class="compare"><thead><tr>'
               f'<th scope="col">配置</th>{head}</tr></thead><tbody>{body}</tbody></table></div></section>')
tools = ''
if cfg.get('tools'):
    body = ''.join(f'<tr><th scope="row">{esc(t[0])}</th><td>{t[1]}</td></tr>' for t in cfg['tools'])
    tools = (f'<section class="wide"><h2>ほかのツールのやり方</h2><div class="scroll"><table class="compare">'
             f'<tbody>{body}</tbody></table></div></section>')
li = lambda items: ''.join(f'<li>{x}</li>' for x in items)
facts = ''.join(f'<div><dt>{esc(k)}</dt><dd>{v}</dd></div>' for k, v in cfg.get('facts', []))
toc = ''.join(f'<a href="#{esc(c["id"])}">{esc(c.get("short", c["title"]))}</a>' for c in cfg['cases'])
sheet_toc = ''.join(f'<a href="#sheet-{esc(k)}">{esc(d["label"])}</a>' for k, d in designs.items())

page = f'''<title>{esc(cfg["title"])}</title>
<link rel="preconnect" href="https://fonts.googleapis.com"><link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans+JP:wght@400;500;700&display=swap">
<style>
/* Layout: one reading column for the brief and the decision, then full-width figures: the same part of one map
   under every design, one canvas per record whose placement is switched in place, and one sheet per design.
   Tokens follow the S4 pairs pages; figures always sit on white paper because they are drawn for it. */
:root {{
  --bg: #f6f7f9; --panel: #ffffff; --fg: #1d2330; --muted: #5b6475; --line: #d9dee7;
  --accent: #1f5fbf; --accent-soft: #e6eefb; --ok: #0f7a4a; --bad: #b3261e; --paper: #ffffff;
  --sans: "IBM Plex Sans JP", "Hiragino Sans", "Noto Sans JP", system-ui, sans-serif;
  --mono: "IBM Plex Mono", ui-monospace, "SFMono-Regular", Menlo, monospace;
}}
@media (prefers-color-scheme: dark) {{ :root:not([data-theme="light"]) {{
  --bg: #12151b; --panel: #1a1f27; --fg: #e5e9f0; --muted: #a2abbb; --line: #2c3340;
  --accent: #7fb0ff; --accent-soft: #1d2b44; --ok: #5fd39a; --bad: #ff8a80; --paper: #ffffff; color-scheme: dark; }} }}
:root[data-theme="dark"] {{
  --bg: #12151b; --panel: #1a1f27; --fg: #e5e9f0; --muted: #a2abbb; --line: #2c3340;
  --accent: #7fb0ff; --accent-soft: #1d2b44; --ok: #5fd39a; --bad: #ff8a80; --paper: #ffffff; color-scheme: dark; }}
body {{ background: var(--bg); color: var(--fg); font-family: var(--sans); line-height: 1.7; }}
.wrap {{ max-width: 1760px; margin: 0 auto; padding-inline: 20px; padding-block: 28px 64px; display: grid; gap: 32px; }}
.brief {{ max-width: 74ch; display: grid; gap: 18px; }}
.wide {{ display: grid; gap: 10px; min-width: 0; }}
h1 {{ font-size: 1.6rem; line-height: 1.3; margin: 0; text-wrap: balance; }}
h2 {{ font-size: 1.05rem; margin: 0 0 6px; letter-spacing: .02em; }}
h3 {{ font-size: 1.15rem; margin: 0; line-height: 1.4; }}
.eyebrow {{ font-family: var(--mono); font-size: .78rem; color: var(--muted); letter-spacing: .06em; text-transform: uppercase; }}
p, ul, ol, dl {{ margin: 0; }} ul, ol {{ padding-left: 1.3em; display: grid; gap: 6px; }}
a {{ color: var(--accent); }}
code {{ font-family: var(--mono); font-size: .9em; }}
.ask {{ background: var(--accent-soft); border: 1px solid var(--accent); border-radius: 8px; padding: 14px 18px; }}
.ask ol {{ gap: 10px; }}
dl.designs {{ display: grid; gap: 8px; }}
dl.designs dt {{ font-weight: 500; }}
dl.designs dd {{ margin: 0; color: var(--muted); font-size: .92rem; }}
dt b, figcaption b, .notes b, table.compare th[scope=row] b {{ font-family: var(--mono); font-weight: 500; color: var(--accent); margin-right: 8px; }}
dl.facts {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(220px, 1fr)); gap: 10px 24px; }}
dl.facts div {{ min-width: 0; }} dl.facts dt {{ font-size: .78rem; color: var(--muted); }}
dl.facts dd {{ margin: 0; font-family: var(--mono); font-size: .85rem; overflow-wrap: anywhere; }}
.scroll {{ overflow-x: auto; }}
table {{ border-collapse: collapse; font-size: .88rem; }}
th, td {{ text-align: left; vertical-align: top; padding: 6px 12px; border-bottom: 1px solid var(--line); }}
thead th {{ font-size: .78rem; font-weight: 500; color: var(--muted); white-space: nowrap; }}
table.compare {{ min-width: 1280px; }} table.compare th[scope=row] {{ min-width: 13em; font-weight: 500; }}
table.compare thead th {{ white-space: normal; min-width: 5.5em; }} table.compare td:nth-last-child(-n+2) {{ min-width: 13em; }}
table.compare tr.pick {{ background: var(--accent-soft); }}
table.metrics td {{ font-family: var(--mono); font-variant-numeric: tabular-nums; text-align: right; white-space: nowrap; }}
table.metrics thead th:not(:first-child) {{ text-align: right; }}
table.metrics td.ok {{ color: var(--ok); }} table.metrics td.bad {{ color: var(--bad); font-weight: 500; }}
table.metrics tr.on {{ background: var(--accent-soft); }}
table.metrics th button {{ font: inherit; font-family: var(--mono); color: var(--accent); background: none; border: 0; padding: 0; cursor: pointer; text-decoration: underline; }}
nav.toc {{ display: flex; flex-wrap: wrap; gap: 6px; }}
nav.toc a {{ font-size: .82rem; text-decoration: none; border: 1px solid var(--line); border-radius: 4px; padding: 2px 8px; background: var(--panel); }}
.cases, .sheets {{ display: grid; gap: 56px; }}
.case {{ display: grid; gap: 12px; min-width: 0; padding-top: 16px; border-top: 2px solid var(--line); }}
.case header {{ display: grid; gap: 4px; max-width: 90ch; }}
.meta {{ font-family: var(--mono); font-size: .8rem; color: var(--muted); }}
.bar {{ display: flex; flex-wrap: wrap; gap: 8px 16px; align-items: center; position: sticky; top: env(safe-area-inset-top, 0px); z-index: 2; background: var(--bg); padding-block: 6px; }}
.tabs {{ display: flex; flex-wrap: wrap; gap: 4px; }}
.tabs button {{ font: inherit; display: grid; text-align: left; line-height: 1.25; padding: 5px 12px; border: 1px solid var(--line); border-radius: 6px; background: var(--panel); color: var(--fg); cursor: pointer; }}
.tabs button b {{ font-family: var(--mono); font-weight: 500; }} .tabs button span {{ font-size: .74rem; color: var(--muted); }}
.tabs button[aria-checked=true] {{ border-color: var(--accent); background: var(--accent-soft); box-shadow: inset 0 -3px 0 var(--accent); }}
.hold {{ font: inherit; font-size: .82rem; padding: 5px 12px; border: 1px dashed var(--muted); border-radius: 6px; background: transparent; color: var(--muted); cursor: pointer; touch-action: none; user-select: none; }}
.hold.down {{ background: var(--accent-soft); color: var(--fg); }}
button:focus-visible, a:focus-visible, img:focus-visible, .stack:focus-visible {{ outline: 2px solid var(--accent); outline-offset: 2px; }}
.frame {{ overflow-x: auto; background: var(--paper); border: 1px solid var(--line); border-radius: 6px; padding: 6px; }}
.stack {{ position: relative; width: min(100%, var(--w)); cursor: zoom-in; }}
.stack.actual {{ width: var(--w); cursor: zoom-out; }}
.stack img {{ position: absolute; display: block; max-width: none; height: auto; }}
[hidden] {{ display: none !important; }}
img.zoom {{ display: block; width: min(100%, var(--w)); height: auto; cursor: zoom-in; }}
img.zoom.actual {{ width: var(--w); max-width: none; cursor: zoom-out; }}
.notes {{ max-width: 90ch; font-size: .92rem; min-height: 1.7em; }}
.same {{ max-width: 90ch; font-size: .85rem; color: var(--muted); margin: 0; }}
.rec {{ display: block; margin-top: 2px; color: var(--muted); font-size: .92rem; }}
table.compare .ok {{ color: var(--ok); }} table.compare .bad {{ color: var(--bad); }}
.cropset {{ display: grid; gap: 8px; min-width: 0; }} .cropset > p {{ font-size: .92rem; max-width: 90ch; }}
.crops {{ display: grid; grid-template-columns: repeat(auto-fill, minmax(min(100%, 340px), 1fr)); gap: 14px; }}
.crops figure, figure.sheet {{ margin: 0; display: grid; gap: 4px; min-width: 0; align-content: start; }}
.crops .frame {{ padding: 0; }} .crops img {{ display: block; width: 100%; height: auto; }}
.crops figcaption {{ font-size: .85rem; }} .crops figcaption span {{ display: block; color: var(--muted); font-size: .8rem; }}
.crops.whole {{ grid-template-columns: repeat(auto-fill, minmax(min(100%, 520px), 1fr)); }}
.crops img.zoom {{ width: 100%; }} .crops img.zoom.actual {{ width: var(--w); max-width: none; }}
.crops .frame {{ overflow-x: auto; }}
.variant {{ display: grid; gap: 10px; min-width: 0; padding-top: 16px; border-top: 2px solid var(--line); }} .variant > p {{ max-width: 90ch; font-size: .92rem; }}
figure.sheet figcaption {{ font-weight: 500; }} figure.sheet figcaption span {{ display: block; font-weight: 400; font-size: .9rem; color: var(--muted); max-width: 90ch; }}
.hint {{ font-size: .85rem; color: var(--muted); max-width: 90ch; }}
footer {{ font-size: .85rem; color: var(--muted); max-width: 74ch; }}
ul.plain {{ max-width: 90ch; font-size: .92rem; color: var(--muted); margin: 6px 0 0; padding-left: 1.2em; }}
</style>
<div class="wrap">
<header class="brief">
  <div class="eyebrow">{esc(cfg["eyebrow"])}</div>
  <h1>{esc(cfg["heading"])}</h1>
  <p>{cfg["intro"]}</p>
</header>
{f'<section class="wide"><h2>{esc(cfg["hero"]["title"])}</h2>{hero}</section>' if hero else ''}
<div class="brief">
  {f'<section><h2>推奨と進め方</h2><ul>{li(cfg["recommendation"])}</ul></section>' if cfg.get("recommendation") else ''}
  <section class="ask"><h2>ご判断をお願いしたいこと</h2><ol>{li(cfg["questions"])}</ol></section>
  <section><h2>調査でわかったこと</h2><ul>{li(cfg["findings"])}</ul></section>
  <section><h2>候補</h2><dl class="designs">{design_list}</dl></section>
  {f'<dl class="facts">{facts}</dl>' if facts else ''}
</div>
{d11}
{f'<section class="wide" id="docs"><h2>{esc(cfg["docs_title"])}</h2><p class="hint">{cfg["docs_intro"]}</p><ul class="plain">{li(cfg["docs_unchanged"])}</ul></section><div class="sheets">{docs}</div>' if docs else ''}
{compare}
<section class="wide"><h2>ケースごとに見る</h2>
  <p class="hint">各ケースは同じレコード・同じ文字サイズです。ボタンか ← → キーで配置を切り替えると、地図の中心を合わせたまま入れ替わります。図をクリックすると原寸表示になります。その下は、同じ部分を原寸で切り出して全候補を並べたものです。</p>
  <nav class="toc" aria-label="ケース一覧">{toc}</nav></section>
<main class="cases">
{chr(10).join(case_html)}
</main>
{f'<section class="wide"><h2>候補ごとに見る</h2><p class="hint">1 つの候補を全ケースに適用した一覧です。クリックで原寸表示になります。</p><nav class="toc" aria-label="候補一覧">{sheet_toc}</nav></section><div class="sheets">{sheets}</div>' if sheets else ''}
{f'<section class="wide"><h2>設定で変えられるところ</h2><p class="hint">同じレコードで、1 つの候補の設定だけを変えた図です。クリックで原寸表示になります。</p></section><div class="sheets">{variants}</div>' if variants else ''}
{tools}
<footer>{cfg.get("footer", "")}</footer>
</div>
<script>
document.querySelectorAll('.case').forEach((el) => {{
  const tabs = [...el.querySelectorAll('.tabs button')];
  const stack = el.querySelector('.stack');
  let current = el.dataset.first;
  const show = (k) => {{
    el.querySelectorAll('.stack img, .notes p').forEach((n) => {{ n.hidden = n.dataset.k !== k; }});
    el.querySelectorAll('table.metrics tbody tr').forEach((r) => r.classList.toggle('on', r.dataset.k === k));
  }};
  const select = (k, focus) => {{
    current = k;
    tabs.forEach((t) => {{
      const on = t.dataset.k === k;
      t.setAttribute('aria-checked', on ? 'true' : 'false');
      t.tabIndex = on ? 0 : -1;
      if (on && focus) t.focus();
    }});
    show(k);
  }};
  const warm = () => el.querySelectorAll('.stack img').forEach((img) => {{ img.loading = 'eager'; }});
  el.addEventListener('pointerenter', warm, {{ once: true }});
  el.addEventListener('focusin', warm, {{ once: true }});
  tabs.forEach((t, i) => {{
    t.addEventListener('click', () => select(t.dataset.k, false));
    t.addEventListener('keydown', (e) => {{
      const step = (e.key === 'ArrowRight' || e.key === 'ArrowDown') ? 1 : (e.key === 'ArrowLeft' || e.key === 'ArrowUp') ? -1 : 0;
      if (!step) return;
      e.preventDefault();
      select(tabs[(i + step + tabs.length) % tabs.length].dataset.k, true);
    }});
  }});
  el.querySelectorAll('table.metrics button').forEach((b) => b.addEventListener('click', () => select(b.dataset.k, false)));
  const hold = el.querySelector('.hold');
  const base = tabs[0].dataset.k;
  const press = (down) => {{ hold.classList.toggle('down', down); show(down ? base : current); }};
  hold.addEventListener('pointerdown', (e) => {{ hold.setPointerCapture(e.pointerId); press(true); }});
  ['pointerup', 'pointercancel', 'blur'].forEach((n) => hold.addEventListener(n, () => press(false)));
  hold.addEventListener('keydown', (e) => {{ if (e.key === ' ' || e.key === 'Enter') {{ e.preventDefault(); press(true); }} }});
  hold.addEventListener('keyup', () => press(false));
  stack.addEventListener('click', () => stack.classList.toggle('actual'));
}});
document.querySelectorAll('img.zoom').forEach((img) => {{
  img.tabIndex = 0;
  const toggle = () => img.classList.toggle('actual');
  img.addEventListener('click', toggle);
  img.addEventListener('keydown', (e) => {{ if (e.key === 'Enter' || e.key === ' ') {{ e.preventDefault(); toggle(); }} }});
}});
</script>
'''
(out / 'index.html').write_text(page)
files = list((out / 'img').iterdir())
print(out / 'index.html', len(cfg['cases']), 'cases', len(files), 'files', round(sum(f.stat().st_size for f in files) / 1e6, 1), 'MB')
