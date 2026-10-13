"""Build the second Owner page: the answers to the Owner's three points on H2 (LP-D11, 2026-10-12).

Usage: python make_d11_page.py <page.json> <out_dir> <main_page_index.html>
The style block is taken from the built main page so that both pages share one look; texts and picture groups
come from d11_texts.py (through fill_texts.py -> page.json). Figures go to <out_dir>/img/ (lossless WebP).
"""
import html, json, re, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import figs
import d11_texts

cfg = json.loads(Path(sys.argv[1]).read_text())
out = Path(sys.argv[2]); (out / 'img').mkdir(parents=True, exist_ok=True)
main_html = Path(sys.argv[3]).read_text(encoding='utf-8')
head = main_html[main_html.index('<link rel="preconnect"'):main_html.index('</style>') + len('</style>')]
esc = html.escape
dirs = cfg['render_dirs']
used = set()


def save(im, name):
    im.save(out / 'img' / name, format='WEBP', lossless=True, method=4)
    used.add(name)
    return name


def picture_groups(groups, prefix):
    """One case under several placements or settings: whole maps on one shared canvas, or one map region at 1:1."""
    text = ''
    for vi, group in enumerate(groups, 1):
        loaded_v = [(v, figs.load(figs.find_render(dirs, group['case'], v['render']))) for v in group['panels']]
        if group.get('zoom'):
            images = [figs.crop_zoom(pv, group['box'], group['zoom']) for _, pv in loaded_v]
        elif group.get('box'):
            images = [figs.crop_local(pv, group['box']) for _, pv in loaded_v]
        else:
            images = figs.aligned([pv for _, pv in loaded_v], 1.0)
        if group.get('mark') and group.get('box'):
            for (v, pv), im in zip(loaded_v, images):
                if v['render'] not in group.get('mark_skip', ()):
                    figs.draw_marks(im, pv, group['box'], group.get('zoom', 1), group['mark'])
        cells = ''
        for (v, pv), im in zip(loaded_v, images):
            name = save(im, f'{prefix}{vi}--{v["render"]}.webp')
            cells += (f'<figure><div class="frame"><img class="zoom" style="--w:{im.width}px" src="img/{esc(name)}" width="{im.width}" '
                      f'height="{im.height}" loading="lazy" alt="{esc(v["label"])} {esc(v["name"])}"></div>'
                      f'<figcaption><b>{esc(v["label"])}</b>{esc(v["name"])}<span>{esc(figs.summary_ja(pv))}</span></figcaption></figure>')
        text += (f'<section class="variant" id="{prefix}{vi}"><h3>{esc(group["title"])}</h3><p>{group.get("note", "")}</p>{f'<div style="overflow-x:auto;max-width:calc(100vw - 48px)">{group["extra"]}</div>' if group.get("extra") else ""}'
                 f'<div class="crops{"" if group.get("box") else " whole"}">{cells}</div></section>')
    return text


parts = ''
toc = ''
for part in cfg['d11_parts']:
    letter = part['letter']                       # fixed per part: image names do not change when parts are reordered
    groups = picture_groups(cfg.get(f'd11_{part["key"]}', []), f'd11{letter}')
    points = ''.join(f'<li>{p}</li>' for p in part['points'])
    toc += f'<a href="#d11-{part["key"]}">{esc(part["title"])}</a>'
    parts += (f'<section class="wide" id="d11-{part["key"]}"><h2>{esc(part["title"])}</h2><ul class="plain">{points}</ul></section>'
              + (f'<div class="sheets">{groups}</div>' if groups else ''))
questions = ''.join(f'<li>{q}</li>' for q in cfg.get('d11_questions', []))
page = f'''<title>Label Placement H2 Follow-up</title>
{head}
<div class="wrap">
<header class="brief">
  <div class="eyebrow">gbdraw · Circular のラベル配置 · H2 へのご指摘 2026-10-12</div>
  <h1>{esc(cfg["d11_title"])}</h1>
  <p>{cfg["d11_intro"]} 候補の比較と H2 までの経緯は<a href="{esc(d11_texts.MAIN_URL)}">元のページ</a>にあります。図をクリックすると原寸になります。</p>
  <nav class="toc" aria-label="ご指摘の一覧">{toc}</nav>
</header>
<div class="brief">
  {f'<section class="ask"><h2>ご判断をお願いしたいこと</h2><ol>{questions}</ol></section>' if questions else ''}
</div>
{parts}
<footer>{cfg.get("footer", "")}</footer>
</div>
<script>
document.querySelectorAll('img.zoom').forEach((img) => {{
  img.tabIndex = 0;
  const toggle = () => img.classList.toggle('actual');
  img.addEventListener('click', toggle);
  img.addEventListener('keydown', (e) => {{ if (e.key === 'Enter' || e.key === ' ') {{ e.preventDefault(); toggle(); }} }});
}});
</script>
'''
(out / 'index.html').write_text(page, encoding='utf-8')
for old in (out / 'img').iterdir():
    if old.name not in used:
        old.unlink()
size = sum(p.stat().st_size for p in (out / 'img').iterdir())
print(f'{out}/index.html {len(used)} images {size / 1e6:.1f} MB, html {len(page.encode()) / 1e3:.0f} KB')
