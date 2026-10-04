#!/usr/bin/env python
"""Evidence-only BGC Similarity alignment flows (S00 step 6). Usage:
python bgc_align_flows.py dev:1 dev:2 main:1   (each run = fresh browser context)"""
import json, os, re, sys, time
from pathlib import Path
from playwright.sync_api import sync_playwright

BASE = Path(os.environ.get('S00_BASELINE_DIR', '/home/kawato/gbdraw-baselines/gui-remediation-20260928'))
OUT = BASE / 'evidence/bgc-align'
FIXTURE = 'gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw-session.json'
APPS = {'dev': ('http://127.0.0.1:4302/gbdraw/web/index.html', BASE / 'dev-57cef3ba' / FIXTURE),
        'main': ('http://127.0.0.1:4301/gbdraw/web/index.html', BASE / 'main-4556e04e' / FIXTURE)}
GENES = ['livA', 'livE', 'racM', 'racL']
DIALOG = 'Select alignment anchors'

INIT_JS = r"""(() => {
  const act = { constructions: [], posts: {} };
  window.__BGC_WORKER__ = act; window.__BGC_REJECTIONS__ = [];
  window.addEventListener('unhandledrejection', e => window.__BGC_REJECTIONS__.push(String(e.reason?.message || e.reason)));
  const N = window.Worker, post = N.prototype.postMessage;
  N.prototype.postMessage = function (m, t) {
    if (this.__bgc) { const k = m?.type === 'helper' ? 'helper:' + m.operation : String(m?.type); act.posts[k] = (act.posts[k] || 0) + 1; }
    return t === undefined ? post.call(this, m) : post.call(this, m, t);
  };
  window.Worker = new Proxy(N, { construct(target, args) {
    const w = Reflect.construct(target, args, target), url = String(args[0] || '');
    act.constructions.push(url.split('/').pop().split('?')[0]);
    if (url.includes('diagram-generation-worker')) w.__bgc = true;
    return w;
  } });
})();"""

SNAP_JS = r"""async (genes) => {
  const { state } = await import('./js/state.js');
  let committed = null;
  try { committed = (await import('./js/services/config.js')).getCommittedCanonicalSession?.() || null; } catch (e) {}
  const app = window.__GBDRAW_APP__, h = window.__GBDRAW_HISTORY__;
  const svg = document.querySelector('.gbdraw-preview-surface svg');
  const clone = v => v == null ? null : JSON.parse(JSON.stringify(v));
  const plain = s => String(s ?? '').replace(/<[^>]+>/g, '');
  const text = sel => document.querySelector(sel)?.innerText?.trim() ?? null;
  const arrow = (svgId) => {
    const els = svg && svgId ? [...svg.querySelectorAll('[data-gbdraw-feature-id]')]
      .filter(e => e.getAttribute('data-gbdraw-feature-id') === svgId) : [];
    const pts = [];
    for (const el of els) {
      const d = el.getAttribute('d') || '';
      if (!d || /[^MLZz0-9eE.,\s+-]/.test(d)) continue;
      const m = el.getScreenCTM(), n = (d.match(/-?\d*\.?\d+(?:[eE][-+]?\d+)?/g) || []).map(Number);
      for (let i = 0; i + 1 < n.length; i += 2) { const p = new DOMPoint(n[i], n[i + 1]).matrixTransform(m); pts.push([p.x, p.y]); }
    }
    if (!pts.length) return { arrow: els.length ? 'unparsed' : 'absent', parts: els.length };
    const u = pts.filter((p, i) => pts.findIndex(q => Math.hypot(q[0] - p[0], q[1] - p[1]) < 0.01) === i);
    const xs = u.map(p => p[0]), maxX = Math.max(...xs), minX = Math.min(...xs), tol = 0.5;
    const atMax = u.filter(p => maxX - p[0] <= tol).length, atMin = u.filter(p => p[0] - minX <= tol).length;
    return { arrow: atMax === 1 && atMin >= 2 ? '→' : atMin === 1 && atMax >= 2 ? '←' : 'no-head',
      parts: els.map(e => e.getAttribute('data-gbdraw-feature-part')), screenX: [+minX.toFixed(1), +maxX.toFixed(1)] };
  };
  const ef = state.extractedFeatures?.value || [];
  const geneOf = f => f?.gene || f?.qualifiers?.gene?.[0] || f?.qualifiers?.gene || '';
  const view = f => f ? { gene: geneOf(f), recordKey: f.recordKey ?? f.record_key, biologicalFeatureId: f.biologicalFeatureId,
    svgId: f.svg_id, start: f.start, end: f.end, sourceStrand: f.strand, orthogroupId: f.orthogroupId ?? null, ...arrow(f.svg_id) } : null;
  const groups = svg ? [...svg.querySelectorAll('[data-gbdraw-composition-role="primary"][data-record-key]')]
    .filter(e => !e.hasAttribute('data-gbdraw-definition-part')) : [];
  const seqs = state.linearSeqs || [];
  const records = (committed?.renderRequest?.records || []).map((r, i) => {
    const s = seqs.find(x => String(x.uid) === r.recordKey) || {}, g = groups.find(x => x.getAttribute('data-record-key') === r.recordKey)
      || svg?.querySelector(`g[data-gbdraw-record-index="${i}"][transform]`);
    return { recordKey: r.recordKey, label: plain(r.presentation?.label), reverseComplement: r.presentation?.reverseComplement ?? null,
      region: clone(r.region), linearSeq: { region_start: s.region_start ?? null, region_end: s.region_end ?? null, region_reverse: s.region_reverse ?? null },
      svgGroup: g ? { transform: g.getAttribute('transform'), translationX: g.getAttribute('data-record-translation-x') } : null };
  });
  const plan = clone(state.similarityAlignmentPlan?.value);
  const find = a => ef.find(f => a && f.recordKey === a.recordKey && f.biologicalFeatureId === a.biologicalFeatureId);
  const sha = async s => [...new Uint8Array(await crypto.subtle.digest('SHA-256', new TextEncoder().encode(s)))]
    .map(b => b.toString(16).padStart(2, '0')).join('').slice(0, 16);
  const content = state.results?.value?.[0]?.content;
  return {
    records,
    features: Object.fromEntries(genes.map(g => [g, view(ef.find(f => geneOf(f) === g))])),
    plan: plan ? { groupId: plan.groupId, reference: plan.reference ? { ...plan.reference, gene: geneOf(find(plan.reference)) } : null,
      records: (plan.records || []).map(r => ({ recordKey: r.recordKey, status: r.status, rationale: r.rationale ?? r.reason ?? null,
        anchor: r.anchor ? view(find(r.anchor)) : null })) } : null,
    legacySelectedAlignmentGroup: state.selectedOrthogroupAlignmentFeature?.value ?? null,
    receipt: clone(state.similarityAlignmentResetReceipt?.value) ? 'present' : null,
    result: { count: state.results?.value?.length || 0, contentSha256_16: typeof content === 'string' ? await sha(content) : String(typeof content),
      contentLength: typeof content === 'string' ? content.length : null, previewSvgSha256_16: svg ? await sha(svg.outerHTML) : null },
    history: h ? { undo: h.getUndoCount(), redo: h.getRedoCount(), undoLabel: h.undoLabel?.() ?? null,
      redoLabel: h.redoLabel?.() ?? null, revision: h.revision?.value ?? null } : null,
    ui: { status: app.similarityAlignmentStatus ?? null, dialogOpen: app.similarityAlignmentDialogOpen ?? null,
      summary: text('[data-similarity-alignment-summary]'), statusText: text('[data-similarity-alignment-status]'),
      notice: text('[data-similarity-alignment-notice]'), alignmentError: app.similarityAlignmentError?.summary ?? null,
      errorLog: app.errorLog ? (app.errorLog.summary || app.errorLog.message || String(app.errorLog)) : null,
      generationFeedback: text('[data-generation-application-feedback]') },
    worker: clone(window.__BGC_WORKER__), rejections: [...window.__BGC_REJECTIONS__]
  };
}"""

REVIEW_JS = r"""async () => {
  const { state } = await import('./js/state.js');
  const d = document.querySelector('[data-similarity-alignment-dialog]');
  if (!d) return null;
  const t = el => el ? el.innerText.trim() : null, app = window.__GBDRAW_APP__;
  const gene = a => { const f = (state.extractedFeatures.value || []).find(x => x.recordKey === a?.recordKey && x.biologicalFeatureId === a?.biologicalFeatureId); return f?.gene || null; };
  const draft = app.similarityAlignmentDraft;
  return {
    fullText: t(d), count: t(d.querySelector('[data-similarity-alignment-count]')),
    reference: t(d.querySelector('[data-similarity-alignment-reference]')),
    directionMode: [...d.querySelectorAll('input[name="alignment-direction-mode"]')].find(i => i.checked)?.value || null,
    directionPreview: [...d.querySelectorAll('[data-similarity-alignment-direction]')].map(t),
    applyReason: t(d.querySelector('#similarity-alignment-apply-reason')),
    rows: [...d.querySelectorAll('fieldset[data-alignment-record-key]')].map(fs => {
      const key = fs.getAttribute('data-alignment-record-key'), row = draft?.rows?.find(r => r.recordKey === key);
      return { recordKey: key, legend: t(fs.querySelector('legend')), text: t(fs),
        radios: [...fs.querySelectorAll('input[type=radio]')].map(i => ({ ariaLabel: i.getAttribute('aria-label'), checked: i.checked })),
        draftChoice: row?.choice ?? null, recommendedKey: row?.recommendedKey ?? null,
        candidates: (row?.candidates || []).map(c => ({ key: c.key, label: c.label, gene: gene(c.anchor), featureIdentifier: c.featureIdentifier,
          directEvidence: c.directEvidence, representative: c.representative })) };
    })
  };
}"""

SETTLE_JS = r"""async () => {
  const { state } = await import('./js/state.js');
  const a = window.__GBDRAW_APP__, h = window.__GBDRAW_HISTORY__, c = state.results?.value?.[0]?.content;
  return { busy: Boolean(a.similarityAlignmentBusy || state.processing?.value || state.featureExtractionPending?.value || a.sessionImportPending),
    fp: [state.results?.value?.length, typeof c === 'string' ? c.length : 0, h?.revision?.value, JSON.stringify(window.__BGC_WORKER__.posts)].join('|') };
}"""


def settle(page, quiet=1.5, timeout=300):
    prev, since, t0 = None, 0, time.time()
    while time.time() - t0 < timeout:
        s = page.evaluate(SETTLE_JS)
        if s['busy']:
            prev = None
        elif s['fp'] != prev:
            prev, since = s['fp'], time.time()
        elif time.time() - since >= quiet:
            return round(time.time() - t0, 1)
        page.wait_for_timeout(250)
    raise TimeoutError('app did not settle')


def classify(prev, cur):
    out = []
    for r in cur['records']:
        p = next((x for x in prev['records'] if x['recordKey'] == r['recordKey']), None)
        a, b = bool(p and p['reverseComplement']), bool(r['reverseComplement'])
        out.append({'recordKey': r['recordKey'], 'before': a, 'after': b,
                    'change': 'NEW reversal' if b and not a else 'un-reversed' if a and not b
                    else 'retained reversed' if a else 'retained forward'})
    return out


class Run:
    def __init__(self, browser, app, flow):
        self.app, self.flow, self.tag = app, flow, f'{app}-flow{flow}'
        self.url, self.fixture = APPS[app]
        self.ctx = browser.new_context(viewport={'width': 1440, 'height': 900})
        self.ctx.add_init_script(INIT_JS)
        self.page = self.ctx.new_page()
        self.log = {'app': app, 'flow': flow, 'url': self.url, 'fixture': str(self.fixture), 'steps': [],
                    'pageErrors': [], 'consoleErrors': [], 'deviations': [], 'aids': []}
        self.page.on('pageerror', lambda e: self.log['pageErrors'].append(str(e)))
        self.page.on('console', lambda m: m.type == 'error' and self.log['consoleErrors'].append(m.text))
        self.shot_n = 0

    def shot(self, name, locator=None):
        self.shot_n += 1
        path = OUT / f'{self.tag}-{self.shot_n:02d}-{name}.png'
        (locator or self.page).screenshot(path=str(path))
        return path.name

    def close_popup(self):
        popup = self.page.locator('.feature-popup[role="dialog"]')
        if popup.count() and popup.first.is_visible():
            self.page.get_by_role('button', name='Close feature popup').click()

    def result_shot(self, name):
        svg = self.page.evaluate("() => document.querySelector('.gbdraw-preview-surface svg')?.outerHTML || ''")
        tab = self.ctx.new_page()
        tab.set_content(f'<body style="margin:0;background:#fff">{svg}</body>')
        tab.evaluate("() => { const s = document.querySelector('svg'); const b = s.viewBox.baseVal;"
                     " s.setAttribute('width', '1440'); s.setAttribute('height', String(b && b.width ? 1440 * b.height / b.width : 900)); }")
        self.shot_n += 1
        path = OUT / f'{self.tag}-{self.shot_n:02d}-{name}.png'
        tab.screenshot(path=str(path), full_page=True)
        tab.close()
        return path.name

    def step(self, name, extra=None):
        self.close_popup()
        snap = self.page.evaluate(SNAP_JS, GENES)
        prev = self.log['steps'][-1]['snapshot'] if self.log['steps'] else None
        entry = {'step': name, **(extra or {}), 'snapshot': snap}
        if prev:
            entry['orientationChange'] = classify(prev, snap)
            entry['resultChanged'] = prev['result']['contentSha256_16'] != snap['result']['contentSha256_16']
            entry['previewChanged'] = prev['result']['previewSvgSha256_16'] != snap['result']['previewSvgSha256_16']
            entry['workerPostsDelta'] = {k: v - prev['worker']['posts'].get(k, 0) for k, v in snap['worker']['posts'].items()
                                         if v != prev['worker']['posts'].get(k, 0)}
        entry['screenshots'] = [self.shot(f'{name}-preview'), self.result_shot(f'{name}-result-fit')]
        (OUT / f'{self.tag}-{name}-preview.svg').write_text(self.page.evaluate(
            "() => document.querySelector('.gbdraw-preview-surface svg')?.outerHTML || ''"))
        self.log['steps'].append(entry)
        print(self.tag, name, json.dumps({k: entry.get(k) for k in ('orientationChange', 'resultChanged', 'workerPostsDelta')}), flush=True)

    def load(self):
        self.page.goto(self.url, wait_until='domcontentloaded')
        self.page.wait_for_function('() => window.__GBDRAW_APP__')
        with self.page.expect_file_chooser() as fc:
            self.page.get_by_role('button', name='Load Session').click()
        fc.value.set_files(str(self.fixture))
        self.page.wait_for_function('() => window.__GBDRAW_APP__.results?.length > 0', timeout=300000)
        self.step('00-load', {'settleSeconds': settle(self.page)})

    def open_popup(self, gene):
        popup = self.page.locator('.feature-popup[role="dialog"]')
        self.close_popup()
        svg_id = self.page.evaluate('async g => { const { state } = await import("./js/state.js");'
                                    ' return state.extractedFeatures.value.find(f => f.gene === g)?.svg_id; }', gene)
        feat = self.page.locator(f'.gbdraw-preview-surface svg [data-gbdraw-feature-id="{svg_id}"]').first
        try:
            feat.click(timeout=10000)
        except Exception as error:
            self.log['deviations'].append(f'{gene}: real click failed ({str(error)[:160]}); used dispatchEvent click')
            feat.dispatch_event('click')
        popup.first.wait_for(state='visible', timeout=30000)
        clicked = self.page.evaluate('() => { const c = window.__GBDRAW_APP__.clickedFeature; '
                                     'return { gene: c?.feat?.gene, svgId: c?.feat?.svg_id, orthogroupId: c?.orthogroupId }; }')
        assert clicked['gene'] == gene, clicked
        return popup.first, clicked

    def capture_review(self, label):
        dialog = self.page.get_by_role('dialog', name=DIALOG)
        review = self.page.evaluate(REVIEW_JS)
        review['screenshots'] = [self.shot(f'{label}-review', dialog)]
        for row in review['rows']:
            if row['draftChoice'] is None:
                dialog.locator(f'fieldset[data-alignment-record-key="{row["recordKey"]}"]').scroll_into_view_if_needed()
                review['screenshots'].append(self.shot(f'{label}-review-{row["recordKey"]}', dialog))
        review['displayedInternalIds'] = sorted(set(re.findall(
            r'\bh_[a-z0-9]{20,}\b|\bf[0-9a-f]{8}\b|\brecord-\d+\b|Internal ID: [^\n]+', review['fullText'])))
        return review

    def resolve_unselected(self, label, rows, direction):
        dialog = self.page.get_by_role('dialog', name=DIALOG)
        for row in rows:
            if row['draftChoice'] is not None:
                continue
            idx = next((i for i, c in enumerate(row['candidates']) if c['gene'] == 'racM'), None)
            planned = idx is not None
            if idx is None:
                idx = next((i for i, c in enumerate(row['candidates']) if c['key'] == row['recommendedKey']), 0)
            fs = dialog.locator(f'fieldset[data-alignment-record-key="{row["recordKey"]}"]')
            fs.locator('input[type=radio]').nth(idx).check()
            self.log['aids'].append({'step': label, 'recordKey': row['recordKey'], 'selected': row['candidates'][idx],
                                     'plannedBaselineAid': planned, 'direction': direction})

    def apply(self, label):
        dialog = self.page.get_by_role('dialog', name=DIALOG)
        retries = []
        for _ in range(3):
            dialog.get_by_role('button', name='Apply', exact=True).click()
            self.page.wait_for_function('() => !window.__GBDRAW_APP__.similarityAlignmentBusy', timeout=300000)
            if not dialog.count():
                return retries
            err = self.page.evaluate('() => window.__GBDRAW_APP__.similarityAlignmentError?.summary || ""')
            retries.append(err)
            self.shot(f'{label}-apply-error', dialog)
            if 'Apply again' not in err:
                return retries
        return retries

    def align(self, gene, label, review_direction=None):
        popup, clicked = self.open_popup(gene)
        extra = {'clickedFeature': clicked, 'popupScreenshot': self.shot(f'{label}-popup', popup)}
        if self.app == 'main':
            popup.locator('button[title="Align by this similarity group"]').click()
            extra['control'] = 'legacy popup Align (runAnalysis)'
            extra['settleSeconds'] = settle(self.page)
            return self.step(label, extra)
        name = 'Review alignment options…' if review_direction else 'Align…'
        extra['control'] = f'popup {name}'
        popup.get_by_role('button', name=name, exact=True).click()
        self.page.wait_for_function('() => !window.__GBDRAW_APP__.similarityAlignmentBusy', timeout=300000)
        dialog = self.page.get_by_role('dialog', name=DIALOG)
        extra['reviewOpened'] = dialog.count() > 0
        if extra['reviewOpened']:
            dialog.wait_for(state='visible')
            review = self.capture_review(label)
            extra['review'] = review
            direction = review_direction or 'Keep current directions'
            dialog.get_by_role('radio', name=direction, exact=True).check()
            self.resolve_unselected(label, review['rows'], direction)
            extra['reviewBeforeApply'] = self.capture_review(f'{label}-chosen')
            extra['applyRetries'] = self.apply(label)
        extra['settleSeconds'] = settle(self.page)
        self.step(label, extra)

    def save(self):
        self.log['finalRejections'] = self.page.evaluate('() => window.__BGC_REJECTIONS__')
        (OUT / f'{self.tag}.json').write_text(json.dumps(self.log, indent=1, ensure_ascii=False))
        self.ctx.close()


def main(specs):
    with sync_playwright() as p:
        browser = p.chromium.launch()
        print('chromium', browser.version, flush=True)
        for spec in specs:
            app, flow = spec.split(':')
            run = Run(browser, app, int(flow))
            try:
                run.load()
                run.align('livA', '01-livA', 'All selected arrows right' if flow == '2' else None)
                run.align('livE', '02-livE')
            except Exception as error:
                run.log['fatal'] = repr(error)[:2000]
                run.shot('fatal')
                print(run.tag, 'FATAL', repr(error)[:500], flush=True)
            finally:
                run.save()
        browser.close()


if __name__ == '__main__':
    main(sys.argv[1:] or ['dev:1', 'dev:2', 'main:1'])
