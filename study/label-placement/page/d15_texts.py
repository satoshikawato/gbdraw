"""Texts and picture groups for the Owner's answers of 2026-10-12 (LP-D15). Called by fill_texts.py after d11_texts.

Three new parts in front of the older ones, each with its own image letter: `end` (e; where a leader ends, prototype
J3), `back` (f; examples for the `back` measure), `blocks` (g; examples for prototype L). `end` is added only when
the J3 renders exist (work/J4/*__J_end.*).
"""
import json
from pathlib import Path

WORK = Path(__file__).resolve().parents[1] / 'work'
T01 = 'doc_tcli01_human_mitochondrion'
NAMES = {
    'p1_f14': 'P1 14 pt', 'p1_f8': 'P1 8 pt', 'p1_open_f14': 'P1 14 pt（目盛りなし）', 'lambda_f14': 'λ 14 pt', 'lambda_open_f14': 'λ 14 pt（目盛りなし）',
    'n4_f14': 'N4 14 pt', 'epifagus_f14': 'Epifagus の色素体 14 pt', 'tobacco_f14': 'タバコ葉緑体 14 pt', 'tobacco_f8': 'タバコ葉緑体 8 pt',
    'mg1655_rna_f8': '大腸菌 MG1655 の RNA 8 pt', 'plasmid153_f8': 'プラスミド 153 kb 8 pt', 'plasmid153_f14': 'プラスミド 153 kb 14 pt',
    't5_f8': 'T5 8 pt', 't5_f14': 'T5 14 pt', 'gal_hmmt_atskew': 'ギャラリーの mtDNA（AT skew、20 pt）', 'gal_hmmt_basic': 'ギャラリーの mtDNA',
    'gal_tobacco_chloroplast_hz': 'ギャラリーの葉緑体', 'mito_f14': 'mtDNA 14 pt', 'mito_open_f14': 'mtDNA 14 pt（目盛りなし）', 'pbr322_f14': 'pBR322 14 pt',
    'pbr322_f10': 'pBR322 10 pt', 'pbr322_open_f14': 'pBR322 14 pt（目盛りなし）', 'plasmid37_f14': 'プラスミド 37 kb 14 pt',
    'plasmid37_open_f14': 'プラスミド 37 kb 14 pt（目盛りなし）', T01: 'チュートリアルのヒト mtDNA',
    'doc_tcli09_precomputed_circular_rings': 'チュートリアルの mtDNA（リング付き）', 'doc_tcli03_mitochondrial_features_highlighted': 'チュートリアルの mtDNA（Highlight）',
}
SAME = {'gal_hmmtdna_basic_circular': 'gal_hmmt_basic', 'gal_hmmtdna_atskew': 'gal_hmmt_atskew'}      # the same map under two ids


def _load(sub, case, design):
    return json.loads((WORK / sub / f'{case}__{design}.metrics.json').read_text())


def _cases(sub, design):
    return sorted(p.name.split('__')[0] for p in (WORK / sub).glob(f'*__{design}.metrics.json'))


def _back_over(labels_json, px=1.0):
    """Labels of one render whose box comes more than `px` nearer to the ring than the end of their leader."""
    import sys
    sys.path.insert(0, str(WORK.parent / 'scripts'))
    import lp_metrics
    labs = [l for l in json.loads(Path(labels_json).read_text())['labels'] if not l.get('is_embedded') and l.get('aabb_local') and l.get('leader_points')]
    back = lp_metrics.back_px(labs, [tuple(l['aabb_local']) for l in labs], [[tuple(q) for q in l['leader_points']] for l in labs], [False] * len(labs))
    return sum(1 for value, _over, _text in back if value > px)


def back_part(V):
    m = lambda sub, case, design: _load(sub, case, design)
    h = lambda sub, case, design: f"{m(sub, case, design)['canvas_px'][1]:,}"
    ln = lambda sub, case, design: round(m(sub, case, design)['leader_len_sum'])
    over1 = sum(_back_over(WORK / 'J' / f'{c}__J_sym2.labels.json') for c in _cases('J', 'J_sym2'))
    p0, p1, p2 = ('back', 'p1_f14', 'J_sym2_b0'), ('J', 'p1_f14', 'J_sym2'), ('back', 'p1_f14', 'J_sym2_b17')
    part = {'key': 'back', 'letter': 'f', 'title': 'ご質問 2 の例: 「ラベルがリング側へ戻る」のはどういう場合か',
            'points': [
        '<b>平たく言うと</b>: ラベルの枠のどこかが、そのリーダー線の終点よりもリング（円）に近い所にあることです。'
        '下の図では、線の終点を通る円（リングと同じ中心）を<b style="color:#e27000">橙の破線</b>で描きました。枠が破線より内側（リング側）に入っている分が「戻り」です。'
        '橙の枠は戻りが 1 px を超えるラベル、橙の点は線の終点です。',
        '<b>どういう場合に起きるか</b>: 文字が、線の終点から 12 時（または 6 時）に近い側へ向かって伸びるラベルで起きます。文字は水平で、リングは円なので、'
        '円の曲がりの分だけ、文字の一部が終点よりリングに近くなります。'
        '12 時・6 時のすぐ近くでは数 px です（J2 で最大 8.7 px、H2・H3 でも最大 13 px ありました。1、2 組目の図）。'
        '12 時から離れた所（11 時や 1 時）で長い文字を 12 時の方向へ伸ばすと、数十 px になります（G の Phd-like antitoxin は 24.5 px、T5 14 pt では最大 68 px、今の配置は最大 94 px。3、4 組目の図）。'
        '最初にご指摘いただいたのは、この大きいほうです。',
        '<b>お聞きしているのは、小さいほうをどこまで許すかです。</b>大きいほうはどの案でも禁止のままです。'
        f'小さいほうも禁止する（上限 0 px）と、12 時・6 時の近くでも文字は終点から外側へしか伸ばせないので、ラベルが 1 段ずつ高く積まれます。'
        f'P1 14 pt の高さは {h(*p1)} → {h(*p0)} px、線の合計は {ln(*p1):,} → {ln(*p0):,} px（{(ln(*p0) / ln(*p1) - 1) * 100:.0f} % 増）。'
        f'N4 14 pt は {h("J", "n4_f14", "J_sym2")} → {h("back", "n4_f14", "J_sym2_b0")} px、λ 14 pt は {h("J", "lambda_f14", "J_sym2")} → {h("back", "lambda_f14", "J_sym2_b0")} px（5 組目の図）。',
        f'<b>上限は「枠の高さの半分」で足ります。</b>J2 の既定がこれで、14 pt なら約 8 px です。上限を枠の高さ（約 17 px）まで広げても、P1 14 pt の高さは {h(*p1)} → {h(*p2)} px、'
        f'線の合計は {ln(*p1):,} → {ln(*p2):,} px で、ほとんど変わりません。前回のページでは「枠の高さまで許す」とお聞きしましたが、半分で同じ結果が得られるので、推奨を「枠の高さの半分まで」に改めます。'
        f'J2 で戻りが 1 px を超えるラベルは 29 図の 2,267 個のうち {over1} 個です。',
    ]}
    groups = [
        {'title': 'P1 14 pt の 12 時付近を 3 倍に: 戻りを禁止、枠の高さの半分まで（J2）、枠の高さまで', 'case': 'p1_f14', 'box': [-150, -604, 180, -494], 'zoom': 3, 'mark': ('back',),
         'note': 'Lpa-like transcriptional activator と terminase 2 個です。戻りを禁止すると（左）、3 個とも 12 時の線より右に置かれ、Lpa-like transcriptional activator は 1 段上に積まれます。'
                 'J2（中）は、文字の一部が破線の円より最大 8.6 px 内側（リング側）に入る代わりに、3 個が低く収まります。右は上限を枠の高さにした場合で、最大 16.7 px 入ります。',
         'panels': [V('J_sym2_b0', '戻りを禁止', '上限 0 px'), V('J_sym2', 'J2', '上限は枠の高さの半分'), V('J_sym2_b17', '枠の高さまで', '上限 17 px')]},
        {'title': 'λ 14 pt の 12 時付近を 3 倍に: 戻りを禁止した場合と J2', 'case': 'lambda_f14', 'box': [-285, -528, 175, -458], 'zoom': 3, 'mark': ('back',),
         'note': 'Rz-like spanin（11 時の feature）は、J2 では feature のすぐ上に置かれ、文字が 12 時の方向へ伸びます。文字の右端は終点より 8.6 px リングに近くなります。'
                 '戻りを禁止すると、このラベルは左の列へ回り、線が長くなります。',
         'panels': [V('J_sym2_b0', '戻りを禁止', '上限 0 px'), V('J_sym2', 'J2', '上限は枠の高さの半分')]},
        {'title': '大きく戻る例（どの案でも禁止）: P1 14 pt の 11 時付近を 3 倍に（試作 L の 2 つの設定）', 'case': 'p1_f14', 'box': [-320, -520, -40, -340], 'zoom': 3, 'mark': ('back',),
         'note': '左は、試作 L で「文字がリング側へ伸びるのも許す」設定にした場合です。UmuD-like protein は線の終点から 12 時の方向へ伸び、文字の右端が終点より 29.3 px リングに近くなります（Doc-like toxin は 19.4 px）。'
                 '右は推奨の設定で、同じ 2 個は線から外向き（左）に伸び、戻りは最大 8.6 px に収まります。',
         'panels': [V('L', 'L（文字がリング側へ伸びるのも許す）', '枠の高さを超える戻りが 2 個'), V('L_f0o1', 'L（推奨の設定）', '枠の高さを超える戻りは 0 個')]},
        {'title': '最初にご指摘いただいた形: P1 14 pt の 11 時付近を 3 倍に（G と J2）', 'case': 'p1_f14', 'box': [-330, -452, -20, -332], 'zoom': 3, 'mark': ('back',),
         'note': 'G の Phd-like antitoxin は線の終点から 12 時の方向へ伸び、文字の右端が終点より 24.5 px リングに近くなります（Doc-like toxin は 21.0 px）。'
                 'どちらも枠の高さ（14 px）を超えるので、新しい測り方でも「戻っている」に数えます。J2 は、この場所のラベルをすべて線から外向きに伸ばします。',
         'panels': [V('G', 'G', '枠の高さを超える戻りが 2 個'), V('J_sym2', 'J2', '0 個')]},
        {'title': '戻りを禁止した場合の全体: N4 14 pt（ラベル 69 個）', 'case': 'n4_f14',
         'note': f"高さは戻りを禁止すると {h('back', 'n4_f14', 'J_sym2_b0')} px、J2 は {h('J', 'n4_f14', 'J_sym2')} px です。12 時の塔が高くなります。",
         'panels': [V('J_sym2_b0', '戻りを禁止', '上限 0 px'), V('J_sym2', 'J2', '上限は枠の高さの半分')]},
    ]
    return part, groups


def blocks_part(V):
    cases = [c for c in _cases('L3', 'L_f0o1') if c not in SAME]
    rows = []
    for c in cases:
        m = _load('L3', c, 'L_f0o1'); blocks = (m['proto']['notes'].get('blocks') or {})
        top, bottom = (blocks.get('top') or {}).get('moved', 0), (blocks.get('bottom') or {}).get('moved', 0)
        rows.append((c, top, bottom, _load('J', c, 'J_sym2')['canvas_px'][1], m['canvas_px'][1], (m['proto']['notes'].get('behind') or {}).get('leaders', 0), m['labels_external'], m['seconds']))
    made = [r for r in rows if r[1] or r[2]]
    none = [r for r in rows if not (r[1] or r[2])]
    where = lambda r: '・'.join(x for x in (f'上 {r[1]} 個' if r[1] else '', f'下 {r[2]} 個' if r[2] else '') if x)
    made_txt = '、'.join(f'{NAMES[r[0]]}（{where(r)}）' for r in sorted(made, key=lambda r: -(r[1] + r[2])))
    none_txt = '、'.join(NAMES[r[0]] for r in none)
    table = ('<table class="cmp"><thead><tr><th>図</th><th>ラベル</th><th>行に組んだラベル</th><th>高さ J2 → L (px)</th><th>背面を通る線</th></tr></thead><tbody>'
             + ''.join(f'<tr><td>{NAMES[r[0]]}</td><td>{r[6]}</td><td>{where(r)}</td><td>{r[3]:,} → {r[4]:,}</td><td>{r[5]} 本</td></tr>' for r in sorted(made, key=lambda r: r[4] - r[3]))
             + '</tbody></table>')
    slow = sorted(rows, key=lambda r: -r[7])[:3]
    part = {'key': 'blocks', 'letter': 'g', 'title': 'ご質問 3 の例: 12 時・6 時を行で組む L は、どの図で何が変わるか',
            'points': [
        '<b>平たく言うと</b>: 12 時・6 時の近くに feature が多い図では、H2 も J2 も 1 行にラベルを左右 1 個ずつしか置けないので、ラベルが上下に塔のように積まれます。'
        'L はそこを「1 行に 2〜4 個」の行に組み直します。ラベルは 1 個も省きません。その代わり、奥の行のラベルの線が、手前の行のラベルの枠の背面を通ります。'
        '印を付けた拡大図では、背面を通る線を<b style="color:#b026c4">紫</b>で描き直し、通られる枠を紫で囲みました。',
        f'<b>行に組まれる図（29 図のうち {len(made) + sum(1 for k, v in SAME.items() if any(r[0] == v for r in made))} 図。ギャラリーの mtDNA は同じ図が 2 つずつあるので、下の一覧では 1 つにまとめています）</b>: {made_txt}。ラベルが 50 個前後より多く、12 時か 6 時の近くに集まっている図です。',
        f'<b>行に組まれない図（{len(none) + sum(1 for k, v in SAME.items() if any(r[0] == v for r in none))} 図）</b>: {none_txt}。'
        'まばらな図は塔ができないので組みません。T5 14 pt とプラスミド 153 kb 14 pt は逆に混みすぎていて、行に組んでも高さが縮まないので、組みません（高さを決めているのが 12 時・6 時ではなく左右の列です）。'
        'これらの図では L は H2 に「線の終点は枠の辺のどこでも」を足しただけの配置になります（最後の図）。',
        '<b>行に組まれる図の変化</b>（J2 との比較。L の左右の列は H2 のままなので、J2 の対称な扇形はまだ入っていません）:<div style="overflow-x:auto;max-width:calc(100vw - 48px)">' + table + '</div>',
        '<b>見た目の課題</b>（前回と同じ）: 同じ行で隣り合うラベルが 3 px しか離れないので、ひと続きの言葉に見える所があります。3 本以上の線が 1 つのラベルの背面を通る所は線を追いにくくなります'
        '（P1 14 pt の Lpa-like transcriptional activator）。T5 8 pt のように 80 本以上が背面を通る図は、線を追うのが難しくなります。',
        '<b>試作の状態</b>: このページの L は、探索の時間切れを外して描き直しました。同じ入力で 2 回描いて、混んだ 5 図ともバイト単位で同じ結果です（再現できます）。'
        f'時間は長く、{NAMES[slow[0][0]]} で {slow[0][7]:.0f} 秒、{NAMES[slow[1][0]]} で {slow[1][7]:.0f} 秒、{NAMES[slow[2][0]]} で {slow[2][7]:.0f} 秒かかります（H2 は 2〜3 秒）。製品にするには探索の作り直しが要ります。',
    ]}
    J2, L = V('J_sym2', 'J2', '1 行に左右 1 個ずつ'), V('L_f0o1', 'L', '12 時・6 時を行に組む')
    r = {x[0]: x for x in rows}
    hh = lambda c: f'高さ {r[c][3]:,} → {r[c][4]:,} px、背面を通る線 {r[c][5]} 本。'
    Lm = V('L_f0o1', 'L', '背面を通る線を紫で')
    groups = [
        {'title': 'P1 14 pt の全体（ラベル 111 個）: J2 と L', 'case': 'p1_f14', 'panels': [J2, L],
         'note': hh('p1_f14') + f' 上の {r["p1_f14"][1]} 個、下の {r["p1_f14"][2]} 個を行に組んでいます。'},
        {'title': 'P1 14 pt の 12 時付近を 2 倍に', 'case': 'p1_f14', 'panels': [J2, L], 'box': [-560, -620, 560, -300], 'zoom': 2,
         'note': 'J2 では 20 個あまりが 1 行に 1〜2 個ずつ積まれます。L は同じラベルを 3〜4 個の行に並べます。'},
        {'title': 'その中央を 3 倍に: 描かれたままの L', 'case': 'p1_f14', 'panels': [L], 'box': [-215, -478, 215, -318], 'zoom': 3,
         'note': '線がラベルの背面を通る所です。枠は白の半透明で塗ってあり、背面の線は薄く見えます。'},
        {'title': '同じ場所に印を付けたもの', 'case': 'p1_f14', 'panels': [Lm], 'box': [-215, -478, 215, -318], 'zoom': 3, 'mark': ('behind',),
         'note': '紫の線が、紫で囲んだ枠の背面を通っています（印は図の上に重ねて描いたので、ここでは線が文字の手前に見えます）。'
                 'Lpa-like transcriptional activator の背面を 6 本、ImcB を 2 本、terminase を 2 本が通ります。'},
        {'title': 'P1 14 pt の 6 時付近を 2 倍に', 'case': 'p1_f14', 'panels': [J2, L], 'box': [-560, 300, 560, 620], 'zoom': 2,
         'note': '下側の 14 個を行に組んでいます。'},
        {'title': '同じ 6 時付近の L に印を付けたもの', 'case': 'p1_f14', 'panels': [Lm], 'box': [-560, 300, 560, 620], 'zoom': 2, 'mark': ('behind',),
         'note': '紫の線が、紫で囲んだ枠の背面を通っています。'},
        {'title': 'λ 14 pt の全体（ラベル 71 個）: J2 と L', 'case': 'lambda_f14', 'panels': [J2, L],
         'note': hh('lambda_f14') + ' 行に組まれるのは上だけです。'},
        {'title': 'λ 14 pt の 12 時付近を 2 倍に', 'case': 'lambda_f14', 'panels': [J2, L], 'box': [-360, -540, 360, -300], 'zoom': 2,
         'note': '背面を通る線は 5 本です。'},
        {'title': '同じ 12 時付近の L に印を付けたもの', 'case': 'lambda_f14', 'panels': [Lm], 'box': [-360, -540, 360, -300], 'zoom': 2, 'mark': ('behind',),
         'note': '紫の線が、紫で囲んだ枠の背面を通っています。'},
        {'title': 'N4 14 pt の全体（ラベル 69 個）: J2 と L', 'case': 'n4_f14', 'panels': [J2, L],
         'note': hh('n4_f14') + ' 12 時の塔がなくなります。'},
        {'title': '行に組まれない図: pBR322 14 pt の全体（H2、J2、L）', 'case': 'pbr322_f14',
         'note': f"12 時・6 時にラベルが集まらないので、L は行を作りません。高さは H2 {_load('I', 'pbr322_f14', 'I')['canvas_px'][1]} px、J2 {r['pbr322_f14'][3]} px、L {r['pbr322_f14'][4]} px。"
                 'H2 との違いは、線の終点が枠の辺のどこでもよい規則によるものだけです。',
         'panels': [V('I', 'H2', ''), J2, L]},
    ]
    return part, groups


END_KEYS = ('leaders_end_far_half', 'leaders_end_off_middle', 'leader_ends_unclear', 'leader_end_gap_ink_over5')


def _lop_mod():
    """The J2 agent's lopsided-cluster count (work/J/lop.py; definition in logs/design-J2.md)."""
    import sys
    sys.path.insert(0, str(WORK / 'J'))
    import lop
    return lop


def _end_sum(sub, design, cases):
    """Totals of the leader-end measures and of the costs over `cases` for one render set."""
    import sys
    sys.path.insert(0, str(WORK.parent / 'scripts'))
    import fan_order
    tot = dict.fromkeys(END_KEYS + ('leaders_straight', 'leaders_drawn', 'leader_len_sum', 'height', 'lopsided', 'clusters', 'mixed'), 0)
    tot['gap_max'], tot['clean'] = 0.0, []
    for case in cases:
        folder = sub if (WORK / sub / f'{case}__{design}.metrics.json').exists() else 'docs'
        m = _load(folder, case, design)
        for key in END_KEYS + ('leaders_straight', 'leaders_drawn', 'leader_len_sum'):
            tot[key] += m.get(key) or 0
        tot['height'] += m['canvas_px'][1]
        tot['gap_max'] = max(tot['gap_max'], m.get('leader_end_gap_ink_max') or 0.0)
        labels = str(WORK / folder / f'{case}__{design}.labels.json')
        tot['lopsided'] += sum(1 for c in _lop_mod().clusters(labels) if c['lop'])
        n, bad = fan_order.count(labels)
        tot['clusters'] += n; tot['mixed'] += len(bad)
        if not bad:
            tot['clean'].append(case)
    return tot


def end_part(V):
    cases = _cases('J', 'J_sym2')
    have = lambda design: all((WORK / 'J4' / f'{c}__{design}.metrics.json').exists() for c in cases)
    after = 'J_fan' if have('J_fan') else 'J_mid' if have('J_mid') else 'J_end'
    fan = after == 'J_fan'
    h2, h3, j2, j3 = _end_sum('I', 'I', cases), _end_sum('J', 'J', cases), _end_sum('J', 'J_sym2', cases), _end_sum('J4', after, cases)
    n = j2['leaders_drawn']
    row = lambda name, t: (f"<tr><td>{name}</td><td>{t['leaders_end_far_half']}</td><td>{t['leaders_end_off_middle']}</td><td>{t['leader_ends_unclear']}</td>"
                           f"<td>{t['leader_end_gap_ink_over5']:,}</td><td>{t['gap_max']:.1f}</td><td>{t['mixed']} / {t['clusters']}</td><td>{t['leaders_straight']:,}</td>"
                           f"<td>{round(t['leader_len_sum']):,}</td></tr>")
    table = ('<table class="cmp"><thead><tr><th>配置</th><th>反対側の半分で終わる線</th><th>真横から来て、端の真ん中で終わらない線</th><th>どのラベルの線か紛らわしい終点</th>'
             '<th>終点から文字まで 5 px 超</th><th>その最大 (px)</th><th>折れ線が直線に挟まれた塊 / 塊</th><th>直線</th><th>線の合計 (px)</th></tr></thead><tbody>'
             + row('H2', h2) + row('H3', h3) + row('J2', j2) + row('<b>J3</b>', j3) + '</tbody></table>')
    left = []
    for case in cases:
        m = _load('J4', case, after)
        if any(m.get(k) for k in END_KEYS if k != 'leaders_end_off_middle'):
            left.append(f"{NAMES.get(SAME.get(case, case), case)}（{m['leaders_end_far_half']} / {m['leader_ends_unclear']} / {m['leader_end_gap_ink_over5']}）")
    clean = len(cases) - len(left)
    t2, t3 = _load('J', T01, 'J_sym2'), _load('J4', T01, after)
    part = {'key': 'end', 'letter': 'e', 'title': 'ご指摘（線の終点、塊の中の並び）への対応: 試作 J3',
            'points': [
        '<b>ご指摘</b>: (a) 上の feature から来る線は、枠の上の頂点から辺の真ん中までの間で終わるべきで、下の頂点で終わるのはおかしい（tRNA-Trp、tRNA-Asn）。tRNA-Cys の線は終点が枠から遠い。'
        '(b) tRNA-Ala を右へずらせば tRNA-Asn の線は直線でよい。(c) tRNA-Trp は枠を少し下げて、枠の左端の真ん中から線を出す。'
        '(d) 塊の一番外が折れ線で、真ん中は直線になるのが望ましい。外側が直線で、中に折れ線が混ざり、また直線が出てくるのはよくない。',
        f'<b>まず数えました</b>（29 図の線 {n:,} 本。下の表）。線が来る向きから見て枠の辺の<b style="color:#d61e1e">反対側の半分</b>で終わる線は、J2 に {j2["leaders_end_far_half"]} 本ありました。'
        f'終点が自分の文字よりほかのラベルの文字に近い（<b style="color:#b026c4">どのラベルの線か紛らわしい</b>）線は {j2["leader_ends_unclear"]} 本、'
        f'終点から自分の文字まで 5 px を超える線は {j2["leader_end_gap_ink_over5"]:,} 本（最大 {j2["gap_max"]:.1f} px）です。'
        f'折れ線が直線に挟まれている塊は {j2["clusters"]} 個のうち {j2["mixed"]} 個でした。1 つの図だけの問題ではありませんでした。',
        '<b>J3 の規則</b>: (1) 線は、枠の辺のうち線が来る側の半分で終わる。(2) 終点と文字の間隔を、どのラベルでも同じにする（約 3〜4.5 px）。'
        '(3) ほぼ真横から枠の左右の端に着く線は、端の真ん中で終わる。'
        + ('(4) 折れ線のある塊は「直線の扇形」に組み直す: どの線も feature からまっすぐ伸ばし、枠の左右の端の真ん中で終える。上のラベルの枠が次の線を塞がないよう、ラベルを自分の線に沿って外へずらす。'
           '組み直しは、塊が元の範囲に収まり、直線が増え、直線と折れ線が交互にならない場合だけ採用します。' if fan else
           '(4) 塊の中の並び（直線は真ん中に 1 続き、折れるのは外側だけ）は、まだ試作に入っていません。'),
        f'<b>結果</b>: 反対側の半分で終わる線は {j2["leaders_end_far_half"]} → {j3["leaders_end_far_half"]} 本、紛らわしい終点は {j2["leader_ends_unclear"]} → {j3["leader_ends_unclear"]} 本、'
        f'文字まで 5 px を超える線は {j2["leader_end_gap_ink_over5"]:,} → {j3["leader_end_gap_ink_over5"]} 本、'
        f'<b style="color:#146edc">真横から来て端の真ん中で終わらない線</b>は {j2["leaders_end_off_middle"]} → {j3["leaders_end_off_middle"]} 本になりました。29 図のうち {clean} 図は最初の 3 つとも 0 です。'
        f'直線は {j2["leaders_straight"]:,} → {j3["leaders_straight"]:,} 本です。ご指摘の tRNA 5 個は 5 本とも直線で平行になり、どれも枠の左端の真ん中で終わります（1 組目の図の右端）。'
        '重なりなどの必須の条件は 29 図すべてで満たしています。',
        f'<b>(d) はまだ一部しかできていません。</b>折れ線が直線に挟まれている塊は、29 図の合計では {j2["mixed"]} → {j3["mixed"]} 個（{j3["clusters"]} 個中）で、'
        + ('減っていません。' if j3['mixed'] >= j2['mixed'] else '少し減っただけです。')
        + f'この数が 0 の図は {len(j2["clean"])} → {len(j3["clean"])} 図です（ヒト mtDNA の各図、pBR322、プラスミド 37 kb など、ラベルがまばらな図）。'
        'λ、Epifagus、葉緑体と、混んだ図では残ります。feature どうしが近い塊（線の出口が約 11 px より近い）では、上のラベルの枠が次の直線を塞ぐので、直線の扇形にできません。'
        'そこでは「真ん中を直線、外側へ行くほど大きく折る」扇形が要りますが、これはまだ作っていません（次の作業です）。',
        f'<b>代わりに失うもの</b>: 線の合計は J2 より {100 * (j3["leader_len_sum"] / j2["leader_len_sum"] - 1):.1f} % 長くなります。'
        '直線の扇形では、ラベルの列がリングから 5〜23 px 離れます（ヒト mtDNA の tRNA 5 個で線が +65 px）。'
        f'左右が対称でない塊は {j2["lopsided"]} → {j3["lopsided"]} 個です。図の高さの合計は {j2["height"]:,} → {j3["height"]:,} px。'
        '計算時間は、混んだ図で 20〜65 秒かかります（J2 は約 10 秒）。製品に入れるときは探索に上限を設けます。',
        (f'<b>終点の規則の残り</b>（反対側 / 紛らわしい / 5 px 超）: {"、".join(left)}。' if left else '<b>終点の規則の残り</b>: 3 つの数は 29 図すべてで 0 です。')
        + '混んだ図の 12 時・6 時で起きます。線の向きが 70 度までに制限された扇形の中では、規則どおりに置ける場所が残らないラベルがあり、そのラベルは規則をゆるめて置いています（最後の組の図）。',
        '前回お見せした「線に沿ってラベルをリング側へ寄せる」（shift）は、J3 の上では線の合計が 1 % 短くなるだけで、終点の数は変わりません。既定にはしません。',
    ]}
    AFTER = V(after, 'J3', '直線の扇形まで' if fan else '終点の規則')
    marks = ('far', 'offmid', 'unclear')
    legend = ('<b style="color:#d61e1e">赤</b>は反対側の半分で終わる線、<b style="color:#146edc">青</b>は真横から来て端の真ん中で終わらない線、'
              '<b style="color:#b026c4">紫</b>の丸はどのラベルの線か紛らわしい終点です。')
    pair = [V('J_sym2', 'J2', '前回の版'), AFTER]
    steps = ([V('J_sym2', 'J2', '前回の版'), V('J_end', '途中の版 1', '規則 1、2 のみ（(a)(b) まで）'), V('J_mid', '途中の版 2', '規則 3 も（(c) まで）'), AFTER]
             if fan else pair)
    groups = [
        {'title': 'ご指摘の箇所: ヒト mtDNA のチュートリアル図の tRNA 5 個を 3 倍に（印なし）', 'case': T01, 'box': [290, 140, 480, 310], 'zoom': 3,
         'note': '上から tRNA-Trp、tRNA-Ala、tRNA-Asn、tRNA-Cys、tRNA-Tyr です。途中の版 1 は、外側の 2 本が折れて真ん中の 3 本が直線ですが、tRNA-Trp の線が枠の左下の頂点で終わります（(c) のご指摘）。'
                 '途中の版 2 は tRNA-Trp を直しましたが、直線と折れ線が交互になりました（(d) のご指摘）。右端の J3 は 5 本とも直線で、枠の左端の真ん中で終わります。'
                 'ラベルの列は途中の版より 5〜23 px 右へ寄ります。',
         'panels': steps},
        {'title': '同じ場所に印を付けたもの', 'case': T01, 'box': [290, 140, 480, 310], 'zoom': 3, 'mark': marks, 'note': legend, 'panels': steps},
        {'title': 'mtDNA 14 pt の 4 時付近を 3 倍に', 'case': 'mito_f14', 'box': [240, 150, 440, 300], 'zoom': 3, 'mark': marks, 'note': legend, 'panels': pair},
        {'title': 'pBR322 14 pt の 4 時付近を 3 倍に', 'case': 'pbr322_f14', 'box': [225, 245, 375, 335], 'zoom': 3, 'mark': marks, 'note': legend, 'panels': pair},
        {'title': 'λ 14 pt の 9 時付近を 3 倍に（直線と折れ線が交互に残る所）', 'case': 'lambda_f14', 'box': [-613, -58, -383, 112], 'zoom': 3, 'mark': marks,
         'note': legend + ' この塊は feature どうしが近く、直線の扇形にできませんでした。', 'panels': pair},
        {'title': 'ギャラリーの葉緑体の 3 時付近を 3 倍に', 'case': 'gal_tobacco_chloroplast_hz', 'box': [300, -117, 530, 53], 'zoom': 3, 'mark': marks, 'note': legend, 'panels': pair},
        {'title': 'ヒト mtDNA のチュートリアル図の全体と、29 図の合計', 'case': T01, 'note': '表の数は少ないほど良く、直線は多いほど良い数です。', 'extra': table, 'panels': pair},
        {'title': '残っている所: P1 14 pt の 12 時付近を 2 倍に', 'case': 'p1_f14', 'box': [-230, -620, 230, -450], 'zoom': 2, 'mark': marks,
         'note': legend + ' 12 時の扇形の中では、規則どおりに置けないラベルが残ります。', 'panels': pair},
    ]
    return part, groups


def span_part(V):
    part = {'key': 'span', 'letter': 'h', 'title': 'トラックをまたぐ underlay（ご決定: SPAN-C）',
            'points': [
        '<b>ご決定</b>（2026-10-12、「SPAN-Cでいきましょう」）: アノテーションのトラック（overlay）に、もう 1 つのトラックを <code>to_slot=&lt;トラックの id&gt;</code> で指定すると、'
        '帯が 2 つのトラックの内側の端から外側の端まで広がります。Circular と Linear で同じ意味です。',
        '<b>進め方</b>: まず描画エンジン（CLI、Python API、Session ファイル）。GUI は、読み込んだ Session の指定をそのまま保って描きます。'
        'GUI で編集する欄は、開いている型付けの PR（#1022、#1023）が dev に入った後に足します。0.14.0 のタグに間に合えば 0.14.0、間に合わなければ 0.14.x です。',
        '<b>いまの状態</b>: エンジンの PR を作っています。開く前に別のエージェントが差分を確認します。'
        '目標の見た目は、下の「葉緑体の図の Region annotation」の節の「くさび + 輪」です。',
        '<b>ギャラリーとチュートリアルの葉緑体の図</b>は、新しいラベル配置が決まった時に一度で差し替えます（二度作り直さないための段取りで、こちらで決めました）。',
    ]}
    return part, []


QUESTIONS = [
    '<b>線の終点と塊の並びの規則（J3）はこの方向でよいですか</b>（先頭の節）。線は枠の辺のうち線が来る側の半分で終わる、終点と文字の間隔はどのラベルでも同じ、'
    'ほぼ真横から来る線は枠の端の真ん中で終わる、折れ線のある塊は直線の扇形に組み直す、の 4 つです。'
    '<span class="rec">推奨: はい。代わりに、線の合計が J2 より約 7 % 長くなり、直線の扇形ではラベルの列がリングから 5〜23 px 離れます。'
    'feature どうしが近い塊では直線の扇形にできないので、そこは「真ん中を直線、外側ほど大きく折る」扇形を次に作ります。'
    'tRNA 5 個のように全部を直線にできる塊でも、外側を折った扇形のほうがよい場合はお知らせください。</span>',
    '<b>ほぼ真横から来る線の規則は、12 時・6 時（枠の上下の辺に着く線）にも同じように当てはめますか</b>。tRNA-Trp へのご指摘は左右の端についてでした。'
    '<span class="rec">推奨: 上下の辺では「辺の真ん中」ではなく、いまの「線が来る側の半分」のままにします。12 時・6 時では文字が横に長く、真ん中に合わせると塔が高くなるためです。</span>',
    '<b>「ラベルがリング側へ戻る」のは、枠の高さの半分までなら許してよいですか</b>（2 番目の節の図。12 時・6 時の近くで、文字の一部が線の終点より数 px リングに近くなる場合です）。'
    '<span class="rec">推奨: はい（14 pt で約 8 px まで）。禁止すると P1 14 pt は 74 px、N4 14 pt は 214 px 高くなります。枠の高さまで広げても得るものはありません。'
    '11 時や 1 時で文字が 12 時の方向へ数十 px 伸びる形は、どの案でも禁止のままです。</span>',
    '<b>12 時・6 時を「1 行に複数のラベル」で組む L を取り入れてよいですか</b>（3 番目の節の図）。行に組むと、内側の行のラベルの線が、外側の行のラベルの背面を通ります。'
    '<span class="rec">推奨: はい。行に組まれるのは混んだ図だけで（29 図のうち 14 図）、P1 14 pt は 86 px、N4 14 pt は 91 px、T5 8 pt は 152 px 低くなります。'
    'ラベルが少ない図は J2 と同じです。背面を通る線の描き方（下地のみか、下線つきか。「ご指摘 3」の節の図）は図を見てお決めください。</span>',
    f'<b>元のページの 7 つの質問</b>（<a href="https://claude.ai/artifact/2PH4kKKkmqTF4RQv1nLeyG">元のページ</a>の「ご判断をお願いしたいこと」）もまだ伺っていません: '
    '入りきらないときに何を譲るか、リングの内側を使うか、0.14.x に「今の見た目のまま重ねない」を先に入れるか、追加の表現、Auto の文字サイズ、'
    '「hypothetical protein」を既定で隠すか、です（1 つ目の「H2 でよいか」は、上の質問に置き換わりました）。',
]


def apply(cfg, V):
    new_parts = []
    builds = (end_part, back_part, blocks_part, span_part) if (WORK / 'J4' / f'{T01}__J_end.metrics.json').exists() else (back_part, blocks_part, span_part)
    for build in builds:
        part, groups = build(V)
        new_parts.append(part); cfg[f"d11_{part['key']}"] = groups
    cfg['d11_parts'] = new_parts + cfg['d11_parts']
    cfg['d11_questions'] = QUESTIONS
    cfg['d11_title'] = 'H2 へのご指摘と、その後のご指摘への対応（2026-10-12）'
    cfg['d11_intro'] = ('<b>今回足したもの</b>（先頭の 4 節）: 線の終点の規則を直した試作 J3、「ラベルがリング側へ戻る」の例、12 時・6 時を行で組む L の例、'
                        'トラックをまたぐ underlay の進め方です。その下は前回までの内容（J2、H3、Region annotation、G / K / L）で、変えていません。'
                        'ラベル配置はどれも調査用の試作で、製品のコードは変えていません。')
