"""Texts and picture groups for the Owner's feedback on H2 (LP-D11, 2026-10-12). Called at the end of fill_texts.py.

Three parts, each its own list so that adding one does not renumber the images of another:
d11_straight (point 1, prototype J), d11_region (point 2, Region annotation variants), d11_over (point 3, prototype K),
d11_sym (LP-D12 to LP-D14, prototype J2; shown first). Each part has a fixed `letter` for its image names.
"""
import json
from pathlib import Path

WORK = Path(__file__).resolve().parents[1] / 'work'
UL, JUNCTION = [-430, -290, -110, 70], [-500, -30, -180, 330]
CHL = 'chl_region'
MAIN_URL = 'https://claude.ai/artifact/2PH4kKKkmqTF4RQv1nLeyG'
D11_URL = 'https://claude.ai/artifact/MkpZPJ1zRXpP3DAByEtqjR'          # the second page (LP-D11)
T01, T03 = 'doc_tcli01_human_mitochondrion', 'doc_tcli03_mitochondrial_features_highlighted'


def _m(cfg, case, design):
    import figs
    stem = figs.find_render(cfg['render_dirs'], case, design)
    assert stem, (case, design)
    return json.loads(Path(str(stem) + '.metrics.json').read_text())


def _st(cfg, case, design):
    m = _m(cfg, case, design)
    return f"直線 {m['leaders_straight']}/{m['leaders_drawn']} 本、曲がりの平均 {m['leader_turn_deg_mean']:.1f}°、線の合計 {round(m['leader_len_sum']):,} px"


def _three(cfg, case):
    return f"今の配置: {_st(cfg, case, 'cur')}。H2: {_st(cfg, case, 'I')}。H3: {_st(cfg, case, 'J')}。"


def _two(cfg, case, a='J', b='J_sym2'):
    return f"H3: {_st(cfg, case, a)}。J2: {_st(cfg, case, b)}。"


def _sum(design, key):
    return sum(json.loads(q.read_text())[key] for q in sorted((WORK / 'J').glob(f'*__{design}.metrics.json')))


def _lop(design, case='*'):
    """Lopsided clusters as the J2 agent's work/J/lop.py counts them (definition in logs/design-J2.md)."""
    import sys
    sys.path.insert(0, str(WORK / 'J'))
    import lop
    return sum(1 for q in sorted((WORK / 'J').glob(f'{case}__{design}.labels.json')) for c in lop.clusters(str(q)) if c['lop'])


def _sum3(key):
    return sum(json.loads(q.read_text())[key] for q in sorted((WORK / 'J3').glob('*__J_sym3.metrics.json')))


def _off_corner(design, tol=3.7):
    """Leaders that end more than `tol` px away from every corner of their label box."""
    import math
    n = 0
    for q in sorted((WORK / 'J').glob(f'*__{design}.labels.json')):
        for lab in json.loads(q.read_text())['labels']:
            box, lead = lab.get('aabb_local'), lab.get('leader_points')
            if lab.get('is_embedded') or not box or not lead:
                continue
            n += min(math.dist(lead[-1], (x, y)) for x in (box[0], box[2]) for y in (box[1], box[3])) > tol
    return n


def _tower(cfg, case, design):
    """How far the labels reach above and below the ring (px beyond the outer edge of the feature band)."""
    import figs
    stem = figs.find_render(cfg['render_dirs'], case, design)
    ext = [lab for lab in json.loads(Path(str(stem) + '.labels.json').read_text())['labels'] if not lab.get('is_embedded') and lab.get('aabb_local')]
    ring = max(float(lab.get('feature_outer_radius_px', 0)) for lab in ext)
    return round(-min(lab['aabb_local'][1] for lab in ext) - ring), round(max(lab['aabb_local'][3] for lab in ext) - ring)


def _lrow(cfg, case, design):
    m = _m(cfg, case, design)
    up, down = _tower(cfg, case, design)
    return f"高さ {m['canvas_px'][1]:,} px（リングの上へ {up} px、下へ {down} px）、線の合計 {round(m['leader_len_sum']):,} px"


def apply(cfg, V):
    for sub in ('region_page', 'J', 'K', 'L3', 'L2', 'J3', 'back', 'J4'):   # L3 before L2: L without wall-clock limits (LP-D15 session)
        d = str(WORK / sub)
        if d not in cfg['render_dirs']:
            cfg['render_dirs'].append(d)
    l2 = [d for d in cfg['render_dirs'] if d.rstrip('/').endswith('/L2')]
    l3 = [d for d in cfg['render_dirs'] if d.rstrip('/').endswith('/L3')]
    if l2 and l3 and cfg['render_dirs'].index(l3[0]) > cfg['render_dirs'].index(l2[0]):   # page.json keeps the old order: the unlimited renders win
        cfg['render_dirs'].remove(l3[0]); cfg['render_dirs'].insert(cfg['render_dirs'].index(l2[0]), l3[0])
    if D11_URL:
        cfg['intro'] += (f' <b>2026-10-12 追記 2</b>: H2 へのご指摘 3 点（線の角度、葉緑体の図の Region annotation、P1 の上下の伸び）への'
                         f'対応は、図が多いので<a href="{D11_URL}">別のページ</a>にまとめました。')
    cfg['d11_url'] = D11_URL
    cfg['d11_title'] = 'H2 へのご指摘 3 点と、その後のご指摘（2026-10-12）'
    cfg['d11_intro'] = ('1. リーダー線の角度も小さくする（直線にできる所は直線に）。2. 葉緑体の図で、内側ラベルと Region annotation の間が空く。'
                        '3. P1 14 pt は H2 でも上下に伸びる（線の上にラベルが被るのを許す案）。'
                        'どれも調査用の試作と描き分けで、製品のコードは変えていません。'
                        ' <b>追記</b>: その後にいただいたご指摘（塊の中の線を対称な扇形に、線の終点は枠の辺のどこでも、'
                        '直線は頂点まで延ばさない）への対応を先頭の節に、12 時・6 時を行で組む試作 L を「ご指摘 3」の節の後半に足しました。')
    st_j, st_2, len_j, len_2 = _sum('J', 'leaders_straight'), _sum('J_sym2', 'leaders_straight'), _sum('J', 'leader_len_sum'), _sum('J_sym2', 'leader_len_sum')
    drawn = _sum('J_sym2', 'leaders_drawn')
    len_1 = _sum('J_sym', 'leader_len_sum')
    back_j = max(json.loads(q.read_text())['label_back_px_max'] for q in (WORK / 'J').glob('*__J.metrics.json'))
    back_2 = max(json.loads(q.read_text())['label_back_px_max'] for q in (WORK / 'J').glob('*__J_sym2.metrics.json'))
    assert _sum('J_sym2', 'labels_back_over_height') == 0 and _sum('J', 'labels_back_over_height') == 0
    hgt = lambda case, design: _m(cfg, case, design)['canvas_px'][1]
    sym_part = {'key': 'sym', 'letter': 'd', 'title': 'その後のご指摘: 対称な扇形、線の終点は枠の辺のどこでも（試作 J2）',
                'points': [
        '<b>J2 は H3 に次の 3 つを足した版です。</b>(1) 塊の中で 1 本だけ傾く並びを、両端が逆向きに傾く扇形にそろえる。(2) リーダー線の終点は、ラベルの枠の辺のどこでもよい。'
        '(3) 直線は頂点まで延ばさず、枠の辺に届いた所で止める。ほかの規則（間隔 3 px、交差なし、省かない、折れは 110° 以上）は H3 のままで、29 の図と設定違い 6 枚のすべてで守れています。',
        f'<b>対称な扇形</b>: 「偏った塊」（上下に接するラベルの並びで、線の曲がる向きが片方に寄っているもの）は 29 図の合計で {_lop("J")} → {_lop("J_sym2")} 個。'
        'ご指摘の tRNA 5 個の並び（ヒト mtDNA の図）と pBR322 の BsgI / BspEI は、両端が逆向きに傾き、間は直線です。',
        f'<b>代償</b>: 直線の線（曲がり 3° 未満）は {drawn:,} 本のうち {st_j:,} → {st_2:,} 本に減ります。扇形では 1 本ではなく 2 本が傾くことと、辺の規則では「短くて少し傾いた線」が「長い直線」より選ばれることが理由です。'
        '曲がりの重み（0.8 のまま）を上げれば直線は戻りますが、線は長くなります。重みはまだ調整していません。',
        f'<b>辺の規則の効果</b>: ラベルが feature のすぐ横や真上に寄るので、線の合計が {round(len_j):,} → {round(len_2):,} px（{(1 - len_2 / len_j) * 100:.1f} % 減）になります。'
        f'キャンバスの高さは T5 8 pt で {hgt("t5_f8", "J"):,} → {hgt("t5_f8", "J_sym2"):,} px、P1 14 pt で {hgt("p1_f14", "J"):,} → {hgt("p1_f14", "J_sym2"):,} px。'
        f'角から 3.7 px より離れた所で終わる線は {drawn:,} 本のうち {_off_corner("J_sym2")} 本です。',
        '<b>「ラベルがリング側へ戻らない」の測り方を変えました。</b>これまでは「線が届く側の端より、反対の端のほうがリングに近いラベル」を数えていましたが、線が辺の途中に届くとこの数え方は当てはまりません。'
        '新しい測り方は「枠のうちリングにいちばん近い点が、線の終点より何 px リング側にあるか」で、これが枠の高さ（14 pt で 12〜17 px）を超えたら「戻っている」と数えます。'
        f'H3 は最大 {back_j} px、J2 は最大 {back_2} px で、どちらも 0 個です（G は 23 図で 28 個、今の配置は 57 個）。古い数え方だと J2 は {_sum("J_sym2", "labels_toward_axis")} 個になります。'
        '下の「角の規則と辺の規則」の図でご確認ください。',
        f'<b>3 つめ（頂点まで延ばさない）の効き方は小さいです。</b>ご指摘の tRNA-Ala は、頂点へ向かう 2.8 px の横棒がなくなり、線が左の辺へまっすぐ届くようになりました。'
        f'ただし枠そのものは動いていません（上の tRNA-Trp との間隔 3 px と、リングとの間隔で止まっています）。全体でも線の合計は {round(len_1):,} → {round(len_2):,} px（{(1 - len_2 / len_1) * 100:.1f} % 減）です。'
        '<b>追加で調べた結果</b>: tRNA-Ala の枠はこれ以上リングへ寄せられません（枠の角がリングとの間隔 3 px にちょうど接しています。線を辺の中ほどに届かせるには枠をリングから離すことになり、線はかえって長くなります）。'
        f'一方、隣の tRNA-Trp と tRNA-Asn は寄せられます。線の向きを変えずに枠をリング側へ寄せる処理を足すと、tRNA-Trp が 12 px、tRNA-Asn が 9 px 寄ります（最初の図の右端）。'
        f'29 図すべてで規則を守れ、線の合計は {round(len_2):,} → {round(_sum3("leader_len_sum")):,} px、直線の本数は変わりません。効き方は小さいですが害がないので、実装の案に入れます。',
        f'<b>効かない所</b>: 混んだ図で 12 時・6 時に積み上がる塊は対称になりません（偏った塊は P1 14 pt で {_lop("J", "p1_f14")} → {_lop("J_sym2", "p1_f14")} 個、T5 8 pt で {_lop("J", "t5_f8")} → {_lop("J_sym2", "t5_f8")} 個）。'
        'ここは 12 時・6 時を「1 行に複数ラベル」で組む試作 L で扱います（ご指摘 3 の節）。',
        '<b>まだ試していないこと</b>: 線の出る位置を feature の中で滑らせること、斜めの直線。',
    ]}
    P1, lm = 'p1_f14', lambda case, design: _m(cfg, case, design)
    hh = lambda case, design: f"{lm(case, design)['canvas_px'][1]:,}"
    hall = lambda case, sub, design: f"{json.loads((WORK / sub / f'{case}__{design}.metrics.json').read_text())['canvas_px'][1]:,}"
    over_all = sum(json.loads(q.read_text())['labels_back_over_height'] for q in (WORK / 'L').glob('*__L.metrics.json'))
    back_all = max(json.loads(q.read_text())['label_back_px_max'] for q in (WORK / 'L').glob('*__L.metrics.json'))
    assert lm(P1, 'L_f0o1')['labels_back_over_height'] == 0
    LPOINTS = [
        f'<b>試作 L（12 時・6 時を「1 行に複数のラベル」で組む）の結果</b>: P1 14 pt の高さは H2 {hh(P1, "I")} px、K {hh(P1, "K")} px、J2 {hh(P1, "J_sym2")} px に対して L は {hh(P1, "L_f0o1")} px です。'
        f'リングの上に出る高さは {_tower(cfg, P1, "I")[0]} → {_tower(cfg, P1, "L_f0o1")[0]} px、下は {_tower(cfg, P1, "I")[1]} → {_tower(cfg, P1, "L_f0o1")[1]} px。'
        f'線の合計は {round(lm(P1, "I")["leader_len_sum"]):,} → {round(lm(P1, "L_f0o1")["leader_len_sum"]):,} px、いちばん長い線は {round(lm(P1, "I")["leader_len_max"])} → {round(lm(P1, "L_f0o1")["leader_len_max"])} px。'
        'ラベルは 1 個も省かず、文字どうし・ほかの物との間隔 3 px、交差なし、折れ 110° 以上も守れています。線の終点は J2 と同じく枠の辺のどこでもよい規則です。',
        f'<b>L が緩める規則は 2 つで、それぞれ切り替えられます。</b>(a) 線がほかのラベルの背面を通る（K と同じ上限: 1 本の線が通るラベルは 3 個まで。P1 14 pt で {lm(P1, "L_f0o1")["leader_through_text"]} か所）。'
        f'(b) 極のまとまりの中では、文字がリング側へ戻る向きにも伸びる。P1 14 pt の高さは、どちらも無しで {hh(P1, "L_f0o0")} px、(b) だけで {hh(P1, "L_f1o0")} px、(a) だけで {hh(P1, "L_f0o1")} px、両方で {hh(P1, "L")} px。'
        '効くのは (a) です。',
        f'<b>推奨は (a) だけを許す設定です</b>（下の図で「L」と書いたもの）。(b) を足しても P1 14 pt で {lm(P1, "L_f0o1")["canvas_px"][1] - lm(P1, "L")["canvas_px"][1]} px しか縮まず、'
        f'「リング側へ戻らない」を破るラベルが出ます（枠の高さを超えて戻るラベルが 29 図の合計で {over_all} 個、最大 {back_all} px）。(a) だけなら 0 個です。',
        f'<b>ほかの図</b>（H2 → L）: N4 14 pt {hh("n4_f14", "I")} → {hh("n4_f14", "L_f0o1")} px、λ 14 pt {hh("lambda_f14", "I")} → {hh("lambda_f14", "L_f0o1")} px、P1 8 pt {hh("p1_f8", "I")} → {hh("p1_f8", "L_f0o1")} px、'
        f'プラスミド 153 kb 8 pt {hh("plasmid153_f8", "I")} → {hh("plasmid153_f8", "L_f0o1")} px、T5 8 pt {hh("t5_f8", "I")} → {hh("t5_f8", "L_f0o1")} px（K は {hh("t5_f8", "K")} px）。'
        'まばらな図は極のまとまりを作らず、H2 に辺の規則を足しただけの配置になります。'
        f'14 pt の混んだ図は L の得意な形ではありません（両方を許す設定で、T5 14 pt は {hall("t5_f14", "I", "I")} → {hall("t5_f14", "L", "L")} px で K の {hall("t5_f14", "K", "K")} px と同程度、'
        f'プラスミド 153 kb 14 pt は {hall("plasmid153_f14", "I", "I")} → {hall("plasmid153_f14", "L", "L")} px で縮みません）。高さを決めているのが極ではなく左右の列だからです。',
        '<b>見た目の課題</b>: 同じ行で隣り合うラベルが 3 px しか離れないので、ひと続きの言葉に見える所があります（「Coi ImcB site-specific integrase」など）。3 本以上の線が 1 つのラベルの背面を通る所（Lpa-like transcriptional activator）は線を追いにくくなります。'
        '下線つきの描き方（下の図）にすると、どの線がどのラベルかは分かりやすくなります。行の中の間隔を広げる設定はまだ試していません。',
        '<b>試作としての限界</b>: (1) L は H2 を土台にしていて、J2 の対称な扇形とはまだ一緒になっていません。(2) 試作は探索を実行時間で打ち切る作りで、混んだ図では機械の負荷によって配置が少し変わっていました。'
        '<b>2026-10-12 夜に、このページの L の図（既定の「線が背面を通るのを許す」設定）を、打ち切りを外して描いたものに差し替えました</b>（同じ入力で 2 回描いて 5 図ともバイト単位で一致。プラスミド 153 kb 8 pt は 1,107 → 1,086 px に変わり、ほかは同じ高さです）。'
        '(3) 打ち切りを外すと大きい図では 50〜130 秒かかります（H2 は 2〜3 秒）。製品にするときは、打ち切りを回数で決めることと、探索の作り直しが要ります。',
    ]
    cfg['d11_parts'] = [
        sym_part,
        {'key': 'straight', 'letter': 'a', 'title': 'ご指摘 1: リーダー線の角度（試作 H3）',
         'points': [
             '<b>H3 は H2 に「線の曲がり」を 2 つめの目的として足した版です。</b>置き場所の良さを「線の長さ + 0.8 px × 曲がりの角度」で比べ、'
             'feature の真上（放射方向の延長）にラベルの端が来る場所を必ず候補に入れます。ほかの規則（間隔 3 px、交差なし、省かない、リング側へ戻らない、折れは 110° 以上）は H2 のままで、29 の図すべてで守れています。',
             '<b>直線の線</b>（曲がり 3° 未満）は 29 の図の合計で 2,267 本のうち H2 552 本 → H3 1,310 本。線の合計は 4.4 % 長くなります（まばらな図は 0〜3 %、最大で 12 %）。キャンバスはほぼ同じです。',
             '<b>まばらな図はほぼ全部が直線になります</b>: ヒト mtDNA のチュートリアル図 24 本中 23 本、pBR322 34 本中 33 本、ギャラリーの mtDNA 32 本中 31 本。'
             '残るのは、tRNA が 5 個並ぶ所の最後の 1 個（12〜14°）と、隣のラベルの文字に阻まれる 1 個です。',
             '<b>混んだ図は直線にしきれません</b>: P1 14 pt は 111 本中 62 本（H2 は 18 本）、T5 8 pt は 194 本中 55 本（H2 は 16 本）。12 時・6 時に積まれるラベルと 3 時・9 時の階段では、順番と間隔を守るために曲がりが残ります。',
             '<b>まだ試していないこと</b>: 斜めの直線（feature の中心からずらして線を出す、feature の中で線の出る位置を滑らせる）。試作の枠組みに feature の両端の位置を渡す変更が要るので、今回は入れていません。入れれば、長い feature のラベルをさらに直線にできます。',
             '<b>重みの選び方</b>: 重みを上げるほど直線が増えますが、0.4 までで大半が直線になり、1 を超えると線が長くなるほうが目立ちます（下の「重みを変えると」の図）。既定は 0.8 にしました。',
         ]},
        {'key': 'region', 'letter': 'b', 'title': 'ご指摘 2: 葉緑体の図の Region annotation（14 通りを描き分け）',
         'points': [
             'ラベルの配置はどの案でも H2 のままです（内側 67 個、外側 52 個、省略 0）。変えたのは Region annotation と GC content の輪だけです。',
             '<b>今の gbdraw で描ける案</b>: 詰める、underlay、OGDRAW 風、2 本のトラックに同じ色を敷く、underlay と名前の細い輪。'
             'Region の表（<code>mark</code> と色の列）とトラックの指定を変えるだけです。underlay は既存の <code>mark=highlight</code> で、feature 帯の背面に敷かれます。',
             '<b>新しいコードが要る案</b>: トラックをまたぐ 1 枚の扇形（下の該当図は SVG を加工した見本です）。'
             '今の underlay は 1 本のトラックの幅にしか敷けません。「このトラックからこのトラックまで」という幅の指定を足せば描けます（描画の順序とラベルの配置は変えずに済みます）。',
             '<b>推奨</b>: underlay と名前の細い輪。空白がなくなり、LSC・IRb・SSC・IRa の名前も残り、既存の機能だけで描けます。'
             '名前を凡例だけにしてよければ underlay 単独（濃さ 40 %）がいちばん静かです。',
             '「詰める」だけの案は、輪をいちばん深い内側ラベル（rRNA の名前）に合わせることになり、ほかの場所では 20〜45 px の空白が残ります。'
             'これ以上外へ出すと内側ラベルが 8〜10 個、外側へ押し出されます。',
             '採用する場合は、ギャラリー「Quadripartite structure of a chloroplast genome」とチュートリアル「Make an annotated map of the tobacco chloroplast genome」の'
             '表・コマンド・図を差し替えます。',
         ]},
        {'key': 'over', 'letter': 'c', 'title': 'ご指摘 3: 12 時・6 時にラベルが積み上がる件（G との違い、試作 K、試作 L）',
         'points': [
             '<b>G が 12 時・6 時をうまく処理して見える理由</b>（P1 14 pt で数えました）: (1) 置けないラベルを省いています。上の 35° 以内にある 24 個のうち 8 個、下の 23 個のうち 7 個、全体で 111 個のうち 26 個。'
             '(2) ほかの物との間隔を守っていません。feature・リングの線・目盛りに重なるラベルが 21 個、3 px 未満が 29 個。'
             '(3) リングの内側の空き（目盛りの間）にも置き、文字が線の左右どちらへも伸びます（リング側へ戻るラベルが 5 個）。'
             'その結果、1 行に複数のラベルが並び、リングより上に出るラベルは 6 個で済んでいます（H2 は 23 個、H3 は 20 個）。',
             '<b>H2・H3 が積み上がる理由</b>: 省かず、3 px を守り、文字は線から外向きにだけ伸ばすので、12 時・6 時では 1 行に左右 1 個ずつしか置けません。'
             '111 個を約 19 px 間隔で 2 列に並べると最低でも約 1,050 px の高さが要ります（リングは約 830 px、H2 は 1,230 px、H3 は 1,171 px）。探し方を良くしても、この前提では縮みません。',
             '<b>G から取り入れられること</b>: 1 行に複数のラベルを置くことです。そのためには、(a) 線がほかのラベルの背面を通るのを許す（試作 K）、'
             '(b) 12 時・6 時の近くでは文字が線の左右どちらへ伸びてもよいことにする、のどちらか、または両方が要ります。省略と間隔の違反は取り入れません。',
             '<b>試作 K（線がラベルの背面を通るのを許す）の結果</b>: 規則は 50 枚すべてで守れましたが、縮み方は小さいです。P1 14 pt の高さは 1,230 px → 1,165 px（強い設定で 1,139 px）、'
             'T5 8 pt は 1,678 → 1,490 px。K は「列に 1 行 1 個」の並べ方を残したまま、列の端のラベルだけを空きへ動かす作りで、線は交差できないので動かせるのは極に近い数個だけでした。'
             '線の合計は短くなります（P1 14 pt で 8,975 → 7,939 px）。まばらな 21 図は H2 と同じ配置のままです。',
             '<b>K の見た目</b>（下の図）: 下地つきの下線がいちばん読めます（下線が線の終点とラベルを結びます）。濃さ 100 % の下地は線がそこで止まって見えます。'
             'feature の色の文字は、薄い青が白地で読みにくく、ご懸念のとおりでした。線が 3 個以上のラベルの背面を通る所（T5 8 pt）は、どの線がどのラベルか追いにくくなります。',
             LPOINTS[0], LPOINTS[1], LPOINTS[2], LPOINTS[3], LPOINTS[4], LPOINTS[5],
         ]},
    ]
    real, edit = '今の gbdraw で描けます', 'SVG を加工した見本（新しいコードが要ります）'
    cfg['d11_region'] = [
        {'title': '全体図で 5 案', 'case': CHL,
         'note': '左から、今の図（輪は <code>r=0.65</code>、いちばん近いラベルまで 44 px）、詰める、underlay、underlay と名前の細い輪、トラックをまたぐ underlay と輪。'
                 '最後の 1 枚だけ SVG を加工した見本です。図をクリックすると原寸になります。',
         'panels': [V('base', '今の図', 'ブラケットの輪（r=0.65）'), V('tight_gc_follow', '詰める', '輪を r=0.75 へ、GC content も外へ'),
                    V('underlay_o40', 'underlay', 'feature 帯の背面に色、名前は凡例'), V('underlay_ring', 'underlay + 輪', '背面の色と、名前を載せた細い輪'),
                    V('wedge_ring', 'またぐ underlay + 輪', edit)]},
        {'title': '詰める: 輪を内側ラベルのすぐ内側へ（左上を原寸で）', 'case': CHL, 'box': UL,
         'note': '輪を <code>r=0.65, w=20px</code> から <code>r=0.75, w=12px</code> に変えると、いちばん深いラベルまでの間が 44 px から 6 px になります。'
                 'GC content をそのままにすると、空白が輪の内側へ移るだけです。GC content も外へ出して幅を広げる（<code>r=0.62, w=0.10</code>）と詰まって見えます。'
                 + real + '。',
         'panels': [V('base', '今の図', 'r=0.65'), V('tight_gc_stay', '輪だけ外へ', 'GC content はそのまま'),
                    V('tight_gc_follow', '輪と GC を外へ', 'GC の幅 0.08 → 0.10'), V('tight_gc_follow_thin', '輪と GC を外へ', 'GC の幅は 0.08 のまま')]},
        {'title': 'underlay: 領域を feature 帯の背面に敷く（左上を原寸で）', 'case': CHL, 'box': UL,
         'note': 'Region の表の <code>mark</code> を <code>bracket</code> から <code>highlight</code> に変え、単一コピー領域と IR を 2 色で塗ります。濃さ 25 %、40 %、55 %、IR にハッチ。'
                 'feature の色や輪郭、ラベルはどの濃さでも読めます。空いた所へ GC content を広げています（<code>r=0.68, w=0.12</code>）。'
                 '領域の名前は図から消えて凡例だけになります（highlight に名前を出す設定は、名前が feature 帯の外側に出て外側ラベルと重なるので使えませんでした）。' + real + '。',
         'panels': [V('underlay_o25', '25 %', '静かだが IR の違いが薄い'), V('underlay_o40', '40 %', '釣り合いがよい'),
                    V('underlay_o55', '55 %', 'やや重い'), V('underlay_hatch', 'ハッチ', 'IR を斜線で（色は 1 色）')]},
        {'title': 'OGDRAW 風: IR を太い線で（左上を原寸で）', 'case': CHL, 'box': UL,
         'note': '左は詰めた輪の中で IR だけを太く（9 px、単一コピーは 2 px）し、名前を残した案。右は IR を軸の上の太い線（14 px）として feature の背面に描いた案で、名前はなく凡例だけです。'
                 + real + '。',
         'panels': [V('og_ring', '輪の中で IR を太く', '名前あり。図の中でいちばん重い線になる'), V('og_axis', '軸の上に IR の太線', '名前なし、IR だけを示す')]},
        {'title': 'トラックをまたぐ: 同じ色を 2 本のトラックに敷く、1 枚の扇形にする（左上を原寸で）', 'case': CHL, 'box': UL,
         'note': '1 枚目は同じ領域の色を feature 帯と GC content の帯の両方に敷いたものです（間のラベルの帯は白いまま）。' + real + '。'
                 '2 枚目は GC content の帯から feature 帯までを 1 枚の扇形にしたものです（濃さ 20 %）。これは SVG を加工した見本で、新しいコードが要ります。内側ラベル 67 個は色の上に載りますが、黒い文字は読めます。'
                 '3、4 枚目は名前を載せた細い輪を足した版です。',
         'panels': [V('across_repeat', '2 本に同じ色', real), V('across_wedge', '1 枚の扇形', edit),
                    V('underlay_ring', 'underlay + 輪', real), V('wedge_ring', '扇形 + 輪', edit)]},
        {'title': '同じ案を SSC と IRb の境界で（原寸）', 'case': CHL, 'box': JUNCTION,
         'note': '領域の境界の両側にラベルがある場所です。境界がどう見えるかを比べます。',
         'panels': [V('base', '今の図', 'ブラケットの輪'), V('tight_gc_follow', '詰める', '輪と GC を外へ'), V('underlay_o40', 'underlay', '40 %'),
                    V('og_axis', '軸の上に IR の太線', ''), V('underlay_ring', 'underlay + 輪', real), V('wedge_ring', '扇形 + 輪', edit)]},
    ]
    NOW = [V('cur', 'Current', '今の配置'), V('I', 'H2', '今のページの推奨'), V('J', 'H3', '直線を優先')]
    HJ = [V('I', 'H2', '今のページの推奨'), V('J', 'H3', '直線を優先')]
    TW = [V('J_tw0', '重み 0', 'H2 と同じ'), V('J_tw04', '重み 0.4', ''), V('J', '重み 0.8', 'H3 の既定'), V('J_tw2', '重み 2', ''), V('J_tw5', '重み 5', '')]
    tw_note = lambda case: '。'.join(f"{v['label']}: {_st(cfg, case, v['render'])}" for v in TW) + '。'
    cfg['d11_straight'] = [
        {'title': 'ご指摘の箇所: ヒト mtDNA のチュートリアル図の 12 時付近を原寸で', 'case': T01, 'panels': NOW, 'box': [-330, -490, 300, -300],
         'note': 'tRNA-Pro、tRNA-Thr、tRNA-Phe、tRNA-Glu、tRNA-Val の線が、H3 では 1 本の直線になります。図全体では ' + _three(cfg, T01)},
        {'title': '同じ図の 6 時付近を原寸で', 'case': T01, 'panels': NOW, 'box': [-300, 280, 300, 490],
         'note': 'ATP8、tRNA-Lys、tRNA-Asp、tRNA-Ser の付近です。'},
        {'title': '同じ図の 3 時〜5 時を原寸で', 'case': T01, 'panels': NOW, 'box': [240, -10, 540, 330],
         'note': 'tRNA が 5 個並ぶ所です。H3 でも最後の 1 個（tRNA-Tyr）は 12〜14° 曲がります。'},
        {'title': 'ヒト mtDNA のチュートリアル図の全体', 'case': T01, 'panels': NOW,
         'note': _three(cfg, T01)},
        {'title': 'お手本にされた箇所: チュートリアル「Highlight」の完成図の 6 時付近を原寸で', 'case': T03, 'panels': NOW, 'box': [-260, 330, 160, 480],
         'note': 'ATP8 は H2 でも直線でした。H3 でも同じ場所です。' + _three(cfg, T03)},
        {'title': '葉緑体の図の左上を原寸で（内側ラベルと外側ラベル）', 'case': 'gal_tobacco_chloroplast_hz', 'panels': HJ, 'box': UL,
         'note': f"H2: {_st(cfg, 'gal_tobacco_chloroplast_hz', 'I')}。H3: {_st(cfg, 'gal_tobacco_chloroplast_hz', 'J')}。"},
        {'title': 'λ 14 pt の左側を原寸で（ラベル 71 個）', 'case': 'lambda_f14', 'panels': NOW, 'box': [-700, -330, -200, 330],
         'note': _three(cfg, 'lambda_f14') + ' ラベルは feature の放射方向の延長に付くので、リングから少し離れる所があります。'},
        {'title': 'P1 8 pt の全体（ラベル 107 個）', 'case': 'p1_f8', 'panels': NOW, 'note': _three(cfg, 'p1_f8')},
        {'title': 'P1 14 pt の全体（ラベル 111 個）', 'case': 'p1_f14', 'panels': NOW,
         'note': _three(cfg, 'p1_f14') + ' 上下に積まれるラベルは H3 でも残ります（ご指摘 3 で扱います）。'},
        {'title': 'T5 8 pt の右側を原寸で（ラベル 194 個）', 'case': 't5_f8', 'panels': HJ, 'box': [330, -330, 830, 330],
         'note': f"H2: {_st(cfg, 't5_f8', 'I')}。H3: {_st(cfg, 't5_f8', 'J')}。ここまで混むと、直線にできる線は 3 割ほどです。"},
        {'title': '重みを変えると: ヒト mtDNA のチュートリアル図の 12 時付近', 'case': T01, 'panels': TW, 'box': [-330, -490, 300, -300],
         'note': tw_note(T01)},
        {'title': '重みを変えると: λ 14 pt の左側', 'case': 'lambda_f14', 'panels': TW, 'box': [-700, -330, -200, 330],
         'note': tw_note('lambda_f14')},
        {'title': '重みを変えると: P1 14 pt の全体', 'case': 'p1_f14', 'panels': TW, 'note': tw_note('p1_f14')},
    ]
    HJ2 = [V('J', 'H3', '直線を優先'), V('J_sym2', 'J2', '対称な扇形、辺の規則')]
    STEPS = [V('J', 'H3', '直線を優先'), V('J_sym', 'J2（1 回目）', '対称な扇形、辺の規則'), V('J_sym2', 'J2', '頂点まで延ばさない'),
             V('J_sym3', 'J2 + 寄せる', '枠をリング側へ寄せる')]
    EDGE = [V('J_sym_e0', '角の規則', '線は文字の端の角に届く'), V('J_sym2', '辺の規則（J2）', '線は枠の辺のどこでも')]
    cfg['d11_sym'] = [
        {'title': 'ご指摘の箇所: ヒト mtDNA のチュートリアル図の tRNA 5 個を 3 倍に拡大', 'case': T01, 'panels': STEPS, 'box': [290, 140, 480, 310], 'zoom': 3,
         'note': 'H3 では tRNA-Tyr だけが傾いていました。J2 では tRNA-Trp と tRNA-Tyr が逆向きに 15.9° ずつ傾き、間の 3 本は直線です。'
                 '1 回目と 3 枚目の違いは tRNA-Ala の線の終点です。1 回目は枠の左上の頂点へ向かう短い横棒があり、3 枚目は線が左の辺へそのまま届きます（差は 2〜3 px）。'
                 '右端は、線の向きを変えずに枠をリング側へ寄せた版です。tRNA-Trp と tRNA-Asn が寄り、tRNA-Ala は動きません。'},
        {'title': '同じ図の 3 時〜5 時を原寸で', 'case': T01, 'panels': HJ2, 'box': [240, -10, 540, 330], 'note': _two(cfg, T01)},
        {'title': 'pBR322 14 pt の BsgI / BspEI を 3 倍に拡大', 'case': 'pbr322_f14', 'panels': HJ2, 'box': [225, 245, 375, 335], 'zoom': 3,
         'note': 'H3 では BspEI 1664 だけが 15.7° 傾いていました。J2 では 2 本が左右に 5.3° ずつ開きます。' + _two(cfg, 'pbr322_f14')},
        {'title': 'mtDNA 14 pt（目盛りが外側の図）の tRNA 5 個を 2 倍に拡大', 'case': 'mito_f14', 'panels': HJ2, 'box': [240, 150, 440, 300], 'zoom': 2,
         'note': _two(cfg, 'mito_f14') + ' 扇形の両端が傾くので、いちばん大きい曲がりは 12.0° から 17.5° になります。'},
        {'title': 'λ 14 pt の左側を原寸で（ラベル 71 個）', 'case': 'lambda_f14', 'panels': HJ2, 'box': [-700, -330, -200, 330],
         'note': _two(cfg, 'lambda_f14') + ' 枠の右の辺の途中に届く線が増えます。4 個以上のラベルが同じ向きに傾く並び（偏った塊 2 個）は残っています。'},
        {'title': '葉緑体の図の左上を原寸で（内側ラベルと外側ラベル）', 'case': 'gal_tobacco_chloroplast_hz', 'panels': HJ2, 'box': UL,
         'note': _two(cfg, 'gal_tobacco_chloroplast_hz') + ' 内側ラベルが feature のすぐ横に寄ります。'},
        {'title': '葉緑体の図の SSC と IRb の境界を原寸で', 'case': 'gal_tobacco_chloroplast_hz', 'panels': HJ2, 'box': JUNCTION,
         'note': 'ndhH、ndhA などの内側ラベルの線が短くなります。'},
        {'title': '角の規則と辺の規則: P1 14 pt の 12 時付近を原寸で', 'case': 'p1_f14', 'panels': EDGE, 'box': [-560, -620, 560, -300],
         'note': '辺の規則では、線が枠の下の辺の途中に届くラベルがあります（C1.100 など）。文字は線の左右にまたがりますが、枠はリングとほぼ平行で、リング側へ戻る量は最大 8.6 px（Lpa-like transcriptional activator）です。'
                 '古い数え方ではこの図の 4 個（C1.100、terminase 2 個、Lpa-like transcriptional activator）が「戻っている」に数えられます。'},
        {'title': '角の規則と辺の規則: λ 14 pt の 12 時付近を原寸で', 'case': 'lambda_f14', 'panels': EDGE, 'box': [-330, -530, 330, -300],
         'note': 'この図でリング側へ戻る量がいちばん大きいのは Rz-like spanin の 8.6 px です。'},
        {'title': 'P1 14 pt の全体（ラベル 111 個）', 'case': 'p1_f14', 'panels': HJ2,
         'note': _two(cfg, 'p1_f14') + ' 上下に積まれるラベルは J2 でも残ります（試作 L で扱います）。'},
        {'title': 'P1 8 pt の全体（ラベル 107 個）', 'case': 'p1_f8', 'panels': HJ2, 'note': _two(cfg, 'p1_f8')},
        {'title': 'T5 8 pt の全体（ラベル 194 個）', 'case': 't5_f8', 'panels': HJ2,
         'note': _two(cfg, 't5_f8') + ' 辺の規則で列が詰まり、キャンバスの高さが 1 割ほど減ります。'},
    ]
    LREC, LBOTH = V('L_f0o1', 'L', '線が背面を通るのを許す'), V('L', 'L（両方）', '文字がリング側へ戻るのも許す')
    GHJ = [V('G', 'G', '26 個を省略、間隔は守らない'), V('I', 'H2', '全部置く'), V('J', 'H3', '全部置く、直線を優先')]
    cfg['d11_over'] = [
        {'title': 'P1 14 pt の 12 時付近を原寸で: G、H2、H3', 'case': 'p1_f14', 'panels': GHJ, 'box': [-560, -620, 560, -300],
         'note': 'G は上の 24 個のうち 16 個だけを置いています（HNH endonuclease、PdcA、terminase 2 個、Lpa-like transcriptional activator、UpfN、UpfM、PmgV などがありません）。'
                 'Coi、ImcB、Mat はリングの内側です。「Doc-like toxin」や「site-specific integrase」は feature やリングの線に接しています。'},
        {'title': 'P1 14 pt の 6 時付近を原寸で: G、H2、H3', 'case': 'p1_f14', 'panels': GHJ, 'box': [-560, 300, 560, 620],
         'note': 'G は下の 23 個のうち 16 個を置いています。'},
        {'title': '試作 K: P1 14 pt の全体（H2、K の 3 段階）', 'case': 'p1_f14',
         'note': 'K は線がラベルの背面を通るのを許した版です。高さは H2 1,230 px、弱い設定 1,211 px、既定 1,165 px、強い設定 1,139 px。背面を通る線は 6 本、11 本、16 本。',
         'panels': [V('I', 'H2', '線は背面を通らない'), V('K_mild', 'K 弱', '1 本の線が通れるラベルは 1 個まで'), V('K', 'K 既定', '3 個まで'), V('K_strong', 'K 強', '6 個まで')]},
        {'title': '試作 K の見た目: P1 14 pt の 12 時付近を原寸で', 'case': 'p1_f14', 'box': [-560, -620, 560, -300],
         'note': '配置は同じで、描き方だけを変えています。下地は白で、濃さ 60 %、既定、100 %。「線を切る」はラベルの所で線を中断します（LoVis4u の描き方に近い）。'
                 '下線はラベルの下に線を引き、リーダー線がその端につながります。',
         'panels': [V('K_nobg', '下地なし', '線が文字の上を通る'), V('K_bg06', '下地 60 %', ''), V('K', '下地（既定）', ''), V('K_bg10', '下地 100 %', '線が止まって見える'),
                    V('K_gap', '線を切る', ''), V('K_ul', '下線 + 下地', 'いちばん読める'), V('K_ulnb', '下線のみ', ''), V('K_tcol', 'feature の色の文字', '薄い色が読みにくい')]},
        {'title': '試作 L: P1 14 pt の全体（H2、K、J2、L）', 'case': P1,
         'note': f"H2: {_lrow(cfg, P1, 'I')}。K: {_lrow(cfg, P1, 'K')}。J2: {_lrow(cfg, P1, 'J_sym2')}。L: {_lrow(cfg, P1, 'L_f0o1')}。L（両方）: {_lrow(cfg, P1, 'L')}。",
         'panels': [V('I', 'H2', '塔が残る'), V('K', 'K', '線が背面を通る'), V('J_sym2', 'J2', '対称な扇形、辺の規則'), LREC, LBOTH]},
        {'title': '試作 L: P1 14 pt の 12 時付近を原寸で（G、H2、J2、L）', 'case': P1, 'box': [-560, -620, 560, -300],
         'note': 'L は上の 24 個を全部置いたまま、3〜4 個のラベルが並ぶ行にまとめます。G と違って省略はなく、feature やリングの線に触れるラベルもありません。'
                 'Lpa-like transcriptional activator の背面を 6 本の線が通っています。',
         'panels': [V('G', 'G', '26 個を省略、間隔は守らない'), V('I', 'H2', '全部置く'), V('J_sym2', 'J2', '全部置く'), LREC]},
        {'title': '試作 L の 4 通り: P1 14 pt の 12 時付近を原寸で', 'case': P1, 'box': [-560, -620, 560, -300],
         'note': '(a) 線がラベルの背面を通る、(b) 文字がリング側へ戻る向きにも伸びる、の有り無しです。',
         'panels': [V('L_f0o0', 'どちらも無し', f"高さ {hh(P1, 'L_f0o0')} px"), V('L_f1o0', '(b) だけ', f"高さ {hh(P1, 'L_f1o0')} px"),
                    V('L_f0o1', '(a) だけ（推奨）', f"高さ {hh(P1, 'L_f0o1')} px"), V('L', '両方', f"高さ {hh(P1, 'L')} px")]},
        {'title': '試作 L: P1 14 pt の 6 時付近を原寸で（G、H2、J2、L）', 'case': P1, 'box': [-560, 300, 560, 620],
         'note': 'G は下の 23 個のうち 16 個だけを置いています。L は 23 個を全部置いたまま、複数のラベルが並ぶ行にまとめます。',
         'panels': [V('G', 'G', '26 個を省略、間隔は守らない'), V('I', 'H2', '全部置く'), V('J_sym2', 'J2', '全部置く'), LREC]},
        {'title': '試作 L の見た目: 下線つき（P1 14 pt の 12 時付近を原寸で）', 'case': P1, 'box': [-560, -620, 560, -300],
         'note': '配置は同じで、描き方だけを変えています。下線はラベルの下に線を引き、リーダー線がその端につながります（K の節の「下線 + 下地」と同じ描き方）。',
         'panels': [V('L_f0o1', '下地のみ', '背面を通る線は下地で薄くなる'), V('L_f0o1_ul', '下線 + 下地', '')]},
        {'title': '試作 L: N4 14 pt の全体（ラベル 69 個）', 'case': 'n4_f14',
         'note': f"H2: {_lrow(cfg, 'n4_f14', 'I')}。L: {_lrow(cfg, 'n4_f14', 'L_f0o1')}。",
         'panels': [V('I', 'H2', ''), LREC]},
        {'title': '試作 L: プラスミド 153 kb 8 pt の全体（ラベル 179 個）', 'case': 'plasmid153_f8',
         'note': f"H2: {_lrow(cfg, 'plasmid153_f8', 'I')}。L: {_lrow(cfg, 'plasmid153_f8', 'L_f0o1')}。背面を通る線が {lm('plasmid153_f8', 'L_f0o1')['leader_through_text']} か所あります。",
         'panels': [V('I', 'H2', ''), LREC]},
        {'title': '試作 L: T5 8 pt の全体（ラベル 194 個）', 'case': 't5_f8',
         'note': f"H2: {_lrow(cfg, 't5_f8', 'I')}。K: {_lrow(cfg, 't5_f8', 'K')}。L: {_lrow(cfg, 't5_f8', 'L_f0o1')}。"
                 f"背面を通る線が {lm('t5_f8', 'L_f0o1')['leader_through_text']} か所あり、ここまで混むと線は追いにくくなります。",
         'panels': [V('I', 'H2', ''), V('K', 'K', ''), LREC]},
    ]
    cfg['d11_questions'] = [
        '<b>J2（H3 + 対称な扇形 + 辺の規則）を新しい配置の既定にしてよいですか</b>。H2 との違いは、線の曲がりも数えること、塊を対称な扇形にそろえること、線の終点を枠の辺のどこでもよくしたことです。'
        '<span class="rec">推奨: はい。ご指摘の箇所は対称になり、線の合計は H3 より 5 % 短くなります。直線の本数は H3 より 1 割ほど減るので、曲がりの重みを実装のときに調整します。</span>',
        '<b>「ラベルがリング側へ戻らない」を、枠の高さまでは許す測り方にしてよいですか</b>（先頭の節の「角の規則と辺の規則」の図）。'
        '<span class="rec">推奨: はい。H2・H3 でも最大 13 px は起きていた量で、J2 は最大 8.7 px です。12 時・6 時の近くでは文字が線の左右にまたがれるようになります。</span>',
        '<b>12 時・6 時を「1 行に複数のラベル」で組む L を取り入れてよいですか</b>。そのためには、線がほかのラベルの背面を通ることを許す必要があります（ご指摘 3 の節の図）。'
        '<span class="rec">推奨: はい。P1 14 pt の高さが 1,230 px から 1,063 px になります。文字がリング側へ戻る向きに伸びること（L のもう 1 つの緩和）は、差が 16 px しかないので許しません。'
        '背面を通る線の描き方（下地のみか、下線つきか）は図を見てお決めください。</span>',
        '<b>トラックをまたぐ underlay と輪（葉緑体の図の Region annotation）をどのリリースに入れますか</b>。ご希望は承りました。新しい機能（アノテーションのトラックに「どのトラックからどのトラックまで」を指定する）になります。'
        '<span class="rec">推奨: 新しいラベル配置と同じ 0.15.0 で、別の機能として計画します。ギャラリーとチュートリアルの葉緑体の図はそのときに差し替えます。</span>',
    ]
