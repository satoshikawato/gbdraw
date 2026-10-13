"""Japanese texts of the Owner page, kept apart from the page structure: fill_texts.py rewrites page.json in place.

Run it before make_label_page.py. Facts come from logs/design-*.md, logs/part2-breaking-points.md,
logs/survey-tools.md and logs/review-figures.md.
"""
import json
import sys
from pathlib import Path

PAGE = Path(__file__).with_name('page.json')
cfg = json.loads(PAGE.read_text())

cfg['eyebrow'] = 'gbdraw · Circular のラベル配置 · 提案 2026-10-12'
cfg['heading'] = 'Circular のラベル配置: 今の配置と候補を同じ記録で比べる'

DESC = {
    'cur': '外側のラベルは全部ひとつの楕円の上に置き、文字は水平、引き出し線は折れ線。動かせるのは楕円上の角度だけです。'
           '入りきらないときの手当てがなく、残った重なりをそのまま描きます。',
    'A': '置き方は今のまま。入りきらないときだけ、文字を少し小さくし、ラベルの帯を少し広げて、重ならなくなる最初の組み合わせを使います。'
         '今すでに収まっている図は 1 バイトも変わりません。A1（文字だけ縮める）、A2（帯だけ広げる）、A3（ラベルを省く）は下の「設定で変えられるところ」に並べました。',
    'B': '内側の列は今の楕円のまま。入りきらないラベルを外側の列へ送り、内側のラベルの行間を通る細い線でつなぎます。'
         '何も省きません。代わりに横幅が増えます。',
    'D': '全部のラベルを左右 2 本の列にし、列をリングに沿わせます。時計回りの順番を保つので線は交差しません。'
         '何も省きません。代わりに縦が伸び、12 時と 6 時に平行線の「櫛」ができます。',
    'E': '近くに固まったラベルだけを、そろえた縦リストにまとめます（プラスミド地図の書き方）。離れたラベルは今と同じ見た目です。'
         '線は 1 行ずつ（fan）か、リスト全体に 1 本とブラケット（bracket）。既定では入りきらない分を省き、設定で省かずに縦を伸ばせます。',
    'R': 'すでにある <code>--label_placement radial</code>。文字を放射状に置き、重ならなくなるまで半径を広げます。'
         '重なりも交差も出ませんが、キャンバスが大きくなり、文字の向きがばらばらになります。',
}

cfg['metric_columns'] = [
    ['labels_external', '置いたラベル'], ['dropped', '省いたラベル'], ['font_size', '文字 pt'],
    ['text_overprint_pairs', '文字の重なり'], ['label_gap_min_px', 'ラベル間の最小すき間 px'],
    ['labels_toward_axis', 'リング側へ戻るラベル'],
    ['labels_on_objects', '目盛りなどに重なるラベル'], ['labels_near_objects', '目盛りなどから 3 px 未満のラベル'],
    ['leader_crossing_pairs', '線の交差'], ['leader_through_text', '文字を横切る線'],
    ['leader_bends_lt110', '折れが 110° 未満の線'],
    ['leader_len_sum', '線の合計 px'], ['leader_len_max', '線の最大長 px'], ['canvas_px', 'キャンバス px'],
]

cfg['findings'] = [
    '今の配置は「全部のラベルをひとつの楕円に並べる」方式で、<b>入る数に上限があります</b>。14 pt では楕円 1 周に 90〜110 個ほど。'
    '上限は記録の長さではなくラベルの数で決まります（N4 70 kb・69 個は 14 pt で問題なし、P1 95 kb・111 個は重なる）。',
    '<b>入りきらないときの手当てがありません</b>。どの最適化も「いちばんましな状態」で終わるので、残った重なりをそのまま描きます。'
    'P1 の左上の重なりはさらに単純な不具合で、最後のラベルが 269° に張り付き、直前の 1 個としか重なりを調べないため、1 個おきに同じ場所へ積まれます（OV-456 の原因）。',
    '重ならずに描ける最大の文字: mtDNA 26 pt、pVNTG4 18 pt、N4 16 pt、P1 12 pt、タバコ葉緑体 10 pt、153 kb プラスミド 8 pt。'
    'T5（194 個）は既定の 8 pt でも線が地図を横切ります（交差 2,085 組）。',
    'すでにある放射状配置（R）は重なりも交差も出しませんが、キャンバスが P1 14 pt で 1735 × 1580、T5 14 pt で 2288 × 2530 になります。',
    'ほかのツールで 100〜200 個のラベルに耐える方法を公開しているのは CGView.js / Proksee と Circos だけで、どちらも'
    '「線を伸ばして逃がす」と「入らないものは描かない」の組み合わせです。SnapGene、Geneious、MacVector も入らないラベルは隠します。',
    '簡単な地図（14 pt で 70 個くらいまで）では、今の配置がすでにきれいで、どの候補もほぼ同じ絵になります。'
    '差が出るのは入りきらない地図と、ラベルが 1 か所に固まる地図です。',
    '195 個を 14 pt で描く（T5、153 kb プラスミド）と、どの候補でも掲載できる図になりません。文字を小さくするか、ラベルを減らすしかありません。',
]

cfg['tools'] = [
    ['Geneious Prime', 'ラベルの位置は Inside / Outside / Inside or Outside から選べ、「Hide excessive labels」で多すぎるラベルを隠せます'
     '（<a href="https://manual.geneious.com/en/latest/Sequences.html">マニュアル</a>）。配置の方法そのものは公開されていません。'
     'いただいたスクリーンショットからは、feature のすぐ近くにリングの内外どちらにも置く、短い直線の引き出し線、密な所は斜めの階段、同じ名前の隣どうしをまとめる、文字を feature の色にする、と読み取れます。'],
    ['CGView.js / Proksee', 'ラベルを優先度順に並べ、上位 50 個は線を伸ばすか傾けて重ならない位置に置き、残りは重なるなら描きません。'
     '隣り合って重なるラベルを「島」にまとめて扇形に開きます（<a href="https://github.com/stothard-group/cgview-js/blob/main/src/LabelPlacementAngled.js">ソース</a>）。'],
    ['Circos', '同じ角度のまま外側へ積み、「snuggle」でリングに沿って少しずらします。それでも入らないラベルは描きません'
     '（<a href="https://circos.ca/documentation/tutorials/2d_tracks/text_2/">ドキュメント</a>）。'],
    ['SnapGene', '外側と内側のラベルを併用し、混んでいて入らないラベルは隠します。対処は文字を小さくするかウィンドウを広げること'
     '（<a href="https://www.snapgene.com/updates/snapgene-2-8-0-release-notes">リリースノート</a>）。'],
    ['CLC Genomics Workbench', '制限酵素サイトのラベルに「Radial」（切断点のすぐ近く）と「Stacked」（円の左右にそろえて全部並べる）の 2 通りがあります'
     '（<a href="https://resources.qiagenbioinformatics.com/manuals/clcgenomicsworkbench/700/Dynamic_restriction_sites.html">マニュアル</a>）。候補 D の「まっすぐな列」に当たります。'],
    ['SimVector', 'MCS の酵素を縦 1 列のリストにする「Bracket View」があります'
     '（<a href="https://premierbiosoft.com/plasmid_maps/VectorNTI.html">製品ページ</a>）。候補 E のブラケットに当たります。'],
    ['OGDRAW、DNAPlotter、GenomeDiagram', '重なりを避ける仕組みは見当たりません。名前を短くする、手で直す、が案内されています。'],
]

cfg['facts'] = [
    ['条件', 'どの図も <code>--separate_strands --track_type tuckin --labels out --legend right</code>。同じ記録・同じ文字サイズ・同じキャンバス基準で、ラベルの置き方だけを変えています。'],
    ['試作', '候補はすべて調査用の試作です（ローカルのブランチ <code>study/label-placement</code>、push していません）。製品コードは変えていません。'],
    ['図の見方', '各ケースの図は地図の中心をそろえて重ねてあります。タブで切り替えると、地図は動かずラベルだけが入れ替わります。図をクリックすると原寸になります。'],
]

NOTES = {
    'p1_f14': {
        'note': '最初にご指摘いただいたケースです。',
        'notes': {
            'cur': '左上で文字が重なって読めません（20 組）。楕円の上に置く場所がなくなり、最後の 12 個が同じ角度に積まれています。',
            'A': '見た目は今のまま、文字 12.5 pt・帯 1.15 倍で全 111 個が収まります。ラベルが外へ出るぶん引き出し線が長くなります。',
            'B': '全部置けます。左側は外の列へ向かう線が束になり、12 時の上に階段ができます。',
            'D': '今の見た目に近いまま全部置けます。高さが 19 % 増え、12 時の上に縦線の櫛ができます。',
            'E': '密な所だけ縦リストになります。既定では 12 個を省略。省かない設定だと高さが 1576 px になります。',
            'R': '重なりも交差もありませんが、キャンバスは 1735 × 1580。地図が小さく見えます。',
        }},
    'pbr322_f14': {
        'note': '制限酵素サイトが固まるプラスミドの典型です（サイトは Bio.Restriction で計算した 1 か所切断の酵素）。',
        'notes': {
            'cur': '重なりはありません。固まった所は扇形に開きます。',
            'A': 'A の判定は「重なり」より厳しく（文字を横切る線も不可）、13.5 pt に下がります。重なりだけで判定すれば今のままです。',
            'B': '今とほぼ同じです。',
            'D': '上のラベルが横へ引っ張られ、長い水平線が出ます。小さい地図では今より落ちます。',
            'E': '固まったサイトが 2〜4 行のそろったリストになります。プラスミド地図でよく見る書き方です。',
            'R': '放射状。小さい地図ではキャンバスの増え方は小さめです。',
        }},
    'mito_f14': {
        'notes': {
            'cur': '重なりなし。今のままで問題ありません。',
            'B': '今と同じです（1 段で足ります）。',
            'D': '列がリングに沿うだけで、今とほぼ同じです。',
            'E': 'tRNA が続く所（Trp・Ala・Asn・Cys・Tyr）が縦リストになります。',
        }},
    'mg1655_rna_f8': {
        'note': 'ラベルが数か所に固まる例です（rRNA オペロンとその周りの tRNA）。',
        'notes': {
            'cur': '重なりはありませんが、固まりごとに扇形に開き、どの線がどのラベルか追いにくくなります。',
            'A': 'A の厳しい判定では 6.5 pt まで下がります。重なりだけで判定すれば今のままです。',
            'B': '今とほぼ同じ見た目です。',
            'D': '全部が左右の列に並び、線が長くなります。',
            'E': 'オペロンごとに 1 つのリストになり、線は 1 本とブラケットだけ。この候補がいちばん合うケースです。',
        }},
    'tobacco_f14': {
        'notes': {
            'cur': '左上で文字が積み重なり（50 組）、線が 1,062 組交差します。',
            'A': '文字 11.5 pt・帯 1.25 倍で全 140 個が収まります（1785 × 1223）。',
            'B': '内側に短い名前、外側の列に長い名前。全部置けますが 1867 × 1252 に広がります。',
            'D': '全部置けて線も短いままですが、高さが 43 % 増えます（1555 × 1417）。',
            'E': '既定では 22 個を省略して 1533 × 1322。',
            'R': '2427 × 2206。Gallery のタバコ葉緑体はこの配置です。',
        }},
    't5_f14': {
        'note': '195 個を 14 pt で描く極端な例です。どの候補でも掲載できる図にはなりません。',
        'notes': {
            'cur': '170 組が重なり、線が 2,532 組交差します。',
            'A': '文字 8.5 pt・帯 1.5 倍まで下げてようやく収まります（1886 × 1448）。',
            'B': '4 段。2290 × 1414 に広がり、右半分は線の束で灰色の帯になります。',
            'D': '全部置けますが 1610 × 1965。上下の櫛が大きくなります。',
            'E': '既定では 85 個を省略。省かないと高さ 3081 px。',
            'R': '2288 × 2530。',
        }},
    't5_f8': {
        'note': '既定（Auto）の 8 pt です。重なりはありませんが、今の配置では線が地図を横切ります（OV-472）。',
        'notes': {
            'cur': '重なりは 0 ですが、線が地図の上を横切り（交差 2,085 組）、ラベルが feature から最大 629 px 離れます。',
            'A': '文字 6.5 pt・帯 1.2 倍。線は地図を横切らなくなりますが、文字が小さくなります。',
            'B': '右半分にラベルが偏るため、右側だけ 2 段になり、右上に大きな階段ができます。',
            'D': '左右 2 本の高い列（1359 × 1344）。交差はありません。',
            'E': '59 個を省略して 1329 × 1139。',
            'R': '1627 × 1668。',
        }},
    'plasmid153_f14': {
        'note': '179 個を 14 pt で描く極端な例です。',
        'notes': {
            'cur': 'ほぼ全周で文字が重なります（143 組）。',
            'A': '文字 10.5 pt・帯 1.35 倍で収まります（1862 × 1315）。',
            'B': '左右とも 2 段、2203 × 1352。',
            'D': '1650 × 1857。',
            'E': '既定では 64 個を省略。',
        }},
    'n4_f14': {
        'notes': {
            'cur': '重なりなし。今のままで問題ありません。',
            'B': '12 時を越えたラベルが地図の上に高い階段になります。今より落ちます。',
            'D': '今とほぼ同じ、少し縦に伸びます。',
            'E': '3 個を省略。今より良くはなりません。',
        }},
    'lambda_f14': {
        'notes': {
            'cur': '重なりなし。線の交差が 2 組あります。',
            'A': 'A の厳しい判定では 13.5 pt・帯 1.05 倍になります。重なりだけで判定すれば今のままです。',
        }},
    'plasmid37_f14': {
        'notes': {
            'cur': '重なりなし。今のままで問題ありません。',
            'D': '上のラベルが横へ引っ張られ、長い水平線が出ます。',
            'E': '3 時の固まりが縦リストになります。',
        }},
}

DESC['C'] = ('楕円をやめ、ラベルをリングの外の空いている場所へ、feature からまっすぐな線 1 本で届く範囲に置きます。'
             '混んだ所から順に、線が短く、ほかのラベルや線に当たらない最初の場所を取ります。線は交差させません。'
             '簡単な地図ではきれいな扇形になりますが、混んだ地図では置き場所を失ったラベルを省きます。')
DESC['G'] = ('Geneious Prime のように、ラベルを feature のすぐ横に置きます。リングの外だけでなく、何も描かれていなければ内側も使います。'
             '線は短い直線で、交差させません。ほかのラベル、線、feature、目盛り、GC トラック、中央のタイトルから 3 px 以上離します。'
             '中央が空いた地図ではいちばんコンパクトで、線の合計も最小です。混んだ地図や内側にトラックがある地図では、置けないラベルを省きます'
             '（外へ逃がす仕組みがまだありません）。')

DESC['I'] = ('H を、いただいた 4 つの指摘で直した版です。ラベルは 1 個ずつ、線がいちばん短くなる所（3 時・9 時では feature の横、12 時・6 時では真上・真下）から置き始め、'
             '順番を保ったまま、空いていなければ列に沿ってずらします。線は feature からリングの縁まで放射方向に出て、そこで 1 回だけ、110° 以上のゆるい角度で折れます。'
             'ラベルと線は、ほかのラベルだけでなく、目盛り・目盛りの文字・トラック・比較リングからも 3 px 以上離します（リングの内側でも外側でも）。'
             '内側ラベル（<code>--labels both</code>）も同じ規則で置きます。何も省かず、線は交差せず、文字はリングと反対の向きへ伸びます。')
DESC['H'] = ('<b>H2 の前の版です</b>（目盛りの文字に重なることがあり、12 時・6 時で線が横へ伸び、折れが直角に近い）。'
             'D と同じく、ラベルの順番を保った列に置くので、何も省かず、線も交差しません。違いは置く場所で、ラベルをリングの線から 3 px の所まで寄せ、'
             'feature のすぐ横に置きます（混んでいない所の線は 5〜10 px）。リングの内側に何も描かれていなければ、線の合計が短くなるラベルだけを内側の列へ移します。'
             '文字は、線がつく端からリングと反対の向きへ伸ばします（外側は外へ、内側は中心へ）。入りきらない分は列が上下に伸び、12 時と 6 時では縦の線になります。')
if not any(d['key'] == 'H' for d in cfg['designs']):
    cfg['designs'].insert(1, {'key': 'H', 'label': 'H', 'name': '順番を保って、リングの両岸に', 'en': 'ordered, both sides of the ring'})
if not any(d['key'] == 'I' for d in cfg['designs']):
    cfg['designs'].insert(1, {'key': 'I', 'label': 'H2', 'name': 'feature のすぐ横に、ゆるい折れ線で（推奨）', 'en': 'beside the feature, gentle elbows'})
for d in cfg['designs']:
    if d['key'] == 'H':
        d['name'] = '順番を保って、リングの両岸に（H2 の前の版）'
for sub in ('H', 'I'):
    sub_dir = str(Path(__file__).resolve().parents[1] / 'work' / sub)
    if sub_dir not in cfg['render_dirs']:
        cfg['render_dirs'].append(sub_dir)
for key in ('labels_on_objects', 'labels_near_objects', 'leader_bends_lt110'):
    if key not in cfg.setdefault('zero_is_good', []):
        cfg['zero_is_good'].append(key)

for d in cfg['designs']:
    d['desc'] = DESC.get(d['key'], d.get('desc', ''))
for case in cfg['cases']:
    extra = NOTES.get(case['id'])
    if extra:
        if 'note' in extra:
            case['note'] = extra['note']
        notes = dict(case.get('notes', {}))
        notes.update(extra['notes'])
        case['notes'] = notes


import figs  # noqa: E402  (crop boxes come from the drawn labels of today's render)


def dense_box(case, w=560, h=400, design='cur'):
    panel = figs.load(figs.find_render(cfg['render_dirs'], case, design))
    return list(figs.densest_box(panel['labels'], w, h))


V = lambda render, label, name: {'render': render, 'label': label, 'name': name}
A_SET = [V('cur', 'Current', '今の配置'), V('A', 'A', '文字と帯を一緒に（推奨）'), V('A1', 'A1', '文字だけ縮める'),
         V('A2', 'A2', '帯だけ広げる'), V('A3', 'A3', 'ラベルを省く')]
E_LEAD = [V('cur', 'Current', '今の配置'), V('E_fan', 'E fan', '1 行に 1 本の線'), V('E_bracket', 'E bracket', 'リストに 1 本の線とブラケット')]
cfg['variants'] = [
    {'title': 'A: 入りきらないとき、何を譲るか（P1 14 pt、左上を原寸で）',
     'note': '同じ 111 個のラベルです。A は文字 12.5 pt・帯 1.15 倍、A1 は 8 pt、A2 は 14 pt のまま帯 1.25 倍、A3 は 14 pt のまま 12 個を省略。',
     'case': 'p1_f14', 'panels': A_SET, 'box': dense_box('p1_f14')},
    {'title': 'A: 同じ 4 通りを地図全体で（P1 14 pt）',
     'note': 'キャンバスは A 1773 × 1130、A1 1441 × 982、A2 1928 × 1226、A3 1636 × 992（今は 1706 × 990）。',
     'case': 'p1_f14', 'panels': A_SET},
    {'title': 'A: T5 の既定 8 pt では',
     'note': '文字だけ縮めても（A1、6 pt）線は地図を横切ったままです。帯を広げる（A2）か、両方を少しずつ動かす（A、6.5 pt）と解消します。A3 は 42 個を省略。',
     'case': 't5_f8', 'panels': A_SET},
    {'title': 'E: 線の引き方（pBR322、原寸）',
     'note': 'fan は 1 行ずつ feature まで線を引きます。bracket は固まり全体にブラケットを掛け、線は 1 本だけ。静かですが、どの行がどのサイトかは並び順で読むことになります。',
     'case': 'pbr322_f14', 'panels': E_LEAD, 'box': dense_box('pbr322_f14')},
    {'title': 'E: 線の引き方（MG1655 の tRNA・rRNA、原寸）',
     'note': 'rRNA オペロンのように 5〜8 個が 1 か所に固まる所では、bracket のほうが読みやすくなります。',
     'case': 'mg1655_rna_f8', 'panels': E_LEAD, 'box': dense_box('mg1655_rna_f8')},
    {'title': 'E: 入りきらないとき、省くか伸ばすか（P1 14 pt）',
     'note': '既定は 12 個を省いて 1596 × 1315。省かない設定では全 111 個を置いて 1596 × 1576。',
     'case': 'p1_f14', 'panels': [V('cur', 'Current', '今の配置'), V('E', 'E', '入らない分を省く（既定）'), V('E_grow', 'E grow', '省かずに縦を伸ばす')]},
    {'title': 'D: リングに沿う列か、まっすぐな列か（pBR322）',
     'note': 'まっすぐな列（CLC の「Stacked」と同じ考え方）は表のように整いますが、線が長くなります（最長 530 px）。',
     'case': 'pbr322_f14', 'panels': [V('cur', 'Current', '今の配置'), V('D', 'D', 'リングに沿う列（既定）'), V('D_flat', 'D flat', 'まっすぐな列')]},
]



def load_metrics(case, design):
    stem = figs.find_render(cfg['render_dirs'], case, design)
    return json.loads(Path(str(stem) + '.metrics.json').read_text()) if stem else None


def placed_note(case, design, tail=''):
    """One sentence of facts read from the render's metrics, so the numbers cannot drift from the figure."""
    m = load_metrics(case, design)
    if not m:
        return None
    proto = m.get('proto') or {}
    dropped = len(proto.get('dropped') or [])
    inside = (proto.get('notes') or {}).get('inside')
    text = f"{m['labels_external'] + dropped} 個のうち {m['labels_external']} 個を配置"
    text += f"、{dropped} 個を省略" if dropped else '（全部）'
    if design in ('G', 'H', 'I') and inside:
        text += f"、{inside} 個はリングの内側"
    if m.get('leader_crossing_pairs'):
        text += f"、交差 {m['leader_crossing_pairs']} 組"
    w, h = m['canvas_px']
    return f"{text}。キャンバス {w} × {h}。{tail}"


TAILS = {
    ('p1_f14', 'C'): '入りきらない分は地図の上下に積まれます。',
    ('p1_f14', 'G'): '内側は GC トラックでふさがっていて、ほとんど使えません。',
    ('pbr322_f14', 'G'): 'ラベルがリングに寄り添い、今より小さく収まります。',
    ('mito_f14', 'G'): '線の合計は 1,351 px（今は 2,898 px）。',
    ('mg1655_rna_f8', 'G'): '固まりが大きすぎて置き場所が足りません。この配置が苦手な形です。',
    ('plasmid153_f8', 'C'): '今は全部置けているケースです。',
    ('plasmid153_f8', 'G'): '今は全部置けているケースです。',
    ('lambda_f14', 'G'): '内側がふさがっているため、外側だけの結果とほぼ同じです。',
    ('pbr322_open_f14', 'G'): '内側に入るのは ClaI、HindIII、BstZ17I の 3 個です。',
    ('p1_open_f14', 'G'): '内側を使わない設定では省略が 28 個に増えます。',
    ('p1_f14', 'H'): '線の合計 10,564 px（今は 12,731 px で重なり 20 組、D は 12,526 px）。内側は GC トラックでふさがっているので全部外側です。12 時と 6 時の縦線は残ります。',
    ('p1_open_f14', 'H'): '線の合計 6,699 px（D 12,526 px。G は 9 個を省いて 4,622 px）。',
    ('t5_f8', 'H'): '今は線の交差が 2,085 組。H は交差 0 で、線の合計は 36,908 px（今は 67,238 px）。高さは増えます。',
    ('t5_f14', 'H'): '全部置けますが、縦に 2 倍になります。この記録を 14 pt で描くこと自体に無理があります。',
    ('plasmid153_f14', 'H'): '全部置けますが、縦に 1.9 倍になります。この記録を 14 pt で描くこと自体に無理があります。',
    ('tobacco_f14', 'H'): '今は重なり 50 組・交差 1,062 組。H は全部置けて、線の合計 20,019 px（今は 37,840 px）。高さは増えます。',
    ('p1_f14', 'I'): '線の合計 8,975 px（今は 12,731 px で重なり 20 組、H は 10,564 px）。いちばん長い線は 260 px（H 413 px）。12 時と 6 時に積まれるラベルは残ります。',
    ('p1_open_f14', 'I'): '線の合計 5,107 px（H 6,699 px。G は 9 個を省いて 4,622 px）。',
    ('t5_f8', 'I'): '交差 0、線の合計 36,803 px（今は 67,238 px、H 36,908 px）。<b>H より縦に伸びます</b>（H は 1318 × 1354）。片側にラベルが偏った記録が苦手です。',
    ('t5_f14', 'I'): '全部置けますが、<b>H よりさらに縦に伸びます</b>（H は 1566 × 1985）。この記録を 14 pt で描くこと自体に無理があります。',
    ('plasmid153_f14', 'I'): '全部置けますが、縦に 1.9 倍になります（H とほぼ同じ）。この記録を 14 pt で描くこと自体に無理があります。',
    ('tobacco_f14', 'I'): '今は重なり 50 組・交差 1,062 組。H2 は全部置けて、線の合計 17,411 px（今は 37,840 px、H 20,019 px）。高さは増えます。',
}
OPEN_NOTE = ('内側の GC トラックを消した地図です（<code>--no-gc --no-skew</code>）。Geneious の図と同じく中央が空いているので、'
             '候補 G、H、H2 が内側を使えます。ほかの候補は外側しか使わないため、GC トラックありの図と同じ配置です。')
for case in cfg['cases']:
    notes = dict(case.get('notes', {}))
    for design in ('C', 'G', 'F', 'H', 'I'):
        tail = TAILS.get((case['id'], design), '')
        if design in ('H', 'I') and not tail:
            h, cur = load_metrics(case['id'], design), load_metrics(case['id'], 'cur')
            tail = f"線の合計 {round(h['leader_len_sum']):,} px（今は {round(cur['leader_len_sum']):,} px）。" if h and cur else ''
        note = placed_note(case['id'], design, tail)
        if note:
            notes[design] = note
    if notes:
        case['notes'] = notes
    if case['id'].endswith('_open_f14'):
        case['note'] = OPEN_NOTE

G_SET = [V('cur', 'Current', '今の配置'), V('G', 'G', '内側も使う・間隔 3 px（既定）'), V('G_in0', 'G in0', '外側だけ'),
         V('G_clear6', 'G clear6', '間隔 6 px'), V('G_cross', 'G cross', '線の交差を許す'),
         V('G_merge', 'G merge', '同じ名前の隣どうしをまとめる'), V('G_color', 'G color', '文字を feature の色に')]
cfg['variants'] = [
    {'title': 'G: 内側を使うか、間隔をどれだけ取るか（P1、内側のトラックなし、原寸）',
     'note': '既定は 111 個のうち 102 個を配置（25 個が内側）。外側だけにすると省略は 9 個から 28 個に増えます。間隔を 6 px にすると省略 12 個。'
             '交差を許すと省略は 4 個に減りますが、交差が 14 組出ます。まとめる設定は「名前 (×2)」と書きます。色付きの文字は薄い色で読みにくくなります。',
     'case': 'p1_open_f14', 'panels': G_SET, 'box': dense_box('p1_open_f14', 640, 460, 'G')},
    {'title': 'G: 同じ設定を λ で（内側のトラックなし、原寸）',
     'note': '既定は 71 個のうち 68 個を配置（9 個が内側）。外側だけでは省略 8 個、間隔 6 px では省略 4 個。',
     'case': 'lambda_open_f14', 'panels': G_SET, 'box': dense_box('lambda_open_f14', 640, 460, 'G')},
    {'title': 'G: pBR322 では（内側のトラックなし、地図全体）',
     'note': '小さいプラスミドではほとんどのラベルがリングのすぐ外に付き、固まった所だけ内側に入ります。外側だけにすると 3 個を省略。',
     'case': 'pbr322_open_f14', 'panels': G_SET[:4]},
] + cfg['variants']


# ---- texts that depend on every candidate's results -------------------------------------------------------------
for d in cfg['designs']:
    if d['key'] == 'F':
        d['desc'] = ('「重ならない・間隔を保つ・交差しない・省かない」を守ったうえで、リード線の合計が最小になる位置を厳密に解きます。'
                     '厳密に解けたのは、順番を保って 1 本の楕円に並べる形だけでした（交差しないことを証明できる形）。'
                     '入りきらないときは楕円を 0.05 倍ずつ広げます。楕円に縛られない版は交差を消しきれず、比較用にだけ載せています。')

cfg['intro'] = ('P1 を 14 pt で描くとラベルが重なる件（OV-456）から始めた調査の結果です。23 の記録で、今の配置と 10 の候補を'
                '同じ記録・同じ文字サイズで描きました。製品のコードはまだ何も変えていません。途中でいただいた条件'
                '（楕円に縛らない、ほかのオブジェクトと最低限の間隔を取る、リード線を交差させない、ラベルをリング側へ戻さない、ラベルを省かない、線の合計を短く）で'
                '各候補を比べています。<b>2026-10-12 追記</b>: 候補 H へのご指摘（12 時・6 時で線が横へ伸びる、目盛りにラベルが重なる、線の折れが鋭い、'
                '内側ラベルの線が長い）を受けて、直した版 <b>H2</b> を追加し、推奨を H2 に替えました。'
                '<a href="#docs">ご指摘の箇所の H と H2 の比較</a>から先にご覧ください。')
cfg['hero']['title'] = '同じ場所を原寸で: P1 14 pt の左上（今の配置で文字が重なる所）'

WITH_TRACKS = [c['id'] for c in cfg['cases'] if not c['id'].endswith('_open_f14')]


def across(design):
    rows = [m for m in (load_metrics(c, design) for c in WITH_TRACKS) if m]
    gaps = [m['label_gap_min_px'] for m in rows if m.get('label_gap_min_px') is not None]
    crossing = [m['leader_crossing_pairs'] for m in rows]
    drops = [len((m.get('proto') or {}).get('dropped') or []) for m in rows]
    return rows, (min(gaps) if gaps else None), crossing, drops


def turned(design, cases=None):
    """Records and labels where a label runs back toward the ring (LP-D09)."""
    vals = [m.get('labels_toward_axis') or 0 for m in (load_metrics(c, design) for c in (cases or WITH_TRACKS)) if m]
    return sum(1 for v in vals if v), sum(vals)


def compare_row(design, good, weak, pick=False):
    rows, gap, crossing, drops = across(design)
    if gap is None:
        gap_cell = '文字が斜めのため未計測'
    elif gap < 0:
        gap_cell = f'<span class="bad">重なる（{-gap:g} px 食い込む）</span>'
    elif gap < 2.95:
        gap_cell = f'<span class="bad">{gap:g} px まで詰まる</span>'
    else:
        gap_cell = f'<span class="ok">{gap:g} px 以上</span>'
    hit = sum(1 for v in crossing if v)
    cross_cell = '<span class="ok">0</span>' if not hit else f'<span class="bad">{hit} 記録で発生（最大 {max(crossing):,} 組）</span>'
    lost = sum(1 for v in drops if v)
    drop_cell = '<span class="ok">0</span>' if not lost else f'<span class="bad">{lost} 記録で計 {sum(drops):,} 個</span>'
    back_n, back_sum = turned(design)
    back_cell = '<span class="ok">0</span>' if not back_n else f'<span class="bad">{back_n} 記録で計 {back_sum} 個</span>'
    every = [m for m in (load_metrics(c['id'], design) for c in cfg['cases']) if m]  # all records, as in the recommendation
    on = sum(m.get('labels_on_objects') or 0 for m in every)
    near = sum(m.get('labels_near_objects') or 0 for m in every)
    if on:
        obj_cell = f'<span class="bad">重なる {on} 個、3 px 未満 {near} 個</span>'
    elif near:
        obj_cell = f'<span class="bad">3 px 未満 {near} 個</span>'
    else:
        obj_cell = '<span class="ok">0</span>'
    bends = sum(m.get('leader_bends_lt110') or 0 for m in every)
    leaders = sum(m.get('labels_external') or 0 for m in every)
    bend_cell = '<span class="ok">0</span>' if not bends else f'{bends:,} / {leaders:,} 本'
    p1 = load_metrics('p1_f14', design)
    p1_drop = len((p1.get('proto') or {}).get('dropped') or [])
    size = '/'.join(f'{v:g}' for v in p1.get('font_size', []))
    sum_cell = f"{round(p1['leader_len_sum']):,} px" + (f"（{p1['labels_external']} 個分）" if p1_drop else '') + ('' if size == '14' else f"（{size} pt）")
    return {'design': design, 'pick': pick,
            'cells': [gap_cell, cross_cell, back_cell, obj_cell, bend_cell, drop_cell, sum_cell, ' × '.join(str(v) for v in p1['canvas_px']), good, weak]}


cfg['compare'] = {
    'columns': ['ラベル間の最小すき間', '線の交差', 'リング側へ戻るラベル', '目盛り・トラックに重なる、近すぎるラベル（23 記録）', '折れが 110° 未満の線（23 記録）', '省いたラベル', '線の合計（P1 14 pt）', 'キャンバス（P1 14 pt）', '向いている地図', '弱い所'],
    'rows': [
        compare_row('cur', '14 pt で 70 個くらいまで', '入りきらないと重なる。線が交差することがある。'),
        compare_row('I', '全部のラベルを feature のすぐ横に置きたいとき。中央が空いていれば内側も使う。内側ラベルも同じ規則で置く', '12 時と 6 時にラベルが積まれる（H より短い）。片側に偏ったとても混んだ記録では H より縦に伸びる。', pick=True),
        compare_row('H', 'H2 の前の版', '目盛りの文字に近すぎる。12 時・6 時で線が横へ伸びる。折れが直角に近い。'),
        compare_row('A', '今の見た目を変えたくないとき', '文字が小さくなるか、キャンバスが広がる。すき間と交差の考え方は今のまま。'),
        compare_row('B', '名前の長さがばらつく 100〜180 個', '外の列へ向かう線が束になる。横に広がる。'),
        compare_row('C', 'まばらな 100 個くらいまで', '混むとラベルを省く。'),
        compare_row('D', '全部のラベルを確実に置きたいとき', '縦に伸びる。12 時と 6 時に平行線の櫛。外側しか使わない。'),
        compare_row('E', 'ラベルが数か所に固まる地図（制限酵素サイト、rRNA オペロン）', '既定では省く。省かない設定では D と同じく縦に伸びる。'),
        compare_row('G', '中央が空いた 100 個くらいまで', '混むと省く。内側にトラックがあると内側を使えない。'),
        compare_row('F', '今きれいに描けている地図に、交差なしの保証を付けたいとき', '混むと楕円ごと大きくなり、線が長くなる。'),
        compare_row('R', 'キャンバスが大きくてよいとき', 'キャンバスが最大 2.5 倍。文字の向きがそろわない。'),
    ],
}

ALL_CASES = [c['id'] for c in cfg['cases']]


def h_facts():
    """Numbers of the recommendation, read from the renders."""
    shortest = 0
    clean = []
    for case in ALL_CASES:
        sums = {d: load_metrics(case, d)['leader_len_sum'] for d in 'BDFHI'}
        shortest += min(sums, key=sums.get) == 'I'
        cur = load_metrics(case, 'cur')
        if not cur['text_overprint_pairs']:
            clean.append(1 - sums['I'] / cur['leader_len_sum'])
    for d in 'BDFHI':
        for case in ALL_CASES:
            m = load_metrics(case, d)
            assert not (m.get('proto') or {}).get('dropped') and not m['leader_crossing_pairs'] and not m.get('labels_toward_axis'), (case, d)
    return shortest, round(min(clean) * 100), round(max(clean) * 100)


H_SHORTEST, H_MIN, H_MAX = h_facts()
P1H, P1O = load_metrics('p1_f14', 'H'), load_metrics('p1_open_f14', 'H')
P1I, P1IO = load_metrics('p1_f14', 'I'), load_metrics('p1_open_f14', 'I')


def objects_all(design):
    rows = [load_metrics(c, design) for c in ALL_CASES]
    return (sum(m.get('labels_on_objects') or 0 for m in rows), sum(m.get('labels_near_objects') or 0 for m in rows),
            sum(m.get('leader_bends_lt110') or 0 for m in rows), sum(m.get('labels_external') or 0 for m in rows))


H_OBJ, I_OBJ, CUR_OBJ = objects_all('H'), objects_all('I'), objects_all('cur')
assert I_OBJ[:3] == (0, 0, 0), I_OBJ
assert H_SHORTEST == len(ALL_CASES), H_SHORTEST
cfg['recommendation'] = [
    '<b>進む方向</b>: 候補 <b>H2</b> を新しい配置として勧めます（H を、いただいた指摘で直した版）。条件（最低限の間隔・線が交差しない・ラベルがリング側へ戻らない・省かない）を '
    f'{len(ALL_CASES)} 記録すべてで満たし、順番を保つ候補（B・D・F・H・H2）の中で、全 {H_SHORTEST} 記録で線の合計がいちばん短くなりました。'
    f'今の配置と比べても、今重ならずに描けている記録で {H_MIN}〜{H_MAX} % 短くなります（λ 14 pt: 8,171 px → 4,091 px）。'
    f"P1 14 pt は 111 個を全部置いて {round(P1I['leader_len_sum']):,} px です（今は 12,731 px で重なり 20 組、H は {round(P1H['leader_len_sum']):,} px）。",
    '<b>H から直した所</b>: (1) 12 時・6 時の近くのラベルを、円のすぐ外で feature の真上・真下に置きます（ヒト mtDNA の図の ATP8: 線が 77 px → 24 px）。'
    f'(2) 目盛り・目盛りの文字・トラック・比較リングから 3 px 以上離します。{len(ALL_CASES)} 記録で、H は目盛りなどに重なるラベルが {H_OBJ[0]} 個、3 px 未満が {H_OBJ[1]} 個ありましたが、H2 はどちらも 0 です。'
    f'(3) 線の折れは 1 回だけで、110° 以上です。折れが 110° 未満の線は、H が {H_OBJ[3]:,} 本のうち {H_OBJ[2]:,} 本、今の配置が {CUR_OBJ[2]:,} 本、H2 は 0 本です。'
    '(4) 内側ラベル（<code>--labels both</code>）も同じ規則で置きます。葉緑体の図では feature 側の線が 47 px → 11 px になりました。',
    '<b>訂正</b>: 前の版のこのページは、H について「重なりはありません」と書いていました。数えていたのはラベル同士の重なりだけで、'
    '目盛りの文字との重なりを見ていませんでした（チュートリアル「Highlight」の図で、H は「ATP8」を「9 kbp」の上に置きます）。'
    '今の版は、ラベルのない同じ図と重ね合わせて、描かれているすべての物との間隔を数えています。',
    '<b>見た目</b>: 簡単な地図では、ラベルが feature のすぐ横に付きます（pBR322、λ、mtDNA）。中央が空いた地図では内側も使い、G（Geneious 風）に近い絵になります。'
    f"P1 の内側トラックなしでは 111 個のうち {P1IO['proto']['notes']['inside']} 個が内側で、線の合計は {round(P1IO['leader_len_sum']):,} px です（H 6,699 px、D 12,526 px。G は 9 個を省いて 4,622 px）。"
    '内側に GC トラックがある既定の図では、外側ラベルは内側に入りません。',
    '<b>H2 に残る弱点</b>: (a) ラベルが片側に偏ったとても混んだ記録は、H より縦に伸びます（T5 8 pt: H 1318 × 1354 → H2 1328 × 1678、T5 14 pt: 1566 × 1985 → 1620 × 2630）。'
    '折れを 110° 以上に保つため、12 時・6 時をまたいで反対側の列へ逃がせるラベルが減るからです。'
    f"(b) 12 時と 6 時の近くにラベルが積まれるのは残ります（P1 14 pt で上に 24 個、下に 18 個。いちばん長い線は {round(P1I['leader_len_max'])} px、H は {round(P1H['leader_len_max'])} px）。"
    '(c) 内側のレーンにある feature は、リングの縁までの線が長くなります（ヒト mtDNA の tRNA-Pro: 77 px のうち 62 px がリングの中）。'
    '(d) ほかのラベルの線との間隔は最小 2.4 px で、3 px を保証していません。(e) 葉緑体の図では、内側に頼んだ 69 個のうち 2 個（tRNA-His、rps14）が入りきらず外側に出ます。'
    '(f) 目盛りの間隔を手で指定した図は、試作では未対応です。この 6 点は実装計画で扱います。',
    '<b>線の合計について</b>: 目的関数（合計・二乗和・最大値）を変えても結果はほとんど変わりませんでした（F）。線が短くなるのは、feature の近くに置き場所があるときです。'
    'H2 が短いのは、リングのすぐ横と、12 時・6 時では真上・真下を使うからです。',
    '<b>入りきらないとき</b>: 既定では省かず、列を上下に伸ばします（キャンバスが大きくなる）。省くのは使う人が選んだときだけにして、省いた数を知らせます。',
    '<b>段階</b>: (1) 0.14.x で「今の見た目のまま重ねない」（A）と 269° の不具合の修正。今収まっている図は 1 バイトも変わりません。'
    '(2) 0.15.0 で H2 を <code>--label_placement</code> の新しい値として追加し、Gallery で確かめてから既定にするかを決める。'
    '(3) そのあとで、E のブラケット付きリストを「固まったラベルの書き方」として追加。',
    '<b>採らない候補</b>: B（線の束が読みにくい）、C と G（混むとラベルを省き、ラベルがリング側へ戻る。G の置き方は H2 が引き継ぎます）、D と H（H2 が同じ保証で線が短い）、'
    'F の楕円版（混むと楕円ごと大きくなる）。',
]

cfg['questions'] = [
    '<b>新しい配置は H2 でよいですか</b>（feature のすぐ横、12 時・6 時では真上・真下、折れは 110° 以上、目盛りやトラックから 3 px、交差なし・リング側へ戻らない・省かない）。'
    '<span class="rec">推奨: はい。片側に偏った混んだ記録が縦に伸びる点（T5）は、実装計画で H の「反対側の列へ逃がす」方法と組み合わせて直します。</span>',
    '<b>入りきらないときに何を譲るか</b>: (a) キャンバスを伸ばす、(b) 文字を小さくする、(c) ラベルを省く。'
    '<span class="rec">推奨: 既定は (a)。文字サイズが Auto のときだけ (b) も使う。(c) は明示的に選んだときだけで、省いた数を表示する。</span>',
    '<b>リングの内側</b>を既定で使いますか。内側のラベルは、リング側の端に線をつけ、中心へ向かって伸ばしています'
    '（「リング側へ戻らない」をこう読みました。地図の中心から遠ざかる向きという意味でしたら、内側は使えなくなります）。'
    '<span class="rec">推奨: 使う。何も描かれていない所だけで、線の合計が短くなるラベルだけ。GC トラックがある既定の図では内側に入りません。</span>',
    '<b>0.14.x に「今の見た目のまま重ねない」（A）を先に入れますか</b>。重なりが残ったときだけ働きます。'
    '<span class="rec">推奨: 入れる。文字サイズを明示したときは文字を変えず、帯だけ広げる（A2）。Auto のときは文字と帯を一緒に動かす（A）。'
    '変えたときはログと Web に知らせを出す。</span>',
    '<b>追加の表現</b>: 同じ名前の隣どうしをまとめる（「hypothetical protein (×2)」）、文字を feature の色にする、ブラケット付きの縦リスト（E）。'
    '<span class="rec">推奨: まとめるのはオプションとして入れる。色は見送り（薄い色が読めない）。ブラケットは 0.15.0 のあとで。</span>',
    '<b>Auto の文字サイズ</b>を、記録の長さ（今は 75 kb が境界）ではなく「ラベルが収まる最大の文字」で決めますか。'
    '<span class="rec">推奨: 新しい配置と一緒に 0.15.0 で。長さでは決まりません（N4 70 kb・69 個は 14 pt で収まり、P1 95 kb・111 個は収まらない）。'
    '100 kb への引き上げはそれまで保留のままにします。</span>',
    '<b>既定で「hypothetical protein」などを隠しますか</b>。今は CLI も Web も隠しません（設定ファイルの語リストは使われていません）。'
    '<span class="rec">推奨: 今のまま隠さない。T5 は 194 個のうち 90 個がこの種の名前で、Web の Blacklist を使えば今の配置でもきれいに描けます。</span>',
]

_g, _c, _cur = turned('G', ALL_CASES), turned('C', ALL_CASES), turned('cur', ALL_CASES)
assert all(turned(d, ALL_CASES) == (0, 0) for d in 'BDFHI'), 'B, D, F, H, H2 must not turn labels toward the ring'
BACK_FINDING = ('<b>ラベルがリング側へ戻る</b>（線がつく端より、反対の端のほうがリングに近い）のは、G が ' + f'{len(ALL_CASES)} 記録のうち {_g[0]} 記録・計 {_g[1]} 個、'
                f'C が {_c[0]} 記録・{_c[1]} 個、今の配置が {_cur[0]} 記録・{_cur[1]} 個（T5 8 pt で 19 個）でした。'
                '順番を保つ B・D・F・H・H2 は全記録で 0 です。リングの内側のラベルは、中心へ向かって伸びるものを可として数えています。')
OBJ_FINDING = ('<b>目盛りやトラックとの間隔</b>を、ラベルのない同じ図と重ねて数えました。' + f'{len(ALL_CASES)} 記録で、H は目盛りなどに重なるラベルが {H_OBJ[0]} 個、3 px 未満が {H_OBJ[1]} 個。'
               f"G は重なるラベルが {objects_all('G')[0]} 個。今の配置は、ラベルは 0 個ですが、線がトラックや目盛りの上を通ります"
               f"（T5 8 pt で {load_metrics('t5_f8', 'cur')['leaders_on_objects']} 本、T5 14 pt で {load_metrics('t5_f14', 'cur')['leaders_on_objects']} 本）。H2 はラベルも線も 0 です。"
               f"<b>線の折れ</b>が 110° 未満の線は、今の配置 {CUR_OBJ[2]:,} 本、H {H_OBJ[2]:,} 本、D {objects_all('D')[2]:,} 本、H2 0 本（{H_OBJ[3]:,} 本のうち）。")
cfg['findings'] = [
    '今の配置は「全部のラベルをひとつの楕円に並べる」方式で、<b>入る数に上限があります</b>（14 pt で 1 周に 90〜110 個ほど）。'
    '上限は記録の長さではなくラベルの数で決まります。',
    '<b>入りきらないときの手当てがありません</b>。残った重なりをそのまま描きます。P1 の左上はさらに単純な不具合で、'
    '最後のラベルが 269° に張り付き、直前の 1 個としか重なりを調べないため、1 個おきに同じ場所へ積まれます。',
    '今の配置は重なりを避けるだけで、<b>ラベル間のすき間は 0.1〜0.4 px まで詰まります</b>（λ、mtDNA、T5 8 pt）。試作の候補は 3 px を保っています。',
    '重ならずに描ける最大の文字: mtDNA 26 pt、pVNTG4 18 pt、N4 16 pt、P1 12 pt、タバコ葉緑体 10 pt、153 kb プラスミド 8 pt。'
    'T5（194 個）は既定の 8 pt でも線が地図を横切ります（交差 2,085 組）。',
    '<b>簡単な地図</b>（14 pt で 70 個くらいまで）では今の配置がすでにきれいで、順番を保つ候補はどれもほぼ同じ絵になります。',
    '<b>自由に置く候補（C・G）は混んだ所でラベルを省きます</b>。交差を禁じると省略が増えます（P1 14 pt で C 24 個、G 26 個）。'
    '順番を保った列（H、H2）なら、同じ「feature のすぐ横」を省略なしで描けます。',
    BACK_FINDING,
    OBJ_FINDING,
    '<b>195 個を 14 pt で描く</b>（T5、153 kb プラスミド）と、どの候補でも掲載できる図になりません。文字を小さくするか、ラベルを減らすしかありません。',
    'ほかのツールで 100〜200 個に耐える方法を公開しているのは CGView.js / Proksee と Circos だけで、どちらも'
    '「線を伸ばして逃がす」と「入らないものは描かない」の組み合わせです。SnapGene、Geneious、MacVector も入らないラベルは隠します。',
]


def have(case, panels):
    return [p for p in panels if figs.find_render(cfg['render_dirs'], case, p['render'])]


F_P1 = [V('cur', 'Current', '今の配置'), V('F', 'F', '楕円の上で厳密に（既定）'), V('F_c6', 'F clear6', '間隔 6 px'),
        V('F_free', 'F free', '楕円に縛らない版（交差が残る）'), V('D', 'D', '左右の列（比較）')]
F_LAMBDA = [V('cur', 'Current', '今の配置'), V('F', 'F', '合計を最小に'), V('F_max', 'F max', '最大値を最小に'),
            V('F_free', 'F free', '楕円に縛らない版（交差が残る）'), V('D', 'D', '左右の列（比較）')]
BACK_SET = [V('cur', 'Current', '今の配置'), V('C', 'C', '外側に自由に'), V('G', 'G', 'Geneious 風'), V('D', 'D', '左右の列'), V('H', 'H', '順番を保って両岸に')]
H_SET = [V('cur', 'Current', '今の配置'), V('D', 'D', '左右の列'), V('G', 'G', 'Geneious 風（省略あり）'), V('H', 'H', '内側も使う・間隔 3 px（既定）'),
         V('H_in0', 'H in0', '外側だけ'), V('H_clear6', 'H clear6', '間隔 6 px')]
H_PLAIN = [V('cur', 'Current', '今の配置'), V('D', 'D', '左右の列'), V('G', 'G', 'Geneious 風（省略あり）'), V('H', 'H', '順番を保って両岸に')]


def hsum(case, render):
    return f"{round(load_metrics(case, render)['leader_len_sum']):,} px"


def canvas(case, render):
    return ' × '.join(str(v) for v in load_metrics(case, render)['canvas_px'])


H2_PLAIN = [V('cur', 'Current', '今の配置'), V('H', 'H', '前の版'), V('I', 'H2', '直した版（推奨）')]
H2_OPEN = [V('cur', 'Current', '今の配置'), V('G', 'G', 'Geneious 風（省略あり）'), V('H', 'H', '前の版'), V('I', 'H2', '直した版（推奨）')]
H2_SET = [V('I', 'H2', '内側も使う・間隔 3 px（既定）'), V('I_in0', 'H2 in0', '外側だけ'), V('I_clear6', 'H2 clear6', '間隔 6 px')]
cfg['variants'] = [
    {'title': 'H2: 内側に GC トラックがある既定の図（P1 14 pt、地図全体）',
     'note': f"111 個を全部置きます。線の合計は今 12,731 px（重なり 20 組）、H {hsum('p1_f14', 'H')}、H2 {hsum('p1_f14', 'I')}。"
             f"キャンバスは今 1706 × 990、H {canvas('p1_f14', 'H')}、H2 {canvas('p1_f14', 'I')}。",
     'case': 'p1_f14', 'panels': have('p1_f14', H2_PLAIN)},
    {'title': 'H2: 同じ図の 12 時付近を原寸で（P1 14 pt）',
     'note': 'H は線が直角に折れて真上へ伸びます。H2 は折れが 110° 以上で、ラベルが feature の真上に近づきます。積まれるラベルの数はほぼ同じです。',
     'case': 'p1_f14', 'panels': have('p1_f14', H2_PLAIN[1:]), 'box': [-340, -640, 300, -280]},
    {'title': 'H2: 中央が空いた地図で（P1、内側のトラックなし、地図全体）',
     'note': f"H2 は 111 個のうち {P1IO['proto']['notes']['inside']} 個が内側で、線の合計 {hsum('p1_open_f14', 'I')}（H は {P1O['proto']['notes']['inside']} 個が内側で {hsum('p1_open_f14', 'H')}）。"
             f"G は 9 個を省いて 4,622 px。キャンバスは H {canvas('p1_open_f14', 'H')}、H2 {canvas('p1_open_f14', 'I')}。",
     'case': 'p1_open_f14', 'panels': have('p1_open_f14', H2_OPEN)},
    {'title': 'H2: 内側を使うか、間隔をどれだけ取るか（P1、内側のトラックなし、原寸）',
     'note': f"既定 {hsum('p1_open_f14', 'I')}、外側だけ {hsum('p1_open_f14', 'I_in0')}、間隔 6 px {hsum('p1_open_f14', 'I_clear6')}。間隔は、ラベル同士にも、目盛りやトラックにも同じ値を使います。",
     'case': 'p1_open_f14', 'panels': have('p1_open_f14', H2_SET), 'box': dense_box('p1_open_f14', 640, 460, 'G')},
    {'title': 'H2: 簡単な地図では（λ 14 pt、地図全体）',
     'note': f"線の合計は今 8,171 px、H {hsum('lambda_f14', 'H')}、H2 {hsum('lambda_f14', 'I')}。",
     'case': 'lambda_f14', 'panels': have('lambda_f14', H2_PLAIN)},
    {'title': 'H2: 小さいプラスミドでは（pBR322 と制限酵素サイト 33 個、地図全体）',
     'note': f"線の合計は今 3,050 px、H {hsum('pbr322_f14', 'H')}、H2 {hsum('pbr322_f14', 'I')}。",
     'case': 'pbr322_f14', 'panels': have('pbr322_f14', H2_PLAIN)},
    {'title': 'H2 の弱点: 片側に偏ったとても混んだ記録（T5 の既定 8 pt、194 個）',
     'note': f"交差は今 2,085 組、H と H2 は 0。線の合計は今 67,238 px、H {hsum('t5_f8', 'H')}、H2 {hsum('t5_f8', 'I')}。"
             f"<b>キャンバスは H {canvas('t5_f8', 'H')}、H2 {canvas('t5_f8', 'I')} で、H2 のほうが縦に伸びます</b>。"
             'T5 は右半分に 140 個、左半分に 55 個と偏っていて、H は 12 時・6 時をまたいで左の列へ逃がしますが、H2 は折れを 110° 以上に保つためほとんど逃がせません。',
     'case': 't5_f8', 'panels': have('t5_f8', H2_PLAIN)},
    {'title': 'H: 中央が空いた地図で（P1、内側のトラックなし、原寸）',
     'note': f"H は 111 個を全部置きます。既定は {P1O['proto']['notes']['inside']} 個が内側で線の合計 {hsum('p1_open_f14', 'H')}、外側だけにすると {hsum('p1_open_f14', 'H_in0')}、"
             f"間隔 6 px では {hsum('p1_open_f14', 'H_clear6')}。G は 102 個を置いて 9 個を省略（4,622 px）、D は 12,526 px、今は 12,731 px（重なり 20 組）。",
     'case': 'p1_open_f14', 'panels': have('p1_open_f14', H_SET), 'box': dense_box('p1_open_f14', 640, 460, 'G')},
    {'title': 'H: 同じ比較を地図全体で（P1、内側のトラックなし）',
     'note': 'キャンバスは今 1706 × 990、D 1624 × 1179、G 1551 × 904、H 1574 × 1081。',
     'case': 'p1_open_f14', 'panels': have('p1_open_f14', H_PLAIN)},
    {'title': 'H: 内側に GC トラックがある既定の図（P1 14 pt、地図全体）',
     'note': f"内側はふさがっているので、H は全部外側です。線の合計は今 12,731 px（重なり 20 組）、D 12,526 px、H {hsum('p1_f14', 'H')}。"
             'キャンバスは今 1706 × 990、D 1624 × 1179、H 1574 × 1224。G は 26 個を省略します。',
     'case': 'p1_f14', 'panels': have('p1_f14', H_PLAIN)},
    {'title': 'H: 簡単な地図では（λ 14 pt、地図全体）',
     'note': f"線の合計は今 8,171 px、D 8,256 px、H {hsum('lambda_f14', 'H')}。G は 8 個を省略します。",
     'case': 'lambda_f14', 'panels': have('lambda_f14', H_PLAIN)},
    {'title': 'H: 小さいプラスミドでは（pBR322 と制限酵素サイト 33 個、地図全体）',
     'note': f"線の合計は今 3,050 px、D 3,936 px、H {hsum('pbr322_f14', 'H')}。G は 1 個を省略します（970 px）。",
     'case': 'pbr322_f14', 'panels': have('pbr322_f14', H_PLAIN)},
    {'title': 'H: λ で内側を使うと（内側のトラックなし、原寸）',
     'note': f"71 個を全部置きます。既定は 10 個が内側で {hsum('lambda_open_f14', 'H')}、外側だけでは {hsum('lambda_open_f14', 'H_in0')}、間隔 6 px では {hsum('lambda_open_f14', 'H_clear6')}。"
             'G は 3 個を省略します（2,973 px）。',
     'case': 'lambda_open_f14', 'panels': have('lambda_open_f14', H_SET), 'box': dense_box('lambda_open_f14', 640, 460, 'G')},
    {'title': 'リング側へ戻るラベル（P1 14 pt の 12 時付近、原寸）',
     'note': 'G では「Phd-like antitoxin」「Doc-like toxin」が、線のつく左端から 12 時の方向へ伸びて、先がリングに近づきます'
             '（内側の「Coi」「ImcB」も同じ）。今の配置でも、269° に張り付いた 3 個が同じ形になります。D と H は、線がつく端からリングと反対の向きへ文字を伸ばすので 0 個です。',
     'case': 'p1_f14', 'panels': have('p1_f14', BACK_SET), 'box': [-340, -600, 240, -290]},
    {'title': 'F: 線の合計を最小にする配置（P1 14 pt、地図全体）',
     'note': '線の合計: 今 12,731 px（重なり 20 組）、F 27,594 px（楕円 1.25 倍）、間隔 6 px では 35,699 px、楕円に縛らない版 11,953 px'
             '（交差 1 組と文字を横切る線 1 本が残る）、D 12,526 px。',
     'case': 'p1_f14', 'panels': have('p1_f14', F_P1)},
    {'title': 'F: 目的関数を変えても（λ 14 pt、地図全体）',
     'note': '合計を最小にする版と、最大値を最小にする版はほぼ同じ絵です（いちばん長い線が 235.0 px → 233.6 px）。',
     'case': 'lambda_f14', 'panels': have('lambda_f14', F_LAMBDA)},
    {'title': 'F: T5 の既定 8 pt（今は交差 2,085 組）',
     'note': 'F は交差 0 にできますが、楕円が 1.5 倍になります（1784 × 1358、合計 70,473 px）。D は 40,309 px。',
     'case': 't5_f8', 'panels': have('t5_f8', [F_P1[0], F_P1[1], F_P1[3], F_P1[4]])},
] + cfg['variants']

cfg['footer'] = ('図はすべて調査用の試作で描いたものです（ローカルのブランチ <code>study/label-placement</code>）。'
                 '数字は各図の配置データから計算しています。線の合計は feature からラベルまでの折れ線全体の長さです。'
                 '省いたラベルがある候補では、置いたラベルの分だけの合計です。')

# ---- current Gallery and tutorial drawings under H (Owner 2026-10-12) --------------------------------------------
DOCS_DIR = str(Path(__file__).resolve().parents[1] / 'work' / 'docs')
if DOCS_DIR not in cfg['render_dirs']:
    cfg['render_dirs'].append(DOCS_DIR)


def pair_note(case):
    """Facts of one Current / H / H2 group, read from the three renders."""
    a, h, b = load_metrics(case, 'cur'), load_metrics(case, 'H'), load_metrics(case, 'I')
    for key in ('text_overlap_pairs', 'leader_crossing_pairs', 'leader_through_text', 'labels_on_objects', 'labels_near_objects', 'leader_bends_lt110'):
        assert not b[key], (case, key, b[key])
    return (f"外側のラベル {b['labels_external']} 個。線の合計は今 {round(a['leader_len_sum']):,} px、H {round(h['leader_len_sum']):,} px、H2 {round(b['leader_len_sum']):,} px。"
            f"キャンバスは今 {a['canvas_px'][0]} × {a['canvas_px'][1]}、H2 {b['canvas_px'][0]} × {b['canvas_px'][1]}。")


def leader_px(case, design, name, folder):
    """Whole leader length of one outside label, from the render's label list."""
    import math
    import lp_metrics
    data = json.loads((Path(__file__).resolve().parents[1] / 'work' / folder / f'{case}__{design}.labels.json').read_text())
    got = [sum(math.dist(pts[i], pts[i + 1]) for i in range(len(pts) - 1))
           for lab in data['labels'] if lab['label_text'] == name and not lab.get('is_embedded')
           for pts in [lp_metrics.label_leader(lab)]]
    assert len(got) == 1, (case, design, name, got)
    return round(got[0])


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
T01, T03 = 'doc_tcli01_human_mitochondrion', 'doc_tcli03_mitochondrial_features_highlighted'


def three(case, name):
    return f"{name} 今 {leader_px(case, 'cur', name, 'docs')} px、H {leader_px(case, 'H', name, 'docs')} px、H2 {leader_px(case, 'I', name, 'I')} px"


NOW_H = [V('cur', 'Current', '今の図'), V('H', 'H', '前の版（H）'), V('I', 'H2', '新しい配置（H2）')]
TOB = [V('cur', 'Current', '今の図（radial）'), V('hz', 'Current 水平', '今の配置で radial を外すと'), V('H', 'H', '前の版（H、外側だけ）'),
       V('I', 'H2', '新しい配置（H2、内側ラベルも）')]
_tob = {k: load_metrics('gal_tobacco_chloroplast', k) for k in ('cur', 'hz', 'H', 'I')}
_tob_notes = _tob['I']['proto']['notes']
_t01 = {k: load_metrics(T01, k) for k in ('cur', 'H', 'I')}
cfg['docs_title'] = 'ご指摘の 4 点と、今のギャラリーとチュートリアルの図を新しい配置（H2）で描いた結果'
cfg['docs_intro'] = ('最初の 5 組は、H へのご指摘の箇所を原寸で切り出したものです（今の配置、前の版 H、直した版 H2）。'
                     'そのあとに、dev（2aa6fbf0）のギャラリーとチュートリアルに載っている Circular の図を、掲載されているコマンドのまま 3 通りで描いた全体図を並べました。'
                     'どの組も同じデータ・同じ設定で、違うのはラベルの配置だけです。図をクリックすると原寸になります。'
                     'ラベルのある図は 6 種類でした（ギャラリー 3、チュートリアル 3。同じコマンドの図は 1 つにまとめています）。')
cfg['docs_unchanged'] = [
    '<b>変わらない図</b>（外側ラベルがないので、H2 でも SVG が 1 バイトも変わりません）: チュートリアル「Highlight selected features on the human mitochondrial map」の最初の図、'
    '「Plot read depth, GC content, and GC skew on a genome map」の 2 枚。ギャラリー「Every replicon of a multi-replicon bacterial genome on one canvas」もラベルなしです'
    '（入力の GenBank がリポジトリにないため、これだけ描いていません）。',
    '<b>対象外</b>: Linear の図（ギャラリー 6、チュートリアル 4）は今回の変更の範囲外です。チュートリアルの GUI スクリーンショットは撮り直していません。'
    '中の図は、同じページの CLI の図と同じ設定のものを下に載せています。',
    '<b>内側ラベル（<code>--labels both</code>）</b>: H2 は内側ラベルも同じ規則で置きます。前の版 H は外側ラベルだけを置き、内側は今の配置のままでした。',
]
cfg['docs'] = [
    {'title': '指摘 1（線が横へ伸びる）: ヒト mtDNA のチュートリアル図の 6 時付近を原寸で',
     'note': 'H2 はラベルを円のすぐ外、feature の真下に置きます。線の長さ: ' + '。'.join(three(T01, n) for n in ('ATP8', 'tRNA-Lys', 'tRNA-Asp')) + '。',
     'case': T01, 'panels': NOW_H, 'box': [-300, 280, 300, 490]},
    {'title': '指摘 1 の続き: 同じ図の 12 時付近を原寸で',
     'note': three(T01, 'tRNA-Thr') + '。' + three(T01, 'tRNA-Pro') + '。tRNA-Pro は短くなっていません。feature が内側のレーンにあり、線 77 px のうち 62 px がリングの中を通るためです'
             '（リングの外に出てからは 15 px）。',
     'case': T01, 'panels': NOW_H, 'box': [-330, -490, 60, -300]},
    {'title': '指摘 2（目盛りとの重なり）: チュートリアル「Highlight」の完成図の 6 時付近を原寸で',
     'note': 'H は「ATP8」を目盛りの文字「9 kbp」の上に置いていました。H2 は「9 kbp」と「8 kbp」の間に置き、どちらからも 3 px 以上離します。' + three(T03, 'ATP8') + '。',
     'case': T03, 'panels': NOW_H, 'box': [-260, 330, 160, 480]},
    {'title': '指摘 3（線の折れ）: ヒト mtDNA のチュートリアル図の 3 時〜5 時を原寸で',
     'note': f"この図で折れが 110° 未満の線は、今 {_t01['cur']['leader_bends_lt110']} 本、H {_t01['H']['leader_bends_lt110']} 本、H2 {_t01['I']['leader_bends_lt110']} 本（24 本のうち）。"
             'H2 の線は、feature からリングの縁まで放射方向に出て、そこで 1 回だけ折れます（110° 以上）。',
     'case': T01, 'panels': NOW_H, 'box': [240, -10, 540, 330]},
    {'title': '指摘 4（内側ラベルの線）: 葉緑体の図の左上を原寸で',
     'note': '今の配置は、内側ラベルの feature 側の線が 47 px で固定です（内側ラベルのオフセットを変えても同じでした）。H2 は平均 12 px で、線全体は平均 83 px → 21 px。'
             f"内側に頼んだ {_tob_notes['asked_inside']} 個のうち {_tob_notes['inside']} 個を内側に置き、入りきらない {len(_tob_notes['overflow'])} 個（{'、'.join(_tob_notes['overflow'])}）は外側に出しています。"
             '今の配置（水平）では内側ラベル 3 個（5S、4.5S、23S ribosomal RNA）がほかの物に重なりますが、H2 は 0 個です。',
     'case': 'gal_tobacco_chloroplast', 'panels': [TOB[1], TOB[3]], 'box': [-430, -290, -110, 70]},
    {'title': 'ギャラリー: A first circular genome map (human mitochondrial genome)',
     'note': pair_note('gal_hmmtdna_basic_circular'),
     'case': 'gal_hmmtdna_basic_circular', 'panels': NOW_H},
    {'title': 'ギャラリー: Strand composition of the human mitochondrial genome (AT skew)',
     'note': pair_note('gal_hmmtdna_atskew') + 'ラベルは 20 pt です。',
     'case': 'gal_hmmtdna_atskew', 'panels': NOW_H},
    {'title': 'ギャラリー: Quadripartite structure of a chloroplast genome (Nicotiana tabacum) / チュートリアル: Make an annotated map of the tobacco chloroplast genome',
     'note': 'この 2 つは同じコマンドで、同じ図です。今の図は <code>--labels both --label_placement radial</code>（文字を放射状に置く）です。'
             'H2 は水平の文字なので、比較のために「今の配置で radial を外した図」も並べました。'
             f"キャンバスは {_tob['cur']['canvas_px'][0]} × {_tob['cur']['canvas_px'][1]}（radial）、"
             f"{_tob['hz']['canvas_px'][0]} × {_tob['hz']['canvas_px'][1]}（水平）、{_tob['I']['canvas_px'][0]} × {_tob['I']['canvas_px'][1]}（H2）。"
             f"線の合計は水平 {round(_tob['hz']['leader_len_sum']):,} px、H2 {round(_tob['I']['leader_len_sum']):,} px（内側ラベルを含む）。",
     'case': 'gal_tobacco_chloroplast', 'panels': TOB},
    {'title': '同じ葉緑体の図の右上を原寸で',
     'note': '外側ラベルと内側ラベルの違いを原寸で比べます。',
     'case': 'gal_tobacco_chloroplast', 'panels': TOB, 'box': [40, -470, 560, -90]},
    {'title': 'チュートリアル: Draw a labeled circular map of the human mitochondrial genome / Save an interactive mitochondrial figure and reopen it later',
     'note': 'この 2 ページの図は同じコマンドです。' + pair_note(T01),
     'case': T01, 'panels': NOW_H},
    {'title': 'チュートリアル: Show how three animal mitochondrial genomes match the human mitochondrial genome',
     'note': pair_note('doc_tcli09_precomputed_circular_rings') + '外側に比較リングが 3 本あり、ラベルはその外に付きます。',
     'case': 'doc_tcli09_precomputed_circular_rings', 'panels': NOW_H},
    {'title': 'チュートリアル: Highlight selected features on the human mitochondrial map（完成図）',
     'note': pair_note(T03) + 'ほとんどのラベルは feature の中に描かれていて、外に出るのは 2 個だけです。目盛りと目盛りの文字はリングの外側にあります。',
     'case': T03, 'panels': NOW_H},
]

import d11_texts  # LP-D11: the Owner's three points on H2
d11_texts.apply(cfg, V)
import d15_texts  # LP-D15: where a leader ends, examples for `back` and for L
d15_texts.apply(cfg, V)

PAGE.write_text(json.dumps(cfg, ensure_ascii=False, indent=1))
print('texts written:', sum(bool(d['desc']) for d in cfg['designs']), 'design descriptions,',
      sum(bool(c.get('notes')) for c in cfg['cases']), 'cases with notes')
