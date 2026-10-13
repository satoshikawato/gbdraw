import json,math,sys
sys.path.insert(0,'/home/kawato/gbdraw-baselines/label-placement-20261011/work/J')
import lop
names=('Trp','Ala','Asn','Cys','Tyr','Leu','Ser','His','Ile','Gln','Met')
for l in json.load(open(sys.argv[1]))['labels']:
    if l.get('leader_points') and not l.get('is_embedded') and l['label_text'][5:] in names and l['label_text'].startswith('tRNA-') and abs(l['leader_points'][-1][0])>200 and abs(l['leader_points'][-1][1])<300:
        p=l['leader_points']; b=l['aabb_local']
        print(f"{l['label_text']:9s} y0={b[1]:6.1f} turn={lop.turns(p):6.1f} r={math.hypot(*p[-1]):6.1f} last={math.dist(p[-2],p[-1]):5.1f} {l.get('port')}")
