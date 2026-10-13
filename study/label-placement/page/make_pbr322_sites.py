"""Write page/inputs/pBR322_sites.gb: pBR322 (J01749.1) with its genes, origin, and restriction sites.

A plasmid map as in the Owner's references: the CDS and rep_origin of the NCBI record are kept, the other
features are dropped, and every enzyme of ENZYMES that cuts the plasmid once is added as a misc_feature over its
recognition site, whose /note is "<enzyme> <first base of the site>" (the numbering of the classic maps).
Run: /home/kawato/micromamba/bin/python3.13 -I make_pbr322_sites.py
"""
from pathlib import Path
from Bio import Restriction, SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature

# Enzymes of the classic pBR322 maps (New England Biolabs names); only the single cutters are drawn.
ENZYMES = """EcoRI ClaI HindIII EcoRV NheI BamHI SphI SalI EagI NruI BspMI BsmI StyI AvaI MscI BspEI PvuII BstZ17I
NdeI Tth111I AflIII AlwNI SapI PstI PvuI ScaI SspI AatII AseI BsaI AhdI XmnI BsgI PpuMI BsrBI DrdI BsmBI EarI
BstAPI PflMI BseRI XbaI KpnI XhoI SacI NcoI SmaI NotI BglII""".split()

HERE = Path('/home/kawato/gbdraw-baselines/label-placement-20261011/page')
rec = SeqIO.read(HERE / 'inputs-ncbi' / 'J01749.1.gb', 'genbank')
# Keep the full-length source only: the record has a second source for the pSC101-derived part, and the map
# title would be taken from it.
keep = [f for f in rec.features if f.type in ('CDS', 'rep_origin') or (f.type == 'source' and len(f) == len(rec))]
names = {'tetracycline resistance protein': 'tet', 'beta-lactamase': 'bla', 'ROP protein': 'rop'}
for f in keep:
    if f.type == 'CDS':
        f.qualifiers['product'] = [names.get(f.qualifiers['product'][0], f.qualifiers['product'][0])]
    if f.type == 'rep_origin':
        f.qualifiers['note'] = ['ori']
seq = str(rec.seq).upper()
n = len(seq)
sites = []
for name in ENZYMES:
    enzyme = getattr(Restriction, name)
    if len(enzyme.search(rec.seq, linear=False)) != 1:
        continue
    starts = sorted({m.start() % n for m in enzyme.compsite.finditer(seq + seq[:enzyme.size - 1])})
    if len(starts) == 1:
        sites.append((starts[0], name, enzyme.size))
for start, name, size in sorted(sites):
    keep.append(SeqFeature(FeatureLocation(start, min(start + size, n), strand=1), type='misc_feature',
                           qualifiers={'note': [f'{name} {start + 1}']}))
rec.features = sorted(keep, key=lambda f: (int(f.location.start), f.type != 'source'))
out = HERE / 'inputs' / 'pBR322_sites.gb'
SeqIO.write(rec, out, 'genbank')
print(out, len(sites), 'sites:', ', '.join(f'{e} {s + 1}' for s, e, _ in sorted(sites)))
