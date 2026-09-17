"""Small, deliberately direct reference for unit fields and ordered lookups.

Groups start in input order, then each group is independently sorted. Alias
ambiguity uses complete target sets, independent of the production reduction.
No production resolver, rank, alias, or grouping helper is used.
"""
from __future__ import annotations

import logging

from gbdraw.analysis.collinearity_units import CollinearityUnit, CollinearityUnitIndex
from gbdraw.exceptions import ValidationError


def locus(protein):
    xrefs = [str(x).strip() for x in (protein.db_xref or ())]
    gene_ids = [x for x in xrefs if x.startswith('GeneID:') and x != 'GeneID:']
    choices = [protein.gene_parent_id, protein.locus_tag, protein.gene_id, *gene_ids]
    return next((str(x).strip() for x in choices if str(x or '').strip()), None)


def position(protein):
    return int(protein.start), int(protein.end), int(protein.feature_index), str(protein.protein_id)


def build(extraction, *, records=None, mode='auto'):
    normalized = str(mode).strip().lower() or 'auto'
    if normalized not in ('auto', 'cds', 'locus'):
        raise ValidationError('collinear_unit_mode must be one of: auto, cds, locus')
    ids = ([str(r.id) for r in records] if records is not None else
           [str(ps[0].record_id) if ps else '' for ps in extraction.proteins_by_record])
    result = CollinearityUnitIndex([], {}, {}, [], [])
    fallback, collapsed = [], []
    for ri, proteins in enumerate(extraction.proteins_by_record):
        if normalized == 'locus':
            missing = [p for p in sorted(proteins, key=position) if locus(p) is None]
            if missing:
                examples = ', '.join(f'{p.record_id}:{p.protein_id}' for p in missing[:5])
                raise ValidationError("collinear_unit_mode='locus' requires stable locus identifiers for all CDS proteins; "
                                      f'missing examples: {examples}')
        groups = []
        for p in proteins:
            key = locus(p)
            matches = [g for g in groups if g[0] == 'locus' and g[1] == key]
            if normalized != 'cds' and key is not None:
                if matches:
                    matches[0][2].append(p)
                else:
                    groups.append(('locus', key, [p]))
            else:
                groups.append(('cds', key, [p]))
        # Establish stable group insertion ties from the global position order.
        ordered = sorted(proteins, key=position)
        groups.sort(key=lambda g: min(next(i for i, p in enumerate(ordered) if p is member)
                                      for member in g[2]))
        groups.sort(key=lambda g: tuple(min(column) for column in zip(*(position(p) for p in g[2]))))
        record_units = []
        for order, (kind, key, members) in enumerate(groups):
            members = sorted(members, key=position)
            rep = sorted(members, key=lambda p: (-int(p.protein_length), -bool(p.source_protein_id),
                         -max(0, int(p.end)-int(p.start)), int(p.feature_index), str(p.protein_id)))[0]
            uid = f'gbd_r{ri+1:04d}_unit{len(result.unit_by_id)+1:06d}'
            name = str(key or rep.label or rep.protein_id)
            alias_values = [uid, rep.protein_id, rep.source_protein_id, rep.feature_svg_id, key, name]
            for p in members:
                alias_values.extend([p.protein_id, p.source_protein_id, p.feature_svg_id,
                                     p.locus_tag, p.gene_id, p.old_locus_tag, p.gene])
            aliases = tuple(dict.fromkeys(str(a).strip() for a in alias_values if str(a or '').strip()))
            known = set(p.strand for p in members) & {-1, 1}
            unit = CollinearityUnit(uid, kind, ri, ids[ri] if ri < len(ids) else rep.record_id,
                       order, rep.protein_id, str(rep.feature_svg_id or ''),
                       min(int(p.start) for p in members), max(int(p.end) for p in members),
                       next(iter(known)) if len(known) == 1 else None, key, name,
                       tuple(p.protein_id for p in members), aliases)
            record_units.append(unit)
            result.unit_by_id[uid] = unit
            for p in members:
                result.unit_by_protein_id[str(p.protein_id)] = unit
        result.units_by_record.append(record_units)
        targets = {}
        for u in record_units:
            for a in u.aliases:
                targets.setdefault(a, set()).add(u.unit_id)
        result.aliases_by_record.append({a: next(iter(us)) for a, us in targets.items() if len(us) == 1})
        result.ambiguous_aliases_by_record.append({a for a, us in targets.items() if len(us) > 1})
        if normalized == 'auto':
            label = ids[ri] if ri < len(ids) else str(ri+1)
            n = sum(u.unit_kind == 'cds' for u in record_units)
            c = sum(len(u.cds_members)-1 for u in record_units if u.unit_kind == 'locus')
            if n:
                fallback.append(f'{label}: {n}')
            if c:
                collapsed.append(f'{label}: {c}')
    if fallback:
        logging.getLogger(__name__).warning(
            'WARNING: collinear_unit_mode auto used CDS units for CDS proteins without stable locus IDs (%s).',
            ', '.join(fallback[:5]))
    if collapsed:
        logging.getLogger(__name__).info(
            'INFO: Collapsed multiple CDS proteins into locus collinearity units (%s).', ', '.join(collapsed[:5]))
    return result
