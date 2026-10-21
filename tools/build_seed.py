#!/usr/bin/env python3
"""Build the review-database seed of the translation review page.

Usage: build_seed.py [--versions VERSIONS_JSON] SNAPSHOT_DIR OUT_DIR ROWS_DIR [ROWS_DIR ...]
  SNAPSHOT_DIR  <dir>/rev-<loc>/<chunk>.json, as ArtifactData `list` with `out_dir`
                writes them: {"items": {<key id>: {t, s, c, by, at}}}; "__*" fields
                added by the tool are dropped.
  OUT_DIR       receives rev-<loc>/<chunk>.json (whole documents) and batch-<loc>.json
                (ArtifactData `batch` entries; file_path relative to OUT_DIR).
  ROWS_DIR      rows-<loc>.json = {<chunk>: {<key id>: row}}, merged over the snapshot.
                Several folders are applied in the order given; a later row replaces an
                earlier one for the same key. An item whose "by" is non-empty was edited
                in the page and is never overwritten; it is reported instead.
  --versions    {"rev-<loc>": {<chunk>: version}}: pins each write with if_version when
                the target database already holds the document. Omit for a new database.
"""
import json
import sys
from pathlib import Path

LOCALES = ['ja', 'zh-Hans', 'zh-Hant']


def main() -> None:
    args = sys.argv[1:]
    versions = {}
    if args[:1] == ['--versions'] and len(args) > 1:
        versions = json.loads(Path(args[1]).read_text(encoding='utf-8'))
        args = args[2:]
    if len(args) < 3:
        sys.exit(__doc__)
    snap, out, *rows_dirs = (Path(a) for a in args)
    total = 0
    for loc in LOCALES:
        coll = f'rev-{loc}'
        layers = [json.loads(f.read_text(encoding='utf-8')) for f in (d / f'rows-{loc}.json' for d in rows_dirs) if f.exists()]
        (out / coll).mkdir(parents=True, exist_ok=True)
        writes, kept = [], []
        chunks = sorted({p.stem for p in (snap / coll).glob('*.json')}.union(*layers))
        for chunk in chunks:
            src = snap / coll / f'{chunk}.json'
            items = json.loads(src.read_text(encoding='utf-8')).get('items', {}) if src.exists() else {}
            for rows in layers:
                for key_id, row in rows.get(chunk, {}).items():
                    if items.get(key_id, {}).get('by'):
                        kept.append((chunk, key_id))
                        continue
                    items[key_id] = row
            doc = out / coll / f'{chunk}.json'
            doc.write_text(json.dumps({'items': items}, ensure_ascii=False), encoding='utf-8')
            write = {'op': 'set', 'collection': coll, 'doc_id': chunk, 'file_path': f'{coll}/{chunk}.json'}
            if chunk in versions.get(coll, {}):
                write['if_version'] = versions[coll][chunk]
            writes.append(write)
            total += len(items)
        (out / f'batch-{loc}.json').write_text(json.dumps(writes, indent=1), encoding='utf-8')
        print(f'{coll}: {len(writes)} documents; page edits kept: {kept or "none"}')
    print(f'items: {total}')


if __name__ == '__main__':
    main()
