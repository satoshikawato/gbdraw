#!/usr/bin/env python3
"""Build the review-database seed of the translation review page.

Usage: build_seed.py SNAPSHOT_DIR ROWS_DIR OUT_DIR [VERSIONS_JSON]
  SNAPSHOT_DIR  <dir>/rev-<loc>/<chunk>.json, as ArtifactData `list` with `out_dir`
                writes them: {"items": {<key id>: {t, s, c, by, at}}}; "__*" fields
                added by the tool are dropped.
  ROWS_DIR      rows-<loc>.json = {<chunk>: {<key id>: row}}, merged over the snapshot.
                An item whose "by" is non-empty was edited in the page and is never
                overwritten; it is reported instead.
  OUT_DIR       receives rev-<loc>/<chunk>.json (whole documents) and batch-<loc>.json
                (ArtifactData `batch` entries; file_path relative to OUT_DIR).
  VERSIONS_JSON {"rev-<loc>": {<chunk>: version}}: pins each write with if_version when
                the target database already holds the document. Omit for a new database.
"""
import json
import sys
from pathlib import Path

LOCALES = ['ja', 'zh-Hans', 'zh-Hant']


def main() -> None:
    if len(sys.argv) not in (4, 5):
        sys.exit(__doc__)
    snap, rows_dir, out = (Path(a) for a in sys.argv[1:4])
    versions = json.loads(Path(sys.argv[4]).read_text(encoding='utf-8')) if len(sys.argv) == 5 else {}
    total = 0
    for loc in LOCALES:
        coll = f'rev-{loc}'
        rows_file = rows_dir / f'rows-{loc}.json'
        rows = json.loads(rows_file.read_text(encoding='utf-8')) if rows_file.exists() else {}
        (out / coll).mkdir(parents=True, exist_ok=True)
        writes, kept = [], []
        chunks = sorted({p.stem for p in (snap / coll).glob('*.json')} | set(rows))
        for chunk in chunks:
            src = snap / coll / f'{chunk}.json'
            items = json.loads(src.read_text(encoding='utf-8')).get('items', {}) if src.exists() else {}
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
