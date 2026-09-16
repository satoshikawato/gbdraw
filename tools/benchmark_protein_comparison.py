#!/usr/bin/env python3
"""Reproduce protein-comparison baselines without changing production code.

Run each source tree in a fresh process. Timing (one warmup, three samples),
profiling/counters, and tracemalloc are separate invocations. Inputs and complete
ordered stage results are hashed; --artifacts also saves the semantic oracles.
The benchmark records current behavior, including errors, not Product authority.
"""
from __future__ import annotations

import argparse
import ast
import inspect
import base64
from collections import Counter
from contextlib import ExitStack
import cProfile
from dataclasses import fields, is_dataclass
import gc
import gzip
import hashlib
import importlib
import importlib.metadata
import io
import json
import math
from pathlib import Path
import platform
import pstats
import random
import resource
import statistics
import subprocess
import sys
import tempfile
import time
import tracemalloc
from unittest.mock import patch

NAME = "gbdraw-protein-comparison"
SEED = 20260915
POLICY = {"warmups": 1, "samples": 3, "regressionPct": 10.0,
          "maxNoisePct": 5.0, "noise": "100 * MAD / median",
          "noisyDecision": "inconclusive; no automatic repeat",
          "semanticGate": "exact ordered stage digest, including types and errors"}
GALLERIES = ("hepatoplasmataceae_collinear", "hepatoplasmataceae_orthogroup",
             "vibrio-harveyi-group-collinear")
CASES = ("hsp-edges", "hsp-1", "hsp-1000", "hsp-many", "dense-24", "support-edges",
         "sparse-2", "sparse-200", "sparse-400", "sparse-800",
         "giant-200", "giant-800", "unrelated-200", "unrelated-800",
         "cache-49", "cache-64", "cache-81", "path-8", "path-12", "path-16", "path-24", "path-32", "path-56",
         "merge-300", "merge-600", "merge-1200", "merge-edges", "manifest",
         "units-multicds", "gallery-collinear", "gallery-collinear-off", "gallery-collinear-vibrio",
         "gallery-collinear-vibrio-off", "gallery-orthogroup", "gallery-orthogroup-unbounded", "render-gallery")


def digest(data):
    return hashlib.sha256(data).hexdigest()


def json_bytes(value):
    return json.dumps(value, ensure_ascii=False, allow_nan=False,
                      separators=(",", ":")).encode()


def canonical(value):
    """Preserve ordering, tuple/list distinctions, dtypes, and nonfinite values."""
    import numpy as np
    import pandas as pd
    if isinstance(value, pd.DataFrame):
        return {"type": "DataFrame", "columns": list(value.columns),
                "dtypes": [str(x) for x in value.dtypes],
                "index": canonical(value.index.tolist()),
                "rows": canonical(list(value.itertuples(index=False, name=None)))}
    if isinstance(value, np.generic):
        return canonical(value.item())
    if value is pd.NA:
        return {"type": "pd.NA"}
    if type(value).__name__ == "OrthologPathCollection":
        return {"type": "OrthologPathCollection", "value": canonical(value.to_payload())}
    if is_dataclass(value):
        return {"type": type(value).__name__,
                "fields": [[f.name, canonical(getattr(value, f.name))] for f in fields(value)]}
    if isinstance(value, dict):
        return {"type": "dict", "items": [[canonical(k), canonical(v)] for k, v in value.items()]}
    if isinstance(value, (tuple, list)):
        return {"type": type(value).__name__, "items": [canonical(x) for x in value]}
    if isinstance(value, (set, frozenset)):
        items = [canonical(x) for x in value]
        return {"type": type(value).__name__, "items": sorted(items, key=json_bytes)}
    if isinstance(value, bytes):
        return {"type": "bytes", "base64": base64.b64encode(value).decode()}
    if isinstance(value, float) and not math.isfinite(value):
        return {"type": "float", "value": str(value)}
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    raise TypeError(f"Unsupported semantic value: {type(value).__name__}")


def guard_imports(root):
    for name, module in list(sys.modules.items()):
        if name == "gbdraw" or name.startswith("gbdraw."):
            filename = getattr(module, "__file__", None)
            if not filename or not Path(filename).resolve().is_relative_to(root / "gbdraw"):
                raise ValueError(f"foreign gbdraw import: {name} from {filename}; use a fresh process")


def load_source(root):
    root = root.resolve()
    if not (root / "gbdraw/analysis/protein_colinearity.py").is_file():
        raise ValueError(f"source root does not contain gbdraw: {root}")
    guard_imports(root)
    sys.path.insert(0, str(root))
    pc = importlib.import_module("gbdraw.analysis.protein_colinearity")
    cc = importlib.import_module("gbdraw.analysis.collinearity")
    guard_imports(root)
    return pc, cc


def source_info(root):
    def git(*args):
        return subprocess.check_output(["git", "-C", str(root), *args], text=True).strip()
    paths = git("ls-files", "--cached", "--others", "--exclude-standard",
                "gbdraw", "tools", "tests").splitlines()
    changed = git("status", "--porcelain", "--untracked-files=all")
    branch = git("branch", "--show-current")
    return {"root": str(root), "head": git("rev-parse", "HEAD"),
            "branch": branch,
            "upstream": git("for-each-ref", "--format=%(upstream:short)", "refs/heads/" + branch) if branch else "",
            "status": changed, "diffSha256": digest(subprocess.check_output(
                ["git", "-C", str(root), "diff", "--binary", "HEAD"])),
            "filesSha256": {p: digest((root / p).read_bytes()) for p in paths if (root / p).is_file()},
            "runnerSha256": digest(Path(__file__).read_bytes())}


def helpers(root):
    text = (root / "gbdraw/web/js/app/python-helpers.js").read_text()
    if not text.startswith("export const PYTHON_HELPERS = `") or text.count("`") != 2:
        raise ValueError("Python helper embedding changed; review the extraction boundary")
    namespace = {}
    exec(compile(text.split("`", 1)[1].rsplit("`", 1)[0],
                 str(root / "gbdraw/web/js/app/python-helpers.js"), "exec"), namespace)
    return namespace


def hit(q, s, **kw):
    return dict(query=q, subject=s, identity=90.0, alignment_length=100,
                mismatches=0, gap_opens=0, qstart=1, qend=100, sstart=1,
                send=100, evalue=1e-30, bitscore=200.0) | kw


def protein(pc, pid, record, index=0, length=100):
    return pc.CdsProtein(protein_id=pid, record_index=record, feature_index=index,
                         record_id=f"record_{record}", start=index * 330,
                         end=index * 330 + length * 3, strand=1, label=pid,
                         protein_length=length, sequence="M" * length)


def frame(rows):
    import pandas as pd
    from gbdraw.io.comparisons import COMPARISON_COLUMNS
    return pd.DataFrame(rows, columns=COMPARISON_COLUMNS)


def safe_call(fn):
    try:
        return fn()
    except (ValueError, TypeError, OverflowError, KeyError) as exc:
        return {"observedException": type(exc).__name__, "message": str(exc)}


def synthetic(pc, name, seed):
    size = int(name.rsplit("-", 1)[1]) if name.rsplit("-", 1)[1].isdigit() else 0
    rng = random.Random(seed)
    if name.startswith("hsp-"):
        pm = {p: protein(pc, p, i, length=1000) for i, p in enumerate(("q", "s"))}
        if name == "hsp-many":
            size = 10
        if size:
            pm = {f"{p}{i}": protein(pc, f"{p}{i}", r, i, 1000)
                  for i in range(size) for r, p in enumerate(("q", "s"))}
            rows = [hit(f"q{i}", f"s{i}", qstart=(j % 5) * 150 + 1, qend=(j % 5) * 150 + 300,
                        sstart=(j % 5) * 150 + 1, send=(j % 5) * 150 + 300,
                        alignment_length=300, bitscore=500 - j % 3) for i in range(size)
                    for j in range(1000 if name == "hsp-many" else 3)]
            rng.shuffle(rows)
            tables = {"multi_hsp": frame(rows)}
        else:
            cases = {
                "overlap_duplicate_tie": [hit("q", "s", qend=300, send=300, alignment_length=300)] * 2
                    + [hit("q", "s", qstart=200, qend=500, sstart=200, send=500)],
                "disjoint_reverse_clamp": [hit("q", "s", qstart=300, qend=1, sstart=300, send=1),
                    hit("q", "s", qstart=800, qend=1200, sstart=-30, send=100)],
                "unknown_missing": [hit("unknown", "s"), hit(None, "s"), hit("q", None), hit("q", "s")],
                "bad_coordinates": [hit("q", "s", qstart="bad", qend=float("nan"),
                                         sstart=float("inf")), hit("q", "s", bitscore=150)],
                "pair_row_order": [hit("s", "q"), hit("q", "s"), hit("s", "q", identity=88)],
                "nan_score": [hit("q", "s", bitscore=float("nan"))],
                "infinite_score": [hit("q", "s", bitscore=float("inf"))],
                "infinite_length": [hit("q", "s", alignment_length=float("inf"))],
                "empty": [],
            }
            tables = {k: frame(v) for k, v in cases.items()}
        return pm, tables
    if name.startswith("path-"):
        if size >= 24 and not hasattr(pc, "OrthologPathCollection"):
            raise ValueError("R>=24 is compact-only; legacy enumeration is prohibited")
        pm = {f"p{i}": protein(pc, f"p{i}", i) for i in range(size)}
        return pm, {(i, j): frame([hit(f"p{i}", f"p{j}")])
                    for i in range(size) for j in range(size)}
    if name.startswith("sparse-"):
        pm = {f"{p}{i}": protein(pc, f"{p}{i}", r, i)
              for i in range(size) for r, p in enumerate(("q", "s", "u"))}
        return pm, {(0, 1): frame([hit(f"q{i}", f"s{i}") for i in range(size)]),
                    (1, 0): frame([hit(f"s{i}", f"q{i}") for i in range(size)])}
    if name.startswith("giant-"):
        # Two star cores: many members, only twenty sparse incoming attachments.
        # Two records bound the legacy path depth; no large all-record DAG.
        pm = {}
        rows = {(0, 1): [], (1, 0): [], (0, 0): []}
        for group in range(2):
            subject = f"s{group}"
            pm[subject] = protein(pc, subject, 1, group)
            for i in range(size):
                query = f"q{group}_{i}"
                pm[query] = protein(pc, query, 0, group * size + i)
                rows[0, 1].append(hit(query, subject))
                rows[1, 0].append(hit(subject, query))
            for i in range(10):
                pid = f"u{group}_{i}"
                pm[pid] = protein(pc, pid, 0, 2 * size + group * 10 + i)
                rows[0, 0].append(hit(f"q{group}_{i}", pid, bitscore=150))
        return pm, {pair: frame(values) for pair, values in rows.items()}
    if name.startswith("dense-"):
        pm = {f"{p}{i}": protein(pc, f"{p}{i}", r, i)
              for i in range(size) for r, p in enumerate(("q", "s", "u"))}
        rows = {(i, j): [] for i, j in ((0, 1), (1, 0), (2, 0), (2, 1))}
        for q in pm.values():
            for s in pm.values():
                if (q.record_index, s.record_index) not in rows:
                    continue
                rows[q.record_index, s.record_index].append(hit(q.protein_id, s.protein_id,
                    bitscore=250 if q.record_index != 2 and q.feature_index == s.feature_index else 100 + rng.randrange(20)))
        return pm, {k: frame(v) for k, v in rows.items()}
    if name == "support-edges":
        # Two cross-record cores, incoming-only support, record-local evidence,
        # domain-only evidence, best/second ties, and two competing unassigneds.
        pm = {p: protein(pc, p, r, i) for i, (p, r) in enumerate(
            [("q0", 0), ("s0", 1), ("q1", 0), ("s1", 1),
             ("u0", 2), ("u1", 2), ("local", 0), ("domain", 2)])}
        rows = [("q0", "s0", 200, 100), ("s0", "q0", 200, 100),
                ("q1", "s1", 200, 100), ("s1", "q1", 200, 100),
                ("q0", "u0", 150, 100), ("s0", "u0", 150, 100),
                ("u1", "q0", 160, 100), ("u1", "q1", 160, 100),
                ("q0", "local", 180, 100), ("local", "q0", 180, 100),
                ("domain", "q1", 200, 10), ("q1", "domain", 200, 10)]
        tables = {}
        for q, s, score, length in rows:
            tables.setdefault((pm[q].record_index, pm[s].record_index), []).append(
                hit(q, s, bitscore=score, alignment_length=length, qend=length, send=length))
        return pm, {k: frame(v) for k, v in tables.items()}
    raise ValueError(name)


def gallery_input(root, name, pc):
    from Bio import SeqIO
    path = root / f"gbdraw/web/gallery/sessions/{name}.gbdraw-session.json.gz"
    session = json.loads(gzip.decompress(path.read_bytes()))
    specs = session["renderRequest"]["records"]
    records, source_hashes = [], {}
    for spec in specs:
        rid = spec["source"]["resourceId"]
        data = base64.b64decode(session["resources"][rid]["data"], validate=True)
        source_hashes[rid] = digest(data)
        # These saved recipes use one already-resolved GenBank record per resource.
        # Reject different recipes instead of silently selecting their first record.
        if spec["selector"] is not None or spec.get("region") is not None:
            raise ValueError("Gallery recipe requires selector/region planning")
        parsed = list(SeqIO.parse(io.StringIO(data.decode()), "genbank"))
        if len(parsed) != 1:
            raise ValueError("Gallery source cardinality changed; use the record planner")
        records.extend(parsed)
    keys = [x["recordKey"] for x in specs]
    extraction = pc.extract_protein_identity_manifest(records, record_instance_keys=keys)
    manifest = pc.validate_protein_identity_manifest(session["proteinIdentityManifest"])
    rebuilt = extraction.identity_manifest.to_dict()
    for key in keys:
        saved = session["proteinIdentityManifest"]["recordInstances"][key]
        fresh = rebuilt["recordInstances"][key]
        for field in ("runtimeIds", "runtimeBindingHash"):
            if fresh[field] != saved[field]:
                raise ValueError(f"Gallery {name}: mismatched {key}.{field}")
        analysis = session["proteinIdentityManifest"]["recordAnalyses"][saved["recordAnalysisId"]]
        fresh_analysis = rebuilt["recordAnalyses"][fresh["recordAnalysisId"]]
        if analysis["proteinSetHash"] != fresh_analysis["proteinSetHash"]:
            raise ValueError(f"Gallery {name}: protein sequence identity changed")
    entries = session["losatCache"]["entries"]
    raw = {}
    for entry in entries:
        pc.validate_protein_raw_entry_references(entry, manifest)
        pair = (keys.index(entry["queryRecordInstanceKey"]), keys.index(entry["subjectRecordInstanceKey"]))
        if pair in raw:
            raise ValueError(f"Duplicate Gallery direction {pair}")
        raw[pair] = entry["text"]
    inventory = {"session": str(path.relative_to(root)), "sessionSha256": digest(path.read_bytes()),
                 "savedSessionVersion": session["version"], "savedRequestSchema": session["renderRequest"]["schema"],
                 "sourceSha256": source_hashes, "sourceCount": len(source_hashes),
                 "recordIds": [r.id for r in records], "recordKeys": keys,
                 "proteinCount": len(extraction.protein_map), "directionalTableCount": len(raw),
                 "rawSha256": [[list(k), digest(v.encode())] for k, v in raw.items()],
                 "rawRows": sum(len(v.splitlines()) for v in raw.values()),
                 "manifestSha256": digest(json_bytes(session["proteinIdentityManifest"])),
                 "runtimeAndProteinSequenceIdentityMatch": True,
                 "savedSearchSettings": [{k: v for k, v in e.items() if k not in ("text", "filename", "key")}
                                         for e in entries],
                 "savedLosatConfig": session["config"]["losat"],
                 "rawSearchTime": None, "historicalSourceInvocations": None}
    return session, extraction, records, raw, inventory


def manifest_case(root, pc):
    session, _, _, _, inventory = gallery_input(root, GALLERIES[2], pc)
    ns = helpers(root)
    payload = session["proteinIdentityManifest"]
    pairs = [{"queryRecordInstanceKey": e["queryRecordInstanceKey"],
              "subjectRecordInstanceKey": e["subjectRecordInstanceKey"],
              "expectedOptions": {k: e[k] for k in ("program", "args", "outfmt", "searchContext") if k in e}}
             for e in session["losatCache"]["entries"]]
    text, pair_text = json.dumps(payload), json.dumps(pairs)
    def batch():
        result = json.loads(ns["build_protein_losat_cache_keys_json"](text, pair_text))
        if "error" in result:
            raise ValueError(result["error"])
        return result["keys"]
    typed = pc.validate_protein_identity_manifest(payload)
    expected = [pc.build_protein_losat_cache_key(pc.build_protein_losat_pair_identity(
        typed, query_record_instance_key=p["queryRecordInstanceKey"],
        subject_record_instance_key=p["subjectRecordInstanceKey"]), **{
            "program": p["expectedOptions"]["program"], "outfmt": p["expectedOptions"]["outfmt"],
            "args": p["expectedOptions"]["args"], "search_context": p["expectedOptions"].get("searchContext")}) for p in pairs]
    if batch() != expected or expected != [e["key"] for e in session["losatCache"]["entries"]]:
        raise ValueError("batched / typed / saved key mismatch")
    by_direction = {(p["queryRecordInstanceKey"], p["subjectRecordInstanceKey"]): key
                    for p, key in zip(pairs, expected)}
    if any(q != s and by_direction.get((s, q)) == key for (q, s), key in by_direction.items()):
        raise ValueError("query/subject direction lost from cache identity")
    # Probe order, direction and searched-DB identity with the same actual helper.
    altered = json.loads(pair_text)
    altered[0]["expectedOptions"]["searchContext"] = digest(b"benchmark-different-database")
    changed = json.loads(ns["build_protein_losat_cache_keys_json"](text, json.dumps(altered)))["keys"]
    reversed_keys = json.loads(ns["build_protein_losat_cache_keys_json"](text, json.dumps(pairs[::-1])))["keys"]
    tampered = json.loads(text)
    instance = next(iter(tampered["recordInstances"].values()))
    first_id = next(iter(instance["runtimeIds"]))
    instance["runtimeIds"][first_id] = "h_invalid"
    rejected = "error" in json.loads(ns["build_protein_losat_cache_keys_json"](json.dumps(tampered), pair_text))
    if changed[0] == expected[0] or changed[1:] != expected[1:] or reversed_keys != expected[::-1] or not rejected:
        raise ValueError("batch identity invariant failed")
    inventory.update({"batchKeysEqualTypedAndSaved": True, "orderAndSearchContextVerified": True,
                      "tamperedRuntimeIdRejected": True, "manifestJsonBytes": len(text.encode()),
                      "pairJsonBytes": len(pair_text.encode()), "keys": len(pairs),
                      "warmInputBytesPerHelperCall": len(text.encode()) + len(pair_text.encode()),
                      "browserPayload": {"identityManifest": payload, "pairs": pairs},
                      "expectedKeys": expected})
    return inventory, {"batch_identity": batch}


def merge_case(cc, name):
    from dataclasses import replace
    def anchor(q, s, reverse=False):
        return cc.CollinearityAnchor(
            query_protein_id=f"q{q}", subject_protein_id=f"s{s}",
            query_record_index=0, subject_record_index=1, query_order=q, subject_order=s,
            query_start=q*10+1, query_end=q*10+9, subject_start=s*10+1, subject_end=s*10+9,
            identity=90., evalue=1e-30, bitscore=200., alignment_length=100,
            query_feature_svg_id=f"fq{q}", subject_feature_svg_id=f"fs{s}", source="benchmark",
            query_unit_id=f"qu{q}", subject_unit_id=f"su{s}", query_unit_kind="cds", subject_unit_kind="cds",
            query_locus_id=None, subject_locus_id=None, query_display_name=f"q{q}", subject_display_name=f"s{s}",
            query_strand=1, subject_strand=-1 if reverse else 1)
    if name == "merge-edges":
        anchors = [anchor(0, 0), anchor(1, 1), anchor(2, 3), anchor(4, 4), anchor(5, 5)]
        params = cc.LosslessCollinearityParameters(max_unit_gap=3, max_diagonal_drift=0, merge_orientation="strand")
        operations = {f"max_conflicts_{i}": lambda i=i: cc.cluster_lossless_collinearity_anchors(
            anchors, params=replace(params, max_conflicts=i)) for i in (0, 1)}
        operations["reverse"] = lambda: cc.cluster_lossless_collinearity_anchors(
            [anchor(i, 8-i, True) for i in range(9)], params=params)
        return {"anchors": canonical(anchors), "parameters": canonical(params)}, operations
    size = int(name.split("-")[1])
    anchors = [anchor(i, i if i % 3 != 2 else size * 2 + i) for i in range(size)]
    params = cc.LosslessCollinearityParameters(max_unit_gap=2, max_diagonal_drift=0, max_conflicts=0)
    return {"anchors": canonical(anchors), "parameters": canonical(params)}, {
        "cluster": lambda: cc.cluster_lossless_collinearity_anchors(anchors, params=params)}



def isolated_merge_operations(cc, operations):
    """Capture real merge inputs once; preparation is outside every measured run."""
    isolated = {}
    original = cc._merge_lossless_clusters
    for stage, fn in operations.items():
        if stage not in {"cluster", "post_search", "max_conflicts_0", "max_conflicts_1", "reverse"}:
            continue
        inputs = []
        def capture(blocks, *, anchors, params):
            inputs.append((tuple(blocks), tuple(anchors), params))
            return original(blocks, anchors=anchors, params=params)
        with patch.object(cc, "_merge_lossless_clusters", capture):
            fn()
        def replay(inputs=inputs):
            return tuple(cc._merge_lossless_clusters(blocks, anchors=anchors, params=params)
                         for blocks, anchors, params in inputs)
        isolated[stage + "_merge"] = replay
    return isolated

def cache_case(root, pc, name):
    size = int(name.split("-")[1])
    ns = helpers(root)
    raw = "q\ts\t90\t100\t0\t0\t1\t100\t1\t100\t1e-30\t200\n"
    def run():
        ns["_WEB_LOSATP_FILTERED_HIT_CACHE"].clear()
        ns["_WEB_LOSATP_CONVERTED_PAYLOAD_CACHE"].clear()
        ns["_WEB_LOSATP_CACHE_ORDER"].clear()
        def parsed():
            return pc.filter_protein_hits_by_thresholds(pc.parse_losatp_outfmt6(raw),
                bitscore=50, evalue=0.01, identity=0, alignment_length=0)
        for i in range(size):
            ns["_web_losatp_cache_set"]("filtered", str(i), parsed())
        ns["_web_losatp_cache_set"]("converted", "previous-result", "{}")
        hits = misses = 0
        outputs = []
        for i in range(size):
            cached = ns["_web_losatp_cache_get"]("filtered", str(i))
            if cached is None:
                misses += 1
                cached = parsed()
                ns["_web_losatp_cache_set"]("filtered", str(i), cached)
            else:
                hits += 1
            outputs.append(cached)
        retained = ns["_WEB_LOSATP_FILTERED_HIT_CACHE"]
        return {"semantic": outputs, "diagnostics": {"warmHits": hits, "warmMisses": misses, "coldParses": size, "warmParses": misses,
                "retainedFilteredTables": len(retained),
                "retainedDataFrameBytes": sum(int(df.memory_usage(index=True, deep=True).sum()) for df in retained.values()),
                "convertedPayloads": len(ns["_WEB_LOSATP_CONVERTED_PAYLOAD_CACHE"])}}
    return {"tableCount": size, "rawSha256": digest(raw.encode()),
            "boundary": "real Python helper LRU functions + native parse/filter; no raw search or JS derived hit"}, {"cache_cycle": run}


def build_case(root, pc, cc, name, seed, path_representation="graph", unit_stages=False):
    if name == "units-multicds":
        from dataclasses import replace
        from gbdraw.analysis.collinearity_units import build_collinearity_unit_index
        rng = random.Random(seed)
        rows = []
        for r in range(3):
            row = [replace(protein(pc, f"p{r}-{i}", r, i), locus_tag=f"L{i % 20}",
                           gene=f"shared{i % 3}", strand=-1 if i % 2 else 1)
                   for i in range(400)]
            rng.shuffle(row)
            rows.append(row)
        extraction = pc.ProteinExtractionResult(rows, {p.protein_id: p for row in rows for p in row})
        return {"seed": seed, "inputSha256": digest(json_bytes(canonical(extraction))),
                "proteinCount": 1200, "unitMode": "auto", "recordCount": 3,
                "boundary": "prepared extraction -> complete ordered unit index"}, {
                    "unit_index": lambda: build_collinearity_unit_index(extraction, mode="auto")}
    if name.startswith("unrelated-"):
        import pandas as pd
        # Membership-stage isolation: adding groups changes neither evidence nor
        # unassigned proteins. Anchors are an already prepared stage input.
        size = int(name.split("-")[1])
        pm = {f"m{i}": protein(pc, f"m{i}", 0, i) for i in range(size)}
        pm.update({f"u{i}": protein(pc, f"u{i}", 1, i) for i in range(200)})
        groups = {f"og_{i}": {f"m{i}"} for i in range(size)}
        mapping = {pid: group for group, members in groups.items() for pid in members}
        best = {(f"m{i}", f"u{i}"): next(pd.DataFrame([hit(f"m{i}", f"u{i}", normalized_score=1.0, min_coverage=1.0)]).itertuples(index=False))
                for i in range(200)}
        def support_scan():
            index = pc._index_core_support_evidence(best, mapping) if hasattr(pc, "_index_core_support_evidence") else None
            result = {}
            for i in range(200):
                pid = f"u{i}"
                candidates = []
                for group, evidence in (index.get(pid, {}) if index is not None else groups).items():
                    candidate = (pc._build_core_support_candidate(pid, group, evidence, {}, pm)
                                 if index is not None else
                                 pc._build_core_support_candidate(pid, group, sorted(evidence), best, {}, pm))
                    if candidate is not None:
                        candidates.append(candidate)
                result[pid] = sorted(candidates, key=pc._core_support_sort_key)
            return result
        return {"groups": size, "unassigned": 200, "evidence": len(best), "seed": seed,
                "inputSha256": digest(json_bytes(canonical((pm, groups, best)))),
                "boundary": "support candidate stage with precomputed membership; includes index construction"}, {"support_scan": support_scan}
    if name == "render-gallery":
        from gbdraw.api import load_session_document, materialize_session, render_session
        paths = [root / f"gbdraw/web/gallery/sessions/{g}.gbdraw-session.json.gz" for g in GALLERIES[:2]]
        documents = [load_session_document(p) for p in paths]
        def replay(document):
            with tempfile.TemporaryDirectory(prefix="gbdraw-comparison-render-") as directory:
                with materialize_session(document, output_directory=directory) as materialized:
                    result = render_session(materialized)
                    return [p.read_bytes() for p in result.output_paths]
        return {"sessionSha256": [digest(p.read_bytes()) for p in paths],
                "boundary": "saved comparison typed replay + materialization + SVG export; no raw search",
                "settings": "complete saved Gallery recipes, independent of the post-search benchmark settings"}, {
                    g: lambda document=d: replay(document) for g, d in zip(GALLERIES, documents)}
    if name == "manifest":
        return manifest_case(root, pc)
    if name.startswith("cache-"):
        return cache_case(root, pc, name)
    if name.startswith("merge-"):
        return merge_case(cc, name)
    if name.startswith("gallery-"):
        gallery = GALLERIES[2 if "-vibrio" in name else 0 if name.startswith("gallery-collinear") else 1]
        _, extraction, records, raw, inventory = gallery_input(root, gallery, pc)
        settings = {"bitscore": 50, "evalue": 0.01, "identity": 0, "alignment_length": 0}
        def parse():
            return {k: pc.parse_losatp_outfmt6(v) for k, v in raw.items()}
        parsed = parse()
        def filter_hits():
            return {k: pc.filter_protein_hits_by_thresholds(v, **settings) for k, v in parsed.items()}
        tables = filter_hits()
        member_max_hits = None if name.endswith("-unbounded") else 5
        inventory.update({"benchmarkThresholds": settings, "memberMaxHits": member_max_hits,
                          "inference": not name.endswith("-off"), "scope": "adjacent" if name.startswith("gallery-collinear") else "all",
                          "blockParameters": canonical(cc.LosslessCollinearityParameters()) if name.startswith("gallery-collinear") else None,
                          "unitMode": "auto", "edgeMode": "rbh", "maxRelatedEdges": 2,
                          "filteredRows": sum(len(x) for x in tables.values()),
                          "inputDataFrameDeepBytes": sum(int(x.memory_usage(deep=True).sum()) for x in tables.values())})
        def aggregate():
            return {k: pc._aggregate_hsps_by_protein_pair(v, extraction.protein_map) for k, v in tables.items()}
        if name.startswith("gallery-collinear"):
            def analyze():
                return cc.build_orthogroup_collinearity_blocks_from_hits(tables, extraction,
                            records=records, orthogroup_member_max_hits=member_max_hits, infer_orthogroups=not name.endswith("-off"), search_scope="adjacent")
        else:
            def analyze():
                return pc.select_rbh_orthogroup_edges_from_directional_hits(tables, extraction.protein_map,
                            record_count=len(records), orthogroup_member_max_hits=member_max_hits)
        from gbdraw.session_request_codec import encode_canonical_typed_resource
        from gbdraw.web_support.orthogroup_metadata import serialize_orthogroups_payload
        result = analyze()
        # Stage inputs are prepared once outside measurement. The post_search stage
        # itself includes its own normalization, inference, projection and blocks.
        operations = {"parse": parse, "filter": filter_hits, "hsp_aggregate": aggregate,
            "post_search": analyze,
            "metadata": lambda: serialize_orthogroups_payload(result.orthogroups, records=records),
            "display_tables": (lambda: cc.convert_collinearity_blocks_to_pair_comparisons(result, records=records))
                if name.startswith("gallery-collinear") else
                (lambda: [pc.convert_pair_protein_hits_to_genomic_links(
                    table, extraction.protein_map, extraction.protein_map, result.orthogroups)
                    for table in result.adjacent_display_edges_by_pair.values()]),
            "typed_serialization": lambda: encode_canonical_typed_resource("result", result if name.startswith("gallery-collinear") else result.orthogroups)}
        if unit_stages and name.startswith("gallery-collinear"):
            from gbdraw.analysis.collinearity_units import build_collinearity_unit_index
            operations["unit_index"] = lambda: build_collinearity_unit_index(extraction, records=records, mode="auto")
        return inventory, operations
    pm, tables = synthetic(pc, name, seed)
    inventory = {"seed": seed, "generator": name, "proteinCount": len(pm),
                 "inputSha256": digest(json_bytes(canonical((pm, tables)))),
                 "tableCount": len(tables), "rows": sum(len(x) for x in tables.values()),
                 "inputDataFrameDeepBytes": sum(int(x.memory_usage(deep=True).sum()) for x in tables.values()),
                 "memberMaxHits": None, "inference": True}
    if name.startswith("hsp-"):
        return inventory, {k: lambda v=v: safe_call(lambda: pc._aggregate_hsps_by_protein_pair(v, pm)) for k, v in tables.items()}
    if name.startswith("path-"):
        size = int(name.split("-")[1])
        inventory.update({"expectedPaths": 2**(size-2), "expectedEdges": size*(size-1)//2})
    def select():
        return pc.select_rbh_orthogroup_edges_from_directional_hits(tables, pm,
            record_count=max(p.record_index for p in pm.values())+1, orthogroup_member_max_hits=None,
            **({"path_representation": path_representation} if path_representation != "graph" else {}))
    if name.startswith("path-") and hasattr(pc, "OrthogroupGraphResult"):
        from gbdraw.session_request_codec import encode_canonical_typed_resource
        from gbdraw.web_support.orthogroup_metadata import serialize_orthogroups_payload
        result = select().orthogroups
        operations = {
            "selector": select,
            "metadata": lambda: serialize_orthogroups_payload(result),
            "typed_serialization": lambda: encode_canonical_typed_resource("result", result),
        }
        if path_representation == "graph":
            operations["path_graph"] = lambda: pc._build_ortholog_path_indexes(result.ortholog_edges_by_orthogroup_id, pm)
        return inventory, operations
    return inventory, {"selector": select}


def path_browser_inputs(root, pc, cc, names=("sparse-2", "path-24", "gallery-collinear", "gallery-collinear-off", "gallery-orthogroup")):
    """Use the same saved evidence and native production helper as the browser."""
    from dataclasses import asdict
    inputs = []
    for name in names:
        if name.startswith("gallery-"):
            gallery = GALLERIES[2 if "-vibrio" in name else 0 if name.startswith("gallery-collinear") else 1]
            _, extraction, records, raw, inventory = gallery_input(root, gallery, pc)
            pm = extraction.protein_map
            lengths = [len(r.seq) for r in records]
            ids = [r.id for r in records]
        else:
            pm, tables = synthetic(pc, name, SEED)
            raw = {pair: table.to_csv(sep="\t", header=False, index=False, lineterminator="\n")
                   for pair, table in tables.items()}
            count = max(p.record_index for p in pm.values()) + 1
            lengths = [max((p.end for p in pm.values() if p.record_index == i), default=300) for i in range(count)]
            ids = [f"record_{i}" for i in range(count)]
            inventory = {"generator": name, "inputSha256": digest(json_bytes(canonical((pm, tables))))}
        records_payload = [{"recordIndex": i, "recordId": ids[i],
                            "proteinCacheKey": f"{name}-record-{i}",
                            "proteinMap": {pid: asdict(p) for pid, p in pm.items() if p.record_index == i},
                            "viewTransform": {"length": lengths[i], "reverse": False}}
                           for i in range(len(ids))]
        offset, texts, pairs = 0, [], []
        for (q, t), text in raw.items():
            size = len(text.encode())
            pairs.append({"pairIndex": min(q, t), "queryIndex": q, "subjectIndex": t,
                          "cacheKey": f"{name}-raw-{q}-{t}", "displayPair": t == q+1,
                          "rawTsvOffset": offset, "rawTsvBytes": size})
            texts.append(text)
            offset += size
        parameters = {"mode": "collinear" if name.startswith("gallery-collinear") else "orthogroup",
                      "bitscore": 50, "evalue": "1e-2", "identity": 0, "alignmentLength": 0,
                      "orthogroupMemberMaxHits": 5 if name.startswith("gallery-") else None,
                      "collinearInferOrthogroups": not name.endswith("-off")}
        inputs.append({"name": name, "inventory": inventory, "parameters": parameters,
                       "pairsText": json.dumps({"records": records_payload, "pairs": pairs}, separators=(",", ":")),
                       "rawText": "".join(texts)})
    return inputs


def operation_probe(pc, cc, fn):
    import pandas as pd
    from gbdraw.web_support import orthogroup_metadata as metadata
    from gbdraw.analysis import collinearity_units as cu

    counters = Counter()
    profile = cProfile.Profile()
    aggregate_depth = 0
    with ExitStack() as stack:
        # Scope pandas observations to the actual aggregation call, excluding
        # fixture preparation, normalization and the semantic serializer.
        for name in ("itertuples", "groupby", "copy", "_constructor_from_mgr"):
            original = getattr(pd.DataFrame, name)
            def pandas_call(df, *args, _name=name, _original=original, **kwargs):
                observed = aggregate_depth > 0
                if observed:
                    counters["hsp.pandas." + _name + ".calls"] += 1
                result = _original(df, *args, **kwargs)
                if observed and _name == "itertuples":
                    def counted_rows():
                        for row in result:
                            counters["hsp.itertuples.rows"] += 1
                            yield row
                    return counted_rows()
                return result
            stack.enter_context(patch.object(pd.DataFrame, name, pandas_call))
        # Observe removed unit work only in this module; never instrument timing.
        unit_source = inspect.getsource(cu)
        member_sort_lines = {n.lineno for n in ast.walk(ast.parse(unit_source))
                             if isinstance(n, ast.Call) and isinstance(n.func, ast.Name)
                             and n.func.id == "sorted" and n.args
                             and isinstance(n.args[0], ast.Name) and n.args[0].id == "members"}
        for key in ("instances", "memberSortCalls", "memberSortReferences", "explicitSetCalls", "aliasVisits"):
            counters["units." + key] = 0
        original_unit = cu.CollinearityUnit
        def make_unit(*args, **kwargs):
            counters["units.instances"] += 1
            return original_unit(*args, **kwargs)
        def unit_sorted(values, *args, **kwargs):
            counters["units.sortCalls"] += 1
            if sys._getframe(1).f_lineno in member_sort_lines:
                counters["units.memberSortCalls"] += 1
                counters["units.memberSortReferences"] += len(values)
            return sorted(values, *args, **kwargs)
        def unit_set(*args):
            counters["units.explicitSetCalls"] += 1
            return set(*args)
        original_aliases = cu._unit_aliases
        def unit_aliases(**kwargs):
            result = original_aliases(**kwargs)
            counters["units.aliasVisits"] += len(result)
            return result
        stack.enter_context(patch.object(cu, "CollinearityUnit", make_unit))
        stack.enter_context(patch.object(cu, "sorted", unit_sorted, create=True))
        stack.enter_context(patch.object(cu, "set", unit_set, create=True))
        stack.enter_context(patch.object(cu, "_unit_aliases", unit_aliases))
        merge_depth = 0
        conflict_depth = 0
        bisect_depth = 0
        class CountedAnchors:
            def __init__(self, values):
                self.values = values
            def __len__(self):
                return len(self.values)
            def __getitem__(self, index):
                counters["merge.queryIndexReads" if bisect_depth else "merge.conflictCandidateVisits"] += 1
                return self.values[index]
            def __iter__(self):
                for anchor in self.values:
                    counters["merge.conflictCandidateVisits"] += 1
                    yield anchor
        original_merge = cc._merge_lossless_clusters
        def merge_probe(*args, **kwargs):
            nonlocal merge_depth
            merge_depth += 1
            counters["merge.calls"] += 1
            try:
                return original_merge(*args, **kwargs)
            finally:
                merge_depth -= 1
        stack.enter_context(patch.object(cc, "_merge_lossless_clusters", merge_probe))
        original_can_merge = cc._lossless_clusters_can_merge
        def can_merge_probe(left, right, **kwargs):
            counters["merge.boundaryTests"] += 1
            left_path = left.anchors if hasattr(left, "anchors") else left
            orientation = left.orientation if hasattr(left, "orientation") else kwargs["orientation"]
            if orientation == right.orientation and left_path and right.anchors:
                counters["merge.endpointReads"] += 2
            accepted = original_can_merge(left, right, **kwargs)
            if accepted:
                counters["merge.accepted"] += 1
                counters["merge.rightAnchorReferencesJoined"] += len(right.anchors)
            return accepted
        stack.enter_context(patch.object(cc, "_lossless_clusters_can_merge", can_merge_probe))
        original_conflicts = cc._lossless_conflicts_between_clusters
        def conflict_probe(left, right, anchors, **kwargs):
            nonlocal conflict_depth
            conflict_depth += 1
            counters["merge.conflictCalls"] += 1
            if (left.anchors if hasattr(left, "anchors") else left) and (right.anchors if hasattr(right, "anchors") else right):
                counters["merge.endpointReads"] += 2
            try:
                return original_conflicts(left, right, CountedAnchors(anchors), **kwargs)
            finally:
                conflict_depth -= 1
        stack.enter_context(patch.object(cc, "_lossless_conflicts_between_clusters", conflict_probe))
        for name in ("bisect_left", "bisect_right"):
            if not hasattr(cc, name):
                continue
            original = getattr(cc, name)
            def bisect_probe(*args, _original=original, **kwargs):
                nonlocal bisect_depth
                bisect_depth += 1
                try:
                    return _original(*args, **kwargs)
                finally:
                    bisect_depth -= 1
            stack.enter_context(patch.object(cc, name, bisect_probe))
        for name in ("_path_sorted_anchors", "_lossless_block_from_anchors"):
            original = getattr(cc, name)
            def merge_work(*args, _name=name, _original=original, **kwargs):
                if merge_depth:
                    kind = "pathSort" if _name == "_path_sorted_anchors" else "materialization"
                    values = args[0] if args else kwargs["anchors"]
                    counters["merge." + kind + "Calls"] += 1
                    counters["merge." + kind + "AnchorReferences"] += len(values)
                    if conflict_depth and kind == "pathSort":
                        counters["merge.conflictEndpointSortCalls"] += 1
                return _original(*args, **kwargs)
            stack.enter_context(patch.object(cc, name, merge_work))
        if hasattr(metadata, "_index_rbh_groups"):
            original_rbh_index = metadata._index_rbh_groups
            def rbh_index(groups):
                counters["metadata.rbhIndexBuilds"] += 1
                counters["metadata.rbhIndexInputGroups"] += len(groups)
                counters["metadata.rbhIndexInputMembers"] += sum(map(len, groups.values()))
                started = time.perf_counter_ns()
                ids, index = original_rbh_index(groups)
                counters["metadata.rbhIndexBuildNanoseconds"] += time.perf_counter_ns() - started
                counters["metadata.rbhIndexReferences"] += sum(map(len, index.values()))
                counters["metadata.rbhIndexOwnedContainerBytes"] = max(
                    counters["metadata.rbhIndexOwnedContainerBytes"],
                    sys.getsizeof(ids) + sys.getsizeof(index) + sum(map(sys.getsizeof, index.values())))
                return ids, index
            stack.enter_context(patch.object(metadata, "_index_rbh_groups", rbh_index))
        original_text = metadata._text
        def metadata_text(value):
            counters["metadata.textNormalizations"] += 1
            return original_text(value)
        stack.enter_context(patch.object(metadata, "_text", metadata_text))
        # The old metadata owner materializes both edge tuples for every lookup.
        def edge_amount(args, kwargs):
            if hasattr(pc, "_index_orthogroup_edges"):
                return 0
            groups, gid = args[:2]
            return (len(groups.ortholog_edges_by_orthogroup_id.get(gid, ())) +
                    len(groups.related_edges_by_orthogroup_id.get(gid, ()))) if groups and gid else 0
        original_anchor_metadata = cc._orthogroup_edge_metadata_for_anchor
        def anchor_metadata(*args, **kwargs):
            counters["metadata.anchorLookups"] += 1
            if not hasattr(pc, "_index_orthogroup_edges"):
                counters["metadata.anchorEdgeReferencesMaterialized"] += edge_amount((args[3], args[2]), {})
            return original_anchor_metadata(*args, **kwargs)
        stack.enter_context(patch.object(cc, "_orthogroup_edge_metadata_for_anchor", anchor_metadata))
        original_best = pc._best_evidence_between_protein_and_members
        def best_visits(*args, **kwargs):
            legacy = "same_record" in kwargs
            counters["support.memberVisits" if legacy else "support.evidenceVisits"] += len(args[1])
            counters["support.reductionCalls"] += 1
            return original_best(*args, **kwargs)
        stack.enter_context(patch.object(pc, "_best_evidence_between_protein_and_members", best_visits))
        for name in ("_index_core_support_evidence", "_index_orthogroup_edges", "_orthogroup_member_counts"):
            if not hasattr(pc, name):
                continue
            original = getattr(pc, name)
            def index_build(*args, _name=name, _original=original, **kwargs):
                counters[_name + ".calls"] += 1
                if _name == "_index_orthogroup_edges":
                    groups, gid = args[:2]
                    if args[3][1] == 0:
                        counters[_name + ".inputItems"] += len(groups.ortholog_edges_by_orthogroup_id.get(gid, ())) + len(groups.related_edges_by_orthogroup_id.get(gid, ()))
                else:
                    counters[_name + ".inputItems"] += len(args[0])
                started = time.perf_counter_ns()
                index = _original(*args, **kwargs)
                counters[_name + ".buildNanoseconds"] += time.perf_counter_ns() - started
                if _name == "_index_orthogroup_edges":
                    counters[_name + ".edgeVisits"] += index[1] - args[3][1]
                    return index
                if _name == "_orthogroup_member_counts":
                    counters[_name + ".retainedEntries"] += len(index)
                    counters[_name + ".maxOwnedContainerBytes"] = max(counters[_name + ".maxOwnedContainerBytes"], sys.getsizeof(index))
                if _name == "_index_core_support_evidence":
                    buckets = [values for groups in index.values() for values in groups.values()]
                    counters[_name + ".retainedEvidenceReferences"] += sum(map(len, buckets))
                    # Owned containers only: existing row/string objects excluded.
                    owned = sys.getsizeof(index) + sum(sys.getsizeof(groups) for groups in index.values())
                    owned += sum(sys.getsizeof(values) + sum(sys.getsizeof(item) for item in values) for values in buckets)
                    counters[_name + ".maxOwnedContainerBytes"] = max(counters[_name + ".maxOwnedContainerBytes"], owned)
                return index
            stack.enter_context(patch.object(pc, name, index_build))
            # Collinearity imports shared metadata functions directly.
            if hasattr(cc, name):
                stack.enter_context(patch.object(cc, name, index_build))
        for module, name, amount in (
            (pc, "validate_protein_identity_manifest", None),
            (pc, "_aggregate_hsps_by_protein_pair", lambda a, k: len(a[0])),
            (pc, "parse_losatp_outfmt6", None),
            (pc, "_raw_hsp_representative_rank", None),
            (pc, "_build_core_support_candidate", None),
            (pc, "_edge_metadata_for_protein_pair", edge_amount),
            (pc, "_build_ortholog_path_indexes" if hasattr(pc, "_build_ortholog_path_indexes") else "_build_ortholog_paths", None),
            (cc, "_lossless_conflicts_between_clusters", lambda a, k: len(a[2])),
        ):
            original = getattr(module, name)
            def counted(*args, _name=name, _original=original, _amount=amount, **kwargs):
                nonlocal aggregate_depth
                counters[_name + ".calls"] += 1
                if _amount:
                    counters[_name + ".inputItems"] += _amount(args, kwargs)
                aggregation = _name == "_aggregate_hsps_by_protein_pair"
                aggregate_depth += int(aggregation)
                try:
                    indexes = args[4] if _name == "_edge_metadata_for_protein_pair" and len(args) > 4 else None
                    state = indexes.get(args[1]) if indexes is not None else None
                    previous = len(state[0]) if state else 0
                    result = _original(*args, **kwargs)
                    state = indexes.get(args[1]) if indexes is not None else None
                    if state:
                        indexed = state[0]
                        counters["_index_orthogroup_edges.retainedEntries"] += len(indexed) - previous
                        owned = sys.getsizeof(state) + sys.getsizeof(indexed) + sum(map(sys.getsizeof, indexed))
                        counters["_index_orthogroup_edges.maxOwnedContainerBytes"] = max(counters["_index_orthogroup_edges.maxOwnedContainerBytes"], owned)
                    return result
                finally:
                    aggregate_depth -= int(aggregation)
            stack.enter_context(patch.object(module, name, counted))
            if module is pc and getattr(cc, name, None) is original:
                stack.enter_context(patch.object(cc, name, counted))
        if hasattr(pc, "OrthologPathCollection"):
            for name in ("iter_paths", "_path"):
                original = getattr(pc.OrthologPathCollection, name)
                counters["paths." + name + ".calls"] = 0
                def paths_call(*args, _name=name, _original=original, **kwargs):
                    counters["paths." + _name + ".calls"] += 1
                    return _original(*args, **kwargs)
                stack.enter_context(patch.object(pc.OrthologPathCollection, name, paths_call))
            original_adapter = pc.materialize_ortholog_paths
            counters["paths.materialize.calls"] = 0
            def adapter(*args, **kwargs):
                counters["paths.materialize.calls"] += 1
                return original_adapter(*args, **kwargs)
            stack.enter_context(patch.object(pc, "materialize_ortholog_paths", adapter))
        result = profile.runcall(fn)
    stats = pstats.Stats(profile)
    # Profile details are diagnostic; none of these instrumented times enters the timing gate.
    calls = []
    for (filename, line, name), (primitive, total, own, cumulative, _) in stats.stats.items():
        if "/gbdraw/" in filename or name in ("itertuples", "groupby"):
            calls.append({"file": filename.split("/gbdraw/", 1)[-1], "line": line, "name": name,
                          "calls": total, "primitiveCalls": primitive, "ownSeconds": own,
                          "cumulativeSeconds": cumulative})
    return result, dict(counters), sorted(calls, key=lambda x: -x["cumulativeSeconds"])[:35]


def summarize_result(value):
    groups = getattr(value, "orthogroups", None)
    summary = {}
    if isinstance(groups, dict):
        groups = value
    if groups is not None:
        summary = {"groups": len(groups.orthogroups),
                   "members": sum(len(x) for x in groups.orthogroups.values()),
                   "edges": sum(len(x) for x in groups.ortholog_edges_by_orthogroup_id.values()),
                   "paths": (sum(x.count for x in groups.path_indexes_by_orthogroup_id.values())
                             if hasattr(groups, "path_indexes_by_orthogroup_id") else
                             sum(len(x) for x in groups.ortholog_paths_by_orthogroup_id.values()))}
        if hasattr(groups, "path_indexes_by_orthogroup_id"):
            summary.update(nodes=sum(len(x.nodes) for x in groups.path_indexes_by_orthogroup_id.values()),
                           transitions=sum(len(x.transitions) for x in groups.path_indexes_by_orthogroup_id.values()))
    if hasattr(value, "blocks"):
        summary.update({"blocks": len(value.blocks), "anchors": sum(len(b.anchors) for b in value.blocks)})
    if isinstance(value, dict) and ("warmHits" in value or "observedException" in value):
        summary.update(value)
    return summary


def measure_stage(pc, cc, fn, args, artifact):
    samples = []
    hashes = []
    counts = None
    profile = None
    warmups = args.warmups if args.measure == "timing" else 0
    n = args.samples if args.measure == "timing" else 1
    for _ in range(warmups):
        fn()
    for _ in range(n):
        gc.collect()
        if args.measure == "memory":
            tracemalloc.start()
            try:
                result = fn()
                retained, peak = tracemalloc.get_traced_memory()
            finally:
                tracemalloc.stop()
            samples.append(peak)
        elif args.measure == "probe":
            result, counts, profile = operation_probe(pc, cc, fn)
        else:
            start = time.perf_counter_ns()
            result = fn()
            samples.append((time.perf_counter_ns()-start)/1e6)
        diagnostics = result.get("diagnostics") if isinstance(result, dict) else None
        data = json_bytes(canonical(result["semantic"] if diagnostics is not None else result))
        hashes.append(digest(data))
        summary = diagnostics if diagnostics is not None else summarize_result(result)
        output_bytes = (len(result) if isinstance(result, bytes) else
                        len(result.encode()) if isinstance(result, str) else None)
        del result
    if len(set(hashes)) != 1:
        raise ValueError("semantic results changed between identical samples")
    if artifact:
        artifact.parent.mkdir(parents=True, exist_ok=True)
        artifact.write_bytes(gzip.compress(data, mtime=0))
    report = {"semanticSha256": hashes[0], "semanticBytes": len(data), "summary": summary,
              "samples": samples, "unit": {"timing": "ms", "memory": "tracemalloc bytes", "probe": "counts/profile"}[args.measure],
              "sampleSemanticSha256": hashes}
    if output_bytes is not None:
        report["outputBytes"] = output_bytes
    if samples:
        median = statistics.median(samples)
        mad = statistics.median(abs(x-median) for x in samples)
        report.update({"median": median, "mad": mad, "noisePct": 100 * mad/median if median else 0})
    if args.measure == "memory":
        report["retainedBytes"] = retained
    if counts is not None:
        report.update({"operations": counts, "profile": profile})
    return report


def compare_reports(base, current, *, semantic_only=False):
    for report in (base, current):
        if report.get("benchmark") != NAME or report.get("schema") != 1 or not report.get("cases"):
            raise ValueError("invalid benchmark report")
        if report.get("exitCode", 0) != 0:
            raise ValueError("cannot compare a failed benchmark run")
    for key in ("benchmark", "schema", "settings", "dependencies"):
        if base[key] != current[key]:
            raise ValueError(f"incomparable reports: {key}")
    if not semantic_only and base["measurement"] != current["measurement"]:
        raise ValueError("incomparable reports: measurement")
    if set(base["cases"]) != set(current["cases"]):
        raise ValueError("incomparable case inventory")
    policy = base["settings"]["policy"]
    rows = {}
    for name, bcase in base["cases"].items():
        ccase = current["cases"][name]
        if bcase["fixture"] != ccase["fixture"] or bcase["stages"].keys() != ccase["stages"].keys():
            raise ValueError(f"incomparable inputs/settings/stages: {name}")
        for stage, b in bcase["stages"].items():
            c = ccase["stages"][stage]
            equal = b["semanticSha256"] == c["semanticSha256"]
            decision = "pass" if equal else "semantic_mismatch"
            row = {"semanticEqual": equal}
            if base["measurement"] == "timing" and not semantic_only:
                for measurement in (b, c):
                    if not measurement["samples"] or any(not math.isfinite(x) or x <= 0 for x in measurement["samples"]):
                        raise ValueError("invalid timing samples")
                delta = 100 * (c["median"] / b["median"] - 1)
                row["timeDeltaPct"] = delta
                if equal:
                    if min(len(b["samples"]), len(c["samples"])) < policy["samples"] or max(b["noisePct"], c["noisePct"]) > policy["maxNoisePct"]:
                        decision = "inconclusive"
                    elif delta > policy["regressionPct"]:
                        decision = "regression"
            row["decision"] = decision
            rows[f"{name}/{stage}"] = row
    failed = any(r["decision"] in ("semantic_mismatch", "regression") for r in rows.values())
    inconclusive = any(r["decision"] == "inconclusive" for r in rows.values())
    code = 1 if failed else 2 if inconclusive else 0
    return {"benchmark": NAME, "baselineHead": base["source"]["head"],
            "currentHead": current["source"]["head"], "semanticOnly": semantic_only,
            "comparisons": rows, "exitCode": code}, code


def read_report(path):
    data = path.read_bytes()
    return json.loads(gzip.decompress(data) if path.suffix == ".gz" else data)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run")
    run.add_argument("--source-root", type=Path, required=True)
    run.add_argument("--cases", nargs="+", choices=CASES, default=list(CASES))
    run.add_argument("--measure", choices=("timing", "probe", "memory"), default="timing")
    run.add_argument("--seed", type=int, default=SEED)
    run.add_argument("--warmups", type=int, default=POLICY["warmups"])
    run.add_argument("--samples", type=int, default=POLICY["samples"])
    run.add_argument("--path-representation", choices=("graph", "exhaustive"), default="graph",
                     help="Same-source explicit-output control, restricted to R<=16 path cases")
    run.add_argument("--unit-stages", action="store_true", help="Add prepared extraction -> full unit index for Collinear cases")
    run.add_argument("--merge-stages", action="store_true", help="Add isolated merge stages with captured production inputs")
    run.add_argument("--stages", nargs="+", help="Measure only named stages; preparation remains outside stage timing")
    run.add_argument("--artifacts", type=Path)
    run.add_argument("--output", type=Path, required=True)
    comp = sub.add_parser("compare")
    comp.add_argument("--baseline", type=Path, required=True)
    comp.add_argument("--current", type=Path, required=True)
    comp.add_argument("--output", type=Path, required=True)
    comp.add_argument("--semantics-only", action="store_true", help="compare scientific output across timing/profile/memory runs")
    browser = sub.add_parser("browser", help="measure the production diagram Worker helper with Python Playwright")
    browser.add_argument("--source-root", type=Path, required=True)
    browser.add_argument("--samples", type=int, default=POLICY["samples"])
    browser.add_argument("--output", type=Path, required=True)
    paths_browser = sub.add_parser("path-browser", help="S06 production conversion Worker timing/memory")
    paths_browser.add_argument("--source-root", type=Path, required=True)
    paths_browser.add_argument("--samples", type=int, default=POLICY["samples"])
    paths_browser.add_argument("--measure", choices=("timing", "memory"), default="timing")
    paths_browser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.command == "compare":
        report, code = compare_reports(read_report(args.baseline), read_report(args.current), semantic_only=args.semantics_only)
    elif args.command == "path-browser":
        from protein_comparison_browser import run_path_browser
        root = args.source_root.resolve()
        pc, cc = load_source(root)
        report = {"benchmark": NAME, "schema": 1, "command": sys.argv,
                  "source": source_info(root), "measurement": args.measure,
                  "result": run_path_browser(root, path_browser_inputs(root, pc, cc), args.samples, args.measure)}
        code = 0
    elif args.command == "browser":
        from protein_comparison_browser import run_browser
        if args.samples < 1:
            parser.error("samples must be positive")
        root = args.source_root.resolve()
        pc, _ = load_source(root)
        fixture, _ = manifest_case(root, pc)
        payload = fixture.pop("browserPayload")
        expected = fixture.pop("expectedKeys")
        report = {"benchmark": NAME, "schema": 1, "command": sys.argv,
                  "source": source_info(root), "fixture": fixture, "measurement": "browser",
                  "result": run_browser(root, payload, expected, args.samples), "exitCode": 0}
        code = 0
    else:
        if args.samples < 1 or args.warmups < 0:
            parser.error("samples must be positive and warmups nonnegative")
        if args.path_representation == "exhaustive" and any(name not in {"path-8", "path-12", "path-16"} for name in args.cases):
            parser.error("Explicit-output controls are restricted to R=8/12/16; never expand R>=24")
        root = args.source_root.resolve()
        pc, cc = load_source(root)
        report = {"benchmark": NAME, "schema": 1, "command": sys.argv,
                  "source": source_info(root), "measurement": args.measure,
                  "settings": {"policy": POLICY, "seed": args.seed, "warmups": args.warmups, "samples": args.samples},
                  "dependencies": {"python": platform.python_version(), **{p: importlib.metadata.version(p) for p in ("pandas", "numpy", "biopython", "svgwrite")}},
                  "environment": {"platform": platform.platform(), "pythonExecutable": sys.executable,
                                  "cpu": platform.processor(), "rssScope": "whole process including fixture preparation; KiB on Linux"},
                  "pathTheoryOnly": [{"records": r, "pathsDecimal": str(2**(r-2)), "executed": False} for r in (24, 32, 56)],
                  "cases": {}}
        for name in args.cases:
            print(f"{args.measure}: {name}", file=sys.stderr, flush=True)
            fixture, operations = build_case(root, pc, cc, name, args.seed, args.path_representation, args.unit_stages)
            if args.merge_stages:
                operations.update(isolated_merge_operations(cc, operations))
            if args.path_representation != "graph":
                fixture["pathRepresentation"] = args.path_representation
            if args.stages:
                operations = {stage: fn for stage, fn in operations.items() if stage in args.stages}
                if not operations:
                    parser.error(f"No selected stages in {name}")
            # A large transport fixture is written once, not embedded in every report.
            browser_payload = fixture.pop("browserPayload", None)
            if browser_payload and args.artifacts:
                args.artifacts.mkdir(parents=True, exist_ok=True)
                (args.artifacts / "manifest-browser.json").write_bytes(json_bytes(browser_payload))
            fixture.pop("expectedKeys", None)
            result = {"fixture": fixture, "stages": {}}
            for stage, fn in operations.items():
                print(f"  {name}/{stage}", file=sys.stderr, flush=True)
                artifact = args.artifacts / f"{name}.{stage}.json.gz" if args.artifacts else None
                result["stages"][stage] = measure_stage(pc, cc, fn, args, artifact)
            if name.startswith("path-"):
                summary = result["stages"]["selector"]["summary"]
                if summary["paths"] != fixture["expectedPaths"] or summary["edges"] != fixture["expectedEdges"]:
                    raise ValueError("complete ordered record DAG path formula mismatch")
            report["cases"][name] = result
            guard_imports(root)
        report["processMaxRssKiB"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        report["exitCode"] = code = 0
    args.output.parent.mkdir(parents=True, exist_ok=True)
    data = (json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False) + "\n").encode()
    args.output.write_bytes(gzip.compress(data, mtime=0) if args.output.suffix == ".gz" else data)
    return code


if __name__ == "__main__":
    raise SystemExit(main())
