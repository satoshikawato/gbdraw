#!/usr/bin/env python3
"""Reproduce protein-comparison baselines without changing production code.

Run each source tree in a fresh process. Timing (one warmup, seven samples),
profiling/counters, and tracemalloc are separate invocations. Inputs and complete
ordered stage results are hashed; --artifacts also saves the semantic oracles.
The benchmark records current behavior, including errors, not Product authority.
"""
from __future__ import annotations

import argparse
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
POLICY = {"warmups": 1, "samples": 7, "regressionPct": 10.0,
          "maxNoisePct": 5.0, "noise": "100 * MAD / median",
          "noisyDecision": "inconclusive; repeat both trees, never waive",
          "semanticGate": "exact ordered stage digest, including types and errors"}
GALLERIES = ("hepatoplasmataceae_collinear", "hepatoplasmataceae_orthogroup",
             "vibrio-harveyi-group-collinear")
CASES = ("hsp-edges", "hsp-1", "hsp-1000", "hsp-many", "dense-24", "support-edges",
         "sparse-200", "sparse-400", "sparse-800",
         "cache-49", "cache-64", "cache-81", "path-8", "path-12", "path-16",
         "merge-300", "merge-600", "merge-1200", "merge-edges", "manifest",
         "gallery-collinear", "gallery-orthogroup", "render-gallery")


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
        pm = {f"p{i}": protein(pc, f"p{i}", i) for i in range(size)}
        return pm, {(i, j): frame([hit(f"p{i}", f"p{j}")])
                    for i in range(size) for j in range(size)}
    if name.startswith("sparse-"):
        pm = {f"{p}{i}": protein(pc, f"{p}{i}", r, i)
              for i in range(size) for r, p in enumerate(("q", "s", "u"))}
        return pm, {(0, 1): frame([hit(f"q{i}", f"s{i}") for i in range(size)]),
                    (1, 0): frame([hit(f"s{i}", f"q{i}") for i in range(size)])}
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


def build_case(root, pc, cc, name, seed):
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
        gallery = GALLERIES[0 if name.endswith("collinear") else 1]
        _, extraction, records, raw, inventory = gallery_input(root, gallery, pc)
        settings = {"bitscore": 50, "evalue": 0.01, "identity": 0, "alignment_length": 0}
        def parse():
            return {k: pc.parse_losatp_outfmt6(v) for k, v in raw.items()}
        parsed = parse()
        def filter_hits():
            return {k: pc.filter_protein_hits_by_thresholds(v, **settings) for k, v in parsed.items()}
        tables = filter_hits()
        inventory.update({"benchmarkThresholds": settings, "memberMaxHits": 5,
                          "inference": True, "scope": "adjacent" if name.endswith("collinear") else "all",
                          "blockParameters": canonical(cc.LosslessCollinearityParameters()) if name.endswith("collinear") else None,
                          "unitMode": "auto", "edgeMode": "rbh", "maxRelatedEdges": 2,
                          "filteredRows": sum(len(x) for x in tables.values()),
                          "inputDataFrameDeepBytes": sum(int(x.memory_usage(deep=True).sum()) for x in tables.values())})
        def aggregate():
            return {k: pc._aggregate_hsps_by_protein_pair(v, extraction.protein_map) for k, v in tables.items()}
        if name.endswith("collinear"):
            def analyze():
                return cc.build_orthogroup_collinearity_blocks_from_hits(tables, extraction,
                            records=records, orthogroup_member_max_hits=5, infer_orthogroups=True, search_scope="adjacent")
        else:
            def analyze():
                return pc.select_rbh_orthogroup_edges_from_directional_hits(tables, extraction.protein_map,
                            record_count=len(records), orthogroup_member_max_hits=5)
        from gbdraw.session_request_codec import encode_canonical_typed_resource
        from gbdraw.web_support.orthogroup_metadata import serialize_orthogroups_payload
        result = analyze()
        # Stage inputs are prepared once outside measurement. The post_search stage
        # itself includes its own normalization, inference, projection and blocks.
        return inventory, {"parse": parse, "filter": filter_hits, "hsp_aggregate": aggregate,
            "post_search": analyze,
            "metadata": lambda: serialize_orthogroups_payload(result.orthogroups, records=records),
            "typed_serialization": lambda: encode_canonical_typed_resource("result", result if name.endswith("collinear") else result.orthogroups)}
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
            record_count=max(p.record_index for p in pm.values())+1, orthogroup_member_max_hits=None)
    return inventory, {"selector": select}


def operation_probe(pc, cc, fn):
    import pandas as pd

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
        for module, name, amount in (
            (pc, "validate_protein_identity_manifest", None),
            (pc, "_aggregate_hsps_by_protein_pair", lambda a, k: len(a[0])),
            (pc, "parse_losatp_outfmt6", None),
            (pc, "_raw_hsp_representative_rank", None),
            (pc, "_build_core_support_candidate", None),
            (pc, "_build_ortholog_paths", None),
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
                    return _original(*args, **kwargs)
                finally:
                    aggregate_depth -= int(aggregation)
            stack.enter_context(patch.object(module, name, counted))
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
    if groups is not None:
        summary = {"groups": len(groups.orthogroups),
                   "members": sum(len(x) for x in groups.orthogroups.values()),
                   "edges": sum(len(x) for x in groups.ortholog_edges_by_orthogroup_id.values()),
                   "paths": sum(len(x) for x in groups.ortholog_paths_by_orthogroup_id.values())}
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
                _, peak = tracemalloc.get_traced_memory()
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
        del result
    if len(set(hashes)) != 1:
        raise ValueError("semantic results changed between identical samples")
    if artifact:
        artifact.parent.mkdir(parents=True, exist_ok=True)
        artifact.write_bytes(gzip.compress(data, mtime=0))
    report = {"semanticSha256": hashes[0], "semanticBytes": len(data), "summary": summary,
              "samples": samples, "unit": {"timing": "ms", "memory": "tracemalloc bytes", "probe": "counts/profile"}[args.measure],
              "sampleSemanticSha256": hashes}
    if samples:
        median = statistics.median(samples)
        report.update({"median": median, "noisePct": 100 * statistics.median(abs(x-median) for x in samples)/median if median else 0})
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
                    if min(len(b["samples"]), len(c["samples"])) < POLICY["samples"] or max(b["noisePct"], c["noisePct"]) > POLICY["maxNoisePct"]:
                        decision = "inconclusive"
                    elif delta > POLICY["regressionPct"]:
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
    args = parser.parse_args(argv)
    if args.command == "compare":
        report, code = compare_reports(read_report(args.baseline), read_report(args.current), semantic_only=args.semantics_only)
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
            fixture, operations = build_case(root, pc, cc, name, args.seed)
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
