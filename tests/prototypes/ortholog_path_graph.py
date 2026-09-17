"""S02 test-only DAG path index. Never imported by production.

This is a feasibility oracle, not an approved API, writer, or cycle fallback.
Ordered paths are the current protein-deduplicated source-to-sink language.
Only explicit path retrieval constructs OrthologPath objects.
"""

from __future__ import annotations

from collections import deque
from dataclasses import replace
import re

from gbdraw.analysis.protein_colinearity import (
    OrthologPath,
    _edge_id,
    _protein_sort_key,
)
from gbdraw.exceptions import ValidationError


class DagPathIndex:
    """Count, rank, unrank, shared membership, and first-edge rank on a DAG."""

    def __init__(self, group_id, edges, protein_map):
        self.group_id = group_id
        self.edges = tuple(edges)
        self.keys = {pid: _protein_sort_key(p) for pid, p in protein_map.items()}
        if any(pid != p.protein_id for pid, p in protein_map.items()):
            raise ValidationError("DAG index requires canonical protein-map keys")
        self.chosen = {}
        for edge in edges:
            u, v = edge.query_protein_id, edge.subject_protein_id
            if edge.edge_kind not in {"rbh", "coortholog"} or u not in self.keys or v not in self.keys:
                continue
            # Minimizing each parallel edge minimizes the entire edge tuple.
            pair = (u, v)
            eid = _edge_id(edge)
            self.chosen[pair] = min(self.chosen.get(pair, eid), eid)
        nodes = {pid for pair in self.chosen for pid in pair}
        self.children = {pid: [] for pid in nodes}
        parents = {pid: [] for pid in nodes}
        for u, v in self.chosen:
            self.children[u].append(v)
            parents[v].append(u)
        for children in self.children.values():
            children.sort(key=self.keys.__getitem__)
        self.starts = sorted((v for v in nodes if not parents[v]), key=self.keys.__getitem__)
        indegree = {v: len(parents[v]) for v in nodes}
        queue = deque(self.starts)
        order = []
        while queue:
            u = queue.popleft()
            order.append(u)
            for v in self.children[u]:
                indegree[v] -= 1
                if indegree[v] == 0:
                    queue.append(v)
        if len(order) != len(nodes):
            raise ValidationError("DAG path index cannot represent a cycle")
        self.suffix = {}
        self.offset = {}
        for u in reversed(order):
            total = 0
            for v in self.children[u]:
                self.offset[u, v] = total
                total += self.suffix[v]
            self.suffix[u] = total if self.children[u] else 1
        self.count = sum(self.suffix[u] for u in self.starts)
        self.prefix = dict.fromkeys(nodes, 0)
        self.first_prefix_rank = {}
        start_offset = 0
        self.start_offsets = {}
        for u in self.starts:
            self.prefix[u] = 1
            self.first_prefix_rank[u] = start_offset
            self.start_offsets[u] = start_offset
            start_offset += self.suffix[u]
        for u in order:
            for v in self.children[u]:
                self.prefix[v] += self.prefix[u]
                rank = self.first_prefix_rank[u] + self.offset[u, v]
                self.first_prefix_rank[v] = min(self.first_prefix_rank.get(v, rank), rank)
        self.containing_count = {v: self.prefix[v] * self.suffix[v] for v in nodes}
        self.first_edge_rank = {
            eid: self.first_prefix_rank[u] + self.offset[u, v] + 1
            for (u, v), eid in self.chosen.items()
        }
        self.materialized_paths = 0

    def path_at(self, rank: int) -> OrthologPath:
        """One-based rank, with exactly the legacy sorted ID and shared tuple."""
        if type(rank) is not int:
            raise TypeError("path rank must be an int, excluding bool")
        if not 1 <= rank <= self.count:
            raise IndexError("path rank is outside the collection")
        offset = rank - 1
        proteins = []
        candidates = self.starts
        while candidates:
            for node in candidates:
                if offset < self.suffix[node]:
                    proteins.append(node)
                    candidates = self.children[node]
                    break
                offset -= self.suffix[node]
            else:
                raise AssertionError("invalid suffix count")
        self.materialized_paths += 1
        return OrthologPath(
            orthogroup_id=self.group_id,
            path_id=f"{self.group_id}.path_{rank}",
            protein_ids=tuple(proteins),
            edge_ids=tuple(self.chosen[u, v] for u, v in zip(proteins, proteins[1:])),
            shared_protein_ids=tuple(sorted(
                (v for v in proteins if self.containing_count[v] > 1),
                key=self.keys.__getitem__,
            )),
        )

    def rank_of(self, protein_ids: tuple[str, ...]) -> int:
        if type(protein_ids) is not tuple or not all(type(v) is str for v in protein_ids):
            raise TypeError("protein_ids must be a tuple of strings")
        if len(protein_ids) < 2 or protein_ids[0] not in self.start_offsets:
            raise ValidationError("not a complete source-to-sink path")
        if any((u, v) not in self.chosen for u, v in zip(protein_ids, protein_ids[1:])):
            raise ValidationError("path contains an unknown transition")
        if self.children[protein_ids[-1]]:
            raise ValidationError("path does not end at a sink")
        return 1 + self.start_offsets[protein_ids[0]] + sum(
            self.offset[u, v] for u, v in zip(protein_ids, protein_ids[1:])
        )

    def path_by_id(self, path_id: str) -> OrthologPath:
        if type(path_id) is not str:
            raise TypeError("path ID must be a string")
        match = re.fullmatch(re.escape(self.group_id) + r"\.path_([1-9][0-9]*)", path_id)
        if match is None:
            raise ValidationError("malformed or wrong-group path ID")
        return self.path_at(int(match[1]))

    def iter_paths(self):
        for rank in range(1, self.count + 1):
            yield self.path_at(rank)

    def updated_edges(self):
        return tuple(
            replace(edge, path_id=f"{self.group_id}.path_{self.first_edge_rank[_edge_id(edge)]}")
            if _edge_id(edge) in self.first_edge_rank else edge
            for edge in self.edges
        )

    def summary(self):
        """Test-only JSON transport probe; not a candidate production writer."""
        return {
            "count": str(self.count),
            "containingCounts": {v: str(n) for v, n in sorted(self.containing_count.items())},
            "firstEdgeRanks": {eid: str(n) for eid, n in sorted(self.first_edge_rank.items())},
            "nodes": len(self.children),
            "transitions": len(self.chosen),
            "materializedPaths": self.materialized_paths,
        }
