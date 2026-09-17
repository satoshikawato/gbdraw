"""Lossless ordered ortholog paths, without expansion on the normal path.

The DAG denotes protein-deduplicated source-to-sink paths. Parallel evidence
chooses the smallest edge ID; all evidence remains in the enclosing result.
Explicit collections preserve supplied legacy corpora, including non-DAG paths.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass, field
import re
from typing import Iterator, Literal, Mapping, Sequence

from gbdraw.exceptions import ValidationError

PATH_ORDERING = "protein-key-edge-tuple-v1"
ProteinOrderKey = tuple[int, int, int, str]


@dataclass(frozen=True)
class OrthologPath:
    """One explicitly retrieved path, with the original tuple contract."""

    orthogroup_id: str
    path_id: str
    protein_ids: tuple[str, ...]
    edge_ids: tuple[str, ...]
    shared_protein_ids: tuple[str, ...] = ()


def ortholog_edge_id(edge) -> str:
    return (
        f"{edge.orthogroup_id}:"
        f"{edge.query_record_index}:{edge.query_protein_id}->"
        f"{edge.subject_record_index}:{edge.subject_protein_id}:"
        f"{edge.edge_kind}"
    )


@dataclass(frozen=True)
class OrthologPathCollection:
    """An immutable DAG or an exact explicit corpus; never implicitly iterable.

    ``count`` is an arbitrary-precision integer (there is deliberately no len).
    ``path_at`` uses one-based ranks. ``iter_paths`` is the explicit exhaustive
    boundary and costs at least the size of its output. Construction sorts nodes
    and successors, then performs DAG dynamic programming with big integers.
    """

    orthogroup_id: str
    kind: Literal["dag", "explicit"]
    nodes: tuple[tuple[str, ProteinOrderKey], ...] = ()
    transitions: tuple[tuple[str, str, str], ...] = ()
    paths: tuple[OrthologPath, ...] = ()
    _keys: dict = field(init=False, repr=False, compare=False)
    _children: dict = field(init=False, repr=False, compare=False)
    _chosen: dict = field(init=False, repr=False, compare=False)
    _suffix: dict = field(init=False, repr=False, compare=False)
    _offset: dict = field(init=False, repr=False, compare=False)
    _starts: dict = field(init=False, repr=False, compare=False)
    _containing: dict = field(init=False, repr=False, compare=False)
    _first: dict = field(init=False, repr=False, compare=False)
    _count: int = field(init=False, repr=False, compare=False)
    _tied_keys: bool = field(init=False, repr=False, compare=False)

    def __post_init__(self):
        if not isinstance(self.orthogroup_id, str):
            raise ValidationError("Path collection requires a string orthogroup ID")
        if self.kind == "explicit":
            if self.nodes or self.transitions or type(self.paths) is not tuple:
                raise ValidationError("Explicit paths must be a tuple without a DAG")
            if not all(isinstance(p, OrthologPath) for p in self.paths):
                raise ValidationError("Explicit corpus requires OrthologPath values")
            # Keep arbitrary IDs, repeated paths, cycles and supplied shared data.
            containing = {}
            for path in self.paths:
                for pid in set(path.protein_ids):
                    containing[pid] = containing.get(pid, 0) + 1
            state = ({}, {}, {}, {}, {}, {}, containing, {}, len(self.paths))
        elif self.kind == "dag":
            if self.paths:
                raise ValidationError("DAG collections cannot contain explicit paths")
            state = self._build_index()
        else:
            raise ValidationError("Unknown ortholog path representation")
        for name, value in zip(
            ("_keys", "_children", "_chosen", "_suffix", "_offset", "_starts",
             "_containing", "_first", "_count"), state
        ):
            object.__setattr__(self, name, value)
        object.__setattr__(self, "_tied_keys", len(set(self._keys.values())) != len(self._keys))
        if self._tied_keys:
            self._index_tied_first_paths()

    def _build_index(self):
        keys = {}
        for pid, key in self.nodes:
            if (type(pid) is not str or pid in keys or len(key) != 4
                    or any(type(v) is not int for v in key[:3]) or type(key[3]) is not str):
                raise ValidationError("DAG requires unique node IDs and protein order keys")
            keys[pid] = key
        if tuple(keys) != tuple(sorted(keys, key=lambda pid: (keys[pid], pid))):
            raise ValidationError("DAG nodes are not in canonical order")
        children = {pid: [] for pid in keys}
        indegree = dict.fromkeys(keys, 0)
        chosen = {}
        for u, v, eid in self.transitions:
            if (u not in keys or v not in keys or (u, v) in chosen
                    or type(eid) is not str or not eid):
                raise ValidationError("Invalid or duplicate DAG transition")
            chosen[u, v] = eid
            children[u].append(v)
            indegree[v] += 1
        if len(set(chosen.values())) != len(chosen):
            raise ValidationError("DAG transition edge IDs must be unique")
        if any(not children[u] and not indegree[u] for u in keys):
            raise ValidationError("DAG paths cannot contain isolated nodes")
        if self.transitions != tuple(sorted(self.transitions, key=lambda e: (keys[e[0]], keys[e[1]], e[2]))):
            raise ValidationError("DAG transitions are not in canonical order")
        sources = [u for u in keys if not indegree[u]]
        queue = deque(sources)
        order = []
        while queue:
            u = queue.popleft()
            order.append(u)
            for v in children[u]:
                indegree[v] -= 1
                if not indegree[v]:
                    queue.append(v)
        if len(order) != len(keys):
            raise ValidationError("Ortholog path DAG contains a cycle")
        suffix, offsets = {}, {}
        for u in reversed(order):
            total = 0
            for v in children[u]:
                offsets[u, v] = total
                total += suffix[v]
            suffix[u] = total if children[u] else 1
        starts, prefix, earliest = {}, dict.fromkeys(keys, 0), {}
        count = 0
        for u in sources:
            starts[u] = earliest[u] = count
            prefix[u] = 1
            count += suffix[u]
        for u in order:
            for v in children[u]:
                prefix[v] += prefix[u]
                rank = earliest[u] + offsets[u, v]
                earliest[v] = min(earliest.get(v, rank), rank)
        containing = {u: prefix[u] * suffix[u] for u in keys}
        first = {eid: earliest[u] + offsets[u, v] + 1 for (u, v), eid in chosen.items()}
        return keys, children, chosen, suffix, offsets, starts, containing, first, count

    @classmethod
    def from_edges(cls, group_id: str, edges: Sequence, protein_map: Mapping):
        chosen = {}
        for edge in edges:
            u, v = edge.query_protein_id, edge.subject_protein_id
            if edge.edge_kind not in {"rbh", "coortholog"} or u not in protein_map or v not in protein_map:
                continue
            eid = ortholog_edge_id(edge)
            chosen[u, v] = min(chosen.get((u, v), eid), eid)
        keys = {}
        for pid in {pid for pair in chosen for pid in pair}:
            protein = protein_map[pid]
            keys[pid] = (int(protein.record_index), int(protein.start), int(protein.end), str(protein.protein_id))
        return cls(
            group_id, "dag",
            tuple(sorted(keys.items(), key=lambda item: (item[1], item[0]))),
            tuple(sorted(((u, v, eid) for (u, v), eid in chosen.items()),
                         key=lambda e: (keys[e[0]], keys[e[1]], e[2]))),
        )

    @property
    def count(self) -> int:
        return self._count

    def containing_count(self, protein_id: str) -> int:
        return self._containing[protein_id]

    def first_path_id(self, edge) -> str | None:
        rank = self._first.get(ortholog_edge_id(edge))
        return f"{self.orthogroup_id}.path_{rank}" if rank is not None else edge.path_id

    def validate_edges(self, edges: Sequence) -> None:
        """Check persisted graph selection and first-path IDs against evidence."""
        if self.kind == "explicit":
            return
        chosen = {}
        for edge in edges:
            u, v = edge.query_protein_id, edge.subject_protein_id
            if edge.edge_kind in {"rbh", "coortholog"}:
                if u not in self._keys or v not in self._keys:
                    raise ValidationError("DAG omits an evidence endpoint")
                eid = ortholog_edge_id(edge)
                chosen[u, v] = min(chosen.get((u, v), eid), eid)
            if edge.path_id != self.first_path_id(edge):
                raise ValidationError("DAG evidence has an inconsistent first-path ID")
        if chosen != self._chosen:
            raise ValidationError("DAG transitions do not match selected evidence")

    def _path(self, rank: int, proteins: Sequence[str], edge_ids: Sequence[str]) -> OrthologPath:
        return OrthologPath(
            self.orthogroup_id, f"{self.orthogroup_id}.path_{rank}",
            tuple(proteins), tuple(edge_ids),
            tuple(sorted((v for v in proteins if self._containing[v] > 1), key=self._keys.__getitem__)),
        )

    def _label_frontier(self, frontier):
        """Weighted next nodes; merging prefixes must retain their multiplicity."""
        following = {}
        for u, weight in frontier.items():
            for v in self._children[u]:
                following[v] = following.get(v, 0) + weight
        return following

    def _matching_suffixes(self, labels):
        matches = [{} for _ in labels]
        for position in range(len(labels) - 1, -1, -1):
            for u, key in self._keys.items():
                if key == labels[position]:
                    matches[position][u] = (
                        int(not self._children[u]) if position == len(labels) - 1 else
                        sum(matches[position + 1].get(v, 0) for v in self._children[u])
                    )
        return matches

    def _tied_path(self, rank):
        # Legacy sorting compares the entire protein-key sequence before any
        # edge ID. Caller-built maps can have equal keys, unlike extracted CDSs.
        remaining, labels = rank - 1, []
        frontier = dict.fromkeys(self._starts, 1)
        while frontier:
            for key in sorted({self._keys[u] for u in frontier}):
                selected = {u: w for u, w in frontier.items() if self._keys[u] == key}
                count = sum(w * self._suffix[u] for u, w in selected.items())
                if remaining < count:
                    break
                remaining -= count
            labels.append(key)
            terminal = sum(w for u, w in selected.items() if not self._children[u])
            if remaining < terminal:
                break
            remaining -= terminal
            frontier = self._label_frontier(selected)
        matches = self._matching_suffixes(labels)
        sources = [u for u in self._starts if matches[0].get(u)]
        proteins, edges = [], []
        for position in range(1, len(labels)):
            choices = sorted((self._chosen[u, v], u, v) for u in sources
                             for v in self._children[u] if matches[position].get(v))
            for eid, u, v in choices:
                if remaining < matches[position][v]:
                    break
                remaining -= matches[position][v]
            if not proteins:
                proteins.append(u)
            proteins.append(v)
            edges.append(eid)
            sources = [v]
        return self._path(rank, proteins, edges)

    def _tied_rank(self, proteins):
        labels = [self._keys[u] for u in proteins]
        rank, frontier = 1, dict.fromkeys(self._starts, 1)
        for position, key in enumerate(labels):
            rank += sum(w * self._suffix[u] for u, w in frontier.items() if self._keys[u] < key)
            frontier = {u: w for u, w in frontier.items() if self._keys[u] == key}
            if position < len(labels) - 1:
                rank += sum(w for u, w in frontier.items() if not self._children[u])
                frontier = self._label_frontier(frontier)
        matches = self._matching_suffixes(labels)
        sources = [u for u in self._starts if matches[0].get(u)]
        for position, (u, v) in enumerate(zip(proteins, proteins[1:]), 1):
            target = self._chosen[u, v]
            rank += sum(matches[position].get(y, 0) for x in sources for y in self._children[x]
                        if self._chosen[x, y] < target)
            sources = [v]
        return rank

    def _index_tied_first_paths(self):
        # _suffix insertion order is reverse topological. Find the least full
        # (key sequence, edge sequence) through each edge, then rank just that
        # path. This remains polynomial and never expands the path corpus.
        suffix = {}
        for u in self._suffix:
            suffix[u] = min(
                (((self._keys[u],) + suffix[v][0], (self._chosen[u, v],) + suffix[v][1],
                  (u,) + suffix[v][2]) for v in self._children[u]),
                default=((self._keys[u],), (), (u,)),
            )
        for (u, v), eid in self._chosen.items():
            best = {u: ((self._keys[u],) + suffix[v][0], (eid,) + suffix[v][1], (u,) + suffix[v][2])}
            for node in self._suffix:
                if node != u:
                    choices = [((self._keys[node],) + best[child][0],
                                (self._chosen[node, child],) + best[child][1], (node,) + best[child][2])
                               for child in self._children[node] if child in best]
                    if choices:
                        best[node] = min(choices)
            self._first[eid] = self._tied_rank(min(best[s] for s in self._starts if s in best)[2])

    def path_at(self, rank: int) -> OrthologPath:
        if type(rank) is not int:
            raise TypeError("Path rank must be an int, excluding bool")
        if not 1 <= rank <= self.count:
            raise IndexError("Path rank is outside the collection")
        if self.kind == "explicit":
            return self.paths[rank - 1]
        if self._tied_keys:
            return self._tied_path(rank)
        offset, proteins, edges = rank - 1, [], []
        candidates = self._starts
        while candidates:
            for node in candidates:
                if offset < self._suffix[node]:
                    if proteins:
                        edges.append(self._chosen[proteins[-1], node])
                    proteins.append(node)
                    candidates = self._children[node]
                    break
                offset -= self._suffix[node]
        return self._path(rank, proteins, edges)

    def path_by_id(self, path_id: str) -> OrthologPath:
        if type(path_id) is not str:
            raise TypeError("Path ID must be a string")
        if self.kind == "explicit":
            matches = [p for p in self.paths if p.path_id == path_id]
            if not matches:
                raise KeyError(path_id)
            if len(matches) != 1:
                raise ValidationError("Ambiguous explicit path ID")
            return matches[0]
        match = re.fullmatch(re.escape(self.orthogroup_id) + r"\.path_([1-9][0-9]*)", path_id)
        if match is None:
            raise ValidationError("Malformed or wrong-group path ID")
        return self.path_at(int(match[1]))

    def rank_of(self, protein_ids: tuple[str, ...]) -> int:
        if type(protein_ids) is not tuple or not all(type(v) is str for v in protein_ids):
            raise TypeError("protein_ids must be a tuple of strings")
        if self.kind == "explicit":
            matches = [rank for rank, p in enumerate(self.paths, 1) if p.protein_ids == protein_ids]
            if len(matches) != 1:
                raise ValidationError("Missing or ambiguous explicit path sequence")
            return matches[0]
        if len(protein_ids) < 2 or protein_ids[0] not in self._starts:
            raise ValidationError("Not a complete source-to-sink path")
        pairs = list(zip(protein_ids, protein_ids[1:]))
        if any(pair not in self._chosen for pair in pairs) or self._children[protein_ids[-1]]:
            raise ValidationError("Not a complete source-to-sink path")
        if self._tied_keys:
            return self._tied_rank(protein_ids)
        return 1 + self._starts[protein_ids[0]] + sum(self._offset[pair] for pair in pairs)

    def iter_paths(self) -> Iterator[OrthologPath]:
        """Explicit expansion with a reusable DFS stack, in legacy order."""
        if self.kind == "explicit":
            yield from self.paths
            return
        if self._tied_keys:
            for rank in range(1, self.count + 1):
                yield self.path_at(rank)
            return
        rank = 0
        for source in self._starts:
            proteins, edges = [source], []
            stack = [iter(self._children[source])]
            while stack:
                node = next(stack[-1], None)
                if node is None:
                    stack.pop()
                    proteins.pop()
                    if edges:
                        edges.pop()
                    continue
                edges.append(self._chosen[proteins[-1], node])
                proteins.append(node)
                stack.append(iter(self._children[node]))
                if not self._children[node]:
                    rank += 1
                    yield self._path(rank, proteins, edges)

    def to_payload(self) -> dict:
        """Persist only the lossless representation; no DP state or expansion."""
        body = {"orthogroupId": self.orthogroup_id, "kind": self.kind, "count": str(self.count)}
        if self.kind == "explicit":
            body["paths"] = self.paths
        else:
            body.update(ordering=PATH_ORDERING,
                        nodes=[{"proteinId": pid, "orderKey": key} for pid, key in self.nodes],
                        transitions=[{"queryProteinId": u, "subjectProteinId": v, "edgeId": eid}
                                     for u, v, eid in self.transitions])
        return body

    @classmethod
    def from_payload(cls, body: dict):
        common = {"orthogroupId", "kind", "count"}
        if body.get("kind") == "explicit":
            if set(body) != common | {"paths"}:
                raise ValidationError("Invalid explicit path payload fields")
            result = cls(body["orthogroupId"], "explicit", paths=tuple(body["paths"]))
        elif body.get("kind") == "dag":
            if set(body) != common | {"ordering", "nodes", "transitions"} or body["ordering"] != PATH_ORDERING:
                raise ValidationError("Invalid DAG path payload fields or ordering")
            if (any(set(n) != {"proteinId", "orderKey"} for n in body["nodes"])
                    or any(set(e) != {"queryProteinId", "subjectProteinId", "edgeId"} for e in body["transitions"])):
                raise ValidationError("Invalid DAG node or transition fields")
            result = cls(body["orthogroupId"], "dag",
                         tuple((n["proteinId"], tuple(n["orderKey"])) for n in body["nodes"]),
                         tuple((e["queryProteinId"], e["subjectProteinId"], e["edgeId"]) for e in body["transitions"]))
        else:
            raise ValidationError("Unknown path payload kind")
        if type(body["count"]) is not str or body["count"] != str(result.count):
            raise ValidationError("Path count must be the exact canonical decimal string")
        return result
