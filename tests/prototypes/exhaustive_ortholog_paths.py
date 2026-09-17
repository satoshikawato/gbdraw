"""Frozen S04 exhaustive oracle. Test-only; never use for R >= 24."""
from dataclasses import replace
from typing import Mapping, Sequence
from gbdraw.analysis.protein_colinearity import CdsProtein, OrthologEdge, OrthologPath


def _protein_sort_key(protein):
    return (int(protein.record_index), int(protein.start), int(protein.end), str(protein.protein_id))


def _edge_id(edge):
    return (
        f"{edge.orthogroup_id}:"
        f"{edge.query_record_index}:{edge.query_protein_id}->"
        f"{edge.subject_record_index}:{edge.subject_protein_id}:"
        f"{edge.edge_kind}"
    )


def _path_sort_key(
    protein_ids: Sequence[str],
    protein_map: Mapping[str, CdsProtein],
) -> tuple[tuple[int, int, int, str], ...]:
    return tuple(_protein_sort_key(protein_map[protein_id]) for protein_id in protein_ids)


def _build_ortholog_paths(
    edges_by_group: Mapping[str, Sequence[OrthologEdge]],
    protein_map: Mapping[str, CdsProtein],
) -> tuple[dict[str, tuple[OrthologEdge, ...]], dict[str, tuple[OrthologPath, ...]]]:
    updated_edges_by_group: dict[str, tuple[OrthologEdge, ...]] = {}
    paths_by_group: dict[str, tuple[OrthologPath, ...]] = {}
    for group_id, edges in edges_by_group.items():
        path_edges = [
            edge
            for edge in edges
            if edge.edge_kind in {"rbh", "coortholog"}
            and edge.query_protein_id in protein_map
            and edge.subject_protein_id in protein_map
        ]
        if not path_edges:
            updated_edges_by_group[group_id] = tuple(edges)
            paths_by_group[group_id] = ()
            continue
        outgoing: dict[str, list[OrthologEdge]] = {}
        incoming: dict[str, list[OrthologEdge]] = {}
        for edge in path_edges:
            outgoing.setdefault(edge.query_protein_id, []).append(edge)
            incoming.setdefault(edge.subject_protein_id, []).append(edge)
        for edge_list in outgoing.values():
            edge_list.sort(
                key=lambda edge: (
                    edge.subject_record_index,
                    _protein_sort_key(protein_map[edge.subject_protein_id]),
                    _edge_id(edge),
                )
            )
        nodes = set(outgoing).union(incoming)
        start_nodes = [
            node
            for node in nodes
            if node not in incoming
        ] or list(nodes)
        start_nodes.sort(key=lambda protein_id: _protein_sort_key(protein_map[protein_id]))

        raw_paths: list[tuple[tuple[str, ...], tuple[str, ...]]] = []

        def walk(node: str, protein_path: tuple[str, ...], edge_path: tuple[str, ...]) -> None:
            next_edges = outgoing.get(node, [])
            if not next_edges:
                if edge_path:
                    raw_paths.append((protein_path, edge_path))
                return
            for edge in next_edges:
                if edge.subject_protein_id in protein_path:
                    if edge_path:
                        raw_paths.append((protein_path, edge_path))
                    continue
                walk(
                    edge.subject_protein_id,
                    (*protein_path, edge.subject_protein_id),
                    (*edge_path, _edge_id(edge)),
                )

        for start_node in start_nodes:
            walk(start_node, (start_node,), ())

        deduped: dict[tuple[str, ...], tuple[str, ...]] = {}
        for protein_path, edge_path in raw_paths:
            current = deduped.get(protein_path)
            if current is None or edge_path < current:
                deduped[protein_path] = edge_path
        sorted_paths = sorted(
            deduped.items(),
            key=lambda item: (_path_sort_key(item[0], protein_map), item[1]),
        )
        protein_path_counts: dict[str, int] = {}
        for protein_path, _edge_path in sorted_paths:
            for protein_id in set(protein_path):
                protein_path_counts[protein_id] = protein_path_counts.get(protein_id, 0) + 1

        edge_path_id: dict[str, str] = {}
        paths: list[OrthologPath] = []
        for path_index, (protein_path, edge_path) in enumerate(sorted_paths, start=1):
            path_id = f"{group_id}.path_{path_index}"
            for edge_id in edge_path:
                edge_path_id.setdefault(edge_id, path_id)
            shared_protein_ids = tuple(
                sorted(
                    (
                        protein_id
                        for protein_id in protein_path
                        if protein_path_counts.get(protein_id, 0) > 1
                    ),
                    key=lambda protein_id: _protein_sort_key(protein_map[protein_id]),
                )
            )
            paths.append(
                OrthologPath(
                    orthogroup_id=group_id,
                    path_id=path_id,
                    protein_ids=tuple(protein_path),
                    edge_ids=tuple(edge_path),
                    shared_protein_ids=shared_protein_ids,
                )
            )

        updated_edges_by_group[group_id] = tuple(
            replace(edge, path_id=edge_path_id.get(_edge_id(edge), edge.path_id))
            for edge in edges
        )
        paths_by_group[group_id] = tuple(paths)
    return updated_edges_by_group, paths_by_group
