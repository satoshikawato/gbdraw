"""Skeleton (S0 probe)."""

from __future__ import annotations


def read_feature_sequences_json(request_json: str, resource_paths_json: str, output_directory: str, features_json: str) -> str:
    raise NotImplementedError


def read_match_sequences_json(request_json: str, resource_paths_json: str, output_directory: str, spans_json: str) -> str:
    raise NotImplementedError


def read_embed_sequences_json(request_json: str, resource_paths_json: str, output_directory: str, record_keys_json: str) -> str:
    raise NotImplementedError
