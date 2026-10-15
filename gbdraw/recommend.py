"""A first drawing proposed from the contents of the input records.

``gbdraw auto`` and the Web app can hand a GenBank file to :func:`propose`
and get back a complete, explained set of drawing settings: the mode, the
grouping, the palette and color rules, which feature types to draw, whether
and how to label, the strand layout, the quantitative tracks, the circular
track preset and the legend position. Every setting carries the reason it
was chosen, so the proposal can be shown before it is drawn and edited
afterwards.

The module reads records only. It never draws: a :class:`Proposal` is turned
into mode-command arguments (:meth:`Proposal.cli_args`) that the ``circular``
and ``linear`` commands render as usual, so the proposal is also a
reproducible command line.
"""

from __future__ import annotations

import json
import re
from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from importlib import resources
from pathlib import Path
from typing import Literal

from Bio.SeqRecord import SeqRecord

from gbdraw.core.record_metadata import infer_record_source_metadata
from gbdraw.features.shapes import DEFAULT_FEATURE_RENDERINGS
from gbdraw.mode_profiles import DEFAULT_FEATURE_TYPES, DiagramMode

GenomeKind = Literal[
    "mitochondrion",
    "plastid",
    "plasmid",
    "phage",
    "virus_dna",
    "virus_rna",
    "gene_cluster",
    "bacterial_chromosome",
    "contig",
    "eukaryotic_region",
    "unknown",
]
Grouping = Literal["single", "grid"]
ColorRule = tuple[str, str, str, str, str]
"""One specific-color rule: feature type, qualifier key, regex, color, caption."""

# Label density: how many CDS and RNA genes a drawing can label.
LABEL_ALL_MAX = 150
"""Up to this many labelable features, every one of them is labelled."""
LABEL_NAMED_MAX = 400
"""Up to this many, only features with a ``gene`` name and the rRNA/tRNA genes."""
GENE_NAME_COVERAGE_MIN = 0.6
"""Fraction of CDS with a ``gene`` qualifier above which gene names make the labels."""
OVERLAP_FRACTION_MIN = 0.2
"""Fraction of drawn features nested in another one above which lanes are split."""
HYPOTHETICAL_FRACTION_MIN = 0.05
"""Fraction of hypothetical CDS from which their labels are left out."""
HYPOTHETICAL_LABEL_BLACKLIST = ("hypothetical", "uncharacterized", "unknown")
"""The label keywords left out then (the packaged default, which ``--label_blacklist`` replaces)."""
COMPLETE_GENOME_MIN_LENGTH = 500_000
"""A prokaryotic record at least this long is treated as a chromosome."""
_RNA_TYPES = ("rRNA", "tRNA", "tmRNA", "ncRNA", "misc_RNA")
_HYPOTHETICAL = re.compile(r"hypothetical|uncharacteri[sz]ed|unknown", re.IGNORECASE)
_MITOCHONDRION = re.compile(r"\bmitochondri(on|al|a)\b", re.IGNORECASE)
_PLASTID = re.compile(r"\b(chloroplast|plastid|plastome)\b", re.IGNORECASE)
_PLASMID = re.compile(r"\bplasmid\b", re.IGNORECASE)
_PHAGE = re.compile(r"\bphage\b", re.IGNORECASE)
_GENE_CLUSTER = re.compile(r"gene cluster|biosynthetic|\bBGC\d{7}\b", re.IGNORECASE)
_COMPLETE = re.compile(r"complete (genome|sequence|plastome|chromosome)", re.IGNORECASE)
_TOOL_SIGNATURES: tuple[tuple[str, re.Pattern[str]], ...] = (
    ("bakta", re.compile(r"\bbakta\b", re.IGNORECASE)),
    ("pharokka", re.compile(r"\bpharokka\b", re.IGNORECASE)),
    ("prokka", re.compile(r"\bprokka\b", re.IGNORECASE)),
    ("dfast", re.compile(r"\bdfast\b", re.IGNORECASE)),
    ("pgap", re.compile(r"NCBI Prokaryotic Genome Annotation Pipeline|\bPGAP\b", re.IGNORECASE)),
    ("antismash", re.compile(r"\bantismash\b", re.IGNORECASE)),
    ("mitos", re.compile(r"\bmitos\b", re.IGNORECASE)),
)

# The auto color scheme. Fills carry black labels, so they sit at OKLCH
# lightness 0.55-0.90 with moderate chroma; the quantitative rings are muted so
# the features stay in front; the categories of one rule set were stepped in
# lightness as well as hue so that neighbors stay apart under deuteranopia and
# protanopia (OKLab distance >= 8 for the four mitochondrial complexes and the
# four antiSMASH kinds). ``-d`` rows override the palette; ``-t`` rows color
# the CDS that match a qualifier regex, searched case-insensitively.
BASE_DEFAULT_COLORS: tuple[tuple[str, str], ...] = (
    ("CDS", "#73b1e6"),
    ("rRNA", "#6ed9d2"),
    ("tRNA", "#e5d994"),
    ("tmRNA", "#bdd594"),
    ("ncRNA", "#a9e8c6"),
    ("misc_RNA", "#bddbb6"),
    ("repeat_region", "#e4e4e4"),
    ("misc_feature", "#dedede"),
    ("gene", "#73b1e6"),
    ("D-loop", "#e2d6c2"),
    ("rep_origin", "#625245"),
    ("oriT", "#625245"),
    ("mobile_element", "#d78e8d"),
    ("mat_peptide", "#c5b3ea"),
    ("gc_content", "#c8bca9"),
    ("skew_high", "#89cdc3"),
    ("skew_low", "#b89dd9"),
)
"""The palette overrides every proposal carries: soft blue CDS, sand GC ring, teal and lilac skew."""
# Overrides of the base for the kinds with a rule set: the CDS a rule does not
# match turn neutral so the categories stand out, and the plastome rings and
# RNA genes move away from its many greens and yellows.
_KIND_DEFAULT_COLORS: Mapping[str, tuple[tuple[str, str], ...]] = {
    "mitochondrion": (("CDS", "#c3d3e4"),),
    "plastid": (
        ("CDS", "#ded6c9"),
        ("rRNA", "#fda485"),
        ("tRNA", "#c19c7f"),
        ("repeat_region", "#ebebeb"),
        ("gc_content", "#c8b39c"),
        ("skew_high", "#96bbe3"),
        ("skew_low", "#dda5ba"),
    ),
    "antismash": (("CDS", "#e3d5c2"),),
}
PLASTOME_RULES: tuple[ColorRule, ...] = (
    ("CDS", "gene", r"^psa[A-Z]$", "#0b764d", "Photosystem I"),
    ("CDS", "gene", r"^psb[A-Z]$", "#6db453", "Photosystem II"),
    ("CDS", "gene", r"^pet[A-Z]$", "#1c989e", "Cytochrome b/f complex"),
    ("CDS", "gene", r"^atp[A-Z]$", "#c4bb5a", "ATP synthase"),
    ("CDS", "gene", r"^ndh[A-Z]$", "#f8dc86", "NADH dehydrogenase"),
    ("CDS", "gene", r"^rbcL$", "#79dbb2", "RubisCO large subunit"),
    ("CDS", "gene", r"^rpo[A-Z]\d?$", "#c04448", "RNA polymerase"),
    ("CDS", "gene", r"^rpl\d+$", "#406eb3", "Ribosomal proteins (LSU)"),
    ("CDS", "gene", r"^rps\d+$", "#90bff3", "Ribosomal proteins (SSU)"),
    ("CDS", "gene", r"^(clpP|matK)$", "#8d6cc2", "clpP, matK"),
    ("CDS", "gene", r"^(accD|ccsA|cemA|infA)$", "#e08bbc", "Other genes"),
    ("CDS", "gene", r"^ycf\d+$", "#cacaca", "Hypothetical reading frames (ycf)"),
)
MITOCHONDRION_RULES: tuple[ColorRule, ...] = (
    ("CDS", "gene", r"^(nd|nad)\d+l?$", "#6da9e7", "Complex I (NADH dehydrogenase)"),
    ("CDS", "gene", r"^(cytb|cob)$", "#fab45f", "Complex III (cytochrome b)"),
    ("CDS", "gene", r"^cox?\d$", "#c75374", "Complex IV (cytochrome c oxidase)"),
    ("CDS", "gene", r"^atp\d$", "#42a878", "ATP synthase"),
)
ANTISMASH_RULES: tuple[ColorRule, ...] = (
    ("CDS", "gene_kind", r"^biosynthetic$", "#c13a46", "Core biosynthetic"),
    ("CDS", "gene_kind", r"^biosynthetic-additional$", "#e88d94", "Additional biosynthetic"),
    ("CDS", "gene_kind", r"^transport$", "#5497d9", "Transport"),
    ("CDS", "gene_kind", r"^regulatory$", "#5dad70", "Regulatory"),
    ("CDS", "gene_kind", r"^other$", "#b7b7b7", "Other"),
)


@dataclass(frozen=True)
class RecordFacts:
    """What one record says about itself, before any drawing decision."""

    record_id: str
    length: int
    topology: str
    molecule_type: str
    organism: str
    strain: str
    organelle: str | None
    plasmid: str | None
    taxonomy: tuple[str, ...]
    definition: str
    annotation_tool: str
    feature_counts: Mapping[str, int]
    cds_qualifier_coverage: Mapping[str, float]
    hypothetical_fraction: float
    plus_strand_fraction: float
    gc_percent: float
    complete: bool
    kind: GenomeKind

    @property
    def labelable(self) -> int:
        """The CDS and RNA genes, which get labels when labels are drawn."""
        return sum(self.feature_counts.get(kind, 0) for kind in ("CDS", *_RNA_TYPES))


@dataclass(frozen=True)
class Decision:
    """One proposed setting and why it was chosen."""

    setting: str
    value: object
    reason: str


@dataclass
class Proposal:
    """A complete first drawing for a set of records, with its reasons."""

    facts: tuple[RecordFacts, ...]
    mode: DiagramMode = "circular"
    grouping: Grouping = "single"
    feature_types: tuple[str, ...] = DEFAULT_FEATURE_TYPES
    palette: str = "default"
    default_colors: tuple[tuple[str, str], ...] = BASE_DEFAULT_COLORS
    color_rules_name: str | None = None
    color_rules: tuple[ColorRule, ...] = ()
    label_scope: str = "none"
    qualifier_priority: tuple[tuple[str, str], ...] = ()
    label_whitelist: tuple[tuple[str, str, str], ...] = ()
    label_blacklist: tuple[str, ...] = ()
    separate_strands: bool = False
    resolve_overlaps: bool = False
    show_gc: bool = False
    show_skew: bool = False
    track_type: str = "tuckin"
    legend: str = "left"
    topology_override: str | None = None
    species: str | None = None
    decisions: list[Decision] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    def decide(self, setting: str, value: object, reason: str) -> None:
        setattr(self, setting, value)
        self.decisions.append(Decision(setting, value, reason))

    def tables(self) -> dict[str, tuple[tuple[str, ...], ...]]:
        """The header-less TSV tables the proposal needs, by table name."""
        tables: dict[str, tuple[tuple[str, ...], ...]] = {}
        if self.default_colors:
            tables["default_colors"] = self.default_colors
        if self.color_rules:
            tables["color_table"] = self.color_rules
        if self.qualifier_priority:
            tables["qualifier_priority"] = self.qualifier_priority
        if self.label_whitelist:
            tables["label_whitelist"] = self.label_whitelist
        return tables

    def write_tables(self, directory: Path, prefix: str) -> dict[str, Path]:
        """Write the tables next to the output as ``<prefix>.auto.<table>.tsv``."""
        written: dict[str, Path] = {}
        for name, rows in self.tables().items():
            path = directory / f"{prefix}.auto.{name}.tsv"
            path.write_text("".join("\t".join(row) + "\n" for row in rows), encoding="utf-8")
            written[name] = path
        return written

    def cli_args(self, tables: Mapping[str, Path | str] | None = None) -> list[str]:
        """The mode-command arguments after ``gbdraw <mode> --gbk FILE``.

        ``tables`` maps a table name from :meth:`tables` to the path it was
        written to; a table without a path is named ``<table>.tsv``.
        """
        paths = {name: str((tables or {}).get(name, f"{name}.tsv")) for name in self.tables()}
        args: list[str] = []
        if self.grouping == "grid":
            args.append("--multi_record_canvas")
        if self.palette != "default":
            args += ["-p", self.palette]
        if "default_colors" in paths:
            args += ["-d", paths["default_colors"]]
        if "color_table" in paths:
            args += ["-t", paths["color_table"]]
        if tuple(self.feature_types) != DEFAULT_FEATURE_TYPES:
            args += ["-k", ",".join(self.feature_types)]
        if self.mode == "circular":
            if self.label_scope != "none":
                args += ["--labels", "out" if self.label_scope == "outer" else self.label_scope]
            if self.track_type != "tuckin":
                args += ["--track_type", self.track_type]
            args += ["--gc" if self.show_gc else "--no-gc", "--skew" if self.show_skew else "--no-skew"]
        else:
            if self.label_scope != "none":
                args += ["--show_labels", self.label_scope]
            if self.show_gc:
                args.append("--gc")
            if self.show_skew:
                args.append("--skew")
        if self.separate_strands:
            args.append("--separate_strands")
        if self.resolve_overlaps:
            args.append("--resolve_overlaps")
        if "qualifier_priority" in paths:
            args += ["--qualifier_priority", paths["qualifier_priority"]]
        if "label_whitelist" in paths:
            args += ["--label_whitelist", paths["label_whitelist"]]
        elif self.label_blacklist:
            args += ["--label_blacklist", ",".join(self.label_blacklist)]
        if self.topology_override:
            args += ["--record_topology", self.topology_override]
        if self.species:
            args += ["--species", self.species]
        args += ["-l", self.legend]
        return args

    def describe(self) -> str:
        """The proposal as text: the facts it rests on, each decision, and the notes."""
        head = self.facts[0]
        kinds = Counter(facts.kind for facts in self.facts)
        counts = Counter[str]()
        for facts in self.facts:
            counts.update(facts.feature_counts)
        census = ", ".join(
            f"{kind} {count}"
            for kind, count in counts.most_common()
            if kind not in {"source", "gene"}
        )
        lines = [
            f"Input: {len(self.facts)} record(s), {sum(f.length for f in self.facts):,} bp, "
            + ", ".join(f"{kind} x{count}" for kind, count in kinds.most_common()),
            f"Organism: {head.organism or '(none)'} {head.strain}".rstrip()
            + f" [{' > '.join(head.taxonomy[:3]) or 'no taxonomy'}]; annotation: {head.annotation_tool}",
            f"Features: {census or 'none'} (hypothetical {head.hypothetical_fraction:.0%}, "
            f"gene names on {head.cds_qualifier_coverage.get('gene', 0.0):.0%} of CDS)",
            "",
            "Proposed settings:",
        ]
        width = max((len(decision.setting) for decision in self.decisions), default=0)
        for decision in self.decisions:
            lines.append(f"  {decision.setting:<{width}}  {decision.value!s:<22}  {decision.reason}")
        lines.extend(f"  note: {note}" for note in self.notes)
        return "\n".join(lines)

    def to_dict(self) -> dict[str, object]:
        """A JSON-compatible view for the Web app or scripts."""
        return {
            "mode": self.mode,
            "grouping": self.grouping,
            "featureTypes": list(self.feature_types),
            "palette": self.palette,
            "defaultColors": [list(row) for row in self.default_colors],
            "colorRulesName": self.color_rules_name,
            "colorRules": [list(rule) for rule in self.color_rules],
            "labelScope": self.label_scope,
            "qualifierPriority": [list(row) for row in self.qualifier_priority],
            "labelWhitelist": [list(row) for row in self.label_whitelist],
            "labelBlacklist": list(self.label_blacklist),
            "separateStrands": self.separate_strands,
            "resolveOverlaps": self.resolve_overlaps,
            "showGc": self.show_gc,
            "showSkew": self.show_skew,
            "trackType": self.track_type,
            "legend": self.legend,
            "topologyOverride": self.topology_override,
            "species": self.species,
            "decisions": [
                {"setting": d.setting, "value": d.value, "reason": d.reason} for d in self.decisions
            ],
            "notes": list(self.notes),
            "records": [
                {
                    "recordId": f.record_id,
                    "length": f.length,
                    "topology": f.topology,
                    "kind": f.kind,
                    "organism": f.organism,
                    "annotationTool": f.annotation_tool,
                    "featureCounts": dict(f.feature_counts),
                }
                for f in self.facts
            ],
        }

    def to_json(self) -> str:
        return json.dumps(self.to_dict(), indent=2, default=str)


def _qualifier_text(feature: object, key: str) -> list[str]:
    qualifiers = getattr(feature, "qualifiers", None) or {}
    return [str(value) for value in qualifiers.get(key, ())]


def _detect_annotation_tool(record: SeqRecord, cds_features: Sequence[object]) -> str:
    """Name the annotation tool from the qualifiers it writes, else from the comment."""
    cds_qualifiers = Counter[str]()
    for feature in cds_features:
        cds_qualifiers.update((getattr(feature, "qualifiers", None) or {}).keys())
    if cds_qualifiers["phrog"] or cds_qualifiers["vfdb_short_name"]:
        return "pharokka"
    if cds_qualifiers["gene_kind"] or cds_qualifiers["gene_functions"]:
        return "antismash"
    if any(
        value.startswith("COG:")
        for feature in cds_features[:200]
        for value in _qualifier_text(feature, "note")
    ):
        return "bakta"
    annotations = record.annotations or {}
    haystack = " ".join(
        (
            str(annotations.get("comment", "") or ""),
            json.dumps(annotations.get("structured_comment") or {}, default=str),
            " ".join(value for feature in cds_features[:50] for value in _qualifier_text(feature, "inference")),
        )
    )
    for name, pattern in _TOOL_SIGNATURES:
        if pattern.search(haystack):
            return name
    if "RefSeq" in (annotations.get("keywords") or ()):
        return "refseq"
    return "unknown"


def _gc_percent(record: SeqRecord) -> float:
    try:
        sequence = str(record.seq).upper()
    except Exception:  # an undefined (CON) sequence has no letters
        return 0.0
    if not sequence:
        return 0.0
    return 100.0 * (sequence.count("G") + sequence.count("C")) / len(sequence)


def _classify(
    *,
    definition: str,
    organism: str,
    organelle: str | None,
    plasmid: str | None,
    taxonomy: Iterable[str],
    molecule_type: str,
    topology: str,
    length: int,
    complete: bool,
    annotation_tool: str,
) -> GenomeKind:
    lineage = {name.lower() for name in taxonomy}
    organelle_text = (organelle or "").lower()
    if "mitochondri" in organelle_text or _MITOCHONDRION.search(definition):
        return "mitochondrion"
    if any(word in organelle_text for word in ("plastid", "chloroplast")) or _PLASTID.search(definition):
        return "plastid"
    if plasmid or _PLASMID.search(definition):
        return "plasmid"
    if "viruses" in lineage or "viroids" in lineage or "virus" in organism.lower():
        if "rna" in molecule_type.lower() or "riboviria" in lineage:
            return "virus_rna"
        if _PHAGE.search(f"{definition} {organism}") or lineage & {"caudoviricetes", "duplodnaviria", "uroviricota"}:
            return "phage"
        return "virus_dna"
    if annotation_tool == "antismash" or _GENE_CLUSTER.search(definition):
        return "gene_cluster"
    if lineage & {"bacteria", "archaea"}:
        if complete or topology == "circular" or length >= COMPLETE_GENOME_MIN_LENGTH:
            return "bacterial_chromosome"
        return "contig"
    if "eukaryota" in lineage:
        return "eukaryotic_region"
    if topology == "circular" and length >= COMPLETE_GENOME_MIN_LENGTH:
        return "bacterial_chromosome"
    return "unknown"


def inspect_record(record: SeqRecord) -> RecordFacts:
    """Read the facts of one record and classify it."""
    annotations = record.annotations or {}
    metadata = infer_record_source_metadata(record)
    features = list(record.features)
    counts = Counter[str](str(feature.type) for feature in features)
    cds = [feature for feature in features if feature.type == "CDS"]
    qualifier_counts = Counter[str]()
    for feature in cds:
        qualifier_counts.update((feature.qualifiers or {}).keys())
    coverage = {key: count / len(cds) for key, count in qualifier_counts.items()} if cds else {}
    hypothetical = (
        sum(1 for feature in cds if any(_HYPOTHETICAL.search(p) for p in _qualifier_text(feature, "product")))
        / len(cds)
        if cds
        else 0.0
    )
    stranded = [
        feature
        for feature in features
        if feature.type in {"CDS", "rRNA", "tRNA"} and getattr(feature, "location", None) is not None
    ]
    plus = sum(1 for feature in stranded if (feature.location.strand or 1) > 0)
    plasmid = next(
        (values[0] for feature in features if feature.type == "source" for values in [_qualifier_text(feature, "plasmid")] if values),
        None,
    )
    definition = str(record.description or "")
    taxonomy = tuple(str(name) for name in annotations.get("taxonomy", ()) or ())
    topology = str(annotations.get("topology") or "unknown")
    molecule_type = str(annotations.get("molecule_type") or "")
    complete = bool(_COMPLETE.search(definition))
    tool = _detect_annotation_tool(record, cds)
    return RecordFacts(
        record_id=str(record.id),
        length=len(record),
        topology=topology,
        molecule_type=molecule_type,
        organism=metadata.organism,
        strain=metadata.strain,
        organelle=metadata.organelle,
        plasmid=plasmid,
        taxonomy=taxonomy,
        definition=definition,
        annotation_tool=tool,
        feature_counts=dict(counts),
        cds_qualifier_coverage=coverage,
        hypothetical_fraction=hypothetical,
        plus_strand_fraction=plus / len(stranded) if stranded else 0.5,
        gc_percent=_gc_percent(record),
        complete=complete,
        kind=_classify(
            definition=definition,
            organism=metadata.organism,
            organelle=metadata.organelle,
            plasmid=plasmid,
            taxonomy=taxonomy,
            molecule_type=molecule_type,
            topology=topology,
            length=len(record),
            complete=complete,
            annotation_tool=tool,
        ),
    )


def _overlap_fraction(records: Sequence[SeqRecord], types: set[str]) -> float:
    """The fraction of drawn features with a part that lies over half inside another feature.

    Location parts, not whole locations, are compared: a feature joined across
    the origin of a circular record would otherwise span the whole record.
    """
    total = nested = 0
    for record in records:
        parts: list[tuple[int, int, int]] = []
        drawn = 0
        for index, feature in enumerate(record.features):
            location = getattr(feature, "location", None)
            if feature.type not in types or location is None:
                continue
            drawn += 1
            parts.extend((int(part.start), int(part.end), index) for part in location.parts)
        parts.sort()
        total += drawn
        nested_features: set[int] = set()
        for index, (start, end, owner) in enumerate(parts):
            length = max(1, end - start)
            for other_start, other_end, other_owner in parts[max(0, index - 32) : index + 32]:
                if other_owner == owner:
                    continue
                shared = min(end, other_end) - max(start, other_start)
                if shared > 0.5 * min(length, max(1, other_end - other_start)):
                    nested_features.add(owner)
                    break
        nested += len(nested_features)
    return nested / total if total else 0.0


def _preset_rules(name: str) -> tuple[ColorRule, ...]:
    """The rows of a color table that ships with the Web app (Bakta COG, Pharokka).

    The browser wheel does not carry the preset files; there the rules are empty
    and the Web app applies its own copy of the table.
    """
    presets = resources.files("gbdraw").joinpath("web").joinpath("presets")
    try:
        text = presets.joinpath(f"{name}_color_table.txt").read_text(encoding="utf-8")
    except OSError:
        return ()
    rows: list[ColorRule] = []
    for line in text.splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 5 and not line.startswith("#"):
            rows.append((parts[0], parts[1], parts[2], parts[3], parts[4]))
    return tuple(rows)


def propose(records: Sequence[SeqRecord]) -> Proposal:
    """Propose a first drawing of ``records`` from what they contain."""
    if not records:
        raise ValueError("propose needs at least one record")
    facts = tuple(inspect_record(record) for record in records)
    head = facts[0]
    kinds = {f.kind for f in facts}
    circular_kinds = {"mitochondrion", "plastid", "plasmid", "bacterial_chromosome"}
    count = len(facts)
    proposal = Proposal(facts=facts)

    # Mode and grouping.
    if count == 1:
        if head.topology == "circular":
            proposal.decide("mode", "circular", "the LOCUS line says circular")
        elif head.kind in circular_kinds and head.complete:
            proposal.decide(
                "mode", "circular", f"the LOCUS line says linear, but a complete {head.kind} genome is circular"
            )
            proposal.decide("topology_override", "circular", "treat the record as circular (--record_topology)")
        else:
            proposal.decide("mode", "linear", f"a {head.kind} record with a linear or unknown topology")
    elif {f.topology for f in facts} == {"circular"} and kinds <= circular_kinds | {"phage", "virus_dna"}:
        proposal.decide("mode", "circular", f"all {count} records are circular (chromosomes and plasmids)")
        proposal.decide("grouping", "grid", "one canvas with every record (--multi_record_canvas)")
    else:
        proposal.decide("mode", "linear", f"{count} records, not all circular, stacked as rows")
        if count > 30:
            proposal.notes.append(
                f"{count} records is a lot for one figure; a records table can select the longest ones"
            )

    # Feature types.
    counts = Counter[str]()
    for f in facts:
        counts.update(f.feature_counts)
    types = list(DEFAULT_FEATURE_TYPES)
    if counts["CDS"] == 0 and counts["gene"] > 0:
        types.append("gene")
        proposal.notes.append("no CDS features; the gene features are drawn instead")
    for kind, reason in (
        ("rep_origin", "replication origins are annotated"),
        ("D-loop", "the mitochondrial control region is annotated"),
        ("mobile_element", "mobile elements are annotated"),
        ("oriT", "oriT is annotated"),
    ):
        if counts[kind]:
            types.append(kind)
            proposal.decisions.append(Decision("feature_types", f"+{kind}", reason))
    mature = counts["mat_peptide"]
    if mature and mature >= max(3, counts["CDS"] // 2):
        types.append("mat_peptide")
        proposal.decisions.append(
            Decision(
                "feature_types",
                "+mat_peptide",
                f"{mature} mature peptides against {counts['CDS']} CDS: a polyprotein genome",
            )
        )
    elif mature:
        proposal.notes.append(f"{mature} mat_peptide feature(s) are not drawn; add mat_peptide to -k to draw them")
    if 0 < counts["regulatory"] <= 30:
        proposal.notes.append(f"{counts['regulatory']} regulatory feature(s) can be drawn by adding regulatory to -k")
    proposal.feature_types = tuple(types)
    # Underlays (repeat_region by default) are backdrops that genes lie on by design.
    overlap = _overlap_fraction(
        records, {kind for kind in types if DEFAULT_FEATURE_RENDERINGS.get(kind) != "underlay"}
    )
    if overlap >= OVERLAP_FRACTION_MIN:
        proposal.decide(
            "resolve_overlaps", True, f"{overlap:.0%} of the drawn features lie inside another one (--resolve_overlaps)"
        )

    # Labels.
    labelable = sum(f.labelable for f in facts)
    gene_coverage = head.cds_qualifier_coverage.get("gene", 0.0)
    if proposal.mode == "circular":
        if proposal.grouping == "single" and labelable <= LABEL_ALL_MAX:
            proposal.decide("label_scope", "outer", f"{labelable} labels fit outside the circle")
            proposal.decide("track_type", "middle", "the middle preset reads best with outer labels")
        elif proposal.grouping == "single" and labelable <= LABEL_NAMED_MAX:
            proposal.decide("label_scope", "outer", f"{labelable} candidates: only named genes and RNA genes are labelled")
            proposal.label_whitelist = (("CDS", "gene", ".+"), ("rRNA", "product", ".+"), ("tRNA", "product", ".+"))
        else:
            proposal.decide("label_scope", "none", f"{labelable} labels would not be readable")
    elif count == 1 and labelable <= LABEL_ALL_MAX:
        proposal.decide("label_scope", "all", f"one record with {labelable} labels")
    elif count > 1 and labelable <= LABEL_ALL_MAX * count:
        proposal.decide("label_scope", "first", f"{count} records: only the first row is labelled")
    else:
        proposal.decide("label_scope", "none", f"{labelable // count} labels per record would not be readable")
    if proposal.label_scope != "none":
        if gene_coverage >= GENE_NAME_COVERAGE_MIN:
            priority = "gene,product"
            reason = f"gene names on {gene_coverage:.0%} of CDS are shorter than products"
        else:
            priority = "product,gene,locus_tag"
            reason = f"gene names on only {gene_coverage:.0%} of CDS, so products label the CDS"
        rows = [("CDS", priority)]
        if "mat_peptide" in proposal.feature_types:
            rows.append(("mat_peptide", "product,gene"))
        proposal.decide("qualifier_priority", priority, reason)
        proposal.qualifier_priority = tuple(rows)
        if not proposal.label_whitelist and head.hypothetical_fraction >= HYPOTHETICAL_FRACTION_MIN:
            proposal.label_blacklist = HYPOTHETICAL_LABEL_BLACKLIST
            proposal.decisions.append(
                Decision(
                    "label_blacklist",
                    ",".join(HYPOTHETICAL_LABEL_BLACKLIST),
                    f"{head.hypothetical_fraction:.0%} of CDS are hypothetical; their labels say nothing",
                )
            )

    # Colors: the base scheme, then the rule set the annotation allows.
    tool = head.annotation_tool
    preset = _preset_rules(tool) if tool in {"bakta", "pharokka"} else ()
    scheme = "base"
    if tool == "bakta" and preset:
        proposal.color_rules_name, proposal.color_rules = "bakta_cog", preset
        proposal.decide("colors", "bakta COG", "Bakta COG categories (note COG:X) color the CDS")
    elif tool == "pharokka" and preset:
        proposal.color_rules_name, proposal.color_rules = "pharokka", preset
        proposal.decide("colors", "pharokka", "Pharokka functions (function, phrog) color the CDS")
    elif tool == "antismash":
        scheme = "antismash"
        proposal.color_rules_name, proposal.color_rules = "antismash_gene_kind", ANTISMASH_RULES
        proposal.decide("colors", "antiSMASH gene_kind", "core and additional biosynthetic, transport, regulatory and other CDS")
    elif head.kind == "plastid" and gene_coverage >= 0.5:
        scheme = "plastid"
        proposal.color_rules_name, proposal.color_rules = "plastome_function", PLASTOME_RULES
        proposal.decide("colors", "plastome families", "psa, psb, pet, atp, ndh, rbcL, rpo, rpl, rps, clpP/matK and ycf gene names")
    elif head.kind == "mitochondrion" and gene_coverage >= 0.5:
        scheme = "mitochondrion"
        proposal.color_rules_name, proposal.color_rules = "mitochondrion_complexes", MITOCHONDRION_RULES
        proposal.decide("colors", "respiratory complexes", "nad, cob, cox and atp gene names")
    else:
        proposal.decide("colors", "base scheme", "no functional annotation to color by; soft blue CDS")
    colors = dict(BASE_DEFAULT_COLORS)
    colors.update(_KIND_DEFAULT_COLORS.get(scheme, ()))
    proposal.default_colors = tuple(colors.items())

    # Strands, GC content, skew.
    plus = head.plus_strand_fraction
    if plus >= 0.95 or plus <= 0.05:
        proposal.decide("separate_strands", False, f"{max(plus, 1 - plus):.0%} of the genes are on one strand")
    elif head.kind == "gene_cluster":
        proposal.decide("separate_strands", False, "a gene cluster is more compact in one lane")
    else:
        proposal.decide("separate_strands", True, f"genes on both strands ({plus:.0%} on the plus strand)")
    if proposal.mode == "circular":
        proposal.decide("show_gc", True, "a circular diagram shows the GC content ring")
        if head.kind in {"bacterial_chromosome", "plasmid"}:
            proposal.decide("show_skew", True, "the GC skew sign change marks the replication origin and terminus")
        else:
            proposal.decide("show_skew", True, "the circular default")
        if head.kind == "mitochondrion" and head.gc_percent < 45:
            proposal.notes.append(f"GC {head.gc_percent:.1f}%: an AT skew (-n AT) may read better")
    else:
        proposal.decide("show_gc", False, "the linear default (add --gc for the GC content track)")
        proposal.decide("show_skew", False, "the linear default")

    # Legend and title.
    captions = len({rule[4] for rule in proposal.color_rules})
    if proposal.mode == "circular":
        proposal.decide("legend", "left" if captions <= 12 else "right", "the legend beside the circle")
    else:
        proposal.decide("legend", "bottom", "the legend under the rows")
    if not head.organism or head.organism.lower() in {"synthetic construct", "unknown", "unidentified"}:
        proposal.decide("species", head.record_id, "no organism name, so the record ID titles the drawing")
    if head.feature_counts.get("CDS", 0) == 0 and "gene" not in proposal.feature_types:
        proposal.notes.append("no CDS or gene features: the drawing shows the quantitative tracks only")
    if any(f.length > 20_000_000 for f in facts):
        proposal.notes.append("over 20 Mb: a region (--region) draws a readable part of it")
    return proposal


__all__ = [
    "ANTISMASH_RULES",
    "ColorRule",
    "Decision",
    "GenomeKind",
    "Grouping",
    "HYPOTHETICAL_FRACTION_MIN",
    "HYPOTHETICAL_LABEL_BLACKLIST",
    "LABEL_ALL_MAX",
    "LABEL_NAMED_MAX",
    "MITOCHONDRION_RULES",
    "PLASTOME_RULES",
    "Proposal",
    "RecordFacts",
    "inspect_record",
    "propose",
]
