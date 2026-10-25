from __future__ import annotations

import re
from collections.abc import Iterable, Iterator
from xml.sax.saxutils import escape
from typing import TYPE_CHECKING, NamedTuple, TypeVar

from Bio.SeqFeature import ExactPosition, SeqFeature
from Bio.SeqRecord import SeqRecord

from gbdraw.core.text import parse_mixed_content_text

if TYPE_CHECKING:
    import numpy as np
    from numpy.typing import NDArray


_COORD_BASE_KEY = "gbdraw_coord_base"
_COORD_STEP_KEY = "gbdraw_coord_step"
_SOURCE_FEATURE_INDEX_ATTR = "_gbdraw_source_feature_index"
_SOURCE_FEATURE_PARTS_ATTR = "_gbdraw_source_feature_location_parts"
_SOURCE_FEATURE_ANCHOR_PROFILE_ATTR = "_gbdraw_source_feature_anchor_profile"
_SOURCE_FEATURE_HASH_ATTR = "_gbdraw_source_feature_hash"
_SourceBase = TypeVar("_SourceBase", int, "NDArray[np.int64]")


class _SourceFeatureAnchorProfile(NamedTuple):
    precision: str
    operator: str
    part_order: str
    strand: str


class _SourceFeatureHashStamp(NamedTuple):
    """The source hash of a feature, computed once for the record it was loaded in.

    The hash also depends on the feature's type and location, so the stamp keeps
    the type and the location object it was computed for; a reader accepts it
    only for that same type and location object.
    """

    record_id: str | None
    feature_type: str
    location: object
    feature_hash: str


def _iter_source_features(features: Iterable[SeqFeature] | None) -> Iterator[SeqFeature]:
    """Yield top-level and nested source features in their original order."""
    for feature in features or ():
        yield feature
        yield from _iter_source_features(getattr(feature, "sub_features", None))


def _feature_source_index_map(features: Iterable[SeqFeature] | None) -> dict[int, int]:
    """Map nested feature objects to their flattened source-order ordinals."""
    return {id(feature): index for index, feature in enumerate(_iter_source_features(features))}


def _source_feature_index(feature: object) -> int | None:
    """Return a pre-transform source ordinal attached to a feature, if any."""

    try:
        index = int(getattr(feature, _SOURCE_FEATURE_INDEX_ATTR))
    except (AttributeError, TypeError, ValueError):
        return None
    return index if index >= 0 else None


def _source_feature_location_parts(
    feature: object,
) -> tuple[tuple[int, int, int | None], ...] | None:
    """Return pre-transform biological location parts, if attached."""

    value = getattr(feature, _SOURCE_FEATURE_PARTS_ATTR, None)
    if not isinstance(value, tuple):
        return None
    parts: list[tuple[int, int, int | None]] = []
    for part in value:
        if not isinstance(part, tuple) or len(part) != 3:
            return None
        try:
            start = int(part[0])
            end = int(part[1])
        except (TypeError, ValueError):
            return None
        strand = part[2]
        if strand not in {-1, 1}:
            strand = None
        parts.append((start, end, strand))
    return tuple(parts) or None


def _location_anchor_profile(
    location: object,
    *,
    coord_step: int,
) -> _SourceFeatureAnchorProfile:
    """Normalize only source-location facts needed by Web anchor resolution."""

    raw_parts = list(getattr(location, "parts", None) or [location])
    parts = [part for part in raw_parts if part is not None]
    exact = bool(parts) and all(
        isinstance(getattr(part, "start", None), ExactPosition)
        and isinstance(getattr(part, "end", None), ExactPosition)
        for part in parts
    )
    if len(parts) == 1:
        operator = "single"
    else:
        raw_operator = str(getattr(location, "operator", "") or "")
        operator = raw_operator if raw_operator in {"join", "order"} else "unknown"

    direction = 1 if int(coord_step) >= 0 else -1
    strands = [
        int(part.strand) * direction if getattr(part, "strand", None) in {-1, 1} else None
        for part in parts
    ]
    if strands and all(strand == 1 for strand in strands):
        strand = "+"
    elif strands and all(strand == -1 for strand in strands):
        strand = "-"
    elif strands and all(strand is None for strand in strands):
        strand = "unstranded"
    else:
        strand = "mixed"

    if operator == "single" and strand in {"+", "-"}:
        part_order = "biological"
    elif operator == "single" and strand == "unstranded":
        part_order = "source-forward"
    elif operator == "join" and strand in {"+", "-"}:
        part_order = "biological"
    elif operator == "join" and strand == "unstranded":
        part_order = "source-forward"
    else:
        part_order = "ambiguous"

    return _SourceFeatureAnchorProfile(
        precision="exact" if exact else "fuzzy",
        operator=operator,
        part_order=part_order,
        strand=strand,
    )


def _source_feature_anchor_profile(
    feature: object,
    *,
    coord_step: int = 1,
) -> _SourceFeatureAnchorProfile:
    """Return pre-transform source facts, computing them before first copy."""

    stored = getattr(feature, _SOURCE_FEATURE_ANCHOR_PROFILE_ATTR, None)
    if isinstance(stored, _SourceFeatureAnchorProfile):
        return stored
    return _location_anchor_profile(
        getattr(feature, "location", None),
        coord_step=coord_step,
    )


def _stamp_source_feature_hash(
    feature: object, record_id: str | None, feature_hash: str | None,
) -> None:
    """Store the source hash computed when the feature's source was loaded."""

    if feature_hash is None:
        if hasattr(feature, _SOURCE_FEATURE_HASH_ATTR):
            delattr(feature, _SOURCE_FEATURE_HASH_ATTR)
        return
    setattr(
        feature,
        _SOURCE_FEATURE_HASH_ATTR,
        _SourceFeatureHashStamp(
            record_id, str(getattr(feature, "type", "") or ""), getattr(feature, "location", None), feature_hash,
        ),
    )


def _feature_hash_stamp(feature: object) -> _SourceFeatureHashStamp | None:
    """The stamp when it was computed for the feature's current type and location object."""

    stamp = getattr(feature, _SOURCE_FEATURE_HASH_ATTR, None)
    if (
        isinstance(stamp, _SourceFeatureHashStamp)
        and stamp.feature_type == str(getattr(feature, "type", "") or "")
        and stamp.location is getattr(feature, "location", None)
    ):
        return stamp
    return None


def _stamped_source_feature_hash(feature: object, record_id: str | None) -> str | None:
    """The stored source hash when it was computed for ``record_id`` and the feature as it is, else None."""

    stamp = _feature_hash_stamp(feature)
    return stamp.feature_hash if stamp is not None and stamp.record_id == record_id else None


def _mapped_feature_location_parts(
    feature: object,
    *,
    coord_base: int,
    coord_step: int,
) -> tuple[tuple[int, int, int | None], ...]:
    """Map a feature's current parts into its biological source coordinates."""

    location = getattr(feature, "location", None)
    raw_parts = list(getattr(location, "parts", None) or [location])
    parts: list[tuple[int, int, int | None]] = []
    for part in raw_parts:
        if part is None:
            continue
        try:
            start, end = _absolute_display_interval(
                int(part.start),
                int(part.end),
                coord_base,
                coord_step,
            )
        except (AttributeError, TypeError, ValueError):
            continue
        strand = getattr(part, "strand", None)
        if strand in {-1, 1}:
            strand = int(strand) * (1 if int(coord_step) >= 0 else -1)
        else:
            strand = None
        parts.append((start, end, strand))
    return tuple(parts)


def _copy_source_feature_identity(
    source: object,
    target: object,
    *,
    fallback_index: int,
    coord_base: int,
    coord_step: int,
) -> None:
    """Carry source ordinal, biological parts and source hash across a transformation."""

    index = _source_feature_index(source)
    setattr(
        target,
        _SOURCE_FEATURE_INDEX_ATTR,
        int(fallback_index) if index is None else index,
    )
    parts = _source_feature_location_parts(source)
    if parts is None:
        parts = _mapped_feature_location_parts(
            source,
            coord_base=coord_base,
            coord_step=coord_step,
        )
    if parts:
        setattr(target, _SOURCE_FEATURE_PARTS_ATTR, parts)
    setattr(
        target,
        _SOURCE_FEATURE_ANCHOR_PROFILE_ATTR,
        _source_feature_anchor_profile(source, coord_step=coord_step),
    )
    stamp = _feature_hash_stamp(source)
    if stamp is not None:
        _stamp_source_feature_hash(target, stamp.record_id, stamp.feature_hash)


def _read_coord_map(record: object) -> tuple[int, int]:
    annotations = getattr(record, "annotations", None) or {}
    try:
        base = int(annotations.get(_COORD_BASE_KEY, 1))
    except (TypeError, ValueError):
        base = 1
    try:
        step = int(annotations.get(_COORD_STEP_KEY, 1))
    except (TypeError, ValueError):
        step = 1
    if step == 0:
        step = 1
    return base, (1 if step > 0 else -1)


def _source_to_local_index(
    source: _SourceBase, coord_base: int, coord_step: int
) -> _SourceBase:
    """Return the 0-based local index of a 1-based source base.

    This inverts the record coordinate map ``source = base + step * local``.
    A numpy integer array of source bases is mapped elementwise. The result is
    outside ``0..len(record) - 1`` when the base is outside a cropped record.
    """

    return (source - coord_base) * coord_step


def _record_source_start(record: SeqRecord) -> int:
    """The 1-based source position of a record's first base on its source strand.

    It is where the record's search frame starts (PD-OI-073): the crop start of
    a cropped record, also after a reverse complement, and 1 for a whole record.
    """

    coord_base, coord_step = _read_coord_map(record)
    return coord_base if coord_step == 1 else coord_base - len(record.seq) + 1


def _write_coord_map(record: SeqRecord, *, base: int, step: int) -> None:
    if getattr(record, "annotations", None) is None:
        record.annotations = {}
    record.annotations[_COORD_BASE_KEY] = int(base)
    record.annotations[_COORD_STEP_KEY] = 1 if int(step) >= 0 else -1


def _absolute_display_interval(
    start: int,
    end: int,
    coord_base: int,
    coord_step: int,
) -> tuple[int, int]:
    if end <= start:
        coord = coord_base + (coord_step * start)
        return coord - 1, coord
    first_coord = coord_base + (coord_step * start)
    last_coord = coord_base + (coord_step * (end - 1))
    return min(first_coord, last_coord) - 1, max(first_coord, last_coord)


class RecordSourceMetadata(NamedTuple):
    organism: str
    strain: str
    replicon: str | None
    organelle: str | None
    # The first word of the record's lineage (GenBank ORGANISM lines), if any.
    lineage_top: str | None = None


def infer_record_source_metadata(record: SeqRecord) -> RecordSourceMetadata:
    """Extract source-feature metadata used in definition labels."""
    annotations = getattr(record, "annotations", None) or {}
    organism = str(annotations.get("organism", "") or "").strip()
    taxonomy = annotations.get("taxonomy") or []
    lineage_top = (str(taxonomy[0]).strip() or None) if taxonomy else None
    strain = ""
    replicon: str | None = None
    organelle: str | None = None

    for feature in getattr(record, "features", []):
        if getattr(feature, "type", None) != "source":
            continue

        qualifiers = getattr(feature, "qualifiers", {}) or {}
        if "organism" in qualifiers and qualifiers["organism"]:
            organism = str(qualifiers["organism"][0]).strip()
        if "isolate" in qualifiers and qualifiers["isolate"]:
            strain = str(qualifiers["isolate"][0]).strip()
        elif "strain" in qualifiers and qualifiers["strain"]:
            strain = str(qualifiers["strain"][0]).strip()

        if "chromosome" in qualifiers and qualifiers["chromosome"]:
            replicon = f"Chromosome {str(qualifiers['chromosome'][0]).strip()}"
        elif "plasmid" in qualifiers and qualifiers["plasmid"]:
            replicon = str(qualifiers["plasmid"][0]).strip()

        if "organelle" in qualifiers and qualifiers["organelle"]:
            organelle = str(qualifiers["organelle"][0]).strip()

    return RecordSourceMetadata(
        organism=organism,
        strain=strain,
        replicon=replicon,
        organelle=organelle,
        lineage_top=lineage_top,
    )


class OrganismLines(NamedTuple):
    """A record's name as a taxon line and a strain line (text with ``<i>``)."""

    taxon: str
    strain: str


_NON_ORGANISM_NAMES = frozenset(
    {
        "synthetic construct",
        "artificial sequence",
        "unidentified",
        "unidentified organism",
        "unknown",
        "unknown organism",
        "vector",
        "cloning vector",
        "unidentified cloning vector",
        "expression vector",
    }
)
_CANDIDATUS_PREFIX = re.compile(r"(?:candidatus|ca\.)\s+", re.IGNORECASE)
_NON_ORGANISM_START = re.compile(
    r"(?:synthetic|artificial|cloning|expression|shuttle|transformation|binary|gene\s+trap|vector|"
    r"unidentified|unknown|unclassified|uncultured|environmental)\b",
    re.IGNORECASE,
)
_SEQUENCE_WORDS = frozenset(
    {"dna", "rna", "molecule", "vector", "construct", "sequence", "sequences", "episomal", "replicon",
     "transposon", "integron", "plasmid"}
)
# A virus word or an ICTV rank ending at the end of a word or before "-" or "["
# ("virus-associated"), never inside a Latin epithet ("Treponema phagedenis").
_VIRUS_WORD = re.compile(
    r"(?:virus|viruses|viroid|viroids|viriform|satellite|satellites|phage|phages|virinae|viridae|virales|"
    r"viricetes|viricetidae|viricota|viricotina|virineae|virites|virae|vira|viria|virid|virids)(?![a-z])",
    re.IGNORECASE,
)
_VIRUS_ACRONYM = re.compile(r"[A-Z]{0,3}NPV")
_ICTV_GENUS = re.compile(r"[A-Z][a-z]*(?:virus|viroid|viriform|satellite)")
_CELLULAR_LINEAGES = frozenset({"Bacteria", "Archaea", "Eukaryota"})
_VIRUS_LINEAGE = "Viruses"
_HYBRID = re.compile(r"(?:^|\s)[x×](?:\s|$)|×")
_GENUS = re.compile(r"[A-Z][a-z]+")
_BRACKETED_GENUS = re.compile(r"\[([A-Z][a-z]+)\]")
# Family, subfamily and class endings: such a name is no genus, and ICZN sets it upright.
_ABOVE_GENUS = re.compile(r"[A-Z][a-z]*(?:idae|inae|aceae|ineae|phyceae|mycetes|mycota|opsida)")
_EPITHET = re.compile(r"[a-z]{2,}(?:-[a-z]+)*")
# Rank words before a Latin infraspecific epithet; upright, part of the taxon line.
_LATIN_RANK_WORDS = frozenset({"subsp.", "subsp", "ssp.", "var.", "subvar.", "f.", "forma", "nothosubsp.", "nothovar."})
# Infrasubspecific ranks and designations: they start the strain line.
_STRAIN_LINE_WORDS = frozenset(
    {"str.", "str", "strain", "substr.", "substr", "serovar", "sv.", "ser.", "serotype", "serogroup", "pv.",
     "pathovar", "bv.", "biovar", "biotype", "genomovar", "gv.", "morphovar", "chemovar", "phagovar", "ecotype",
     "cv.", "cultivar", "type", "subtype", "genotype", "isolate", "breed", "race", "pathotype", "phylotype",
     "ribotype", "morphotype", "variant", "mutant", "assemblage"}
)
_NOT_EPITHETS = frozenset(
    {"sp", "spp", "bacterium", "bacteria", "archaeon", "endosymbiont", "symbiont", "symbiote", "phytoplasma",
     "group", "complex", "clade", "subgroup", "cluster", "lineage", "genomosp", "genosp", "genomospecies",
     "genospecies", "x", "hybrid", "vector", "agent", "organism", "eukaryote", "cyanobacterium",
     "proteobacterium", "of", "and", "like", "cf", "aff", "nr", "gen", "nov", "sect", "subsect", "subgen",
     "environmental", "uncultured", "unidentified", "unclassified", "metagenome", "mixed", "culture", "sequence",
     "plasmid", "chromosome", "segment", "satellite", "virus", "phage", "clone", "oral", "cytoform",
     "genogroup", "form", "incertae", "sedis", "ex", "sensu", "hybrids", "icg", "ectomycorrhiza", "archeaon",
     "subcluster", "from", "in", "on", "at", "near", "with", "the", "for", "by", "to", "taxon", "related",
     "species", "subspecies", "unnamed", "fungal", "fungus", "contaminant", "mycorrhizal"}
    | {word.rstrip(".") for word in _STRAIN_LINE_WORDS}
)
# The last part of a word ("S-endosymbiont") that names an organism kind; the
# Latin words before it name a host.
_ORGANISM_NOUN = re.compile(r"(?:phytoplasma|biont|symbiote|bacterium|archaeon|archaeote|archaeum|parasite)$")
# NCBI status and homonym suffixes ("(SeqCode)", "(nom. nud.)", "(in: starfish)").
_STATUS_SUFFIX = re.compile(r"\((?:SeqCode|nom\. [a-z.]+|nomen [a-z]+|in: [^()]+|ex .*\b\d{4})\)")
_DESIGNATION_PART = re.compile(r"[-_/:=;,()\[\]]+")
# A lowercase word after an unranked epithet ("tomato", "sp."): it ends the taxon there.
_LOWERCASE_WORD = re.compile(r"[a-z]+(?:-[a-z]+)*\.?")
# A strain that is one lowercase word ("pneumoniae") names an epithet, not a boundary.
_LOWERCASE_STRAIN = re.compile(r"[a-z]{3,}")


def _letters(text: str) -> str:
    return re.sub(r"[^0-9a-z]+", "", text.lower())


def _names_strain(tokens: list[str], strain: str) -> int | None:
    """Index of the first word where the name spells the strain, else None.

    Letters and digits only: "So ce 56" is named by "So ce56", "Q6390" by
    "Marseille-Q6390", but "7" is not named by "O157:H7".
    """
    target = _letters(strain)
    if not target:
        return None
    for index, token in enumerate(tokens):
        if target in {_letters(part) for part in _DESIGNATION_PART.split(token) if part}:
            return index
        joined = ""
        for later in tokens[index:]:
            joined += _letters(later)
            if joined == target:
                return index
            if not target.startswith(joined):
                break
    return None


def _strain_line(rest: list[str], strain: str, name_tokens: list[str]) -> str:
    """The words after the taxon, then the strain unless the name holds it (OV-1830)."""
    words = " ".join(rest)
    if not _letters(strain) or _names_strain(name_tokens, strain) is not None:
        return words
    parts = [part for part in re.split(r"\s*;\s*", strain) if part]
    missing = [part for part in parts if _names_strain(name_tokens, part) is None]
    return f"{words} {'; '.join(missing)}".strip()


def _is_epithet(token: str) -> bool:
    # The last part also counts: "funestus-like" is not an epithet.
    last = token.rsplit("-", 1)[-1]
    return bool(_EPITHET.fullmatch(token)) and last not in _NOT_EPITHETS and not _ORGANISM_NOUN.search(last)


def _parenthetical_end(tokens: list[str], start: int) -> int:
    if start >= len(tokens) or not tokens[start].startswith("("):
        return start
    depth = 0
    for end in range(start, len(tokens)):
        depth += tokens[end].count("(") - tokens[end].count(")")
        if depth <= 0:
            return end + 1
    return start


def _name_suffix(words: list[str]) -> str | None:
    """A parenthetical that belongs to the taxon line, as text with ``<i>``.

    NCBI's status suffixes are upright; a host name ("Buchnera aphidicola
    (Aphis fabae)") is italic. Anything else ("(ATCC 35973)") is strain text.
    """
    text = " ".join(words)
    if _STATUS_SUFFIX.fullmatch(text):
        return text
    inner = text[1:-1].split() if text.startswith("(") and text.endswith(")") else []
    if len(inner) >= 2 and _GENUS.fullmatch(inner[0]):
        if _is_epithet(inner[1]):
            return f"(<i>{inner[0]} {inner[1]}</i>{''.join(' ' + word for word in inner[2:])})"
        if inner[1] in {"sp.", "cf."}:
            return f"(<i>{inner[0]}</i> {' '.join(inner[1:])})"
    return None


def _taxon_words(tokens: list[str], strain_at: int) -> tuple[list[str], int] | None:
    """The taxon line's parts and the number of words they take, or None to keep the name whole.

    The taxon is a genus, then "sp." or a species epithet, then Latin
    infraspecific ranks (subsp., var., f., f. sp.) or one unranked zoological
    epithet, then NCBI's status or host suffixes.
    """
    bracketed = _BRACKETED_GENUS.fullmatch(tokens[0])
    genus = bracketed.group(1) if bracketed else tokens[0]
    if not (bracketed or _GENUS.fullmatch(genus)) or _ABOVE_GENUS.fullmatch(genus):
        return None
    parts = [f"[<i>{genus}</i>]" if bracketed else f"<i>{genus}</i>"]
    if len(tokens) == 1:
        return parts, 1
    if tokens[1] == "sp.":
        end = _parenthetical_end(tokens, 2)
        homonym = end > 2 and tokens[2] == "(in:"
        return parts + ["sp."] + (tokens[2:end] if homonym else []), end if homonym else 2
    if strain_at <= 1 or not _is_epithet(tokens[1]):
        return None
    parts.append(f"<i>{tokens[1]}</i>")
    index = 2
    while index < strain_at and tokens[index] in _LATIN_RANK_WORDS:
        rank = "f. sp." if tokens[index] == "f." and index + 1 < len(tokens) and tokens[index + 1] == "sp." else tokens[index]
        epithet_at = index + len(rank.split())
        if epithet_at >= strain_at or not _is_epithet(tokens[epithet_at]):
            return None
        parts.append(f"{rank} <i>{tokens[epithet_at]}</i>")
        index = epithet_at + 1
    if (
        index == 2
        and index < strain_at
        and _is_epithet(tokens[index])
        and (index + 1 in (len(tokens), strain_at) or not _LOWERCASE_WORD.fullmatch(tokens[index + 1]))
    ):
        parts.append(f"<i>{tokens[index]}</i>")
        index += 1
    while (end := _parenthetical_end(tokens, index)) > index and (suffix := _name_suffix(tokens[index:end])):
        parts.append(suffix)
        index = end
    rest = tokens[index:]
    if rest and index < strain_at and rest[0][:1].islower() and rest[0] not in _STRAIN_LINE_WORDS:
        return None
    if any(_ORGANISM_NOUN.search(word.lower().strip("()'\",;").split("-")[-1]) for word in rest):
        return None
    return parts, index


def split_organism_lines(metadata: RecordSourceMetadata) -> OrganismLines:
    """Split the source /organism and strain into a taxon line and a strain line.

    The taxon is italic after the codes: genus and epithets italic, rank
    words upright, an italic "Ca." before an upright Candidatus name
    (OV-1720), a virus name upright and an ICTV binomial species italic. A
    name that the rule cannot read as a taxon stays whole and upright on the
    taxon line. The strain line never repeats what the name already holds.
    """
    organism = " ".join(str(metadata.organism or "").split())
    strain = " ".join(str(metadata.strain or "").split())
    prefix = _CANDIDATUS_PREFIX.match(organism)
    name = organism[prefix.end():] if prefix else organism
    tokens = name.split()
    whole = organism.split()
    # Bakta writes "." when it was given no organism.
    if not re.search(r"[A-Za-z]", name) or name.lower() in _NON_ORGANISM_NAMES:
        return OrganismLines("", strain)
    keep_whole = OrganismLines(organism, _strain_line([], strain, whole))
    # The lines are markup ("<i>"); a name with markup characters keeps today's plain text.
    if re.search(r"[&<>]", organism + strain):
        return keep_whole
    if _is_virus(tokens, metadata.lineage_top):
        # An ICTV binomial ("Lambdavirus lambda") is italic; a virus name with a
        # number or letter ("Enterovirus J", "Cypovirus 1") is no species name.
        if prefix is None and len(tokens) == 2 and _ICTV_GENUS.fullmatch(tokens[0]) and _EPITHET.fullmatch(tokens[1]):
            return OrganismLines(f"<i>{organism}</i>", keep_whole.strain)
        return keep_whole
    taxon = _taxon_of(name, tokens, strain)
    if taxon is None:
        return keep_whole
    parts, used = taxon
    line = " ".join(parts).replace("</i> <i>", " ")
    if prefix:
        line = "<i>Ca.</i> " + re.sub(r"</?i>", "", line)
    return OrganismLines(line, _strain_line(tokens[used:], strain, whole))


def _is_virus(tokens: list[str], lineage: str | None) -> bool:
    return lineage == _VIRUS_LINEAGE or (
        lineage not in _CELLULAR_LINEAGES
        and any(_VIRUS_WORD.search(token) or _VIRUS_ACRONYM.fullmatch(token) for token in tokens)
    )


def _taxon_of(name: str, tokens: list[str], strain: str) -> tuple[list[str], int] | None:
    """The taxon of a name that is no virus (``_taxon_words``), or None to keep it whole."""
    if (
        _NON_ORGANISM_START.match(name)
        or _HYBRID.search(name)
        or any(token.lower().strip("()'\",;") in _SEQUENCE_WORDS for token in tokens)
    ):
        return None
    strain_at = _names_strain(tokens, strain)
    # A strain that spells the genus, an epithet ("pneumoniae") or the whole
    # name marks no boundary; the strain line still leaves it out.
    if strain_at is None or strain_at < 2 or _LOWERCASE_STRAIN.fullmatch(strain):
        strain_at = len(tokens)
    return _taxon_words(tokens, strain_at)


def _markup(parts: list[tuple[str, bool]]) -> str:
    return "".join(f"<i>{escape(text)}</i>" if italic else escape(text) for text, italic in parts if text)


def organism_name_lines(metadata: RecordSourceMetadata, typed_taxon: str = "", typed_strain: str = "") -> OrganismLines:
    """The two name lines a Circular definition draws, as text with ``<i>``.

    Without typed text they are the Auto lines (``split_organism_lines``). A
    typed Strain keeps both lines as typed, with the Auto taxon line when no
    line 1 is typed. A typed line 1 without a typed Strain is read with the
    Auto rule, the record's strain and lineage: when it reads a taxon and
    words remain, the typed text is cut there into two lines, its words and
    italics kept as typed. Otherwise the typed text is line 1 and the Auto
    strain line is line 2, unless the typed text already names the record's
    strain.
    """
    auto = split_organism_lines(metadata)
    if typed_strain:
        return OrganismLines(typed_taxon or auto.taxon, typed_strain)
    if not typed_taxon:
        return auto
    spans = [(str(part["text"] or ""), bool(part["italic"])) for part in parse_mixed_content_text(typed_taxon)]
    plain = "".join(text for text, _italic in spans)
    strain = " ".join(str(metadata.strain or "").split())
    words = plain.split()
    prefix = _CANDIDATUS_PREFIX.match(" ".join(words))
    prefix_words = len(prefix.group(0).split()) if prefix else 0
    name_tokens = words[prefix_words:]
    name = " ".join(name_tokens)
    readable = (
        re.search(r"[A-Za-z]", name)
        and name.lower() not in _NON_ORGANISM_NAMES
        and not _is_virus(name_tokens, metadata.lineage_top)
        and not _is_virus(str(metadata.organism or "").split(), metadata.lineage_top)
    )
    taxon = _taxon_of(name, name_tokens, strain) if readable else None
    if taxon is not None and prefix_words + taxon[1] < len(words):
        # Cut the spans between the last taxon word and the next word.
        cut = prefix_words + taxon[1]
        bounds = [match.span() for match in re.finditer(r"\S+", plain)]
        line_end, rest_start = bounds[cut - 1][1], bounds[cut][0]
        line: list[tuple[str, bool]] = []
        rest: list[tuple[str, bool]] = []
        offset = 0
        for text, italic in spans:
            span_start, offset = offset, offset + len(text)
            line.append((text[: max(0, line_end - span_start)], italic))
            rest.append((text[max(0, rest_start - span_start):], italic))
        return OrganismLines(_markup(line), _markup(rest))
    named = _names_strain(words, strain) is not None
    return OrganismLines(typed_taxon, "" if named else auto.strain)


def organism_name_rule_data() -> dict[str, object]:
    """The word lists and patterns of ``split_organism_lines``, as data for the Web.

    tools/generate_web_contracts.py writes them into the browser, where
    ``formatInferredOrganismStrain`` (services/record-discovery.js) runs the
    same steps without the Worker; tools/twins.json registers that twin. A
    pattern is its source and its flags ("i" for IGNORECASE).
    """

    def pattern(compiled: re.Pattern[str]) -> dict[str, str]:
        return {"source": compiled.pattern, "flags": "i" if compiled.flags & re.IGNORECASE else ""}

    return {
        "nonOrganismNames": sorted(_NON_ORGANISM_NAMES),
        "candidatusPrefix": pattern(_CANDIDATUS_PREFIX),
        "nonOrganismStart": pattern(_NON_ORGANISM_START),
        "sequenceWords": sorted(_SEQUENCE_WORDS),
        "virusWord": pattern(_VIRUS_WORD),
        "virusAcronym": pattern(_VIRUS_ACRONYM),
        "ictvGenus": pattern(_ICTV_GENUS),
        "cellularLineages": sorted(_CELLULAR_LINEAGES),
        "virusLineage": _VIRUS_LINEAGE,
        "hybrid": pattern(_HYBRID),
        "genus": pattern(_GENUS),
        "bracketedGenus": pattern(_BRACKETED_GENUS),
        "aboveGenus": pattern(_ABOVE_GENUS),
        "epithet": pattern(_EPITHET),
        "latinRankWords": sorted(_LATIN_RANK_WORDS),
        "strainLineWords": sorted(_STRAIN_LINE_WORDS),
        "notEpithets": sorted(_NOT_EPITHETS),
        "organismNoun": pattern(_ORGANISM_NOUN),
        "statusSuffix": pattern(_STATUS_SUFFIX),
        "designationPart": pattern(_DESIGNATION_PART),
        "lowercaseWord": pattern(_LOWERCASE_WORD),
        "lowercaseStrain": pattern(_LOWERCASE_STRAIN),
    }


def format_inferred_definition(metadata: RecordSourceMetadata) -> str:
    """The record's name on one line: the taxon line, then the strain line."""
    lines = split_organism_lines(metadata)
    return " ".join(line for line in lines if line)


def format_replicon_label(metadata: RecordSourceMetadata) -> str:
    """Name a replicon, falling back to the existing capitalized organelle label."""
    return str(metadata.replicon or "").strip() or str(metadata.organelle or "").strip().capitalize()


__all__ = [
    "OrganismLines",
    "RecordSourceMetadata",
    "format_inferred_definition",
    "format_replicon_label",
    "infer_record_source_metadata",
    "organism_name_lines",
    "split_organism_lines",
]
