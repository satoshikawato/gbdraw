from __future__ import annotations

import re
from collections.abc import Iterable, Iterator
from xml.sax.saxutils import escape
from typing import TYPE_CHECKING, Literal, NamedTuple, TypedDict, TypeVar

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
    # The record's lineage (GenBank ORGANISM lines), highest taxon first.
    lineage: tuple[str, ...] = ()


def infer_record_source_metadata(record: SeqRecord) -> RecordSourceMetadata:
    """Extract source-feature metadata used in definition labels."""
    annotations = getattr(record, "annotations", None) or {}
    organism = str(annotations.get("organism", "") or "").strip()
    taxonomy = annotations.get("taxonomy") or []
    lineage = tuple(word for word in (str(taxon).strip() for taxon in taxonomy) if word)
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
        lineage=lineage,
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
# What a Candidatus name is drawn after: an italic "Ca." (OV-1720); the name is upright.
_CANDIDATUS_MARK = "<i>Ca.</i>"
# A name needs a letter: Bakta writes "." when it was given no organism.
_HAS_LETTER = re.compile(r"[A-Za-z]")
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
# The nomenclature code of a name, from its record's lineage; "any" when the
# lineage names none (other eukaryotes, or no lineage, as for GFF3 input).
_NomenclatureCode = Literal["zoological", "botanical", "prokaryotic", "any"]
_ANY_CODE: _NomenclatureCode = "any"
# Lineage words that name the code.
_CODE_LINEAGES: dict[str, _NomenclatureCode] = {
    "Bacteria": "prokaryotic", "Archaea": "prokaryotic", "Metazoa": "zoological", "Viridiplantae": "botanical",
    "Fungi": "botanical", "Rhodophyta": "botanical", "Ochrophyta": "botanical",
}
# Endings of names above genus under each code: such a first word is no genus.
# "any" takes every ending, so a name of unknown code stays whole and upright.
# The botanical list holds the ICN order ending (-ales) too, so "Helotiales sp.
# EF0021" stays whole (OV-1834). The division ending -phyta is left out: plant
# and insect genera end in it (Xerophyta, Leucophyta, Symphyta; OV-1835).
_ABOVE_GENUS_ENDINGS: dict[_NomenclatureCode, tuple[str, ...]] = {
    "zoological": ("idae", "inae"),
    "botanical": ("aceae", "ineae", "inae", "idae", "opsida", "phyceae", "mycetes", "mycota", "ales"),
    "prokaryotic": ("aceae", "ales", "ineae", "idae"),
}
_ABOVE_GENUS_ENDINGS[_ANY_CODE] = tuple(dict.fromkeys(e for endings in _ABOVE_GENUS_ENDINGS.values() for e in endings))
_ABOVE_GENUS = {
    code: re.compile(rf"[A-Z][a-z]*(?:{'|'.join(endings)})") for code, endings in _ABOVE_GENUS_ENDINGS.items()
}
# The codes that set names above genus in italics (the ICNP and the ICN; the
# ICZN does not), with the endings that make a placeholder's first word such a
# name ("Lachnospiraceae bacterium A4"); -ota is the ICNP phylum ending.
_ITALIC_ABOVE_GENUS_ENDINGS: dict[_NomenclatureCode, tuple[str, ...]] = {
    "prokaryotic": (*_ABOVE_GENUS_ENDINGS["prokaryotic"], "ota"),
    "botanical": _ABOVE_GENUS_ENDINGS["botanical"],
}
_ITALIC_ABOVE_GENUS = {
    code: re.compile(rf"[A-Z][a-z]*(?:{'|'.join(endings)})") for code, endings in _ITALIC_ABOVE_GENUS_ENDINGS.items()
}
# The word before a placeholder's higher-taxon name ("uncultured Lachnospiraceae bacterium").
_PLACEHOLDER_PREFIX = "uncultured"
# Hybrid signs, kept as the record writes them and upright (ICN Rec. H.3A.2
# allows the letter x); the multiplication sign may also stand joined to the
# next word ("Mentha ×piperita").
_HYBRID_SIGNS = frozenset({"x", "×"})
_JOINED_HYBRID_SIGN = "×"
# A word after a sign that starts a parent of a formula, or a nothogenus: a
# capital or a bracketed genus.
_PARENT_START = re.compile(r"\[?[A-Z]")
# A hyphenated genus (ICN Art. 20.3: Pseudo-nitzschia, Sciuro-hypnum), read only
# under the botanical code: "Legionella-like" and "Gram-positive" start
# placeholders elsewhere. Neither part may be a word of _NOT_EPITHETS.
_HYPHEN_GENUS = re.compile(r"([A-Z][a-z]+)-([a-z]{3,})")
_HYPHEN_GENUS_CODE: _NomenclatureCode = "botanical"
# A letter-hyphen epithet ("w-albus", "c-aureum") is an epithet too.
_EPITHET = re.compile(r"(?:[a-z]-)?[a-z]{2,}(?:-[a-z]+)*")
# The word after a genus that names no species ("Bacillus sp."), and the words
# after a host genus that name none ("(Aphis cf. fabae)").
_UNNAMED_SPECIES = "sp."
_UNNAMED_HOST_WORDS = frozenset({_UNNAMED_SPECIES, "cf."})
# The word that starts NCBI's homonym suffix after "sp." ("Bacillus sp. (in: Bacteria)").
_HOMONYM_START = "(in:"
# The two-word rank of a forma specialis ("Fusarium oxysporum f. sp. lycopersici").
_FORMA_SPECIALIS = "f. sp."
# Characters stripped from a word before it is looked up in a word list.
_WORD_STRIP = "()'\",;"
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
     "species", "subspecies", "unnamed", "fungal", "fungus", "contaminant", "mycorrhizal", "non"}
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
# The separator of the strain's parts ("K-12; MG1655").
_STRAIN_SEPARATOR = re.compile(r"\s*;\s*")


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
    parts = [part for part in _STRAIN_SEPARATOR.split(strain) if part]
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
        if inner[1] in _UNNAMED_HOST_WORDS:
            return f"(<i>{inner[0]}</i> {' '.join(inner[1:])})"
    return None


def _nomenclature_code(lineage: tuple[str, ...]) -> _NomenclatureCode:
    return next((_CODE_LINEAGES[word] for word in lineage if word in _CODE_LINEAGES), _ANY_CODE)


def _taxon_words(tokens: list[str], strain_at: int, code: _NomenclatureCode) -> tuple[list[str], int] | None:
    """The taxon line's parts and the number of words they take, or None to keep the name whole.

    The taxon is a genus, then "sp." or a species epithet, then Latin
    infraspecific ranks (subsp., var., f., f. sp.) or one unranked zoological
    epithet, then NCBI's status or host suffixes.
    """
    bracketed = _BRACKETED_GENUS.fullmatch(tokens[0])
    genus = bracketed.group(1) if bracketed else tokens[0]
    hyphen = _HYPHEN_GENUS.fullmatch(genus) if code == _HYPHEN_GENUS_CODE else None
    if (
        not (bracketed or _GENUS.fullmatch(genus) or (hyphen and not {hyphen.group(1).lower(), hyphen.group(2)} & _NOT_EPITHETS))
        or _ABOVE_GENUS[code].fullmatch(genus)
    ):
        return None
    parts = [f"[<i>{genus}</i>]" if bracketed else f"<i>{genus}</i>"]
    if len(tokens) == 1:
        return parts, 1
    if tokens[1] == _UNNAMED_SPECIES:
        end = _parenthetical_end(tokens, 2)
        homonym = end > 2 and tokens[2] == _HOMONYM_START
        return parts + [_UNNAMED_SPECIES] + (tokens[2:end] if homonym else []), end if homonym else 2
    if strain_at <= 1 or not _is_epithet(tokens[1]):
        return None
    parts.append(f"<i>{tokens[1]}</i>")
    index = 2
    while index < strain_at and tokens[index] in _LATIN_RANK_WORDS:
        rank = _FORMA_SPECIALIS if " ".join(tokens[index:index + 2]) == _FORMA_SPECIALIS else tokens[index]
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
    if any(_ORGANISM_NOUN.search(word.lower().strip(_WORD_STRIP).split("-")[-1]) for word in rest):
        return None
    return parts, index


def split_organism_lines(metadata: RecordSourceMetadata) -> OrganismLines:
    """Split the source /organism and strain into a taxon line and a strain line.

    The taxon is italic after the codes: genus and epithets italic, rank
    words upright, an italic "Ca." before an upright Candidatus name
    (OV-1720), a virus name upright and an ICTV binomial species italic. A
    name that the rule cannot read as a taxon stays whole on the taxon line,
    upright except for a hybrid's parents and a placeholder's higher-taxon
    name (``_whole_name_line``); a Candidatus name kept whole keeps its italic
    "Ca.". Other text with markup characters stays as written. The strain line
    never repeats what the name already holds.
    """
    strain = _single_spaced(metadata.strain)
    auto = _auto_taxon_line(metadata)
    if auto is None:
        return OrganismLines("", strain)
    strain_line = _strain_line(auto.rest, strain, _single_spaced(metadata.organism).split())
    return OrganismLines(auto.line, escape(strain_line) if auto.escaped else strain_line)


class _AutoTaxonLine(NamedTuple):
    """The Auto rule's line 1, and the /organism words it leaves to line 2."""

    line: str
    rest: list[str]
    # The line escapes markup characters, so line 2 must escape them too.
    escaped: bool = False


def _single_spaced(text: str | None) -> str:
    return " ".join(str(text or "").split())


def _auto_taxon_line(metadata: RecordSourceMetadata) -> _AutoTaxonLine | None:
    """Line 1 of ``split_organism_lines``, or None when the record names no organism."""
    organism = _single_spaced(metadata.organism)
    strain = _single_spaced(metadata.strain)
    prefix = _CANDIDATUS_PREFIX.match(organism)
    name = organism[prefix.end():] if prefix else organism
    tokens = name.split()
    if not _HAS_LETTER.search(name) or name.lower() in _NON_ORGANISM_NAMES:
        return None
    # The lines are markup ("<i>"), so a name with markup characters is kept
    # whole: a Candidatus name keeps its italic "Ca." and escapes the rest
    # (PD-OI-107 keeps the Ca. form); any other name stays plain text as written.
    if re.search(r"[&<>]", organism + strain):
        return _AutoTaxonLine(f"{_CANDIDATUS_MARK} {escape(name)}", [], True) if prefix else _AutoTaxonLine(organism, [])
    # A Candidatus name kept whole: an italic "Ca.", then the upright name.
    candidatus = _AutoTaxonLine(f"{_CANDIDATUS_MARK} {name}", []) if prefix else None
    if _is_virus(tokens, metadata.lineage):
        # An ICTV binomial ("Lambdavirus lambda") is italic; a virus name with a
        # number or letter ("Enterovirus J", "Cypovirus 1") is no species name.
        if prefix is None and len(tokens) == 2 and _ICTV_GENUS.fullmatch(tokens[0]) and _EPITHET.fullmatch(tokens[1]):
            return _AutoTaxonLine(f"<i>{organism}</i>", [])
        return candidatus or _AutoTaxonLine(organism, [])
    taxon = _taxon_of(name, tokens, strain, metadata.lineage)
    if taxon is None:
        return candidatus or _AutoTaxonLine(_whole_name_line(name, tokens, metadata.lineage), [])
    parts, used = taxon
    line = " ".join(parts).replace("</i> <i>", " ")
    if prefix:
        line = f"{_CANDIDATUS_MARK} " + re.sub(r"</?i>", "", line)
    return _AutoTaxonLine(line, tokens[used:])


def _strain_line_under(typed_line: str, metadata: RecordSourceMetadata) -> str:
    """Line 2 under a typed line 1: the Auto strain line read against the typed text.

    The typed text replaces the Auto line 1, so a strain part that only the
    /organism named is drawn here. Nothing is drawn when the typed text names
    the record's strain (PD-OI-107); otherwise the words the Auto rule puts
    after the taxon, then each strain part that neither the typed text nor
    those words name.
    """
    words = "".join(str(part["text"] or "") for part in parse_mixed_content_text(typed_line)).split()
    strain = _single_spaced(metadata.strain)
    if _names_strain(words, strain) is not None:
        return ""
    auto = _auto_taxon_line(metadata)
    rest = auto.rest if auto else []
    return _strain_line(rest, strain, words + rest)


def _hybrid_line(tokens: list[str], code: _NomenclatureCode) -> str | None:
    """A hybrid name or formula as text with ``<i>``, or None when a part is no taxon.

    A sign before a capitalized word starts a parent (a formula, "x
    Triticosecale"); a sign before an epithet belongs to the nothotaxon
    ("Citrus x limon", "subsp. x tunetana"). Each parent is read like a name
    without a strain (``_taxon_words``); its other words stay upright, and so
    does every sign, which must stand before a word of a parent's name
    ("Phytophthora sp. x amnicola" stays upright).
    """
    parents: list[list[tuple[str, str]]] = []  # each parent: (sign before the word, word)
    sign = ""
    for token in tokens:
        if token in _HYBRID_SIGNS:
            if sign:
                return None
            sign = token + " "
            continue
        word = token
        if token.startswith(_JOINED_HYBRID_SIGN) and len(token) > 1:
            if sign:
                return None
            sign, word = _JOINED_HYBRID_SIGN, token[1:]
        if sign and _PARENT_START.match(word):
            parents.append([])
        elif sign and not (parents and _is_epithet(word)):
            return None
        elif not parents:
            parents.append([])
        parents[-1].append((sign, word))
        sign = ""
    if sign or not parents:
        return None
    words: list[str] = []
    for parent in parents:
        plain = [word for _sign, word in parent]
        reading = _taxon_words(plain, len(plain), code)
        if reading is None:
            return None
        parts, used = reading
        marked = " ".join(parts).split(" ")
        # One marked word per word taken (the zip below relies on it), and no
        # sign before the parent's strain or other words.
        if len(marked) != used or any(sign for sign, _word in parent[used:]):
            return None
        words += [sign + text for (sign, _word), text in zip(parent, marked + plain[used:])]
    return " ".join(words).replace("</i> <i>", " ")


def _whole_name_line(name: str, tokens: list[str], lineage: tuple[str, ...]) -> str:
    """A name the rule keeps whole, as text with ``<i>``.

    A hybrid is read part by part (``_hybrid_line``). A placeholder's first
    word, or the word after "uncultured", is italic when the lineage's code
    sets names above genus in italics, the word has that code's ending, and
    the record's lineage holds the word as written (no misspelling, no
    Candidatus name written without its prefix). Any other name stays upright.
    """
    code = _nomenclature_code(lineage)
    if _HYBRID.search(name):
        return _hybrid_line(tokens, code) or name
    higher = _ITALIC_ABOVE_GENUS.get(code)
    at = 1 if len(tokens) > 1 and tokens[0] == _PLACEHOLDER_PREFIX else 0
    if higher is None or not higher.fullmatch(tokens[at]) or tokens[at] not in lineage:
        return name
    return " ".join([*tokens[:at], f"<i>{tokens[at]}</i>", *tokens[at + 1:]])


def _is_virus(tokens: list[str], lineage: tuple[str, ...]) -> bool:
    top = lineage[0] if lineage else None
    return top == _VIRUS_LINEAGE or (
        top not in _CELLULAR_LINEAGES
        and any(_VIRUS_WORD.search(token) or _VIRUS_ACRONYM.fullmatch(token) for token in tokens)
    )


def _taxon_of(name: str, tokens: list[str], strain: str, lineage: tuple[str, ...]) -> tuple[list[str], int] | None:
    """The taxon of a name that is no virus (``_taxon_words``), or None to keep it whole."""
    if (
        _NON_ORGANISM_START.match(name)
        or _HYBRID.search(name)
        or any(token.lower().strip(_WORD_STRIP) in _SEQUENCE_WORDS for token in tokens)
    ):
        return None
    strain_at = _names_strain(tokens, strain)
    # A strain that spells the genus, an epithet ("pneumoniae") or the whole
    # name marks no boundary; the strain line still leaves it out.
    if strain_at is None or strain_at < 2 or _LOWERCASE_STRAIN.fullmatch(strain):
        strain_at = len(tokens)
    return _taxon_words(tokens, strain_at, _nomenclature_code(lineage))


def _markup(parts: list[tuple[str, bool]]) -> str:
    return "".join(f"<i>{escape(text)}</i>" if italic else escape(text) for text, italic in parts if text)


def organism_name_lines(
    metadata: RecordSourceMetadata, typed_taxon: str = "", typed_strain: str = "", *, record_label: str = ""
) -> OrganismLines:
    """The two name lines a Circular definition draws, as text with ``<i>``.

    Without typed text they are the Auto lines (``split_organism_lines``). A
    typed Record label is line 1 as typed and wins over a typed Species;
    line 2 is then the typed strain or the Auto strain line (PD-OI-012:
    explicit values override their lines, empty values keep the inferred
    output). A typed Strain keeps both lines as typed, with the Auto taxon
    line when no line 1 is typed. A typed Species (``typed_taxon``) without a
    typed Strain is read with the Auto rule, the record's strain and lineage:
    when it reads a taxon and words remain, the typed text is cut there into
    two lines, its words and italics kept as typed (DF-D09). Otherwise the
    typed text is line 1. Under a typed line 1 and no typed line 2, line 2 is
    the Auto strain line read against the typed text (``_strain_line_under``).
    """
    if record_label:
        return OrganismLines(record_label, typed_strain or _strain_line_under(record_label, metadata))
    if typed_strain:
        return OrganismLines(typed_taxon or split_organism_lines(metadata).taxon, typed_strain)
    if not typed_taxon:
        return split_organism_lines(metadata)
    spans = [(str(part["text"] or ""), bool(part["italic"])) for part in parse_mixed_content_text(typed_taxon)]
    plain = "".join(text for text, _italic in spans)
    strain = _single_spaced(metadata.strain)
    words = plain.split()
    prefix = _CANDIDATUS_PREFIX.match(" ".join(words))
    prefix_words = len(prefix.group(0).split()) if prefix else 0
    name_tokens = words[prefix_words:]
    name = " ".join(name_tokens)
    readable = (
        _HAS_LETTER.search(name)
        and name.lower() not in _NON_ORGANISM_NAMES
        and not _is_virus(name_tokens, metadata.lineage)
        and not _is_virus(str(metadata.organism or "").split(), metadata.lineage)
    )
    taxon = _taxon_of(name, name_tokens, strain, metadata.lineage) if readable else None
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
    return OrganismLines(typed_taxon, _strain_line_under(typed_taxon, metadata))


class NamePatternData(TypedDict):
    """A pattern of the name rule as the Web reads it: its source, and "i" for IGNORECASE."""

    source: str
    flags: Literal["i", ""]


class OrganismNameRuleData(TypedDict):
    """The word lists, words and patterns of ``split_organism_lines``; every field is required."""

    nonOrganismNames: list[str]
    candidatusPrefix: NamePatternData
    candidatusMark: str
    hasLetter: NamePatternData
    nonOrganismStart: NamePatternData
    sequenceWords: list[str]
    virusWord: NamePatternData
    virusAcronym: NamePatternData
    ictvGenus: NamePatternData
    cellularLineages: list[str]
    virusLineage: str
    hybrid: NamePatternData
    genus: NamePatternData
    bracketedGenus: NamePatternData
    hyphenGenus: NamePatternData
    hyphenGenusCode: _NomenclatureCode
    codeLineages: dict[str, _NomenclatureCode]
    anyCode: _NomenclatureCode
    aboveGenus: dict[_NomenclatureCode, NamePatternData]
    italicAboveGenus: dict[_NomenclatureCode, NamePatternData]
    placeholderPrefix: str
    hybridSigns: list[str]
    joinedHybridSign: str
    parentStart: NamePatternData
    epithet: NamePatternData
    unnamedSpecies: str
    unnamedHostWords: list[str]
    homonymStart: str
    formaSpecialis: str
    wordStrip: str
    latinRankWords: list[str]
    strainLineWords: list[str]
    notEpithets: list[str]
    organismNoun: NamePatternData
    statusSuffix: NamePatternData
    designationPart: NamePatternData
    lowercaseWord: NamePatternData
    lowercaseStrain: NamePatternData
    strainSeparator: NamePatternData


def organism_name_rule_data() -> OrganismNameRuleData:
    """The word lists, words and patterns of ``split_organism_lines``, as data for the Web.

    tools/generate_web_contracts.py writes them into the browser, where
    ``formatInferredOrganismStrain`` (services/record-discovery.js) runs the
    same steps without the Worker; tools/twins.json registers that twin.
    """

    def pattern(compiled: re.Pattern[str]) -> NamePatternData:
        return {"source": compiled.pattern, "flags": "i" if compiled.flags & re.IGNORECASE else ""}

    return {
        "nonOrganismNames": sorted(_NON_ORGANISM_NAMES),
        "candidatusPrefix": pattern(_CANDIDATUS_PREFIX),
        "candidatusMark": _CANDIDATUS_MARK,
        "hasLetter": pattern(_HAS_LETTER),
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
        "hyphenGenus": pattern(_HYPHEN_GENUS),
        "hyphenGenusCode": _HYPHEN_GENUS_CODE,
        "codeLineages": dict(_CODE_LINEAGES),
        "anyCode": _ANY_CODE,
        "aboveGenus": {code: pattern(compiled) for code, compiled in _ABOVE_GENUS.items()},
        "italicAboveGenus": {code: pattern(compiled) for code, compiled in _ITALIC_ABOVE_GENUS.items()},
        "placeholderPrefix": _PLACEHOLDER_PREFIX,
        "hybridSigns": sorted(_HYBRID_SIGNS),
        "joinedHybridSign": _JOINED_HYBRID_SIGN,
        "parentStart": pattern(_PARENT_START),
        "epithet": pattern(_EPITHET),
        "unnamedSpecies": _UNNAMED_SPECIES,
        "unnamedHostWords": sorted(_UNNAMED_HOST_WORDS),
        "homonymStart": _HOMONYM_START,
        "formaSpecialis": _FORMA_SPECIALIS,
        "wordStrip": _WORD_STRIP,
        "latinRankWords": sorted(_LATIN_RANK_WORDS),
        "strainLineWords": sorted(_STRAIN_LINE_WORDS),
        "notEpithets": sorted(_NOT_EPITHETS),
        "organismNoun": pattern(_ORGANISM_NOUN),
        "statusSuffix": pattern(_STATUS_SUFFIX),
        "designationPart": pattern(_DESIGNATION_PART),
        "lowercaseWord": pattern(_LOWERCASE_WORD),
        "lowercaseStrain": pattern(_LOWERCASE_STRAIN),
        "strainSeparator": pattern(_STRAIN_SEPARATOR),
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
