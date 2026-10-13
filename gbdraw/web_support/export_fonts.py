"""The Web export's fonts at the Worker boundary (D-05).

``gbdraw.render.font_embedding`` owns both rules: the ``@font-face`` CSS an
exported SVG, Interactive SVG or PNG carries, and the TrueType faces a PDF
draws its text with. These functions only translate them to JSON for the
browser, which inserts the CSS and registers the faces with jsPDF.
"""

from __future__ import annotations

import base64
import json

from gbdraw.exceptions import GbdrawError, ValidationError
from gbdraw.render.font_embedding import PdfFontRun, pdf_fonts, svg_font_embedding


def svg_font_embedding_json(svg_text: str) -> str:
    """``{css, refusedFamilies}`` for the custom fonts ``svg_text`` draws with.

    A font that cannot be embedded fails the export with the
    ``FONT_EMBEDDING`` reason, whose guidance names the opt-out (D-09).
    """
    try:
        embedding = svg_font_embedding(str(svg_text))
    except GbdrawError:
        raise
    except Exception as exc:
        raise ValidationError(
            f"Could not embed the fonts in the SVG: {exc}",
            diagnostic={"code": "INPUT_INVALID", "reason": "FONT_EMBEDDING"},
        ) from exc
    return json.dumps({"css": embedding.css, "refusedFamilies": list(embedding.refused_families)})


def pdf_fonts_json(runs_json: str, known_json: str = "[]") -> str:
    """``{faces, fonts, substitutions}`` for the text runs of a PDF.

    ``runs_json``: ``[{fontFamily, bold, italic}]``; ``known_json``: the font
    files the browser already holds. ``fonts`` maps each other file to its
    base64 TrueType bytes.
    """
    runs: list[PdfFontRun] = [
        (str(run["fontFamily"]), bool(run["bold"]), bool(run["italic"])) for run in json.loads(str(runs_json))
    ]
    known = frozenset(str(name) for name in json.loads(str(known_json)))
    chosen = pdf_fonts(runs, known=known)
    return json.dumps({
        "faces": [{"family": face.family, "file": face.file} for face in chosen.faces],
        "fonts": {name: base64.b64encode(data).decode("ascii") for name, data in chosen.fonts.items()},
        "substitutions": [
            {"family": item.family, "substitute": item.substitute, "reason": item.reason}
            for item in chosen.substitutions
        ],
    })


__all__ = ["pdf_fonts_json", "svg_font_embedding_json"]
