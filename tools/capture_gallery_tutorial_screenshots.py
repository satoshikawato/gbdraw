#!/usr/bin/env python3
from __future__ import annotations

import argparse
import contextlib
import http.server
import io
import json
import math
import re
import socketserver
import sys
import threading
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, TypedDict

try:
    from PIL import Image, ImageChops
except ImportError as exc:  # pragma: no cover - dependency guard
    raise SystemExit(
        "Pillow is required. Install the dev dependencies with `pip install -e .[dev]`."
    ) from exc

if TYPE_CHECKING:
    from playwright.sync_api import ElementHandle


REPO_ROOT = Path(__file__).resolve().parents[1]
WEB_ROOT = REPO_ROOT / "gbdraw" / "web"
GALLERY_ROOT = WEB_ROOT / "gallery"
GALLERY_TUTORIAL_ROOT = GALLERY_ROOT / "tutorials"
GALLERY_MEDIA_ROOT = GALLERY_ROOT / "media"

READY_SECTIONS = (
    "manualSteps",
    "postGenerationEdits",
    "colorRules",
)
OPERATION_REQUIRED_SECTIONS = (
    "manualSteps",
    "postGenerationEdits",
)

GENERIC_OPERATION_BODY_RE = re.compile(r"^\s*Confirm\b", re.IGNORECASE)
GENERIC_OPERATION_MEDIA_RE = re.compile(
    r"operation screenshot showing|"
    r"in the web app session view\.?$|"
    r"in the gallery page view\.?$|"
    r"in the generated preview view\.?$|"
    r"in the popup view\.?$|"
    r"operation screenshot",
    re.IGNORECASE,
)
GALLERY_MEDIA_REFERENCE_RE = re.compile(r"^\./media/([^/]+)/")

DEFAULT_MAX_IMAGE_WIDTH = 4200
DEFAULT_MAX_IMAGE_HEIGHT = 4200
DEFAULT_MAX_FILE_SIZE_KB = 2400
DEFAULT_WEBP_QUALITY = 95
DEFAULT_DEVICE_SCALE_FACTOR = 3
DEFAULT_VIEWPORT = {"width": 1280, "height": 900}
DEFAULT_GENERATE_TIMEOUT_MS = 1_800_000
CAPTURE_ACTION_TYPES = frozenset(
    {
        "click",
        "fill",
        "select",
        "check",
        "uncheck",
        "waitForSelector",
        "waitForFunction",
        "evaluate",
        "wait",
        "resetSettings",
        "generate",
        "upload",
    }
)
CAPTURE_SELECT_MENU_SELECTOR = ".gbdraw-capture-select-menu"
CAPTURE_SELECT_CONTROL_SELECTOR = ".gbdraw-capture-select-control"
# A drawn openSelect menu shows at most this many rows, scrolled so the
# selected option is in view, like a native select list.
OPEN_SELECT_MAX_VISIBLE_OPTIONS = 17
SESSION_FILE_CHOOSER_SELECTOR = "[data-session-file-picker]"
SESSION_FILE_CHOOSER_NAME = "Files in this Session"
CAPTURE_POSITION_PICK = re.compile(r"\bchildren\s*(?:\?\.)?\s*\[")
CAPTURE_RECORD_READS = (
    "inspectCircularSourceRecords",
    "refreshLinearRecordSelectors",
    "Inspect source records",
    "Inspect records",
)
# A crop edge never cuts a line of text or a control a reader sees (OV-3410): an edge that
# crosses one moves to the nearer position that crosses none, outward at most CROP_EDGE_REACH
# CSS px from where the crop put it or inward within the crop's padding, never into what the
# crop shows; it stops CROP_EDGE_MARGIN px clear of what it passes. CROP_EDGE_TOLERANCE only
# absorbs rounding: a glyph or a border 1 CSS px over the edge is a visible cut.
CROP_EDGE_REACH = 160.0
CROP_EDGE_MARGIN = 4.0
CROP_EDGE_TOLERANCE = 0.1
CROP_EDGE_ROUNDS = 8
# The frame a capture draws around its highlightSelector never covers text or a control a reader
# sees, the highlighted element's own text and the value its field shows included (OV-3492): an outline HIGHLIGHT_FRAME_WIDTH px wide,
# drawn HIGHLIGHT_FRAME_OFFSET px outside the element, with a glow HIGHLIGHT_FRAME_GLOW px wide
# (HTML elements only), moves in (down to just inside the element's edge) where text outside sits
# closer, and stays HIGHLIGHT_FRAME_CLEARANCE px clear of it; text inside the element does not move
# it, and text across the element's edge that it cannot pass within its reach fails the capture, as
# does a drawn frame that paints a pixel of text (highlight_frame_style,
# assert_highlight_frames_cover_no_text). An open-select overlay frames its control and its
# selected option inside them, in the same colour and width (highlight_frame_css).
HIGHLIGHT_FRAME_COLOR = "#f59e0b"
HIGHLIGHT_FRAME_GLOW_COLOR = "rgba(245, 158, 11, 0.18)"
HIGHLIGHT_FRAME_WIDTH = 3.0
HIGHLIGHT_FRAME_OFFSET = 3.0
HIGHLIGHT_FRAME_GLOW = 4.0
HIGHLIGHT_FRAME_CLEARANCE = 1.0


class QuietGalleryRequestHandler(http.server.SimpleHTTPRequestHandler):
    def end_headers(self) -> None:
        self.send_header("Cross-Origin-Opener-Policy", "same-origin")
        self.send_header("Cross-Origin-Embedder-Policy", "require-corp")
        self.send_header("Cross-Origin-Resource-Policy", "same-origin")
        super().end_headers()

    def log_message(self, format: str, *args: object) -> None:
        return


@dataclass(frozen=True)
class TutorialContext:
    example_id: str
    tutorial_path: Path
    section: str
    step_index: int
    step_title: str
    operation_index: int | None = None
    operation_title: str = ""

    @property
    def label(self) -> str:
        parts = [
            self.example_id,
            self.section,
            f"step {self.step_index}",
        ]
        if self.step_title:
            parts.append(self.step_title)
        if self.operation_index is not None:
            parts.append(f"operation {self.operation_index}")
        if self.operation_title:
            parts.append(self.operation_title)
        return " / ".join(parts)


@dataclass
class ValidationResult:
    errors: list[str]
    warnings: list[str]
    media_count: int = 0
    operation_count: int = 0
    operation_media_count: int = 0


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def as_array(value: Any) -> list[Any]:
    return value if isinstance(value, list) else []


def as_text(value: Any) -> str:
    return value.strip() if isinstance(value, str) else ""


def tutorial_status(sample: dict[str, Any]) -> str:
    status = as_text(sample.get("tutorialStatus")).lower()
    if status in {"ready", "draft", "planned"}:
        return status
    return "ready" if sample.get("tutorial") else "planned"


def load_ready_examples(example_id: str | None = None) -> list[dict[str, Any]]:
    examples = read_json(GALLERY_ROOT / "examples.json")
    if not isinstance(examples, list):
        raise RuntimeError("gallery/examples.json must contain an array.")

    ready = [sample for sample in examples if isinstance(sample, dict) and tutorial_status(sample) == "ready"]
    if example_id:
        ready = [sample for sample in ready if sample.get("id") == example_id]
        if not ready:
            raise RuntimeError(f"No ready Gallery example found for {example_id}.")
    return ready


def resolve_gallery_reference(reference: str, base_dir: Path = GALLERY_ROOT) -> Path | None:
    if not reference or "://" in reference:
        return None
    if reference.startswith("/"):
        return (WEB_ROOT / reference.lstrip("/")).resolve()
    return (base_dir / reference).resolve()


def media_entries(media: Any) -> list[dict[str, Any]]:
    if media is None:
        return []
    raw_entries = media if isinstance(media, list) else [media]
    entries: list[dict[str, Any]] = []
    for entry in raw_entries:
        if isinstance(entry, str):
            entries.append({"src": entry})
        elif isinstance(entry, dict):
            entries.append(entry)
    return entries


def iter_step_contexts(sample: dict[str, Any], tutorial: dict[str, Any]):
    tutorial_path = resolve_gallery_reference(as_text(sample.get("tutorial")))
    if tutorial_path is None:
        raise RuntimeError(f"Tutorial path for {sample.get('id')} is not local.")

    for section in READY_SECTIONS:
        for step_index, step in enumerate(as_array(tutorial.get(section)), start=1):
            if isinstance(step, str):
                step_object: dict[str, Any] = {"body": step}
            elif isinstance(step, dict):
                step_object = step
            else:
                continue
            context = TutorialContext(
                example_id=as_text(sample.get("id")),
                tutorial_path=tutorial_path,
                section=section,
                step_index=step_index,
                step_title=as_text(step_object.get("title")),
            )
            yield context, step_object


def iter_operation_contexts(sample: dict[str, Any], tutorial: dict[str, Any]):
    for step_context, step in iter_step_contexts(sample, tutorial):
        for operation_index, operation in enumerate(as_array(step.get("operations")), start=1):
            if isinstance(operation, str):
                operation_object: dict[str, Any] = {"body": operation}
            elif isinstance(operation, dict):
                operation_object = operation
            else:
                continue
            yield (
                TutorialContext(
                    example_id=step_context.example_id,
                    tutorial_path=step_context.tutorial_path,
                    section=step_context.section,
                    step_index=step_context.step_index,
                    step_title=step_context.step_title,
                    operation_index=operation_index,
                    operation_title=as_text(operation_object.get("title")),
                ),
                operation_object,
            )


def media_is_video(entry: dict[str, Any], source: str) -> bool:
    media_type = as_text(entry.get("type")).lower()
    return media_type == "video" or source.lower().split("?", 1)[0].endswith((".mp4", ".webm", ".ogg"))


def add_media_validation(
    result: ValidationResult,
    context: TutorialContext,
    entry: dict[str, Any],
    *,
    operation_media: bool,
    max_width: int,
    max_height: int,
    max_file_size_kb: int,
) -> None:
    source = as_text(entry.get("src")) or as_text(entry.get("href"))
    if not source:
        result.errors.append(f"{context.label}: media entry is missing src.")
        return

    result.media_count += 1
    if operation_media:
        result.operation_media_count += 1

    is_video = media_is_video(entry, source)
    if not is_video and not as_text(entry.get("alt")):
        result.errors.append(f"{context.label}: {source} is missing alt text.")
    if operation_media and not as_text(entry.get("caption")):
        result.errors.append(f"{context.label}: {source} is missing an operation caption.")

    path = resolve_gallery_reference(source)
    if path is None:
        result.warnings.append(f"{context.label}: external media was not checked: {source}")
        return
    if not path.exists():
        result.errors.append(f"{context.label}: media file does not exist: {source}")
        return
    if path.is_dir():
        result.errors.append(f"{context.label}: media source resolves to a directory: {source}")
        return

    file_size_kb = path.stat().st_size / 1024
    if file_size_kb > max_file_size_kb:
        message = f"{context.label}: {source} is {file_size_kb:.0f} KiB, above {max_file_size_kb} KiB."
        if operation_media:
            result.errors.append(message)
        else:
            result.warnings.append(message)

    if is_video:
        return

    try:
        with Image.open(path) as image:
            width, height = image.size
    except Exception as exc:
        result.errors.append(f"{context.label}: could not read image dimensions for {source}: {exc}")
        return

    too_large = width > max_width or height > max_height
    full_screen_like = width >= 5400 or (width >= 4200 and height >= 2700) or (width >= 3600 and height >= 3000)
    if too_large or full_screen_like:
        message = (
            f"{context.label}: {source} dimensions are {width}x{height}; "
            f"limit is {max_width}x{max_height}."
        )
        if operation_media:
            result.errors.append(message)
        else:
            result.warnings.append(message)


def add_operation_text_validation(
    result: ValidationResult,
    context: TutorialContext,
    operation: dict[str, Any],
) -> None:
    body = as_text(operation.get("body"))
    if GENERIC_OPERATION_BODY_RE.search(body):
        result.errors.append(
            f"{context.label}: operation body uses generic Confirm wording."
        )

    for entry in media_entries(operation.get("media")):
        for field in ("alt", "caption"):
            text = as_text(entry.get(field))
            if GENERIC_OPERATION_MEDIA_RE.search(text):
                result.errors.append(
                    f"{context.label}: media {field} uses generic screenshot wording."
                )


def add_capture_contract_validation(
    result: ValidationResult,
    context: TutorialContext,
    sample: dict[str, Any],
    operation: dict[str, Any],
) -> None:
    data_dependent = operation.get("dataDependent") is True
    generic_media = operation.get("genericMedia") is True
    capture = operation.get("capture")

    if "dataDependent" in operation and not isinstance(operation.get("dataDependent"), bool):
        result.errors.append(f"{context.label}: dataDependent must be a boolean.")
    if "genericMedia" in operation and not isinstance(operation.get("genericMedia"), bool):
        result.errors.append(f"{context.label}: genericMedia must be a boolean.")
    if data_dependent and generic_media:
        result.errors.append(
            f"{context.label}: data-dependent media cannot be declared generic."
        )

    for entry in media_entries(operation.get("media")):
        source = as_text(entry.get("src")) or as_text(entry.get("href"))
        match = GALLERY_MEDIA_REFERENCE_RE.match(source)
        if not match or match.group(1) == context.example_id:
            continue
        owner = match.group(1)
        if data_dependent:
            result.errors.append(
                f"{context.label}: data-dependent media must be owned by "
                f"{context.example_id}, not {owner}: {source}"
            )
        elif not generic_media:
            result.errors.append(
                f"{context.label}: cross-example media requires genericMedia=true: {source}"
            )

    if not isinstance(capture, dict):
        if data_dependent:
            result.errors.append(
                f"{context.label}: data-dependent media requires capture metadata."
            )
        return

    assert_app_state = capture.get("assertAppState")
    visible_controls = capture.get("visibleControls")
    visible_text = capture.get("visibleText")

    if "replayThrough" in capture:
        replay_through = capture.get("replayThrough")
        if (
            not isinstance(replay_through, int)
            or isinstance(replay_through, bool)
            or context.section != "manualSteps"
            or not 1 <= replay_through <= context.step_index
        ):
            result.errors.append(
                f"{context.label}: capture.replayThrough must name a manual step from 1 to this step."
            )
        if as_text(capture.get("source")).lower() not in {"webapp", "web-app", "app"}:
            result.errors.append(f"{context.label}: capture.replayThrough needs a webapp capture.")

    if "freshPage" in capture:
        if capture.get("freshPage") is not True or as_text(capture.get("source")).lower() not in {"webapp", "web-app", "app"}:
            result.errors.append(f"{context.label}: capture.freshPage must be true on a webapp capture.")
        if "replayThrough" in capture or as_text(capture.get("session")):
            result.errors.append(f"{context.label}: capture.freshPage excludes session and replayThrough.")

    if assert_app_state is not None and not isinstance(assert_app_state, dict):
        result.errors.append(f"{context.label}: capture.assertAppState must be an object.")
    elif isinstance(assert_app_state, dict):
        for state_path in assert_app_state:
            if not as_text(state_path):
                result.errors.append(
                    f"{context.label}: capture.assertAppState paths must be non-empty strings."
                )

    if visible_controls is not None and not isinstance(visible_controls, list):
        result.errors.append(f"{context.label}: capture.visibleControls must be an array.")
    elif isinstance(visible_controls, list):
        for index, spec in enumerate(visible_controls, start=1):
            if not isinstance(spec, dict):
                result.errors.append(
                    f"{context.label}: visibleControls entry {index} must be an object."
                )
                continue
            if not (as_text(spec.get("label")) or as_text(spec.get("selector"))):
                result.errors.append(
                    f"{context.label}: visibleControls entry {index} needs label or selector."
                )
            if (
                data_dependent
                and "value" not in spec
                and "checked" not in spec
                and "pressed" not in spec
            ):
                result.errors.append(
                    f"{context.label}: visibleControls entry {index} needs value, checked, or pressed."
                )

    if visible_text is not None and not isinstance(visible_text, list):
        result.errors.append(f"{context.label}: capture.visibleText must be an array.")
    elif isinstance(visible_text, list):
        for index, text in enumerate(visible_text, start=1):
            if not as_text(text):
                result.errors.append(
                    f"{context.label}: visibleText entry {index} must be a non-empty string."
                )

    if not data_dependent:
        return

    source = as_text(capture.get("source")).lower()
    if source not in {"webapp", "web-app", "app"}:
        result.errors.append(
            f"{context.label}: data-dependent capture source must be webapp."
        )

    capture_session_ref = as_text(capture.get("session"))
    sample_session_ref = as_text(sample.get("session"))
    if not capture_session_ref:
        result.errors.append(
            f"{context.label}: data-dependent capture must name its session explicitly."
        )
    elif not sample_session_ref:
        result.errors.append(
            f"{context.label}: Gallery example is missing its session reference."
        )
    else:
        capture_session = resolve_gallery_reference(capture_session_ref)
        sample_session = resolve_gallery_reference(sample_session_ref)
        if capture_session is None or not capture_session.is_file():
            result.errors.append(
                f"{context.label}: capture session does not exist: {capture_session_ref}"
            )
        elif capture_session != sample_session:
            result.errors.append(
                f"{context.label}: capture session must match the Gallery example session."
            )

    if not isinstance(assert_app_state, dict) or not assert_app_state:
        result.errors.append(
            f"{context.label}: data-dependent capture needs non-empty assertAppState."
        )
    if not (
        isinstance(visible_controls, list)
        and visible_controls
        or isinstance(visible_text, list)
        and visible_text
    ):
        result.errors.append(
            f"{context.label}: data-dependent capture must declare visibleControls or visibleText."
        )


def add_step_text_validation(
    result: ValidationResult,
    context: TutorialContext,
    step: dict[str, Any],
) -> None:
    """No step tells the reader to use the uploader's Session file chooser, which captures hide;
    no capture finds its card by text or an element by position, or does what the tool does for
    every capture: read the records, hide the Settings pane under a popup."""

    if SESSION_FILE_CHOOSER_NAME.lower() in json.dumps(step, ensure_ascii=False).lower():
        result.errors.append(
            f"{context.label}: names {SESSION_FILE_CHOOSER_NAME}, which tutorial captures hide: "
            "the reader uploads the files and has no Session."
        )
    captures = [
        operation["capture"]
        for operation in as_array(step.get("operations"))
        if isinstance(operation, dict) and isinstance(operation.get("capture"), dict)
    ]
    open_selects = [
        spec
        for capture in captures
        for spec in as_array(capture.get("openSelect"))
        + ([capture["openSelect"]] if isinstance(capture.get("openSelect"), dict) else [])
        if isinstance(spec, dict)
    ]
    if any("cardTitle" in capture for capture in captures):
        result.errors.append(
            f"{context.label}: a capture finds its card by text (cardTitle), which matches the first card "
            "that names it anywhere; give selector \".card[data-settings-section='<key>']\"."
        )
    if any("cardTitle" in spec or "selectIndex" in spec for spec in open_selects):
        result.errors.append(
            f"{context.label}: an openSelect finds its select by card text and position (cardTitle, "
            "selectIndex); give its selector, for example \"select[aria-label='Palette']\"."
        )
    actions = [json.dumps(action, ensure_ascii=False) for capture in captures for action in as_array(capture.get("actions"))]
    if any(CAPTURE_POSITION_PICK.search(action) for action in actions):
        result.errors.append(
            f"{context.label}: a capture action finds an element by position, which a new sibling "
            "shifts; find it by its id, data attribute or accessible name."
        )
    if any(".settings-pane" in action for action in actions):
        result.errors.append(
            f"{context.label}: a capture action changes the Settings pane; the tool hides it under "
            "a popup for every capture (hide_settings_pane_under_popup)."
        )
    if any(read in action for action in actions for read in CAPTURE_RECORD_READS):
        result.errors.append(
            f"{context.label}: a capture action reads the records, which the tool does for every "
            "capture after the replay (inspect_shown_records)."
        )


def add_step_replay_validation(
    result: ValidationResult,
    context: TutorialContext,
    step: dict[str, Any],
) -> None:
    if "replay" not in step:
        return
    replay = step.get("replay")
    if context.section != "manualSteps" or not isinstance(replay, list):
        result.errors.append(f"{context.label}: replay must be an array on a manual step.")
        return
    for index, action in enumerate(replay, start=1):
        action_type = as_text(action.get("type")) if isinstance(action, dict) else ""
        if action_type not in CAPTURE_ACTION_TYPES or action_type == "resetSettings":
            result.errors.append(
                f"{context.label}: replay action {index} has an unsupported type: {action_type or action!r}"
            )
        elif action_type == "upload":
            try:
                resolve_upload_file(action)
            except ValueError as error:
                result.errors.append(f"{context.label}: replay action {index} {error}")


def resolve_upload_file(action: dict[str, Any]) -> Path:
    """The file an `upload` action chooses: a file under gbdraw/web, named by its path from there.

    The Gallery serves that folder, so the path is the one a reader downloads (`gallery/files/...`).
    """

    reference = as_text(action.get("path"))
    web_root = WEB_ROOT.resolve()
    path = (web_root / reference).resolve()
    if not as_text(action.get("selector")):
        raise ValueError("upload needs the selector of the button that opens the file chooser.")
    if not reference or not path.is_relative_to(web_root) or not path.is_file():
        raise ValueError(f"upload path must name a file under gbdraw/web: {reference or '(none)'}")
    return path


def validate_tutorial_media(
    samples: list[dict[str, Any]],
    *,
    max_width: int,
    max_height: int,
    max_file_size_kb: int,
) -> ValidationResult:
    result = ValidationResult(errors=[], warnings=[])

    for sample in samples:
        tutorial_ref = as_text(sample.get("tutorial"))
        tutorial_path = resolve_gallery_reference(tutorial_ref)
        if tutorial_path is None or not tutorial_path.exists():
            result.errors.append(f"{sample.get('id')}: tutorial JSON does not exist: {tutorial_ref}")
            continue
        tutorial = read_json(tutorial_path)

        example_operation_count = 0
        for context, step in iter_step_contexts(sample, tutorial):
            if (
                context.section in OPERATION_REQUIRED_SECTIONS
                and not as_array(step.get("operations"))
                and not media_entries(step.get("media"))
            ):
                result.errors.append(f"{context.label}: step is missing operations.")
            add_step_text_validation(result, context, step)
            add_step_replay_validation(result, context, step)
            for entry in media_entries(step.get("media")):
                add_media_validation(
                    result,
                    context,
                    entry,
                    operation_media=False,
                    max_width=max_width,
                    max_height=max_height,
                    max_file_size_kb=max_file_size_kb,
                )

        for context, operation in iter_operation_contexts(sample, tutorial):
            result.operation_count += 1
            example_operation_count += 1
            add_operation_text_validation(result, context, operation)
            add_capture_contract_validation(result, context, sample, operation)
            operation_media = media_entries(operation.get("media"))
            if not operation_media:
                result.errors.append(f"{context.label}: operation is missing media.")
                continue
            for entry in operation_media:
                add_media_validation(
                    result,
                    context,
                    entry,
                    operation_media=True,
                    max_width=max_width,
                    max_height=max_height,
                    max_file_size_kb=max_file_size_kb,
                )

        if example_operation_count == 0:
            result.warnings.append(
                f"{sample.get('id')}: no operation entries yet. Migration is still incomplete for this tutorial."
            )

    return result


@contextlib.contextmanager
def serve_web_root():
    def handler(*args, **kwargs):
        return QuietGalleryRequestHandler(
            *args,
            directory=str(WEB_ROOT),
            **kwargs,
        )

    with socketserver.TCPServer(("127.0.0.1", 0), handler) as httpd:
        port = httpd.server_address[1]
        thread = threading.Thread(target=httpd.serve_forever, daemon=True)
        thread.start()
        try:
            yield f"http://127.0.0.1:{port}"
        finally:
            httpd.shutdown()
            thread.join(timeout=5)


def output_path_for_operation(operation: dict[str, Any]) -> Path | None:
    entries = media_entries(operation.get("media"))
    if not entries:
        return None
    source = as_text(entries[0].get("src")) or as_text(entries[0].get("href"))
    path = resolve_gallery_reference(source)
    if path is None:
        return None
    return path


def import_playwright():
    try:
        from playwright.sync_api import sync_playwright
    except ImportError as exc:  # pragma: no cover - dependency guard
        raise RuntimeError(
            "Python Playwright is required. Install it with `python -m pip install playwright` "
            "and run `python -m playwright install chromium`."
        ) from exc
    return sync_playwright


def wait_for_app_shell(page) -> None:
    page.wait_for_function(
        """
        () => {
          const app = window.__GBDRAW_APP__;
          if (!app) return false;
          return Object.keys(app.paletteDefinitions || {}).length > 0;
        }
        """,
        timeout=120000,
    )


def load_web_app_session(page, session_path: Path) -> None:
    page.locator('input[accept*=".gbdraw"]').set_input_files(str(session_path))
    # A loaded Session shows the notice `Session loaded: <file name>`.
    try:
        page.locator("[data-app-notice]").filter(has_text="Session loaded").wait_for(
            state="visible", timeout=180_000
        )
    except Exception as error:
        raise RuntimeError(f"Could not load capture session: {error}") from error
    if page.locator("[data-app-notice] li").count():
        raise RuntimeError(
            "Could not load capture session: "
            f"{page.locator('[data-app-notice]').inner_text()}"
        )
    page.wait_for_function(
        """
        () => {
          const app = window.__GBDRAW_APP__;
          return Boolean(
            app &&
            Array.isArray(app.results) &&
            app.results.length > 0
          );
        }
        """,
        timeout=120000,
    )


# Whether the shown drawing's input records are read, as a reader's upload reads them: a
# Circular drawing's Source records count, a Linear card's Record select. A drawing without
# input files has nothing to read. A loaded Session defers the read ("Records not inspected").
SHOWN_RECORDS_READ_SCRIPT = """
() => {
  const app = window.__GBDRAW_APP__;
  if (app.mode === 'circular') {
    const discovery = app.circularRecordDiscoveryState;
    return !discovery.hasInput || discovery.status === 'ready';
  }
  return app.linearSeqs
    .filter((seq) => app.shownInputComplete('linear', seq))
    .every((seq) => !app.linearRecordSelectorDisabled(seq));
}
"""


def dismiss_session_notice(page) -> None:
    """Close the `Session loaded` notice with its Dismiss button (OV-3369).

    A tutorial reader uploaded the files and has no Session, so no capture shows the notice. It
    stands above the Result: a popup that a capture places at a fixed height met the Preview
    toolbar the notice moved down. Reset Settings, which a replay starts with, closes it too.
    """

    notice = page.locator("[data-app-notice]")
    notice.locator(".app-notice-dismiss").click()
    notice.wait_for(state="hidden")


def inspect_shown_records(page) -> None:
    """Read the shown drawing's records, as the reader's upload did (OV-3139).

    A tutorial reader uploads the files, and an upload reads their records ("1 source record
    inspected", the Record select naming the record). A loaded Session leaves them unread
    until **Inspect source records** (a Linear card's **Inspect records**), which changes
    neither the Result nor Undo history; this presses that action through the app.
    """

    page.evaluate(
        """
        async () => {
          const app = window.__GBDRAW_APP__;
          await (app.mode === 'circular' ? app.inspectCircularSourceRecords() : app.refreshLinearRecordSelectors());
        }
        """
    )
    try:
        page.wait_for_function(SHOWN_RECORDS_READ_SCRIPT, timeout=DEFAULT_GENERATE_TIMEOUT_MS)
    except Exception as error:
        raise RuntimeError(f"Could not read the records of the capture session: {error}") from error
    wait_for_record_scroll_anchor(page)


# A record read or a Record choice that opens the Circular single-record panel puts the Settings
# pane back where it was scrolled once the panel has painted. The capture awaits the app's own
# restore, so it never undoes a scroll the capture makes next (OV-3547, HmmtDNA basic 3-1: the
# Layout card ended up below the page). A page without the app has none.
RECORD_SCROLL_ANCHOR_SCRIPT = (
    "async () => { const app = window.__GBDRAW_APP__; if (app) await app.waitForCircularRecordScrollAnchor(); }"
)


def wait_for_record_scroll_anchor(page) -> None:
    page.evaluate(RECORD_SCROLL_ANCHOR_SCRIPT)


def assert_capture_records_read(page) -> None:
    """A webapp crop shows the records read, as the reader who uploaded the files sees them."""

    if not page.evaluate(SHOWN_RECORDS_READ_SCRIPT):
        raise RuntimeError(
            "capture shows records not inspected: the reader who uploaded the files sees them read."
        )


def apply_capture_action(page, action: dict[str, Any]) -> None:
    action_type = as_text(action.get("type")) or as_text(action.get("action"))
    selector = as_text(action.get("selector"))
    value = action.get("value")

    if action_type == "click":
        click_options = {"force": True} if action.get("force") is True else {}
        page.locator(selector).first.click(**click_options)
    elif action_type == "fill":
        page.locator(selector).fill("" if value is None else str(value))
    elif action_type == "select":
        page.locator(selector).select_option(str(value))
    elif action_type == "check":
        page.locator(selector).check()
    elif action_type == "uncheck":
        page.locator(selector).uncheck()
    elif action_type == "waitForSelector":
        page.locator(selector).first.wait_for(state=as_text(action.get("state")) or "visible")
    elif action_type == "waitForFunction":
        script = as_text(action.get("script"))
        timeout = action.get("timeout")
        if script:
            page.wait_for_function(script, timeout=int(timeout) if timeout else None)
    elif action_type == "evaluate":
        script = as_text(action.get("script"))
        if script:
            page.evaluate(script)
    elif action_type == "wait":
        page.wait_for_timeout(int(action.get("ms") or 250))
    elif action_type == "resetSettings":
        page.get_by_role("button", name="Reset Settings", exact=True).click()
        confirm = page.get_by_role("dialog", name="Reset settings?")
        confirm.get_by_role("button", name="Reset all drawings", exact=True).click()
        confirm.wait_for(state="hidden")
        # Reset redraws the shown Result ("Updating diagram"); the mode buttons, Undo, and
        # other History-backed controls stay disabled until that update ends.
        page.wait_for_function(
            "() => { const app = window.__GBDRAW_APP__; return Boolean(app) && !app.processing && !app.labelReflowProcessing; }",
            timeout=DEFAULT_GENERATE_TIMEOUT_MS,
        )
        # History can stay busy for many seconds after that (OV-1700); every replay starts
        # with a drawing-tab click, so wait until the drawing tabs are enabled.
        page.wait_for_function(
            "() => [...document.querySelectorAll('.drawing-bar [role=\"tab\"]')].every((tab) => !tab.disabled)",
            timeout=DEFAULT_GENERATE_TIMEOUT_MS,
        )
    elif action_type == "generate":
        generate_and_wait(page, int(action.get("timeout") or DEFAULT_GENERATE_TIMEOUT_MS))
    elif action_type == "upload":
        path = resolve_upload_file(action)
        with page.expect_file_chooser() as chooser:
            page.locator(selector).first.click()
        chooser.value.set_files(str(path))
    else:
        raise RuntimeError(f"Unsupported capture action: {action_type}")
    # The next action measures what this one moved: the Settings pane the app puts back after a
    # Record choice, a zoom or pan that moves the preview in a 0.2 s transition, where a box read
    # before it ends is the old one (OV-2280).
    wait_for_record_scroll_anchor(page)
    wait_for_css_transitions(page)


def generate_and_wait(page, timeout_ms: int) -> None:
    previous_run = page.evaluate("() => String(window.__GBDRAW_APP__?.lastRunInfo?.startedAtIso || '')")
    page.get_by_role("button", name="Generate Diagram", exact=True).click()
    page.wait_for_function(
        """
        (previous) => {
          const app = window.__GBDRAW_APP__;
          if (!app || app.processing) return false;
          const run = String(app.lastRunInfo?.startedAtIso || '');
          return Boolean(app.errorLog) || Boolean(run && run !== previous && (app.results || []).length > 0);
        }
        """,
        arg=previous_run,
        timeout=timeout_ms,
    )
    error = page.evaluate("() => window.__GBDRAW_APP__?.errorLog?.summary || window.__GBDRAW_APP__?.errorLog || null")
    if error:
        raise RuntimeError(f"Generate Diagram failed: {error}")


def replay_tutorial_steps(page, sample: dict[str, Any], through: int) -> None:
    """Reset Settings, then perform the reader's actions of manual steps 1..through."""

    tutorial_path = resolve_gallery_reference(as_text(sample.get("tutorial")))
    if tutorial_path is None:
        raise RuntimeError(f"Tutorial path for {sample.get('id')} is not local.")
    steps = as_array(read_json(tutorial_path).get("manualSteps"))[:through]
    apply_capture_action(page, {"type": "resetSettings"})
    for step in steps:
        for action in as_array(step.get("replay") if isinstance(step, dict) else None):
            apply_capture_action(page, action)


def open_capture_details(page, titles: Any) -> None:
    detail_titles = [title for title in as_array(titles) if as_text(title)]
    if not detail_titles:
        return

    page.evaluate(
        """
        (titles) => {
          const normalize = (value) => String(value || '').replace(/\\s+/g, ' ').trim().toLowerCase();
          const wanted = titles.map(normalize);
          document.querySelectorAll('details').forEach((details) => {
            const summary = details.querySelector('summary');
            const label = normalize(summary?.textContent || '');
            if (wanted.some((title) => label.includes(title))) details.open = true;
          });
        }
        """,
        detail_titles,
    )


def apply_open_select_overlays(page, open_select: Any) -> None:
    specs = as_array(open_select) if isinstance(open_select, list) else ([open_select] if isinstance(open_select, dict) else [])

    page.evaluate(
        """
        ({ specs, maxVisibleOptions }) => {
          const normalize = (value) => String(value || '').replace(/\\s+/g, ' ').trim();
          const normalizeKey = (value) => normalize(value).toLowerCase();

          document.querySelectorAll('.gbdraw-capture-select-menu').forEach((node) => node.remove());
          document.querySelectorAll('.gbdraw-capture-select-control').forEach((node) => {
            node.classList.remove('gbdraw-capture-select-control');
          });
          if (!specs.length) return;
          let style = document.getElementById('gbdraw-capture-select-style');
          if (!style) {
            style = document.createElement('style');
            style.id = 'gbdraw-capture-select-style';
            style.textContent = `
              .gbdraw-capture-select-control { ${CONTROL_FRAME} }
              .gbdraw-capture-select-menu {
                position: absolute;
                z-index: 2147483647;
                box-sizing: border-box;
                padding: 4px;
                border: 1px solid #d1d5db;
                border-radius: 6px;
                background: #ffffff;
                box-shadow: 0 12px 28px rgba(15, 23, 42, 0.18), 0 2px 8px rgba(15, 23, 42, 0.12);
                color: #334155;
                font-family: Inter, "Noto Sans JP", system-ui, sans-serif;
                pointer-events: none;
              }
              .gbdraw-capture-select-option {
                min-height: 24px;
                display: flex;
                align-items: center;
                padding: 3px 8px;
                border-radius: 4px;
                font-size: 12px;
                line-height: 1.25;
                white-space: nowrap;
              }
              .gbdraw-capture-select-option.is-selected {
                ${OPTION_FRAME}
                background: #fff7ed;
                color: #1d4ed8;
                font-weight: 800;
              }
            `;
            document.head.appendChild(style);
          }

          const selectNearLabel = (labelText, nth = 0) => {
            const labelKey = normalizeKey(labelText);
            const matches = [];
            document.querySelectorAll('label').forEach((label) => {
              if (!normalizeKey(label.textContent).includes(labelKey)) return;
              if (label.control?.tagName === 'SELECT') {
                matches.push(label.control);
                return;
              }
              const local = label.parentElement?.querySelector('select');
              if (local) {
                matches.push(local);
                return;
              }
              const container = label.closest('div, section, details, .card');
              const nearby = container?.querySelector('select');
              if (nearby) matches.push(nearby);
            });
            return matches[nth] || null;
          };

          const resolveSelect = (spec) => {
            if (spec.selector) {
              const selected = document.querySelector(spec.selector);
              if (selected?.tagName === 'SELECT') return selected;
            }
            if (spec.label) return selectNearLabel(spec.label, Number(spec.nth || 0));
            return null;
          };

          specs.forEach((spec, index) => {
            const select = resolveSelect(spec);
            if (!select) {
              throw new Error(`Could not find select for openSelect ${spec.label || spec.selector || index}.`);
            }
            if (spec.hideFollowingSiblings === true) {
              let sibling = select.parentElement?.nextElementSibling || null;
              while (sibling) {
                sibling.style.visibility = 'hidden';
                sibling = sibling.nextElementSibling;
              }
            }
            if (index === 0) select.scrollIntoView({ block: 'center', inline: 'nearest' });
          });

          document.querySelectorAll('.app-main, .settings-pane, .settings-scroll').forEach((element) => {
            element.scrollLeft = 0;
          });
          document.documentElement.scrollLeft = 0;
          document.body.scrollLeft = 0;
          window.scrollTo({ left: 0, top: window.scrollY });

          specs.forEach((spec, index) => {
            const select = resolveSelect(spec);
            if (!select) {
              throw new Error(`Could not find select for openSelect ${spec.label || spec.selector || index}.`);
            }

            const optionLabels = spec.options && Array.isArray(spec.options) && spec.options.length
              ? spec.options.map(normalize)
              : Array.from(select.options).map((option) => normalize(option.textContent));
            const selectedLabel = normalize(spec.selectedLabel);
            if (!optionLabels.some((label) => normalizeKey(label) === normalizeKey(selectedLabel))) {
              throw new Error(`Selected label "${selectedLabel}" is not in openSelect options: ${optionLabels.join(', ')}`);
            }

            const selectedOption = Array.from(select.options).find(
              (option) => normalizeKey(option.textContent) === normalizeKey(selectedLabel)
            );
            if (selectedOption && select.value !== selectedOption.value) {
              select.value = selectedOption.value;
              select.dispatchEvent(new Event('change', { bubbles: true }));
            }

            select.classList.add('gbdraw-capture-select-control');
            const rect = select.getBoundingClientRect();
            const menu = document.createElement('div');
            menu.className = 'gbdraw-capture-select-menu';
            menu.dataset.captureOpenSelect = String(index);
            const minWidth = Number(spec.minWidth || 0);
            menu.style.minWidth = `${Math.max(rect.width, minWidth, 120)}px`;

            optionLabels.forEach((label) => {
              const option = document.createElement('div');
              option.className = 'gbdraw-capture-select-option';
              if (normalizeKey(label) === normalizeKey(selectedLabel)) option.classList.add('is-selected');
              option.textContent = label;
              menu.appendChild(option);
            });

            document.body.appendChild(menu);
            const rows = Array.from(menu.children);
            if (rows.length > maxVisibleOptions) {
              const rowHeight = rows[0].offsetHeight;
              const menuStyle = getComputedStyle(menu);
              const chrome = ['paddingTop', 'paddingBottom', 'borderTopWidth', 'borderBottomWidth']
                .reduce((sum, key) => sum + parseFloat(menuStyle[key] || '0'), 0);
              menu.style.maxHeight = `${rowHeight * maxVisibleOptions + chrome}px`;
              menu.style.overflowX = 'hidden';
              menu.style.overflowY = 'auto';
              const selectedIndex = rows.findIndex((row) => row.classList.contains('is-selected'));
              const firstRow = Math.max(0, selectedIndex - Math.floor(maxVisibleOptions / 2));
              menu.scrollTop = firstRow * rowHeight;
            }
            const placement = normalizeKey(spec.placement || 'below');
            const gap = Number(spec.gap ?? 4);
            if (placement === 'right') {
              menu.style.left = `${window.scrollX + rect.right + gap}px`;
              menu.style.top = `${window.scrollY + rect.top}px`;
            } else if (placement === 'left') {
              menu.style.left = `${window.scrollX + rect.left - menu.offsetWidth - gap}px`;
              menu.style.top = `${window.scrollY + rect.top}px`;
            } else if (placement === 'above') {
              menu.style.left = `${window.scrollX + rect.left}px`;
              menu.style.top = `${window.scrollY + rect.top - menu.offsetHeight - gap}px`;
            } else {
              menu.style.left = `${window.scrollX + rect.left}px`;
              menu.style.top = `${window.scrollY + rect.bottom + gap}px`;
            }
          });
        }
        """.replace("${CONTROL_FRAME}", highlight_frame_css(-HIGHLIGHT_FRAME_WIDTH, HIGHLIGHT_FRAME_GLOW, inset=True))
        .replace("${OPTION_FRAME}", highlight_frame_css(-2.0, 0.0, inset=True)),
        {"specs": specs, "maxVisibleOptions": OPEN_SELECT_MAX_VISIBLE_OPTIONS},
    )
    if not specs:
        return
    page.wait_for_timeout(100)


def resolve_capture_target(page, capture: dict[str, Any]):
    selector = as_text(capture.get("selector"))
    if selector:
        return page.locator(selector).first

    return None


def normalize_crop_padding(value: Any, default: int = 12) -> dict[str, float]:
    if isinstance(value, dict):
        return {
            "top": float(value.get("top", default)),
            "right": float(value.get("right", default)),
            "bottom": float(value.get("bottom", default)),
            "left": float(value.get("left", default)),
        }
    if value is None:
        amount = float(default)
    else:
        amount = float(value)
    return {"top": amount, "right": amount, "bottom": amount, "left": amount}


def page_scroll_dimensions(page) -> dict[str, float]:
    dimensions = page.evaluate(
        """
        () => ({
          width: Math.max(document.documentElement.scrollWidth, document.body?.scrollWidth || 0),
          height: Math.max(document.documentElement.scrollHeight, document.body?.scrollHeight || 0),
        })
        """
    )
    return {"width": float(dimensions["width"]), "height": float(dimensions["height"])}


def padded_clip_from_boxes(
    boxes: list[dict[str, Any]],
    padding: dict[str, float],
    dimensions: dict[str, float],
) -> dict[str, float]:
    if not boxes:
        raise RuntimeError("Cannot compute screenshot clip without visible boxes.")

    min_x = min(float(box["x"]) for box in boxes)
    min_y = min(float(box["y"]) for box in boxes)
    max_x = max(float(box["x"]) + float(box["width"]) for box in boxes)
    max_y = max(float(box["y"]) + float(box["height"]) for box in boxes)
    # A crop trims its padding at the page edge, never the element it shows (OV-2280).
    if min_x < -1 or min_y < -1 or max_x > dimensions["width"] + 1 or max_y > dimensions["height"] + 1:
        raise RuntimeError(
            f"The crop would cut its element: it spans x {min_x:.0f}..{max_x:.0f}, "
            f"y {min_y:.0f}..{max_y:.0f} on a {dimensions['width']:.0f} x {dimensions['height']:.0f} page."
        )

    x = max(0.0, min_x - padding["left"])
    y = max(0.0, min_y - padding["top"])
    right = min(dimensions["width"], max_x + padding["right"])
    bottom = min(dimensions["height"], max_y + padding["bottom"])
    return {
        "x": x,
        "y": y,
        "width": max(1.0, right - x),
        "height": max(1.0, bottom - y),
    }


def capture_selector_clip(page, target_locator, padding: dict[str, float]) -> dict[str, float]:
    box = target_locator.evaluate(
        """
        (element) => {
          const rect = element.getBoundingClientRect();
          return {
            x: rect.left + window.scrollX,
            y: rect.top + window.scrollY,
            width: rect.width,
            height: rect.height,
          };
        }
        """
    )
    return padded_clip_from_boxes([box], padding, page_scroll_dimensions(page))


def capture_open_select_clip(
    page,
    padding: dict[str, float],
    include_selectors: list[str] | None = None,
) -> dict[str, float]:
    include_selectors = include_selectors or []
    boxes = page.evaluate(
        """
        (includeSelectors) => {
          const selectors = [
            '.gbdraw-capture-select-control',
            '.gbdraw-capture-select-menu',
            ...includeSelectors,
          ];
          return selectors.flatMap((selector) => Array.from(document.querySelectorAll(selector))).map((element) => {
            const rect = element.getBoundingClientRect();
            return {
              x: rect.left + window.scrollX,
              y: rect.top + window.scrollY,
              width: rect.width,
              height: rect.height,
            };
          }).filter((box) => box.width > 0 && box.height > 0);
        }
        """,
        include_selectors,
    )
    if not boxes:
        raise RuntimeError("openSelect crop requested, but no capture select overlay exists.")
    return padded_clip_from_boxes(boxes, padding, page_scroll_dimensions(page))


def capture_multi_selector_clip(
    page,
    selectors: list[str],
    padding: dict[str, float],
) -> dict[str, float]:
    boxes = page.evaluate(
        """
        (selectors) => selectors.flatMap((selector) => {
          try {
            return Array.from(document.querySelectorAll(selector));
          } catch {
            return [];
          }
        }).map((element) => {
          const rect = element.getBoundingClientRect();
          return {
            x: rect.left + window.scrollX,
            y: rect.top + window.scrollY,
            width: rect.width,
            height: rect.height,
          };
        }).filter((box) => box.width > 0 && box.height > 0)
        """,
        selectors,
    )
    if not boxes:
        raise RuntimeError("clipSelectors crop requested, but no visible selector boxes were found.")
    return padded_clip_from_boxes(boxes, padding, page_scroll_dimensions(page))


# The lines of text and the controls drawn in a region of the page, in page coordinates, each
# trimmed to the part a reader sees: a hit test at points over it must find it, not a popup or
# menu above it, the pane that scrolls it out of view, its own clip (a screen-reader line), or
# nothing. For the measure, hit tests see the measured text and every painted layer the page
# leaves without pointer events (pointer-events), and the parked pointer hovers nothing.
class TextBox(TypedDict):
    """A line of text or a control a reader sees, or the seen part of one, in page CSS px
    (VISIBLE_TEXT_BOXES_SCRIPT). The edge rule, the crop guard and the highlight frame read it."""

    x: float
    y: float
    width: float
    height: float
    text: str


VISIBLE_TEXT_BOXES_SCRIPT = """
({ region, exclude }) => {
  const normalize = (value) => String(value || '').replace(/\\s+/g, ' ').trim();
  const left = region.x - window.scrollX;
  const top = region.y - window.scrollY;
  const right = left + region.width;
  const bottom = top + region.height;
  const shown = (element) => element.checkVisibility({ opacityProperty: true, visibilityProperty: true });
  // What covers a point is what hides it (M1, L1): the first layer of the hit list whose paint is
  // opaque there, passing over layers that paint nothing there (a clear full-cover overlay, an svg
  // root, a group, the inside of a border-only frame) and layers whose paint lets the text show
  // through (a tint such as a selection marquee, a layer at an opacity below 1).
  const SHAPES = new Set(['path', 'rect', 'circle', 'ellipse', 'line', 'polyline', 'polygon', 'text',
    'tspan', 'textPath', 'image', 'use', 'foreignObject']);
  const REPLACED = new Set(['IMG', 'CANVAS', 'VIDEO', 'IFRAME', 'INPUT', 'SELECT', 'TEXTAREA', 'BUTTON']);
  const SIDES = ['Top', 'Right', 'Bottom', 'Left'];
  // The alpha of a computed colour: `transparent`, `rgba(r, g, b, a)` or a colour function with `/ a`.
  const alphaOf = (color) => {
    if (color === 'transparent') return 0;
    const slash = color.lastIndexOf('/');
    if (slash >= 0) {
      const value = color.slice(slash + 1).trim();
      return value.endsWith('%') ? parseFloat(value) / 100 : parseFloat(value);
    }
    return color.startsWith('rgba(') ? parseFloat(color.slice(color.lastIndexOf(',') + 1)) : 1;
  };
  // A layer paints when it draws anything a reader sees: a shape, a replaced element, a background
  // or a border (a control's pill, a layer to make hittable for the measure).
  const paintsCache = new Map();
  const paints = (element) => {
    if (paintsCache.has(element)) return paintsCache.get(element);
    const style = getComputedStyle(element);
    let result;
    if (style.opacity === '0') result = false;
    else if (element instanceof SVGElement) result = SHAPES.has(element.tagName);
    else result = REPLACED.has(element.tagName) || alphaOf(style.backgroundColor) > 0
      || style.backgroundImage !== 'none'
      || SIDES.some((side) => parseFloat(style[`border${side}Width`]) > 0 && alphaOf(style[`border${side}Color`]) > 0);
    paintsCache.set(element, result);
    return result;
  };
  // Where an HTML layer's paint is opaque: nowhere when it or an ancestor is drawn at an opacity
  // below 1; its whole box for a replaced element or a background image; the box its background
  // clips to for a background of alpha 1; the strip of each side whose border has alpha 1.
  const opaqueCache = new Map();
  const opaqueParts = (element) => {
    if (opaqueCache.has(element)) return opaqueCache.get(element);
    const parts = [];
    const box = element.getBoundingClientRect();
    let solid = true;
    for (let layer = element; layer && solid; layer = layer.parentElement) solid = getComputedStyle(layer).opacity === '1';
    if (solid) {
      const style = getComputedStyle(element);
      const rect = box;
      const [top, right, bottom, left] = SIDES.map((side) => parseFloat(style[`border${side}Width`]) || 0);
      if (REPLACED.has(element.tagName) || style.backgroundImage !== 'none') parts.push(rect);
      else if (alphaOf(style.backgroundColor) === 1) {
        const inset = style.backgroundClip === 'border-box' ? [0, 0, 0, 0] : [top, right, bottom, left].map(
          (border, index) => border + (style.backgroundClip === 'content-box'
            ? parseFloat(style[`padding${SIDES[index]}`]) || 0 : 0));
        parts.push(new DOMRect(rect.left + inset[3], rect.top + inset[0],
          rect.width - inset[1] - inset[3], rect.height - inset[0] - inset[2]));
      }
      const strips = [
        new DOMRect(rect.left, rect.top, rect.width, top), new DOMRect(rect.right - right, rect.top, right, rect.height),
        new DOMRect(rect.left, rect.bottom - bottom, rect.width, bottom), new DOMRect(rect.left, rect.top, left, rect.height),
      ];
      SIDES.forEach((side, index) => {
        if (strips[index].width > 0 && strips[index].height > 0 && alphaOf(style[`border${side}Color`]) === 1) {
          parts.push(strips[index]);
        }
      });
    }
    opaqueCache.set(element, { box, parts });
    return opaqueCache.get(element);
  };
  // The hit test found `layer` at the point, which it may have snapped to the pixel grid (a point
  // 0.4 px left of a drawer hits the drawer): the point counts where it lies within the layer's box.
  const hides = (layer, x, y) => {
    if (layer instanceof SVGElement) return paints(layer);
    const { box, parts } = opaqueParts(layer);
    const px = Math.min(Math.max(x, box.left), box.right - 1e-6);
    const py = Math.min(Math.max(y, box.top), box.bottom - 1e-6);
    return parts.some((rect) => px >= rect.left && px < rect.right && py >= rect.top && py < rect.bottom);
  };
  // The boxes that clip an owner, its own included (a scrolled list, an svg, a screen-reader
  // line in a 1 px box or under `clip: rect(0 0 0 0)`): its text outside them is not seen. A
  // clip-path counts by its element's box; within the viewport the hit test below follows its shape.
  const clipsCache = new Map();
  const clipsOf = (owner) => {
    if (!clipsCache.has(owner)) {
      const clips = [];
      for (let element = owner; element && element !== document.body; element = element.parentElement) {
        const style = getComputedStyle(element);
        const rect = element.getBoundingClientRect();
        if (style.overflowX !== 'visible' || style.overflowY !== 'visible' || style.clipPath !== 'none') clips.push(rect);
        if (style.clip.startsWith('rect(') && ['absolute', 'fixed'].includes(style.position)) {
          const [clipTop, clipRight, clipBottom, clipLeft] = style.clip.slice(5, -1).split(/[,\\s]+/)
            .map((value, index) => (value === 'auto' ? [0, rect.width, rect.height, 0][index] : parseFloat(value)));
          clips.push(new DOMRect(rect.left + clipLeft, rect.top + clipTop,
            Math.max(0, clipRight - clipLeft), Math.max(0, clipBottom - clipTop)));
        }
      }
      clipsCache.set(owner, clips);
    }
    return clipsCache.get(owner);
  };
  // A point is seen where the owner itself is hit before anything that paints: where its text
  // paints, not where only an ancestor's box lies (a clipped line, the corner of a rotated label).
  // Outside the viewport no hit test answers: count the point as seen.
  const sees = (owner, x, y) => {
    if (clipsOf(owner).some((rect) => x < rect.left || x >= rect.right || y < rect.top || y >= rect.bottom)) return false;
    if (x < 0 || y < 0 || x >= window.innerWidth || y >= window.innerHeight) return true;
    for (const hit of document.elementsFromPoint(x, y)) {
      if (hit === owner || owner.contains(hit)) return true;
      if (hides(hit, x, y)) return false;
    }
    return false;
  };
  const candidates = [];
  const measured = new Set();
  const add = (owner, rect, text) => {
    if (rect.width > 0 && rect.height > 0 && rect.right > left && rect.left < right
        && rect.bottom > top && rect.top < bottom) {
      candidates.push({ owner, rect, text });
      measured.add(owner);
    }
  };
  // A text node's lines, in HTML and in SVG (an SVG label's box holds its rotated glyphs).
  const walker = document.createTreeWalker(document.body, NodeFilter.SHOW_TEXT);
  for (let node = walker.nextNode(); node; node = walker.nextNode()) {
    const owner = node.parentElement;
    const text = normalize(node.nodeValue);
    if (!text || !owner || owner.closest('script, style, template, noscript, select, textarea, title')) continue;
    if (!shown(owner)) continue;
    const range = document.createRange();
    range.selectNodeContents(node);
    for (const rect of range.getClientRects()) add(owner, rect, text);
  }
  const CONTROLS = 'input:not([type="hidden"]), select, textarea, button, '
    + '[role="button"], [role="checkbox"], [role="radio"], [role="switch"], [role="tab"]';
  // A control's visible shape may be drawn by a wrapper (M2: the pill of a Palette clock note
  // around its help button): the nearest ancestor within WRAP px that paints a border or a
  // background and holds no other control.
  const WRAP = 8;
  const shapeOf = (control) => {
    const own = control.getBoundingClientRect();
    let element = control.parentElement;
    for (let depth = 0; element && depth < 4; depth += 1, element = element.parentElement) {
      const rect = element.getBoundingClientRect();
      if (rect.left < own.left - WRAP || rect.top < own.top - WRAP || rect.right > own.right + WRAP
          || rect.bottom > own.bottom + WRAP || element.querySelectorAll(CONTROLS).length > 1) break;
      if (!(element instanceof SVGElement) && paints(element)) return { owner: element, rect: new DOMRect(
        Math.min(own.left, rect.left), Math.min(own.top, rect.top),
        Math.max(own.right, rect.right) - Math.min(own.left, rect.left),
        Math.max(own.bottom, rect.bottom) - Math.min(own.top, rect.top)) };
    }
    return { owner: control, rect: own };
  };
  // The controls that show their value as text: a field's value, a list's chosen option.
  const VALUE_TEXT = 'input:not([type="checkbox"], [type="radio"], [type="range"], [type="color"], [type="image"]), '
    + 'select, textarea';
  for (const control of document.querySelectorAll(CONTROLS)) {
    // A drawn group that opens settings (role="button") is figure: its labels count as text.
    if (control instanceof SVGElement || !shown(control)) continue;
    const name = normalize(control.getAttribute('aria-label') || control.textContent || control.value);
    const label = `${control.tagName.toLowerCase()}: ${name}`;
    const { owner, rect } = shapeOf(control);
    if (control !== exclude && owner !== exclude) {
      add(owner, rect, label);
      continue;
    }
    // The excluded element's own shape (the control, or the pill that draws it) is what the frame
    // goes around; the value text it shows is text like any other: its content box (the border
    // box less border and padding) stays inside the frame.
    if (!control.matches(VALUE_TEXT)) continue;
    const style = getComputedStyle(control);
    const [padTop, padRight, padBottom, padLeft] = SIDES.map((side) => (parseFloat(style[`border${side}Width`]) || 0)
      + (parseFloat(style[`padding${side}`]) || 0));
    const own = control.getBoundingClientRect();
    add(control, new DOMRect(own.left + padLeft, own.top + padTop, own.width - padLeft - padRight,
      own.height - padTop - padBottom), label);
  }
  // The part of each box a reader sees, as rectangles that keep its shape: a label half under
  // a popup is an L, and one bounding box around it would cross edges the label does not reach.
  // Each run of seen cells down a column, and along a row, becomes a box (neighbouring runs
  // that match merge), so an edge between two cells of one run still crosses that run's box.
  const runs = (cells) => {
    const found = [];
    let first = -1;
    cells.forEach((value, index) => {
      if (value && first < 0) first = index;
      if (first >= 0 && (!value || index === cells.length - 1)) {
        found.push([first, value ? index : index - 1]);
        first = -1;
      }
    });
    return found;
  };
  const strips = (lines, across) => {
    const merged = [];
    lines.forEach((line, index) => {
      for (const [from, to] of runs(line)) {
        const open = merged.find((strip) => strip.from === from && strip.to === to && strip.last === index - 1);
        if (open) open.last = index;
        else merged.push({ first: index, last: index, from, to });
      }
    });
    return merged.map((strip) => (across
      ? { minColumn: strip.first, maxColumn: strip.last, minRow: strip.from, maxRow: strip.to }
      : { minRow: strip.first, maxRow: strip.last, minColumn: strip.from, maxColumn: strip.to }));
  };
  const whole = (rect, text) => [{
    x: rect.left + window.scrollX, y: rect.top + window.scrollY, width: rect.width, height: rect.height, text,
  }];
  // An owner the page leaves without pointer events (a figure label) is made hittable, so the
  // hit list says what lies above it, and so is every shown HTML layer without pointer events
  // that paints (a help tip, a hover summary, a capture's select menu), so it covers what lies
  // under it (L1). Nothing else changes its pointer events.
  const covers = [...document.body.querySelectorAll('*')].filter((element) => element instanceof HTMLElement
    && getComputedStyle(element).pointerEvents === 'none' && paints(element) && shown(element));
  for (const element of [...measured, ...covers]) element.setAttribute('data-gbdraw-measured', '');
  const style = document.createElement('style');
  style.textContent = '[data-gbdraw-measured] { pointer-events: auto !important; }';
  document.head.appendChild(style);
  try {
  return candidates.flatMap(({ owner, rect, text }) => {
    // A box seen at the nine points of a coarse grid is seen whole. A box partly covered is
    // sampled in cells of about 2 CSS px; a cell counts as seen when its centre or a corner is
    // (L1), so a seen part keeps the cells its cover only starts in.
    // Points stay half a CSS px inside the box, where its own hit test answers.
    const at = (fx, fy) => sees(owner, rect.left + Math.min(Math.max(fx * rect.width, 0.5), rect.width - 0.5),
      rect.top + Math.min(Math.max(fy * rect.height, 0.5), rect.height - 0.5));
    if ([0, 0.5, 1].every((fx) => [0, 0.5, 1].every((fy) => at(fx, fy)))) return whole(rect, text);
    const columns = Math.min(200, Math.max(3, Math.ceil(rect.width / 2)));
    const rows = Math.min(40, Math.max(3, Math.ceil(rect.height / 2)));
    const corner = Array.from({ length: rows + 1 }, (_, row) => Array.from({ length: columns + 1 },
      (__, column) => at(column / columns, row / rows)));
    const seen = Array.from({ length: rows }, (_, row) => Array.from({ length: columns }, (__, column) => (
      corner[row][column] || corner[row][column + 1] || corner[row + 1][column] || corner[row + 1][column + 1]
      || at((column + 0.5) / columns, (row + 0.5) / rows))));
    const cellWidth = rect.width / columns;
    const cellHeight = rect.height / rows;
    const byColumn = Array.from({ length: columns }, (_, column) => seen.map((line) => line[column]));
    const boxes = new Map();
    for (const part of [...strips(seen, false), ...strips(byColumn, true)]) {
      const box = {
        x: rect.left + part.minColumn * cellWidth + window.scrollX,
        y: rect.top + part.minRow * cellHeight + window.scrollY,
        width: (part.maxColumn - part.minColumn + 1) * cellWidth,
        height: (part.maxRow - part.minRow + 1) * cellHeight,
        text,
      };
      boxes.set(`${box.x},${box.y},${box.width},${box.height}`, box);
    }
    return [...boxes.values()];
  });
  } finally {
    style.remove();
    for (const element of [...measured, ...covers]) element.removeAttribute('data-gbdraw-measured');
  }
}
"""


def visible_text_boxes(page, region: dict[str, float], exclude: ElementHandle | None = None) -> list[TextBox]:
    """The lines of text and the controls a reader sees in ``region`` (VISIBLE_TEXT_BOXES_SCRIPT),
    less the shape of the control ``exclude`` (an element handle: a highlighted field, whose own box
    and the pill that draws it are what its frame goes around) when one is given. Text inside it counts,
    and so does the value a field, list or text area shows: its content box stands for it."""

    return page.evaluate(VISIBLE_TEXT_BOXES_SCRIPT, {"region": region, "exclude": exclude})


def park_pointer(page) -> None:
    """Move the pointer off the page (OV-3411), so a crop shows no hover that the last click of
    the replay or of the capture's actions left (BGC 7-2's Skip row)."""

    page.mouse.move(-1, -1)


def _crop_sides(clip: dict[str, float]) -> dict[str, float]:
    return {
        "left": clip["x"],
        "top": clip["y"],
        "right": clip["x"] + clip["width"],
        "bottom": clip["y"] + clip["height"],
    }


def _edge_crossings(
    sides: dict[str, float], boxes: list[TextBox], dimensions: dict[str, float], side: str
) -> list[TextBox]:
    """The boxes that ``side`` of the crop crosses; an edge at the page edge crosses none."""

    tolerance = CROP_EDGE_TOLERANCE
    vertical = side in {"left", "right"}
    position = sides[side]
    if position <= 0 or position >= dimensions["width" if vertical else "height"]:
        return []
    start, size, across_start, across_size = ("x", "width", "y", "height") if vertical else ("y", "height", "x", "width")
    low, high = (sides["top"], sides["bottom"]) if vertical else (sides["left"], sides["right"])
    return [
        box
        for box in boxes
        if float(box[start]) < position - tolerance
        and float(box[start]) + float(box[size]) > position + tolerance
        and float(box[across_start]) < high - tolerance
        and float(box[across_start]) + float(box[across_size]) > low + tolerance
    ]


def _cleared_edge(
    sides: dict[str, float],
    boxes: list[TextBox],
    dimensions: dict[str, float],
    side: str,
    origin: float,
    limit: float,
) -> float:
    """Where ``side`` crosses no box: the nearer of the first clear position outward, at most
    CROP_EDGE_REACH from ``origin``, and the first one inward, never past ``limit`` (the
    subject); outward on a tie; where it is when neither is free."""

    start, size = ("x", "width") if side in {"left", "right"} else ("y", "height")
    extent = dimensions["width" if side in {"left", "right"} else "height"]
    outward = -1.0 if side in {"left", "top"} else 1.0
    cleared: list[float] = []
    for direction in (outward, -outward):
        trial = dict(sides)
        while crossings := _edge_crossings(trial, boxes, dimensions, side):
            near = [float(box[start]) for box in crossings]
            far = [float(box[start]) + float(box[size]) for box in crossings]
            position = min(near) - CROP_EDGE_MARGIN if direction < 0 else max(far) + CROP_EDGE_MARGIN
            position = min(max(position, 0.0), extent)
            beyond_reach = abs(position - origin) > CROP_EDGE_REACH
            into_subject = (position - limit) * outward < 0
            if (direction == outward and beyond_reach) or (direction != outward and into_subject):
                break
            trial[side] = position
        else:
            cleared.append(trial[side])
    return min(cleared, key=lambda position: abs(position - sides[side]), default=sides[side])


def clear_crop_edges(
    clip: dict[str, float],
    boxes: list[TextBox],
    padding: dict[str, float],
    dimensions: dict[str, float],
) -> dict[str, float]:
    """The crop with each edge moved off the text and controls it cuts (OV-3410).

    An edge moves inward only within its padding, so the crop keeps all of its subject.
    """

    origin = _crop_sides(clip)
    limits = {
        "left": origin["left"] + padding["left"],
        "top": origin["top"] + padding["top"],
        "right": origin["right"] - padding["right"],
        "bottom": origin["bottom"] - padding["bottom"],
    }
    sides = dict(origin)
    for _ in range(CROP_EDGE_ROUNDS):
        moved = False
        for side in ("top", "bottom", "left", "right"):
            position = _cleared_edge(sides, boxes, dimensions, side, origin[side], limits[side])
            moved = moved or position != sides[side]
            sides[side] = position
        if not moved:
            break
    return {
        "x": sides["left"],
        "y": sides["top"],
        "width": sides["right"] - sides["left"],
        "height": sides["bottom"] - sides["top"],
    }


class EdgeMove(TypedDict):
    """How far clear_crop_edges moved one crop edge (CSS px, outward positive) and the text it would cut."""

    by: float
    off: list[str]


def crop_edge_moves(
    clip: dict[str, float], cleared: dict[str, float], boxes: list[TextBox], dimensions: dict[str, float]
) -> dict[str, EdgeMove]:
    """Each edge that clear_crop_edges moved: by how many CSS px (outward positive) and the text
    of the boxes it would cut where it stood, in the final crop (the capture log names them)."""

    origin, sides = _crop_sides(clip), _crop_sides(cleared)
    moves: dict[str, EdgeMove] = {}
    for side in ("top", "bottom", "left", "right"):
        if sides[side] == origin[side]:
            continue
        outward = -1.0 if side in {"left", "top"} else 1.0
        crossed = _edge_crossings({**sides, side: origin[side]}, boxes, dimensions, side)
        moves[side] = {
            "by": round((sides[side] - origin[side]) * outward, 1),
            "off": sorted({box["text"][:80] for box in crossed}),
        }
    return moves


def crop_edge_cuts(
    clip: dict[str, float], boxes: list[TextBox], dimensions: dict[str, float]
) -> list[tuple[str, TextBox]]:
    """Each edge of the crop with each box of text or control it crosses."""

    sides = _crop_sides(clip)
    return [
        (side, box)
        for side in ("top", "bottom", "left", "right")
        for box in _edge_crossings(sides, boxes, dimensions, side)
    ]


def assert_crop_edges_cut_no_text(
    clip: dict[str, float], boxes: list[TextBox], dimensions: dict[str, float]
) -> None:
    """Refuse a crop whose edge still cuts text or a control a reader sees (OV-3410)."""

    cuts = crop_edge_cuts(clip, boxes, dimensions)
    if cuts:
        raise RuntimeError(
            "the crop cuts text: "
            + "; ".join(f"{side} edge cuts {box['text'][:60]!r}" for side, box in cuts)
        )


class FrameStyle(TypedDict):
    """A highlight frame: its outline offset from the element's edge (negative: inside it) and the
    width of its glow outside the element, in CSS px."""

    offset: float
    glow: float


# Below rounding error: a frame that ends where a box's clearance starts does not meet it.
_FRAME_EPSILON = 1e-6


def _frame_span(element: dict[str, float], box: TextBox) -> tuple[float, float]:
    """The distances from ``element``'s edge (outward positive, inward negative) at which a contour
    around the element meets ``box``: from where it first touches the box to where it has passed it."""

    left, top = float(element["x"]), float(element["y"])
    right, bottom = left + float(element["width"]), top + float(element["height"])
    box_left, box_top = float(box["x"]), float(box["y"])
    box_right, box_bottom = box_left + float(box["width"]), box_top + float(box["height"])
    return (
        max(box_left - right, left - box_right, box_top - bottom, top - box_bottom),
        max(left - box_left, box_right - right, top - box_top, box_bottom - bottom),
    )


def _frame_bands(style: FrameStyle, width: float = HIGHLIGHT_FRAME_WIDTH) -> list[tuple[float, float]]:
    """The bands a frame paints, as distances from its element's edge: the outline and the glow."""

    bands = [(style["offset"], style["offset"] + width)]
    if style["glow"] > 0:
        bands.append((0.0, style["glow"]))
    return bands


def highlight_frame_cuts(
    element: dict[str, float],
    style: FrameStyle,
    boxes: list[TextBox],
    margin: float,
    width: float = HIGHLIGHT_FRAME_WIDTH,
) -> list[TextBox]:
    """The boxes of text that a frame around ``element`` (its outline and its glow) comes within
    ``margin`` px of (OV-3492)."""

    cuts: list[TextBox] = []
    for box in boxes:
        near, far = _frame_span(element, box)
        if any(inner < far + margin - _FRAME_EPSILON and outer > near - margin + _FRAME_EPSILON
               for inner, outer in _frame_bands(style, width)):
            cuts.append(box)
    return cuts


def highlight_frame_style(element: dict[str, float], boxes: list[TextBox], *, glow: bool = True) -> FrameStyle | None:
    """The frame nearest its full size (an outline HIGHLIGHT_FRAME_OFFSET px out, a glow
    HIGHLIGHT_FRAME_GLOW px wide) that stays HIGHLIGHT_FRAME_CLEARANCE px clear of every box of
    text, the element's own included (``boxes`` holds no box of the element's own control shape, only the
    content box of the value it shows), or
    as clear as the element's own edge allows, or None when no outline position within its reach
    meets none (OV-3492).

    The outline moves in, down to the element's own edge inside it, where text outside the element
    sits closer; it keeps a box outside the element outside, and a box inside the element (a label
    on a figure path, the element's own caption) inside. A box across the element's edge is passed on its outer side when the
    outline still lies within its reach; the glow narrows, or goes, where text sits in it. ``glow``
    False: the element draws no glow (an SVG element), so none is planned or kept clear."""

    clearance = HIGHLIGHT_FRAME_CLEARANCE
    spans = []
    for box in boxes:
        near, far = _frame_span(element, box)
        if near >= 0:  # outside the element: kept outside the frame
            spans.append((max(near - clearance, 0.0), float("inf")))
        elif far <= 0:  # inside the element: kept inside the frame
            spans.append((near - clearance, min(far + clearance, 0.0)))
        else:  # across the element's edge
            spans.append((near - clearance, far + clearance))

    def clear(inner: float, outer: float) -> bool:
        return not any(inner < high - _FRAME_EPSILON and outer > low + _FRAME_EPSILON for low, high in spans)

    width = HIGHLIGHT_FRAME_WIDTH
    trials = [HIGHLIGHT_FRAME_OFFSET, *(low - width for low, _ in spans), *(high for _, high in spans)]
    offsets = [offset for offset in trials if -width <= offset <= HIGHLIGHT_FRAME_OFFSET and clear(offset, offset + width)]
    if not offsets:
        return None
    glows = [HIGHLIGHT_FRAME_GLOW, *(low for low, _ in spans if 0 < low < HIGHLIGHT_FRAME_GLOW)] if glow else []
    width_of_glow = max((value for value in glows if clear(0.0, value)), default=0.0)
    # Rounded to 1/100 px: box edges in fractional px must not give an offset of -1e-14 px, which
    # Chrome floors to -1 px (ATskew 7-1); the clearance absorbs the rounding.
    return {"offset": round(max(offsets), 2) + 0.0, "glow": round(width_of_glow, 2) + 0.0}


def highlight_frame_css(offset: float, glow: float, *, inset: bool = False) -> str:
    """The declarations of every frame a capture draws: its colour and width, the outline
    ``offset`` (negative: inside the element) and a glow ``glow`` px wide, outside the element or
    ``inset`` inside it."""

    return (
        f"outline: {HIGHLIGHT_FRAME_WIDTH:g}px solid {HIGHLIGHT_FRAME_COLOR} !important; "
        f"outline-offset: {offset:g}px !important; "
        f"box-shadow: {'inset ' if inset else ''}0 0 0 {glow:g}px {HIGHLIGHT_FRAME_GLOW_COLOR} !important;"
    )


HIGHLIGHT_ELEMENT_BOX_SCRIPT = """
(element) => {
  const rect = element.getBoundingClientRect();
  return { x: rect.left + window.scrollX, y: rect.top + window.scrollY, width: rect.width, height: rect.height };
}
"""

# How far around a highlighted SVG element Chromium may draw its outline: around the element's
# stroke box as it approximates it (a miter join widens it by up to 4 stroke widths), not around
# its getBoundingClientRect (OV-3492).
HIGHLIGHT_OUTLINE_SEARCH = 64.0
# Two screenshots that differ only in the colour of a frame show the pixels it paints: a frame
# of the same geometry in a clear colour leaves the page's layers and their rasterization as they
# are (adding or removing an outline re-rasterizes the shapes around an SVG element).
# A 1 px probe outline, drawn clear and then in a colour no capture uses:
_PROBE_CLEAR_CSS = "outline: 1px solid rgba(0, 0, 0, 0) !important; outline-offset: 0px !important;"
_PROBE_CSS = "outline: 1px solid rgb(255, 0, 255) !important; outline-offset: 0px !important;"
# The drawn frame, made clear:
_CLEAR_FRAME_SCRIPT = """(element) => {
  const shadow = getComputedStyle(element).boxShadow;
  element.style.setProperty('outline-color', 'rgba(0, 0, 0, 0)', 'important');
  if (shadow !== 'none') element.style.setProperty('box-shadow', shadow.replace(/rgba?\\([^)]*\\)/g, 'rgba(0, 0, 0, 0)'), 'important');
}"""
_ADD_STYLE_SCRIPT = "(element, css) => { element.style.cssText += `; ${css}`; }"
# While a frame is drawn clear and in colour, the element runs no transition or animation (a field
# or a button eases every property over 150 ms), so each screenshot shows the style it was given.
_STILL_CSS = "transition: none !important; animation: none !important;"
# Back to the saved style without a transition: first the saved values with none, then the saved
# transitions, once the browser has taken those values.
_RESTORE_STILL_SCRIPT = """([element, css, still]) => {
  element.style.cssText = `${css}; ${still}`;
  getComputedStyle(element).outlineColor;
  element.style.cssText = css;
}"""
_SAVE_STYLE_SCRIPT = "(element) => element.style.cssText"
_RESTORE_STYLE_SCRIPT = "(element, css) => { element.style.cssText = css; }"


def _grown(box: dict[str, float], by: float) -> dict[str, float]:
    x, y = max(0.0, float(box["x"]) - by), max(0.0, float(box["y"]) - by)
    return {
        "x": x,
        "y": y,
        "width": float(box["x"]) + float(box["width"]) + by - x,
        "height": float(box["y"]) + float(box["height"]) + by - y,
    }


def _painted_by(page, handle: ElementHandle, clip: dict[str, float], clear: str, drawn: str | None) -> tuple[Image.Image, float]:
    """The pixels of ``clip`` that a frame paints: those that differ between a screenshot with the
    frame clear (``clear``: declarations to add, or a script that clears the drawn frame) and one
    with it drawn (``drawn``: declarations to add, or None for the frame as it is). Returns a mask
    in device px from the top left of ``clip`` (non-zero where painted) and the device px per CSS
    px. The element's style is restored."""

    page_size = page_scroll_dimensions(page)
    clip = {
        "x": clip["x"],
        "y": clip["y"],
        "width": min(clip["width"], page_size["width"] - clip["x"]),
        "height": min(clip["height"], page_size["height"] - clip["y"]),
    }

    def shot() -> Image.Image:
        with Image.open(io.BytesIO(page.screenshot(type="png", clip=clip))) as image:
            return image.convert("RGB")

    saved = handle.evaluate(_SAVE_STYLE_SCRIPT)
    try:
        handle.evaluate(_ADD_STYLE_SCRIPT, _STILL_CSS)
        if clear.startswith("("):
            handle.evaluate(clear)
        else:
            handle.evaluate(_ADD_STYLE_SCRIPT, clear)
        unpainted = shot()
        handle.evaluate(_RESTORE_STYLE_SCRIPT, f"{saved}; {_STILL_CSS}")
        if drawn is not None:
            handle.evaluate(_ADD_STYLE_SCRIPT, drawn)
        painted = shot()
    finally:
        page.evaluate(_RESTORE_STILL_SCRIPT, [handle, saved, _STILL_CSS])
    red, green, blue = ImageChops.difference(unpainted, painted).split()
    return ImageChops.lighter(ImageChops.lighter(red, green), blue), painted.width / float(clip["width"])


def highlight_outline_box(page, handle: ElementHandle, element: dict[str, float], svg: bool) -> dict[str, float]:
    """The box the frame's outline is drawn around: an HTML element's border box (``element``); for
    an SVG element, the box Chromium draws a 1 px outline around, found by comparing screenshots
    with and without it. Its glow (a box shadow) is drawn only on HTML elements."""

    if not svg:
        return element
    clip = _grown(element, HIGHLIGHT_OUTLINE_SEARCH)
    mask, scale = _painted_by(page, handle, clip, _PROBE_CLEAR_CSS, _PROBE_CSS)
    painted = mask.getbbox()
    if painted is None:
        return element
    left, top, right, bottom = painted
    return {
        "x": clip["x"] + left / scale + 1.0,
        "y": clip["y"] + top / scale + 1.0,
        "width": (right - left) / scale - 2.0,
        "height": (bottom - top) / scale - 2.0,
    }


_FRAME_REACH = max(HIGHLIGHT_FRAME_OFFSET + HIGHLIGHT_FRAME_WIDTH, HIGHLIGHT_FRAME_GLOW) + HIGHLIGHT_FRAME_CLEARANCE


def _quoted_texts(boxes: list[TextBox]) -> str:
    return ", ".join(repr(text[:60]) for text in sorted({box["text"] for box in boxes}))


def draw_highlight_frames(page, selector: str) -> None:
    """Frame each element ``selector`` matches, clear of the text a reader sees around and inside it
    (highlight_frame_style). Run it on the page as it is captured (after the pointer
    park, the scroll resets, the hidden panes and the open-select overlays), so the frame is placed
    on what the crop shows. A frame that no position within its reach keeps off text fails the
    capture (OV-3492)."""

    for handle in page.query_selector_all(selector):
        svg = bool(handle.evaluate("(element) => element instanceof SVGElement"))
        box = highlight_outline_box(page, handle, handle.evaluate(HIGHLIGHT_ELEMENT_BOX_SCRIPT), svg)
        boxes = visible_text_boxes(page, _grown(box, _FRAME_REACH), handle)
        style = highlight_frame_style(box, boxes, glow=not svg)
        if style is None:
            # Every position the outline can take, as one band from the box's edge inside to its reach.
            positions: FrameStyle = {"offset": -HIGHLIGHT_FRAME_WIDTH, "glow": 0.0}
            reach = HIGHLIGHT_FRAME_OFFSET + 2 * HIGHLIGHT_FRAME_WIDTH
            raise RuntimeError(
                f"no highlight frame around {selector!r} within its reach stays off text: "
                + _quoted_texts(highlight_frame_cuts(box, positions, boxes, HIGHLIGHT_FRAME_CLEARANCE, reach))
            )
        handle.evaluate(
            "(element, css) => { element.style.cssText += `; ${css}`; }",
            highlight_frame_css(style["offset"], style["glow"]) + " border-radius: 4px !important;",
        )


def assert_highlight_frames_cover_no_text(page, selector: str) -> None:
    """Refuse a capture whose drawn highlight frame paints a pixel of text or a control a reader
    sees, the element's own text and the value its field shows included (OV-3492): the pixels the frame paints are those that
    differ between the frame drawn clear and drawn in its colour. Run it on the page as it is captured, after the CSS transitions."""

    for handle in page.query_selector_all(selector):
        around = _grown(handle.evaluate(HIGHLIGHT_ELEMENT_BOX_SCRIPT), HIGHLIGHT_OUTLINE_SEARCH + _FRAME_REACH)
        boxes = visible_text_boxes(page, around, handle)
        if not boxes:
            continue
        mask, scale = _painted_by(page, handle, around, _CLEAR_FRAME_SCRIPT, None)
        cuts = []
        for box in boxes:
            left = max(0, math.floor((float(box["x"]) - around["x"]) * scale))
            top = max(0, math.floor((float(box["y"]) - around["y"]) * scale))
            right = min(mask.width, math.ceil((float(box["x"]) + float(box["width"]) - around["x"]) * scale))
            bottom = min(mask.height, math.ceil((float(box["y"]) + float(box["height"]) - around["y"]) * scale))
            if right > left and bottom > top and mask.crop((left, top, right, bottom)).getbbox() is not None:
                cuts.append(box)
        if cuts:
            raise RuntimeError(f"the highlight frame around {selector!r} covers text: {_quoted_texts(cuts)}")


def assert_capture_app_state(page, capture: dict[str, Any]) -> None:
    expectations = capture.get("assertAppState")
    if not isinstance(expectations, dict) or not expectations:
        return

    mismatches = page.evaluate(
        """
        (expectations) => {
          const resolvePath = (root, path) => String(path).split('.').reduce(
            (value, segment) => value == null ? undefined : value[segment],
            root,
          );
          const equal = (actual, expected) => {
            if (Object.is(actual, expected)) return true;
            try {
              return JSON.stringify(actual) === JSON.stringify(expected);
            } catch {
              return false;
            }
          };
          const describe = (value) => {
            if (value === undefined) return '<undefined>';
            if (value instanceof File) return { name: value.name, type: value.type };
            try {
              return JSON.parse(JSON.stringify(value));
            } catch {
              return String(value);
            }
          };
          const app = window.__GBDRAW_APP__;
          if (!app) return [{ path: '<app>', expected: 'mounted app', actual: '<undefined>' }];
          return Object.entries(expectations).flatMap(([path, expected]) => {
            const actual = resolvePath(app, path);
            return equal(actual, expected) ? [] : [{ path, expected, actual: describe(actual) }];
          });
        }
        """,
        expectations,
    )
    if mismatches:
        details = "; ".join(
            f"{item['path']} expected {item['expected']!r}, got {item['actual']!r}"
            for item in mismatches
        )
        raise RuntimeError(f"capture assertAppState failed: {details}")


def assert_capture_visible_contract(
    page,
    capture: dict[str, Any],
    clip: dict[str, float],
) -> None:
    controls = capture.get("visibleControls")
    text_entries = capture.get("visibleText")
    controls = controls if isinstance(controls, list) else []
    text_entries = text_entries if isinstance(text_entries, list) else []
    if not controls and not text_entries:
        return

    failures = page.evaluate(
        """
        ({ controls, textEntries, clip }) => {
          const normalize = (value) => String(value ?? '').replace(/\\s+/g, ' ').trim();
          const key = (value) => normalize(value).toLowerCase();
          const pageRect = (element) => {
            const rect = element.getBoundingClientRect();
            return {
              x: rect.left + window.scrollX,
              y: rect.top + window.scrollY,
              width: rect.width,
              height: rect.height,
            };
          };
          const isVisible = (element) => {
            if (!element) return false;
            const style = getComputedStyle(element);
            const rect = element.getBoundingClientRect();
            return style.display !== 'none'
              && style.visibility !== 'hidden'
              && Number(style.opacity || 1) > 0
              && rect.width > 0
              && rect.height > 0;
          };
          const insideClip = (rect) => {
            const tolerance = 1;
            return rect.x >= clip.x - tolerance
              && rect.y >= clip.y - tolerance
              && rect.x + rect.width <= clip.x + clip.width + tolerance
              && rect.y + rect.height <= clip.y + clip.height + tolerance;
          };
          const controlForLabel = (label) => {
            if (label.control) return label.control;
            const nested = label.querySelector('input, select, textarea, button');
            if (nested) return nested;
            let container = label.parentElement;
            for (let depth = 0; container && depth < 3; depth += 1) {
              const nearby = container.querySelector('input, select, textarea, button');
              if (nearby) return nearby;
              container = container.parentElement;
            }
            return null;
          };
          const resolveControl = (spec) => {
            if (spec.selector) {
              const matches = Array.from(document.querySelectorAll(spec.selector));
              const element = matches[Number(spec.nth || 0)] || null;
              return { element, label: null };
            }
            const labelKey = key(spec.label);
            const labels = Array.from(document.querySelectorAll('label')).filter(
              (label) => key(label.textContent).includes(labelKey),
            );
            const label = labels[Number(spec.nth || 0)] || null;
            return { element: label ? controlForLabel(label) : null, label };
          };

          const failures = [];
          controls.forEach((spec, index) => {
            const name = spec.label || spec.selector || `control ${index + 1}`;
            let resolved;
            try {
              resolved = resolveControl(spec);
            } catch (error) {
              failures.push(`${name}: invalid selector (${error.message})`);
              return;
            }
            const { element, label } = resolved;
            if (!isVisible(element)) {
              failures.push(`${name}: control is not visible`);
              return;
            }
            const controlRect = pageRect(element);
            if (!insideClip(controlRect)) {
              failures.push(
                `${name}: control is outside the final crop `
                + `(control=${JSON.stringify(controlRect)}, clip=${JSON.stringify(clip)})`,
              );
            }
            const labelRect = label ? pageRect(label) : null;
            if (label && (!isVisible(label) || !insideClip(labelRect))) {
              failures.push(
                `${name}: label is outside the final crop `
                + `(label=${JSON.stringify(labelRect)}, clip=${JSON.stringify(clip)})`,
              );
            }
            if (Object.prototype.hasOwnProperty.call(spec, 'checked')) {
              const actual = Boolean(element.checked);
              if (actual !== Boolean(spec.checked)) {
                failures.push(`${name}: expected checked=${Boolean(spec.checked)}, got ${actual}`);
              }
            }
            if (Object.prototype.hasOwnProperty.call(spec, 'pressed')) {
              const actual = element.getAttribute('aria-pressed') === 'true';
              if (actual !== Boolean(spec.pressed)) {
                failures.push(`${name}: expected pressed=${Boolean(spec.pressed)}, got ${actual}`);
              }
            }
            if (Object.prototype.hasOwnProperty.call(spec, 'value')) {
              const expected = normalize(spec.value);
              const selected = element.tagName === 'SELECT'
                ? normalize(element.options[element.selectedIndex]?.textContent)
                : '';
              const actual = normalize(element.value);
              if (expected !== actual && expected !== selected) {
                failures.push(`${name}: expected value "${expected}", got "${selected || actual}"`);
              }
            }
          });

          textEntries.forEach((expectedText, index) => {
            const expectedKey = key(expectedText);
            const candidates = Array.from(document.querySelectorAll('body *'))
              .filter((element) => !['SCRIPT', 'STYLE'].includes(element.tagName))
              .filter((element) => isVisible(element) && key(element.textContent).includes(expectedKey))
              .map((element) => ({ element, rect: pageRect(element) }))
              .sort((left, right) => (
                left.rect.width * left.rect.height - right.rect.width * right.rect.height
              ));
            const match = candidates.find(({ rect }) => insideClip(rect));
            if (!match) {
              failures.push(`visibleText ${index + 1}: "${normalize(expectedText)}" is not visible in the final crop`);
            }
          });
          return failures;
        }
        """,
        {
            "controls": controls,
            "textEntries": text_entries,
            "clip": clip,
        },
    )
    if failures:
        raise RuntimeError("capture visibility contract failed: " + "; ".join(failures))


def prepare_capture_page(page, base_url: str, sample: dict[str, Any], capture: dict[str, Any]) -> None:
    source = as_text(capture.get("source")).lower() or "gallery"
    viewport = capture.get("viewport") if isinstance(capture.get("viewport"), dict) else DEFAULT_VIEWPORT
    page.set_viewport_size(
        {
            "width": int(viewport.get("width") or DEFAULT_VIEWPORT["width"]),
            "height": int(viewport.get("height") or DEFAULT_VIEWPORT["height"]),
        }
    )

    if source == "gallery":
        page.goto(f"{base_url}/gallery/#{sample['id']}", wait_until="domcontentloaded")
        tab = as_text(capture.get("tab")) or "tutorial"
        if tab != "preview":
            page.locator(f"[data-gallery-tab='{tab}']").click()
        panel_selector = {
            "preview": "#preview-panel",
            "tutorial": "#tutorial-content",
            "command": "#command-panel",
            "files": "#files-content",
        }.get(tab, "#tutorial-content")
        page.locator(panel_selector).wait_for(state="visible")
    elif source in {"webapp", "web-app", "app"}:
        page.goto(f"{base_url}/", wait_until="domcontentloaded")
        wait_for_app_shell(page)
        # capture.freshPage: a reader's new page, before any example is loaded.
        session_ref = "" if capture.get("freshPage") is True else (
            as_text(capture.get("session")) or as_text(sample.get("session"))
        )
        session_path = resolve_gallery_reference(session_ref)
        if session_path is not None and session_path.exists():
            load_web_app_session(page, session_path)
            dismiss_session_notice(page)
            # A loaded Session folds its Linear File cards. A tutorial reader uploaded
            # the files, so every capture shows the cards open.
            page.evaluate(
                "() => document.querySelectorAll('[data-linear-source-fold]')"
                ".forEach((fold) => { fold.open = true; })"
            )
            # The uploader offers the files of a loaded Session (**Files in this Session**). A
            # tutorial reader uploaded the files and has no Session, so no capture shows it
            # (OV-3330); `add_step_text_validation` keeps every step from naming it.
            page.add_style_tag(content=f"{SESSION_FILE_CHOOSER_SELECTOR} {{ display: none !important; }}")
        replay_through = capture.get("replayThrough")
        if isinstance(replay_through, int):
            replay_tutorial_steps(page, sample, replay_through)
        # After the replay, so the read meets the settings of the reader's step, as an upload
        # does: a read that finds one record opens the single-record panel only when those
        # settings allow it.
        if session_path is not None and session_path.exists():
            inspect_shown_records(page)
    else:
        raise RuntimeError(f"Unsupported capture source: {source}")

    open_capture_details(page, capture.get("openDetails"))

    for action in as_array(capture.get("actions")):
        if isinstance(action, dict):
            apply_capture_action(page, action)

    wait_selector = as_text(capture.get("waitForSelector"))
    if wait_selector:
        page.locator(wait_selector).wait_for(state="visible")



def hide_settings_pane_under_popup(page) -> None:
    """A popup a capture opens may sit over the Settings pane: hide the pane, so the crop
    shows no settings around the popup. The pane keeps its place, so nothing moves."""

    page.evaluate(
        """
        () => {
          const popup = document.querySelector('.feature-popup');
          if (!popup?.checkVisibility()) return;
          document.querySelector('.settings-pane')?.style.setProperty('visibility', 'hidden', 'important');
        }
        """
    )


def hide_page_around_settings_pane(page) -> None:
    """An open-select crop shows a settings control and its menu, which may reach past the
    Settings pane: hide the rest of the page (header, Result, the footer's text), so the crop's
    edge cuts no figure, button or text there (OV-3395). Everything keeps its place, and the
    footer itself stays: the pane's Generate bar is anchored to it and would hide with it."""

    page.evaluate(
        """
        () => {
          document.querySelectorAll('.app-header, #result-pane, [data-app-footer] > *').forEach((element) => {
            element.style.setProperty('visibility', 'hidden', 'important');
          });
        }
        """
    )


def wait_for_css_transitions(page) -> None:
    """Let running CSS transitions end, so a crop never shows a control mid-change.

    The caret of a card or section turns over 200 ms after it opens. A transition inside a
    closed <details> never settles, because the browser skips that content, and the reader
    never sees it: wait only for the transitions the page renders (OV-2450).
    """

    unfinished = page.evaluate(
        """
        () => {
          const running = new Set(
            document.getAnimations().filter(
              (animation) => animation instanceof CSSTransition && animation.effect?.target?.checkVisibility(),
            ),
          );
          const settled = [...running].map((animation) => animation.finished
            .catch(() => undefined)
            .then(() => running.delete(animation)));
          const describe = () => [...running].map((animation) => {
            const target = animation.effect.target;
            return `${target.tagName.toLowerCase()}${target.id ? `#${target.id}` : ''} ${animation.transitionProperty}`;
          });
          return Promise.race([
            Promise.all(settled).then(() => []),
            new Promise((resolve) => setTimeout(() => resolve(describe()), 5000)),
          ]);
        }
        """
    )
    if unfinished:
        raise RuntimeError(f"CSS transitions did not finish in 5 s: {', '.join(unfinished)}")


def capture_operation(page, sample: dict[str, Any], operation: dict[str, Any], base_url: str, quality: int) -> Path:
    capture = operation.get("capture")
    if not isinstance(capture, dict):
        raise RuntimeError("operation is missing capture metadata.")

    target = output_path_for_operation(operation)
    if target is None:
        raise RuntimeError("operation media must point at a local output path.")

    prepare_capture_page(page, base_url, sample, capture)
    park_pointer(page)
    target_locator = resolve_capture_target(page, capture)
    if target_locator is not None:
        target_locator.scroll_into_view_if_needed()
        page.wait_for_timeout(100)
    page.evaluate(
        """
        () => {
          document.querySelectorAll('.app-main, .settings-pane, .settings-scroll').forEach((element) => {
            element.scrollLeft = 0;
          });
          document.documentElement.scrollLeft = 0;
          document.body.scrollLeft = 0;
          window.scrollTo({ left: 0, top: window.scrollY });
        }
        """
    )
    hide_settings_pane_under_popup(page)
    if as_text(capture.get("crop")) == "openSelect":
        hide_page_around_settings_pane(page)
    apply_open_select_overlays(page, capture.get("openSelect"))
    highlight_selector = as_text(capture.get("highlightSelector"))
    if highlight_selector:
        draw_highlight_frames(page, highlight_selector)
    wait_for_css_transitions(page)
    if highlight_selector:
        assert_highlight_frames_cover_no_text(page, highlight_selector)

    clip_selectors = [as_text(selector) for selector in as_array(capture.get("clipSelectors")) if as_text(selector)]

    no_padding = normalize_crop_padding(0, default=0)
    if as_text(capture.get("crop")) == "openSelect":
        padding = normalize_crop_padding(capture.get("cropPadding"), default=14)
        clip = capture_open_select_clip(page, padding, clip_selectors)
    elif clip_selectors:
        padding = normalize_crop_padding(capture.get("cropPadding"), default=24)
        clip = capture_multi_selector_clip(page, clip_selectors, padding)
    elif target_locator is not None:
        padding = (
            no_padding
            if capture.get("cropPadding") is None
            else normalize_crop_padding(capture.get("cropPadding"), default=12)
        )
        clip = capture_selector_clip(page, target_locator, padding)
    else:
        box = capture.get("boundingBox")
        if not isinstance(box, dict):
            raise RuntimeError("capture metadata needs selector, crop=openSelect, or boundingBox.")
        padding = no_padding
        clip = {
            "x": float(box["x"]),
            "y": float(box["y"]),
            "width": float(box["width"]),
            "height": float(box["height"]),
        }

    dimensions = page_scroll_dimensions(page)
    around = CROP_EDGE_REACH + CROP_EDGE_MARGIN
    boxes = visible_text_boxes(
        page,
        {
            "x": clip["x"] - around,
            "y": clip["y"] - around,
            "width": clip["width"] + 2 * around,
            "height": clip["height"] + 2 * around,
        },
    )
    cleared = clear_crop_edges(clip, boxes, padding, dimensions)
    moves = crop_edge_moves(clip, cleared, boxes, dimensions)
    if moves:
        print(f"crop edges moved off text (CSS px, outward positive): {json.dumps(moves)}", flush=True)
    clip = cleared
    assert_crop_edges_cut_no_text(clip, boxes, dimensions)
    if as_text(capture.get("source")).lower() in {"webapp", "web-app", "app"}:
        assert_capture_records_read(page)
    assert_capture_app_state(page, capture)
    assert_capture_visible_contract(page, capture, clip)
    png_bytes = page.screenshot(type="png", clip=clip)

    target.parent.mkdir(parents=True, exist_ok=True)
    with Image.open(io.BytesIO(png_bytes)) as image:
        image.convert("RGB").save(target, "WEBP", quality=quality, method=6)
    return target


def capture_examples(
    samples: list[dict[str, Any]],
    *,
    operation_filter: str,
    quality: int,
    device_scale_factor: float,
    keep_going: bool,
) -> tuple[int, list[str]]:
    capturable: list[tuple[dict[str, Any], TutorialContext, dict[str, Any]]] = []
    for sample in samples:
        tutorial_path = resolve_gallery_reference(as_text(sample.get("tutorial")))
        if tutorial_path is None:
            continue
        tutorial = read_json(tutorial_path)
        for context, operation in iter_operation_contexts(sample, tutorial):
            target = output_path_for_operation(operation)
            target_name = target.name if target is not None else ""
            if operation_filter and operation_filter not in target_name and operation_filter not in context.label:
                continue
            if isinstance(operation.get("capture"), dict):
                capturable.append((sample, context, operation))

    if not capturable:
        print("No operations with capture metadata were found.")
        return 0, []

    sync_playwright = import_playwright()
    with serve_web_root() as base_url, sync_playwright() as playwright:
        try:
            browser = playwright.chromium.launch(headless=True)
        except Exception as exc:  # pragma: no cover - environment dependent
            raise RuntimeError(
                "Could not launch Playwright Chromium. In sandboxed agent environments, rerun this "
                f"command with sandbox escalation. Underlying error: {exc}"
            ) from exc
        page = browser.new_page(device_scale_factor=device_scale_factor)
        captured_count = 0
        failures: list[str] = []
        try:
            for sample, context, operation in capturable:
                try:
                    target = capture_operation(page, sample, operation, base_url, quality)
                    rel = target.relative_to(REPO_ROOT)
                    with Image.open(target) as image:
                        print(f"captured {context.label}: {rel} {image.size[0]}x{image.size[1]}", flush=True)
                    captured_count += 1
                except Exception as exc:
                    message = f"failed {context.label}: {exc}"
                    if not keep_going:
                        raise
                    print(message, file=sys.stderr, flush=True)
                    failures.append(message)
        finally:
            browser.close()
    if failures:
        print(f"capture completed with {len(failures)} failure(s).", file=sys.stderr)
    return captured_count, failures


def print_validation(result: ValidationResult, *, strict: bool) -> int:
    for warning in result.warnings:
        print(f"warning: {warning}", file=sys.stderr)
    for error in result.errors:
        print(f"error: {error}", file=sys.stderr)

    print(
        "checked "
        f"{result.media_count} media entries, "
        f"{result.operation_count} operations, "
        f"{result.operation_media_count} operation media entries"
    )
    if result.errors or (strict and result.warnings):
        return 1
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Capture and validate operation-level Gallery tutorial screenshots.",
    )
    target = parser.add_mutually_exclusive_group(required=True)
    target.add_argument("--example", help="Capture or check one ready Gallery example ID.")
    target.add_argument("--all", action="store_true", help="Capture or check all ready Gallery examples.")
    parser.add_argument("--check", action="store_true", help="Validate tutorial media instead of capturing screenshots.")
    parser.add_argument("--strict", action="store_true", help="Treat validation warnings as failures.")
    parser.add_argument("--operation", default="", help="Capture only operations whose label or output filename contains this text.")
    parser.add_argument("--quality", type=int, default=DEFAULT_WEBP_QUALITY, help="WebP quality for captured screenshots.")
    parser.add_argument(
        "--device-scale-factor",
        type=float,
        default=DEFAULT_DEVICE_SCALE_FACTOR,
        help="Browser device scale factor for captured screenshots.",
    )
    parser.add_argument("--keep-going", action="store_true", help="Continue capturing after an operation fails.")
    parser.add_argument("--max-width", type=int, default=DEFAULT_MAX_IMAGE_WIDTH, help="Maximum accepted image width.")
    parser.add_argument("--max-height", type=int, default=DEFAULT_MAX_IMAGE_HEIGHT, help="Maximum accepted image height.")
    parser.add_argument(
        "--max-file-size-kb",
        type=int,
        default=DEFAULT_MAX_FILE_SIZE_KB,
        help="Maximum accepted media file size in KiB.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    samples = load_ready_examples(None if args.all else args.example)

    if args.check:
        result = validate_tutorial_media(
            samples,
            max_width=args.max_width,
            max_height=args.max_height,
            max_file_size_kb=args.max_file_size_kb,
        )
        return print_validation(result, strict=args.strict)

    _, failures = capture_examples(
        samples,
        operation_filter=args.operation,
        quality=args.quality,
        device_scale_factor=args.device_scale_factor,
        keep_going=args.keep_going,
    )
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
