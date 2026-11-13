#!/usr/bin/env python3
from __future__ import annotations

import argparse
import contextlib
import http.server
import io
import json
import re
import socketserver
import sys
import tempfile
import threading
from dataclasses import dataclass
from pathlib import Path
from typing import Any

try:
    from PIL import Image
except ImportError as exc:  # pragma: no cover - dependency guard
    raise SystemExit(
        "Pillow is required. Install the dev dependencies with `pip install -e .[dev]`."
    ) from exc


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from gbdraw.session_container import read_session_container  # noqa: E402

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
        "drag",
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
                session_upload_name(action) if "sessionFile" in action else resolve_upload_file(action)
            except ValueError as error:
                result.errors.append(f"{context.label}: replay action {index} {error}")
        elif action_type == "drag" and (problem := drag_action_problem(action)):
            result.errors.append(f"{context.label}: replay action {index} {problem}")


def drag_action_problem(action: dict[str, Any]) -> str | None:
    """What a `drag` lacks: a record ID, two of its bases (1-based, in its file), and a window index if any."""

    def whole(value: Any, low: int) -> bool:
        return isinstance(value, int) and not isinstance(value, bool) and value >= low

    if not as_text(action.get("record")):
        return "drag needs the record ID it is made on."
    if not (whole(action.get("from"), 1) and whole(action.get("to"), 1)):
        return "drag needs `from` and `to`, bases of the record counted from 1."
    if "window" in action and not whole(action.get("window"), 0):
        return "drag `window` must be the index (from 0) of one window of the record."
    return None


def add_created_drawings_validation(result: ValidationResult, example_id: str, tutorial: dict[str, Any]) -> None:
    """`createdDrawings`: each drawing a manual step's replay creates, by its drawing ID and that step.

    The Gallery Session holds them; the replay starts without them (`replay_tutorial_steps`).
    """

    steps = as_array(tutorial.get("manualSteps"))
    seen: set[str] = set()
    for index, entry in enumerate(as_array(tutorial.get("createdDrawings")), start=1):
        label = f"{example_id}: createdDrawings {index}"
        drawing_id = as_text(entry.get("id")) if isinstance(entry, dict) else ""
        step = entry.get("step") if isinstance(entry, dict) else None
        if not drawing_id:
            result.errors.append(f"{label} needs the drawing ID of the drawing a step creates.")
        elif drawing_id in seen:
            result.errors.append(f"{label} names {drawing_id} twice.")
        elif not (isinstance(step, int) and not isinstance(step, bool) and 1 <= step <= len(steps)):
            result.errors.append(f"{label} must name the manual step (1 to {len(steps)}) that creates {drawing_id}.")
        elif not (isinstance(steps[step - 1], dict) and as_array(steps[step - 1].get("replay"))):
            result.errors.append(f"{label} names manual step {step}, which has no replay to create {drawing_id}.")
        seen.add(drawing_id)


def session_upload_name(action: dict[str, Any]) -> str:
    """The file an `upload` action chooses from the capture's Gallery Session, by the name it was uploaded as.

    A reader downloads that file from Files (an input the Session holds, such as a genome too large for
    gbdraw/web); the replay writes the Session's copy under that name and chooses it.
    """

    name = as_text(action.get("sessionFile"))
    if not as_text(action.get("selector")):
        raise ValueError("upload needs the selector of the button that opens the file chooser.")
    if "path" in action:
        raise ValueError("upload names a file by `path` or by `sessionFile`, not both.")
    if not name or Path(name).name != name:
        raise ValueError(f"upload sessionFile must be the file name an input of the Session was uploaded as: {name or '(none)'}")
    return name


def session_input_file_bytes(session_path: Path, name: str) -> bytes:
    """The bytes of the input file of a Session container that was uploaded as `name`."""

    container = read_session_container(session_path)
    web_files = container.document.get("webFiles")
    original_names = web_files.get("resourceOriginalNames") if isinstance(web_files, dict) else None
    resource_ids = sorted(
        resource_id for resource_id, original in (original_names or {}).items() if original == name
    )
    if not resource_ids:
        raise RuntimeError(f"upload: {session_path.name} holds no input file named {name}.")
    with container.resource(resource_ids[0]).open_bytes() as stream:
        return stream.read()


@contextlib.contextmanager
def upload_file(action: dict[str, Any], session_path: Path | None):
    """The local file an `upload` chooses: under gbdraw/web (`path`), or the Session's copy (`sessionFile`)."""

    if "sessionFile" not in action:
        yield resolve_upload_file(action)
        return
    name = session_upload_name(action)
    if session_path is None:
        raise RuntimeError(f"upload: {name} comes from the capture's Session, and this capture loads none.")
    with tempfile.TemporaryDirectory(prefix="gbdraw-capture-upload-") as folder:
        path = Path(folder) / name
        path.write_bytes(session_input_file_bytes(session_path, name))
        yield path


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

        add_created_drawings_validation(result, as_text(sample.get("id")), tutorial)
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


def assert_capture_records_read(page) -> None:
    """A webapp crop shows the records read, as the reader who uploaded the files sees them."""

    if not page.evaluate(SHOWN_RECORDS_READ_SCRIPT):
        raise RuntimeError(
            "capture shows records not inspected: the reader who uploaded the files sees them read."
        )


def apply_capture_action(page, action: dict[str, Any], session_path: Path | None = None) -> None:
    action_type = as_text(action.get("type")) or as_text(action.get("action"))
    selector = as_text(action.get("selector"))
    value = action.get("value")

    if action_type == "click":
        # A click that starts long work (Create new drawing) names how long the page may keep it waiting.
        click_options: dict[str, Any] = {"force": True} if action.get("force") is True else {}
        if action.get("timeout") is not None:
            click_options["timeout"] = int(action["timeout"])
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
        with upload_file(action, session_path) as path:
            with page.expect_file_chooser() as chooser:
                page.locator(selector).first.click()
            chooser.value.set_files(str(path))
    elif action_type == "drag":
        drag_between_bases(page, action)
    else:
        raise RuntimeError(f"Unsupported capture action: {action_type}")
    # The next action measures what this one moved: a zoom or pan moves the preview in a
    # 0.2 s transition, and a box read before it ends is the old one (OV-2280).
    wait_for_css_transitions(page)


# The records of the shown Result whose axis Python wrote, read by the preview's own reader.
SHOWN_AXIS_RECORDS_SCRIPT = """
async () => {
  const { state } = await import('./js/state.js');
  const { shownAxisRecords } = await import('./js/services/region-axis.js');
  return shownAxisRecords(state);
}
"""

# Points of a record's axis frame on the screen. `sourceBaseAt`, the preview's map from a drawn
# base to the file's, checks each display index `drag_frame_points` found.
AXIS_POINTS_TO_SCREEN_SCRIPT = """
async ({ recordIndex, points }) => {
  const { state } = await import('./js/state.js');
  const { shownAxisRecords, regionAxisFrame, sourceBaseAt } = await import('./js/services/region-axis.js');
  const record = shownAxisRecords(state).find((entry) => entry.recordIndex === recordIndex);
  const frame = regionAxisFrame(state.svgContainer.value?.querySelector('svg') ?? document, recordIndex);
  if (!record || !frame) return { error: `the shown Result draws no record ${recordIndex}` };
  const wrong = points.find((point) => sourceBaseAt(record.axis.source, point.index) !== point.base);
  if (wrong) return { error: `display index ${wrong.index} shows base ${sourceBaseAt(record.axis.source, wrong.index)}, not ${wrong.base}` };
  const matrix = frame.getScreenCTM();
  if (!matrix) return { error: 'the record is not on the screen' };
  const screen = points.map((point) => new DOMPoint(point.x, point.y).matrixTransform(matrix));
  if (screen.some(({ x, y }) => x < 0 || y < 0 || x > window.innerWidth || y > window.innerHeight)) {
    return { error: 'a drag end is outside the viewport' };
  }
  return { points: screen.map(({ x, y }) => ({ x, y })) };
}
"""


def drag_frame_points(records: list[dict[str, Any]], action: dict[str, Any]) -> dict[str, Any]:
    """Where a `drag` presses and releases, in the record's axis frame.

    `records` are the shown Result's axis records (`shownAxisRecords` in services/region-axis.js).
    The drag names a record by its ID (and a window by its index when the record is drawn in
    several) and two of its bases, 1-based in the input file. Each end is the middle of the drawn
    base on its Linear row, half way across the feature tracks, so the preview's
    `linearRowPosition` reads it back as that base; a drag stays on one row.
    """

    record_id = as_text(action.get("record"))
    matches = [record for record in records if record.get("recordId") == record_id]
    if "window" in action:
        matches = [record for record in matches if record.get("windowIndex") == action.get("window")]
    if not matches:
        raise ValueError(f"drag: the shown Result draws no record {record_id} (window {action.get('window', '-')})")
    if len(matches) > 1:
        raise ValueError(f"drag: {record_id} is drawn in {len(matches)} windows; name one with `window`.")
    record = matches[0]
    axis = record["axis"]
    if axis.get("kind") != "linear":
        raise ValueError("drag: only a Linear record's axis is supported.")
    source = axis["source"]
    points: list[dict[str, Any]] = []
    rows: set[int] = set()
    for base in (int(action["from"]), int(action["to"])):
        # The inverse of `sourceBaseAt` (mod(anchorBase - 1 + step * index, length) + 1);
        # the page checks it with `sourceBaseAt` itself.
        index = (source["step"] * (base - source["anchorBase"])) % source["sourceLength"]
        row = next(
            (number for number, row in enumerate(axis["rows"]) if row["displayStart"] <= index < row["displayEnd"]),
            None,
        )
        if row is None:
            raise ValueError(f"drag: base {base} of {record_id} is not drawn.")
        drawn = axis["rows"][row]
        fraction = (index + 0.5 - drawn["displayStart"]) / (drawn["displayEnd"] - drawn["displayStart"])
        points.append(
            {
                "base": base,
                "index": index,
                "x": drawn["xStartPx"] + fraction * (drawn["xEndPx"] - drawn["xStartPx"]),
                "y": (drawn["bandTopPx"] + drawn["bandBottomPx"]) / 2,
            }
        )
        rows.add(row)
    if len(rows) > 1:
        raise ValueError(f"drag: bases {action['from']} and {action['to']} of {record_id} are not on one row.")
    return {"recordIndex": record["recordIndex"], "points": points}


def drag_between_bases(page, action: dict[str, Any]) -> None:
    """A reader's drag with Select region on, from one base of a record to another (real mouse events)."""

    mapped = drag_frame_points(page.evaluate(SHOWN_AXIS_RECORDS_SCRIPT), action)
    screen = page.evaluate(AXIS_POINTS_TO_SCREEN_SCRIPT, mapped)
    if not isinstance(screen, dict) or "points" not in screen:
        raise RuntimeError(f"drag: {screen.get('error') if isinstance(screen, dict) else screen}")
    start, end = screen["points"]
    page.mouse.move(start["x"], start["y"])
    page.mouse.down()
    page.mouse.move(end["x"], end["y"], steps=8)
    page.mouse.up()


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


def replay_tutorial_steps(page, sample: dict[str, Any], through: int, session_path: Path | None = None) -> None:
    """Start as the reader did, then perform the reader's actions of manual steps 1..through.

    The start: the Session without the drawings its steps create (`createdDrawings`, by drawing
    ID; a replay through the creating step creates each again), and Reset Settings.
    """

    tutorial_path = resolve_gallery_reference(as_text(sample.get("tutorial")))
    if tutorial_path is None:
        raise RuntimeError(f"Tutorial path for {sample.get('id')} is not local.")
    tutorial = read_json(tutorial_path)
    steps = as_array(tutorial.get("manualSteps"))[:through]
    for entry in as_array(tutorial.get("createdDrawings")):
        remove_drawing(page, as_text(entry.get("id")))
    apply_capture_action(page, {"type": "resetSettings"})
    for step in steps:
        for action in as_array(step.get("replay") if isinstance(step, dict) else None):
            apply_capture_action(page, action, session_path)


def remove_drawing(page, drawing_id: str) -> None:
    """Delete a drawing by its ID through its tab's Drawing menu, confirming when the app asks."""

    tab = page.locator(f'.drawing-bar [role="tab"][data-drawing-id="{drawing_id}"]')
    page.wait_for_function(
        "() => [...document.querySelectorAll('.drawing-bar [role=\"tab\"]')].every((tab) => !tab.disabled)",
        timeout=DEFAULT_GENERATE_TIMEOUT_MS,
    )
    tab.click()
    menu = page.locator(".drawing-bar-menu")
    menu.locator("summary").click()
    menu.get_by_role("button", name="Delete drawing", exact=True).click()
    # The app asks only for a drawing with a Result or edits.
    confirm = page.get_by_role("dialog").filter(has=page.get_by_role("button", name="Delete drawing", exact=True))
    page.wait_for_function(
        """(id) => !document.querySelector(`.drawing-bar [role="tab"][data-drawing-id="${id}"]`)
          || [...document.querySelectorAll('[role="dialog"]')].some((dialog) => dialog.textContent.includes('Keep drawing'))""",
        arg=drawing_id,
        timeout=DEFAULT_GENERATE_TIMEOUT_MS,
    )
    if confirm.count():
        confirm.get_by_role("button", name="Delete drawing", exact=True).click()
    tab.wait_for(state="detached")
    wait_for_css_transitions(page)


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
              .gbdraw-capture-select-control {
                outline: 3px solid #f59e0b !important;
                outline-offset: -3px !important;
                box-shadow: inset 0 0 0 4px rgba(245, 158, 11, 0.18) !important;
              }
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
                outline: 3px solid #f59e0b;
                outline-offset: -2px;
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
        """,
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
            replay_tutorial_steps(page, sample, replay_through, session_path)
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

    highlight_selector = as_text(capture.get("highlightSelector"))
    if highlight_selector:
        page.add_style_tag(
            content=(
                f"{highlight_selector} {{ "
                "outline: 3px solid #f59e0b !important; "
                "outline-offset: 3px !important; "
                "box-shadow: 0 0 0 4px rgba(245, 158, 11, 0.18) !important; "
                "border-radius: 4px !important; "
                "}}"
            )
        )
        page.wait_for_timeout(100)


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
    never sees it: wait only for the transitions the page renders (OV-2450). A busy page (a
    dialog of many rows opening while a large File is read) may need more than a second or two.
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
            new Promise((resolve) => setTimeout(() => resolve(describe()), 15000)),
          ]);
        }
        """
    )
    if unfinished:
        raise RuntimeError(f"CSS transitions did not finish in 15 s: {', '.join(unfinished)}")


def capture_operation(page, sample: dict[str, Any], operation: dict[str, Any], base_url: str, quality: int) -> Path:
    capture = operation.get("capture")
    if not isinstance(capture, dict):
        raise RuntimeError("operation is missing capture metadata.")

    target = output_path_for_operation(operation)
    if target is None:
        raise RuntimeError("operation media must point at a local output path.")

    prepare_capture_page(page, base_url, sample, capture)
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
    wait_for_css_transitions(page)

    clip_selectors = [as_text(selector) for selector in as_array(capture.get("clipSelectors")) if as_text(selector)]

    if as_text(capture.get("crop")) == "openSelect":
        clip = capture_open_select_clip(
            page,
            normalize_crop_padding(capture.get("cropPadding"), default=14),
            clip_selectors,
        )
    elif clip_selectors:
        clip = capture_multi_selector_clip(
            page,
            clip_selectors,
            normalize_crop_padding(capture.get("cropPadding"), default=24),
        )
    elif target_locator is not None and capture.get("cropPadding") is not None:
        clip = capture_selector_clip(
            page,
            target_locator,
            normalize_crop_padding(capture.get("cropPadding"), default=12),
        )
    elif target_locator is not None:
        clip = capture_selector_clip(
            page,
            target_locator,
            normalize_crop_padding(0, default=0),
        )
    else:
        box = capture.get("boundingBox")
        if not isinstance(box, dict):
            raise RuntimeError("capture metadata needs selector, crop=openSelect, or boundingBox.")
        clip = {
            "x": float(box["x"]),
            "y": float(box["y"]),
            "width": float(box["width"]),
            "height": float(box["height"]),
        }

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
