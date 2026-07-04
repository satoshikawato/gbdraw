#!/usr/bin/env python3
from __future__ import annotations

import argparse
import contextlib
import http.server
import io
import json
import socketserver
import sys
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
WEB_ROOT = REPO_ROOT / "gbdraw" / "web"
GALLERY_ROOT = WEB_ROOT / "gallery"
GALLERY_TUTORIAL_ROOT = GALLERY_ROOT / "tutorials"
GALLERY_MEDIA_ROOT = GALLERY_ROOT / "media"

READY_SECTIONS = (
    "quickReproduce",
    "manualSteps",
    "postGenerationEdits",
    "colorRules",
)
OPERATION_REQUIRED_SECTIONS = (
    "quickReproduce",
    "manualSteps",
    "postGenerationEdits",
)

DEFAULT_MAX_IMAGE_WIDTH = 1400
DEFAULT_MAX_IMAGE_HEIGHT = 1400
DEFAULT_MAX_FILE_SIZE_KB = 650
DEFAULT_WEBP_QUALITY = 86
DEFAULT_VIEWPORT = {"width": 1280, "height": 900}


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
    full_screen_like = width >= 1800 or (width >= 1400 and height >= 900) or (width >= 1200 and height >= 1000)
    if too_large or full_screen_like:
        message = (
            f"{context.label}: {source} dimensions are {width}x{height}; "
            f"limit is {max_width}x{max_height}."
        )
        if operation_media:
            result.errors.append(message)
        else:
            result.warnings.append(message)


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
            if context.section in OPERATION_REQUIRED_SECTIONS and not as_array(step.get("operations")):
                result.errors.append(f"{context.label}: step is missing operations.")
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
    handler = lambda *args, **kwargs: QuietGalleryRequestHandler(
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


def wait_for_web_app_ready(page) -> None:
    page.wait_for_function(
        """
        () => {
          const app = window.__GBDRAW_APP__;
          if (!app) return false;
          const status = String(app.loadingStatus || '');
          return app.pyodideReady === true || status.startsWith('Startup Error:');
        }
        """,
        timeout=120000,
    )


def load_web_app_session(page, session_path: Path) -> None:
    page.locator('input[accept=".json"]').set_input_files(str(session_path))
    page.wait_for_function(
        """
        () => {
          const app = window.__GBDRAW_APP__;
          return Boolean(app && Array.isArray(app.results) && app.results.length > 0);
        }
        """,
        timeout=120000,
    )


def apply_capture_action(page, action: dict[str, Any]) -> None:
    action_type = as_text(action.get("type")) or as_text(action.get("action"))
    selector = as_text(action.get("selector"))
    value = action.get("value")

    if action_type == "click":
        page.locator(selector).click()
    elif action_type == "fill":
        page.locator(selector).fill("" if value is None else str(value))
    elif action_type == "select":
        page.locator(selector).select_option(str(value))
    elif action_type == "check":
        page.locator(selector).check()
    elif action_type == "uncheck":
        page.locator(selector).uncheck()
    elif action_type == "waitForSelector":
        page.locator(selector).wait_for(state=as_text(action.get("state")) or "visible")
    elif action_type == "evaluate":
        script = as_text(action.get("script"))
        if script:
            page.evaluate(script)
    elif action_type == "wait":
        page.wait_for_timeout(int(action.get("ms") or 250))
    else:
        raise RuntimeError(f"Unsupported capture action: {action_type}")


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
        wait_for_web_app_ready(page)
        session_ref = as_text(capture.get("session")) or as_text(sample.get("session"))
        session_path = resolve_gallery_reference(session_ref)
        if session_path is not None and session_path.exists():
            load_web_app_session(page, session_path)
    else:
        raise RuntimeError(f"Unsupported capture source: {source}")

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


def capture_operation(page, sample: dict[str, Any], operation: dict[str, Any], base_url: str, quality: int) -> Path:
    capture = operation.get("capture")
    if not isinstance(capture, dict):
        raise RuntimeError("operation is missing capture metadata.")

    target = output_path_for_operation(operation)
    if target is None:
        raise RuntimeError("operation media must point at a local output path.")

    prepare_capture_page(page, base_url, sample, capture)

    selector = as_text(capture.get("selector"))
    if selector:
        png_bytes = page.locator(selector).screenshot(type="png")
    else:
        box = capture.get("boundingBox")
        if not isinstance(box, dict):
            raise RuntimeError("capture metadata needs selector or boundingBox.")
        png_bytes = page.screenshot(
            type="png",
            clip={
                "x": float(box["x"]),
                "y": float(box["y"]),
                "width": float(box["width"]),
                "height": float(box["height"]),
            },
        )

    target.parent.mkdir(parents=True, exist_ok=True)
    with Image.open(io.BytesIO(png_bytes)) as image:
        image.convert("RGB").save(target, "WEBP", quality=quality, method=6)
    return target


def capture_examples(samples: list[dict[str, Any]], *, operation_filter: str, quality: int) -> int:
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
        return 0

    sync_playwright = import_playwright()
    with serve_web_root() as base_url, sync_playwright() as playwright:
        try:
            browser = playwright.chromium.launch(headless=True)
        except Exception as exc:  # pragma: no cover - environment dependent
            raise RuntimeError(
                "Could not launch Playwright Chromium. In sandboxed agent environments, rerun this "
                f"command with sandbox escalation. Underlying error: {exc}"
            ) from exc
        page = browser.new_page()
        try:
            for sample, context, operation in capturable:
                target = capture_operation(page, sample, operation, base_url, quality)
                rel = target.relative_to(REPO_ROOT)
                with Image.open(target) as image:
                    print(f"captured {context.label}: {rel} {image.size[0]}x{image.size[1]}")
        finally:
            browser.close()
    return len(capturable)


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

    capture_examples(samples, operation_filter=args.operation, quality=args.quality)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
