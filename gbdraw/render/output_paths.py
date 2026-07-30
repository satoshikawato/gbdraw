"""Filesystem validation shared by diagram and session writers."""

from __future__ import annotations

import errno
import os
import shutil
from pathlib import Path
from typing import Iterable

from gbdraw.exceptions import ValidationError

_WINDOWS_DEVICE_NAMES = frozenset(
    {
        "CON",
        "PRN",
        "AUX",
        "NUL",
        "CONIN$",
        "CONOUT$",
        *(f"COM{suffix}" for suffix in "123456789¹²³"),
        *(f"LPT{suffix}" for suffix in "123456789¹²³"),
    }
)
_WINDOWS_RESERVED_CHARACTERS = frozenset('"*:<>?\\|')


def is_windows_reserved_filename_component(value: object) -> bool:
    """Return whether one filename component is reserved on Windows."""

    component = str(value)
    if not component or component[-1] in {" ", "."}:
        return True
    if any(
        ord(character) < 32
        or character in _WINDOWS_RESERVED_CHARACTERS
        for character in component
    ):
        return True
    device_name = component.split(".", 1)[0].rstrip(" ").upper()
    return device_name in _WINDOWS_DEVICE_NAMES


def preflight_output_paths(
    paths: Iterable[str | Path],
    *,
    overwrite: bool,
) -> tuple[Path, ...]:
    """Validate a complete output set without changing the filesystem."""

    output_paths = tuple(Path(path) for path in paths)
    invalid_targets: list[Path] = []
    invalid_parents: list[Path] = []
    existing: list[Path] = []
    for path in output_paths:
        raw_path = str(path)
        if any(
            ord(character) < 32 or ord(character) == 127
            for character in raw_path
        ):
            raise ValidationError(
                f"Output path must not contain ASCII control characters: {raw_path!r}."
            )
        reserved_component = next(
            (
                component
                for component in path.parts
                if component not in {path.anchor, ".", ".."}
                and is_windows_reserved_filename_component(component)
            ),
            None,
        )
        if reserved_component is not None:
            raise ValidationError(
                "Output path contains a Windows-reserved filename component: "
                f"{reserved_component!r}."
            )
        try:
            path_is_symlink = path.is_symlink()
            path_exists = path.exists()
            if (
                (path_is_symlink and not path_exists)
                or (path_exists and not path.is_file())
            ):
                invalid_targets.append(path)

            existing_ancestor = path.parent
            while (
                not existing_ancestor.exists()
                and not existing_ancestor.is_symlink()
                and existing_ancestor.parent != existing_ancestor
            ):
                existing_ancestor = existing_ancestor.parent
            if (
                existing_ancestor.is_symlink()
                and not existing_ancestor.exists()
            ) or (
                existing_ancestor.exists()
                and not existing_ancestor.is_dir()
            ):
                invalid_parents.append(existing_ancestor)
            if path_exists:
                existing.append(path)
        except (OSError, ValueError) as exc:
            raise ValidationError(
                f"Could not inspect output path: {path}."
            ) from exc

    if invalid_targets:
        raise ValidationError(
            "Output target(s) are not replaceable files: "
            + ", ".join(str(path) for path in invalid_targets)
            + "."
        )
    if invalid_parents:
        raise ValidationError(
            "Output parent path(s) are not directories: "
            + ", ".join(str(path) for path in invalid_parents)
            + "."
        )
    if existing and not overwrite:
        raise ValidationError(
            "Output file(s) already exist: "
            + ", ".join(str(path) for path in existing)
            + ". Use overwrite=True to replace."
        )
    return output_paths


def _copy_staged_file_exclusively(
    staged_path: Path,
    output_path: Path,
) -> None:
    """Fallback no-replace commit for filesystems without hard links."""

    output_fd: int | None = None
    try:
        output_fd = os.open(
            output_path,
            os.O_WRONLY | os.O_CREAT | os.O_EXCL,
            0o600,
        )
        with staged_path.open("rb") as source:
            output_file = os.fdopen(output_fd, "wb")
            output_fd = None
            with output_file:
                shutil.copyfileobj(source, output_file)
                output_file.flush()
                os.fsync(output_file.fileno())
    except Exception:
        if output_fd is not None:
            try:
                os.close(output_fd)
            except OSError:
                pass
        raise
    staged_path.unlink()


def commit_staged_output_file(
    staged_path: str | Path,
    output_path: str | Path,
    *,
    overwrite: bool,
) -> None:
    """Commit a same-directory staged file without following the final path."""

    staged = Path(staged_path)
    output = Path(output_path)
    if overwrite:
        staged.replace(output)
        return
    if os.name == "nt":
        staged.rename(output)
        return
    try:
        os.link(staged, output, follow_symlinks=False)
    except NotImplementedError:
        _copy_staged_file_exclusively(staged, output)
        return
    except OSError as exc:
        fallback_errors = {
            errno.ENOSYS,
            errno.EPERM,
            getattr(errno, "EOPNOTSUPP", errno.EPERM),
        }
        if exc.errno not in fallback_errors:
            raise
        _copy_staged_file_exclusively(staged, output)
        return
    staged.unlink()


__all__ = [
    "commit_staged_output_file",
    "is_windows_reserved_filename_component",
    "preflight_output_paths",
]
