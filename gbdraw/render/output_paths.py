"""Filesystem validation shared by diagram and session writers."""

from __future__ import annotations

import errno
import os
from contextlib import contextmanager, suppress
from pathlib import Path
from tempfile import mkdtemp
from typing import Iterable, Iterator, Mapping

from gbdraw.exceptions import ExportError, ValidationError

# How each surface spells the permission to replace an output.
PYTHON_OVERWRITE_OPTION = "overwrite=True"
CLI_OVERWRITE_OPTION = "--overwrite"


class OutputExistsError(ValidationError):
    """Outputs exist and the caller did not allow replacing them.

    The one owner of the refusal's wording. ``str()`` names Python's
    ``overwrite=True``; another surface words it with :meth:`message` and its
    own option (the CLI prints ``--overwrite``).
    """

    def __init__(
        self,
        paths: Iterable[str | Path],
        *,
        subject: str = "Output file",
        written: Iterable[str | Path] = (),
    ) -> None:
        self.paths = tuple(str(path) for path in paths)
        self.subject = subject
        self.written = tuple(str(path) for path in written)
        super().__init__(self.message(PYTHON_OVERWRITE_OPTION))

    def message(self, overwrite_option: str) -> str:
        one = len(self.paths) == 1
        text = (
            f"{self.subject}{'' if one else 's'} already exist{'s' if one else ''}: "
            f"{', '.join(self.paths)}. Use {overwrite_option} to replace {'it' if one else 'them'}."
        )
        if self.written:
            text += f" Already written: {', '.join(self.written)}."
        return text

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
        raise OutputExistsError(existing)
    return output_paths


def _reserve_and_replace(staged_path: Path, output_path: Path) -> None:
    """No-replace commit without hard links: reserve the target, then rename onto it.

    The exclusive create refuses an existing target, as a link would, and the
    rename moves the staged file itself, so the output keeps the staged file's
    mode, as a link does. A failed rename removes the reservation again.
    """

    os.close(os.open(output_path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600))
    try:
        os.replace(staged_path, output_path)
    except BaseException:
        with suppress(OSError):
            os.unlink(output_path)
        raise


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
    if not hasattr(os, "link"):
        # Pyodide's file system, where the Web writes its outputs, has no os.link.
        _reserve_and_replace(staged, output)
        return
    try:
        os.link(staged, output, follow_symlinks=False)
    except NotImplementedError:
        _reserve_and_replace(staged, output)
        return
    except OSError as exc:
        fallback_errors = {
            errno.ENOSYS,
            errno.EPERM,
            getattr(errno, "EOPNOTSUPP", errno.EPERM),
        }
        if exc.errno not in fallback_errors:
            raise
        _reserve_and_replace(staged, output)
        return
    staged.unlink()


@contextmanager
def staged_output_files(output_paths: Iterable[str | Path]) -> Iterator[dict[Path, Path]]:
    """Map each output to a staged path in its own directory, for one output set.

    The staged files live in a private directory beside their outputs, so a
    commit is a rename or link on one filesystem and each file gets the usual
    permissions. Writers fill every staged file first and then call
    :func:`commit_staged_outputs`; on exit, whatever was not committed is
    removed, so a failed generation leaves no output behind.
    """

    directories: dict[Path, Path] = {}
    staged: dict[Path, Path] = {}
    try:
        for output_path in (Path(path) for path in output_paths):
            parent = output_path.parent
            if parent not in directories:
                try:
                    directories[parent] = Path(mkdtemp(prefix=".gbdraw-staged-", dir=parent))
                except OSError as exc:
                    raise ExportError(f"Could not stage output files in {parent}: {exc}.") from exc
            staged[output_path] = directories[parent] / output_path.name
        yield staged
    finally:
        # Only the staged files can be left: no recursive removal is needed.
        for staged_path in staged.values():
            with suppress(OSError):
                staged_path.unlink(missing_ok=True)
        for directory in directories.values():
            with suppress(OSError):
                directory.rmdir()


def commit_staged_outputs(staged: Mapping[Path, Path], *, overwrite: bool) -> None:
    """Commit the staged files of :func:`staged_output_files` in order.

    A target that appeared after the preflight stops a no-replace commit
    (``ValidationError``) and any other filesystem failure raises
    ``ExportError``; both name the outputs committed before it.
    """

    committed: list[str] = []
    for output_path, staged_path in staged.items():
        try:
            commit_staged_output_file(staged_path, output_path, overwrite=overwrite)
        except OSError as exc:
            if isinstance(exc, FileExistsError):
                raise OutputExistsError((output_path,), written=committed) from exc
            already = f" Already written: {', '.join(committed)}." if committed else ""
            raise ExportError(f"Could not write output file {output_path}: {exc}.{already}") from exc
        committed.append(str(output_path))


__all__ = [
    "commit_staged_output_file",
    "commit_staged_outputs",
    "is_windows_reserved_filename_component",
    "CLI_OVERWRITE_OPTION",
    "OutputExistsError",
    "PYTHON_OVERWRITE_OPTION",
    "preflight_output_paths",
    "staged_output_files",
]
