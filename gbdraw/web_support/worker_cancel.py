"""The Web Worker's cancel request.

The page writes SIGINT into the Worker's Pyodide interrupt buffer when the user
cancels; Pyodide delivers it to ``on_signal``, which only records the request.
The prepared-input cache stops the request at its next checkpoint
(``gbdraw.api.prepared.cancel_checkpoint``), a point between builds, so no
build, cache commit, font table or module body is left half done. The Worker's
drain step clears the request before it tells the page it is ready again.
"""

from __future__ import annotations

from types import FrameType


class WorkerCancelRequest:
    """Whether the page has canceled the Worker's current requests."""

    __slots__ = ("_requested",)

    def __init__(self) -> None:
        self._requested = False

    def on_signal(self, signum: int, frame: FrameType | None) -> None:
        """The SIGINT handler: record the request; never raise here."""

        self._requested = True

    def check(self) -> None:
        """Raise ``KeyboardInterrupt`` at a checkpoint while a cancel is requested."""

        if self._requested:
            raise KeyboardInterrupt

    def clear(self) -> None:
        """The drain step: later requests run again."""

        self._requested = False
