"""Runtime read firewall for candidate-blind Experiment 22 Phase M."""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path


def _resolve(value: str) -> str:
    try:
        return str(Path(value).resolve())
    except (OSError, TypeError, ValueError):
        return str(value)


_forbidden = tuple(
    _resolve(item)
    for item in json.loads(os.environ.get("E22_FORBIDDEN_READ_ROOTS", "[]"))
)
_audit_log = os.environ.get("E22_FIREWALL_AUDIT_LOG")


def audit(event: str, args: tuple) -> None:
    if event != "open" or not args:
        return
    path = args[0]
    if not isinstance(path, (str, bytes, os.PathLike)):
        return
    resolved = _resolve(os.fsdecode(path))
    denied = any(resolved == root or resolved.startswith(root + os.sep) for root in _forbidden)
    if denied:
        if _audit_log:
            with Path(_audit_log).open("a", encoding="utf-8") as handle:
                handle.write(json.dumps({"event": event, "path": resolved, "decision": "DENY"}) + "\n")
        raise PermissionError(f"Experiment 22 Phase M firewall denied read: {resolved}")


sys.addaudithook(audit)

