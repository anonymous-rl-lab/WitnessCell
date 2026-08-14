#!/usr/bin/env python3
"""Require an exact vX.Y.Z tag matching project.version."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import tomllib


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tag", required=True)
    parser.add_argument("--pyproject", type=Path, default=Path("pyproject.toml"))
    arguments = parser.parse_args()
    match = re.fullmatch(r"v(0|[1-9]\d*)\.(0|[1-9]\d*)\.(0|[1-9]\d*)", arguments.tag)
    if match is None:
        raise SystemExit("release tag must match vX.Y.Z exactly")
    with arguments.pyproject.open("rb") as handle:
        version = str(tomllib.load(handle)["project"]["version"])
    if arguments.tag != f"v{version}":
        raise SystemExit(
            f"release tag {arguments.tag!r} does not match project version {version!r}"
        )
    print(f"RELEASE TAG VERIFIED: {arguments.tag}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
