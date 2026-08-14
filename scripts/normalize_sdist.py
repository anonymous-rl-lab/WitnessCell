#!/usr/bin/env python3
"""Normalize an sdist tarball for byte-reproducible release builds."""

from __future__ import annotations

import argparse
import gzip
import os
import tarfile
from io import BytesIO
from pathlib import Path, PurePosixPath
from tempfile import NamedTemporaryFile


def _epoch(argument: int | None) -> int:
    raw = str(argument) if argument is not None else os.environ.get("SOURCE_DATE_EPOCH")
    if raw is None:
        raise ValueError("provide --epoch or set SOURCE_DATE_EPOCH")
    value = int(raw)
    if value < 0:
        raise ValueError("SOURCE_DATE_EPOCH must be non-negative")
    return value


def normalize(path: Path, epoch: int) -> None:
    records: list[tuple[tarfile.TarInfo, bytes | None]] = []
    with tarfile.open(path, "r:gz") as source:
        for member in source.getmembers():
            name = PurePosixPath(member.name)
            if name.is_absolute() or ".." in name.parts:
                raise ValueError(f"unsafe sdist member: {member.name}")
            if not (member.isdir() or member.isfile()):
                raise ValueError(f"unsupported sdist member type: {member.name}")
            payload = None
            if member.isfile():
                extracted = source.extractfile(member)
                if extracted is None:
                    raise ValueError(f"cannot read sdist member: {member.name}")
                payload = extracted.read()
            records.append((member, payload))

    with NamedTemporaryFile(
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        temporary = Path(handle.name)
        with gzip.GzipFile(
            filename="",
            mode="wb",
            fileobj=handle,
            compresslevel=9,
            mtime=epoch,
        ) as compressed:
            with tarfile.open(
                fileobj=compressed,
                mode="w",
                format=tarfile.PAX_FORMAT,
            ) as destination:
                for original, payload in sorted(records, key=lambda item: item[0].name):
                    member = tarfile.TarInfo(original.name)
                    member.type = tarfile.DIRTYPE if original.isdir() else tarfile.REGTYPE
                    member.mode = 0o755 if original.isdir() else 0o644
                    member.mtime = epoch
                    member.uid = 0
                    member.gid = 0
                    member.uname = ""
                    member.gname = ""
                    if payload is None:
                        member.size = 0
                        destination.addfile(member)
                    else:
                        member.size = len(payload)
                        destination.addfile(member, BytesIO(payload))
    os.replace(temporary, path)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("sdist", type=Path)
    parser.add_argument("--epoch", type=int)
    arguments = parser.parse_args()
    normalize(arguments.sdist, _epoch(arguments.epoch))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
