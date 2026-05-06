# -*- coding: utf-8 -*-

"""
Copyright [2009-2026] EMBL-European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

import logging
import os
import tempfile
import typing as ty
import weakref
from pathlib import Path

from sqlitedict import SqliteDict

LOGGER = logging.getLogger(__name__)


def _iter_string_attrs(item: ty.Any) -> ty.Iterator[str]:
    cls_slots = getattr(type(item), "__slots__", None)
    if cls_slots is not None:
        for name in cls_slots:
            v = getattr(item, name, None)
            if isinstance(v, str):
                yield v
        return
    instance_dict = getattr(item, "__dict__", None)
    if instance_dict is not None:
        for v in instance_dict.values():
            if isinstance(v, str):
                yield v


def _estimate_size(key: str, value: ty.Any) -> int:
    """
    Cheap proxy for "bytes of payload" — sums string lengths of the key and any
    string attributes reachable on the value. Not RAM-accurate (Python object
    overhead is significant), but monotonic and tunable: set the threshold below
    your real memory budget and tune from there.
    """
    total = len(key)
    if isinstance(value, list):
        for item in value:
            for s in _iter_string_attrs(item):
                total += len(s)
    return total


class SpillDict:
    """
    Mapping that starts in memory and spills to a SqliteDict on disk once the
    estimated payload size crosses a threshold. After spill, all reads and
    writes go to disk; the in-memory dict is migrated and cleared.

    Only implements the slice of the Mapping protocol that ena.Context uses.
    """

    def __init__(self, threshold_bytes: int, spill_path: ty.Optional[Path] = None):
        self._mem: ty.Dict[str, ty.Any] = {}
        self._disk: ty.Optional[SqliteDict] = None
        self._size = 0
        self._threshold = threshold_bytes
        self._spill_path = spill_path
        self._owns_spill_file = spill_path is None
        self._finalizer: ty.Optional[weakref.finalize] = None

    def __setitem__(self, key, value):
        if self._disk is not None:
            self._disk[key] = value
            return
        self._mem[key] = value
        self._size += _estimate_size(key, value)
        if self._size >= self._threshold:
            self._spill()

    def __getitem__(self, key):
        if self._disk is not None:
            return self._disk[key]
        return self._mem[key]

    def __contains__(self, key):
        if self._disk is not None:
            return key in self._disk
        return key in self._mem

    def __len__(self):
        if self._disk is not None:
            return len(self._disk)
        return len(self._mem)

    def commit(self):
        if self._disk is not None:
            self._disk.commit()

    @property
    def spilled(self) -> bool:
        return self._disk is not None

    def close(self):
        if self._finalizer is not None and self._finalizer.alive:
            self._finalizer()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False

    @staticmethod
    def _cleanup(disk: ty.Optional[SqliteDict], path: ty.Optional[Path], owns: bool):
        if disk is not None:
            try:
                disk.close()
            except Exception:
                LOGGER.exception("Error closing SqliteDict at %s", path)
        if owns and path is not None:
            try:
                os.unlink(path)
            except FileNotFoundError:
                pass
            except OSError:
                LOGGER.exception("Error removing SpillDict temp file %s", path)

    def _spill(self):
        if self._spill_path is None:
            fd, name = tempfile.mkstemp(prefix="ena-dr-")
            os.close(fd)
            self._spill_path = Path(name)
        LOGGER.info(
            "SpillDict crossing threshold (%d entries, ~%d payload bytes); "
            "spilling to %s",
            len(self._mem),
            self._size,
            self._spill_path,
        )
        self._disk = SqliteDict(filename=str(self._spill_path))
        for k, v in self._mem.items():
            self._disk[k] = v
        self._disk.commit()
        self._mem.clear()
        self._finalizer = weakref.finalize(
            self, self._cleanup, self._disk, self._spill_path, self._owns_spill_file
        )
