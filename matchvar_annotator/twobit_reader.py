# -*- coding: utf-8 -*-
"""Minimal UCSC .2bit sequence reader (pure Python)."""
from __future__ import annotations

import struct
from pathlib import Path
from typing import Dict, Optional, Tuple, Union

PathLike = Union[str, Path]


class TwoBitFile:
    """Read DNA from a UCSC 2bit file (forward genomic strand)."""

    def __init__(self, path: PathLike):
        self.path = Path(path)
        self._fh = self.path.open("rb")
        header = self._fh.read(16)
        signature = struct.unpack(">I", header[:4])[0]
        if signature == 0x1A412743:
            self._endian = ">"
        elif signature == 0x4327411A:
            self._endian = "<"
        else:
            raise ValueError(f"not a valid 2bit file: {self.path}")
        version, seq_count, _reserved = struct.unpack(f"{self._endian}III", header[4:16])
        if version != 0:
            raise ValueError(f"unsupported 2bit version {version}: {self.path}")

        self._index: Dict[str, int] = {}
        for _ in range(seq_count):
            name_size = struct.unpack(f"{self._endian}B", self._fh.read(1))[0]
            name = self._fh.read(name_size).decode("ascii")
            offset = struct.unpack(f"{self._endian}I", self._fh.read(4))[0]
            self._index[name] = offset

    def close(self) -> None:
        if self._fh:
            self._fh.close()
            self._fh = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()

    def resolve_chrom(self, chrom: str) -> Optional[str]:
        if chrom in self._index:
            return chrom
        alt = chrom[3:] if chrom.lower().startswith("chr") else f"chr{chrom}"
        if alt in self._index:
            return alt
        if chrom.upper() in ("M", "MT", "CHRM", "CHRMT"):
            for cand in ("chrM", "chrMT", "M", "MT"):
                if cand in self._index:
                    return cand
        return None

    def chroms(self) -> Dict[str, int]:
        out: Dict[str, int] = {}
        for name, offset in self._index.items():
            dna_size, _n_blocks = self._read_record_header(offset)[:2]
            out[name] = dna_size
        return out

    def _read_record_header(self, offset: int) -> Tuple[int, int, int]:
        self._fh.seek(offset)
        dna_size = struct.unpack(f"{self._endian}I", self._fh.read(4))[0]
        n_blocks = struct.unpack(f"{self._endian}I", self._fh.read(4))[0]
        self._fh.seek(8 * n_blocks, 1)
        mask_blocks = struct.unpack(f"{self._endian}I", self._fh.read(4))[0]
        self._fh.seek(8 * mask_blocks, 1)
        self._fh.read(4)
        packed_offset = self._fh.tell()
        return dna_size, n_blocks, packed_offset

    def fetch(self, chrom: str, start: int, end: int) -> str:
        """Fetch [start, end) 0-based half-open sequence (uppercase)."""
        resolved = self.resolve_chrom(chrom)
        if resolved is None:
            raise KeyError(f"chromosome not in 2bit: {chrom}")
        chrom = resolved
        if end < start:
            raise ValueError("end must be >= start")
        if start < 0:
            raise ValueError("start must be >= 0")

        dna_size, _n_blocks, packed_offset = self._read_record_header(self._index[chrom])
        if end > dna_size:
            raise ValueError(f"interval exceeds chromosome length {chrom}:{dna_size}")

        byte_start = start // 4
        byte_end = (end + 3) // 4
        self._fh.seek(packed_offset + byte_start)
        raw = self._fh.read(byte_end - byte_start)
        bases = []
        table = ("T", "C", "A", "G")
        for byte in raw:
            bases.append(table[(byte >> 6) & 3])
            bases.append(table[(byte >> 4) & 3])
            bases.append(table[(byte >> 2) & 3])
            bases.append(table[byte & 3])
        slice_start = start % 4
        return "".join(bases[slice_start:slice_start + (end - start)])
