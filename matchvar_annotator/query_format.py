# -*- coding: utf-8 -*-
"""Query-file formats: VCF (CHROM POS REF ALT) and MATCHVAR 5-column input.

VCF alleles are left-padded. Annotation uses 1-based closed intervals with
``-`` for an empty allele (insertion Ref / deletion Alt). That mapping lives
here so table annotation can read a VCF directly.
"""
from __future__ import annotations

import gzip
import os
import re
from typing import Dict, Iterator, List, Optional, Sequence, Tuple

_DNA = re.compile(r"^[ACGTN]+$", re.I)
_INT = re.compile(r"^-?\d+$")
_VCF_ALT_OK = re.compile(r"^[ACGTN*]+$", re.I)


def _is_int(text: str) -> bool:
    return bool(_INT.match(str(text or "").strip()))


def _is_dna(text: str) -> bool:
    return bool(_DNA.match(str(text or "").strip()))


def _is_allele(text: str) -> bool:
    t = str(text or "").strip()
    return t in ("-", ".", "*") or bool(_DNA.match(t))


def normalize_vcf_alleles(pos: int, ref: str, alt: str) -> Tuple[int, int, str, str]:
    """Map VCF POS/REF/ALT to 1-based closed query alleles.

    Insertion: ``start == end``, Ref ``-``, Alt = inserted bases (after POS).
    Deletion: Start/End span deleted bases, Alt ``-``.
    SNV: Start == End == POS.
    """
    start = int(pos)
    ref_s = (ref or "").upper()
    alt_s = (alt or "").upper()
    if not ref_s or not alt_s or ref_s in (".", "-") or alt_s.startswith("<"):
        end = start + max(len(ref_s), 1) - 1
        return start, end, ref_s or "-", alt_s or "-"

    if len(ref_s) == 1 and len(alt_s) == 1:
        return start, start, ref_s, alt_s

    # Trim identical prefix (VCF left anchor)
    i = 0
    while i < len(ref_s) and i < len(alt_s) and ref_s[i] == alt_s[i]:
        i += 1
    # Trim identical suffix
    j = 0
    while (
        i + j < len(ref_s)
        and i + j < len(alt_s)
        and ref_s[-(j + 1)] == alt_s[-(j + 1)]
    ):
        j += 1
    core_ref = ref_s[i:len(ref_s) - j] if j else ref_s[i:]
    core_alt = alt_s[i:len(alt_s) - j] if j else alt_s[i:]
    new_start = start + i

    if not core_ref and core_alt:
        # insertion after the last remaining ref base (the anchor)
        ins_pos = new_start - 1 if i > 0 else new_start
        if ins_pos < 1:
            ins_pos = 1
        return ins_pos, ins_pos, "-", core_alt
    if core_ref and not core_alt:
        return new_start, new_start + len(core_ref) - 1, core_ref, "-"
    if core_ref and core_alt:
        return new_start, new_start + len(core_ref) - 1, core_ref, core_alt
    return start, start, ref_s, alt_s


def sniff_query_format(path: str) -> str:
    """Return ``vcf`` or ``mvinput``."""
    name = os.path.basename(path or "").lower()
    if name.endswith(".vcf") or name.endswith(".vcf.gz"):
        return "vcf"
    try:
        with _open_text(path) as fh:
            for raw in fh:
                line = raw.strip()
                if not line:
                    continue
                if line.startswith("##fileformat=VCF") or line.upper().startswith("#CHROM"):
                    return "vcf"
                if line.startswith("#"):
                    continue
                parts = _split_fields(line)
                if len(parts) == 4 and _is_int(parts[1]) and _is_dna(parts[2]) and _looks_vcf_alt(parts[3]):
                    return "vcf"
                if len(parts) >= 5 and not _is_int(parts[2]) and _is_dna(parts[3]) and _looks_vcf_alt(parts[4]):
                    return "vcf"
                return "mvinput"
    except OSError:
        return "mvinput"
    return "mvinput"


def _looks_vcf_alt(text: str) -> bool:
    tok = (text or "").strip()
    if not tok or tok == ".":
        return False
    return all(_VCF_ALT_OK.match(a) or a.startswith("<") for a in tok.split(","))


def _split_fields(line: str) -> List[str]:
    if "\t" in line:
        return line.rstrip("\n\r").split("\t")
    return line.strip().split()


def open_query_text(path: str):
    """Open a query file (plain text or ``.gz``) as text."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def _open_text(path: str):
    return open_query_text(path)


def parse_query_line(line: str, fmt: str = "auto") -> List[Dict[str, object]]:
    """Parse one query line into one or more variant dicts."""
    line = (line or "").strip()
    if not line or line.startswith("#"):
        return []
    parts = _split_fields(line)
    if fmt == "auto":
        fmt = _guess_line_format(parts)
    if fmt == "vcf":
        return _parse_vcf_fields(parts, line)
    return _parse_mvinput_fields(parts, line)


def _guess_line_format(parts: Sequence[str]) -> str:
    if len(parts) == 4 and _is_int(parts[1]) and _is_dna(parts[2]):
        return "vcf"
    if len(parts) >= 5 and _is_int(parts[1]) and not _is_int(parts[2]) and _is_dna(parts[3]):
        return "vcf"
    return "mvinput"


def _parse_mvinput_fields(parts: Sequence[str], line: str) -> List[Dict[str, object]]:
    if len(parts) < 5:
        return []
    try:
        start = int(parts[1])
        end = int(parts[2])
    except (TypeError, ValueError):
        return []
    rec = {
        "chrom": parts[0],
        "start": start,
        "end": end,
        "ref": parts[3] if parts[3] not in (".",) else "-",
        "alt": parts[4] if parts[4] not in (".",) else "-",
        "original_line": line,
        "extra_fields": list(parts[5:]),
        "vcf_pos": None,
        "vcf_ref": None,
        "vcf_alt": None,
    }
    rec["query_line"] = f"{rec['chrom']}\t{rec['start']}\t{rec['end']}\t{rec['ref']}\t{rec['alt']}"
    return [rec]


def _parse_vcf_fields(parts: Sequence[str], line: str) -> List[Dict[str, object]]:
    if len(parts) < 4:
        return []
    chrom = parts[0]
    try:
        pos = int(parts[1])
    except (TypeError, ValueError):
        return []
    if len(parts) == 4:
        ref, alt_field = parts[2], parts[3]
        extra: List[str] = []
    else:
        # CHROM POS ID REF ALT [QUAL FILTER INFO ...]
        ref, alt_field = parts[3], parts[4]
        extra = [parts[2]] + list(parts[5:])
    out: List[Dict[str, object]] = []
    for alt in (alt_field or "").split(","):
        alt = alt.strip()
        if not alt or alt == ".":
            continue
        if alt.startswith("<"):
            continue
        start, end, nref, nalt = normalize_vcf_alleles(pos, ref, alt)
        out.append(
            {
                "chrom": chrom,
                "start": start,
                "end": end,
                "ref": nref,
                "alt": nalt,
                "original_line": line,
                "query_line": f"{chrom}\t{start}\t{end}\t{nref}\t{nalt}",
                "extra_fields": extra,
                "vcf_pos": pos,
                "vcf_ref": ref.upper(),
                "vcf_alt": alt.upper(),
            }
        )
    return out


def iter_query_records(path: str, fmt: Optional[str] = None) -> Iterator[Dict[str, object]]:
    kind = fmt or sniff_query_format(path)
    with _open_text(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            for rec in parse_query_line(line, fmt=kind):
                yield rec
