# -*- coding: utf-8 -*-
"""Shared ground-truth pathogenicity labels for pipeline / metrics / visualization."""

from __future__ import annotations

import re
from typing import Iterable, Optional

PATHOGENIC_FUNC_RE = re.compile(
    r"splicing|stopgain|stoploss|frameshift|nonsynonymous|nonsense|missense",
    re.IGNORECASE,
)

PATHOGENIC_TYPES = frozenset({
    "SPLICING",
    "SPLICE_SITE",
    "FRAMESHIFT",
    "NONSENSE",
    "STOPLOSS",
    "STOPGAIN",
})

_TYPE_RE = re.compile(r"(?:^|;)TYPE=([^;]+)")
_FRM_EQ_RE = re.compile(r"(?:^|;)FRAMESHIFT=([^;]+)")
_FRM_FLAG_RE = re.compile(r"(?:^|;)FRAMESHIFT(?:;|$)")

FUNC_LABEL_COLUMNS = (
    "Function.refGene",
    "Func.refGene",
    "ExonicEffect.refGene",
    "ExonicFunc.refGene",
    "Function.ensGene",
    "Func.ensGene",
    "ExonicEffect.ensGene",
    "ExonicFunc.ensGene",
    "Function.knownGene",
    "Func.knownGene",
    "ExonicEffect.knownGene",
    "ExonicFunc.knownGene",
)


def frameshift_from_info(info: str) -> bool:
    """True if INFO has FRAMESHIFT=true or a bare FRAMESHIFT flag."""
    text = info or ""
    match = _FRM_EQ_RE.search(text)
    if match:
        return match.group(1).lower().strip() in ("true", "1", "yes")
    return bool(_FRM_FLAG_RE.search(text))


def label_from_vcf_info(info: str) -> int:
    """Return 1 for pathogenic VCF INFO, else 0."""
    text = info or ""
    type_match = _TYPE_RE.search(text)
    variant_type = type_match.group(1).upper().strip() if type_match else ""
    return int(variant_type in PATHOGENIC_TYPES or frameshift_from_info(text))


def func_is_pathogenic(text: str) -> bool:
    return bool(PATHOGENIC_FUNC_RE.search(str(text or "")))


def first_present_column(columns: Iterable[str], candidates: Iterable[str]) -> Optional[str]:
    colset = set(columns)
    for name in candidates:
        if name in colset:
            return name
    return None
