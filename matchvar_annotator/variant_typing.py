# -*- coding: utf-8 -*-
"""MATCHVAR 变异类型命名（ExonicEffect / VarType）。

规则：
- 外显子：标准类型（missense / synonymous / nonsense / stoploss / frameshift / inframe / splicing）
- 内含子：SNV 或 INDEL
- UTR / 上游 / 下游 / 基因间：RV（Regulatory variation）
- 其它区域（ncRNA 等）：按等位基因长度归为 SNV / INDEL

多注释时与 Gene / cHGVS 条目数对齐，逗号分隔。
"""

from __future__ import annotations

import re
from typing import List, Optional, Sequence, Tuple


def split_annotation_parts(text: str) -> List[str]:
    """按逗号/分号拆分注释字段（忽略空段）。"""
    if not text or str(text).strip() in ("", ".", "NA", "N/A", "None", "nan"):
        return []
    parts = re.split(r"[;,]", str(text))
    return [p.strip() for p in parts if p.strip()]


def snv_or_indel(ref: str, alt: str) -> str:
    """根据 Ref/Alt 判定 SNV 或 INDEL。"""
    r = str(ref or "").replace("-", "").replace(".", "").replace("*", "")
    a = str(alt or "").replace("-", "").replace(".", "").replace("*", "")
    if len(r) == 1 and len(a) == 1 and r and a:
        return "SNV"
    return "INDEL"


def _region_of_part(chgvs_part: str, function_part: str) -> str:
    """返回 exonic / intronic / utr / splicing / other。

    upstream / downstream / intergenic 按调控区处理，归入 utr（→ RV）。
    """
    c = (chgvs_part or "").lower()
    f = (function_part or "").strip().lower()
    blob = f"{c} {f}"

    if "splicing" in blob or "r.spl" in c:
        return "splicing"
    # UTR 及近基因/基因间调控区（上游、下游、intergenic）一律按调控变异
    if (
        "utr5" in blob
        or "utr3" in blob
        or re.search(r"\butr\b", blob)
        or "upstream" in blob
        or "downstream" in blob
        or "intergenic" in blob
        or f in ("upstream", "downstream", "intergenic")
        or re.search(r"c\.\-\d", c)
        or re.search(r"c\.\*\d", c)
    ):
        return "utr"
    if "intronic" in blob or "intron" in c or f == "intronic":
        return "intronic"
    # c.123+45 / c.123-45 为内含子偏移
    if re.search(r"c\.\d+[+-]\d+", c):
        return "intronic"
    if f == "exonic" or ":exon" in c or re.search(r":exon\d+", c):
        return "exonic"
    # c.HGVS 在外显子上常见且不含 intron/utr
    if re.search(r":c\.", c) and "intron" not in c and "utr" not in c and "intronic" not in c:
        # 仍可能是近剪接；无 intron 标记时按外显子
        if "exonic" in f or f in ("", "exonic"):
            return "exonic"
    if "ncrna" in f or f == "ncrna":
        return "other"
    return "other"


def _exonic_type_from_text(exonic_effect: str, chgvs_part: str, p_part: str = "") -> str:
    """外显子标准变异类型。"""
    blob = f"{exonic_effect or ''} {chgvs_part or ''} {p_part or ''}".lower()

    if "splicing" in blob or "r.spl" in blob:
        return "splicing"
    if "stoploss" in blob:
        return "stoploss"
    if "stopgain" in blob or "nonsense" in blob or re.search(r"p\.\w+ter\b", blob):
        return "nonsense"
    # nonframeshift 必须先于 frameshift（后者是前者的子串）
    if "nonframeshift" in blob or "inframe" in blob:
        return "inframe"
    if re.search(r"p\.[A-Za-z]{3}\d+_[A-Za-z]{3}\d+del", blob):
        return "inframe"
    if "frameshift" in blob or "fs*" in blob or re.search(r"p\.\w+fs", blob):
        return "frameshift"
    # 必须先判 nonsynonymous / missense，再判 synonymous（子串包含关系）
    if "nonsynonymous" in blob or "missense" in blob:
        return "missense"
    if re.search(r"(?<!non)synonymous", blob) or "p.(=)" in blob or re.search(r"p\.\(=\)", blob) or "p.=" in blob:
        return "synonymous"
    # 从 p.HGVS 推断
    if re.search(r"p\.\(?=", chgvs_part or "") or "p.(=)" in (chgvs_part or ""):
        return "synonymous"
    if re.search(r"p\.[A-Za-z]{3}\d+[A-Za-z]{3}", chgvs_part or p_part or ""):
        return "missense"
    if "duplication" in blob and "frameshift" not in blob:
        return "inframe"
    return "missense"  # 外显子默认


def classify_one(
    *,
    function_part: str,
    chgvs_part: str,
    ref: str,
    alt: str,
    exonic_effect_hint: str = "",
    p_hgvs_part: str = "",
) -> str:
    """对单条注释返回统一类型名。"""
    region = _region_of_part(chgvs_part, function_part)
    if region == "utr":
        return "RV"
    if region == "intronic":
        return snv_or_indel(ref, alt)
    if region == "splicing":
        return "splicing"
    if region == "exonic":
        return _exonic_type_from_text(exonic_effect_hint, chgvs_part, p_hgvs_part)
    # ncRNA 等
    return snv_or_indel(ref, alt)


def _pad_parts(parts: Sequence[str], n: int) -> List[str]:
    if n <= 0:
        return []
    items = list(parts)
    if not items:
        return [""] * n
    if len(items) >= n:
        return items[:n]
    # 不足时用最后一个补齐
    return items + [items[-1]] * (n - len(items))


def classify_row(
    *,
    function: str,
    gene: str,
    c_hgvs: str,
    ref: str,
    alt: str,
    exonic_effect: str = "",
    p_hgvs: str = "",
) -> Tuple[str, str]:
    """
    返回 (ExonicEffect, VarType)，二者内容相同，逗号分隔。

    条目数优先对齐 cHGVS，其次 Gene，再次 Function。
    """
    c_parts = split_annotation_parts(c_hgvs)
    g_parts = split_annotation_parts(gene)
    f_parts = split_annotation_parts(function)
    e_parts = split_annotation_parts(exonic_effect)
    p_parts = split_annotation_parts(p_hgvs)

    n = len(c_parts) or len(g_parts) or len(f_parts) or 1

    c_parts = _pad_parts(c_parts, n)
    f_parts = _pad_parts(f_parts, n)
    e_parts = _pad_parts(e_parts, n)
    p_parts = _pad_parts(p_parts, n)

    types: List[str] = []
    for i in range(n):
        types.append(
            classify_one(
                function_part=f_parts[i] if i < len(f_parts) else "",
                chgvs_part=c_parts[i] if i < len(c_parts) else "",
                ref=ref,
                alt=alt,
                exonic_effect_hint=e_parts[i] if i < len(e_parts) else (exonic_effect or ""),
                p_hgvs_part=p_parts[i] if i < len(p_parts) else (p_hgvs or ""),
            )
        )

    joined = ",".join(types)
    return joined, joined
