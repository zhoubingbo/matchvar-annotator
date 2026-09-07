# -*- coding: utf-8 -*-
"""
MATCHVAR 输出列名约定。

基因注释协议（refGene / ensGene / knownGene）固定 7 列：
  Function / Gene / cHGVS / MANE Select / ExonicEffect / VarType / pHGVS

坐标列之后追加 gHGVS（基因组 HGVS）。
"""
from __future__ import annotations

import re
from typing import Dict, List, Optional, Tuple

GENE_DBTYPES = ("refGene", "ensGene", "knownGene", "ncbiRefSeq", "gencode")

# 旧版列名 → MATCHVAR（兼容读历史结果；平台主链路已直接用 MATCHVAR 列名）
LEGACY_TO_MATCHVAR = {
    "Func": "Function",
    "GeneDetail": "cHGVS",
    "ExonicFunc": "ExonicEffect",
    "AAChange": "pHGVS",
}


def remap_legacy_column(name: str, protocol: Optional[str] = None) -> str:
    """将旧版列名映射为 MATCHVAR 列名。

    支持 ``Func.refGene`` / ``Func_refGene``；未知列原样返回。
    """
    if not name:
        return name
    sep = "." if "." in name else ("_" if "_" in name else "")
    if sep:
        head, _, tail = name.partition(sep)
        mapped = LEGACY_TO_MATCHVAR.get(head)
        if mapped:
            return f"{mapped}{sep}{tail}"
        return name
    mapped = LEGACY_TO_MATCHVAR.get(name)
    if mapped and protocol:
        return f"{mapped}.{protocol}"
    return mapped or name


def gene_column_names(protocol: str, underline: bool = False) -> List[str]:
    """返回某基因协议的 MATCHVAR 列名（顺序固定）。"""
    sep = "_" if underline else "."
    return [
        f"Function{sep}{protocol}",
        f"Gene{sep}{protocol}",
        f"cHGVS{sep}{protocol}",
        f"MANE Select{sep}{protocol}",
        f"ExonicEffect{sep}{protocol}",
        f"VarType{sep}{protocol}",
        f"pHGVS{sep}{protocol}",
    ]


def gene_column_keys(protocol: str, underline: bool = False) -> dict:
    """语义键 → 实际列名。"""
    names = gene_column_names(protocol, underline=underline)
    return {
        "function": names[0],
        "gene": names[1],
        "c_hgvs": names[2],
        "mane_select": names[3],
        "exonic_effect": names[4],
        "vartype": names[5],
        "p_hgvs": names[6],
    }


def annotation_headers_for_gene_dbs() -> dict:
    """供 ANNOTATION_HEADERS 使用的基因库展开定义。"""
    return {db: gene_column_names(db, underline=False) for db in GENE_DBTYPES}


def format_g_hgvs(chrom: str, start, end, ref: str, alt: str) -> str:
    """
    由基因组坐标生成 g.HGVS。

    期望输入为 1-based closed 区间 + VCF 风格 REF/ALT（indel 可含左侧锚碱基）。
    """
    try:
        chrom_s = str(chrom or "").strip()
        if not chrom_s:
            return "."
        # 保留 chr 前缀以与结果表 Chr 列一致
        if not chrom_s.lower().startswith("chr") and chrom_s.upper() != "MT":
            if chrom_s.isdigit() or chrom_s.upper() in ("X", "Y", "M", "MT"):
                chrom_s = f"chr{chrom_s}"

        start_i = int(float(str(start)))
        end_i = int(float(str(end)))
        ref_s = str(ref or "").replace("-", "").replace(".", "").replace("*", "").upper()
        alt_s = str(alt or "").replace("-", "").replace(".", "").replace("*", "").upper()

        # SNV
        if len(ref_s) == 1 and len(alt_s) == 1:
            return f"{chrom_s}:g.{start_i}{ref_s}>{alt_s}"

        # VCF-padded deletion: REF longer, shares left base
        if len(ref_s) > len(alt_s) and ref_s.startswith(alt_s) and len(alt_s) >= 1:
            deleted = ref_s[len(alt_s) :]
            del_start = start_i + len(alt_s)
            del_end = del_start + len(deleted) - 1
            if len(deleted) == 1:
                return f"{chrom_s}:g.{del_start}del{deleted}"
            return f"{chrom_s}:g.{del_start}_{del_end}del{deleted}"

        # VCF-padded insertion: ALT longer, shares left base
        if len(alt_s) > len(ref_s) and alt_s.startswith(ref_s) and len(ref_s) >= 1:
            inserted = alt_s[len(ref_s) :]
            pos = start_i + len(ref_s) - 1
            return f"{chrom_s}:g.{pos}_{pos + 1}ins{inserted}"

        # delins / substitution of unequal length without clean padding
        if ref_s and alt_s:
            if start_i == end_i and len(ref_s) == 1:
                return f"{chrom_s}:g.{start_i}delins{alt_s}"
            return f"{chrom_s}:g.{start_i}_{end_i}delins{alt_s}"

        if ref_s and not alt_s:
            if start_i == end_i:
                return f"{chrom_s}:g.{start_i}del{ref_s}"
            return f"{chrom_s}:g.{start_i}_{end_i}del{ref_s}"

        if alt_s and not ref_s:
            return f"{chrom_s}:g.{start_i}_{start_i + 1}ins{alt_s}"

        return f"{chrom_s}:g.{start_i}_{end_i}"
    except Exception:
        return "."


# NC_ accession (without version) → chr name
_NC_TO_CHR = {n: f"chr{n}" for n in range(1, 23)}
_NC_TO_CHR[23] = "chrX"
_NC_TO_CHR[24] = "chrY"

_G_PREFIX = re.compile(
    r"^(?:(?P<acc>NC_\d+\.\d+|chr[0-9XYMTmt]+|[0-9XYMTmt]+):)?g\.(?P<body>.+)$",
    re.I,
)
_G_SNV = re.compile(r"^(?P<pos>\d+)(?P<ref>[ACGTN]+)>(?P<alt>[ACGTN]+)$", re.I)
_G_DELINS = re.compile(
    r"^(?P<a>\d+)(?:_(?P<b>\d+))?del(?P<delseq>[ACGTN]*)ins(?P<ins>[ACGTN]+)$",
    re.I,
)
_G_DELINS2 = re.compile(
    r"^(?P<a>\d+)(?:_(?P<b>\d+))?delins(?P<ins>[ACGTN]+)$",
    re.I,
)
_G_DEL = re.compile(r"^(?P<a>\d+)(?:_(?P<b>\d+))?del(?P<delseq>[ACGTN]*)$", re.I)
_G_INS = re.compile(r"^(?P<a>\d+)_(?P<b>\d+)ins(?P<ins>[ACGTN]+)$", re.I)
_G_DUP = re.compile(
    r"^(?P<a>\d+)(?:_(?P<b>\d+))?dup(?P<dupseq>[ACGTN]*)$",
    re.I,
)


def _normalize_g_chrom(acc: Optional[str]) -> str:
    raw = (acc or "").strip()
    if not raw:
        return ""
    if raw.upper().startswith("NC_"):
        try:
            num = int(raw.split("_", 1)[1].split(".", 1)[0])
        except (IndexError, ValueError):
            return raw
        return _NC_TO_CHR.get(num, raw)
    if raw.upper() in ("M", "MT"):
        return "chrM"
    if not raw.lower().startswith("chr"):
        return f"chr{raw}"
    return raw if raw.lower().startswith("chr") else f"chr{raw}"


def parse_g_hgvs(hgvs: str) -> Dict[str, object]:
    """Parse genomic HGVS into MATCHVAR / ANNOVAR 1-based closed alleles.

    Returns dict: chrom, start, end, ref, alt, g_hgvs.
    Insertion: start=end=left flank, ref='-', alt=inserted sequence.
    Deletion: start/end span deleted bases, ref=deleted sequence (may be empty
    if the HGVS omitted it), alt='-'.
    """
    text = (hgvs or "").strip().replace(" ", "")
    if not text:
        raise ValueError("empty g.HGVS")
    m = _G_PREFIX.match(text)
    if not m:
        raise ValueError(f"not a g.HGVS string: {hgvs}")
    chrom = _normalize_g_chrom(m.group("acc"))
    body = m.group("body")

    def _out(start: int, end: int, ref: str, alt: str) -> Dict[str, object]:
        if not chrom:
            raise ValueError(f"g.HGVS missing chromosome: {hgvs}")
        return {
            "chrom": chrom,
            "start": int(start),
            "end": int(end),
            "ref": (ref or "-").upper() if ref else "-",
            "alt": (alt or "-").upper() if alt else "-",
            "g_hgvs": text if ":" in text else f"{chrom}:{text}",
        }

    m = _G_SNV.match(body)
    if m:
        pos = int(m.group("pos"))
        return _out(pos, pos, m.group("ref"), m.group("alt"))

    m = _G_DELINS.match(body) or _G_DELINS2.match(body)
    if m:
        a = int(m.group("a"))
        b = int(m.group("b") or a)
        deleted = (m.groupdict().get("delseq") or "")
        ins = m.group("ins") or ""
        return _out(a, b, deleted or "-", ins)

    m = _G_INS.match(body)
    if m:
        a = int(m.group("a"))
        ins = m.group("ins") or ""
        return _out(a, a, "-", ins)

    m = _G_DEL.match(body)
    if m:
        a = int(m.group("a"))
        b = int(m.group("b") or a)
        deleted = m.group("delseq") or ""
        return _out(a, b, deleted or "-", "-")

    m = _G_DUP.match(body)
    if m:
        a = int(m.group("a"))
        b = int(m.group("b") or a)
        dupseq = (m.group("dupseq") or "").upper()
        if dupseq:
            return _out(b, b, "-", dupseq)
        parsed = _out(a, b, "-", "-")
        parsed["_dup"] = True
        return parsed

    raise ValueError(f"unsupported g.HGVS body: {hgvs}")


def fill_g_hgvs_alleles(parsed: Dict[str, object], fetch_1based) -> Dict[str, object]:
    """Fill omitted REF/dup sequence using a 1-based closed genome fetch(chrom, start, end)."""
    ref = str(parsed.get("ref") or "")
    alt = str(parsed.get("alt") or "")
    chrom = str(parsed["chrom"])
    start = int(parsed["start"])
    end = int(parsed["end"])
    if parsed.get("_dup"):
        seq = fetch_1based(chrom, start, end).upper()
        parsed["ref"] = "-"
        parsed["alt"] = seq
        parsed["start"] = end
        parsed["end"] = end
        parsed.pop("_dup", None)
        return parsed
    if ref in ("", "-") and alt not in ("", "-") and start == end:
        # insertion — no REF to fetch
        parsed["ref"] = "-"
        return parsed
    if ref in ("", "-") and alt == "-":
        parsed["ref"] = fetch_1based(chrom, start, end).upper()
        parsed["alt"] = "-"
        return parsed
    return parsed


def to_mvinput_line(parsed: Dict[str, object], keep_chr: bool = False) -> str:
    """MATCHVAR 5-column line: chrom start end ref alt."""
    chrom = str(parsed["chrom"])
    if not keep_chr and chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    ref = str(parsed.get("ref") or "-") or "-"
    alt = str(parsed.get("alt") or "-") or "-"
    return f"{chrom}\t{parsed['start']}\t{parsed['end']}\t{ref}\t{alt}"


def pick_gene_field(row, field: str, preferred_db: str = "refGene") -> str:
    """
    从结果/Series 中按 MATCHVAR 列名取字段；兼容历史旧列名。
    field: function | gene | c_hgvs | exonic_effect | vartype | p_hgvs
    """
    dbs = [preferred_db] + [d for d in GENE_DBTYPES if d != preferred_db]
    key_map = {
        "function": ("Function", "Func"),
        "gene": ("Gene", "Gene"),
        "c_hgvs": ("cHGVS", "GeneDetail"),
        "mane_select": ("MANE Select", "MANE Select"),
        "exonic_effect": ("ExonicEffect", "ExonicFunc"),
        "vartype": ("VarType", "VarType"),
        "p_hgvs": ("pHGVS", "AAChange"),
    }
    new_prefix, old_prefix = key_map[field]
    for db in dbs:
        for prefix in (new_prefix, old_prefix):
            for sep in (".", "_"):
                col = f"{prefix}{sep}{db}"
                try:
                    val = row.get(col) if hasattr(row, "get") else row[col]
                except Exception:
                    val = None
                if val is not None and str(val) not in ("", "nan", "None", "NA", "N/A", "."):
                    return str(val)
    return ""
