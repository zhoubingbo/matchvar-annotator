# -*- coding: utf-8 -*-
"""MANE 转录本映射加载与匹配（供 table_matchvar / annotate_variation 共用）。

优先读取 resources/mane 下最新的 MANE *.gtf / *.gtf.gz；
键为基因符号（gene），同时保留 gene_id 别名。
同一基因可包含 MANE Select + MANE Plus Clinical 多条转录本。
"""

from __future__ import annotations

import gzip
import logging
import os
import re
from typing import Any, Dict, Iterable, Optional, Set, TextIO

logger = logging.getLogger(__name__)

_ATTR_RE = {
    "gene_id": re.compile(r'gene_id\s+"([^"]+)"'),
    "gene": re.compile(r'(?:^|;)\s*gene\s+"([^"]+)"'),
    "transcript_id": re.compile(r'transcript_id\s+"([^"]+)"'),
    "ensembl": re.compile(r'db_xref\s+"Ensembl:(ENST[^"]+)"'),
    "tag": re.compile(r'tag\s+"([^"]+)"'),
}


def project_resources_dir() -> str:
    env = os.environ.get("MATCHVAR_RESOURCES_DIR")
    if env:
        return os.path.abspath(env)
    pkg_dir = os.path.dirname(os.path.abspath(__file__))
    candidates = [
        os.path.join(os.path.dirname(pkg_dir), "resources"),
        os.path.join(os.getcwd(), "resources"),
    ]
    for path in candidates:
        if os.path.isdir(path):
            return path
    return candidates[0]


def resolve_mane_file(
    explicit: Optional[str] = None,
    resources_dir: Optional[str] = None,
) -> Optional[str]:
    """解析可用的 MANE 映射文件路径。"""
    if explicit and os.path.isfile(explicit):
        return os.path.abspath(explicit)

    resources = resources_dir or project_resources_dir()
    search_dirs = [
        os.path.join(resources, "geneinfo", "hg38"),
        os.path.join(resources, "mane", "hg38"),
        os.path.join(resources, "mane"),
    ]
    candidates = []
    for mane_dir in search_dirs:
        if not os.path.isdir(mane_dir):
            continue
        for name in os.listdir(mane_dir):
            lower = name.lower()
            if "mane" not in lower:
                continue
            if lower.endswith(".gtf") or lower.endswith(".gtf.gz"):
                full = os.path.join(mane_dir, name)
                if os.path.isfile(full) or os.path.islink(full):
                    candidates.append(full)
    if candidates:
        candidates.sort(key=lambda p: os.path.getmtime(p) if os.path.exists(p) else 0, reverse=True)
        return candidates[0]

    legacy = os.path.join(resources, "humandb", "mane_transcript.txt")
    if os.path.isfile(legacy):
        return legacy
    legacy2 = os.path.join(resources, "humandb", "hg19", "mane_transcript.txt")
    if os.path.isfile(legacy2):
        return legacy2
    for bb_name in ("mane.bb", "hg19_mane.bb", "hg38_mane.bb"):
        for folder in (
            os.path.join(resources, "humandb"),
            resources,
        ):
            cand = os.path.join(folder, bb_name)
            if os.path.isfile(cand):
                return cand
    return None


def _open_text(path: str) -> TextIO:
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def base_transcript_id(transcript_id: Optional[str]) -> Optional[str]:
    if not transcript_id:
        return None
    tid = str(transcript_id).strip()
    if not tid or tid in (".", "NA", "None"):
        return None
    return tid.split(".")[0] if "." in tid else tid


def _empty_entry() -> Dict[str, Any]:
    return {
        "refseq": "",
        "ensembl": "",
        "base_refseq": "",
        "base_ensembl": "",
        "tag": "",
        "base_ids": set(),
        "transcripts": [],
    }


def _add_transcript(
    store: Dict[str, Dict[str, Any]],
    gene_key: str,
    refseq: str,
    ensembl: Optional[str],
    tag: str,
) -> None:
    if not gene_key or not refseq:
        return
    entry = store.setdefault(gene_key, _empty_entry())
    base_ref = base_transcript_id(refseq)
    base_ens = base_transcript_id(ensembl)
    tx = {
        "refseq": refseq,
        "ensembl": ensembl or "",
        "base_refseq": base_ref or "",
        "base_ensembl": base_ens or "",
        "tag": tag or "",
    }
    # 去重
    for existing in entry["transcripts"]:
        if existing.get("base_refseq") == base_ref and existing.get("base_ensembl") == base_ens:
            return
    entry["transcripts"].append(tx)
    if base_ref:
        entry["base_ids"].add(base_ref)
    if base_ens:
        entry["base_ids"].add(base_ens)

    # 主记录优先 MANE Select
    is_select = "MANE Select" in (tag or "")
    current_is_select = "MANE Select" in (entry.get("tag") or "")
    if not entry["refseq"] or (is_select and not current_is_select):
        entry["refseq"] = refseq
        entry["ensembl"] = ensembl or ""
        entry["base_refseq"] = base_ref or ""
        entry["base_ensembl"] = base_ens or ""
        entry["tag"] = tag or ""


def _parse_gtf_line(line: str, store: Dict[str, Dict[str, Any]]) -> None:
    parts = line.split("\t")
    if len(parts) >= 9:
        feature = (parts[2] or "").strip()
        # 正式 GTF：仅 transcript / mRNA 行
        if feature not in ("transcript", "mRNA"):
            return
        attributes = parts[8]
    elif "transcript_id" in line:
        # 兼容仅含 attributes 的简化行
        attributes = line
    else:
        return

    gene_symbol = None
    gene_id = None
    transcript_id = None
    ensembl_id = None
    tag = ""

    m = _ATTR_RE["gene"].search(attributes)
    if m:
        gene_symbol = m.group(1)
    m = _ATTR_RE["gene_id"].search(attributes)
    if m:
        gene_id = m.group(1)
    m = _ATTR_RE["transcript_id"].search(attributes)
    if m:
        transcript_id = m.group(1)
    m = _ATTR_RE["ensembl"].search(attributes)
    if m:
        ensembl_id = m.group(1)
    m = _ATTR_RE["tag"].search(attributes)
    if m:
        tag = m.group(1)

    if not transcript_id:
        return
    # 只要 MANE Select / Plus Clinical（若无 tag 也收录，兼容非官方文件）
    if tag and ("MANE" not in tag):
        return

    primary_key = gene_symbol or gene_id
    if not primary_key:
        return
    _add_transcript(store, primary_key, transcript_id, ensembl_id, tag)
    # gene_id 与 gene 不同时建立别名，便于兼容旧逻辑
    if gene_id and gene_id != primary_key:
        # 共享同一条目，保证后续 Select 优先更新同步到别名
        store[gene_id] = store[primary_key]


def _parse_three_col_line(line: str, store: Dict[str, Dict[str, Any]]) -> None:
    parts = line.split("\t")
    if len(parts) < 2:
        return
    gene_id = parts[0].strip()
    refseq = parts[1].strip()
    ensembl = parts[2].strip() if len(parts) >= 3 else ""
    if gene_id and refseq:
        _add_transcript(store, gene_id, refseq, ensembl or None, "MANE")


_MANE_STATUS_RANK = {
    "MANE Select": 2,
    "MANE Plus Clinical": 1,
}


def _load_mane_from_bigbed(path: str) -> Dict[str, Dict[str, Any]]:
    """Parse UCSC mane.bb (bigGenePred extras: gene symbol, NM, status)."""
    try:
        import pyBigWig
    except ImportError:
        logger.warning("mane.bb found but pyBigWig is not installed; skip MANE mapping")
        return {}
    bw = pyBigWig.open(path)
    if bw is None:
        logger.warning("cannot open mane.bb: %s", path)
        return {}
    store: Dict[str, Dict[str, Any]] = {}
    try:
        chroms = bw.chroms() or {}
        for chrom, length in chroms.items():
            entries = bw.entries(chrom, 0, int(length), withString=True) or []
            for _start, _end, rest in entries:
                parts = (rest or "").split("\t")
                if len(parts) < 19:
                    continue
                status = parts[-1].strip()
                if status not in _MANE_STATUS_RANK:
                    continue
                gene = (parts[15] if len(parts) > 15 else "").strip() or (parts[9] if len(parts) > 9 else "").strip()
                ens = (parts[0] or "").strip()
                nm = (parts[18] if len(parts) > 18 else "").strip()
                if not nm.upper().startswith(("NM_", "NR_")):
                    continue
                tag = status
                if gene:
                    _add_transcript(store, gene, nm, ens or None, tag)
                if ens:
                    _add_transcript(store, ens.split(".")[0], nm, ens, tag)
    finally:
        bw.close()
    logger.info("已从 mane.bb 加载 MANE 映射: file=%s genes=%s", path, len(store))
    return store


def load_mane_transcripts(
    mane_file: Optional[str] = None,
    resources_dir: Optional[str] = None,
) -> Dict[str, Dict[str, Any]]:
    """加载 gene → MANE 转录本映射。"""
    path = resolve_mane_file(mane_file, resources_dir)
    if not path:
        logger.warning("未找到 MANE 转录本文件（resources/mane、humandb/mane_transcript.txt 或 mane.bb）")
        return {}

    if path.lower().endswith(".bb"):
        return _load_mane_from_bigbed(path)

    store: Dict[str, Dict[str, Any]] = {}
    try:
        with _open_text(path) as fh:
            for raw in fh:
                line = raw.strip()
                if not line or line.startswith("#"):
                    continue
                # GTF（含 gzip）或含 transcript_id 的 attributes 行
                if "\t" in line and (
                    line.count("\t") >= 8 or "transcript_id" in line
                ):
                    cols = line.split("\t")
                    if len(cols) >= 9:
                        _parse_gtf_line(line, store)
                    elif len(cols) >= 2 and not cols[0].startswith("chr"):
                        _parse_three_col_line(line, store)
                    else:
                        _parse_gtf_line(line, store)
                elif "\t" in line:
                    _parse_three_col_line(line, store)
                elif "transcript_id" in line:
                    _parse_gtf_line(f".\t.\ttranscript\t.\t.\t.\t.\t.\t{line}", store)

        n_primary = sum(
            1 for gene, info in store.items()
            if info.get("refseq") and info.get("transcripts")
            and any(tx.get("refseq") == info.get("refseq") for tx in info.get("transcripts") or [])
        )
        # 更稳妥：以拥有 transcripts 且 gene 出现在自身 transcripts 所属键中计数
        n_with_tx = sum(1 for info in store.values() if info.get("transcripts"))
        logger.info(
            "已加载 MANE 转录本映射: file=%s entries=%s with_transcripts=%s",
            path, len(store), n_with_tx,
        )
        return store
    except Exception as exc:
        logger.error("加载 MANE 文件失败 (%s): %s", path, exc)
        return {}


def mane_base_ids(mane_info: Any) -> Set[str]:
    """从一条 gene 的 MANE 记录取出所有可匹配的 base transcript ID。"""
    ids: Set[str] = set()
    if not mane_info:
        return ids
    if isinstance(mane_info, str):
        bid = base_transcript_id(mane_info)
        if bid:
            ids.add(bid)
        return ids
    if isinstance(mane_info, dict):
        for key in ("base_refseq", "base_ensembl", "refseq", "ensembl"):
            bid = base_transcript_id(mane_info.get(key))
            if bid:
                ids.add(bid)
        for bid in mane_info.get("base_ids") or []:
            if bid:
                ids.add(str(bid))
        for tx in mane_info.get("transcripts") or []:
            if isinstance(tx, dict):
                for key in ("base_refseq", "base_ensembl", "refseq", "ensembl"):
                    bid = base_transcript_id(tx.get(key))
                    if bid:
                        ids.add(bid)
    return ids


def all_mane_base_ids(mane_transcripts: Dict[str, Dict[str, Any]]) -> Set[str]:
    """汇总全部 MANE 转录本 base ID（跨基因），用于多基因重叠位点过滤。"""
    ids: Set[str] = set()
    for info in (mane_transcripts or {}).values():
        ids |= mane_base_ids(info)
    return ids


def filter_parts_by_global_mane_ids(
    parts: Iterable[str],
    allowed_ids: Set[str],
    joiner: str = ";",
    fallback_to_original: bool = True,
) -> str:
    """仅保留转录本 ID 落在 allowed_ids 中的注释片段。

    无 MANE 匹配时：默认回退保留原片段（fallback_to_original=True）。
    """
    original = joiner.join(p.strip() for p in parts if (p or "").strip())
    if not allowed_ids:
        return original
    kept = []
    for part in parts:
        p = (part or "").strip()
        if not p:
            continue
        # 构造假 mane_info 以便复用 token 抽取
        fake = {"base_ids": allowed_ids}
        if is_mane_transcript_token(p, fake):
            kept.append(p)
    if kept:
        return joiner.join(kept)
    return original if fallback_to_original else ""


def is_mane_transcript_token(token: str, mane_info: Any) -> bool:
    """判断注释片段（如 NM_xxx:exon:c.xx 或 GENE:NM_xxx:...）是否属于 MANE 转录本。"""
    if not token or not mane_info:
        return False
    allowed = mane_base_ids(mane_info)
    if not allowed:
        return False
    # 取出片段中所有看起来像转录本 ID 的 token
    candidates = re.findall(r"\b((?:NM|NR|XM|XR|ENST)_[0-9]+(?:\.[0-9]+)?)\b", token)
    if not candidates and ":" in token:
        # 回退：取第一个冒号前
        head = token.split(":")[0].strip()
        # GENE:NM_xxx 形式取第二段
        if re.match(r"^[A-Za-z0-9_-]+$", head) and token.count(":") >= 1:
            parts = token.split(":")
            for p in parts[:3]:
                if re.match(r"^(?:NM|NR|XM|XR|ENST)_", p):
                    candidates.append(p)
                    break
        else:
            candidates.append(head)
    for c in candidates:
        bid = base_transcript_id(c)
        if bid and bid in allowed:
            return True
    return False


def filter_annotation_parts(
    parts: Iterable[str],
    mane_info: Any,
    joiner: str = ";",
    fallback_to_original: bool = True,
) -> str:
    """过滤注释片段列表，仅保留 MANE 转录本。

    无匹配时：默认回退保留原片段（fallback_to_original=True）。
    """
    original = joiner.join(p.strip() for p in parts if (p or "").strip())
    kept = []
    for part in parts:
        p = (part or "").strip()
        if not p:
            continue
        if is_mane_transcript_token(p, mane_info):
            kept.append(p)
    if kept:
        return joiner.join(kept)
    return original if fallback_to_original else ""


def write_mane_tsv(mane_transcripts: Dict[str, Dict[str, Any]], out_path: str) -> None:
    """写出 gene / refseq / ensembl 三列临时映射（每转录本一行）。"""
    with open(out_path, "w", encoding="utf-8") as fh:
        written = set()
        for gene, info in mane_transcripts.items():
            if not isinstance(info, dict):
                fh.write(f"{gene}\t{info}\t\n")
                continue
            txs = info.get("transcripts") or []
            if not txs and info.get("refseq"):
                txs = [{
                    "refseq": info.get("refseq", ""),
                    "ensembl": info.get("ensembl", ""),
                }]
            for tx in txs:
                refseq = (tx.get("refseq") or "").strip()
                ensembl = (tx.get("ensembl") or "").strip()
                key = (gene, refseq, ensembl)
                if not refseq or key in written:
                    continue
                written.add(key)
                fh.write(f"{gene}\t{refseq}\t{ensembl}\n")
