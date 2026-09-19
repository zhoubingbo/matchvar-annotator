#!/usr/bin/env python3
"""从 dbNSFP 原始表生成 MATCHVAR filter 库。

原始文件示例：``dbNSFP5.3a.gz``（约 50GB，TSV，首行列名以 ``#chr`` 开头）。
输出：``hg19_dbNSFP53a.txt``（或 ``hg19_dbnsfp531a_rankscore.txt``），格式为::

  #Chr  Start  End  Ref  Alt  [选定注释列...]

坐标默认取 ``hg19_chr`` / ``hg19_pos(1-based)``（``--buildver hg19``），
Ref/Alt 取 ``ref`` / ``alt``。SNV 的 Start=End=pos。
染色体无 ``chr`` 前缀时自动补上（与 ``vcf2clinvar_db.py`` 一致）。

流式读写，不把整表载入内存，适合大文件。

示例::

  # 仅导出 *_rankscore（推荐给 Fusion-ML / protocol dbnsfp531a_rankscore）
  python matchvar/dbnsfp2matchvar_db.py /data/dbNSFP5.3a.gz \\
      -o /path/to/humandb/hg19_dbnsfp531a_rankscore.txt \\
      --rankscore-only --buildver hg19

  # 指定列
  python matchvar/dbnsfp2matchvar_db.py dbNSFP5.3a.gz \\
      -o hg19_dbNSFP53a.txt \\
      --columns REVEL_rankscore,CADD_raw_rankscore,SIFT_converted_rankscore

  # 冒烟：用仓库内小样例
  python matchvar/dbnsfp2matchvar_db.py data/test_dbNSFP5.3a.txt \\
      -o /tmp/hg19_dbNSFP53a.txt --rankscore-only
"""

from __future__ import annotations

import argparse
import gzip
import logging
import sys
import time
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Sequence, TextIO, Tuple

from tqdm import tqdm

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

CORE_HEADER = ["#Chr", "Start", "End", "Ref", "Alt"]

# dbNSFP 列名（随版本略有差异；用 normalize 做模糊匹配）
CHR_CANDIDATES = {
    "hg19": ("hg19_chr",),
    "hg18": ("hg18_chr",),
    "hg38": ("#chr", "chr"),
}
POS_CANDIDATES = {
    "hg19": ("hg19_pos(1-based)", "hg19_pos"),
    "hg18": ("hg18_pos(1-based)", "hg18_pos"),
    "hg38": ("pos(1-based)", "pos"),
}
REF_CANDIDATES = ("ref",)
ALT_CANDIDATES = ("alt",)


def input_bytes_tell(handle: TextIO) -> Optional[Callable[[], int]]:
    """返回读取磁盘上已消耗字节数的函数（gzip 用压缩偏移，便于对文件大小做进度条）。"""
    raw = getattr(handle, "buffer", None)
    if raw is not None and getattr(raw, "fileobj", None) is not None:
        return raw.fileobj.tell
    if raw is not None and hasattr(raw, "tell"):
        return raw.tell
    if hasattr(handle, "tell"):
        return handle.tell
    return None


def open_text(path: Path, mode: str = "rt") -> TextIO:
    if "b" in mode:
        raise ValueError("仅支持文本模式")
    if str(path).endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8", errors="replace")  # type: ignore[return-value]
    return open(path, mode, encoding="utf-8", errors="replace")


def normalize_header_name(name: str) -> str:
    return name.strip().lstrip("\ufeff")


def build_header_index(header_fields: Sequence[str]) -> Dict[str, int]:
    return {normalize_header_name(name): i for i, name in enumerate(header_fields)}


def resolve_column(index: Dict[str, int], candidates: Sequence[str], label: str) -> str:
    for name in candidates:
        if name in index:
            return name
        # 兼容去掉括号的写法
        for key in index:
            if key.replace(" ", "") == name.replace(" ", ""):
                return key
    raise KeyError(f"找不到列 {label}，候选: {candidates}；实际列示例: {list(index)[:12]}")


def normalize_chrom(chrom: str, *, ensure_chr: bool = True, strip_chr: bool = False, mito: str = "M") -> str:
    value = str(chrom).strip()
    if not value or value == ".":
        return value
    bare = value[3:] if value.lower().startswith("chr") else value
    if bare.upper() in {"M", "MT", "MITO"}:
        bare = mito
    if strip_chr:
        return bare
    if ensure_chr and not bare.lower().startswith("chr"):
        return f"chr{bare}"
    return value if value.lower().startswith("chr") else bare


def select_annotation_columns(
    header_fields: Sequence[str],
    *,
    columns: Optional[Sequence[str]] = None,
    rankscore_only: bool = False,
    column_pattern_suffix: Optional[str] = None,
) -> List[str]:
    """选择输出在 Ref/Alt 之后的注释列."""
    names = [normalize_header_name(c) for c in header_fields]
    # 排除坐标相关列，避免重复输出
    skip = {
        "#chr",
        "chr",
        "pos(1-based)",
        "pos",
        "ref",
        "alt",
        "hg19_chr",
        "hg19_pos(1-based)",
        "hg19_pos",
        "hg18_chr",
        "hg18_pos(1-based)",
        "hg18_pos",
        "hg38_chr",
        "hg38_pos(1-based)",
    }

    if columns:
        missing = [c for c in columns if c not in names]
        if missing:
            raise KeyError(f"--columns 中以下列不存在: {missing[:20]}")
        return list(columns)

    if rankscore_only or column_pattern_suffix:
        suffix = column_pattern_suffix or "rankscore"
        picked = [
            c
            for c in names
            if c not in skip and (c.endswith(suffix) or c.endswith(f"_{suffix}"))
        ]
        if not picked:
            raise ValueError(f"未匹配到以 {suffix} 结尾的列")
        return picked

    # 默认：全部非坐标列（文件会很大；更推荐 --rankscore-only）
    return [c for c in names if c not in skip]


def convert(
    input_path: Path,
    outfile: Path,
    *,
    buildver: str = "hg19",
    columns: Optional[Sequence[str]] = None,
    rankscore_only: bool = False,
    ensure_chr: bool = True,
    strip_chr: bool = False,
    mito: str = "M",
    missing: str = ".",
    progress_every: int = 1_000_000,
    max_rows: Optional[int] = None,
) -> None:
    buildver = buildver.lower()
    outfile.parent.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    written = 0
    skipped = 0
    scanned = 0
    total_bytes = input_path.stat().st_size
    show_bar = progress_every > 0

    with open_text(input_path, "rt") as fin, open(outfile, "w", encoding="utf-8", buffering=1024 * 1024) as fout:
        header_line = fin.readline()
        if not header_line:
            raise ValueError(f"空文件: {input_path}")
        header_fields = header_line.rstrip("\n").split("\t")
        index = build_header_index(header_fields)

        chr_col = resolve_column(index, CHR_CANDIDATES.get(buildver, CHR_CANDIDATES["hg38"]), "chr")
        pos_col = resolve_column(index, POS_CANDIDATES.get(buildver, POS_CANDIDATES["hg38"]), "pos")
        ref_col = resolve_column(index, REF_CANDIDATES, "ref")
        alt_col = resolve_column(index, ALT_CANDIDATES, "alt")

        anno_cols = select_annotation_columns(
            header_fields,
            columns=columns,
            rankscore_only=rankscore_only,
        )
        anno_idx = [index[c] for c in anno_cols]
        i_chr, i_pos, i_ref, i_alt = index[chr_col], index[pos_col], index[ref_col], index[alt_col]

        out_header = CORE_HEADER + anno_cols
        fout.write("\t".join(out_header) + "\n")
        logger.info(
            "输入=%s | buildver=%s | 坐标列=%s/%s | 注释列=%d | 输出=%s",
            input_path,
            buildver,
            chr_col,
            pos_col,
            len(anno_cols),
            outfile,
        )
        logger.info("注释列示例: %s%s", anno_cols[:8], " ..." if len(anno_cols) > 8 else "")

        tell = input_bytes_tell(fin)
        seen_bytes = tell() if tell else 0
        bar = tqdm(
            total=total_bytes if tell else None,
            desc="dbNSFP",
            unit="B" if tell else "行",
            unit_scale=True,
            unit_divisor=1024,
            file=sys.stderr,
            mininterval=0.5,
            disable=not show_bar,
        )
        if tell and show_bar:
            bar.update(min(seen_bytes, total_bytes))

        try:
            for line in fin:
                if max_rows is not None and written >= max_rows:
                    break
                scanned += 1
                line = line.rstrip("\n")
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) <= max(i_chr, i_pos, i_ref, i_alt, *(anno_idx or [0])):
                    skipped += 1
                    continue

                chrom_raw = parts[i_chr].strip()
                pos_raw = parts[i_pos].strip()
                ref = parts[i_ref].strip()
                alt = parts[i_alt].strip()
                if chrom_raw in {"", "."} or pos_raw in {"", "."} or ref in {"", "."} or alt in {"", "."}:
                    skipped += 1
                    continue
                try:
                    pos_i = int(float(pos_raw))
                except ValueError:
                    skipped += 1
                    continue

                chrom = normalize_chrom(
                    chrom_raw, ensure_chr=ensure_chr, strip_chr=strip_chr, mito=mito
                )
                # dbNSFP 每行已是单等位基因 SNV/替换；Start=End=pos
                row = [chrom, str(pos_i), str(pos_i), ref, alt]
                for j in anno_idx:
                    val = parts[j].strip() if j < len(parts) else ""
                    row.append(val if val != "" else missing)
                fout.write("\t".join(row) + "\n")
                written += 1

                if show_bar and scanned % progress_every == 0:
                    if tell:
                        now = tell()
                        delta = now - seen_bytes
                        if delta > 0:
                            bar.update(delta)
                            seen_bytes = now
                    else:
                        bar.update(progress_every)
                    bar.set_postfix(written=written, skipped=skipped, refresh=False)
        finally:
            if show_bar and tell:
                now = tell()
                delta = now - seen_bytes
                if delta > 0:
                    bar.update(delta)
            if show_bar:
                bar.set_postfix(written=written, skipped=skipped, refresh=False)
                bar.close()

    elapsed = time.time() - t0
    logger.info(
        "完成: scanned=%d written=%d skipped=%d → %s | 耗时 %.1f 分钟",
        scanned,
        written,
        skipped,
        outfile,
        elapsed / 60.0,
    )


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    examples = (
        "示例:\n"
        "1) 50GB dbNSFP5.3a.gz → 仅 rankscore 库（推荐）:\n"
        "   python matchvar/dbnsfp2matchvar_db.py /data/dbNSFP5.3a.gz \\\n"
        "       -o /path/humandb/hg19_dbnsfp531a_rankscore.txt \\\n"
        "       --rankscore-only --buildver hg19\n\n"
        "2) 全量列（体积很大，慎用）:\n"
        "   python matchvar/dbnsfp2matchvar_db.py dbNSFP5.3a.gz -o hg19_dbNSFP53a.txt\n\n"
        "3) 指定列:\n"
        "   python matchvar/dbnsfp2matchvar_db.py dbNSFP5.3a.gz -o out.txt \\\n"
        "       --columns REVEL_rankscore,CADD_raw_rankscore,MetaSVM_rankscore\n\n"
        "4) 用仓库样例冒烟:\n"
        "   python matchvar/dbnsfp2matchvar_db.py data/test_dbNSFP5.3a.txt \\\n"
        "       -o /tmp/hg19_dbNSFP53a.txt --rankscore-only\n\n"
        "生成后可用于:\n"
        "   --protocol refGene,clinvar,dbnsfp531a_rankscore\n"
        "大库建议再执行 scripts/build_tabix_indexes.py 做 bgzip+tabix。\n"
    )
    parser = argparse.ArgumentParser(
        description="从 dbNSFP5.3a 生成 MATCHVAR 用 hg19_dbNSFP53a / dbnsfp531a 库",
        epilog=examples,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("input", type=str, help="dbNSFP 原始 TSV / .gz（如 dbNSFP5.3a.gz）")
    parser.add_argument(
        "-o",
        "--outfile",
        type=str,
        default="hg19_dbNSFP53a.txt",
        help="输出库路径（默认 hg19_dbNSFP53a.txt）",
    )
    parser.add_argument(
        "--buildver",
        type=str,
        default="hg19",
        choices=["hg19", "hg18", "hg38"],
        help="坐标版本：决定使用 hg19_chr/pos 还是 #chr/pos（默认 hg19）",
    )
    parser.add_argument(
        "--columns",
        type=str,
        default=None,
        help="要保留的注释列，逗号分隔（不含坐标列）",
    )
    parser.add_argument(
        "--rankscore-only",
        action="store_true",
        help="只导出列名以 rankscore 结尾的列（推荐）",
    )
    parser.add_argument(
        "--strip-chr-prefix",
        action="store_true",
        help="去掉 chr 前缀（默认会自动补上 chr）",
    )
    parser.add_argument("--mito", type=str, default="M", help="线粒体染色体名（默认 M）")
    parser.add_argument("--missing", type=str, default=".", help="空值填充（默认 .）")
    parser.add_argument(
        "--progress-every",
        type=int,
        default=20_000,
        help="进度条刷新间隔（行数，默认 20000；0=关闭进度条）",
    )
    parser.add_argument(
        "--max-rows",
        type=int,
        default=None,
        help="最多写出多少数据行（调试用）",
    )
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    input_path = Path(args.input).expanduser()
    if not input_path.exists():
        logger.error("输入不存在: %s", input_path)
        return 1
    outfile = Path(args.outfile).expanduser()
    columns = None
    if args.columns:
        columns = [c.strip() for c in args.columns.split(",") if c.strip()]

    try:
        convert(
            input_path,
            outfile,
            buildver=args.buildver,
            columns=columns,
            rankscore_only=bool(args.rankscore_only),
            ensure_chr=not bool(args.strip_chr_prefix),
            strip_chr=bool(args.strip_chr_prefix),
            mito=args.mito,
            missing=args.missing,
            progress_every=int(args.progress_every),
            max_rows=args.max_rows,
        )
    except Exception as exc:  # noqa: BLE001
        logger.error("转换失败: %s", exc)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
