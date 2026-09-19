#!/usr/bin/env python3
"""从 ClinVar（或同类）VCF 生成 MATCHVAR filter 库文件。

输出格式与 humandb 中 ``{buildver}_clinvar.txt`` 一致：
  #Chr  Start  End  Ref  Alt  [INFO字段...]

前 5 列的坐标 / 等位基因归一化逻辑对齐 ``matchvar/convert2matchvar.py``
（left-normalize + adjustStartEndRefAlt），以便与 MATCHVAR 注释键一致。

示例::

  python matchvar/vcf2clinvar_db.py data/clinvar.vcf \\
      -o resources/humandb/hg19_clinvar.txt \\
      --info-fields CLNSIG,CLNREVSTAT,CLNALLELEID,GENEINFO

  # 注释时：
  #   --protocol refGene,clinvar,dbnsfp531a_rankscore
  #   humandb 目录下放置 hg19_clinvar.txt（或 hg19_clinvar.txt.gz + .tbi）
"""

from __future__ import annotations

import argparse
import gzip
import logging
import re
import sys
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, TextIO, Tuple

from tqdm import tqdm

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

CORE_HEADER = ["#Chr", "Start", "End", "Ref", "Alt"]
DEFAULT_INFO_FIELDS = ["CLNSIG", "CLNREVSTAT", "CLNALLELEID", "CLNDN", "GENEINFO"]


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


def open_text(path: Path, mode: str = "rt"):
    """按后缀打开普通或 gzip 文本."""
    if "b" in mode:
        raise ValueError("仅支持文本模式")
    if str(path).endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8", errors="replace")
    return open(path, mode, encoding="utf-8", errors="replace")


def adjust_start_end_ref_alt(
    start: int, end: int, ref: str, alt: str
) -> Tuple[int, int, str, str]:
    """与 convert2matchvar.adjust_start_end_ref_alt 一致."""
    _start, _end, _ref, _alt = start, end, ref, alt
    while len(_ref) > 0 and len(_alt) > 0 and _ref[-1] == _alt[-1]:
        _ref = _ref[:-1]
        _alt = _alt[:-1]
        _end -= 1
        if not _ref:
            _ref = "-"
            _start -= 1
            break
        if not _alt:
            _alt = "-"
            break

    while len(_ref) > 0 and len(_alt) > 0 and _ref[0] == _alt[0]:
        _ref = _ref[1:]
        _alt = _alt[1:]
        _start += 1
        if not _ref:
            _ref = "-"
            _start -= 1
            break
        if not _alt:
            _alt = "-"
            break

    if _alt == "-" and _ref != "-":
        _end = _start + len(_ref) - 1
    elif _ref == "-" and _alt != "-":
        _end = _start
    elif _ref != "-" and _alt != "-":
        _end = _start + len(_ref) - 1
    return _start, _end, _ref, _alt


def left_normalize(pos: int, ref: str, alt: str) -> Tuple[int, int, str, str]:
    """与 convert2matchvar.left_normalize 一致：VCF → MATCHVAR Chr/Start/End/Ref/Alt."""
    start = int(pos)
    if len(ref) == 1 and len(alt) == 1:
        return start, start, ref, alt

    if len(ref) < len(alt):
        head = alt[: len(ref)]
        if head == ref:
            new_start = start + len(ref) - 1
            new_end = new_start
            return adjust_start_end_ref_alt(new_start, new_end, "-", alt[len(ref) :])
        return adjust_start_end_ref_alt(start, start + len(ref) - 1, ref, alt)

    if len(ref) > len(alt):
        head = ref[: len(alt)]
        if head == alt:
            new_start = start + len(alt)
            new_end = start + len(ref) - 1
            return adjust_start_end_ref_alt(new_start, new_end, ref[len(alt) :], "-")
        return adjust_start_end_ref_alt(start, start + len(ref) - 1, ref, alt)

    head = ref[: len(ref) - 1]
    if alt.startswith(head):
        new_start = start + len(ref) - 1
        new_ref = ref[-1] if ref else "-"
        new_alt = alt[-1] if alt else "-"
        return adjust_start_end_ref_alt(new_start, new_start, new_ref, new_alt)
    return adjust_start_end_ref_alt(start, start + len(ref) - 1, ref, alt)


def normalize_chrom(
    chrom: str,
    *,
    ensure_chr: bool = True,
    strip_chr: bool = False,
    mito: str = "M",
) -> str:
    """规范化染色体名。

    - ``ensure_chr=True``（默认）：没有 ``chr`` 前缀时自动补上（如 ``1`` → ``chr1``）。
    - ``strip_chr=True``：去掉 ``chr`` 前缀（与 MATCHVAR 内部规范化一致）。
    """
    value = str(chrom).strip()
    if not value:
        return value

    bare = value[3:] if value.lower().startswith("chr") else value
    if bare.upper() in {"M", "MT", "MITO"}:
        bare = mito

    if strip_chr:
        return bare
    if ensure_chr and not bare.lower().startswith("chr"):
        return f"chr{bare}"
    return value if value.lower().startswith("chr") else bare

def parse_info(info: str) -> Dict[str, str]:
    """解析 VCF INFO：KEY=VALUE;KEY2=VALUE2。无 '=' 的 token 记为 _BARE."""
    result: Dict[str, str] = {}
    if not info or info == ".":
        return result
    for token in info.split(";"):
        token = token.strip()
        if not token:
            continue
        if "=" in token:
            key, value = token.split("=", 1)
            result[key] = value
        else:
            # 兼容非标准 ClinVar 子集：INFO 整段为 Uncertain_significance
            result.setdefault("_BARE", token)
            result[token] = "1"
    return result


def info_value_for_allele(
    info_map: Dict[str, str],
    field: str,
    allele_index: int,
    n_alts: int,
) -> str:
    """取 INFO 字段；若为逗号分隔且长度=等位基因数，取对应等位基因的值."""
    if field == "_BARE":
        return info_map.get("_BARE", ".")
    raw = info_map.get(field)
    if raw is None or raw == "":
        return "."
    if "," in raw and n_alts > 1:
        parts = raw.split(",")
        if len(parts) == n_alts:
            return parts[allele_index] if allele_index < len(parts) else "."
    return raw


def _tick_bar(bar, tell, seen_bytes: int, written: int, skipped: int) -> int:
    if tell is not None:
        now = tell()
        delta = now - seen_bytes
        if delta > 0:
            bar.update(delta)
            seen_bytes = now
    else:
        bar.update(1)
    bar.set_postfix(written=written, skipped=skipped, refresh=False)
    return seen_bytes


def convert(
    vcf_path: Path,
    outfile: Path,
    info_fields: Sequence[str],
    *,
    strip_chr: bool = False,
    ensure_chr: bool = True,
    mito: str = "M",
    keep_indel_ref: bool = False,
    include_id: bool = False,
    id_column: str = "CLNALLELEID",
    missing: str = ".",
    pass_only: bool = False,
) -> None:
    header = list(CORE_HEADER)
    if include_id:
        header.append(id_column)
    header.extend(list(info_fields))

    outfile.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    skipped = 0
    scanned = 0
    refresh_every = 2_000
    total_bytes = vcf_path.stat().st_size

    with open_text(vcf_path, "rt") as fin, open(outfile, "w", encoding="utf-8") as fout:
        fout.write("\t".join(header) + "\n")
        tell = input_bytes_tell(fin)
        seen_bytes = tell() if tell else 0
        bar = tqdm(
            total=total_bytes if tell else None,
            desc="ClinVar",
            unit="B" if tell else "行",
            unit_scale=True,
            unit_divisor=1024,
            file=sys.stderr,
            mininterval=0.5,
        )
        if tell:
            bar.update(min(seen_bytes, total_bytes))
        try:
            for line_no, line in enumerate(fin, 1):
                scanned += 1
                line = line.rstrip("\n")
                if not line or line.startswith("#"):
                    if scanned % refresh_every == 0:
                        seen_bytes = _tick_bar(bar, tell, seen_bytes, written, skipped)
                    continue
                cols = line.split("\t")
                if len(cols) < 8:
                    skipped += 1
                    continue
                chrom, pos, vid, ref, alt, _qual, filt, info = cols[:8]
                if pass_only and filt not in ("PASS", ".", ""):
                    skipped += 1
                    continue

                try:
                    pos_i = int(pos)
                except ValueError:
                    logger.warning("第 %d 行 POS 非法: %s", line_no, pos)
                    skipped += 1
                    continue

                chrom_out = normalize_chrom(
                    chrom, ensure_chr=ensure_chr, strip_chr=strip_chr, mito=mito
                )
                info_map = parse_info(info)
                alts = [a for a in alt.split(",") if a]
                if not alts:
                    skipped += 1
                    continue

                for i_alt, alt_i in enumerate(alts):
                    if alt_i.upper() == "<DEL>":
                        end_m = re.search(r"\bEND=(\d+)\b", info)
                        if end_m:
                            new_start, new_end, new_ref, new_alt = (
                                pos_i,
                                int(end_m.group(1)),
                                "0",
                                "-",
                            )
                        else:
                            new_start, new_end = pos_i, pos_i + len(ref) - 1
                            new_ref, new_alt = ref, "-"
                    elif alt_i.upper() in {"<DUP>", "<INV>", "<INS>"}:
                        end_m = re.search(r"\bEND=(\d+)\b", info)
                        if end_m:
                            new_start, new_end, new_ref, new_alt = (
                                pos_i,
                                int(end_m.group(1)),
                                "0",
                                "0",
                            )
                        else:
                            new_start, new_end = pos_i, pos_i + len(ref) - 1
                            new_ref, new_alt = ref, "-"
                    elif keep_indel_ref:
                        new_start, new_end = pos_i, pos_i + len(ref) - 1
                        new_ref, new_alt = ref, alt_i
                    else:
                        new_start, new_end, new_ref, new_alt = left_normalize(pos_i, ref, alt_i)

                    row = [
                        chrom_out,
                        str(new_start),
                        str(new_end),
                        new_ref,
                        new_alt,
                    ]
                    if include_id:
                        allele_id = (
                            info_map.get("ALLELEID")
                            or info_map.get("CLNALLELEID")
                            or (vid if vid not in {".", ""} else missing)
                        )
                        if "," in str(allele_id) and len(alts) > 1:
                            parts = str(allele_id).split(",")
                            allele_id = parts[i_alt] if i_alt < len(parts) else missing
                        row.append(str(allele_id) if allele_id not in {None, ""} else missing)

                    for field in info_fields:
                        val = info_value_for_allele(info_map, field, i_alt, len(alts))
                        # 非标准 VCF：INFO 整段为意义字符串且请求 CLNSIG
                        if (val == "." or val == "") and field == "CLNSIG" and "_BARE" in info_map:
                            val = info_map["_BARE"]
                        row.append(val if val not in {None, ""} else missing)

                    fout.write("\t".join(row) + "\n")
                    written += 1

                if scanned % refresh_every == 0:
                    seen_bytes = _tick_bar(bar, tell, seen_bytes, written, skipped)
        finally:
            if tell:
                now = tell()
                delta = now - seen_bytes
                if delta > 0:
                    bar.update(delta)
            bar.set_postfix(written=written, skipped=skipped, refresh=False)
            bar.close()

    logger.info(
        "完成: 写入 %d 行 → %s（跳过 %d 行）| INFO 列: %s",
        written,
        outfile,
        skipped,
        ",".join(info_fields) if info_fields else "(无)",
    )


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    examples = (
        "示例:\n"
        "1) 标准 ClinVar VCF → hg19_clinvar.txt（默认提取 CLNSIG 等）:\n"
        "   python matchvar/vcf2clinvar_db.py clinvar.vcf \\\n"
        "       -o resources/humandb/hg19_clinvar.txt\n\n"
        "2) 只提取 CLNSIG / CLNREVSTAT:\n"
        "   python matchvar/vcf2clinvar_db.py clinvar.vcf -o hg19_clinvar.txt \\\n"
        "       --info-fields CLNSIG,CLNREVSTAT\n\n"
        "3) 去掉 chr 前缀，并写入 ID 列:\n"
        "   python matchvar/vcf2clinvar_db.py clinvar.vcf -o hg19_clinvar.txt \\\n"
        "       --strip-chr-prefix --include-id\n\n"
        "生成后用于 MATCHVAR:\n"
        "   python run_pipeline.py --annotate-only --variant-input query.vcf \\\n"
        "       --humandb resources/humandb --protocol refGene,clinvar --buildver hg19\n"
    )
    parser = argparse.ArgumentParser(
        description="从 ClinVar VCF 生成 MATCHVAR 用 clinvar 库（hg19_clinvar.txt）",
        epilog=examples,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "vcf",
        type=str,
        help="输入 ClinVar VCF / VCF.gz",
    )
    parser.add_argument(
        "-o",
        "--outfile",
        type=str,
        default="hg19_clinvar.txt",
        help="输出库文件路径（默认: hg19_clinvar.txt）",
    )
    parser.add_argument(
        "--info-fields",
        type=str,
        default=",".join(DEFAULT_INFO_FIELDS),
        help=(
            "从 INFO 提取的字段，逗号分隔；"
            f"默认: {','.join(DEFAULT_INFO_FIELDS)}。"
            "对应 INFO 中 KEY=value（如 CLNSIG=Likely_benign）"
        ),
    )
    parser.add_argument(
        "--buildver",
        type=str,
        default=None,
        help="若指定且 -o 为默认名，则输出为 {buildver}_clinvar.txt",
    )
    parser.add_argument(
        "--strip-chr-prefix",
        action="store_true",
        help="去掉染色体 chr 前缀（默认会补全为 chr1 / chrM 等形式）",
    )
    parser.add_argument(
        "--keep-chr-prefix",
        action="store_true",
        help="兼容旧参数：与默认行为相同（无 chr 时自动补上）",
    )
    parser.add_argument(
        "--mito",
        type=str,
        default="M",
        help="线粒体染色体输出名（默认 M）",
    )
    parser.add_argument(
        "--keep-indel-ref",
        action="store_true",
        help="不做 indel left-normalize（对齐 convert2matchvar -keepindelref）",
    )
    parser.add_argument(
        "--include-id",
        action="store_true",
        help="在 Ref/Alt 后增加 ID 列（优先 INFO.ALLELEID / CLNALLELEID，否则用 VCF ID）",
    )
    parser.add_argument(
        "--id-column",
        type=str,
        default="CLNALLELEID",
        help="--include-id 时的列名（默认 CLNALLELEID）",
    )
    parser.add_argument(
        "--missing",
        type=str,
        default=".",
        help="缺失 INFO 填充符（默认 .）",
    )
    parser.add_argument(
        "--pass-only",
        action="store_true",
        help="仅保留 FILTER 为 PASS 或 . 的记录",
    )
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    vcf_path = Path(args.vcf).expanduser()
    if not vcf_path.exists():
        logger.error("VCF 不存在: %s", vcf_path)
        return 1

    outfile = Path(args.outfile).expanduser()
    if args.buildver and args.outfile == "hg19_clinvar.txt":
        outfile = Path(f"{args.buildver}_clinvar.txt")

    info_fields = [f.strip() for f in str(args.info_fields).split(",") if f.strip()]
    strip_chr = bool(args.strip_chr_prefix)
    ensure_chr = not strip_chr  # 默认补全 chr；仅 --strip-chr-prefix 时不补

    try:
        convert(
            vcf_path,
            outfile,
            info_fields,
            strip_chr=strip_chr,
            ensure_chr=ensure_chr,
            mito=args.mito,
            keep_indel_ref=args.keep_indel_ref,
            include_id=args.include_id,
            id_column=args.id_column,
            missing=args.missing,
            pass_only=args.pass_only,
        )
    except Exception as exc:  # noqa: BLE001
        logger.error("转换失败: %s", exc)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
