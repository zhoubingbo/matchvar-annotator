# -*- coding: utf-8 -*-
"""Convert g.HGVS strings to MATCHVAR input and run table annotation."""
from __future__ import annotations

import argparse
import logging
import os
import sys
import tempfile
from typing import Iterable, List, Optional, Sequence, Union

logger = logging.getLogger(__name__)

try:
    from .column_names import fill_g_hgvs_alleles, parse_g_hgvs, to_mvinput_line
    from .resource_files import GenomeSequence, discover_gene_protocols
    from .table_matchvar import TableAnnotator
except ImportError:  # pragma: no cover - script entry
    from column_names import fill_g_hgvs_alleles, parse_g_hgvs, to_mvinput_line
    from resource_files import GenomeSequence, discover_gene_protocols
    from table_matchvar import TableAnnotator


def read_ghgvs_lines(source: Union[str, Iterable[str]]) -> List[str]:
    """Read g.HGVS strings from a list, a file path, or '-' (stdin)."""
    if isinstance(source, str):
        if source == "-":
            lines = sys.stdin.read().splitlines()
        elif os.path.isfile(source):
            with open(source, "r", encoding="utf-8") as fh:
                lines = fh.readlines()
        else:
            # Treat a single HGVS token as input
            return [source.strip()] if source.strip() else []
    else:
        lines = list(source)
    out: List[str] = []
    for raw in lines:
        line = (raw or "").strip()
        if not line or line.startswith("#"):
            continue
        # allow TSV: first column is g.HGVS
        token = line.split("\t", 1)[0].strip()
        if token:
            out.append(token)
    return out


def ghgvs_to_mvinput_lines(
    hgvs_list: Sequence[str],
    genome: Optional[GenomeSequence] = None,
    keep_chr: bool = False,
) -> List[str]:
    lines: List[str] = []
    fetch = None
    if genome is not None and genome.ready:
        fetch = genome.fetch_1based
    for hgvs in hgvs_list:
        parsed = parse_g_hgvs(hgvs)
        ref = str(parsed.get("ref") or "")
        alt = str(parsed.get("alt") or "")
        needs_seq = (ref in ("", "-") and alt == "-") or (
            "dup" in str(parsed.get("g_hgvs") or "").lower() and alt in ("", "-")
        )
        if needs_seq:
            if fetch is None:
                raise ValueError(
                    f"{hgvs} omits the allele sequence; provide a genome FASTA or 2bit"
                )
            parsed = fill_g_hgvs_alleles(parsed, fetch)
        lines.append(to_mvinput_line(parsed, keep_chr=keep_chr))
    return lines


def write_mvinput(lines: Sequence[str], path: str) -> str:
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        for line in lines:
            fh.write(line.rstrip("\n") + "\n")
    return path


def annotate_g_hgvs(
    hgvs: Union[str, Sequence[str]],
    dbloc: str,
    buildver: str = "hg19",
    outfile: Optional[str] = None,
    protocol: Optional[str] = None,
    operation: Optional[str] = None,
    polish: bool = True,
    use_mane_transcript: bool = False,
    mane_file: Optional[str] = None,
    thread: Optional[int] = None,
    keep_temp: bool = False,
):
    """Convert g.HGVS → .mvinput → annotation table.

    Returns a pandas DataFrame when pandas can read the TSV; otherwise the
    output path string.
    """
    if isinstance(hgvs, str):
        hgvs_list = read_ghgvs_lines(hgvs)
    else:
        hgvs_list = [h.strip() for h in hgvs if (h or "").strip()]
    if not hgvs_list:
        raise ValueError("no g.HGVS strings to annotate")

    dbloc = os.path.abspath(dbloc)
    genome = None
    try:
        genome = GenomeSequence.from_dbloc(dbloc, buildver)
    except FileNotFoundError:
        genome = None

    mv_lines = ghgvs_to_mvinput_lines(hgvs_list, genome=genome, keep_chr=False)

    if outfile:
        out_prefix = os.path.abspath(outfile)
        workdir = os.path.dirname(out_prefix) or os.getcwd()
    else:
        workdir = tempfile.mkdtemp(prefix="matchvar_ghgvs_")
        out_prefix = os.path.join(workdir, "annot")

    mv_path = out_prefix + ".mvinput"
    write_mvinput(mv_lines, mv_path)

    if not protocol:
        discovered = discover_gene_protocols(dbloc, buildver)
        if not discovered:
            raise FileNotFoundError(
                f"no gene/region databases found in {dbloc} "
                "(need refGene.txt, ncbiRefSeq.bb, gencode*.bb, or cytoBand)"
            )
        protocol = ",".join(p for p, _ in discovered)
        operation = ",".join(op for _, op in discovered)
    elif not operation:
        ops = []
        for p in protocol.split(","):
            p = p.strip()
            if p in ("refGene", "ensGene", "knownGene", "ncbiRefSeq", "gencode"):
                ops.append("g")
            elif p == "cytoBand":
                ops.append("r")
            else:
                ops.append("f")
        operation = ",".join(ops)

    annotator = TableAnnotator(
        queryfile=mv_path,
        dbloc=dbloc,
        outfile=out_prefix,
        buildver=buildver,
        protocol=protocol,
        operation=operation,
        nopolish=not polish,
        polish=polish,
        use_mane_transcript=use_mane_transcript,
        mane_file=mane_file,
        thread=thread,
        remove=not keep_temp,
    )
    annotator.run_annotation()
    tsv = f"{out_prefix}.{buildver}_multianno.tsv"
    if genome is not None:
        genome.close()
    try:
        import pandas as pd

        return pd.read_csv(tsv, sep="\t", dtype=str)
    except Exception:
        return tsv


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Convert g.HGVS to a MATCHVAR annotation table",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  matchvar-ghgvs variants.txt --resources-dir resources --output result
  matchvar-ghgvs --hgvs 'chr5:g.131008083_131008084insTGACAGTTGTTTGCACAG' \\
      --humandb resources/humandb --output fnip1
  matchvar-ghgvs variants.ghgvs -p refGene,cytoBand -op g,r --buildver hg19
  matchvar-ghgvs --hgvs 'chr10:g.114925413_114925426del' \\
      --humandb resources/humandb/hg19 --use-mane-transcript -o tcf7l2
""",
    )
    parser.add_argument(
        "input",
        nargs="?",
        help="File of g.HGVS (one per line), a single g.HGVS string, or - for stdin",
    )
    parser.add_argument(
        "--hgvs",
        action="append",
        default=[],
        help="g.HGVS string (repeatable). Used in addition to INPUT",
    )
    parser.add_argument(
        "--humandb",
        help="humandb directory (genePred / bigBed / 2bit / cytoBand)",
    )
    parser.add_argument(
        "--resources-dir", "-r",
        help="resources directory that contains humandb/",
    )
    parser.add_argument(
        "--buildver", "--genome-version", "-g",
        dest="buildver",
        default="hg19",
        choices=["hg19", "hg38"],
    )
    parser.add_argument("--output", "-o", default="ghgvs_annot", help="output prefix")
    parser.add_argument("--protocol", "-p", help="comma-separated protocols (default: auto-detect)")
    parser.add_argument("--operation", "-op", help="comma-separated operations matching --protocol")
    parser.add_argument("--threads", "-t", type=int, default=None)
    parser.add_argument(
        "--use-mane-transcript",
        action="store_true",
        help=(
            "Keep only MANE Select / MANE Plus Clinical transcripts in "
            "cHGVS/pHGVS (auto-loads mane.bb, MANE GTF, or mane_transcript.txt)"
        ),
    )
    parser.add_argument(
        "--mane-file",
        help="MANE mapping file (GTF, mane.bb, or gene/refseq/ensembl TSV). Implies --use-mane-transcript",
    )
    parser.add_argument("--nopolish", action="store_true")
    parser.add_argument("--keep-temp", action="store_true")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )

    items: List[str] = list(args.hgvs or [])
    if args.input:
        items.extend(read_ghgvs_lines(args.input))
    if not items:
        parser.error("provide INPUT or --hgvs")

    if args.humandb:
        dbloc = os.path.abspath(args.humandb)
    elif args.resources_dir:
        dbloc = os.path.join(os.path.abspath(args.resources_dir), "humandb")
    else:
        env = os.environ.get("MATCHVAR_RESOURCES_DIR")
        if env:
            dbloc = os.path.join(os.path.abspath(env), "humandb")
        else:
            here = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "resources", "humandb"))
            dbloc = here if os.path.isdir(here) else os.path.abspath("resources/humandb")
    if not os.path.isdir(dbloc):
        logger.error("humandb directory does not exist: %s", dbloc)
        return 1

    try:
        result = annotate_g_hgvs(
            items,
            dbloc=dbloc,
            buildver=args.buildver,
            outfile=args.output,
            protocol=args.protocol,
            operation=args.operation,
            polish=not args.nopolish,
            use_mane_transcript=bool(args.use_mane_transcript or args.mane_file),
            mane_file=args.mane_file,
            thread=args.threads,
            keep_temp=args.keep_temp,
        )
    except Exception as exc:
        logger.error("%s", exc)
        if args.verbose:
            raise
        return 1

    tsv = f"{os.path.abspath(args.output)}.{args.buildver}_multianno.tsv"
    n = len(result) if hasattr(result, "__len__") and not isinstance(result, str) else "?"
    logger.info("annotation written: %s (%s rows)", tsv, n)
    return 0


if __name__ == "__main__":
    sys.exit(main())
