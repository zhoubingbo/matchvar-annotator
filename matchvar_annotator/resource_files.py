# -*- coding: utf-8 -*-
"""Discover and materialize humandb resources (genePred txt, bigBed, 2bit, FASTA).

Recognized extra files in ``humandb`` (in addition to ``{buildver}_refGene.txt``
and a genome FASTA):

- ``ncbiRefSeq.bb`` / ``gencodeV49lift37.bb`` / ``gencode*.bb``
- ``mane.bb``
- ``{buildver}.2bit`` / ``hg19.2bit``
- ``{buildver}.chrom.sizes``
- ``cytoBand.txt.gz``
"""
from __future__ import annotations

import gzip
import logging
import os
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple, Union

logger = logging.getLogger(__name__)

PathLike = Union[str, Path]

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
_PRIMARY_CHROM_RE = re.compile(r"^(?:chr)?(?:[0-9]{1,2}|[XYM]|MT)$", re.I)

GENE_PROTOCOLS = ("refGene", "ncbiRefSeq", "ensGene", "gencode", "knownGene")
GENEPRED_LIKE = {"refGene", "ncbiRefSeq", "ensGene", "gencode"}


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def chrom_aliases(chrom: str) -> List[str]:
    raw = (chrom or "").strip()
    if not raw:
        return []
    if raw.lower().startswith("chr"):
        rest = raw[3:]
        aliases = [raw, rest]
        if rest.upper() in ("M", "MT"):
            aliases.extend(["chrM", "chrMT", "M", "MT"])
        return list(dict.fromkeys(aliases))
    aliases = [raw, f"chr{raw}"]
    if raw.upper() in ("M", "MT"):
        aliases.extend(["chrM", "chrMT", "M", "MT"])
    return list(dict.fromkeys(aliases))


def strip_chr(chrom: str) -> str:
    c = (chrom or "").strip()
    if c.lower().startswith("chr"):
        c = c[3:]
    if c.upper() in ("M", "MT"):
        return "M"
    return c


def is_primary_chrom(chrom: str) -> bool:
    return bool(_PRIMARY_CHROM_RE.match((chrom or "").strip()))


def _first_existing(paths: Iterable[PathLike]) -> Optional[str]:
    for p in paths:
        if p and os.path.isfile(str(p)):
            return os.path.abspath(str(p))
    return None


def list_dir_files(dbloc: str) -> List[str]:
    try:
        return os.listdir(dbloc)
    except OSError:
        return []


def resolve_genome_2bit(dbloc: str, buildver: str = "hg19") -> Optional[str]:
    dbloc = os.path.abspath(dbloc)
    parent = os.path.dirname(dbloc)
    names = [
        f"{buildver}.2bit",
        "hg19.2bit" if buildver == "hg19" else "hg38.2bit",
        f"{buildver}.2bit",
    ]
    candidates: List[str] = []
    for name in names:
        candidates.append(os.path.join(dbloc, name))
        candidates.append(os.path.join(parent, name))
        candidates.append(os.path.join(parent, "reference", buildver, name))
    return _first_existing(candidates)


def resolve_genome_fasta(dbloc: str, buildver: str = "hg19") -> Optional[str]:
    dbloc = os.path.abspath(dbloc)
    parent = os.path.dirname(dbloc)
    names = [
        f"{buildver}.fa",
        f"{buildver}.fasta",
        "Homo_sapiens_assembly37.fasta" if buildver == "hg19" else "Homo_sapiens_assembly38.fasta",
        "hg19.fa" if buildver == "hg19" else "hg38.fa",
        "GRCh37.fa" if buildver == "hg19" else "GRCh38.fa",
    ]
    candidates: List[str] = []
    for name in names:
        candidates.append(os.path.join(dbloc, name))
        candidates.append(os.path.join(parent, name))
    return _first_existing(candidates)


def resolve_chrom_sizes(dbloc: str, buildver: str = "hg19") -> Optional[str]:
    dbloc = os.path.abspath(dbloc)
    parent = os.path.dirname(dbloc)
    names = [f"{buildver}.chrom.sizes", "hg19.chrom.sizes" if buildver == "hg19" else "hg38.chrom.sizes"]
    candidates: List[str] = []
    for name in names:
        candidates.append(os.path.join(dbloc, name))
        candidates.append(os.path.join(parent, name))
    return _first_existing(candidates)


def resolve_cytoband(dbloc: str, buildver: str = "hg19") -> Optional[str]:
    dbloc = os.path.abspath(dbloc)
    return _first_existing(
        [
            os.path.join(dbloc, f"{buildver}_cytoBand.txt"),
            os.path.join(dbloc, f"{buildver}_cytoBand.txt.gz"),
            os.path.join(dbloc, "cytoBand.txt"),
            os.path.join(dbloc, "cytoBand.txt.gz"),
        ]
    )


def resolve_mane_bb(dbloc: str, buildver: str = "hg19") -> Optional[str]:
    dbloc = os.path.abspath(dbloc)
    parent = os.path.dirname(dbloc)
    return _first_existing(
        [
            os.path.join(dbloc, "mane.bb"),
            os.path.join(dbloc, f"{buildver}_mane.bb"),
            os.path.join(dbloc, f"mane.{buildver}.bb"),
            os.path.join(parent, "mane.bb"),
        ]
    )


def _gencode_bb_names(buildver: str) -> List[str]:
    if buildver == "hg19":
        return ["gencodeV49lift37.bb", "gencodeV49.bb"]
    return ["gencodeV49.bb", "gencodeV49lift37.bb"]


def resolve_gene_bigbed(dbloc: str, protocol: str, buildver: str = "hg19") -> Optional[str]:
    """Return a gene-model bigBed for this protocol, if present."""
    dbloc = os.path.abspath(dbloc)
    names: List[str] = []
    proto = (protocol or "").strip()
    if proto in ("refGene", "ncbiRefSeq"):
        names.extend(["ncbiRefSeq.bb", f"{buildver}_ncbiRefSeq.bb", f"ncbiRefSeq.{buildver}.bb"])
    if proto in ("ensGene", "gencode"):
        names.extend(_gencode_bb_names(buildver))
        names.extend([f"{buildver}_gencode.bb", "gencode.bb"])
    found = _first_existing(os.path.join(dbloc, n) for n in names)
    if found:
        return found
    if proto in ("ensGene", "gencode") and os.path.isdir(dbloc):
        for fn in sorted(list_dir_files(dbloc)):
            if fn.lower().startswith("gencode") and fn.lower().endswith(".bb"):
                return os.path.join(dbloc, fn)
    return None


def resolve_gene_txt(dbloc: str, protocol: str, buildver: str = "hg19") -> Optional[str]:
    dbloc = os.path.abspath(dbloc)
    stems = [protocol, f"{buildver}_{protocol}"]
    if protocol == "refGene":
        stems.extend(["ncbiRefSeq", f"{buildver}_ncbiRefSeq"])
    if protocol == "ensGene":
        stems.extend(["gencode", f"{buildver}_gencode"])
    for stem in stems:
        for suf in (".txt", ".txt.gz"):
            cand = os.path.join(dbloc, f"{stem}{suf}")
            if os.path.isfile(cand):
                return cand
    return None


def resolve_mrna_fasta(dbloc: str, protocol: str, buildver: str = "hg19") -> Optional[str]:
    dbloc = os.path.abspath(dbloc)
    stems = [protocol, f"{buildver}_{protocol}"]
    if protocol == "refGene":
        stems.extend(["ncbiRefSeq", f"{buildver}_ncbiRefSeq"])
    if protocol == "ensGene":
        stems.extend(["gencode", f"{buildver}_gencode"])
    for stem in stems:
        for suf in (".fa", ".fasta"):
            cand = os.path.join(dbloc, f"{stem}Mrna{suf}")
            if os.path.isfile(cand):
                return cand
    return None


class GenomeSequence:
    """Fetch genomic DNA from FASTA (pyfaidx) or 2bit."""

    def __init__(self, fasta_path: Optional[str] = None, twobit_path: Optional[str] = None):
        self.fasta_path = fasta_path
        self.twobit_path = twobit_path
        self._fasta = None
        self._twobit = None
        if fasta_path:
            try:
                from pyfaidx import Fasta

                try:
                    self._fasta = Fasta(fasta_path, as_raw=True, sequence_always_upper=True)
                except TypeError:
                    self._fasta = Fasta(fasta_path, as_raw=True)
            except Exception as exc:
                logger.warning("failed to open FASTA %s: %s", fasta_path, exc)
                self._fasta = None
        if twobit_path and self._fasta is None:
            from .twobit_reader import TwoBitFile

            self._twobit = TwoBitFile(twobit_path)

    @classmethod
    def from_dbloc(cls, dbloc: str, buildver: str = "hg19") -> "GenomeSequence":
        fasta = resolve_genome_fasta(dbloc, buildver)
        twobit = resolve_genome_2bit(dbloc, buildver)
        if not fasta and not twobit:
            raise FileNotFoundError(
                f"no genome FASTA or 2bit under {dbloc} (expected Homo_sapiens_assembly37.fasta or {buildver}.2bit)"
            )
        return cls(fasta_path=fasta, twobit_path=twobit)

    @property
    def ready(self) -> bool:
        return self._fasta is not None or self._twobit is not None

    def close(self) -> None:
        if self._twobit is not None:
            self._twobit.close()
            self._twobit = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()

    def _fasta_chrom(self, chrom: str) -> Optional[str]:
        if self._fasta is None:
            return None
        for alias in chrom_aliases(chrom):
            if alias in self._fasta:
                return alias
        return None

    def fetch(self, chrom: str, start0: int, end0: int) -> str:
        """0-based half-open genomic sequence, uppercase."""
        if self._fasta is not None:
            key = self._fasta_chrom(chrom)
            if key is None:
                raise KeyError(f"chromosome not in FASTA: {chrom}")
            return str(self._fasta[key][start0:end0]).upper()
        if self._twobit is not None:
            return self._twobit.fetch(chrom, start0, end0).upper()
        raise RuntimeError("no genome sequence backend")

    def fetch_1based(self, chrom: str, start: int, end: int) -> str:
        """1-based closed interval."""
        return self.fetch(chrom, start - 1, end)


def _open_pybigwig(path: str):
    try:
        import pyBigWig
    except ImportError as exc:
        raise RuntimeError(
            "reading bigBed (*.bb) requires pyBigWig; install with: pip install pyBigWig"
        ) from exc
    bw = pyBigWig.open(path)
    if bw is None:
        raise RuntimeError(f"cannot open bigBed: {path}")
    return bw


def _parse_int_list(text: str) -> List[int]:
    return [int(x) for x in (text or "").rstrip(",").split(",") if x]


def parse_bigbed_gene_entry(
    chrom: str, start: int, end: int, rest: str
) -> Optional[Dict[str, object]]:
    """Parse BED12 / bigGenePred extra fields into a genePred-like dict (0-based)."""
    parts = rest.strip().split("\t") if rest else []
    if len(parts) < 8:
        return None
    name = parts[0].strip()
    if not name:
        return None
    strand = parts[2].strip() if len(parts) > 2 else "+"
    if strand not in {"+", "-"}:
        strand = "+"
    try:
        thick_start = int(parts[3])
        thick_end = int(parts[4])
        block_count = int(parts[6])
        block_sizes = _parse_int_list(parts[7])
        block_starts = _parse_int_list(parts[8]) if len(parts) > 8 else []
    except (ValueError, IndexError):
        return None
    if len(block_sizes) != block_count or len(block_starts) != block_count:
        return None
    exon_starts = [start + rel for rel in block_starts]
    exon_ends = [s + size for s, size in zip(exon_starts, block_sizes)]
    gene = ""
    if len(parts) > 15 and parts[15].strip():
        gene = parts[15].strip()
    elif len(parts) > 9:
        gene = parts[9].strip()
    cds_stat_s = parts[10].strip() if len(parts) > 10 else "cmpl"
    cds_stat_e = parts[11].strip() if len(parts) > 11 else "cmpl"
    return {
        "name": name,
        "chrom": chrom,
        "strand": strand,
        "txStart": start,
        "txEnd": end,
        "cdsStart": thick_start,
        "cdsEnd": thick_end,
        "exonCount": block_count,
        "exonStarts": exon_starts,
        "exonEnds": exon_ends,
        "name2": gene or name,
        "cdsStartStat": cds_stat_s if cds_stat_s in {"none", "unk", "incmpl", "cmpl"} else "cmpl",
        "cdsEndStat": cds_stat_e if cds_stat_e in {"none", "unk", "incmpl", "cmpl"} else "cmpl",
    }


def _gene_pred_line(rec: Dict[str, object]) -> str:
    starts = ",".join(str(x) for x in rec["exonStarts"]) + ","
    ends = ",".join(str(x) for x in rec["exonEnds"]) + ","
    return "\t".join(
        [
            "0",
            str(rec["name"]),
            str(rec["chrom"]),
            str(rec["strand"]),
            str(rec["txStart"]),
            str(rec["txEnd"]),
            str(rec["cdsStart"]),
            str(rec["cdsEnd"]),
            str(rec["exonCount"]),
            starts,
            ends,
            "0",
            str(rec.get("name2") or rec["name"]),
            str(rec.get("cdsStartStat") or "cmpl"),
            str(rec.get("cdsEndStat") or "cmpl"),
        ]
    )


def materialize_bigbed_genepred(bb_path: str, out_path: Optional[str] = None) -> str:
    """Convert a gene bigBed to UCSC refGene-style genePred text (cached beside the bb)."""
    bb_path = os.path.abspath(bb_path)
    if out_path is None:
        out_path = bb_path + ".genePred.txt"
    if os.path.isfile(out_path) and os.path.getmtime(out_path) >= os.path.getmtime(bb_path):
        return out_path

    bw = _open_pybigwig(bb_path)
    records: List[Dict[str, object]] = []
    try:
        chroms = bw.chroms() or {}
        ordered = sorted(chroms.items(), key=lambda kv: (0 if is_primary_chrom(kv[0]) else 1, kv[0]))
        for chrom, length in ordered:
            entries = bw.entries(chrom, 0, int(length), withString=True) or []
            for start, end, rest in entries:
                rec = parse_bigbed_gene_entry(chrom, int(start), int(end), rest or "")
                if rec:
                    records.append(rec)
    finally:
        bw.close()
    if not records:
        raise RuntimeError(f"no transcripts parsed from {bb_path}")

    tmp = out_path + ".tmp"
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with open(tmp, "w", encoding="utf-8") as fh:
        for rec in records:
            fh.write(_gene_pred_line(rec) + "\n")
    os.replace(tmp, out_path)
    logger.info("materialized %s transcripts from %s → %s", len(records), bb_path, out_path)
    return out_path


def ensure_gene_pred(dbloc: str, protocol: str, buildver: str = "hg19") -> Optional[str]:
    """Return a genePred .txt path for annotation (classic file or materialized bigBed)."""
    proto = protocol
    if proto.lower() in ("gene", "refgene"):
        proto = "refGene"
    elif proto.lower() == "ensgene":
        proto = "ensGene"
    elif proto.lower() == "knowngene":
        proto = "knownGene"
    elif proto.lower() == "ncbirefseq":
        proto = "ncbiRefSeq"
    elif proto.lower() == "gencode":
        proto = "gencode"

    txt = resolve_gene_txt(dbloc, proto, buildver)
    if txt and not txt.endswith(".bb"):
        # Prefer real genePred text over a later bb fallback, except when the
        # hit is a different protocol's file used only as last resort.
        if proto in ("refGene", "ensGene", "knownGene"):
            stem = os.path.basename(txt)
            if proto.lower() in stem.lower() or f"{buildver}_{proto}" in stem:
                return txt
            if proto == "refGene" and "refGene" in stem:
                return txt
            if proto == "ensGene" and "ensGene" in stem:
                return txt
        else:
            return txt
        # refGene.txt missing but ncbiRefSeq.txt present — use it.
        if proto == "refGene":
            return txt

    bb = resolve_gene_bigbed(dbloc, proto, buildver)
    if bb:
        try:
            return materialize_bigbed_genepred(bb)
        except Exception as exc:
            logger.warning("failed to materialize %s: %s", bb, exc)
    return txt


def load_chrom_sizes(path: str) -> Dict[str, int]:
    sizes: Dict[str, int] = {}
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    sizes[parts[0]] = int(parts[1])
                except ValueError:
                    continue
    return sizes


def parse_evf_transcript_ids(evf_path: str) -> Set[str]:
    ids: Set[str] = set()
    if not evf_path or not os.path.isfile(evf_path):
        return ids
    pat = re.compile(r"\b((?:NM|NR|XM|XR|ENST)[_:]?[0-9]+(?:\.[0-9]+)?)\b", re.I)
    with open(evf_path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            for m in pat.finditer(line):
                ids.add(m.group(1))
                ids.add(m.group(1).split(".")[0])
    return ids


def _load_genepred_records(gene_file: str) -> Dict[str, Dict[str, object]]:
    recs: Dict[str, Dict[str, object]] = {}
    opener = gzip.open if gene_file.endswith(".gz") else open
    with opener(gene_file, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            # bin, name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds, score, name2, ...
            if len(parts) < 11:
                continue
            has_bin = parts[0].isdigit() and len(parts) > 3 and parts[3] in ("+", "-")
            if has_bin and len(parts) >= 11:
                name, chrom, strand = parts[1], parts[2], parts[3]
                tx_s, tx_e = int(parts[4]), int(parts[5])
                cds_s, cds_e = int(parts[6]), int(parts[7])
                starts = _parse_int_list(parts[9])
                ends = _parse_int_list(parts[10])
                gene = parts[12] if len(parts) > 12 else name
            else:
                name, chrom, strand = parts[0], parts[1], parts[2]
                tx_s, tx_e = int(parts[3]), int(parts[4])
                cds_s, cds_e = int(parts[5]), int(parts[6])
                starts = _parse_int_list(parts[8])
                ends = _parse_int_list(parts[9])
                gene = parts[11] if len(parts) > 11 else name
            rec = {
                "name": name,
                "chrom": chrom,
                "strand": strand,
                "txStart": tx_s,
                "txEnd": tx_e,
                "cdsStart": cds_s,
                "cdsEnd": cds_e,
                "exonStarts": starts,
                "exonEnds": ends,
                "name2": gene,
            }
            recs[name] = rec
            recs.setdefault(name.split(".")[0], rec)
    return recs


def splice_transcript_seq(rec: Dict[str, object], genome: GenomeSequence) -> str:
    starts: Sequence[int] = rec["exonStarts"]  # type: ignore[assignment]
    ends: Sequence[int] = rec["exonEnds"]  # type: ignore[assignment]
    chrom = str(rec["chrom"])
    parts = [genome.fetch(chrom, int(s), int(e)) for s, e in zip(starts, ends)]
    seq = "".join(parts)
    if rec.get("strand") == "-":
        seq = reverse_complement(seq)
    return seq


def extract_mrna_fasta(
    gene_file: str,
    genome: GenomeSequence,
    outfile: str,
    transcript_ids: Optional[Iterable[str]] = None,
) -> str:
    """Write spliced transcript FASTA for polish (header = transcript id)."""
    recs = _load_genepred_records(gene_file)
    wanted: Optional[Set[str]] = None
    if transcript_ids is not None:
        wanted = set()
        for tid in transcript_ids:
            wanted.add(tid)
            wanted.add(tid.split(".")[0])
        if not wanted:
            raise RuntimeError("no transcript IDs to splice from EVF")
    written = 0
    os.makedirs(os.path.dirname(os.path.abspath(outfile)) or ".", exist_ok=True)
    with open(outfile, "w", encoding="utf-8") as fh:
        seen: Set[str] = set()
        for key, rec in recs.items():
            name = str(rec["name"])
            base = name.split(".")[0]
            if name in seen:
                continue
            if wanted is not None and name not in wanted and base not in wanted:
                continue
            seen.add(name)
            try:
                seq = splice_transcript_seq(rec, genome)
            except Exception as exc:
                logger.debug("skip transcript %s: %s", name, exc)
                continue
            if not seq:
                continue
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")
            written += 1
    logger.info("wrote %s spliced transcripts to %s", written, outfile)
    return outfile


def discover_gene_protocols(dbloc: str, buildver: str = "hg19") -> List[Tuple[str, str]]:
    """Return (protocol, operation) pairs available in humandb."""
    out: List[Tuple[str, str]] = []
    ref_txt = resolve_gene_txt(dbloc, "refGene", buildver)
    ref_bb = resolve_gene_bigbed(dbloc, "ncbiRefSeq", buildver)
    if ref_txt and "refgene" in os.path.basename(ref_txt).lower():
        out.append(("refGene", "g"))
        if ref_bb:
            out.append(("ncbiRefSeq", "g"))
    elif ref_bb:
        out.append(("ncbiRefSeq", "g"))
    elif ref_txt:
        out.append(("refGene", "g"))

    ens_txt = resolve_gene_txt(dbloc, "ensGene", buildver)
    gencode_bb = resolve_gene_bigbed(dbloc, "gencode", buildver)
    if ens_txt and "ensgene" in os.path.basename(ens_txt).lower():
        out.append(("ensGene", "g"))
        if gencode_bb:
            out.append(("gencode", "g"))
    elif gencode_bb:
        out.append(("gencode", "g"))
    elif ens_txt:
        out.append(("ensGene", "g"))

    if resolve_gene_txt(dbloc, "knownGene", buildver):
        out.append(("knownGene", "g"))
    if resolve_cytoband(dbloc, buildver):
        out.append(("cytoBand", "r"))

    seen = set()
    uniq: List[Tuple[str, str]] = []
    for p, op in out:
        if p not in seen:
            seen.add(p)
            uniq.append((p, op))
    return uniq


def resolve_region_file(dbloc: str, protocol: str, buildver: str = "hg19") -> Optional[str]:
    if protocol == "cytoBand":
        return resolve_cytoband(dbloc, buildver)
    return _first_existing(
        [
            os.path.join(dbloc, f"{buildver}_{protocol}.txt"),
            os.path.join(dbloc, f"{buildver}_{protocol}.txt.gz"),
            os.path.join(dbloc, f"{protocol}.txt"),
            os.path.join(dbloc, f"{protocol}.txt.gz"),
        ]
    )
