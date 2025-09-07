#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build bgzip compression and Tabix indexing for the large database under resources/humandb.

 - Auto-discover files matching {buildver}_*.txt
 - Produce {file}.gz and {file}.gz.tbi (skip if already exists)
 - Column mapping defaults: chrom=1st, start=2nd, end=2nd (for SNV). Use --end-col for interval databases
 - Header/comment lines starting with '#' are handled implicitly (no explicit meta argument)

Examples:
  python utils/build_tabix_indexes.py --humandb resources/humandb --buildver hg19 --threads 8
  python utils/build_tabix_indexes.py --humandb resources/humandb --buildver hg19 --end-col 3
  python utils/build_tabix_indexes.py --humandb resources/humandb --buildver hg19 --threads 8 --min-size-gb 10
  # Use improved indexing script with forced pre-sorting
  python utils/build_tabix_indexes.py --humandb resources/humandb --buildver hg19 --pre-sort --sort-buffer-mb 4096 --threads 8
"""

import os
import sys
import argparse
import glob
import traceback
from concurrent.futures import ThreadPoolExecutor, as_completed
import tempfile
import shutil
import subprocess

try:
    import pysam
except Exception as e:
    print("Missing dependency 'pysam'. Please install: pip install pysam", file=sys.stderr)
    raise


def is_already_indexed(txt_path: str, force_rebuild: bool = False) -> bool:
    if force_rebuild:
        return False
    gz = txt_path + ".gz"
    tbi = gz + ".tbi"
    return os.path.exists(gz) and os.path.exists(tbi)


def verify_tabix_index(gz_path: str, test_queries: list = None) -> bool:
    """
    Verify that a Tabix index is working correctly by testing queries.
    
    Args:
        gz_path: Path to the compressed file
        test_queries: List of (chrom, start, end) tuples to test
        
    Returns:
        True if index is working, False otherwise
    """
    try:
        tbx = pysam.TabixFile(gz_path)
        
        # Test basic functionality
        contigs = tbx.contigs
        if not contigs:
            print(f"Warning: No contigs found in {gz_path}")
            return False
        
        # Test specific queries if provided
        if test_queries:
            for chrom, start, end in test_queries:
                try:
                    results = list(tbx.fetch(chrom, start, end))
                    print(f"Query {chrom}:{start}-{end}: {len(results)} results")
                except Exception as e:
                    print(f"Query {chrom}:{start}-{end} failed: {e}")
                    return False
        
        tbx.close()
        return True
        
    except Exception as e:
        print(f"Tabix verification failed for {gz_path}: {e}")
        return False


def diagnose_index_issues(gz_path: str) -> dict:
    """
    Diagnose potential issues with a Tabix index.
    
    Args:
        gz_path: Path to the compressed file
        
    Returns:
        Dictionary with diagnosis results
    """
    diagnosis = {
        'file_exists': os.path.exists(gz_path),
        'tbi_exists': os.path.exists(gz_path + '.tbi'),
        'file_size': 0,
        'tbi_size': 0,
        'can_open': False,
        'contigs': [],
        'sample_queries': []
    }
    
    if not diagnosis['file_exists']:
        return diagnosis
    
    diagnosis['file_size'] = os.path.getsize(gz_path)
    if diagnosis['tbi_exists']:
        diagnosis['tbi_size'] = os.path.getsize(gz_path + '.tbi')
    
    try:
        tbx = pysam.TabixFile(gz_path)
        diagnosis['can_open'] = True
        diagnosis['contigs'] = list(tbx.contigs)
        
        # Test a few sample queries
        for contig in diagnosis['contigs'][:3]:  # Test first 3 contigs
            try:
                results = list(tbx.fetch(contig, 1, 1000))
                diagnosis['sample_queries'].append({
                    'contig': contig,
                    'query': f"{contig}:1-1000",
                    'results': len(results),
                    'success': True
                })
            except Exception as e:
                diagnosis['sample_queries'].append({
                    'contig': contig,
                    'query': f"{contig}:1-1000",
                    'error': str(e),
                    'success': False
                })
        
        tbx.close()
        
    except Exception as e:
        diagnosis['error'] = str(e)
    
    return diagnosis


def _detect_end_col(txt_path: str, default_end_col: int = 2) -> int:
    """Heuristically detect if file is an interval database: if any row has col2 != col3,
    treat it as interval (end column = 3), otherwise return 2.
    """
    try:
        with open(txt_path, 'r', encoding='utf-8', errors='replace') as f:
            checked = 0
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) < 3:
                    break
                checked += 1
                if parts[1] != parts[2]:
                    return 3
                if checked >= 100:
                    break
    except Exception:
        return default_end_col
    return 2


def _external_presort(txt_path: str, sorted_path: str, tempdir: str, buffer_mb: int) -> bool:
    """Use GNU coreutils sort to pre-sort by chr then start (numeric), preserving header lines (#)."""
    try:
        os.makedirs(tempdir, exist_ok=True)
    except Exception:
        pass
    
    # Check if GNU sort is available
    try:
        subprocess.run(['sort', '--version'], capture_output=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        print(f"GNU sort not available, cannot pre-sort: {txt_path}", file=sys.stderr)
        return False
    
    # Compose shell command: write headers first, then sorted body
    # Use parallel sort if available, with proper error handling
    cmd = (
        f"(grep '^#' '{txt_path}' > '{sorted_path}' && "
        f"grep -v '^#' '{txt_path}' | sort -S {buffer_mb}M -T '{tempdir}' -k1,1 -k2,2n --parallel=4 >> '{sorted_path}')"
    )
    
    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=3600)  # 1 hour timeout
        if result.returncode != 0:
            print(f"Sort command failed: {result.stderr}", file=sys.stderr)
            return False
        return os.path.exists(sorted_path) and os.path.getsize(sorted_path) > 0
    except subprocess.TimeoutExpired:
        print(f"Sort command timed out for: {txt_path}", file=sys.stderr)
        return False
    except Exception as e:
        print(f"Sort command error: {e}", file=sys.stderr)
        return False


def build_index_for_file(txt_path: str, end_col: int, skip_header: int,
                         pre_sort: bool, sort_buffer_mb: int, sort_tempdir: str, keep_sorted: bool) -> bool:
    gz = txt_path + ".gz"
    tbi = gz + ".tbi"
    try:
        # Select column mapping based on filename
        seq_col, start_col, end_col0 = _select_tabix_columns(txt_path, end_col)

        # Optional pre-sort
        src_path = txt_path
        sorted_path = txt_path + ".sorted"
        if pre_sort:
            tmpdir = sort_tempdir or tempfile.gettempdir()
            ok = _external_presort(txt_path, sorted_path, tmpdir, sort_buffer_mb)
            if not ok:
                print(f"Failed to pre-sort (external sort): {txt_path}", file=sys.stderr)
            else:
                src_path = sorted_path

        # Compress (use original output name for .gz, even if sorted source is used)
        pysam.tabix_compress(src_path, gz, force=True)
        # Build index: rely on default handling of comment lines ('#'); do not pass meta for compatibility
        pysam.tabix_index(
            gz,
            seq_col=seq_col,
            start_col=start_col,
            end_col=end_col0,
            preset=None,
            force=True
        )
        
        # Verify the index was created correctly
        if not verify_tabix_index(gz):
            print(f"Warning: Index verification failed for {txt_path}", file=sys.stderr)
            return False
        
        # cleanup sorted
        if src_path == sorted_path and not keep_sorted:
            try:
                os.remove(sorted_path)
            except Exception:
                pass
        return True
    except Exception as e:
        print(f"Failed to build index: {txt_path} -> {e}", file=sys.stderr)
        traceback.print_exc()
        # Cleanup partial outputs
        try:
            if os.path.exists(gz):
                os.remove(gz)
            if os.path.exists(tbi):
                os.remove(tbi)
        except Exception:
            pass
        # Retry once with pre-sort if unsorted positions caused failure
        msg = str(e)
        print(f"Error message: {msg}", file=sys.stderr)
        need_retry = ('Unsorted' in msg or 'unsorted' in msg or 'hts_idx_push' in msg or 'building of index' in msg)
        print(f"Need retry: {need_retry}, pre_sort: {pre_sort}", file=sys.stderr)
        if not pre_sort and need_retry:
            try:
                print(f"Retrying with external pre-sort due to unsorted positions: {txt_path}", file=sys.stderr)
                tmpdir = sort_tempdir or tempfile.gettempdir()
                sorted_path = txt_path + ".sorted"
                if _external_presort(txt_path, sorted_path, tmpdir, sort_buffer_mb):
                    # Clean up any partial files from first attempt
                    try:
                        if os.path.exists(gz):
                            os.remove(gz)
                        if os.path.exists(tbi):
                            os.remove(tbi)
                    except Exception:
                        pass
                    # Try again compress+index using sorted source
                    pysam.tabix_compress(sorted_path, gz, force=True)
                    pysam.tabix_index(
                        gz,
                        seq_col=seq_col,
                        start_col=start_col,
                        end_col=end_col0,
                        preset=None,
                        force=True
                    )
                    if not keep_sorted:
                        try:
                            os.remove(sorted_path)
                        except Exception:
                            pass
                    print(f"Successfully built index after pre-sort: {txt_path}", file=sys.stderr)
                    return True
                else:
                    print(f"Pre-sort failed, cannot retry: {txt_path}", file=sys.stderr)
            except Exception as e2:
                print(f"Retry after pre-sort failed: {txt_path} -> {e2}", file=sys.stderr)
                traceback.print_exc()
        return False


def discover_txt_files(humandb: str, buildver: str) -> list:
    pattern = os.path.join(humandb, f"{buildver}_*.txt") if buildver else os.path.join(humandb, "*.txt")
    return sorted(glob.glob(pattern))


def discover_gz_files(humandb: str, buildver: str) -> list:
    """Discover compressed files for diagnosis"""
    pattern = os.path.join(humandb, f"{buildver}_*.txt.gz") if buildver else os.path.join(humandb, "*.txt.gz")
    return sorted(glob.glob(pattern))


def _select_tabix_columns(txt_path: str, cli_end_col: int):
    """
    Choose column mapping by filename, return 0-based (seq_col, start_col, end_col).
    - *_ensGene.txt / *_knownGene.txt / *_refGene.txt: UCSC gene table format
      1:bin 2:name 3:chrom 4:strand 5:txStart 6:txEnd ... => seq=3rd(2), start=5th(4), end=6th(5)
    - *cytoBand.txt: 1:chr 2:start 3:end => (0,1,2)
    - *clinvar.txt: 1:Chr 2:Start 3:End 4:Ref 5:Alt ... => (0,1,2)
    - Others default: 1:chr 2:start 3:end
      If CLI provides --end-col, use it; otherwise try _detect_end_col()
    """
    base = os.path.basename(txt_path)
    name = base.lower()
    if any(k in name for k in ["ensgene", "knowngene", "refgene"]):
        return 2, 4, 5
    if "cytoband" in name:
        return 0, 1, 2
    if "clinvar" in name:
        # ClinVar format: #Chr Start End Ref Alt ...
        return 0, 1, 2
    # default
    end_col_1based = cli_end_col if cli_end_col else _detect_end_col(txt_path, 2)
    end_col_0based = max(0, end_col_1based - 1)
    return 0, 1, end_col_0based


def main():
    parser = argparse.ArgumentParser(description="Build bgzip+tabix indexes for humandb txt files")
    parser.add_argument("--humandb", default=os.path.join("resources", "humandb"), help="Path to humandb directory")
    parser.add_argument("--buildver", default="hg19", help="Genome build prefix, e.g., hg19/hg38. Empty to process all")
    parser.add_argument("--end-col", type=int, default=2, help="End column (1-based). Use 2 for SNV; interval DBs may use 3")
    parser.add_argument("--skip-header", type=int, default=1, help="Skip header lines starting with '#': 1=yes, 0=no")
    parser.add_argument("--threads", type=int, default=4, help="Parallel compression threads")
    parser.add_argument("--min-size-gb", type=float, default=5.0, help="Only index TXT files with size ≥ this threshold (GB). Default 5GB")
    parser.add_argument("--pre-sort", action="store_true", help="Force external pre-sort before indexing")
    parser.add_argument("--sort-buffer-mb", type=int, default=2048, help="Buffer size for GNU sort (MB)")
    parser.add_argument("--sort-tempdir", type=str, default="", help="Temp directory for GNU sort (default: system temp)")
    parser.add_argument("--keep-sorted", action="store_true", help="Keep the intermediate .sorted file")
    parser.add_argument("--verify", action="store_true", help="Verify existing indexes after building")
    parser.add_argument("--force-rebuild", action="store_true", help="Force rebuild even if index already exists")
    parser.add_argument("--diagnose", action="store_true", help="Diagnose existing indexes without rebuilding")
    args = parser.parse_args()

    humandb = os.path.abspath(args.humandb)
    if not os.path.isdir(humandb):
        print(f"humandb directory not found: {humandb}", file=sys.stderr)
        sys.exit(1)

    txt_files = discover_txt_files(humandb, args.buildver)
    if not txt_files:
        print(f"No matching TXT files found: {humandb}/{args.buildver}_*.txt")
        return

    # Handle diagnosis mode
    if args.diagnose:
        print("Diagnosing existing Tabix indexes...\n")
        gz_files = discover_gz_files(humandb, args.buildver)
        if not gz_files:
            print(f"No compressed files found: {humandb}/{args.buildver}_*.txt.gz")
            return
        
        for gz_file in gz_files:
            print(f"Diagnosing {os.path.basename(gz_file)}:")
            diagnosis = diagnose_index_issues(gz_file)
            print(f"  File exists: {diagnosis['file_exists']}")
            print(f"  TBI exists: {diagnosis['tbi_exists']}")
            print(f"  File size: {diagnosis['file_size']:,} bytes")
            print(f"  TBI size: {diagnosis['tbi_size']:,} bytes")
            print(f"  Can open: {diagnosis['can_open']}")
            print(f"  Contigs: {diagnosis['contigs']}")
            if 'error' in diagnosis:
                print(f"  Error: {diagnosis['error']}")
            for query in diagnosis['sample_queries']:
                if query['success']:
                    print(f"  Query {query['query']}: {query['results']} results")
                else:
                    print(f"  Query {query['query']}: FAILED - {query['error']}")
            print()
        return

    print(f"Found {len(txt_files)} TXT files, filtering (≥{args.min_size_gb}GB and non-gene tables)...\n")

    # Skip gene-definition tables and small files
    min_bytes = int(args.min_size_gb * (1024**3))
    def need_index(path: str) -> bool:
        base = os.path.basename(path).lower()
        if any(k in base for k in ["ensgene", "knowngene", "refgene"]):
            return False
        try:
            sz = os.path.getsize(path)
        except Exception:
            return False
        return sz >= min_bytes

    candidates = [p for p in txt_files if need_index(p)]
    skipped_small_or_gene = len(txt_files) - len(candidates)
    if skipped_small_or_gene:
        print(f"Skipped by policy (below threshold or gene table): {skipped_small_or_gene}. To process: {len(candidates)}.\n")

    todo = [f for f in candidates if not is_already_indexed(f, args.force_rebuild)]
    skipped = len(candidates) - len(todo)
    if skipped and not args.force_rebuild:
        print(f"Already indexed, skipped {skipped} files. To index: {len(todo)}.\n")
    elif args.force_rebuild:
        print(f"Force rebuild enabled, will rebuild {len(todo)} files.\n")

    success = 0
    if todo:
        with ThreadPoolExecutor(max_workers=max(1, args.threads)) as ex:
            fut_map = {ex.submit(build_index_for_file, f, args.end_col, args.skip_header,
                                 args.pre_sort, args.sort_buffer_mb, args.sort_tempdir, args.keep_sorted): f for f in todo}
            for fut in as_completed(fut_map):
                f = fut_map[fut]
                ok = False
                try:
                    ok = fut.result()
                except Exception as e:
                    print(f"Index build exception: {f} -> {e}", file=sys.stderr)
                if ok:
                    success += 1
                    print(f"Done: {os.path.basename(f)}.gz(.tbi)")
                else:
                    print(f"Failed: {os.path.basename(f)}", file=sys.stderr)

    print("\nIndex build finished:")
    print(f"  Total: {len(txt_files)}")
    print(f"  Skipped (already indexed): {skipped}")
    print(f"  Succeeded: {success}")
    print(f"  Failed: {len(todo) - success}")
    
    # Verify indexes if requested
    if args.verify and success > 0:
        print("\nVerifying indexes...")
        verified_count = 0
        for txt_file in txt_files:
            gz_file = txt_file + ".gz"
            if os.path.exists(gz_file):
                if verify_tabix_index(gz_file):
                    verified_count += 1
                    print(f"✅ {os.path.basename(gz_file)}")
                else:
                    print(f"❌ {os.path.basename(gz_file)}")
        print(f"Verified {verified_count}/{len(txt_files)} indexes")


if __name__ == "__main__":
    main()
