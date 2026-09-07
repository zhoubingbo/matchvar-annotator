#!/usr/bin/env python3
"""
MATCHVAR编码变化分析工具
分析DNA水平变异对蛋白质序列的影响
"""

import os
import sys
import argparse
import logging
import re
from typing import Dict, List, Tuple, Optional

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# 密码子表
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TCT': 'S', 'TCC': 'S', 'TAT': 'Y', 'TAC': 'Y',
    'TGT': 'C', 'TGC': 'C', 'TTA': 'L', 'TCA': 'S', 'TAA': '*', 'TGA': '*',
    'TTG': 'L', 'TCG': 'S', 'TAG': '*', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L',
    'CCT': 'P', 'CCC': 'P', 'CAT': 'H', 'CAC': 'H', 'CGT': 'R', 'CGC': 'R',
    'CTA': 'L', 'CTG': 'L', 'CCA': 'P', 'CCG': 'P', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ACT': 'T', 'ACC': 'T',
    'AAT': 'N', 'AAC': 'N', 'AGT': 'S', 'AGC': 'S', 'ATA': 'I', 'ACA': 'T',
    'AAA': 'K', 'AGA': 'R', 'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GCT': 'A', 'GCC': 'A', 'GAT': 'D', 'GAC': 'D',
    'GGT': 'G', 'GGC': 'G', 'GTA': 'V', 'GTG': 'V', 'GCA': 'A', 'GCG': 'A',
    'GAA': 'E', 'GAG': 'E', 'GGA': 'G', 'GGG': 'G'
}

# 线粒体密码子表
CODON_TABLE_MT = {
    'TTT': 'F', 'TTC': 'F', 'TCT': 'S', 'TCC': 'S', 'TAT': 'Y', 'TAC': 'Y',
    'TGT': 'C', 'TGC': 'C', 'TTA': 'L', 'TCA': 'S', 'TAA': '*', 'TGA': 'W',
    'TTG': 'L', 'TCG': 'S', 'TAG': '*', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L',
    'CCT': 'P', 'CCC': 'P', 'CAT': 'H', 'CAC': 'H', 'CGT': 'R', 'CGC': 'R',
    'CTA': 'L', 'CTG': 'L', 'CCA': 'P', 'CCG': 'P', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ACT': 'T', 'ACC': 'T',
    'AAT': 'N', 'AAC': 'N', 'AGT': 'S', 'AGC': 'S', 'ATA': 'M', 'ACA': 'T',
    'AAA': 'K', 'AGA': '*', 'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': '*',
    'GTT': 'V', 'GTC': 'V', 'GCT': 'A', 'GCC': 'A', 'GAT': 'D', 'GAC': 'D',
    'GGT': 'G', 'GGC': 'G', 'GTA': 'V', 'GTG': 'V', 'GCA': 'A', 'GCG': 'A',
    'GAA': 'E', 'GAG': 'E', 'GGA': 'G', 'GGG': 'G'
}

# 氨基酸单字母到三字母映射（HGVS p. 表示用三字母）
AA_ONE_TO_THREE = {
    'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
    'Q': 'Gln', 'E': 'Glu', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
    'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
    'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val',
    '*': 'Ter', 'X': 'Xaa'
}

def reverse_complement(seq: str) -> str:
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(comp.get(b, 'N') for b in reversed(seq.upper()))

# ===================== 通用 c.HGVS -> p. 引擎辅助函数 =====================

def translate_protein(dna: str, chrom: str) -> str:
    dna = dna.upper()
    protein = []
    table = CODON_TABLE_MT if chrom.upper() in ['M', 'MT', 'CHRM', 'CHRMT'] else CODON_TABLE
    for i in range(0, len(dna) - 2, 3):
        codon = dna[i:i+3]
        aa = table.get(codon, 'X')
        protein.append(aa)
        if aa == '*':
            break
    return ''.join(protein)

def three_letter(aa: str) -> str:
    return AA_ONE_TO_THREE.get(aa, 'Xaa')

def diff_proteins(wt: str, mut: str) -> Tuple[int, str, str, int]:
    """返回 (pos1, aa1_diff, aa2_diff, pos2)
    pos1/pos2 为1-based位置窗口，aa1_diff/aa2_diff为差异段（可能为空）"""
    i = 0
    L1, L2 = len(wt), len(mut)
    while i < L1 and i < L2 and wt[i] == mut[i]:
        i += 1
    if i == L1 and i == L2:
        return 0, '', '', 0  # 完全相同
    # 从尾部修剪公共后缀（仅在非frameshift情况下处理更合理，这里统一先取差异段）
    j1, j2 = L1 - 1, L2 - 1
    while j1 >= i and j2 >= i and wt[j1] == mut[j2]:
        j1 -= 1
        j2 -= 1
    aa1_diff = wt[i:j1+1]
    aa2_diff = mut[i:j2+1]
    return i + 1, aa1_diff, aa2_diff, (j1 + 1)  # 1-based 起始，pos2 使用wt端结束位置

def format_p_hgvs_from_diff(wt: str, mut: str, chrom: str, effect_hint: Optional[str] = None) -> Tuple[str, str]:
    """根据蛋白差异生成 p. 注释与 effect。

    整码 indel/dup/delins/错义/终止丢失逻辑对齐 test/local_hgvs_frameshift.format_inframe_protein_hgvs；
    移码仍输出 MATCHVAR 惯用的 fs*N。
    translate_protein 会把终止符 '*' 留在串内；整码比较时先去掉，与 TransVar 肽段一致。
    """
    if wt == mut:
        return 'p.(=)', 'synonymous_SNV'

    # frameshift：保留 '*'，用于计算 fs*N
    if effect_hint == 'frameshift':
        start_idx = 0
        while start_idx < len(wt) and start_idx < len(mut) and wt[start_idx] == mut[start_idx]:
            start_idx += 1
        stop_idx = mut.find('*', start_idx)
        if stop_idx == -1:
            fs_len = max(0, len(mut) - start_idx)
        else:
            fs_len = stop_idx - start_idx + 1  # 含终止密码子本身
        ref_aa = wt[start_idx] if start_idx < len(wt) and wt[start_idx] != '*' else '?'
        new_aa = mut[start_idx] if start_idx < len(mut) and mut[start_idx] != '*' else '?'
        pos1 = start_idx + 1
        p = f'p.{three_letter(ref_aa)}{pos1}{three_letter(new_aa)}fs*{fs_len}'
        return p, 'frameshift'

    # 整码路径：去掉终止符再比较（对齐 TransVar _translate_peptide）
    wt = wt.rstrip('*')
    mut = mut.rstrip('*')
    if wt == mut:
        return 'p.(=)', 'synonymous_SNV'

    start = 0
    while start < len(wt) and start < len(mut) and wt[start] == mut[start]:
        start += 1

    # 终止丢失 / C 端延伸：原肽已走完，突变肽仍继续
    if start >= len(wt) and start < len(mut):
        stop_pos = len(wt) + 1  # 原终止密码子位置（1-based）
        first_aa = mut[start]
        ext_ter = (len(mut) - start) + 1
        return (
            f'p.Ter{stop_pos}{three_letter(first_aa)}extTer{ext_ter}',
            'stoploss',
        )

    # 提前终止：突变肽更短且无公共后缀可剪到非空 mid（新终止截断）
    if start < len(wt) and start >= len(mut):
        return f'p.{three_letter(wt[start])}{start + 1}Ter', 'stopgain'

    end_wt = len(wt)
    end_mut = len(mut)
    while end_wt > start and end_mut > start and wt[end_wt - 1] == mut[end_mut - 1]:
        end_wt -= 1
        end_mut -= 1

    wt_mid = wt[start:end_wt]
    mut_mid = mut[start:end_mut]
    pos = start + 1  # 1-based

    # 纯插入（优先识别串联重复 → p.…dup）
    if not wt_mid and mut_mid:
        n_ins = len(mut_mid)
        if start >= n_ins and wt[start - n_ins:start] == mut_mid:
            left_pos = pos - n_ins
            if n_ins == 1:
                return f'p.{three_letter(mut_mid[0])}{left_pos}dup', 'nonframeshift_duplication'
            return (
                f'p.{three_letter(mut_mid[0])}{left_pos}_'
                f'{three_letter(mut_mid[-1])}{pos - 1}dup',
                'nonframeshift_duplication',
            )
        ins = ''.join(three_letter(x) for x in mut_mid)
        if start == 0:
            if not wt:
                return 'p.?', 'unknown'
            return (
                f'p.{three_letter(wt[0])}1delins{ins}{three_letter(wt[0])}',
                'nonframeshift_insertion',
            )
        left_aa = wt[start - 1]
        right_aa = wt[start]
        left_pos = pos - 1
        return (
            f'p.{three_letter(left_aa)}{left_pos}_{three_letter(right_aa)}{pos}ins{ins}',
            'nonframeshift_insertion',
        )

    # 纯缺失
    if wt_mid and not mut_mid:
        n_del = len(wt_mid)
        if n_del == 1:
            return f'p.{three_letter(wt_mid[0])}{pos}del', 'nonframeshift_deletion'
        end_pos = pos + n_del - 1
        return (
            f'p.{three_letter(wt_mid[0])}{pos}_'
            f'{three_letter(wt_mid[-1])}{end_pos}del',
            'nonframeshift_deletion',
        )

    # 替换（含单氨基酸错义 / delins）
    if wt_mid and mut_mid:
        if len(wt_mid) == 1 and len(mut_mid) == 1:
            if wt_mid == mut_mid:
                return f'p.{three_letter(wt_mid[0])}{pos}=', 'synonymous_SNV'
            return (
                f'p.{three_letter(wt_mid[0])}{pos}{three_letter(mut_mid[0])}',
                'nonsynonymous_SNV',
            )
        alt = ''.join(three_letter(x) for x in mut_mid)
        if len(wt_mid) == 1:
            return (
                f'p.{three_letter(wt_mid[0])}{pos}delins{alt}',
                'inframe_substitution',
            )
        end_pos = pos + len(wt_mid) - 1
        return (
            f'p.{three_letter(wt_mid[0])}{pos}_'
            f'{three_letter(wt_mid[-1])}{end_pos}delins{alt}',
            'inframe_substitution',
        )

    return 'p.?', 'unknown'

def apply_c_hgvs_to_cds(cds: str, cchange: str) -> Tuple[str, Optional[str]]:
    """将 c.HGVS 应用于 CDS，返回(突变后CDS, effect_hint)。effect_hint 可为 frameshift/nonframeshift/None"""
    s = cds
    c = cchange.strip()
    # 统一大写
    c = c.replace('c.', '')
    # 1) SNV: 180A>G 或 C180G
    m = re.match(r'^(\d+)([ACGTN])>([ACGTN])$', c, re.IGNORECASE)
    if m:
        pos = int(m.group(1))
        alt = m.group(3).upper()
        idx = pos - 1
        if 0 <= idx < len(s):
            return s[:idx] + alt + s[idx+1:], None
        return s, None
    m = re.match(r'^([ACGTN])(\d+)([ACGTN])$', c, re.IGNORECASE)
    if m:
        pos = int(m.group(2))
        alt = m.group(3).upper()
        idx = pos - 1
        if 0 <= idx < len(s):
            return s[:idx] + alt + s[idx+1:], None
        return s, None
    # 2) del / delSEQ / N_Mdel / N_MdelSEQ
    m = re.match(r'^(\d+)del([ACGTN]*)$', c, re.IGNORECASE)
    if m:
        start = int(m.group(1))
        end = start
        l = len(m.group(2) or '1')  # 未给出长度则按1
        end = start + l - 1
        return s[:start-1] + s[end:], 'frameshift' if l % 3 != 0 else 'nonframeshift'
    m = re.match(r'^(\d+)_(\d+)del([ACGTN]*)$', c, re.IGNORECASE)
    if m:
        start = int(m.group(1)); end = int(m.group(2))
        l = end - start + 1
        return s[:start-1] + s[end:], 'frameshift' if l % 3 != 0 else 'nonframeshift'
    # 3) ins: N_MinsSEQ
    m = re.match(r'^(\d+)_(\d+)ins([ACGTN]+)$', c, re.IGNORECASE)
    if m:
        left = int(m.group(1)); right = int(m.group(2))
        seq = m.group(3).upper()
        # HGVS 插入在两个碱基“之间”，此处按 left 位置后插入
        return s[:left] + seq + s[left:], 'frameshift' if len(seq) % 3 != 0 else 'nonframeshift'
    # 4) dup: Ndup[SEQ] / N_Mdup[SEQ]（重复；无 SEQ 时从 CDS 截取）
    m = re.match(r'^(\d+)dup([ACGTN]*)$', c, re.IGNORECASE)
    if m:
        pos = int(m.group(1))
        seq = (m.group(2) or s[pos - 1:pos]).upper()
        return s[:pos] + seq + s[pos:], 'dup' if len(seq) % 3 == 0 else 'frameshift'
    m = re.match(r'^(\d+)_(\d+)dup([ACGTN]*)$', c, re.IGNORECASE)
    if m:
        start = int(m.group(1)); end = int(m.group(2))
        seq = (m.group(3) or s[start - 1:end]).upper()
        # 在区间末端后插入重复序列（转录 3′ 侧）
        return s[:end] + seq + s[end:], 'dup' if len(seq) % 3 == 0 else 'frameshift'
    # 5) delins: NdelinsSEQ / N_MdelinsSEQ
    m = re.match(r'^(\d+)delins([ACGTN]+)$', c, re.IGNORECASE)
    if m:
        pos = int(m.group(1)); seq = m.group(2).upper()
        # 删除1个并插入seq
        return s[:pos-1] + seq + s[pos:], 'frameshift' if (1 - len(seq)) % 3 != 0 else 'nonframeshift'
    m = re.match(r'^(\d+)_(\d+)delins([ACGTN]+)$', c, re.IGNORECASE)
    if m:
        start = int(m.group(1)); end = int(m.group(2)); seq = m.group(3).upper()
        l = end - start + 1
        return s[:start-1] + seq + s[end:], 'frameshift' if ((l - len(seq)) % 3) != 0 else 'nonframeshift'
    # 6) 块替换：N_MSEQ（无delins关键词）
    m = re.match(r'^(\d+)_(\d+)([ACGTN]+)$', c, re.IGNORECASE)
    if m:
        start = int(m.group(1)); end = int(m.group(2)); seq = m.group(3).upper()
        l = end - start + 1
        return s[:start-1] + seq + s[end:], 'frameshift' if ((l - len(seq)) % 3) != 0 else 'nonframeshift'
    # 未能解析
    return s, None

def normalize_c_hgvs(cdot: str) -> str:
    """规范化一些非标准/冗余的 c.HGVS 表达，尽量贴近 HGVS 书写。
    - N_NdelinsSEQ -> NdelinsSEQ
    - del-ins -> delins
    - NBase>- -> Ndel（删除单碱基时可省略具体碱基）
    """
    cd = cdot
    cd = cd.replace('del-ins', 'delins').replace('del- ins', 'delins')
    m = re.match(r'^(c\.)?(\d+)_\2delins([ACGTN]+)$', cd, re.IGNORECASE)
    if m:
        cd = f"c.{m.group(2)}delins{m.group(3)}"
    # 单碱基删除 c.7202T>- 或 c.7202A>- → c.7202del
    m = re.match(r'^(c\.)?(\d+)[ACGTN]?>\-$', cd, re.IGNORECASE)
    if m:
        cd = f"c.{m.group(2)}del"
    # 若缺少前缀c.，补上
    if not cd.startswith('c.'):
        cd = 'c.' + cd
    return cd


def _format_dup_cdot(start: int, end: int) -> str:
    if start == end:
        return f"c.{start}dup"
    return f"c.{start}_{end}dup"


def _ins_to_dup_if_tandem(dna: str, left: int, ins: str) -> Optional[str]:
    """若插入序列等于其 5′ 侧邻接序列 → HGVS 优先记为 dup（3′ 端已对齐后调用）。"""
    k = len(ins)
    if k < 1 or left < k:
        return None
    if dna[left - k:left].upper() != ins.upper():
        return None
    return _format_dup_cdot(left - k + 1, left)


def normalize_coding_indel_hgvs(cdot: str, coding_dna: str) -> str:
    """按 HGVS 对 CDS 内 indel 做 3′ 规则归一化，并优先将串联插入写为 dup。

    - 仅处理纯编码区数字坐标（不含 c.*N / c.N+d / c.-N）
    - 3′ = CDS 坐标增大方向（转录本 3′）
    - 插入在 3′ 对齐后若与紧邻 5′ 序列相同 → ``c.a_bdup``
    """
    if not cdot or not coding_dna:
        return cdot
    cd = normalize_c_hgvs(cdot)
    body = cd[2:] if cd.lower().startswith('c.') else cd
    # 剪接/UTR 偏移不在此归一化
    if re.search(r'[+*\-]', body):
        return cd
    dna = coding_dna.upper()
    n = len(dna)

    # --- insertion: c.A_BinsSEQ ---
    m = re.match(r'^(\d+)_(\d+)ins([ACGTN]+)$', body, re.IGNORECASE)
    if m:
        left, right = int(m.group(1)), int(m.group(2))
        ins = m.group(3).upper()
        if right != left + 1 or left < 1 or left > n:
            return cd
        # 3′ roll：若插入点 3′ 碱基 == 插入序列首碱基，则右移并旋转
        from collections import deque
        buf = deque(ins)
        pos = left  # 插在 pos 之后；3′ 侧碱基 0-based 下标为 pos
        while pos < n and buf and dna[pos] == buf[0]:
            buf.append(buf.popleft())
            pos += 1
        ins2 = ''.join(buf)
        left2 = pos
        dup = _ins_to_dup_if_tandem(dna, left2, ins2)
        if dup:
            return dup
        return f"c.{left2}_{left2 + 1}ins{ins2}"

    # --- deletion: c.Adel / c.A_Bdel ---
    m = re.match(r'^(\d+)(?:_(\d+))?del([ACGTN]*)$', body, re.IGNORECASE)
    if m:
        start = int(m.group(1))
        end = int(m.group(2) or start)
        delseq = (m.group(3) or '').upper()
        if start < 1 or end < start or end > n:
            return cd
        obs = dna[start - 1:end]
        if delseq and delseq != obs:
            return cd  # 与本地 CDS 不一致时不强行滚动
        delseq = obs
        # 3′ roll deletion
        while end < n and delseq and dna[end] == delseq[0]:
            delseq = delseq[1:] + dna[end]
            start += 1
            end += 1
        if start == end:
            return f"c.{start}del{delseq}"
        return f"c.{start}_{end}del{delseq}"

    # --- duplication: c.Adup / c.A_Bdup ---
    m = re.match(r'^(\d+)(?:_(\d+))?dup([ACGTN]*)$', body, re.IGNORECASE)
    if m:
        start = int(m.group(1))
        end = int(m.group(2) or start)
        dupseq = (m.group(3) or '').upper()
        if start < 1 or end < start or end > n:
            return cd
        obs = dna[start - 1:end]
        if dupseq and dupseq != obs:
            return cd
        # dup ≡ 在 end 之后插入 obs；按插入做 3′ 对齐后再写回 dup
        from collections import deque
        buf = deque(obs)
        pos = end
        while pos < n and buf and dna[pos] == buf[0]:
            buf.append(buf.popleft())
            pos += 1
        ins2 = ''.join(buf)
        left2 = pos
        dup = _ins_to_dup_if_tandem(dna, left2, ins2)
        return dup if dup else _format_dup_cdot(left2 - len(ins2) + 1, left2)

    return cd


def _is_inframe_len(n: int) -> bool:
    """净长度变化是否为 3 的倍数（整码）。"""
    return n % 3 == 0


def detect_duplication_cdot(cdot: str, coding_dna: str) -> Optional[Tuple[str, Tuple[int, int], str]]:
    """兼容旧接口：委托 HGVS 3′ 归一化；若结果为 dup 则解析区间。"""
    norm = normalize_coding_indel_hgvs(cdot, coding_dna)
    m = re.match(r'^c\.(\d+)(?:_(\d+))?dup([ACGTN]*)$', norm, re.IGNORECASE)
    if not m:
        return None
    a = int(m.group(1))
    b = int(m.group(2) or a)
    seq = (m.group(3) or coding_dna[a - 1:b]).upper()
    return norm, (a, b), seq

class CodingChange:
    """编码变化分析器"""
    
    def __init__(self, evffile: str, genefile: str, fastafile: str, **kwargs):
        self.evffile = evffile
        self.genefile = genefile
        self.fastafile = fastafile
        # 平台 polish 会传：includesnp / alltranscript / newevf / outfile
        # -includesnp：处理 SNV；默认仅 indel（与 MATCHVAR coding_change 一致）。平台 polish 会传 True。
        self.includesnp = kwargs.get('includesnp', False)
        self.alltranscript = kwargs.get('alltranscript', False)
        self.newevf = kwargs.get('newevf')
        self.outfile = kwargs.get('outfile')
        # 独立 CLI 可选：onlyAltering 跳过同义；mrnaseq/codingseq 预留序列输出
        self.mrnaseq = kwargs.get('mrnaseq', False)
        self.onlyAltering = kwargs.get('onlyAltering', False)
        self.codingseq = kwargs.get('codingseq', False)
        self.tolerate = kwargs.get('tolerate', False)
        self.verbose = kwargs.get('verbose', False)
        self.man = kwargs.get('man', False)
        
        # 内部变量
        self.queue = []
        self.need_trans = set()
        self.flagged_transcript = set()
        self.mrnastart = {}
        self.mrnaend = {}
        self.mrna_sequences = {}
        self.newevf_p: Dict[Tuple[str, str], str] = {}
        self.newevf_function: Dict[Tuple[str, str], str] = {}
        # 记录修正后的 c.HGVS（用于回写 EVF 覆盖错误的 c.）
        self.newevf_c: Dict[Tuple[str, str], str] = {}
        # 基因/转录本元信息：用于 g→c 推导（负链互补等）
        self._gene_meta: Dict[str, Dict[str, object]] = {}
        
        # 验证参数
        self._validate_arguments()
    
    def _validate_arguments(self):
        """验证参数"""
        if self.codingseq and not self.mrnaseq:
            raise ValueError("Error in argument: --mrnaseq is required when --codingseq is specified")
        
        if self.newevf and not self.alltranscript:
            raise ValueError("Error in argument: --alltranscript arguments are required when you specify -newevf")
    
    def run_analysis(self):
        """运行编码变化分析"""
        try:
            # 读取外显子变异功能文件
            self._read_evf_file()
            
            # 读取基因定义文件
            self._read_gene_file()
            
            # 读取FASTA文件
            self._read_fasta_file()
            
            # 处理队列中的变异
            self._process_variants()
            
            # 输出结果
            self._write_output()
            
            logger.info("编码变化分析完成")
            
        except Exception as e:
            logger.error(f"编码变化分析失败: {e}")
            raise
    
    def _read_evf_file(self):
        """读取外显子变异功能文件"""
        logger.info(f"读取外显子变异功能文件: {self.evffile}")
        
        try:
            with open(self.evffile, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    
                    if not line.startswith('line'):
                        logger.warning(f"跳过无效记录: {line}")
                        continue
                    
                    parts = line.split('\t')
                    if len(parts) < 6:
                        continue
                    
                    line_num = parts[0]
                    function = parts[1]
                    annotation = parts[2]
                    chrom = parts[3]
                    start = int(parts[4])
                    end = int(parts[5])
                    ref = parts[6] if len(parts) > 6 else ''
                    alt = parts[7] if len(parts) > 7 else ''
                    
                    # 跳过影响整个基因的变异
                    if 'wholegene' in annotation:
                        continue
                    
                    # 跳过标记为unknown的变异
                    if 'unknown' in annotation.lower():
                        continue

                    # 无 -includesnp 时跳过 SNV（长度相等的单碱基替换）
                    if not self.includesnp:
                        ref_n = '' if ref in ('-', '.', '*') else ref
                        alt_n = '' if alt in ('-', '.', '*') else alt
                        if len(ref_n) == 1 and len(alt_n) == 1:
                            continue
                    
                    # 解析注释中的转录本信息（含c.与可选p.）
                    items = self._parse_annotation_transcripts(annotation)
                    for it in items:
                        transcript = it['transcript']
                        self.need_trans.add(transcript)
                        self.queue.append({
                            'line_num': line_num,
                            'function': function,
                            'annotation': annotation,
                            'chrom': chrom,
                            'start': start,
                            'end': end,
                            'ref': ref,
                            'alt': alt,
                            'transcript': transcript,
                            'cchange': it.get('c', ''),
                            'p_from_ann': it.get('p', '')
                        })
        
        except Exception as e:
            logger.error(f"读取外显子变异功能文件失败: {e}")
            raise
    
    def _parse_annotation_transcripts(self, annotation: str) -> List[Dict[str, str]]:
        """解析注释中的转录本信息，返回 [{transcript, c, p}]
        兼容两种格式：
        1) GENE:TRANSCRIPT:exonX:c.xxx[:p.yyy]
        2) TRANSCRIPT:exonX:c.xxx[:p.yyy]
        """
        items: List[Dict[str, str]] = []
        patterns = [
            # 带基因名：宽松捕获c.，忽略并吞掉任意p.后缀到逗号结束
            r'^[^\t,]*?([\w\-\.\@\/]+?):([\w\.\-]+?):exon\d+:(c\.[^:,\s]+)(?::p\.[^,]+)?$',
            # 仅转录本
            r'^[^\t,]*?([\w\.\-]+):exon\d+:(c\.[^:,\s]+)(?::p\.[^,]+)?$'
        ]
        for item in [x for x in annotation.split(',') if x]:
            item = item.strip()
            m = None
            # 先尝试基因+转录本形式
            m = re.match(patterns[0], item)
            if m:
                items.append({'transcript': m.group(2), 'c': m.group(3).replace('c.', ''), 'p': ''})
                continue
            # 再尝试仅转录本形式
            m = re.match(patterns[1], item)
            if m:
                items.append({'transcript': m.group(1), 'c': m.group(2).replace('c.', ''), 'p': ''})
        return items
    
    def _read_gene_file(self):
        """读取基因定义文件"""
        logger.info(f"读取基因定义文件: {self.genefile}")
        
        try:
            with open(self.genefile, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    
                    parts = line.split('\t')
                    if len(parts) < 15:
                        continue
                    
                    name = parts[1]
                    base_name = name.split('.')[0]
                    chrom = parts[2]
                    strand = parts[3]
                    txStart = int(parts[4])
                    txEnd = int(parts[5])
                    cdsStart = int(parts[6])
                    cdsEnd = int(parts[7])
                    exonStarts = parts[9]
                    exonEnds = parts[10]

                    try:
                        # 转为列表与1-based坐标
                        exon_start_list = [int(x) for x in exonStarts.rstrip(',').split(',') if x]
                        exon_end_list = [int(x) for x in exonEnds.rstrip(',').split(',') if x]
                        exon_start_1 = [x + 1 for x in exon_start_list]
                        exon_end_1 = [x for x in exon_end_list]  # 作为1-based闭区端
                        txStart1 = txStart + 1
                        cdsStart1 = cdsStart + 1

                        mrna_start = None
                        mrna_end = None
                        if strand == '+':
                            intron = 0
                            for i in range(len(exon_start_1)):
                                if i > 0:
                                    intron += (exon_start_1[i] - exon_end_1[i-1] - 1)
                                if exon_start_1[i] <= cdsStart1 <= exon_end_1[i]:
                                    mrna_start = cdsStart1 - txStart1 + 1 - intron
                                if exon_start_1[i] <= cdsEnd <= exon_end_1[i]:
                                    mrna_end = cdsEnd - txStart1 + 1 - intron
                        else:
                            intron = 0
                            for i in range(len(exon_start_1)-1, -1, -1):
                                if i < len(exon_start_1) - 1:
                                    intron += (exon_start_1[i+1] - exon_end_1[i] - 1)
                                if exon_start_1[i] <= cdsEnd <= exon_end_1[i]:
                                    mrna_start = txEnd - cdsEnd + 1 - intron
                                if exon_start_1[i] <= cdsStart1 <= exon_end_1[i]:
                                    mrna_end = txEnd - cdsStart1 + 1 - intron

                        # 同时为含版本号与无版本号的转录本ID建立映射，避免键不匹配
                        ms = int(mrna_start) if mrna_start else None
                        me = int(mrna_end) if mrna_end else None
                        self.mrnastart[name] = ms
                        self.mrnaend[name] = me
                        self.mrnastart[base_name] = ms
                        self.mrnaend[base_name] = me

                        # 预计算编码片段（1-based闭区间），用于 g→c 位置映射
                        # UCSC genePred：正负链均为 cdsStart < cdsEnd（0-based 半开）
                        coding_segments = []
                        cdsStart1 = cdsStart + 1
                        cdsEnd1 = cdsEnd  # 0-based exclusive end ≡ 1-based inclusive last CDS base
                        for s0, e0 in zip(exon_start_list, exon_end_list):
                            s1 = s0 + 1
                            e1 = e0
                            cs = max(s1, cdsStart1)
                            ce = min(e1, cdsEnd1)
                            if cs <= ce:
                                coding_segments.append((cs, ce))
                        # 统一按转录本5'→3'方向排列
                        if strand == '+':
                            coding_segments.sort(key=lambda x: x[0])
                        else:
                            # 负链时，转录本5'端在基因组高坐标，按起点降序排列
                            coding_segments.sort(key=lambda x: x[0], reverse=True)

                        meta = {
                            'chrom': chrom,
                            'strand': strand,
                            'coding_segments': coding_segments,
                            'cds_len': sum(ce - cs + 1 for cs, ce in coding_segments)
                        }
                        self._gene_meta[name] = meta
                        self._gene_meta[base_name] = meta
                    except Exception as e:
                        logger.warning(f"计算CDS mRNA坐标失败 {name}: {e}")
                        self.mrnastart[name] = None
                        self.mrnaend[name] = None
                        self.mrnastart[base_name] = None
                        self.mrnaend[base_name] = None
        
        except Exception as e:
            logger.error(f"读取基因定义文件失败: {e}")
            raise
    
    def _read_fasta_file(self):
        """读取FASTA文件"""
        logger.info(f"读取FASTA文件: {self.fastafile}")
        
        try:
            current_transcript = None
            current_sequence = ""
            
            with open(self.fastafile, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    
                    if line.startswith('>'):
                        # 保存前一个转录本的序列
                        if current_transcript and current_sequence:
                            self.mrna_sequences[current_transcript] = current_sequence
                        
                        # 开始新的转录本
                        # 仅取ID的第一个token，避免尾随描述导致ID不匹配
                        header = line[1:].strip()
                        token = header.split()[0]
                        # 去掉版本号以匹配注释中的转录本ID（NM_xxx.1 -> NM_xxx）
                        current_transcript = token.split('.')[0]
                        current_sequence = ""
                    else:
                        current_sequence += line
            
            # 保存最后一个转录本的序列
            if current_transcript and current_sequence:
                self.mrna_sequences[current_transcript] = current_sequence
        
        except Exception as e:
            logger.error(f"读取FASTA文件失败: {e}")
            raise
    
    def _process_variants(self):
        """处理变异队列"""
        logger.info(f"处理 {len(self.queue)} 个变异")
        
        for variant in self.queue:
            try:
                self._process_single_variant(variant)
            except Exception as e:
                logger.warning(f"处理变异失败 {variant['line_num']}: {e}")
                continue
    
    def _process_single_variant(self, variant: Dict):
        """处理单个变异"""
        # 去掉版本号以匹配FASTA键
        transcript = variant['transcript'].split('.')[0]
        
        # 检查转录本是否存在
        if transcript not in self.mrna_sequences:
            logger.warning(f"转录本 {transcript} 不存在于FASTA文件中")
            return
        
        # 获取mRNA序列
        mrna_seq = self.mrna_sequences[transcript]
        
        # 获取CDS位置
        if transcript not in self.mrnastart or transcript not in self.mrnaend:
            logger.warning(f"转录本 {transcript} 缺少CDS位置信息")
            return
        
        cds_start = self.mrnastart[transcript]
        cds_end = self.mrnaend[transcript]
        
        # 分析变异对蛋白质的影响
        protein_change = self._analyze_protein_change(variant, mrna_seq, cds_start, cds_end)

        # -onlyAltering：跳过同义（不写回 EVF）
        if self.onlyAltering and isinstance(protein_change, dict):
            eff = (protein_change.get('effect') or '').lower()
            if eff.startswith('synonymous'):
                return
        
        # 更新变异信息
        variant['protein_change'] = protein_change
        # 记录用于回写EVF的p.与功能
        try:
            line_num = variant.get('line_num', '')
            transcript = variant.get('transcript', '')
            p_hgvs = protein_change.get('p_hgvs', 'p.?') if isinstance(protein_change, dict) else 'p.?'
            effect = protein_change.get('effect', '') if isinstance(protein_change, dict) else ''
            if line_num and transcript:
                self.newevf_p[(line_num, transcript)] = p_hgvs.replace('p.', '')  # 存三字母，不带前缀，回写时加p.
                self.newevf_function[(line_num, transcript)] = effect
                # 若本次计算修正了 c.，也记录以便回写
                if 'corrected_cdot' in variant:
                    self.newevf_c[(line_num, transcript)] = variant['corrected_cdot']
        except Exception:
            pass
    
    def _analyze_protein_change(self, variant: Dict, mrna_seq: str, cds_start: int, cds_end: int) -> Dict:
        """分析蛋白质变化"""
        try:
            # 根据染色体选择密码子表
            chrom = variant.get('chrom', '')
            if chrom.upper() in ['M', 'MT', 'CHRM', 'CHRMT']:
                codon_table = CODON_TABLE_MT
                if self.verbose:
                    logger.info(f"使用线粒体密码子表处理染色体: {chrom}")
            else:
                codon_table = CODON_TABLE
            
            # 优先使用（或修正）c.HGVS 计算 p.，必要时根据基因组坐标推导 c.（SNV）
            cchange = variant.get('cchange', '')
            coding_dna = mrna_seq[cds_start-1:cds_end]
            if not coding_dna or len(coding_dna) < 3:
                return {'type': 'unknown', 'effect': 'unknown', 'p_hgvs': 'p.?'}

            def _valid_cdot(cdot: str) -> bool:
                """宽松校验 c.HGVS，支持 SNV/ins/del/delins/dup（含区间）。"""
                if not cdot:
                    return False
                cd = normalize_c_hgvs(cdot)
                cd = cd.replace('c.', '')
                patterns = [
                    r'^(\d+)[ACGTN]?>[ACGTN]$',                          # SNV: 180A>G or C180G（normalize 后）
                    r'^(\d+)del([ACGTN]*)$',                              # del / delSEQ
                    r'^(\d+)_([\d]+)del([ACGTN]*)$',                     # N_Mdel / N_MdelSEQ
                    r'^(\d+)_([\d]+)ins([ACGTN]+)$',                     # N_MinsSEQ
                    r'^(\d+)delins([ACGTN]+)$',                           # NdelinsSEQ
                    r'^(\d+)_([\d]+)delins([ACGTN]+)$',                  # N_MdelinsSEQ
                    r'^(\d+)dup([ACGTN]*)$',                              # Ndup[SEQ]
                    r'^(\d+)_([\d]+)dup([ACGTN]*)$'                      # N_Mdup[SEQ]
                ]
                for p in patterns:
                    if re.match(p, cd, re.IGNORECASE):
                        return True
                return False

            cdot_to_use = None
            if cchange and _valid_cdot(cchange):
                # 原始 c. 合法，但需要与基因组推导的 c. 交叉校验（位置/参考碱基），不一致则改用推导值
                orig_norm = normalize_c_hgvs(cchange)
                cdot_to_use = orig_norm
                if len(variant.get('ref', '')) == 1 and len(variant.get('alt', '')) == 1:
                    inferred = self._infer_c_from_genome_snv(variant)
                    if inferred and _valid_cdot(inferred):
                        # 对比位置和参考碱基
                        mo = re.match(r'^c\.(\d+)([ACGTN])>([ACGTN])$', orig_norm, re.IGNORECASE)
                        mi = re.match(r'^(?:c\.)?(\d+)([ACGTN])>([ACGTN])$', inferred, re.IGNORECASE)
                        if mo and mi:
                            pos_o, ref_o, alt_o = int(mo.group(1)), mo.group(2), mo.group(3)
                            pos_i, ref_i, alt_i = int(mi.group(1)), mi.group(2), mi.group(3)
                            # 若参考碱基不匹配或位置差异明显，则以推导为准
                            if ref_o.upper() != ref_i.upper() or pos_o != pos_i:
                                cdot_to_use = f"c.{pos_i}{ref_i}>{alt_i}"
                                variant['corrected_cdot'] = cdot_to_use
            else:
                if len(variant.get('ref', '')) == 1 and len(variant.get('alt', '')) == 1:
                    inferred = self._infer_c_from_genome_snv(variant)
                    if inferred and _valid_cdot(inferred):
                        # 规范化为带前缀的形式（避免在f-string中直接使用带反斜杠的正则）
                        m_norm = re.match(r'^(?:c\.)?(\d+[ACGTN]>[ACGTN])$', inferred, re.IGNORECASE)
                        if m_norm:
                            cdot_to_use = f"c.{m_norm.group(1)}"
                            variant['corrected_cdot'] = cdot_to_use
                else:
                    # INDEL：从基因组推导初步 c.，随后尝试 dup 归并
                    inferred_indel = self._infer_c_from_genome_indel(variant)
                    if inferred_indel:
                        cdot_to_use = normalize_c_hgvs(inferred_indel)
                        variant['corrected_cdot'] = cdot_to_use

            if cdot_to_use:
                # HGVS 3′ 规则归一化（重复区滚动）+ 串联插入优先写 dup
                try:
                    norm_c = normalize_coding_indel_hgvs(cdot_to_use, coding_dna)
                    if norm_c:
                        if norm_c != cdot_to_use:
                            cdot_to_use = norm_c
                        variant['corrected_cdot'] = cdot_to_use
                except Exception:
                    pass
                # intronic 插入（c.N+/-a_N+/-binsSEQ）不在此处归并，由 annotate_variation 负责；
                # 但若 EVF 给出 intronic ins，后续 p. 计算会忽略，不会生成错误的蛋白注释。
                mutated_cds, effect_hint = apply_c_hgvs_to_cds(coding_dna, cdot_to_use)
                if re.search(r'dup', cdot_to_use or '', re.I) and effect_hint != 'frameshift':
                    effect_hint = 'dup'
                wt_protein = translate_protein(coding_dna, chrom)
                mut_protein = translate_protein(mutated_cds, chrom)
                # SNV 强制短路：若新密码子为终止子，则直接输出 stopgain（避免被diff误判为删除）
                m_snv = re.match(r'^c\.(\d+)[ACGTN]?>[ACGTN]$', cdot_to_use, re.IGNORECASE)
                if m_snv:
                    pos_nt = int(m_snv.group(1))
                    codon_idx = (pos_nt - 1) // 3
                    wt_codon = coding_dna[codon_idx*3: codon_idx*3+3]
                    mut_codon = mutated_cds[codon_idx*3: codon_idx*3+3]
                    # 须用当前染色体对应表（核 / 线粒体），勿写死 CODON_TABLE
                    wt_aa = codon_table.get(wt_codon.upper(), 'X')
                    mut_aa = codon_table.get(mut_codon.upper(), 'X')
                    if mut_aa == '*':
                        prot_pos = codon_idx + 1
                        p_hgvs = f"p.{three_letter(wt_aa)}{prot_pos}Ter"
                        effect = 'stopgain'
                        return {'type': effect, 'effect': effect, 'p_hgvs': p_hgvs}
                p_hgvs, eff = format_p_hgvs_from_diff(wt_protein, mut_protein, chrom, effect_hint)
                effect = eff
                # 规范 frameshift 命名：明确 fs*X 位点
                if effect_hint == 'frameshift':
                    # 若 format 已给出 fs*X，保持；否则计算新终止子位置
                    if 'fs*' not in (p_hgvs or ''):
                        try:
                            # 从差异起点向后找 '*' 终止子（pos1 为 1-based）
                            pos1, _, _, _ = diff_proteins(wt_protein, mut_protein)
                            start_idx = max(0, pos1 - 1)
                            tail = mut_protein[start_idx:]
                            stop_idx = tail.find('*')
                            fs_len = stop_idx + 1 if stop_idx >= 0 else max(0, len(mut_protein) - start_idx)
                            aa_new = (
                                mut_protein[start_idx]
                                if start_idx < len(mut_protein) and mut_protein[start_idx] != '*'
                                else '?'
                            )
                            p_hgvs = (
                                f"p.{three_letter(wt_protein[start_idx])}{pos1}"
                                f"{three_letter(aa_new)}fs*{fs_len}"
                            )
                        except Exception:
                            pass
                    effect = 'frameshift'
                return {'type': effect.replace('_', ' ') if effect else 'unknown', 'effect': effect or 'unknown', 'p_hgvs': p_hgvs}

            # 回退：无可用 c. 时，用基因模型从基因组推导 c. 再算 p.（勿用 start-cds_start 当 mRNA 坐标）
            inferred = self._infer_c_from_genome_snv(variant) or self._infer_c_from_genome_indel(variant)
            if inferred:
                try:
                    cdot = normalize_coding_indel_hgvs(normalize_c_hgvs(inferred), coding_dna)
                    mutated_cds, effect_hint = apply_c_hgvs_to_cds(coding_dna, cdot)
                    wt_protein = translate_protein(coding_dna, chrom)
                    mut_protein = translate_protein(mutated_cds, chrom)
                    p_hgvs, eff = format_p_hgvs_from_diff(wt_protein, mut_protein, chrom, effect_hint)
                    variant['corrected_cdot'] = cdot
                    return {
                        'type': (eff or 'unknown').replace('_', ' '),
                        'effect': eff or 'unknown',
                        'p_hgvs': p_hgvs,
                    }
                except Exception:
                    pass
            cds_pos0 = self._genome_to_cds_index0(variant)
            ref = variant['ref']
            alt = variant['alt']
            if cds_pos0 is None:
                return {'type': 'unknown', 'effect': 'unknown', 'p_hgvs': 'p.?'}
            if alt in ['-', '*', '.'] and len(ref) >= 1:
                return self._analyze_deletion(variant, coding_dna, cds_pos0, codon_table)
            if len(ref) == len(alt):
                if len(ref) == 1:
                    return self._analyze_substitution(variant, coding_dna, cds_pos0, codon_table)
                return self._analyze_block_substitution(variant, coding_dna, cds_pos0, codon_table)
            if len(ref) > len(alt):
                return self._analyze_deletion(variant, coding_dna, cds_pos0, codon_table)
            return self._analyze_insertion(variant, coding_dna, cds_pos0, codon_table)
        
        except Exception as e:
            if self.tolerate:
                logger.warning(f"分析蛋白质变化失败: {e}")
                return {'type': 'unknown', 'effect': 'unknown'}
            else:
                raise

    def _genome_to_cds_index0(self, variant: Dict) -> Optional[int]:
        """基因组坐标 → CDS 0-based 下标（用于粗略回退）；失败返回 None。"""
        transcript = variant.get('transcript')
        try:
            gpos = int(variant.get('start'))
        except Exception:
            return None
        if not transcript:
            return None
        meta = self._gene_meta.get(transcript) or self._gene_meta.get(str(transcript).split('.')[0])
        if not meta:
            return None
        chrom = variant.get('chrom')
        if meta.get('chrom') and str(meta['chrom']) != str(chrom):
            return None
        strand = meta['strand']
        traversed = 0
        for (s, e) in meta['coding_segments']:
            if s <= gpos <= e:
                if strand == '+':
                    offset = gpos - s + 1
                else:
                    offset = e - gpos + 1
                return traversed + offset - 1  # 0-based
            traversed += (e - s + 1)
        return None

    def _infer_c_from_genome_snv(self, variant: Dict) -> Optional[str]:
        """基于基因组坐标与基因模型，将SNV推导为 c.HGVS（考虑正负链与外显子拼接）。"""
        transcript = variant.get('transcript')
        chrom = variant.get('chrom')
        try:
            gpos = int(variant.get('start'))
        except Exception:
            return None
        ref = (variant.get('ref') or '').upper()
        alt = (variant.get('alt') or '').upper()
        meta = self._gene_meta.get(transcript) or self._gene_meta.get(
            str(transcript).split('.')[0] if transcript else ''
        )
        if not meta:
            return None
        if meta.get('chrom') and str(meta['chrom']) != str(chrom):
            return None
        strand = meta['strand']
        segments: List[Tuple[int, int]] = meta['coding_segments']
        # g→c 映射（按转录本方向累计，负链片段内偏移需用 e-gpos+1）
        traversed = 0
        pos_in_cds = None
        for (s, e) in segments:
            if s <= gpos <= e:
                if strand == '+':
                    offset = gpos - s + 1
                else:
                    offset = e - gpos + 1
                pos_in_cds = traversed + offset
                break
            traversed += (e - s + 1)
        if pos_in_cds is None:
            return None
        if strand == '-':
            ref = reverse_complement(ref)
            alt = reverse_complement(alt)
        return f"c.{pos_in_cds}{ref}>{alt}"

    def _infer_c_from_genome_indel(self, variant: Dict) -> Optional[str]:
        """基于基因组坐标与基因模型，将INDEL推导为 c.HGVS（初步表达：ins/del/delins），
        后续会由 detect_duplication_cdot/normalize 进一步归并（如 dup）和左对齐。
        """
        transcript = variant.get('transcript')
        chrom = variant.get('chrom')
        try:
            gpos = int(variant.get('start'))
        except Exception:
            return None
        raw_ref = variant.get('ref') or ''
        raw_alt = variant.get('alt') or ''
        # MATCHVAR 占位符：'-'/'.'/ '*' 表示空等位基因
        ref_g = '' if raw_ref in ('-', '.', '*') else raw_ref.upper()
        alt_g = '' if raw_alt in ('-', '.', '*') else raw_alt.upper()
        meta = self._gene_meta.get(transcript) or self._gene_meta.get(
            str(transcript).split('.')[0] if transcript else ''
        )
        if not meta:
            return None
        if meta.get('chrom') and str(meta['chrom']) != str(chrom):
            return None
        strand = meta['strand']
        segments: List[Tuple[int, int]] = meta['coding_segments']
        # g→c：取插入/缺失发生处的 cDNA 基座位置（负链片段内偏移用 e-gpos+1）
        traversed = 0
        pos_in_cds = None
        for (s, e) in segments:
            if s <= gpos <= e:
                if strand == '+':
                    offset = gpos - s + 1
                else:
                    offset = e - gpos + 1
                pos_in_cds = traversed + offset
                break
            traversed += (e - s + 1)
        if pos_in_cds is None:
            return None
        # 规范 indel：插入用 c.N_(N+1)ins（切勿 delins，否则 apply_c_hgvs 会误删 1bp）
        if len(alt_g) > len(ref_g):
            ins_seq = alt_g[len(ref_g):] if ref_g else alt_g
            if strand == '-':
                ins_seq = reverse_complement(ins_seq)
            if strand == '+':
                return f"c.{pos_in_cds}_{pos_in_cds + 1}ins{ins_seq}"
            # 负链：基因组左对齐插入点对应转录 3′ 侧翼碱基
            return f"c.{pos_in_cds - 1}_{pos_in_cds}ins{ins_seq}"
        # 缺失
        if len(ref_g) > len(alt_g):
            del_seq = ref_g[len(alt_g):] if alt_g else ref_g
            if strand == '-':
                del_seq = reverse_complement(del_seq)
            # 负链时，pos_in_cds 为区间右端，需向左回溯长度
            if strand == '+':
                start = pos_in_cds
                end = pos_in_cds + len(del_seq) - 1
            else:
                start = pos_in_cds - len(del_seq) + 1
                end = pos_in_cds
            if len(del_seq) == 1:
                return f"c.{start}del{del_seq}"
            return f"c.{start}_{end}del{del_seq}"
        # delins（等长替换）
        if len(ref_g) == len(alt_g) and len(ref_g) > 1:
            seq = alt_g
            if strand == '-':
                seq = reverse_complement(seq)
            if strand == '+':
                start = pos_in_cds
                end = pos_in_cds + len(ref_g) - 1
            else:
                start = pos_in_cds - len(ref_g) + 1
                end = pos_in_cds
            return f"c.{start}_{end}delins{seq}"
        return None
    
    # 主路径走 c.HGVS→p.；以下为无有效 c. / 推导失败时的粗略回退
    def _analyze_substitution(self, variant: Dict, coding_dna: str, cds_pos0: int, codon_table: Dict) -> Dict:
        """分析替换变异。cds_pos0 为 CDS 0-based 下标；等位基因按转录本方向。"""
        ref = (variant.get('ref') or '').upper()
        alt = (variant.get('alt') or '').upper()
        transcript = variant.get('transcript')
        strand = '+'
        if transcript and transcript in self._gene_meta:
            strand = self._gene_meta[transcript].get('strand', '+') or '+'
        if strand == '-':
            ref = reverse_complement(ref)
            alt = reverse_complement(alt)
        if cds_pos0 < 0 or cds_pos0 >= len(coding_dna):
            return {'type': 'unknown', 'effect': 'unknown', 'p_hgvs': 'p.?'}
        if ref and coding_dna[cds_pos0].upper() != ref:
            # 与本地 CDS 不一致时仍按坐标替换 alt（兼容占位/模糊 REF）
            pass
        codon_start = (cds_pos0 // 3) * 3
        original_codon = coding_dna[codon_start:codon_start + 3].upper()
        if len(original_codon) < 3:
            return {'type': 'unknown', 'effect': 'unknown', 'p_hgvs': 'p.?'}
        pos_in_codon = cds_pos0 % 3
        new_codon_list = list(original_codon)
        new_codon_list[pos_in_codon] = alt if alt else 'N'
        new_codon = ''.join(new_codon_list)
        original_aa = codon_table.get(original_codon, 'X')
        new_aa = codon_table.get(new_codon, 'X')
        protein_pos = (cds_pos0 // 3) + 1
        if original_aa == new_aa:
            p_hgvs = f"p.{three_letter(original_aa)}{protein_pos}=" if original_aa not in ('X',) else "p.?"
            return {'type': 'synonymous', 'effect': 'synonymous_SNV', 'p_hgvs': p_hgvs}
        if new_aa == '*':
            return {
                'type': 'stopgain', 'effect': 'stopgain',
                'p_hgvs': f"p.{three_letter(original_aa)}{protein_pos}Ter",
            }
        if original_aa == '*':
            return {
                'type': 'stoploss', 'effect': 'stoploss',
                'p_hgvs': f"p.Ter{protein_pos}{three_letter(new_aa)}",
            }
        return {
            'type': 'nonsynonymous', 'effect': 'nonsynonymous_SNV',
            'p_hgvs': f"p.{three_letter(original_aa)}{protein_pos}{three_letter(new_aa)}",
        }

    def _analyze_block_substitution(self, variant: Dict, coding_dna: str, cds_pos0: int, codon_table: Dict) -> Dict:
        """分析块替换变异（粗略：按长度判整码/移码）。"""
        ref = variant.get('ref') or ''
        alt = variant.get('alt') or ''
        delta = abs(len(ref) - len(alt))
        if _is_inframe_len(len(ref)) and _is_inframe_len(len(alt)) and _is_inframe_len(delta):
            return {'type': 'in-frame_substitution', 'effect': 'inframe_substitution', 'p_hgvs': 'p.?'}
        return {'type': 'frameshift_substitution', 'effect': 'frameshift_substitution', 'p_hgvs': 'p.?'}

    def _analyze_deletion(self, variant: Dict, coding_dna: str, cds_pos0: int, codon_table: Dict) -> Dict:
        """分析删除变异（粗略）。"""
        ref = variant.get('ref') or ''
        alt = variant.get('alt') or ''
        if alt in ('-', '.', '*'):
            alt = ''
        del_len = len(ref) - len(alt) if len(ref) >= len(alt) else len(ref)
        if _is_inframe_len(del_len):
            return {'type': 'in-frame_deletion', 'effect': 'nonframeshift_deletion', 'p_hgvs': 'p.?'}
        return {'type': 'frameshift_deletion', 'effect': 'frameshift_deletion', 'p_hgvs': 'p.?'}

    def _analyze_insertion(self, variant: Dict, coding_dna: str, cds_pos0: int, codon_table: Dict) -> Dict:
        """分析插入变异（粗略）。"""
        ref = variant.get('ref') or ''
        alt = variant.get('alt') or ''
        if ref in ('-', '.', '*'):
            ref = ''
        ins_len = len(alt) - len(ref) if len(alt) >= len(ref) else len(alt)
        if _is_inframe_len(ins_len):
            return {'type': 'in-frame_insertion', 'effect': 'nonframeshift_insertion', 'p_hgvs': 'p.?'}
        return {'type': 'frameshift_insertion', 'effect': 'frameshift_insertion', 'p_hgvs': 'p.?'}
    
    def _write_output(self):
        """输出结果"""
        if self.outfile:
            output_file = self.outfile
        else:
            output_file = f"{self.evffile}.coding_change"
        
        try:
            # 写辅助输出
            with open(output_file, 'w', encoding='utf-8') as f:
                for variant in self.queue:
                    if 'protein_change' in variant:
                        line = f"{variant['line_num']}\t{variant['function']}\t{variant['annotation']}\t"
                        line += f"{variant['chrom']}\t{variant['start']}\t{variant['end']}\t"
                        line += f"{variant['ref']}\t{variant['alt']}\t"
                        pc = variant['protein_change']
                        line += f"{pc.get('type','')}\t{pc.get('effect','')}\t{pc.get('p_hgvs','p.?')}\n"
                        f.write(line)

            # 如果指定 -newevf，则重写exonic_variant_function（逐条注释精确回写c.与p.，并必要时修正功能类别）
            if self.newevf:
                try:
                    with open(self.evffile, 'r', encoding='utf-8') as evf_in, open(self.newevf, 'w', encoding='utf-8') as nf:
                        for raw in evf_in:
                            line = raw.strip()
                            if not line:
                                continue
                            parts = line.split('\t')
                            if len(parts) < 8 or not parts[0].startswith('line'):
                                nf.write(raw)
                                continue
                            line_id = parts[0]
                            exonic_func = parts[1]
                            annot_field = parts[2].rstrip(',')
                            chrom, start, end, ref, alt = parts[3], parts[4], parts[5], parts[6], parts[7]

                            # 按逗号分割每个转录本注释，逐个回写c./p.
                            items = [x for x in annot_field.split(',') if x]
                            new_items: List[str] = []
                            for item in items:
                                # 尝试匹配两种格式
                                # 1) GENE:TRANSCRIPT:exonX:c.xxx[:p.yyy]
                                m1 = re.match(
                                    r'^([\w\-\.\@\/]+?):([\w\.\-]+?):(exon\d+):(c\.[^:\s,]+)(:p\.[^:\s,]+)?$',
                                    item,
                                )
                                # 2) TRANSCRIPT:exonX:c.xxx[:p.yyy]
                                m2 = re.match(
                                    r'^([\w\.\-]+):(exon\d+):(c\.[^:\s,]+)(:p\.[^:\s,]+)?$',
                                    item,
                                )
                                if m1:
                                    gene_name = m1.group(1)
                                    transcript = m1.group(2)
                                    exon_tag = m1.group(3)
                                    cdot = m1.group(4)
                                    p_new = self.newevf_p.get((line_id, transcript))
                                    c_new = self.newevf_c.get((line_id, transcript))
                                    cdot = c_new if c_new else cdot
                                    if p_new:
                                        item = f"{gene_name}:{transcript}:{exon_tag}:{cdot}:p.{p_new}"
                                elif m2:
                                    transcript = m2.group(1)
                                    exon_tag = m2.group(2)
                                    cdot = m2.group(3)
                                    p_new = self.newevf_p.get((line_id, transcript))
                                    c_new = self.newevf_c.get((line_id, transcript))
                                    cdot = c_new if c_new else cdot
                                    if p_new:
                                        item = f"{transcript}:{exon_tag}:{cdot}:p.{p_new}"
                                new_items.append(item)

                            new_annot = ','.join(new_items) + (',' if new_items else '')
                            # 若可据功能修正 exonic_func（参考 Perl 逻辑）
                            for item in items:
                                m = re.match(r'^([\w\-\.\@\/]+?):([\w\.\-]+?):exon\d+:', item) or re.match(r'^([\w\.\-]+):exon\d+:', item)
                                if not m:
                                    continue
                                t = m.group(2) if m.lastindex and m.lastindex >= 2 else m.group(1)
                                eff = self.newevf_function.get((line_id, t))
                                if not eff:
                                    continue
                                eff_l = eff.lower()
                                ef_l = exonic_func.lower()
                                if eff_l.startswith('synonymous'):
                                    if ef_l.startswith('frameshift'):
                                        exonic_func = 'non' + exonic_func
                                    else:
                                        exonic_func = 'synonymous SNV'
                                elif 'stopgain' in eff_l:
                                    if exonic_func != 'startloss':
                                        exonic_func = 'stopgain'
                                elif 'stoploss' in eff_l:
                                    if exonic_func not in ('startloss', 'stopgain'):
                                        exonic_func = 'stoploss'
                                elif 'startloss' in eff_l:
                                    exonic_func = 'startloss'
                                elif 'duplication' in eff_l:
                                    # HGVS 将串联插入规范为 dup 后，同步修正外显子功能类别
                                    exonic_func = (
                                        'frameshift_duplication'
                                        if ('frameshift' in eff_l and 'nonframeshift' not in eff_l)
                                        else 'nonframeshift_duplication'
                                    )
                                break

                            nf.write('\t'.join([line_id, exonic_func, new_annot, chrom, start, end, ref, alt]) + '\n')
                except Exception as e:
                    logger.warning(f"重写exonic_variant_function失败: {e}")
        
        except Exception as e:
            logger.error(f"写入输出文件失败: {e}")
            raise

def main():
    """主函数"""
    examples = (
        "示例:\n"
        "1) 基于EVF/基因/FASTA运行编码变化分析，并生成修正后的EVF：\n"
        "   python utils/matchvar/coding_change.py \\\n+        result.refGene.exonic_variant_function \\\n+        /Users/James/PycharmProjects/Variant_Data_Simulation_2.0/resources/humandb/hg19_refGene.txt \\\n+        /Users/James/PycharmProjects/Variant_Data_Simulation_2.0/resources/humandb/hg19_refGeneMrna.fa \\\n+        -alltranscript -newevf result.refGene.exonic_variant_function.fixed\n\n"
        "2) 仅输出辅助结果到指定文件：\n"
        "   python utils/matchvar/coding_change.py evf.txt hg19_refGene.txt hg19_refGeneMrna.fa -outfile out.coding_change\n\n"
        "3) 容忍错误并开启详细日志：\n"
        "   python utils/matchvar/coding_change.py evf.txt hg19_refGene.txt hg19_refGeneMrna.fa -tolerate -v\n"
    )
    parser = argparse.ArgumentParser(
        description='MATCHVAR编码变化分析工具',
        epilog=examples,
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('evffile', help='外显子变异功能文件')
    parser.add_argument('genefile', help='基因定义文件')
    parser.add_argument('fastafile', help='FASTA文件')
    parser.add_argument('-includesnp', action='store_true', help='包含SNP（平台 polish 会开启）')
    parser.add_argument('-alltranscript', action='store_true', help='所有转录本')
    parser.add_argument('-newevf', help='新的外显子变异功能文件')
    parser.add_argument('-outfile', help='输出文件')
    # 独立 CLI 可选（平台 polish 通常不传）
    parser.add_argument('-mrnaseq', action='store_true', help='mRNA序列')
    parser.add_argument('-onlyAltering', action='store_true', help='仅输出改变氨基酸的变异')
    parser.add_argument('-codingseq', action='store_true', help='编码序列')
    parser.add_argument('-tolerate', action='store_true', help='容忍错误')
    parser.add_argument('-verbose', '-v', action='store_true', help='详细输出')
    parser.add_argument('-man', '-m', action='store_true', help='显示手册')
    
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args()
    
    # 创建编码变化分析器
    analyzer = CodingChange(
        evffile=args.evffile,
        genefile=args.genefile,
        fastafile=args.fastafile,
        includesnp=args.includesnp,
        mrnaseq=args.mrnaseq,
        onlyAltering=args.onlyAltering,
        codingseq=args.codingseq,
        alltranscript=args.alltranscript,
        newevf=args.newevf,
        outfile=args.outfile,
        tolerate=args.tolerate,
        verbose=args.verbose,
        man=args.man
    )
    
    # 运行分析
    try:
        analyzer.run_analysis()
        logger.info("编码变化分析成功完成")
    except Exception as e:
        logger.error(f"编码变化分析过程中发生错误: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main() 