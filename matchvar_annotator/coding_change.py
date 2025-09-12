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
    """根据蛋白差异生成p.注释与effect。优先处理frameshift/stopgain等情况。"""
    pos1, aa1, aa2, pos2 = diff_proteins(wt, mut)
    if pos1 == 0:
        # 无氨基酸差异：HGVS 推荐 p.(=)
        return 'p.(=)', 'synonymous_SNV'
    # 若遇到'*'，判断stopgain/stoploss
    # 简化：当aa2含'*'且位置在差异段内，视为stopgain
    if '*' in aa2:
        return f'p.{three_letter(wt[pos1-1])}{pos1}Ter', 'stopgain'
    if '*' in wt[pos1-1:pos2]:  # 野生型差异窗口含终止
        return f'p.Ter{pos1}{three_letter(mut[pos1-1])}', 'stoploss'
    # frameshift 专用：若上层判定为frameshift，则计算至新终止子的距离
    if effect_hint == 'frameshift':
        # 终止位置：从差异起点在 mutated 蛋白中寻找 '*'
        start_idx = max(0, pos1 - 1)
        stop_idx = mut.find('*', start_idx)
        if stop_idx == -1:
            fs_len = max(0, len(mut) - start_idx)
        else:
            fs_len = stop_idx - start_idx
        new_aa = mut[start_idx] if start_idx < len(mut) and mut[start_idx] != '*' else '?'
        p = f'p.{three_letter(wt[pos1-1])}{pos1}{three_letter(new_aa)}fs*{fs_len}'
        return p, 'frameshift'
    # 同义变化（单点且AA相同）已在上层过滤；此处处理常见情况
    if len(aa1) == 1 and len(aa2) == 1:
        if aa1 == aa2:
            return f'p.{three_letter(aa1)}{pos1}{three_letter(aa2)}', 'synonymous_SNV'
        else:
            return f'p.{three_letter(aa1)}{pos1}{three_letter(aa2)}', 'nonsynonymous_SNV'
    # 指定为dup的非移码插入：优先输出 dup 语法
    if effect_hint == 'dup':
        # 以差异窗口为界，aa2 为被复制的氨基酸序列
        if len(aa2) >= 1:
            left = wt[pos1-1] if pos1-1 < len(wt) else ''
            right = wt[pos2-1] if pos2-1 < len(wt) and pos2 >= pos1 else left
            dup_seq = ''.join(three_letter(x) for x in aa2)
            if len(aa2) == 1:
                return f'p.{three_letter(aa2[0])}{pos1}dup', 'nonframeshift_duplication'
            else:
                return f'p.{three_letter(left)}{pos1}_{three_letter(right)}{pos2}dup{dup_seq}', 'nonframeshift_duplication'
    # 插入
    if len(mut) > len(wt) and len(aa2) > len(aa1):
        left = wt[pos1-1] if pos1-1 < len(wt) else ''
        right = wt[pos2-1] if pos2-1 < len(wt) and pos2 >= pos1 else left
        ins_seq = ''.join(three_letter(x) for x in aa2)
        return f'p.{three_letter(left)}{pos1}_{three_letter(right)}{pos2}ins{ins_seq}', 'nonframeshift_insertion'
    # 删除
    if len(mut) < len(wt) and len(aa1) > len(aa2):
        if len(aa1) == 1:
            return f'p.{three_letter(aa1)}{pos1}del', 'nonframeshift_deletion'
        else:
            return f'p.{three_letter(aa1[0])}{pos1}_{three_letter(aa1[-1])}{pos2}del', 'nonframeshift_deletion'
    # 替换/复合
    return f'p.{three_letter(aa1[0])}{pos1}_{three_letter(aa1[-1])}{pos2}delins' + ''.join(three_letter(x) for x in aa2), 'inframe_substitution'

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
    # 4) dup: NdupSEQ / N_MdupSEQ（重复）
    m = re.match(r'^(\d+)dup([ACGTN]+)$', c, re.IGNORECASE)
    if m:
        pos = int(m.group(1)); seq = m.group(2).upper()
        return s[:pos] + seq + s[pos:], 'dup' if len(seq) % 3 == 0 else 'frameshift'
    m = re.match(r'^(\d+)_([\d]+)dup([ACGTN]+)$', c, re.IGNORECASE)
    if m:
        start = int(m.group(1)); end = int(m.group(2)); seq = m.group(3).upper()
        # 在区间末端后插入重复序列
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

def _ceil_div3(n: int) -> int:
    return (n + 2) // 3

def _floor_div3(n: int) -> int:
    return n // 3

def detect_duplication_cdot(cdot: str, coding_dna: str) -> Optional[Tuple[str, Tuple[int, int], str]]:
    """检测 c.HGVS 是否可归并为 dup，并返回 (dup_cdot, (n1, n2), nt_seq)。
    规则（改进版）：
    - 对 c.N_MinsSEQ：若 coding_dna[(N-len(seq)) : N] == seq，则视为复制了 [N-len(seq)+1, N] 区间 → c.(N-len+1)_(N)dupSEQ
    - 对 c.N_MdelinsSEQ：若删除长度为0或插入序列等于左侧相邻序列，同上（作为 dup）。
    - 增加更严格的重复检测条件，避免将简单插入错误识别为重复
    """
    if not cdot:
        return None
    cd = normalize_c_hgvs(cdot)
    cd_nopre = cd.replace('c.', '')
    
    # ins 情形（尝试严格匹配，失败则尝试左移或修剪尾碱基以获得最简重复）
    m = re.match(r'^(\d+)_(\d+)ins([ACGTN]+)$', cd_nopre, re.IGNORECASE)
    if m:
        left = int(m.group(1))  # 插入点左侧核苷酸位置（1-based）
        seq = m.group(3).upper(); k = len(seq)
        
        # 增加重复检测的严格条件
        # 1. 序列长度必须大于等于3个碱基
        if k < 3:
            return None
            
        # 2. 检查是否为简单的重复模式（如AAAA, TTTT等）
        if len(set(seq)) == 1:
            # 单碱基重复，需要更严格的匹配
            if left >= k * 2 and coding_dna[left - k*2:left] == seq + seq:
                return f"c.{left - k*2 + 1}_{left}dup{seq}", (left - k*2 + 1, left), seq
            return None
        
        # 3. 精确匹配：左侧k碱基等于插入序列
        if left >= k and coding_dna[left - k:left] == seq:
            # 额外检查：确保不是简单的插入
            # 如果插入序列在基因组中频繁出现，可能是简单插入
            seq_count = coding_dna.count(seq)
            if seq_count > 5:  # 如果序列在CDS中出现超过5次，可能是简单插入
                return None
            return f"c.{left - k + 1}_{left}dup{seq}", (left - k + 1, left), seq
        
        # 4. 尝试修剪尾部（常见多写一碱基）：seq[:-1]
        if k > 1 and left >= (k - 1) and coding_dna[left - (k - 1):left] == seq[:-1]:
            seq2 = seq[:-1]; k2 = k - 1
            # 同样检查修剪后的序列
            if len(set(seq2)) == 1 and k2 < 3:
                return None
            return f"c.{left - k2 + 1}_{left}dup{seq2}", (left - k2 + 1, left), seq2
        
        # 5. 尝试左移一位匹配
        if left - 1 >= k and coding_dna[left - 1 - k:left - 1] == seq:
            # 检查左移后的匹配是否合理
            if len(set(seq)) == 1 and k < 3:
                return None
            return f"c.{left - 1 - k + 1}_{left - 1}dup{seq}", (left - k, left - 1), seq
        return None
    
    # delins 情形（尽量识别成等效dup）
    m = re.match(r'^(\d+)_(\d+)delins([ACGTN]+)$', cd_nopre, re.IGNORECASE)
    if m:
        start = int(m.group(1)); end = int(m.group(2)); seq = m.group(3).upper(); k = len(seq)
        anchor = start  # 插入位置视为 start 左侧
        
        # 增加重复检测的严格条件
        if k < 3:
            return None
            
        if len(set(seq)) == 1 and k < 3:
            return None
            
        if anchor >= k and coding_dna[anchor - k:anchor] == seq:
            # 检查序列在CDS中的出现频率
            seq_count = coding_dna.count(seq)
            if seq_count > 5:
                return None
            return f"c.{anchor - k + 1}_{anchor}dup{seq}", (anchor - k + 1, anchor), seq
        
        if k > 1 and anchor >= (k - 1) and coding_dna[anchor - (k - 1):anchor] == seq[:-1]:
            seq2 = seq[:-1]; k2 = k - 1
            if len(set(seq2)) == 1 and k2 < 3:
                return None
            return f"c.{anchor - k2 + 1}_{anchor}dup{seq2}", (anchor - k2 + 1, anchor), seq2
        
        if anchor - 1 >= k and coding_dna[anchor - 1 - k:anchor - 1] == seq:
            if len(set(seq)) == 1 and k < 3:
                return None
            return f"c.{anchor - k}_{anchor - 1}dup{seq}", (anchor - k, anchor - 1), seq
    return None

class CodingChange:
    """编码变化分析器"""
    
    def __init__(self, evffile: str, genefile: str, fastafile: str, **kwargs):
        self.evffile = evffile
        self.genefile = genefile
        self.fastafile = fastafile
        self.includesnp = kwargs.get('includesnp', False)
        self.mrnaseq = kwargs.get('mrnaseq', False)
        self.onlyAltering = kwargs.get('onlyAltering', False)
        self.codingseq = kwargs.get('codingseq', False)
        self.alltranscript = kwargs.get('alltranscript', False)
        self.newevf = kwargs.get('newevf')
        self.outfile = kwargs.get('outfile')
        self.tolerate = kwargs.get('tolerate', False)
        
        # 添加新的重要参数
        self.verbose = kwargs.get('verbose', False)  # 详细输出
        self.man = kwargs.get('man', False)  # 手册
        
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
                        coding_segments = []
                        cdsStart1 = cdsStart + 1
                        cdsEnd1 = cdsEnd
                        for s0, e0 in zip(exon_start_list, exon_end_list):
                            s1 = s0 + 1
                            e1 = e0
                            if strand == '-':
                                # 对于负链基因: cdsStart > cdsEnd
                                cs = max(s1, cdsEnd1)  # 使用cdsEnd作为下界
                                ce = min(e1, cdsStart1)  # 使用cdsStart作为上界
                            else:
                                # 对于正链基因: cdsStart < cdsEnd
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
                    r'^(\d+)dup([ACGTN]+)$',                              # NdupSEQ
                    r'^(\d+)_([\d]+)dup([ACGTN]+)$'                      # N_MdupSEQ
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
                # dup 归并：若可将 delins/ins 归并为 dup，则替换 c. 并记录 corrected_cdot
                try:
                    dup_try = detect_duplication_cdot(cdot_to_use, coding_dna)
                    if dup_try:
                        cdot_to_use = dup_try[0]
                        variant['corrected_cdot'] = cdot_to_use
                except Exception:
                    pass
                # intronic 插入（c.N+/-a_N+/-binsSEQ）不在此处归并，由 annotate_variation 负责；
                # 但若 EVF 给出 intronic ins，后续 p. 计算会忽略，不会生成错误的蛋白注释。
                mutated_cds, effect_hint = apply_c_hgvs_to_cds(coding_dna, cdot_to_use)
                wt_protein = translate_protein(coding_dna, chrom)
                mut_protein = translate_protein(mutated_cds, chrom)
                # SNV 强制短路：若新密码子为终止子，则直接输出 stopgain（避免被diff误判为删除）
                m_snv = re.match(r'^c\.(\d+)[ACGTN]?>[ACGTN]$', cdot_to_use, re.IGNORECASE)
                if m_snv:
                    pos_nt = int(m_snv.group(1))
                    codon_idx = (pos_nt - 1) // 3
                    wt_codon = coding_dna[codon_idx*3: codon_idx*3+3]
                    mut_codon = mutated_cds[codon_idx*3: codon_idx*3+3]
                    wt_aa = CODON_TABLE.get(wt_codon.upper(), 'X')
                    mut_aa = CODON_TABLE.get(mut_codon.upper(), 'X')
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
                            # 从差异起点向后找 '*' 终止子
                            pos1, _, _, _ = diff_proteins(wt_protein, mut_protein)
                            tail = mut_protein[pos1:]
                            stop_idx = tail.find('*')
                            fs_len = stop_idx + 1 if stop_idx >= 0 else 0
                            aa_new = mut_protein[pos1-1] if pos1-1 < len(mut_protein) and mut_protein[pos1-1] != '*' else '?'
                            p_hgvs = f"p.{three_letter(wt_protein[pos1-1])}{pos1}{three_letter(aa_new)}fs*{fs_len}"
                        except Exception:
                            pass
                    effect = 'frameshift'
                return {'type': effect.replace('_', ' ') if effect else 'unknown', 'effect': effect or 'unknown', 'p_hgvs': p_hgvs}

            # 回退：使用粗略的基于基因组坐标的方法
            mrna_pos = variant['start'] - cds_start
            ref = variant['ref']
            alt = variant['alt']
            if alt in ['-', '*'] and len(ref) >= 1:
                return self._analyze_deletion(variant, mrna_seq, mrna_pos, codon_table)
            if len(ref) == len(alt):
                if len(ref) == 1:
                    return self._analyze_substitution(variant, mrna_seq, mrna_pos, codon_table)
                else:
                    return self._analyze_block_substitution(variant, mrna_seq, mrna_pos, codon_table)
            elif len(ref) > len(alt):
                return self._analyze_deletion(variant, mrna_seq, mrna_pos, codon_table)
            else:
                return self._analyze_insertion(variant, mrna_seq, mrna_pos, codon_table)
        
        except Exception as e:
            if self.tolerate:
                logger.warning(f"分析蛋白质变化失败: {e}")
                return {'type': 'unknown', 'effect': 'unknown'}
            else:
                raise

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
        if not transcript or transcript not in self._gene_meta:
            return None
        meta = self._gene_meta[transcript]
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
        ref_g = (variant.get('ref') or '').upper()
        alt_g = (variant.get('alt') or '').upper()
        if not transcript or transcript not in self._gene_meta:
            return None
        meta = self._gene_meta[transcript]
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
        # 规范 indel 两端：以MATCHVAR风格，假设ref/alt共享左侧锚碱基，取差异部分
        # 插入
        if len(alt_g) > len(ref_g):
            ins_seq = alt_g[len(ref_g):]
            if strand == '-':
                ins_seq = reverse_complement(ins_seq)
            
            # 检查原始REF是否为"-"，如果是则优先使用delins格式
            original_ref = variant.get('ref', '')
            if original_ref == '-' or original_ref == '.' or original_ref == '*':
                # REF为"-"的插入，使用delins格式，单个位置
                return f"c.{pos_in_cds}delins{ins_seq}"
            else:
                # 常规插入格式
                if strand == '+':
                    return f"c.{pos_in_cds}_{pos_in_cds+1}ins{ins_seq}"
                else:
                    return f"c.{pos_in_cds-1}_{pos_in_cds}ins{ins_seq}"
        # 缺失
        if len(ref_g) > len(alt_g):
            del_seq = ref_g[len(alt_g):]
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
    
    def _analyze_substitution(self, variant: Dict, mrna_seq: str, mrna_pos: int, codon_table: Dict) -> Dict:
        """分析替换变异"""
        ref = variant['ref']
        alt = variant['alt']
        
        # 获取原始密码子
        codon_start = (mrna_pos // 3) * 3
        original_codon = mrna_seq[codon_start:codon_start + 3]
        
        # 计算变异在密码子中的位置
        pos_in_codon = mrna_pos % 3
        
        # 构建新密码子
        new_codon = list(original_codon)
        new_codon[pos_in_codon] = alt
        new_codon = ''.join(new_codon)
        
        # 翻译密码子
        original_aa = codon_table.get(original_codon, 'X')
        new_aa = codon_table.get(new_codon, 'X')
        
        protein_pos = (mrna_pos // 3) + 1
        # 构建p.注释
        def three(aa: str) -> str:
            return AA_ONE_TO_THREE.get(aa, 'Xaa')
        
        if original_aa == new_aa:
            p_hgvs = f"p.{three(original_aa)}{protein_pos}{three(new_aa)}" if original_aa != 'X' else "p.?"
            return {'type': 'synonymous', 'effect': 'synonymous_SNV', 'p_hgvs': p_hgvs}
        # 终止子产生/消失
        if new_aa == '*':
            p_hgvs = f"p.{three(original_aa)}{protein_pos}Ter"
            return {'type': 'stopgain', 'effect': 'stopgain', 'p_hgvs': p_hgvs}
        if original_aa == '*':
            p_hgvs = f"p.Ter{protein_pos}{three(new_aa)}"
            return {'type': 'stoploss', 'effect': 'stoploss', 'p_hgvs': p_hgvs}
        # 普通错义
        p_hgvs = f"p.{three(original_aa)}{protein_pos}{three(new_aa)}"
        return {'type': 'nonsynonymous', 'effect': 'nonsynonymous_SNV', 'p_hgvs': p_hgvs}
    
    def _analyze_block_substitution(self, variant: Dict, mrna_seq: str, mrna_pos: int, codon_table: Dict) -> Dict:
        """分析块替换变异"""
        ref = variant['ref']
        alt = variant['alt']
        
        # 检查是否为3的倍数
        if len(ref) % 3 == 0 and len(alt) % 3 == 0:
            return {'type': 'in-frame_substitution', 'effect': 'inframe_substitution'}
        else:
            return {'type': 'frameshift_substitution', 'effect': 'frameshift_substitution'}
    
    def _analyze_deletion(self, variant: Dict, mrna_seq: str, mrna_pos: int, codon_table: Dict) -> Dict:
        """分析删除变异"""
        ref = variant['ref']
        
        # 检查是否为3的倍数
        if len(ref) % 3 == 0:
            return {'type': 'in-frame_deletion', 'effect': 'nonframeshift_deletion', 'p_hgvs': 'p.?'}
        else:
            return {'type': 'frameshift_deletion', 'effect': 'frameshift_deletion', 'p_hgvs': 'p.?'}
    
    def _analyze_insertion(self, variant: Dict, mrna_seq: str, mrna_pos: int, codon_table: Dict) -> Dict:
        """分析插入变异"""
        alt = variant['alt']
        
        # 检查是否为3的倍数
        if len(alt) % 3 == 0:
            return {'type': 'in-frame_insertion', 'effect': 'nonframeshift_insertion', 'p_hgvs': 'p.?'}
        else:
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
                                m1 = re.match(r'^([\w\-\.\@\/]+?):([\w\.\-]+?):(exon\d+):(c\.[\w\->]+)(:p\.[\w\*]+)?$', item)
                                # 2) TRANSCRIPT:exonX:c.xxx[:p.yyy]
                                m2 = re.match(r'^([\w\.\-]+):(exon\d+):(c\.[\w\->]+)(:p\.[\w\*]+)?$', item)
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
    parser.add_argument('-includesnp', action='store_true', help='包含SNP')
    parser.add_argument('-mrnaseq', action='store_true', help='mRNA序列')
    parser.add_argument('-onlyAltering', action='store_true', help='仅改变序列')
    parser.add_argument('-codingseq', action='store_true', help='编码序列')
    parser.add_argument('-alltranscript', action='store_true', help='所有转录本')
    parser.add_argument('-newevf', help='新的外显子变异功能文件')
    parser.add_argument('-outfile', help='输出文件')
    parser.add_argument('-tolerate', action='store_true', help='容忍错误')
    
    # 添加新的重要参数
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