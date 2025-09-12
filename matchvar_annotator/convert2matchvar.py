#!/usr/bin/env python3
"""
MATCHVAR数据格式转换工具：将各种格式的变异文件转换为MATCHVAR输入格式
"""

import os
import sys
import argparse
import logging
import re
from typing import Dict, List, Tuple, Optional, Any

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# IUPAC密码子表
IUPAC = {
    'R': 'AG', 'Y': 'CT', 'S': 'CG', 'W': 'AT', 'K': 'GT', 'M': 'AC',
    'A': 'AA', 'C': 'CC', 'G': 'GG', 'T': 'TT', 'B': 'CGT', 'D': 'AGT',
    'H': 'ACT', 'V': 'ACG', 'N': 'ACGT', '.': '-', '-': '-'
}

IUPAC_REV = {v: k for k, v in IUPAC.items()}

class Convert2Matchvar:
    """MATCHVAR数据格式转换工具：将各种格式的变异文件转换为MATCHVAR输入格式"""
    def __init__(self, variantfile: str, **kwargs):
        self.variantfile = variantfile
        self.outfile = kwargs.get('outfile')
        self.format = kwargs.get('format', 'pileup')
        self.includeinfo = kwargs.get('includeinfo', False)
        self.snpqual = kwargs.get('snpqual', 20)
        self.snppvalue = kwargs.get('snppvalue', 1)
        self.coverage = kwargs.get('coverage', 0)
        self.maxcoverage = kwargs.get('maxcoverage')
        self.chr = kwargs.get('chr')
        self.chrmt = kwargs.get('chrmt', 'M')
        self.altcov = kwargs.get('altcov')
        self.allelicfrac = kwargs.get('allelicfrac', False)
        self.fraction = kwargs.get('fraction', 0)
        self.species = kwargs.get('species', False)
        self.filter = kwargs.get('filter')
        self.confraction = kwargs.get('confraction', 0)
        self.allallele = kwargs.get('allallele', False)
        self.withzyg = kwargs.get('withzyg', False)
        self.comment = kwargs.get('comment', False)
        self.allsample = kwargs.get('allsample', False)
        self.genoqual = kwargs.get('genoqual')
        self.varqual = kwargs.get('varqual')
        self.dbsnpfile = kwargs.get('dbsnpfile')
        self.withfreq = kwargs.get('withfreq', False)
        self.withfilter = kwargs.get('withfilter', False)
        self.seqdir = kwargs.get('seqdir')
        self.seqfile = kwargs.get('seqfile')
        self.inssize = kwargs.get('inssize')
        self.delsize = kwargs.get('delsize')
        self.subsize = kwargs.get('subsize')
        self.genefile = kwargs.get('genefile')
        self.splicing_threshold = kwargs.get('splicing_threshold')
        self.context = kwargs.get('context', False)
        self.avsnpfile = kwargs.get('avsnpfile')
        self.keepindelref = kwargs.get('keepindelref', False)
        
        # 添加新的重要参数
        self.verbose = kwargs.get('verbose', False)  # 详细输出
        self.man = kwargs.get('man', False)  # 手册
        
        # 处理参数
        self._process_arguments()
    
    def _process_arguments(self):
        """处理参数"""
        if not self.format:
            self.format = 'pileup'
            logger.info("NOTICE: the default --format argument is set as 'pileup'")
        
        if self.format == 'vcf':
            self.format = 'vcf4'
        
        if self.allsample:
            if not self.withfreq and not self.outfile:
                raise ValueError("Error in argument: please specify --outfile when --allsample is specified (unless -withfreq is set)")
            if self.format != 'vcf4':
                raise ValueError("Error in argument: the --allsample argument is supported only if --format is 'vcf4'")
            if not self.withfreq:
                logger.info(f"NOTICE: output files will be written to {self.outfile}.<samplename>.mvinput")
        
        # 验证参数
        if self.snpqual and self.format not in ['pileup', 'vcf4old']:
            raise ValueError("Error in argument: the --snpqual is supported only for the 'pileup' or 'vcf4old' format")
        
        if self.snppvalue and self.format != 'gff3-solid':
            raise ValueError("Error in argument: the --snppvalue is supported only for the 'gff3-solid' format")
        
        if not self.snpqual and self.format == 'pileup':
            self.snpqual = 20
            logger.info("NOTICE: the default --snpqual argument for pileup format is set as 20")
        
        if not self.snppvalue:
            self.snppvalue = 1
        
        if not self.coverage:
            self.coverage = 0

        # vcf4 模式下，不强制使用 snpqual 过滤；若未提供则设为0（不过滤）
        if self.format == 'vcf4' and self.snpqual is None:
            self.snpqual = 0.0
        
        if self.fraction is not None:
            if self.format not in ['pileup', 'vcf4']:
                raise ValueError("Error in argument: the '--fraction' argument is supported for the pileup or vcf4 format only")
            if self.format == 'vcf4old':
                logger.info("NOTICE: the --fraction argument works ONLY on indels for vcf4old format")
            if not (0 <= self.fraction <= 1):
                raise ValueError("Error in argument: the --fraction argument must be between 0 and 1 inclusive")
        else:
            self.fraction = 0
        
        if self.withfilter and self.format not in ['vcf4', 'vcf4old']:
            raise ValueError("Error in argument: the '-withfilter' argument is supported for the vcf4 or vcf4old format only")
        
        if self.confraction is not None:
            if self.format == 'vcf4old':
                logger.info("NOTICE: the --confraction argument works ONLY on indels for vcf4old format")
            if not (0 <= self.confraction <= 1):
                raise ValueError("Error in argument: the --confraction argument must be between 0 and 1 inclusive")
        else:
            self.confraction = 0
        
        if self.altcov is not None:
            if self.format != 'pileup':
                raise ValueError("Error in argument: the '--altcov' argument is supported for the '--format pileup' only")
            if self.altcov >= self.coverage:
                raise ValueError("Error in argument: the --altcov argument must be less than --coverage")
            if self.altcov <= 0:
                raise ValueError("Error in argument: the --altcov argument must be a positive integer")
        
        if self.species and self.format != 'gff3-solid':
            raise ValueError("Error in argument: the '--species' argument is only necessary for the '--format gff3-solid'")
    
    def convert(self):
        """执行转换"""
        try:
            with open(self.variantfile, 'r', encoding='utf-8') as f:
                if self.format == 'vcf4':
                    self._convert_vcf4(f)
                elif self.format == 'pileup':
                    self._convert_pileup(f)
                elif self.format == 'gff3-solid':
                    self._convert_gff3_solid(f)
                else:
                    raise ValueError(f"Unsupported format: {self.format}")
        except Exception as e:
            logger.error(f"Error during conversion: {e}")
            raise
    
    def _convert_vcf4(self, input_file):
        """VCF4 → MATCHVAR"""
        logger.info("Converting VCF4 format...")

        header_lines: List[str] = []
        sample_names: List[str] = []

        # 读取头
        for line in input_file:
            line = line.rstrip('\n')
            if line.startswith('##'):
                header_lines.append(line)
                continue
            if line.startswith('#'):
                parts = line.split('\t')
                if len(parts) > 9:
                    sample_names = parts[9:]
                header_lines.append(line)
                break

        # 逐行处理
        for raw in input_file:
            raw = raw.rstrip('\n')
            if not raw or raw.startswith('#'):
                continue
            cols = raw.split('\t')
            if len(cols) < 8:
                continue

            chrom, start, vid, ref, alt, qual, filt, info = cols[:8]
            fmt = cols[8] if len(cols) > 8 else ''
            samples = cols[9:] if len(cols) > 9 else ['NULL']

            # 质量与过滤
            try:
                if self.snpqual is not None and qual not in ['.', 'NA', 'NaN', 'nan', '']:
                    if float(qual) < float(self.snpqual):
                        continue
            except Exception:
                pass
            if self.withfilter and filt not in ('PASS', '.'):
                continue

            # 多等位
            alts = alt.split(',') if alt else ['.']

            # 解析 FORMAT 位置
            gtpos = dppos = gqpos = None
            if fmt:
                fmt_keys = fmt.split(':')
                for idx, key in enumerate(fmt_keys):
                    if key == 'GT':
                        gtpos = idx
                    elif key == 'DP' and dppos is None:
                        dppos = idx
                    elif key == 'GQ':
                        gqpos = idx
                    elif key == 'NR' and dppos is None:
                        dppos = idx
                    elif key == 'AD' and dppos is None:
                        dppos = idx

            # 进一步调整坐标和等位基因，参考convert2annovar.pl中的adjustStartEndRefAlt
            def adjust_start_end_ref_alt(_start: int, _end: int, _ref: str, _alt: str) -> Tuple[int, int, str, str]:
                # 从末尾开始移除相同的碱基
                while len(_ref) > 0 and len(_alt) > 0 and _ref[-1] == _alt[-1]:
                    _ref = _ref[:-1]
                    _alt = _alt[:-1]
                    _end -= 1
                    if not _ref:
                        _ref = '-'
                        _start -= 1  # 插入位置需要调整
                        break
                    if not _alt:
                        _alt = '-'
                        break
                
                # 从开头开始移除相同的碱基
                while len(_ref) > 0 and len(_alt) > 0 and _ref[0] == _alt[0]:
                    _ref = _ref[1:]
                    _alt = _alt[1:]
                    _start += 1
                    if not _ref:
                        _ref = '-'
                        _start -= 1  # 插入位置需要调整
                        break
                    if not _alt:
                        _alt = '-'
                        break
                
                # 确保删除变异的终止位置正确
                if _alt == '-' and _ref != '-':
                    # 删除变异：终止位置应该是起始位置 + 删除长度 - 1
                    _end = _start + len(_ref) - 1
                elif _ref == '-' and _alt != '-':
                    # 插入变异：起始和终止位置相同
                    _end = _start
                elif _ref != '-' and _alt != '-':
                    # 替换变异：终止位置应该是起始位置 + 参考序列长度 - 1
                    _end = _start + len(_ref) - 1
                
                return _start, _end, _ref, _alt
            
            # 归一化 indel 坐标：返回 (newstart,newend,newref,newalt)
            def left_normalize(_start: int, _ref: str, _alt: str) -> Tuple[int, int, str, str]:
                _start = int(_start)
                if len(_ref) == 1 and len(_alt) == 1:
                    return _start, _start, _ref, _alt
                
                # 插入或块替换
                if len(_ref) < len(_alt):
                    head = _alt[:len(_ref)]
                    if head == _ref:
                        new_start = _start + len(_ref) - 1
                        new_end = _start + len(_ref) - 1
                        new_ref = '-'
                        new_alt = _alt[len(_ref):]
                        # 进一步调整
                        return adjust_start_end_ref_alt(new_start, new_end, new_ref, new_alt)
                    new_start = _start
                    new_end = _start + len(_ref) - 1
                    new_ref = _ref
                    new_alt = _alt
                    # 进一步调整
                    return adjust_start_end_ref_alt(new_start, new_end, new_ref, new_alt)
                
                # 删除或块替换
                elif len(_ref) > len(_alt):
                    head = _ref[:len(_alt)]
                    if head == _alt:
                        new_start = _start + len(_alt)
                        new_end = _start + len(_ref) - 1
                        new_ref = _ref[len(_alt):]
                        new_alt = '-'
                        # 进一步调整
                        return adjust_start_end_ref_alt(new_start, new_end, new_ref, new_alt)
                    new_start = _start
                    new_end = _start + len(_ref) - 1
                    new_ref = _ref
                    new_alt = _alt
                    # 进一步调整
                    return adjust_start_end_ref_alt(new_start, new_end, new_ref, new_alt)
                
                # 相同长度的块替换
                else:
                    # 特殊处理：检查是否可以通过移除最后一个字符来简化变异
                    head = _ref[:len(_ref) - 1]
                    if _alt.startswith(head):
                        new_start = _start + len(_ref) - 1
                        new_end = _start + len(_ref) - 1
                        # 移除最后一个字符
                        new_ref = _ref[-1] if _ref else '-'
                        new_alt = _alt[-1] if _alt else '-'
                        # 进一步调整
                        return adjust_start_end_ref_alt(new_start, new_end, new_ref, new_alt)
                    else:
                        new_start = _start
                        new_end = _start + len(_ref) - 1
                        new_ref = _ref
                        new_alt = _alt
                        # 进一步调整
                        return adjust_start_end_ref_alt(new_start, new_end, new_ref, new_alt)

            # 统计 withfreq 需要的计数
            def compute_sample_metrics(sample_str: str) -> Tuple[Optional[int], Optional[int], Optional[int]]:
                rd = gq = None
                dp_val = None
                if sample_str and sample_str != 'NULL' and fmt:
                    arr = sample_str.split(':')
                    if dppos is not None and dppos < len(arr):
                        try:
                            if ',' in arr[dppos]:
                                dp_val = sum(int(x) for x in arr[dppos].split(',') if x.isdigit())
                            else:
                                dp_val = int(arr[dppos]) if arr[dppos].isdigit() else None
                        except Exception:
                            dp_val = None
                    if gqpos is not None and gqpos < len(arr):
                        try:
                            gq = int(arr[gqpos]) if arr[gqpos].isdigit() else None
                        except Exception:
                            gq = None
                if dp_val is None:
                    m = re.search(r'\bDP=(\d+)\b', info)
                    if m:
                        dp_val = int(m.group(1))
                if dp_val is None:
                    m = re.search(r'\bDP4=(\d+),(\d+),(\d+),(\d+)\b', info)
                    if m:
                        dp_val = sum(int(x) for x in m.groups())
                rd = dp_val
                return rd, gq, dp_val

            # allele 循环
            for i_alt, alt_i in enumerate(alts):
                # 统计频率
                count_alt = 0
                count_all = 0
                # 逐样本统计
                for s in samples:
                    if s == 'NULL' or not fmt or gtpos is None:
                        continue
                    arr = s.split(':')
                    if gtpos >= len(arr):
                        continue
                    gt = arr[gtpos]
                    a1 = a2 = None
                    if re.match(r'^[\d\.]+[/|][\d\.]+$', gt):
                        a1, a2 = re.split(r'[/|]', gt)
                    elif re.match(r'^[\d\.]+$', gt):
                        a1 = a2 = gt
                    else:
                        continue
                    
                    # 处理 . 和 0 的情况
                    if a1 != '.':
                        count_all += 1
                    if a2 != '.':
                        count_all += 1
                    
                    # 将 . 转换为 0
                    if a1 == '.':
                        a1 = '0'
                    if a2 == '.':
                        a2 = '0'
                    
                    # 统计等位基因
                    if a1 == '0' and a2 == '0':  # ref/ref call or unknown call
                        # 不增加 count_alt
                        pass
                    elif a1 == a2:  # homozygous
                        if a1.isdigit() and int(a1) == i_alt + 1:
                            count_alt += 2
                    else:  # heterozygous
                        if a1.isdigit() and int(a1) == i_alt + 1:
                            count_alt += 1
                        if a2.isdigit() and int(a2) == i_alt + 1:
                            count_alt += 1

                # 处理特殊变异类型
                if alt_i.upper() == '<DEL>':  # 删除变异
                    # 如果INFO字段中有END信息，使用它
                    end_match = re.search(r'\bEND=(\d+)\b', info)
                    if end_match:
                        newstart, newend = int(start), int(end_match.group(1))
                        newref, newalt = '0', '-'  # 与 convert2annovar.pl 保持一致
                    else:
                        # 删除变异：起始位置是start，终止位置是start + len(ref) - 1
                        newstart, newend = int(start), int(start) + len(ref) - 1
                        newref, newalt = ref, '-'
                elif alt_i.upper() in ('<DUP>', '<INV>', '<INS>'):  # 复制、倒位、插入变异
                    # 如果INFO字段中有END信息，使用它
                    end_match = re.search(r'\bEND=(\d+)\b', info)
                    if end_match:
                        newstart, newend = int(start), int(end_match.group(1))
                        newref, newalt = '0', '0'  # 与 convert2annovar.pl 保持一致
                    else:
                        newstart, newend = int(start), int(start) + len(ref) - 1
                        newref, newalt = ref, '-'
                else:
                    # 归一化坐标
                    if self.keepindelref:
                        # 如果 keepindelref 为真，不进行坐标归一化
                        newstart, newend = int(start), int(start) + len(ref) - 1
                        newref, newalt = ref, alt_i
                    else:
                        newstart, newend, newref, newalt = left_normalize(start, ref, alt_i)

                if self.includeinfo:
                    if self.withzyg:
                        # 合子性/质量/深度三列
                        # 简化：仅按首样本给出 het/hom
                        zyg = '.'
                        if samples and fmt and gtpos is not None:
                            gt0 = samples[0].split(':')[gtpos] if gtpos < len(samples[0].split(':')) else ''
                            if re.match(r'^(\d+)[/|](\d+)$', gt0):
                                g1, g2 = re.match(r'^(\d+)[/|](\d+)$', gt0).groups()
                                zyg = 'hom' if g1 == g2 and g1 == str(i_alt + 1) else ('het' if g1 == str(i_alt + 1) or g2 == str(i_alt + 1) else '.')
                        rd0, gq0, _ = compute_sample_metrics(samples[0])
                        out = [chrom, str(newstart), str(newend), newref, newalt, zyg, qual, str(rd0) if rd0 is not None else '.']
                        print('\t'.join(out + [raw]))
                    elif self.withfreq:
                        # 频率、质量、深度 + 原始行
                        freq = (count_alt / count_all) if count_all else '.'
                        # 近似深度用 INFO.DP 或首样本DP
                        rd0, _, _ = compute_sample_metrics(samples[0])
                        out = [chrom, str(newstart), str(newend), newref, newalt, f"{freq if isinstance(freq,str) else (f'{freq:.4g}')}", qual, str(rd0) if rd0 is not None else '.']
                        print('\t'.join(out + [raw]))
                    else:
                        # 附带 VCF 关键列
                        extra = [chrom, start, vid, ref, alt_i, qual, filt, info, fmt, samples[0] if samples else 'NULL']
                        print('\t'.join([chrom, str(newstart), str(newend), newref, newalt] + extra))
                else:
                    if self.withfreq:
                        freq = (count_alt / count_all) if count_all else '.'
                        rd0, _, _ = compute_sample_metrics(samples[0])
                        out = [chrom, str(newstart), str(newend), newref, newalt, f"{freq if isinstance(freq,str) else (f'{freq:.4g}')}", qual, str(rd0) if rd0 is not None else '.']
                        print('\t'.join(out))
                    else:
                        # 默认输出合子性/质量/深度
                        zyg = '.'
                        if samples and fmt and gtpos is not None:
                            gt0 = samples[0].split(':')[gtpos] if gtpos < len(samples[0].split(':')) else ''
                            if re.match(r'^(\d+)[/|](\d+)$', gt0):
                                g1, g2 = re.match(r'^(\d+)[/|](\d+)$', gt0).groups()
                                zyg = 'hom' if g1 == g2 and g1 == str(i_alt + 1) else ('het' if g1 == str(i_alt + 1) or g2 == str(i_alt + 1) else '.')
                        rd0, _, _ = compute_sample_metrics(samples[0])
                        out = [chrom, str(newstart), str(newend), newref, newalt]
                        if self.withfilter:
                            out += [zyg, filt, qual, str(rd0) if rd0 is not None else '.']
                        else:
                            out += [zyg, qual, str(rd0) if rd0 is not None else '.']
                        print('\t'.join(out))
    
    def _convert_pileup(self, input_file):
        """转换pileup格式"""
        logger.info("Converting pileup format...")
        
        for line in input_file:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            parts = line.split('\t')
            if len(parts) < 6:
                continue
            
            chrom = parts[0]
            pos = parts[1]
            ref = parts[2]
            coverage = int(parts[3])
            bases = parts[4]
            qualities = parts[5]
            
            # 跳过不符合覆盖度要求的变异
            if coverage < self.coverage:
                continue
            
            if self.maxcoverage and coverage > self.maxcoverage:
                continue
            
            # 解析碱基和频率
            variants = self._parse_pileup_bases(bases, ref)
            
            for variant in variants:
                alt, freq = variant
                
                # 检查频率要求
                if freq < self.fraction:
                    continue
                
                # 检查替代等位基因覆盖度
                if self.altcov and freq * coverage < self.altcov:
                    continue
                
                # 构建MATCHVAR输出行
                output_line = f"{chrom}\t{pos}\t{pos}\t{ref}\t{alt}"
                print(output_line)
    
    def _parse_pileup_bases(self, bases: str, ref: str) -> List[Tuple[str, float]]:
        """解析pileup格式的碱基字符串"""
        variants = []
        base_counts = {}
        
        # 统计碱基频率
        for base in bases.upper():
            if base in 'ACGT':
                base_counts[base] = base_counts.get(base, 0) + 1
        
        total = sum(base_counts.values())
        if total == 0:
            return variants
        
        # 计算频率并过滤
        for base, count in base_counts.items():
            if base != ref:
                freq = count / total
                variants.append((base, freq))
        
        return variants
    
    def _convert_gff3_solid(self, input_file):
        """转换GFF3-SOLID格式"""
        logger.info("Converting GFF3-SOLID format...")
        
        for line in input_file:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            parts = line.split('\t')
            if len(parts) < 9:
                continue
            
            chrom = parts[0]
            source = parts[1]
            feature = parts[2]
            start = parts[3]
            end = parts[4]
            score = parts[5]
            strand = parts[6]
            frame = parts[7]
            attributes = parts[8]
            
            # 解析属性
            attr_dict = self._parse_gff_attributes(attributes)
            
            # 检查P值
            if 'pvalue' in attr_dict:
                pvalue = float(attr_dict['pvalue'])
                if pvalue > self.snppvalue:
                    continue
            
            # 提取变异信息
            if 'ref' in attr_dict and 'alt' in attr_dict:
                ref = attr_dict['ref']
                alt = attr_dict['alt']
                
                # 构建MATCHVAR输出行
                output_line = f"{chrom}\t{start}\t{end}\t{ref}\t{alt}"
                print(output_line)
    
    def _parse_gff_attributes(self, attributes: str) -> Dict[str, str]:
        """解析GFF属性字符串"""
        attr_dict = {}
        
        for attr in attributes.split(';'):
            if '=' in attr:
                key, value = attr.split('=', 1)
                attr_dict[key.strip()] = value.strip()
        
        return attr_dict

def main():
    """主函数"""
    examples = (
        "示例:\n"
        "1) 从VCF转换为MATCHVAR（附带原始信息列）:\n"
        "   python utils/matchvar/convert2matchvar.py \\\n+        /Users/James/PycharmProjects/Variant_Data_Simulation_2.0/resources/202511.family.vcf \\\n+        -format vcf4 -includeinfo > out.mvinput\n\n"
        "2) 从VCF转换并输出频率列（withfreq）:\n"
        "   python utils/matchvar/convert2matchvar.py input.vcf -format vcf4 -withfreq > out.mvinput\n\n"
        "3) Pileup格式按覆盖度与频率阈值过滤:\n"
        "   python utils/matchvar/convert2matchvar.py input.pileup -format pileup -coverage 10 -fraction 0.05 > out.mvinput\n\n"
        "4) GFF3-SOLID示例（按P值阈值）:\n"
        "   python utils/matchvar/convert2matchvar.py input.gff3 -format gff3-solid -snppvalue 0.01 > out.mvinput\n\n"
        "5) VCF多样本，按样本输出（需 -outfile 前缀）:\n"
        "   python utils/matchvar/convert2matchvar.py cohort.vcf -format vcf4 -allsample -withfreq -outfile result\n\n"
        "提示: 可与 table_matchvar.py 配合使用：\n"
        "   python utils/matchvar/table_matchvar.py out.mvinput resources/humandb -outfile anno -buildver hg19 -protocol refGene -operation g\n"
    )
    parser = argparse.ArgumentParser(
        description='MATCHVAR数据格式转换工具',
        epilog=examples,
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('variantfile', help='输入变异文件')
    parser.add_argument('-outfile', help='输出文件')
    parser.add_argument('-format', default='pileup', help='输入格式')
    parser.add_argument('-includeinfo', action='store_true', help='包含信息字段')
    parser.add_argument('-snpqual', type=float, help='SNP质量阈值')
    parser.add_argument('-snppvalue', type=float, help='SNP P值阈值')
    parser.add_argument('-coverage', type=int, help='覆盖度阈值')
    parser.add_argument('-maxcoverage', type=int, help='最大覆盖度')
    parser.add_argument('-chr', help='染色体名称')
    parser.add_argument('-chrmt', default='M', help='线粒体染色体名称')
    parser.add_argument('-altcov', type=int, help='替代等位基因覆盖度')
    parser.add_argument('-allelicfrac', action='store_true', help='等位基因频率')
    parser.add_argument('-fraction', type=float, help='频率阈值')
    parser.add_argument('-species', action='store_true', help='物种信息')
    parser.add_argument('-filter', help='过滤条件')
    parser.add_argument('-confraction', type=float, help='置信频率')
    parser.add_argument('-allallele', action='store_true', help='所有等位基因')
    parser.add_argument('-withzyg', action='store_true', help='包含合子性')
    parser.add_argument('-comment', action='store_true', help='包含注释')
    parser.add_argument('-allsample', action='store_true', help='所有样本')
    parser.add_argument('-genoqual', type=float, help='基因型质量')
    parser.add_argument('-varqual', type=float, help='变异质量')
    parser.add_argument('-dbsnpfile', help='dbSNP文件')
    parser.add_argument('-withfreq', action='store_true', help='包含频率')
    parser.add_argument('-withfilter', action='store_true', help='包含过滤信息')
    parser.add_argument('-seqdir', help='序列目录')
    parser.add_argument('-seqfile', help='序列文件')
    parser.add_argument('-inssize', type=int, help='插入大小')
    parser.add_argument('-delsize', type=int, help='删除大小')
    parser.add_argument('-subsize', type=int, help='替换大小')
    parser.add_argument('-genefile', help='基因文件')
    parser.add_argument('-splicing_threshold', type=int, help='剪接阈值')
    parser.add_argument('-context', action='store_true', help='上下文')
    parser.add_argument('-avsnpfile', help='AV SNP文件')
    parser.add_argument('-keepindelref', action='store_true', help='保留INDEL参考')
    
    # 添加新的重要参数
    parser.add_argument('-verbose', '-v', action='store_true', help='详细输出')
    parser.add_argument('-man', '-m', action='store_true', help='显示手册')
    
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args()
    
    # 创建转换器
    converter = Convert2Matchvar(
        variantfile=args.variantfile,
        outfile=args.outfile,
        format=args.format,
        includeinfo=args.includeinfo,
        snpqual=args.snpqual,
        snppvalue=args.snppvalue,
        coverage=args.coverage,
        maxcoverage=args.maxcoverage,
        chr=args.chr,
        chrmt=args.chrmt,
        altcov=args.altcov,
        allelicfrac=args.allelicfrac,
        fraction=args.fraction,
        species=args.species,
        filter=args.filter,
        confraction=args.confraction,
        allallele=args.allallele,
        withzyg=args.withzyg,
        comment=args.comment,
        allsample=args.allsample,
        genoqual=args.genoqual,
        varqual=args.varqual,
        dbsnpfile=args.dbsnpfile,
        withfreq=args.withfreq,
        withfilter=args.withfilter,
        seqdir=args.seqdir,
        seqfile=args.seqfile,
        inssize=args.inssize,
        delsize=args.delsize,
        subsize=args.subsize,
        genefile=args.genefile,
        splicing_threshold=args.splicing_threshold,
        context=args.context,
        avsnpfile=args.avsnpfile,
        keepindelref=args.keepindelref,
        verbose=args.verbose,
        man=args.man
    )
    
    # 执行转换
    try:
        converter.convert()
        logger.info("转换成功完成")
    except Exception as e:
        logger.error(f"转换过程中发生错误: {e}")
        sys.exit(1)

def parse_vcf_file(vcf_content: str) -> List[Dict[str, Any]]:
    """
    解析VCF文件内容并返回变异列表
    
    Args:
        vcf_content: VCF文件内容字符串
        
    Returns:
        变异列表，每个变异包含chromosome, position, reference, alternate字段
    """
    variants = []
    lines = vcf_content.strip().split('\n')
    
    logger.info(f"解析VCF内容，总行数: {len(lines)}")
    
    for i, line in enumerate(lines):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
            
        parts = line.split('\t')
        if len(parts) < 5:
            logger.warning(f"VCF第 {i+1} 行字段数不足，跳过: {line}")
            continue
            
        try:
            chromosome = parts[0]
            position = int(parts[1])
            reference = parts[3]
            # 处理多个替代等位基因，只取第一个
            alternate = parts[4].split(',')[0] if parts[4] != '.' else parts[3]
            
            # 验证碱基序列
            if not _validate_bases(reference, alternate):
                logger.warning(f"VCF第 {i+1} 行无效的碱基序列，跳过: ref={reference}, alt={alternate}")
                continue
                
            variant = {
                'chromosome': chromosome,
                'position': position,
                'reference': reference,
                'alternate': alternate
            }
            
            variants.append(variant)
            
        except (ValueError, IndexError) as e:
            logger.warning(f"VCF第 {i+1} 行解析错误，跳过: {e}, 行内容: {line}")
            continue
    
    logger.info(f"VCF总共解析到 {len(variants)} 个变异")
    return variants


def parse_bed_file(bed_content: str) -> List[Dict[str, Any]]:
    """
    解析BED文件内容并返回变异列表
    
    Args:
        bed_content: BED文件内容字符串
        
    Returns:
        变异列表，每个变异包含chromosome, position, reference, alternate字段
    """
    variants = []
    lines = bed_content.strip().split('\n')
    
    logger.info(f"解析BED内容，总行数: {len(lines)}")
    
    for i, line in enumerate(lines):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
            
        parts = line.split('\t')
        if len(parts) < 3:
            logger.warning(f"BED第 {i+1} 行字段数不足，跳过: {line}")
            continue
            
        try:
            chromosome = parts[0]
            start = int(parts[1])  # BED是0-based
            end = int(parts[2])
            
            # BED文件格式有多种变体
            if len(parts) >= 6:
                # 如果有第6列，可能包含序列信息
                reference = parts[3] if len(parts) > 3 and parts[3] != '.' else 'N'
                alternate = parts[4] if len(parts) > 4 and parts[4] != '.' else 'A'
            elif len(parts) >= 4:
                # 如果有第4列，可能是变异名称或序列信息
                name = parts[3]
                if '>' in name:  # 格式如 "A>T"
                    ref_alt = name.split('>')
                    reference = ref_alt[0] if len(ref_alt) > 0 else 'N'
                    alternate = ref_alt[1] if len(ref_alt) > 1 else 'A'
                else:
                    reference = 'N'  # 默认参考碱基
                    alternate = 'A'  # 默认替代碱基
            else:
                # 标准3列BED文件，使用默认碱基
                reference = 'N'
                alternate = 'A'
            
            # 验证碱基序列
            if not _validate_bases(reference, alternate):
                logger.warning(f"BED第 {i+1} 行无效的碱基序列，使用默认值: ref={reference}, alt={alternate}")
                reference = 'N'
                alternate = 'A'
            
            # BED是0-based，转换为1-based位置
            position = start + 1
            
            variant = {
                'chromosome': chromosome,
                'position': position,
                'reference': reference,
                'alternate': alternate
            }
            
            variants.append(variant)
            
        except (ValueError, IndexError) as e:
            logger.warning(f"BED第 {i+1} 行解析错误，跳过: {e}, 行内容: {line}")
            continue
    
    logger.info(f"BED总共解析到 {len(variants)} 个变异")
    return variants


def _validate_bases(reference: str, alternate: str) -> bool:
    """
    验证碱基序列是否有效
    
    Args:
        reference: 参考碱基序列
        alternate: 替代碱基序列
        
    Returns:
        是否有效
    """
    if not reference or not alternate:
        return False
        
    # 允许的碱基字符（包括IUPAC代码）
    valid_bases = set('ACGTRYSWKMBDHVN.-')
    
    # 检查参考碱基
    if not all(base.upper() in valid_bases for base in reference):
        return False
        
    # 检查替代碱基
    if not all(base.upper() in valid_bases for base in alternate):
        return False
        
    return True


if __name__ == "__main__":
    main() 