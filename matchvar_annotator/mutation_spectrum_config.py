#!/usr/bin/env python3
"""
突变谱模型配置文件
包含MuRaL模型参数和gnomAD约束分数配置
"""

import os
from typing import Dict, List, Optional
from dataclasses import dataclass

@dataclass
class MuRaLConfig:
    """MuRaL模型配置"""
    # 模型文件路径
    model_file: Optional[str] = None
    
    # 默认突变率参数（基于Alexandrov et al. 2013, Nature）
    base_mutation_rates: Dict[str, float] = None
    
    # CpG位点增强因子
    cpg_enhancement_factor: float = 10.0
    
    # TpC上下文增强因子
    tpc_enhancement_factor: float = 2.0
    
    # 最小突变概率阈值
    min_mutation_probability: float = 1e-6
    
    def __post_init__(self):
        if self.base_mutation_rates is None:
            self.base_mutation_rates = {
                'C>T': 0.3,   # CpG位点高突变率
                'T>C': 0.1,
                'C>A': 0.1,
                'T>A': 0.05,
                'C>G': 0.05,
                'T>G': 0.05,
                'A>T': 0.1,
                'A>C': 0.1,
                'A>G': 0.1,
                'G>T': 0.1,
                'G>C': 0.1,
                'G>A': 0.1
            }

@dataclass
class GnomADConfig:
    """gnomAD约束分数配置"""
    # 约束分数文件路径
    constraint_file: Optional[str] = None
    
    # 约束阈值 (LOEUF < threshold 表示约束)
    constraint_threshold: float = 0.35
    
    # 高度约束阈值
    high_constraint_threshold: float = 0.1
    
    # 默认约束分数（基于gnomAD v2.1.1）
    default_constraint_scores: Dict[str, float] = None
    
    def __post_init__(self):
        if self.default_constraint_scores is None:
            self.default_constraint_scores = {
                # 极高约束基因 (LOEUF < 0.1) - 基于真实gnomAD v2.1.1数据
                'NIPBL': 0.032,  # 极高约束
                'RB1': 0.10,     # 极高约束 (默认值)
                
                # 高度约束基因 (LOEUF < 0.35)
                'APC': 0.161,    # 高度约束
                'MSH2': 0.334,   # 高度约束
                'STK11': 0.245,  # 高度约束
                
                # 中等约束基因 (LOEUF 0.35-0.65)
                'TP53': 0.469,   # 中等约束
                'BRCA2': 0.635,  # 中等约束
                'MLH1': 0.575,   # 中等约束
                'MSH6': 0.498,   # 中等约束
                'PTEN': 0.507,   # 中等约束
                'CDH1': 0.430,   # 中等约束
                'FANCM': 0.593,  # 中等约束
                
                # 低约束基因 (LOEUF > 0.65)
                'BRCA1': 0.915,  # 低约束
                'VHL': 0.927,    # 低约束
                'PMS2': 1.266,   # 低约束
                'ATM': 0.710,    # 低约束
                'CHEK2': 1.530,  # 低约束
                'PALB2': 1.006,  # 低约束
                'BARD1': 1.358,   # 低约束
                'BRIP1': 0.786,  # 低约束
                'RAD51C': 1.491,  # 低约束
                'RAD51D': 1.218, # 低约束
                'NBN': 1.010,    # 低约束
                'MRE11A': 0.839, # 低约束
                'RAD50': 0.873,  # 低约束
                'FANCC': 1.043,  # 低约束
                'FANCG': 1.091,  # 低约束
                'FANCA': 1.366,  # 低约束
            }

@dataclass
class MutationSpectrumConfig:
    """突变谱生成器总配置"""
    # MuRaL模型配置
    mural_config: MuRaLConfig
    
    # gnomAD配置
    gnomad_config: GnomADConfig
    
    # 变异生成参数
    max_variants: int = 1000
    min_probability: float = 0.001
    use_constraint_filter: bool = True
    
    # 混合生成参数
    spectrum_ratio: float = 0.7
    
    # 输出配置
    output_format: str = 'vcf'  # 'vcf', 'bed', 'json'
    include_probability: bool = True
    include_context: bool = True
    include_constraint: bool = True
    
    def __init__(self, 
                 mural_config: Optional[MuRaLConfig] = None,
                 gnomad_config: Optional[GnomADConfig] = None,
                 **kwargs):
        self.mural_config = mural_config or MuRaLConfig()
        self.gnomad_config = gnomad_config or GnomADConfig()
        
        # 设置其他参数
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)

# 预定义配置
DEFAULT_CONFIG = MutationSpectrumConfig()

# 高精度配置（更多变异，更严格过滤）
HIGH_PRECISION_CONFIG = MutationSpectrumConfig(
    max_variants=5000,
    min_probability=0.0001,
    use_constraint_filter=True,
    spectrum_ratio=0.9,
    include_probability=True,
    include_context=True,
    include_constraint=True
)

# 快速配置（较少变异，宽松过滤）
FAST_CONFIG = MutationSpectrumConfig(
    max_variants=500,
    min_probability=0.01,
    use_constraint_filter=False,
    spectrum_ratio=0.5,
    include_probability=False,
    include_context=False,
    include_constraint=False
)

# 研究配置（平衡设置）
RESEARCH_CONFIG = MutationSpectrumConfig(
    max_variants=2000,
    min_probability=0.001,
    use_constraint_filter=True,
    spectrum_ratio=0.8,
    include_probability=True,
    include_context=True,
    include_constraint=True
)

def get_config(config_name: str = 'default') -> MutationSpectrumConfig:
    """
    获取预定义配置
    
    Args:
        config_name: 配置名称 ('default', 'high_precision', 'fast', 'research')
        
    Returns:
        配置对象
    """
    configs = {
        'default': DEFAULT_CONFIG,
        'high_precision': HIGH_PRECISION_CONFIG,
        'fast': FAST_CONFIG,
        'research': RESEARCH_CONFIG
    }
    
    if config_name not in configs:
        raise ValueError(f"未知的配置名称: {config_name}. 可用配置: {list(configs.keys())}")
    
    return configs[config_name]

def create_custom_config(**kwargs) -> MutationSpectrumConfig:
    """
    创建自定义配置
    
    Args:
        **kwargs: 配置参数
        
    Returns:
        自定义配置对象
    """
    return MutationSpectrumConfig(**kwargs)

# 配置验证函数
def validate_config(config: MutationSpectrumConfig) -> List[str]:
    """
    验证配置的有效性
    
    Args:
        config: 配置对象
        
    Returns:
        错误信息列表
    """
    errors = []
    
    # 验证变异数量
    if config.max_variants <= 0:
        errors.append("max_variants 必须大于 0")
    
    # 验证概率阈值
    if not 0 <= config.min_probability <= 1:
        errors.append("min_probability 必须在 0-1 之间")
    
    # 验证比例
    if not 0 <= config.spectrum_ratio <= 1:
        errors.append("spectrum_ratio 必须在 0-1 之间")
    
    # 验证约束阈值
    if config.gnomad_config.constraint_threshold <= 0:
        errors.append("constraint_threshold 必须大于 0")
    
    # 验证文件路径
    if config.mural_config.model_file and not os.path.exists(config.mural_config.model_file):
        errors.append(f"MuRaL模型文件不存在: {config.mural_config.model_file}")
    
    if config.gnomad_config.constraint_file and not os.path.exists(config.gnomad_config.constraint_file):
        errors.append(f"gnomAD约束文件不存在: {config.gnomad_config.constraint_file}")
    
    return errors

# 使用示例
if __name__ == "__main__":
    print("突变谱模型配置文件")
    print("可用配置:")
    
    configs = ['default', 'high_precision', 'fast', 'research']
    for config_name in configs:
        config = get_config(config_name)
        print(f"- {config_name}: max_variants={config.max_variants}, "
              f"min_probability={config.min_probability}, "
              f"spectrum_ratio={config.spectrum_ratio}")
    
    # 验证配置
    print("\n配置验证:")
    for config_name in configs:
        config = get_config(config_name)
        errors = validate_config(config)
        if errors:
            print(f"{config_name}: 错误 - {errors}")
        else:
            print(f"{config_name}: 有效")
