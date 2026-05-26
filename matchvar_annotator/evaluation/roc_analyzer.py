#!/usr/bin/env python3
"""
ROC曲线分析器
实现模拟变异与ClinVar数据库比较的ROC曲线分析
"""

import os
import logging
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Optional, Any, Tuple
from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score
from sklearn.metrics import confusion_matrix, classification_report
import json
from datetime import datetime

# 设置中文字体支持
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

logger = logging.getLogger(__name__)

class ROCAnalyzer:
    """
    ROC曲线分析器
    实现模拟变异与ClinVar数据库比较的ROC曲线分析
    """
    
    def __init__(self, output_dir: str = "datasimulator/results"):
        """
        初始化ROC分析器
        
        Args:
            output_dir: 输出目录
        """
        self.output_dir = output_dir
        os.makedirs(f"{output_dir}/roc_analysis", exist_ok=True)
        os.makedirs(f"{output_dir}/plots", exist_ok=True)
        
        # 设置绘图样式
        plt.style.use('seaborn-v0_8')
        sns.set_palette("husl")
        
        logger.info(f"ROC分析器初始化完成，输出目录: {output_dir}")
    
    def perform_roc_analysis(self, 
                            comparison_result: Dict[str, Any],
                            gene_name: str,
                            save_results: bool = True,
                            use_corrected_scores: bool = False) -> Dict[str, Any]:
        """
        执行ROC分析
        
        Args:
            comparison_result: 比较结果
            gene_name: 基因名称
            save_results: 是否保存结果
            use_corrected_scores: 是否使用矫正后得分
            
        Returns:
            ROC分析结果
        """
        logger.info(f"开始对基因 {gene_name} 进行ROC分析 (使用矫正得分: {use_corrected_scores})")
        
        if not comparison_result.get('matches'):
            logger.warning("没有匹配数据，无法进行ROC分析")
            return {'error': 'No matches available for ROC analysis'}
        
        # 准备数据
        roc_data = self._prepare_roc_data(comparison_result, use_corrected_scores)
        
        if not roc_data:
            return {'error': 'No valid data for ROC analysis'}
        
        # 计算ROC指标
        roc_metrics = self._calculate_roc_metrics(roc_data)
        
        # 计算PR指标
        pr_metrics = self._calculate_pr_metrics(roc_data)
        
        # 计算混淆矩阵
        confusion_metrics = self._calculate_confusion_metrics(roc_data)
        
        # 合并结果
        analysis_result = {
            'gene_name': gene_name,
            'timestamp': datetime.now().isoformat(),
            'total_samples': len(roc_data),
            'roc_metrics': roc_metrics,
            'pr_metrics': pr_metrics,
            'confusion_metrics': confusion_metrics,
            'data_summary': self._summarize_data(roc_data),
            'roc_data': roc_data  # 添加ROC数据到返回结果中
        }
        
        # 保存结果
        if save_results:
            self._save_roc_results(analysis_result, gene_name)
        
        logger.info("ROC分析完成")
        return analysis_result
    
    def _prepare_roc_data(self, comparison_result: Dict[str, Any], use_corrected_scores: bool = False) -> List[Dict]:
        """
        准备ROC分析数据
        基于ClinVar临床意义和Total_Score进行ROC分析
        
        Args:
            comparison_result: 比较结果
            use_corrected_scores: 是否使用矫正后得分
        """
        matches = comparison_result.get('matches', [])
        if not matches:
            return []
        
        roc_data = []
        
        for match in matches:
            # 提取特征
            features = {
                'position': match.get('simulated_position', 0),
                'ref_base': match.get('simulated_ref', ''),
                'alt_base': match.get('simulated_alt', ''),
                'variant_type': match.get('simulated_type', 'SNV'),
                'hgvs': match.get('simulated_hgvs', ''),
                'clinical_significance': match.get('clinical_significance', 'Unknown'),
                'significance': match.get('significance', 'Unknown'),
                'match_type': match.get('match_type', 'Unknown'),
                'match_score': match.get('match_score', 0.0),
                'hgvs_consistent': match.get('hgvs_consistent', False),
                'base_match': match.get('base_match', False),
                'partial_match': match.get('partial_match', False),
                'gene_name': match.get('gene_name', 'Unknown'),  # 添加基因信息
                'total_score': match.get('total_score', 0.0),  # 添加Total_Score
                'corrected_score': match.get('corrected_score', None)  # 添加矫正后得分
            }
            
            # 使用ClinVar临床意义作为真实标签（三分类）
            # 0 = Benign, 1 = Uncertain significance, 2 = Pathogenic
            clinical_significance = match.get('clinical_significance', 'Unknown')
            clinical_class = self._get_clinical_class(clinical_significance)
            is_pathogenic = self._is_pathogenic(clinical_significance)

            # 计算预测分数（优先使用矫正后得分）
            if use_corrected_scores and features.get('corrected_score') is not None:
                prediction_score = float(features['corrected_score'])
                logger.debug(f"使用矫正后得分: {prediction_score}")
            else:
                prediction_score = self._calculate_prediction_score(features, match)
                logger.debug(f"使用原始得分: {prediction_score}")

            roc_data.append({
                'features': features,
                'true_label': clinical_class,  # 三分类标签
                'prediction_score': prediction_score,
                'clinical_significance': clinical_significance,
                'match_type': match.get('match_type', 'Unknown'),
                'hgvs_consistent': match.get('hgvs_consistent', False),
                'hgvs_simulated': match.get('simulated_hgvs', ''),
                'hgvs_clinvar': match.get('clinvar_hgvs', ''),
                'base_match': match.get('base_match', False),
                'partial_match': match.get('partial_match', False),
                'gene_name': match.get('gene_name', 'Unknown'),  # 添加基因信息
                'total_score': match.get('total_score', 0.0),  # 添加Total_Score
                'is_pathogenic': is_pathogenic,  # 添加致病性标签
                'clinical_class': clinical_class  # 添加三分类标签
            })
        
        return roc_data
    
    def _is_pathogenic(self, clinical_significance: str) -> bool:
        """
        判断临床意义是否为致病性
        """
        pathogenic_terms = [
            'Pathogenic', 'Likely pathogenic', 'Pathogenic/Likely pathogenic'
        ]
        
        return clinical_significance in pathogenic_terms
    
    def _get_clinical_class(self, clinical_significance: str) -> int:
        """
        获取临床意义的三分类标签
        
        Args:
            clinical_significance: 临床意义字符串
            
        Returns:
            分类标签: 0=Benign, 1=Uncertain significance, 2=Pathogenic
        """
        if not clinical_significance:
            return 1  # 默认为Uncertain significance
            
        # Pathogenic (标签2)
        pathogenic_terms = [
            'Pathogenic', 'Likely pathogenic', 'Pathogenic/Likely pathogenic'
        ]
        
        if clinical_significance in pathogenic_terms:
            return 2
                
        # Benign (标签0)
        benign_terms = [
            'Benign', 'Likely benign', 'Benign/Likely benign'
        ]
        
        if clinical_significance in benign_terms:
            return 0
                
        # 默认为Uncertain significance (标签1)
        return 1
    
    def _calculate_prediction_score(self, features: Dict, match: Dict) -> float:
        """
        计算预测分数
        使用Total_Score作为预测分数，归一化到0-1范围
        """
        # 使用Total_Score作为预测分数
        total_score = match.get('total_score', 0.0)
        
        # 如果Total_Score不存在，尝试从features中获取
        if total_score == 0.0:
            total_score = features.get('total_score', 0.0)
        
        # 归一化Total_Score到0-1范围（假设Total_Score范围是0-10）
        # 如果Total_Score已经是0-1范围，直接使用
        if total_score <= 1.0:
            normalized_score = total_score
        else:
            # 假设Total_Score范围是0-10，归一化到0-1
            normalized_score = total_score / 10.0
        
        # 确保分数在0-1范围内
        return max(0.0, min(1.0, normalized_score))
    
    def _calculate_roc_metrics(self, roc_data: List[Dict]) -> Dict[str, Any]:
        """
        计算ROC指标（支持三分类数据，自动转换为二分类）
        """
        if not roc_data:
            return {'error': 'No data for ROC calculation'}
        
        # 提取真实标签和预测分数
        y_true = np.array([item['true_label'] for item in roc_data])
        y_scores = np.array([item['prediction_score'] for item in roc_data])
        
        # 检查数据是否有效
        unique_labels = np.unique(y_true)
        if len(unique_labels) < 2:
            # 如果只有一个类别，返回默认值
            return {
                'fpr': [0.0, 1.0],
                'tpr': [0.0, 1.0],
                'thresholds': [float('inf'), 0.0],
                'auc': 0.5,  # 随机分类器的AUC
                'optimal_thresholds': {
                    'youden': {'threshold': 0.5, 'fpr': 1.0, 'tpr': 1.0, 'youden_index': 0.0},
                    'f1': {'threshold': 0.5, 'fpr': 1.0, 'tpr': 1.0, 'f1_score': 0.0}
                },
                'warning': 'Only one class present in data'
            }
        
        # 检查是否为三分类数据，如果是则过滤并转换为二分类
        if len(unique_labels) > 2:
            logger.info(f"检测到三分类数据，将过滤为只使用Pathogenic和Benign进行二分类分析")
            # 过滤数据：只保留Pathogenic(2)和Benign(0)的数据
            binary_mask = (y_true == 0) | (y_true == 2)  # 只保留Benign和Pathogenic
            y_true_filtered = y_true[binary_mask]
            y_scores_filtered = y_scores[binary_mask]
            
            logger.info(f"原始数据: 总样本数 {len(y_true)}")
            logger.info(f"过滤后数据: 总样本数 {len(y_true_filtered)}")
            logger.info(f"过滤后类别分布: {dict(zip(*np.unique(y_true_filtered, return_counts=True)))}")
            
            if len(y_true_filtered) == 0:
                logger.warning("过滤后没有Pathogenic/Benign数据，使用默认值")
                return {
                    'fpr': [0.0, 1.0],
                    'tpr': [0.0, 1.0],
                    'thresholds': [float('inf'), 0.0],
                    'auc': 0.5,
                    'optimal_thresholds': {
                        'youden': {'threshold': 0.5, 'fpr': 1.0, 'tpr': 1.0, 'youden_index': 0.0},
                        'f1': {'threshold': 0.5, 'fpr': 1.0, 'tpr': 1.0, 'f1_score': 0.0}
                        },
                    'warning': 'No Pathogenic/Benign data after filtering'
                }
            
            # 将过滤后的数据转换为二分类：0=Benign, 1=Pathogenic
            y_true_binary = (y_true_filtered == 2).astype(int)  # 只有Pathogenic(2)为1，Benign(0)为0
            y_scores = y_scores_filtered  # 使用过滤后的分数
            logger.info(f"二分类转换: Benign({np.sum(y_true_binary == 0)}) vs Pathogenic({np.sum(y_true_binary == 1)})")
        else:
            y_true_binary = y_true
            logger.info(f"使用原始二分类数据: 类别分布 {dict(zip(*np.unique(y_true_binary, return_counts=True)))}")
        
        try:
            # 计算ROC曲线
            fpr, tpr, thresholds = roc_curve(y_true_binary, y_scores)
            roc_auc = auc(fpr, tpr)
            
            # 检查AUC是否有效
            if np.isnan(roc_auc) or np.isinf(roc_auc):
                roc_auc = 0.5
            
            # 计算多种最佳阈值方法
            optimal_thresholds = self._calculate_optimal_thresholds(y_true_binary, y_scores, fpr, tpr, thresholds)
            
            return {
                'fpr': fpr.tolist(),
                'tpr': tpr.tolist(),
                'thresholds': thresholds.tolist(),
                'auc': float(roc_auc),
                'optimal_thresholds': optimal_thresholds
            }
            
        except Exception as e:
            logger.warning(f"ROC calculation failed: {e}")
            return {
                'fpr': [0.0, 1.0],
                'tpr': [0.0, 1.0],
                'thresholds': [float('inf'), 0.0],
                'auc': 0.5,
                'optimal_thresholds': {
                    'youden': {'threshold': 0.5, 'fpr': 1.0, 'tpr': 1.0, 'youden_index': 0.0},
                    'f1': {'threshold': 0.5, 'fpr': 1.0, 'tpr': 1.0, 'f1_score': 0.0}
                },
                'error': str(e)
            }
    
    def _calculate_pr_metrics(self, roc_data: List[Dict]) -> Dict[str, Any]:
        """
        计算PR指标
        """
        if not roc_data:
            return {'error': 'No data for PR calculation'}
        
        # 提取真实标签和预测分数
        y_true = np.array([item['true_label'] for item in roc_data])
        y_scores = np.array([item['prediction_score'] for item in roc_data])
        
        # 检查是否为三分类数据，如果是则过滤并转换为二分类
        unique_labels = np.unique(y_true)
        if len(unique_labels) > 2:
            # 过滤数据：只保留Pathogenic(2)和Benign(0)的数据
            binary_mask = (y_true == 0) | (y_true == 2)  # 只保留Benign和Pathogenic
            y_true_filtered = y_true[binary_mask]
            y_scores_filtered = y_scores[binary_mask]
            
            if len(y_true_filtered) == 0:
                return {
                    'precision': [1.0, 1.0],
                    'recall': [1.0, 0.0],
                    'thresholds': [0.0],
                    'auc': 1.0 if sum(y_true) > 0 else 0.0,
                    'warning': 'No Pathogenic/Benign data after filtering'
                }
            
            # 将过滤后的数据转换为二分类：0=Benign, 1=Pathogenic
            y_true_binary = (y_true_filtered == 2).astype(int)
            y_scores = y_scores_filtered
        else:
            y_true_binary = y_true
        
        try:
            # 计算PR曲线
            precision, recall, pr_thresholds = precision_recall_curve(y_true_binary, y_scores)
            pr_auc = average_precision_score(y_true_binary, y_scores)
            
            # 检查AUC是否有效
            if np.isnan(pr_auc) or np.isinf(pr_auc):
                pr_auc = 1.0 if sum(y_true) > 0 else 0.0
            
            return {
                'precision': precision.tolist(),
                'recall': recall.tolist(),
                'thresholds': pr_thresholds.tolist(),
                'auc': float(pr_auc)
            }
            
        except Exception as e:
            logger.warning(f"PR calculation failed: {e}")
            return {
                'precision': [1.0, 1.0],
                'recall': [1.0, 0.0],
                'thresholds': [0.0],
                'auc': 1.0 if sum(y_true) > 0 else 0.0,
                'error': str(e)
            }
    
    def _calculate_optimal_thresholds(self, y_true: np.ndarray, y_scores: np.ndarray, 
                                    fpr: np.ndarray, tpr: np.ndarray, thresholds: np.ndarray) -> Dict[str, Any]:
        """
        计算多种最佳阈值方法
        
        Args:
            y_true: 真实标签
            y_scores: 预测分数
            fpr: 假正率
            tpr: 真正率
            thresholds: 阈值
            
        Returns:
            包含多种最佳阈值的字典
        """
        from sklearn.metrics import precision_score, recall_score, f1_score
        
        optimal_thresholds = {}
        
        try:
            # 1. Youden指数方法 (TPR - FPR的最大值)
            youden_indices = tpr - fpr
            youden_idx = np.argmax(youden_indices)
            optimal_thresholds['youden'] = {
                'threshold': float(thresholds[youden_idx]),
                'fpr': float(fpr[youden_idx]),
                'tpr': float(tpr[youden_idx]),
                'youden_index': float(youden_indices[youden_idx])
            }
            
            # 2. F1分数方法
            f1_scores = []
            for threshold in thresholds:
                y_pred = (y_scores >= threshold).astype(int)
                if len(np.unique(y_pred)) > 1:  # 确保有正负样本
                    f1 = f1_score(y_true, y_pred, zero_division=0)
                else:
                    f1 = 0.0
                f1_scores.append(f1)
            
            f1_idx = np.argmax(f1_scores)
            optimal_thresholds['f1'] = {
                'threshold': float(thresholds[f1_idx]),
                'fpr': float(fpr[f1_idx]),
                'tpr': float(tpr[f1_idx]),
                'f1_score': float(f1_scores[f1_idx])
            }
            
            # 3. 精确率-召回率平衡方法 (Precision = Recall的点)
            precision_scores = []
            recall_scores = []
            for threshold in thresholds:
                y_pred = (y_scores >= threshold).astype(int)
                if len(np.unique(y_pred)) > 1:
                    precision = precision_score(y_true, y_pred, zero_division=0)
                    recall = recall_score(y_true, y_pred, zero_division=0)
                else:
                    precision = recall = 0.0
                precision_scores.append(precision)
                recall_scores.append(recall)
            
            # 找到精确率和召回率最接近的点
            precision_recall_diff = np.abs(np.array(precision_scores) - np.array(recall_scores))
            pr_idx = np.argmin(precision_recall_diff)
            optimal_thresholds['precision_recall'] = {
                'threshold': float(thresholds[pr_idx]),
                'fpr': float(fpr[pr_idx]),
                'tpr': float(tpr[pr_idx]),
                'precision': float(precision_scores[pr_idx]),
                'recall': float(recall_scores[pr_idx])
            }
            
            # 4. 几何平均值方法 (sqrt(TPR * (1-FPR)))
            geometric_means = np.sqrt(tpr * (1 - fpr))
            gm_idx = np.argmax(geometric_means)
            optimal_thresholds['geometric_mean'] = {
                'threshold': float(thresholds[gm_idx]),
                'fpr': float(fpr[gm_idx]),
                'tpr': float(tpr[gm_idx]),
                'geometric_mean': float(geometric_means[gm_idx])
            }
            
        except Exception as e:
            logger.warning(f"Optimal threshold calculation failed: {e}")
            # 返回默认值
            optimal_thresholds = {
                'youden': {'threshold': 0.5, 'fpr': 1.0, 'tpr': 1.0, 'youden_index': 0.0},
                'f1': {'threshold': 0.5, 'fpr': 1.0, 'tpr': 1.0, 'f1_score': 0.0},
                'precision_recall': {'threshold': 0.5, 'fpr': 1.0, 'tpr': 1.0, 'precision': 0.0, 'recall': 1.0},
                'geometric_mean': {'threshold': 0.5, 'fpr': 1.0, 'tpr': 1.0, 'geometric_mean': 0.0}
            }
        
        return optimal_thresholds
    
    def _calculate_confusion_metrics(self, roc_data: List[Dict]) -> Dict[str, Any]:
        """
        计算混淆矩阵指标
        """
        if not roc_data:
            return {'error': 'No data for confusion matrix calculation'}
        
        # 使用0.5作为阈值
        threshold = 0.5
        y_true = np.array([item['true_label'] for item in roc_data])
        y_scores = np.array([item['prediction_score'] for item in roc_data])
        
        # 检查是否为三分类数据，如果是则过滤并转换为二分类
        unique_labels = np.unique(y_true)
        if len(unique_labels) > 2:
            # 过滤数据：只保留Pathogenic(2)和Benign(0)的数据
            binary_mask = (y_true == 0) | (y_true == 2)  # 只保留Benign和Pathogenic
            y_true_filtered = y_true[binary_mask]
            y_scores_filtered = y_scores[binary_mask]
            
            if len(y_true_filtered) == 0:
                return {
                    'confusion_matrix': [[0, 0], [0, 0]],
                    'classification_report': {},
                    'threshold_used': threshold,
                    'warning': 'No Pathogenic/Benign data after filtering'
                }
            
            # 将过滤后的数据转换为二分类：0=Benign, 1=Pathogenic
            y_true_binary = (y_true_filtered == 2).astype(int)
            y_scores = y_scores_filtered
        else:
            y_true_binary = y_true
        
        # 使用过滤后的数据计算预测标签
        y_pred = np.array([1 if score >= threshold else 0 for score in y_scores])
        
        # 计算混淆矩阵
        cm = confusion_matrix(y_true_binary, y_pred)
        
        # 计算分类报告
        report = classification_report(y_true_binary, y_pred, output_dict=True)
        
        return {
            'confusion_matrix': cm.tolist(),
            'classification_report': report,
            'threshold_used': threshold
        }
    
    def _summarize_data(self, roc_data: List[Dict]) -> Dict[str, Any]:
        """
        总结数据概况
        """
        if not roc_data:
            return {'error': 'No data to summarize'}
        
        total_samples = len(roc_data)
        positive_samples = sum(1 for item in roc_data if item['true_label'] == 1)
        negative_samples = total_samples - positive_samples
        
        # 按临床意义分组
        significance_counts = {}
        for item in roc_data:
            sig = item['clinical_significance']
            significance_counts[sig] = significance_counts.get(sig, 0) + 1
        
        # 按匹配类型分组
        match_type_counts = {}
        for item in roc_data:
            match_type = item['match_type']
            match_type_counts[match_type] = match_type_counts.get(match_type, 0) + 1
        
        # 按基因分组（多基因分析）
        gene_counts = {}
        for item in roc_data:
            gene = item.get('gene_name', 'Unknown')
            gene_counts[gene] = gene_counts.get(gene, 0) + 1
        
        return {
            'total_samples': total_samples,
            'positive_samples': positive_samples,
            'negative_samples': negative_samples,
            'positive_rate': positive_samples / total_samples if total_samples > 0 else 0,
            'significance_distribution': significance_counts,
            'match_type_distribution': match_type_counts,
            'gene_distribution': gene_counts  # 添加基因分布
        }
    
    def _save_roc_results(self, analysis_result: Dict[str, Any], gene_name: str):
        """
        保存ROC分析结果
        """
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # 保存JSON结果
        json_file = f"{self.output_dir}/roc_analysis/{gene_name}_roc_analysis_{timestamp}.json"
        with open(json_file, 'w', encoding='utf-8') as f:
            json.dump(analysis_result, f, ensure_ascii=False, indent=2)
        
        logger.info(f"ROC分析结果已保存到: {json_file}")
    
    def create_roc_plots(self, 
                        analysis_result: Dict[str, Any],
                        gene_name: str) -> Dict[str, str]:
        """
        创建ROC分析图表
        
        Args:
            analysis_result: ROC分析结果
            gene_name: 基因名称
            
        Returns:
            图表文件路径字典
        """
        logger.info(f"为基因 {gene_name} 创建ROC分析图表")
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        plot_files = {}
        
        try:
            # 1. ROC曲线图
            if 'roc_metrics' in analysis_result and 'error' not in analysis_result['roc_metrics']:
                plot_files['roc_curve'] = self._plot_roc_curve(
                    analysis_result['roc_metrics'], gene_name, timestamp
                )
            
            # 2. PR曲线图
            if 'pr_metrics' in analysis_result and 'error' not in analysis_result['pr_metrics']:
                plot_files['pr_curve'] = self._plot_pr_curve(
                    analysis_result['pr_metrics'], gene_name, timestamp
                )
            
            # 3. 混淆矩阵热图
            if 'confusion_metrics' in analysis_result and 'error' not in analysis_result['confusion_metrics']:
                plot_files['confusion_matrix'] = self._plot_confusion_matrix(
                    analysis_result['confusion_metrics'], gene_name, timestamp
                )
            
            # 4. 综合性能图
            plot_files['performance_summary'] = self._plot_performance_summary(
                analysis_result, gene_name, timestamp
            )
            
            logger.info(f"成功创建 {len(plot_files)} 个ROC分析图表")
            
        except Exception as e:
            logger.error(f"创建ROC分析图表失败: {e}")
        
        return plot_files
    
    def _plot_roc_curve(self, roc_metrics: Dict, gene_name: str, timestamp: str) -> str:
        """
        绘制ROC曲线
        """
        plt.figure(figsize=(10, 8))
        
        fpr = roc_metrics['fpr']
        tpr = roc_metrics['tpr']
        auc_score = roc_metrics['auc']
        
        # 绘制ROC曲线
        plt.plot(fpr, tpr, color='darkorange', lw=2, 
                label=f'ROC Curve (AUC = {auc_score:.3f})')
        
        # 绘制对角线
        plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--', 
                label='Random Classifier')
        
        # 标记最佳阈值点（使用Youden指数方法作为主要显示）
        optimal_thresholds = roc_metrics.get('optimal_thresholds', {})
        if 'youden' in optimal_thresholds:
            youden = optimal_thresholds['youden']
            plt.plot(youden['fpr'], youden['tpr'], 'ro', markersize=8, 
                    label=f'Optimal Threshold (FPR={youden["fpr"]:.3f}, TPR={youden["tpr"]:.3f})')
        
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate (FPR)', fontsize=12)
        plt.ylabel('True Positive Rate (TPR)', fontsize=12)
        plt.title(f'{gene_name} ROC Curve Analysis', fontsize=16, fontweight='bold')
        plt.legend(loc="lower right")
        plt.grid(alpha=0.3)
        
        filename = f"{self.output_dir}/plots/{gene_name}_roc_curve_{timestamp}.png"
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close()
        
        return filename
    
    def _plot_pr_curve(self, pr_metrics: Dict, gene_name: str, timestamp: str) -> str:
        """
        绘制PR曲线
        """
        plt.figure(figsize=(10, 8))
        
        precision = pr_metrics['precision']
        recall = pr_metrics['recall']
        auc_score = pr_metrics['auc']
        
        # 绘制PR曲线
        plt.plot(recall, precision, color='darkgreen', lw=2, 
                label=f'PR Curve (AUC = {auc_score:.3f})')
        
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('Recall', fontsize=12)
        plt.ylabel('Precision', fontsize=12)
        plt.title(f'{gene_name} PR Curve Analysis', fontsize=16, fontweight='bold')
        plt.legend(loc="lower left")
        plt.grid(alpha=0.3)
        
        filename = f"{self.output_dir}/plots/{gene_name}_pr_curve_{timestamp}.png"
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close()
        
        return filename
    
    def _plot_confusion_matrix(self, confusion_metrics: Dict, gene_name: str, timestamp: str) -> str:
        """
        绘制混淆矩阵热图
        """
        plt.figure(figsize=(8, 6))
        
        cm = np.array(confusion_metrics['confusion_matrix'])
        
        # 创建热图
        sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', 
                   xticklabels=['Predicted Negative', 'Predicted Positive'],
                   yticklabels=['Actual Negative', 'Actual Positive'])
        
        plt.title(f'{gene_name} Confusion Matrix', fontsize=16, fontweight='bold')
        plt.xlabel('Predicted Label', fontsize=12)
        plt.ylabel('Actual Label', fontsize=12)
        
        filename = f"{self.output_dir}/plots/{gene_name}_confusion_matrix_{timestamp}.png"
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close()
        
        return filename
    
    def _plot_performance_summary(self, analysis_result: Dict, gene_name: str, timestamp: str) -> str:
        """
        绘制综合性能摘要图
        """
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
        
        # 1. AUC比较
        roc_auc = analysis_result.get('roc_metrics', {}).get('auc', 0)
        pr_auc = analysis_result.get('pr_metrics', {}).get('auc', 0)
        
        metrics = ['ROC AUC', 'PR AUC']
        values = [roc_auc, pr_auc]
        colors = ['skyblue', 'lightcoral']
        
        bars1 = ax1.bar(metrics, values, color=colors, edgecolor='black')
        ax1.set_title('AUC Performance Comparison', fontweight='bold')
        ax1.set_ylabel('AUC Value')
        ax1.set_ylim(0, 1)
        
        # 添加数值标签
        for bar, value in zip(bars1, values):
            ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01, 
                    f'{value:.3f}', ha='center', va='bottom', fontweight='bold')
        
        # 2. 数据分布
        data_summary = analysis_result.get('data_summary', {})
        positive_rate = data_summary.get('positive_rate', 0)
        
        ax2.pie([positive_rate, 1-positive_rate], 
               labels=['Positive Samples', 'Negative Samples'], 
               colors=['lightcoral', 'lightblue'],
               autopct='%1.1f%%', startangle=90)
        ax2.set_title('Sample Distribution', fontweight='bold')
        
        # 3. 临床意义分布
        sig_dist = data_summary.get('significance_distribution', {})
        if sig_dist:
            categories = list(sig_dist.keys())
            counts = list(sig_dist.values())
            
            bars3 = ax3.bar(categories, counts, color='lightgreen', edgecolor='black')
            ax3.set_title('Clinical Significance Distribution', fontweight='bold')
            ax3.set_ylabel('Count')
            ax3.tick_params(axis='x', rotation=45)
            
            # 添加数值标签
            for bar, count in zip(bars3, counts):
                ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1, 
                        str(count), ha='center', va='bottom', fontweight='bold')
        
        # 4. 匹配类型分布
        match_dist = data_summary.get('match_type_distribution', {})
        if match_dist:
            categories = list(match_dist.keys())
            counts = list(match_dist.values())
            
            bars4 = ax4.bar(categories, counts, color='orange', edgecolor='black')
            ax4.set_title('Match Type Distribution', fontweight='bold')
            ax4.set_ylabel('Count')
            
            # 添加数值标签
            for bar, count in zip(bars4, counts):
                ax4.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1, 
                        str(count), ha='center', va='bottom', fontweight='bold')
        
        plt.suptitle(f'{gene_name} ROC Analysis Comprehensive Performance Summary', fontsize=16, fontweight='bold')
        plt.tight_layout()
        
        filename = f"{self.output_dir}/plots/{gene_name}_performance_summary_{timestamp}.png"
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close()
        
        return filename
    
    def perform_multiclass_roc_analysis(self, 
                                       comparison_result: Dict[str, Any],
                                       gene_name: str,
                                       save_results: bool = True) -> Dict[str, Any]:
        """
        执行多分类ROC分析（三分类：Benign, Uncertain significance, Pathogenic）
        
        Args:
            comparison_result: ClinVar比较结果
            gene_name: 基因名称
            save_results: 是否保存结果
            
        Returns:
            多分类ROC分析结果
        """
        logger.info(f"开始对基因 {gene_name} 进行多分类ROC分析")
        
        try:
            # 准备ROC数据
            roc_data = self._prepare_roc_data(comparison_result)
            
            if not roc_data:
                logger.warning(f"基因 {gene_name} 没有可用的ROC数据")
                return {'error': 'No ROC data available'}
            
            # 提取数据
            true_labels = np.array([item['true_label'] for item in roc_data])
            prediction_scores = np.array([item['prediction_score'] for item in roc_data])
            clinical_significances = [item['clinical_significance'] for item in roc_data]
            
            # 统计各类别数量
            class_counts = np.bincount(true_labels)
            class_names = ['Benign', 'Uncertain significance', 'Pathogenic']
            
            logger.info(f"类别分布: {dict(zip(class_names, class_counts))}")
            
            # 计算多分类ROC曲线
            from sklearn.metrics import roc_curve, auc
            from sklearn.preprocessing import label_binarize
            
            # 二值化标签
            y_bin = label_binarize(true_labels, classes=[0, 1, 2])
            n_classes = y_bin.shape[1]
            
            # 计算每个类别的ROC曲线
            fpr = dict()
            tpr = dict()
            roc_auc = dict()
            
            for i in range(n_classes):
                fpr[i], tpr[i], _ = roc_curve(y_bin[:, i], prediction_scores)
                roc_auc[i] = auc(fpr[i], tpr[i])
            
            # 计算微平均和宏平均
            from sklearn.metrics import roc_curve, auc
            from itertools import cycle
            
            # 微平均
            fpr_micro, tpr_micro, _ = roc_curve(y_bin.ravel(), np.tile(prediction_scores, n_classes))
            roc_auc_micro = auc(fpr_micro, tpr_micro)
            
            # 宏平均
            all_fpr = np.unique(np.concatenate([fpr[i] for i in range(n_classes)]))
            mean_tpr = np.zeros_like(all_fpr)
            for i in range(n_classes):
                mean_tpr += np.interp(all_fpr, fpr[i], tpr[i])
            mean_tpr /= n_classes
            roc_auc_macro = auc(all_fpr, mean_tpr)
            
            # 创建结果
            result = {
                'gene_name': gene_name,
                'analysis_type': 'multiclass_roc',
                'class_names': class_names,
                'class_counts': class_counts.tolist(),
                'roc_metrics': {
                    'fpr': {str(i): fpr[i].tolist() for i in range(n_classes)},
                    'tpr': {str(i): tpr[i].tolist() for i in range(n_classes)},
                    'auc': {class_names[i]: roc_auc[i] for i in range(n_classes)},
                    'micro_auc': roc_auc_micro,
                    'macro_auc': roc_auc_macro
                },
                'total_samples': len(roc_data),
                'timestamp': datetime.now().isoformat()
            }
            
            # 保存结果
            if save_results:
                result_file = f"{self.output_dir}/roc_analysis/{gene_name}_multiclass_roc_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
                with open(result_file, 'w', encoding='utf-8') as f:
                    json.dump(result, f, ensure_ascii=False, indent=2)
                logger.info(f"多分类ROC分析结果已保存到: {result_file}")
            
            logger.info(f"多分类ROC分析完成")
            logger.info(f"各类别AUC: {dict(zip(class_names, [roc_auc[i] for i in range(n_classes)]))}")
            logger.info(f"微平均AUC: {roc_auc_micro:.3f}")
            logger.info(f"宏平均AUC: {roc_auc_macro:.3f}")
            
            return result
            
        except Exception as e:
            logger.error(f"多分类ROC分析失败: {e}")
            return {'error': str(e)}
    
    def create_multiclass_roc_plots(self, 
                                   roc_result: Dict[str, Any], 
                                   gene_name: str) -> List[str]:
        """
        创建多分类ROC曲线图表
        
        Args:
            roc_result: ROC分析结果
            gene_name: 基因名称
            
        Returns:
            生成的图表文件列表
        """
        logger.info(f"为基因 {gene_name} 创建多分类ROC分析图表")
        
        try:
            if 'error' in roc_result:
                logger.error(f"无法创建图表: {roc_result['error']}")
                return []
            
            roc_metrics = roc_result.get('roc_metrics', {})
            class_names = roc_result.get('class_names', ['Benign', 'Uncertain significance', 'Pathogenic'])
            
            if not roc_metrics:
                logger.warning("没有可用的ROC指标数据")
                return []
            
            # 创建图表
            fig, axes = plt.subplots(2, 2, figsize=(15, 12))
            fig.suptitle(f'Multi-class ROC Analysis for {gene_name}', fontsize=16, fontweight='bold')
            
            # 1. 各类别ROC曲线
            ax1 = axes[0, 0]
            colors = ['blue', 'orange', 'green']
            
            for i, class_name in enumerate(class_names):
                if str(i) in roc_metrics['fpr'] and str(i) in roc_metrics['tpr']:
                    fpr = roc_metrics['fpr'][str(i)]
                    tpr = roc_metrics['tpr'][str(i)]
                    auc_score = roc_metrics['auc'][class_name]
                    
                    ax1.plot(fpr, tpr, color=colors[i], linewidth=2,
                            label=f'{class_name} (AUC = {auc_score:.3f})')
            
            ax1.plot([0, 1], [0, 1], 'k--', alpha=0.5, linewidth=1, label='Random')
            ax1.set_xlabel('False Positive Rate', fontweight='bold')
            ax1.set_ylabel('True Positive Rate', fontweight='bold')
            ax1.set_title('Multi-class ROC Curves', fontweight='bold')
            ax1.legend()
            ax1.grid(alpha=0.3)
            ax1.set_xlim([0, 1])
            ax1.set_ylim([0, 1.05])
            
            # 2. 微平均ROC曲线
            ax2 = axes[0, 1]
            if 'micro_auc' in roc_metrics:
                fpr_micro = roc_metrics['fpr'].get('micro', [])
                tpr_micro = roc_metrics['tpr'].get('micro', [])
                
                if fpr_micro and tpr_micro:
                    ax2.plot(fpr_micro, tpr_micro, color='red', linewidth=2,
                            label=f'Micro-average (AUC = {roc_metrics["micro_auc"]:.3f})')
            
            ax2.plot([0, 1], [0, 1], 'k--', alpha=0.5, linewidth=1, label='Random')
            ax2.set_xlabel('False Positive Rate', fontweight='bold')
            ax2.set_ylabel('True Positive Rate', fontweight='bold')
            ax2.set_title('Micro-average ROC Curve', fontweight='bold')
            ax2.legend()
            ax2.grid(alpha=0.3)
            ax2.set_xlim([0, 1])
            ax2.set_ylim([0, 1.05])
            
            # 3. 宏平均ROC曲线
            ax3 = axes[1, 0]
            if 'macro_auc' in roc_metrics:
                fpr_macro = roc_metrics['fpr'].get('macro', [])
                tpr_macro = roc_metrics['tpr'].get('macro', [])
                
                if fpr_macro and tpr_macro:
                    ax3.plot(fpr_macro, tpr_macro, color='purple', linewidth=2,
                            label=f'Macro-average (AUC = {roc_metrics["macro_auc"]:.3f})')
            
            ax3.plot([0, 1], [0, 1], 'k--', alpha=0.5, linewidth=1, label='Random')
            ax3.set_xlabel('False Positive Rate', fontweight='bold')
            ax3.set_ylabel('True Positive Rate', fontweight='bold')
            ax3.set_title('Macro-average ROC Curve', fontweight='bold')
            ax3.legend()
            ax3.grid(alpha=0.3)
            ax3.set_xlim([0, 1])
            ax3.set_ylim([0, 1.05])
            
            # 4. AUC比较柱状图
            ax4 = axes[1, 1]
            auc_scores = [roc_metrics['auc'][class_name] for class_name in class_names]
            bars = ax4.bar(class_names, auc_scores, color=colors, alpha=0.7)
            ax4.set_ylabel('AUC Score', fontweight='bold')
            ax4.set_title('AUC Scores by Class', fontweight='bold')
            ax4.set_ylim([0, 1])
            
            # 添加数值标签
            for bar, score in zip(bars, auc_scores):
                ax4.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                        f'{score:.3f}', ha='center', va='bottom', fontweight='bold')
            
            # 添加微平均和宏平均
            ax4.axhline(y=roc_metrics.get('micro_auc', 0), color='red', linestyle='--', 
                       label=f'Micro-avg: {roc_metrics.get("micro_auc", 0):.3f}')
            ax4.axhline(y=roc_metrics.get('macro_auc', 0), color='purple', linestyle='--', 
                       label=f'Macro-avg: {roc_metrics.get("macro_auc", 0):.3f}')
            ax4.legend()
            ax4.grid(alpha=0.3)
            
            plt.tight_layout()
            
            # 保存图表
            plot_file = f"{self.output_dir}/plots/{gene_name}_multiclass_roc_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
            plt.savefig(plot_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            logger.info(f"成功创建多分类ROC分析图表: {plot_file}")
            return [plot_file]
            
        except Exception as e:
            logger.error(f"创建多分类ROC图表失败: {e}")
            return []

# 使用示例
if __name__ == "__main__":
    # 设置日志
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    
    # 创建ROC分析器
    roc_analyzer = ROCAnalyzer()
    
    # 模拟比较结果
    sample_comparison = {
        'matches': [
            {
                'simulated_position': 100,
                'simulated_ref': 'A',
                'simulated_alt': 'G',
                'clinical_significance': 'Pathogenic',
                'match_type': 'exact'
            },
            {
                'simulated_position': 200,
                'simulated_ref': 'C',
                'simulated_alt': 'T',
                'clinical_significance': 'Benign',
                'match_type': 'position_only'
            }
        ]
    }
    
    # 执行ROC分析
    result = roc_analyzer.perform_roc_analysis(sample_comparison, "BRCA1")
    print("ROC分析结果:", json.dumps(result, indent=2, ensure_ascii=False))
    
    # 创建图表
    plots = roc_analyzer.create_roc_plots(result, "BRCA1")
    print("生成的图表:", plots)
