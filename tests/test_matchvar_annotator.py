#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MATCHVAR注释器测试
"""

import pytest
import os
import tempfile
import pandas as pd
from unittest.mock import Mock, patch

from matchvar_annotator import MatchvarRunner


class TestMatchvarRunner:
    """测试MatchvarRunner类"""
    
    def setup_method(self):
        """测试前准备"""
        self.temp_dir = tempfile.mkdtemp()
        self.runner = MatchvarRunner(
            resources_dir=self.temp_dir,
            genome_version="hg19",
            thread_count=2
        )
    
    def teardown_method(self):
        """测试后清理"""
        import shutil
        shutil.rmtree(self.temp_dir, ignore_errors=True)
    
    def test_init(self):
        """测试初始化"""
        assert self.runner.genome_version == "hg19"
        assert self.runner.thread_count == 2
        assert os.path.exists(self.runner.resources_dir)
        assert os.path.exists(self.runner.path_humandb)
        assert os.path.exists(self.runner.input_files)
    
    def test_get_protocol_categories(self):
        """测试获取协议分类"""
        categories = self.runner.get_protocol_categories()
        assert isinstance(categories, dict)
        assert 'gene_info' in categories
        assert 'database' in categories
        assert 'prediction' in categories
    
    def test_validate_protocols(self):
        """测试协议验证"""
        # 有效协议
        valid_protocols = ['refGene', 'exac03', 'avsift']
        is_valid, errors = self.runner.validate_protocols(valid_protocols)
        assert is_valid
        assert len(errors) == 0
        
        # 无效协议
        invalid_protocols = ['refGene', 'invalid_protocol']
        is_valid, errors = self.runner.validate_protocols(invalid_protocols)
        assert not is_valid
        assert len(errors) > 0
    
    def test_add_custom_protocol(self):
        """测试添加自定义协议"""
        success = self.runner.add_custom_protocol(
            name="test_protocol",
            operation="f",
            description="Test protocol",
            category="test"
        )
        assert success
        
        # 重复添加应该失败
        success = self.runner.add_custom_protocol(
            name="test_protocol",
            operation="f", 
            description="Test protocol",
            category="test"
        )
        assert not success
    
    def test_remove_protocol(self):
        """测试移除协议"""
        # 添加协议
        self.runner.add_custom_protocol(
            name="test_protocol",
            operation="f",
            description="Test protocol", 
            category="test"
        )
        
        # 移除协议
        success = self.runner.remove_protocol("test_protocol")
        assert success
        
        # 移除不存在的协议应该失败
        success = self.runner.remove_protocol("nonexistent_protocol")
        assert not success
    
    @patch('subprocess.run')
    def test_build_matchvar_command(self, mock_subprocess):
        """测试构建MATCHVAR命令"""
        protocols = ['refGene', 'exac03']
        command = self.runner.build_matchvar_command(
            input_file="test.vcf",
            protocols=protocols,
            buildver="hg19",
            output_prefix="test_output"
        )
        
        assert isinstance(command, str)
        assert "refGene" in command
        assert "exac03" in command
        assert "test.vcf" in command
        assert "test_output" in command
    
    def test_get_custom_databases(self):
        """测试获取自定义数据库"""
        custom_dbs = self.runner.get_custom_databases()
        assert isinstance(custom_dbs, dict)
        # 由于是空目录，应该没有自定义数据库
        assert len(custom_dbs) == 0


class TestPackageImport:
    """测试包导入"""
    
    def test_import_package(self):
        """测试包导入"""
        import matchvar_annotator
        assert hasattr(matchvar_annotator, 'MatchvarRunner')
        assert hasattr(matchvar_annotator, 'TableAnnotator')
        assert hasattr(matchvar_annotator, 'Convert2Matchvar')
        assert hasattr(matchvar_annotator, 'CodingChange')
    
    def test_package_version(self):
        """测试包版本"""
        import matchvar_annotator
        assert hasattr(matchvar_annotator, '__version__')
        assert matchvar_annotator.__version__ is not None
    
    def test_package_metadata(self):
        """测试包元数据"""
        import matchvar_annotator
        assert hasattr(matchvar_annotator, '__author__')
        assert hasattr(matchvar_annotator, '__email__')
        assert hasattr(matchvar_annotator, '__description__')


if __name__ == '__main__':
    pytest.main([__file__])
