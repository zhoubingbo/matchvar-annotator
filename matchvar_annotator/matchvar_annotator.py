import subprocess
import pandas as pd
import os
import sys
import locale
import psutil
from typing import List, Dict, Optional, Tuple
from datetime import datetime


def _detect_python_executable() -> str:
    """Detect a suitable Python interpreter to invoke child scripts.
    Priority:
    1) Environment variable PYTHON_EXECUTABLE
    2) Project local .venv (Windows/Linux/macOS)
    3) Current interpreter (sys.executable)
    """
    # 1) Explicit env var
    env_py = os.environ.get('PYTHON_EXECUTABLE')
    if env_py and os.path.exists(env_py):
        return env_py

    # 2) Local .venv
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
    venv_unix = os.path.join(project_root, '.venv', 'bin', 'python')
    venv_win = os.path.join(project_root, '.venv', 'Scripts', 'python.exe')
    if os.name == 'nt' and os.path.exists(venv_win):
        return venv_win
    if os.path.exists(venv_unix):
        return venv_unix

    # 3) Fallback to current interpreter
    return sys.executable

PYTHON_EXECUTABLE = _detect_python_executable()

def get_system_encoding():
    """Get the system encoding, ensuring cross-platform compatibility"""
    try:
        # Try to get the system default encoding
        system_encoding = locale.getpreferredencoding()
        # If the encoding is not available, use UTF-8
        if system_encoding.lower() in ['cp1252', 'gbk', 'gb2312']:
            return 'utf-8'
        return system_encoding
    except:
        return 'utf-8'

def run_subprocess_safe(command: str, **kwargs):
    """Run subprocess safely, handle encoding issues, and support real-time output"""
    # Set default encoding parameters
    default_kwargs = {
        'shell': True,
        'text': True,
        'encoding': get_system_encoding(),
        'errors': 'replace'  # Replace with a placeholder when encountering undecodable characters
    }
    
    # Check if real-time output is needed
    real_time_output = kwargs.pop('real_time_output', True)
    
    if real_time_output:
        # Real-time output mode: do not capture output, display directly
        default_kwargs.update({
            'stdout': None,  # Directly output to the console
            'stderr': None,  # Directly output to the console
        })
    else:
        # Capture output mode: used for scenarios that need to return results
        default_kwargs.update({
            'capture_output': True,
        })
    
    # Merge other parameters
    default_kwargs.update(kwargs)
    
    try:
        return subprocess.run(command, **default_kwargs)
    except UnicodeDecodeError:
        # If the default encoding fails, try UTF-8
        default_kwargs['encoding'] = 'utf-8'
        return subprocess.run(command, **default_kwargs)

class MatchvarRunner:
    """Improved MATCHVAR annotator, supporting dynamic protocol parameters"""
    
    # Predefined protocol configurations - simplified to three main categories
    PROTOCOL_CONFIGS = {
        # Gene information related
        'refGene': {
            'operation': 'g',
            'description': 'RefSeq gene annotation',
            'category': 'gene_info'
        },
        'ensGene': {
            'operation': 'g', 
            'description': 'Ensembl gene annotation',
            'category': 'gene_info'
        },
        'knownGene': {
            'operation': 'g',
            'description': 'UCSC known gene annotation', 
            'category': 'gene_info'
        },
        
        # 频率数据库
        'exac03': {
            'operation': 'f',
            'description': 'ExAC exome frequency data',
            'category': 'database'
        },
        'gnomad211_genome': {
            'operation': 'f',
            'description': 'gnomAD genome frequency data',
            'category': 'database'
        },
        'esp6500siv2_all': {
            'operation': 'f',
            'description': 'ESP6500 population frequency data',
            'category': 'database'
        },
        '1000g2015aug_all': {
            'operation': 'f',
            'description': '1000 genome frequency data',
            'category': 'database'
        },
        'cytoBand': {
            'operation': 'r',
            'description': 'Cytoband information',
            'category': 'database'
        },
        'clinvar': {
            'operation': 'f',
            'description': 'ClinVar clinical variant database',
            'category': 'database'
        },
        'cosmic70': {
            'operation': 'f',
            'description': 'COSMIC cancer variant database',
            'category': 'database'
        },
        
        # 功能预测软件
        'avsift': {
            'operation': 'f',
            'description': 'SIFT function prediction',
            'category': 'prediction'
        },
        'dbnsfp42a': {
            'operation': 'f',
            'description': 'dbNSFP functional prediction',
            'category': 'prediction'
        },
        'revel': {
            'operation': 'f',
            'description': 'REVEL pathogenicity prediction',
            'category': 'prediction'
        },
        'cadd13gt10': {
            'operation': 'f',
            'description': 'CADD deleterious prediction',
            'category': 'prediction'
        },
        'AlphaMissense': {
            'operation': 'f',
            'description': 'AlphaMissense pathogenicity prediction',
            'category': 'prediction'
        }
    }

    def __init__(self, resources_dir: str = None, input_files_dir: str = None, genome_version: str = 'hg19', thread_count: int = 4):
        """Initialize MATCHVAR runner with custom resources directory and genome version support"""
        if resources_dir:
            self.resources_dir = resources_dir
        else:
            # 获取项目根目录
            current_dir = os.path.dirname(os.path.abspath(__file__))
            project_root = os.path.dirname(os.path.dirname(current_dir))  # 从utils/matchvar回到项目根目录
            self.resources_dir = os.path.join(project_root, 'resources')
        
        self.path_humandb = os.path.join(self.resources_dir, 'humandb')
        
        # 使用传入的输入文件目录，如果没有传入则使用默认值
        self.input_files = input_files_dir or os.path.join(self.resources_dir, 'input_files')
        
        # MATCHVAR文件在utils/matchvar目录中，不在resources目录中
        self.matchvar_path = os.path.dirname(__file__)  # 当前文件就在utils/matchvar目录中
        self.python_table_matchvar = os.path.join(self.matchvar_path, 'table_matchvar.py')
        
        # 设置基因组版本和线程数
        self.genome_version = genome_version
        self.thread_count = int(thread_count) if thread_count else 4
        print(f"初始化MATCHVAR注释器，基因组版本: {self.genome_version}, 线程数: {self.thread_count}")
        
        # 确保目录存在
        os.makedirs(self.path_humandb, exist_ok=True)
        os.makedirs(self.input_files, exist_ok=True)
        
        # 检查自定义数据库目录
        self.custom_db_dir = os.path.join(self.resources_dir, 'custom_databases')
        os.makedirs(self.custom_db_dir, exist_ok=True)

    def get_custom_databases(self) -> Dict[str, Dict]:
        """Get custom databases configuration based on genome version"""
        custom_databases = {}
        
        # 检查custom_databases目录中的自定义数据库
        if os.path.exists(self.custom_db_dir):
            for filename in os.listdir(self.custom_db_dir):
                # 根据基因组版本选择对应的数据库文件
                expected_prefix = f"{self.genome_version}_"
                if filename.startswith(expected_prefix) and filename.endswith('.txt'):
                    # 检查是否是自定义数据库（不在标准协议配置中）
                    db_name = filename[len(expected_prefix):-4]  # 去掉基因组版本前缀和 '.txt'
                    if db_name not in self.PROTOCOL_CONFIGS:
                        file_path = os.path.join(self.custom_db_dir, filename)
                        
                        # 验证文件格式是否正确（至少包含Chr, Start, End, Ref, Alt列）
                        try:
                            with open(file_path, 'r', encoding='utf-8') as f:
                                first_line = f.readline().strip()
                                if first_line.startswith('#'):
                                    # 跳过注释行，读取第二行
                                    second_line = f.readline().strip()
                                    if second_line:
                                        # 检查是否有足够的列（至少5列：Chr, Start, End, Ref, Alt）
                                        cols = second_line.split('\t')
                                        if len(cols) >= 5:
                                            custom_databases[db_name] = {
                                                'operation': 'f',  # 默认使用频率操作
                                                'description': f'Custom Database: {db_name} ({self.genome_version})',
                                                'category': 'custom_database',
                                                'is_custom': True,
                                                'filename': filename,
                                                'file_path': file_path
                                            }
                        except Exception as e:
                            print(f"验证自定义数据库文件 {filename} 时出错: {e}")
                            continue
        
        # 也检查humandb目录中的自定义数据库（向后兼容）
        if os.path.exists(self.path_humandb):
            for filename in os.listdir(self.path_humandb):
                # 根据基因组版本选择对应的数据库文件
                expected_prefix = f"{self.genome_version}_"
                if filename.startswith(expected_prefix) and filename.endswith('.txt'):
                    # 检查是否是自定义数据库（不在标准协议配置中）
                    db_name = filename[len(expected_prefix):-4]  # 去掉基因组版本前缀和 '.txt'
                    if db_name not in self.PROTOCOL_CONFIGS and db_name not in custom_databases:
                        file_path = os.path.join(self.path_humandb, filename)
                        
                        # 验证文件格式是否正确（至少包含Chr, Start, End, Ref, Alt列）
                        try:
                            with open(file_path, 'r', encoding='utf-8') as f:
                                first_line = f.readline().strip()
                                if first_line.startswith('#'):
                                    # 跳过注释行，读取第二行
                                    second_line = f.readline().strip()
                                    if second_line:
                                        # 检查是否有足够的列（至少5列：Chr, Start, End, Ref, Alt）
                                        cols = second_line.split('\t')
                                        if len(cols) >= 5:
                                            custom_databases[db_name] = {
                                                'operation': 'f',  # 默认使用频率操作
                                                'description': f'Custom Database: {db_name} ({self.genome_version})',
                                                'category': 'custom_database',
                                                'is_custom': True,
                                                'filename': filename,
                                                'file_path': file_path
                                            }
                        except Exception as e:
                            print(f"验证自定义数据库文件 {filename} 时出错: {e}")
                            continue
        
        print(f"找到 {len(custom_databases)} 个自定义数据库文件 (基因组版本: {self.genome_version})")
        return custom_databases

    def get_enhanced_protocol_configs(self) -> Dict[str, Dict]:
        """Get enhanced protocol configurations including custom databases"""
        protocols = self.PROTOCOL_CONFIGS.copy()
        custom_databases = self.get_custom_databases()
        protocols.update(custom_databases)
        return protocols

    def validate_protocols(self, protocols: List[str]) -> Tuple[bool, List[str]]:
        """Verify if the protocol list is valid, including custom databases"""
        errors = []
        enhanced_configs = self.get_enhanced_protocol_configs()
        
        for protocol in protocols:
            if protocol not in enhanced_configs:
                errors.append(f"Unknown protocol: {protocol}")
        
        return len(errors) == 0, errors

    def build_matchvar_command(self, input_file: str, protocols: List[str], 
                            buildver: str = 'hg19', output_prefix: str = 'matchvar_result',
                            additional_args: Dict = None) -> str:
        """
        Build the MATCHVAR command with support for custom databases and MANE transcript filtering
        """
        # Verify protocols
        is_valid, errors = self.validate_protocols(protocols)
        if not is_valid:
            raise ValueError(f"Protocol验证失败: {'; '.join(errors)}")
        
        # Get enhanced configurations
        enhanced_configs = self.get_enhanced_protocol_configs()
        
        # Build protocol parameters
        protocol_str = ','.join(protocols)
        
        # Build operation parameters
        operations = []
        for protocol in protocols:
            operations.append(enhanced_configs[protocol]['operation'])
        operation_str = ','.join(operations)
        
        # Build argument parameters (for most protocols, argument is empty)
        arguments = ['' for _ in protocols]  # 使用空字符串作为默认参数
        argument_str = ','.join(arguments)
        
        # Ensure all paths are absolute paths
        python_table_matchvar_abs = os.path.abspath(self.python_table_matchvar)
        input_file_abs = os.path.abspath(os.path.join(self.input_files, input_file))
        humandb_abs = os.path.abspath(self.path_humandb)
        
        cmd_parts = [
            PYTHON_EXECUTABLE,
            python_table_matchvar_abs,
            input_file_abs,
            humandb_abs,
            '-outfile', output_prefix,
            '-buildver', buildver,
            '-protocol', protocol_str,
            '-operation', operation_str
        ]
        
        # 添加argument参数（如果非空）
        if argument_str and argument_str != '.':
            cmd_parts.extend(['-argument', argument_str])
        
        # Add MANE transcript filtering if requested
        print(f"=== matchvar_annotator.py 调试 ===")
        print(f"additional_args: {additional_args}")
        print(f"use_mane_transcript: {additional_args.get('use_mane_transcript') if additional_args else None}")
        if additional_args and additional_args.get('use_mane_transcript'):
            mane_file = os.path.join(humandb_abs, 'mane_transcript.txt')
            if os.path.exists(mane_file):
                cmd_parts.extend(['-mane_file', mane_file])
                cmd_parts.extend(['-use_mane_transcript'])
                print(f"✅ 添加MANE转录本过滤: {mane_file}")
            else:
                print(f"❌ 警告: MANE转录本文件不存在: {mane_file}")
        else:
            print(f"ℹ️  未启用MANE转录本过滤")
        
        # Add otherinfo parameter for complete annotation information
        cmd_parts.append('-otherinfo')
        print("添加-otherinfo参数以获取完整注释信息")
        
        # Add thread count parameter
        cmd_parts.extend(['-thread', str(self.thread_count)])
        print(f"添加线程数参数: {self.thread_count}")
        
        # 检查是否有自定义数据库，如果有，需要特殊处理
        custom_db_protocols = []
        for protocol in protocols:
            if enhanced_configs[protocol].get('is_custom'):
                custom_db_protocols.append(protocol)
        
        if custom_db_protocols:
            print(f"发现自定义数据库协议: {custom_db_protocols}")
            # 将自定义数据库文件复制到humandb目录中，以便MATCHVAR能够处理
            self._prepare_custom_databases_for_matchvar(custom_db_protocols, enhanced_configs)
        
        # Add additional arguments if provided
        if additional_args:
            for key, value in additional_args.items():
                if isinstance(value, bool):
                    if value:
                        cmd_parts.extend([f'-{key}'])
                else:
                    cmd_parts.extend([f'-{key}', str(value)])
        
        return ' '.join(cmd_parts)
    
    def _prepare_custom_databases_for_matchvar(self, custom_db_protocols: List[str], enhanced_configs: Dict[str, Dict]):
        """
        为MATCHVAR准备自定义数据库文件
        
        Args:
            custom_db_protocols: 自定义数据库协议列表
            enhanced_configs: 增强的协议配置
        """
        try:
            for protocol in custom_db_protocols:
                config = enhanced_configs[protocol]
                source_file = config['file_path']
                target_file = os.path.join(self.path_humandb, config['filename'])
                
                # 如果目标文件不存在或源文件更新，则复制文件
                if not os.path.exists(target_file) or os.path.getmtime(source_file) > os.path.getmtime(target_file):
                    import shutil
                    shutil.copy2(source_file, target_file)
                    print(f"已将自定义数据库 {protocol} 复制到 humandb 目录: {target_file}")
                else:
                    print(f"自定义数据库 {protocol} 已存在于 humandb 目录中")
                    
        except Exception as e:
            print(f"准备自定义数据库时出错: {e}")
            import traceback
            traceback.print_exc()
    
    def run_matchvar(self, input_file: str, protocols: List[str] = None, 
                   buildver: str = 'hg19', output_prefix: str = 'matchvar_result',
                   additional_args: Dict = None, timeout: int = 600) -> Optional[pd.DataFrame]:
        """
        Run MATCHVAR annotation
        
        Args:
            input_file: input file name
            protocols: protocol list, if None, use default configuration
            buildver: genome version
            output_prefix: output file prefix
            additional_args: additional arguments (including use_mane_transcript)
            
        Returns:
            Annotation result DataFrame, return None if failed
        """
        # Use default protocols if not specified
        if protocols is None:
            protocols = ['refGene', 'ensGene', 'knownGene', 'cytoBand', 'exac03', 'avsift', 'dbnsfp42a', 
                       'dbscsnv11', 'gnomad211_genome', 'esp6500siv2_all', 
                       'revel', 'cadd13gt10', 'AlphaMissense']
        
        # Log MANE transcript setting
        if additional_args and additional_args.get('use_mane_transcript'):
            print("启用MANE转录本过滤模式")
        else:
            print("使用所有转录本注释模式")
        
        try:
            # Check the input file type
            input_path = os.path.join(self.input_files, input_file)
            is_vcf = input_file.lower().endswith('.vcf')
            
            # 清理可能存在的.orig文件，避免重命名冲突
            self._cleanup_orig_files(output_prefix)
            
            # Convert VCF to MATCHVAR format if needed
            if is_vcf:
                print(f"Detected VCF file: {input_file}")
                converted_file = self._convert_vcf_to_matchvar_input(input_file)
                if not converted_file:
                    print(f"VCF conversion failed for {input_file}")
                    return None
                input_file = converted_file
                print(f"Using converted file: {input_file}")
            
            # Build MATCHVAR command
            command = self.build_matchvar_command(input_file, protocols, buildver, output_prefix, additional_args)
            print(f"Executing MATCHVAR command: {command}")
            
            # Execute MATCHVAR with process management
            result = self._execute_with_process_management(command, timeout=timeout)
            
            if result.returncode == 0:
                print("MATCHVAR execution completed successfully")
                try:
                    result_df = self._process_matchvar_results(output_prefix, protocols)
                    if result_df is not None:
                        print(f"成功处理MATCHVAR结果，包含 {len(result_df)} 行数据")
                    else:
                        print("MATCHVAR结果处理失败，返回None")
                    return result_df
                except Exception as e:
                    print(f"处理MATCHVAR结果时出错: {e}")
                    import traceback
                    traceback.print_exc()
                    return None
            else:
                print(f"MATCHVAR execution failed with return code: {result.returncode}")
                return None
                
        except TimeoutError as e:
            print(f"MATCHVAR execution timeout: {e}")
            return None
        except Exception as e:
            print(f"Exception occurred during MATCHVAR execution: {e}")
            return None
    
    def _cleanup_orig_files(self, output_prefix: str):
        """
        清理可能存在的.orig文件，避免重命名冲突
        
        Args:
            output_prefix: 输出文件前缀
        """
        try:
            # 查找所有可能的.orig文件
            for protocol in self.PROTOCOL_CONFIGS.keys():
                orig_file = os.path.join(self.input_files, f"{output_prefix}.{protocol}.exonic_variant_function.orig")
                if os.path.exists(orig_file):
                    print(f"清理已存在的.orig文件: {orig_file}")
                    os.remove(orig_file)
                    
            # 也清理其他可能的.orig文件
            for filename in os.listdir(self.input_files):
                if filename.endswith('.orig') and output_prefix in filename:
                    orig_file = os.path.join(self.input_files, filename)
                    print(f"清理已存在的.orig文件: {orig_file}")
                    os.remove(orig_file)
                    
        except Exception as e:
            print(f"清理.orig文件时出错: {e}")
            # 不抛出异常，继续执行
    
    def _execute_with_process_management(self, command: str, timeout: int = 600):
        """
        执行命令并进行进程管理
        
        Args:
            command: 要执行的命令
            timeout: 超时时间（秒），默认10分钟
            
        Returns:
            执行结果
        """
        process = None
        try:
            # 启动进程 - 直接输出到控制台，保持原来的输出格式
            process = subprocess.Popen(
                command,
                shell=True,
                stdout=None,  # 直接输出到控制台
                stderr=None,  # 直接输出到控制台
                text=True,
                encoding=get_system_encoding(),
                errors='replace'
            )
            
            print(f"MATCHVAR进程已启动，PID: {process.pid}")
            
            # 等待进程完成，带超时控制
            try:
                process.wait(timeout=timeout)
                print(f"MATCHVAR进程正常完成，返回码: {process.returncode}")
            except subprocess.TimeoutExpired:
                print(f"MATCHVAR进程超时（{timeout}秒），强制终止")
                process.terminate()
                try:
                    process.wait(timeout=10)  # 给进程10秒时间优雅退出
                except subprocess.TimeoutExpired:
                    print(f"强制杀死MATCHVAR进程: PID {process.pid}")
                    process.kill()
                    process.wait()
                raise TimeoutError(f"MATCHVAR执行超时（{timeout}秒）")
            
            # 创建结果对象 - 由于输出直接到控制台，这里创建空的结果
            result = subprocess.CompletedProcess(
                args=command,
                returncode=process.returncode,
                stdout="",  # 输出已经直接显示在控制台
                stderr=""   # 输出已经直接显示在控制台
            )
            
            print(f"MATCHVAR进程已完成，返回码: {result.returncode}")
            return result
            
        except Exception as e:
            print(f"执行MATCHVAR命令时出错: {e}")
            # 创建错误结果
            return subprocess.CompletedProcess(
                args=command,
                returncode=-1,
                stdout="",
                stderr=str(e)
            )
        finally:
            # 确保进程被清理
            if process and process.poll() is None:
                print(f"强制终止MATCHVAR进程: PID {process.pid}")
                try:
                    process.terminate()
                    process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    print(f"强制杀死MATCHVAR进程: PID {process.pid}")
                    process.kill()
                except Exception as e:
                    print(f"清理进程时出错: {e}")
            
            # 清理可能的子进程
            self._cleanup_child_processes()
    
    def _cleanup_child_processes(self):
        """
        清理可能的子进程
        """
        try:
            current_pid = os.getpid()
            current_process = psutil.Process(current_pid)
            
            # 查找当前进程的子进程
            children = current_process.children(recursive=True)
            
            for child in children:
                try:
                    # 检查是否是MATCHVAR相关的进程
                    if any(keyword in child.name().lower() for keyword in ['python', 'matchvar', 'perl']):
                        print(f"清理子进程: PID {child.pid}, 名称: {child.name()}")
                        child.terminate()
                        
                        # 等待进程终止
                        try:
                            child.wait(timeout=5)
                        except psutil.TimeoutExpired:
                            print(f"强制杀死子进程: PID {child.pid}")
                            child.kill()
                            
                except (psutil.NoSuchProcess, psutil.AccessDenied) as e:
                    print(f"无法清理子进程: {e}")
                    
        except Exception as e:
            print(f"清理子进程时出错: {e}")
    
    def _process_matchvar_results(self, output_prefix: str, protocols: List[str]) -> Optional[pd.DataFrame]:
        """
        处理MATCHVAR执行结果，包括自定义数据库的后处理
        
        Args:
            output_prefix: 输出文件前缀
            protocols: 使用的协议列表
            
        Returns:
            注释结果DataFrame
        """
        try:
            # 查找输出文件
            output_files = []
            
            # 根据基因组版本查找TSV格式的输出文件
            tsv_file = f"{output_prefix}.{self.genome_version}_multianno.tsv"
            tsv_path = os.path.join(self.input_files, tsv_file)
            if os.path.exists(tsv_path):
                output_files.append(tsv_path)
            
            # 如果TSV文件不存在，查找当前目录
            if not output_files:
                tsv_path = tsv_file
                if os.path.exists(tsv_path):
                    output_files.append(tsv_path)
            
            if output_files:
                output_path = output_files[0]
                file_size = os.path.getsize(output_path)
                print(f"找到MATCHVAR输出文件: {output_path}, 文件大小: {file_size} bytes")
                
                # 读取结果文件
                try:
                    # 使用更严格的TSV读取方式
                    result_df = pd.read_csv(output_path, sep='\t', low_memory=False,
                                          na_values=[], keep_default_na=False,
                                          encoding='utf-8', quoting=3)  # quoting=3 means no quoting
                    
                    # 验证基本列是否存在
                    required_columns = ['Chr', 'Start', 'End', 'Ref', 'Alt']
                    missing_columns = [col for col in required_columns if col not in result_df.columns]
                    
                    if missing_columns:
                        print(f"警告：TSV文件缺少必需的列: {missing_columns}")
                        print(f"实际列名: {list(result_df.columns)[:10]}")
                        return None
                    
                    # 验证数据完整性：检查前几行数据
                    for idx in range(min(3, len(result_df))):
                        row = result_df.iloc[idx]
                        chr_val = str(row.get('Chr', ''))
                        alt_val = str(row.get('Alt', ''))
                        
                        # 检查Alt字段是否包含功能注释而不是碱基
                        if alt_val in ['exonic', 'intronic', 'intergenic', 'utr5', 'utr3', 'splicing']:
                            print(f"检测到列错位问题：第{idx+1}行 Alt={alt_val}，这不是碱基序列")
                            print(f"行数据: Chr={chr_val}, Alt={alt_val}")
                            # 尝试重新解析文件
                            print("尝试使用不同的解析方式...")
                            try:
                                result_df = pd.read_csv(output_path, sep='\t', low_memory=False,
                                                      na_values=[], keep_default_na=False,
                                                      encoding='utf-8-sig', engine='python')
                                # 重新验证
                                test_row = result_df.iloc[0]
                                test_alt = str(test_row.get('Alt', ''))
                                if test_alt in ['exonic', 'intronic', 'intergenic', 'utr5', 'utr3', 'splicing']:
                                    print("重新解析仍有问题，返回None")
                                    return None
                                else:
                                    print("重新解析成功")
                                    break
                            except Exception as e:
                                print(f"重新解析失败: {e}")
                                return None
                    
                    # 将所有的NaN值替换为空字符串
                    result_df = result_df.fillna('')
                    print(f"成功读取TSV文件，包含 {len(result_df)} 行数据")
                    
                    # 显示前几行的关键列以供调试
                    if len(result_df) > 0:
                        sample_data = result_df[['Chr', 'Start', 'End', 'Ref', 'Alt']].head(3)
                        print("前3行基本数据:")
                        for idx, row in sample_data.iterrows():
                            print(f"  行{idx+1}: Chr={row['Chr']}, Start={row['Start']}, End={row['End']}, Ref={row['Ref']}, Alt={row['Alt']}")
                    
                    # 处理自定义数据库注释
                    try:
                        result_df = self._add_custom_database_annotations(result_df, protocols)
                        print("自定义数据库注释处理完成")
                    except Exception as e:
                        print(f"处理自定义数据库注释时出错: {e}")
                        import traceback
                        traceback.print_exc()
                        # 即使自定义数据库处理失败，也返回原始结果
                        print("继续使用原始MATCHVAR结果")
                    
                    return result_df
                except Exception as e:
                    print(f"读取TSV文件失败: {e}")
                    return None
            else:
                print(f"未找到MATCHVAR输出文件: {tsv_file}")
                return None
                
        except Exception as e:
            print(f"处理MATCHVAR结果时出错: {e}")
            return None
    
    def _add_custom_database_annotations(self, result_df: pd.DataFrame, protocols: List[str]) -> pd.DataFrame:
        """
        为结果DataFrame添加自定义数据库注释
        
        Args:
            result_df: MATCHVAR结果DataFrame
            protocols: 使用的协议列表
            
        Returns:
            添加了自定义数据库注释的DataFrame
        """
        try:
            enhanced_configs = self.get_enhanced_protocol_configs()
            
            # 查找自定义数据库协议
            custom_db_protocols = []
            for protocol in protocols:
                if enhanced_configs[protocol].get('is_custom'):
                    custom_db_protocols.append(protocol)
            
            if not custom_db_protocols:
                return result_df
            
            print(f"开始处理自定义数据库注释: {custom_db_protocols}")
            
            # 为每个自定义数据库添加注释
            for custom_db_name in custom_db_protocols:
                custom_db_config = enhanced_configs[custom_db_name]
                custom_db_path = custom_db_config['file_path']
                
                # 读取自定义数据库文件
                custom_annotations = self._load_custom_database(custom_db_path)
                
                if custom_annotations:
                    # 为结果DataFrame添加自定义数据库列
                    result_df = self._merge_custom_annotations(result_df, custom_annotations, custom_db_name)
                    print(f"成功添加自定义数据库 {custom_db_name} 的注释")
                else:
                    print(f"无法加载自定义数据库 {custom_db_name}")
            
            return result_df
            
        except Exception as e:
            print(f"添加自定义数据库注释时出错: {e}")
            return result_df
    
    def _load_custom_database(self, db_path: str) -> Dict[str, Dict]:
        """
        加载自定义数据库文件
        
        Args:
            db_path: 自定义数据库文件路径
            
        Returns:
            自定义数据库注释字典
        """
        try:
            annotations = {}
            
            with open(db_path, 'r', encoding='utf-8') as f:
                lines = f.readlines()
                
                # 跳过注释行
                data_lines = [line.strip() for line in lines if not line.strip().startswith('#')]
                
                for line in data_lines:
                    if not line:
                        continue
                    
                    parts = line.split('\t')
                    if len(parts) >= 5:
                        # 创建唯一键：chr:start:end:ref:alt
                        chr_pos = parts[0]
                        start = parts[1]
                        end = parts[2]
                        ref = parts[3]
                        alt = parts[4]
                        
                        key = f"{chr_pos}:{start}:{end}:{ref}:{alt}"
                        
                        # 存储所有列的数据
                        annotation_data = {
                            'chr': chr_pos,
                            'start': start,
                            'end': end,
                            'ref': ref,
                            'alt': alt
                        }
                        
                        # 添加额外的列（如果有的话）
                        for i, value in enumerate(parts[5:], 5):
                            annotation_data[f'col_{i}'] = value
                        
                        annotations[key] = annotation_data
            
            print(f"加载自定义数据库，包含 {len(annotations)} 条记录")
            return annotations
            
        except Exception as e:
            print(f"加载自定义数据库失败: {e}")
            return {}
    
    def _merge_custom_annotations(self, result_df: pd.DataFrame, custom_annotations: Dict[str, Dict], db_name: str) -> pd.DataFrame:
        """
        将自定义数据库注释合并到结果DataFrame中
        
        Args:
            result_df: 结果DataFrame
            custom_annotations: 自定义数据库注释
            db_name: 数据库名称
            
        Returns:
            合并后的DataFrame
        """
        try:
            # 创建新列来存储自定义数据库的注释
            result_df[f'{db_name}_found'] = ''
            result_df[f'{db_name}_data'] = ''
            
            # 为每一行查找匹配的自定义数据库记录
            for idx, row in result_df.iterrows():
                # 创建键来匹配自定义数据库
                chr_pos = str(row.get('Chr', ''))
                start = str(row.get('Start', ''))
                end = str(row.get('End', ''))
                ref = str(row.get('Ref', ''))
                alt = str(row.get('Alt', ''))
                
                key = f"{chr_pos}:{start}:{end}:{ref}:{alt}"
                
                if key in custom_annotations:
                    result_df.at[idx, f'{db_name}_found'] = 'Yes'
                    # 将自定义数据库的数据转换为字符串
                    custom_data = custom_annotations[key]
                    data_str = ';'.join([f"{k}={v}" for k, v in custom_data.items() if k not in ['chr', 'start', 'end', 'ref', 'alt']])
                    result_df.at[idx, f'{db_name}_data'] = data_str
                else:
                    result_df.at[idx, f'{db_name}_found'] = 'No'
            
            return result_df
            
        except Exception as e:
            print(f"合并自定义数据库注释时出错: {e}")
            return result_df
    
    def _convert_vcf_to_matchvar_input(self, vcf_file: str) -> Optional[str]:
        """
        Convert VCF file to MATCHVAR input format
        
        Args:
            vcf_file: VCF file name
            
        Returns:
            Converted file name, return None if failed
        """
        try:
            # Get the path of the convert2matchvar.py script
            convert2matchvar_script = os.path.join(self.matchvar_path, 'convert2matchvar.py')
            
            if not os.path.exists(convert2matchvar_script):
                print(f"convert2matchvar.py does not exist: {convert2matchvar_script}")
                return None
            
            # Generate the output file name
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            output_file = f"converted_{timestamp}.mvinput"
            output_path = os.path.join(self.input_files, output_file)
            
            # Build the conversion command
            input_path = os.path.join(self.input_files, vcf_file)
            command = f"{PYTHON_EXECUTABLE} {convert2matchvar_script} -includeinfo -allsample -withfreq -format vcf4 {input_path} > {output_path}"
            
            print(f"Executing VCF conversion command: {command}")
            
            # Execute the conversion, using more secure encoding processing
            result = subprocess.run(
                command, 
                shell=True, 
                capture_output=True, 
                text=True, 
                encoding='utf-8',
                errors='replace'  # Use replace mode to handle encoding errors
            )
            
            if result.returncode == 0 and os.path.exists(output_path):
                file_size = os.path.getsize(output_path)
                print(f"VCF conversion completed, output file: {output_path}, file size: {file_size} bytes")
                
                # Verify the conversion result
                with open(output_path, 'r', encoding='utf-8', errors='replace') as f:
                    content = f.read()
                    lines = content.strip().split('\n')
                    print(f"The converted file contains {len(lines)} rows of data")
                    if lines:
                        print(f"First row example: {lines[0]}")
                
                return output_file
            else:
                print(f"VCF conversion failed, return code: {result.returncode}")
                if result.stderr:
                    print(f"Error information: {result.stderr}")
                return None
                
        except Exception as e:
            print(f"Exception occurred during VCF conversion: {e}")
            return None
    
    def run_matchvar_with_custom_protocols(self, input_file: str, 
                                        gene_info: List[str] = None,
                                        databases: List[str] = None,
                                        predictions: List[str] = None,
                                        buildver: str = 'hg19',
                                        output_prefix: str = 'matchvar_result') -> Optional[pd.DataFrame]:
        """
        Run MATCHVAR using custom categorized protocols
        
        Args:
            input_file: input file name
            gene_info: gene information protocols
            databases: database protocols  
            predictions: prediction software protocols
            buildver: genome version
            output_prefix: output file prefix
            
        Returns:
            Annotation result DataFrame
        """
        # Build the protocol list
        protocols = []
        
        # Add gene information
        if gene_info:
            protocols.extend(gene_info)
        
        # Add database
        if databases:
            protocols.extend(databases)
        
        # Add prediction software
        if predictions:
            protocols.extend(predictions)
        
        # If no protocol is specified, use default configuration
        if not protocols:
            protocols = ['refGene', 'cytoBand', 'exac03', 'avsift', 'dbnsfp42a', 
                       'dbscsnv11', 'gnomad211_genome', 'esp6500siv2_all', 
                       'revel', 'cadd13gt10', 'AlphaMissense']
        
        return self.run_matchvar(input_file, protocols, buildver, output_prefix)
    
    def get_protocol_categories(self) -> Dict[str, List[str]]:
        """Get protocols grouped by category"""
        categories = {}
        for protocol, config in self.PROTOCOL_CONFIGS.items():
            category = config['category']
            if category not in categories:
                categories[category] = []
            categories[category].append(protocol)
        return categories
    
    def add_custom_protocol(self, name: str, operation: str, description: str, 
                           category: str) -> bool:
        """
        Add custom protocol configuration
        
        Args:
            name: protocol name
            operation: operation type (g/r/f)
            description: description
            category: category
            
        Returns:
            Whether the addition is successful
        """
        if name in self.PROTOCOL_CONFIGS:
            print(f"Protocol {name} already exists")
            return False
        
        self.PROTOCOL_CONFIGS[name] = {
            'operation': operation,
            'description': description,
            'category': category
        }
        
        print(f"Successfully added custom protocol: {name}")
        return True
    
    def remove_protocol(self, name: str) -> bool:
        """
        Remove protocol configuration
        
        Args:
            name: protocol name
            
        Returns:
            Whether the removal is successful
        """
        if name in self.PROTOCOL_CONFIGS:
            del self.PROTOCOL_CONFIGS[name]
            print(f"Successfully removed protocol: {name}")
            return True
        else:
            print(f"Protocol {name} does not exist")
            return False