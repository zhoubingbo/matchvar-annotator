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
        
        # Frequency databases
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
        
        # Functional prediction software
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
        self.resources_dir = resources_dir or os.path.join(os.path.dirname(__file__), '..', 'resources')
        self.path_humandb = os.path.join(self.resources_dir, 'humandb')
        
        # Use the incoming input file directory, if not provided, use the default value
        self.input_files = input_files_dir or os.path.join(self.resources_dir, 'input_files')
        
        # MATCHVAR file is in the utils/matchvar directory, not in the resources directory
        self.matchvar_path = os.path.dirname(__file__)  # The current file is in the utils/matchvar directory
        self.python_table_matchvar = os.path.join(self.matchvar_path, 'table_matchvar.py')
        
        # Set genome version and thread count
        self.genome_version = genome_version
        self.thread_count = int(thread_count) if thread_count else 4
        print(f"Initialize MATCHVAR annotator, genome version: {self.genome_version}, thread count: {self.thread_count}")
        
        # Ensure directory exists
        os.makedirs(self.path_humandb, exist_ok=True)
        os.makedirs(self.input_files, exist_ok=True)
        
        # Check custom database directory
        self.custom_db_dir = os.path.join(self.resources_dir, 'custom_databases')
        os.makedirs(self.custom_db_dir, exist_ok=True)

    def get_custom_databases(self) -> Dict[str, Dict]:
        """Get custom databases configuration based on genome version"""
        custom_databases = {}
        
        # Check custom database directory in humandb directory
        if os.path.exists(self.path_humandb):
            for filename in os.listdir(self.path_humandb):
                # Select database file according to genome version
                expected_prefix = f"{self.genome_version}_"
                if filename.startswith(expected_prefix) and filename.endswith('.txt'):
                    # Check if it is a custom database (not in standard protocol configuration)
                    db_name = filename[len(expected_prefix):-4]  # Remove genome version prefix and '.txt'
                    if db_name not in self.PROTOCOL_CONFIGS:
                        file_path = os.path.join(self.path_humandb, filename)
                        
                        # Verify file format is correct (at least contains Chr, Start, End, Ref, Alt columns)
                        try:
                            with open(file_path, 'r', encoding='utf-8') as f:
                                first_line = f.readline().strip()
                                if first_line.startswith('#'):
                                    # Skip comment line, read second line
                                    second_line = f.readline().strip()
                                    if second_line:
                                        # Check if there are enough columns (at least 5 columns: Chr, Start, End, Ref, Alt)
                                        cols = second_line.split('\t')
                                        if len(cols) >= 5:
                                            custom_databases[db_name] = {
                                                'operation': 'f',  # Default use frequency operation
                                                'description': f'Custom Database: {db_name} ({self.genome_version})',
                                                'category': 'custom_database',
                                                'is_custom': True,
                                                'filename': filename,
                                                'file_path': file_path
                                            }
                        except Exception as e:
                            print(f"Error validating custom database file {filename}: {e}")
                            continue
        
        print(f"Found {len(custom_databases)} custom database files (genome version: {self.genome_version})")
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
            raise ValueError(f"Protocol validation failed: {'; '.join(errors)}")
        
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
        arguments = ['' for _ in protocols]  # Use empty string as default parameter
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
        
        # Add argument parameter (if not empty)
        if argument_str and argument_str != '.':
            cmd_parts.extend(['-argument', argument_str])
        
        # Add MANE transcript filtering if requested
        if additional_args and additional_args.get('use_mane_transcript'):
            mane_file = os.path.join(humandb_abs, 'mane_transcript.txt')
            if os.path.exists(mane_file):
                cmd_parts.extend(['-mane_file', mane_file])
                print(f"Add MANE transcript filtering: {mane_file}")
            else:
                print(f"Warning: MANE transcript file does not exist: {mane_file}")
        
        # Add otherinfo parameter for complete annotation information
        cmd_parts.append('-otherinfo')
        print("Add -otherinfo parameter to get complete annotation information")
        
        # Add thread count parameter
        cmd_parts.extend(['-thread', str(self.thread_count)])
        print(f"Add thread count parameter: {self.thread_count}")
        
        # Check if there are custom databases, if there are, we need to handle them specially
        custom_db_protocols = []
        for protocol in protocols:
            if enhanced_configs[protocol].get('is_custom'):
                custom_db_protocols.append(protocol)
        
        if custom_db_protocols:
            print(f"Found custom database protocol: {custom_db_protocols}")
            # For custom databases, we need to ensure they are handled correctly
            # Currently MATCHVAR may not support using custom databases as standard protocols
            # We may need to use other methods, such as converting custom databases to standard format
            print("Warning: Custom databases may not be handled correctly because of MATCHVAR's limitations")
        
        # Add additional arguments if provided
        if additional_args:
            for key, value in additional_args.items():
                if isinstance(value, bool):
                    if value:
                        cmd_parts.extend([f'-{key}'])
                else:
                    cmd_parts.extend([f'-{key}', str(value)])
        
        return ' '.join(cmd_parts)
    
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
            print("Enable MANE transcript filtering mode")
        else:
            print("Use all transcript annotation mode")
        
        try:
            # Check the input file type
            input_path = os.path.join(self.input_files, input_file)
            is_vcf = input_file.lower().endswith('.vcf')
            
            # Clean up possible .orig files to avoid renaming conflicts
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
                        print(f"Successfully processed MATCHVAR results, containing {len(result_df)} rows")
                    else:
                        print("MATCHVAR results processing failed, returning None")
                    return result_df
                except Exception as e:
                    print(f"Error processing MATCHVAR results: {e}")
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
        Clean up possible .orig files to avoid renaming conflicts
        
        Args:
            output_prefix: output file prefix
        """
        try:
            # Find all possible .orig files
            for protocol in self.PROTOCOL_CONFIGS.keys():
                orig_file = os.path.join(self.input_files, f"{output_prefix}.{protocol}.exonic_variant_function.orig")
                if os.path.exists(orig_file):
                    print(f"Clean up existing .orig file: {orig_file}")
                    os.remove(orig_file)
                    
            # Also clean up other possible .orig files
            for filename in os.listdir(self.input_files):
                if filename.endswith('.orig') and output_prefix in filename:
                    orig_file = os.path.join(self.input_files, filename)
                    print(f"Clean up existing .orig file: {orig_file}")
                    os.remove(orig_file)
                    
        except Exception as e:
            print(f"Error cleaning .orig file: {e}")
            # Do not throw exception, continue execution
    
    def _execute_with_process_management(self, command: str, timeout: int = 600):
        """
        Execute command and manage processes
        
        Args:
            command: command to execute
            timeout: timeout time (seconds), default 10 minutes
            
        Returns:
            Execution result
        """
        process = None
        try:
            # Start process - directly output to console, keep original output format
            process = subprocess.Popen(
                command,
                shell=True,
                stdout=None,  # Directly output to console
                stderr=None,  # Directly output to console
                text=True,
                encoding=get_system_encoding(),
                errors='replace'
            )
            
            print(f"MATCHVAR process started, PID: {process.pid}")
            
            # Wait for process to complete, with timeout control
            try:
                process.wait(timeout=timeout)
                print(f"MATCHVAR process completed normally, return code: {process.returncode}")
            except subprocess.TimeoutExpired:
                print(f"MATCHVAR process timed out ({timeout} seconds), force termination")
                process.terminate()
                try:
                    process.wait(timeout=10)  # Give process 10 seconds to gracefully exit
                except subprocess.TimeoutExpired:
                    print(f"Force kill MATCHVAR process: PID {process.pid}")
                    process.kill()
                    process.wait()
                raise TimeoutError(f"MATCHVAR execution timed out ({timeout} seconds)")
            
            # Create result object - since output is directly to console, create empty result here
            result = subprocess.CompletedProcess(
                args=command,
                returncode=process.returncode,
                stdout="",  # Output is already directly displayed on console
                stderr=""   # Output is already directly displayed on console
            )
            
            print(f"MATCHVAR process completed, return code: {result.returncode}")
            return result
            
        except Exception as e:
            print(f"Error executing MATCHVAR command: {e}")
            # Create error result
            return subprocess.CompletedProcess(
                args=command,
                returncode=-1,
                stdout="",
                stderr=str(e)
            )
        finally:
            # Ensure process is cleaned up
            if process and process.poll() is None:
                print(f"Force terminate MATCHVAR process: PID {process.pid}")
                try:
                    process.terminate()
                    process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    print(f"Force kill MATCHVAR process: PID {process.pid}")
                    process.kill()
                except Exception as e:
                    print(f"Error cleaning process: {e}")
            
            # Clean up possible child processes
            self._cleanup_child_processes()
    
    def _cleanup_child_processes(self):
        """
        Clean up possible child processes
        """
        try:
            current_pid = os.getpid()
            current_process = psutil.Process(current_pid)
            
            # Find current process's children
            children = current_process.children(recursive=True)
            
            for child in children:
                try:
                    # Check if it is a MATCHVAR related process
                    if any(keyword in child.name().lower() for keyword in ['python', 'matchvar', 'perl']):
                        print(f"Clean up child process: PID {child.pid}, name: {child.name()}")
                        child.terminate()
                        
                        # Wait for process to terminate
                        try:
                            child.wait(timeout=5)
                        except psutil.TimeoutExpired:
                            print(f"Force kill child process: PID {child.pid}")
                            child.kill()
                            
                except (psutil.NoSuchProcess, psutil.AccessDenied) as e:
                    print(f"Cannot clean up child process: {e}")
                    
        except Exception as e:
            print(f"Error cleaning child process: {e}")
    
    def _process_matchvar_results(self, output_prefix: str, protocols: List[str]) -> Optional[pd.DataFrame]:
        """
        Process MATCHVAR execution results, including custom database post-processing
        
        Args:
            output_prefix: output file prefix
            protocols: used protocol list
            
        Returns:
            Annotation result DataFrame
        """
        try:
            # Find output files
            output_files = []
            
            # Find TSV format output files according to genome version
            tsv_file = f"{output_prefix}.{self.genome_version}_multianno.tsv"
            tsv_path = os.path.join(self.input_files, tsv_file)
            if os.path.exists(tsv_path):
                output_files.append(tsv_path)
            
            # If TSV file does not exist, find current directory
            if not output_files:
                tsv_path = tsv_file
                if os.path.exists(tsv_path):
                    output_files.append(tsv_path)
            
            if output_files:
                output_path = output_files[0]
                file_size = os.path.getsize(output_path)
                print(f"Find MATCHVAR output file: {output_path}, file size: {file_size} bytes")
                
                # Read result file
                try:
                    # Use more strict TSV reading method
                    result_df = pd.read_csv(output_path, sep='\t', low_memory=False,
                                          na_values=[], keep_default_na=False,
                                          encoding='utf-8', quoting=3)  # quoting=3 means no quoting
                    
                    # Verify if basic columns exist
                    required_columns = ['Chr', 'Start', 'End', 'Ref', 'Alt']
                    missing_columns = [col for col in required_columns if col not in result_df.columns]
                    
                    if missing_columns:
                        print(f"Warning: TSV file missing required columns: {missing_columns}")
                        print(f"Actual column names: {list(result_df.columns)[:10]}")
                        return None
                    
                    # Verify data completeness: check first few rows
                    for idx in range(min(3, len(result_df))):
                        row = result_df.iloc[idx]
                        chr_val = str(row.get('Chr', ''))
                        alt_val = str(row.get('Alt', ''))
                        
                        # Check if Alt field contains functional annotation instead of base sequence
                        if alt_val in ['exonic', 'intronic', 'intergenic', 'utr5', 'utr3', 'splicing']:
                            print(f"Detected column misalignment problem: Row {idx+1} Alt={alt_val}, not base sequence")
                            print(f"Row data: Chr={chr_val}, Alt={alt_val}")
                            # Try to reparse file
                            print("Trying different parsing method...")
                            try:
                                result_df = pd.read_csv(output_path, sep='\t', low_memory=False,
                                                      na_values=[], keep_default_na=False,
                                                      encoding='utf-8-sig', engine='python')
                                # Re-verify
                                test_row = result_df.iloc[0]
                                test_alt = str(test_row.get('Alt', ''))
                                if test_alt in ['exonic', 'intronic', 'intergenic', 'utr5', 'utr3', 'splicing']:
                                    print("Re-parsing still has problems, returning None")
                                    return None
                                else:
                                    print("Re-parsing successful")
                                    break
                            except Exception as e:
                                print(f"Re-parsing failed: {e}")
                                return None
                    
                    # Replace all NaN values with empty string
                    result_df = result_df.fillna('')
                    print(f"Successfully read TSV file, containing {len(result_df)} rows of data")
                    
                    # Display first few rows of key columns for debugging
                    if len(result_df) > 0:
                        sample_data = result_df[['Chr', 'Start', 'End', 'Ref', 'Alt']].head(3)
                        print("First 3 rows of basic data:")
                        for idx, row in sample_data.iterrows():
                            print(f"  Row {idx+1}: Chr={row['Chr']}, Start={row['Start']}, End={row['End']}, Ref={row['Ref']}, Alt={row['Alt']}")
                    
                    # Process custom database annotations
                    try:
                        result_df = self._add_custom_database_annotations(result_df, protocols)
                        print("Custom database annotations processed successfully")
                    except Exception as e:
                        print(f"Error processing custom database annotations: {e}")
                        import traceback
                        traceback.print_exc()
                        # Even if custom database processing fails, return the original result
                        print("Continue using original MATCHVAR results")
                    
                    return result_df
                except Exception as e:
                    print(f"Failed to read TSV file: {e}")
                    return None
            else:
                print(f"MATCHVAR output file not found: {tsv_file}")
                return None
                
        except Exception as e:
            print(f"Error processing MATCHVAR results: {e}")
            return None
    
    def _add_custom_database_annotations(self, result_df: pd.DataFrame, protocols: List[str]) -> pd.DataFrame:
        """
        Add custom database annotations to the result DataFrame
        
        Args:
            result_df: MATCHVAR result DataFrame
            protocols: used protocol list
            
        Returns:
            DataFrame with custom database annotations
        """
        try:
            enhanced_configs = self.get_enhanced_protocol_configs()
            
            # Find custom database protocols
            custom_db_protocols = []
            for protocol in protocols:
                if enhanced_configs[protocol].get('is_custom'):
                    custom_db_protocols.append(protocol)
            
            if not custom_db_protocols:
                return result_df
            
            print(f"Start processing custom database annotations: {custom_db_protocols}")
            
            # Add annotations for each custom database
            for custom_db_name in custom_db_protocols:
                custom_db_config = enhanced_configs[custom_db_name]
                custom_db_path = custom_db_config['file_path']
                
                # Read custom database file
                custom_annotations = self._load_custom_database(custom_db_path)
                
                if custom_annotations:
                    # Add custom database columns to the result DataFrame
                    result_df = self._merge_custom_annotations(result_df, custom_annotations, custom_db_name)
                    print(f"Successfully added custom database {custom_db_name} annotations")
                else:
                    print(f"Cannot load custom database {custom_db_name}")
            
            return result_df
            
        except Exception as e:
            print(f"Error adding custom database annotations: {e}")
            return result_df
    
    def _load_custom_database(self, db_path: str) -> Dict[str, Dict]:
        """
        Load custom database file
        
        Args:
            db_path: custom database file path
            
        Returns:
            Custom database annotations dictionary
        """
        try:
            annotations = {}
            
            with open(db_path, 'r', encoding='utf-8') as f:
                lines = f.readlines()
                
                # Skip comment lines
                data_lines = [line.strip() for line in lines if not line.strip().startswith('#')]
                
                for line in data_lines:
                    if not line:
                        continue
                    
                    parts = line.split('\t')
                    if len(parts) >= 5:
                        # Create unique key: chr:start:end:ref:alt
                        chr_pos = parts[0]
                        start = parts[1]
                        end = parts[2]
                        ref = parts[3]
                        alt = parts[4]
                        
                        key = f"{chr_pos}:{start}:{end}:{ref}:{alt}"
                        
                        # Store data for all columns
                        annotation_data = {
                            'chr': chr_pos,
                            'start': start,
                            'end': end,
                            'ref': ref,
                            'alt': alt
                        }
                        
                        # Add additional columns (if any)
                        for i, value in enumerate(parts[5:], 5):
                            annotation_data[f'col_{i}'] = value
                        
                        annotations[key] = annotation_data
            
            print(f"Loaded custom database, containing {len(annotations)} records")
            return annotations
            
        except Exception as e:
            print(f"Failed to load custom database: {e}")
            return {}
    
    def _merge_custom_annotations(self, result_df: pd.DataFrame, custom_annotations: Dict[str, Dict], db_name: str) -> pd.DataFrame:
        """
        Merge custom database annotations into the result DataFrame
        
        Args:
            result_df: result DataFrame
            custom_annotations: custom database annotations
            db_name: database name
            
        Returns:
            Merged DataFrame
        """
        try:
            # Create new columns to store custom database annotations
            result_df[f'{db_name}_found'] = ''
            result_df[f'{db_name}_data'] = ''
            
            # Find matching custom database records for each row
            for idx, row in result_df.iterrows():
                # Create key to match custom database
                chr_pos = str(row.get('Chr', ''))
                start = str(row.get('Start', ''))
                end = str(row.get('End', ''))
                ref = str(row.get('Ref', ''))
                alt = str(row.get('Alt', ''))
                
                key = f"{chr_pos}:{start}:{end}:{ref}:{alt}"
                
                if key in custom_annotations:
                    result_df.at[idx, f'{db_name}_found'] = 'Yes'
                    # Convert custom database data to string
                    custom_data = custom_annotations[key]
                    data_str = ';'.join([f"{k}={v}" for k, v in custom_data.items() if k not in ['chr', 'start', 'end', 'ref', 'alt']])
                    result_df.at[idx, f'{db_name}_data'] = data_str
                else:
                    result_df.at[idx, f'{db_name}_found'] = 'No'
            
            return result_df
            
        except Exception as e:
            print(f"Error merging custom database annotations: {e}")
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