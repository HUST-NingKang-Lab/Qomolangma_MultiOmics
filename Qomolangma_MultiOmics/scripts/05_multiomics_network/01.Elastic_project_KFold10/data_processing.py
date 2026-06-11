import os
import pandas as pd
import numpy as np
from pathlib import Path
import logging
from typing import Tuple, List
from sklearn.model_selection import train_test_split
from utils import create_directories

def validate_data_quality(df: pd.DataFrame) -> Tuple[bool, str]:
    """验证数据质量"""
    # 检查数据形状
    if df.empty:
        return False, "数据为空"
    
    # 检查缺失值比例
    missing_ratio = df.isnull().sum().sum() / (df.shape[0] * df.shape[1])
    if missing_ratio > 0.5:
        return False, f"缺失值比例过高: {missing_ratio:.2%}"
    
    # 检查数值数据
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    if len(numeric_cols) == 0:
        return False, "没有数值型数据"
    
    # 检查无穷值
    if np.any(np.isinf(df[numeric_cols].values)):
        return False, "数据包含无穷值"
    
    return True, "数据质量良好"

def split_train_test(data_path: str, output_dir: str, test_size: float = 0.2, 
                    random_state: int = 42) -> Tuple[str, str, pd.DataFrame, pd.DataFrame]:
    """划分训练集和测试集"""
    logging.info(f"正在读取数据文件: {data_path}")
    
    try:
        # 读取数据，不指定index_col，让pandas自动生成整数索引
        df = pd.read_csv(data_path)
        logging.info(f"数据形状: {df.shape}")
        
        # 验证数据质量
        is_valid, message = validate_data_quality(df)
        if not is_valid:
            raise ValueError(f"数据质量检查失败: {message}")
        
        # 随机划分训练集和测试集
        train_indices, test_indices = train_test_split(
            df.index,  # 使用自动生成的整数索引
            test_size=test_size, 
            random_state=random_state
        )
        
        train_df = df.loc[train_indices].copy()
        test_df = df.loc[test_indices].copy()
        
        logging.info(f"训练集样本数: {len(train_df)}")
        logging.info(f"测试集样本数: {len(test_df)}")
        
        # 根据输入文件名生成输出文件名
        input_filename = Path(data_path).stem
        train_filename = f"{input_filename}_train.csv"
        test_filename = f"{input_filename}_test.csv"
        
        # 创建子目录
        sub_output_dir = os.path.join(output_dir, input_filename)
        create_directories([sub_output_dir])
        
        # 保存训练集和测试集
        train_path = os.path.join(sub_output_dir, train_filename)
        test_path = os.path.join(sub_output_dir, test_filename)
        
        # 保存时保留原始索引（如果有需要），或者重置索引
        train_df.to_csv(train_path, index=False)  # 不保存索引
        test_df.to_csv(test_path, index=False)    # 不保存索引
        
        logging.info(f"训练集已保存至: {train_path}")
        logging.info(f"测试集已保存至: {test_path}")
        
        return train_path, test_path, train_df, test_df
        
    except Exception as e:
        logging.error(f"数据划分失败: {e}")
        raise

def load_metadata(metadata_path: str, input_filename: str) -> Tuple[List[str], List[str], pd.DataFrame]:
    """加载metadata文件并根据文件名提取特征"""
    logging.info(f"正在读取metadata文件: {metadata_path}")
    
    try:
        metadata = pd.read_csv(metadata_path)
        
        # 验证metadata必需列
        required_cols = ['Source', 'Group', 'Feature']
        missing_cols = [col for col in required_cols if col not in metadata.columns]
        if missing_cols:
            raise ValueError(f"metadata文件缺少必需列: {missing_cols}")
        
        y_source = None
        x_source = None
        
        # 检测因变量类型 - 改进正则匹配
        if "metabolin" in input_filename.lower():
            y_source = "Metabolin_lg"
        elif "lipid" in input_filename.lower():
            y_source = "Lipid_lg"
        
        # 检测自变量类型 - 改进正则匹配
        if "gut" in input_filename.lower():
            x_source = "gut"
        elif "oral" in input_filename.lower():
            x_source = "oral"
        elif "skin" in input_filename.lower():
            x_source = "skin"
        
        # 验证提取的信息
        if y_source is None or x_source is None:
            raise ValueError(f"无法从文件名 '{input_filename}' 中提取有效的自变量或因变量信息")
        
        # 提取自变量和因变量特征名
        X_features = metadata[
            ((metadata['Source'] == x_source) | (metadata['Source'] == 'people')) & 
            (metadata['Group'] == 'X')
        ]['Feature'].tolist()
        
        y_features = metadata[
            (metadata['Source'] == y_source) & 
            (metadata['Group'] == 'Y')
        ]['Feature'].tolist()
        
        logging.info(f"自变量来源: {x_source}, 特征数: {len(X_features)}")
        logging.info(f"因变量来源: {y_source}, 特征数: {len(y_features)}")

        return X_features, y_features, metadata
        
    except Exception as e:
        logging.error(f"metadata加载失败: {e}")
        raise