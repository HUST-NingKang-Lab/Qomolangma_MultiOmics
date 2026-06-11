import os
import sys
import pickle
import logging
from pathlib import Path
from typing import List, Dict, Optional

def setup_logging(log_file: str = None) -> None:
    """设置日志配置"""
    log_format = '%(asctime)s - %(levelname)s - %(message)s'
    if log_file:
        logging.basicConfig(
            level=logging.INFO,
            format=log_format,
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )
    else:
        logging.basicConfig(level=logging.INFO, format=log_format)

def create_directories(paths: List[str]) -> None:
    """创建必要的目录"""
    for path in paths:
        Path(path).mkdir(parents=True, exist_ok=True)
        logging.info(f"创建/确认目录: {path}")

def load_cache(cache_path: str) -> Optional[Dict]:
    """加载缓存文件"""
    if os.path.exists(cache_path):
        try:
            with open(cache_path, 'rb') as f:
                return pickle.load(f)
        except Exception as e:
            logging.warning(f"缓存文件加载失败 {cache_path}: {e}")
    return None

def save_cache(data: Dict, cache_path: str) -> None:
    """保存缓存文件"""
    try:
        with open(cache_path, 'wb') as f:
            pickle.dump(data, f)
        logging.info(f"缓存已保存至: {cache_path}")
    except Exception as e:
        logging.error(f"缓存保存失败 {cache_path}: {e}")