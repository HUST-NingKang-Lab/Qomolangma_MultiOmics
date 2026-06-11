"""
数据预处理模块
将微生物丰度数据转换为适合贝叶斯建模的格式
"""

import pandas as pd
import numpy as np
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# ==================== 配置 ====================
OUTPUT_DIR = "data/processed"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("开始数据预处理...")
print("=" * 60)

# ==================== 数据导入 ====================
# 年龄信息
metadata = pd.read_csv(
     "data/metadata/Matedata_Information.csv"
)
pair_subject_age = metadata[['people', 'age_group']].drop_duplicates()
pair_subject_age.columns = ['subject_id', 'age_group']

# 研究对象和时间点
people = [f"S{i}" for i in range(1, 14)]
time12 = ['Ta', 'Tb'] + [f'T{i}' for i in range(1, 11)]

# 时间映射字典
time_mapping = {
    'T1': 1, 'T2': 2, 'Ta': 3, 'T3': 4, 'Tb': 5, 'T4': 6,
    'T5': 7, 'T6': 8, 'T7': 9, 'T8': 10, 'T9': 11, 'T10': 12
}

# 读取差异分析结果
Lipid_DE = pd.read_csv("data/processed/Lipid_DE_results.csv").query("qval < 0.05")
Metabolin_DE = pd.read_csv("data/processed/Metabolin_DE_results.csv").query("qval < 0.05")

# 读取丰度表

def load_abundance(file_path, de_results, people, time12):
    """加载并筛选丰度表"""
    df = pd.read_csv(file_path)
    sig_microbes = de_results['original_name'].unique()
    cols = ['subject_id', 'time'] + [col for col in sig_microbes if col in df.columns]
    df = df[cols]
    df = df[df['subject_id'].isin(people) & df['time'].isin(time12)]
    return df

Lipid_abs = load_abundance("data/processed/Lipid_log10.csv", Lipid_DE, people, time12)
Metabolin_abs = load_abundance("data/processed/Metabolin_log10.csv", Metabolin_DE, people, time12)

print(f"加载的微生物数: Lipid={len(Lipid_abs.columns)-2}, Metabolin={len(Metabolin_abs.columns)-2}")

# ==================== 数据预处理函数 ====================
def preprocess_data(abs_data, site_name, time_mapping, pair_subject_age):
    """
    预处理微生物丰度数据
    
    参数:
        abs_data: 丰度表 (宽格式)
        site_name: 部位名称 ('Metabolin', 'Lipid')
        time_mapping: 时间点映射字典
        pair_subject_age: 年龄分组信息
    
    返回:
        长格式数据框，包含log转换后的丰度
    """
    # 获取微生物列
    microbe_cols = [col for col in abs_data.columns 
                   if col not in ['subject_id', 'time']]
    

    df = abs_data.copy()
    
    # 转换为长格式
    df_long = df.melt(
        id_vars=['subject_id', 'time'],
        value_vars=microbe_cols,
        var_name='Metabolite',
        value_name='log_abundance'
    )

    # 对同一个 subject_id + time + microbe 取平均（如果有多次采样）
    df_long = (
        df_long
        .groupby(['subject_id', 'time', 'Metabolite'], as_index=False)
        ['log_abundance']
        .mean()   # 也可以用 .median()，更稳健
    )

    # 添加时间数值
    df_long['time_numeric'] = df_long['time'].map(time_mapping)
    df_long['site'] = site_name
    
    # 合并年龄信息
    df_long = df_long.merge(pair_subject_age, on='subject_id', how='left')
    df_long = df_long.dropna(subset=['age_group'])
    
    # 添加年龄编码
    df_long['age_numeric'] = (df_long['age_group'] == 'old').astype(int)
    
    # 数据类型优化
    df_long['subject_id'] = df_long['subject_id'].astype('category')
    df_long['Metabolite'] = df_long['Metabolite'].astype('category')
    df_long['site'] = df_long['site'].astype('category')
    df_long['age_group'] = df_long['age_group'].astype('category')
    
    return df_long

# ==================== 处理三个部位 ====================
print("\n处理各部位数据...")
Lipid_processed = preprocess_data(Lipid_abs, 'Lipid', time_mapping, pair_subject_age)
Metabolin_processed = preprocess_data(Metabolin_abs, 'Metabolin', time_mapping, pair_subject_age)

# 合并所有数据
all_data = pd.concat([Lipid_processed, Metabolin_processed], 
                     ignore_index=True)

# ==================== 数据质量检查 ====================
print("\n数据质量检查:")
print(f"总观测数: {len(all_data):,}")
print(f"显著代谢物数: {all_data['Metabolite'].nunique()}")
print(f"个体数: {all_data['subject_id'].nunique()}")
print(f"\n年龄分布:")
print(all_data['age_group'].value_counts())
print(f"\n部位分布:")
print(all_data['site'].value_counts())

# 检查缺失值
missing_summary = all_data.isnull().sum()
if missing_summary.sum() > 0:
    print("\n警告: 发现缺失值")
    print(missing_summary[missing_summary > 0])

# 检查无穷值
inf_check = np.isinf(all_data['log_abundance']).sum()
if inf_check > 0:
    print(f"\n警告: 发现{inf_check}个无穷值")

# ==================== 保存预处理数据 ====================
output_path =  "data/processed/preprocessed_data_for_counterfactual.csv"
all_data.to_csv(output_path, index=False)
print(f"\n预处理数据已保存至: {output_path}")

# 保存元数据摘要
summary_stats = all_data.groupby(['site', 'age_group']).agg({
    'subject_id': 'nunique',
    'Metabolite': 'nunique',
    'log_abundance': ['count', 'mean', 'std']
}).round(3)

summary_stats.to_csv( "data/processed/preprocessing_summary.csv")
print("摘要统计已保存")

# ==================== 可视化数据分布 ====================
try:
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    plt.style.use('seaborn-v0_8-whitegrid')
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    
    for idx, site in enumerate(['Lipid', 'Metabolin']):
        site_data = all_data[all_data['site'] == site]
        axes[idx].hist(site_data['log_abundance'], bins=50, alpha=0.7, edgecolor='black')
        axes[idx].set_title(f'{site.capitalize()} - Log Abundance', fontsize=12, fontweight='bold')
        axes[idx].set_xlabel('Log Abundance')
        axes[idx].set_ylabel('Frequency')
        axes[idx].grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "data_distribution.png", dpi=300, bbox_inches='tight')
    print("数据分布图已保存")
    
except ImportError:
    print("matplotlib未安装，跳过可视化")

print("\n" + "=" * 60)
print("数据预处理完成！")
print("=" * 60)