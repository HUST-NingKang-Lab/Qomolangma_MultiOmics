"""
最终版：筛选真正具有年龄差异的活跃微生物
条件：
1. proportion_divergent > 0
2. 至少一个时间点 mean(p_adj) < 0.05
输出表格：site | MAG | Microbe | mean_plasticity
"""

import pandas as pd
from pathlib import Path

# ==================== 路径设置 ====================
INPUT_DIR = "data/processed"
OUTPUT_DIR = "data/processed"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ====================

# ==================== 加载数据 ====================
print("正在加载数据...")
microbe_summary = pd.read_csv(INPUT_DIR / "counterfactual_summary_by_microbe_annotated.csv")
full_data = pd.read_csv(INPUT_DIR / "counterfactual_predictions_annotated.csv")

print(f"汇总表行数: {len(microbe_summary)}")
print(f"完整预测表行数: {len(full_data)}")

# ==================== 筛选活跃微生物 ====================
active_microbes = []

for site in ['gut', 'oral', 'skin']:
    print(f"\n处理 {site.upper()}...")
    
    # 1. 从汇总表筛选 proportion_divergent > 0 的微生物
    candidates = microbe_summary[
        (microbe_summary['site'] == site) &
        (microbe_summary['proportion_divergent'] > 0)
    ][['microbe', 'mean_plasticity']].drop_duplicates()
    
    if len(candidates) == 0:
        print(f"  {site}: 无 proportion_divergent > 0 的微生物")
        continue
    
    print(f"  候选微生物数量（proportion_divergent > 0）: {len(candidates)}")
    
    # 2. 对每个候选微生物，检查是否至少有一个时间点 p_adj < 0.05
    significant_microbes = []
    
    for microbe in candidates['microbe']:
        subset = full_data[
            (full_data['site'] == site) &
            (full_data['microbe'] == microbe)
        ]
        
        if len(subset) == 0:
            continue
            
        # 每个时间点平均 p_adj
        p_per_time = subset.groupby('time')['p_adj'].mean()
        
        # 判断是否有至少一个时间点显著
        if (p_per_time < 0.05).any():
            mean_plas = candidates[candidates['microbe'] == microbe]['mean_plasticity'].iloc[0]
            significant_microbes.append({
                'site': site,
                'microbe': microbe,
                'mean_plasticity': mean_plas
            })
    
    print(f"  真正显著（≥1 time point p_adj<0.05）: {len(significant_microbes)} 个")
    
    active_microbes.extend(significant_microbes)

# 转为 DataFrame
active_df = pd.DataFrame(active_microbes)

if len(active_df) == 0:
    print("\n警告：未找到任何满足条件的活跃微生物！")
else:
    print(f"\n总共筛选出 {len(active_df)} 个真正具有年龄差异的活跃微生物")

# ==================== 合并 MAG 信息（如果有） ====================
# 尝试从原始表格中提取 MAG 列（常见列名：MAG, genome, bin 等）
mag_column = None
possible_mag_cols = ['MAG', 'genome', 'bin', 'contig', 'assembly']
for col in possible_mag_cols:
    if col in full_data.columns:
        mag_column = col
        break

if mag_column:
    print(f"检测到 MAG 列: {mag_column}")
    mag_map = full_data[['microbe', mag_column]].drop_duplicates()
    active_df = active_df.merge(mag_map, on='microbe', how='left')
    # 重命名
    active_df = active_df.rename(columns={mag_column: 'MAG'})
    # 调整列顺序
    active_df = active_df[['site', 'MAG', 'microbe', 'mean_plasticity']]
else:
    print("未检测到 MAG 列，将使用 microbe 作为标识")
    active_df['MAG'] = active_df['microbe']  # 占位
    active_df = active_df[['site', 'MAG', 'microbe', 'mean_plasticity']]

# 重命名列为你想要的
active_df.columns = ['site', 'MAG', 'Microbe', 'mean_plasticity']

# 排序：按 |mean_plasticity| 降序
active_df = active_df.reindex(active_df['mean_plasticity'].abs().sort_values(ascending=False).index).reset_index(drop=True)

# ==================== 保存结果 ====================
output_file = OUTPUT_DIR / "active_age_divergent_microbes.csv"
active_df.to_csv(output_file, index=False)
print(f"\n已保存活跃微生物表格：")
print(f"   {output_file}")
print(f"   共 {len(active_df)} 行数据")

# 同时保存一个带显著性排序的版本（可选）
top_by_site = []
for site in ['gut', 'oral', 'skin']:
    subset = active_df[active_df['site'] == site].copy()
    subset = subset.reindex(subset['mean_plasticity'].abs().sort_values(ascending=False).index)
    top_by_site.append(subset.head(20))  # 每个site前20

pd.concat(top_by_site).to_csv(OUTPUT_DIR / "active_age_divergent_microbes_top20_per_site.csv", index=False)

# ==================== 打印统计 ====================
print("\n" + "="*60)
print("最终统计（真正活跃的年龄差异微生物）：")
print(active_df['site'].value_counts().sort_index())
print("\n前10名（按 |plasticity| 最大）：")
print(active_df.head(10)[['site', 'Microbe', 'mean_plasticity']].to_string(index=False))

print("\n所有文件已保存至：")
print(OUTPUT_DIR)
print("\nDone！这批微生物可以自信地称为：Age-associated dynamically active microbes")