"""
分析活跃菌方向性分布：elderly-high vs young-high vs Normal
输入：active_age_divergent_microbes.csv
输出：三个site的饼图（oral/gut/skin）
"""

import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# ==================== 路径设置 ====================
BASE_DIR = Path("/public/home/xiaokechen/01.HuaDa_Qomolangma")
ACTIVE_DIR = BASE_DIR / "data/20250901/06.贝叶斯论证具有年龄差异的时间差异菌/05.活跃菌"

active_file = "data/processed/active_age_divergent_microbes.csv"
output_dir = ACTIVE_DIR
output_dir.mkdir(parents=True, exist_ok=True)

# 总菌数（你提供的真实背景总数）
TOTAL_COUNTS = {
    'oral': 103,
    'gut' : 22,
    'skin': 6
}

COLORS = {
    'oral': '#4DBBD5',   # 青蓝
    'gut' : '#E64B35',   # 红橙
    'skin': '#00A087'    # 墨绿
}

DIRECTION_COLORS = {
    'elderly-high active bacteria': '#C1666B',  # 偏红（老人富集）
    'young-high active bacteria'   : '#4281A4',  # 偏蓝（年轻人富集）
    'Normal'                      : '#DDDDDD'   # 灰色
}

# ==================== 加载并聚合活跃菌 ====================
print("正在读取活跃菌表格...")
df = pd.read_csv(active_file)

print(f"原始活跃菌数量: {len(df)}")
print("按 site + MAG 去重并聚合 mean_plasticity...")

# 关键：同一个 MAG 在多个样本中可能出现，取平均 plasticity
df_agg = df.groupby(['site', 'MAG'], as_index=False)['mean_plasticity'].mean()

print(f"聚合后唯一 MAG 数量: {len(df_agg)}")
print(df_agg['site'].value_counts().sort_index())

# ==================== 分类并补足 Normal ====================
results = []

for site in ['oral', 'gut', 'skin']:
    site_df = df_agg[df_agg['site'] == site].copy()
    total = TOTAL_COUNTS[site]
    
    # 分类
    elderly = len(site_df[site_df['mean_plasticity'] > 0.5])
    young   = len(site_df[site_df['mean_plasticity'] < -0.5])
    normal  = total - elderly - young
    
    if normal < 0:
        print(f"警告: {site} 活跃菌总数超过背景总数！将 Normal 设为 0")
        normal = 0
    
    results.append({
        'site': site,
        'elderly-high active bacteria': elderly,
        'young-high active bacteria'   : young,
        'Normal'                      : normal
    })

pie_df = pd.DataFrame(results).set_index('site')
print("\n各部位活跃菌方向性分布：")
print(pie_df)

# ==================== 绘图：三个饼图 ====================
fig, axes = plt.subplots(1, 3, figsize=(15, 5), dpi=300)
axes = axes.flatten()

for idx, site in enumerate(['oral', 'gut', 'skin']):
    ax = axes[idx]
    data = pie_df.loc[site]
    labels = [f"{k}\n({v})" for k, v in data.items()]
    sizes = data.values
    colors = [DIRECTION_COLORS[label.split('\n')[0]] for label in labels]
    
    wedges, texts, autotexts = ax.pie(
        sizes,
        labels=labels,
        colors=colors,
        autopct=lambda pct: f'{pct:.1f}%\n({int(pct/100*sum(sizes))})' if pct > 0 else '',
        startangle=90,
        textprops={'fontsize': 12, 'fontweight': 'bold'},
        wedgeprops=dict(width=0.4, edgecolor='white', linewidth=2)
    )
    
    # 美化文字颜色
    for autotext in autotexts:
        autotext.set_color('black')
        autotext.set_fontweight('bold')
        autotext.set_fontsize(11)
    
    ax.set_title(f'{site.capitalize()}\n(Total active: {TOTAL_COUNTS[site]})',
                 fontsize=16, fontweight='bold', pad=20)

plt.suptitle('Directionality of Age-Associated Dynamically Active Microbes Across Body Sites',
             fontsize=20, fontweight='bold', y=1.05)

plt.tight_layout()

# 保存
plt.savefig("Fig_active_bacteria_direction_pie_charts.pdf", bbox_inches='tight', dpi=400)
plt.savefig("Fig_active_bacteria_direction_pie_charts.png", bbox_inches='tight', dpi=400)
plt.close()

print("\n饼图已保存：")
print("   Fig_active_bacteria_direction_pie_charts.pdf")
print("   Fig_active_bacteria_direction_pie_charts.png")

# 同时保存统计表
pie_df.to_csv(output_dir / "summary_active_bacteria_direction_by_site.csv")
print("   summary_active_bacteria_direction_by_site.csv")

print("\n" + "="*60)
print("最终结果汇总：")
for site in ['oral', 'gut', 'skin']:
    e = pie_df.loc[site, 'elderly-high active bacteria']
    y = pie_df.loc[site, 'young-high active bacteria']
    n = pie_df.loc[site, 'Normal']
    print(f"{site.capitalize():5}: Elderly-high = {e:2d} | Young-high = {y:2d} | Normal = {n:2d}")

print("\n分析完成！")