import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
import matplotlib as mpl
# 读取数据
df = pd.read_csv('data/processed/123_filtered_data_micro.csv')

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
# 定义颜色
microbe_colors = ['#529AC6', '#F49600', '#29A15C', '#D8D3E7', '#9B9C2D', '#efe551', '#C0DFA2', '#F4849F']
metabolite_colors = ['#D8D3E7', '#efe551', '#C0DFA2']
phenotype_colors = ['#F6AF3A', '#EB716B', '#8982BB', '#A977A6', '#9B9C2D', '#D3E2E7', '#C0DFA2', '#88CEE6']

# 处理菌名称：取第5个分号后面的所有内容
def extract_microbe_name(name):
    if pd.isna(name):
        return 'Unknown'
    parts = str(name).split(';')
    if len(parts) > 5:
        return ';'.join(parts[5:]).strip()
    else:
        return str(name).strip()

# 应用处理函数
df.iloc[:, 0] = df.iloc[:, 0].apply(extract_microbe_name)

# 获取唯一的身体部位
body_sites = df.iloc[:, 9].unique()[:3]  # 取前3个身体部位

# 获取列名
microbe_col = df.columns[0]
metabolite_col = df.columns[1]
phenotype_col = df.columns[2]

# 创建一个空的图形
fig = go.Figure()

gap = 0.05  # 桑基图之间的间隔，可以根据需要调整这个值
chart_height = (1 - 2 * gap) / 3  # 每个桑基图的高度
# 为每个身体部位创建桑基图
for i, site in enumerate(body_sites):
    # 过滤当前身体部位的数据
    site_data = df[df.iloc[:, 9] == site]
    
    if site_data.empty:
        continue
    
    # 获取唯一的菌、代谢、表型
    microbes = site_data.iloc[:, 0].unique()
    metabolites = site_data.iloc[:, 1].unique()
    phenotypes = site_data.iloc[:, 2].unique()
    
    # 创建节点列表（菌 + 代谢 + 表型）
    all_nodes = list(microbes) + list(metabolites) + list(phenotypes)
    node_indices = {node: idx for idx, node in enumerate(all_nodes)}
    
    # 创建链接数据
    source = []
    target = []
    value = []
    link_colors = []  # 存储链接颜色
    
    # 统计菌到代谢的关系
    microbe_metabolite_counts = site_data.groupby([site_data.columns[0], site_data.columns[1]]).size()
    
    for (microbe, metabolite), count in microbe_metabolite_counts.items():
        source_idx = node_indices[microbe]
        target_idx = node_indices[metabolite]
        source.append(source_idx)
        target.append(target_idx)
        value.append(count)
        
        # 确定源节点颜色并添加透明度
        if microbe in microbes:
            color_idx = list(microbes).index(microbe) % len(microbe_colors)
            color = microbe_colors[color_idx]
        elif microbe in metabolites:
            color_idx = list(metabolites).index(microbe) % len(metabolite_colors)
            color = metabolite_colors[color_idx]
        else:
            color_idx = list(phenotypes).index(microbe) % len(phenotype_colors)
            color = phenotype_colors[color_idx]
        
        # 添加透明度 (50%)
        link_colors.append(f"rgba{tuple(int(color.lstrip('#')[j:j+2], 16) for j in (0, 2, 4)) + (0.3,)}")
    
    # 统计代谢到表型的关系
    metabolite_phenotype_counts = site_data.groupby([site_data.columns[1], site_data.columns[2]]).size()
    
    for (metabolite, phenotype), count in metabolite_phenotype_counts.items():
        source_idx = node_indices[metabolite]
        target_idx = node_indices[phenotype]
        source.append(source_idx)
        target.append(target_idx)
        value.append(count)
        
        # 确定源节点颜色并添加透明度
        if metabolite in microbes:
            color_idx = list(microbes).index(metabolite) % len(microbe_colors)
            color = microbe_colors[color_idx]
        elif metabolite in metabolites:
            color_idx = list(metabolites).index(metabolite) % len(metabolite_colors)
            color = metabolite_colors[color_idx]
        else:
            color_idx = list(phenotypes).index(metabolite) % len(phenotype_colors)
            color = phenotype_colors[color_idx]
        
        # 添加透明度 (50%)
        link_colors.append(f"rgba{tuple(int(color.lstrip('#')[j:j+2], 16) for j in (0, 2, 4)) + (0.3,)}")
    
    # 创建节点颜色
    node_colors = []
    for node in all_nodes:
        if node in microbes:
            # 循环使用菌的颜色
            color_idx = list(microbes).index(node) % len(microbe_colors)
            node_colors.append(microbe_colors[color_idx])
        elif node in metabolites:
            # 循环使用代谢的颜色
            color_idx = list(metabolites).index(node) % len(metabolite_colors)
            node_colors.append(metabolite_colors[color_idx])
        else:
            # 循环使用表型的颜色
            color_idx = list(phenotypes).index(node) % len(phenotype_colors)
            node_colors.append(phenotype_colors[color_idx])
    
    # 计算垂直位置
    y_start = 1 - (i+1) * chart_height - i * gap
    y_end = 1 - i * chart_height - i * gap
    
    # 创建桑基图
    sankey = go.Sankey(
        valueformat=".0f",
        valuesuffix="",
        # 定义节点
        node=dict(
            pad=15,
            thickness=20,
            line=dict(width=0),  # 去除黑色边框
            label=all_nodes,
            color=node_colors,
            hovertemplate='%{label}<extra></extra>'
        ),
        # 定义链接
        link=dict(
            source=source,
            target=target,
            value=value,
            color=link_colors,  # 使用源节点颜色并设置透明度
            hovertemplate='源: %{source.label}<br>目标: %{target.label}<br>数量: %{value}<extra></extra>'
        ),
        domain=dict(x=[0, 1], y=[y_start, y_end]),  # 垂直分布
        name=f"{site}"  # 添加名称以便在悬停时显示
    )
    
    fig.add_trace(sankey)
    
    # 添加桑基图标题
    fig.add_annotation(
        x=0.5,
        y=y_end + 0.05,  # 在桑基图上方
        xref="paper",
        yref="paper",
        text=f"{site}",
        showarrow=False,
        font=dict(size=19, color="black"),
        align="center"
    )
    
    # 添加列名
    # 菌列名
    fig.add_annotation(
        x=0.0,
        y=y_end + 0.01,
        xref="paper",
        yref="paper",
        text=microbe_col,
        showarrow=False,
        font=dict(size=12, color="black"),
        align="center"
    )
    
    # 代谢列名
    fig.add_annotation(
        x=0.5,
        y=y_end + 0.01,
        xref="paper",
        yref="paper",
        text=metabolite_col,
        showarrow=False,
        font=dict(size=12, color="black"),
        align="center"
    )
    
    # 表型列名
    fig.add_annotation(
        x=1,
        y=y_end + 0.02,
        xref="paper",
        yref="paper",
        text=phenotype_col,
        showarrow=False,
        font=dict(size=12, color="black"),
        align="center"
    )

# 更新布局
fig.update_layout(
    title_text="",
    title_x=0.5,
    font_size=12,
    height=1200,  # 增加高度以适应垂直排列
    width=1000,
    showlegend=False
)

# 保存为PNG文件
pio.write_image(fig, "sankey_diagrams_vertical.pdf", width=1000, height=1200, scale=2)
# fig.write_image("sankey_diagrams_vertical.png", width=1000, height=1200, scale=2)
# fig.write_image(fig, "sankey_diagrams_vertical.pdf", width=1000, height=1200, scale=2)
print("桑基图已生成并保存为 sankey_diagrams_vertical.png")