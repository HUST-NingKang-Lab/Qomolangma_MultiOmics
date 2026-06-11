import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import warnings
from matplotlib.patches import FancyArrowPatch
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.colors as mcolors

warnings.filterwarnings('ignore')

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

def load_data():
    """加载所有数据文件"""
    # 请确保目录下存在这些文件
    df_filtered = pd.read_csv('data/processed/123_filtered_data_with_microbiome.csv')
    gut_species = pd.read_csv('data/processed/gut_species.csv')
    oral_species = pd.read_csv('data/processed/oral_species.csv')
    skin_species = pd.read_csv('data/processed/skin_species.csv')
    lipid_data = pd.read_csv('data/processed/Lipid_log10.csv')
    metabolin_data = pd.read_csv('data/processed/Metabolin_log10.csv')
    cognitive_data = pd.read_csv('data/processed/cognitive.csv')
    physiology_data = pd.read_csv('data/processed/physiology.csv')
    return df_filtered, gut_species, oral_species, skin_species, lipid_data, metabolin_data, cognitive_data, physiology_data

def preprocess_data(df):
    """预处理数据"""
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    if len(numeric_cols) <= 1: 
        for col in df.columns:
            if col != 'subject_id':
                df[col] = pd.to_numeric(df[col], errors='coerce')
        numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    df_avg = df.groupby('subject_id')[numeric_cols].mean().reset_index()
    return df_avg

def filter_data_by_body_site(df_filtered, body_site):
    if body_site == 'all': return df_filtered
    else: return df_filtered[df_filtered.iloc[:, 9] == body_site]

def get_species_data(species_name, source, species_dict):
    if source == 'gut': return species_dict['gut']
    elif source == 'oral': return species_dict['oral']
    elif source == 'skin': return species_dict['skin']
    else: return None

def get_metabolite_data(metabolite_name, source, metabolite_dict):
    if source == 'Lipid': return metabolite_dict['lipid']
    elif source == 'Metabolite': return metabolite_dict['metabolite']
    else: return None

def get_phenotype_data(phenotype_name, source, phenotype_dict):
    if source == 'cognitive': return phenotype_dict['cognitive']
    elif source == 'physiology': return phenotype_dict['physiology']
    else: return None

def extract_microbe_label(full_name, df_filtered):
    matching_rows = df_filtered[df_filtered.iloc[:, 0] == full_name]
    if not matching_rows.empty and len(matching_rows.iloc[0]) > 12:
        col_13_value = matching_rows.iloc[0, 12]
        if pd.notna(col_13_value) and isinstance(col_13_value, str):
            parts = col_13_value.split(';')
            if len(parts) > 5: return ';'.join(parts[5:])
            else: return col_13_value
        else: return full_name
    else: return full_name

def calculate_average_values(df_filtered, species_dict, metabolite_dict, phenotype_dict, subject_filter=None):
    """计算节点平均值"""
    node_data = {}
    for idx, row in df_filtered.iterrows():
        microbe = row.iloc[0] 
        metabolite = row.iloc[1] 
        phenotype = row.iloc[2] 
        microbe_source = row.iloc[9] if len(row) > 9 else 'gut' 
        metabolite_source = row.iloc[10] if len(row) > 10 else 'Lipid' 
        phenotype_source = row.iloc[11] if len(row) > 11 else 'cognitive' 
        
        microbe_label = extract_microbe_label(microbe, df_filtered)
        
        # Microbe
        species_df = get_species_data(microbe, microbe_source, species_dict)
        microbe_avg = 0.1
        if species_df is not None:
            if subject_filter: species_df = species_df[species_df['subject_id'].isin(subject_filter)]
            if microbe in species_df.columns:
                microbe_avg = species_df.groupby('subject_id')[microbe].mean().mean() or 0.1

        # Metabolite
        metabolite_df = get_metabolite_data(metabolite, metabolite_source, metabolite_dict)
        metabolite_avg = 0.1
        if metabolite_df is not None:
            if subject_filter: metabolite_df = metabolite_df[metabolite_df['subject_id'].isin(subject_filter)]
            if metabolite in metabolite_df.columns:
                metabolite_avg = metabolite_df.groupby('subject_id')[metabolite].mean().mean() or 0.1
        
        # Phenotype
        phenotype_df = get_phenotype_data(phenotype, phenotype_source, phenotype_dict)
        phenotype_avg = 0.1
        if phenotype_df is not None:
            if subject_filter: phenotype_df = phenotype_df[phenotype_df['subject_id'].isin(subject_filter)]
            if phenotype in phenotype_df.columns:
                phenotype_avg = phenotype_df.groupby('subject_id')[phenotype].mean().mean() or 0.1
        
        node_data[microbe] = {'type': 'microbe', 'value': microbe_avg, 'label': microbe_label}
        node_data[metabolite] = {'type': 'metabolite', 'value': metabolite_avg, 'label': metabolite}
        node_data[phenotype] = {'type': 'phenotype', 'value': phenotype_avg, 'label': phenotype}
    
    return node_data

def create_network(df_filtered, node_data):
    """创建网络图"""
    G = nx.Graph()
    for node, data in node_data.items():
        G.add_node(node, **data)
    
    for idx, row in df_filtered.iterrows():
        microbe = row.iloc[0]
        metabolite = row.iloc[1]
        phenotype = row.iloc[2]
        
        try: r1 = float(row['Actual_Spearman_f12'])
        except: r1 = 0.0
        try: r2 = float(row['Actual_Spearman_f23'])
        except: r2 = 0.0
        
        G.add_edge(microbe, metabolite, weight=r1, type='microbe_metabolite')
        G.add_edge(metabolite, phenotype, weight=r2, type='metabolite_phenotype')
    return G

def get_edge_style(correlation):
    """
    样式计算
    |r| < 0.4: 灰色
    """
    abs_r = abs(correlation)
    width = max(abs_r * 4, 0.5)
    
    if abs_r < 0.4:
        return '#CCCCCC', 1.0
    
    norm_val = (abs_r - 0.4) / (1.0 - 0.4)
    norm_val = max(0, min(1, norm_val))
    color_intensity = 0.4 + (0.6 * norm_val) 
    
    if correlation > 0:
        return cm.Reds(color_intensity), width
    else:
        return cm.Blues(color_intensity), width

def draw_curved_edges(G, pos, ax, edges):
    """绘制弧形边"""
    for u, v in edges:
        u_pos = pos[u]
        v_pos = pos[v]
        correlation = G[u][v]['weight']
        color, linewidth = get_edge_style(correlation)
        
        arrow = FancyArrowPatch(
            u_pos, v_pos,
            arrowstyle='-',
            connectionstyle=f"arc3,rad=0.1",
            color=color,
            linewidth=linewidth,
            alpha=0.9,
            zorder=1
        )
        ax.add_patch(arrow)

def generate_randomized_concentric_layout(G):
    """
    随机同心布局
    - Metabolite: 随机分布在中心圆 (半径 0 ~ 0.3)
    - Microbe & Phenotype: 随机分布在外环 (半径 0.6 ~ 1.0)
    
    注意：没有设置 seed，所以每次调用都会产生不同的随机位置
    """
    pos = {}
    
    # 区分节点类型
    metabolite_nodes = [n for n, d in G.nodes(data=True) if d['type'] == 'metabolite']
    other_nodes = [n for n, d in G.nodes(data=True) if d['type'] != 'metabolite']
    
    # 1. 放置中心节点 (Metabolites) - 随机
    for node in metabolite_nodes:
        # 随机半径 [0, 0.3]
        r = np.random.uniform(0, 0.3)
        # 随机角度 [0, 2pi]
        theta = np.random.uniform(0, 2 * np.pi)
        pos[node] = np.array([r * np.cos(theta), r * np.sin(theta)])
        
    # 2. 放置外围节点 (Microbes & Phenotypes) - 随机
    for node in other_nodes:
        # 随机半径 [0.6, 1.0] - 留出0.3-0.6的空白带，让连线更清楚
        r = np.random.uniform(0.6, 1.0)
        theta = np.random.uniform(0, 2 * np.pi)
        pos[node] = np.array([r * np.cos(theta), r * np.sin(theta)])
        
    return pos

def plot_network(G, ax, title):
    """绘制网络图"""
    if len(G.nodes()) == 0:
        ax.text(0.5, 0.5, "无数据", ha='center', va='center', transform=ax.transAxes, fontsize=12)
        ax.set_title(title, fontsize=10)
        ax.axis('off')
        return ax
    
    # 使用修改后的随机同心布局
    pos = generate_randomized_concentric_layout(G)
    
    microbe_nodes = [n for n, d in G.nodes(data=True) if d['type'] == 'microbe']
    metabolite_nodes = [n for n, d in G.nodes(data=True) if d['type'] == 'metabolite']
    phenotype_nodes = [n for n, d in G.nodes(data=True) if d['type'] == 'phenotype']
    
    edges = G.edges()
    
    # 计算节点大小
    microbe_values = [G.nodes[n]['value'] for n in microbe_nodes]
    metabolite_values = [G.nodes[n]['value'] for n in metabolite_nodes]
    phenotype_values = [G.nodes[n]['value'] for n in phenotype_nodes]
    
    microbe_max = max(microbe_values) if microbe_values else 1
    metabolite_max = max(metabolite_values) if metabolite_values else 1
    phenotype_max = max(phenotype_values) if phenotype_values else 1
    
    microbe_sizes = [max(np.log(G.nodes[n]['value'] + 1) * 500 / np.log(microbe_max + 1), 50) for n in microbe_nodes]
    metabolite_sizes = [max(np.log(G.nodes[n]['value'] + 1) * 500 / np.log(metabolite_max + 1), 50) for n in metabolite_nodes]
    phenotype_sizes = [max(np.log(G.nodes[n]['value'] + 1) * 300 / np.log(phenotype_max + 1), 50) for n in phenotype_nodes]
    
    # 绘图
    draw_curved_edges(G, pos, ax, edges)
    
    # 绘制节点
    if microbe_nodes:
        nx.draw_networkx_nodes(G, pos, nodelist=microbe_nodes, node_color='#E29135',
                              node_size=microbe_sizes, alpha=1.0, edgecolors='white', linewidths=0.5, ax=ax)
    if metabolite_nodes:
        nx.draw_networkx_nodes(G, pos, nodelist=metabolite_nodes, node_color='#B8DBB3',
                              node_size=metabolite_sizes, alpha=1.0, edgecolors='white', linewidths=0.5, ax=ax)
    if phenotype_nodes:
        nx.draw_networkx_nodes(G, pos, nodelist=phenotype_nodes, node_color='#719AAC',
                              node_size=phenotype_sizes, alpha=1.0, edgecolors='white', linewidths=0.5, ax=ax)
    
    # 标签
    labels = {}
    for node in G.nodes():
        node_data = G.nodes[node]
        if 'label' in node_data:
            label = node_data['label']
            if len(label) > 30: label = label[:30] + '...'
            labels[node] = label
        else:
            labels[node] = node
    
    nx.draw_networkx_labels(G, pos, labels, font_size=6, font_family='sans-serif', ax=ax)
    
    ax.set_title(title, fontsize=12, pad=5)
    ax.axis('off')
    return ax

def create_legend(ax):
    """图例"""
    ax.axis('off')
    
    legend_elements = [
        mpl.lines.Line2D([0], [0], marker='o', color='w', markerfacecolor='#E29135', 
                         markersize=8, label='Microbiome (Outer)'),
        mpl.lines.Line2D([0], [0], marker='o', color='w', markerfacecolor='#B8DBB3', 
                         markersize=8, label='Metabolite (Center)'),
        mpl.lines.Line2D([0], [0], marker='o', color='w', markerfacecolor='#719AAC', 
                         markersize=8, label='Phenotype (Outer)')
    ]
    ax.legend(handles=legend_elements, loc='upper center', bbox_to_anchor=(0.5, 1.0), 
              frameon=False, fontsize=10, title="Node Types")
    
    cax = ax.inset_axes([0.2, 0.1, 0.6, 0.05])
    
    colors = []
    blues = plt.cm.Blues(np.linspace(1, 0.2, 60))
    colors.extend(blues)
    grays = [mcolors.to_rgba('#CCCCCC')] * 80
    colors.extend(grays)
    reds = plt.cm.Reds(np.linspace(0.2, 1, 60))
    colors.extend(reds)
    
    custom_cmap = mcolors.ListedColormap(colors)
    x = np.linspace(-1, 1, 200)
    
    cax.imshow([x], aspect='auto', cmap=custom_cmap, extent=[-1, 1, 0, 1])
    cax.set_yticks([])
    cax.set_xticks([-1, -0.4, 0, 0.4, 1])
    cax.set_xticklabels(['-1', '-0.4', '0', '0.4', '1'], fontsize=8)
    cax.set_title("Spearman Correlation", fontsize=10)
    
    ax.text(0.5, 0.25, "Blue: Negative (-)\nRed: Positive (+)\nGray: |r| < 0.4", 
            ha='center', va='top', fontsize=9, transform=ax.transAxes)

def main():
    print("正在加载数据...")
    df_filtered, gut_species, oral_species, skin_species, lipid_data, metabolin_data, cognitive_data, physiology_data = load_data()
    
    print("预处理数据...")
    gut_species_avg = preprocess_data(gut_species)
    oral_species_avg = preprocess_data(oral_species)
    skin_species_avg = preprocess_data(skin_species)
    lipid_data_avg = preprocess_data(lipid_data)
    metabolin_data_avg = preprocess_data(metabolin_data)
    cognitive_data_avg = preprocess_data(cognitive_data)
    physiology_data_avg = preprocess_data(physiology_data)
    
    species_dict = {'gut': gut_species_avg, 'oral': oral_species_avg, 'skin': skin_species_avg}
    metabolite_dict = {'lipid': lipid_data_avg, 'metabolite': metabolin_data_avg}
    phenotype_dict = {'cognitive': cognitive_data_avg, 'physiology': physiology_data_avg}
    
    elderly_subjects = ['S1', 'S10']
    young_subjects = ['S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8', 'S9', 'S11']
    female_subjects = ['S6', 'S7', 'S8', 'S11']
    male_subjects = ['S1', 'S2', 'S3', 'S4', 'S5', 'S9', 'S10']
    body_sites = ['gut', 'oral', 'skin']
    
    plots_info = [
        (None, "All"),
        (elderly_subjects, "Old"),
        (young_subjects, "Young"),
        (male_subjects, "Male"),
        (female_subjects, "Female")
    ]
    
    for body_site in body_sites:
        print(f"\n正在处理身体部位: {body_site}")
        df_body_site = filter_data_by_body_site(df_filtered, body_site)
        if len(df_body_site) == 0: continue
        
        fig, axes = plt.subplots(2, 3, figsize=(15, 10)) 
        axes = axes.flatten()
        networks = []
        
        for i, (subject_filter, title) in enumerate(plots_info):
            print(f"  正在处理分组: {title}")
            try:
                node_data = calculate_average_values(df_body_site, species_dict, metabolite_dict, phenotype_dict, subject_filter)
                G = create_network(df_body_site, node_data)
                networks.append((G, title))
            except Exception as e:
                print(f"  出错: {e}")
                networks.append((nx.Graph(), title))
        
        for i, (G, title) in enumerate(networks):
            plot_network(G, axes[i], title)
        
        create_legend(axes[5])
        filename = f"{body_site}_network_concentric_random.pdf"
        fig.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"已保存: {filename}")
        plt.show()

if __name__ == "__main__":
    main()