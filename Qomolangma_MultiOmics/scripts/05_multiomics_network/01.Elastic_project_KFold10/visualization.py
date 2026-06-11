import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import logging
import os

plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'SimHei']
plt.rcParams['axes.unicode_minus'] = False

def create_visualizations(stable_features_df: pd.DataFrame, metabolite_name: str, 
                         figure_dir: str) -> None:
    """创建可视化图表"""
    try:
        # 绘图：特征选择频率
        plt.figure(figsize=(12, 8))
        top_n = min(20, len(stable_features_df))
        if top_n > 0:
            sns.barplot(
                data=stable_features_df.head(top_n), 
                x='Selection_Frequency', y='Feature'
            )
            plt.title(f'{metabolite_name} - Top {top_n} Stable Features')
            plt.xlabel('Selection Frequency')
            plt.tight_layout()
            plt.savefig(
                os.path.join(figure_dir, 'feature_selection_frequency.png'),
                dpi=300, bbox_inches='tight'
            )
            plt.close()
        
        # 绘图：权重分布
        if len(stable_features_df) > 0:
            plt.figure(figsize=(12, 8))
            weights_df = stable_features_df[stable_features_df['Average_Weight'] != 0].copy()
            if len(weights_df) > 0:
                top_n_weights = min(20, len(weights_df))
                weights_df_sorted = weights_df.reindex(
                    weights_df['Average_Weight'].abs().sort_values(ascending=False).index
                )
                
                colors = ['red' if x < 0 else 'blue' for x in weights_df_sorted.head(top_n_weights)['Average_Weight']]
                
                plt.barh(range(top_n_weights), weights_df_sorted.head(top_n_weights)['Average_Weight'], color=colors)
                plt.yticks(range(top_n_weights), weights_df_sorted.head(top_n_weights)['Feature'])
                plt.xlabel('Average Weight')
                plt.title(f'{metabolite_name} - Feature Weights (Top {top_n_weights})')
                plt.tight_layout()
                plt.savefig(
                    os.path.join(figure_dir, 'feature_weights.png'),
                    dpi=300, bbox_inches='tight'
                )
                plt.close()
    except Exception as e:
        logging.warning(f"图表创建失败 {metabolite_name}: {e}")

def create_summary_visualizations(summary_df: pd.DataFrame, feature_metabolite_counts: pd.DataFrame, 
                                figure_dir: str) -> None:
    """创建汇总可视化图表"""
    try:
        # 1. 代谢物稳定特征数量分布
        plt.figure(figsize=(12, 8))
        top_metabolites = summary_df.head(20)
        sns.barplot(data=top_metabolites, x='Stable_Feature_Count', y='Metabolite')
        plt.title('各代谢物的稳定特征数量 (Top 20)')
        plt.xlabel('稳定特征数量')
        plt.tight_layout()
        plt.savefig(os.path.join(figure_dir, 'metabolite_feature_counts.png'), dpi=300, bbox_inches='tight')
        plt.close()
        
        # 2. 最常被选择的特征
        plt.figure(figsize=(12, 8))
        top_features = feature_metabolite_counts.head(20)
        sns.barplot(data=top_features, x='Metabolite_Count', y='Feature')
        plt.title('在最多代谢物中被选择的特征 (Top 20)')
        plt.xlabel('被选择的代谢物数量')
        plt.tight_layout()
        plt.savefig(os.path.join(figure_dir, 'top_selected_features.png'), dpi=300, bbox_inches='tight')
        plt.close()
        
        # 3. 稳定特征数量分布直方图
        plt.figure(figsize=(10, 6))
        plt.hist(summary_df['Stable_Feature_Count'], bins=20, alpha=0.7, edgecolor='black')
        plt.title('稳定特征数量分布')
        plt.xlabel('稳定特征数量')
        plt.ylabel('代谢物数量')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(figure_dir, 'feature_count_distribution.png'), dpi=300, bbox_inches='tight')
        plt.close()
        
    except Exception as e:
        logging.error(f"汇总图表创建失败: {e}")