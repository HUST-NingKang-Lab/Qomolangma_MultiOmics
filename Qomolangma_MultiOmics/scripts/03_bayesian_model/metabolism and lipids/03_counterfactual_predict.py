"""
反事实预测（手动计算最终版）
完全不依赖model.predict()，直接从后验参数计算
最稳定、最可控的方案
"""

import pandas as pd
import numpy as np
from pathlib import Path
import pickle
import warnings
warnings.filterwarnings('ignore')

import arviz as az
from patsy import dmatrix
from joblib import Parallel, delayed
import time

# ==================== 配置 ====================
INPUT_DIR = "data/processed"
MODEL_DIR = "scripts/03_bayesian_model/model"
OUTPUT_DIR = "data/processed"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

N_JOBS = 16
N_POSTERIOR_SAMPLES = 500

print("=" * 60)
print("反事实预测（手动计算最终版）")
print("=" * 60)

# ==================== 加载数据 ====================
print("\n加载数据...")
all_data = pd.read_csv(INPUT_DIR / "preprocessed_data_for_counterfactual.csv")

print("加载模型元数据...")
with open(MODEL_DIR / "young_only_models_metadata.pkl", 'rb') as f:
    young_models = pickle.load(f)

successful_models = [m for m in young_models if m['success']]
print(f"成功加载 {len(successful_models)} 个年轻人模型")

old_data = all_data[all_data['age_group'] == 'old'].copy()
young_data = all_data[all_data['age_group'] == 'youth'].copy()

print(f"\n老年人数据: {len(old_data):,} 观测点")
print(f"年轻人数据: {len(young_data):,} 观测点")

# ==================== 手动预测函数 ====================
def manual_predict_from_posterior(model_metadata, old_subjects_data):
    """
    从后验参数直接计算预测
    
    基于诊断结果，参数结构为:
    - Intercept: (chains, draws) = (2, 1000)
    - bs(time_scaled, df=3): (chains, draws, 3) = (2, 1000, 3)
    - sigma: (chains, draws) = (2, 1000)
    
    模型: log_abundance ~ Intercept + sum(β_i * bs_i(time))
    """
    site = model_metadata['site']
    Metabolite = model_metadata['Metabolite']
    
    # 筛选老年人数据
    old_subset = old_subjects_data[
        (old_subjects_data['site'] == site) & 
        (old_subjects_data['Metabolite'] == Metabolite)
    ].copy()
    
    if len(old_subset) == 0:
        return None
    
    try:
        # === 步骤1: 加载后验参数 ===
        idata_path = model_metadata['idata_path']
        idata = az.from_netcdf(idata_path)
        posterior = idata.posterior
        available_params = list(posterior.data_vars.keys())
        
        # === 步骤2: 标准化老年人的时间 ===
        time_mean = model_metadata['time_mean']
        time_std = model_metadata['time_std']
        old_subset['time_scaled'] = (old_subset['time_numeric'] - time_mean) / time_std
        
        # === 步骤3: 构建样条基函数矩阵 ===
        # 使用patsy生成与Bambi一致的样条基
        X_spline = dmatrix(
            "bs(time_scaled, df=3) - 1",  # -1移除截距
            data=old_subset,
            return_type='dataframe'
        )
        
        # === 步骤4: 提取后验参数（基于诊断结果）===
        # 截距
        if 'Intercept' not in available_params:
            print(f"警告: {site}-{Metabolite} 缺少Intercept")
            return None
        
        intercept_samples = posterior['Intercept'].values.flatten()
        
        # **样条系数：确认是 'bs(time_scaled, df=3)' 格式，shape=(chains, draws, 3)**
        spline_key = 'bs(time_scaled, df=3)'
        
        if spline_key not in available_params:
            print(f"警告: {site}-{Metabolite} 缺少样条参数")
            return None
        
        spline_data = posterior[spline_key].values  # shape: (2, 1000, 3)
        
        # 提取3个样条系数，每个是 (chains * draws) 的一维数组
        spline_coef_samples = []
        for i in range(3):
            coef = spline_data[..., i].flatten()  # 提取第i个系数
            spline_coef_samples.append(coef)
        
        spline_coef_samples = np.array(spline_coef_samples)  # shape: (3, 2000)
        
        # Sigma（观测误差）
        if 'sigma' not in available_params:
            print(f"提示: {site}-{Metabolite} 未找到sigma，使用默认值1")
            sigma_samples = np.ones_like(intercept_samples)
        else:
            sigma_samples = posterior['sigma'].values.flatten()
        
        # 限制样本数
        n_samples = min(N_POSTERIOR_SAMPLES, len(intercept_samples))
        intercept_samples = intercept_samples[:n_samples]
        spline_coef_samples = spline_coef_samples[:, :n_samples]
        sigma_samples = sigma_samples[:n_samples]
        
        # === 步骤5: 计算后验预测分布 ===
        n_obs = len(old_subset)
        predictions = np.zeros((n_samples, n_obs))
        
        for s in range(n_samples):
            # 固定效应预测（不包含随机效应）
            mu = intercept_samples[s] + X_spline.values @ spline_coef_samples[:, s]
            
            # 后验预测分布（加入观测误差）
            y_pred = np.random.normal(mu, sigma_samples[s], size=n_obs)
            
            predictions[s, :] = y_pred
        
        # === 步骤6: 计算后验统计量 ===
        counter_mean = predictions.mean(axis=0)
        counter_lower = np.percentile(predictions, 2.5, axis=0)
        counter_upper = np.percentile(predictions, 97.5, axis=0)
        counter_sd = predictions.std(axis=0)
        
        # === 步骤7: 构建结果DataFrame ===
        results = pd.DataFrame({
            'site': site,
            'Metabolite': Metabolite,
            'subject_id': old_subset['subject_id'].values,
            'time': old_subset['time'].values,
            'time_numeric': old_subset['time_numeric'].values,
            'observed_log_abundance': old_subset['log_abundance'].values,
            'counterfactual_mean': counter_mean,
            'counterfactual_lower': counter_lower,
            'counterfactual_upper': counter_upper,
            'counterfactual_sd': counter_sd
        })
        
        # === 步骤8: 计算年龄可塑性指数 ===
        results['plasticity_index'] = (
            results['observed_log_abundance'] - results['counterfactual_mean']
        )
        
        # 判断显著偏离
        results['is_divergent'] = (
            (results['observed_log_abundance'] < results['counterfactual_lower']) |
            (results['observed_log_abundance'] > results['counterfactual_upper'])
        )
        
        # 偏离方向
        results['divergence_direction'] = 'none'
        results.loc[results['observed_log_abundance'] < results['counterfactual_lower'], 
                   'divergence_direction'] = 'lower'
        results.loc[results['observed_log_abundance'] > results['counterfactual_upper'], 
                   'divergence_direction'] = 'upper'
        
        # 零膨胀标记
        results['young_model_predicts_zero'] = counter_mean < -15
        results['observed_is_zero'] = results['observed_log_abundance'] < -15
        results['zero_consistency'] = (
            results['young_model_predicts_zero'] == results['observed_is_zero']
        )
        
        # === 步骤9: 计算p值 ===
        p_values = []
        for i in range(n_obs):
            obs = results.iloc[i]['observed_log_abundance']
            samples = predictions[:, i]
            
            # 双侧p值
            p_lower = np.mean(samples >= obs)
            p_upper = np.mean(samples <= obs)
            p_value = 2 * min(p_lower, p_upper)
            p_values.append(p_value)
        
        results['p_value'] = p_values
        
        return results
        
    except Exception as e:
        error_msg = str(e)
        # 只打印关键错误
        if len(error_msg) > 100:
            error_msg = error_msg[:100] + "..."
        print(f"失败: {site}-{Metabolite}: {error_msg}")
        return None

# ==================== 并行计算 ====================
print(f"\n开始并行计算 (使用 {N_JOBS} 个进程)...")
start_time = time.time()

results = Parallel(n_jobs=N_JOBS, verbose=10)(
    delayed(manual_predict_from_posterior)(
        model_metadata, old_data
    )
    for model_metadata in successful_models
)

# 过滤并合并
successful_predictions = [r for r in results if r is not None]
failed_predictions = len(results) - len(successful_predictions)

elapsed = time.time() - start_time

print(f"\n并行计算完成，用时: {elapsed/60:.2f} 分钟")
print(f"成功: {len(successful_predictions)}/{len(results)} ({len(successful_predictions)/len(results)*100:.1f}%)")
print(f"失败: {failed_predictions}")

if len(successful_predictions) == 0:
    print("\n错误：所有反事实预测都失败了！")
    print("请检查idata文件和后验参数命名")
    exit(1)

counterfactual_df = pd.concat(successful_predictions, ignore_index=True)

# ==================== FDR校正 ====================
print("\n进行FDR校正...")
from statsmodels.stats.multitest import multipletests

_, p_adj, _, _ = multipletests(
    counterfactual_df['p_value'], 
    alpha=0.05, 
    method='fdr_bh'
)
counterfactual_df['p_adj'] = p_adj
counterfactual_df['is_divergent_adj'] = p_adj < 0.05

# ==================== 结果汇总 ====================
print("\n" + "=" * 60)
print("反事实分析完成")
print("=" * 60)

print(f"\n基本统计:")
print(f"  预测的微生物数: {counterfactual_df['Metabolite'].nunique()}")
print(f"  老年人观测点总数: {len(counterfactual_df):,}")
print(f"  显著偏离 (FDR<0.05): {counterfactual_df['is_divergent_adj'].sum()} ({counterfactual_df['is_divergent_adj'].mean()*100:.1f}%)")

print(f"\n可塑性指数统计:")
print(f"  均值: {counterfactual_df['plasticity_index'].mean():.4f}")
print(f"  中位数: {counterfactual_df['plasticity_index'].median():.4f}")
print(f"  标准差: {counterfactual_df['plasticity_index'].std():.4f}")
print(f"  范围: [{counterfactual_df['plasticity_index'].min():.2f}, {counterfactual_df['plasticity_index'].max():.2f}]")

# 各部位统计
print("\n各部位偏离统计:")
site_summary = counterfactual_df.groupby('site').agg({
    'Metabolite': 'nunique',
    'plasticity_index': ['count', 'mean', 'std'],
    'is_divergent_adj': ['sum', 'mean']
}).round(3)
site_summary.columns = ['n_Metabolites', 'n_obs', 'mean_plasticity', 'sd_plasticity', 'n_divergent', 'prop_divergent']
print(site_summary)

# 零值一致性
print("\n零值一致性分析:")
zero_stats = counterfactual_df.groupby('site').agg({
    'zero_consistency': 'mean',
    'young_model_predicts_zero': 'mean',
    'observed_is_zero': 'mean'
}).round(3)
zero_stats.columns = ['一致性', '年轻预测零', '老年观测零']
print(zero_stats)

# ==================== 保存结果 ====================
print("\n保存结果...")

counterfactual_df.to_csv(
    OUTPUT_DIR / "counterfactual_predictions.csv",
    index=False
)

# 微生物水平汇总
Metabolite_summary = counterfactual_df.groupby(['site', 'Metabolite']).agg({
    'time_numeric': 'count',
    'plasticity_index': ['mean', 'std', 'min', 'max'],
    'is_divergent': 'sum',
    'is_divergent_adj': 'sum',
    'p_value': 'mean',
    'p_adj': 'mean',
    'zero_consistency': 'mean'
}).round(4)

Metabolite_summary.columns = [
    'n_timepoints', 
    'mean_plasticity', 'sd_plasticity', 'min_plasticity', 'max_plasticity',
    'n_divergent_raw', 'n_divergent_adj',
    'mean_p', 'mean_p_adj',
    'zero_consistency'
]

Metabolite_summary['proportion_divergent'] = (
    Metabolite_summary['n_divergent_adj'] / Metabolite_summary['n_timepoints']
).round(3)

Metabolite_summary = Metabolite_summary.sort_values('mean_plasticity', key=abs, ascending=False)
Metabolite_summary.to_csv(OUTPUT_DIR / "counterfactual_summary_by_Metabolite.csv")

# 时间点汇总
time_summary = counterfactual_df.groupby(['site', 'time_numeric']).agg({
    'plasticity_index': ['mean', lambda x: x.std() / np.sqrt(len(x))],
    'is_divergent_adj': 'mean',
    'zero_consistency': 'mean'
}).round(4)
time_summary.columns = ['mean_plasticity', 'se_plasticity', 'prop_divergent', 'zero_consistency']
time_summary.to_csv(OUTPUT_DIR / "counterfactual_summary_by_time.csv")

# ==================== 主要发现 ====================
print("\n" + "=" * 60)
print("主要发现")
print("=" * 60)

print("\n老年人显著高于年轻人预期 (Top 5):")
top_positive = Metabolite_summary.nlargest(5, 'mean_plasticity')[
    ['mean_plasticity', 'n_divergent_adj', 'proportion_divergent']
]
print(top_positive)

print("\n老年人显著低于年轻人预期 (Top 5):")
top_negative = Metabolite_summary.nsmallest(5, 'mean_plasticity')[
    ['mean_plasticity', 'n_divergent_adj', 'proportion_divergent']
]
print(top_negative)

print("\n结果文件:")
print(f"  主要结果: counterfactual_predictions.csv")
print(f"  微生物汇总: counterfactual_summary_by_Metabolite.csv")
print(f"  时间汇总: counterfactual_summary_by_time.csv")
print("\n✓ 完成！")