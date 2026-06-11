"""
细菌反事实预测（Microbe 级 FDR + 统一 proportion_divergent = 0.3）
最终分类：time-differential / elderly-active / youth-active / delayed-response
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
from statsmodels.stats.multitest import multipletests

# ==================== 配置 ====================
INPUT_DIR = "data/processed"
MODEL_DIR = "scripts/03_bayesian_model/model"
OUTPUT_DIR = "data/processed"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

N_JOBS = 16
N_POSTERIOR_SAMPLES = 500

# ===== 统一分析阈值（与代谢物一致）=====
PROPORTION_THRESHOLD = 0.3
DIRECTION_CONSISTENCY = 0.7
EFFECT_SIZE = 0.3
EARLY_PROP_THRESHOLD = 0.2
LATE_PROP_THRESHOLD = 0.5

print("=" * 60)
print("细菌反事实预测（Microbe 级 FDR + 统一 0.3 阈值）")
print("=" * 60)

# ==================== 加载数据 ====================
all_data = pd.read_csv(INPUT_DIR / "preprocessed_data_for_counterfactual.csv")

with open(MODEL_DIR / "young_only_models_metadata.pkl", 'rb') as f:
    young_models = pickle.load(f)

successful_models = [m for m in young_models if m['success']]
old_data = all_data[all_data['age_group'] == 'old'].copy()

# ==================== 手动预测函数 ====================
def manual_predict_from_posterior(model_metadata, old_subjects_data):

    site = model_metadata['site']
    microbe = model_metadata['microbe']

    old_subset = old_subjects_data[
        (old_subjects_data['site'] == site) &
        (old_subjects_data['microbe'] == microbe)
    ].copy()

    if len(old_subset) == 0:
        return None

    try:
        idata = az.from_netcdf(model_metadata['idata_path'])
        posterior = idata.posterior
        available_params = list(posterior.data_vars.keys())

        # 标准化时间
        time_mean = model_metadata['time_mean']
        time_std = model_metadata['time_std']
        old_subset['time_scaled'] = (old_subset['time_numeric'] - time_mean) / time_std

        X_spline = dmatrix("bs(time_scaled, df=3) - 1",
                           data=old_subset, return_type='dataframe')

        intercept_samples = posterior['Intercept'].values.flatten()
        spline_data = posterior['bs(time_scaled, df=3)'].values

        spline_coef_samples = np.array([
            spline_data[..., i].flatten() for i in range(3)
        ])

        sigma_samples = posterior['sigma'].values.flatten()

        n_samples = min(N_POSTERIOR_SAMPLES, len(intercept_samples))
        intercept_samples = intercept_samples[:n_samples]
        spline_coef_samples = spline_coef_samples[:, :n_samples]
        sigma_samples = sigma_samples[:n_samples]

        n_obs = len(old_subset)
        predictions = np.zeros((n_samples, n_obs))

        for s in range(n_samples):
            mu = intercept_samples[s] + X_spline.values @ spline_coef_samples[:, s]
            predictions[s, :] = np.random.normal(mu, sigma_samples[s], n_obs)

        counter_mean = predictions.mean(axis=0)
        counter_lower = np.percentile(predictions, 2.5, axis=0)
        counter_upper = np.percentile(predictions, 97.5, axis=0)

        results = pd.DataFrame({
            'site': site,
            'microbe': microbe,
            'time_numeric': old_subset['time_numeric'].values,
            'observed_log_abundance': old_subset['log_abundance'].values,
            'counterfactual_mean': counter_mean,
            'counterfactual_lower': counter_lower,
            'counterfactual_upper': counter_upper
        })

        results['plasticity_index'] = results['observed_log_abundance'] - results['counterfactual_mean']
        results['is_divergent'] = (
            (results['observed_log_abundance'] < counter_lower) |
            (results['observed_log_abundance'] > counter_upper)
        )

        results['divergence_direction'] = 'none'
        results.loc[results['observed_log_abundance'] > counter_upper, 'divergence_direction'] = 'upper'
        results.loc[results['observed_log_abundance'] < counter_lower, 'divergence_direction'] = 'lower'

        results['p_value'] = [
            2 * min(
                np.mean(predictions[:, i] >= results.iloc[i]['observed_log_abundance']),
                np.mean(predictions[:, i] <= results.iloc[i]['observed_log_abundance'])
            ) for i in range(n_obs)
        ]

        return results

    except:
        return None

# ==================== 并行预测 ====================
results = Parallel(n_jobs=N_JOBS)(
    delayed(manual_predict_from_posterior)(m, old_data)
    for m in successful_models
)

counterfactual_df = pd.concat([r for r in results if r is not None], ignore_index=True)

# ==================== Microbe 级 FDR ====================
microbe_p = counterfactual_df.groupby(['site','microbe'])['p_value'].min().reset_index()
_, p_adj, _, _ = multipletests(microbe_p['p_value'], method='fdr_bh')
microbe_p['p_adj_microbe'] = p_adj
microbe_p['is_divergent_microbe'] = p_adj < 0.05

counterfactual_df = counterfactual_df.merge(
    microbe_p, on=['site','microbe'], how='left'
)

# ==================== 汇总 & 分类 ====================
summary = counterfactual_df.groupby(['site','microbe']).agg(
    n_timepoints = ('time_numeric','count'),
    mean_plasticity=('plasticity_index','mean'),
    n_divergent=('is_divergent','sum'),
    dirs=('divergence_direction',list),
    p_adj_microbe=('p_adj_microbe','first'),
    is_divergent_microbe=('is_divergent_microbe','first')
).reset_index()

summary['proportion_divergent'] = summary['n_divergent'] / summary['n_timepoints']

# 方向一致性
def direction_stats(lst):
    d = [x for x in lst if x!='none']
    if len(d)==0: return 0,0
    return d.count('upper')/len(d), d.count('lower')/len(d)

summary[['p_up','p_down']] = summary.apply(
    lambda r: pd.Series(direction_stats(r['dirs'])), axis=1
)

summary['is_time_differential'] = summary['proportion_divergent'] >= PROPORTION_THRESHOLD

summary['is_elderly_active'] = (
    summary['is_divergent_microbe'] &
    summary['is_time_differential'] &
    (summary['mean_plasticity'] > 0) &
    (summary['p_up'] >= DIRECTION_CONSISTENCY) &
    (summary['mean_plasticity'].abs() >= EFFECT_SIZE)
)

summary['is_youth_active'] = (
    summary['is_divergent_microbe'] &
    summary['is_time_differential'] &
    (summary['mean_plasticity'] < 0) &
    (summary['p_down'] >= DIRECTION_CONSISTENCY) &
    (summary['mean_plasticity'].abs() >= EFFECT_SIZE)
)

# 延迟响应
phase = counterfactual_df.copy()
q1, q2 = phase['time_numeric'].quantile([1/3, 2/3])

phase['phase'] = 'mid'
phase.loc[phase['time_numeric']<=q1,'phase']='early'
phase.loc[phase['time_numeric']>q2,'phase']='late'

phase_prop = phase.groupby(['site','microbe','phase'])['is_divergent'].mean().unstack()
summary = summary.merge(phase_prop, on=['site','microbe'], how='left')

summary['is_delayed_response'] = (
    (summary['early'] < EARLY_PROP_THRESHOLD) &
    (summary['late'] >= LATE_PROP_THRESHOLD)
)

# ==================== 保存结果 ====================
counterfactual_df.to_csv(
    OUTPUT_DIR / "counterfactual_predictions_with_microbeFDR.csv", index=False)

summary.to_csv(
    OUTPUT_DIR / "counterfactual_summary_by_microbe_unified0.3.csv", index=False)

summary.loc[summary['is_time_differential'],
            ['site','microbe']].to_csv(
    OUTPUT_DIR / "time_differential_microbes_unified0.3.csv", index=False)

summary.loc[summary['is_elderly_active'],
            ['site','microbe']].to_csv(
    OUTPUT_DIR / "elderly_active_microbes_unified0.3.csv", index=False)

summary.loc[summary['is_youth_active'],
            ['site','microbe']].to_csv(
    OUTPUT_DIR / "youth_active_microbes_unified0.3.csv", index=False)

summary.loc[summary['is_delayed_response'],
            ['site','microbe']].to_csv(
    OUTPUT_DIR / "delayed_response_microbes_unified0.3.csv", index=False)

print("\n✅ 细菌反事实分析 + 统一阈值分类 完成！")
