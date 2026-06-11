"""
贝叶斯分层模型训练（稳健版）
增加重试机制和更好的错误处理
"""

import pandas as pd
import numpy as np
from pathlib import Path
import pickle
import warnings
warnings.filterwarnings('ignore')

# 贝叶斯建模库
import bambi as bmb
import pymc as pm
import arviz as az
from joblib import Parallel, delayed
import time
import os

# **关键：设置PyTensor避免文件锁冲突**
os.environ['PYTENSOR_FLAGS'] = 'base_compiledir=/tmp/pytensor_bayesian_cache'

# ==================== 配置 ====================
INPUT_DIR = "data/processed"
OUTPUT_DIR = "data/processed"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

N_JOBS = 12  # 减少并行数避免资源竞争
N_CHAINS = 2
N_SAMPLES = 1000
N_WARMUP = 400
MAX_RETRIES = 2  # 失败后重试次数

print("=" * 60)
print("开始训练年轻人专属贝叶斯模型（稳健版）")
print("=" * 60)

# ==================== 加载数据 ====================
all_data = pd.read_csv(INPUT_DIR / "preprocessed_data_for_counterfactual.csv")
young_only_data = all_data[all_data['age_group'] == 'youth'].copy()

print(f"\n年轻人数据统计:")
print(f"  观测点数: {len(young_only_data):,}")
print(f"  个体数: {young_only_data['subject_id'].nunique()}")
print(f"  微生物数: {young_only_data['microbe'].nunique()}")

combinations = young_only_data[['site', 'microbe']].drop_duplicates().reset_index(drop=True)
print(f"\n需要拟合 {len(combinations)} 个年轻人专属模型")

# ==================== 带重试的拟合函数 ====================
def fit_young_only_model_with_retry(idx, site, microbe, young_data, output_dir, max_retries=2):
    """带重试机制的模型拟合"""
    
    for attempt in range(max_retries + 1):
        try:
            result = fit_young_only_model_core(idx, site, microbe, young_data, output_dir, attempt)
            if result['success']:
                return result
            # 如果是数据问题（非技术错误），不重试
            if result.get('failure_reason') in ['insufficient_data', 'excessive_zeros']:
                return result
        except Exception as e:
            if attempt == max_retries:
                return {
                    'index': idx,
                    'site': site,
                    'microbe': microbe,
                    'success': False,
                    'error': f"Failed after {max_retries+1} attempts: {str(e)}",
                    'failure_reason': 'max_retries_exceeded'
                }
            time.sleep(1)  # 短暂等待后重试
    
    return result

def fit_young_only_model_core(idx, site, microbe, young_data, output_dir, attempt=0):
    """核心拟合逻辑"""
    
    # 筛选数据
    model_data = young_data[
        (young_data['site'] == site) & 
        (young_data['microbe'] == microbe)
    ].copy()
    
    # 数据检查
    if len(model_data) < 10:
        return {
            'index': idx,
            'site': site,
            'microbe': microbe,
            'success': False,
            'n_obs': len(model_data),
            'error': 'Insufficient data (<10 observations)',
            'failure_reason': 'insufficient_data'
        }
    
    zero_proportion = (model_data['log_abundance'] < -15).mean()
    if zero_proportion > 0.95:
        return {
            'index': idx,
            'site': site,
            'microbe': microbe,
            'success': False,
            'n_obs': len(model_data),
            'zero_proportion': float(zero_proportion),
            'error': f'Excessive zeros ({zero_proportion*100:.1f}%)',
            'failure_reason': 'excessive_zeros'
        }
    
    # 标准化时间
    time_mean = model_data['time_numeric'].mean()
    time_std = model_data['time_numeric'].std()
    model_data['time_scaled'] = (model_data['time_numeric'] - time_mean) / time_std
    
    # 构建模型
    formula = "log_abundance ~ bs(time_scaled, df=3) + (1 + time_scaled | subject_id)"
    
    # **简化版：使用默认先验（更稳定）**
    model = bmb.Model(
        formula,
        data=model_data,
        family='gaussian'
    )
    
    # **只为关键参数设置弱信息先验**
    model.set_priors({
        "Intercept": bmb.Prior("Normal", mu=0, sigma=5),
        "sigma": bmb.Prior("Exponential", lam=1)
    })
    
    # 拟合
    idata = model.fit(
        draws=N_SAMPLES,
        tune=N_WARMUP,
        chains=N_CHAINS,
        cores=1,
        random_seed=123 + idx + attempt * 1000,  # 重试时使用不同种子
        target_accept=0.90,
        max_treedepth=12,
        progressbar=False,
        return_inferencedata=True
    )
    
    # 收敛性检查
    rhat = az.rhat(idata)
    max_rhat = float(rhat.to_array().max().values)
    convergence_status = 'good' if max_rhat < 1.05 else 'acceptable' if max_rhat < 1.1 else 'poor'
    
    # 保存
    save_path = Path(output_dir) / f"idata_{idx:06d}_{site}_{microbe}.nc"
    idata.to_netcdf(save_path)
    
    return {
        'index': idx,
        'site': site,
        'microbe': microbe,
        'success': True,
        'n_obs': len(model_data),
        'max_rhat': max_rhat,
        'convergence_status': convergence_status,
        'zero_proportion': float(zero_proportion),
        'idata_path': str(save_path),
        'formula': formula,
        'time_mean': float(time_mean),
        'time_std': float(time_std),
        'n_subjects': model_data['subject_id'].nunique(),
        'attempt': attempt,
        'error': None,
        'failure_reason': None
    }

# ==================== 并行拟合 ====================
print(f"\n开始并行拟合 (使用 {N_JOBS} 个进程，最多重试{MAX_RETRIES}次)...")
start_time = time.time()

results = Parallel(n_jobs=N_JOBS, verbose=10, timeout=None)(
    delayed(fit_young_only_model_with_retry)(
        idx, 
        row['site'], 
        row['microbe'], 
        young_only_data,
        OUTPUT_DIR,
        MAX_RETRIES
    )
    for idx, row in combinations.iterrows()
)

elapsed = time.time() - start_time

# ==================== 结果汇总 ====================
successful_models = [r for r in results if r['success']]
failed_models = [r for r in results if not r['success']]

print("\n" + "=" * 60)
print("模型拟合完成")
print("=" * 60)
print(f"总用时: {elapsed/60:.2f} 分钟")
print(f"平均每个模型: {elapsed/len(results):.2f} 秒")
print(f"成功: {len(successful_models)}/{len(results)} ({len(successful_models)/len(results)*100:.1f}%)")
print(f"失败: {len(failed_models)}")

# 重试统计
if successful_models:
    retry_counts = [r.get('attempt', 0) for r in successful_models]
    print(f"\n重试统计（成功模型）:")
    print(f"  首次成功: {retry_counts.count(0)}")
    if max(retry_counts) > 0:
        print(f"  重试1次后成功: {retry_counts.count(1)}")
        print(f"  重试2次后成功: {retry_counts.count(2)}")

# 失败原因
if failed_models:
    from collections import Counter
    reasons = Counter(r.get('failure_reason', 'unknown') for r in failed_models)
    print(f"\n失败原因分析:")
    for reason, count in reasons.most_common():
        print(f"  {reason}: {count}")

# 收敛性统计
if successful_models:
    from collections import Counter
    conv_counts = Counter(r['convergence_status'] for r in successful_models)
    print(f"\n收敛性统计:")
    for status, count in conv_counts.items():
        print(f"  {status}: {count}")
    
    zero_props = [r['zero_proportion'] for r in successful_models]
    print(f"\n零值比例统计:")
    print(f"  平均: {np.mean(zero_props)*100:.2f}%")
    print(f"  中位数: {np.median(zero_props)*100:.2f}%")
    print(f"  范围: {np.min(zero_props)*100:.2f}% - {np.max(zero_props)*100:.2f}%")

# ==================== 保存结果 ====================
print("\n保存模型元数据...")

with open(OUTPUT_DIR / "young_only_models_metadata.pkl", 'wb') as f:
    pickle.dump(results, f, protocol=pickle.HIGHEST_PROTOCOL)

summary_df = pd.DataFrame([
    {
        'index': r['index'],
        'site': r['site'],
        'microbe': r['microbe'],
        'n_obs': r['n_obs'],
        'max_rhat': r.get('max_rhat', np.nan),
        'convergence_status': r.get('convergence_status', 'failed'),
        'success': r['success'],
        'attempt': r.get('attempt', 0),
        'idata_path': r.get('idata_path', ''),
        'failure_reason': r.get('failure_reason', ''),
        'error': str(r.get('error', ''))[:200]
    }
    for r in results
])
summary_df.to_csv(OUTPUT_DIR / "young_models_summary.csv", index=False)

if failed_models:
    failed_df = pd.DataFrame([
        {
            'site': m['site'], 
            'microbe': m['microbe'], 
            'failure_reason': m.get('failure_reason', 'unknown'),
            'error': str(m['error'])[:500]
        }
        for m in failed_models
    ])
    failed_df.to_csv(OUTPUT_DIR / "young_models_failed.csv", index=False)

print(f"\n✅ 完成！文件已保存至: {OUTPUT_DIR}")
print(f"   成功率: {len(successful_models)/len(results)*100:.1f}%")
print(f"   收敛良好: {sum(1 for r in successful_models if r['convergence_status']=='good')}")