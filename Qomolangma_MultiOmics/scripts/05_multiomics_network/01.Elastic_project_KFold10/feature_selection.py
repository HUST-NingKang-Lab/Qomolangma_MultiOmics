import os
import pandas as pd
import numpy as np
import logging
from typing import Tuple, List, Dict, Optional
from sklearn.linear_model import ElasticNet
from sklearn.model_selection import KFold, GridSearchCV
from sklearn.preprocessing import StandardScaler
from utils import create_directories, load_cache, save_cache

def preprocess_data(X: pd.DataFrame, y: pd.Series) -> Tuple[pd.DataFrame, pd.Series, StandardScaler]:
    """数据预处理：清理和标准化"""
    # 移除缺失值
    valid_indices = ~(X.isnull().any(axis=1) | y.isnull())
    X_clean = X.loc[valid_indices].copy()
    y_clean = y.loc[valid_indices].copy()
    
    # 移除方差为0的特征
    var_threshold = 1e-8
    feature_vars = X_clean.var()
    valid_features = feature_vars[feature_vars > var_threshold].index
    X_clean = X_clean[valid_features]
    
    # 标准化
    scaler = StandardScaler()
    X_scaled = pd.DataFrame(
        scaler.fit_transform(X_clean),
        columns=X_clean.columns,
        index=X_clean.index
    )
    
    return X_scaled, y_clean, scaler

def elasticnet_feature_selection_cv(X: pd.DataFrame, y: pd.Series, 
                                   cv_folds: int = 10, random_state: int = 42) -> Tuple[List, List, List]:
    """
    使用交叉验证进行ElasticNet特征选择
    """
    # 优化的参数网格
    alphas = np.logspace(-4, 1, 20)
    l1_ratios = np.linspace(0.1, 0.9, 5)
    
    # 存储每折选择的特征
    fold_selected_features = []
    fold_weights = []
    fold_best_params = []
    fold_scores = []
    
    # 交叉验证
    kf = KFold(n_splits=cv_folds, shuffle=True, random_state=random_state)
    
    for fold, (train_idx, val_idx) in enumerate(kf.split(X)):
        try:
            X_train_fold = X.iloc[train_idx]
            y_train_fold = y.iloc[train_idx]
            X_val_fold = X.iloc[val_idx]
            y_val_fold = y.iloc[val_idx]
            
            # 预处理
            X_train_processed, y_train_processed, scaler = preprocess_data(X_train_fold, y_train_fold)
            X_val_processed = pd.DataFrame(
                scaler.transform(X_val_fold[X_train_processed.columns]),
                columns=X_train_processed.columns,
                index=X_val_fold.index
            )
            
            # 检查处理后的数据
            if X_train_processed.empty or len(X_train_processed) < 5:
                logging.warning(f"第{fold+1}折: 有效数据太少，跳过")
                continue
            
            # ElasticNet网格搜索
            elastic_net = ElasticNet(
                max_iter=10000,
                tol=1e-4,
                random_state=random_state,
                selection='cyclic'
            )
            
            param_grid = {
                'alpha': alphas,
                'l1_ratio': l1_ratios
            }
            
            grid_search = GridSearchCV(
                elastic_net, param_grid, cv=5,
                scoring='r2', n_jobs=-1,
                error_score='raise'
            )
            
            grid_search.fit(X_train_processed, y_train_processed)
            
            # 获取最佳参数和模型
            best_params = grid_search.best_params_
            best_model = grid_search.best_estimator_
            
            # 验证集性能
            val_score = best_model.score(X_val_processed, y_val_fold)
            
            fold_best_params.append(best_params)
            fold_scores.append(val_score)
            
            # 记录非零权重的特征
            feature_mask = np.abs(best_model.coef_) > 1e-6
            selected_features = X_train_processed.columns[feature_mask].tolist()
            fold_selected_features.append(selected_features)
            
            # 记录权重
            weights_dict = dict(zip(X_train_processed.columns, best_model.coef_))
            fold_weights.append(weights_dict)
            
            logging.debug(f"第{fold+1}折完成，选择特征数: {len(selected_features)}, 验证R²: {val_score:.4f}")
            
        except Exception as e:
            logging.warning(f"第{fold+1}折处理失败: {e}")
            continue
    
    if len(fold_selected_features) == 0:
        logging.error("所有折都处理失败")
        return [], [], []
    
    logging.info(f"成功完成 {len(fold_selected_features)}/{cv_folds} 折，平均验证R²: {np.mean(fold_scores):.4f}")
    
    return fold_selected_features, fold_weights, fold_best_params

def get_stable_features(fold_selected_features: List[List[str]], 
                       feature_names: List[str], threshold: float = 0.8) -> Tuple[List[str], Dict[str, int]]:
    """获取稳定特征"""
    feature_counts = {}
    total_folds = len(fold_selected_features)
    
    # 统计每个特征在多少折中被选择
    for features in fold_selected_features:
        for feature in features:
            feature_counts[feature] = feature_counts.get(feature, 0) + 1
    
    # 选择稳定特征
    min_count = threshold * total_folds
    stable_features = [
        feature for feature, count in feature_counts.items() 
        if count >= min_count
    ]
    
    logging.info(f"稳定特征阈值: {threshold:.1%}, 最小出现次数: {min_count:.1f}")
    logging.info(f"稳定特征数量: {len(stable_features)}")
    
    return stable_features, feature_counts

def calculate_average_weights(fold_weights: List[Dict], stable_features: List[str]) -> Dict[str, float]:
    """计算稳定特征的平均权重"""
    average_weights = {}
    
    for feature in stable_features:
        weights = [weights_dict.get(feature, 0) for weights_dict in fold_weights]
        # 只计算非零权重的平均值
        non_zero_weights = [w for w in weights if abs(w) > 1e-8]
        if non_zero_weights:
            average_weights[feature] = np.mean(non_zero_weights)
        else:
            average_weights[feature] = 0
    
    return average_weights

def process_single_metabolite(X_train: pd.DataFrame, y_metabolite: pd.Series, 
                            metabolite_name: str, output_dir: str, figure_dir: str, 
                            input_filename: str) -> Optional[Dict]:
    """处理单个代谢物的特征筛选"""
    from visualization import create_visualizations
    logging.info(f"处理代谢物: {metabolite_name}")
    
    # 创建代谢物专用目录
    metabolite_data_dir = os.path.join(output_dir, input_filename, metabolite_name)
    metabolite_figure_dir = os.path.join(figure_dir, input_filename, metabolite_name)
    create_directories([metabolite_data_dir, metabolite_figure_dir])
    
    # 检查缓存
    cache_path = os.path.join(metabolite_data_dir, "feature_selection_cache.pkl")
    cached_results = load_cache(cache_path)
    
    if cached_results is None:
        # 数据质量检查
        valid_ratio = (~y_metabolite.isnull()).sum() / len(y_metabolite)
        if valid_ratio < 0.5:
            logging.warning(f"代谢物 {metabolite_name} 有效数据比例过低: {valid_ratio:.2%}")
            return None
        
        # 进行特征选择
        try:
            fold_selected_features, fold_weights, fold_best_params = elasticnet_feature_selection_cv(
                X_train, y_metabolite
            )
            
            if not fold_selected_features:
                logging.warning(f"代谢物 {metabolite_name} 特征选择失败")
                return None
            
            # 获取稳定特征
            stable_features, feature_counts = get_stable_features(
                fold_selected_features, X_train.columns.tolist()
            )
            
            # 计算平均权重
            average_weights = calculate_average_weights(fold_weights, stable_features)
            
            # 缓存结果
            cached_results = {
                'fold_selected_features': fold_selected_features,
                'fold_weights': fold_weights,
                'fold_best_params': fold_best_params,
                'stable_features': stable_features,
                'feature_counts': feature_counts,
                'average_weights': average_weights
            }
            save_cache(cached_results, cache_path)
            
        except Exception as e:
            logging.error(f"代谢物 {metabolite_name} 处理失败: {e}")
            return None
    else:
        logging.info(f"从缓存加载结果: {cache_path}")
        stable_features = cached_results['stable_features']
        feature_counts = cached_results['feature_counts']
        average_weights = cached_results['average_weights']
    
    if len(stable_features) == 0:
        logging.warning(f"代谢物 {metabolite_name} 没有找到稳定特征")
        return None
    
    # 保存结果
    try:
        # 保存稳定特征
        stable_features_df = pd.DataFrame({
            'Feature': stable_features,
            'Average_Weight': [average_weights[f] for f in stable_features],
            'Selection_Frequency': [feature_counts[f]/len(cached_results['fold_selected_features']) for f in stable_features]
        })
        stable_features_df = stable_features_df.sort_values('Selection_Frequency', ascending=False)
        stable_features_df.to_csv(
            os.path.join(metabolite_data_dir, 'stable_features.csv'), index=False
        )
        
        # 保存所有特征的选择频率
        total_folds = len(cached_results['fold_selected_features'])
        all_features_freq = pd.DataFrame([
            {'Feature': feature, 'Selection_Count': count, 'Selection_Frequency': count/total_folds}
            for feature, count in feature_counts.items()
        ])
        all_features_freq = all_features_freq.sort_values('Selection_Count', ascending=False)
        all_features_freq.to_csv(
            os.path.join(metabolite_data_dir, 'all_features_frequency.csv'), index=False
        )
        
        # 保存最佳参数统计
        if cached_results['fold_best_params']:
            param_stats = pd.DataFrame(cached_results['fold_best_params'])
            param_summary = pd.DataFrame({
                'Parameter': ['alpha', 'l1_ratio'],
                'Mean': [param_stats['alpha'].mean(), param_stats['l1_ratio'].mean()],
                'Std': [param_stats['alpha'].std(), param_stats['l1_ratio'].std()],
                'Min': [param_stats['alpha'].min(), param_stats['l1_ratio'].min()],
                'Max': [param_stats['alpha'].max(), param_stats['l1_ratio'].max()]
            })
            param_summary.to_csv(
                os.path.join(metabolite_data_dir, 'parameter_summary.csv'), index=False
            )
        
        # 创建可视化
        create_visualizations(stable_features_df, metabolite_name, metabolite_figure_dir)
        
    except Exception as e:
        logging.error(f"结果保存失败 {metabolite_name}: {e}")
        return None
    
    return {
        'metabolite': metabolite_name,
        'stable_feature_count': len(stable_features),
        'stable_features': stable_features,
        'feature_counts': feature_counts,
        'average_weights': average_weights
    }