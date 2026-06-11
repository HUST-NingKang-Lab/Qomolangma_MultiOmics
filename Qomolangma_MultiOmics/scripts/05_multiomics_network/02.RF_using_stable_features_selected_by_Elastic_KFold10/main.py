import os
import pandas as pd
import sys
import numpy as np
import logging
from pathlib import Path
from typing import Tuple, List, Dict, Optional
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold, GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
from scipy.stats import spearmanr
import pickle
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'SimHei']
plt.rcParams['axes.unicode_minus'] = False

def setup_logging(log_file: str) -> None:
    """设置日志配置"""
    log_format = '%(asctime)s - %(levelname)s - %(message)s'
    logging.basicConfig(
        level=logging.INFO,
        format=log_format,
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )

def create_directories(paths: List[str]) -> None:
    """创建必要的目录"""
    for path in paths:
        Path(path).mkdir(parents=True, exist_ok=True)
        logging.info(f"创建/确认目录: {path}")

def load_cache(cache_path: str) -> Optional[Dict]:
    """加载缓存文件"""
    if os.path.exists(cache_path):
        try:
            with open(cache_path, 'rb') as f:
                return pickle.load(f)
        except Exception as e:
            logging.warning(f"缓存文件加载失败 {cache_path}: {e}")
    return None

def save_cache(data: Dict, cache_path: str) -> None:
    """保存缓存文件"""
    try:
        with open(cache_path, 'wb') as f:
            pickle.dump(data, f)
        logging.info(f"缓存已保存至: {cache_path}")
    except Exception as e:
        logging.error(f"缓存保存失败 {cache_path}: {e}")

def preprocess_data(X: pd.DataFrame, y: pd.Series, scaler: Optional[StandardScaler] = None) -> Tuple[pd.DataFrame, pd.Series, StandardScaler]:
    """数据预处理：清理和标准化"""
    valid_indices = ~(X.isnull().any(axis=1) | y.isnull())
    X_clean = X.loc[valid_indices].copy()
    y_clean = y.loc[valid_indices].copy()
    
    var_threshold = 1e-8
    feature_vars = X_clean.var()
    valid_features = feature_vars[feature_vars > var_threshold].index
    X_clean = X_clean[valid_features]
    
    if scaler is None:
        scaler = StandardScaler()
        X_scaled = pd.DataFrame(
            scaler.fit_transform(X_clean),
            columns=X_clean.columns,
            index=X_clean.index
        )
    else:
        X_scaled = pd.DataFrame(
            scaler.transform(X_clean),
            columns=X_clean.columns,
            index=X_clean.index
        )
    
    return X_scaled, y_clean, scaler

def train_random_forest_model(X: pd.DataFrame, y: pd.Series, cv_folds: int = 10, 
                             random_state: int = 42) -> Tuple[RandomForestRegressor, Dict, List[Dict], StandardScaler]:
    """训练随机森林模型并进行交叉验证"""
    param_grid = {
        'n_estimators': [100, 200],
        'max_depth': [10, 20, None],
        'min_samples_split': [2, 5],
        'min_samples_leaf': [1, 2]
    }
    
    metrics = {'r2': [], 'mse': [], 'mae': [], 'spearman': []}
    feature_importances = []
    
    kf = KFold(n_splits=cv_folds, shuffle=True, random_state=random_state)
    
    for fold, (train_idx, val_idx) in enumerate(kf.split(X)):
        X_train_fold = X.iloc[train_idx]
        y_train_fold = y.iloc[train_idx]
        X_val_fold = X.iloc[val_idx]
        y_val_fold = y.iloc[val_idx]
        
        X_train_processed, y_train_processed, scaler = preprocess_data(X_train_fold, y_train_fold)
        X_val_processed = pd.DataFrame(
            scaler.transform(X_val_fold[X_train_processed.columns]),
            columns=X_train_processed.columns,
            index=X_val_fold.index
        )
        
        if X_train_processed.empty or len(X_train_processed) < 5:
            logging.warning(f"第{fold+1}折: 有效数据太少，跳过")
            continue
            
        rf = RandomForestRegressor(random_state=random_state)
        
        grid_search = GridSearchCV(
            rf, param_grid, cv=cv_folds, scoring='r2', n_jobs=-1, error_score='raise'
        )
        grid_search.fit(X_train_processed, y_train_processed)
        
        best_model = grid_search.best_estimator_
        y_pred = best_model.predict(X_val_processed)
        
        metrics['r2'].append(r2_score(y_val_fold, y_pred))
        metrics['mse'].append(mean_squared_error(y_val_fold, y_pred))
        metrics['mae'].append(mean_absolute_error(y_val_fold, y_pred))
        metrics['spearman'].append(spearmanr(y_val_fold, y_pred)[0])
        
        feature_importances.append(dict(zip(X_train_processed.columns, best_model.feature_importances_)))
    
    avg_metrics = {
        'r2': np.mean(metrics['r2']) if metrics['r2'] else np.nan,
        'mse': np.mean(metrics['mse']) if metrics['mse'] else np.nan,
        'mae': np.mean(metrics['mae']) if metrics['mae'] else np.nan,
        'spearman': np.mean(metrics['spearman']) if metrics['spearman'] else np.nan
    }
    
    return grid_search.best_estimator_, avg_metrics, feature_importances, scaler

def predict_and_evaluate(model: RandomForestRegressor, X: pd.DataFrame, y: pd.Series, 
                        scaler: StandardScaler) -> Dict:
    """对数据集进行预测并计算评价指标"""
    X_processed, y_processed, _ = preprocess_data(X, y, scaler)
    if X_processed.empty:
        return {'r2': np.nan, 'mse': np.nan, 'mae': np.nan, 'spearman': np.nan}
    
    y_pred = model.predict(X_processed)
    return {
        'r2': r2_score(y_processed, y_pred),
        'mse': mean_squared_error(y_processed, y_pred),
        'mae': mean_absolute_error(y_processed, y_pred),
        'spearman': spearmanr(y_processed, y_pred)[0]
    }

def calculate_average_importances(importances: List[Dict], features: List[str]) -> Dict:
    """计算特征的平均重要性"""
    avg_importances = {}
    for feature in features:
        imps = [imp_dict.get(feature, 0) for imp_dict in importances]
        non_zero_imps = [i for i in imps if i > 1e-8]
        avg_importances[feature] = np.mean(non_zero_imps) if non_zero_imps else 0
    return avg_importances

def process_metabolite_prediction(X_train: pd.DataFrame, y_train: pd.Series, 
                                 X_test: pd.DataFrame, y_test: pd.Series, 
                                 metabolite_name: str, stable_features: List[str],
                                 output_dir: str, figure_dir: str, input_filename: str) -> Dict:
    """处理单个代谢物的预测任务"""
    metabolite_data_dir = os.path.join(output_dir, input_filename, metabolite_name)
    metabolite_figure_dir = os.path.join(figure_dir, input_filename, metabolite_name)
    create_directories([metabolite_data_dir, metabolite_figure_dir])
    
    cache_path = os.path.join(metabolite_data_dir, 'prediction_results.pkl')
    cached_results = load_cache(cache_path)
    
    if cached_results is None:
        X_train_subset = X_train[stable_features]
        X_test_subset = X_test[stable_features]
        
        valid_ratio = (~y_train.isnull()).sum() / len(y_train)
        if valid_ratio < 0.5:
            logging.warning(f"代谢物 {metabolite_name} 有效数据比例过低: {valid_ratio:.2%}")
            return None
            
        model, cv_metrics, feature_importances, scaler = train_random_forest_model(X_train_subset, y_train)
        
        train_metrics = predict_and_evaluate(model, X_train_subset, y_train, scaler)
        test_metrics = predict_and_evaluate(model, X_test_subset, y_test, scaler)
        
        avg_importances = calculate_average_importances(feature_importances, stable_features)
        
        model_path = os.path.join(metabolite_data_dir, 'random_forest_model.pkl')
        with open(model_path, 'wb') as f:
            pickle.dump({'model': model, 'scaler': scaler}, f)
        logging.info(f"模型已保存至: {model_path}")
        
        results = {
            'cv_metrics': cv_metrics,
            'train_metrics': train_metrics,
            'test_metrics': test_metrics,
            'avg_importances': avg_importances
        }
        save_cache(results, cache_path)
        
        X_test_processed, y_test_processed, _ = preprocess_data(X_test_subset, y_test, scaler)
        if not X_test_processed.empty:
            y_pred = model.predict(X_test_processed)
            plt.figure(figsize=(10, 6))
            plt.scatter(y_test_processed, y_pred, alpha=0.5)
            plt.plot([y_test_processed.min(), y_test_processed.max()], 
                    [y_test_processed.min(), y_test_processed.max()], 'r--', lw=2)
            plt.xlabel('实际值')
            plt.ylabel('预测值')
            plt.title(f'{metabolite_name} - 测试集预测结果')
            plt.tight_layout()
            plt.savefig(os.path.join(metabolite_figure_dir, 'prediction_scatter.png'), dpi=300)
            plt.close()
            
            importance_df = pd.DataFrame({
                'Feature': list(avg_importances.keys()),
                'Importance': list(avg_importances.values())
            }).sort_values('Importance', ascending=False)
            
            plt.figure(figsize=(12, 8))
            top_n = min(20, len(importance_df))
            sns.barplot(data=importance_df.head(top_n), x='Importance', y='Feature')
            plt.title(f'{metabolite_name} - 特征重要性 (Top {top_n})')
            plt.xlabel('平均特征重要性')
            plt.tight_layout()
            plt.savefig(os.path.join(metabolite_figure_dir, 'feature_importance.png'), dpi=300)
            plt.close()
    else:
        logging.info(f"从缓存加载预测结果: {cache_path}")
        cv_metrics = cached_results['cv_metrics']
        train_metrics = cached_results['train_metrics']
        test_metrics = cached_results['test_metrics']
        avg_importances = cached_results['avg_importances']
    
    return {
        'metabolite': metabolite_name,
        'cv_metrics': cv_metrics,
        'train_metrics': train_metrics,
        'test_metrics': test_metrics,
        'avg_importances': avg_importances
    }

def load_metadata(metadata_path: str, input_filename: str) -> Tuple[List[str], List[str], pd.DataFrame]:
    """加载metadata文件并根据文件名提取特征"""
    logging.info(f"正在读取metadata文件: {metadata_path}")
    
    try:
        metadata = pd.read_csv(metadata_path)
        
        required_cols = ['Source', 'Group', 'Feature']
        missing_cols = [col for col in required_cols if col not in metadata.columns]
        if missing_cols:
            raise ValueError(f"metadata文件缺少必需列: {missing_cols}")
        
        y_source = None
        x_source = None
        
        if "metabolin" in input_filename.lower():
            y_source = "Metabolin_lg"
        elif "lipid" in input_filename.lower():
            y_source = "Lipid_lg"
        
        if "gut" in input_filename.lower():
            x_source = "gut"
        elif "oral" in input_filename.lower():
            x_source = "oral"
        elif "skin" in input_filename.lower():
            x_source = "skin"
        
        if y_source is None or x_source is None:
            raise ValueError(f"无法从文件名 '{input_filename}' 中提取有效的自变量或因变量信息")
        
        X_features = metadata[
            ((metadata['Source'] == x_source) | (metadata['Source'] == 'people')) & 
            (metadata['Group'] == 'X')
        ]['Feature'].tolist()
        
        y_features = metadata[
            (metadata['Source'] == y_source) & 
            (metadata['Group'] == 'Y')
        ]['Feature'].tolist()
        
        logging.info(f"自变量来源: {x_source}, 特征数: {len(X_features)}")
        logging.info(f"因变量来源: {y_source}, 特征数: {len(y_features)}")

        return X_features, y_features, metadata
        
    except Exception as e:
        logging.error(f"metadata加载失败: {e}")
        raise

def main(args):
    # 设置参数
    input_file = args.input_file
    metadata_file = args.metadata_file
    output_dir = args.output_dir
    figure_dir = args.figure_dir
    test_size = args.test_size
    random_state = args.random_state
    cv_folds = args.cv_folds
    
    # 设置日志
    log_dir = os.path.join(output_dir, 'logs')
    create_directories([log_dir])
    setup_logging(os.path.join(log_dir, 'random_forest_prediction.log'))
    
    # 创建输出目录
    create_directories([output_dir, figure_dir])
    
    logging.info("="*80)
    logging.info("随机森林预测分析开始")
    logging.info("="*80)
    
    try:
        # 加载特征筛选步骤中的训练集和测试集
        train_test_dir = args.train_test_dir
        input_filename = Path(input_file).stem
        train_path = os.path.join(train_test_dir, input_filename, f'{input_filename}_train.csv')
        test_path = os.path.join(train_test_dir, input_filename, f'{input_filename}_test.csv')
        
        if not os.path.exists(train_path) or not os.path.exists(test_path):
            raise FileNotFoundError(f"训练集或测试集文件不存在: {train_path}, {test_path}")
        
        train_df = pd.read_csv(train_path)
        test_df = pd.read_csv(test_path)
        logging.info(f"已加载训练集: {train_path}, 样本数: {len(train_df)}")
        logging.info(f"已加载测试集: {test_path}, 样本数: {len(test_df)}")
        
        # 加载metadata
        X_features, y_features, metadata = load_metadata(metadata_file, input_filename)
        
        # 提取存在的特征列
        existing_X_features = [f for f in X_features if f in train_df.columns]
        existing_y_features = [f for f in y_features if f in train_df.columns]
        
        X_train = train_df[existing_X_features]
        y_train = train_df[existing_y_features]
        X_test = test_df[existing_X_features]
        y_test = test_df[existing_y_features]
        
        # 处理每个代谢物
        results = []
        all_importances = []
        
        for metabolite in tqdm(existing_y_features, desc="Processing metabolites"):
            # 加载稳定特征
            feature_selection_dir = os.path.join(
                args.feature_selection_dir, input_filename, metabolite
            )
            stable_features_path = os.path.join(feature_selection_dir, 'stable_features.csv')
            
            if not os.path.exists(stable_features_path):
                logging.warning(f"代谢物 {metabolite} 的稳定特征文件不存在，跳过")
                continue
                
            stable_features_df = pd.read_csv(stable_features_path)
            stable_features = stable_features_df['Feature'].tolist()
            
            if not stable_features:
                logging.warning(f"代谢物 {metabolite} 没有稳定特征，跳过")
                continue
                
            result = process_metabolite_prediction(
                X_train, y_train[metabolite], X_test, y_test[metabolite],
                metabolite, stable_features, output_dir, figure_dir, input_filename
            )
            
            if result is not None:
                results.append({
                    'Metabolite': metabolite,
                    'CV_R2': result['cv_metrics']['r2'],
                    'CV_MSE': result['cv_metrics']['mse'],
                    'CV_MAE': result['cv_metrics']['mae'],
                    'CV_Spearman': result['cv_metrics']['spearman'],
                    'Train_R2': result['train_metrics']['r2'],
                    'Train_MSE': result['train_metrics']['mse'],
                    'Train_MAE': result['train_metrics']['mae'],
                    'Train_Spearman': result['train_metrics']['spearman'],
                    'Test_R2': result['test_metrics']['r2'],
                    'Test_MSE': result['test_metrics']['mse'],
                    'Test_MAE': result['test_metrics']['mae'],
                    'Test_Spearman': result['test_metrics']['spearman']
                })
                
                for feature, importance in result['avg_importances'].items():
                    all_importances.append({
                        'Metabolite': metabolite,
                        'Feature': feature,
                        'Average_Importance': importance
                    })
        
        # 保存统计表
        metrics_df = pd.DataFrame(results)
        metrics_df.to_csv(os.path.join(output_dir, f'{input_filename}_prediction_metrics.csv'), index=False)
        logging.info(f"预测评价指标已保存至: {os.path.join(output_dir, f'{input_filename}_prediction_metrics.csv')}")
        
        importances_df = pd.DataFrame(all_importances)
        importances_df.to_csv(os.path.join(output_dir, f'{input_filename}_feature_importances.csv'), index=False)
        logging.info(f"特征重要性已保存至: {os.path.join(output_dir, f'{input_filename}_feature_importances.csv')}")
        
        # 创建汇总可视化
        plt.figure(figsize=(12, 6))
        metrics_melted = pd.melt(
            metrics_df, 
            id_vars=['Metabolite'], 
            value_vars=['CV_R2', 'Train_R2', 'Test_R2', 'CV_Spearman', 'Train_Spearman', 'Test_Spearman'],
            var_name='Metric', 
            value_name='Value'
        )
        sns.boxplot(x='Metric', y='Value', data=metrics_melted)
        plt.title('R²和Spearman相关系数分布 - 交叉验证 vs 训练集 vs 测试集')
        plt.savefig(os.path.join(figure_dir, 'metrics_distribution.png'), dpi=300)
        plt.close()
        
        logging.info("="*80)
        logging.info("随机森林预测分析完成")
        logging.info("="*80)
        
    except Exception as e:
        logging.error(f"程序执行失败: {e}")
        sys.exit(1)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Random Forest Prediction Script')
    parser.add_argument('--input_file', type=str, 
                        default="data/processed/gut_clr_Lipid_lg_split_001.csv",
                        help='Path to input CSV file')
    parser.add_argument('--metadata_file', type=str,
                        default="data/metadata/Matedata_Information.csv",
                        help='Path to metadata CSV file')
    parser.add_argument('--output_dir', type=str,
                        default="03.RF_output",
                        help='Output directory for results')
    parser.add_argument('--figure_dir', type=str,
                        default="03.RF_output",
                        help='Output directory for figures')
    parser.add_argument('--train_test_dir', type=str,
                        default="01.train_test_dataset",
                        help='Directory containing train and test datasets')
    parser.add_argument('--feature_selection_dir', type=str,
                        default="01.ElasticNet特征筛选",
                        help='Directory containing feature selection results')
    parser.add_argument('--test_size', type=float, default=0.2, 
                        help='Test set size ratio (default: 0.2)')
    parser.add_argument('--random_state', type=int, default=42, 
                        help='Random seed (default: 42)')
    parser.add_argument('--cv_folds', type=int, default=10, 
                        help='Number of cross-validation folds (default: 10)')
    
    args = parser.parse_args()
    main(args)