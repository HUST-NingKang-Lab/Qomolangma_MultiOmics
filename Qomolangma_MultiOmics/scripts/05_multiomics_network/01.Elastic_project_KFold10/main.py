import argparse
import sys
import logging
import os
from pathlib import Path
from tqdm import tqdm
import pandas as pd
from utils import setup_logging, create_directories
from data_processing import split_train_test, load_metadata, validate_data_quality
from feature_selection import process_single_metabolite
from visualization import create_summary_visualizations

def main():
    # 设置命令行参数
    parser = argparse.ArgumentParser(description='ElasticNet特征筛选脚本（优化版）')
    parser.add_argument('--input_file', '-i', 
                       default='data/processed/gut_clr_Lipid_lg_split_001.csv',
                       help='输入CSV文件路径')
    parser.add_argument('--metadata_file', '-m',
                       default='data/metadata/Matedata_Information.csv',
                       help='metadata文件路径')
    parser.add_argument('--output_dir', '-o',
                       default='01.ElasticNet_Kfold10',
                       help='输出目录')
    parser.add_argument('--figure_dir', '-f',
                       default='01.ElasticNet_Kfold10',
                       help='图片输出目录')
    parser.add_argument('--test_size', '-t', type=float, default=0.2,
                       help='测试集比例 (默认: 0.2)')
    parser.add_argument('--cv_folds', '-c', type=int, default=10,
                       help='交叉验证折数 (默认: 10)')
    parser.add_argument('--random_state', '-r', type=int, default=42,
                       help='随机种子 (默认: 42)')
    parser.add_argument('--threshold', type=float, default=0.8,
                       help='稳定特征选择阈值 (默认: 0.8)')
    
    args = parser.parse_args()
    
    # 设置日志
    log_dir = os.path.join(args.output_dir, 'logs')
    create_directories([log_dir])
    log_file = os.path.join(log_dir, 'elasticnet_feature_selection.log')
    setup_logging(log_file)
    
    # 输入验证
    if not os.path.exists(args.input_file):
        logging.error(f"输入文件不存在: {args.input_file}")
        sys.exit(1)
    
    if not os.path.exists(args.metadata_file):
        logging.error(f"metadata文件不存在: {args.metadata_file}")
        sys.exit(1)
    
    if not 0 < args.test_size < 1:
        logging.error(f"测试集比例必须在0和1之间: {args.test_size}")
        sys.exit(1)
    
    # 获取输入文件名
    input_filename = Path(args.input_file).stem
    
    # 创建输出目录
    train_test_dir = os.path.join(args.output_dir, '01.train_test_dataset')
    feature_selection_dir = os.path.join(args.output_dir, '02.ElasticNet特征筛选')
    figure_dir = os.path.join(args.figure_dir, '02.ElasticNet特征筛选')
    
    create_directories([train_test_dir, feature_selection_dir, figure_dir])
    
    logging.info("="*80)
    logging.info("ElasticNet特征筛选分析开始")
    logging.info("="*80)
    
    try:
        # 1. 划分训练集和测试集
        logging.info("1. 划分训练集和测试集...")
        train_path, test_path, train_df, test_df = split_train_test(
            args.input_file, 
            train_test_dir, 
            test_size=args.test_size, 
            random_state=args.random_state
        )
        
        # 2. 加载metadata
        logging.info("2. 加载metadata...")
        X_features, y_features, metadata = load_metadata(args.metadata_file, input_filename)

        # 提取存在的特征列
        existing_X_features = [f for f in X_features if f in train_df.columns]
        existing_y_features = [f for f in y_features if f in train_df.columns]

        # 提取训练数据
        X_train = train_df[existing_X_features]
        y_train = train_df[existing_y_features]
        
        # 3. 对每个代谢物进行特征筛选
        logging.info("3. 开始处理代谢物...")
        results = []
        feature_metabolite_counts = {}
        
        for metabolite in tqdm(existing_y_features, desc="Processing metabolites"):
            result = process_single_metabolite(
                X_train=X_train,
                y_metabolite=y_train[metabolite],
                metabolite_name=metabolite,
                output_dir=feature_selection_dir,
                figure_dir=figure_dir,
                input_filename=input_filename
            )
            
            if result is not None:
                results.append({
                    'Metabolite': result['metabolite'],
                    'Stable_Feature_Count': result['stable_feature_count']
                })
                
                # 统计特征出现在多少代谢物中
                for feature in result['stable_features']:
                    feature_metabolite_counts[feature] = feature_metabolite_counts.get(feature, 0) + 1
        
        if not results:
            logging.error("没有成功处理任何代谢物")
            sys.exit(1)
        
        # 4. 生成汇总结果
        logging.info("4. 生成汇总结果...")
        summary_df = pd.DataFrame(results)
        summary_df = summary_df.sort_values('Stable_Feature_Count', ascending=False)
        
        # 保存汇总结果
        summary_path = os.path.join(feature_selection_dir, f'{input_filename}_summary.csv')
        summary_df.to_csv(summary_path, index=False)
        logging.info(f"汇总结果已保存至: {summary_path}")
        
        # 处理特征-代谢物计数
        feature_metabolite_df = pd.DataFrame([
            {'Feature': feature, 'Metabolite_Count': count}
            for feature, count in feature_metabolite_counts.items()
        ])
        feature_metabolite_df = feature_metabolite_df.sort_values('Metabolite_Count', ascending=False)
        
        feature_counts_path = os.path.join(feature_selection_dir, f'{input_filename}_feature_metabolite_counts.csv')
        feature_metabolite_df.to_csv(feature_counts_path, index=False)
        logging.info(f"特征-代谢物计数已保存至: {feature_counts_path}")
        
        # 5. 创建汇总可视化
        logging.info("5. 创建汇总可视化...")
        create_summary_visualizations(summary_df, feature_metabolite_df, figure_dir)
        
        logging.info("="*80)
        logging.info("ElasticNet特征筛选分析完成")
        logging.info("="*80)
        
    except Exception as e:
        logging.error(f"程序执行失败: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()