# step1: 创建运行脚本
input_csv="/public/home/zhuxue/Qomolangma_Human/raw_data/sample_path.csv"
output_script="/public/home/zhuxue/Qomolangma_Human/02.MAG/05.MAG_abs/code/coverm.sh"

mkdir -p "/public/home/zhuxue/Qomolangma_Human/02.MAG/05.MAG_abs/code"

> "$output_script"

# 使用关联数组存储每个sample的路径
declare -A sample_fq1
declare -A sample_fq2

# 读取CSV文件并处理
while IFS=',' read -r sample path || [[ -n "$sample" ]]; do
    # 去除空格和回车符
    sample=$(echo "$sample" | tr -d '[:space:]')
    path=$(echo "$path" | tr -d '[:space:]')
    
    # 跳过空行和标题行
    if [[ -z "$sample" ]] || [[ "$sample" == "Sample" ]] || [[ "$sample" == "sample" ]]; then
        continue
    fi
    
    # 根据文件名判断是_1还是_2（支持.fastq.gz和.fq.gz）
    if [[ "$path" == *"_1.fastq.gz" ]] || [[ "$path" == *"_1.fq.gz" ]]; then
        sample_fq1["$sample"]="$path"
    elif [[ "$path" == *"_2.fastq.gz" ]] || [[ "$path" == *"_2.fq.gz" ]]; then
        sample_fq2["$sample"]="$path"
    fi
    
done < "$input_csv"

# 生成coverm命令
for sample in "${!sample_fq1[@]}"; do
    if [[ -n "${sample_fq1[$sample]}" ]] && [[ -n "${sample_fq2[$sample]}" ]]; then
        fq1="${sample_fq1[$sample]}"
        fq2="${sample_fq2[$sample]}"
        out="/public/home/zhuxue/Qomolangma_Human/02.MAG/05.MAG_abs/coverm_out/${sample}"
        
        # 生成完整的coverm命令（确保在一行）
        echo "coverm genome -1 $fq1 -2 $fq2 -d /public/home/zhuxue/Qomolangma_Human/02.MAG/04.MAG_dRep/02.drep_out/dereplicated_genomes -x fa --min-read-percent-identity 95 --min-read-aligned-percent 75 --min-covered-fraction 0 --proper-pairs-only --methods count covered_bases covered_fraction relative_abundance rpkm tpm -t 64 -o $out" >> "$output_script"

    else
        echo "警告: 样本 $sample 缺少配对的fastq文件" >&2
    fi
done

# step2： 拆分运行脚本
cd /public/home/zhuxue/Qomolangma_Human/02.MAG/05.MAG_abs/code
split -l 25 coverm.sh --numeric-suffixes=1 --suffix-length=2 split_code/code_

# step3: 获得TPM丰度表
import os
import pandas as pd
from pathlib import Path
# 定义文件路径
directory = "/public/home/zhuxue/Qomolangma_Human/02.MAG/05.MAG_abs/coverm_out"
# 初始化一个空的 DataFrame 用于存储合并结果
merged_df = None

# 遍历目录下的所有文件
for file_path in Path(directory).glob("*"):
    if file_path.is_file():
        # 读取文件，假设文件是制表符分隔的
        df = pd.read_csv(file_path, sep="\t")
        # 提取第一列（Genome）和最后一列（TPM）
        temp_df = df.iloc[:, [0, -1]].copy()
        # 将最后一列的列名改为文件名（不含路径和扩展名）
        filename = file_path.stem
        temp_df.columns = ["Genome", filename]
        # 如果是第一个文件，初始化 merged_df
        if merged_df is None:
            merged_df = temp_df
        else:
            # 合并数据，基于 Genome 列
            merged_df = merged_df.merge(temp_df, on="Genome", how="outer")

# 保存合并后的结果到 CSV 文件
output_file = "/public/home/zhuxue/Qomolangma_Human/02.MAG/05.MAG_abs/merged_coverm_data.csv"
merged_df.to_csv(output_file, index=False)

print(f"Merged data saved to {output_file}")