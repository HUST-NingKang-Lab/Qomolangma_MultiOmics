# step1: 创建运行脚本
input_csv="/public/home/zhuxue/Qomolangma_Human/raw_data/sample_path.csv"
output_script="/public/home/zhuxue/Qomolangma_Human/02.MAG/05.MAG_abs/code/coverm.sh"

mkdir -p "/public/home/zhuxue/Qomolangma_Human/02.MAG/05.MAG_abs/code/"

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
        out="/public/home/zhuxue/Qomolangma_Human/02.MAG/05.MAG_abs/${sample}"
        
        # 生成完整的coverm命令（确保在一行）
        echo "coverm genome -1 $fq1 -2 $fq2 -d /public/home/zhuxue/Qomolangma_Human/02.MAG/04.MAG_dRep/02.drep_out/dereplicated_genomes -x fa --min-read-percent-identity 95 --min-read-aligned-percent 75 --proper-pairs-only -m tpm -t 64 -o $out" >> "$output_script"
    else
        echo "警告: 样本 $sample 缺少配对的fastq文件" >&2
    fi
done

# step2： 拆分运行脚本
cd /public/home/zhuxue/Qomolangma_Human/02.MAG/05.MAG_abs/code
split -l 350 coverm.sh --numeric-suffixes=1 --suffix-length=2 split_code/code_

# step3: 检查运行失败的样本
find /public/home/zhuxue/Qomolangma_Human/02.MAG/05.MAG_abs \
-type f -empty \
-not -path "/public/home/zhuxue/Qomolangma_Human/02.MAG/05.MAG_abs/code/error_li.txt" \
-printf "%f\n" > /public/home/zhuxue/Qomolangma_Human/02.MAG/05.MAG_abs/code/error_li.txt

# step4: 筛选出报错的脚本并重新运行
cd /public/home/zhuxue/Qomolangma_Human/02.MAG/05.MAG_abs/code

# 逐行读取error_li.txt中的值
while IFS= read -r pattern; do
    # 跳过空行
    [ -z "$pattern" ] && continue
    # 使用grep查找包含pattern的行，并追加到输出文件
    grep "$pattern" coverm.sh >> error_coverm.sh
done < error_li.txt

# step6: 拆分重新运行的脚本
cd /public/home/zhuxue/Qomolangma_Human/02.MAG/05.MAG_abs/code
split -l 10 error_coverm.sh --numeric-suffixes=1 --suffix-length=2 split_code_2/code_

source /public/home/zhuxue/tools/miniconda3/bin/activate
conda activate coverm_env