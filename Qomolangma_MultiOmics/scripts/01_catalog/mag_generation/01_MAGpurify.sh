# step1：激活magpurify环境
source /public/home/zhuxue/miniconda3/bin/activate
conda activate magpurify-2.1.2

# step2：去除MAG中所有非标准的碱基，生成cleaned.fa
awk '/^>/ {print; next} {gsub(/[^ATCGatcg]/, ""); print}' QML0000-240606-02B_bin.10.strict.fa > cleaned.fa

## 执行三次，分别对gut/oral/skin
src_dir="/public/home/zhuxue/Qomolangma_Human/02.MAG/01.MAG_raw/gut"
dest_dir="/public/home/zhuxue/Qomolangma_Human/02.MAG/02.MAG_rm_contamination/01.MAG_rm_nonstandard/gut"

# 遍历源目录下所有的 .fa 文件
for file in "$src_dir"/*.fa; do
    filename=$(basename "$file")  # 提取文件名
    dest_file="$dest_dir/$filename"  # 设置目标文件路径

    # 处理每个 .fa 文件并输出到目标目录
    awk '/^>/ {print; next} {gsub(/[^ATCGatcg]/, ""); print}' "$file" > "$dest_file"
done

# step3：依次运行前5个模块，生成的结果文件放在magpurify_output
# 线程：32（注意3、4模块不能选择线程）

#!/bin/bash

# 输入参数
INPUT_FASTA="$1"
OUTPUT_DIR="$2"
DB_PATH="/public/home/zhuxue/Qomolangma_Human/database/MAGpurify-db/MAGpurify-db-v1.0"
THREADS=32

# 检查输入文件是否存在
if [ ! -f "$INPUT_FASTA" ]; then
    echo "错误: 输入文件 '$INPUT_FASTA' 不存在"
    exit 1
fi

# 获取基本文件名
FILENAME=$(basename -- "$INPUT_FASTA")
BASENAME="${FILENAME%.*}"

# 创建输出目录
mkdir -p "$OUTPUT_DIR" || { echo "无法创建输出目录: $OUTPUT_DIR"; exit 1; }

# 创建日志文件
LOG_FILE="$OUTPUT_DIR/${BASENAME}_process.log"
echo "MAGpurify处理日志 - $(date)" > "$LOG_FILE"
echo "输入文件: $INPUT_FASTA" >> "$LOG_FILE"

# 函数：执行命令并记录
run_command() {
    local step="$1"
    local cmd="$2"
    echo "" >> "$LOG_FILE"
    echo "===== 步骤 $step: $(date) =====" >> "$LOG_FILE"
    echo "执行: $cmd" >> "$LOG_FILE"
    
    # 执行命令并捕获输出
    eval "$cmd" >> "$LOG_FILE" 2>&1
    
    # 检查退出状态
    if [ $? -ne 0 ]; then
        echo "错误: 步骤 $step 执行失败" | tee -a "$LOG_FILE"
        exit 1
    else
        echo "步骤 $step 完成: $(date)" >> "$LOG_FILE"
    fi
}

# 按顺序执行命令
echo "开始处理 $INPUT_FASTA..."
{
    run_command 1 "magpurify phylo-markers --db '$DB_PATH' --threads $THREADS '$INPUT_FASTA' '$OUTPUT_DIR/'"
    run_command 2 "magpurify clade-markers --db '$DB_PATH' --threads $THREADS '$INPUT_FASTA' '$OUTPUT_DIR/'"
    run_command 3 "magpurify tetra-freq '$INPUT_FASTA' '$OUTPUT_DIR/'"
    run_command 4 "magpurify gc-content '$INPUT_FASTA' '$OUTPUT_DIR/'"
    run_command 5 "magpurify known-contam --db '$DB_PATH' --threads $THREADS '$INPUT_FASTA' '$OUTPUT_DIR/'"
} | tee -a "$LOG_FILE"

# 最终检查和处理
if [ $? -eq 0 ]; then
    echo "" >> "$LOG_FILE"
    echo "===== 所有步骤成功完成 =====" >> "$LOG_FILE"
    echo "输出目录: $OUTPUT_DIR" >> "$LOG_FILE"
    
    # 生成结果摘要
    SUMMARY_FILE="$OUTPUT_DIR/${BASENAME}_summary.txt"
    echo "MAGpurify处理摘要" > "$SUMMARY_FILE"
    echo "输入文件: $INPUT_FASTA" >> "$SUMMARY_FILE"
    echo "开始时间: $(head -1 "$LOG_FILE" | cut -d'-' -f2-)" >> "$SUMMARY_FILE"
    echo "结束时间: $(tail -1 "$LOG_FILE" | cut -d'-' -f2-)" >> "$SUMMARY_FILE"
    
    echo "" >> "$SUMMARY_FILE"
    echo "生成的输出文件:" >> "$SUMMARY_FILE"
    find "$OUTPUT_DIR" -type f -name "${BASENAME}*" -exec ls -lh {} \; >> "$SUMMARY_FILE"
    
    echo "" >> "$SUMMARY_FILE"
    echo "各步骤执行结果:" >> "$SUMMARY_FILE"
    grep "步骤 .*:.*" "$LOG_FILE" >> "$SUMMARY_FILE"
    
    echo ""
    echo "处理成功完成!"
    echo "详细日志: $LOG_FILE"
    echo "结果摘要: $SUMMARY_FILE"
else
    echo ""
    echo "处理过程中出错，请检查日志文件: $LOG_FILE"
    exit 1
fi

# 使用：sh MAG_magpurify_test.sh cleaned.fa output_dir

# step4: 批量生成命令
# cd /public/home/zhuxue/Qomolangma_Human/02.MAG/02.MAG_rm_contamination/code
# 【执行三次 gut/oral/skin】
input_dir="/public/home/zhuxue/Qomolangma_Human/02.MAG/02.MAG_rm_contamination/01.MAG_rm_nonstandard/skin"
output_dir="/public/home/zhuxue/Qomolangma_Human/02.MAG/02.MAG_rm_contamination/02.magpurify_out_skin"

output_script="magpurify_code_skin.sh"
> "$output_script"  # 清空或新建脚本文件

# 遍历输入目录下的所有文件
for file in "$input_dir"/*; do
  if [[ -f "$file" ]]; then
    filename=$(basename "$file" .fa)
    echo "sh magpurify.sh \"$file\" \"$output_dir/$filename\"" >> "$output_script"
  fi
done

# step4:使用clean-bin对magpurify_output中的结果进行识别汇总，删去cleaned.fa中的污染序列，并生成去除污染序列后的新MAG，final_cleaned.fa
# 批量生成命令
INPUT_BASE="/public/home/zhuxue/Qomolangma_Human/02.MAG/02.MAG_rm_contamination/01.MAG_rm_nonstandard"
MID_BASE="/public/home/zhuxue/Qomolangma_Human/02.MAG/02.MAG_rm_contamination"
OUTPUT_BASE="/public/home/zhuxue/Qomolangma_Human/02.MAG/02.MAG_rm_contamination/03.magpurify_out"

OUTPUT_FILE="/public/home/zhuxue/Qomolangma_Human/02.MAG/02.MAG_rm_contamination/code/clean_bin.sh"

> "$OUTPUT_FILE"

# Process each folder (gut, skin, oral)
for folder in gut skin oral; do
    # Find all .fa files in the current folder
    for fa_file in "$INPUT_BASE/$folder"/*.fa; do
        if [ -f "$fa_file" ]; then
            # Get the filename without path and .fa extension
            filename=$(basename "$fa_file" .fa)
            
            # Construct second argument path
            mid_path="${MID_BASE}/02.magpurify_out_${folder}/${filename}"
            
            # Construct third argument path
            out_path="${OUTPUT_BASE}/${folder}/${filename}.fa"
            
            # Write the magpurify command to the output file
            echo "magpurify clean-bin \"$fa_file\" \"$mid_path\" \"$out_path\"" >> "$OUTPUT_FILE"
        fi
    done
donecd /public/home/zhuxue/Qomolangma_Human/02.MAG/02.MAG_rm_contamination/code
split -l 350 clean_bin.sh --numeric-suffixes=1 --suffix-length=2 split_code_2/code_magpurify clean-bin cleaned.fasta magpurify_output final_cleaned.fa

# magpurify clean-bin "/public/home/zhuxue/Qomolangma_Human/02.MAG/02.MAG_rm_contamination/01.MAG_rm_nonstandard/gut/QML0000-240606-02B_bin.10.strict.fa" "/public/home/zhuxue/Qomolangma_Human/02.MAG/02.MAG_rm_contamination/02.magpurify_out_gut/QML0000-240606-02B_bin.10.strict" "/public/home/zhuxue/Qomolangma_Human/02.MAG/02.MAG_rm_contamination/03.magpurify_out/gut/QML0000-240606-02B_bin.10.strict.fa"