# step1:生成批量执行文件
# 1、每个MAG对应的antismash命令
#!/bin/bash

# 定义输入和输出路径
input_dir="/public/home/zhuxue/Qomolangma_Human/02.MAG/04.MAG_dRep/01.raw_MAG_rm_S14_purify_Inferior_quality_rm_Tc"
output_file="antismash.sh"

# 确保输出目录存在
mkdir -p "$output_dir"

# 创建输出文件（如果已存在则覆盖）
> "$output_file"

# 遍历所有FA文件
for fa_file in "$input_dir"/*.fa; do
    # 获取不带路径的文件名
    filename=$(basename -- "$fa_file")
    
    # 提取basename（去除.fa扩展名）
    basename="${filename%.fa}"
    
    # 生成antismash命令
    command="antismash --cb-general --cb-knownclusters --cb-subclusters --cc-mibig --asf --allow-long-headers "
    command+="--output-dir ${output_dir}/${basename} "
    command+="--minlength 5000 --cpus 64 "
    command+="\"$fa_file\" "
    command+="--genefinding-tool prodigal-m --taxon bacteria "
    command+="--output-basename $basename"
    
    # 将命令写入文件
    echo "$command" >> "$output_file"
done

echo "已生成 $(wc -l < "$output_file") 条antismash命令到 $output_file"
echo "每个FA文件对应的输出目录已创建在: $output_dir"

# 2、增加创建对应文件夹命令
#!/bin/bash

# 设置输入和输出文件
input_file="antismash.sh"
output_file="antismash_with_mkdir.sh"
log_file="process.log"

# 创建目录结构的主路径
antismash_out_dir="/public/home/zhuxue/Qomolangma_Human/02.MAG/06.MAG_antismash/antismash_out"
mkdir -p "$antismash_out_dir"

# 处理每行命令
counter=0
while IFS= read -r line; do
    ((counter++))
    echo "处理命令 $counter/$total_commands" | tee -a "$log_file"
    
    # 提取--output-basename的值（使用复杂模式匹配）
    output_basename=$(echo "$line" | grep -Po -- '--output-basename\s+\K[\w\.,;:@=+()-]+[^\s"]*' | head -1)
    
    # 验证提取的值
    if [ -z "$output_basename" ]; then
        echo "错误: 未找到有效的--output-basename值" | tee -a "$log_file"
        echo "原始命令: $line" | tee -a "$log_file"
        continue
    fi
    
    # 创建目录命令
    mkdir_cmd="mkdir -p \"${antismash_out_dir}/${output_basename}\""
    
    # 添加到输出文件
    echo "$mkdir_cmd" >> "$output_file"
    echo "$line" >> "$output_file"
    echo "" >> "$output_file"  # 添加空行分隔
    
    # 创建目录（可选）
    mkdir -p "${antismash_out_dir}/${output_basename}"
    
done < "$input_file"

echo "处理完成: $(date)" | tee -a "$log_file"
echo "生成文件: $output_file 包含 $counter 条增强命令" | tee -a "$log_file"

# step2:拆分执行文件
#!/bin/bash

# 设置输入文件和输出目录
input_file="/public/home/zhuxue/Qomolangma_Human/02.MAG/06.MAG_antismash/code/antismash_with_mkdir.sh"
output_dir="/public/home/zhuxue/Qomolangma_Human/02.MAG/06.MAG_antismash/code/split"

# 计算总行数和每个文件的行数
total_lines=$(wc -l < "$input_file")
lines_per_file=$(( (total_lines + 48) / 49 ))  # 向上取整确保覆盖所有行

# 使用split命令进行拆分
split -d -a 3 -l "$lines_per_file" "$input_file" "${output_dir}/command_part_"

# 重命名文件以添加.sh扩展名
for file in "$output_dir"/command_part_*; do
    mv "$file" "${file}.sh"
done

# 报告结果
file_count=$(ls "$output_dir" | wc -l)
echo "拆分完成!"

# step3:批量运行（49个执行文件）
source /public/home/zhuxue/miniconda3/bin/activate
conda activate antiSMASH

# 设置工作目录
split_dir="/public/home/zhuxue/Qomolangma_Human/02.MAG/06.MAG_antismash/code/split"
log_dir="/public/home/zhuxue/Qomolangma_Human/02.MAG/06.MAG_antismash/code/slurm_logs"
mkdir -p "$log_dir"

# 获取脚本文件列表
script_files=($(find "$split_dir" -maxdepth 1 -name "*.sh" | sort))
total_files=${#script_files[@]}

# 获取当前任务对应的脚本
script_index=$((SLURM_ARRAY_TASK_ID - 1))
script="${script_files[$script_index]}"
script_name=$(basename "$script")

# 执行脚本
echo "[$(date +'%F %T')] 开始执行: $script_name"
bash "$script"
exit_status=$?
echo "[$(date +'%F %T')] 执行完成: $script_name (状态: $exit_status)"

exit $exit_status”