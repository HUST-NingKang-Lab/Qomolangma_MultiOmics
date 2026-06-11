# 下载参考的MAG数据
# HumGut_gut：https://arken.nmbu.no/~larssn/humgut/ （2025-6-12， 17G）

# eHOMD_oral: https://www.ehomd.org/download/dld_genome_table_all/browser
# 【口腔的三条数据下载不下来】
# 【GCF_000499705：被NCBI管理员删除】
# 【GCF_030218885：被新版本替代】
# 【GCF_030223965：被新版本替代】

# HSMG_skin：
# 1、删除文件“data_download_links.txt”中的第1列，再将文件拆分成20个文件
# 处理空格分隔的文件
awk '{
    # 保存第二列（$2）之后的所有内容，包括空格
    rest = $0;
    sub(/^[[:space:]]*[^[:space:]]+[[:space:]]*/, "", rest);
    print "wget -c " rest
}' data_download_links.txt > processed_links.txt

# 拆分文件（每个800行）
split -d -a 2 -l 800 processed_links.txt download_part_ --additional-suffix=".txt"

# 清理中间文件
rm processed_links.txt

# 2、开始下载
#!/bin/bash

# 设置源文件和目标路径
download_dir="/public/home/zhuxue/Qomolangma_Human/database/compare_MAG/skin/download_path"
output_dir="/public/home/zhuxue/Qomolangma_Human/database/compare_MAG/skin/raw_MAG"

# 创建目标目录
mkdir -p "$output_dir"

# 创建下载函数
download_files() {
    local part_file="$1"
    local part_number="${part_file##*/download_part_}"
    part_number="${part_number%.txt}"
    local target_dir="${output_dir}/part_${part_number}"
    
    # 创建分区目录
    mkdir -p "$target_dir"
    
    echo "开始下载part_${part_number} (800个文件)..."
    cd "$target_dir" || return 1
    
    # 执行下载命令并记录日志
    while read -r cmd; do
        echo "正在下载: ${cmd##*/}"
        # 重试逻辑（最多尝试3次）
        for i in {1..3}; do
            $cmd && break
            echo "重试($i/3): ${cmd##*/}"
            sleep $((i*5))
        done
    done < "$part_file"
    
    echo "完成下载part_${part_number}"
}

# 导出函数以便在并行环境中使用
export -f download_files
export output_dir

# 并行执行所有下载
max_parallel=8  # 控制最大同时下载进程数
find "$download_dir" -name "download_part_*.txt" | parallel -j "$max_parallel" --bar download_files {}

echo "所有下载任务已完成！"
echo "文件保存在: $output_dir"### （1）将skin的数据整合在一起
cp /public/home/zhuxue/Qomolangma_Human/database/compare_MAG/HSMG_skin/raw_MAG01/batch_*/*.fa.gz /public/home/zhuxue/Qomolangma_Human/database/compare_MAG/HSMG_skin/raw_MAG

### （2）对skin的数据去冗余
source /public/home/zhuxue/miniconda3/bin/activate
conda activate drep
INPUT="/public/home/zhuxue/Qomolangma_Human/database/compare_MAG/HSMG_skin/raw_MAG/*.fa"
OUTPUT="/public/home/zhuxue/Qomolangma_Human/database/compare_MAG/HSMG_skin/drep_MAG"
dRep dereplicate $OUTPUT -nc 0.3 -sa 0.95 -p 64 -l 0 --skip_plots --ignoreGenomeQuality -g $INPUT

### （3）和skin MAG进行去冗余比较
cp /public/home/zhuxue/Qomolangma_Human/database/compare_MAG/HSMG_skin/drep_MAG/dereplicated_genomes/*.fa /public/home/zhuxue/Qomolangma_Human/02.MAG/08.MAG_compare/HSMG_skin/mix_MAG/
cp /public/home/zhuxue/Qomolangma_Human/02.MAG/04.MAG_dRep/02.drep_out/dereplicated_genomes/*.fa /public/home/zhuxue/Qomolangma_Human/02.MAG/08.MAG_compare/HSMG_skin/mix_MAG/

source /public/home/zhuxue/miniconda3/bin/activate
conda activate drep
INPUT="/public/home/zhuxue/Qomolangma_Human/02.MAG/08.MAG_compare/HSMG_skin/mix_MAG/*.fa"
OUTPUT="/public/home/zhuxue/Qomolangma_Human/02.MAG/08.MAG_compare/HSMG_skin/drep_out"
dRep dereplicate $OUTPUT -nc 0.3 -sa 0.95 -p 64 -l 0 --skip_plots --ignoreGenomeQuality -g $INPUT

# eHOMD_oral：
# 1、提取文件“browser”第10列，并对原下载链接进行修改：
# 提供路径为：
#  ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/008/185/GCA_000008185.1_ASM818v1
# 下载路径为：
#  ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/008/185/GCA_000008185.1_ASM818v1/GCA_000008185.1_ASM818v1_genomic.fna.gz
# 并在每一行前面加上“wget -c”，生成download.sh
awk -F'/' '{
    # 获取basename（最后一个斜杠后的内容）
    basename = $NF
    # 构建下载命令
    print "wget -c " $0 "/" basename "_genomic.fna.gz"
}' url.txt > download.sh

# 2、拆分为3个可执行文件
total_lines=$(wc -l < download.sh)

# 计算每个文件应包含的行数（向上取整）
lines_per_file=$(( (total_lines + 2) / 3 ))

# 使用split命令分割文件
split -l $lines_per_file -d -a 1 download.sh "/public/home/zhuxue/Qomolangma_Human/database/compare_MAG/oral/split/part_"

# 为所有分割文件添加.sh扩展名并确保可执行
for file in /public/home/zhuxue/Qomolangma_Human/database/compare_MAG/oral/split/part_*; do
    mv "$file" "${file}.sh"
    chmod +x "${file}.sh"
done

# Gut
# 对oral的数据去冗余
source /public/home/zhuxue/miniconda3/bin/activate
conda activate drep
INPUT="/public/home/zhuxue/Qomolangma_Human/database/compare_MAG/eHOMD_oral/raw_MAG/*.fna"
OUTPUT="/public/home/zhuxue/Qomolangma_Human/database/compare_MAG/eHOMD_oral/drep_MAG"
dRep dereplicate $OUTPUT -nc 0.3 -sa 0.95 -p 64 -l 0 --skip_plots --ignoreGenomeQuality -g $INPUT

### （3）和oral MAG进行去冗余比较
cp /public/home/zhuxue/Qomolangma_Human/database/compare_MAG/eHOMD_oral/drep_MAG/dereplicated_genomes/*.fna /public/home/zhuxue/Qomolangma_Human/02.MAG/08.MAG_compare/eHOMD_oral/mix_MAG
cp /public/home/zhuxue/Qomolangma_Human/02.MAG/04.MAG_dRep/02.drep_out/dereplicated_genomes/*.fa /public/home/zhuxue/Qomolangma_Human/02.MAG/08.MAG_compare/eHOMD_oral/mix_MAG

source /public/home/zhuxue/miniconda3/bin/activate
conda activate drep
INPUT="/public/home/zhuxue/Qomolangma_Human/02.MAG/08.MAG_compare/eHOMD_oral/mix_MAG/*a"
OUTPUT="/public/home/zhuxue/Qomolangma_Human/02.MAG/08.MAG_compare/eHOMD_oral/drep_out"
dRep dereplicate $OUTPUT -nc 0.3 -sa 0.95 -p 64 -l 0 --skip_plots --ignoreGenomeQuality -g $INPUT# 正式处理

## (1) oral
cd /public/home/zhuxue/Qomolangma_Human/database/compare_MAG/oral
cp raw_MAG/*/*.fna.gz drep_MAG/temp/
# 对oral的数据去冗余
source /public/home/zhuxue/miniconda3/bin/activate
conda activate drep
INPUT="/public/home/zhuxue/Qomolangma_Human/database/compare_MAG/eHOMD_oral/raw_MAG/*.fna"
OUTPUT="/public/home/zhuxue/Qomolangma_Human/database/compare_MAG/eHOMD_oral/drep_MAG"
dRep dereplicate $OUTPUT -nc 0.3 -sa 0.95 -p 64 -l 0 --skip_plots --ignoreGenomeQuality -g $INPUT