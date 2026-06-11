### （1）把文件平均分成分成多分，分批次运行
src_dir="/public/home/zhuxue/Qomolangma_Human/02.MAG/02.MAG_rm_contamination/03.magpurify_out/skin"
dst_root="/public/home/zhuxue/Qomolangma_Human/02.MAG/03.MAG_checkM/MAG_tmp/skin"

mkdir -p "$dst_root"

# 获取所有文件
files=($(find "$src_dir" -maxdepth 1 -type f))
total=${#files[@]}

# 计算每份的文件数量（向上取整分5份）
group_size=$(( (total + 4) / 5 ))  # 向上取整：(total + 份數 -1 ) / 份數

# 随机打乱文件顺序
shuffled=($(printf "%s\n" "${files[@]}" | shuf))

# 分5组处理
for i in $(seq 1 5); do
    subdir=$(printf "%s/skin_%02d" "$dst_root" "$i")
    mkdir -p "$subdir"
    start=$(( (i - 1) * group_size ))
    end=$(( start + group_size - 1 ))

    for j in $(seq $start $end); do
        [[ $j -ge $total ]] && break
        cp "${shuffled[$j]}" "$subdir/"
    done
done

### （2）运行CheckM
source /public/home/zhuxue/miniconda3/bin/activate
conda activate metawrap_1.3.2

cd /public/home/zhuxue/Qomolangma_Human/02.MAG/03.MAG_checkM

INPUTDIR="MAG_tmp/skin/skin_01"
TMPDIR="CheckM_temp_dir/skin_01"  # 这个路径无法自己创建
OUTPUT="CheckM_results_out/skin_01_checkm.tsv"
WorkDIR="Check_out/skin_01"

checkm lineage_wf -t 64 --pplacer_threads 64 -x fa --tmpdir "$TMPDIR" -f "$OUTPUT" --tab_table "$INPUTDIR" "$WorkDIR"

### （2）合并结果
cd /public/home/zhuxue/Qomolangma_Human/02.MAG/03.MAG_checkM/CheckM_results_out/summary
awk 'FNR==1 && NR!=1 {next} 1' ../gut_* > gut_bins_checkM.tsv
awk 'FNR==1 && NR!=1 {next} 1' ../oral_* > oral_bins_checkM.tsv
awk 'FNR==1 && NR!=1 {next} 1' ../skin_* > skin_bins_checkM.tsv

awk 'FNR==1 && NR!=1 {next} 1' *checkM.tsv > all_bins_checkM.tsv

### （3）为drep的输入处理结果
cd /public/home/zhuxue/Qomolangma_Human/02.MAG/04.MAG_dRep/temp

use_raw_MAG_li="use_raw_MAG_li.txt"
all_bins_checkM="all_bins_checkM.tsv"
output_file="use_bins_checkM.tsv"

# 提取表头 + 筛选数据行
(head -n 1 "$all_bins_checkM" && tail -n +2 "$all_bins_checkM" | \
awk 'NR==FNR {a[$1]; next} $1 in a' "$use_raw_MAG_li" - ) > "$output_file"

# 获得低质量的checkM数据
awk -F'\t' 'NR == 1 || ($12 < 50 || $13 > 10)' use_bins_checkM.tsv > Inferior_quality_bins_checkM.tsv
awk -F'\t' 'NR > 1 {print $1}' Inferior_quality_bins_checkM.tsv > Inferior_quality_MAG_li.txt


# 从drep的输入文件夹，删掉这些文件
target_dir="/public/home/zhuxue/Qomolangma_Human/02.MAG/04.MAG_dRep/01.raw_MAG_rm_S14_purify_Inferior_quality"
prefix_file="/public/home/zhuxue/Qomolangma_Human/02.MAG/03.MAG_checkM/CheckM_bu/CheckM_results_out/Inferior_quality_MAG_li.txt"
# 逐行读取文件中的前缀，并删除对应的 .fa 文件
while IFS= read -r prefix
do
    rm -f "${target_dir}/${prefix}.fa"  # 删除匹配的 .fa 文件
done < "$prefix_file"

echo "删除完成。"