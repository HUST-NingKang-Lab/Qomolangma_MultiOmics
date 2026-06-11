# step1：安装与配置BiG-SCAPE
# 1、 下载Bigscape(1.1.5, python=3.9）
conda create -n bigscape python=3.8
conda activate bigscape

git clone https://github.com/medema-group/BiG-SCAPE.git
cd BiG-SCAPE
pip install -r requirements.txt

# 2、 下载Pfam数据库（38.0     2025.08.26）
# wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam38.0/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm  # 生成索引文件

# 3、下载fasttree （2.1.11）
conda install -c bioconda fasttree

# 4、下载mafft
conda install -c bioconda mafft


# step2：处理数据
# 运行antiSMASH后，每个MAG会生成独立目录，其中包含regionGBK文件，该文件包含BGC的基因注释信息，将所有MAG的regionGBK文件汇总到同一目录，确保文件名无冲突
mkdir gbk_files
#antismash_out文件夹下是所有的MAG的BGC文件夹，每个文件夹下有多个gbk文件
ls /public/home/zhuxue/Qomolangma_Human/02.MAG/06.MAG_antismash/antismash_out/ >MAG_list
for i in `cat MAG_list`
do
  echo $i
  #把*region*.gbk提取出来，如果没有region就跳过，如果是cluster就修改一下
  find /public/home/zhuxue/Qomolangma_Human/02.MAG/06.MAG_antismash/antismash_out/$i -name "*region*.gbk" >tmp_list
  for j in `cat tmp_list`
  do
    #重命名成genome_name_region的形式并拷贝到一个目录
    cp $j gbk_files/${i}_$(basename $j)
  done
done

# step3：参数设置
python3 bigscape.py \
  -i /public/home/zhuxue/Qomolangma_Human/02.MAG/07.MAG_bigscape/gbk_files \       #输入gbk文件夹
  -o /public/home/zhuxue/Qomolangma_Human/02.MAG/07.MAG_bigscape/bigscape_out \    # 输出目录
  --pfam_dir /public/home/zhuxue/Qomolangma_Human/database/Pfam \
  -c 64 \                   # 使用64个CPU核心
  --mix \
  --mibig \                  # 整合MIBiG v3.1
  --mode auto