### step1: 下载数据集
wget -c ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
wget -c ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin
gunzip Rfam.cm.gz
cmpress Rfam.cm

### step2: batch command
cd /public/home/zhuxue/Qomolangma_Human/02.MAG/09.MAG_infernal/code
export Rfam="/public/home/zhuxue/MT_zky/database/Rfam"
export Binning="/public/home/zhuxue/Qomolangma_Human/02.MAG/04.MAG_dRep/01.raw_MAG_rm_S14_purify_Inferior_quality_rm_Tc"
export RNA_detect="/public/home/zhuxue/Qomolangma_Human/02.MAG/09.MAG_infernal/infernal_out"
for i in $Binning/*.fa;do
filename=$(basename "$i")
echo cmscan --cpu 64 --cut_ga --rfam --noali --fmt 2 --oskip --tblout $RNA_detect/$filename $Rfam/Rfam.cm $Binning/$filename >> infernal.sh
done

### step3: split code
cd /public/home/zhuxue/Qomolangma_Human/02.MAG/09.MAG_infernal/code
split -l 300 infernal.sh --numeric-suffixes=1 --suffix-length=2 split_code/code_

source /public/home/zhuxue/miniconda3/bin/activate
conda activate metawrap_1.3.2