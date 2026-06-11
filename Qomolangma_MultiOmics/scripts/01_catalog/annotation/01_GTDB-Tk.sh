# step1: 下载226版本的GTDB
# https://data.ace.uq.edu.au/public/gtdb/data/releases/release226/226.0/auxillary_files/gtdbtk_package/
# step2: 建树
source /public/home/zhuxue/miniconda3/bin/activate
conda activate gtdbtk-2.4.1
conda env config vars set GTDBTK_DATA_PATH="/public/home/zhuxue/Qomolangma_Human/database/GTDB226/release226"

cd /public/home/zhuxue/Qomolangma_Human/02.MAG/07.MAG_gtdb
gtdbtk classify_wf --genome_dir /public/home/zhuxue/Qomolangma_Human/02.MAG/04.MAG_dRep/02.drep_out/dereplicated_genomes --skip_ani_screen --out_dir ann --extension fa --cpus 64 --force --tmpdir tmp --scratch_dir scratch

cd /public/home/zhuxue/Qomolangma_Human/02.MAG/07.MAG_gtdb
gtdbtk infer --msa_file /public/home/zhuxue/Qomolangma_Human/02.MAG/07.MAG_gtdb/ann/align/gtdbtk.bac120.user_msa.fasta.gz --out_dir infer/infer_user_bac --cpus 64 --prefix user_bac120 --tmpdir ./tmp/user_bac120
gtdbtk infer --msa_file /public/home/zhuxue/Qomolangma_Human/02.MAG/07.MAG_gtdb/ann/align/gtdbtk.ar53.user_msa.fasta.gz --out_dir infer/infer_user_ar --cpus 64 --prefix user_ar53 --tmpdir ./tmp/user_ar53
gtdbtk infer --msa_file /public/home/zhuxue/Qomolangma_Human/02.MAG/07.MAG_gtdb/ann/align/gtdbtk.bac120.msa.fasta.gz --out_dir infer/infer_all_bac --cpus 60 --prefix all_bac120 --tmpdir ./tmp/all_bac120
gtdbtk infer --msa_file /public/home/zhuxue/Qomolangma_Human/02.MAG/07.MAG_gtdb/ann/align/gtdbtk.ar53.msa.fasta.gz --out_dir infer/infer_all_ar --cpus 60 --prefix all_ar53 --tmpdir ./tmp/all_ar53