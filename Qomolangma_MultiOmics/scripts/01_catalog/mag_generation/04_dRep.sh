source /public/home/zhuxue/miniconda3/bin/activate
conda activate drep
INPUT="/public/home/zhuxue/Qomolangma_Human/02.MAG/04.MAG_dRep/01.raw_MAG_rm_S14_purify_Inferior_quality/*.fa"
OUTPUT="/public/home/zhuxue/Qomolangma_Human/02.MAG/04.MAG_dRep/02.drep_out"
dRep dereplicate $OUTPUT -nc 0.3 -sa 0.95 -p 64 -l 0 --skip_plots --ignoreGenomeQuality -g $INPUT

# 去除了低质量的，所以直接dRep，不用GenomeQuality