#python3.11 run.py \
#  --bam /storage/home/jhkim/Projects/NIPT/GCX-NIPT-260121/Results/WGS.cfDNA.ControlSamples/SRR38035992/bam/SRR38035992.recal.filtered.bam \
#  --bed /storage/home/jhkim/Projects/NIPT/GCX-NIPT-260121/Workplace/scripts/genes.bed \
#  --out /storage/home/jhkim/Projects/NIPT/GCX-NIPT-260121/Results/SRR38035992.tsv \
#  --short-max 150 \
#  --mono-min 151 \
#  --mono-max  220 \
#  --min-mapq 30

python3.11 check_fragment_count.py \
  --bam /storage/home/jhkim/Projects/NIPT/GCX-NIPT-260121/Results/WGS.cfDNA.ControlSamples/SRR38035992/bam/SRR38035992.recal.filtered.bam \
  --bed /storage/home/jhkim/Projects/NIPT/GCX-NIPT-260121/Resources/mrd_ic.bed \
  --out /storage/home/jhkim/Projects/NIPT/GCX-NIPT-260121/Results/SRR38035992.tsv \
  --short-max 150 \
  --mono-min 151 \
  --mono-max  220 \
  --min-mapq 30
