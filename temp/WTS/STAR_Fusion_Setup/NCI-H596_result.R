


base_dir="/storage/home/kangsm/analysis/gcx.wts.pipeline.fusionGenesDetection/NCI-H596_RNAseq_test_run"

r1=openxlsx::read.xlsx(sprintf("%s/NCI-H596_rep1/fusionGenes_spliceVariants/NCI-H596_rep1.FusionGenes.SpliceVariants.CancerIntrons.Summary.xlsx",base_dir), sheet=5)
r2=openxlsx::read.xlsx(sprintf("%s/NCI-H596_rep2/fusionGenes_spliceVariants/NCI-H596_rep2.FusionGenes.SpliceVariants.CancerIntrons.Summary.xlsx",base_dir), sheet=5)
r3=openxlsx::read.xlsx(sprintf("%s/NCI-H596_rep3/fusionGenes_spliceVariants/NCI-H596_rep3.FusionGenes.SpliceVariants.CancerIntrons.Summary.xlsx",base_dir), sheet=5)

vl = list(
    rep1=r1$intron,
    rep2=r2$intron,
    rep3=r3$intron
)
vennListName=names(vl)
VN = ggVennDiagram::Venn(vl)

vl2 = list(
    rep1=unique(r1$genes),
    rep2=unique(r2$genes),
    rep3=unique(r3$genes)
)
vennListName=names(vl2)
ggVennDiagram(vl2)
VN = ggVennDiagram::Venn(vl2)
