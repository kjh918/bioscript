
options(stringsAsFactors=FALSE)   
suppressPackageStartupMessages(library("openxlsx"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("dplyr"))

INDEXING = function(X) { sapply(unique(X), function(y) list(y)) }

sampleList = read.delim('sample.list2', header=F)[,1]

optitype.ResDir = 'optitype_results'
hla_la.ResDir   = 'HLA-LA-result'

Optitype.Res = ldply(lapply( INDEXING(sampleList), function(sl) {
    opti.res = read.delim(sprintf('%s/%s_result.tsv',optitype.ResDir, sl))[,-1]
    return(opti.res)
}), .id='sampleid')

HLA.LA.Res = lapply( INDEXING(sampleList), function(sl) {
    hla.res = read.delim(sprintf('%s/%s/R1_bestguess_G.txt', hla_la.ResDir, sl))
    hla.res = hla.res %>% select(c('Locus','Chromosome','Allele')) %>% mutate( sampleid = sl)
    return(hla.res)
})

HLA.LA.Res2 = lapply(HLA.LA.Res, function(y) {
    
    y2 = y %>% mutate( tag = paste0(Locus, Chromosome), optitype="" )
    for( i in 1:6 )
    {
        y2[i,'optitype'] = Optitype.Res[which(Optitype.Res$sampleid %in% y2$sampleid), y2[i, 'tag'] ]
    } 
    y2 = y2 %>% select(c('sampleid','Locus','Chromosome','Allele','optitype')) %>% 
        rename( ID=sampleid, 'HLA-LA'=Allele, OptiType=optitype)
    return(y2)
} )
openxlsx::write.xlsx(HLA.LA.Res2, 'OrganoidScience_HLA_Typing_Result.xlsx')



HLARes = read.delim(sprintf('%s/R1_bestguess_G.txt', hla_la.ResDir))

RES = HLARes %>% select(c('Locus','Chromosome','Allele')) %>% mutate( digit = paste0(Locus, Chromosome))
RES$OptiType=""

for( n in 1:6 ){
    RES[n, 'OptiType'] = OptitypeRes[ 4 ,which( colnames(OptitypeRes) == RES[n,'digit'])]
}




