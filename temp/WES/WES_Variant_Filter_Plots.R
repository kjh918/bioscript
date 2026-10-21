
library(ggplot2)

vcf_idepth1 = read.delim('3-100X.mutect2.single.idepth')
vcf_idepth3 = read.delim('3-100X.mutect2.single.filtered.only.idepth')



vcf_smdepth1 = read.delim('3-100X.mutect2.single.ldepth.mean') %>% mutate(STATUS='Before Filtering')
vcf_smdepth3 = read.delim('3-100X.mutect2.single.filtered.only.ldepth.mean') %>% mutate(STATUS='After Filtering')


rbind(vcf_smdepth1, vcf_smdepth3) %>% ggplot( aes(MEAN_DEPTH)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + xlim(0,100) + facet_wrap(~STATUS)

ggplot(vcf_smdepth1, aes(MEAN_DEPTH)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + xlim(0,100)

ggplot(vcf_smdepth3, aes(MEAN_DEPTH)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + xlim(0,100)

summary(vd$MEAN_DEPTH)


ggplot(mf@data, aes(tumor_f)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + xlim(0,100)
