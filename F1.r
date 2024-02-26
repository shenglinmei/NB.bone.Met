


library(scProcess)
library(ggpubr)
library('lib.r')


scon = read.RDS('dat.conos.rds')

# Marker gene expression
gs=c('ISL1','PRPH','TH','PNMT','GAP43','EPAS1','RGS5','HAND2','MYCN')
b = lapply(sn(gs),function(x) scon$plotGraph(gene=x,title=x,size=0.1,raster = TRUE))
fig=cowplot::plot_grid(plotlist=b, ncol=3, nrow=3)
fig


# cell annotations
anof = readRDS('cell.ano.rds')
anof = as.factor(cells)
anoM.pal <- setNames(sample(rainbow(length(levels(anof)))),levels(anof))


a1 = scon$plotGraph(groups = anof,size=0.05,raster = TRUE,alpha=0.07,plot.na=F,palette=anoM.palf,font.size = c(3.7,4.1))
a1
a2 = scon$plotGraph(groups = anof,size=0.05,raster = TRUE,alpha=0.07,plot.na=F,palette=anoM.palf,font.size = c(3.7,4.1))
a2

ggsave('no.met.overview.pdf',a1,height=4,width=4)


ggsave('met.overview.pdf',a2,height=4,width=4)


a1 = scon$plotGraph(groups = anof,size=0.05,raster = TRUE,alpha=0.07,plot.na=F,palette=anoM.palf,font.size = c(3.7,4.1))
a1
ggsave('F1.overview.pdf',a1,height=4,width=4)


# sample groups

Met= c("NB35-BM1","NB33-BM3",'NB34_BM','BM10','NBBM2','NBBM4','NB34_Bone')
# 
sample.group[sample.group %in% Met]='metastatic'
sample.group[sample.group %in% c("NB21-BM2R","NB37-BM5" ,"BM6", "BM7",'NBBM1','NBBM3','NBBM5','NBBM6')]='non-metastatic'
sample.group['NBBM6']='non-metastatic'
sample.group

getConditionPerCell=
function(sample.per.cell,sample.groups) {
    sample.per.cell %>%
      {setNames(as.character(sample.groups[as.character(.)]), names(.))} %>%
      as.factor()
}

ssamp = scon$getDatasetPerCell() %>% Toch()
stype = getConditionPerCell(ssamp,as.factor(sample.group)[unique(ssamp)])
table(stype)

a1=scon$plotGraph(groups=stype,size=0.01,alpha=0.1,raster = TRUE,show.legend=TRUE,mark.groups=F)+theme(legend.position=c(0.16,0.89))
a1

ggsave('F1.fraction.pdf',a1,height=4,width=4)


# sample fractions
sample.group = factor(sample.group,levels=c('non-metastatic','metastatic'))
cname = names(anof)
cname = setdiff(cname,rm)
ano2 = data.frame(Cell = anof[cname], SampleType = ssamp[cname])
tmp2 <- acast(ano2, Cell ~ SampleType, fun.aggregate = length)
head(tmp2)
tmp3 <- (sweep(tmp2, 2, colSums(tmp2), FUN = "/"))
tmp4 <- melt(tmp3)
head(tmp4)
names(tmp4) <- c("annot", "sample", "pc.of.sample")
head(tmp4)
tmp4$Group = sample.group[Toch(tmp4$sample)]

library(ggpubr)
library(ggplot2)
ylab = 'fraction of total cells'
#tmp4 = tmp4[tmp4$cell %in% c('macrophage','Proliferating T','CTL','Treg','Progenitors','Naive T','Thelper'),]
 p <- ggplot(na.omit(tmp4),aes(x=annot,y=pc.of.sample,dodge=Group,fill=Group))+geom_boxplot(notch=FALSE,outlier.shape=NA)  +  geom_point(position = position_jitterdodge(jitter.width=0.1),color=adjustcolor(1,alpha=0.3),pch=19,size=0.5)+theme_classic()+
      theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y = element_text(angle = 90, hjust = 0.5))  +xlab("") +ylab(ylab)+theme(legend.position="top") #+
    #  scale_fill_manual(values=fraction.palette1)
p=p+scale_fill_manual(values=fraction.palette1)+ylim(c(min(tmp4$pc.of.sample),max(tmp4$pc.of.sample)*1.4))
p=p+ stat_compare_means(label = "p.signif")#,method='t.test'
p


# Tumor markers


gs=c('PHOX2B','CDK4','MDK','HAND2','KCNQ2','MDK')

b = lapply(sn(gs),function(x) scon$plotGraph(gene=x,title=x,size=0.1,raster = TRUE))
fig=cowplot::plot_grid(plotlist=b, ncol=3, nrow=3)
fig