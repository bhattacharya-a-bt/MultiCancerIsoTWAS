require(data.table)
require(wesanderson)
require(dplyr)
require(vroom)
setwd('/Users/abhattacharya3/OneDrive - Inside MD Anderson/Fall 2024/PanCan IsoTWAS/Results/')

pli_v4 = fread('/Users/abhattacharya3/Downloads/gnomad.v4.1.constraint_metrics.tsv')
pli_gene = pli_v4 %>%
  group_by(gene_id) %>%
  summarize(pli = max(lof.pLI))
colnames(pli_gene) = c('Gene','pLI')

shet = readxl::read_excel('41588_2024_1820_MOESM4_ESM.xlsx',sheet = 2)
shet = data.frame(Gene = shet$ensg,
                  shet = shet$post_mean)

isotwas = as.data.frame(vroom('isoTWAS_FineMap_PanCan_090524.tsv'))[,1:17]
isotwas = isotwas[!duplicated(isotwas),]

twas = as.data.frame(fread('TWAS_FineMap_PanCan_090524.tsv'))
twas = twas[!duplicated(twas),]


isotwas = isotwas[,-grep('pLI', colnames(isotwas))] 
isotwas = merge(isotwas,pli_gene,by='Gene',all.x = T)
isotwas = merge(isotwas,shet,by='Gene',all.x=T)
isotwas$GTA = paste(isotwas$Gene,isotwas$Indication,sep=':')
isotwas$TTA = paste(isotwas$Transcript,isotwas$Indication,sep=':')
print(length(unique(isotwas$GTA)))
print(length(unique(isotwas$TTA)))
print(length(unique(isotwas$GTA[isotwas$pLI > .9])))
print(length(unique(isotwas$GTA[isotwas$shet > .1])))

#isotwas = subset(isotwas,in_cred_set == TRUE)
print(length(unique(isotwas$Gene)))
print(length(unique(isotwas$Gene[isotwas$pLI > .9])))
print(length(unique(isotwas$Gene[isotwas$shet > .1])))
fwrite(isotwas,'isoTWAS_FineMap_PanCan_090524_update.tsv',sep='\t',
       quote=F,col.names=T,row.names=F)

prop.test(x = c(length(unique(isotwas$GTA[isotwas$shet > .1])),
                sum(shet$shet > .1)),n = c(length(unique(isotwas$GTA)),nrow(shet)))


twas = twas[!duplicated(twas),]
twas$GTA = paste(twas$Gene,twas$Indication,sep=':')
twas = twas[,-grep('pLI', colnames(twas))] 
twas = merge(twas,pli_gene,by='Gene',all.x = T)
twas = merge(twas,shet,by='Gene',all.x = T)
print(length(unique(twas$GTA)))
print(length(unique(twas$GTA[twas$pLI > .9])))
print(length(unique(twas$GTA[twas$shet > .1])))

#twas = subset(twas,in_cred_set == TRUE)
print(length(unique(twas$Gene)))
print(length(unique(twas$Gene[twas$pLI > .9])))
print(length(unique(twas$Gene[twas$shet > .1])))


prop.test(x = c(length(unique(twas$Gene[twas$shet > .1])),
                sum(shet$shet > .1)),n = c(length(unique(twas$Gene)),nrow(shet)))


prop.test(x = c(length(unique(isotwas$GTA[isotwas$shet > .1])),
                length(unique(twas$GTA[twas$shet > .1]))),
          n = c(length(unique(isotwas$GTA)),
                length(unique(twas$GTA))))

fwrite(isotwas,'TWAS_FineMap_PanCan_090524_update.tsv',sep='\t',
       quote=F,col.names=T,row.names=F)


ceres = as.data.frame(
  fread('/Users/abhattacharya3/Downloads/CRISPR_(Project_Score,_CERES)_subsetted.csv')
)


### Figure 2
# Number of Genes
require(dplyr)
require(ggplot2)
summary_twas = twas %>%
  group_by(Indication) %>%
  summarise(Number = length(unique(GTA)),
            Number_pLI = length(unique(GTA[shet > .1])))
summary_twas$Method = 'TWAS'


summary_isotwas = isotwas %>%
  group_by(Indication) %>%
  summarise(Number = length(unique(GTA)),
            Number_pLI = length(unique(GTA[shet > .1])))
summary_isotwas$Method = 'isoTWAS'

isotwas_sum = sum(summary_isotwas$Number)
print(isotwas_sum)
twas_sum = sum(summary_twas$Number)
print(twas_sum)
isotwas_sumpli = sum(summary_isotwas$Number_pLI)
print(isotwas_sumpli)
twas_sumpli = sum(summary_twas$Number_pLI)
print(twas_sumpli)

all_num = rbind(summary_isotwas,summary_twas)
all_num$Indication = factor(all_num$Indication,
                            levels = c('Breast (Overall)',
                                       'Breast (ER+)',
                                       'Breast (ER-)',
                                       'Colorectal',
                                       'Lung',
                                       'Lung (ADE)',
                                       'Lung (SQC)',
                                       'Ovarian (overall)',
                                       'Ovarian (serous)',
                                       'Prostate (overall)',
                                       'Prostate (advanced)',
                                       'Endometrial'))
levels(all_num$Indication) = c('BRCA',
                               'BRCA (ER+)',
                               'BRCA (ER-)',
                               'CRC',
                               'Lung',
                               'LUAD',
                               'LUSC',
                               'OV',
                               'OV (ser)',
                               'PRCA',
                               'PRCA (adv)',
                               'UCEC')
all_num$Method = factor(all_num$Method,
                        levels = c('TWAS','isoTWAS'))
fwrite(all_num,
       '/Users/abhattacharya3/OneDrive - Inside MD Anderson/Fall 2024/PanCan IsoTWAS/Results/ASHG/NumGene_noFM_090524.tsv',
       quote=F,
       row.names = F,
       col.names = T,
       sep='\t')

SigGene = ggplot(data = all_num,
                 aes(x = Indication,
                     y = Number,
                     fill = Method)) +
  geom_col(position = position_dodge2(width = .9)) +
  theme_minimal() +
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=8),
        plot.title = element_text(size = 8),
        legend.title=element_text(size=8),
        legend.text=element_text(size=8),
        strip.text = element_text(size=8),
        panel.spacing=unit(1, "lines"),
        panel.border = element_rect(color = NA,
                                    fill = NA, linewidth = .1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.key.size = unit(0.2, "cm"),
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust=1),
        legend.position.inside = c(0.5, 0.9)) +
  scale_fill_manual(values = wes_palette("Darjeeling1")) +
  xlab('') +
  ylab('# signficant genes') +
  guides(
    fill  = guide_legend(position = "inside")
  )


df_p = data.frame(Indication = unique(all_num$Indication),
                  Method = 'isoTWAS',
                  P = 0)



ddd = subset(all_num,grepl('BRCA',Indication))
prop.test(x = c(sum(ddd$Number_pLI[1:3]),
                            sum(ddd$Number_pLI[4:6])),
                      n =  c(sum(ddd$Number[1:3]),
                             sum(ddd$Number[4:6])),correct = T)$p.value  

for (i in 1:nrow(df_p)){
  
  ddd = subset(all_num,Indication == df_p$Indication[i])
  df_p$P[i] = prop.test(x = ddd$Number_pLI,
                        n = ddd$Number,correct = T)$p.value   
  
}
df_p$FDR = p.adjust(df_p$P,'BH')
df_p$Label = ''
df_p$Label = ifelse(round(df_p$FDR,2) < .05 &
                      df_p$FDR >= 0.01,'*',df_p$Label)
df_p$Label = ifelse(df_p$FDR < .01 &
                      df_p$FDR >= 0.005,'**',df_p$Label)
df_p$Label = ifelse(df_p$FDR < .005,'***',df_p$Label)
plicomp =  merge(all_num,df_p,
                 by = c('Indication','Method'),
                 all.x = T)
plicomp = plicomp[,c(colnames(all_num),
                     'Label')]
plicomp$Label[is.na(plicomp$Label)] = ''

pLIGene = ggplot(data = plicomp,
                 aes(x = Indication,
                     y = Number_pLI/Number,
                     fill = Method)) +
  geom_col(position = position_dodge2(width = .9)) +
  theme_minimal() +
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=8),
        plot.title = element_text(size = 8),
        legend.title=element_text(size=8),
        legend.text=element_text(size=8),
        strip.text = element_text(size=8),
        panel.spacing=unit(1, "lines"),
        panel.border = element_rect(color = NA,
                                    fill = NA, linewidth = .1),
        legend.position='right',
        legend.box = "vertical",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.key.size = unit(0.2, "cm"),
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust=1)) +
  scale_fill_manual(values = wes_palette("Darjeeling1")) +
  xlab('') +
  ylab(expression('% genes with'~s[het] > 0.1)) +
  geom_text(aes(label = Label),
            size = 3,angle = 90) + 
  geom_hline(yintercept = mean(shet$shet > .1),
             linetype = 2,
             color = 'black')

require(readxl)
require(ggplot2)
require(wesanderson)

comparison_df = 
  read_xlsx('/Users/abhattacharya3/OneDrive - Inside MD Anderson/Fall 2024/PanCan IsoTWAS/Manuscript/Tables/SupplementalTable_SummaryNumbers.xlsx')
comparison_df$Method = factor(comparison_df$Method,
                              levels = c('TWAS','isoTWAS'))
comparison_df$Indication = factor(comparison_df$Indication,
                                  levels = c('BRCA',
                                             'BRCA (ER+)',
                                             'BRCA (ER-)',
                                             'CRC',
                                             'Lung',
                                             'LUAD',
                                             'LUSC',
                                             'OV',
                                             'OV (ser)',
                                             'PRCA',
                                             'PRCA (adv)',
                                             'UCEC'))


effN = ggplot(data = comparison_df[complete.cases(comparison_df),],
              aes(x = Indication,
                  y = `Percent increase in effective sample size`)) +
  geom_point(size = 2) +
  geom_linerange(aes(ymin = `Percent increase in effective sample size` - 1.96*`Jackknife SE`,
                     ymax = `Percent increase in effective sample size` + 1.96*`Jackknife SE`)) +
  theme_minimal() +
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=8),
        plot.title = element_text(size = 8),
        legend.title=element_text(size=8),
        legend.text=element_text(size=8),
        strip.text = element_text(size=8),
        panel.spacing=unit(1, "lines"),
        panel.border = element_rect(color = NA,
                                    fill = NA, linewidth = .1),
        legend.position='right',
        legend.box = "vertical",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.key.size = unit(0.2, "cm"),
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust=1)) +
  xlab('') +
  ylab(expression('% increase in'~N[eff]))

require(cowplot)
firstrow = plot_grid(SigGene,
                     effN,
                     pLIGene+ guides(fill = 'none'),
                     rel_widths = c(1,.8,1),
                     labels = c('A','B','C'),nrow=1)
# ggsave(filename = '/Users/abhattacharya3/OneDrive - Inside MD Anderson/Fall 2024/PanCan IsoTWAS/Manuscript/Figures/Fig2.pdf',
#        plot = firstrow,
#        width = 6.5,
#        height = 3)







require(data.table)
require(wesanderson)
require(dplyr)
require(vroom)
setwd('/Users/abhattacharya3/OneDrive - Inside MD Anderson/Fall 2024/PanCan IsoTWAS/Results/')

pli_v4 = fread('/Users/abhattacharya3/Downloads/gnomad.v4.1.constraint_metrics.tsv')
pli_gene = pli_v4 %>%
  group_by(gene_id) %>%
  summarize(pli = max(lof.pLI))
colnames(pli_gene) = c('Gene','pLI')

shet = readxl::read_excel('41588_2024_1820_MOESM4_ESM.xlsx',sheet = 2)
shet = data.frame(Gene = shet$ensg,
                  shet = shet$post_mean)

isotwas = as.data.frame(vroom('isoTWAS_FineMap_PanCan_090524.tsv'))[,1:17]
isotwas = isotwas[!duplicated(isotwas),]

twas = as.data.frame(fread('TWAS_FineMap_PanCan_090524.tsv'))
twas = twas[!duplicated(twas),]


isotwas = isotwas[,-grep('pLI', colnames(isotwas))] 
isotwas = merge(isotwas,pli_gene,by='Gene',all.x = T)
isotwas = merge(isotwas,shet,by='Gene',all.x=T)
isotwas$GTA = paste(isotwas$Gene,isotwas$Indication,sep=':')
isotwas$TTA = paste(isotwas$Transcript,isotwas$Indication,sep=':')
print(length(unique(isotwas$GTA)))
print(length(unique(isotwas$TTA)))
print(length(unique(isotwas$GTA[isotwas$pLI > .9])))
print(length(unique(isotwas$GTA[isotwas$shet > .1])))

#isotwas = subset(isotwas,in_cred_set == TRUE)
print(length(unique(isotwas$Gene)))
print(length(unique(isotwas$Gene[isotwas$pLI > .9])))
print(length(unique(isotwas$Gene[isotwas$shet > .1])))
fwrite(isotwas,'isoTWAS_FineMap_PanCan_090524_update.tsv',sep='\t',
       quote=F,col.names=T,row.names=F)

prop.test(x = c(length(unique(isotwas$Gene[isotwas$shet > .1])),
                sum(shet$shet > .1)),n = c(length(unique(isotwas$Gene)),nrow(shet)))


twas = twas[!duplicated(twas),]
twas$GTA = paste(twas$Gene,twas$Indication,sep=':')
twas = twas[,-grep('pLI', colnames(twas))] 
twas = merge(twas,pli_gene,by='Gene',all.x = T)
twas = merge(twas,shet,by='Gene',all.x = T)
print(length(unique(twas$GTA)))
print(length(unique(twas$GTA[twas$pLI > .9])))
print(length(unique(twas$GTA[twas$shet > .1])))

#twas = subset(twas,in_cred_set == TRUE)
print(length(unique(twas$Gene)))
print(length(unique(twas$Gene[twas$pLI > .9])))
print(length(unique(twas$Gene[twas$shet > .1])))

isotwas = isotwas[order(isotwas$Indication,isotwas$`Screening Adjusted P`),]
isotwas$GTA = paste(isotwas$Gene,isotwas$Indication,sep=':')
isotwas = isotwas[!duplicated(isotwas$GTA),]
gene_freq = as.data.frame(table(isotwas$Gene))
colnames(gene_freq) = c('Gene','Freq')
gene_freq = gene_freq[order(gene_freq$Freq,decreasing = T),]
isotwas_multi = subset(isotwas,Gene %in% gene_freq$Gene[gene_freq$Freq >= 5])

inTWAS = unique(isotwas_multi$HGNC[isotwas_multi$HGNC %in%
                              twas$HGNC[twas$in_cred_set == TRUE]])


df_hm = as.data.frame(expand.grid(
  Indication = unique(isotwas$Indication),
  Gene = unique(isotwas_multi$HGNC)
))

map_name_abbrev = data.frame(Indication = unique(isotwas$Indication),
                             Cancer = c('BRCA (ER+)',
                                                 'BRCA (ER-)',
                                                 'BRCA',
                                                 'CRC',
                                                 'UCEC',
                                                 'Lung',
                                                 'LUAD',
                                                 'LUSC',
                                                 'OV',
                                                 'OV (ser)',
                                                 'PRCA (adv)',
                                                 'PRCA'))

manifest_df = data.frame(File = c('0_AdvProstate.subtypes.AUTO_Finalcleaned_MAF0.01.gz',
                                  '0_all.meta.gwas.bcfr.icogs.oncoarray.overall.AUTO_Finalcleaned_MAF0.01.gz',
                                  '0_Cimba_braca1_Breast_ERneg.AUTO_Finalcleaned_MAF0.01.gz',
                                  '0_Cimba_Ovarian_serous.AUTO_Finalcleaned_MAF0.01.gz',
                                  '0_CrossCancer_colorectal.AUTO_Finalcleaned_MAF0.01.gz',
                                  '0_CrossCancer_ENDOMETRIAL.AUTO_Finalcleaned_MAF0.01.gz',
                                  '0_CrossCancer_LUNG_1000G.AUTO_Finalcleaned_MAF0.01.gz',
                                  '0_CrossCancer_LUNGadenocarcinoma_1000G.AUTO_Finalcleaned_MAF0.01.gz',
                                  '0_CrossCancer_LUNGsquamousCell_1000G.AUTO_Finalcleaned_MAF0.01.gz',
                                  '0_ERposBreastCancer.updMarch2016.AUTO_Finalcleaned_MAF0.01.gz',
                                  '0_ovarian.European_overall.AUTO_Finalcleaned_MAF0.01.gz',
                                  '0_ProstateCancer_Overall_Meta_All_20160411.AUTO_Finalcleaned_MAF0.01.gz'),
                         Cancer = c('PRCA (adv)',
                                    'BRCA',
                                    'BRCA (ER-)',
                                    'OV (ser)',
                                    'CRC',
                                    'UCEC',
                                    'Lung',
                                    'LUAD',
                                    'LUSC',
                                    'BRCA (ER+)',
                                    'OV',
                                    'PRCA'))

df_hm = merge(merge(df_hm,map_name_abbrev,by='Indication'),
              manifest_df,by='Cancer')
df_hm$File = file.path('isoTWAS',
                       sapply(strsplit(df_hm$File,
                                       '_MAF'),
                              function(x) x[1]),
                       'ScreenConfirm_isoTWAS.tsv')


df_hm$Z = 0
df_hm$Asterisk = ''
for (i in 1:nrow(df_hm)){
  print(i)
  fff = fread(df_hm$File[i])
  fff = subset(fff,HGNC == df_hm$Gene[i])
  a = sign(fff$Z[which.max(abs(fff$Z))])
  df_hm$Z[i] = a * abs(qnorm(ACAT::ACAT(unique(fff$Screen.P))))
  df_hm$Asterisk[i] = ifelse(paste(df_hm$Gene[i],
                                   df_hm$Indication[i],
                                   sep = ':') %in%
                               paste(isotwas$HGNC,
                                     isotwas$Indication,
                                     sep = ':'),
                             '*','')
  
}

hm_mat = matrix(ncol = length(unique(df_hm$Cancer)),
                nrow = length(unique(isotwas_multi$HGNC)))
colnames(hm_mat) = unique(df_hm$Cancer)
rownames(hm_mat) = unique(isotwas_multi$HGNC)

for (i in colnames(hm_mat)){
  for (j in rownames(hm_mat)){
    
    hm_mat[j,i] = subset(df_hm, Cancer == i &
                           Gene == j)$Z
    
  }
}
hm_mat[is.infinite(abs(hm_mat))] = 0

require(ggdendroplot)
rowclust = hclust(dist((hm_mat)))
colclust = hclust(dist(t(hm_mat)))
hm = hmReady(sign(hm_mat) * sqrt(abs(hm_mat)),colclus = colclust, rowclus = rowclust)
hm = merge(hm,data.frame(rowid = df_hm$Gene,
                         variable = df_hm$Cancer,
                         Asterisk = df_hm$Asterisk),
           by = c('rowid','variable'))


hmplot <- ggplot(data = hm,
                 aes(x=x,y=y,label=Asterisk)) + 
  geom_tile(data=hm, aes(x=x, y=y, fill=value)) +
  geom_text(nudge_y = -.3)
hmplot <- hmplot + scale_fill_gradient2(low = 'blue',
                                        mid = 'white',
                                        high = 'red',
                                        midpoint = 0)
hmplot = hmplot + geom_dendro(colclust, ylim=c(nrow(hm_mat) + .5, 
                                               nrow(hm_mat) + 2)) +
  geom_dendro(rowclust, xlim=c(ncol(hm_mat) + .5,
                               ncol(hm_mat) + 1.5),pointing = 'side') +
  theme_hm()+                                                          
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=8),
        plot.title = element_text(size = 8),
        legend.title=element_text(size=8),
        legend.text=element_text(size=8),
        strip.text = element_text(size=8),
        panel.spacing=unit(1, "lines"),
        panel.border = element_rect(color = NA,
                                    fill = NA, size = .1),
        legend.position='bottom',
        legend.box = "horizontal",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(face='italic',
                                   color=ifelse(unique(hm$rowid) %in% inTWAS,
                                                wes_palette('Darjeeling1')[3],
                                                wes_palette('Darjeeling1')[2])[rowclust$order]),
        legend.key.size = unit(0.4, "cm")) +
  xlab('') +
  ylab('') +
  labs(fill = expression("Z\n(sqrt scale)"))

require(cowplot)
require(ggpubr)
lz_leg = ggpubr::as_ggplot(ggpubr::get_legend(hmplot))
ppp = plot_grid(hmplot + guides(fill = 'none'),
                lz_leg,
                rel_heights = c(1,.08),
                nrow=2)
ggsave(filename = '/Users/abhattacharya3/OneDrive - Inside MD Anderson/Talks/DCCP Grand Rounds - March 2024/Heatmap.pdf',
       plot = (hmplot + guides(fill = 'none')),
       width = 6,
       height = 8)

require(enrichR)
listEnrichrSites()
setEnrichrSite("Enrichr")
dbs <- listEnrichrDbs()
dbs_query <- c("GO_Biological_Process_2023", "GO_Cellular_Component_2023", 
         "GO_Molecular_Function_2023",
         "KEGG_2021_Human",
         "Reactome_2022",
         "ChEA_2022")
enriched <- enrichr(unique(isotwas_multi$HGNC[!isotwas_multi$HGNC %in% 
                                              inTWAS]), dbs_query)
df_multi = data.frame(Term = NA,
                      Overlap = NA,
                      P.value = NA,
                      Adjusted.P.value = NA,
                      Odds.Ratio = NA,
                      Combined.Score = NA,
                      Genes = NA,
                      Ontology = NA)
for (i in 1:length(dbs_query)){
  print(i)
  ff = enriched[[i]]
  ff$Ontology = names(enriched)[i]
  if (i == 6){
    
    df_multi = rbind(df_multi,
                     subset(ff,
                            Adjusted.P.value < 
                              .1)[1:10,colnames(df_multi)])
  } else {
    df_multi = rbind(df_multi,
                     subset(ff,
                            round(Adjusted.P.value,2) <= 
                              .1)[1:10,colnames(df_multi)])
  }
  
}
df_multi = df_multi[complete.cases(df_multi),]

df_multi$Ontology = as.factor(df_multi$Ontology)
levels(df_multi$Ontology) = c('ChEA','BP','CC')
# df_multi$Term = sapply(strsplit(df_multi$Term,
#                                 '[(]'),
#                        function(x) x[1])


df_multi = df_multi[order(df_multi$Ontology,
                          df_multi$Odds.Ratio),]
df_multi$Term = factor(df_multi$Term,
                       levels = unique(df_multi$Term))

levels(df_multi$Term) = c('RUNX (mouse bone)',
                          'KLF4 (mouse MESC)',
                          'YY1 (human kidney embryo)',
                          'ESR1 (human T47D)',
                          'LDB1 (mouse BM-HSCs)',
                          'GATA4 (mouse HL-1)',
                          'CRX (mouse retina)',
                          'AF4 (human SEM)',
                          'BRD4 (human BCBL1)',
                          'ESR1 (human MCF-7)',
                          'Intracellular monoatomic cation homeostasis',
                          'Cell response to nutrient levels',
                          'Lipid biosynthesis',
                          'Fatty acid biosynthesis',
                          'Response to metal ion',
                          'Response to peptide hormone',
                          'Pos. regulation of fat cell differentiation',
                          'Membrane protein proteolysis',
                          'Membrane protein ectodomian proteolysis',
                          'Neg. reg. of CD8+ T-cell activation',
                          'Endosome membrane',
                          'Early endosome membrane',
                          'Phagocytic vesicle',
                          'Phagocytic vesicle membrane',
                          'Autolysosome',
                          'HFE-transferrin receptor complex',
                          'MHC class I protein complex')

df_multi = df_multi[order(df_multi$Ontology,
                          df_multi$Odds.Ratio),]
df_multi$Term = factor(as.character(df_multi$Term),
                       levels = unique(df_multi$Term))
                          
ora_multi = ggplot(df_multi[complete.cases(df_multi),], 
       aes(x=Term, 
           y=log2(Odds.Ratio))) +
  geom_segment(aes(x=Term, 
                   xend=Term, 
                   y=0, 
                   yend=log2(Odds.Ratio)), 
                color="grey") +
  geom_point(aes(color = Ontology),
             size = 1.2) +
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=8),
        plot.title = element_text(size = 8),
        legend.title=element_text(size=8),
        legend.text=element_text(size=8),
        strip.text = element_text(size=8),
        panel.spacing=unit(1, "lines"),
        panel.border = element_rect(color = NA,
                                    fill = NA, size = .1),
        legend.position=,
        legend.box = "vertical",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.key.size = unit(0.2, "cm"),
        strip.text.y = element_text(angle = 0)) +
  scale_color_manual(values = c(
    "#6A3D9A", # purple
    "#FF7F00", # orange
    "black", "gold1",
    "skyblue2", "#FB9A99", # lt pink
    "palegreen2",
    "#CAB2D6", # lt purple
    "#FDBF6F", # lt orange
    "gray70", "khaki2",
    "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
    "darkturquoise", "green1", "yellow4", "yellow3",
    "darkorange4", "brown"
  )) +
  facet_grid(Ontology~.,space='free',scales='free') +
  coord_flip() +
  guides(color = 'none') +
  ylab(expression(log[2]~"Enrichment Ratio")) +
  xlab('')
ggsave(filename = '/Users/abhattacharya3/OneDrive - Inside MD Anderson/Talks/DCCP Grand Rounds - March 2024/ORA_talk.pdf',
       plot = ora_multi,
       width = 8,
       height = 8)


firstrow = plot_grid(plot_grid(SigGene,
                     effN,
                     pLIGene+ guides(fill = 'none'),
                     rel_widths = c(1,.7,1),
                     labels = c('A','B','C'),nrow=1),
                     plot_grid(ppp,ora_multi,
                               rel_widths = c(.6,1),
                               labels = c('D','E')),
                     rel_heights = c(1,2.2),ncol=1,
                     label_size = 12)
ggsave(filename = '/Users/abhattacharya3/OneDrive - Inside MD Anderson/Fall 2024/PanCan IsoTWAS/Manuscript/Figures/Fig2.pdf',
       plot = firstrow,
       width = 7.5,
       height = 7.5)






multi_isogenes = unique(isotwas_multi$Gene)
isotwas_multi = isotwas_multi[order(isotwas_multi$Gene),]
isotwas_multi = subset(isotwas_multi,!Gene %in% twas$Gene)

ceres = as.data.frame(
  fread('/Users/abhattacharya3/Downloads/CRISPR_GeCKO_(DepMap,_CERES)_subsetted.csv')
)
multi_ceres = colMeans(ceres[,colnames(ceres) %in% 
                               c(unique(isotwas_multi$HGNC))])


