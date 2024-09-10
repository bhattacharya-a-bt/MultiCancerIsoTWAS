require(data.table)
require(topr)
ldblocks = 
  fread('/Users/abhattacharya3/OneDrive - Inside MD Anderson/Fall 2024/PanCan IsoTWAS/Results/pyrho_EUR_LD_blocks.bed')
for (f in 2:nrow(ldblocks)){
  
  if (ldblocks$chr[f] == ldblocks$chr[f-1]){
    if (ldblocks$start[f] != ldblocks$end[f-1]){
      
      ldblocks = rbind(ldblocks,
                      data.frame(chr = ldblocks$chr[f],
                                 start = ldblocks$end[f-1],
                                 end = ldblocks$start[f]))
      
    }
  }
  
}

manifest_df = data.frame(File = c('0_AdvProstate.subtypes.AUTO_Finalcleaned_MAF0.01.gz',
                                  '0_all.meta.gwas.bcfr.icogs.oncoarray.overall.AUTO_Finalcleaned_MAF0.01.gz',
                                  '0_Cimba_braca1_Breast_ERneg.AUTO_Finalcleaned_MAF0.01.gz',
                                  '0_Cimba_Ovarian_serous.AUTO_Finalcleaned_MAF0.01.gz',
                                  'Colorectal_Fernandez.tsv.gz',
                                  '0_CrossCancer_ENDOMETRIAL.AUTO_Finalcleaned_MAF0.01.gz',
                                  'Lung_Byun.tsv.gz',
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
manifest_df$Total = 
  file.path('/Users/abhattacharya3/OneDrive - Inside MD Anderson/Fall 2024/PanCan IsoTWAS/Results/GWAS/',
            manifest_df$File)

all_top_snps = data.frame(CHROM = c(),
                          POS = c(),
                          ID = c(),
                          REF = c(),
                          ALT = c(),
                          P = c(),
                          OR = c(),
                          AF = c(),
                          TWAS = c(),
                          isoTWAS = c(),
                          POS_hg38 = c(),
                          LDBlock = c(),
                          Method = c(),
                          Cancer = c(),
                          LDBlock_Method = c())
for (i in 1:nrow(manifest_df)){
  
  print(manifest_df$Cancer[i])
  gwas = fread(manifest_df$Total[i])
  if (manifest_df$Cancer[i] == 'UCEC'){
    gwas$MAF = gwas$Rsq_ave
  }
  
  if (i %in% c(5,7)){
    gwas$Beta = gwas$BETA
    gwas$Position = gwas$BP
    gwas$rsid = gwas$SNP
    gwas$MAF = NA
  }
    
  gwas_df = data.frame(CHROM = paste0('chr',gwas$CHR),
                       POS = gwas$Position,
                       ID = gwas$rsid,
                       REF = gwas$A2,
                       ALT = gwas$A1,
                       P = gwas$P,
                       OR = exp(gwas$Beta),
                       AF = gwas$MAF)
  
  top_snps = get_lead_snps(gwas_df,verbose = T,keep_chr = F)
  prefix = strsplit(manifest_df$File[i],'_MAF')[[1]][1]
  
  twas = fread(
    '/Users/abhattacharya3/OneDrive - Inside MD Anderson/Fall 2024/PanCan IsoTWAS/Results/TWAS_FineMap_PanCan_090524.tsv'
  )
  twas = twas[!duplicated(twas),]
  
  twas = subset(twas,in_cred_set == TRUE)
  twas$Indication = as.factor(twas$Indication)
  levels(twas$Indication) = c('BRCA (ER-)',
                              'BRCA (ER+)',
                              'BRCA',
                              'CRC',
                              'UCEC',
                              'Lung',
                              'LUAD',
                              'LUSC',
                              'OV',
                              'OV (ser)',
                              'PRCA (adv)',
                              'PRCA')
  twas = subset(twas,Indication == manifest_df$Cancer[i])
  
  
  isotwas = as.data.frame(vroom('isoTWAS_FineMap_PanCan_090524.tsv'))[,1:17]
  isotwas = isotwas[!duplicated(isotwas),]
  isotwas = subset(isotwas,in_cred_set == TRUE)
  isotwas$Indication = as.factor(isotwas$Indication)
  levels(isotwas$Indication) = c('BRCA (ER-)',
                              'BRCA (ER+)',
                              'BRCA',
                              'CRC',
                              'UCEC',
                              'Lung',
                              'LUAD',
                              'LUSC',
                              'OV',
                              'OV (ser)',
                              'PRCA (adv)',
                              'PRCA')
  isotwas = subset(isotwas,Indication == manifest_df$Cancer[i])
  
  top_snps$TWAS = F
  top_snps$isoTWAS = F
  top_snps$ID[is.na(top_snps$ID)] = paste(top_snps$CHROM[is.na(top_snps$ID)],
                                          top_snps$POS[is.na(top_snps$ID)],
                                          sep=':')
  
  if (!(i %in% c(5,7))){
  library(rtracklayer)
  path = '/Users/abhattacharya3/OneDrive - Inside MD Anderson/Fall 2024/PanCan IsoTWAS/Results/hg19ToHg38.over.chain'
  ch = import.chain(path)
  df.snp <- data.frame(chr=paste0('chr',top_snps$CHROM), 
                       start=top_snps$POS, 
                       end=top_snps$POS, 
                       score=1:nrow(top_snps),
                       id = top_snps$ID)
  
  gr  = makeGRangesFromDataFrame(df.snp, 
                                 ignore.strand=TRUE,
                                 keep.extra.columns=TRUE)
  cur19 = as.data.frame(unlist(liftOver(gr, ch)))
  
  genelocs = data.frame(ID = cur19$id,
                        chr = sapply(strsplit(as.character(cur19$seqnames),'r'),
                                     function(x) x[2]),
                        POS_hg38 = cur19$start,
                        right = cur19$end)
  top_snps = merge(top_snps,genelocs[,c('ID','POS_hg38')],by='ID')
  } else {
    top_snps$POS_hg38 = top_snps$POS
  }
  top_snps$LDBlock = ''
  for (t in 1:nrow(top_snps)){
    
    ld_this = subset(ldblocks,chr == paste0('chr',top_snps$CHROM[t]))
    ld_this = subset(ld_this,start < top_snps$POS[t] & end > top_snps$POS[t])
    if (nrow(ld_this) == 1){
      top_snps$LDBlock[t] = paste(ld_this$chr[1],
                                  ld_this$start[1],
                                  ld_this$end[1],
                                  sep=':')
    } else {
      top_snps$LDBlock[t] = top_snps$ID[t]
    }
    
  }
  
  
  for (t in 1:nrow(top_snps)){
    
    ld_this = subset(ldblocks,chr == paste0('chr',top_snps$CHROM[t]))
    ld_this = subset(ld_this,start < top_snps$POS[t] & end > top_snps$POS[t])
    if (nrow(ld_this) == 1){
      top_snps$LDBlock[t] = paste(ld_this$chr[1],
                                  ld_this$start[1],
                                  ld_this$end[1],
                                  sep=':')
    } else {
      top_snps$LDBlock[t] = top_snps$ID[t]
    }
    
  }
  
  print('Counting GWAS tags')
  for (q in 1:nrow(top_snps)){
    
    twas_sub = subset(twas,Chromosome == top_snps$CHROM[q])
    pos = top_snps$POS_hg38[q]
    twas_inGWAS = apply(twas_sub,1,function(x){
      
      return(ifelse(pos < as.numeric(x[6]) + 1e6 & pos > as.numeric(x[6]) - 1e6,
                    1,
                    0))
      
    })
    
    isotwas_sub = subset(isotwas,Chromosome == top_snps$CHROM[q])
    pos = top_snps$POS[q]
    isotwas_inGWAS = apply(isotwas_sub,1,function(x){
      
      return(ifelse(pos < as.numeric(x[6]) + 1e6 & pos > as.numeric(x[6]) - 1e6,
                    1,
                    0))
      
    })
    
    top_snps$TWAS[q] = (sum(twas_inGWAS) > 0)
    top_snps$isoTWAS[q] = (sum(isotwas_inGWAS) > 0)
    
  }
  top_snps = as.data.frame(top_snps)
  top_snps$Method = ''
  top_snps$Method = ifelse(top_snps$isoTWAS == TRUE & top_snps$TWAS == TRUE,
                           'Both',
                           top_snps$Method)
  top_snps$Method = ifelse(top_snps$isoTWAS == TRUE & top_snps$TWAS == FALSE,
                           'isoTWAS',
                           top_snps$Method)
  top_snps$Method = ifelse(top_snps$isoTWAS == FALSE & top_snps$TWAS == TRUE,
                           'TWAS',
                           top_snps$Method)
  top_snps$Method = ifelse(top_snps$isoTWAS == FALSE & top_snps$TWAS == FALSE,
                           'Neither',
                           top_snps$Method)
  top_snps$Cancer = manifest_df$Cancer[i]
  top_snps = top_snps[order(as.numeric(top_snps$CHROM),
                            as.numeric(top_snps$POS_hg38)),]
  top_snps$LDBlock_Method = paste(top_snps$LDBlock,top_snps$Method,sep='_')
  top_snps = top_snps[!duplicated(top_snps$LDBlock_Method),]
  all_top_snps = rbind(all_top_snps,
                       top_snps)
  
}

require(dplyr)
require(ggplot2)
summary_gwas = all_top_snps %>%
  group_by(Cancer) %>%
  summarise(Both = sum(Method == 'Both'),
            TWAS = sum(Method == 'TWAS'),
            isoTWAS = sum(Method == 'isoTWAS'),
            Total = sum(Method == 'Neither'))

summary_gwas$Total[summary_gwas$Cancer=='CRC'] =
  205 - summary_gwas$Both[summary_gwas$Cancer=='CRC'] - summary_gwas$TWAS[summary_gwas$Cancer=='CRC']- summary_gwas$isoTWAS[summary_gwas$Cancer=='CRC']

fwrite(summary_gwas,
       '/Users/abhattacharya3/OneDrive - Inside MD Anderson/Fall 2024/PanCan IsoTWAS/Results/ASHG/NumGWASLoci_afterFM.tsv',
       quote=F,
       row.names = F,
       col.names = T,
       sep='\t')
summary_gwas = fread('/Users/abhattacharya3/OneDrive - Inside MD Anderson/Fall 2024/PanCan IsoTWAS/Results/ASHG/NumGWASLoci_afterFM.tsv')
 
both_sum = sum(summary_gwas$Both)
twas_sum = sum(summary_gwas$TWAS)
isotwas_sum = sum(summary_gwas$isoTWAS)
total_sum = sum(summary_gwas[,c(2:5)])

summary_melt = reshape2::melt(summary_gwas,
                              measure.vars =
                                c('Both','isoTWAS',
                                  'TWAS','Total'))
colnames(summary_melt) = c('Cancer','Method','Number')
summary_melt$Method = factor(summary_melt$Method,
                             c('Total',
                               'isoTWAS',
                               'TWAS',
                               'Both'))
levels(summary_melt$Method) = c('Total','Only tagged\nby isoTWAS',
                                'Only tagged\nby TWAS','Both')
summary_melt$Cancer = factor(summary_melt$Cancer,
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

require(wesanderson)
gwasloci = ggplot(data = summary_melt,
                  aes(x = Cancer,
                      y = Number,
                      fill = Method)) +
  geom_col(color = 'black') +
  theme_minimal() +
  theme(axis.text=element_text(size=8.5),
        axis.title=element_text(size=9),
        plot.title = element_text(size = 9),
        legend.title=element_text(size=9),
        legend.text=element_text(size=8),
        strip.text = element_text(size=9),
        panel.spacing=unit(1, "lines"),
        panel.border = element_rect(color = NA,
                                    fill = NA, linewidth = .1),
        legend.position='bottom',
        legend.box = "horizontal",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.key.size = unit(0.2, "cm"),
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust=1),
        legend.position.inside = c(0.5, 0.8)) +
  scale_fill_manual(values = c('grey',
                               wes_palette("Darjeeling1")[c(2,1,3)])) +
  xlab('') +
  ylab('# tagged independent\nGWAS loci') +
  labs(fill = 'GWAS\nloci')
a = ggpubr::as_ggplot(ggpubr::get_legend(gwasloci))
ggsave(filename = '/Users/abhattacharya3/OneDrive - Inside MD Anderson/Talks/DCCP Grand Rounds - March 2024/GWASCancerPlot.pdf',
       plot = gwasloci,
       width = 5.56,
       height = 5.4)

require(readxl)
mesc = read_xlsx('/Users/abhattacharya3/OneDrive - Inside MD Anderson/Fall 2024/PanCan IsoTWAS/Manuscript/Tables/SupplementalTable_MESC.xlsx')
mesc$Label = ''
mesc$Label = ifelse(round(mesc$`FDR-adjusted P`,2) < .05 &
                      mesc$`FDR-adjusted P` >= 0.01,'*',mesc$Label)
mesc$Label = ifelse(mesc$`FDR-adjusted P` < .01 &
                      mesc$`FDR-adjusted P` >= 0.005,'**',mesc$Label)
mesc$Label = ifelse(mesc$`FDR-adjusted P` < .005,'***',mesc$Label)
mesc$Cancer = factor(mesc$Cancer,
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

pp=ggplot(data = mesc,
          aes(x = Cancer,
              y = `h2med/h2`,
              fill = Feature)) +
  geom_bar(stat = "identity", 
           alpha = 0.7, 
           position = position_dodge(width = .9),
           color = 'black') +
  geom_linerange(aes(ymin = `h2med/h2` - `h2med/h2 (SE)`,
                     ymax = `h2med/h2` + `h2med/h2 (SE)`),
                 position = position_dodge2(width = .9),
                 size = .5) +
  xlab('') +
  ylab(expression(h[med]^2/h^2)) +
  theme(axis.text=element_text(size=8.5),
        axis.title=element_text(size=9),
        plot.title = element_text(size = 9),
        legend.title=element_text(size=9),
        legend.text=element_text(size=8.5),
        strip.text = element_text(size=9),
        panel.spacing=unit(1, "lines"),
        panel.border = element_rect(color = NA,
                                    fill = NA, size = .1),
        legend.position=,
        legend.box = "horizontal",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.key.size = unit(0.2, "cm"),
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust=1)) +
  scale_fill_manual(values = wes_palette("Darjeeling1")) +
  geom_text(aes(y = `h2med/h2` + `h2med/h2 (SE)`,
                label = Label),
            position = position_dodge2(width = .9),
            hjust = -.2,
            vjust = .8,
            angle = 90) +
  labs(fill = 'Feature') +
  ylim(c(-.25,.5))

all_top_snps = all_top_snps[order(all_top_snps$CHROM,
                                  all_top_snps$POS),]


require(karyoploteR)
both = toGRanges(data.frame(
  chr = paste0('chr',all_top_snps$CHROM[all_top_snps$Method == 'Both']),
  start = all_top_snps$POS_hg38[all_top_snps$Method == 'Both'],
  end = all_top_snps$POS_hg38[all_top_snps$Method == 'Both']
))

twas = toGRanges(data.frame(
  chr = paste0('chr',all_top_snps$CHROM[all_top_snps$Method == 'TWAS']),
  start = all_top_snps$POS_hg38[all_top_snps$Method == 'TWAS'],
  end = all_top_snps$POS_hg38[all_top_snps$Method == 'TWAS']
))

isotwas = toGRanges(data.frame(
  chr = paste0('chr',all_top_snps$CHROM[all_top_snps$Method == 'isoTWAS']),
  start = all_top_snps$POS_hg38[all_top_snps$Method == 'isoTWAS'],
  end = all_top_snps$POS_hg38[all_top_snps$Method == 'isoTWAS']
))

neither = toGRanges(data.frame(
  chr = paste0('chr',all_top_snps$CHROM[all_top_snps$Method == 'Neither']),
  start = all_top_snps$POS_hg38[all_top_snps$Method == 'Neither'],
  end = all_top_snps$POS_hg38[all_top_snps$Method == 'Neither']
))

ggsave(filename = '/Users/abhattacharya3/OneDrive - Inside MD Anderson/Talks/DCCP Grand Rounds - March 2024/MESCCancerPlot.pdf',
       plot = pp,
       width = 6,
       height = 2.5)

clpp_all = readRDS('/Users/abhattacharya3/OneDrive - Inside MD Anderson/Fall 2024/PanCan IsoTWAS/Results/QTL/CLPP/results_eCAVIAR.RDS')



ddd = rbind(data.frame(Indication = clpp_all$isotwas_df$Cancer,
                       Gene = clpp_all$isotwas_df$Gene,
                       HGNC = clpp_all$isotwas_df$HGNC,
                       Transcript = clpp_all$isotwas_df$Isoform,
                       CLPP = clpp_all$isotwas_df$CLPP_isoform,
                       QTL_P = clpp_all$isotwas_df$isoQTL_P,
                       GWAS_P = clpp_all$isotwas_df$GWAS_P,
                       Method = 'isoTWAS'),
            data.frame(Indication = clpp_all$twas_df$Cancer,
                       Gene = clpp_all$twas_df$Gene,
                       HGNC = clpp_all$twas_df$HGNC,
                       Transcript = '',
                       CLPP = clpp_all$twas_df$CLPP_gene,
                       QTL_P = clpp_all$twas_df$eQTL_P,
                       GWAS_P = clpp_all$twas_df$GWAS_P,
                       Method = 'TWAS'))
ddd$Indication = as.factor(ddd$Indication)
levels(ddd$Indication) = c('BRCA (ER-)',
                           'BRCA (ER+)',
                           'BRCA',
                           'CRC',
                           'UCEC',
                           'Lung',
                           'LUAD',
                           'LUSC',
                           'OV',
                           'PRCA (adv)',
                           'PRCA')
ddd$Indication = factor(as.character(ddd$Indication),
                        levels = c('BRCA',
                                   'BRCA (ER+)',
                                   'BRCA (ER-)',
                                   'CRC',
                                   'Lung',
                                   'LUAD',
                                   'LUSC',
                                   'OV',
                                   'PRCA',
                                   'PRCA (adv)',
                                   'UCEC'))
ddd$Method = factor(ddd$Method,
                    levels = c('TWAS','isoTWAS'))
levels(ddd$Method) = c('Gene','Isoform')
fwrite(ddd,
       '/Users/abhattacharya3/OneDrive - Inside MD Anderson/Fall 2024/PanCan IsoTWAS/Results/CLPP_SuppTable.tsv',
       sep='\t',
       col.names = T,
       row.names=F,
       quote=F)

ddd = subset(ddd,GWAS_P < 5e-8 & QTL_P < 1e-6)
ddd_sum = ddd %>%
  group_by(Indication,Method) %>%
  summarise(Min = quantile(CLPP,.1),
            Low = quantile(CLPP,.25),
            Median = quantile(CLPP,.5),
            High = quantile(CLPP,.75),
            Max = quantile(CLPP,.9),
            Prop = mean(CLPP > .01))

viol=ggplot(data = ddd_sum,
          aes(x = Indication,
              fill = Method)) +
  geom_boxplot(aes(middle = Median,
                   lower = Low,
                   upper = High,
                   ymin = Min,
                   ymax = Max),
               stat = 'identity',
               width = .5,
               color = 'black') +
  xlab('') +
  ylab('CLPP') +
  guides(fill = guide_legend(position='inside')) +
  theme(axis.text=element_text(size=8.5),
        axis.title=element_text(size=9),
        plot.title = element_text(size = 9),
        legend.title=element_text(size=9),
        legend.text=element_text(size=8.5),
        strip.text = element_text(size=9),
        panel.spacing=unit(1, "lines"),
        panel.border = element_rect(color = NA,
                                    fill = NA, size = .1),
        legend.position.inside = c(0.2, 0.75),
        legend.box = "horizontal",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.key.size = unit(0.2, "cm"),
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust=1)) +
  scale_fill_manual(values = wes_palette("Darjeeling1")) +
  scale_y_sqrt() + labs(fill = 'Feature')


prop=ggplot(data = ddd_sum,
            aes(x = Indication,
                y = Prop,
                fill = Method)) +
  geom_col(position = position_dodge2(width = .9)) +
  xlab('') +
  ylab('Proportion of loci with CLPP > 0.01') +
  guides(fill = guide_legend(position='inside')) +
  theme(axis.text=element_text(size=8.5),
        axis.title=element_text(size=9),
        plot.title = element_text(size = 9),
        legend.title=element_text(size=9),
        legend.text=element_text(size=8.5),
        strip.text = element_text(size=9),
        panel.spacing=unit(1, "lines"),
        panel.border = element_rect(color = NA,
                                    fill = NA, size = .1),
        legend.position.inside = c(0.1, 0.9),
        legend.box = "horizontal",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.key.size = unit(0.2, "cm"),
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust=1)) +
  scale_fill_manual(values = wes_palette("Darjeeling1")) +
  labs(fill = 'Feature')
ggsave(filename = '/Users/abhattacharya3/OneDrive - Inside MD Anderson/Fall 2024/PanCan IsoTWAS/Manuscript/Figures/clpp_supp.pdf',
       plot = prop,
       width = 6.5,
       height = 6)

require(cowplot)
secondrow = plot_grid(
  plot_grid(
    plot_grid(
      gwasloci + guides(fill = 'none'),
      a,
      rel_heights = c(1,.1),
      nrow=2,
      labels = 'A',
      label_size = 12),
    viol,
    rel_widths = c(1,1),
    nrow=1,
    labels=c('','B'),
    label_size = 12),
  pp,
  rel_heights = c(1,1),
  labels = c('','C'),
  nrow=2,
  label_size = 12
  )
ggsave(filename = '/Users/abhattacharya3/OneDrive - Inside MD Anderson/Fall 2024/PanCan IsoTWAS/Manuscript/Figures/Fig3.pdf',
       plot = secondrow,
       width = 6.5,
       height = 6)


setwd('/Users/abhattacharya3/OneDrive - Inside MD Anderson/Fall 2024/PanCan IsoTWAS/Results/QTL/CLPP/eCAVIAR/')
a = list.files()
for (i in a){
  
  clpp_all = readRDS(i)
  out_df = data.frame(Cancer = clpp_all$Cancer,
                      Tissue = clpp_all$Tissue,
                      Gene = clpp_all$Gene,
                      Transcript = clpp_all$Transcript,
                      SNP = clpp_all$eQTL$SNP,
                      isoQTL_P = clpp_all$isoQTL$P,
                      eQTL_P = clpp_all$eQTL$P)
  
  GWAS = data.frame(SNP = clpp_all$GWAS$rsid,
                    GWAS_P = clpp_all$GWAS$P)
  out_df = merge(out_df,GWAS,by='SNP')
  
  colnames(clpp_all$Coloc_Gene)[1] =
    colnames(clpp_all$Coloc_Tx)[1] = 'SNP'
  
  out_df = merge(merge(out_df,clpp_all$Coloc_Gene,by='SNP'),
                 clpp_all$Coloc_Tx,by='SNP')
  out_df = merge(out_df,clpp_all$LD,by='SNP')
  colnames(out_df)[9:12] = c('Probability eQTL in causal set',
                             'eQTL CLPP',
                             'Probability isoQTL in causal set',
                             'isoQTL CLPP')
  fwrite(out_df,
         '/Users/abhattacharya3/OneDrive - Inside MD Anderson/Fall 2024/PanCan IsoTWAS/Manuscript/RawResults/QTLResults_SuppData3.tsv',
         sep = '\t',
         quote=F,
         append=T,
         row.names=F)
  
}


