rm(list=ls())
setwd('/Users/abhattacharya3/OneDrive - Inside MD Anderson/Fall 2024/PanCan IsoTWAS/Results/QTL/')
fff = list.files()
fff = fff[!grepl('.pdf',fff)]
fff = fff[grepl('ER-',fff)]


require(vroom)
regulome = vroom('/Users/abhattacharya3/Downloads/ENCFF250UJY.tsv')
colnames(regulome)[4] = 'SNP'

inGene = data.table::fread('/Users/abhattacharya3/OneDrive - Inside MD Anderson/Fall 2024/PanCan IsoTWAS/Results/isoTWAS/InGeneCandidates.tsv')

inGene$TissueAll = as.factor(inGene$Tissue)
levels(inGene$TissueAll)

bim = data.table::fread('/Users/abhattacharya3/OneDrive - Inside MD Anderson/Fall 2024/PanCan IsoTWAS/Results/isoTWAS/GTEx.WGS.838.passOnly.geno0.05.hwe0.00001.dbsnp.SNPsOnly.NoAmbig.LDREF.bim')
colnames(bim) = c('Chromosome','SNP','V3','Position','A1','A2')

require(data.table)
require(ggplot2)
require(bigsnpr)
snps = snp_attach(snp_readBed2('/Users/abhattacharya3/OneDrive - Inside MD Anderson/Fall 2024/PanCan IsoTWAS/Results/isoTWAS/GTEx.WGS.838.passOnly.geno0.05.hwe0.00001.dbsnp.SNPsOnly.NoAmbig.LDREF.bed',
                               backingfile = tempfile()))

require(LDlinkR)
require(susieR)
require(coloc)

f = 'LAMC1_ColonTransverse.RDS'
all_res = readRDS(f)







df_plot = rbind(
  data.frame(Feature = 'LAMC1',
             SNP = all_res$eQTL$SNP,
             P = all_res$eQTL$P),
  data.frame(Feature = all_res$isoQTL$Feature,
             SNP = all_res$isoQTL$SNP,
             P = all_res$isoQTL$P)
)
df_plot = merge(bim[,c('SNP','Position')],
                df_plot,by='SNP')

df_plot = rbind(df_plot,data.frame(Feature = all_res$Cancer,
                                   SNP = all_res$GWAS$rsid,
                                   P = all_res$GWAS$P,
                                   Position = all_res$GWAS$Position_hg38))

inGene_this = subset(inGene,HGNC == all_res$Gene)

pli = unique(inGene_this$pLI)
pli = pli[!is.na(pli)]

df_plot = subset(df_plot,Position < inGene_this$End[1] + 1e4 &
                   Position > inGene_this$Start[1] - 1e4)

snps_this = snp_attach(subset(snps,
                              ind.col = which(snps$map$marker.ID %in% df_plot$SNP),
                              backingfile = tempfile()))


sss = subset(df_plot,Feature == 'ENST00000466964.1')
lead_snp = sss$SNP[which.min(sss$P)]

my_proxies <- LDproxy(snp = lead_snp, 
                      pop = "CEU", 
                      r2d = "r2", 
                      token = 'f199984caa85',
                      win_size = 3e5)
ld_df = my_proxies
ld_df$SNP = my_proxies$RS_Number
ld_df$LD_all = ld_df$R2

df_plot = merge(df_plot,ld_df,by='SNP')
df_plot = merge(df_plot,regulome,by='SNP')


df_plot$Label = ifelse(df_plot$SNP == lead_snp,
                       lead_snp,'')
df_plot$LDBlock = cut(df_plot$LD_all, 
                      breaks = c(0,.2,.4,.6,.8,1),
                      include.lowest=TRUE)
df_plot$Sig = ifelse(df_plot$Feature == all_res$Cancer,
                     -log10(5e-8),
                     -log10(1e-6))

df_plot = df_plot[!is.na(df_plot$SNP),]




require(ggrepel)
top='ENST00000466964.1'

df_plot$Feature = factor(as.character(df_plot$Feature),
                         levels = c(all_res$Cancer,
                                    'LAMC1',
                                    top,
                                    all_res$Transcript[all_res$Transcript != top]))

gwas_this = data.frame(CHROM = 'chr1',
                       POS = df_plot$Position[df_plot$Feature == 'Colorectal'],
                       ID =  df_plot$SNP[df_plot$Feature == 'Colorectal'],
                       P = df_plot$P[df_plot$Feature == 'Colorectal'],
                       R2 = df_plot$LD_all[df_plot$Feature == 'Colorectal'])

df_plot$Strip = df_plot$Feature
levels(df_plot$Strip) = c('Colorectal',
                          'LAMC1 (no TWAS model)',
                          'ENST00000466964.1 (Z = -6.0, effect)',
                          'ENST00000258341.5 (Z = -6.0)',
                          'ENST00000479499.1 (Z = -0.5)',
                          'ENST00000478064.1 (Z = 5.5)',
                          'ENST00000495918.1 (Z = 8.1)')
locuszoom = ggplot(data = df_plot,
                   aes(x = Position,
                       y = -log10(P),
                       label = Label,
                       color = LDBlock)) +
  facet_wrap(~Strip,ncol = 1,scales = 'free_y') +
  geom_hline(aes(yintercept = Sig),
             linetype = 2,
             color = 'black',
             size = .3) +
  geom_vline(xintercept = c(inGene_this$Start[1],
                            inGene_this$End[1]),
             linetype = 2,
             color = 'black',
             size = .3) +
  geom_point(aes(size = ifelse(Label != '',
                               'not lead','lead'),
                 shape = ifelse(Label != '',
                                'not lead','lead'))) +
  theme_minimal() +
  theme(axis.text=element_text(size=6),
        axis.title=element_text(size=8),
        plot.title = element_text(size = 8),
        legend.title=element_text(size=8),
        legend.text=element_text(size=7),
        strip.text = element_text(size=8),
        panel.spacing=unit(1, "lines"),
        panel.border = element_rect(color = NA,
                                    fill = NA, size = .1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(size = 5),
        axis.line = element_line(colour = "black"),
        legend.margin = margin(0, 0, 0, 0),
        legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"),
        legend.direction = 'horizontal',legend.key.spacing = unit(0,'mm'),
        legend.location = 'bottom') +
  scale_size_discrete("level", range=c(.8,2)) +
  guides(size = 'none',
         shape = 'none') +
  xlab(paste0('Position on Chromosome ',all_res$GWAS$CHR[1])) +
  ylab(expression(-log[10]~"P-value")) +
  labs(color = 'LD') +
  scale_color_manual(values = c('navyblue',
                                'lightblue',
                                'darkgreen',
                                'orange',
                                'red'))
lz_leg = ggpubr::as_ggplot(ggpubr::get_legend(locuszoom))


library(magrittr)
library(dplyr)
library(ggplot2)
library(ggtranscript)
setwd('/Users/abhattacharya3/OneDrive - Inside MD Anderson/Fall 2024/PanCan IsoTWAS/Results/QTL/')

gtf <- rtracklayer::import('gencode.v38.annotation.gtf')
gtf$gene_id = sapply(strsplit(gtf$gene_id,'[.]'),
                     function(x) x[1])

gtf = subset(gtf, 
             gene_id %in% c('ENSG00000135862'))
tx_list = unique(df_plot$Feature)

df_annotation = data.frame(gene_name = c('LAMC1'),
                           N_tx = NA,
                           N_tx_modeled = NA,
                           identifiedTx = c('ENSG00000135862'))
for (i in 1:nrow(df_annotation)){
  
  gtf_this = as.data.frame(subset(gtf,
                                  gene_name == df_annotation$gene_name[i] &
                                    type == 'transcript'))
  df_annotation$N_tx[i] = nrow(gtf_this)
  df_annotation$N_tx_modeled[i] = sum(gtf_this$transcript_id %in%
                                        tx_list)
  
}

df_annotation$StripName = 
  paste0(df_annotation$gene_name,' (',
         df_annotation$N_tx, ' total isoforms, ',
         df_annotation$N_tx_modeled,' in isoTWAS model)')

gtf <- subset(gtf,transcript_id %in% gtf$transcript_id) %>% dplyr::as_tibble()
gtf$Color = as.factor(ifelse(gtf$transcript_id %in% df_plot$Feature,
                             'Yes','No'))
gtf$Color = relevel(gtf$Color,ref='Yes')
gtf = merge(gtf,df_annotation,by='gene_name')
gtf_exons <- gtf %>% 
  dplyr::filter(type == "exon")

g = 'LAMC1'
ggg = subset(gtf_exons,gene_name == g)

ggg$transcript_id = as.factor(ggg$transcript_id)
ggg$transcript_id = relevel(ggg$transcript_id,
                            ref = top)
ggg = subset(ggg,transcript_id %in% df_plot$Feature)
ggg$transcript_id = droplevels(ggg$transcript_id)
ggg$Color = ifelse(ggg$transcript_id == top,
                   'Yes','No')

tx_plot = ggg %>%
  ggplot(aes(
    xstart = start,
    xend = end,
    y = transcript_id
  )) +
  geom_range(
    aes(fill = Color)
  ) +
  geom_intron(
    data = to_intron(ggg, "transcript_id"),
    aes(strand = strand)
  ) +
  #facet_wrap(~StripName,scales = 'free',ncol = 2) + 
  theme_minimal() +
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=8),
        plot.title = element_text(size = 10),
        legend.title=element_text(size=8),
        legend.text=element_text(size=7),
        strip.text = element_text(size=8),
        panel.spacing=unit(1, "lines"),
        panel.border = element_rect(color = NA,
                                    fill = NA, size = .1),
        legend.position="bottom",
        legend.box = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  xlab(paste0('Position on Chromosome ',all_res$GWAS$CHR[1])) +
  ylab('Isoform') +
  labs(fill = 'Prioritized by isoTWAS') + guides(fill='none') +
  geom_vline(xintercept = df_plot$Position[df_plot$Feature == 'Colorectal' & 
                                             df_plot$probability_score > .9 &
                                             df_plot$Position > inGene_this$Start[1] &
                                             df_plot$Position < inGene_this$End[1] & 
                                             df_plot$LD_all > .8],
             linetype = 1,
             color = 'grey',
             size = .2) +
  geom_vline(xintercept = df_plot$Position[df_plot$Feature == top & 
                                             df_plot$LD_all > .8 &
                                             df_plot$Position > inGene_this$Start[1]-100 &
                                             df_plot$Position < inGene_this$End[1]+100 &
                                             df_plot$probability_score > .9],
             linetype = 2,
             color = 'red',
             size = .3) +
  scale_x_continuous(breaks = c(183100000,183125000,183145000),limits = c(183100000,
                                                                          183150000))

df_annotate = data.frame(Annotation = rep(c('GWAS SNP w/ P < 5e-8',paste0('isoQTL in LD = 1 with ',lead_snp)),
                                          each = 20),
                         Number = c(rnorm(20),rnorm(20)))
ddd = ggplot(data = df_annotate,
             aes(x = Annotation,
                 y = Number,
                 color = Annotation,
                 linetype = Annotation)) +
  geom_line(linetype = 2) +
  scale_color_manual(values = c('grey','red')) +
  scale_linetype_manual(values = c(1,2)) +
  theme_minimal() +
  theme(axis.text=element_text(size=7.5),
        axis.title=element_text(size=8),
        plot.title = element_text(size = 8),
        legend.title=element_text(size=8),
        legend.text=element_text(size=7),
        strip.text = element_text(size=8),
        panel.spacing=unit(1, "lines"),
        panel.border = element_rect(color = NA,
                                    fill = NA, size = .1),
        legend.position="bottom",
        legend.box = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
a = ggpubr::as_ggplot(ggpubr::get_legend(ddd))


forest_gwas = subset(all_res$GWAS,rsid %in% 
                       df_plot$SNP[df_plot$Feature == top & 
                                     df_plot$LD_all == 1 &
                                     df_plot$Position > inGene_this$Start[1]-100 &
                                     df_plot$Position < inGene_this$End[1]+100 &
                                     grepl('missen',df_plot$Function)])
forest_gwas$Feature = 'Colorectal'
forest_gwas = forest_gwas[,c('rsid','Feature','Beta','SE','Position_hg38')]
colnames(forest_gwas) = c('SNP','Feature','Beta','SE','Position')

forest_eqtl = subset(all_res$eQTL,SNP %in% 
                       df_plot$SNP[df_plot$Feature == top & 
                                     df_plot$LD_all == 1 &
                                     df_plot$Position > inGene_this$Start[1]-100 &
                                     df_plot$Position < inGene_this$End[1]+100 &
                                     grepl('missen',df_plot$Function)])

forest_eqtl = forest_eqtl[,c('SNP','Feature','Beta','SE')]
forest_eqtl = merge(forest_eqtl,bim[,c('SNP','Position')],
                    by='SNP')

forest_isoqtl = subset(all_res$isoQTL,SNP %in% 
                         df_plot$SNP[df_plot$Feature == top & 
                                       df_plot$LD_all == 1 &
                                       df_plot$Position > inGene_this$Start[1]-100 &
                                       df_plot$Position < inGene_this$End[1]+100 &
                                       grepl('missen',df_plot$Function)])

forest_isoqtl = forest_isoqtl[,c('SNP','Feature','Beta','SE')]
forest_isoqtl = merge(forest_isoqtl,bim[,c('SNP','Position')],
                      by='SNP')

forest_df = rbind(forest_gwas,
                  forest_eqtl,
                  forest_isoqtl)
forest_df$Feature = factor(forest_df$Feature,
                           c('Colorectal',
                             'ENSG00000135862',
                             top,
                             'ENST00000258341.5',
                             'ENST00000479499.1',
                             'ENST00000478064.1',
                             'ENST00000495918.1'))
levels(forest_df$Feature)[2] = 'LAMC1'
forest_df$Lower = forest_df$Beta - abs(qnorm(1e-6))*forest_df$SE
forest_df$Upper = forest_df$Beta + abs(qnorm(1e-6))*forest_df$SE
forest_df$SNP = paste0(forest_df$SNP,' (missense SNP)')

fplot = ggplot(data = forest_df,
               aes(x = Feature,
                   y = Beta,
                   color = Feature)) +
  geom_point(size=2) +
  theme_minimal() +
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=8),
        plot.title = element_text(size = 10),
        legend.title=element_text(size=8),
        legend.text=element_text(size=7),
        strip.text = element_text(size=8),
        panel.spacing=unit(1, "lines"),
        panel.border = element_rect(color = NA,
                                    fill = NA, size = .1),
        legend.position="bottom",
        legend.box = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  guides(color = 'none',size='none') +
  xlab('') +
  ylab('Effect') +
  geom_errorbar(aes(ymin = Lower,
                    ymax = Upper),
                width = .2) +
  scale_color_manual(values = c('black',
                                RColorBrewer::brewer.pal(n=3,'Set1')[1],
                                scales::hue_pal()(4)[3],
                                rep(scales::hue_pal()(4)[1],
                                    5))) +
  geom_hline(yintercept = 0,
             color = 'grey',
             linetype = 2) +
  coord_flip() + facet_wrap(~SNP,ncol=1)

require(cowplot)
total = plot_grid(plot_grid(locuszoom + guides(color = 'none'),
                            lz_leg,
                            nrow = 2,
                            rel_heights = c(1,.05)),
                  plot_grid(tx_plot,
                            a,
                            fplot,
                            rel_heights = c(.9,.1,1),
                            labels=c('B','','C'),
                            nrow=3),
                  nrow = 1,
                  rel_widths = c(.75,1),
                  labels = c('A',''))

ggsave(plot = total,
       filename = '/Users/abhattacharya3/OneDrive - Inside MD Anderson/Fall 2024/PanCan IsoTWAS/Manuscript/Figures/Fig_LAMC1.pdf',
       height = 7,
       width = 7)

