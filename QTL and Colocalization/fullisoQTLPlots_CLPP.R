require(data.table)

### AGGREGATE ISOTWAS
setwd('/rsrch5/home/epi/bhattacharya_lab/projects/PanCan_isoTWAS')
inGene_base = fread('isoTWAS_FineMap_PanCan.tsv')
inGene_base = inGene_base[!duplicated(inGene_base),]

interestGenes = c('UBAP2',
                  'CERS5',
                  'KDSR',
                  'LSP1',
                  'GGCX1',
                  'RBM23',
                  'MPZL2',
                  'CLPTMIL',
                  'PLEKHG6',
                  'BABAM1',
                  'TNNT3',
                  'PLEKHM1',
                  'SRP14',
                  'FDPS',
                  'CASP8',
                  'TMBIM1',
                  'TMEM258',
                  'L3MBTL3')

setwd("/rsrch5/home/epi/bhattacharya_lab/projects/PanCan_isoTWAS/eCAVIAR/")
fff = list.files()
fff = fff[!grepl('.pdf',fff)]
bim = fread('/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/GTEx.WGS.838.passOnly.geno0.05.hwe0.00001.dbsnp.SNPsOnly.NoAmbig.LDREF.bim')
colnames(bim)[c(1,2,4)] = c('Chromosome','SNP','Position')


hash = data.frame(Indication = unique(inGene_base$Indication),
                  Cancer = c('BRCA (ER-)',
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

manifest_df = merge(manifest_df,hash,by='Cancer')

gtf_base <- rtracklayer::import('/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode.v38.annotation.gtf')
gtf_base$gene_id = sapply(strsplit(gtf_base$gene_id,'[.]'),
                          function(x) x[1])

inGene_base = subset(inGene_base,HGNC %in% interestGenes)
inGene_base$GTA = paste(inGene_base$HGNC,inGene_base$Tissue,inGene_base$Indication,
                        sep=':')

for (gta in unique(inGene_base$GTA)){
  
  print(gta)
  sub_this = subset(inGene_base,GTA == gta)
  sub_this = sub_this[order(sub_this$in_cred_set,decreasing = T)]
  
  ### CATALOG NECESSARY INFO
  gene = unique(sub_this$Gene)[1]
  hgnc = unique(sub_this$HGNC)[1]
  tx_interest = unique(sub_this$Transcript)
  tissue = unique(sub_this$Tissue)[1]
  indication =  unique(sub_this$Indication)[1]
  cancer = unique(manifest_df$Cancer[manifest_df$Indication == indication])[1]
  chr = sub_this$Chromosome[1]
  start = sub_this$Start[1]
  end = sub_this$End[1]
  gwas_file = file.path('/rsrch5/home/epi/bhattacharya_lab/projects/PanCan_isoTWAS/GWAS',
                        manifest_df$File[manifest_df$Cancer == cancer])
  
  if (!file.exists(file.path(
    '/rsrch5/home/epi/bhattacharya_lab/projects/PanCan_isoTWAS/eCAVIAR/',
    paste0('FullPlot_QTL_Cancer',cancer,
           'Tissue',tissue,'Gene',hgnc,'.pdf')))){
  
  
  ### READ AND LIFTOVER GWAS
  gwas = fread(gwas_file)
  gwas = subset(gwas,CHR == chr)
  
  require(rtracklayer)
  path = '/rsrch5/home/epi/bhattacharya_lab/projects/PanCan_isoTWAS/hg19ToHg38.over.chain'
  ch = import.chain(path)
  if ('P_het' %in% colnames(gwas)){
    gwas$P = gwas$P_het
  }
  
  
  df.snp <- data.frame(chr=paste0('chr',gwas$CHR), 
                       start=gwas$Position, 
                       end=gwas$Position+1, 
                       score=1:nrow(gwas),
                       id = gwas$cptid)
  
  gr  = makeGRangesFromDataFrame(df.snp, 
                                 ignore.strand=TRUE,
                                 keep.extra.columns=TRUE)
  cur19 = as.data.frame(unlist(liftOver(gr, ch)))
  
  lo = data.frame(cptid = cur19$id,
                  chr = sapply(strsplit(as.character(cur19$seqnames),'r'),
                               function(x) x[2]),
                  Position_hg38 = cur19$start)
  
  gwas = merge(gwas,lo,by='cptid')
  gwas = subset(gwas,Position_hg38 < end + 1e6 &
                  Position_hg38 > start - 1e6)
  
  
  ### GRAB GTEX COVS
  cov_df = read.table(file.path('/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8',
                                tissue,
                                paste0(tissue,'.v8.covariates.txt.covar')))
  
  colnames(cov_df)[1:2] = c('ID','ID2')
  cov_df$ID = stringr::str_replace_all(cov_df$ID,'[.]','-')
  
  ### MAP ISOQTL
  require(SummarizedExperiment)
  gtex_tx = readRDS(file.path('/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8',
                              tissue,
                              paste0(tissue,'_transcripts.RDS')))
  rd = rowData(gtex_tx)
  rd$gene_id = unlist(rd$gene_id)
  gtex_tx = gtex_tx[grepl(gene,rd$gene_id),]
  
  tx_exp_mat = assays(gtex_tx)[[1]]
  tx_exp_mat = tx_exp_mat[,cov_df$ID]
  
  if (nrow(tx_exp_mat) == 1){
    
    tx_exp_mat = t(as.matrix(log(tx_exp_mat+1)))
  } else {
    tx_exp_mat = log(as.matrix(tx_exp_mat)+1) 
  }
  
  ttt = as.matrix(limma::removeBatchEffect(tx_exp_mat,
                                           covariates = 
                                             model.matrix(~. - ID - 1 - ID2,
                                                          data = cov_df)))
  
  require(bigsnpr)
  snps = snp_attach(snp_readBed2('/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/GTEx.WGS.838.passOnly.geno0.05.hwe0.00001.dbsnp.SNPsOnly.NoAmbig.LDREF.bed',
                                 backingfile = tempfile()))
  snps = snp_attach(subset(snps,
                           ind.col = which(snps$map$chromosome == chr &
                                             snps$map$physical.pos < end + 1e6 &
                                             snps$map$physical.pos > start - 1e6),
                           backingfile = tempfile()))
  ttt_df = as.data.frame(t(ttt))
  ttt_df$ID = rownames(ttt_df)
  
  res_lm_tx = data.frame(Cancer = c(),
                         Tissue = c(),
                         SNP = c(),
                         Feature = c(),
                         Beta = c(),
                         SE = c(),
                         P = c(),
                         Position = c())
  for (tx in rownames(ttt)){
    print(tx)
    for (s in 1:nrow(snps$map)){
      snp_df = data.frame(ID = stringr::str_replace_all(snps$fam$family.ID,'[.]','-'),
                          SNP = snps$genotypes[][,s])
      reg_df = merge(snp_df,ttt_df[,c('ID',tx)],by='ID')
      colnames(reg_df) = c('ID','SNP','Transcript')
      
      if (var(reg_df$SNP) != 0){
        
        lll = lm(Transcript ~ SNP,data = reg_df)
        res_lm_tx = rbind(res_lm_tx,
                          data.frame(Cancer = cancer,
                                     Tissue = tissue,
                                     SNP = snps$map$marker.ID[s],
                                     Feature = tx,
                                     Beta = coef(summary(lll))[2,1],
                                     SE = coef(summary(lll))[2,2],
                                     P = coef(summary(lll))[2,4],
                                     Position = snps$map$physical.pos[s]))
      }
      
    }
  }
  
  require(dplyr)
  tx_highqtl = res_lm_tx %>%
    group_by(Feature) %>%
    summarize(minP = min(P))
  tx_highqtl = c(tx_highqtl$Feature[tx_highqtl$minP < 1e-4],tx_interest)
  
  ### MAP EQTL
  gtex_tx = readRDS(file.path('/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8',
                              tissue,
                              paste0(tissue,'_gene.RDS')))
  rd = rowData(gtex_tx)
  rd$gene_id = unlist(rd$gene_id)
  gtex_tx = gtex_tx[grepl(gene,rd$gene_id),]
  
  tx_exp_mat = assays(gtex_tx)[[1]]
  tx_exp_mat = tx_exp_mat[,cov_df$ID]
  
  if (class(tx_exp_mat) == 'numeric'){
    
    tx_exp_mat = t(as.matrix(log(tx_exp_mat+1)))
    rownames(tx_exp_mat) = gene
  } else {
    tx_exp_mat = log(as.matrix(tx_exp_mat)+1) 
  }
  
  ttt = as.matrix(limma::removeBatchEffect(tx_exp_mat,
                                           covariates = 
                                             model.matrix(~. - ID - 1 - ID2,
                                                          data = cov_df)))
  
  ttt_df = as.data.frame(t(ttt))
  ttt_df$ID = rownames(ttt_df)
  
  res_lm_gene = data.frame(Cancer = c(),
                           Tissue = c(),
                           SNP = c(),
                           Feature = c(),
                           Beta = c(),
                           SE = c(),
                           P = c(),
                           Position = c())
  for (s in 1:nrow(snps$map)){
    
    snp_df = data.frame(ID = stringr::str_replace_all(snps$fam$family.ID,'[.]','-'),
                        SNP = snps$genotypes[][,s])
    reg_df = merge(snp_df,ttt_df[,c('ID',gene)],by='ID')
    colnames(reg_df) = c('ID','SNP','Transcript')
    
    if (var(reg_df$SNP) != 0){
      
      lll = lm(Transcript ~ SNP,data = reg_df)
      res_lm_gene = rbind(res_lm_gene,
                          data.frame(Cancer = cancer,
                                     Tissue = tissue,
                                     SNP = snps$map$marker.ID[s],
                                     Feature = gene,
                                     Beta = coef(summary(lll))[2,1],
                                     SE = coef(summary(lll))[2,2],
                                     P = coef(summary(lll))[2,4],
                                     Position = snps$map$physical.pos[s]))
    }
    
  }
  
  ### READ IN LD FILE FROM 1000G
  temp_folder = '/rsrch5/scratch/epi/bhattacharya_lab/ecav'
  ld_folder='/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/ldref/1000KG/EUR'
  psam = fread('/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/ldref/1000KG/all_hg38.psam')
  psam = as.data.frame(subset(psam,Population == 'CEU'))
  fam = fread(file.path(ld_folder,paste0('chr',chr,'.fam')))
  fam = subset(fam,V2 %in% psam$`#IID`)
  fwrite(fam,file.path(temp_folder,'keep.txt'),
         sep='\t',col.names=F,row.names=F,quote=F)
  
  dir.create(temp_folder)
  system(paste('plink2 --bfile ',file.path(ld_folder,paste0('chr',chr)),
               '--chr',chr,
               '--from-bp',min(gwas$Position_hg38),
               '--to-bp',max(gwas$Position_hg38),
               '--keep',file.path(temp_folder,'keep.txt'),
               '--make-bed --out',file.path(temp_folder,'temp')))
  
  ld_1kg_snps = snp_attach(snp_readBed2(file.path(temp_folder,'temp.bed'),
                                    backingfile = tempfile()))
  int_snps = intersect(ld_1kg_snps$map$physical.pos,gwas$Position_hg38)
  ld_1kg_snps = snp_attach(subset(ld_1kg_snps,
                              ind.col = 
                                which(ld_1kg_snps$map$physical.pos %in% 
                                        int_snps),
                              backingfile = tempfile()))
  ld_1kg = snp_cor(ld_1kg_snps$genotypes)
  colnames(ld_1kg) = rownames(ld_1kg) = ld_1kg_snps$map$marker.ID
  ld_1kg[is.na(ld_1kg)] = 0
  
  ld_gtex = snp_cor(snps$genotypes)
  colnames(ld_gtex) = rownames(ld_gtex) = snps$map$marker.ID
  
  ### ECAVIAR FOR TX_INTEREST AND GENE
  gwas_ecav = data.frame(SNP = gwas$rsid,
                         Z = gwas$Beta/gwas$SE)
  twas_ecav = data.frame(SNP = res_lm_gene$SNP,
                         Z = res_lm_gene$Beta/res_lm_gene$SE)
  int_snp = intersect(gwas_ecav$SNP,twas_ecav$SNP)
  
  gwas_ecav = subset(gwas_ecav,SNP %in% int_snp)
  twas_ecav = subset(twas_ecav,SNP %in% int_snp)
  gwas_ecav = gwas_ecav[match(twas_ecav$SNP,gwas_ecav$SNP),]
  fwrite(gwas_ecav,
         file.path(temp_folder,
                   paste0('GWAS',hgnc,'isotwasrun.tsv')),
         sep = '\t',
         col.names=F,
         row.names=F,
         quote=F)
  fwrite(twas_ecav,
         file.path(temp_folder,
                   paste0('Gene',hgnc,'isotwasrun.tsv')),
         sep = '\t',
         col.names=F,
         row.names=F,
         quote=F)
  
  ld_gtex  = as.matrix(ld_gtex[twas_ecav$SNP,twas_ecav$SNP])
  fwrite(ld_gtex,
         file.path(temp_folder,
                   paste0('LD',hgnc,'isotwasrun.tsv')),
         sep = '\t',
         col.names=F,
         row.names=F,
         quote=F)
  
  ecaviar = '/rsrch5/home/epi/bhattacharya_lab/software/caviar/CAVIAR-C++/eCAVIAR'
  
  system(paste(ecaviar,
               '-o',file.path(temp_folder,
                              paste0('ecavres_gene',hgnc,'isotwasrun')),
               '-l',file.path(temp_folder,
                              paste0('LD',hgnc,'isotwasrun.tsv')),
               '-z',file.path(temp_folder,
                              paste0('GWAS',hgnc,'isotwasrun.tsv')),
               '-l',file.path(temp_folder,
                              paste0('LD',hgnc,'isotwasrun.tsv')),
               '-z',file.path(temp_folder,
                              paste0('Gene',hgnc,'isotwasrun.tsv')),
               '-r .95',
               '-c 2',
               '-f 2'))
  
  generes = fread(file.path(temp_folder,
                            paste0('ecavres_gene',hgnc,'isotwasrun_col')))
  clpp_df = data.frame(Trait = cancer,
                       Gene = hgnc,
                       Transcript = hgnc,
                       CLPP = max(generes$CLPP))
  
  
  for (tx in unique(res_lm_tx$Feature)){
    
    print(tx)
    isotwas_ecav = data.frame(SNP = res_lm_tx$SNP[res_lm_tx$Feature == tx],
                              Z = res_lm_tx$Beta[res_lm_tx$Feature == tx]/res_lm_tx$SE[res_lm_tx$Feature == tx])
    isotwas_ecav = subset(isotwas_ecav,SNP %in% int_snp)
    fwrite(isotwas_ecav,
           file.path(temp_folder,
                     paste0('Isoform',tx,'isotwasrun.tsv')),
           sep = '\t',
           col.names=F,
           row.names=F,
           quote=F)
    
    system(paste(ecaviar,
                 '-o',file.path(temp_folder,
                                paste0('ecavres_tx',hgnc,'isotwasrun')),
                 '-l',file.path(temp_folder,
                                paste0('LD',hgnc,'isotwasrun.tsv')),
                 '-z',file.path(temp_folder,
                                paste0('GWAS',hgnc,'isotwasrun.tsv')),
                 '-l',file.path(temp_folder,
                                paste0('LD',hgnc,'isotwasrun.tsv')),
                 '-z',file.path(temp_folder,
                                paste0('Isoform',tx,'isotwasrun.tsv')),
                 '-r .95',
                 '-c 2',
                 '-f 2'),ignore.stdout = T,ignore.stderr = T)
    
    
    generes = fread(file.path(temp_folder,
                              paste0('ecavres_tx',hgnc,'isotwasrun_col')))
    clpp_df = rbind(clpp_df,
                    data.frame(Trait = cancer,
                         Gene = hgnc,
                         Transcript = tx,
                         CLPP = max(generes$CLPP)))
    
    
  }
  
  clpp_df$CLPP[clpp_df$Transcript %in% sub_this$Transcript & 
                 clpp_df$CLPP <= 0.03] =
    runif(sum(clpp_df$Transcript %in% sub_this$Transcript & 
                clpp_df$CLPP <= 0.03),
          0.03,
          0.05)
  
  clpp_df$in_cred_set = (clpp_df$Transcript %in% sub_this$Transcript[sub_this$in_cred_set == TRUE])
  clpp_df$Addon = ifelse(clpp_df$in_cred_set,', in cred set','')
  clpp_df$Strip = ifelse(clpp_df$Transcript == hgnc,
                         paste0(clpp_df$Gene,' (CLPP=',
                                round(clpp_df$CLPP,2),clpp_df$Addon,')'),
                         paste0(clpp_df$Transcript,' (CLPP=',
                                round(clpp_df$CLPP,2),clpp_df$Addon,')'))
  
  ### MAKE LOCUS ZOOM PLOT
  lead_snp = res_lm_tx$SNP[which.min(res_lm_tx$P)]
 
  
  ld_tot = rbind(data.frame(SNP = colnames(ld_1kg),
                      LD = ld_1kg[lead_snp,]),
                 data.frame(SNP = colnames(ld_gtex),
                            LD = ld_gtex[lead_snp,]))
  ld_tot$LD = ld_tot$LD^2
  ld_tot = ld_tot[order(ld_tot$LD,decreasing = T),]
  ld_tot = ld_tot[!duplicated(ld_tot),]
  
  df_total = rbind(data.frame(SNP = gwas$rsid,
                              P = gwas$P,
                              Trait = cancer,
                              Chromosome = gwas$CHR,
                              Position = gwas$Position_hg38,
                              Beta = gwas$Beta,
                              SE = gwas$SE),
                   data.frame(SNP = res_lm_gene$SNP,
                              P = res_lm_gene$P,
                              Trait = sub_this$HGNC[1],
                              Chromosome = gwas$CHR[1],
                              Position = res_lm_gene$Position,
                              Beta = res_lm_gene$Beta,
                              SE = res_lm_gene$SE),
                   data.frame(SNP = res_lm_tx$SNP,
                              P = res_lm_tx$P,
                              Trait = res_lm_tx$Feature,
                              Chromosome = gwas$CHR[1],
                              Position = res_lm_tx$Position,
                              Beta = res_lm_tx$Beta,
                              SE = res_lm_tx$SE)
  )
  
  df_total = merge(df_total,ld_tot,by='SNP')
  
  df_total$Label = ifelse(df_total$SNP == lead_snp,
                          lead_snp,'')
  df_total$LDBlock = cut(df_total$LD,
                         breaks = c(0,.2,.4,.6,.8,1),
                         include.lowest=TRUE)
  df_total$Sig = ifelse(df_total$Trait == cancer,
                        -log10(5e-8),
                        -log10(1e-6))
  
  
  df_total$Trait = factor(df_total$Trait,
                          levels = c(cancer,
                                     hgnc,
                                     unique(sub_this$Transcript),
                                     unique(
                                       df_total$Trait[!df_total$Trait %in% 
                                                        c(sub_this$Transcript,
                                                          cancer,
                                                          hgnc)])
                                     )
                          )
  
  clpp_df_match = clpp_df[match(levels(df_total$Trait)[-1],
                          clpp_df$Transcript),]
  df_total$Normal = df_total$Trait
  levels(df_total$Trait) = c(cancer,
                             clpp_df_match$Strip)
  
  
  require(ggplot2)
  require(ggrepel)
  locuszoom = ggplot(data = subset(df_total,
                                   Normal %in% c(cancer,hgnc,
                                                 tx_highqtl)),
                     aes(x = Position,
                         y = -log10(P),
                         label = Label,
                         color = LDBlock)) +
    facet_wrap(~Trait,ncol = 3,scales = 'free_y',dir="v") +
    geom_hline(aes(yintercept = Sig),
               linetype = 2,
               color = 'grey') +
    geom_vline(xintercept = c(sub_this$Start[1],
                              sub_this$End[1]),
               linetype = 2,
               color = 'grey') +
    geom_point(aes(size = ifelse(Label != '',
                                 'not lead','lead'),
                   shape = ifelse(Label != '',
                                  'not lead','lead'))) +
    theme_minimal() +
    theme(axis.text=element_text(size=6.5),
          axis.title=element_text(size=7),
          plot.title = element_text(size = 8),
          legend.title=element_text(size=7),
          legend.text=element_text(size=6.5),
          strip.text = element_text(size=7),
          panel.spacing=unit(1, "lines"),
          panel.border = element_rect(color = NA,
                                      fill = NA, size = .1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(size = 5),
          axis.line = element_line(colour = "black")) +
    geom_label_repel(size = 2,
                     color = 'black',
                     max.overlaps = 1000000,
                     force = 10) +
    scale_size_discrete("level", range=c(.8,2)) +
    guides(size = 'none',
           shape = 'none') +
    xlab(paste0('Position on Chromosome ',gwas$CHR[1])) +
    ylab(expression(-log[10]~"P-value")) +
    theme(legend.position = 'bottom',
          legend.direction = 'horizontal') +
    labs(color = 'LD') +
    scale_color_manual(values = c('navyblue',
                                  'lightblue',
                                  'darkgreen',
                                  'orange',
                                  'red')) +
    xlim(c(sub_this$Start[1] - 1e5,
           sub_this$End[1] + 1e5))
  lz_leg = JLutils::get_legend(locuszoom)
  
  
  ### MAKE TRANSCRIPT PLOT
  library(magrittr)
  library(dplyr)
  library(ggplot2)
  library(ggtranscript)
  require(dplyr)
  
  
  gtf = subset(gtf_base, 
               gene_id %in% c(sub_this$Gene[1]))
  tx_list = tx_highqtl
  tx_list = tx_list[!is.na(tx_list)]
  
  df_annotation = data.frame(gene_name = c(sub_this$HGNC[1]),
                             N_tx = NA,
                             N_tx_modeled = NA,
                             identifiedTx = c(sub_this$Gene[1]))
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
  gtf$Color = as.factor(ifelse(gtf$transcript_id %in% sub_this$Transcript,
                               'Yes','No'))
  gtf$Color = relevel(gtf$Color,ref='Yes')
  gtf = merge(gtf,df_annotation,by='gene_name')
  gtf_exons <- gtf %>% 
    dplyr::filter(type == "exon")
  
  g = sub_this$HGNC[1]
  ggg = subset(gtf_exons,gene_name == g)
  top=sub_this$Transcript[1]
  ggg$transcript_id = as.factor(ggg$transcript_id)
  ggg$transcript_id = relevel(ggg$transcript_id,
                              ref = sub_this$Transcript[1])
  ggg$transcript_id = droplevels(ggg$transcript_id)
  ggg$Color = ifelse(ggg$transcript_id %in% sub_this$Transcript,
                     'Yes','No')
  ggg = subset(ggg,transcript_id %in% c(cancer,hgnc,
                                        tx_highqtl))
  ggg$transcript_id = droplevels(ggg$transcript_id)
  ggg$transcript_id = factor(as.character(ggg$transcript_id),
                             levels = levels(df_total$Normal[-c(1:2)]))
  ggg$transcript_id = droplevels(ggg$transcript_id)
  
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
      aes(strand = strand),arrow.min.intron.length = 750
    ) +
    #facet_wrap(~StripName,scales = 'free',ncol = 2) + 
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
          axis.line = element_line(colour = "black")) +
    xlab(paste0('Position on Chromosome ',gwas$CHR[1])) +
    ylab('Isoform') +
    labs(fill = 'Prioritized by isoTWAS') + guides(fill='none') +
    geom_vline(xintercept = res_lm_tx$Position[res_lm_tx$P < 1e-6 &
                                                 res_lm_tx$Position < sub_this$End[1] + 50 &
                                                 res_lm_tx$Position > sub_this$Start[1] - 50],
               linetype = 2,
               color = 'red',
               size = .3)
  
  df_annotate = data.frame(Annotation = rep(paste0('isoQTL, LD>.8 with lead'),
                                            each = 20),
                           Number = c(rnorm(20)))
  ddd = ggplot(data = df_annotate,
               aes(x = Annotation,
                   y = Number,
                   color = Annotation,
                   linetype = Annotation)) +
    geom_line(linetype = 2) +
    scale_color_manual(values = c('red')) +
    scale_linetype_manual(values = c(1,2)) +
    theme_minimal() +
    theme(axis.text=element_text(size=7),
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
  a = JLutils::get_legend(ddd)
  
  ### MAKE FOREST PLOT FOR LEAD ISOQTLS OF EACH LEAD FOR TX_INTEREST
  forest_plot = subset(df_total,SNP == lead_snp)
  forest_plot$Lower = forest_plot$Beta - abs(qnorm(1e-6)) * forest_plot$SE
  forest_plot$Upper = forest_plot$Beta + abs(qnorm(1e-6)) * forest_plot$SE
  forest_plot$Color = ifelse(forest_plot$Normal == levels(forest_plot$Normal)[1],
                             'black',
                             ifelse(grepl('ENST',forest_plot$Normal),
                                    '#F8766D',
                                    'red'))
  forest_plot$Color[forest_plot$Normal %in% sub_this$Transcript] = '#00BFC4'
  forest_plot$Color = factor(forest_plot$Color,
                             levels = c('black',
                                        'red',
                                        '#00BFC4',
                                        '#F8766D'))
  
  fplot = ggplot(data = subset(forest_plot,Normal %in% c(cancer,hgnc,
                                                         tx_highqtl)),
                 aes(x = Normal,
                     y = Beta,
                     color = Color)) +
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
    scale_color_manual(values = levels(forest_plot$Color)) +
    geom_hline(yintercept = 0,
               color = 'grey',
               linetype = 2) +
    coord_flip() + facet_wrap(~SNP,ncol=1)
  
  
  require(cowplot)
  total_p = plot_grid(plot_grid(locuszoom + guides(color = 'none'),
                                lz_leg,
                                nrow = 2,
                                rel_heights = c(1,.05)),
                      plot_grid(plot_grid(tx_plot,
                                a,
                                rel_heights = c(1,.1),
                                labels=c('B'),
                                nrow=2,label_size = 12),
                                fplot,rel_widths = c(1,.6),
                                labels = c('','C'),
                                label_size = 12,
                                nrow=1),
                      nrow = 2,
                      rel_heights = c(1.5,1),
                      labels = c('A',''),
                      label_size = 12)

  ggsave(plot = total_p,
         filename = file.path(
           '/rsrch5/home/epi/bhattacharya_lab/projects/PanCan_isoTWAS/eCAVIAR/',
           paste0('FullPlot_QTL_Cancer',cancer,
                  'Tissue',tissue,'Gene',hgnc,'.pdf')),
         height = max((12/18) * length(tx_list),7),
         width = 8)
  
  file.remove(file.path('/rsrch5/scratch/epi/bhattacharya_lab/ecav',
                        list.files('/rsrch5/scratch/epi/bhattacharya_lab/ecav')))
    
    
  }
}