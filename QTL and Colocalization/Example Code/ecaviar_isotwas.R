library("optparse")
option_list <- list(
  make_option(c("-i", "--index"), action="store_true", default=TRUE,
              help="index",type='numeric')
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
index = opt$index


require(data.table)
require(ggplot2)
require(vroom)


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
                         Cancer = c('Prostate (advanced)',
                                    'Breast (Overall)',
                                    'Breast (ER-)',
                                    'Ovarian (serou)',
                                    'Colorectal',
                                    'Endometrial',
                                    'Lung',
                                    'Lung (ADE)',
                                    'Lung (SQC)',
                                    'Breast (ER+)',
                                    'Ovarian (overall)',
                                    'Prostate (overall)'))
manifest_df$File = file.path('/rsrch5/home/epi/bhattacharya_lab/projects/PanCan_isoTWAS/GWAS',
                             manifest_df$File)

setwd('/rsrch5/home/epi/bhattacharya_lab/projects/PanCan_isoTWAS')
inGene_base = fread('isoTWAS_FineMap_PanCan.tsv')
inGene_base = inGene_base[!duplicated(inGene_base),]


chunk = ceiling(nrow(inGene_base)/300)
start = ((index-1) * chunk) + 1
end = min(c(((index) * chunk),nrow(inGene_base)))

for (index in start:end){

inGene = inGene_base[index,]

gene = inGene$Gene[1]
tx = inGene$Transcript[1]
tissue = inGene$Tissue[1]
cancer = inGene$Indication[1]
chr = inGene$Chromosome[1]
start = inGene$Start[1]
end = inGene$End[1]
gwas_file = manifest_df$File[manifest_df$Cancer == cancer]

print(index)
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

cov_df = read.table(file.path('/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8',
                              tissue,
                              paste0(tissue,'.v8.covariates.txt.covar')))

colnames(cov_df)[1:2] = c('ID','ID2')
cov_df$ID = stringr::str_replace_all(cov_df$ID,'[.]','-')

require(SummarizedExperiment)

### MAP ISOQTL
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
                       P = c())
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
                                 P = coef(summary(lll))[2,4]))
  }
  
}

tx_vec = reg_df$Transcript


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

res_lm_gene = data.frame(Cancer = c(),
                         Tissue = c(),
                         SNP = c(),
                         Feature = c(),
                         Beta = c(),
                         SE = c(),
                         P = c())
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
                                   P = coef(summary(lll))[2,4]))
  }
  
}

ld = snp_cor(snps$genotypes)
rownames(ld) = colnames(ld) = snps$map$marker.ID
ld_df = data.frame(SNP = rownames(ld),
                   LD = ld[res_lm_tx$SNP[which.min(res_lm_tx$P)],])

interest = 'No'
if (min(res_lm_tx$P) < 1e-6){
  interest = 'Yes'
}

df_snp = data.frame(ID = stringr::str_replace_all(snps$fam$family.ID,'[.]','-'),
                    SNP = snps$genotypes[][,
                                           which(snps$map$marker.ID == 
                                                   res_lm_tx$SNP[which.min(res_lm_tx$P)])])
df_snp = df_snp[match(reg_df$ID,df_snp$ID),]

### Prep eCAVIAR files
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

ld_gt  = as.matrix(ld[twas_ecav$SNP,twas_ecav$SNP])
fwrite(ld,
       file.path(temp_folder,
                 paste0('LD',hgnc,'isotwasrun.tsv')),
       sep = '\t',
       col.names=F,
       row.names=F,
       quote=F)

ecaviar = '/rsrch5/home/epi/bhattacharya_lab/software/caviar/CAVIAR-C++/eCAVIAR'


for (tx in unique(res_lm_tx$Feature)){
  
  isotwas_ecav = data.frame(SNP = res_lm_tx[res_lm_tx$Feature == tx]$SNP,
                            Z = res_lm_tx[res_lm_tx$Feature == tx]$Beta/res_lm_tx[res_lm_tx$Feature == tx]$SE)
  isotwas_ecav = subset(isotwas_ecav,SNP %in% int_snp)
  fwrite(isotwas_ecav,
         file.path(temp_folder,
                   paste0('Isoform',tx,'isotwasrun.tsv')),
         sep = '\t',
         col.names=F,
         row.names=F,
         quote=F)
  
}





system(paste(ecaviar,
             '-o',file.path(temp_folder,
                            paste0('ecavres_gene',index,'isotwasrun')),
             '-l',file.path(temp_folder,
                            paste0('LD',index,'isotwasrun.tsv')),
             '-z',file.path(temp_folder,
                            paste0('GWAS',index,'isotwasrun.tsv')),
             '-l',file.path(temp_folder,
                            paste0('LD',index,'isotwasrun.tsv')),
             '-z',file.path(temp_folder,
                            paste0('Gene',index,'isotwasrun.tsv')),
             '-r .95',
             '-c 2',
             '-f 2'))

clpp_res_gene = fread(file.path(temp_folder,
                           paste0('ecavres_gene',index,'isotwasrun_col')))


system(paste(ecaviar,
             '-o',file.path(temp_folder,
                            paste0('ecavres_isoform',index,'isotwasrun')),
             '-l',file.path(temp_folder,
                            paste0('LD',index,'isotwasrun.tsv')),
             '-z',file.path(temp_folder,
                            paste0('GWAS',index,'isotwasrun.tsv')),
             '-l',file.path(temp_folder,
                            paste0('LD',index,'isotwasrun.tsv')),
             '-z',file.path(temp_folder,
                            paste0('Isoform',index,'isotwasrun.tsv')),
             '-r .95',
             '-c 2',
             '-f 2'))

clpp_res_isoform = fread(file.path(temp_folder,
                                paste0('ecavres_isoform',index,'isotwasrun_col')))



saveRDS(list(GWAS = gwas,
             eQTL = res_lm_gene,
             isoQTL = res_lm_tx,
             Gene = gene,
             Transcript = tx,
             Cancer = cancer,
             Tissue = tissue,
             LD = ld_df,
             SNPVec = df_snp$SNP,
             TxVec = tx_vec,
             GeneVec = reg_df$Transcript,
             Coloc_Gene = clpp_res_gene,
             Coloc_Tx = clpp_res_isoform),
        file.path('/rsrch5/home/epi/bhattacharya_lab/projects/PanCan_isoTWAS/eCAVIAR/',
                  paste0('QTL_Tissue',tissue,'_',
                         'Cancer',cancer,'_',
                         'Transcript',tx,'_',
                         'Gene',gene,'.RDS')))
}