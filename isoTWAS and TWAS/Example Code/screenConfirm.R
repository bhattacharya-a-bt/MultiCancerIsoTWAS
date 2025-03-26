setwd('/Users/abhattacharya3/OneDrive - Inside MD Anderson/Fall 2024/PanCan IsoTWAS/Results/')

require(data.table)
require(isotwas)

require(biomaRt)
ensembl <- useEnsembl(biomart = "ensembl", 
                      dataset = "hsapiens_gene_ensembl")

manifest = as.data.frame(
  fread('twas_cancer_tissue_full.tsv',header=T)
)

manifest$Cancer = sapply(
  strsplit(manifest$Cancer,'_MAF'),
  function(x) x[1]
)  
gnomad = 
    fread('gnomad_plof.txt')
gnomad$HGNC = gnomad$gene

for (cancer in unique(manifest$Cancer)){
 
  setwd('/Users/abhattacharya3/OneDrive - Inside MD Anderson/Fall 2024/PanCan IsoTWAS/Results/')
  setwd('isoTWAS/')
  setwd(cancer)
  print(cancer)
  
  
  ddd = list.dirs()
  manifest_this = subset(manifest,Cancer == cancer)
  ddd = ddd[!ddd %in% paste0('./',manifest_this$Tissue)]
  unlink(ddd,recursive = T)
  
  ddd = list.dirs()
  ddd = ddd[ddd != '.']
  for (d in ddd){
    setwd(d)
    print(d)
    
    isotwas_res = vroom::vroom(list.files()[1],show_col_types = F)
    bm = getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol',
                              'chromosome_name','start_position','end_position',
                              'gene_biotype'),
               filters = 'ensembl_gene_id',
               values = unique(isotwas_res$Gene), 
               mart = ensembl)
    colnames(bm) = c('Gene','HGNC','Chromosome','Start','End','Biotype')
    
    isotwas_res = isotwas_res[complete.cases(isotwas_res),]
    isotwas_res = isotwas_res[order(abs(isotwas_res$Z),decreasing = T),]
    isotwas_res = isotwas_res[!duplicated(isotwas_res$Transcript) & 
                                abs(isotwas_res$Z) < Inf,]
    colnames(isotwas_res)[8] = 'R2'
    isotwas_res = merge(bm,isotwas_res,by='Gene')
    
    www = which(isotwas_res$Chromosome == 6 &
                  isotwas_res$Start < 35e6 &
                  isotwas_res$End > 27e6)
    isotwas_res = isotwas_res[-www,]
    
    gene = data.frame(Gene = unique(isotwas_res$Gene),
                      HGNC = NA,
                      Chromosome = NA,
                      Start = NA,
                      End = NA,
                      Biotype = NA,
                      Screen.P = NA)
    
    require(tidyverse)
    gene = isotwas_res %>%
      group_by(Gene) %>%
      summarise(HGNC = unique(HGNC),
                Chromosome = unique(Chromosome),
                Start = unique(Start),
                End = unique(End),
                Biotype = unique(Biotype),
                Screen.P = isotwas::p_screen(P))
    
    alpha1=.05
    G = nrow(gene)
    gene$Screen.P.Adjusted = p.adjust(gene$Screen.P,method = 'fdr')
    R = length(unique(gene$Gene[gene$Screen.P.Adjusted < alpha1]))
    alpha2 = (R*alpha1)/G
    isoform_new = as.data.frame(matrix(nrow = 0,
                                       ncol = ncol(isotwas_res)+2))
    colnames(isoform_new) = c(colnames(isotwas_res),'Screen.P','Confirmation.P')
    gene = gene[order(gene$Screen.P),]
    ttt = merge(isotwas_res,
                gene[,c('Gene','Screen.P',
                        'Screen.P.Adjusted')])
    isoform_new = ttt %>%
      group_by(Gene) %>%
      summarise(Transcript = Transcript,
                Confirmation.P = isotwas::p_confirm(P,alpha = alpha2))
    isoform_new = merge(isoform_new,ttt,by=c('Gene','Transcript'))
    isoform_new$Confirmation.P = ifelse(isoform_new$Screen.P.Adjusted < 0.05,
                                        isoform_new$Confirmation.P,
                                        1)
    isoform_new$Tissue = strsplit(d,'./')[[1]][2]
    isoform_new = isoform_new[!duplicated(isoform_new),]
    isoform_sig = subset(isoform_new, 
                         Screen.P < alpha1 &
                           Confirmation.P < alpha2 &
                           permute.P < 0.05)
    setwd('/Users/abhattacharya3/OneDrive - Inside MD Anderson/Fall 2024/PanCan IsoTWAS/Results/')
    setwd('isoTWAS/')
    setwd(cancer)
    
    fwrite(isoform_new,
           'ScreenConfirm_isoTWAS.tsv',
           append=T,
           row.names=F,
           quote=F)
    
    if (nrow(isoform_sig) > 1){
      isoform_sig$Chromosome = as.numeric(isoform_sig$Chromosome)
      
      isoform_sig = isoform_sig[order(isoform_sig$Chromosome,
                                      isoform_sig$Start),]
      isoform_sig$Indication = manifest$Indication[manifest$Cancer == cancer][1]
      isoform_sig = isoform_sig[,c('Indication','Gene','HGNC','Tissue','Chromosome',
                                   'Start','End','Biotype','Transcript',
                                   'Z','P','permute.P','topSNP','topSNP.P',
                                   'Screen.P.Adjusted','Confirmation.P','R2')]
      colnames(isoform_sig) = c('Indication','Gene','HGNC','Tissue','Chromosome',
                                'Start','End','Biotype','Transcript',
                                'Z','P','Permutation P','Top GWAS SNP','Top GWAS P',
                                'Screening Adjusted P','Confirmation P','R2')
      
      isoform_sig = merge(isoform_sig,
                          gnomad[,c('HGNC','pLI')],
                          by = 'HGNC',
                          all.x=T)
      
      
      fwrite(isoform_sig,
             'SignificantAssociations_isoTWAS_noFineMap.tsv',
             append=T,
             row.names=F,
             quote=F)
    }
    
    
  }
  
  
  setwd('/Users/abhattacharya3/OneDrive - Inside MD Anderson/Fall 2024/PanCan IsoTWAS/Results/')
}
