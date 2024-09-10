require(data.table)
require(SummarizedExperiment)
setwd('/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8')
tissues = list.files()
tissues = tissues[!grepl('.bed',tissues)]
tissues = tissues[!grepl('.fam',tissues)]
tissues = tissues[!grepl('.bim',tissues)]
tissues = tissues[!grepl('pca',tissues)]

for (t in tissues){
  
  setwd(t)
  
  fff = list.files()
  
  cov_file = fff[grepl('covar',fff)]
  cov_df = fread(cov_file)
  colnames(cov_df)[1:2] = c('ID','ID2')
  cov_df$ID = stringr::str_replace_all(cov_df$ID,'[.]','-')
  
  tx = readRDS(fff[grepl('transcripts.RDS',fff)])
  tx = tx[,cov_df$ID]
  rd = as.data.frame(rowData(tx))
  tx_df = data.frame(GENE = as.character(rowData(tx)$tx_name),
                     CHR = sapply(strsplit(as.character(rowData(tx)$seqnames),
                                    'r'),function(x) x[2]),
                     GENE_COORD = as.numeric(rowData(tx)$start))
  tx_mat = log(as.matrix(assays(tx)[[1]])+1)
  tx_mat = t(as.matrix(limma::removeBatchEffect(tx_mat,
                                       covariates = 
                                         model.matrix(~. - ID - 1 - ID2,
                                                      data = cov_df))))
  tx_df = cbind(tx_df,t(tx_mat))
  colnames(tx_df) = c('GENE','CHR','GENE_COORD',stringr::str_replace_all(colnames(tx),'-','.'))
  fwrite(tx_df,
         paste0('MESC_tx_exp.tsv'),
         sep='\t',
         col.names=T,
         row.names=F,
         quote=F)
  rm(tx,tx_df,tx_mat)
  
  tx = readRDS(fff[grepl('gene.RDS',fff)])
  tx = tx[,cov_df$ID]
  rd$GENE = unlist(rd$gene_id)
  rd$CHR = sapply(strsplit(as.character(rd$seqnames),
                           'r'),function(x) x[2])
  rd$GENE_COORD = as.numeric(rd$start)
  tx_df = data.frame(GENE = as.character(rowData(tx)$gene_id))
  tx_df = merge(tx_df,rd[,c('GENE','CHR','GENE_COORD')],by='GENE')
  tx_df = tx_df[!duplicated(tx_df$GENE),]
  tx = tx[tx_df$GENE,]
  tx_mat = log(as.matrix(assays(tx)[[1]])+1)
  tx_mat = t(as.matrix(limma::removeBatchEffect(tx_mat,
                                                covariates = model.matrix(~. - ID - 1 - ID2,data = cov_df))))
  tx_df = cbind(tx_df,t(tx_mat))
  colnames(tx_df) = c('GENE','CHR','GENE_COORD',stringr::str_replace_all(colnames(tx),'-','.'))
  fwrite(tx_df,
         paste0('MESC_gene_exp.tsv'),
         sep='\t',
         col.names=T,
         row.names=F,
         quote=F)
  rm(tx_df,tx_mat)
  
  fam_keep = data.frame(V1=stringr::str_replace_all(colnames(tx),'-','.'),
                        V2=stringr::str_replace_all(colnames(tx),'-','.'),
                        V3=-9)
  rm(tx)
  fwrite(fam_keep,
         paste0('keep.fam'),
         sep='\t',
         col.names=T,
         row.names=F,
         quote=F)
  
  setwd('/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8')
  
}




for (t in tissues){
  
  setwd('/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8')
  setwd(t)
  print(t)
  tx = fread('MESC_tx_exp.tsv')
  
  imputeMean = function(x){
    x[is.na(x)] = mean(x,na.rm=T)
  }
  
  rowvars = apply(tx[,4:ncol(tx)],1,var)
  print(sum(rowvars == 0))
  tx = tx[rowvars > 0,]
  if (any(rowMeans(is.na(tx[,4:ncol(tx)]))>0)){
    print(t)
    tx[,4:ncol(tx)] = 
      apply(tx[,4:ncol(tx)],1,imputeMean)
    
    fwrite(tx,
           paste0('MESC_gene_exp.tsv'),
           sep='\t',
           col.names=T,
           row.names=F,
           quote=F)
  }
  
  tx = fread('MESC_gene_exp.tsv')
  
  rowvars = apply(tx[,4:ncol(tx)],1,var)
  print(sum(rowvars == 0))
  tx = tx[rowvars > 0,]
  if (any(rowMeans(is.na(tx[,4:ncol(tx)]))>0)){
    print(t)
    tx[,4:ncol(tx)] = 
      apply(tx[,4:ncol(tx)],1,imputeMean)
    print(sum(rowvars == 0))
    tx = tx[rowvars > 0,]
    fwrite(tx,
           paste0('MESC_gene_exp.tsv'),
           sep='\t',
           col.names=T,
           row.names=F,
           quote=F)
  }
  
  setwd('/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8')
  
}