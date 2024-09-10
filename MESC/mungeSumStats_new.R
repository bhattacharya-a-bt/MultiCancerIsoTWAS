library("optparse")
option_list <- list(
  make_option(c("-i", "--index"), action="store_true", default=TRUE,
              help="index",type='numeric')
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
index = opt$index

library(MungeSumstats)

setwd('/rsrch5/home/epi/bhattacharya_lab/data/munged_GWAS')

fff = c('GCST90244168_buildGRCh38.tsv.gz',
        'GCST90244167_buildGRCh38.tsv.gz',
        'GCST90129505_buildGRCh37.tsv',
        'GCST90134661_buildGRCh37.tsv')
df = data.frame(files = fff,
                cancer = c('OvarianSerous_Dareng',
                           'Ovarian_Dareng',
                           'Colorectal_Fernandez',
                           'Lung_Byun'))

i = index
    
  file_this = df$files[i]
  cancer = df$cancer[i]
  
  effective_n <- function(ncase, ncontrol)
  {
    return(2 / (1/ncase + 1/ncontrol))
  }
  
  dir.create('NewGWAS')
  a = data.table::fread(file_this)
  
  if (i == 1){
    
    
    a$N_Cases = 15588
    a$N_Controls = 105724
    
    
    a$CHR = a$chromosome
    a$Position = a$base_pair_location
    a$Beta = a$beta
    a$SE = a$standard_error
    a$N = ceiling(effective_n(a$N_Cases,a$N_Controls))
    a$P = 2*pnorm(-abs(a$Beta/a$SE))
    
    a = a[,c('SNP','CHR','Position','Beta','SE','P','N','N_Cases','N_Controls')]
    
  }
  
  if (i == 2){
    
    
    a$N_Cases = 23394
    a$N_Controls = 105724
    
    a$CHR = a$chromosome
    a$Position = a$base_pair_location
    a$Beta = a$beta
    a$SE = a$standard_error
    a$N = ceiling(effective_n(a$N_Cases,a$N_Controls))
    a$P = 2*pnorm(-abs(a$Beta/a$SE))
    
    a = a[,c('SNP','CHR','Position','Beta','SE','P','N','N_Cases','N_Controls')]
  }
  
  if (i == 3){
    
    a$SNP = a$rs_id
    a$CHR = a$chromosome
    a$Position = a$base_pair_location
    a$Beta = a$beta
    a$SE = a$standard_error
    a$N_Cases = 100204
    a$N_Controls = 154587
    a$N = ceiling(effective_n(a$N_Cases,a$N_Controls))
    a$P = a$p_value
    
    a = a[,c('SNP','CHR','Position','Beta','SE','P','N','N_Cases','N_Controls')]
    
  }
  
  if (i == 4){
    
    a$SNP = a$variant_id
    a$CHR = a$chromosome
    a$Position = a$base_pair
    a$Beta = log(a$odds_ratio)
    a$SE = a$standard_error
    a$N_Cases = 61047
    a$N_Controls = 947237
    a$N = ceiling(effective_n(a$N_Cases,a$N_Controls))
    a$P = a$p_value
    
    a = a[,c('SNP','CHR','Position','Beta','SE','P','N','N_Cases','N_Controls')]
    
  }
  
  data.table::fwrite(a,file.path('NewGWAS',file_this),
         sep='\t',col.names=T,row.names=F)
  
  reformatted <- 
    MungeSumstats::format_sumstats(path=file.path('NewGWAS',file_this),
                                   ref_genome = c('GRCh38','GRCh38',
                                                  'GRCh37','GRCh37')[i],
                                   convert_ref_genome = 'GRCh38',
                                   save_format = 'LDSC',
                                   save_path = paste0(cancer,'.tsv.gz'),
                                   drop_indels = T,
                                   snp_ids_are_rs_ids = c(F,F,T,T)[i],
                                   nThread = 6,compute_n = 'ldsc',
                                   local_chain = '/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/liftover/hg19ToHg38.over.chain')
  
  
# setwd('/rsrch5/home/epi/bhattacharya_lab/projects/PanCan_isoTWAS/GWAS')
# require(data.table)
# 
# fff = list.files()
# fff = fff[grepl('.tsv.gz',fff)]
# 
# for (f in fff){
#   
#   a = fread(f)
#   if ('Neff' %in% colnames(a)){
#   a$N = a$Neff}
#   a = a[,-ncol(a),with=F]
#   fwrite(a,f,sep='\t',col.names=T,row.names=F,quote=F)
#   
# }
# 
