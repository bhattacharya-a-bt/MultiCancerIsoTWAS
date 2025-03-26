library("optparse")
option_list <- list(
  make_option(c("-i", "--index"), action="store_true", default=TRUE,
              help="index",type='numeric')
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

##ceiling = 43
##number = 1000

index = opt$index
#number = opt$number


manifest = as.data.frame(
  fread('/u/project/pasaniuc/yhc1998/yhc/isoTWAS/twas_cancer_tissue_full.tsv',header=T)
)
cancer = manifest$Cancer[index]
tissue = manifest$Tissue[index]

folder="/u/project/pasaniuc/abtbhatt/GTEx"
list_tissue = list.files(folder)


#' Compute weighted burden test
#'
#' The function takes in a gene expression model in MOSTWAS form
#' and GWAS summary statistics and
#' carries out the weighted burden Z-test for a trait
#'
#' @param mod data.frame, model for a given isoform
#' @param ld matrix, ld reference matrix
#' @param gene character, gene name
#' @param sumStats data frame, GWAS summary statistics
#' @param chr character, colnames in sumStats that keeps the chromosome
#' @param pos character, colnames in sumStats that keeps the position
#' @param a1 character, colnames in sumStats that keeps the ALT allele
#' @param a2 character, colnames in sumStats that keeps the REF allele
#' @param Z character, colnames in sumStats that keeps the Z score
#' @param beta character, colnames in sumStats that keeps the effect size
#' @param se character, colnames in sumStats that keeps the standard error
#' @param R2cutoff numeric, predictive R2 cutoff
#' @param alpha numeric, P-value threshold for permutation testing
#' @param nperms numeric, number of permutations
#' @param usePos logical, use SNP positions vs. SNP ids
#'
#' @return list of results for burden and permutation tests
#'
#' @importFrom boot boot
#' @importFrom stats pnorm
#'
#' @export
burdenTest <- function(mod,
                       ld,
                       gene,
                       sumStats,
                       chr,
                       pos,
                       a1,
                       a2,
                       Z = NULL,
                       beta = NULL,
                       se = NULL,
                       R2cutoff = 0.01,
                       alpha = 2.5e-6,
                       nperms = 1e3,
                       usePos = F){
  
  
  if (all(is.null(c(Z,beta,se)))){
    stop('Please provide a column name for the Z-score or beta and SE.')
  }
  
  if (is.null(Z) | any(is.null(c(beta,se)))){
    stop('Please provide a column name for the Z-score or beta and SE.')
  }
  
  colnames(sumStats)[which(colnames(sumStats) == chr)] = 'Chromosome'
  colnames(sumStats)[which(colnames(sumStats) == pos)] = 'Position'
  colnames(sumStats)[which(colnames(sumStats) == a1)] = 'A1'
  colnames(sumStats)[which(colnames(sumStats) == a2)] = 'A2'
  
  if (!is.null(Z)){
    colnames(sumStats)[which(colnames(sumStats) == Z)] = 'Z'
  }
  
  if (!all(is.null(c(beta,se)))){
    colnames(sumStats)[which(colnames(sumStats) == beta)] = 'Beta'
    colnames(sumStats)[which(colnames(sumStats) == se)] = 'SE'
  }
  
  if (!'Z' %in% colnames(sumStats)){
    sumStats$Z = sumStats$Beta/sumStats$SE
  }
  
  if (mod$R2[1] <= R2cutoff){
    return(paste0('The isoform is not predicted at R2 > ',
                  R2cutoff))
  }
  
  if (usePos){
    
    sumStats$SNP = paste(sumStats$Chromosome,sumStats$Position,sep=':')
    mod$SNP = paste(mod$Chromosome,mod$Position,sep=':')
    
  }
  
  tot = merge(mod,sumStats,by = 'SNP')
  
  if (nrow(tot) == 0){
    return('SNPs not found.')
  }
  
  tot$Z = ifelse(tot$A1.x == tot$A1.y,
                 tot$Z,
                 -1 * tot$Z)
  
  calculateTWAS <- function(effects,
                            Z,
                            LD,
                            indices){
    effects = effects[indices]
    twasZ = as.numeric(effects %*% Z)
    twasr2pred = as.numeric(effects %*% LD %*% effects)
    if (twasr2pred > 0){
      twas = as.numeric(twasZ/sqrt(twasr2pred))
    } else {
      twas = 0
    }
    return(twas)
  }
  
  twasLD = as.numeric(tot$Weight %*% tot$Z) /
    sqrt(as.numeric(tot$Weight) %*% ld[tot$SNP,tot$SNP] %*% as.numeric(tot$Weight))
  twasLD = as.numeric(twasLD)
  P = 2*stats::pnorm(-abs(twasLD))
  
  if (P <= alpha){
    permutationLD = boot::boot(data = tot$Weight,
                               statistic = calculateTWAS,
                               R = nperms,
                               sim = 'permutation',
                               Z = tot$Z,
                               LD = ld[tot$SNP,tot$SNP])
    permute.p = (nperms * mean(abs(permutationLD$t) >
                                 abs(permutationLD$t0)) + 1)/(nperms+1)
  } else {
    permute.p = 1}
  
  return(data.frame(Gene = gene,
                    Transcript = tot$Transcript[1],
                    Z = twasLD,
                    P = 2*pnorm(-abs(twasLD)),
                    permute.P = permute.p,
                    topSNP = tot$SNP[which.max(abs(tot$Z))],
                    topSNP.P = 2*pnorm(-abs(max(abs(tot$Z))))))
}


### GATHER TWAS/ISOTWAS FOLDER
twas_folder=paste0("/u/project/pasaniuc/abtbhatt/GTEx/",tissue,"/TWAS")


# DEFINE GWAS FOLDER/FILE
gwas_folder="/u/project/pasaniuc/ARCHIVE/amajumda/Cancers/analyses1/0_Raw"
fff = list.files(gwas_folder)
fff = fff[grepl('_MAF0.01',fff)]

trait = sapply(strsplit(fff,'_MAF'),
               function(x) x[1])



tr = trait[trait == cancer]
sumstats_file = file.path(gwas_folder,
                          fff[grep(cancer,fff)])

#library(isotwas)
#library(MOSTWAS)
library(data.table)
library(bigsnpr)


print(tr)

#READ IN GWAS SUMMSTATS 
bim = fread('/u/project/pasaniuc/pasaniucdata/GTEXv8_geno/GTEx.WGS.838.passOnly.geno0.05.hwe0.00001.dbsnp.SNPsOnly.NoAmbig.LDREF.bim')
sumstats = fread(sumstats_file)
sumstats$SNP = sumstats$rsid
sumstats = subset(sumstats,SNP %in% bim$V2)
rm(bim)
sumstats$A1 = toupper(sumstats$A1)
sumstats$A2 = toupper(sumstats$A2)


# DEFINE OUTPUT FOR TWAS AND ISOTWAS SUMM STATS
#dir.create(file.path('/u/project/pasaniuc/yhc1998/yhc/isoTWAS/TWAS',
#                     tr,tissue), recursive = T)

outFile_twas = file.path('/u/project/pasaniuc/yhc1998/yhc/isoTWAS/TWAS',
                         tr,tissue,
                         paste0(tr,'_TWAS.tsv'))


fff = list.files(twas_folder)
toDoGenes = sapply(strsplit(fff,'_TWAS'),
                   function(x) x[1])




if (file.exists(outFile_twas)){
  ttt = fread(outFile_twas)
  toDoGenes = toDoGenes[!toDoGenes %in% c(ttt$Gene)]
}



doneFile = 
  paste0('/u/project/pasaniuc/yhc1998/yhc/isoTWAS/TISSUE/',
         tr,"_",tissue,'_done_dev.tsv')

file.remove(doneFile)

if (!exists('ttt')){
  
  
  fwrite(data.frame(Trait = tr,
                    Gene = 'test',
                    Index = 1),
         doneFile,
         sep='\t',
         col.names=T,
         row.names=F,
         quote=F)
  
}

if (!file.exists(doneFile)){
  
  fwrite(data.frame(Trait = tr,
                    Gene = unique(c(ttt$Gene)),
                    Index = 1),
         doneFile,
         sep='\t',
         col.names=T,
         row.names=F,
         quote=F)
  
}











runAssociation = function(gene){
  
  print(gene)
  
  if (nrow(subset(vroom::vroom(doneFile,show_col_types=F),
                  Gene == gene & Trait == tr)) == 0){
    
    if (file.exists(file.path(twas_folder,
                              paste0(gene,'_TWAS.RDS')))){
      twas_model = readRDS(file.path(twas_folder,
                                     paste0(gene,'_TWAS.RDS')))
    } else {
      twas_model = as.data.frame(isotwas_model[1:2,])
      twas_model$R2 = 0
    }
    
    if (any(twas_model$R2 > 0.01)){
      
      
      ### USE PLINK TO CREATE A SMALLER .BED
      
      ### USE BED_COR() to create ld matrix
      snps_ld = snp_attach(
        snp_readBed2(paste0('/u/project/pasaniuc/abtbhatt/LD_out/',
                            gene,
                            '.bed'),
                     backingfile = tempfile()))
      LD_twas = bigsnpr::snp_cor(snps_ld$genotypes)
      colnames(LD_twas) = rownames(LD_twas) = snps_ld$map$marker.ID
      
      
      sumstats.cur = subset(sumstats,SNP %in% unique(twas_model$SNP))

      tot = twas_model[,c('SNP','Chromosome','Position')]
      tot = tot[!duplicated(tot$SNP),]
      sumstats.cur = merge(sumstats.cur,tot,by = 'SNP')
      #sumstats.cur$Beta = sumstats.cur$Z
      twas_model = subset(twas_model,SNP %in% sumstats$SNP)
      twas_model$A1 = twas_model$REF
      twas_model$A2 = twas_model$ALT
      twas_model$Transcript = twas_model$Feature
      if (nrow(twas_model) > 1){twas_model = as.data.frame(apply(twas_model,2,unlist))}
      twas_model$Chromosome = as.numeric(twas_model$Chromosome)
      twas_model$Weight = as.numeric(twas_model$Weight)
      twas_model$R2 = as.numeric(twas_model$R2)
      twas_model$Position = as.numeric(twas_model$Position)
      print(max(twas_model$R2))
      
      ### RUN TWAS ASSOCIATION
      if (nrow(twas_model) > 0){
        gene_df_twas = burdenTest(mod = twas_model,
                                  ld = LD_twas,
                                  gene = gene,
                                  sumStats = sumstats.cur,
                                  chr = 'Chromosome',
                                  pos = 'Position',
                                  a1 = 'A1',
                                  a2 = 'A2',
                                  Z = 'Z',
                                  beta = 'Beta',
                                  se = 'SE',
                                  R2cutoff = .01,
                                  alpha = 1e-3,
                                  nperms = 1e3,
                                  usePos = F)
        if (class(gene_df_twas) == 'data.frame'){
          gene_df_twas$R2 = unlist(twas_model$R2)[1]
          fwrite(gene_df_twas,outFile_twas,
                 append = T, sep = '\t',
                 quote = F, row.names=F)
        }
      }
    
    
    fwrite(data.frame(Trait = tr,
                      Gene = gene,
                      Index = 1),
           doneFile,
           sep='\t',
           row.names=F,
           quote=F,append=T)
  }
  }
}


require(pbmcapply)
tc_assoc = function(s){
  tryCatch(runAssociation(s),
           error = function(e) {
             print(paste0('error with ',s))
           }
  )
}

lapply(toDoGenes,
       tc_assoc)
