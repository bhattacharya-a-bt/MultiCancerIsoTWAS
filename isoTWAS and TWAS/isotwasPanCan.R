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
require(isotwas)

manifest = as.data.frame(
  fread('/u/project/pasaniuc/yhc1998/yhc/isoTWAS/twas_cancer_tissue_full.tsv',header=T)
)
cancer = manifest$Cancer[index]
tissue = manifest$Tissue[index]

folder="/u/project/pasaniuc/abtbhatt/GTEx"
list_tissue = list.files(folder)



### GATHER TWAS/ISOTWAS FOLDER
isotwas_folder=paste0("/u/project/pasaniuc/abtbhatt/GTEx/",tissue,"/TWAS")

# DEFINE LD FILE
ld_folder = "/u/project/pasaniuc/pasaniucdata/GTEXv8_geno/GTEx.WGS.838.passOnly.geno0.05.hwe0.00001.dbsnp.SNPsOnly.NoAmbig.LDREF.bed" 

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
#dir.create(file.path('/u/project/pasaniuc/yhc1998/yhc/isoTWAS/TISSUE/Associations',
#                     tr),recursive = T)

# DEFINE OUTPUT FOR TWAS AND ISOTWAS SUMM STATS
outFile_twas = file.path('/u/project/pasaniuc/yhc1998/yhc/isoTWAS/TWAS',
                         tr,tissue,
                         paste0(tr,'_TWAS.tsv'))
outFile_isotwas = file.path('/u/project/pasaniuc/yhc1998/yhc/isoTWAS/TWAS/isoTWAS',
                            tr,tissue,
                            paste0(tr,'_isoTWAS.tsv'))



fff = list.files(isotwas_folder)
toDoGenes = sapply(strsplit(fff,'_iso'),
                   function(x) x[1])




if (file.exists(outFile_isotwas) &
    file.exists(outFile_twas)){
  ttt = fread(outFile_twas)
  iii = fread(outFile_isotwas)
  toDoGenes = toDoGenes[!toDoGenes %in% c(ttt$Gene,
                                          iii$Gene)]
}



doneFile = 
  paste0('/u/project/pasaniuc/yhc1998/yhc/isoTWAS/TISSUE/',
         tr,"_",tissue,'_done_dev.tsv')

file.remove(doneFile)

if (!exists('ttt') & !exists('iii')){
  
  
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
                    Gene = unique(c(ttt$Gene,
                                    iii$Gene)),
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
    
    ### RUN ISOTWAS
    isotwas_model = readRDS(file.path(isotwas_folder,
                                      paste0(gene,'_isoTWAS.RDS')))
    isotwas_model$R2 = unlist(isotwas_model$R2)
    
    if (file.exists(file.path(twas_folder,
                              paste0(gene,'_TWAS.RDS')))){
      twas_model = readRDS(file.path(twas_folder,
                                     paste0(gene,'_TWAS.RDS')))
    } else {
      twas_model = as.data.frame(isotwas_model[1:2,])
      twas_model$R2 = 0
    }
    
    if (any(twas_model$R2 > 0.01) |
        any(unlist(isotwas_model$R2) > 0.01)){
      
      
      ### USE PLINK TO CREATE A SMALLER .BED
      
      ### USE BED_COR() to create ld matrix
      if (file.exists(file.path(paste("/u/project/pasaniuc/yhc1998/yhc/isoTWAS/LD_matrice/",tissue,"/LD_isotwas/isotwas_",gene, ".bed", sep="")))){
        LD_isotwas=snp_attach(snp_readBed2(paste("/u/project/pasaniuc/yhc1998/yhc/isoTWAS/LD_matrice/",tissue,"/LD_isotwas/isotwas_",gene, ".bed", sep=""),
                                           backingfile = tempfile()))
        snpnames = LD_isotwas$map$marker.ID
        LD_isotwas = as.matrix(snp_cor(LD_isotwas$genotypes))
        colnames(LD_isotwas) = rownames(LD_isotwas) = snpnames
        
      } else {
        LD_twas=snp_attach(snp_readBed2(paste("/u/project/pasaniuc/yhc1998/yhc/isoTWAS/LD_matrice/",tissue,"/LD_isotwas/twas_",gene, ".bed", sep=""),backingfile = tempfile()))
        LD_twas = snp_cor(LD_twas$genotypes)
        snpnames = LD_twas$map$marker.ID
        LD_twas = as.matrix(snp_cor(LD_twas$genotypes))
        colnames(LD_twas) = rownames(LD_twas) = snpnames
      }
      
      sumstats.cur = subset(sumstats,SNP %in% unique(c(twas_model$SNP,
                                                       isotwas_model$SNP)))
      
      
      
      
      tot = rbind(twas_model[,c('SNP','Chromosome','Position')],
                  isotwas_model[,c('SNP','Chromosome','Position')])
      tot = tot[!duplicated(tot$SNP),]
      sumstats.cur = merge(sumstats.cur,tot,by = 'SNP')
      #sumstats.cur$Beta = sumstats.cur$Z
      twas_model = subset(twas_model,SNP %in% sumstats$SNP)
      isotwas_model = subset(isotwas_model,
                             SNP %in% sumstats$SNP)
      twas_model$A1 = twas_model$REF
      twas_model$A2 = twas_model$ALT
      twas_model$Transcript = twas_model$Feature
      if (nrow(twas_model) > 1){twas_model = as.data.frame(apply(twas_model,2,unlist))}
      twas_model$Chromosome = as.numeric(twas_model$Chromosome)
      twas_model$Weight = as.numeric(twas_model$Weight)
      twas_model$R2 = as.numeric(twas_model$R2)
      twas_model$Position = as.numeric(twas_model$Position)
      isotwas_model$A1 = isotwas_model$REF
      isotwas_model$A2 = isotwas_model$ALT
      isotwas_model$Transcript = isotwas_model$Feature
      print(max(twas_model$R2))
      print(max(isotwas_model$R2))
      
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
      
      
      isotwas_model = subset(isotwas_model,R2 >= 0.01)
      for (tx in unique(isotwas_model$Transcript)){
        
        if (nrow(isotwas_model) > 0){
          #print(tx)
          tx_df_isotwas = burdenTest(mod = subset(isotwas_model,
                                                  Transcript == tx),
                                     ld = LD_isotwas,
                                     gene = gene,
                                     sumStats = sumstats.cur,
                                     chr = 'Chromosome',
                                     pos = 'Position',
                                     a1 = 'A1',
                                     a2 = 'A2',
                                     Z = 'Z',
                                     beta = 'Beta',
                                     se = 'SE',
                                     R2cutoff = 0.01,
                                     alpha = 1e-3,
                                     nperms = 1e3,
                                     usePos = F)
          if (class(tx_df_isotwas) == 'data.frame'){
            
            tx_df_isotwas$R2 = subset(isotwas_model,
                                      Transcript == tx)$R2[1]
            colnames(tx_df_isotwas) = c('Gene','Transcript',
                                        'Z','P','permute.P','topSNP',
                                        'topSNP.P','R2')
            
            fwrite(tx_df_isotwas,outFile_isotwas,
                   append = T, sep = '\t',
                   quote = F, row.names=F)
            
          }
        }
        
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
