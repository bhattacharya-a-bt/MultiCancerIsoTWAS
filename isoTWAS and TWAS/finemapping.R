require(vroom)
require(data.table)
require(GenomicRanges)
require(Matrix)
require(bigsnpr)
require(xlsx)
require(isotwas)

setwd('/u/project/pasaniuc/abtbhatt')
all_res = fread('AllResults_isoTWAS_PanCan.tsv')

traits = unique(all_res$File)

for (tr in traits[9:12]){
  print(tr)
  
  if (file.exists(file.path(paste0(tr,'_isoTWAS_FOCUS_Pancan.tsv')))){
    
    file.remove(file.path(paste0(tr,'_isoTWAS_FOCUS_Pancan.tsv')))
    
  }
  
  twas = subset(all_res,File == tr)
  twas$pLI[is.na(twas$pLI)] = 0
  res = twas[order(twas$Chromosome,twas$Start),]
  res$Trait = tr
  res$GTA = paste(res$Gene,res$Tissue,sep=':')
  res$TTA = paste(res$Transcript,res$Tissue,sep=':')
  
  chr.table = table(res$Chromosome)
  chr.un = unique(res$Chromosome)
  keep.gta = c()
  keep.tta = c()
  for (c in chr.un){
    res.cur = subset(res,Chromosome == c)
    res.cur = res.cur[order(res.cur$Start),]
    if (nrow(res.cur) > 1){
      for (i in 1:(nrow(res.cur)-1)){
        if (res.cur$End[i] > res.cur$Start[i+1] - 1e6){
          keep.gta = unique(c(keep.gta,
                              c(res.cur$GTA[c(i,i+1)])))
          keep.tta = unique(c(keep.tta,
                              c(res.cur$TTA[c(i,i+1)])))
        }
      }
    }
  }
  
  new_res = subset(res,!TTA %in% keep.tta)
  if (nrow(new_res) > 0){
    new_res$pip = 1
    new_res$in_cred_set = T
    new_res$Group = paste0('ind',1:nrow(new_res))
    new_res$Overlap = 'No'
    new_res = new_res[,c('Indication',
                         'Gene','HGNC','Transcript',
                         'Chromosome',
                         'Start','End','Tissue','Z','Screening Adjusted P',
                         'Confirmation P',"Permutation P",
                         'Top GWAS SNP','Top GWAS P','pLI',
                         "pip",'in_cred_set','Group','Overlap')]
    fwrite(new_res,
           file.path(paste0(tr,'_isoTWAS_FOCUS_Pancan.tsv')),
           append=T,
           row.names= F,quote=F,sep='\t',na = 'NA')
  }
  
  
  res_base = subset(res,TTA %in% keep.tta)
  res_17q = subset(res_base,Chromosome == 17 & Start >= 42800000 &
                     End <= 46800000)
  
  if (nrow(res_17q) > 0){
    res_17q$pip = 1
    res_17q$in_cred_set = T
    res_17q$Group = paste0('ind',1:nrow(res_17q))
    res_17q$Overlap = 'No'
    res_17q = res_17q[,c('Indication',
                         'Gene','HGNC','Transcript',
                         'Chromosome',
                         'Start','End','Tissue','Z','Screening Adjusted P',
                         'Confirmation P',"Permutation P",
                         'Top GWAS SNP','Top GWAS P','pLI',
                         "pip",'in_cred_set','Group','Overlap')]
    fwrite(res_17q,
           file.path(paste0(tr,'_isoTWAS_FOCUS_Pancan.tsv')),
           append=T,
           row.names= F,quote=F,sep='\t',
           na = 'NA')
    res_base = subset(res_base,
                      !(paste0(res_base$Transcript,
                               res_base$Tissue) %in%
                          paste0(res_17q$Transcript,
                                 res_17q$Tissue)))
  }
  
  
  tissues_list = unique(res_base$Tissue)
  for (tissue in tissues_list){
    print(tissue)
    res = subset(res_base,Tissue == tissue)
    chr.un = as.numeric(unique(res$Chromosome))
  if (length(chr.un) >= 1){
    for (c in chr.un[1:length(chr.un)]){
      print(c)
      this.res.tot = subset(res,Chromosome == c)
      this.res.tot$Group = 1
      ggg = 1
      
      if (nrow(this.res.tot) == 1){
        this.res.tot$pip = 1
        this.res.tot$in_cred_set = T
        this.res.tot$Group = paste0('ind',1:nrow(this.res.tot))
        this.res.tot$Overlap = 'No'
        this.res.tot = this.res.tot[,c('Indication',
                                       'Gene','HGNC','Transcript',
                                       'Chromosome',
                                       'Start','End','Tissue','Z','Screening Adjusted P',
                                       'Confirmation P',"Permutation P",
                                       'Top GWAS SNP','Top GWAS P','pLI',
                                       "pip",'in_cred_set','Group','Overlap')]
        fwrite(this.res.tot,
               file.path(paste0(tr,'_isoTWAS_FOCUS_Pancan.tsv')),
               append=T,
               row.names= F,quote=F,sep='\t',
               na = 'NA')
      } else {
      
      for (i in 1:(nrow(this.res.tot)-1)){
        if (this.res.tot$End[i] <= this.res.tot$Start[i+1] - 1e6){
          ggg = ggg+1
          this.res.tot$Group[(i+1):nrow(this.res.tot)] = ggg
        }
      }
      for (g in unique(this.res.tot$Group)){
        print(paste0('Group ',g))
        this.res = subset(this.res.tot,
                          Group == g)
        if (nrow(this.res) == 1){
          this.res$pip = 1
          this.res$in_cred_set = T
          this.res$Group = paste0('ind',1:nrow(this.res))
          this.res$Overlap = 'No'
          this.res = this.res[,c('Indication',
                                 'Gene','HGNC','Transcript',
                                 'Chromosome',
                                 'Start','End','Tissue','Z','Screening Adjusted P',
                                 'Confirmation P',"Permutation P",
                                 'Top GWAS SNP','Top GWAS P','pLI',
                                 "pip",'in_cred_set','Group','Overlap')]
          fwrite(this.res,
                 file.path(paste0(tr,'_isoTWAS_FOCUS_Pancan.tsv')),
                 append=T,
                 row.names= F,quote=F,sep='\t',
                 na = 'NA')
          
        } else {
        all.snps = c()
        omega = c()
        pos = c()
        gene = c()
        snp.chr = c()
        for (i in 1:nrow(this.res)){
          aaa = readRDS(paste0('/u/project/pasaniuc/abtbhatt/GTEx/',
                               this.res$Tissue[i],
                               '/isoTWAS/',
                               this.res$Gene[i],
                               '_isoTWAS.RDS'))
          aaa = subset(aaa,
                       Feature == this.res$Transcript[i])
          Model = data.frame(SNP = aaa$SNP,
                             Chromosome = aaa$Chromosome,
                             Position = aaa$Position,
                             Effect = aaa$Weight,
                             A1 = aaa$ALT,
                             A2 = aaa$REF)
          Model = subset(Model,Effect!=0)
          Model = Model[!duplicated(Model$SNP),]
          all.snps = c(all.snps,
                       as.character(Model$SNP))
          omega = c(omega,
                    as.numeric(Model$Effect))
          gene = c(gene,
                   rep(this.res$TTA[i],nrow(Model)))
          snp.chr = c(snp.chr,
                      as.numeric(Model$Chromosome))
          pos = c(pos,as.numeric(Model$Position))
        }
        tot.df = data.frame(SNP = all.snps,
                            Gene = gene,
                            Effect = omega,
                            Chromosome = snp.chr)
        model.df = as.data.frame(matrix(nrow = length(unique(all.snps)),
                                        ncol = nrow(this.res)+1))
        colnames(model.df) = c('SNP',this.res$TTA)
        model.df$SNP = as.character(unique(all.snps))
        for (q in 1:nrow(this.res)){
          #print(this.res$TTA[g])
          cur.tot.df = subset(tot.df,Gene == this.res$TTA[q])
          cur.tot.df$SNP = as.character(cur.tot.df$SNP)
          for (i in 1:nrow(model.df)){
            #print(i)
            w = which(cur.tot.df$SNP == model.df$SNP[i])
            model.df[i,q+1] = ifelse(length(w) != 0,
                                     cur.tot.df$Effect[w],
                                     0)
          }
        }
        model.df$Chromosome = c
        for (i in 1:nrow(model.df)){
          rrr = subset(tot.df,SNP == model.df$SNP[i])
          model.df$Chromosome[i] = rrr$Chromosome[1]
        }
        
        min = max(c(1,min(pos)-1e6))
        max = max(pos)+1e6
        
        system(paste0('plink ',
                      ' --bfile /u/project/pasaniuc/pasaniucdata/GTEXv8_geno/GTEx.WGS.838.passOnly.geno0.05.hwe0.00001.dbsnp.SNPsOnly.NoAmbig.LDREF',
                      ' --chr ',c,
                      ' --from-bp ',min,
                      ' --to-bp ',max,
                      ' --make-bed --out /u/scratch/a/abtbhatt/temp'),
               ignore.stdout = T,
               ignore.stderr = T)
        
        snps = snp_attach(snp_readBed2(paste0('/u/scratch/a/abtbhatt/temp.bed'),
                                       backingfile = tempfile()))
        snp.set = subset(snps$map,
                         marker.ID %in% model.df$SNP)
        model.df = model.df[match(snp.set$marker.ID,
                                  model.df$SNP),]
        V = snp_cor(snp_attach(subset(snps,
                                      ind.col =
                                        which(snps$map$marker.ID %in%
                                                model.df$SNP),
                                      backingfile=tempfile()))$genotypes)
        
        Omega = Matrix(as.matrix(model.df[,-c(1,ncol(model.df))]))
        zscores = this.res$Z
        m = length(zscores)
        wcor = estimate_cor(as.matrix(Omega),
                            as.matrix(V),intercept=T)[[1]]
        diag(wcor) = 1
        wcor[is.na(wcor)] = 0
        swld = estimate_cor(as.matrix(Omega),
                            as.matrix(V),
                            intercept=T)[[2]]
        null_res = m * log(1 - 1e-3)
        marginal = m * log(1 - 1e-3)
        comb_list = list()
        for (n in 1:min(3,length(zscores))){
          comb_list = c(comb_list,
                        combn(1:length(zscores),n,simplify=F))
        }
        pips = rep(0,length(zscores))
        zscores = get_resid(zscores,as.matrix(swld),as.matrix(wcor))[[1]]
        for (j in 1:length(comb_list)){
          subset = comb_list[[j]]
          local = bayes_factor(zscores,
                               idx_set = subset,
                               wcor = wcor)
          marginal = log(exp(local) + exp(marginal))
          for (idx in subset){
            if (pips[idx] == 0){
              pips[idx] = local
            } else {
              pips[idx] = log(exp(pips[idx]) + exp(local))
            }
          }
          print(pips)
          print(marginal)
        }
        
        iter=2
        while (any(pips == Inf) & iter >= 1){
          
          zscores = this.res$Z
          m = length(zscores)
          null_res = m * log(1 - 1e-3)
          marginal = m * log(1 - 1e-3)
          comb_list = list()
          for (n in 1:min(iter,
                          length(zscores))){
            comb_list = c(comb_list,
                          combn(1:length(zscores),n,simplify=F))
          }
          
          pips = rep(0,length(zscores))
          for (j in 1:length(comb_list)){
            subset = comb_list[[j]]
            local = isotwas::bayes_factor(zscores,
                                          idx_set = subset,
                                          wcor = wcor)
            marginal = log(exp(local) + exp(marginal))
            for (idx in subset){
              if (pips[idx] == 0){
                pips[idx] = local
              } else {
                pips[idx] = log(exp(pips[idx]) + exp(local))
              }
            }
            print(pips)
            print(marginal)
          }
          iter = iter-1
        }
        
        pips = exp(pips - marginal)
        null_res = exp(null_res - marginal)
        this.res$pip = pips
        this.res = this.res[order(this.res$pip,decreasing = T),]
        npost = this.res$pip/sum(this.res$pip)
        csum = cumsum(npost)
        this.res$in_cred_set = F
        for (i in 1:nrow(this.res)){
          this.res$in_cred_set[i] = T
          if (i > 1){
            if (csum[i] > .9 & csum[i-1] < .9){
              this.res$in_cred_set[i] = T
            }
            if (csum[i] < .9){
              this.res$in_cred_set[i] = T
            }
            if (csum[i] > .9 & csum[i-1] > .9){
              this.res$in_cred_set[i] = F
            }
          }
        }
        
        this.res$Overlap = 'Yes'
        this.res = this.res[,c('Indication',
                               'Gene','HGNC','Transcript',
                               'Chromosome',
                               'Start','End','Tissue','Z','Screening Adjusted P',
                               'Confirmation P',"Permutation P",
                               'Top GWAS SNP','Top GWAS P','pLI',
                               "pip",'in_cred_set','Group','Overlap')]
        fwrite(this.res,
               file.path(paste0(tr,'_isoTWAS_FOCUS_Pancan.tsv')),
               append=T,
               row.names= F,quote=F,
               sep='\t',na = 'NA')
        }
      }
      }
  }
  }
  }
}


file.remove('isoTWAS_FineMap_PanCan.tsv')

for (tr in c(traits[1:10],
             '0_Cimba_Ovarian_serous.AUTO_Finalcleaned',
             '0_ovarian.European_overall.AUTO_Finalcleaned')){
  print(tr)
  all_res = as.data.frame(fread(file.path(paste0(tr,'_isoTWAS_FOCUS_Pancan.tsv'))))
  all_res = all_res[order(all_res$Chromosome,all_res$Start),]
  res = all_res
  
  fwrite(res[,c('Indication',
                           'Gene','HGNC','Transcript',
                           'Chromosome',
                           'Start','End','Tissue','Z','Screening Adjusted P',
                           'Confirmation P',"Permutation P",
                           'Top GWAS SNP','Top GWAS P','pLI',
                           "pip",'in_cred_set')],
           'isoTWAS_FineMap_PanCan.tsv',
           append=T,
           row.names=F,
           quote=F,
           sep='\t')
}

all_fm = fread('isoTWAS_FineMap_PanCan.tsv')
all_fm = all_fm[!duplicated(all_fm),]
all_fm$TTA = paste(all_fm$Transcript,all_fm$Tissue,all_fm$Indication,sep=':')
all_fm = all_fm[!duplicated(all_fm$TTA),]

fwrite(all_fm[,c('Indication',
                 'Gene','HGNC','Transcript',
                 'Chromosome',
                 'Start','End','Tissue','Z','Screening Adjusted P',
                 'Confirmation P',"Permutation P",
                 'Top GWAS SNP','Top GWAS P','pLI',
                 "pip",'in_cred_set')],
       'isoTWAS_FineMap_PanCan_090524.tsv',
       append=T,
       row.names=F,
       quote=F,
       sep='\t')









require(vroom)
require(data.table)
require(GenomicRanges)
require(Matrix)
require(bigsnpr)
require(xlsx)
require(isotwas)

setwd('/u/project/pasaniuc/abtbhatt')
all_res = fread('AllResults_TWAS_PanCan.tsv')
traits = unique(all_res$File)

for (tr in traits[9:12]){
  print(tr)

  if (file.exists(file.path(paste0(tr,'_TWAS_FOCUS_Pancan.tsv')))){

    file.exists(file.path(paste0(tr,'_TWAS_FOCUS_Pancan.tsv')))

  }

  twas = subset(all_res,File == tr)
  twas$pLI[is.na(twas$pLI)] = 0
  res = twas[order(twas$Chromosome,twas$Start),]
  res$Trait = tr
  res$Transcript = res$Gene
  res$GTA = paste(res$Gene,res$Tissue,sep=':')
  res$TTA = paste(res$Transcript,res$Tissue,sep=':')

  chr.table = table(res$Chromosome)
  chr.un = unique(res$Chromosome)
  keep.gta = c()
  keep.tta = c()
  for (c in chr.un){
    res.cur = subset(res,Chromosome == c)
    res.cur = res.cur[order(res.cur$Start),]
    if (nrow(res.cur) > 1){
      for (i in 1:(nrow(res.cur)-1)){
        if (res.cur$End[i] > res.cur$Start[i+1] - 1e6){
          keep.gta = unique(c(keep.gta,
                              c(res.cur$GTA[c(i,i+1)])))
          keep.tta = unique(c(keep.tta,
                              c(res.cur$TTA[c(i,i+1)])))
        }
      }
    }
  }

  new_res = subset(res,!TTA %in% keep.tta)
  if (nrow(new_res) > 0){
    new_res$pip = 1
    new_res$in_cred_set = T
    new_res$Group = paste0('ind',1:nrow(new_res))
    new_res$Overlap = 'No'
    new_res = new_res[,c('Indication',
                         'Gene','HGNC',
                         'Chromosome',
                         'Start','End','Tissue','Z','FDR',"Permutation P",
                         'Top GWAS SNP','Top GWAS P','pLI',
                         "pip",'in_cred_set','Group','Overlap')]
    fwrite(new_res,
           file.path(paste0(tr,'_TWAS_FOCUS_Pancan.tsv')),
           append=T,
           row.names= F,quote=F,sep='\t',
           na = NA)
  }


  res_base = subset(res,TTA %in% keep.tta)
  tissues_list = unique(res_base$Tissue)
  for (tissue in tissues_list){
    print(tissue)
  res = subset(res_base,Tissue == tissue)
  chr.un = as.numeric(unique(res$Chromosome))
  if (length(chr.un) >= 1){
    for (c in chr.un[1:length(chr.un)]){
      print(c)
      this.res.tot = subset(res,Chromosome == c)
      this.res.tot$Group = 1
      ggg = 1

      if (nrow(this.res.tot) == 1){
        this.res.tot$pip = 1
        this.res.tot$in_cred_set = T
        this.res.tot$Group = paste0('ind',1:nrow(this.res.tot))
        this.res.tot$Overlap = 'No'
        this.res.tot = this.res.tot[,c('Indication',
                                       'Gene','HGNC',
                                       'Chromosome',
                                       'Start','End','Tissue','Z','FDR',"Permutation P",
                                       'Top GWAS SNP','Top GWAS P','pLI',
                                       "pip",'in_cred_set','Group','Overlap')]
        fwrite(this.res.tot,
               file.path(paste0(tr,'_TWAS_FOCUS_Pancan.tsv')),
               append=T,
               row.names= F,quote=F,sep='\t',
               na = NA)
      } else {

        for (i in 1:(nrow(this.res.tot)-1)){
          if (this.res.tot$End[i] <= this.res.tot$Start[i+1] - 1e6){
            ggg = ggg+1
            this.res.tot$Group[(i+1):nrow(this.res.tot)] = ggg
          }
        }
        for (g in unique(this.res.tot$Group)){
          print(paste0('Group ',g))
          this.res = subset(this.res.tot,
                            Group == g)
          if (nrow(this.res) == 1){
            this.res$pip = 1
            this.res$in_cred_set = T
            this.res$Group = paste0('ind',1:nrow(this.res))
            this.res$Overlap = 'No'
            this.res = this.res[,c('Indication',
                                   'Gene','HGNC',
                                   'Chromosome',
                                   'Start','End','Tissue','Z','FDR',"Permutation P",
                                   'Top GWAS SNP','Top GWAS P','pLI',
                                   "pip",'in_cred_set','Group','Overlap')]
            fwrite(this.res,
                   file.path(paste0(tr,'_TWAS_FOCUS_Pancan.tsv')),
                   append=T,
                   row.names= F,quote=F,sep='\t',
                   na = NA)

          } else {
            all.snps = c()
            omega = c()
            pos = c()
            gene = c()
            snp.chr = c()
            for (i in 1:nrow(this.res)){
              aaa = readRDS(paste0('/u/project/pasaniuc/abtbhatt/GTEx/',
                                   this.res$Tissue[i],
                                   '/TWAS/',
                                   this.res$Gene[i],
                                   '_TWAS.RDS'))
              aaa = subset(aaa,
                           Feature == this.res$Transcript[i])
              Model = data.frame(SNP = aaa$SNP,
                                 Chromosome = aaa$Chromosome,
                                 Position = aaa$Position,
                                 Effect = aaa$Weight,
                                 A1 = aaa$ALT,
                                 A2 = aaa$REF)
              Model = subset(Model,Effect!=0)
              Model = Model[!duplicated(Model$SNP),]
              all.snps = c(all.snps,
                           as.character(Model$SNP))
              omega = c(omega,
                        as.numeric(Model$Effect))
              gene = c(gene,
                       rep(this.res$TTA[i],nrow(Model)))
              snp.chr = c(snp.chr,
                          as.numeric(Model$Chromosome))
              pos = c(pos,as.numeric(Model$Position))
            }
            tot.df = data.frame(SNP = all.snps,
                                Gene = gene,
                                Effect = omega,
                                Chromosome = snp.chr)
            model.df = as.data.frame(matrix(nrow = length(unique(all.snps)),
                                            ncol = nrow(this.res)+1))
            colnames(model.df) = c('SNP',this.res$TTA)
            model.df$SNP = as.character(unique(all.snps))
            for (q in 1:nrow(this.res)){
              #print(this.res$TTA[g])
              cur.tot.df = subset(tot.df,Gene == this.res$TTA[q])
              cur.tot.df$SNP = as.character(cur.tot.df$SNP)
              for (i in 1:nrow(model.df)){
                #print(i)
                w = which(cur.tot.df$SNP == model.df$SNP[i])
                model.df[i,q+1] = ifelse(length(w) != 0,
                                         cur.tot.df$Effect[w],
                                         0)
              }
            }
            model.df$Chromosome = c
            for (i in 1:nrow(model.df)){
              rrr = subset(tot.df,SNP == model.df$SNP[i])
              model.df$Chromosome[i] = rrr$Chromosome[1]
            }

            min = max(c(1,min(pos)-1e6))
            max = max(pos)+1e6

            system(paste0('plink ',
                          ' --bfile /u/project/pasaniuc/pasaniucdata/GTEXv8_geno/GTEx.WGS.838.passOnly.geno0.05.hwe0.00001.dbsnp.SNPsOnly.NoAmbig.LDREF',
                          ' --chr ',c,
                          ' --from-bp ',min,
                          ' --to-bp ',max,
                          ' --make-bed --out /u/scratch/a/abtbhatt/temptwas'),
                   ignore.stdout = T,
                   ignore.stderr = T)

              snps = snp_attach(snp_readBed2(paste0('/u/scratch/a/abtbhatt/temptwas.bed'),
                                             backingfile = tempfile()))
              snp.set = subset(snps$map,
                               marker.ID %in% model.df$SNP)
              model.df = model.df[match(snp.set$marker.ID,
                                        model.df$SNP),]
              V = snp_cor(snp_attach(subset(snps,
                                            ind.col =
                                              which(snps$map$marker.ID %in%
                                                      model.df$SNP),
                                            backingfile=tempfile()))$genotypes)
  
              Omega = Matrix(as.matrix(model.df[,-c(1,ncol(model.df))]))
              zscores = this.res$Z
              m = length(zscores)
              wcor = estimate_cor(as.matrix(Omega),
                                  as.matrix(V),intercept=T)[[1]]
              diag(wcor) = 1
              wcor[is.na(wcor)] = 0
              swld = estimate_cor(as.matrix(Omega),
                                  as.matrix(V),
                                  intercept=T)[[2]]
              null_res = m * log(1 - 1e-3)
              marginal = m * log(1 - 1e-3)
              comb_list = list()
              for (n in 1:min(3,length(zscores))){
                comb_list = c(comb_list,
                              combn(1:length(zscores),n,simplify=F))
              }
              pips = rep(0,length(zscores))
              zscores = get_resid(zscores,as.matrix(swld),as.matrix(wcor))[[1]]
              for (j in 1:length(comb_list)){
                subset = comb_list[[j]]
                local = bayes_factor(zscores,
                                     idx_set = subset,
                                     wcor = wcor)
                marginal = log(exp(local) + exp(marginal))
                for (idx in subset){
                  if (pips[idx] == 0){
                    pips[idx] = local
                  } else {
                    pips[idx] = log(exp(pips[idx]) + exp(local))
                  }
                }
                print(pips)
                print(marginal)
              }

            iter=2
            while (any(pips == Inf) & iter >= 1){

              zscores = this.res$Z
              m = length(zscores)
              null_res = m * log(1 - 1e-3)
              marginal = m * log(1 - 1e-3)
              comb_list = list()
              for (n in 1:min(iter,
                              length(zscores))){
                comb_list = c(comb_list,
                              combn(1:length(zscores),n,simplify=F))
              }

              pips = rep(0,length(zscores))
              for (j in 1:length(comb_list)){
                subset = comb_list[[j]]
                local = isotwas::bayes_factor(zscores,
                                              idx_set = subset,
                                              wcor = wcor)
                marginal = log(exp(local) + exp(marginal))
                for (idx in subset){
                  if (pips[idx] == 0){
                    pips[idx] = local
                  } else {
                    pips[idx] = log(exp(pips[idx]) + exp(local))
                  }
                }
                print(pips)
                print(marginal)
              }
              iter = iter-1
            }
            
            if (is.infinite(marginal)){
              marginal = max(pips)+1
            }

            pips = exp(pips - marginal)
            null_res = exp(null_res - marginal)
            this.res$pip = pips
            this.res = this.res[order(this.res$pip,decreasing = T),]
            npost = this.res$pip/sum(this.res$pip)
            csum = cumsum(npost)
            this.res$in_cred_set = F
            for (i in 1:nrow(this.res)){
              this.res$in_cred_set[i] = T
              if (i > 1){
                if (csum[i] > .9 & csum[i-1] < .9){
                  this.res$in_cred_set[i] = T
                }
                if (csum[i] < .9){
                  this.res$in_cred_set[i] = T
                }
                if (csum[i] > .9 & csum[i-1] > .9){
                  this.res$in_cred_set[i] = F
                }
              }
            }

            this.res$Overlap = 'Yes'
            this.res = this.res[,c('Indication',
                                   'Gene','HGNC',
                                   'Chromosome',
                                   'Start','End','Tissue','Z','FDR',"Permutation P",
                                   'Top GWAS SNP','Top GWAS P','pLI',
                                   "pip",'in_cred_set','Group','Overlap')]
            fwrite(this.res,
                   file.path(paste0(tr,'_TWAS_FOCUS_Pancan.tsv')),
                   append=T,
                   row.names= F,quote=F,
                   sep='\t',na = NA)
          }
        }
      }
    }
    }
  }

}

file.remove('TWAS_FineMap_PanCan.tsv')

for (tr in traits){
  
  all_res = fread(file.path(paste0(tr,'_TWAS_FOCUS_Pancan.tsv')))
  all_res = all_res[order(all_res$Chromosome,all_res$Start)]
  res = all_res
  
  for (c in unique(res$Chromosome)){
    this.res.tot = subset(res,Chromosome == c)
    this.res.tot$Group = paste0('Chromosome',c,':','Group',1)
    ggg = 1
    if (nrow(this.res.tot) > 1){
    for (i in 1:(nrow(this.res.tot)-1)){
      if (this.res.tot$End[i] <= this.res.tot$Start[i+1] - 1e6){
        ggg = ggg+1
        this.res.tot$Group[(i+1):nrow(this.res.tot)] = 
          paste0('Chromosome',c,':','Group',ggg)
      }
    }
    }
    fwrite(this.res.tot[,c('Indication',
                           'Gene','HGNC',
                           'Chromosome',
                           'Start','End','Tissue','Z','FDR',
                           'Permutation P','Top GWAS SNP',
                           'Top GWAS P','pLI','pip','in_cred_set',
                           'Group')],
           'TWAS_FineMap_PanCan.tsv',
           append=T,
           row.names=F,
           quote=F,
           sep='\t')
  }
}

all_fm = fread('TWAS_FineMap_PanCan.tsv')
all_fm = all_fm[!duplicated(all_fm),]
all_fm$GTA = paste(all_fm$Gene,all_fm$Tissue,all_fm$Indication,sep=':')
all_fm = all_fm[order(all_fm$pip,decreasing=T),]
all_fm = all_fm[!duplicated(all_fm$GTA),]

fwrite(all_fm[,c('Indication',
                       'Gene','HGNC',
                       'Chromosome',
                       'Start','End','Tissue','Z','FDR',
                       'Permutation P','Top GWAS SNP',
                       'Top GWAS P','pLI','pip','in_cred_set',
                       'Group')],
       'TWAS_FineMap_PanCan_090524.tsv',
       append=T,
       row.names=F,
       quote=F,
       sep='\t')
