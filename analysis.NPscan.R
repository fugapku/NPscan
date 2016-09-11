NPscan(species = 'A.nidu_AspGD', out.file = 'A.nidu.xls', remove.glycan=T,
       nNPG= 2, nOG = 1, do.viterbi = T, EM.method ='null',pseudocount = 0.1, 
       dom.freq.cutoff = 2)
NPscan(species = 'A.nidu_AspGD', out.file = 'A.nidu.xls', remove.glycan=T,
       nNPG= 2, nOG = 1, do.viterbi = T, EM.method ='null',pseudocount = 0.1, 
       dom.freq.cutoff = 1)
NPscan(species = 'A.nidu_AspGD', out.file = 'A.nidu.xls', remove.glycan=T,
       nNPG= 2, nOG = 2, do.viterbi = T, EM.method ='null',pseudocount = 0.1, 
       dom.freq.cutoff = 2)
NPscan(species = 'A.nidu_AspGD', out.file = 'A.nidu.xls', remove.glycan=T,
       nNPG= 2, nOG = 2, do.viterbi = T, EM.method ='null',pseudocount = 0.1, 
       dom.freq.cutoff = 1)
NPscan(species = 'A.nidu_AspGD', out.file = 'A.nidu.xls', remove.glycan=T,
       nNPG= 1, nOG = 1, do.viterbi = T, EM.method ='null',pseudocount = 0.1, 
       dom.freq.cutoff = 2)
NPscan(species = 'A.nidu_AspGD', out.file = 'A.nidu.xls', remove.glycan=T,
       nNPG= 1, nOG = 1, do.viterbi = T, EM.method ='null',pseudocount = 0.1, 
       dom.freq.cutoff = 1)

# 20160413
NPscan(species = 'A.nidu_AspGD', out.file = 'A.nidu.xls', remove.glycan=T,init.truth = 'both', init.expand = 0, init.weighted = F,
       eval.truth = 'both',nNPG= 2, nOG = 1, do.viterbi = T, EM.method ='null',pseudocount = 0.1, 
       dom.freq.cutoff = 1) # AUC .895
NPscan(species = 'A.nidu_AspGD', out.file = 'A.nidu.xls', remove.glycan=T,init.truth = 'both', init.expand = 0, init.weighted = F,
       eval.truth = 'both',nNPG= 1, nOG = 1, do.viterbi = T, EM.method ='null',pseudocount = 0.1, 
       dom.freq.cutoff = 1) # AUC = 0.900
NPscan(species = 'A.nidu_AspGD', out.file = 'A.nidu.xls', remove.glycan=T,init.truth = 'both', init.expand = 0, init.weighted = F,
       eval.truth = 'both',nNPG= 2, nOG = 2, do.viterbi = T, EM.method ='null',pseudocount = 0.1, 
       dom.freq.cutoff = 1) # AUC = 0.903
NPscan(species = 'A.nidu_AspGD', out.file = 'A.nidu.xls', remove.glycan=T,init.truth = 'IPR', init.expand = 1, init.weighted = F,
       eval.truth = 'note',nNPG= 2, nOG = 1, do.viterbi = T, EM.method ='null',pseudocount = 0.1, 
       dom.freq.cutoff = 1) # AUC 0.865
NPscan(species = 'A.nidu_AspGD', out.file = 'A.nidu.xls', remove.glycan=T,init.truth = 'IPR', init.expand = 3, init.expand.combine = 'prob', init.weighted = F,
       eval.truth = 'note',nNPG= 2, nOG = 1, do.viterbi = T, EM.method ='null',pseudocount = 0.1, 
       dom.freq.cutoff = 1, predict.extra=T) # AUC 

for (et in c('both')){ # , 'domain'
  for (i in c('IPR')){ # , 'MC6', 'note', 'domain'
    for (w in c(F, T)){
      for (r in c(T)){
        for (ex in c(0,1,2,4)){
          NPscan(species = 'A.nidu_AspGD', out.file = 'A.nidu.xls', remove.glycan=r,
                 init.truth = i, init.expand = ex, init.weighted = w, eval.truth = et,
                 nNPG= 2, nOG = 1, do.viterbi = T, EM.method ='null',pseudocount = 0.1, 
                 dom.freq.cutoff = 1)   
          
        }
      }
    }
  }
}

for (et in c('note')){ # , 'domain'
  for (i in c('IPR')){ # , 'MC6', 'note', 'domain'
    for (w in c(T)){
      for (r in c(T)){
        for (ex in c(1,2,4,0)){
          for (dr in rev(c(3,5,7,9)/10)){
            NPscan(species = 'A.nidu_AspGD', out.file = 'A.nidu.xls', remove.glycan=r,init.decay.rate = dr,
                   init.truth = i, init.expand = ex, init.weighted = w, eval.truth = et,remove.pseudocounts=F, 
                   nNPG= 2, nOG = 1, do.viterbi = T, EM.method ='null',pseudocount = 0.1, 
                   dom.freq.cutoff = 1)   
          }
        }
      }
    }
  }
}

for (et in c('note')){ # , 'domain'
  for (i in c('IPR', 'domain', 'MC29e', 'EC6')){ # , 'MC6', 'note', 'domain'
    for (w in c(F, T)){
      for (r in c(F, T)){
        NPscan(species = 'A.nidu_AspGD', out.file = 'A.nidu.xls', remove.glycan=r,
               init.truth = i, init.expand = 1, init.weighted = w, eval.truth = et,remove.pseudocounts=F, 
               nNPG= 2, nOG = 1, do.viterbi = T, EM.method ='null',pseudocount = 0.1, 
               dom.freq.cutoff = 1)   
      }
    }
  }
}
# for EC6 initialization, remove glycan helps, for IRP/domain/MC29e initialization remove glycan does not help
# best performer using 'note' based SM genes as truth is 
#  --- file perf_itDOMAINie1iwTetNOTEiecPROBidr0.88rcFrpFnpg2og1domfr1pc0.1_viterbi.pdf
# i.e. init.truth = 'domain', remove.glycan = F, init.weighted = F
# the second best is 
# ---- perf_itIPRie1iwTetNOTEiecPROBidr0.88rcFrpFnpg2og1domfr1pc0.1_viterbi.pdf
NPscan(species = 'A.nidu_AspGD', out.file = 'A.nidu.xls', remove.glycan=F,
       init.truth = 'domain', init.expand = 1, init.weighted = F, eval.truth = 'note',remove.pseudocounts=F, 
       nNPG= 2, nOG = 1, do.viterbi = T, EM.method ='null',pseudocount = 0.1, 
       dom.freq.cutoff = 1)

# 20160508
data.root = '/Users/yongli/Dropbox/NPGC/NPGCquery_data'
setwd(data.root)
meta = read.table('NPGCquery_meta.txt',header = T,as.is = T, sep= '\t', row.names = 1)

for (i in setdiff(1:nrow(meta), c(1,20))){
  # s = 'A.nidu_AspGD'
  NPscan(species = rownames(meta)[i],  remove.glycan=T,
         init.truth = 'IPR', init.expand = 1, init.weighted = T, eval.truth = 'IPR',remove.pseudocounts=F, 
         nNPG= 2, nOG = 1, do.viterbi = T, EM.method ='null',pseudocount = 0.1, 
         dom.freq.cutoff = 2, proteinID = meta$proteinID[i])
}
NPscan.postprocess()

### 20160517-28, selected speies
files = c('NPScan_GenePred2_Aspergillus_fumigatus Af293 from AspGDitIPRie1iwTetIPRiecPROBidr0.88rcTrpFnpg2og1domfr2pc0.1_viterbi.csv',
          'NPScan_GenePred2_A.nidu_AspGDitIPRie1iwTetIPRiecPROBidr0.88rcTrpFnpg2og1domfr2pc0.1_viterbi.csv',
          'NPScan_GenePred2_Trichoderma_virensitIPRie1iwTetIPRiecPROBidr0.88rcTrpFnpg2og1domfr2pc0.1_viterbi.csv')
NPscan.postprocess(files,tag = 'WalshL5_20p99enzymeOnly', Walsh.only = T, length.cutoff = 5, length.cutoff.max = 21, p.cutoff = 0.99, remove.TF = T, remove.transporter = T)

root = '/Users/yongli/Dropbox/NPGC/NPGCquery_data'
setwd(root)
RORA.pipeline(root = '/Users/yongli/Dropbox/NPGC/NPGCquery_data', cluster.info.file = 'WalshL5_20p99enzymeOnlycluster_info.xlsx', i.start = 1, use.multiple.model.species = T, skip.existing = T, simplify.model.species = F)
setwd(root)
merge.deepAnno(clusterIDs = c('semiUU174', 'semiUU204', 'semiUU559',
                              'UU1', 'UU10', 'UU29', 'UU48'),
               out.file = 'merged.xlsx') # merge the deepAnno of the selected clusters
  
# system('cp */colored_*.xlsx ./')
# system('cp */pretty_pMap*.xlsx ./')
# system('cp */*blastp.hits ./')
# select.CDS.multiModel(re.cluster.model = 'pretty_pMapc(\\w{2}\\d{4})_(.*)$')

# 20160528
f='Aspfu_A1163_1.filtered_proteins.External_Models.gff3'
f = 'Trichoderma_virens.ASM17099v1.23.gff3'
anno = import.gff(f)
m1 = learn.gff.ID.mapping(unlist.multi(anno@elementMetadata@listData$ID), 
                         parent = unlist.multi(anno@elementMetadata@listData$PARENT), 
                         node.type =  as.character(anno@elementMetadata@listData$type))
z = m1$CDS2gene
m1$CDS2transcript('EHK21906')
m1$CDS2gene('EHK21906')
m1$gene2CDS('gene:TRIVIDRAFT_70851')
m1$gene2transcript('gene:TRIVIDRAFT_70851')
z('EHK21906')

f = 'A_nidulans_FGSC_A4_current_features.gff'
anno = import.gff(f)
m2 = learn.gff.ID.mapping(unlist.multi(anno@elementMetadata@listData$ID), 
                         parent = unlist.multi(anno@elementMetadata@listData$Parent), 
                         node.type =  as.character(anno@elementMetadata@listData$type))
m2$CDS2mRNA('AN0010-P')
m2$CDS2gene('AN0010-P')
m2$gene2CDS('AN0010')
m2$gene2mRNA

### detailed RORA on the selected clusters
root = '/Users/yongli/Dropbox/NPGC/NPGCquery_data'
setwd(root)
RORA.pipeline(root = '/Users/yongli/Dropbox/NPGC/NPGCquery_data', i.all = c(1), 
              cluster.info.file = 'WalshL5_20p99enzymeOnlycluster_info.xlsx',
              use.multiple.model.species = T, RORA.iteration = 2, swiss.db = 'swissprot',
              skip.existing =F, extra.genes = 1, tag = '')
RORA.pipeline(root = '/Users/yongli/Dropbox/NPGC/NPGCquery_data', i.all = 1, 
              cluster.info.file = 'WalshL5_20p99enzymeOnlycluster_info.xlsx',
              use.multiple.model.species = T, RORA.iteration = 3, swiss.db = 'fungiRefSwiss70',
              skip.existing =F, extra.genes = 2, tag = 'fungiSwiss70')
RORA.pipeline(root = '/Users/yongli/Dropbox/NPGC/NPGCquery_data', i.all = c(1, 6, 7, 19, 20, 25, 41), 
                          cluster.info.file = 'WalshL5_20p99enzymeOnlycluster_info.xlsx',
                          use.multiple.model.species = T, RORA.iteration = 3, swiss.db = 'fungiRefSwiss70',
                          skip.existing =F, extra.genes = 2, tag = 'fungiSwiss70')

a = system(paste('ls */pMap*Swiss70.xls'), intern = T)
for (a in system(paste('ls */pMap*Swiss70.xls'), intern = T)){
  select.CDS(paste('./', a, sep=''))
}

### a lot of manual selection to obtain file "UU69seqs.fasta"
# verify, UU69seq
setwd('/Users/yongli/Universe/write/Project_Current/9.O.NPbioinformatics/UU201606')
require(xlsx)
seqs = read.fasta('UU69seqs.fasta')
for (i in 1:nrow(seqs)){
  seq.name = rownames(seqs)[i]
  seq.name = strsplit(seq.name, ':')[[1]]
  tab = read.xlsx2(paste('pretty_pMapc', seq.name[1], 'fungiSwiss70.xlsx', sep=''), sheetIndex = 1)
  i.matched = which(tab$ID == seq.name[2] | gsub('transcript:', '',  tab$ID) == seq.name[2])
  if (length(i.matched)){
    cat('matched', tab$seq[i.matched] == seqs[i,1], ':', seq.name, '\n')
  }else{
    cat('not found', ':', seq.name, '\n')
  }
}

### codon optimze
setwd('/Users/yongli/Universe/write/Project_Current/9.O.NPbioinformatics/UU201606')
codon.optimizer(CDS.list.file='UU69seqs.fasta', N.codon.types.to.change = 6, genetic.table=1, host.species='4932', left.extra='', right.extra='',
                restriction.sites = c(BsaI='GGTCTC', BsaI.rc='GAGACC',
                                      AarI='CACCTGC', AarI.rc = 'GCAGGTG',
                                      polyA8='AAAAAAAA', polyC8 = 'CCCCCCCC',
                                      polyG5='GGGGG', polyT8 = 'TTTTTTTT'), tag = 'BsaIAarIPolymers')

