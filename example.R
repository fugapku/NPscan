################
# NPscan pipeline examples
# covering 
#   1) gene cluster prediction, 
#   2) RORA gene structure prediction/improvement,
#   3) codon optimizaiton
# Yong Fuga Li, 20160911
################
# the meta data file is NPGCquery_meta.txt
# it provides the meta data for each species (or more accurately, each genome) for which the binary input data has been prepared
# the binary data files is located in the data folder

### HMM learning and predictions on a genome names "A.nidu_AspGD"
NPscan(species = 'A.nidu_AspGD', # specify the species, should be one of the species specified in file "NPGCquery_meta.txt"
       out.file = 'A.nidu.xls', # output file name
       nNPG= 2, # number of NP gene hidden states (i.e. number of gene cluster types). 2 is used here one captures the glycan containing clusters, another captures the rest
       nOG = 1, # number of other gene hidden states
       remove.glycan=F, # label glycan related domains, and also remove some other non NP domains during HMM learning. 
       dom.freq.cutoff = 1, # ignore low frequency domains during HMM learning.
       init.truth = 'domain', init.expand = 1, init.weighted = F, eval.truth = 'note',remove.pseudocounts=F, 
       do.viterbi = T, EM.method ='null', pseudocount = 0.1)

### Postprocess all NP cluster prediction files with filenames matching '^NPScan_GenePred2_*'
NPscan.postprocess()

### Postprocess NP gene clusters for 3 selected genomes: mainly to compile a set of clusters into a meta data table "*cluster_info.xlsx"
files = c('NPScan_GenePred2_Aspergillus_fumigatus Af293 from AspGDitIPRie1iwTetIPRiecPROBidr0.88rcTrpFnpg2og1domfr2pc0.1_viterbi.csv',
          'NPScan_GenePred2_A.nidu_AspGDitIPRie1iwTetIPRiecPROBidr0.88rcTrpFnpg2og1domfr2pc0.1_viterbi.csv',
          'NPScan_GenePred2_Trichoderma_virensitIPRie1iwTetIPRiecPROBidr0.88rcTrpFnpg2og1domfr2pc0.1_viterbi.csv')
NPscan.postprocess(files,tag = 'WalshL5_20p99enzymeOnly', Walsh.only = T, length.cutoff = 5, length.cutoff.max = 21, p.cutoff = 0.99, remove.TF = T, remove.transporter = T)

### run RORA pipeline to improve the gene prediction for a set of gene clusters specified in a *cluster_info.xlsx
RORA.pipeline(cluster.info.file = 'WalshL5_20p99enzymeOnlycluster_info.xlsx', # specify the cluster information table
              use.multiple.model.species = F, # use more than one model species for gene predictions, if there are multiple species of equal distance to the species being predicted
              RORA.iteration = 1, # number of iterations for improving the gene structure predictions
              extra.genes = 2, # include 2 extra genes on each side of a gene cluster for RORA gene prediction
              root = '/Users/yongli/Dropbox/NPGC/NPGCquery_data', # the data folder
              skip.existing = T)

### merge the deepAnno of a list of clusters specified by the cluster IDs
merge.deepAnno(clusterIDs = c('semiUU174', 'semiUU204', 'semiUU559',
                              'UU1', 'UU10', 'UU29', 'UU48'),
               out.file = 'merged.xlsx') 
  
### automatically learn the transformation rules between two types of identifiers (e.g. transcript ID to gene ID, or protein ID to gene ID)
m2 = learn.gff.ID.mapping(unlist.multi(anno@elementMetadata@listData$ID), 
                         parent = unlist.multi(anno@elementMetadata@listData$Parent), 
                         node.type =  as.character(anno@elementMetadata@listData$type))
m2$CDS2mRNA('AN0010-P')
m2$CDS2gene('AN0010-P')
m2$gene2CDS('AN0010')
m2$gene2mRNA

### another run of RORA on a set of selected clusters, e.g. clusters that is going into synthesis pipeline
root = '/Users/yongli/Dropbox/NPGC/NPGCquery_data'
setwd(root)
RORA.pipeline(cluster.info.file = 'WalshL5_20p99enzymeOnlycluster_info.xlsx', # specify the cluster information table
              use.multiple.model.species = T, # use more than one model species for gene predictions, if there are multiple species of equal distance to the species being predicted
              RORA.iteration = 3, # number of iterations for improving the gene structure predictions
              swiss.db = 'fungiRefSwiss70', # reference protein database to use, either swissprot or fungiRefSwiss70
              extra.genes = 2, # include 2 extra genes on each side of a gene cluster for RORA gene prediction
              root = '/Users/yongli/Dropbox/NPGC/NPGCquery_data', # the data folder
              skip.existing = F, simplify.model.species = F, tag = 'fungiSwiss70')

### some manual work is required to select the optimal set of gene structures based on RORA output

### codon optimze
setwd('/Users/yongli/Universe/write/Project_Current/9.O.NPbioinformatics/UU201606')
codon.optimizer(CDS.list.file='UU69seqs.fasta', # input fasta file
                N.codon.types.to.change = 6, # only optimize the rarest 6 types of codons
                genetic.table=1, 
                host.species='4932', # yeast host
                left.extra='', right.extra='', # adpator sequences to add to 5' and 3' to each CDS
                restriction.sites = c(BsaI='GGTCTC', BsaI.rc='GAGACC', # restriction.sites specifies the sites to remove from the gene (CDS) sequences
                                      AarI='CACCTGC', AarI.rc = 'GCAGGTG',
                                      polyA8='AAAAAAAA', polyC8 = 'CCCCCCCC',
                                      polyG5='GGGGG', polyT8 = 'TTTTTTTT'), 
                tag = 'BsaIAarIPolymers') # tag to add to the output file name

