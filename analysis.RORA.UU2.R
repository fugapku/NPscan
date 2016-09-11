# 20160525
source('/Users/yongli/Dropbox/Galaxy/R/myAnalysis/myFun.R')
source('/Users/yongli/Dropbox/Galaxy/R/myAnalysis/myFun.NP.R')
setwd('/Users/yongli/Universe/write/Project_Current/9.O.NPbioinformatics/AllFungalGenomes/Aspergillus_Binary')
iprscan.tab.file = 'Pt_K85_scafSeq.augoNovo_mANidulans.faa.tsv'
gff.file = 'Pt_K85_scafSeq.augoNovo_mANidulans.gff'
DNA.file = 'Pt_K85_scafSeq.fasta'
pep.fasta.file = 'Pt_K85_scafSeq.augoNovo_mANidulans.faa'

cluster.deepAnno(gene.ranges = c('g5965.t1', 'g5981.t1', 'Batch2C8'),  gff.file = gff.file, DNA.fasta.file = DNA.file, n.cluster.per.file = 2000, extra.genes = 5,
                 gene.definition = 'transcript',
                 proteinID = 'ID', prot.fasta.file = 'Tvirens_v2.FrozenGeneCatalog_20100318.proteins.fasta.gz', 
                 RORA.iteration = 0, species = 'aspergillus_nidulans', plotLogo=F,multialn.method = 'muscle',RORA.topOnly=T,
                 iprscan.tab.file=iprscan.tab.file, out.file='Batch2C8.xlsx',  geneID2cdsID = function(x){paste(x, '.cds', sep='')})
xlsx.color.NPGC('Batch2C8.xlsx')

meta = read.table('/Users/yongli/Dropbox/NPGC/NPGCquery_data/NPGCquery_meta.txt',header = T,as.is = T, sep= '\t', row.names = 1)


source('/Users/yongli/Dropbox/Galaxy/R/myAnalysis/myFun.R')
source('/Users/yongli/Dropbox/Galaxy/R/myAnalysis/myFun.NP.R')
root = '/Users/yongli/Universe/write/Project_Current/9.O.NPbioinformatics/Problematic_PKS'
cluster.info.file = 'Failed_PKS_clusters.xlsx'
setwd(root); RORA.pipeline(root, cluster.info.file, i.start = 1, use.multiple.model.species = T, skip.existing = T, simplify.model.species = F)
setwd(root)
system('cp */colored_*.xlsx ./')
system('cp */pretty_pMap*.xlsx ./')
system('cp */*blastp.hits ./')
select.CDS.multiModel(re.cluster.model = 'pretty_pMapc(\\w{2}\\d{4})_(.*)$')
