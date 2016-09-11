# test start, stop, and irpart
# pMapcAspor1.semiUU48fungiSwiss70
# gene_1216
# gene_1217
# pretty_pMapcAspni_NRRL3_1.semiUU99fungiSwiss70
# gene_2597
# gene_2598
root = '/Users/yongli/Dropbox/NPGC/NPGCquery_data'
setwd(root)
RORA.pipeline(root = '/Users/yongli/Dropbox/NPGC/NPGCquery_data', i.all = 3:4, 
              cluster.info.file = 'WalshL5_20p99enzymeOnly20160802cluster_info_semiUU.xlsx',
              use.multiple.model.species = T, RORA.iteration = 3, swiss.db = 'fungiRefSwiss70',
              skip.existing =F, simplify.model.species = F, extra.genes = 2, tag = 'test')


