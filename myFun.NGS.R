cluster.chipseq <- function(files = c('H1-H3K4me3.bigwig', 'iPSC-H1-H3K4me3.bigwig',
                                      'HUES3-H3K4me3.bigwig', 'iPSC-HUES3-H3K4me3.bigwig'),
                            labels = c('H1', 'H1i', 'HUES3', 'HUES3i'), main='H3K4me3', log=F, method = c('pearson', 'spearman'), ...){
  # cluster whole genome signals
  # 20140708
  require(rtracklayer)
  method = match.arg(method)
  
  profs = list();
  for (i in 1:length(files)){
    cat(files[i], '\n')
    profs[[labels[i]]] = import(files[i])
    if (log){
      profs[[labels[i]]]$score = log(profs[[labels[i]]]$score)
    }
    if (method == 'spearman'){
      profs[[labels[i]]]$score = rank(profs[[labels[i]]]$score)
    }
  }
  R = cor.GRanges(profs)
  hr <- hclust(as.dist(1-R), method="complete");
  plot(hr, xlab='', sub='', main=main,...)
  return(list(cor=R,hr=hr))
}

cor.GRanges <- function(Glist){
  nm = names(Glist);
  l = length(nm)
  R = matrix(1, nrow = l, ncol = l, dimnames = list(nm, nm));
  for (i in 1:(l-1)){
    for (j in (i+1):l){
      cat(nm[i], ':', nm[j], '\n')
      R[i,j] <- R[j,i] <- cor.GRange.pair(Glist[[i]], Glist[[j]])
    }
  }
  return(R)
}

dist.GRanges <- function(Glist){
  nm = names(Glist);
  l = length(nm)
  R = matrix(1, nrow = l, ncol = l, dimnames = list(nm, nm));
  for (i in 1:(l-1)){
    for (j in (i+1):l){
      cat(nm[i], ':', nm[j], '\n')
      R[i,j] <- R[j,i] <- cor.GRange.pair(Glist[[i]], Glist[[j]])
    }
  }
  return(R)
}

cor.GRange.pair <- function(G1=profs$H1, G2=profs$H1i){
  # compute the pearson correlation coefficient of two Granges
  mu1 = mean.GRanges(G1)
  mu2 = mean.GRanges(G2)
  l = length.GRanges(G1)
  p = dot.GRanges(G1, G2)
  r = (p - l*mu1*mu2)/sqrt((sum.GRanges(G1,function(x){return(x^2)})-l*mu1^2)*
                             (sum.GRanges(G1,function(x){return(x^2)}) -l*mu2^2))
  return(r)
}


sum.GRanges <- function(G1, fun=identity){
  # sum of transformed scores
  i = regexpr('_', G1@seqinfo@seqnames)<1
  i = G1@seqnames %in% G1@seqinfo@seqnames[i]
  return(sum(width(G1[i])*fun(score(G1[i]))))
}

dot.GRanges <- function(G1, G2, fun=identity){
  # dot product
  i = regexpr('_', G1@seqinfo@seqnames)<1
  i = G1@seqnames %in% G1@seqinfo@seqnames[i]
  G1 = G1[i]
  i = regexpr('_', G2@seqinfo@seqnames)<1
  i = G2@seqnames %in% G2@seqinfo@seqnames[i]
  G2 = G2[i]
  G = intersect(G1,G2)
  G1m = subsetByOverlaps(G1, G)
  G2m = subsetByOverlaps(G2, G)
  return(sum(fun(G1m$score)*fun(G2m$score) * width(G1m)))
}

mean.GRanges <- function(G1){
  # i = regexpr('_', G1@seqnames)<1
  i = regexpr('_', G1@seqinfo@seqnames)<1
  i = G1@seqnames %in% G1@seqinfo@seqnames[i]
  return(sum(width(G1[i])*score(G1[i]))/length.GRanges(G1))
}

length.GRanges <- function(G1){
  i = regexpr('_', G1@seqinfo@seqnames)<1
  return(sum(as.numeric(G1@seqinfo@seqlengths)[i]))
}

geneRanges2allGenes <- function(gff.file, gene.ranges, gene.definition = 'CDS', id.type = 'protein_id', format='gff3', do.unique = T){
  # 20151007, Yong Li
  # anno = read.gff3(gff.file, format=format)
  anno = import.gff(gff.file) # 20160502
  idx.gene = (anno$type==gene.definition)
  anno = anno[idx.gene, ]
  anno = anno[!duplicated(anno@elementMetadata[,id.type]),] # only keep one CDS feature per gene
  anno = sort.intervals(anno)
  gene.ranges = sub('\\.\\d$', '', gene.ranges) # remove version numbers
  i.gene = sort(match(gene.ranges, sub('\\.\\d$', '', anno@elementMetadata[,id.type])))
  if (do.unique){
    all.IDs = unique(anno@elementMetadata[min(i.gene):max(i.gene),id.type])    
  }else{
    all.IDs = anno@elementMetadata[min(i.gene):max(i.gene),id.type]
  }
  return(all.IDs)
}
geneRanges2ntRanges <- function(anno, gene.ranges, extra.nt=0){
  # 20141011, Yong Li
  if (is.vector(gene.ranges))
    gene.ranges = matrix(gene.ranges, nrow=1, ncol=length(gene.ranges))
  i.gene = match(gene.ranges[,1:2], anno$ID)
  seqnames = matrix(as.character(anno@seqnames[i.gene]), nrow(gene.ranges),2)
  locs.start = matrix(cbind(anno@ranges@start[i.gene]), nrow(gene.ranges),2)
  locs.end = matrix(locs.start + anno@ranges@width[i.gene], nrow(gene.ranges),2)
  locs = data.frame(seqnames=seqnames[,1], start=locs.start[,1], end=locs.end[,2])
  locs$start = locs$start - extra.nt
  locs$end = locs$end + extra.nt
  return(locs)
}

getDNA.subseq <- function(DNA.fasta.file, locs = NULL, extra.nt = 5000, gene.ranges=NULL, anno=NULL, gff.file = NULL ){
  # YF Li,
  # 201410
  require('Biostrings')
  if (is.null(locs)){
    if (!is.null(gff.file)){
      # anno = read.gff3(gff.file, format = 'gff3')
      anno = import.gff(gff.file) # 20160502
    }
    locs = geneRanges2ntRanges(anno, gene.ranges, extra.nt)
  }
  if (is.vector(locs))
    locs = matrix(locs, nrow = 1, ncol=length(locs)) # 20141113
  locs = sort.by(locs, by = locs[,2])
  # anno$range[match(gene.ranges[,1:2], anno$ID)]
  fa=import(DNA.fasta.file, 'fasta', type='DNA')
  names(fa) <- sub('^([^ ]+) .+$','\\1', names(fa))
  st = sapply(1:nrow(locs), FUN = function(x){max(locs[x,2],1)});  # 20141125
  L = nchar(fa[locs[,1]]) # 20141125
  en = sapply(1:nrow(locs), FUN = function(x){min(locs[x,3],L[x])}); # 20141125
  out = subseq(fa[locs[,1]], st, en)
  # names(out) = paste(c(paste(gene.ranges[,3], gene.ranges[,1], gene.ranges[,2], sep='_'), paste(locs[,1], '[',locs[,2],',',locs[,3],']',sep='')),collapse = '|')
  names(out) = paste(c(paste(locs[,1], '[',st,',',en,']',sep=''), paste(gene.ranges[,3], gene.ranges[,1], gene.ranges[,2], sep='_')),collapse = ' ')
  return(out)
}

gff.subset <- function(gff.file=NULL, locs=NULL, out.file='gff_sub.gff', format = 'gff2', shift=T){
  # 20141116
  # anno = read.gff3(gff.file, format='gff3') 
  anno = import.gff(gff.file) # 20160502
  anno.sub = subsetByOverlaps(anno, GRanges(seqnames = locs[,1], ranges = IRanges(start = locs[,2], end=locs[,3])),ignore.strand = FALSE)
  R = anno.sub@ranges;
  if (shift){
    anno.sub = GRanges(paste(anno.sub@seqnames, '[',locs[,2], ',',locs[,3],']',sep=''), ranges = IRanges(start = R@start - locs[,2] + 1, width = R@width), strand = anno.sub@strand, anno.sub@elementMetadata)
  }
  export(anno.sub, con = out.file, format = format)  
}

translate.fasta.v1 <- function(CDS.file, pep.file){
  # Yong Fuga Li, 20141216
  # 
  CDS = read.fasta(CDS.file);
  for (i in 1:nrow(CDS)){
    CDS[i,'seq'] = as.character(translate(DNAString(CDS[i,'seq']), if.fuzzy.codon = 'X'))
  }
  write.fasta(CDS, out.file = pep.file)
}

translate.fasta <- function(CDS.file, pep.file){
  # Yong Fuga Li, 20141216
  # v2, 20150916
  # readDNAStringSet
  CDS = readDNAStringSet(CDS.file, format = 'fasta');
  pro = translate(CDS, if.fuzzy.codon = 'X')
  export(pro, format = 'fasta', con = pep.file)
}

bam.extract.shift <- function(bam.file, locs, tag, shift=F){
  # extract RNA-seq reads from bam files and change coordinates
  # Yong Fuga Li, 20141117
  require('GenomicAlignments')
  bam.out.file = paste(tag, 'RNAseq.bam', sep='')
  sam.out.file = paste(tag, 'RNAseq.sam', sep='')  
  system(paste('samtools view -o', bam.out.file, '-b', bam.file, paste(locs$seqnames[1],':', locs$start[1], '-', locs$end[1], sep='')))    
  if (shift){
    system(paste('samtools view -H -o', sam.out.file, bam.out.file, sep=' ')) 
    system(paste('samtools view -o', 'tmp.sam', bam.out.file, sep=' '))   
    sam = read.csv('tmp.sam', header = F, sep = '\t', quote = "", as.is = T, comment.char = '')
    sam = sam[sam$V3!='',]
    sam$V3[sam$V3!=''] = paste(sam$V3[sam$V3!=''], '[',locs[,2], ',',locs[,3],']',sep='')
    sam$V4[sam$V4!=''] = as.numeric(sam$V4[sam$V4!='']) - locs[,2] + 1
    write.table(sam, file = sam.out.file, col.names = F, row.names = F, sep='\t', quote=F, append = T)
    system(paste('samtools view -b -o', bam.out.file, ' ', sam.out.file))  
  }
  system(paste('samtools sort -n ', bam.out.file, 'tmp'))
  system(paste('mv tmp.bam sorted_', bam.out.file, sep=''))
  system(paste('samtools index ', bam.out.file, sep=''))
  return(bam.out.file)
}

extra.chr <- function(DNA.fasta.file, chr,out.file='chr.fasta'){
  # extract chromosome sequence
  # YF Li, 20141118
  require('Biostrings')
  fa=import(DNA.fasta.file, 'fasta', type='DNA')
  names(fa) <- sub('^([^ ]+) .+$','\\1', names(fa))
  export(fa[chr], con = out.file, format = 'fasta')
  return(out.file)
}

gff.unshift <- function(gff.shifted.file, gff.out.file = paste('unshift_', gff.shifted.file)){
  # 20141223 - allow empty file
  gff = tryCatch(read.table(gff.shifted.file, header = F, as.is=T), error = function(e){NULL}, finally = NULL)
  if (is.null(gff)){
    system(paste('cp ', gff.shifted.file,  gff.out.file))
    return(gff.out.file)
  }
  gff = read.table(gff.shifted.file, header = F, as.is=T)
  chr.tmp = do.call(rbind, strsplit(gff[,1], split = '\\[|\\]|\\,'))
  gff[,1] = chr.tmp[,1];
  chr.ranges = as.numeric(chr.tmp[,2:3]);
  gff[,4:5] = gff[,4:5] + as.numeric(chr.tmp[,2]) - 1
  write.table(gff, file = gff.out.file, sep = '\t', col.names = F, row.names = F, quote=F)
  return(gff.out.file)
}


bam.stats <- function(bam.file, length.bins = c(0, 300, 1000, Inf), tag = '', mapq.cut = 0){
  # 20160424, Yong Fuga Li
  bam = scanBam(bam.file)[[1]]
  hist(bam$mapq, main = tag, xlab = 'mapping quality')
  n.reads = length(bam$qname);
  # sapply(cbind, bam)
  # bam <- as.data.frame(scanBam(f)[[1]])
  i.keep = !is.na(bam$rname) & !is.na(bam$pos) & bam$mapq >=mapq.cut
  n.mapped = sum(i.keep)
  bam = data.frame(rname = bam$rname, qwidth = bam$qwidth, isize = bam$isize, mapq = bam$mapq) # as.data.frame(bam)[1:length(bam$qname),]
  bam = bam[i.keep, ]
  hist(bam$mapq, main = tag, xlab = 'mapping quality')
  sum(bam$isize[i.keep] > 0, na.rm = T)
  sum(bam$isize[i.keep] == 0, na.rm = T)
  n = hist(log10(bam$isize[bam$isize>0 & bam$isize < 10000]), breaks = 50, plot = F)
  plot(x=n$breaks[1:(length(n$breaks))], y=c(n$count,0)+1, log='y', type='s', lwd=2, xlab = 'length', ylab = 'counts+1', main = paste('all_', tag))
  n = hist(log10(bam$isize[bam$isize>0 & bam$isize < 10000 & (bam$rname %in% 'chrM')]), breaks = 50, plot = F)
  plot(x=n$breaks[1:(length(n$breaks))], y=c(n$count,0)+1, log='y', type='s', lwd=2, xlab = 'length', ylab = 'counts+1', main = paste('mito_', tag))
  n = hist(log(bam$isize[bam$isize>0]), breaks = 50, plot = F)
  plot(x=n$breaks[1:(length(n$breaks))], y=c(n$count,0)+1, log='y', type='s', lwd=2, xlab = 'length', ylab = 'counts+1', main = paste('all_', tag))
  # plot(n$mids, n$count, log='y', type='h', lwd=10, lend=2, xlab = 'length', ylab = 'counts', main = paste('all_', tag))
  
  # plot(n$mids, n$count, log='xy', xlab = 'length', y = 'counts')
  # hist(bam$isize[bam$isize>0 & bam$isize < 10000], main = tag, xlab = 'length', log='y')
  # hist(log(bam$isize[bam$isize>0]), main = tag, xlab = 'log_length')
  # plot(mydata_hist$count, log="y", type='h', lwd=10, lend=2)
  sum(bam$isize>0 & bam$isize < 10000, na.rm = T)
  n.mapped.paired = sum(!is.na(bam$isize))
  n.mapped.paired.sameChr = sum(bam$isize!=0, na.rm = T)
  
  bam.paired = bam[bam$isize>0 & !is.na(bam$isize), c('rname', 'qwidth','isize')]
  length.range = cut(bam.paired$isize, breaks = length.bins);
  len.dist = table(length.range)
  n.chrM.by.len = tapply(bam.paired$rname %in% 'chrM',INDEX = length.range, FUN = mean)
  names(n.chrM.by.len) = paste('chrM%', names(n.chrM.by.len))
  names(len.dist) = paste('counts', names(len.dist))
  i.chrM = bam$rname %in% 'chrM'
  i.chrM.paired.reads = bam.paired$rname %in% 'chrM'
  stats = c("mean read length" = round(mean(bam.paired$qwidth)),
            "counts all" = n.mapped, 
            "counts paired/all" = n.mapped.paired/n.mapped, 
            "counts paired sameChr/paired" = n.mapped.paired.sameChr/n.mapped.paired, 
            len.dist/sum(len.dist), "chrM% all" = mean(i.chrM),
            "chrM% paired" = mean(bam$rname[!is.na(bam$isize)] %in% 'chrM'),
            "chrM% paired sameChr" = mean(i.chrM.paired.reads),
            n.chrM.by.len);
  
  return(list(stats = stats, bam.paired = bam.paired))
}

cf.dna.analysis <- function(bam.files=dir(bam.root, '.*.bam$'), meta.file,
                            IDs = sub('sample(\\d+)_S\\d+.bam','#\\1', bam.files), length.bins, tag = '', mapq.cut = 20,
                            bam.root = '/Volumes/SED/MacBook_Backup/Peidong/Alignment2'){
  # 20160424-26
  # for Peidong's data
  require('rtracklayer')
  require(IRanges)
  require(GenomicRanges)
  require(Rsamtools)
  require(xlsx)
  setwd(bam.root)
  meta = read.xlsx(meta.file,1)
  IDs <- sub('sample(\\d+)_S\\d+.bam','#\\1', bam.files)
  rownames(meta) = meta$Index
  tags = as.character(meta[IDs, 'Sample_cfDNA'])
  
  all.data = c()
  pdf(paste('cfDNA_length_distribution', tag, '.pdf', sep=''), width = 20, height = 20)
  par(mfrow = c(8,5))
  for (i in 1:length(bam.files)){
    f = bam.files[i]
    cat(f,'\n')
    dat = bam.stats(f, length.bins = length.bins, tag = tags[i], mapq.cut = mapq.cut)
    all.data = rbind(all.data, dat$stats)
  }
  dev.off()
  rownames(all.data) = IDs
  
  # setwd('/Users/yongli/Universe/write/Project_Current/1.CFS/cellFreeDNA')
  save(list = 'all.data', file = paste('Peidong_CFS_mitoFrag',tag,'.RData', sep=''))
  # load('Peidong_CFS_mitoFrag.RData')
  
  ######
  # visualization
  pdf(paste('TotalRead',tag,'.pdf', sep=''), height = 3.5,width = 9.5)
  # barplot(all.data[,2], horiz = T, las=2)
  print(barplot.gg(t(all.data[,2,drop=F]), ylab = 'log #reads')+theme(legend.position="none")) # + coord_flip() 
  print(barplot.gg(log10(t(all.data[,2,drop=F])), ylab = 'log #reads')+theme(legend.position="none")) # + coord_flip() 
  dev.off()
  
  pdf(paste('Read.distribution',tag,'.pdf', sep=''), height = 3.2,width = 9.5)
  print(plot.multi(1:nrow(all.data), y.multi = all.data[,3:4], legend.position = c(0.9,0), xlab = 'sample ID', ylab = 'fraction of counts',by.name = '', lwd = 2)+theme(legend.position="top"))# +coord_flip()
  print(barplot.gg(t(all.data[,4+(1:(length(length.bins)-1))]), position = position_stack())+theme(legend.position="top"))# +coord_flip()
  # plot.multi(1:64, y.multi = all.data[,5:7], xlab = 'sample ID', ylab = 'fraction of counts')+coord_flip()
  dev.off()
  
  pdf(paste('mitoFraction',tag,'.pdf', sep=''), height = 3.2, width = 9.5)
  print(plot.multi(1:nrow(all.data), y.multi = all.data[,(4:6)+length(length.bins)], xlab = 'sample ID', ylab = 'fraction of mitochondria', by.name = '',legend.position = c(0.9,0), lwd = 2)+theme(legend.position="top"))# +coord_flip()
  print(plot.multi(1:nrow(all.data), y.multi = all.data[,(7+length(length.bins)):(5+2*length(length.bins))], xlab = 'sample ID', ylab = 'fraction of mitochondria',by.name = '', legend.position = c(0.9,0.1), lwd = 2)+theme(legend.position="top"))# +coord_flip()
  for (j in (7+length(length.bins)):(5+2*length(length.bins))){
    # print(plot.multi(1:nrow(all.data), y.multi = all.data[,j, drop=F], xlab = 'sample ID', ylab = 'fraction of mitochondria',by.name = '', legend.position = c(0.9,0.1), lwd = 2)+theme(legend.position="top"))# +coord_flip()
    print(barplot.gg(t(all.data[,j, drop=F]), position = position_stack())+theme(legend.position="top"))# +coord_flip()
  }
  dev.off()
  
  if (0){
    barplot.gg(all.data[,3:7])
    barplot.gg(t(all.data[,3:7]))
    heat
    plot(all.data[,8:13]) # % mitochondria
    barplot.gg(all.data[,8:10])
    barplot.gg(all.data[,11:13])
  }
}

