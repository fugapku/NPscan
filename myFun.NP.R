# myFun.NP.R
read.antiSMASH <- function(file='fungalAntismashClusters.txt', root=NULL){ #'/Users/yongli/Dropbox/Galaxy/Project_Current/t.NPbioinformatics'){
  
  setwd(root)
  a = read.table(file, header=F, sep='\t', comment='', quote="")
  #   Description  Genbank ID  Cluster num	antismash cluster type	cluster length (bp)	Start	End
  
}

get.runs <- function(seq){
  # 20140429: compute the length of runs for each type of elements in the sequence
  # 20140429: YF Li
  
  is.changed = c(T, seq[2:length(seq)]!= seq[1:(length(seq)-1)], T);
  i = which(is.changed)
  k = diff(i) # run length
  names(k) = seq[i[1:(length(i)-1)]]
  return(k)
}

label.runs <- function(seq=NULL, runs=NULL){
  # Yong Fuga Li, 20140528
  if (is.null(runs)){
    runs = get.runs(seq+1);
    loc.runs = cumsum(runs)
    n = length(seq)    
  }else{
    loc.runs = cumsum(runs)
    # n = sum(runs)
    n = loc.runs[length(loc.runs)]
  }
  labels.runs = zeros(n=n)
  idx.pos = names(loc.runs)=='2'
  labels.runs[loc.runs[idx.pos]] = runs[idx.pos]
  return(labels.runs)
}

label.successes <- function(seq, window.size){
  require('TTR')
  if (window.size>length(seq))
    return(zeros(n=length(seq)))
  a = runSum(seq+0,n=window.size)
  a[is.na(a)] = 0;
  return(a)
}


get.local.max.index <- function(a, tie.resolve=c('first', 'last', 'all')){
  # 20140527-28
  # Yong Fuga Li
  
  tie.resolve = match.arg(tie.resolve)
  b = diff(c(-Inf, a,-Inf))
  s = sign(b)
  idx = which(s!=0)
  idx.max = which(diff(s[idx])==-2)
  #   s[idx[idx.max]]
  #   s[idx[idx.max+1]]
  if (tie.resolve=='first'){
    out  = idx[idx.max]
  }else if (tie.resolve == 'last'){
    out = idx[idx.max+1]-1
  }else if (tie.resolve =='all'){
    out = c()
    for (i in 1:length(idx.max)){
      out = c(out, idx[idx.max[i]]:(idx[idx.max[i]+1]-1))
    }
  }
  return(out)
}

label.successes.local.max <- function(seq, window.size, tie.resolve=c('first', 'last', 'All'), default = 0){
  # get the local max # number successes in sliding windows and set non-local maxes to a default value
  # tie.resolve: resolving ties by taking the first, last, or all
  # 20140527-28
  require('TTR')
  tie.resolve = match.arg(tie.resolve)
  if (window.size>length(seq))
    return(zeros(n=length(seq)))
  a = runSum(seq+0,n=window.size)
  a[is.na(a)] = 0;
  to.keep = get.local.max.index(a, tie.resolve)
  a[setdiff(1:length(a),to.keep)] = default
  return(a)
}

count.runs <- function(runs, max.k = max(runs), types=NULL){
  # count the number of runs of each type of elements
  # 20140429: YF Li
  if (is.null(types))
    types = unique(names(runs))
  C = matrix(0, nrow=length(types), ncol=max.k, dimnames=list(types,1:max.k))
  for (t in types){
    uc = unique.count(runs[names(runs)==t])$counts.unique
    uc = uc[intersect(names(uc), colnames(C))]
    C[t, names(uc)] = uc
  }
  return(C)
}

count.successes <- function(seq, window.size=20, weights=NULL, types=types){
  # counts the elements within sliding windows
  # YF Li, 20140429
  # require('IRanges')
  require('TTR') # 20140527
  warning('weights no implemented yet')
  if (!length(types))
    types = unique(seq)
  C = matrix(0, nrow=length(types), ncol=window.size+1, dimnames=list(types,0:window.size))
  for (t in types){
    # n.success = as.vector(runsum(Rle((seq==t)+0),k=window.size))
    if (window.size>length(seq)){
      n.success = c()    
    }else{
      n.success = runSum((seq==t)+0,n=window.size)[window.size:length(seq)]
    }
    uc = unique.count(n.success)$counts.unique
    C[as.character(t), names(uc)] = as.double(uc)
  }
  return(C)
}

count.successes.local <- function(seq, window.size=20, types){
  # counts the elements within sliding windows
  # YF Li, 20140429
  # require('IRanges')
  require('TTR') # 20140527
  warning('weights no implemented yet')
  if (!length(types))
    types = unique(seq)
  C = matrix(0, nrow=length(types), ncol=window.size+1, dimnames=list(types,0:window.size))
  for (t in types){
    # n.success = as.vector(runsum(Rle((seq==t)+0),k=window.size))
    if (window.size>length(seq)){
      n.success = c()    
    }else{
      n.success = label.successes.local.max((seq==t)+0,window.size)
    }
    uc = unique.count(n.success)$counts.unique
    C[as.character(t), names(uc)] = as.double(uc)
  }
  return(C)  
}

successes.expect <- function(N, n, probs){
  # N: sequence length, n: window size; k: successes; probs: success probability
  # counts the elements within sliding windows
  # YF Li, 20140429
  
  if (length(probs)<2)
    stop('need probability profile, i.e. for more than one elements')
  probs = probs/sum(probs);
  if (any(probs<0))
    stop('Need positive probabilities')
  if (is.null(names(probs)))
    names(probs) = 1:length(probs)
  C = matrix(0, nrow=length(probs), ncol=n+1, dimnames=list(names(probs), 0:n))
  for (t in 1:length(probs)){
    C[t,] = (N-n+1)* dbinom(0:n, size=n, prob=probs[t])
  }
  return(C)   
}

run.expect <- function(L, probs, max.k = L){
  # calculate the expected # of runs of length k for each types of elements
  # L: sequence length
  # probs: probability profile for m elements
  # max.k: max run length to evaluate
  # 20140429: YF Li
  
  if (length(probs)<2)
    stop('need probability profile, i.e. for more than one elements')
  probs = probs/sum(probs);
  if (any(probs<0))
    stop('Need positive probabilities')
  if (is.null(names(probs)))
    names(probs) = 1:length(probs)
  ks = 1:max.k
  C = matrix(0, nrow=length(probs), ncol=length(ks), dimnames=list(names(probs), ks))
  for (k in ks){
    C[,k] = ((L-k>=1)*((L-k-1)*(1-probs)^2 + 2*(1-probs)) + (L==k)) * probs^k
  }
  return(C)
}

plot.fdr <- function(observed, expected, quantile.cutoff = 0.5, reverse=T, do.plot=T, log.scale=F,tag = '', ...){
  # plot FDR curved based on observed distribution and theoretical distribution
  # reverse = T ==> higher score more likely to be true
  # reverse = F ==> lower score more likely to be true
  # Yong Fuga Li
  # 20140428
  # 20140503: quantile.cutoff, quantile of uptail instances in the expected distribution to used for the estimation of FDR
  # 20140527: do.plot
  
  
  quantile.cutoff = min(max(0, quantile.cutoff), 1)
  
  if (is.null(names(observed)))
    names(observed) = 1:length(observed)
  if (is.null(names(expected)))
    names(expected) = 1:length(expected)
  if (reverse){
    observed = rev(observed)
    expected = rev(expected)
  }
  epsilon = 1E-15
  
  n.pos = cumsum(observed - expected)
  fdr = cumsum(expected)/(cumsum(observed)+epsilon)
  quant = cumsum(expected)/sum(expected)
  i.max = min(which(quantile.cutoff<=quant))
  idx = (fdr<=1 & fdr>=0 & observed >= 1); # 20150528: add observed >= 1
  idx[seq2(from=i.max+1, to=length(expected), by=1)] = F
  if (any(idx) & do.plot){
    plot(fdr[idx], n.pos[idx], xlab='False Discovery Rate', ylab='# True gene cluster',...)    
  }
  
  if (do.plot){
    #     dat = rbind(data.frame(score = as.factor(as.numeric(names(observed))), counts=observed, observed='observed'), 
    #                 data.frame(score = as.factor(as.numeric(names(expected))), counts=expected, observed='expected'))
    #     print(barchart(counts~score, data= dat, xlab=tag, groups=observed, 
    #                    equispaced.log=T, scales=list(y = list(log = log.scale)), auto.key=T))  
    dat = rbind(data.frame(score = as.numeric(names(observed)), counts=observed, observed='observed'), 
                data.frame(score = as.numeric(names(expected)), counts=expected, observed='expected'))
    g = ggplot(data=dat) + 
      geom_line(aes(x=score ,y=counts,color=observed))
    print(g)
  }
  return(max(max(c(-Inf,n.pos[idx])), 0))
  
}

distribution.diff <- function(sample=labels.succ.local.all, null.samples=labels.succ.local.all.simus, nbins = NULL, quantile.cutoff = 0.5, reverse=T, do.plot=T, log.scale=F, tag = ''){
  # estimate the total number of true instances in sample with multiple null.samples as reference
  # Yong Fuga Li, 20141220, modified from plot.fdr
  # note: the sample and null.samples can be trancated distribution (e.g. filtered to be postive only), so I do not assume equal sizes of the data
  #       but we do assume the full samples are of the same sizes for all samples
  # input: sample - a vector
  #        null.samples - a list of vectors
  # output: 1) total trues in samples; 2) null distribution of total trues and p-values associated with it. 
  #         3) a plot of the sample distribution against null; 4) a plot of the null distribution of total trues
  quantile.cutoff = min(max(0, quantile.cutoff), 1)
  size.sample = length(sample)
  size.total.null = sum(sapply(null.samples, FUN = length)); size.total = size.total.null + size.sample
  n.sample = length(null.samples)
  if (is.null(nbins))
    nbins = round(sqrt(size.sample))
  R = range(c(unlist(null.samples), sample))
  # R = range(sample)
  breaks = seq(from = R[1], to = R[2], by = (R[2]-R[1])/nbins)
  rep.value = round((breaks[2:(nbins+1)] + breaks[1:nbins])/2,4)
  breaks[1] = breaks[1] - (R[2]-R[1])/nbins * 0.01; breaks[nbins+1] = breaks[nbins+1] + (R[2]-R[1])/nbins * 0.01; 
  get.count <- function(x, breaks){
    observed = unique.count(rep.value[cut(x, breaks = breaks)])$counts.unique
    observed = sort.by(observed, as.numeric(names(observed)))  
    observed = mat.fill.row(observed, rep.value)
    return(observed)
  }
  observed = get.count(sample, breaks)
  expected.all = lapply(null.samples, FUN = function(x){get.count(x, breaks = breaks)})
  expected.merged = get.count(unlist(null.samples), breaks)
  n.pos = plot.fdr(observed, expected.merged/n.sample, reverse=T, main='FDR curve', tag=tag)    
  #     ggplot(data=rbind(data.frame(x=as.numeric(names(observed)), y=observed, data='real genome'), 
  #                       data.frame(x=as.numeric(names(expected.merged)), y=expected.merged, data='null'))) + 
  #       geom_line(aes(x=x,y=y,color=data))
  n.pos.null = vector(mode = 'numeric', length = n.sample)
  if (n.sample>1){
    for (i in 1:n.sample){
      n.pos.null[i] = plot.fdr(expected.all[[i]], (expected.merged-expected.all[[i]])/(n.sample-1+1E-10), reverse=T, do.plot = F)
    }      
  }
  d = hist(n.pos.null, plot=F); plot(runMean(d$breaks,2)[2:length(d$breaks)], d$counts, type = 'l',
                                     xlim = c(min(c(d$breaks), n.pos), max(c(n.pos.null, n.pos))), xlab = paste('#true clusters (bin average):', n.pos), ylab='freq'); abline(v=n.pos,  lty = 2)
  
  ################## value -> p-value and value -> fdr
  score2p.value <- function(x){
    x = x/sum(x); x = rev(cumsum(rev(x)))
    return(x)
  }
  
  x = c(breaks[length(breaks)]+10, sort(unlist(null.samples), decreasing = T), breaks[1], breaks[1]-10)
  pvalue = c(0, (0:size.total.null)/size.total.null,1)
  score2pvalue = approxfun(x, pvalue, method='linear')
  x =  c(breaks[length(breaks)]+10, sort(sample, decreasing = T), breaks[1]-10)
  fdr = score2pvalue(x) * size.total.null/n.sample/c(0, 1:size.sample, size.sample); fdr[1] = 0
  fdr[fdr>1] = 1; fdr = cummax(fdr);
  nTruths = c(0, 1:size.sample, size.sample) - score2pvalue(x) * size.total.null/n.sample 
  nTruths[nTruths<0] = 0; nTruths = cummax(nTruths)
  score2fdr = approxfun(x, fdr, method='linear')
  score2ntrue = approxfun(x, nTruths, method='linear')
  # plot(sample,score2pvalue(sample))
  # plot(sample,score2fdr(sample))
  plot(score2fdr(sample),score2ntrue(sample), xlab='q-value', ylab='#true clusters (monotuned)')
  return(list(n.pos = n.pos, p.value = mean(n.pos.null>=n.pos), score2pvalue=score2pvalue, score2fdr=score2fdr, score2ntrue=score2ntrue))
}

NPGC.clustering <- enzyme.clustering <- function(gff.file, iprscan.tab.file = NULL, chromosome.specific=F, 
                                                 gene.definition = c('gene', 'transcript', 'mRNA'), proteinID = 'ID', 
                                                 annotation.by = c('OR', 'desc', 'domain'), 
                                                 tag = 'A_nidulans_FGSC_A4', window.size = 20, log.scale = F, 
                                                 simu.rep = 5, enzyme.definition = c('ase', 'EC6', 'MC29', 'MC29e'),
                                                 prediction.file='Top.Clusters', min.contig.len=4,
                                                 compare.against =c('simulation','theoretical'),
                                                 p.value.cutoff = 0.005,
                                                 outformat=c('csv', 'tab')){
  # statistical analysis of the enzyme runs in a genome
  # chromosome.specific: estimate chromosome specific enzyme probability estimation
  # simu.rep: simulated gene sequences
  # compare.against: using theoretical model or simulation to estimation null distribution, 20140527
  # Yong Fuga Li, 20140428-29
  # 20141124-25: allow the use of domain annotation instead
  
  # enzyme.definition = match.arg(enzyme.definition)
  compare.against = match.arg(compare.against)
  gene.definition = match.arg(gene.definition)  # 20141125
  outformat = match.arg(outformat)
  annotation.by = match.arg(annotation.by)  # 20141125
  require('rtracklayer')
  require('genomeIntervals')
  require(lattice)
  
  # anno = import(gff.file, format='gff')
  gff.format = sub('^.*\\.([^\\.]*$)', '\\1', gff.file)
  # anno = read.gff3(gff.file, format=gff.format)
  anno = import.gff(gff.file) # 20160502
  # anno.chr = anno[anno$type=='chromosome',]
  # chrs = anno.chr@seqnames
  chrs = as.character(unique(anno@seqnames))
  
  ## keep genes only
  # idx.gene = (anno$type=='gene')
  idx.gene = (anno$type==gene.definition) # 20141125
  anno = anno[idx.gene, ]
  anno = sort.intervals(anno)
  colnames(anno@elementMetadata) = toupper(colnames(anno@elementMetadata)) # 20141125
  if (!is.null(anno$NOTE)){ # 20141125
    desc.fname = 'NOTE'
  }else if (!is.null(anno$DESCRIPTION)){
    desc.fname = 'DESCRIPTION'  
  }else{
    warning('No description or Note field for the annotation of genes')      
    desc.fname = 'NOTE'
    anno$NOTE = ''
  }
  
  # read ipr anno: 20141125
  ipr.anno = iprscan.flat(iprscan.tab.file, na.strings = c('-', 'NA', 'NULL'))
  ipr.anno = mat.fill.row(t(t(ipr.anno)), row.names = anno@elementMetadata[,toupper(proteinID)], default = '')[,1]
  names(ipr.anno) = anno$ID
  if (annotation.by %in% 'desc'){
    annotation.text = as.character(as.vector(anno@elementMetadata[[toupper(desc.fname)]]))    
  }else if(annotation.by %in% 'domain'){
    annotation.text = as.character(as.vector(ipr.anno));
  }else if(annotation.by %in% c('OR')){
    annotation.text = paste(as.character(as.vector(anno@elementMetadata[[toupper(desc.fname)]])), as.character(as.vector(ipr.anno)))    
  }
  
  # is.enzyme.ase = regexpr(pattern='ase[ $]', text = annotation.text, perl=T)>0
  is.enzyme.ase = regexpr(pattern='(?: |^)[^ ]+ase(?: |$)', text = annotation.text, perl=T)>0 # 20140519
  is.enzyme.EC6 = regexpr(pattern='(oxidoreductase|transferase|hydrolase|lyase|isomerase|ligase)', text = annotation.text, perl=T, ignore.case=T) > 0 
  is.enzyme.MC29 = regexpr(pattern='(oxidoreductase|hydrolase|dehydrogenase|synthase|reductase|transferase|methyltransferase|oxidase|synthetase|monooxygenase|isomerase|dehydratase|decarboxylase|deaminase|O\\-methyltransferase|transaminase|hydratase|acetyltransferase|N\\-acetyltransferase|dioxygenase|aminotransferase|O\\-acyltransferase|esterase|N\\-methyltransferase|acyltransferase|aldolase|thiolesterase|O\\-acetyltransferase|cyclase)', text = annotation.text, perl=T, ignore.case=T) > 0 
  is.enzyme.MC29e = regexpr(pattern='(oxidoreductase|hydrolase|dehydrogenase|synthase|reductase|transferase|methyltransferase|oxidase|synthetase|monooxygenase|isomerase|dehydratase|decarboxylase|deaminase|O\\-methyltransferase|transaminase|hydratase|acetyltransferase|N\\-acetyltransferase|dioxygenase|aminotransferase|O\\-acyltransferase|esterase|N\\-methyltransferase|acyltransferase|aldolase|O\\-acetyltransferase|cyclase|catalase|hydroxylase|P450|transporter|transcription factor)', text = annotation.text, perl=T, ignore.case=T) > 0 
  cat('# enzymes by ase:', sum(is.enzyme.ase))
  cat('# enzymes by EC 6 class:', sum(is.enzyme.EC6))
  cat('Some none EC6 enzymes', as.vector(annotation.text[is.enzyme.ase & !is.enzyme.EC6])[1:10])
  
  if (sum(is.enzyme.ase)==0 && sum(is.enzyme.EC6)==0){
    warning('No enzyme annotated in the gff file\n')
    return(NULL)
  }
  
  if (enzyme.definition =='ase'){
    is.enzyme = is.enzyme.ase;    
  }else if (enzyme.definition =='EC6'){
    is.enzyme = is.enzyme.EC6;    
  }else if (enzyme.definition =='MC29'){
    is.enzyme = is.enzyme.MC29
  }else if (enzyme.definition == 'MC29e'){
    is.enzyme = is.enzyme.MC29e
  }else{ # 20141125
    is.enzyme = regexpr(pattern=paste('(', enzyme.definition, ')',sep=''), text = annotation.text, perl=T, ignore.case=T) > 0 
    cat('# enzymes:', sum(is.enzyme.EC6))
  }
  # p.enzyme = sum(is.enzyme)/length(anno)
  
  C.run = list()
  C.run.all = c();
  C.run.exp = list()
  C.run.exp.all = c()
  C.run.simu = list()
  C.run.simu.all = c()
  
  C.success = list()
  C.success.all = c();
  C.success.exp = list()
  C.success.exp.all = c()
  C.success.simu = list()
  C.success.simu.all = c()
  
  C.success.local = list()
  C.success.local.all = c();
  C.success.local.simu = list()
  C.success.local.simu.all = c()
  
  types = unique(is.enzyme+1)
  
  L.gene = list()
  p.enzyme = list()
  
  for (i in 1:length(chrs)){
    chr = as.character(chrs[i])
    # cat('processing', chr,'\n')
    #     if (chr == '1099437636266_N_fischeri_NRRL_181'){
    #       1
    #       1
    #     }
    is.in.chr = as.vector(anno@seqnames==chr)
    L.gene[[chr]] = sum(is.in.chr)# number of genes in this chromosome
    seq = is.enzyme[is.in.chr] 
    if (L.gene[[chr]] < min.contig.len)
      next
    if (chromosome.specific){
      p.enzyme[[chr]] = sum(seq)/L.gene[[chr]];      
    }else{
      p.enzyme[[chr]] = sum(is.enzyme)/length(anno)
    }  
    
    runs = get.runs(seq+1);
    labels.runs = label.runs(runs=runs)
    
    labels.succ = label.successes(seq,window.size)
    labels.succ.local = label.successes.local.max(seq,window.size)
    
    C.run[[chr]] = count.runs(runs, types=types)
    C.run.all = sum.union(C.run.all, C.run[[chr]])
    
    C.success[[chr]] = count.successes(seq+1,window.size=window.size, types=types)
    C.success.all = sum.union(C.success.all, C.success[[chr]])
    
    C.success.local[[chr]] = count.successes.local(seq+1,window.size=window.size, types=types)
    C.success.local.all = sum.union(C.success.local.all, C.success.local[[chr]])  
  }
  max.k = ncol(C.run.all)
  # anno[anno@seqnames==chr,]
  
  for (i in 1:length(chrs)){ # recalculate C.run.exp using the global max.k
    chr = as.character(chrs[i])
    
    if (L.gene[[chr]] < min.contig.len)
      next
    C.run.exp[[chr]] = run.expect(L.gene[[chr]], c(1-p.enzyme[[chr]], p.enzyme[[chr]]), max.k=max.k)
    C.run.exp.all = sum.union(C.run.exp.all, C.run.exp[[chr]]) 
    C.success.exp[[chr]] = successes.expect(L.gene[[chr]], n=window.size, probs=c(1-p.enzyme[[chr]], p.enzyme[[chr]]))
    C.success.exp.all = sum.union(C.success.exp.all, C.success.exp[[chr]])
    
    C.run.simu[[chr]] = c() # one chr in all simulations
    C.success.simu[[chr]] = c() # one chr in all simulations
    C.success.local.simu[[chr]] = c() # one chr in all simulations
  }
  
  npos.success.simus <- npos.success.local.simus <- npos.run.simus <- vector('numeric',simu.rep) # number of positive cluster estimated in each of the simulated samples
  C.run.simu.all1s = list(); # record all simulations
  C.success.simu.all1s = list(); # record all simulations
  C.success.local.simu.all1s = list(); # record all simulations
  for (r in 1:simu.rep){
    cat('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b')    
    cat('iteration', r, '\n')    
    # all chr in one simulation: initialization
    C.run.simu.all1 = c();  
    C.success.simu.all1 = c();
    C.success.local.simu.all1 = c();
    
    for (i in 1:length(chrs)){ # recalculate C.run.exp using the global max.k
      chr = as.character(chrs[i])
      if (L.gene[[chr]] < min.contig.len)
        next
      seq.simu = rbinom(L.gene[[chr]], size=1, prob=p.enzyme[[chr]])
      
      # sumamry of one chr in one simulation
      C.run.simu1 = count.runs(get.runs(seq.simu+1), types=types,max.k=max.k)
      C.success.simu1 = count.successes(seq.simu+1, window.size=window.size, types=types)
      C.success.local.simu1 = count.successes.local(seq.simu+1, window.size=window.size, types=types)
      
      # one chr in all simulations
      C.run.simu[[chr]] = sum.union(C.run.simu[[chr]], C.run.simu1/simu.rep) 
      C.success.simu[[chr]] = sum.union(C.success.simu[[chr]], C.success.simu1/simu.rep)
      C.success.local.simu[[chr]] = sum.union(C.success.local.simu[[chr]], C.success.local.simu1/simu.rep)
      
      # all chr in one simulations
      C.run.simu.all1 = sum.union(C.run.simu.all1, C.run.simu1)
      C.success.simu.all1 = sum.union(C.success.simu.all1, C.success.simu1)
      C.success.local.simu.all1 = sum.union(C.success.local.simu.all1, C.success.local.simu1)
      
      # all chr in all simulations
      C.run.simu.all = sum.union(C.run.simu.all, C.run.simu1/simu.rep) 
      C.success.simu.all = sum.union(C.success.simu.all, C.success.simu1/simu.rep)      
      C.success.local.simu.all = sum.union(C.success.local.simu.all, C.success.local.simu1/simu.rep)      
    }
    
    C.run.simu.all1s[[r]] = C.run.simu.all1; # record all simulations
    C.success.simu.all1s[[r]] = C.success.simu.all1; # record all simulations
    C.success.local.simu.all1s[[r]] = C.success.local.simu.all1; # record all simulations
    
  }
  
  # obtain #pos estimation for simulated data
  for (r in 1:simu.rep){
    
    C.run.simu.all1 = C.run.simu.all1s[[r]];
    C.success.simu.all1 = C.success.simu.all1s[[r]];
    C.success.local.simu.all1 = C.success.local.simu.all1s[[r]];
    
    # number of estimated pos in each simulated sample
    if (compare.against=='simulation'){
      npos.run.simus[r] = plot.fdr(C.run.simu.all1[1,], C.run.simu.all[1,], reverse=T, do.plot=F);
      npos.success.simus[r] = plot.fdr(C.success.simu.all1[1,], C.success.simu.all[1,], reverse=T, do.plot=F);      
      npos.success.local.simus[r] = plot.fdr(C.success.local.simu.all1[1,], C.success.local.simu.all[1,], reverse=T, do.plot=F);            
    }else if (compare.against=='theoretical'){
      npos.run.simus[r] = plot.fdr(C.run.simu.all1[1,], C.run.exp.all[1,], reverse=T, do.plot=F);
      npos.success.simus[r] = plot.fdr(C.success.simu.all1[1,], C.success.exp.all[1,], reverse=T, do.plot=F);            
      npos.success.local.simus[r] = plot.fdr(C.success.local.simu.all1[1,], C.success.local.exp.all[1,], reverse=T, do.plot=F);            
    }else{
      stop('compare.against unknown')
    }
  }
  #   round(C.run.exp.all,3) # simulated and theoretical are very close when simu.rep is large, e.g. 2000
  #   C.run.simu.all
  
  pdf(paste('theoretical.distribution.', tag, '.pdf', sep=''),5,4)
  mapping = list(enzyme = '2', 'non-enzyme' = '1')
  for (n in names(mapping)){
    m = mapping[[n]]
    
    t = plot.fdr(C.run.simu.all[m,], C.run.exp.all[m,],reverse=T, main='FDR curve (k-runs): simulated', tag=paste('#', n, 'runs: simulated', simu.rep));
    if (n == 'enzyme')
      npos.run.simu = t
    t = plot.fdr(C.success.simu.all[m,], C.success.exp.all[m,],reverse=T, main='FDR curve (sliding window): simulated', tag=paste('#', n, 'in', window.size, 'genes: simulated', simu.rep));
    if (n == 'enzyme')
      npos.success.simu = t
  }
  dev.off()
  
  ## simulated gene sequences
  pdf(paste('run.stats.', tag, '.pdf', sep=''),5,4)
  for (n in names(mapping)){
    m = mapping[[n]]
    if (compare.against=='theoretical'){
      t = plot.fdr(C.run.all[m,], C.run.exp.all[m,],reverse=T, main='FDR curve (k-runs) against theoretical: all chr', tag=paste('#', n, 'runs: All'))      
    }else if(compare.against=='simulation'){
      t = plot.fdr(C.run.all[m,], C.run.simu.all[m,],reverse=T, main='FDR curve (k-runs) against simulation: all chr', tag=paste('#', n, 'runs: All'))
    }else{
      stop('compare.against unknown')
    }
    if (n == 'enzyme')
      npos.run = t
    # by chromosome plot
    for (i in 1:length(chrs)){ 
      chr = as.character(chrs[i])
      if (L.gene[[chr]] < min.contig.len)
        next
      plot.fdr(C.run[[chr]][m,], C.run.exp[[chr]][m,],reverse=T, main=paste('FDR curve (k-runs) against theoretical:', chr), tag=paste('#', n, 'runs:', chr))
    }
    
    #     dat = rbind(data.frame(run.length = as.factor(as.numeric(colnames(C.run.all))), counts=C.run.all[m,], observed='observed'), 
    #                 data.frame(run.length = as.factor(as.numeric(colnames(C.run.exp.all))), counts=C.run.exp.all[m,], observed='expected'))
    #     print(barchart(counts~run.length, data= dat, xlab=paste(n, 'run length'), groups=observed, 
    # equispaced.log=T, scales=list(y = list(log = log.scale)), auto.key=T))  
    #   plot(C.exp.all[1,], C.all[1,], ylab='observed', xlab='expected', main='non-enzymes')
    #   abline(0,1, lty='dashed', lwd=2)
  }
  dev.off()
  
  ### binomial model
  pdf(paste('window.stats.', tag, '.pdf', sep=''), 5,4)
  for (n in names(mapping)){
    m = mapping[[n]]
    if (compare.against=='theoretical'){
      t = plot.fdr(C.success.all[m,], C.success.exp.all[m,],reverse=T, main='FDR curve (sliding window) against theoretical: all chr', tag=paste('#', n, 'in', window.size, 'genes: All'))
    }else if (compare.against=='simulation'){
      t = plot.fdr(C.success.all[m,], C.success.simu.all[m,],reverse=T, main='FDR curve (sliding window) against simulation: all chr', tag=paste('#', n, 'in', window.size, 'genes: All'))
    }else{
      stop('compare.against unknown')
    }
    if (n == 'enzyme')
      npos.window = t
    # by chromosome plot
    for (i in 1:length(chrs)){ 
      chr = as.character(chrs[i])
      if (L.gene[[chr]] < min.contig.len)
        next
      plot.fdr(C.success[[chr]][m,], C.success.exp[[chr]][m,],reverse=T, main=paste('FDR curve (sliding window):', chr), tag=paste('#', n, 'in', window.size, 'genes:', chr))
    }
  }
  dev.off()
  
  
  ### binomial model
  pdf(paste('local.window.stats.', tag, '.pdf', sep=''), 5,4)
  for (n in names(mapping)){
    m = mapping[[n]]
    if (compare.against=='theoretical'){
      t = plot.fdr(C.success.local.all[m,], C.success.local.exp.all[m,],reverse=T, main='FDR curve (sliding window) against theoretical: all chr', tag=paste('#', n, 'in', window.size, 'genes: All'))
    }else if (compare.against=='simulation'){
      t = plot.fdr(C.success.local.all[m,], C.success.local.simu.all[m,],reverse=T, main='FDR curve (sliding window) against simulation: all chr', tag=paste('#', n, 'in', window.size, 'genes: All'))
    }else{
      stop('compare.against unknown')
    }
    if (n == 'enzyme')
      npos.window.local = t
  }
  dev.off()
  
  
  #### random distribution of npos
  pdf(paste('npos.null.distribution.', tag, '.pdf', sep=''), 5,4)
  p.runCluster = mean(npos.run.simus>npos.run); 
  hist(npos.run.simus, xlab='#Cluster in simulation', main=paste('# clusters by run length:', round(npos.run,2) , 'p-value:', p.runCluster));abline(v=npos.run, lty=5, col='black')
  p.windowCluster = mean(npos.success.simus>npos.window); 
  hist(npos.success.simus, xlab='#Cluster in simulation', main=paste('# clusters by successes:', round(npos.window,2), 'p-value:', p.windowCluster));abline(v=npos.window, lty=5, col='black')
  p.window.localCluster = mean(npos.success.local.simus>npos.window.local); 
  hist(npos.success.local.simus, xlab='#Cluster in simulation', main=paste('# clusters by successes:', round(npos.window.local,2), 'p-value:', p.window.localCluster));abline(v=npos.window.local, lty=5, col='black')
  dev.off()
  
  
  
  ################ output top predictions
  anno.df = as.data.frame(anno,stringsAsFactors=F)
  for (i in 1:length(anno.df)){
    if (class(anno.df[[i]])!='integer')
      anno.df[[i]] = unlist2(anno.df[[i]])
  }
  
  c2p <- function(x){
    x = x/sum(x); x = rev(cumsum(rev(x)))
    return(x)
  }
  
  # padding missing counts of zeros
  C.run.simu.all = cbind('0' = length(anno) - rowSums(C.run.simu.all), C.run.simu.all)
  C.success.simu.all[,'0'] = length(anno) - rowSums(C.success.simu.all) + C.success.simu.all[,'0']
  p.run = c2p(C.run.simu.all[2,])
  p.succ = c2p(C.success.simu.all[2,])
  p.succ.local = c2p(C.success.local.simu.all[2,])
  
  anno.df[, 'run_len'] = 0
  anno.df[, paste('succ_', window.size, sep='')] = 0
  anno.df[, paste('succ_local', window.size, sep='')] = 0
  anno.df[, 'p.value(run_len)'] = 1
  anno.df[, paste('p.value(succ_', window.size, ')',sep='')] = 1
  anno.df[, paste('p.value(succ_local', window.size, ')',sep='')] = 1
  
  for (i in 1:length(chrs)){
    chr = as.character(chrs[i])
    is.in.chr = as.vector(anno@seqnames==chr)
    L.gene[[chr]] = sum(is.in.chr)# number of genes in this chromosome
    if (L.gene[[chr]] < min.contig.len)
      next
    seq = is.enzyme[is.in.chr] 
    if (chromosome.specific){
      p.enzyme[[chr]] = sum(seq)/L.gene[[chr]];      
    }else{
      p.enzyme[[chr]] = sum(is.enzyme)/length(anno)
    }  
    labels.runs = label.runs(seq=seq)    
    labels.succ = label.successes(seq,window.size)
    labels.succ.local = label.successes.local.max(seq,window.size)
    
    # mark the peaks
    cat(chr)
    anno.df[is.in.chr, 'run_len'] = labels.runs
    anno.df[is.in.chr, paste('succ_', window.size, sep='')] = labels.succ
    anno.df[is.in.chr, paste('succ_local', window.size, sep='')] = labels.succ.local
    anno.df[is.in.chr, 'p.value(run_len)'] = p.run[as.character(labels.runs)]
    anno.df[is.in.chr, paste('p.value(succ_', window.size, ')',sep='')] = p.succ[as.character(labels.succ)]
    anno.df[is.in.chr, paste('p.value(succ_local', window.size, ')',sep='')] = p.succ.local[as.character(labels.succ.local)]
  }
  
  # mark the whole clusters
  run.count = 0;
  anno.df[, 'run_clusters'] = ''
  anno.df[, 'succ_clusters'] = ''
  for (i in which(anno.df[, 'run_len']>0)){
    # cat(i)
    run.count = run.count + 1;
    l = anno.df[i, 'run_len'];
    anno.df[(i-l+1):i, 'run_len'] = rowMax(cbind(anno.df[(i-l+1):i, 'run_len'], anno.df[i, 'run_len']))
    anno.df[(i-l+1):i, 'p.value(run_len)'] = rowMin(cbind(anno.df[(i-l+1):i, 'p.value(run_len)'], anno.df[i, 'p.value(run_len)']))
    anno.df[(i-l+1):i, 'run_clusters'] = paste(anno.df[(i-l+1):i, 'run_clusters'], paste('R', run.count,sep=''))
  }
  
  sl = paste('succ_local', window.size, sep='')
  slp = paste('p.value(succ_local', window.size, ')',sep='')
  l = window.size;
  succ.loc.count = 0;
  for (i in which(anno.df[, sl]>0)){
    succ.loc.count = succ.loc.count+1;
    anno.df[(i-l+1):i, sl] = rowMax(cbind(anno.df[(i-l+1):i, sl], anno.df[i, sl]))
    anno.df[(i-l+1):i, slp] = rowMin(cbind(anno.df[(i-l+1):i, slp], anno.df[i, slp]))
    anno.df[(i-l+1):i, 'succ_clusters'] =  paste(anno.df[(i-l+1):i, 'succ_clusters'], paste('S', succ.loc.count,sep=''))
  }
  
  # select top window and run clusters
  to.output.windows = anno.df[,paste('p.value(succ_local', window.size, ')',sep='')] < p.value.cutoff; 
  to.output.runs = anno.df[,'p.value(run_len)'] < p.value.cutoff;
  
  # how many top clusters are included?
  s.names =  anno.df[to.output.windows, 'succ_clusters']
  s.names = strsplit(paste(s.names,collapse=' '), '\\s+',perl=T)[[1]];
  uc = unique.count(s.names)
  n.clusters.localwindows = sum(uc$counts.unique==window.size)
  
  r.names = anno.df[to.output.runs, 'run_clusters']
  n.clusters.runs = length(unique(r.names))
  
  out.names = c(intersect(c('seqnames', 'start', 'end', 'ID', 'Note', 'orf_classification', 'Gene'),colnames(anno.df)),
                colnames(anno.df)[ncol(anno.df)-8+c(7,1,4,8,3,6)])
  if (outformat=='csv'){
    write.table(anno.df[,out.names], file=paste('cluster.anno.full.', tag, '.csv',sep=''),sep=',', row.names=F)
    write.table(anno.df[to.output.windows,out.names], file=paste('cluster.anno.', tag, '.p', p.value.cutoff, '.NWindowClusters',n.clusters.localwindows, '.csv',sep=''),sep=',', row.names=F)
    write.table(anno.df[to.output.runs,out.names], file=paste('cluster.anno.', tag, '.p', p.value.cutoff, '.NRunClusters',n.clusters.runs, '.csv',sep=''),sep=',', row.names=F)    
  }else if (outformat=='tab'){
    write.table(anno.df[,out.names], file=paste('cluster.anno.full.', tag, '.tab',sep=''),sep='\t', row.names=F, quote = F)
    write.table(anno.df[to.output.windows,out.names], file=paste('cluster.anno.', tag, '.p', p.value.cutoff, '.NWindowClusters',n.clusters.localwindows, '.tab',sep=''),sep='\t', row.names=F, quote = F)
    write.table(anno.df[to.output.runs,out.names], file=paste('cluster.anno.', tag, '.p', p.value.cutoff, '.NRunClusters',n.clusters.runs, '.tab',sep=''),sep='\t', row.names=F, quote = F)        
  }
  # write clean per cluster output, 20140611
  write.NPGC <- function(anno.df, i.new.NPG = to.output.windows, window.size=window.size, method=c('WindowLocal', 'Run'),
                         file.out=paste('cluster.anno.clean', tag, '.p', p.value.cutoff, '.NWindowClusters',n.clusters.localwindows, '.tab',sep='')){
    # 20140613
    is.SM = regexpr(pattern='secondary metab', text = as.character(as.vector(anno.df$Note)), perl=T, ignore.case=T)>0
    is.PKS = regexpr(pattern='polyketide synthase', text = as.character(as.vector(anno.df$Note)), perl=T, ignore.case=T)>0
    
    if (method=='WindowLocal'){
      all.SID = anno.df$succ_clusters[i.new.NPG]
      all.SID = strsplit(paste(all.SID,collapse=' '), '\\s+',perl=T)[[1]];
      uc = unique.count(all.SID)
      cluster.names = names(uc$counts.unique[uc$counts.unique==window.size])      
    }else if (method=='Run'){
      r.names = anno.df[i.new.NPG, 'run_clusters']
      cluster.names = unique(r.names) 
    }
    
    clean.table = matrix('',nrow=length(cluster.names),ncol=8,
                         dimnames=list(cluster.names, c('cluster ID', 'chr', 'coordinate', 'gene range', 'min distance to SM genes', 'closest SM gene(s)', 'p-value', 'cluster gene annotations')));
    n.correct.cluster = 0;
    for (nc in cluster.names){
      if (method=='WindowLocal'){
        i.match = regexpr(paste(nc,'(\\s|$)',sep=''), anno.df$succ_clusters)>0
      }else{
        i.match = regexpr(paste(nc,'(\\s|$)',sep=''), anno.df$run_clusters)>0
      }# mapped$cluster.ID[i.match] = nc
      ## get closest SM
      chr = unique(anno.df$seqnames[i.match])
      loc.SM = t(which(is.SM & anno.df$seqnames==chr))
      loc.cluster = t(t(which(i.match)))
      dist.to.SM = repmat(loc.cluster,1,length(loc.SM)) - repmat(loc.SM, length(loc.cluster),1) 
      min.dist.to.SM = min(abs(dist.to.SM))
      #if (min.dist.to.SM)
      if (!min.dist.to.SM) # 20140720
        n.correct.cluster = n.correct.cluster + 1
      closest.SM = which(abs(dist.to.SM)==min.dist.to.SM,arr.ind=T)
      if (!is.null(closest.SM) && length(closest.SM)>0){
        min.dist.to.SM = paste(dist.to.SM[closest.SM], collapse='...')
        closest.SM = loc.SM[closest.SM[,2]]        
      }
      
      # cluster coordinates
      min.l = min(c(anno.df$start[i.match], anno.df$end[i.match]))
      max.l = max(c(anno.df$start[i.match], anno.df$end[i.match]))
      
      # cluster gene ranges
      first.gene = anno.df$ID[min(which(i.match))]
      last.gene = anno.df$ID[max(which(i.match))]
      
      # cluster all gene annotations;
      cluster.anno = paste(anno.df$ID[i.match], anno.df$Note[i.match], sep='|', collapse='\t')
      matchedSM.anno = paste(anno.df$ID[closest.SM], anno.df$Note[closest.SM], sep='|', collapse='...')
      if (method=='WindowLocal'){
        clean.table[nc, ] = c(nc,chr, paste(min.l, '-', max.l), 
                              paste(first.gene, '-', last.gene), min.dist.to.SM, 
                              matchedSM.anno, min(anno.df[i.match,paste('p.value(succ_local', window.size, ')',sep='')]), cluster.anno)
      }else{
        clean.table[nc, ] = c(nc,chr, paste(min.l, '-', max.l), 
                              paste(first.gene, '-', last.gene), min.dist.to.SM, 
                              matchedSM.anno, min(anno.df[i.match,'p.value(run_len)']), cluster.anno)        
      }
      
    }  
    write(x='#Some of the predicted clusters are overlapping. They may indicate a larger cluster if the clusters significantly overlap (according to the coordiates in column 3).', file=file.out, append=F)
    write(x='#Column 5 gives the distance of the cluster to the closest known secondary metabolite genes', file=file.out, append=T)
    write(x='#Column 5, 0 means known SM genes are within the predicted cluster', file=file.out, append=T)
    write(x='#Column 6 gives the gene names and annotations of the closest SM gene(s)', file=file.out, append=T)
    write(x='#Column 5 and column 6, when there are multiple closest SM genes, they are separated by ...', file=file.out, append=T)
    write(x='#Column 8+ gives the gene names and annotations of the genes in the predicted cluster', file=file.out, append=T)
    write(x=paste('#Estimated No. true NP gene clusters:',npos.window.local), file=file.out, append=T)
    write.table(clean.table, file=file.out,sep='\t', row.names=F, quote = F, append=T)    
    # n.SM.cluster = sum((diff(which(is.SM))>1) | (diff.str(anno.df$seqnames[is.SM])))+1
    # number of known SM gene clusters cannot be determined accurately
    out = c(sum(is.SM),sum(is.PKS), sum(i.new.NPG & is.SM), 
            sum(i.new.NPG & is.PKS), n.correct.cluster);
    names(out) = c('#known SM genes', '#known PKSs',
                   paste('#matched SM genes:', method, sep=''),
                   paste('#matched PKS genes:', method, sep=''),
                   paste('#matched SM clusters:', method, sep=''))
    return(out)
  }
  
  a = write.NPGC(anno.df, i.new.NPG = to.output.windows, window.size=window.size, method='WindowLocal',
                 file.out=paste('cluster.annoCompact.', tag, '.p', p.value.cutoff, '.NWindowClusters',n.clusters.localwindows, '.tab',sep=''))
  b = write.NPGC(anno.df, i.new.NPG = to.output.runs, window.size=window.size,method='Run',
                 file.out=paste('cluster.annoCompact.', tag, '.p', p.value.cutoff, '.NRunClusters',n.clusters.runs, '.tab',sep=''))
  n.unknowns = sum(regexpr(pattern='Protein of unknown function', text = annotation.text, perl=T)>0) # 20140529
  n.genes = length(anno)
  return(list(stats = c('#Pos Run Clusters'=npos.run, 'p Pos Run Clusters'=p.runCluster, 
                        '#Pos WindowLocal Clusters'=npos.window.local, 'p Pos WindowLocal Clusters'=p.window.localCluster,
                        "#Top Run Clusters"=n.clusters.runs, "#Top WindowLocal Clusters"=n.clusters.localwindows,
                        a, b[3:5],
                        '#Protein of unknown function'=n.unknowns,'#genes'=n.genes, 'enzyme prob'=sum(is.enzyme)/length(anno)),
              npos.run.simu=npos.run.simu, npos.success.simu=npos.success.simu, 
              npos.run.simus=npos.run.simus, npos.success.simus=npos.success.simus, n.chr = length(chrs)))
}

express.clustering <- function(gff.file="/Users/yongli/Universe/data/NPgenome/Aspergillus/A_nidulans_FGSC_A4_current_features.gff", geMat, iters = 5){
  # detect spacial clustering behavior of genes expression levels
  # 20140729, YF Li
  
  require(gplots)
  ## read gff
  gff.format = sub('^.*\\.([^\\.]*$)', '\\1', gff.file)
  # anno = read.gff3(gff.file, format=gff.format)
  anno = import.gff(gff.file) # 20160502
  idx.gene = (anno$type=='gene')
  anno = anno[idx.gene, ]
  anno = sort.intervals(anno)
  n = length(anno)
  
  ## filter genes
  idx = !is.na(match(anno$ID, rownames(geMat)))
  IDs = anno$ID[idx]
  geMat = geMat[IDs,]
  
  ## get gene modules
  require("fcanalysis",lib="~/Dropbox/Galaxy/R/lib")
  geMat.n =preprocessICA(geMat,F)
  s = ncol(geMat.n)-1
  ica.spatial = ica.do(geMat.n, iters = iters, nComponents = s)
  
  ## analyze the spacial autocorrelation for each sample and each gene module
  autocorr.all = zeros(n = s)
  names(autocorr.all) = 1:s
  autocorr.all.20 <- autocorr.all.R2 <- autocorr.all.Z <- autocorr.all
  pdf('/Users/yongli/Universe/write/Project_Current/t.NPbioinformatics/Nidulans.SlidingWindow/Autocorr.ICAmodules.pdf',8,6)
  par(mfrow=c(2,5))
  for (i in 1:s){
    lag.max = 60
    a = acf(ica.spatial$S[,i],lag.max = lag.max, main=paste('M', i, ' V%:', round(ica.spatial$power[i],4), sep=''))
    autocorr.all[i] = mean(a$acf[2:(lag.max+1)])
    autocorr.all.R2[i] = sqrt(mean(a$acf[2:(lag.max+1)]^2))
    autocorr.all.Z[i] = sum(1/2*log((1+a$acf[2:(lag.max+1)])/(1-a$acf[2:(lag.max+1)]))*sqrt(nrow(ica.spatial$S)-3))/sqrt(lag.max)
    lag.max = 20;
    a = acf(ica.spatial$S[,i],lag.max = lag.max, main=paste('M', i, ' V%:', round(ica.spatial$power[i],4), sep=''), plot=F)
    autocorr.all.20[i] = mean(a$acf)
  }  
  colnames(ica.spatial$A) = sub('nidulans', '', colnames(ica.spatial$A))
  colnames(ica.spatial$A) = sub('.CEL', '', colnames(ica.spatial$A))
  dev.off()
  
  pdf('/Users/yongli/Universe/write/Project_Current/t.NPbioinformatics/Nidulans.SlidingWindow/Autocorr.MeanExpression.pdf',4,4)
  par(mfrow=c(1,1))
  a = acf(rowMeans(geMat),lag.max = 100, main=paste('average gene expression'))
  dev.off()
  
  pdf('/Users/yongli/Universe/write/Project_Current/t.NPbioinformatics/Nidulans.SlidingWindow/clustering.ICA.pdf',6,8)
  heatmap.quick.geMat(ica.spatial$A, id.type = 'symbol', color = bluered(256),sd.cutoff = 0, margins=c(9,2))
  dev.off()
  autocorr.all = (autocorr.all - min(autocorr.all))/(max(autocorr.all)-min(autocorr.all))
  autocorr.all.20 = (autocorr.all.20 - min(autocorr.all.20))/(max(autocorr.all.20)-min(autocorr.all.20))
  
  ## examine known PKS
  is.SM = regexpr(pattern='secondary metab', text = as.character(as.vector(anno$Note[idx])), perl=T, ignore.case=T)>0
  is.PKS = regexpr(pattern='polyketide synthase', text = as.character(as.vector(anno$Note[idx])), perl=T, ignore.case=T)>0
  
  asso.FET = TFAI.FET(ica.spatial$S, mod.full = cbind(SM=is.SM, PKS=is.PKS))
  asso.lm = TFAI.lm(ica.spatial$S, mod.full = cbind(SM=is.SM, PKS=is.PKS), lm.joint = F, normalize.TF = 'none')
  asso.mu.lm = TFAI.lm(rowMeans(geMat), mod.full = cbind(SM=is.SM, PKS=is.PKS), lm.joint = F, normalize.TF = 'none')
  pdf('/Users/yongli/Universe/write/Project_Current/t.NPbioinformatics/Nidulans.SlidingWindow/Autocorr.vs.knownSM.Enrichment.pdf',4,4)
  plot(autocorr.all, -log10(asso.FET$p.value[1,]), xlab='normalized autocorr', ylab='-log10(p-value) FET', main='all SM genes')
  plot(autocorr.all, -log10(asso.FET$p.value[2,]), xlab='normalized autocorr', ylab='-log10(p-value) FET', main='PKS enzymes')
  plot(autocorr.all, -log10(asso.lm$p.value[1,]), xlab='normalized autocorr', ylab='-log10(p-value)', main='all SM genes')
  plot(autocorr.all, -log10(asso.lm$p.value[2,]), xlab='normalized autocorr', ylab='-log10(p-value)', main='PKS enzymes')
  plot(autocorr.all, asso.lm$statistic[1,], xlab='normalized autocorr', ylab='T statistics', main='all SM genes')
  plot(autocorr.all, asso.lm$statistic[2,], xlab='normalized autocorr', ylab='T statistics', main='PKS enzymes')
  plot(autocorr.all, abs(asso.lm$statistic[1,]), xlab='normalized autocorr', ylab='T statistics', main='all SM genes')
  plot(autocorr.all, abs(asso.lm$statistic[2,]), xlab='normalized autocorr', ylab='T statistics', main='PKS enzymes')
  dev.off()
  cor(autocorr.all, -asso.FET$p.value[1,], method = 'spearman')
  cor(autocorr.all, -asso.lm$p.value[1,], method = 'spearman')
  cor(autocorr.all, abs(asso.lm$statistic[1,]), method = 'spearman')
  
  
  cor(autocorr.all.R2, -asso.FET$p.value[1,], method = 'spearman')
  cor(autocorr.all.R2, -asso.lm$p.value[1,], method = 'spearman')
  cor(autocorr.all.R2, abs(asso.lm$statistic[1,]), method = 'spearman')
  
  cor(autocorr.all.Z, -asso.FET$p.value[1,], method = 'spearman')
  cor(autocorr.all.Z, -asso.lm$p.value[1,], method = 'spearman')
  cor(autocorr.all.Z, abs(asso.lm$statistic[1,]), method = 'spearman')
  
  cor(autocorr.all.20, -asso.FET$p.value[1,], method = 'spearman')
  cor(autocorr.all.20, -asso.lm$p.value[1,], method = 'spearman')
  cor(autocorr.all.20, abs(asso.lm$statistic[1,]), method = 'spearman')
  
  venn(list(module16 = which(ica.spatial$S[,16]>3), Known.NPG=which(is.SM)))
  venn(list(module16 = which(ica.spatial$S[,33]>3), Known.NPG=which(is.SM)))
  venn(list(module16 = which(ica.spatial$S[,34]>3), Known.NPG=which(is.SM)))
  venn(list(module16 = which(ica.spatial$S[,13]>3), Known.NPG=which(is.SM)))
  venn(list(module16 = which(ica.spatial$S[,32]>3), Known.NPG=which(is.SM)))
  pdf('/Users/yongli/Universe/write/Project_Current/t.NPbioinformatics/Nidulans.SlidingWindow/ExpressionLevel.NPvOther.pdf',4,4)
  hist.by(rowMeans(geMat[IDs,]), as.factor(c('No', 'Yes')[is.SM+1]), by.name = 'NPG', xlab='expression')
  hist.by(rowMeans(geMat[IDs,1:36]), as.factor(c('No', 'Yes')[is.SM+1]), by.name = 'NPG', main='liquid medium', xlab='expression')
  hist.by(rowMeans(geMat[IDs,37:44]), as.factor(c('No', 'Yes')[is.SM+1]), by.name = 'NPG', main='solid medium', xlab='expression')
  hist.by(rowMeans(geMat.n[,37:44])- rowMeans(geMat.n[,1:36]), as.factor(c('No', 'Yes')[is.SM+1]), by.name = 'NPG', main='solid/liquid difference', xlab='change')
  dev.off()
  
  ica.spatial$autocorr$R = autocorr.all;
  ica.spatial$autocorr$R2 = autocorr.all.R2;
  ica.spatial$autocorr$Z = autocorr.all.Z;  
  ica.spatial$spatial.cluster.index = autocorr.all;
  ica.spatial$spatial.cluster.method = 'mean'
  
  ica.spatial$anno = anno;
  ica.spatial$autocorr.lag = 60;
  ica.spatial$geMat = geMat
  return(ica.spatial)
}

score.spatial.cluster <- function(ica.spatial, gene.range=c('AN8131', 'AN8137'),
                                  score.type = c('R', 'R2', 'Z'), median.substraction=T, do.plot=T){
  # compute the clustering scores of a given gene range
  # YF Li, 20140731
  
  score.type = match.arg(score.type)
  spatial.cluster.score = ica.spatial$autocorr[[score.type]]
  if (median.substraction){
    spatial.cluster.score = spatial.cluster.score - median(spatial.cluster.score)
  }
  
  is = match(gene.range, rownames(ica.spatial$S))
  is = sort(is)
  
  gene.range = rownames(ica.spatial$S)[is]
  #   s1 = sum(colSums(ica.spatial$S[is[1]:is[2],]^2)*ica.spatial$autocorr)
  #   s2 = sum(colSums(ica.spatial$S[is[1]:is[2],])^2*ica.spatial$autocorr)
  k = is[2]-is[1]+1
  rs.unsigned = (apply(ica.spatial$S2, MARGIN = 2, function(x){y = as.vector(runsum(Rle(x), k)); names(y) = names(x)[1:(length(x)-k+1)]; return(y)}) %*% spatial.cluster.score)[,1]
  rs.signed = (apply(ica.spatial$S, MARGIN = 2, function(x){y = as.vector(runsum(Rle(x), k)); names(y) = names(x)[1:(length(x)-k+1)]; return(y)})^2 %*% spatial.cluster.score)[,1]
  fdr.signed = fdr.symmatric(log2(rs.signed),iterative = F, plot = do.plot)
  fdr.unsigned = fdr.symmatric(log2(rs.unsigned),iterative = F, plot = do.plot)
  names(fdr.signed) <- names(fdr.unsigned) <- names(rs.unsigned)
  return(c(s.unsigned = rs.unsigned[gene.range[1]],
           s.signed = rs.signed[gene.range[1]],
           fdr.unsigned = fdr.unsigned[gene.range[1]],
           fdr.signed = fdr.signed[gene.range[1]]))
}


ica.spatial.prep <- function(ica.spatial, K= 50, center.method='median',
                             score.type = c('R', 'R2', 'Z'), median.substraction=F, do.plot=T){
  require(modeest)
  score.type = match.arg(score.type)
  spatial.cluster.score = ica.spatial$autocorr[[score.type]]
  if (median.substraction){
    spatial.cluster.score = spatial.cluster.score - median(spatial.cluster.score)
    spatial.cluster.score = spatial.cluster.score/mad(spatial.cluster.score)
    spatial.cluster.score[spatial.cluster.score<1] = 0
  }
  ica.spatial$S2 = ica.spatial$S^2;
  ica.spatial$rs.unsigned <- ica.spatial$fdr.unsigned <- ica.spatial$rs.signed <- 
    ica.spatial$fdr.signed <- ica.spatial$rs.LowExpression <- ica.spatial$fdr.LowExpression <- c()
  mu = rowMeans(ica.spatial$geMat) # average expression
  for (k in 1:K){  
    rs.unsigned = (apply(ica.spatial$S2, MARGIN = 2, function(x){return(runsum.2(x,k,addzeros=T))}) %*% spatial.cluster.score)[,1]
    if (mean(rs.unsigned<0)<0.01)
      rs.unsigned = log2(rs.unsigned)
    fdr.unsigned = fdr.symmatric(rs.unsigned,iterative = F, plot = do.plot, center.method=center.method)
    rs.signed = (apply(ica.spatial$S, MARGIN = 2, function(x){return(runsum.2(x,k,addzeros=T))})^2 %*% spatial.cluster.score)[,1]    
    if (mean(rs.signed<0)<0.01)
      rs.signed = log2(rs.signed)
    fdr.signed = fdr.symmatric(rs.signed,iterative = F, plot = do.plot, center.method=center.method)
    
    rs.LowExpression = -runsum.2(mu,k=k,addzeros=T)# low expression score
    fdr.LowExpression = fdr.symmatric(rs.LowExpression,iterative = F, plot = do.plot, center.method=center.method)
    
    names(fdr.LowExpression)  <- names(fdr.signed) <- names(fdr.unsigned) <- names(rs.unsigned)
    ica.spatial$rs.unsigned <- cbind(ica.spatial$rs.unsigned, rs.unsigned)
    ica.spatial$fdr.unsigned <- cbind(ica.spatial$fdr.unsigned, fdr.unsigned);
    ica.spatial$rs.signed <- cbind(ica.spatial$rs.signed, rs.signed)
    ica.spatial$fdr.signed <- cbind(ica.spatial$fdr.signed, fdr.signed);
    ica.spatial$rs.LowExpression <- cbind(ica.spatial$rs.LowExpression, rs.LowExpression)
    ica.spatial$fdr.LowExpression <- cbind(ica.spatial$fdr.LowExpression, fdr.LowExpression)
  }
  rownames(ica.spatial$rs.unsigned)  <- rownames(ica.spatial$fdr.unsigned)  <-
    rownames(ica.spatial$rs.signed)  <- rownames(ica.spatial$fdr.signed)  <-
    rownames(ica.spatial$rs.LowExpression)  <- rownames(ica.spatial$fdr.LowExpression)  <- names(rs.unsigned)
  ica.spatial$mu = mu
  ica.spatial$center.method=center.method;
  ica.spatial$score.type = score.type
  ica.spatial$median.substraction=median.substraction
  return(ica.spatial)
}
score.spatial.cluster.2d <- function(ica.spatial, gene.range=c('AN8131', 'AN8137'), 
                                     cor.method = 'spearman', CS.n.neighbor = 3){
  # compute the clustering scores of all possible windows within a given gene range
  # YF Li, 20140730
  require(modeest)
  is = match(gene.range, rownames(ica.spatial$S))
  is = sort(is)
  gene.range = rownames(ica.spatial$S)[is]
  all.genes = rownames(ica.spatial$S)[is[1]:is[2]]
  #   s1 = sum(colSums(ica.spatial$S[is[1]:is[2],]^2)*ica.spatial$autocorr)
  #   s2 = sum(colSums(ica.spatial$S[is[1]:is[2],])^2*ica.spatial$autocorr)
  K = is[2]-is[1]+1
  s.unsigned.2d <- s.signed.2d <- fdr.unsigned.2d <- fdr.signed.2d <- 
    s.LowExpression.2d <- fdr.LowExpression.2d <- matrix(NA, nrow = K, ncol = K, dimnames = list(all.genes, all.genes)) 
  for (k in 1:K){      
    for (i in is[1]:(is[2]-k+1)){
      s.unsigned.2d[all.genes[i-is[1]+1], all.genes[i+k-is[1]]] = ica.spatial$rs.unsigned[all.genes[i-is[1]+1],k];
      fdr.unsigned.2d[all.genes[i-is[1]+1], all.genes[i+k-is[1]]] = ica.spatial$fdr.unsigned[all.genes[i-is[1]+1],k];
      s.signed.2d[all.genes[i-is[1]+1], all.genes[i+k-is[1]]] = ica.spatial$rs.signed[all.genes[i-is[1]+1],k];
      fdr.signed.2d[all.genes[i-is[1]+1], all.genes[i+k-is[1]]] = ica.spatial$fdr.signed[all.genes[i-is[1]+1],k];
      
      s.LowExpression.2d[all.genes[i-is[1]+1], all.genes[i+k-is[1]]] = ica.spatial$rs.LowExpression[all.genes[i-is[1]+1],k];
      fdr.LowExpression.2d[all.genes[i-is[1]+1], all.genes[i+k-is[1]]] = ica.spatial$fdr.LowExpression[all.genes[i-is[1]+1],k];
    }
  }
  R.ext = cor(t(ica.spatial$geMat[max((is[1]-CS.n.neighbor),1):min(is[2]+CS.n.neighbor,nrow(ica.spatial$geMat)),]), method = cor.method)
  R = R.ext[all.genes, all.genes]
  R.ext[R.ext<0] = 0;
  ai = arrayInd(1:length(R.ext),.dim = dim(R.ext));
  ai = ai[abs(ai[,1]-ai[,2])>CS.n.neighbor | ai[,1]==ai[,2],]
  R.ext[ai] = 0;
  CS = rowSums(R.ext^2)[all.genes]
  
  is.anno = match(gene.range, ica.spatial$anno$ID)
  is.anno =  sort(is.anno) 
  all.genes.anno = ica.spatial$anno$ID[is.anno[1]:is.anno[2]]
  
  return(list(s.unsigned = s.unsigned.2d,
              s.signed = s.signed.2d,
              fdr.unsigned = fdr.unsigned.2d,
              fdr.signed = fdr.signed.2d,
              s.lowExpress = s.LowExpression.2d,
              fdr.lowExpress = fdr.LowExpression.2d,
              cor = R, 
              CS = CS,
              mu = ica.spatial$mu[all.genes],
              sd = rowSds(ica.spatial$geMat[all.genes,]),
              err = rowSds(ica.spatial$geMat[all.genes,])/sqrt(ncol(ica.spatial$geMat)),
              geMat = ica.spatial$geMat[all.genes,],
              center.method=ica.spatial$center.method,
              score.type = ica.spatial$score.type, 
              median.substraction=ica.spatial$median.substraction,
              cor.method = cor.method, CS.n.neighbor = CS.n.neighbor,
              all.gene.geMat = all.genes,
              all.gene.anno = all.genes.anno))
}

plot.spatial.cluster.2d <- function(s2d, col = bluered(256), tag='',
                                    heatmap.clustering=T, no.fdr=F){
  # visualize the expression clustering for all pairwise windows
  require(lattice);
  min.logp = 2
  n.color = 32
  p = s2d$fdr.unsigned;
  p[p==0] = min(p[p!=0 & !is.na(p)])/2
  # get proper color scale so that p-value 0.1 is assigned the middle color
  max.logp = max(max(-log10(p),na.rm = T),min.logp-log10(5))+log10(5);
  l1 = round((length(col)+1)/2);  l2 = length(col);  l0 = (l1*max.logp-l2)/(max.logp-1)
  col0 = col[l0:l2]  
  x4 = levelplot(t(-log10(p)),col.regions = col0, xlab='to gene', ylab='from gene', main='-log10(fdr unsiged cluster score)',
                 scales=list(x=list(rot=90),alternating=1),at=seq(0, max.logp, length.out=n.color))
  p = s2d$fdr.signed;
  p[p==0] = min(p[p!=0 & !is.na(p)])/2
  max.logp = max(max(-log10(p),na.rm = T),min.logp-log10(5))+log10(5);
  l1 = round((length(col)+1)/2);  l2 = length(col);  l0 = (l1*max.logp-l2)/(max.logp-1)
  col0 = col[l0:l2] 
  x5 = levelplot(t(-log10(p)),col.regions = col0, xlab='to gene', ylab='from gene', main='-log10(fdr siged cluster score)',
                 scales=list(x=list(rot=90),alternating=1),at=seq(0, max.logp, length.out=n.color))
  p = s2d$fdr.lowExpress;
  p[p==0] = min(p[p!=0 & !is.na(p)])/2
  max.logp = max(max(-log10(p),na.rm = T),min.logp-log10(5))+log10(5);
  l1 = round((length(col)+1)/2);  l2 = length(col);  l0 = (l1*max.logp-l2)/(max.logp-1)
  col0 = col[l0:l2] 
  x6 = levelplot(t(-log10(p)),col.regions = col0, xlab='to gene', ylab='from gene', main='-log10(fdr low express score)',
                 scales=list(x=list(rot=90),alternating=1),at=seq(0, max.logp, length.out=n.color))
  R = s2d$cor; R[lower.tri(R)] <- NA; diag(R) <- NA
  x3 = levelplot(t(R),col.regions = col, xlab='gene A', ylab='gene B', main='pairwise correlation',
                 scales=list(x=list(rot=90),alternating=1))
  x2 = error.bar(s2d$mu, err = s2d$sd, ylab='average expression', main='');
  x2.2 = barchart(s2d$CS~factor(names(s2d$CS), levels = names(s2d$CS)), ylab=paste('Andersen CS', s2d$cor.method), scales = list(x = list(draw = FALSE)), main= tag);
  f = mat2xyz(s2d$geMat, sym=F)
  ng = nrow(s2d$geMat)
  x1 = xyplot(z~y, group=x, data=f,type='l', scales=list(x=list(rot=90),alternating=1), par.settings = list(superpose.line = list(lwd = 3)), col=greenred(ng), at = seq(1, ng, length = ng), xlab='sample',ylab='expression',
              panel = function(...) { 
                panel.text(1, max(s2d$geMat,na.rm =T), "color maps to gene order", pos=4) 
                panel.xyplot(...) 
              })    
  x1.2 = heatmap.lattice(s2d$geMat, top = F, col.regions = col)   
  x = scale(t(s2d$geMat),center = T,scale = F); m = max(abs(x))
  x1.3 = levelplot(x, scales = list(x = list(rot = 90),alternating=1),xlab='', ylab='',
                   colorkey = F, at = seq(-m, m, length = 32), 
                   aspect = 'fill', col.regions = col)
  if (no.fdr){
    print(x1.3, split=c(1,1,3,1), newpage=T)
    print(x1.2, split=c(3,1,3,1), newpage=F)
    print(x2, split=c(2,2,3,2), newpage=F)
    print(x2.2, split=c(2,1,3,2), newpage=F)
    # print(x3, split=c(3,1,3,2), newpage=F)
    
  }else{
    print(x1.3, split=c(1,1,3,2), newpage=T)
    print(x1.2, split=c(3,1,3,2), newpage=F)
    print(x2, split=c(2,2,3,4), newpage=F)
    print(x2.2, split=c(2,1,3,4), newpage=F)
    # print(x3, split=c(3,1,3,2), newpage=F)
    print(x4, split=c(1,2,3,2), newpage=F)
    print(x5, split=c(2,2,3,2), newpage=F)
    print(x6, split=c(3,2,3,2), newpage=F)
    
  }
  
  k = ncol(s2d$fdr.unsigned)
  #   trellis.focus("panel",column = 1,row=1)
  #   panel.text(cex=1, x=(1:k), y=(1:k), labels=rownames(s2d$fdr.unsigned), xpd=TRUE, srt=0, pos=1)
  #   trellis.unfocus()
  #   print(lattice::levelplot(t(s2d$s.signed),col.regions = bluered(32), xlab='from', ylab='to', main='siged cluster score',
  #                      scales=list(x=list(rot=90))))
  #   print(lattice::levelplot(t(s2d$s.unsigned),col.regions = bluered(32), xlab='from', ylab='to', main='unsiged cluster score',
  #                      scales=list(x=list(rot=90))))
}

read.gff3 <- function(con, format=sub('^.*\\.([^\\.]*$)', '\\1', con),
                      genome = NA, asRangedData = F, colnames = NULL,
                      which = NULL, feature.type = NULL){
  # modified from rtklayer to handle the quotation mark bug: 
  # original version: https://github.com/genome-vendor/r-bioc-rtracklayer/blob/master/R/gff.R
  # 20140502
  # Yong Fuga Li
  require('rtracklayer')
  lines <- readLines(con, warn = FALSE) # unfortunately, not a table
  lines <- lines[nzchar(lines)]
  
  ## strip comments
  notComments <- which(substr(lines, start=1L, stop=1L) != "#")
  lines <- lines[notComments]
  
  ### TODO: handle ontologies (store in RangedData)
  
  ## strip FASTA sequence
  fastaHeaders <- which(substr(lines, start=1L, stop=1L) == ">")
  if (length(fastaHeaders))
    lines <- head(lines, fastaHeaders[1] - 1)
  
  ## construct table
  fields <- c("seqname", "source", "type", "start", "end", "score",
              "strand", "phase", "attributes")
  linesSplit <- strsplit(lines, "\t", fixed=TRUE)
  fieldCounts <- elementLengths(linesSplit)
  if (any(fieldCounts > length(fields)) ||
      any(fieldCounts < (length(fields) - 1)))
    stop("GFF files must have ", length(fields),
         " tab-separated columns")
  haveAttr <- fieldCounts == length(fields)
  data <- unlist(linesSplit[haveAttr], use.names=FALSE)
  if (is.null(data))
    data <- character(0)
  haveAttrMat <- matrix(data, ncol=length(fields), byrow=TRUE)
  data <- unlist(linesSplit[!haveAttr], use.names=FALSE)
  if (is.null(data))
    data <- character(0)
  noAttrMat <- matrix(data, ncol=length(fields)-1L, byrow=TRUE)
  noAttrMat <- cbind(noAttrMat, rep.int("", nrow(noAttrMat)))
  table <- rbind(noAttrMat, haveAttrMat)
  colnames(table) <- fields
  
  if (!is.null(feature.type))
    table <- table[table[,"type"] %in% feature.type,,drop=FALSE]
  
  ## handle missings
  table[table == "."] <- NA_character_
  
  attrCol <- table[,"attributes"]
  if (format=='gff3') {
    table <- table[,setdiff(colnames(table), "attributes"),drop=FALSE]
    table[table[,"strand"] == "?","strand"] <- NA_character_
    is_not_NA <- !is.na(table)
    table[is_not_NA] <- urlDecode(table[is_not_NA])
  }
  table[is.na(table[,"strand"]),"strand"] = '*'
  
  extraCols <- c("source", "type", "score", "strand", "phase")
  if (!is.null(colnames))
    extraCols <- intersect(extraCols, colnames)
  xd <- as(table[,extraCols,drop=FALSE], "DataFrame")
  
  if (!is.null(xd$phase))
    xd$phase <- as.integer(as.character(xd$phase))
  if (!is.null(xd$strand))
    xd$strand <- strand(xd$strand)
  if (!is.null(xd$score))
    suppressWarnings(xd$score <- as.numeric(as.character(xd$score)))
  
  if (is.null(colnames) || length(setdiff(colnames, extraCols))) {
    if (format=='gff1') {
      if (is.null(colnames) || "group" %in% colnames)
        attrList <- list(group = factor(attrCol,
                                        levels=unique(attrCol)))
      else attrList <- list()
    } else {
      attrSplit <- strsplit(attrCol, ";", fixed=TRUE)
      attrs <- unlist(attrSplit, use.names=FALSE)
      lines <- rep.int(seq_len(length(attrSplit)),
                       elementLengths(attrSplit))
      attrs <- sub(" *$", "", sub("^ *", "", attrs))
      if (format=='gff3') {
        equals.pos <- regexpr("=", attrs, fixed=TRUE)
        if (any(equals.pos == -1L))
          stop("Some attributes do not conform to 'tag=value' format")
        tags <- substring(attrs, 1L, equals.pos - 1L)
        vals <- substring(attrs, equals.pos + 1L, nchar(attrs))
      } else { # split on first space (FIXME: not sensitive to quotes)
        tags <- sub(" .*", "", attrs) # strip surrounding quotes
        vals <- sub("^\"([^\"]*)\"$", "\\1",
                    sub("^[^ ]* ", "", attrs))
      }
      if (!is.null(colnames)) {
        keep <- tags %in% colnames
        lines <- lines[keep]
        vals <- vals[keep]
        tags <- urlDecode(tags[keep])
      }
      tags <- factor(tags, levels=unique(tags))
      lineByTag <- split(lines, tags)
      valByTag <- split(vals, tags)
      
      ## FIXME: Parent, Alias, Note, DBxref,
      ## Ontology_term are allowed to have multiple
      ## values. We should probably always return them as a
      ## CharacterList.
      multiTags <- c("Parent", "Alias", "Note", "DBxref",
                     "Ontology_term")
      attrList <- sapply(names(lineByTag), function(tagName) {
        vals <- valByTag[[tagName]]
        if (format=='gff3' &&
            (any(grepl(",", vals, fixed=TRUE)) ||
             tagName %in% multiTags)) {
          vals <- CharacterList(strsplit(vals, ",", fixed=TRUE))
          vals <- relist(urlDecode(unlist(vals)), vals)
          coerced <- suppressWarnings(as(vals, "NumericList"))
          if (!any(any(is.na(coerced))))
            vals <- coerced
          vec <- as(rep.int(list(character()), nrow(table)),
                    class(vals))
        } else {
          coerced <- suppressWarnings(as.numeric(vals))
          if (!any(is.na(coerced)))
            vals <- coerced
          if (format=='gff3')
            vals <- urlDecode(vals)
          vec <- rep.int(NA, nrow(table))
        }
        vec[lineByTag[[tagName]]] <- vals
        vec
      }, simplify = FALSE)
    }
    xd <- DataFrame(xd, attrList)
  }
  
  end <- as.integer(table[,"end"])
  GenomicData(IRanges(as.integer(table[,"start"]), end),
              xd, chrom = table[,"seqname"], genome = genome,
              seqinfo = attr(con, "seqinfo"),
              asRangedData = asRangedData)
}

urlDecode <- function(str)
{
  require('RCurl')
  curlUnescape(str)
}


promoter.statistics <- function(gff.file="A_nidulans_FGSC_A4_current_features.gff", 
                                DNA.fasta.file="A_nidulans_FGSC_A4_current_chromosomes.fasta",
                                window.promoter = c(-4000, 1000), k=8, n.top.motifs=10,
                                tag='window4k1k.8mer'){
  # computer based promoter statiscs: gene-gene distances orientations etc
  # k: k-mer size
  # Yong Fuga Li
  # 20140606
  require(Biostrings)
  require(markovchain)
  require(IRanges)
  require(ggplot2)
  require('TTR') # 20140527
  
  # read gff
  gff.format = sub('^.*\\.([^\\.]*$)', '\\1', gff.file)
  # anno = read.gff3(gff.file, format=gff.format)
  anno = import.gff(gff.file) # 20160502
  idx.gene = (anno$type=='gene')
  anno = anno[idx.gene, ]
  anno = sort.intervals(anno)
  n = length(anno)
  
  # read fasta
  fa = import(DNA.fasta.file,format='fasta')
  
  # are gene orientations independent?
  s4 = substr.stats(anno@strand, anno@seqnames)
  # are NP gene orientations different?
  is.SM = regexpr(pattern='secondary metab', text = as.character(as.vector(anno$Note)), perl=T)>0 # 20140519
  sum(is.SM)
  s4.SM = substr.stats(anno@strand, anno@seqnames, is.SM)
  pdf(paste('geneOrientation.pdf',sep=''),width=5,height=3.5)
  for (i in 1:4){
    da = rbind(data.frame(x=names(s4[[i]]), y=s4[[i]]/sum(s4[[i]])*100,genes='all'),
               data.frame(x=names(s4.SM[[i]]), y=s4.SM[[i]]/sum(s4.SM[[i]])*100,genes='NP'))
    q = ggplot(data=da,  aes(x = x, y=y, by=genes, fill=genes))+geom_bar(stat='identity',position='dodge')+
      xlab('orientation')+ylab('%')+theme(axis.text.x = element_text(angle=90, vjust=1))
    print(q)
  }
  dev.off()
  
  # intergene region lengths: for {+/+, -/-} vs {+/-, -/+} intergenic regions
  inter.dist = get.intergene.dist(anno, cutoff=10000)
  inter.dist.SM = get.intergene.dist(anno[is.SM], cutoff=10000)
  gene.len = end(anno@ranges) -  start(anno@ranges)
  pdf(paste('intergeneDist.pdf', sep=''))
  print(hist.by(inter.dist$dist, as.factor(inter.dist$type), by.name='gene orientations',hist=T,xlab='distance',main='all'))
  print(hist.by(inter.dist$dist, as.factor(inter.dist$type.brief), by.name='gene orientations', hist=F,xlab='distance',main='all'))
  print(hist.by(inter.dist$dist, as.factor(inter.dist$type.brief), by.name='gene orientations', hist=T,xlab='distance',main='all'))
  print(hist.by(inter.dist.SM$dist, as.factor(inter.dist.SM$type), by.name='gene orientations',hist=T,xlab='distance',main='NPGC'))
  print(hist.by(inter.dist.SM$dist, as.factor(inter.dist.SM$type.brief), by.name='gene orientations', hist=F,xlab='distance',main='NPGC'))
  print(hist.by(inter.dist.SM$dist, as.factor(inter.dist.SM$type.brief), by.name='gene orientations', hist=T,xlab='distance',main='NPGC'))
  hist.by(gene.len,by=is.SM,by.name='NP gene',main='gene length')
  dev.off()
  
  
  # motif findings around known NPGCs
  names(fa) = sub(pattern='^([^ ]*) .*$', replacement='\\1',names(fa))
  fa.promoter = get.promoter.seq(fa, anno[is.SM],k=window.promoter);
  SM.mstats = motif.stats(fa.promoter, l=k)
  n.shift = 100; is.SM.rand = (which(is.SM)+n.shift-1)%%n+1 # random promoters
  fa.promoter.rand = get.promoter.seq(fa, anno[is.SM.rand],k=window.promoter);
  n.shift = 191; is.SM.rand = (which(is.SM)+n.shift-1)%%n+1 # random promoters
  fa.promoter.rand2 = get.promoter.seq(fa, anno[is.SM.rand],k=window.promoter);
  SM.mstats.rand = motif.stats(fa.promoter.rand, l=k)
  SM.mstats.rand2 = motif.stats(fa.promoter.rand2, l=k)
  write.fasta(fa.promoter,  paste('A.nidulans.NPG.promoter.',tag,'.fa',sep=''))
  write.fasta(fa.promoter.rand,  paste('A.nidulans.rand.promoter.',tag,'.fa',sep=''))
  
  out = motif.comp(SM.mstats, SM.mstats.rand)
  out2 = motif.comp(SM.mstats, SM.mstats.rand2)
  
  # msets = motif.find(fa) # learn a motif sets
  m.anno = motif.annotate(fa.promoter, msets=out$fitered[1:n.top.motifs,]) # annotate sequences by motif sets
  m.anno.rand = motif.annotate(fa.promoter.rand, msets=out$fitered[1:n.top.motifs,]) # annotate sequences by motif sets
  pdf(paste('motifClust',tag,'.pdf',sep=''),width=4,16)
  heatmap.quick.geMat(t((m.anno$count[,colSums(m.anno$count)>0]>0)+0),centering=F,id.type='symbol', distfun=dist,sd.cutoff=-1, lhei=c(1, 14), margins=c(9,5))
  heatmap.quick.geMat(t(m.anno$loc.average[,colSums(m.anno$count)>0]),centering=F,id.type='symbol', distfun=dist,sd.cutoff=-1, lhei=c(1, 14), margins=c(9,5))
  dev.off()
  pdf(paste('motifLocation',tag,'.pdf',sep=''))
  nbins = floor(sqrt(sum(m.anno$loc.average>0)))
  hist(m.anno$loc.average[m.anno$loc.average>0]+window.promoter[1],xlab='distance to CDS',ylab='#Motifs',breaks=nbins,main='NPGC')
  abline(h=(diff(window.promoter)-k+2)/nbins/2^k*n.top.motifs, col='grey', lty='dashed')
  hist(m.anno.rand$loc.average[m.anno.rand$loc.average>0]+window.promoter[1],xlab='distance to CDS',ylab='#Motifs',breaks=nbins,main='rand')
  abline(h=(diff(window.promoter)-k+2)/nbins/2^k*n.top.motifs, col='grey', lty='dashed')
  dev.off()
  motif.associations = asso.FET(t(m.anno$count>0)+0) # testing associations among the k-mers  
  write.table(motif.associations,file=paste('motif.association.', tag, '.xls',sep=''), col.names=NA, sep='\t')
}

get.intergene.dist <- function(anno,cutoff){
  # 20140607
  n = length(anno)
  intergene.dist = -end(anno@ranges)[1:(n-1)] +   start(anno@ranges)[2:n]
  intergene.type = paste(as.vector(anno@strand)[1:(n-1)],as.vector(anno@strand)[2:(n)], sep='') 
  
  to.keep = (abs(intergene.dist) < cutoff)
  intergene.dist = intergene.dist[to.keep]; intergene.type = intergene.type[to.keep]
  ii = intergene.type=='-+'
  intergene.type2 = intergene.type;
  intergene.type2[intergene.type2 %in% c('++', '--')] = '++,--'
  intergene.type2[intergene.type2 %in% c('+-', '-+')] = '+-,-+'
  return(list(dist=intergene.dist, type= intergene.type, type.brief=intergene.type2))
}

motif.annotate <- function(fa.promoter, msets=out$fitered[1:10,]){
  # 20140610, annotate fasta sequences by a set of motifs
  # msets: k x 4 matrix describing a set of motifs
  # fa.promoter: fasta sequences
  n.seq = length(fa.promoter)
  n.motif = nrow(msets)
  motifs = rownames(msets)
  m.anno.count <- m.anno.loc <- matrix(0,nrow=n.motif,ncol=n.seq,dimnames=list(motifs=rownames(msets), seqs=sapply(strsplit(names(fa.promoter),split='\\|'),FUN=function(x){return(x[1])})))
  for (m in 1:n.motif){
    locs = gregexpr(motifs[m], fa.promoter, ignore.case=T)
    m.anno.loc[m,] <- sapply(locs,FUN=function(x){
      if(sum(x>0)>0)
        return(mean(x[x>0]))
      else
        return(0)
    })
    m.anno.count[m,] <- sapply(locs,FUN=function(x){return(sum(x>0))})
  }
  return(list(count=m.anno.count, 
              loc.average=m.anno.loc))
}

get.promoter.seq <- function(fa, anno, k = c(-2000,500)){
  # Yong Fuga Li 20140606
  require(Biostrings)
  fa.promoter = list()
  strands = as.vector(anno@strand)
  for (i in 1:length(anno)){
    if (anno$type[i]!='gene') # only use gene features
      next
    chr = as.character(anno@seqnames[i]);
    if (strands[i]=='+'){
      pseq = substr(fa[[chr]], max(1,start(anno[i])+k[1]),min(start(anno[i])+k[2], length(fa[[chr]])))
    }else{
      pseq = reverseComplement(substr(fa[[chr]], max(1,end(anno[i])-k[2]),min(end(anno[i])-k[1], length(fa[[chr]]))))
    }
    fa.promoter[[paste(anno$ID[i], '|', anno$Note[i],sep='')]] = as.character(pseq)          
  } 
  return(fa.promoter)
}

substr.stats <- function(s, chr, filter=NULL){
  if (!is.null(filter))
    s = s[filter]
  gs1 = as.vector(s)
  n = length(s)
  gs2 = paste(gs1[1:(n-1)], gs1[2:n], sep='')
  gs3 = paste(gs1[1:(n-2)], gs1[2:(n-1)], gs1[3:(n)], sep='')
  gs4 = paste(gs1[1:(n-3)], gs1[2:(n-2)], gs1[3:(n-1)], gs1[4:(n)], sep='')
  uc1 = unique.count(gs1)$counts.unique
  uc2 = unique.count(gs2)$counts.unique
  uc3 = unique.count(gs3)$counts.unique
  uc4 = unique.count(gs4)$counts.unique
  
  d = list()
  uc.d = list()
  d[['2']] = (gs1[1:(n-1)]!= gs1[2:n])+0 # transitions in window size 2
  uc.d[['2']] = unique.count(d[['2']])$counts.unique
  for (i in 3:10){ # number of transitions in window size i
    d[[paste(i)]] = as.vector(runsum(Rle(d[['2']]),i-1))
    uc.d[[paste(i)]] = unique.count(d[[paste(i)]])$counts.unique
  }
  return(list(uc1, uc2, uc3, uc4, uc.d))
}

motif.stats <- function(fa.promoter, l = 8){
  # fa.promoter: fasta sequences 
  # l: motif length
  # Yong Fuga Li, 20140606
  mstats = c();
  nt.freq = c() # nucleotide frequencies
  for (i in 1:length(fa.promoter)){
    mstats = sum.union(mstats, unique.count.substr(fa.promoter[[i]],l))
    nt.freq = sum.union(nt.freq, unique.count.substr(fa.promoter[[i]],1))
  }
  nt.freq = (nt.freq+1)/sum(nt.freq+1)
  p.motifs = mstats;
  motifs = names(mstats);
  for (i in 1:length(p.motifs)){ # obtain motif probability based on frequency model
    p.motifs = prod(nt.freq[strsplit(motifs[i], '')[[1]]])
  }
  n = sum(mstats)
  p = 1-pbinom(mstats-1, size=n, prob=p.motifs)
  p.neg = pbinom(mstats, size=n, prob=p.motifs)
  out = cbind(motif.count=mstats, p.value=p, p.value.deplete=p.neg)
  out = out[order(out[,2],decreasing=F),]
  cat(nt.freq)
  return(out)
}

motif.comp <- function(SM.mstats, SM.mstats.rand){  
  # compare motif stats results
  mstats = cbind.union(SM.mstats[,1],SM.mstats.rand[,1])
  colnames(mstats) = c('SM', 'rand')
  tot = colSums(mstats)
  p = pbinom(mstats[,2],rowSums(mstats),prob=tot[2]/sum(tot),lower.tail=T)
  q = qvalue.2(p)
  p.deplete = pbinom(mstats[,1],rowSums(mstats),prob=tot[1]/sum(tot),lower.tail=T)
  q.deplete = qvalue.2(p.deplete)
  out = cbind(mstats, fold =mstats[,1]/mstats[,2], p.value=p, p.deplete=p.deplete, q.value=q, q.deplete=q.deplete)
  out.cut = out[q.deplete<0.1 | p<0.001,]
  out.cut = out.cut[order(out.cut[,4],decreasing=F),]
  out = out[order(out[,4], decreasing=F),]
  return(list(all = out, fitered=out.cut))
} 

unique.count.substr <- function(string, l){
  # count the number of each length l unique substrings
  sq = strsplit(string,split='')[[1]]
  n = length(sq)
  subsq = sq[1:(n-l+1)]; # get length l sub-strings
  for (j in seq2(from=2, to=min(l, n), by=1)){
    subsq = paste(subsq,sq[j:(n-l+j)], sep='')
  }
  # count sub strings
  return(unique.count(subsq)$counts.unique)
}

sort.intervals <- function(anno, do.strand=F){
  # sort Genome intervals by seqnames and ranges
  # Yong Fuga Li, 20140606
  # by location
  i = order(anno@ranges)
  anno = anno[i,]
  # by chr
  i = order(as.character(anno@seqnames))
  anno = anno[i,]
  
  # by strand
  if (do.strand){
    i = order(anno@strand)
    anno = anno[i,]
  }
  return(anno)
}

sort.gff <- function(gff.file, format = 'gff3', out.file = sub('.([^\\.]+)$', '_sorted.\\1',gff.file), do.strand=F){
  # 20150916, sort GFF file, used in script analysis.KU2015.RORA.R to handle unsorted features from GenBank
  # Yong Fuga Li
  # anno = read.gff3(gff.file, format = format)
  anno = import.gff(gff.file) # 20160502
  anno = sort.intervals(anno, do.strand = do.strand)
  export(anno, out.file, format = 'gff3', append=F) 
  return(out.file)
}

read.orthologs <- function(desc.file = 'desc.txt', 
                           ortholog.file = 'All_Species_Orthologs_by_Jaccard_clustering.txt', root = '/Users/yongli/Universe/data/NPgenome/Aspergillus'){
  # desc.file: species -- gff file mappings
  # ortholog.file: ortholog groups
  # 20140604
  
  require('rtracklayer')
  require('genomeIntervals')
  require(lattice)
  
  
  # read gff files
  desc = read.table(desc.file, header=T, sep='\t', as.is=T)
  n.species = ncol(ortho)-1; # n
  
  # read orthologs
  ortho = read.table(ortholog.file, header=T, sep='\t', as.is=T)
  conservativity = (rowSums(ortho[,2:ncol(ortho)]!='')-1)/(n.species-1) # in all species ==> 1, in one species ==> 0
  
  # read all genome annotations
  anno.all = list() 
  for (i in 1:length(desc$gff)){
    # i = 11
    species = desc$species[i]
    gff.file = desc$gff[i]
    gff.format = sub('^.*\\.([^\\.]*$)', '\\1', gff.file)
    # anno = read.gff3(gff.file, format=gff.format)
    anno = import.gff(gff.file) # 20160502
    idx.gene = (anno$type=='gene')
    anno = anno[idx.gene, ]
    anno = sort.intervals(anno) # 20140606 sort to sort.intervals
    #     is.enzyme.ase = regexpr(pattern='(?: |^)[^ ]+ase(?: |$)', text = as.character(as.vector(anno$Note)), perl=T)>0 # 20140519
    #     is.enzyme.EC6 = regexpr(pattern='(oxidoreductase|transferase|hydrolase|lyase|isomerase|ligase)', text = as.character(as.vector(anno$Note)), perl=T, ignore.case=T) > 0 
    #     is.enzyme.MC29 = regexpr(pattern='(oxidoreductase|hydrolase|dehydrogenase|synthase|reductase|transferase|methyltransferase|oxidase|synthetase|monooxygenase|isomerase|dehydratase|decarboxylase|deaminase|O\\-methyltransferase|transaminase|hydratase|acetyltransferase|N\\-acetyltransferase|dioxygenase|aminotransferase|O\\-acyltransferase|esterase|N\\-methyltransferase|acyltransferase|aldolase|thiolesterase|O\\-acetyltransferase|cyclase)', text = as.character(as.vector(anno$Note)), perl=T, ignore.case=T) > 0 
    rownames(anno) = anno$ID
    anno.all[[species]] = anno
  }  
  
  # attach single gene conservation levels to each species' annotation
  pdf('conservation.pdf')
  
  for (i in 1:length(desc$gff)){
    # i = 11
    species = desc$species[i]
    gs = ortho[,species] 
    ID2i = 1:length(anno.all[[species]]);
    names(ID2i) = anno.all[[species]]$ID
    rownames(anno.all[[species]]) = anno.all[[species]]$ID
    anno.df = as.data.frame(anno.all[[species]])
    is.SM = regexpr(pattern='secondary metab', text = as.character(as.vector(anno$Note)), perl=T)>0 # 20140519
    
    for (i in 1:length(anno.df)){
      if (class(anno.df[[i]])!='integer')
        anno.df[[i]] = unlist2(anno.df[[i]])
    }
    for (g in 1:length(gs)){
      if (gs[[g]]=='')
        next
      idx = ID2i[strsplit(gs[[g]], split='\\|')[[1]]]
      idx = idx[!is.na(idx)]
      anno.df[idx,'conservativity']= conservativity[g]
    }
    anno.df[is.na(anno.df[,'conservativity']),'conservativity'] = 0;
    hist(anno.df$conservativity, main=species)
    hist.by(anno.df$conservativity,is.SM, by.name='NP gene',hist=T,binwidth=0.05, xlab='conservativity')
    anno.all[[species]] = anno.df
  }  
  
  # neighor context conservation
  
  # output enzyme clusters with the conservation information
  
  # analyze conservation information of known NPGCs
  
  # construct ortholog gene-adjacency graph, nodes --- ortholog groups, edge --- adjacency in each species
  
  # discovery of in-dels & HGT
  
}

ortholog.graph <- function(){
  
}

csv2tab <- function(csv.file){
  csv = read.csv(csv.file, header=T)
  write.table(csv, sub('\\.csv$', '\\.tab',csv.file), quote=F, sep='\t', row.names=F)
}

summerize.cluster <- function(s2d, gene.range = NULL, extra.nt=2500, all.proteins = NULL, anno = NULL, gff.file = NULL, #"/Users/yongli/Universe/data/NPgenome/Aspergillus/A_nidulans_FGSC_A4_current_features.gff",
                              bam.file=NULL,unmapped.bam.file=NULL, 
                              swiss.db = 'swissprot', swiss.fasta.file = paste('/Users/yongli/Universe/data/blastdb/',swiss.db, '.fasta', sep=''), genome.db=NULL,
                              DNA.fasta.file='/Users/yongli/Universe/data/NPgenome/Aspergillus/A_nidulans_FGSC_A4_current_chromosomes.fasta', 
                              prot.fasta.file = 'A_nidulans_FGSC_A4_current_orf_trans_all.fasta', 
                              iprscan.tab.file='A_nidulans_FGSC_A4_iprscan.out.txt',
                              prot.seq = read.fasta(prot.fasta.file, type='AA'), 
                              ipr.anno = iprscan.flat(iprscan.tab.file), multialn.method = 'mafft',
                              intergenic.evidence = T, # 20160805
                              tag=deparse(substitute(s2d)), no.top.hits = 5, RORA.iteration=2, RORA.topOnly=T, plotLogo=F, species=NULL, use.RNAseq.exonpart=T, 
                              minintronlen = 15, maxintronLen = 5000, # to be consistent with tophat2 parameters used for A fumigatus: RNAseq_mapping.sh
                              do.blast=T, do.tblastx=F, extract.bam=!is.null(bam.file),
                              gene.definition = 'gene',
                              geneID2cdsID = geneID2cdsID,# 20141003
                              version = 3 # 20160818, version 3 add start, stop codon, and intergenic region evidences, it assigns different priorities to evidences of different confidence levels
                              # blastp.xml.file = '',# aln.cerevisiae.file='', aln.albicans.file='', aln.NCrassa.file='',  aln.fischeri.file='', aln.self.file='', aln.sp.file='',
){# root = '/Users/yongli/Universe/write/Project_Current/t.NPbioinformatics/Nidulans.SlidingWindow/Annotation'){
  # 20140801
  # 20141003: add automated blast search
  # 20141114: add anno, gene all.genes, add DNA blast searches, and RNA-seq bam reads extraction
  # all.proteins: protein IDs with gene names as names
  # setwd(root)
  # 20141215: modify RORA to RORA.iteration
  #   dir.create(gene.range[3])
  #   setwd(gene.range[3])
  # 20160805: add intergenic.evidence
  
  if (RORA.iteration>0)
    system('VBoxHeadless -s BioLinux7&')
    # system('VBoxManage controlvm BioLinux7 poweroff')
  if (!is.null(all.proteins)){
    genes = intersect2(sub("transcript:","",all.proteins), rownames(prot.seq))
  }else{
    genes = intersect(s2d$all.gene.anno, rownames(prot.seq))  # 20141212 
  }
  if (!is.null(s2d)){
    CS = round(s2d$CS[genes],2); express = round(s2d$mu[genes],2)
  }else{
    CS <- express <- rep('|', times = length(genes)); names(CS) <- names(express) <- genes;
  }
  
  # gene2protein = learn.gff.ID.mapping(genes, )
  if ('anno' %in% colnames(prot.seq)){
    desc = prot.seq[genes, 'anno'];  
  }else{
    desc = '';
  }
  desc = sub('^.*amino acids\\) (.+)$', '\\1', desc)
  # extract gff sub set
  locs = geneRanges2ntRanges(anno, gene.range, extra.nt) 
  gff.sub.file = paste(tag, '.gff', sep='')
  gff.subset(gff.file, locs, out.file=gff.sub.file, format = 'gff3', shift = F)
  # get blast result
  fasta.file = paste(tag, '.fasta', sep='');
  prot.seq = prot.seq[genes, ]; prot.seq[, 'name'] = names(genes); rownames(prot.seq) = names(genes)
  write.fasta(prot.seq, fasta.file)
  if (is.null(ipr.anno) || length(ipr.anno)==0 || is.na(ipr.anno)){
    ipr.anno = vector('character', length = length(genes))
    names(ipr.anno) = genes
  }
  out = cbind('protein seq' = prot.seq[names(genes), 'seq'], name=prot.seq[names(genes), 'name'], length = sapply(prot.seq[names(genes), 'seq'],nchar), Existing.Anno = desc, domains = ipr.anno[genes], CS = CS, express = express)
  
  if (do.blast){
    if (RORA.iteration >0){
      no.top.hits1 = 100L          
    }else{
      no.top.hits1 = no.top.hits
    }
    no.top.hits2 = 100000L
    Sys.setenv(BLASTDB='/Users/yongli/Universe/data/blastdb/')
    blastp.asn.file = paste(tag, '_swissprot.asn', sep='');
    blastp.xml.file = paste(tag, '_swissprot.xml', sep='');
    blastp.hitList = paste(tag, '_swissprot.list', sep='')
    blastp.hitTab = paste(tag, '_swissprot.tab', sep='')
    blastp.hitFasta = paste(tag,'_blastp.fasta', sep='')    
    if (!file.exists(blastp.xml.file) | RORA.iteration>0){
      cat('Blast seaerch', tag)
      system(paste('blastp -query', fasta.file, '-num_threads 6 -db ', swiss.db, '-outfmt 11 -out', blastp.asn.file, '-evalue 1 -max_target_seqs ', no.top.hits1))
      system(paste('blast_formatter -archive', blastp.asn.file, '-outfmt 5 -out', blastp.xml.file, '-max_target_seqs ', no.top.hits))      
    }
    system(paste('blast_formatter -archive', blastp.asn.file, '-outfmt \'6 sseqid\'  -out', blastp.hitList, '-max_target_seqs ', no.top.hits1))
    # system(paste('formatting.pl -idlist ', blastp.hitList, ' -input ', swiss.fasta.file, ' -o ', blastp.hitFasta, sep=''))
    # system(paste('blast_formatter -archive', blastp.asn.file, '-outfmt \'6 qseqid qframe qstart qend evalue qseq sseq sseqid sstart send\'  -out', blastp.hitTab, '-max_target_seqs ', no.top.hits2))
    
    # DNA-blast search   ## get DNA sequence and perform DNA blast
    DNA.seq = getDNA.subseq(DNA.fasta.file, locs = locs)
    DNA.sub.fasta.file = paste(tag, 'DNA_subseq.fasta', sep='')
    blastx.asn.file = paste(tag, 'DNA_subseq.swissprot.asn', sep='')
    blastx.xml.file = paste(tag, 'DNA_subseq.swissprot.xml', sep='')
    blastx.hitList = paste(tag, 'DNA_subseq.swissprot.list', sep='')
    blastx.hitTab = paste(tag, 'DNA_subseq.swissprot.tab', sep='')
    blastx.hitFasta = paste(tag,'_blastx.fasta', sep='')    
    blast.hitList = paste(tag, 'match.list', sep='')
    blast.hitFasta = paste(tag,'match.fasta', sep='')    
    
    blast.AspG.asn.file = paste(tag, 'DNA_subseq.AspGenomes.asn', sep='')
    blast.AspG.xml.file = paste(tag, 'DNA_subseq.AspGenomes.xml', sep='')
    export(DNA.seq, con = DNA.sub.fasta.file, format = 'fasta')
    
    cat('Genome Blast seaerch', tag, '  ',swiss.db,'\n') 
    Sys.setenv(BLASTDB='/Users/yongli/Universe/data/blastdb/')
    if (!file.exists(blastx.xml.file) & !is.null(swiss.db) & RORA.iteration>0){
      #if (!file.exists(blastx.asn.file))
      system(paste('blastx -query', DNA.sub.fasta.file, '-db', swiss.db,  '-num_threads 6 -outfmt 11 -out', blastx.asn.file, '-evalue 1 -max_target_seqs ', no.top.hits2))
      #if (!file.exists(blastx.xml.file))
      system(paste('blast_formatter -archive', blastx.asn.file, '-outfmt 5 -out', blastx.xml.file, '-max_target_seqs ', no.top.hits2))
      # swissSeq = read.fasta(fasta.files = swiss.fasta.file, type = 'AA')
      # swissHits = unique(read.table(blastx.hitList, header=F, as.is=T)[,1]); 
      # system(paste('cdbfasta ',swiss.fasta.file))      
      system(paste('blast_formatter -archive', blastx.asn.file, '-outfmt \'6 sseqid\'  -out', blastx.hitList, '-max_target_seqs ', no.top.hits2))
      # system(paste('formatting.pl -idlist ', blastx.hitList, ' -input ', swiss.fasta.file, ' -o ', blastx.hitFasta, sep=''))
      # system(paste('blast_formatter -archive', blastx.asn.file, '-outfmt \'6 qseqid qframe qstart qend evalue qseq sseq sseqid sstart send\'  -out', blastx.hitTab, '-max_target_seqs ', no.top.hits2))
      # system(paste('rm ', blastx.asn.file)) # 20141125
      system(paste('cat ', blastx.hitList, ' ', blastp.hitList, ' > ', blast.hitList, sep=''))
      system(paste('formatting.pl -idlist ', blast.hitList, ' -input ', swiss.fasta.file, ' -o ', blast.hitFasta, sep=''))
    }
    
    if (do.tblastx){
      cat('Genome Blast seaerch', tag, '  genome.db\n') 
      if (!file.exists(blast.AspG.xml.file) & !is.null(genome.db)){
        #if (!file.exists(blast.AspG.asn.file))
        system(paste('tblastx -query', DNA.sub.fasta.file, '-db', genome.db, '-num_threads 6 -outfmt 11 -out', blast.AspG.asn.file, '-evalue 1 -max_target_seqs ', no.top.hits2))
        #if (!file.exists(blast.AspG.xml.file))
        system(paste('blast_formatter -archive', blast.AspG.asn.file, '-outfmt 5 -out', tblastx.hitList, '-max_target_seqs ', no.top.hits2))      
        # system(paste('rm ', blast.AspG.asn.file)) # 20141125
        
      }      
    }
  }
  
  ######### match predicted proteins with existing protein models and renames predicted genes
  score.file = ''
  for  (iteration in seq2(1,RORA.iteration,1)){
    # extract bam
    if (!is.null(bam.file) & extract.bam){
      bam.out.file = bam.extract.shift(bam.file, locs, tag, shift=F)
    }
    
    # protein evidences: http://bioinf.uni-greifswald.de/bioinf/wiki/pmwiki.php?n=Augustus.IncorporateProteins
    # using protein profiles --proteinprofile=filename: http://bioinf.uni-greifswald.de/augustus/binaries/tutorial/ppx.html
    # RNA-tophat evidence: http://bioinf.uni-greifswald.de/bioinf/wiki/pmwiki.php?n=IncorporatingRNAseq.Tophat
    # EST hits: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1810548/
    # ESTs or assembled RNAseq transcripts: http://bioinf.uni-greifswald.de/bioinf/wiki/pmwiki.php?n=Augustus.IncorporateESTs  
    # Conservation:
    ## old approach AGRIPPA, http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1810548/ 
    if (iteration == 1){
      
      ######### proMap scoring orignal proteins
      pMap = blast2profile.PP(blast.asn.file = blastp.asn.file, 
                              query.gff.file = gff.sub.file, 
                              query.faa.file = fasta.file,
                              DNA.fasta.file = DNA.fasta.file,
                              geneID2cdsID=geneID2cdsID,
                              multialn.method = multialn.method, plot.width = 50, plotLogo =plotLogo,
                              db = swiss.fasta.file)
      nSeq.file = paste('pMap_nSeq_', tag, '_', '', '.faa', sep='')
      nSeq.naive.file = paste('pMap_nSeqNaive_', tag, '_', '', '.faa', sep='')
      cSeq.long.file = paste('pMap_cSeqLong_', tag, '_', '', '.faa', sep='')
      proMap.hint.file = paste(tag, '_proMap2hints.gff', sep='')
      proMap.hint.all.file = paste(tag, '_proMap2hints_all.gff', sep='')
      
      score.file = write.proMap(pMap, nSeq.file = nSeq.file, nSeq.naive.file = nSeq.naive.file, cSeq.long.file = cSeq.long.file, tag = tag, append=F, iteration = '')
      proMap2hints(pMap, gff.file = gff.sub.file, out.file = proMap.hint.all.file, geneID2cdsID=geneID2cdsID, append=F, version = version)
      
      proMap.Mosaichint.file = paste(tag, 'exonerate.nSeq.hints.gff', sep='')
      system(paste('exonerate --model protein2genome --showtargetgff T -q ', nSeq.file, ' -t ', DNA.sub.fasta.file, ' > exonerate.nSeq.out', sep=''))
      system(paste('exonerate2hints.pl --minintronlen=',minintronlen, ' --maxintronlen=', maxintronLen, ' --in=exonerate.nSeq.out --source=P --out=exonerate.hints', sep=''))
      gff.unshift('exonerate.hints', proMap.Mosaichint.file)
      system(paste('cat ', proMap.Mosaichint.file, ' >> ', proMap.hint.all.file, sep=''))
      
      #   proMap.hintNaive.file = paste(tag, 'exonerate.nSeqNaive.hints.gff', sep='')
      #   system(paste('exonerate --model protein2genome --showtargetgff T -q ', nSeq.naive.file, ' -t ', DNA.sub.fasta.file, ' > exonerate.nSeqNaive.out', sep=''))
      #   system(paste('exonerate2hints.pl --minintronlen=',minintronlen, ' --maxintronlen=', maxintronLen, ' --in=exonerate.nSeqNaive.out --source=M --out=exonerate.hints', sep=''))
      #   gff.unshift('exonerate.hints', proMap.hintNaive.file)
      #   
      #   proMap.hintcSeqLong.file = paste(tag, 'exonerate.cSeqLong.hints.gff', sep='')
      #   system(paste('exonerate --model protein2genome --showtargetgff T -q ', cSeq.long.file, ' -t ', DNA.sub.fasta.file, ' > exonerate.cSeqLong.out', sep=''))
      #   system(paste('exonerate2hints.pl --minintronlen=',minintronlen, ' --maxintronlen=', maxintronLen, ' --in=exonerate.cSeqLong.out --source=M --out=exonerate.hints', sep=''))
      #   gff.unshift('exonerate.hints', proMap.hintcSeqLong.file)
      # system(paste('/usr/local/bin/python ~/Universe/code/python/gff2other.py -g', gff.sub.file, '-f', DNA.sub.fasta.file, '-k genbank -s _A_',sep=' '))
      chrseq.file = extra.chr(DNA.fasta.file, locs[,1]) # extract chromosome sequence
      #out.folder = sub('/Users/yongli/Universe/', 'Universe/', getwd())
      out.folder = sub('/Users/yongli/', 'yongli/', getwd())
      
      ############ repeatmasker evidences
      repeatmasker.hint.file = paste(tag, 'rpeatmasker.gff', sep='')
      system(paste('repeatmasker ', DNA.sub.fasta.file, sep=''))
      system(paste('cat ', DNA.sub.fasta.file, '.out | tail -n +3 | perl -ne \'chomp; next if (/^\\s*$/); s/^\\s+//;  @t = split(/\\s+/);print $t[4]."\\t"."repmask\\tnonexonpart\\t".$t[5]."\\t".$t[6]."\\t0\\t.\\t.\\tsrc=RM\\n";\' | sort -n -k 1,1 > ', repeatmasker.hint.file, sep=''))
      gff.unshift(repeatmasker.hint.file, gff.out.file = repeatmasker.hint.file)
      
      ############ denovo predictions
      auguNovoAll.file = paste(tag, '_augoNovoAll.gff', sep='')  
      auguNovoTop.file = paste(tag, '_augoNovoTop.gff', sep='')
      system(paste('sshpass -p abcd ssh fuga@192.168.56.110 \'cd ', out.folder, '; augustus --stopCodonExcludedFromCDS=false --sample=300 --predictionStart=', locs[,2], ' --predictionEnd=', locs[,3], ' --singlestrand=false --species=', species, ' --extrinsicCfgFile=~/',out.folder,'/extrinsic.cfg --alternatives-from-evidence=true  --alternatives-from-sampling=true --minexonintronprob=0.08 --minmeanexonintronprob=0.3 --maxtracks=100 --gff3=on --genemodel=complete ', chrseq.file, ' > ', auguNovoAll.file, '\'', sep=''))
      system(paste('sshpass -p abcd ssh fuga@192.168.56.110 \'cd ', out.folder, '; augustus --stopCodonExcludedFromCDS=false --sample=300 --predictionStart=', locs[,2], ' --predictionEnd=', locs[,3], ' --singlestrand=false --species=', species, ' --extrinsicCfgFile=~/', out.folder, '/extrinsic.cfg --alternatives-from-evidence=false  --alternatives-from-sampling=false --minexonintronprob=0.08 --minmeanexonintronprob=0.3 --maxtracks=100 --gff3=on --genemodel=complete ', chrseq.file, ' > ', auguNovoTop.file, '\'', sep=''))      
      gff.match(gff.file = auguNovoAll.file, gff.reference = gff.sub.file, tag = '', match.by = gene.definition) #, geneID2cdsID=geneID2cdsID); # change gene names 
      gff.match(gff.file = auguNovoTop.file, gff.reference = gff.sub.file, tag = '', match.by = gene.definition) #, geneID2cdsID=geneID2cdsID); # change gene names 
      
      ######### proMap scoring of de novo proteins
      if (RORA.topOnly){
        auguNovo.file = auguNovoTop.file
      }else{
        auguNovo.file = auguNovoAll.file
      }
      system(paste('getAnnoFasta.pl --seqfile=', chrseq.file, ' ', auguNovo.file, sep=''))
      cds.seq.file = sub('.gff', '.codingseq', auguNovo.file);
      fasta.file = sub('.gff', '.aa', auguNovo.file);
      translate.fasta(CDS.file=cds.seq.file, pep.file=fasta.file); #
      blastp.asn.file = sub('.gff', '.asn', auguNovo.file);
      cat('Blast seaerch of Augustus de novo proteins')
      system(paste('blastp -query', fasta.file, '-num_threads 6 -db ', swiss.db,' -outfmt 11 -out', blastp.asn.file, '-evalue 1 -max_target_seqs ', no.top.hits1))
      system(paste('blast_formatter -archive', blastp.asn.file, '-outfmt \'6 sseqid\'  -out', blastp.hitList, '-max_target_seqs ', no.top.hits1))
      system(paste('cat ', blastp.hitList, ' >> ', blast.hitList, sep='')) # add but not replacing 
      
      pMap = blast2profile.PP(blast.asn.file = blastp.asn.file, 
                              query.gff.file = auguNovo.file,
                              query.faa.file = fasta.file,
                              DNA.fasta.file = DNA.fasta.file,
                              geneID2cdsID=function(x){paste(x, '.cds', sep='')},
                              # geneID2cdsID=geneID2cdsID,
                              multialn.method = multialn.method, plot.width = 50, plotLogo =plotLogo, iteration = paste('', sep=''),
                              db = swiss.fasta.file)
      nSeq.file = paste('pMap_nSeq_', tag, '_', 'deNovo', '.faa', sep='')
      nSeq.naive.file = paste('pMap_nSeqNaive_', tag, '_', 'deNovo', '.faa', sep='')
      cSeq.long.file = paste('pMap_cSeqLong_', tag, '_', 'deNovo', '.faa', sep='')
      score.file = write.proMap(pMap, nSeq.file = nSeq.file, nSeq.naive.file = nSeq.naive.file, cSeq.long.file = cSeq.long.file, tag = tag, append=T, iteration = 'deNovo')
      proMap2hints(pMap, gff.file = auguNovo.file, out.file = proMap.hint.file, geneID2cdsID=function(x){paste(x, '.cds', sep='')}, version = version)
      system(paste('cat ', proMap.hint.file, ' >> ', proMap.hint.all.file, sep=''))
      
      proMap.Mosaichint.file = paste(tag, 'exonerate.nSeq.hints','deNovo', '.gff', sep='')
      
      system(paste('exonerate --model protein2genome --showtargetgff T -q ', nSeq.file, ' -t ', DNA.sub.fasta.file, ' > exonerate.nSeq.out', sep=''))
      system(paste('exonerate2hints.pl --minintronlen=',minintronlen, ' --maxintronlen=', maxintronLen, ' --in=exonerate.nSeq.out --source=M --out=exonerate.hints', sep=''))
      gff.unshift('exonerate.hints', proMap.Mosaichint.file)
      system(paste('cat ', proMap.Mosaichint.file, ' >> ', proMap.hint.all.file, sep=''))
    }
    
    # cdbfasta protein.fa
    # cdbfasta genome.fa
    # cat cAfu3g01400_Afu3g01480tblastn.out | allBlastMatches_ncbi-blast.pl > tblastn.matches
    # cat tblastn.matches | perl -e 'while(<>){split; if ($q eq $_[0]){$t .= "\t$_[1]"} else {print "$q$t\n"; $t="\t$_[1]";$q=$_[0];}} print "$q$t\n";' > tblastn.matchlists
    
    ##################### protein hints by exonerate, using all hits proteins
    system(paste('formatting.pl -idlist ', blast.hitList, ' -input ', swiss.fasta.file, ' -o ', blast.hitFasta, sep='')) # prepare all blastx and blastp hits for next rounds of exonerate, 20141216
    exonerate.hint.file = paste(tag, 'exonerate.hints.gff', sep='')
    system(paste('exonerate --model protein2genome --showtargetgff T -q ', blast.hitFasta, ' -t ', DNA.sub.fasta.file, ' > exonerate.out', sep=''))
    system(paste('exonerate2hints.pl --minintronlen=', minintronlen, ' --maxintronlen=', maxintronLen, ' --in=exonerate.out --source=P --out=exonerate.hints', sep=''))
    gff.unshift('exonerate.hints', exonerate.hint.file) 
    
    
    ##################### iteration 1
    all.hints.file = paste(tag, 'all.hints', iteration, sep='')
    if (!is.null(bam.file) & extract.bam){
      if (iteration > 1){
        system(paste('cat ', proMap.hint.all.file, ' ', exonerate.hint.file, ' ', auguHintsAll.file, ' ', auguNovoAll.file, ' all.hints | grep -e \'\tintron\t\' > newIntrons.gff', sep=''))
        system(paste('cat newIntrons.gff | perl -ne \'@array = split(/\\t/, $_);print "$array[0]:$array[3]-$array[4]\\n";\'| sort -u > introns.lst', sep=''))
        system(paste('/Users/yongli/Universe/ubuntu_bin/augustus-3.0.3/scripts/intron2exex.pl --flank=100 --introns=introns.lst --seq=', chrseq.file, ' --exex=exex.fa --map=map.psl', sep=''))
        system(paste('bowtie2-build exex.fa ', tag, '_exex1', sep = ''))
        
        # remapping using unmapped reads
        unmapped.fastq.file = sub('bam', 'fastq', unmapped.bam.file)
        if (!file.exists(unmapped.fastq.file)){
          system(paste('samtools bam2fq -O ', unmapped.bam.file, ' > ', unmapped.fastq.file, sep=''))
        }
        system(paste('bowtie2 --no-unal -p 6 -x ', tag, '_exex1 -U',  unmapped.fastq.file, ' -S bowtieNewIntrons.sam', sep='')) # mapping to the junctions, keep only mapped reads
        # system('samtools view -S -F 4 bowtieNewIntrons1.sam > bowtieNewIntrons.F.sam') # filter to keep mapped reads
        system('samMap.pl bowtieNewIntrons.sam map.psl 100 > bowtie.global.sam')
        system('cat header.txt bowtie.global.sam > bowtie.global.h.sam')
        system('samtools view -bS -o bowtie.global.h.bam bowtie.global.h.sam')
        # join bam files
        system(paste('samtools merge -f both.bam bowtie.global.h.bam ', bam.out.file, sep=''))
        system(paste('samtools sort -n both.bam tmp', iteration, sep=''))
        # system('bam2hints --intronsonly --in=both.ssf.bam --out=hints.2.gff')
        system(paste('filterBam --uniq --in tmp', iteration, '.bam --out tmp', iteration, '_f.bam', sep=''))
        system(paste('samtools sort tmp', iteration, '_f.bam tmp', iteration, '_sf', sep=''))        
      }else{  
        system(paste('filterBam --uniq --in sorted_', bam.out.file, ' --out tmp', iteration, '_f.bam', sep=''))
        system(paste('samtools view -H tmp', iteration, '_f.bam > header.txt', sep=''))
        system(paste('samtools sort tmp', iteration, '_f.bam tmp', iteration, '_sf', sep=''))
      }
      
      hintIntron.file = paste(tag, '_hints_intron',iteration, '.gff', sep='')  
      RNAseq.hint.file = paste(tag, '_RNAseqhints',iteration, '.gff', sep='')    
      # exon parts hints from RNA-seq
      if (use.RNAseq.exonpart){
        system(paste('bam2hints --trunkSS --remove_redundant --minintronlen=',minintronlen, ' --maxintronlen=', maxintronLen, ' --in=tmp_sf.bam --out=', RNAseq.hint.file, sep=''))    
        #DNA.size.file  = paste(DNA.fasta.file, 'chrSize.tab', sep='')
        #system(paste('faSize -detailed -tab ', DNA.fasta.file, ' > ', DNA.size.file, sep=''))
        #system(paste('bam2bigWig bam tmp2_sf ', DNA.size.file, sep=''))
        ## system('bam2wig bam tmp_sf')
        #system('cat tmp_sf.wig | wig2hints.pl --width=10 --margin=10 --minthresh=2 --minscore=4 --prune=0.1 --src=W --type=ep --UCSC=unstranded.track --radius=4.5 --pri=4 --strand="." > hints.ep.gff')
        #system(paste('cat hints.ep.gff ', RNAseq.hint.file,' > hints.tmp', sep=''))
        #system(paste('mv hints.tmp ', RNAseq.hint.file, sep=''))
      }else{
        system(paste('bam2hints --intronsonly --trunkSS --minintronlen=',minintronlen, ' --maxintronlen=', maxintronLen, ' --in=tmp_sf.bam --out=', RNAseq.hint.file, sep=''))        
      }
      # system(paste('cat ', exonerate.hint.file, hintIntron.file, ' > all.hintsIntron', sep=' '))
      system(paste('cat ', repeatmasker.hint.file, ' ', proMap.hint.all.file, ' ', exonerate.hint.file, ' ', RNAseq.hint.file, ' > ', all.hints.file, sep=''))
    }else{
      system(paste('cat ', repeatmasker.hint.file, ' ', proMap.hint.all.file, ' ', exonerate.hint.file, ' > ', all.hints.file, sep=''))      
    }
    
    
    #################### prediction based on all hints combined
    auguHintsAll.file = paste(tag, '_augoHintsAll',iteration, '.gff', sep='')  
    auguHintsTop.file = paste(tag, '_augoHintsTop', iteration, '.gff', sep='')
    # auguHintsIntron.file = paste(tag, '_augoHintsIntron',iteration, '.gff', sep='')  
    system(paste('sshpass -p abcd ssh fuga@192.168.56.110 \'cd ', out.folder, '; augustus --stopCodonExcludedFromCDS=false --sample=300 --predictionStart=', locs[,2], ' --predictionEnd=', locs[,3], ' --singlestrand=false --species=', species, ' --extrinsicCfgFile=~/', out.folder, '/extrinsic.cfg --alternatives-from-evidence=true  --alternatives-from-sampling=true --minexonintronprob=0.08 --minmeanexonintronprob=0.3 --maxtracks=100 --hintsfile=',all.hints.file, ' --allow_hinted_splicesites=atac --introns=on --gff3=on --genemodel=complete ', chrseq.file, ' > ', auguHintsAll.file, '\'', sep=''))
    # system(paste('sshpass -p abcd ssh fuga@192.168.56.110 \'cd ', out.folder, '; augustus --stopCodonExcludedFromCDS=false--sample=300 --predictionStart=', locs[,2], ' --predictionEnd=', locs[,3], ' --singlestrand=false --species=', species, ' --extrinsicCfgFile=extrinsic.cfg --alternatives-from-evidence=true  --alternatives-from-sampling=true --minexonintronprob=0.05 --minmeanexonintronprob=0.3 --maxtracks=100 --hintsfile=all.hintsIntron', ' --allow_hinted_splicesites=atac --introns=on --gff3=on --genemodel=complete ', chrseq.file, ' > ', auguHintsIntron.file, '\'', sep=''))    
    system(paste('sshpass -p abcd ssh fuga@192.168.56.110 \'cd ', out.folder, '; augustus --stopCodonExcludedFromCDS=false --sample=300 --predictionStart=', locs[,2], ' --predictionEnd=', locs[,3], ' --singlestrand=false --species=', species, ' --extrinsicCfgFile=~/', out.folder, '/extrinsic.cfg --alternatives-from-evidence=false  --alternatives-from-sampling=false --minexonintronprob=0.08 --minmeanexonintronprob=0.3 --maxtracks=100 --hintsfile=',all.hints.file, ' --allow_hinted_splicesites=atac --introns=on --gff3=on --genemodel=complete ', chrseq.file, ' > ', auguHintsTop.file, '\'', sep=''))
    gff.match(gff.file = auguHintsAll.file, gff.reference = gff.sub.file, tag = '', match.by = gene.definition) #, geneID2cdsID=geneID2cdsID); # change gene names 
    gff.match(gff.file = auguHintsTop.file, gff.reference = gff.sub.file, tag = '', match.by = gene.definition) #, geneID2cdsID=geneID2cdsID); # change gene names 
    
    #################### proMap scoring and generate new hints
    if (RORA.topOnly){
      auguHints.file = auguHintsTop.file
    }else{
      auguHints.file = auguHintsAll.file
    }
    system(paste('getAnnoFasta.pl --seqfile=', DNA.fasta.file, ' ', auguHints.file, sep=''))
    cds.seq.file = sub('.gff', '.codingseq', auguHints.file);
    fasta.file = sub('.gff', '.aa', auguHints.file);
    translate.fasta(CDS.file=cds.seq.file, pep.file=fasta.file); #
    blastp.asn.file = sub('.gff', '.asn', auguHints.file);
    cat('Blast seaerch of iteration ',iteration, ' predictions', tag)
    system(paste('blastp -query', fasta.file, '-num_threads 6 -db ', swiss.db, ' -outfmt 11 -out', blastp.asn.file, '-evalue 1 -max_target_seqs ', no.top.hits1))
    system(paste('blast_formatter -archive', blastp.asn.file, '-outfmt \'6 sseqid\'  -out', blastp.hitList, '-max_target_seqs ', no.top.hits1))
    system(paste('cat ', blastp.hitList, ' >> ', blast.hitList, sep='')) # add but not replacing 
    
    pMap = blast2profile.PP(blast.asn.file = blastp.asn.file, 
                            query.gff.file = auguHints.file,
                            query.faa.file = fasta.file,
                            DNA.fasta.file = DNA.fasta.file,
                            geneID2cdsID=function(x){paste(x, '.cds', sep='')},
                            # geneID2cdsID=geneID2cdsID,
                            multialn.method = multialn.method, plot.width = 50, plotLogo =plotLogo, iteration = paste('iter', iteration, sep=''),
                            db =  swiss.fasta.file)
    nSeq.file = paste('pMap_nSeq_', tag, '_', iteration, '.faa', sep='')
    nSeq.naive.file = paste('pMap_nSeqNaive_', tag, '_', iteration, '.faa', sep='')
    cSeq.long.file = paste('pMap_cSeqLong_', tag, '_', iteration, '.faa', sep='')
    score.file = write.proMap(pMap, nSeq.file = nSeq.file, nSeq.naive.file = nSeq.naive.file, cSeq.long.file = cSeq.long.file, tag = tag, append=T, iteration = iteration)
    proMap2hints(pMap, gff.file = auguHints.file, out.file = proMap.hint.file, 
                 log.file = paste('log', tag, '.txt', sep=''), geneID2cdsID=function(x){paste(x, '.cds', sep='')}, version = version)
    system(paste('cat ', proMap.hint.file, ' >> ', proMap.hint.all.file, sep=''))
    
    proMap.Mosaichint.file = paste(tag, 'exonerate.nSeq.hints',iteration, '.gff', sep='')
    proMap.hintNaive.file = paste(tag, 'exonerate.nSeqNaive.hints',iteration, '.gff', sep='')
    proMap.hintcSeqLong.file = paste(tag, 'exonerate.cSeqLong.hints',iteration, '.gff', sep='')
    
    system(paste('exonerate --model protein2genome --showtargetgff T -q ', nSeq.file, ' -t ', DNA.sub.fasta.file, ' > exonerate.nSeq.out', sep=''))
    system(paste('exonerate2hints.pl --minintronlen=',minintronlen, ' --maxintronlen=', maxintronLen, ' --in=exonerate.nSeq.out --source=M --out=exonerate.hints', sep=''))
    gff.unshift('exonerate.hints', proMap.Mosaichint.file)
    system(paste('cat ', proMap.Mosaichint.file, ' >> ', proMap.hint.all.file, sep=''))
    #     system(paste('exonerate --model protein2genome --showtargetgff T -q ', nSeq.naive.file, ' -t ', DNA.sub.fasta.file, ' > exonerate.nSeqNaive.out', sep=''))
    #     system(paste('exonerate2hints.pl --minintronlen=',minintronlen, ' --maxintronlen=', maxintronLen, ' --in=exonerate.nSeqNaive.out --source=M --out=exonerate.hints', sep=''))
    #     gff.unshift('exonerate.hints', proMap.hintNaive.file)
    #     
    #     system(paste('exonerate --model protein2genome --showtargetgff T -q ', cSeq.long.file, ' -t ', DNA.sub.fasta.file, ' > exonerate.cSeqLong.out', sep=''))
    #     system(paste('exonerate2hints.pl --minintronlen=',minintronlen, ' --maxintronlen=', maxintronLen, ' --in=exonerate.cSeqLong.out --source=M --out=exonerate.hints', sep=''))
    #     gff.unshift('exonerate.hints', proMap.hintcSeqLong.file)
  }
  
  if (score.file != '' && file.exists(score.file))
    select.CDS(score.file) # chose the top 2 gene models
  # blastp.xml.file =  sub('.gff', '.xml', auguHints2top.file);
  # system(paste('blast_formatter -archive', blastp.asn.file, '-outfmt 5 -out', blastp.xml.file, '-max_target_seqs ', no.top.hits))
  # system(paste('sshpass -p abcd ssh fuga@192.168.56.110 \'poweroff\''))
  ## read alignment results  
  if (blastp.xml.file != '' & file.exists(blastp.xml.file)){
    blast.out = blast.xml.parse(blast.xml = blastp.xml.file, no.top.hits = no.top.hits)
    top.species = sub('^.+\\[([^\\[\\]]+)\\].?$','\\1', blast.out$query$Top_nonself_Hit_def, perl=T)
    names(top.species) = rownames(blast.out$query)
    out = cbind(out, Top_nonself_Hit_species = top.species[names(genes)], blast.out$query[names(genes), c(8,16,19)])  
    add.names = c('cluster',  'RNA-seq reads', 'No.introns',  'intron anno by RNA-seq',  'inron anno by orthologs',
                  'conclusion', 'new.protein.seq');
  }else{
    add.names = c('Top_nonself_Hit_species','Top_nonself_Hit_accession', 'Top_nonself_Hit_identity.percent', 'top.3.hits',
                  'cluster',  'RNA-seq reads', 'No.introns',  'intron anno by RNA-seq',  'inron anno by orthologs',
                  'conclusion', 'new.protein.seq');
  }
  
  out = cbind(out, matrix('|', nrow = nrow(out), 
                          ncol=length(add.names), dimnames = list(names(genes), add.names)))
  out = as.matrix(out)
  out[is.na(out)] = '|';
  invisible(out)
}

select.CDS <- function(score.file = 'pMapcJL1All.xls'){
  # select canidate CDSs based on score and rank and high light them in sorted pMap output file
  # Yong Fuga Li, 20150108
  
  require('xlsx')
  s = read.table(score.file, header = T, sep = '\t')
  # score = as.numeric(s$pHMM.score); score[is.na(score)] = 0
  CDS.rank = regexpr.match('^[^\\.]+(?:\\.t([^\\.]+))?(?:\\.[^\\.]+)?$', s$X, perl=T)[,1]
  CDS.rank[CDS.rank==''] = '1';
  CDS.rank = as.numeric(CDS.rank)
  iter = as.character(s$iteration); iter[iter=='deNovo'] = 0; iter[iter==''] = -1; iter = as.numeric(iter)
  s = sort.by(s, cbind(CDS.rank, iter))
  gene.ID = regexpr.match('^([^\\.]+)(?:\\..+)?$', s$X, perl=T)[,1]
  s = sort.by(s, gene.ID)
  if (1){ # 20160613 - sort by gene locations
    gene.ID = regexpr.match('^([^\\.]+)(?:\\..+)?$', s$X, perl=T)[,1]
    gID.first = which.first.by(gene.ID) 
    gfrom = as.numeric(as.character(s$from))[gID.first]; # 20160613 - sort by gene locations
    names(gfrom) = gene.ID[gID.first];
    gfrom = gfrom[gene.ID];
    s = sort.by(s, gfrom)
  }
  
  score = as.numeric(s$pHMM.score); score[is.na(score)] = 0
  gene.ID = regexpr.match('^([^\\.]+)(?:\\..+)?$', s$X, perl=T)[,1]
  CDS.rank = regexpr.match('^[^\\.]+(?:\\.t([^\\.]+))?(?:\\.[^\\.]+)?$', s$X, perl=T)[,1]
  CDS.rank[CDS.rank==''] = '1'
  CDS.rank = as.numeric(CDS.rank)
  iter = as.character(s$iteration); iter[iter=='deNovo'] = 0; iter[iter==''] = -1; iter = as.numeric(iter)
  
  # suggest the two best candidates 1) 2nd iter top -- iff score + 3 > top score, otherwise chose the one with top score;
  deta.score = 3;
  i.candidate = unlist.dupenames(by(1:nrow(s), INDICES = gene.ID, FUN = function(x){i = which((iter[x] %in% c(max(iter[x]), -1)) & CDS.rank[x]==1); 
  score1 = score[x];
  s.max = max(score1)
  i.max = which(score1==s.max);
  i = unique(c(i[score1[i] + deta.score >= s.max], i.max))
  return(i.select=x[i])})) # candidates: with max score, or close to max && is original/last iteration of CDS prediction
  i.select = unlist.dupenames(by(1:nrow(s), INDICES = gene.ID, FUN = function(x){max.iter = max(iter[x]);
  i = which((iter[x] %in% c(max.iter, -1)) & CDS.rank[x]==1); 
  score1 = score[x];
  s.max = max(score1)
  i.max = which(score1==s.max);
  i = unique(c(i[score1[i] + deta.score >= s.max], i.max));
  i.max.iter = intersect(i, which(iter[x]==max.iter)) # final augustus prediction
  i.original = intersect(i, which(iter[x]==-1)) # original gene
  if (length(i.original)>0)
    i = i.original
  else if (length(i.max.iter)>0)
    i = i.max.iter
  else
    i = i[which.min((CDS.rank[x])[i[(iter[x])[i]==max((iter[x])[i])]])] # the one with highest augustus probability in the last iteration in the candidates
  return(i.select=(x[i])[1])})) # candidates: with max score, or close to max && is original/last iteration of CDS prediction
  i.max = which.max.tie.by(score, by = gene.ID)
  score.improvement.select = sapply(1:length(i.select), function(x){i = (iter == -1 & gene.ID == names(i.select[x])); if (!any(i)) d = 'NA' else d = round(score[i.select[x]] - unique(score[i],3))}) 
  score.improvement.max = sapply(1:length(i.max), function(x){i = (iter == -1 & gene.ID == names(i.max[x])); if (!any(i)) d = 'NA' else d = round(score[i.max[x]] - unique(score[i],3))}) 
  score.improvement.candidate = sapply(1:length(i.candidate), function(x){i = (iter == -1 & gene.ID == names(i.candidate[x])); if (!any(i)) d = 'NA' else d = round(score[i.candidate[x]] - unique(score[i],3))}) 
  ### write.xlsx file with proper highlighting
  # high light 1) the promising and 2) the max scored CDS
  score.file = sub(pattern = '^.+\\/([^\\/]+)$', replacement = '\\1', score.file)
  out.file = paste('pretty_', sub('\\.[^\\.]*$','.xlsx', score.file), sep='')
  names(s)[1] = 'ID'
  s$candidates = ''; s$max.scored = ''; s$selected = ''; 
  s$score.improvement = ''; # 20150123
  s$candidates[i.candidate] = 'Yes'; s$max.scored[i.candidate] = 'Yes'; s$selected[i.select] = 'Yes'; s$score.improvement[i.select] = score.improvement.select; s$score.improvement[i.max] = score.improvement.max; s$score.improvement[i.candidate] = score.improvement.candidate
  s$pHMM.Evalue[s$pHMM.Evalue == Inf] = 1
  write.xlsx2(s, out.file, row.names = F, showNA = F)
  xlsx.color(xlsx.file = out.file, include.header=T, FUN.select = function(x){y = matrix(T, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x)); y},
             font=list(color = NULL, heightInPoints=12, name='Calibri', isItalic=F, isBold=F, isStrikeout=F, underline=NULL), 
             out.file = out.file, na.strings='|')  # change global style
  CDS.groups = cbind(unlist(as.list(by(1:length(gene.ID), gene.ID, FUN = min))), unlist(as.list(by(1:length(gene.ID), gene.ID, FUN = max)))); # group CDS by gene IDs
  xlsx.color(xlsx.file = out.file, row.groups = CDS.groups, out.file=out.file) # frame to indicate the genes
  xlsx.color(xlsx.file = out.file, FUN.select = function(x){y = matrix(F, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x));
  y[i.candidate, ] = T
  return(y)}, fill.color = 'green', out.file = out.file, na.strings='|')  
  xlsx.color(xlsx.file = out.file, FUN.select = function(x){y = matrix(F, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x));
  y[i.max, ] = T
  return(y)}, 
  font=list(color = NULL, isItalic=T, isBold=F, isStrikeout=F, underline=NULL), out.file = out.file, na.strings='|')
  xlsx.color(xlsx.file = out.file, FUN.select = function(x){y = matrix(F, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x));
  y[i.select, ] = T
  return(y)}, 
  font=list(color = 'red', isItalic=F, isBold=T, isStrikeout=F, underline=NULL), out.file = out.file, na.strings='|')
  
  write(paste('green: candidates - either with max phmm score or with scores no less than [max score] - ', deta.score, '
              italic: max scored
              bold red: selected - among the candidates,
              if the original CDS model is no less than [max score] - 3, the original CDS is selected,
              otherwise if the best model in the last iteration of augustus is the is no less than [max score] - 3, that is selected,
              otherwise the highest iteration and highest augustus probability model among the candidates is selected', sep=''), file = 'readme_select_CDS.txt')
}

table.merge <- function(files = m[i,3], extra.columns = m[i,2:3], idx.keep = 1:17, file.format = 'xlsx', out.file = 'KU0011.merge.xls'){
  # Yong Fuga Li, 20151011
  # idx.keep, columns in the orginal table to keep
  require(xlsx)
  require(gdata)
  dat.all = c()
  for (i in 1:length(files)){
    if (file.format == 'xlsx'){
      x = read.xlsx2(files[i], 1)      
    }else{
      x = read.table(files[i], sep='\t', header = T)            
    }
    if (is.null(idx.keep))
      idx.keep = 1:ncol(x)
    extra = repmat(extra.columns[i,,drop=F]);
    colnames(extra) = colnames(extra.columns)
    dat.all = rbind(dat.all, cbind(x[,idx.keep], extra))
  }
  dat.all1 = as.matrix(dat.all)
  rownames(dat.all1) = dat.all1[,1]
  write.table(dat.all1[,2:ncol(dat.all1)], col.names=NA, sep = '\t', file = out.file)
  
}

select.CDS.multiModel <- function(re.cluster.model = 'pretty_pMapc(\\w{2}\\d{4})_(.*).xlsx$', # provide cluster ID and model species
                                  all.files = list.files(pattern = re.cluster.model)) {
  # select canidate CDSs based on score and rank and high light them in sorted pMap output file
  # Yong Fuga Li, 20151011
  # 20160804, add all.files
  
  require('xlsx')
  ## group files
  m = regexpr.match(txt = all.files, pat = re.cluster.model)
  m = cbind(m, all.files)
  m = m[m[,1]!='',,drop=F]
  rownames(m)= m[,1]
  colnames(m) = c('cluster', 'model', 'file')
  # by(m[,2], m[,1], as.character)
  # by(all.files, m[,1], as.character)
  
  ## merge files
  idx = by(1:nrow(m), m[,1], identity)
  for (i in names(idx)){
    # print(i+1)
    out.file = paste('pMapc', i, '.merge.xlx', sep='')
    table.merge(files = m[idx[[i]],3], extra.columns = m[idx[[i]],2:3, drop=F], file.format = 'xlsx', out.file = out.file)
    select.CDS(score.file = out.file)
  }
}

blast.xml.parse <- function(blast.xml = 'AN8127-Alignment.xml', no.top.hits = 3,
                            query.self.min.identity = 0.95,  query.self.len.diff = 10, 
                            query.species = 'Aspergillus nidulans FGSC A4'){
  # ref: https://stat.ethz.ch/pipermail/bioc-sig-sequencing/2010-September/001580.html
  # http://rstudio-pubs-static.s3.amazonaws.com/12097_1352791b169f423f910d93222a4c2d85.html
  # ref: blastSequences
  # query.self.min.identity, query.self.len.diff, and query.species are used to define which hits are itential to query, and hence excluded from the output
  # YF Li, 20140803
  require(XML)
  result <- xmlTreeParse(blast.xml, useInternalNodes=TRUE)
  aa = xmlToList(result)
  
  # read: iteration level
  tags.query = c('Iteration_iter-num', 'Iteration_query-ID', 'Iteration_query-def', 'Iteration_query-len')
  query = list();
  for (t in tags.query)
    query[[t]] = unlist(xpathApply(result, paste("//",t,sep=''), xmlValue))
  query = data.frame(query,stringsAsFactors = F)
  tags.query = make.names(tags.query)
  query[[1]] = as.numeric(query[[1]]);   query[[4]] = as.numeric(query[[4]]); 
  rownames(query) = query$Iteration_query.def
  
  # hit level
  tags.hit = c('Hit_num', 'Hit_id', 'Hit_def', 'Hit_accession', 'Hit_len')
  hit = list();
  for (t in tags.hit)
    hit[[t]] = unlist(xpathApply(result, paste("//",t,sep=''), xmlValue))
  hit = as.data.frame(hit,stringsAsFactors = F)
  tags.hit = make.names(tags.hit)
  hit[[1]] = as.numeric(hit[[1]]);   hit[[5]] = as.numeric(hit[[5]]); 
  
  # hsp level
  tags.hsp = c('Hsp_num', 'Hsp_bit-score','Hsp_score','Hsp_evalue','Hsp_query-from','Hsp_query-to', 
               'Hsp_hit-from','Hsp_hit-to','Hsp_query-frame','Hsp_hit-frame','Hsp_identity','Hsp_positive',
               'Hsp_gaps','Hsp_align-len')
  hsp = list();
  for (t in tags.hsp)
    hsp[[t]] = unlist(xpathApply(result, paste("//",t,sep=''), xmlValue))
  hsp = as.data.frame(hsp,stringsAsFactors = F)
  tags.hsp = make.names(tags.hsp)
  hsp = data.frame(lapply(hsp, FUN = as.numeric))
  
  # expand query
  No.hits = xpathApply(result, "//Iteration_hits", function(x)sum(names(xmlApply(x,xmlName))=='Hit')) # get the number of hits from each qurey
  # No.hits = xpathApply(result, "//Iteration_hits", xmlSize) # xmlSize has a bug that returns 1 for empty node
  No.hits = sapply(No.hits, unlist)
  iQuery4hit = rep(rownames(query), No.hits)
  #   i = hit[['Hit_num']];   i = i[which(c(i[2:length(i)],1) == 1)] # this will messup with a query has 0 hits
  #   iQuery4hit = rep((1:length(i)), i)
  query.hit = as.data.frame(apply(query, MARGIN = 2, function(x) rep(x, times = No.hits)),stringsAsFactors = F)
  query.hit[[1]] = as.numeric(query.hit[[1]]);   query.hit[[4]] = as.numeric(query.hit[[4]]); 
  
  # expand query and hit
  No.hsps = xpathApply(result, "//Hit_hsps", function(x)xmlSize(x)) # get the number of hits from each qurey
  No.hsps = sapply(No.hsps, unlist)
  i = as.numeric(hsp[['Hsp_num']]);   i = i[which(c(i[2:length(i)],1) == 1)]
  iQuery4hsp = rep(iQuery4hit, i);   iHit4hsp = rep((1:length(i)), i)  # all(iQuery4hit[iHit4hsp] == iQuery4hsp) == TRUE
  query.hsp = as.data.frame(apply(query.hit, MARGIN = 2, function(x) rep(x, times = i)),stringsAsFactors = F)
  hit.hsp = as.data.frame(apply(hit, MARGIN = 2, function(x) rep(x, times = i)),stringsAsFactors = F)
  query.hsp[[1]] = as.numeric(query.hsp[[1]]);   query.hsp[[4]] = as.numeric(query.hsp[[4]]); 
  hit.hsp[[1]] = as.numeric(hit.hsp[[1]]);   hit.hsp[[5]] = as.numeric(hit.hsp[[5]]); 
  
  # summerize hsp to hits
  hit.extra =  apply(hsp[,c(2,3,11:14)],MARGIN = 2,FUN = function(x)by(x, INDICES = list(hit=iHit4hsp), FUN = sum))
  colnames(hit.extra) = sub('Hsp', 'Hit', colnames(hit.extra))
  hit = cbind(hit, hit.extra)
  hit$Hit_identity.percent = hit$Hit_identity/hit$Hit_align.len
  hit$Hit_positive.percent = hit$Hit_positive/hit$Hit_align.len
  hit$Hit_gaps.percent = hit$Hit_gaps/hit$Hit_align.len
  
  # sumerize hits to query 
  i.self = (abs(hit$Hit_len - query.hit$Iteration_query.len)< query.self.len.diff | regexpr(query.species, hit$Hit_def)>0)&
    (hit$Hit_identity.percent > query.self.min.identity)
  i.top.hit = by(hit[!i.self,4:13], INDICES = factor(iQuery4hit[!i.self], levels =rownames(query)), 
                 FUN=function(x){i = which.max(x$Hit_bit.score); 
                 return(as.numeric(rownames(x[i, ])))},simplify=T)
  i.top.hit = sapply(i.top.hit, unlist) # note that i.top.hits is the rownames, hence no need to be shifted by i.self
  
  hit$Hit_identity.percent = paste(round(hit$Hit_identity.percent*100,1), '%', sep='')
  hit$Hit_positive.percent = paste(round(hit$Hit_positive.percent*100,1), '%', sep='')
  hit$Hit_gaps.percent = paste(round(hit$Hit_gaps.percent*100,1), '%', sep='')
  query.extra = hit[i.top.hit,]
  rownames(query.extra) = names(i.top.hit)
  
  top.N.hits = by(hit[!i.self,], INDICES = factor(iQuery4hit[!i.self], levels =rownames(query)), 
                  FUN=function(x){i = which.max.n(x$Hit_bit.score, no.top.hits); x = x[i,];
                  top.N.hits = paste(x$Hit_identity.percent, ' | ', x$Hit_accession,' | ', x$Hit_def, sep='', collapse = ' // ')
                  return(top.N.hits)},simplify=T)
  top.N.hits = sapply(top.N.hits, unlist)
  query.extra[[paste('top', no.top.hits, 'hits', sep='.')]] =  top.N.hits
  
  colnames(query.extra) = sub('Hit', 'Top_nonself_Hit', colnames(query.extra))
  query = cbind(query, query.extra[rownames(query),])
  
  blast.out = list(query=query, hit = hit, hsp=hsp, query.hit = query.hit, query.hsp=query.hsp, hit.hsp=hit.hsp,
                   iQuery4hit=iQuery4hit, iQuery4hsp=iQuery4hsp, iHit4hsp=iHit4hsp) 
  return(blast.out)
}

blast.filter <- function(bl, Evalue = 0.1){
  # 20140916, YF Li
  i = bl$hsp$Hsp_evalue < Evalue
  i.hit  = unique(bl$iHit4hsp[i])
  i.query = unique(bl$iQuery4hsp[i])
  bl$query = bl$query[i.query,]
  bl$hit = bl$hit[i.hit,]
  bl$hsp = bl$hsp[i,]
  bl$query.hit = bl$query.hit[i.hit,]
  bl$query.hsp = bl$query.hsp[i,]
  bl$hit.hsp = bl$hit.hsp[i,]
  bl$iQuery4hit = bl$iQuery4hit[i.hit]
  bl$iQuery4hsp = bl$iQuery4hsp[i]
  bl$iHit4hsp = bl$iHit4hsp[i]
  return(bl)
}

blast2coverage <- function(blast.xml, Evalue=0.1, type = c('count', 'bit/length')){
  # 20141125
  # Yong Fuga Li
  
  bl = blast.xml.parse(blast.xml = blast.xml)
  bl.f = blast.filter(bl, Evalue = Evalue)
}

blast2profile <- function(blast.xml='cUUp0_S281DNA_subseq.swissprot.xml', no.top.hits = 10E10, Evalue=0.1, type = c('count', 'bit/length')){
  # 20141125
  # output: 
  #   profile = list(query.DNA, query.AA = matrix(6, n.AA), matches=matrix(nrow=21, ncol=n.aa), insertions = list(locations, inserts.aligned))
  # Yong Fuga Li
  
  require(XML)
  result <- xmlTreeParse(blast.xml, useInternalNodes=TRUE)
  aa = xmlToList(result)
  
  # read: iteration level
  tags.query = c('Iteration_iter-num', 'Iteration_query-ID', 'Iteration_query-def', 'Iteration_query-len')
  query = list();
  for (t in tags.query)
    query[[t]] = unlist(xpathApply(result, paste("//",t,sep=''), xmlValue))
  query = data.frame(query,stringsAsFactors = F)
  tags.query = make.names(tags.query)
  query[[1]] = as.numeric(query[[1]]);   query[[4]] = as.numeric(query[[4]]); 
  rownames(query) = query$Iteration_query.def
  
  # hit level
  tags.hit = c('Hit_num', 'Hit_id', 'Hit_def', 'Hit_accession', 'Hit_len')
  hit = list();
  for (t in tags.hit)
    hit[[t]] = unlist(xpathApply(result, paste("//",t,sep=''), xmlValue))
  hit = as.data.frame(hit,stringsAsFactors = F)
  tags.hit = make.names(tags.hit)
  hit[[1]] = as.numeric(hit[[1]]);   hit[[5]] = as.numeric(hit[[5]]); 
  
  # hsp level
  tags.hsp = c('Hsp_num', 'Hsp_bit-score','Hsp_score','Hsp_evalue','Hsp_query-from','Hsp_query-to', 
               'Hsp_hit-from','Hsp_hit-to','Hsp_query-frame','Hsp_hit-frame','Hsp_identity','Hsp_positive',
               'Hsp_gaps','Hsp_align-len')
  hsp = list();
  for (t in tags.hsp)
    hsp[[t]] = unlist(xpathApply(result, paste("//",t,sep=''), xmlValue))
  hsp = as.data.frame(hsp,stringsAsFactors = F)
  tags.hsp = make.names(tags.hsp)
  hsp = data.frame(lapply(hsp, FUN = as.numeric))
  
  # expand query
  No.hits = xpathApply(result, "//Iteration_hits", function(x)sum(names(xmlApply(x,xmlName))=='Hit')) # get the number of hits from each qurey
  # No.hits = xpathApply(result, "//Iteration_hits", xmlSize) # xmlSize has a bug that returns 1 for empty node
  No.hits = sapply(No.hits, unlist)
  iQuery4hit = rep(rownames(query), No.hits)
  #   i = hit[['Hit_num']];   i = i[which(c(i[2:length(i)],1) == 1)] # this will messup with a query has 0 hits
  #   iQuery4hit = rep((1:length(i)), i)
  query.hit = as.data.frame(apply(query, MARGIN = 2, function(x) rep(x, times = No.hits)),stringsAsFactors = F)
  query.hit[[1]] = as.numeric(query.hit[[1]]);   query.hit[[4]] = as.numeric(query.hit[[4]]); 
  
  # expand query and hit
  No.hsps = xpathApply(result, "//Hit_hsps", function(x)xmlSize(x)) # get the number of hits from each qurey
  No.hsps = sapply(No.hsps, unlist)
  i = as.numeric(hsp[['Hsp_num']]);   i = i[which(c(i[2:length(i)],1) == 1)]
  iQuery4hsp = rep(iQuery4hit, i);   iHit4hsp = rep((1:length(i)), i)  # all(iQuery4hit[iHit4hsp] == iQuery4hsp) == TRUE
  query.hsp = as.data.frame(apply(query.hit, MARGIN = 2, function(x) rep(x, times = i)),stringsAsFactors = F)
  hit.hsp = as.data.frame(apply(hit, MARGIN = 2, function(x) rep(x, times = i)),stringsAsFactors = F)
  query.hsp[[1]] = as.numeric(query.hsp[[1]]);   query.hsp[[4]] = as.numeric(query.hsp[[4]]); 
  hit.hsp[[1]] = as.numeric(hit.hsp[[1]]);   hit.hsp[[5]] = as.numeric(hit.hsp[[5]]); 
  
  # summerize hsp to hits
  hit.extra =  apply(hsp[,c(2,3,11:14)],MARGIN = 2,FUN = function(x)by(x, INDICES = list(hit=iHit4hsp), FUN = sum))
  colnames(hit.extra) = sub('Hsp', 'Hit', colnames(hit.extra))
  hit = cbind(hit, hit.extra)
  hit$Hit_identity.percent = hit$Hit_identity/hit$Hit_align.len
  hit$Hit_positive.percent = hit$Hit_positive/hit$Hit_align.len
  hit$Hit_gaps.percent = hit$Hit_gaps/hit$Hit_align.len
  
  # sumerize hits to query 
  i.self = (abs(hit$Hit_len - query.hit$Iteration_query.len)< query.self.len.diff | regexpr(query.species, hit$Hit_def)>0)&
    (hit$Hit_identity.percent > query.self.min.identity)
  i.top.hit = by(hit[!i.self,4:13], INDICES = factor(iQuery4hit[!i.self], levels =rownames(query)), 
                 FUN=function(x){i = which.max(x$Hit_bit.score); 
                 return(as.numeric(rownames(x[i, ])))},simplify=T)
  i.top.hit = sapply(i.top.hit, unlist) # note that i.top.hits is the rownames, hence no need to be shifted by i.self
  
  hit$Hit_identity.percent = paste(round(hit$Hit_identity.percent*100,1), '%', sep='')
  hit$Hit_positive.percent = paste(round(hit$Hit_positive.percent*100,1), '%', sep='')
  hit$Hit_gaps.percent = paste(round(hit$Hit_gaps.percent*100,1), '%', sep='')
  query.extra = hit[i.top.hit,]
  rownames(query.extra) = names(i.top.hit)
  
  top.N.hits = by(hit[!i.self,], INDICES = factor(iQuery4hit[!i.self], levels =rownames(query)), 
                  FUN=function(x){i = which.max.n(x$Hit_bit.score, no.top.hits); x = x[i,];
                  top.N.hits = paste(x$Hit_identity.percent, ' | ', x$Hit_accession,' | ', x$Hit_def, sep='', collapse = ' // ')
                  return(top.N.hits)},simplify=T)
  top.N.hits = sapply(top.N.hits, unlist)
  query.extra[[paste('top', no.top.hits, 'hits', sep='.')]] =  top.N.hits
  
  colnames(query.extra) = sub('Hit', 'Top_nonself_Hit', colnames(query.extra))
  query = cbind(query, query.extra[rownames(query),])
  
  blast.out = list(query=query, hit = hit, hsp=hsp, query.hit = query.hit, query.hsp=query.hsp, hit.hsp=hit.hsp,
                   iQuery4hit=iQuery4hit, iQuery4hsp=iQuery4hsp, iHit4hsp=iHit4hsp) 
  return(blast.out)
}

cluster.deepAnno <- function(gene.ranges = NULL, gff.file=NULL,
                             geMat=NULL, ica.spatial=NULL, prot.fasta.file=NULL, iprscan.tab.file = NULL, iprscan.table.file = NULL, 
                             bam.file = NULL, unmapped.bam.file=NULL, EST.db = NULL, swiss.db = c('swissprot', 'fungiRefSwiss70'),
                             swiss.fasta.file = paste('/Users/yongli/Universe/data/blastdb/', swiss.db, '.fasta', sep=''),
                             DNA.fasta.file=NULL, genome.db = NULL,
                             in.ID.type = NULL,pat.prot.ID='',
                             prot.seq = read.fasta(prot.fasta.file, pattern = pat.prot.ID, type='AA'), 
                             ipr.anno = iprscan.flat(iprscan.table.file), 
                             gene.definition = c('gene', 'transcript', 'mRNA', 'CDS'), 
                             out.file = NULL,append=F, proteinID = 'ID',
                             geneID2cdsID = NULL,
                             # geneID2cdsID = function(x){paste(x, '-P', sep='')},
                             extra.genes = 0, RORA.iteration=2, RORA.topOnly =T, multialn.method = 'mafft', # mafft is better based on hmmsearch of predicted genes against pHMM models build from blast hits
                             plotLogo=T, species=NULL, do.blast=T, do.tblastx=F, center.method = 'median', score.type = 'R', median.substraction = F, cor.method = 'pearson',
                             n.cluster.per.file = 70, start.from=1,end.to=NULL,
                             extra.nt = 2500, remove.intermediate.files = T, 
                             s2d=NULL, # precomputed s2d
                             version = 3 # 20160818, version 3 add start, stop codon, and intergenic region evidences, it assigns different priorities to evidences of different confidence levels
){
  # deep annotation of multiple predicted clusters
  # YF Li 20140723-0803
  # 20141003: add extra.genes
  # 20141010: add genome file, modify to work without expression data
  # 20141112: modify to work without expression data - add bam.file and EST.db
  require(xlsx)
  require('XLConnect')
  require('Biostrings')
  require(gplots)
  require(rtracklayer)
  root.dir = getwd()
  system('cp /Users/yongli/Universe/write/Project_Current/9.O.NPbioinformatics/extrinsic.cfg ./')
  swiss.db = match.arg(swiss.db); # 20160611
  
  if (is.null(iprscan.table.file))
    iprscan.table.file = iprscan.tab.file  
  gene.definition = match.arg(gene.definition)
  tag = sub('^(.*)\\.xls.*','\\1', out.file)
  if (!is.data.frame(gene.ranges) & !is.matrix(gene.ranges)){
    gene.ranges = matrix(gene.ranges, 1, ncol = length(gene.ranges))
  }
  if (is.null(end.to)){
    end.to=nrow(gene.ranges)  
  }
  if (ncol(gene.ranges)==2){
    gene.ranges = cbind(gene.ranges, paste(gene.ranges[,1], gene.ranges[,2], sep = '_'))
  }else if (ncol(gene.ranges)!=3){
    stop('gene.ranges need 3 or 2 columns')
  }
  
  anno = import.gff(gff.file) # 20160502
  if (is.null(geneID2cdsID)){
    m = learn.gff.ID.mapping(unlist.multi(anno@elementMetadata@listData$ID), 
                             parent = unlist.multi(anno@elementMetadata@listData[[which(tolower(colnames(anno@elementMetadata))=='parent')]]), 
                             node.type =  as.character(anno@elementMetadata@listData$type))
    geneID2cdsID = m[[paste(gene.definition, '2CDS', sep='')]]
  }
  
  if(!is.null(geMat)&!is.null(gff.file)){
    ica.spatial = express.clustering(gff.file, geMat)    
    anno = ica.spatial$anno;    
    is.expressed = !is.na(match(ica.spatial$anno$ID, rownames(ica.spatial$S)))
    names(is.expressed) = ica.spatial$anno$ID
  }else if(!is.null(ica.spatial)){
    anno = ica.spatial$anno;        
    is.expressed = !is.na(match(ica.spatial$anno$ID, rownames(ica.spatial$S)))
    names(is.expressed) = ica.spatial$anno$ID  
  }else if(!is.null(gff.file)){
    # gff.format = sub('^.*\\.([^\\.]*$)', '\\1', gff.file)
    # anno = tryCatch(read.gff3(gff.file, format='gff3'), error = function(e){read.gff3(gff.file, format='gff')}, finally = NULL)
    #anno = read.gff3(gff.file, format=gff.format)
    idx.gene = (anno$type==gene.definition)
    anno = anno[idx.gene, ]
    anno = sort.intervals(anno)   
    # anno$ID = sub('transcript:', '', anno$ID)
    is.expressed = vector('logical', length = length(anno)) | T;
    names(is.expressed) = anno$ID;
  }else{
    stop('Provide ica.spatial or gff.file')
  }
  
  if (!is.null(in.ID.type)){ # 20151001
    gene.ranges[,1:2] = anno$ID[sort(match(gene.ranges[,1:2],  sub('\\.\\d$', '', anno@elementMetadata[,in.ID.type])))]
  }
  
  gene.ranges.original = gene.ranges
  gene.index.ranges = matrix(0, nrow=nrow(gene.ranges), ncol = 2)
  K = 0; # window size
  for (i in start.from:end.to){
    # tweek gene.ranges to find nearest expressed gene
    i.1 <- i.10 <- match(gene.ranges[i,1], names(is.expressed));
    while(i.1 > 1 && (as.character(anno@seqnames[i.1-1])==as.character(anno@seqnames[i.10])) && (!is.expressed[i.1] | i.10 - i.1 < extra.genes))
      i.1 = i.1 - 1;
    gene.ranges[i,1] = names(is.expressed)[i.1]
    i.2 <- i.20 <- match(gene.ranges[i,2], names(is.expressed));
    if (as.character(anno@seqnames[i.10])!=as.character(anno@seqnames[i.20])){
      cat(gene.ranges[i,1], gene.ranges[i,2])
      stop('clusters spans to chromosomes')  
    }
    while(i.2 < length(is.expressed) && (as.character(anno@seqnames[i.2+1])==as.character(anno@seqnames[i.20])) && (!is.expressed[i.2] | i.2 - i.20 < extra.genes))
      i.2 = i.2 + 1;
    gene.ranges[i,2] = names(is.expressed)[i.2]
    # gene.index.ranges = rbind(gene.index.ranges, c(i.1, i.2))
    gene.index.ranges[i,] = c(i.1, i.2) # 20141125
    K = max(K, i.2-i.1+1)
  }
  
  
  if (is.null(s2d) & !is.null(ica.spatial)){
    s2d = list()
    ica.spatial = ica.spatial.prep(ica.spatial, K= K, center.method=center.method,
                                   score.type = score.type, median.substraction=median.substraction, do.plot = F) # precompute scores
    for (i in start.from:end.to){
      s2d[[gene.ranges[i,1]]]= score.spatial.cluster.2d(ica.spatial, gene.range=gene.ranges[i,1:2], cor.method = cor.method)
    }    
  }
  
  # pdf(paste('spatial.cluster.', 'cm_', center.method,'.st_',score.type,'.ms_',median.substraction, '.cor_', cor.method, '.pdf', sep=''),20,12)
  if (!is.null(s2d)){
    for (i in start.from:end.to){
      # cID = paste(gene.ranges[i,1],gene.ranges[i,3], sep='_')
      # cID = paste(gene.ranges.original[i,3],gene.ranges.original[i,1], sep='_')
      cID = paste('c', gene.ranges.original[i,3], sep='')
      fig.file = paste(cID, '.png', sep='')
      png(fig.file, 20,12,units = 'in', res=60); 
      plot.spatial.cluster.2d(s2d[[gene.ranges[i,1]]], tag=paste(gene.ranges.original[i,3],gene.ranges.original[i,1],gene.ranges.original[i,2], sep='_'))
      # plot.spatial.cluster.2d(AN8127, tag=paste(gene.ranges[i,3],gene.ranges[i,1], sep=': '), no.fdr = T)
      dev.off()
      # plot.spatial.cluster.2d(s2d[[i]], tag=cID)  
    }     
    # dev.off()    
  }
  
  for (i in start.from:end.to){
    # cID = paste(gene.ranges[i,1],gene.ranges[i,3], sep='_')
    file.index = floor((i-1)/n.cluster.per.file)
    #cID = paste('c', gene.ranges.original[i,3],'_',gene.ranges.original[i,1], sep='')
    cID = paste('c', gene.ranges.original[i,3], sep='')
    all.genes = names(is.expressed)[gene.index.ranges[i,1]:gene.index.ranges[i,2]] 
    all.proteins = anno@elementMetadata[gene.index.ranges[i,1]:gene.index.ranges[i,2],proteinID]; names(all.proteins) = all.genes;
    tab = summerize.cluster(s2d[[gene.ranges[i,1]]], gene.range =gene.ranges[i,], extra.nt=extra.nt, all.proteins = all.proteins, 
                            swiss.db = swiss.db, swiss.fasta.file = swiss.fasta.file,
                            genome.db=genome.db, anno=anno, gff.file=gff.file, prot.seq = prot.seq, bam.file=bam.file,
                            unmapped.bam.file=unmapped.bam.file, RORA.iteration=RORA.iteration, RORA.topOnly = RORA.topOnly, multialn.method = multialn.method, species=species,plotLogo=plotLogo,
                            DNA.fasta.file=DNA.fasta.file, ipr.anno = ipr.anno, tag = cID, # paste(gene.ranges[i,1], gene.ranges[i,2], sep='_'),  
                            do.blast=do.blast, do.tblastx=do.tblastx, geneID2cdsID=geneID2cdsID, gene.definition=gene.definition, version = version); # blastp.xml.file = NULL,
    setwd(root.dir)
    CDSs = get.CDS(gene.IDs = rownames(tab),
                   gff.file = gff.file, DNA.fasta.file = DNA.fasta.file,
                   geneID2cdsID=geneID2cdsID) # 20141014: retrieve CDS sequences
    for (j in 1:nrow(CDSs)){
      pep = as.character(translate(DNAString(as.character(CDSs[j,1])), if.fuzzy.codon = 'X'));
      if (gsub('\\*', '', pep) != gsub('\\*', '', tab[j,'protein seq'])){
        warning(paste('CDS translation dose not match protein sequences, likely due to ambiguous nucleotide\n', pep ,'\n', tab[j,'protein seq'], '\n'))
      }
    }
    CDSs = cbind(as.matrix(CDSs),nchar(as.character(CDSs[,1])),
                 'coding NT%' = round(nchar(as.character(CDSs[,1]))/CDSs[,3]*100,1),
                 'average exon size' = round(nchar(as.character(CDSs[,1]))/CDSs[,2],1),
                 'average intron size' = round((CDSs[,3]-nchar(as.character(CDSs[,1])))/(CDSs[,2]-1+1E-10),1))
    colnames(CDSs)[1:5] = c('CDS', 'NO.exon', 'CDS_span(nt)', 'CDS_length(nt)', 'coding percentage')
    tab = cbind(CDSs, tab) # 20141014: add CDS sequences
    tab = mat.fill.row(tab, all.genes, default = '|') # 20141014: add non-protein coding genes back    
    extra.info = cbind(as.matrix(anno@ranges[gene.index.ranges[i,1]:gene.index.ranges[i,2]])[,1:2],
                       as.character(anno@strand[gene.index.ranges[i,1]:gene.index.ranges[i,2]]));
    colnames(extra.info) = c('start', 'width', 'strand')
    cluster.boundary = c('', 'Boundary')[1+!is.na(match(rownames(tab), gene.ranges.original))]
    tab = cbind(cluster.boundary = cluster.boundary, extra.info, tab)
    out.file.1 = sub('\\.([^\\.]*$)', c('.xlsx', paste('\\.',file.index, '\\.\\1', sep=''))[(file.index>0)+1],  out.file);
    write.xlsx(tab, out.file.1, sheetName = substr(cID, 1,31),col.names = T, row.names=T,showNA = F,
               append = (((i-1)%%n.cluster.per.file != 0)|append))  
    
    wb <- loadWorkbook(out.file.1, create = TRUE)
    fig.file = paste(cID, '.png', sep='')
    cID = substr(cID, 1,31)
    if (!is.null(s2d)){ # 20141112{
      createName(wb, name = cID, formula = paste(cID, "!$A$", nrow(tab)+6, ':', "$S$", nrow(tab)+63, sep=''))
      addImage(wb, filename = fig.file, name = cID, originalSize = F)  
    }
    setColumnWidth(wb,sheet=cID,column=5,width=256*30)
    # setColumnWidth(wb,sheet=cID,column=10,width=256*30)
    # setColumnWidth(wb,sheet=cID,column=11,width=256*30)
    setColumnWidth(wb,sheet=cID,column=10+5,width=256*30)
    setColumnWidth(wb,sheet=cID,column=10+6,width=256*30)
    setColumnWidth(wb,sheet=cID,column=10+9,width=256*20)
    setColumnWidth(wb,sheet=cID,column=10+12,width=256*25)
    saveWorkbook(wb)
  }
  if (remove.intermediate.files){
    system(paste('rm ', cID, '*', sep = ''))
    fs = setdiff(dir(pattern = paste('.*', cID, sep='')),dir(pattern = paste('pMap', cID, sep='')))
    for (f in fs){
      system(paste('rm ', f))
    }
    system('rm tmp*')
    system('rm exonerate*')
    system('rm hits*')
    system('rm bowtie*')
    
    for (f in c('both.bam', 'exex.fa', 'introns.lst','map.psl',
                'newIntrons.gff', 'header.txt',
                'chr.fasta', 'extrinsic.cfg')){
      system(paste('rm ', f))
    }
  }
  return(s2d)
}

augustus.species.ID <- function(augu.file = '/Users/yongli/Universe/ubuntu_bin/augustus-3.0.3/README.TXT'){
  txt = read.table(augu.file, sep = '\t', quote = '', header = F)
  m = regexpr.match('^([^\\s\\(\\)]*) *\\| (.*)$', txt[regexpr(pattern = '^[^\\)]*\\| .*$', text = txt[,1])>0,1]) # identifiers in () are older versions, ignored
  rownames(m) = m[,2]
  m = m[m[,2] != 'species',1]
  return(m)
}

RORA.pipeline <- function(root = '/Users/yongli/Universe/write/Project_Current/t.NPbioinformatics/KU2015',
                          cluster.info.file = '150925_UbiA_like_terpene_clusters_JGI_boundaries.xlsx',
                          from.id.type = 'protein_id', from.gene.definition = 'CDS',
                          use.multiple.model.species = F, RORA.iteration = 2, # swiss.db = c('swissprot', 'fungiRefSwiss70'),
                          skip.existing =F, plotLogo=F, extra.genes = 1,
                          i.start = 1, i.all = NULL,
                          simplify.model.species = T, tag = '', 
                          GenBank.Only = F,
                          version = 2,# 20160818, version 3, add start, stop codon, and intergenic region evidences, it assigns different priorities to evidences of different confidence levels
                          ...){ 
  # i.start - starting clustering or i.all - index of clusters to included
  # when i.all is provided, it is used instead of i.start
  # 20151006-20151009
  # 20160611, add tag, RORA.iteration, swiss.db, ...
  # YF Li
  require(xlsx)
  
  cluster.info = read.xlsx2(cluster.info.file, 1, stringsAsFactors = F)
  cluster.info = cluster.info[rowSums(cluster.info!='')>0,]
  if (!is.null(i.all))
    cluster.info = cluster.info[i.all,,drop=F]
  else
    cluster.info = cluster.info[i.start:nrow(cluster.info),,drop=F]
  modelSpecies = select.ModelSpecies(unique(cluster.info$species), simplify = simplify.model.species)
  write.table(modelSpecies[[2]], file = paste('model_species_selection',i.start,'.tsv', sep=''), sep='\t')
  cluster.info$First.Protein = sub('\\.\\d$', '', cluster.info$First.Protein)
  cluster.info$Last.Protein = sub('\\.\\d$', '', cluster.info$Last.Protein)
  clusters = cbind(gsub(' ', '', cluster.info$First.Protein), gsub(' ', '', cluster.info$Last.Protein), paste(cluster.info$ClusterID, tag, sep=''))
  auguSpeciesID = augustus.species.ID()
  
  log.file = paste('log', i.start, '.txt', sep='')
  cat('\n\n# log info for the gene prediction tasks\n', file = log.file, append = T)
  cat(date(), file = log.file, append = T)
  
  for (i in which(cluster.info$In.House.Genome!='')){
    dat = cluster.info$In.House.Genome[i]
    dat = gsub("'","", dat)
    a = strsplit(dat, '; ')[[1]]
    b = read.table(text = a, header = F, sep = '=', as.is = T, strip.white = T)
    c = b[,2]
    names(c) = b[,1]
    
    model.species.all = auguSpeciesID[names(auguSpeciesID) %in% modelSpecies[[2]][rownames(modelSpecies[[2]])==cluster.info$species[i],1]]
    if (!use.multiple.model.species){
      model.species.all = model.species.all[1]
    }
    for (model.species in model.species.all){
      setwd(root)
      folder.name = paste(cluster.info$ClusterID[i], model.species, sep='_')
      if (skip.existing & file.exists(folder.name))
        next
      
      dir.create(folder.name)
      setwd(folder.name)
      cluster.deepAnno(gene.ranges = clusters[i,], 
                       species = model.species, RORA.iteration = RORA.iteration, 
                       gff.file = c['gff.file'], DNA.fasta.file = c['DNA.file'], iprscan.tab.file=c['iprscan.tab.file'], 
                       gene.definition = c['gene.definition'] , proteinID = c['proteinID'], prot.fasta.file = c['pep.fasta.file'], 
                       extra.genes = extra.genes, 
                       plotLogo=plotLogo, multialn.method = 'muscle', RORA.topOnly=F,
                       # geneID2cdsID = identity, # 20160528 -- learn geneID2cdsID instead
                       out.file = paste(clusters[i,3], '.xlsx', sep=''),...)
      
      xlsx.color.NPGC(paste(clusters[i,3], '.xlsx', sep=''))
      # system(paste('rm -r ', NCBI.genome.tag, '*', sep=''))
    }
  }
  
  for (i in which(cluster.info$JGI.Genome=='' & cluster.info$GenBank.Genome != '')){
    if (!GenBank.Only)
      break
    model.species.all = auguSpeciesID[names(auguSpeciesID) %in% modelSpecies[[2]][rownames(modelSpecies[[2]])==cluster.info$species[i],1]]
    if (!use.multiple.model.species){
      model.species.all = model.species.all[1]
    }
    for (model.species in model.species.all){
      setwd(root)
      folder.name = paste(cluster.info$ClusterID[i], model.species, sep='_')
      if (skip.existing & file.exists(folder.name))
        next
      
      dir.create(folder.name)
      setwd(folder.name)
      NCBI.genome.tag = sub('_genomic', '', cluster.info$GenBank.Genome[i])
      prot.fasta.file = paste(NCBI.genome.tag, '_protein.faa', sep='')
      gff.file = paste(NCBI.genome.tag, '_genomic.gff', sep='')
      DNA.file = paste(NCBI.genome.tag, '_genomic.fna', sep='')
      iprscan.tab.file = NULL
      download.file(paste('ftp://ftp.ncbi.nlm.nih.gov/genomes/all/', NCBI.genome.tag, '/', prot.fasta.file, '.gz', sep=''),destfile = paste(prot.fasta.file, '.gz', sep=''))
      download.file(paste('ftp://ftp.ncbi.nlm.nih.gov/genomes/all/', NCBI.genome.tag, '/', gff.file, '.gz', sep=''),destfile = paste(gff.file, '.gz', sep=''))
      download.file(paste('ftp://ftp.ncbi.nlm.nih.gov/genomes/all/', NCBI.genome.tag, '/', DNA.file, '.gz', sep=''),destfile = paste(DNA.file, '.gz', sep=''))
      system(paste('gzip -d *.gz -f', sep=''))
      
      files.to.remove = list.files()
      files.to.keep = list.files(pattern = '(\\.xlsx|\\.xls)')
      files.to.remove = setdiff(files.to.remove, c(files.to.keep, 'readme_deepAnno.txt', 'readme_select_CDS.txt'))
      
      setwd(root)
      cat(paste('\n#########\n RORA for', folder.name, '\n', sep=''), file = log.file, append=T)
      cat(paste('protein file', prot.fasta.file), file = log.file, append=T)
      cat(paste('DNA file', DNA.file), file = log.file, append=T)
      cat(paste('gff file', gff.file), file = log.file, append=T)
      setwd(folder.name)
      
      
      
      if (!length(prot.fasta.file) | !length(gff.file) |!length(DNA.file))
        warning('GenBank file Missing')
      
      in.info = c(feature = 'CDS', id.type = 'protein_id')
      locs = gff.id.change(gff.file, in.info = in.info, in.ids = clusters[i,1:2],
                           extra.nt = 2500, out.type = 'nt') # change IDs
      DNA.seq = getDNA.subseq(DNA.fasta.file, locs = locs)
      DNA.sub.file = paste(clusters[i,3], tag, '_Genome_subseq.fa', sep='')
      export(DNA.seq, con = DNA.sub.file, format = 'fasta')
      
      system(paste('sshpass -p abcd ssh fuga@192.168.56.110 \'cd ', out.folder, '; augustus --stopCodonExcludedFromCDS=false --sample=300 --predictionStart=', locs[,2], ' --predictionEnd=', locs[,3], ' --singlestrand=false --species=', species, ' --extrinsicCfgFile=~/',out.folder,'/extrinsic.cfg --alternatives-from-evidence=true  --alternatives-from-sampling=true --minexonintronprob=0.08 --minmeanexonintronprob=0.3 --maxtracks=100 --gff3=on --genemodel=complete ', chrseq.file, ' > ', auguNovoAll.file, '\'', sep=''))
      # sshpass -p abcd ssh fuga@192.168.56.110 'cd NPbioinformatics/TangLab; augustus --sample=300 --singlestrand=false --species=aspergillus_nidulans --alternatives-from-evidence=false  --alternatives-from-sampling=false --minexonintronprob=0.08 --minmeanexonintronprob=0.3 --maxtracks=-1 --protein=on --introns=on --start=on --stop=on --cds=on --gff3=on --genemodel=partial Pt_K85_scafSeq.fasta > Pt_K85_scafSeq.augoNovo_mANidulans.gff'
      
      cluster.deepAnno(gene.ranges = clusters[i,], 
                       species = model.species, RORA.iteration = RORA.iteration, 
                       gff.file = gff.file, DNA.fasta.file = DNA.sub.file, iprscan.tab.file=iprscan.tab.file, 
                       gene.definition = 'gene', proteinID = 'protein_id', prot.fasta.file = prot.fasta.file, extra.genes = extra.genes, 
                       plotLogo=plotLogo, multialn.method = 'muscle', RORA.topOnly=F,
                       # geneID2cdsID = identity, 
                       out.file = paste(clusters[i,3], '.xlsx', sep=''),...)
      
      xlsx.color.NPGC(paste(clusters[i,3], '.xlsx', sep=''))
      system(paste('rm -r ', NCBI.genome.tag, '*', sep=''))
      
      # select.CDS(score.file = paste('pMap', clusters[i,3], '.xlsx', sep=''))        
    }
  }
  
  for (i in which(cluster.info$JGI.Genome!='' & cluster.info$GenBank.Genome == '')){
    model.species.all = auguSpeciesID[names(auguSpeciesID) %in% modelSpecies[[2]][rownames(modelSpecies[[2]])==cluster.info$species[i],1]]
    if (!use.multiple.model.species){
      model.species.all = model.species.all[1]
    }
    for (model.species in model.species.all){
      setwd(root)
      folder.name = paste(cluster.info$ClusterID[i], model.species, sep='_')
      if (skip.existing & file.exists(folder.name))
        next
      
      dir.create(folder.name)
      setwd(folder.name)
      download.file('http://genome.jgi.doe.gov/fungi/fungi.info.html', 'JGI_list.html')
      system(paste('downloadJGIassembly.pl -html JGI_list.html -species ', cluster.info$JGI.Genome[i], sep=''))
      system(paste('gzip -d *.gz -f', sep=''))
      files.to.remove = list.files()
      files.to.keep = list.files(pattern = '(\\.xlsx|\\.xls)')
      files.to.remove = setdiff(files.to.remove, c(files.to.keep, 'readme_deepAnno.txt', 'readme_select_CDS.txt'))
      
      iprscan.tab.file = list.files(pattern = paste(cluster.info$JGI.Genome[i], '.*_IPR.tab', sep=''), ignore.case = T)
      if (!length(iprscan.tab.file))
        iprscan.tab.file = list.files(pattern = '.*.domaininfo.*.tab')
      prot.fasta.file = list.files(pattern = paste(cluster.info$JGI.Genome[i], '.*.aa.fasta', sep=''), ignore.case = T)
      if (!length(prot.fasta.file))
        prot.fasta.file = list.files(pattern = '.*.proteins.fasta')
      gff.file = list.files(pattern = paste(cluster.info$JGI.Genome[i], '.*proteins.*FilteredModels1.gff3', sep=''), ignore.case = T)
      if (!length(gff.file))
        gff.file = list.files(pattern = '.*.gff3')
      DNA.file = list.files(pattern = paste(cluster.info$JGI.Genome[i], '.*_Repeatmasked.fasta', sep=''), ignore.case = T)
      if (!length(DNA.file))
        DNA.file = list.files(pattern = '.*masked.*', ignore.case = T)
      
      setwd(root)
      cat(paste('\n#########\n RORA for', folder.name, '\n', sep=''), file = log.file, append=T)
      cat(paste('\niprscan file', iprscan.tab.file), file = log.file, append=T)
      cat(paste('\nprotein file', prot.fasta.file), file = log.file, append=T)
      cat(paste('\nDNA file', DNA.file), file = log.file, append=T)
      cat(paste('\ngff file', gff.file), file = log.file, append=T)
      setwd(folder.name)
      
      iprscan.tab.file = iprscan.tab.file[1]
      prot.fasta.file = prot.fasta.file[1]
      gff.file = gff.file[1]
      DNA.file = DNA.file[1]
      
      if (!length(iprscan.tab.file) | !length(prot.fasta.file) | !length(gff.file) |!length(DNA.file))
        warning('JGI Missing file')
      
      #       cat('Mapping input IDs from GenBank to JGI\n')
      #       hits = best.blast.hits(from.file = from.fasta.file, 
      #                              from.gff.file = from.gff.file,
      #                              to.file = prot.fasta.file, 
      #                              from.IDs = clusters[i,], id.type = from.id.type)

      cluster.deepAnno(gene.ranges = clusters[i,], # c(hits$sseqid[1], hits$sseqid[length(hits$sseqid)], folder.name), # c('gene11134', 'gene11143', 'KU0001'), 
                       species = model.species, RORA.iteration = RORA.iteration, 
                       gff.file = gff.file, DNA.fasta.file = DNA.file, iprscan.tab.file=iprscan.tab.file, 
                       in.ID.type = 'proteinId',
                       gene.definition = 'gene', proteinID = 'proteinId', prot.fasta.file = prot.fasta.file, extra.genes = extra.genes, 
                       plotLogo=plotLogo, multialn.method = 'muscle', RORA.topOnly=F,
                       # geneID2cdsID = function(x){sub('gene_', 'CDS_', x)}, 
                       out.file = paste(clusters[i,3], '.xlsx', sep=''),...)
      xlsx.color.NPGC(paste(clusters[i,3], '.xlsx', sep=''))
      for (f in files.to.remove){
        system(paste('rm -r ', f, sep=''))
      }
      system(paste('rm -r ', cluster.info$JGI.Genome[i], '*', sep=''))
      # select.CDS(score.file = paste('pMap', clusters[i,3], '.xlsx', sep=''))        
    }
  }
  
  for (i in which(cluster.info$JGI.Genome!='' & cluster.info$GenBank.Genome != '')){ # need ID mapping from NCBI to JGI
    model.species.all = auguSpeciesID[names(auguSpeciesID) %in% modelSpecies[[2]][rownames(modelSpecies[[2]])==cluster.info$species[i],1]]
    if (!use.multiple.model.species){
      model.species.all = model.species.all[1]
    }
    for (model.species in model.species.all){
      setwd(root)
      folder.name = paste(cluster.info$ClusterID[i], model.species, sep='_')
      if (skip.existing & file.exists(folder.name))
        next
      
      dir.create(folder.name)
      setwd(folder.name)
      download.file('http://genome.jgi.doe.gov/fungi/fungi.info.html', 'JGI_list.html')
      system(paste('downloadJGIassembly.pl -html JGI_list.html -species ', cluster.info$JGI.Genome[i], sep=''))
      system(paste('gzip -d *.gz -f', sep=''))
      
      NCBI.genome.tag = sub('_genomic', '', cluster.info$GenBank.Genome[i])
      from.fasta.file = paste(NCBI.genome.tag, '_protein.faa', sep='')
      from.gff.file = paste(NCBI.genome.tag, '_genomic.gff', sep='')
      DNA.file = paste(NCBI.genome.tag, '_genomic.fna', sep='')
      iprscan.tab.file = NULL
      download.file(paste('ftp://ftp.ncbi.nlm.nih.gov/genomes/all/', NCBI.genome.tag, '/', from.fasta.file, '.gz', sep=''),destfile = paste(from.fasta.file, '.gz', sep=''))
      download.file(paste('ftp://ftp.ncbi.nlm.nih.gov/genomes/all/', NCBI.genome.tag, '/', from.gff.file, '.gz', sep=''),destfile = paste(from.gff.file, '.gz', sep=''))
      # download.file(paste('ftp://ftp.ncbi.nlm.nih.gov/genomes/all/', NCBI.genome.tag, '/', DNA.file, '.gz', sep=''),destfile = paste(DNA.file, '.gz', sep=''))
      system(paste('gzip -d *.gz -f', sep=''))
      
      files.to.remove = list.files()
      files.to.keep = list.files(pattern = '(\\.xlsx|\\.xls)')
      files.to.remove = setdiff(files.to.remove, c(files.to.keep, 'readme_deepAnno.txt', 'readme_select_CDS.txt'))
      
      iprscan.tab.file = list.files(pattern = paste(cluster.info$JGI.Genome[i], '.*_IPR.tab', sep=''), ignore.case = T)
      if (!length(iprscan.tab.file))
        iprscan.tab.file = list.files(pattern = '.*.domaininfo.*.tab')
      prot.fasta.file = list.files(pattern = paste(cluster.info$JGI.Genome[i], '.*.aa.fasta', sep=''), ignore.case = T)
      if (!length(prot.fasta.file))
        prot.fasta.file = list.files(pattern = '.*.proteins.fasta')
      gff.file = list.files(pattern = paste(cluster.info$JGI.Genome[i], '.*proteins.*FilteredModels1.gff3', sep=''), ignore.case = T)
      if (!length(gff.file))
        gff.file = list.files(pattern = '.*.gff3')
      DNA.file = list.files(pattern = paste(cluster.info$JGI.Genome[i], '.*_Repeatmasked.fasta', sep=''), ignore.case = T)
      if (!length(DNA.file))
        DNA.file = list.files(pattern = '.*masked.*', ignore.case = T)
      
      setwd(root)
      cat(paste('\n#########\n RORA for', folder.name, '\n', sep=''), file = log.file, append=T)
      cat(paste('\niprscan file', iprscan.tab.file), file = log.file, append=T)
      cat(paste('\nprotein file', prot.fasta.file), file = log.file, append=T)
      cat(paste('\nDNA file', DNA.file), file = log.file, append=T)
      cat(paste('\ngff file', gff.file), file = log.file, append=T)
      setwd(folder.name)
      
      iprscan.tab.file = iprscan.tab.file[1]
      prot.fasta.file = prot.fasta.file[1]
      gff.file = gff.file[1]
      DNA.file = DNA.file[1]
      
      if (!length(iprscan.tab.file) | !length(prot.fasta.file) | !length(gff.file) |!length(DNA.file))
        warning('JGI Missing file')
      
      cat('Mapping input IDs from GenBank to JGI\n')
      hits = best.blast.hits(from.file = from.fasta.file, 
                             from.gff.file = from.gff.file,
                             to.file = prot.fasta.file, 
                             from.IDs = clusters[i,], id.type = from.id.type)
      pat.prot.ID = '^.*\\|.*\\|(.+)\\|(.*)$'; # extract protein names from fasta preambles
      hits$sseqid = sub(pat.prot.ID, '\\1',hits$sseqid)
      sseqid.ordered = geneRanges2allGenes(gff.file, hits$sseqid, id.type = 'proteinId', gene.definition = 'gene')
      
      setwd(root)
      cat(paste(paste(hits$qseqid, collapse = ','), 'mapped to', paste(hits$sseqid, collapse = ',')), file = log.file, append=T)        
      if (length(unique(sseqid.ordered))!=length(unique(hits$sseqid))){
        cat('\nGene number changed after mapping: ', file = log.file, append=T)        
        cat(paste('\n!!!',paste(hits$sseqid, collapse = ','), 'mapped to', paste(sseqid.ordered, collapse = ',')), file = log.file, append=T)        
      }else if (sseqid.ordered[1]!=hits$sseqid[1] | sseqid.ordered[length(sseqid.ordered)] != hits$sseqid[length(hits$sseqid)]){
        cat('\nGene order changed after mapping: ', file = log.file, append=T)        
        cat(paste('\n!!!',paste(hits$sseqid, collapse = ','), 'spans to', paste(sseqid.ordered, collapse = ',')), file = log.file, append=T)        
      }
      setwd(folder.name)
      
      cluster.deepAnno(gene.ranges = c(sseqid.ordered[1], sseqid.ordered[length(sseqid.ordered)], folder.name), # c(hits$sseqid[1], hits$sseqid[length(hits$sseqid)], folder.name), # c('gene11134', 'gene11143', 'KU0001'), 
                       species = model.species, RORA.iteration = RORA.iteration, 
                       gff.file = gff.file, DNA.fasta.file = DNA.file, iprscan.tab.file=iprscan.tab.file, 
                       in.ID.type = 'proteinId',pat.prot.ID= pat.prot.ID,
                       gene.definition = 'gene', proteinID = 'proteinId', prot.fasta.file = prot.fasta.file, extra.genes = extra.genes, 
                       plotLogo=plotLogo, multialn.method = 'muscle', RORA.topOnly=F,
                       # geneID2cdsID = function(x){sub('gene_', 'CDS_', x)}, 
                       out.file = paste(clusters[i,3], '.xlsx', sep=''),...)
      xlsx.color.NPGC(paste(clusters[i,3], '.xlsx', sep=''))
      for (f in files.to.remove){
        system(paste('rm -r ', f, sep=''))
        
      }
      system(paste('rm -r ', cluster.info$JGI.Genome[i], '*', sep=''))
      # select.CDS(score.file = paste('pMap', clusters[i,3], '.xlsx', sep=''))        
    }
  }
  setwd(root)
  system('cp */colored_*.xlsx ./')
  # system('cp */pretty_pMap*.xlsx ./')
  # system('cp */*blastp.hits ./')
  if (1){# only process the folders generated in the current execution of the program
    for (i in cluster.info$ClusterID){
      files = dir(recursive = T, pattern = paste('pretty_pMapc', i, '.xlsx', sep=''))
      files = files[regexpr(paste(i,'_.*\\/', sep=''), files)>0]
      select.CDS.multiModel(all.files = files, 
                            re.cluster.model = '^(.*\\.[^_\\/]*)_([^\\/]*)\\/')
    }
  }else{# process all folders under the current directory
    select.CDS.multiModel(all.files = dir(recursive = T, pattern = paste('pretty_pMapc.*.xlsx', sep='')), 
                          re.cluster.model = '^(.*\\.[^_\\/]*)_([^\\/]*)\\/')
  }
  system('rm ./pMap*.xlx')
  
}

deepAnno.clusters <- clusters.deepAnno <- function(cluster.file = '/Users/yongli/Universe/data/NPgenome/Aspergillus/A_nidulans_FGSC_A4_current/cluster.annoCompact.A_nidulans_FGSC_A4_current.MC29e.simu2000.refSimu.chrNonspecific.w20.p0.005.NWindowClusters98.tab',
                                                   gff.file="/Users/yongli/Universe/data/NPgenome/Aspergillus/A_nidulans_FGSC_A4_current_features.gff",
                                                   DNA.fasta.file='/Users/yongli/Universe/data/NPgenome/Aspergillus/A_nidulans_FGSC_A4_current_chromosomes.fasta',
                                                   prot.fasta.file = "A_nidulans_FGSC_A4_current_orf_trans_all.fasta",
                                                   iprscan.tab.file = 'A_nidulans_FGSC_A4_iprscan.out.txt',
                                                   geMat = NULL, 
                                                   gene.definition = c('gene', 'transcript', 'mRNA'), proteinID = 'ID',
                                                   geneID2cdsID=function(x){paste(x, '-P', sep='')},
                                                   ica.spatial=NULL,n.cluster.per.file=70, start.from=1, end.to=NULL,
                                                   out.file = 'nidulans.deepAnno.all.xlsx',
                                                   RORA.iteration = 2, species = 'aspergillus_nidulans', plotLogo=F,multialn.method = 'muscle',RORA.topOnly=T,
                                                   max.dist.merge = -13, # distance cut off for mergeing clusters
                                                   # negative <=> overlaps, 0<=>next to each other
                                                   extra.genes = 5 # add to each side of the cluster
){ # root = '/Users/yongli/Universe/write/Project_Current/t.NPbioinformatics/Nidulans.SlidingWindow/Annotation'){
  # 20141003-1004, YFLi
  # cluster.file = '/Users/yongli/Universe/data/NPgenome/Aspergillus/A_nidulans_FGSC_A4_current/cluster.annoCompact.A_nidulans_FGSC_A4_current.MC29e.simu2000.refSimu.chrNonspecific.w20.p0.005.NWindowClusters98.tab'
  # gff.file="/Users/yongli/Universe/data/NPgenome/Aspergillus/A_nidulans_FGSC_A4_current_features.gff"
  # prot.fasta.file = "/Users/yongli/Universe/write/Project_Current/t.NPbioinformatics/Nidulans.SlidingWindow/Annotation/A_nidulans_FGSC_A4_current_orf_trans_all.fasta"
  # iprscan.tab.file = '/Users/yongli/Universe/write/Project_Current/t.NPbioinformatics/Nidulans.SlidingWindow/Annotation/A_nidulans_FGSC_A4_iprscan.out.txt'
  #   prot.seq = read.fasta(prot.fasta.file, type='AA')
  #   ipr.anno = iprscan.flat(iprscan.tab.file)
  gene.definition = match.arg(gene.definition);
  if(!is.null(geMat)&!is.null(gff.file)){
    ica.spatial = express.clustering(gff.file, geMat)    
    anno = ica.spatial$anno;    
  }else if(!is.null(ica.spatial)){
    anno = ica.spatial$anno;        
  }else if(!is.null(gff.file)){
    gff.format = sub('^.*\\.([^\\.]*$)', '\\1', gff.file)
    # anno = read.gff3(gff.file, format=gff.format)
    anno = import.gff(gff.file) # 20160502
    idx.gene = (anno$type==gene.definition)
    anno = anno[idx.gene, ]
    anno = sort.intervals(anno)    
  }else{
    stop('Provide ica.spatial or gff.file')
  }  
  # get gene annotation and gene orders
  # anno = read.gff3(gff.file, format='gff3')
  idx.gene = (anno$type==gene.definition)
  anno = anno[idx.gene, ]
  anno = sort.intervals(anno)
  n = length(anno)
  # anno$ID = sub('transcript:', '', anno$ID)
  
  
  # get clusters and merge by distances
  cf = read.csv(cluster.file, header = F, sep='\t', as.is = T, comment.char = '#')
  cf = cf[2:nrow(cf),]
  
  n.zeros = sapply(strsplit(cf$V5, '\\.\\.\\.'), length)
  cf$V5[n.zeros>1] = paste(0,'*', n.zeros[n.zeros>1], sep='') # simply the matched SM gene notation
  gene.ranges = t(sapply(strsplit(cf$V4, ' - '),FUN = 'identity'))  
  locs = matrix(match(gene.ranges, anno$ID),nrow = nrow(gene.ranges), ncol=2)
  gene.ranges = cbind(gene.ranges, cf$V1, cf$V5, sprintf('%.1e', as.numeric(cf$V7)))
  gene.ranges = sort.by(gene.ranges, by = locs[,1])
  locs = sort.by(locs, by = locs[,1])
  nc = nrow(locs)
  cluster.ID = cumsum(c(1, locs[2:nc,1]- locs[1:(nc-1),2] - 1 > max.dist.merge))
  # cbind(gene.ranges, cluster.ID)
  
  # get gene ranges and create cluster names
  s = cbind(by(gene.ranges[,1],INDICES = cluster.ID, FUN = function(x){as.character(x)[1]}), 
            by(gene.ranges[,2],INDICES = cluster.ID, FUN = function(x){as.character(x)[length(x)]}),
            paste(c('UU', 'KN')[((regexpr('\\*', by(gene.ranges[,4],INDICES = cluster.ID, FUN = paste,collapse = '_'))>0)|by(gene.ranges[,4]=='0', INDICES=cluster.ID, FUN=any)) + 1],'p',
                  by(gene.ranges[,5],INDICES = cluster.ID, FUN = function(x){sprintf('%.0e',min(as.numeric(as.character(x))))}), '_',
                  by(gene.ranges[,3],INDICES = cluster.ID, FUN = function(x){paste(as.character(x)[unique(c(1,length(x)))],collapse = '_')}),
                  sep=''))
  s = sort.by(s, by = as.numeric(sub(pattern = '.*p([^ _]*)_S.*', replacement = '\\1', s[,3])))
  s = sort.by(s, by = sub(pattern = '(.*)p[^ _]*_S.*', replacement = '\\1', s[,3]), decreasing=T)
  s[,3] = sub('-0', '_', s[,3])
  s[,3] = sub('-', '_', s[,3])
  s[,3] = sub('p0e\\+00', 'p0', s[,3])
  if (is.null(end.to))
    end.to = nrow(s)
  # s2d = cluster.deepAnno(ica.spatial = ica.spatial, gene.ranges = s[30:31,], prot.seq=prot.seq, ipr.anno = ipr.anno, out.file = out.file, extra.genes=extra.genes, append=F)
  s2d = cluster.deepAnno(ica.spatial = ica.spatial, proteinID = proteinID, gff.file = gff.file,  gene.ranges = s, DNA.fasta.file = DNA.fasta.file, prot.fasta.file = prot.fasta.file, iprscan.table.file = iprscan.tab.file, out.file = out.file, extra.genes=extra.genes, 
                         start.from=start.from, end.to=end.to, n.cluster.per.file=n.cluster.per.file, append=F, geneID2cdsID = geneID2cdsID, gene.definition=gene.definition,
                         RORA.iteration = RORA.iteration, species = species, plotLogo=plotLogo, multialn.method = multialn.method, RORA.topOnly=RORA.topOnly)
  invisible(s2d)
}

get.NCBI.blast <- function(query, db = 'nr', no.hit = 100, filter='L', program='blastp'){
  require('annotate')
  baseUrl <- "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi"
  url0 = paste(baseUrl, '?QUERY=',query,'&DATABASE=',db,'&HITLIST_SIZE=',no.hit, '&FILTER=', filter, '&PROGRAM=', program, '&CMD=Put', sep='')
  post <- htmlTreeParse(url0, useInternalNodes = TRUE)
  x <- post[["string(//comment()[contains(., \"QBlastInfoBegin\")])"]]
  rid <- sub(".*RID = ([[:alnum:]]+).*", "\\1", x)
  rtoe <- as.integer(sub(".*RTOE = ([[:digit:]]+).*", "\\1", x)) * 10
  url1 <- sprintf("%s?RID=%s&FORMAT_TYPE=XML&CMD=Get", baseUrl, 
                  rid)
  message("Waiting for NCBI to process the request")
  result <- .tryParseResult(url1)
  results <- tempfile()
  download.file(url1, destfile = results)
  return(results)
}

landmark.dist.by <- function(idx.query, idx.landmarks, groups, ...){
  # Yong Fuga Li, 20150222
  x = as.list(by(1:length(groups), INDICES = groups, 
                 FUN = function(x){d = landmark.dist(idx.query = idx.query[x], idx.landmarks = idx.landmarks[x]);
                 names(d)=intersect(which(idx.query),x);
                 d}, simplify=T))
  out = c();
  for (z in x){
    out = c(out, z)
  }
  out = sort.by(out, as.numeric(names(out)))
  return(out)
}

dist.landmark <- landmark.dist <- landmark.distances <- function(idx.query, idx.landmarks, query = NULL, landmarks = NULL, sequences = NULL, method = 'min'){
  # idx.query: lenght n logical/indicator vector for the query instances on a sequence of length n
  # idx.landmarks: lenght n logical/indicator vector for the landmarks instances on a sequence of length n
  # query: queries - a subset from sequences
  # landmarks: landmarks - a subset from sequences
  # sequences: a character vectors for the elements in a sequences
  # compute the distances between a set of query instances and a set of predefined landmarks on a linear sequence
  # Yong Fuga Li, 20140820, 20141214
  
  # idx.query = mod$metabolism2nd; idx.landmarks = mod$TF
  if (is.logical(idx.query)){
    loc.q = t(which(idx.query))
    colnames(loc.q) = names(idx.query)[loc.q]
  }else{
    loc.q = t(idx.query)
  }
  if (is.logical(idx.landmarks)){
    loc.l = t(t(which(idx.landmarks)))    
    rownames(loc.l) = names(idx.query)[loc.l]
  }else{
    loc.l = t(t(idx.landmarks))        
  }
  if (length(loc.l)==0)
    return(rep(Inf, length(loc.q)))
  d = repmat(loc.q, length(loc.l),1) - repmat(loc.l,1,length(loc.q))
  d.min = colMin(abs(d))
  return(d.min)
}

dist.landmark.which <- landmark.dist.which <- landmark.distances.which <- function(idx.query, idx.landmarks, query = NULL, 
                                                                                   landmarks = NULL, sequences = NULL, 
                                                                                   method = c('both', 'left', 'right'),
                                                                                   include.equal = T){
  # idx.query: lenght n logical/indicator vector for the query instances on a sequence of length n
  # idx.landmarks: lenght n logical/indicator vector for the landmarks instances on a sequence of length n
  # query: queries - a subset from sequences
  # landmarks: landmarks - a subset from sequences
  # sequences: a character vectors for the elements in a sequences
  # compute the distances between a set of query instances and a set of predefined landmarks on a linear sequence
  # Yong Fuga Li, 20141214
  # 20141219: add method
  # idx.query = mod$metabolism2nd; idx.landmarks = mod$TF
  
  method = match.arg(method);
  if (is.logical(idx.query)){
    loc.q = t(which(idx.query))
  }else{
    loc.q = t(idx.query)
  }
  if (is.logical(idx.landmarks)){
    loc.l = t(t(which(idx.landmarks)))
  }else{
    loc.l = t(t(idx.landmarks)) 
  }
  
  if (length(loc.l)==0)
    return(rep(NA, length(loc.q)))
  d = repmat(loc.q, length(loc.l),1) - repmat(loc.l,1,length(loc.q))
  # d.min = colMin(abs(d))
  if (method=='both'){
    idx = max.col(-t(abs(d)))
  }else if (method == 'left'){
    if (include.equal){
      d[d<0] = Inf;
    }else{
      d[d<=0] =Inf;      
    }
    idx = max.col(-t(d))   
  }else if (method =='right'){
    if (include.equal){
      d[d>0] = -Inf;      
    }else{
      d[d>=0] = -Inf;      
    }
    idx = max.col(t(d))    
  }
  if (is.logical(idx.landmarks)){ # change to the original index if input is indicator vector, 20141218
    idx = loc.l[idx]
  }
  return(idx)
}
blast2profile.DD.P <- blast2profile.tblastx <- function(blast.xml = 'AN8428_blastx_nr.xml', query.file = 'AN8428.fasta',  db = 'nr', by=c('query','db', 'align'),
                                                        root = '/Users/yongli/Universe/write/Project_Current/t.NPbioinformatics/Nidulans.SlidingWindow/AutoAnno'){
  ## obtain the consensus sequence based on blastx
  ## Yong Fuga Li
  ## 20140916
  setwd(root)
  
}

Sys.setenv(WISECONFIGDIR='/Users/yongli/Universe/bin/wise2.4.1/wisecfg/')
veriGene <- blast2profile.DP.P <- blast2profile.blastx <- function(blast.xml = 'AN8428_blastx_nr.xml',
                                                                   query.file = 'AN8428.fasta', 
                                                                   db = 'nr', Evalue = 0.1, 
                                                                   by=c('query','hits', 'align'),
                                                                   max.seq = 100,
                                                                   tag = sub('.xml', '', blast.xml),
                                                                   root = '/Users/yongli/Universe/write/Project_Current/t.NPbioinformatics/Nidulans.SlidingWindow/AutoAnno'){
  ## obtain the consensus sequence based on blastx
  ## Yong Fuga Li
  ## 20140916
  require(annotate)
  require(Biostrings)
  require(lattice)
  
  if (is.null(tag)){
    tag = sub('.xml', '', blast.xml)
    
  }
  
  setwd(root)
  querySeq = read.fasta(fasta.file = query.file, type = 'DNA')
  # parse blast results
  bl = blast.xml.parse(blast.xml = blast.xml)
  bl.f = blast.filter(bl, Evalue = Evalue)
  # retrieve hit sequences from GenBank
  gi = sub('gi\\|([^\\|]*)\\|.*', '\\1', bl.f$hit$Hit_def, perl = T)
  n.seq = min(length(gi), max.seq)
  hitSeq = getSEQS(gi[1:n.seq])
  
  # multiple alignment & HMM bulding
  Evalue.string = sub('\\.', '_', paste(Evalue))
  fa.file =  paste(tag, '_hits', n.seq, '_E', Evalue.string, '.fasta', sep='')
  aln.file =  paste(tag, '_hits', n.seq, '_E', Evalue.string, '_aln.fasta', sep='')
  hmm.file =  paste(tag, '_hits', n.seq, '_E', Evalue.string, '.hmm', sep='')
  hmm2.file =  paste(tag, '_hits', n.seq, '_E', Evalue.string, '.hmm2', sep='')
  genewise.gff.file =  paste(tag, '_hits', n.seq, '_E', Evalue.string, '_wise.gff', sep='')
  genewise.fasta.file =  paste(tag, '_hits', n.seq, '_E', Evalue.string, '_wise.fasta', sep='')
  blastx2.asn.file =  paste(tag, '_hits', n.seq, '_E', Evalue.string, '_wise.asn', sep='')
  blastx2.xml.file =  paste(tag, '_hits', n.seq, '_E', Evalue.string, '_wise.xml', sep='')
  genomescan.fasta.file =  paste(tag, '_hits', n.seq, '_E', Evalue.string, '_GenomeScan_Vertibrate.fasta', sep='')
  blastx3.asn.file =  paste(tag, '_hits', n.seq, '_E', Evalue.string, '_GenomeScan_Vertibrate.asn', sep='')
  blastx3.xml.file =  paste(tag, '_hits', n.seq, '_E', Evalue.string, '_GenomeScan_Vertibrate.xml', sep='')
  
  logo.file =  paste(tag, '_hits', n.seq, '_E', Evalue.string, '_aln_Logo.pdf', sep='')
  write.fasta(hitSeq[,'seq'], out.file = fa.file)
  # system(paste('muscle -in', fa.file, '-out', aln.file))
  system(paste('linsi ', fa.file, '>', aln.file))
  system(paste('hmmbuild', hmm.file, aln.file), ignore.stdout=T, ignore.stderr  = T,intern = F)
  system(paste('hmmconvert -2', hmm.file, '>', hmm2.file))
  
  ### genewise
  system(paste('genewise', hmm2.file, query.file, '-hmmer -both -pep -divide \'\' > ', genewise.fasta.file), ignore.stderr  = T,intern = F)
  system(paste('genewise', hmm2.file, query.file, '-hmmer -both -gff -divide \'\' > ', genewise.gff.file), ignore.stderr  = T,intern = F)
  # system('fasttree hits_aln.fasta > hits.tree')
  if (file.info(genewise.fasta.file)$size > 5){
    system(paste('blastx -subject', genewise.fasta.file, '-query', query.file, '-outfmt 11 -out', blastx2.asn.file, '-evalue 1'))
    system(paste('blast_formatter -archive', blastx2.asn.file, '-outfmt 5 -out', blastx2.xml.file))    
  }
  
  # run genome scan
  
  system(paste('blastx -subject', genomescan.fasta.file, '-query', query.file, '-outfmt 11 -out', blastx3.asn.file, '-evalue 1'))
  system(paste('blast_formatter -archive', blastx3.asn.file, '-outfmt 5 -out', blastx3.xml.file))    
  
  # profile model
  aln = readAAMultipleAlignment(aln.file)
  aln.m = maskGaps(aln)
  prof = consensusMatrix(aln)
  prof.m = consensusMatrix(aln.m)
  pdf(logo.file, 10,6)
  print(seqLogo2(prof))
  print(seqLogo2(prof.m))
  dev.off()
  # 
  # consensusString(aln.m)
  # consensusViews(aln.m)
}

blast2profile.DD.D <- blast2profile.blastn <- function(blast.xml = 'AN8428_blastx_nr.xml', query.file = 'AN8428.fasta', db = 'nr', by=c('query','db', 'align')){
  ## obtain the consensus sequence based on blastx
  ## Yong Fuga Li
  ## 20140916
  
}


proMap.hmm <- function(){ # use pHMM to align, score, and generate hints
  # 20141215
  # cluster all hits by genomic locations
  # MAS of hits in each cluster
  # build pHMM models: hmmbuild -O modifiedMSA.txt cS818AN8444.hmm cS818_augoNovoAN8444.t3.aln
  # align query gene models to the pHMM of the same genomic location: hmmsearch --tblout AN8444.hmmtblout.txt cS818AN8444.hmm AN8444.fa > AN8444.hmm.aln
  # score whole protein, and protein part to generate hints
  
}
proMap <- blastp2profile <- blast2profile.PP <- function(blast.asn.file = 'cAfu3g01400_Afu3g01480_swissprot.asn', 
                                                         query.gff.file = 'cAfu3g01400_Afu3g01480.gff', 
                                                         query.faa.file = 'cAfu3g01400_Afu3g01480.fasta',
                                                         DNA.fasta.file = 'cAfu3g01400_Afu3g01480DNA_subseq.fasta',
                                                         geneID2cdsID=function(x){paste(x, '-P', sep='')},
                                                         remove.identical.hits =T,
                                                         tag = sub('.asn', '', blast.asn.file),
                                                         hints.gff = paste(tag, 'proMap.hints.gff'),
                                                         multialn.method = multialn.method,
                                                         plotLogo = T, plotLogo.noInDel = F,
                                                         plot.width = 50, iteration = '',
                                                         db = '/Users/yongli/Universe/data/blastdb/swissprot.fasta',
                                                         by=c('query','db', 'align'),
                                                         aln.maxiters = 4 # previously 2, now 4, 20160818, ref: http://www.drive5.com/muscle/manual/compromise.html
                                                         ){
  ## obtain the consensus sequence based on blastp
  ## output: visualization, protein scoring (conservation score, aln score, mutation score, ins score, del score), consensus protein fasta, intron hints gff
  ## Yong Fuga Li
  ## 20140916, 
  ## 201412-20141213
  ## version 3, 20160618, add start, stop, intergenic region evidences
  require(annotate)
  require(Biostrings)
  require(lattice)
  require(ggplot2);
  no.top.hits2 = 100000L
  require(gridExtra)
  
  if (is.null(tag))
    tag = sub('.asn', '', blast.asn.file) ;
  
  # get CDS and peptides
  qSeq = read.fasta(fasta.files = query.faa.file, type = 'AA')
  qCDSSeq = get.CDS(gene.IDs = rownames(qSeq),
                    gff.file = query.gff.file, DNA.fasta.file = DNA.fasta.file,
                    geneID2cdsID=geneID2cdsID)
  
  for (j in rownames(qCDSSeq)){
    pep = as.character(translate(DNAString(as.character(qCDSSeq[j,1])),if.fuzzy.codon = 'X'));
    if (gsub('\\*', '', pep) != gsub('\\*', '', qSeq[j,'seq'])){
      cat(pep, '\n')
      cat(qSeq[j, 'seq'], '\n')
      warning('CDS translation dose not match protein sequences')
    }
  }
  exon.sizes = sapply(qCDSSeq$exon.sizes, FUN = function(x)as.numeric(unlist(strsplit(as.character(x), ','))))
  names(exon.sizes) = rownames(qCDSSeq)  
  
  #### parse blast results
  system(paste('blast_formatter -archive', blast.asn.file, '-outfmt \'6 qseqid sseqid evalue pident gaps qstart qend\'  -out hits.txt -max_target_seqs ', no.top.hits2)) # 20160505
  # system(paste('blast_formatter -archive', blast.asn.file, '-outfmt \'6 qseqid sseqid evalue\'  -out hits.txt -max_target_seqs ', no.top.hits2))
  hits = read.table('hits.txt', header = F, sep = '\t', as.is=T)
  if (remove.identical.hits){ # remove the blast hits that is identical to the query, 
    hits = hits[!(hits[[4]]==100 & hits[[5]]==0),]    
  }
  # extract all proteins once
  write(as.character(hits[[2]]), 'hitsList.txt');
  system(paste('formatting.pl -idlist hitsList.txt -input ', db, ' -o hits.fasta', sep=''));
  hits = by(hits[[2]], INDICES = hits[[1]], identity)
  nq = nrow(qSeq)
  #   evidence.proteins = data.frame(cSeq=vector('character', length = nq), # protein level evidences
  #                          score = vector('numeric', length = nq), stringsAsFactors=F)
  # rownames(evidence.proteins) = rownames(qSeq);
  # evidence.introns = data.frame(start=c(), end=c(), coverage=c(), score=c()); # intron level evidences
  scores.all = list();
  
  for (qID in rownames(qCDSSeq)){ # extract hits for invididual proteins
    if (qID == 'gene_7282~gene_7283'){
      cat('Here it is')
      print(1)
      cat('Here it is')
    }
    if (!(qID %in% names(hits))){ # proteins without swissprot htis # 20141215
      scores = list(hmm.global.score=NA,
                    hmm.Evalue = NA, 
                    global.score = NA,
                    local.score = NA,
                    match.score = NA,
                    mean.coverage = 0,
                    total.coverage = 0,
                    nSeq = '',
                    cSeq = '',
                    cSeq.long = '',
                    nSeq.naive = '',
                    all = NA, 
                    new.intron=NA, 
                    intron=NA, 
                    new.exon = NA, 
                    nHits = 0,
                    exon=c(), mutation=c())
      scores$CDS = qCDSSeq[qID, ];
      scores.all[[qID]] = scores;
      next
    }
    write(as.character(hits[[qID]]), 'hits1q.txt');
    system(paste('formatting.pl -idlist hits1q.txt -input hits.fasta -o hits1q.fasta', sep=''));
    qSeq1 = qSeq[qID,'seq']; qIDd = paste(qID, as.character(qCDSSeq[qID,'exon.sizes']), sep='_'); names(qSeq1) = qIDd;
    write.fasta(qSeq1, out.file = 'hits1q.fasta', append=T)
    if (!file.exists(paste(tag, iteration, qID, '.aln', sep=''))){
      cat('Aligning ', qID, '\n')
      if (multialn.method=='muscle'){
        system(paste(multialn.method, ' -maxiters ', aln.maxiters, ' -quiet -in hits1q.fasta > ', tag, iteration, qID, '.aln', sep=''))      
      }else if (multialn.method=='mafft'){
        # system(paste('linsi --thread 6 hits1q.fasta > ', tag, iteration, qID, '.aln', sep=''))        
        system(paste('mafft --auto --thread 6 hits1q.fasta > ', tag, iteration, qID, '.aln', sep=''))        
      }else{
        stop('Unsupported alignment method')
      }      
    }else{
      cat('Using existing file ', paste(tag, iteration, qID, '.aln', sep=''),'\n')
    }
    
    # profile model
    aln = readAAMultipleAlignment(paste(tag, iteration, qID, '.aln', sep=''), 'fasta')
    idx = names(aln@unmasked) != qIDd
    aln.m = tryCatch(maskGaps(aln), error = function(e){aln}, finally = NULL) # 20160806
    # aln.m = maskGaps(aln)
    prof = consensusMatrix(as.matrix(aln)[idx,,drop=F])
    prof.m = consensusMatrix(as.matrix(aln.m)[idx,,drop=F])
    
    # exon locations on the aligned sequence
    idx.alned = strsplit(as.character(aln@unmasked[[qIDd]]), '')[[1]] != '-'
    aa.loc.alned = cumsum(idx.alned)
    intron.loc = cumsum(exon.sizes[[qID]])/3; intron.loc=intron.loc[seq2(1,(length(intron.loc)-1),1)]
    intron.loc.alned = approx(aa.loc.alned, 1:length(aa.loc.alned), xout = intron.loc, method = 'linear')$y    
    nrows = ceiling(length(aa.loc.alned)/plot.width);
    
    if (qID == 'g6901.t1'){
      1      
    }
    ## scoring and proposing new sequence
    #     if (qID == 'Afu3g01440'){
    #       1
    #     }
    scores = score.CDS(aln, qIDd = qIDd, intron.locs = intron.loc)
    s = get.hmmer.global.score(aln, qIDd, paste(tag, iteration, sep='')); scores$hmm.global.score = s[2]; scores$hmm.Evalue = s[1];
    scores$CDS = qCDSSeq[qID, ];
    scores.all[[qID]] = scores;
    ### visualization
    pdf.out.file = paste(tag, iteration, qID, 'aln.pdf', sep='_');
    if (plotLogo &  !file.exists(pdf.out.file)){
      pdf(pdf.out.file, 14, 8/4*nrows)
      seqLogo2(prof, intron.locs = intron.loc.alned, qSeq=as.character(aln@unmasked[[qIDd]]), scores = scores$all[,c('IC', 'global.score', 'coverage')], width = plot.width)
      if (plotLogo.noInDel){
        grid.newpage()
        seqLogo2(prof[,idx.alned], intron.locs = intron.loc, qSeq = qSeq1, scores = scores$all[idx.alned,c('IC', 'global.score', 'coverage')], width = plot.width)        
      }
      # print(seqLogo2(prof.m) + geom_vline(xintercept = inron.loc.alned+0.5))
      dev.off()      
    }
    
    # write evidence file
    # write.evidence(scores, )
    # intron coverage and hint file
    #     evidence.introns = c();
    #     # overall protein scores
    #     evidence.proteins[qID, 'global.score'] = sum(scores$global.score.raw)
    #     evidence.proteins[qID, 'local.score'] = sum(scores$loca.scorel)
    #     evidence.proteins[qID, 'match.score'] = sum(scores$match.score)
    # evidence.proteins[qID, 'cSeq'] = cSeq # consensus sequence
    
    ######## write  hints.gff
    # hints.gff
    
  }
  
  return(scores.all)
  # return(list(proteins=evidence.proteins, exons=evidence.exons, introns=evidence.introns))
  # system(paste('hmmbuild', hmm.file, aln.file), ignore.stdout=T, ignore.stderr  = T,intern = F)
  # system(paste('hmmconvert -2', hmm.file, '>', hmm2.file))
  ### blastx on the consensus seq (not including the query)
  
  ### exonerate (instead of genewise) on the consensus seq
  # system(paste('genewise', hmm2.file, query.file, '-hmmer -both -pep -divide \'\' > ', genewise.fasta.file), ignore.stderr  = T,intern = F)
  # system(paste('genewise', hmm2.file, query.file, '-hmmer -both -gff -divide \'\' > ', genewise.gff.file), ignore.stderr  = T,intern = F)
}

get.hmmer.global.score <- function(aln, qIDd, tag){
  # compute gene - hmm alignment global score
  # Yong Fuga Li, 20141219
  hits.file = paste(tag, qIDd, '_hits.faa', sep = '');
  q.file = paste(tag, qIDd, '.faa', sep = '');
  hmm.file = paste(tag, qIDd, '_hits.hmm', sep = '');
  hmmtblout.file = paste(tag, qIDd, '_hmmsearch.tab', sep = '');
  hmmout.file = paste(tag, qIDd, '_hmmsearch.out', sep = '');
  
  idx = qIDd!=rownames(aln);
  aln.mat = as.matrix(aln)
  aln.seq = paste.row(aln.mat[idx, colSums(aln.mat[idx,,drop=F]!='-')!=0,drop=F], collapse=''); # get hits alignment
  qSeq = paste(aln.mat[qIDd, aln.mat[qIDd,]!= '-'], collapse =''); # get query seq
  names(qSeq) = qIDd
  
  write.fasta(aln.seq, hits.file)
  write.fasta(qSeq, q.file)
  # build pHMM models: hmmbuild -O modifiedMSA.txt cS818AN8444.hmm cS818_augoNovoAN8444.t3.aln
  system(paste('hmmbuild ', hmm.file, hits.file))
  # align query gene models to the pHMM of the same genomic location: hmmsearch --tblout AN8444.hmmtblout.txt cS818AN8444.hmm AN8444.fa > AN8444.hmm.aln
  system(paste('hmmsearch --tblout ', hmmtblout.file, hmm.file, q.file, '> ', hmmout.file))
  temp = tryCatch(read.table(hmmtblout.file, header = F), error = function(e){rbind(c(0,0,0,0,Inf,0))}, finally = NULL)
  if (nrow(temp)>1)
    stop('more than one row in hmmsearch out')
  return(c('Evalue'=temp[1,5],'score'=temp[1,6]))
}

score.CDS <- function(aln, qIDd = qIDd, intron.locs = intron.loc, junction.window = 3,
                      query.weighted.consensus = F, # used query weighted consensus sequence when proposing new exons
                      new.intron.min.size=5, 
                      new.intron.min.junction.coverage=0.3,
                      new.intron.min.junction.score=0,
                      new.exon.min.size = 5,
                      new.exon.old.intron.max.junction.coverage=0.8,
                      new.exon.old.intron.max.junction.score=1,
                      new.exon.old.intron.max.distance = 50,
                      new.exon.min.coverage.stringent = 0.8,
                      new.exon.old.intron.overhaning.exon.maxscore =0,
                      new.exon.min.coverage=0.3){
  # Yong Fuga Li
  # 20141212-14
  # scoring a predicted protein against a set of known proteins (conservation score, aln score, mutation score, ins score, del score)
  # names(aln@unmasked)
  ## version 3, 20160618, add start, stop, intergenic region evidences
  
  idx.alned = strsplit(as.character(aln@unmasked[[qIDd]]), '')[[1]] != '-' # query aligned locations
  aa.loc.alned = cumsum(idx.alned)
  intron.locs = intron.locs #, max(aa.loc.alned)
  intron.loc.alned = approx(aa.loc.alned, 1:length(aa.loc.alned), xout = intron.locs, method = 'linear')$y    
  
  idx = qIDd!=rownames(aln)
  prof = consensusMatrix(as.matrix(aln)[idx,,drop=F])
  if (!('-' %in% rownames(prof))){# fix a bug for alignment without any gaps, 20141223
    prof=rbind('-'=0, prof);
  }
  scores = data.frame(IC = pwm2ic(prof,pseudocount = 1)$IC,
                      coverage = colSums(prof[!(rownames(prof) %in% c('-', '#')),])) # scores along alignment
  scores.intron = data.frame() # scores for the introns
  scores.new.intron = data.frame() # scores for proposed new introns
  scores.exon = data.frame() # scores for the exons
  scores.start = data.frame() # scores for the start codon
  scores.stop = data.frame() # scores for the stop codon
  scores.irpart = data.frame() # scores for the intergenic region part
  
  # calculate alignment score between the query sequences and the rest of tha alignment
  aln = as.matrix(aln)
  aln.local = aln; # local alignment
  for (j in 1:nrow(aln.local)){
    tt = cumsum(aln.local[j,] != '-')
    aln.local[j,tt == 0 | tt == max(tt)] = '*'
  }
  if (exists('BLOSUM62', mode = "matrix"))
    remove('BLOSUM62');
  data(BLOSUM62)
  BLOSUM62 = cbind(rbind(BLOSUM62, '-' = -4), '-' = -4); BLOSUM62['-', '-'] = 0
  BLOSUM62['*',] = 0; BLOSUM62[,'*'] = 0; # N or C term indels
  BLOSUM62.noIndel = BLOSUM62;
  BLOSUM62.noIndel['-',] = 0; BLOSUM62.noIndel[,'-'] = 0; # N or C term indels
  score.mat = mat2xyz(BLOSUM62, sym=F); aapairs = paste(score.mat[,1], score.mat[,2], sep='')
  score.mat = score.mat[,3]; names(score.mat) = aapairs # scoring matrix to vector
  score.mat.noIndel = mat2xyz(BLOSUM62.noIndel, sym = F); aapairs = paste(score.mat.noIndel[,1], score.mat.noIndel[,2], sep='')
  score.mat.noIndel = score.mat.noIndel[,3]; names(score.mat.noIndel) = aapairs # scoring matrix to vector
  scores$global.score = 0; scores$local.score = 0; scores$match.score = 0
  for (i in which(idx)){
    scores$global.score = scores$global.score + score.mat[paste(aln[i,], aln[qIDd,], sep='')]
    scores$local.score = scores$local.score + score.mat[paste(aln.local[i,], aln.local[qIDd,], sep='')]
    scores$match.score = scores$match.score + score.mat.noIndel[paste(aln[i,], aln[qIDd,], sep='')]
  }
  
  nHits = nrow(aln)-1
  scores$global.score.raw = scores$global.score
  scores$global.score = (scores$global.score/(nrow(aln)-1))/max(abs(score.mat))
  scores$IC = scores$IC/max(scores$IC)
  scores$coverage = (scores$coverage)/nHits
  
  scores.intron = data.frame(locs = intron.locs, locs.alned = intron.loc.alned,
                             # coverage = windowMeans(scores$coverage[idx.alned], locs=intron.locs, window.size=2), # alignment evidences at the splicing sites
                             coverage = site.coverage(aln, qIDd, intron.loc.alned, p.indel = 0.5, normalize = T),
                             # match.score = windowMeans(scores$match.score[idx.alned], locs=intron.locs, window.size=2)) # alignment evidences at the splicing sites
                             match.score = windowMeans(scores$match.score, locs=intron.loc.alned, window.size=2)) # 20141218 alignment evidences at the splicing sites 
  
  ##### wrong exon, new intron: old exon coverage = 0, new intron coverage > 30%, match.score > 0
  coverage = round(scores$coverage[idx.alned]);
  scores.new.intron = data.frame(start = which(coverage == 0 & diff(c(Inf, coverage) !=0)), 
                                 end = which(coverage == 0 & diff(c(coverage, Inf) !=0)));
  idx.exon.new = coverage!=0
  intron.locs.new = cumsum(idx.exon.new)[scores.new.intron$start]
  scores.new.intron = cbind(scores.new.intron, 
                            coverage = windowMeans((scores$coverage[idx.alned])[idx.exon.new], locs=intron.locs.new, window.size=2), # alignment evidences at the splicing sites
                            match.score = windowMeans((scores$match.score[idx.alned])[idx.exon.new], locs=intron.locs.new, window.size=2)) # alignment evidences at the splicing sites
  scores.new.intron = scores.new.intron[scores.new.intron$end-scores.new.intron$start+1 >= new.intron.min.size &
                                          scores.new.intron$coverage >= new.intron.min.junction.coverage  &
                                          scores.new.intron$match.score >= new.intron.min.junction.score,]
  # scores$splice.after = # is a splicing site between this aa and the 3' aa
  # scores$mutation = # is the aa location likely containing a mutation
  
  ####### wrong intron, new exon: gap in query sequence, new intron coverage < 30%, match.score < 0
  cSeq.with.Ins = rownames(prof)[max.col(t(prof), ties.method = 'first')]
  prof['-',prof['-',]!=nHits] = 0; # consensus target sequence    
  cSeq = rownames(prof)[max.col(t(prof), ties.method = 'first')]
  nSeq.long <- nSeq <- aln[qIDd,]
  is.potential.new.exon = scores$coverage > new.exon.min.coverage  & nSeq == '-';
  scores.new.exon = data.frame(start = which(is.potential.new.exon & diff(c(Inf, is.potential.new.exon*1)) !=0), 
                               end = which(is.potential.new.exon & diff(c(is.potential.new.exon*1, Inf)) !=0));
  scores.new.exon$dist.to.old.junctions = landmark.dist((scores.new.exon$start+scores.new.exon$end)/2, intron.loc.alned, ncol)
  scores.new.exon$mean.coverage = sapply(seq2(1,nrow(scores.new.exon),1), function(x){mean(scores$coverage[scores.new.exon$start[x]:scores.new.exon$end[x]])})
  # intron.locs.new = cumsum(!is.potential.new.exon)[scores.new.exon$start]
  idx.matched.introns = landmark.dist.which((scores.new.exon$start+scores.new.exon$end)/2, intron.loc.alned) # index of the nearest introns
  
  if (length(intron.loc.alned)){
    scores.new.exon = cbind(scores.new.exon, 
                            nearest.intron.locs = scores.intron[idx.matched.introns,'locs'],
                            nearest.intron.junction.coverage = scores.intron[idx.matched.introns,'coverage'],
                            nearest.intron.junction.score = scores.intron[idx.matched.introns,'match.score'])    
  }else if (nrow(scores.new.exon)){ # intron less cases
    scores.new.exon = cbind(scores.new.exon, 
                            nearest.intron.locs = Inf,
                            nearest.intron.junction.coverage = 1,
                            nearest.intron.junction.score = Inf)        
  }else{
    scores.new.exon = cbind(scores.new.exon, 
                            nearest.intron.locs = vector('numeric',0),
                            nearest.intron.junction.coverage = vector('numeric',0),
                            nearest.intron.junction.score = vector('numeric',0))       
  }
  scores.new.exon = scores.new.exon[((scores.new.exon$nearest.intron.junction.coverage <= new.exon.old.intron.max.junction.coverage  &
                                        scores.new.exon$nearest.intron.junction.score <= new.exon.old.intron.max.junction.score) |
                                       scores.new.exon$mean.coverage >= new.exon.min.coverage.stringent  &
                                       nHits >=3) & scores.new.exon$dist.to.old.junctions <= new.exon.old.intron.max.distance &
                                      scores.new.exon$end-scores.new.exon$start + 1 >=new.exon.min.size,]
  
  ###### propose new protein sequence
  idx.new.intron = unlist(sapply(seq2(1, nrow(scores.new.intron),1), function(x){scores.new.intron$start[x]:scores.new.intron$end[x]}))
  idx.new.exon = unlist(sapply(seq2(1, nrow(scores.new.exon),1), function(x){scores.new.exon$start[x]:scores.new.exon$end[x]}))
  #modify the new sequence based on the new introns and new exons  
  cat(nrow(scores.new.intron), ' new introns for', qIDd, '\n')
  cat(nrow(scores.new.exon), ' new exons for', qIDd, '\n')
  nSeq[which(idx.alned)[idx.new.intron]] = '-'
  nSeq[idx.new.exon] = cSeq[idx.new.exon]
  rbind(cSeq, nSeq)
  
  # scores$seq = aln[qIDd,]
  to.change = nSeq.long == '-' & cSeq.with.Ins !='-'
  nSeq.long[to.change] = cSeq.with.Ins[to.change]
  # scores.intron$coverage = scores.intron$coverage*nHits
  return(list(global.score = sum(scores$global.score.raw),
              local.score = sum(scores$local.score),
              match.score = sum(scores$match.score),
              total.coverage = sum(scores$coverage[aln[qIDd,]!='-']),
              mean.coverage = mean(scores$coverage[aln[qIDd,]!='-']),
              nSeq = paste(nSeq[nSeq!='-'], collapse =''),
              cSeq = paste(cSeq.with.Ins[cSeq.with.Ins!='-'], collapse=''),
              cSeq.long = paste(cSeq[cSeq!='-'], collapse =''), nSeq.naive = paste(nSeq.long[nSeq.long!='-'], collapse =''),
              all = scores, 
              new.intron=scores.new.intron, 
              intron=scores.intron, 
              new.exon = scores.new.exon, 
              nHits = nHits,
              exon=c(), mutation=c()))
  # qSeq = aln@unmasked[[qIDd]]
}

site.coverage <- function(aln, qIDd, intron.loc.alned, p.indel = 0.5, normalize=T){
  # 20141219, Yong Fuga Li
  if (is.null(intron.loc.alned) || !length(intron.loc.alned))
    return(c())
  delta = 1E-10
  idx = qIDd!=rownames(aln)
  qSeq = (aln[qIDd,] != '-'); r = range(which(qSeq)); # 20141231, padding 2 on both sides to indicate the end of sequences
  qSeq = c(2, (aln[qIDd,] != '-') * 1,2); qSeq[c(r[1]-1, r[2]+1)] = 2;  # 20141231, padding 2 on both sides to indicate the end of sequences
  sAln = cbind(0, (aln[idx,,drop=F] != '-')*1, 0) # 20141231, padding 1 on both sides to indicate the end of sequences
  inDels =  t(qSeq != t(sAln))*1;
  #DelDels = t((1-qSeq) * t(1-sAln))
  matches = t(qSeq * t(sAln))
  coverage = 0
  for (i in 1:nrow(sAln)){
    # x.R = landmark.dist.which(intron.loc.alned, inDels[i,]>0, method = 'right', include.equal = F)-1 
    x.R = landmark.dist.which(intron.loc.alned+1, inDels[i,]>0, method = 'right', include.equal = F)-1 # 20141231, adjucting for the 1 padded to the ends - padding 1 on both sides to indicate the end of sequences
    # d.R = x.R -intron.loc.alned; 
    # d.R[d.R<0] = 0
    # d.R = d.R - windowSums(DelDels[i,], locs.s = intron.loc.alned, locs.e = x.R+1)
    # d.R = windowSums(matches[i,], locs.s = intron.loc.alned, locs.e = x.R+1)
    # x.L = landmark.dist.which(intron.loc.alned, inDels[i,2:ncol(inDels)]>0, method = 'left', include.equal = F)+1
    d.R = windowSums(matches[i,], locs.s = intron.loc.alned + 1, locs.e = x.R+1) # 20141231
    x.L = landmark.dist.which(intron.loc.alned+1, inDels[i,]>0, method = 'left', include.equal = F)+1 # 20141231, 
    # d.L = intron.loc.alned - x.L;
    # d.L[d.L<0] = 0
    # d.L = d.L - windowSums(DelDels[i,], locs.s = x.L-1, locs.e = intron.loc.alned)
    d.L = windowSums(matches[i,], locs.s = x.L-1, locs.e = intron.loc.alned+1)
    # rbind(intron.loc.alned, landmark.dist.which(intron.loc.alned+delta, inDels[i,]>0, method = 'left'))
    coverage = coverage + (1-(1-p.indel)^d.R)*(1-(1-p.indel)^d.L)
  }
  if (normalize)
    coverage = coverage/(nrow(aln)-1)
  return(coverage)
}

write.proMap <- function(pMap, score.file = paste('pMap',tag, '.xls', sep=''), nSeq.file = paste('pMap_nSeq',tag, '.faa', sep=''), 
                         nSeq.naive.file = paste('pMap_nSeqNaive',tag, '.faa', sep=''), 
                         cSeq.long.file = paste('pMap_cSeqLong',tag, '.faa', sep=''), tag='', iteration = 0, append = T){
  # write proMap output to files
  # Yong Fuga Li, 20141214
  
  ## protein scores
  scores = c()
  seqs = c()
  for (p in names(pMap)){
    scores = rbind(scores, cbind(pMap[[p]]$CDS, pHMM.score = pMap[[p]]$hmm.global.score, pHMM.Evalue = pMap[[p]]$hmm.Evalue,  
                                 match=pMap[[p]]$match.score, total.coverage=pMap[[p]]$total.coverage, mean.coverage=pMap[[p]]$mean.coverage, 
                                 local=pMap[[p]]$local.score, global=pMap[[p]]$global.score, iteration=iteration))
    seqs = rbind(seqs, c(nSeq = pMap[[p]]$nSeq, cSeq.long=pMap[[p]]$cSeq.long, nSeq.naive = pMap[[p]]$nSeq.naive))
  }
  scores$prot_len = nchar(as.character(scores$seq))/3
  rownames(scores) <- rownames(seqs) <- names(pMap)
  rownames(scores)[scores[,'iteration']!=''] <- paste(names(pMap), scores[,'iteration'], sep='.')[scores[,'iteration']!='']; 
  if (append){
    write.table(scores, append = append, file = score.file, quote = F, sep = '\t', row.names = T, col.names = F)
  }else{
    write.table(scores, append = append, file = score.file, quote = F, sep = '\t', row.names = T, col.names = NA)
  }
  write.fasta(seqs[nchar(seqs[,'nSeq'])>0,'nSeq'], append = F, out.file = nSeq.file)
  write.fasta(seqs[nchar(seqs[,'nSeq.naive'])>0,'nSeq.naive'], append = F, out.file = nSeq.naive.file)
  write.fasta(seqs[nchar(seqs[,'cSeq.long'])>0,'cSeq.long'], append = F, out.file = cSeq.long.file)
  ## proposed protein sequences to fasta format  
  return(score.file)
}

blast2profile.PD.P <- blast2profile.tblastn <- function(blast.xml = 'AN8428_blastx_nr.xml', query.file = 'AN8428.fasta', db = 'nr', by=c('query','db', 'align')){
  ## obtain the consensus sequence based on blastx
  ## Yong Fuga Li
  ## 20140916
  
}

cluster.success.rate <- function(n = 6, alpha.g=0.78, # gene level success rate
                                 par = list(beta.neg = 0.7, # false genes called error
                                            beta.pos=0.1, # true genes called error
                                            gamma.neg=0.25, # make correction among called error from false genes
                                            gamma.pos=0.25, # make correction among called error from true genes
                                            delta=0.8)){ # success rate of correction (only for called errors from false genes)
  # 20140917, success rate in the rescue or abandon approach
  require(gridExtra)
  alpha.c = alpha.g^n # cluster level success rate
  
  abandon.g = (1-alpha.g) * par$beta.neg * (1-par$gamma.neg) + alpha.g * par$beta.pos * (1-par$gamma.pos)
  success.g = alpha.g * (1-par$beta.pos) + (1-alpha.g) * par$beta.neg * par$gamma.neg * par$delta
  fail.g = (1-alpha.g) * (1 - par$beta.neg) + (1-alpha.g) * par$beta.neg * par$gamma.neg * (1-par$delta) +
    alpha.g * par$beta.pos * par$gamma.pos
  p = c(abandon = abandon.g, success=success.g, fail=fail.g)
  abandon.c = 1-(1-abandon.g)^n # % of clusters abandoned
  alpha.g.wRorA = p[2]/(p[2]+p[3])
  alpha.c.wRorA = alpha.g.wRorA ^ n
  
  out = data.frame(success.rate.regular = alpha.c, success.rate.RORA = alpha.c.wRorA, success.rate.Abandon = alpha.c.wRorA, 
                   percent.abandon.RORA = abandon.c, percent.abandon.Abandon = abandon.c, 
                   success.rate.gene = alpha.g, success.rate.gene.RORA = alpha.g.wRorA, success.rate.gene.Abandon = alpha.g.wRorA,
                   row.names = n)
  
  ### only abandon, without rescues
  par$gamma.neg <- par$gamma.pos <- 0
  
  abandon.g = (1-alpha.g) * par$beta.neg * (1-par$gamma.neg) + alpha.g * par$beta.pos * (1-par$gamma.neg)
  success.g = alpha.g * (1-par$beta.pos) + (1-alpha.g) * par$beta.neg * par$gamma.neg * par$delta
  fail.g = (1-alpha.g) * (1 - par$beta.neg) + (1-alpha.g) * par$beta.neg * par$gamma.neg * (1-par$delta) +
    alpha.g * par$beta.pos * par$gamma.pos
  p = c(abandon = abandon.g, success=success.g, fail=fail.g)
  abandon.c = 1-(1-abandon.g)^n # % of clusters abandoned
  alpha.g.wRorA = p[2]/(p[2]+p[3])
  alpha.c.wRorA = alpha.g.wRorA ^ n
  out$success.rate.Abandon = alpha.c.wRorA; 
  out$percent.abandon.Abandon = abandon.c;  
  out$success.rate.gene.Abandon = alpha.g.wRorA;  
  
  out1 = out[1:3]; names(out1) = c('regular', 'RORA', 'RORA\nno rescue');
  xyz = mat2xyz(as.matrix(out1), sym=F)
  q1 = theme.ggplot(ggplot(xyz, aes(x=y, y = z*100, fill = x)) + geom_bar(stat = 'identity', position = 'dodge') + xlab('') + ylab('Cluster success%') + 
                      labs(fill='Cluster Size')) + theme(axis.text.x=element_text(angle=0))
  
  out2 = out[4:5]; names(out2) = c('RORA', 'RORA\nno rescue');
  xyz = mat2xyz(as.matrix(out2), sym=F)
  q2 = theme.ggplot(ggplot(xyz, aes(x=y, y = z*100, fill = x)) + geom_bar(stat = 'identity', position = 'dodge') + xlab('') + ylab('Cluster Abandoned%') + 
                      labs(fill='Cluster Size'), legend.position = 'none') + theme(axis.text.x=element_text(angle=0))
  
  grid.arrange(q1, q2, ncol=2,  widths = c(3,2))
  return(out)
  
}


xlsx.color.NPGC <- color.NPGC.xlsx <- function(xlsx.file = 'nidulans.deepAnno.all.xlsx', out.file=paste('colored_', xlsx.file, sep='')){
  # Yong Fuga Li, 20141004
  write('1. Typically there are 5 extra genes included on each side of a cluster. Below explains the coloring scheme used in the deepAnnotation tables.\n', file = 'readme_deepAnno.txt', append = F)  
  write('2. Black box - cluster boundary: expression > 9 or CS < 0.5\n', file = 'readme_deepAnno.txt', append = T)
  xlsx.color(xlsx.file = xlsx.file, FUN.select = FUN.select.boundary, border.color = 'black', out.file = out.file, na.strings='|')  
  write('3. Green fill - promising expression feature (CS column): expression clustering coefficient > 3\n', file = 'readme_deepAnno.txt', append = T)
  xlsx.color(xlsx.file = out.file, FUN.select = FUN.select.promising, fill.color = 'green', out.file = out.file, na.strings='|')
  write('4. Green box - promising: oxidoreductase (oxidoreductase|P450|oxidase|dehydrogenase|oxygenase|reductase), \n', file = 'readme_deepAnno.txt', append = T)
  xlsx.color(xlsx.file = out.file, FUN.select = FUN.select.oxidoreductase, border.color = 'green', out.file = out.file, na.strings='|')
  write('5. Green fill - promising protein function (domains/annotation/top.5.hits columns): (anabolism) transferase, synthase, synthetase, ligase, \n', file = 'readme_deepAnno.txt', append = T)
  xlsx.color(xlsx.file = out.file, FUN.select = FUN.select.catabolism, fill.color = 'green', out.file = out.file, na.strings='|')  
  write('6. Blue fill - special: llm, laeA, molybdenum containing, \n', file = 'readme_deepAnno.txt', append = T)
  xlsx.color(xlsx.file = out.file, FUN.select = FUN.select.special, fill.color = 'light blue', out.file = out.file, na.strings='|')
  write('7. Purple box - interesting: or length > 800 aa, homology with swissprot proteins high (>75%) or low (<25%), polyketide or alkaloid or terpenoid or terpene or nonribosomal peptide mentioned in domain annotations, swissprot hits, or existing genome annotations, \n', file = 'readme_deepAnno.txt', append = T)
  xlsx.color(xlsx.file = out.file, FUN.select = FUN.select.interesting, border.color = 'purple', out.file = out.file, na.strings='|')
  write('8. Purple fill - potential known clusters: polyketide or alkaloid or terpenoid or terpene or nonribosomal peptide mentioned in domain annotations, swissprot hits, or existing genome annotations, \n', file = 'readme_deepAnno.txt', append = T)
  xlsx.color(xlsx.file = out.file, FUN.select = FUN.select.maybeKU, fill.color = 'purple', out.file = out.file, na.strings='|')  
  write('9. Red box - warning (possible gene structure error): average intron size > 100, or average exon size < 100, \n', file = 'readme_deepAnno.txt', append = T)
  xlsx.color(xlsx.file = out.file, FUN.select = FUN.select.warning, border.color = 'red', out.file = out.file, na.strings='|')  
  write('10. Red fill - boring: human annotated/known SM cluster genes in current genome annotation, \n', file = 'readme_deepAnno.txt', append = T)
  xlsx.color(xlsx.file = out.file, FUN.select = FUN.select.boring, fill.color = 'red', out.file = out.file, na.strings='|')
  # write.xlsx(append = T, file = out.file, 'Green - promising: expression clustering coefficient > 3\nGreen box - promising: P450\nBlue - special: llm, laeA, molybdenum containing\nPurple box - interesting: or length > 1000 aa, homology with swissprot proteins high (>50%) or low (<25%), polyketide or alkaloid or terpenoid or terpene or nonribosomal peptide mentioned in domain annotations, swissprot hits, or existing AspGD annotations\nRed - boring: annotated SM genes in current AspGD annotation')
}


enzyme.clustering.simple <- function(gff.file = gff.file, iprscan.tab.file = iprscan.tab.file, gene.definition = c('gene', 'transcript', 'mRNA'),
                                     enzyme.definition = 'P450', max.window = 6, min.enzymes = 2){
  # find minimum windows of size max.window or smaller that contains the largest number of enzymes
  # both k and n can be ranges
  # Yong FUga Li, 20141007
  
  require(gplots)
  ## read gff
  gff.format = sub('^.*\\.([^\\.]*$)', '\\1', gff.file)
  # anno = read.gff3(gff.file, format=gff.format)
  anno = import.gff(gff.file) # 20160502
  idx.gene = (anno$type==gene.definition)
  anno = anno[idx.gene, ]
  anno = sort.intervals(anno)
  n = length(anno)
  # anno$ID = sub('transcript:', '', anno$ID)
  
  
  ipr.anno = iprscan.flat(iprscan.tab.file)
  
  ipr.anno = mat.fill.row(t(t(ipr.anno)), row.names = anno$ID, default = '')[,1]
  is.enzyme = regexpr(pattern=enzyme.definition, text = as.character(as.vector(anno$Note)), perl=T, ignore.case=T) > 0 | 
    regexpr(pattern=enzyme.definition, text = as.character(as.vector(ipr.anno.1)), perl=T, ignore.case=T) > 0
  window.size = 2:3
  counts = sapply(window.size, FUN = function(x){return(runsum.2(is.enzyme+0, k = x, addzeros = T))})
  rownames(counts) = anno$ID; colnames(counts) = window.size
  max(counts)
  sum(counts>=2)
  tt = unique(which(counts>=2, arr.ind = T)[,1])
  cbind(counts[sort(tt),], sort(tt))
  
  is.cluster = (counts == k.enzyme[1])
  for (i in seq2(2, length(k.enzyme), by = 1))
    is.cluster = is.cluster | (counts == k.enzyme[i])
  
  is.cluster = lapply(apply(k.enzyme, function(x){return(counts == x)}))
  
}

NPGC.query <- NPGC.scan <- function(gff.file=NULL, iprscan.tab.file = NULL, query=list(func = list('P450', 'O-methyltransferase'),
                                                                                       freq = list(2:10, 1:10)), 
                                    window.size=15, out.file = 'Tvirens.xls', window.extend=window.size,
                                    gene.definition = c('gene', 'transcript', 'mRNA'), proteinID = 'ID', 
                                    max.dist.merge = window.size){
  # find a window of size 15 or less that meet the gene function query criteria
  # YF Li
  # 20141028, 20141111
  require('xlsx')
  gff.format = sub('^.*\\.([^\\.]*$)', '\\1', gff.file)
  # anno = read.gff3(gff.file, format=gff.format)
  anno = import.gff(gff.file) # 20160502
  idx.gene = (anno$type==gene.definition)
  colnames(anno@elementMetadata) = toupper(colnames(anno@elementMetadata))
  m = match(anno$PARENT[idx.gene], anno$ID)
  anno$DESCRIPTION[idx.gene] = anno$DESCRIPTION[m]   # transfer annotation from parents to childs
  anno$PARENT = as.character(anno$PARENT); anno$PARENT[is.na(anno$PARENT)] = anno$ID[is.na(anno$PARENT)]  
  anno = anno[idx.gene, ]
  anno = sort.intervals(anno)
  n = length(anno)
  # anno$ID = sub('transcript:', '', anno$ID)
  
  ipr.anno = iprscan.flat(iprscan.tab.file, na.strings = c('-', 'NA', 'NULL'))
  ipr.anno = mat.fill.row(t(t(ipr.anno)), row.names = anno@elementMetadata[,toupper(proteinID)], default = '')[,1]
  names(ipr.anno) = anno$ID
  to.keep = ones(n)
  enzyme.count.all = c()
  is.enzyme.all = c()
  if (!is.null(anno$NOTE)){
    desc.fname = 'NOTE'
  }else if (!is.null(anno$DESCRIPTION)){
    desc.fname = 'DESCRIPTION'  
  }else{
    warning('No description or Note field for the annotation of genes')      
    desc.fname = 'NOTE'
    anno$NOTE = ''
  }
  
  for (i in 1:length(query$func)){
    enzyme.definition = query$func[[i]]
    enzyme.count = query$freq[[i]]
    
    is.enzyme = (regexpr(pattern=enzyme.definition, text = as.character(as.vector(anno@elementMetadata[[desc.fname]])), perl=T, ignore.case=T) > 0);
    is.enzyme[is.na(is.enzyme)] = 0; 
    is.enzyme = is.enzyme | (regexpr(pattern=enzyme.definition, text = as.character(as.vector(ipr.anno)), perl=T, ignore.case=T) > 0)          
    names(is.enzyme) = anno$ID
    
    counts = runsum.by(is.enzyme+0, k = window.size, by = as.character(anno@seqnames))
    to.keep = to.keep & (counts %in% enzyme.count)
    enzyme.count.all = cbind(enzyme.count.all, counts)
    is.enzyme.all = cbind(is.enzyme.all, is.enzyme)
    cat(enzyme.definition, sum(is.enzyme), '\n')
  }
  cat('# of total gene', length(anno))
  colnames(enzyme.count.all) <- colnames(is.enzyme.all) <- query$func
  # sum(to.keep)
  #core.regions = intersect(extend.index(which(to.keep), window.size), which(rowSums(is.enzyme.all)>0))
  core.regions = extend.index(which(to.keep), window.size, sides='down', do.unique = F);
  m = match(core.regions, which(rowSums(is.enzyme.all)>0)) # 20141111
  core.regions =  core.regions[!is.na(m)]
  
  # merge clusters
  gene.ranges = unique(cbind(by(core.regions, names(core.regions), FUN = min), by(core.regions, names(core.regions), FUN = max)), MARGIN = 1, drop=F)
  
  gene.ranges = sort.by(gene.ranges, by = gene.ranges[,1])
  nc = nrow(gene.ranges)
  if (is.null(nc) | !nc){
    cat('\nNumber of clusters: 0')    
    return(NULL)  
  }
  cluster.ID = cumsum(c(1, gene.ranges[seq2(2,nc,1),1]- gene.ranges[seq2(1,nc-1,1),2] - 1 > max.dist.merge))
  gene.ranges = data.frame(from=sapply(by(gene.ranges[,1],INDICES = cluster.ID, FUN = function(x){x[1]}), FUN = 'identity'),
                           to = sapply(by(gene.ranges[,2],INDICES = cluster.ID, FUN = function(x){x[length(x)]}), FUN = 'identity'),
                           name= sapply(by(rownames(gene.ranges),INDICES = cluster.ID, FUN = function(x){paste(c(x[1], x[length(x)]), collapse = '_')}), FUN = 'identity'))
  geneID2clusterID = lapply(1:nrow(gene.ranges), function(x){cbind(as.character(anno$ID)[gene.ranges[x,1]:gene.ranges[x,2]], rep(as.character(gene.ranges[x,3]), gene.ranges[x,2]-gene.ranges[x,1]+1))})
  if (is.list(geneID2clusterID)){
    geneID2clusterID = do.call(rbind, geneID2clusterID)
  }
  gene.ranges = cbind(from=as.character(anno$ID)[gene.ranges[,1]], to=as.character(anno$ID)[gene.ranges[,2]], name=as.character(gene.ranges[,3]))
  cat('\nNumber of clusters: ', nrow(gene.ranges))
  
  #### output
  to.keep.extend = extend.index(core.regions, window.extend, sides='both', do.unique=T)
  to.keep.extend = to.keep.extend[to.keep.extend<=length(anno) & to.keep.extend>=1]
  anno$PARENT[1] ==c()
  is.enzyme.all[] = c('', 'Yes')[is.enzyme.all+1]
  out = cbind(chr = as.character(anno@seqnames)[], gene=anno$ID, 'protein ID' = anno@elementMetadata[,toupper(proteinID)], Existing.Anno = anno@elementMetadata[,toupper(desc.fname)],
              is.enzyme.all, domains = ipr.anno)[to.keep.extend,]
  
  rownames(geneID2clusterID) = geneID2clusterID[,1];
  out = cbind(out, clusterID = mat.fill.row(geneID2clusterID, rownames(out), '')[,2])
  
  write.xlsx(out, out.file)
  return(gene.ranges)
  
}

NPGC.mutliscan <- function(gff.files=NULL, iprscan.tab.files = NULL, prot.fasta.files = NULL,
                           query=list(func = list('P450', 'O-methyltransferase'),
                                      freq = list(2:10, 1:10)), 
                           window.size=15, out.file = 'Tvirens.xls', window.extend=window.size,
                           gene.definition = c('gene', 'transcript', 'mRNA'), proteinID = 'ID', 
                           max.dist.merge = window.size){
  # deepAnno.clusters
  # find a window of size 15 or less that meet the gene function query criteria
  # YF Li
  # Step 1: NPGC.query for multiple genomes using the same query
  # Step 2: Blast search of the identified query against the other genome
  # Step 3: BBH ortholog assignemnet 
  # Step 4: Output all clusters with the ortholog information
  require('xlsx')
  gff.format = sub('^.*\\.([^\\.]*$)', '\\1', gff.file)
  # anno = read.gff3(gff.file, format=gff.format)
  anno = import.gff(gff.file) # 20160502
  idx.gene = (anno$type==gene.definition)
  colnames(anno@elementMetadata) = toupper(colnames(anno@elementMetadata))
  m = match(anno$PARENT[idx.gene], anno$ID)
  anno$DESCRIPTION[idx.gene] = anno$DESCRIPTION[m]   # transfer annotation from parents to childs
  anno$PARENT = as.character(anno$PARENT); anno$PARENT[is.na(anno$PARENT)] = anno$ID[is.na(anno$PARENT)]  
  anno = anno[idx.gene, ]
  anno = sort.intervals(anno)
  n = length(anno)
  #anno$ID = sub('transcript:', '', anno$ID)
  
  ipr.anno = iprscan.flat(iprscan.tab.file, na.strings = c('-', 'NA', 'NULL'))
  ipr.anno = mat.fill.row(t(t(ipr.anno)), row.names = anno@elementMetadata[,toupper(proteinID)], default = '')[,1]
  names(ipr.anno) = anno$ID
  to.keep = ones(n)
  enzyme.count.all = c()
  is.enzyme.all = c()
  if (!is.null(anno$NOTE)){
    desc.fname = 'NOTE'
  }else if (!is.null(anno$DESCRIPTION)){
    desc.fname = 'DESCRIPTION'  
  }else{
    warning('No description or Note field for the annotation of genes')      
    desc.fname = 'NOTE'
    anno$NOTE = ''
  }
  
  for (i in 1:length(query$func)){
    enzyme.definition = query$func[[i]]
    enzyme.count = query$freq[[i]]
    
    is.enzyme = (regexpr(pattern=enzyme.definition, text = as.character(as.vector(anno@elementMetadata[[desc.fname]])), perl=T, ignore.case=T) > 0);
    is.enzyme[is.na(is.enzyme)] = 0; 
    is.enzyme = is.enzyme | (regexpr(pattern=enzyme.definition, text = as.character(as.vector(ipr.anno)), perl=T, ignore.case=T) > 0)          
    names(is.enzyme) = anno$ID
    
    counts = runsum.by(is.enzyme+0, k = window.size, by = as.character(anno@seqnames))
    to.keep = to.keep & (counts %in% enzyme.count)
    enzyme.count.all = cbind(enzyme.count.all, counts)
    is.enzyme.all = cbind(is.enzyme.all, is.enzyme)
    cat(enzyme.definition, sum(is.enzyme), '\n')
  }
  cat('# of total gene', length(anno))
  colnames(enzyme.count.all) <- colnames(is.enzyme.all) <- query$func
  # sum(to.keep)
  #core.regions = intersect(extend.index(which(to.keep), window.size), which(rowSums(is.enzyme.all)>0))
  core.regions = extend.index(which(to.keep), window.size, sides='down', do.unique = F);
  m = match(core.regions, which(rowSums(is.enzyme.all)>0)) # 20141111
  core.regions =  core.regions[!is.na(m)]
  
  # merge clusters
  gene.ranges = unique(cbind(by(core.regions, names(core.regions), FUN = min), by(core.regions, names(core.regions), FUN = max)), MARGIN = 1, drop=F)
  
  gene.ranges = sort.by(gene.ranges, by = gene.ranges[,1])
  nc = nrow(gene.ranges)
  if (is.null(nc) | !nc){
    cat('\nNumber of clusters: 0')    
    return(NULL)  
  }
  cluster.ID = cumsum(c(1, gene.ranges[seq2(2,nc,1),1]- gene.ranges[seq2(1,nc-1,1),2] - 1 > max.dist.merge))
  gene.ranges = data.frame(from=sapply(by(gene.ranges[,1],INDICES = cluster.ID, FUN = function(x){x[1]}), FUN = 'identity'),
                           to = sapply(by(gene.ranges[,2],INDICES = cluster.ID, FUN = function(x){x[length(x)]}), FUN = 'identity'),
                           name= sapply(by(rownames(gene.ranges),INDICES = cluster.ID, FUN = function(x){paste(c(x[1], x[length(x)]), collapse = '_')}), FUN = 'identity'))
  geneID2clusterID = lapply(1:nrow(gene.ranges), function(x){cbind(as.character(anno$ID)[gene.ranges[x,1]:gene.ranges[x,2]], rep(as.character(gene.ranges[x,3]), gene.ranges[x,2]-gene.ranges[x,1]+1))})
  if (is.list(geneID2clusterID)){
    geneID2clusterID = do.call(rbind, geneID2clusterID)
  }
  gene.ranges = cbind(from=as.character(anno$ID)[gene.ranges[,1]], to=as.character(anno$ID)[gene.ranges[,2]], name=as.character(gene.ranges[,3]))
  cat('\nNumber of clusters: ', nrow(gene.ranges))
  
  #### output
  to.keep.extend = extend.index(core.regions, window.extend, sides='both', do.unique=T)
  to.keep.extend = to.keep.extend[to.keep.extend<=length(anno) & to.keep.extend>=1]
  anno$PARENT[1] ==c()
  is.enzyme.all[] = c('', 'Yes')[is.enzyme.all+1]
  out = cbind(chr = as.character(anno@seqnames)[], gene=anno$ID, 'protein ID' = anno@elementMetadata[,toupper(proteinID)], Existing.Anno = anno@elementMetadata[,toupper(desc.fname)],
              is.enzyme.all, domains = ipr.anno)[to.keep.extend,]
  
  rownames(geneID2clusterID) = geneID2clusterID[,1];
  out = cbind(out, clusterID = mat.fill.row(geneID2clusterID, rownames(out), '')[,2])
  
  write.xlsx(out, out.file)
  return(gene.ranges)
  
}
deepAnno.landmarks <- function(landmark.sizes=null, gff.file=NULL, 
                               DNA.fasta.file='/Users/yongli/Universe/data/NPgenome/Aspergillus/A_nidulans_FGSC_A4_current_chromosomes.fasta', 
                               prot.fasta.file=NULL, iprscan.tab.file=NULL, ica.spatial=NULL,  max.dist.merge = 0, 
                               gene.definition = c('gene', 'transcript', 'mRNA'),
                               extra.genes = 20, n.cluster.per.file=40, geMat = NULL,proteinID =c('ID', 'proteinId'),
                               out.file = 'nidulans.deepAnno.llms.xlsx', geneID2cdsID = function(x){sub('gene_(.*)$', 'CDS_\\1', x)}){ #root = '/Users/yongli/Universe/write/Project_Current/t.NPbioinformatics/Nidulans.SlidingWindow/Annotation'){
  # deep annotation around lankmark genes with n.genes on each side
  # YF Li
  # 20141007
  
  prot.seq = read.fasta(prot.fasta.file, type='AA')
  ipr.anno = iprscan.flat(iprscan.tab.file)
  if(!is.null(geMat)&!is.null(gff.file)){
    ica.spatial = express.clustering(gff.file, geMat)    
    anno = ica.spatial$anno;    
  }else if(!is.null(ica.spatial)){
    anno = ica.spatial$anno;        
  }else if(!is.null(gff.file)){
    gff.format = sub('^.*\\.([^\\.]*$)', '\\1', gff.file)
    # anno = read.gff3(gff.file, format=gff.format)
    anno = import.gff(gff.file) # 20160502
    idx.gene = (anno$type==gene.definition)
    anno = anno[idx.gene, ]
    anno = sort.intervals(anno)    
  }else{
    stop('Provide ica.spatial or gff.file')
  }
  # anno$ID = sub('transcript:', '', anno$ID)
  
  
  # get gene annotation and gene orders
  
  # get clusters and merge by distances  
  locs = match(as.character(landmark.sizes[,1]), anno$ID);
  locs = cbind(locs, locs+landmark.sizes[,2]-1)
  gene.ranges = cbind(as.character(landmark.sizes[,1]), anno$ID[locs[,2]])
  gene.ranges = sort.by(gene.ranges, by = locs[,1])
  landmark.sizes = sort.by(landmark.sizes, by = locs[,1])
  locs = sort.by(locs, by = locs[,1])
  nc = nrow(locs)
  cluster.ID = cumsum(c(1, locs[2:nc,1]- locs[1:(nc-1),2] - 1 > max.dist.merge))
  # cbind(gene.ranges, cluster.ID)
  
  # get gene ranges and create cluster names
  if (ncol(landmark.sizes)==3){
    s = cbind(by(as.character(gene.ranges[,1]),INDICES = cluster.ID, FUN = function(x){as.character(x)[1]}), 
              by(as.character(gene.ranges[,2]),INDICES = cluster.ID, FUN = function(x){as.character(x)[length(x)]}),
              by(landmark.sizes[,3], INDICES = cluster.ID, FUN = function(x){paste(as.character(x), collapse = '_')}))
    s = cbind(s[,1:2], paste(s[,3], s[,1], s[,2], sep='_'))    
  }else{
    s = cbind(by(as.character(gene.ranges[,1]),INDICES = cluster.ID, FUN = function(x){as.character(x)[1]}), 
              by(as.character(gene.ranges[,2]),INDICES = cluster.ID, FUN = function(x){as.character(x)[length(x)]}))
    s = cbind(s, paste(s[,1], s[,2], sep='_'))    
  }
  s2d = cluster.deepAnno(ica.spatial = ica.spatial,proteinID =proteinID, gff.file = gff.file, gene.ranges = s, prot.seq=prot.seq, ipr.anno = ipr.anno, out.file = out.file, extra.genes=extra.genes, 
                         DNA.fasta.file = DNA.fasta.file, gene.definition = gene.definition, 
                         n.cluster.per.file=n.cluster.per.file, append=F, geneID2cdsID = geneID2cdsID)
  invisible(s2d)
}

domain.clustering.unsupervised <- function(iprscan.tab.file = iprscan.tab.file, gff.file = gff.file,
                                           window.size = 20){ # enzyme.definition = NULL,
  # find the frequent item set rule of domain annotation combinations among local regions in the genome
  # YF Li
  # 20141007
  
  require(gplots)
  ## read gff
  gff.format = sub('^.*\\.([^\\.]*$)', '\\1', gff.file)
  # anno = read.gff3(gff.file, format=gff.format)
  anno = import.gff(gff.file) # 20160502
  idx.gene = (anno$type=='gene')
  anno = anno[idx.gene, ]
  anno = sort.intervals(anno)
  n = length(anno)
  # read domain annotations
  ipr.anno = iprscan.flat(iprscan.tab.file, out.type = 'itemset')
  ipr.anno.gene = mat.fill.row(t(t(ipr.anno)), row.names = anno$ID, default = '')
  names(ipr.anno.gene) = anno$ID;
  ipr.anno.gene = sapply(ipr.anno.gene, function(x){if (length(x)==1 && x=='') x=c(); return(x)})
  # combine annotation sets to cluster level
  ipr.anno.cluster = sapply(1 : (length(ipr.anno.gene)-window.size+1), FUN = function(x){unique(unlist(ipr.anno.gene[x:(x+window.size-1)]))})
  names(ipr.anno.cluster) = paste(anno$ID[1:(length(ipr.anno.gene)-window.size+1)], anno$ID[window.size:length(ipr.anno.gene)], sep='-')
  
  # create fake gene and clusters
  ipr.anno.gene.rand = relist2(sample(unlist(ipr.anno.gene), replace = F), ipr.anno.gene) # domain level perm
  ipr.anno.gene.perm = sample(ipr.anno.gene) # gene level perm
  ipr.anno.cluster.rand = sapply(1 : (length(ipr.anno.gene)-window.size+1), FUN = function(x){unique(unlist(ipr.anno.gene.perm[x:(x+window.size-1)]))})
  
  # transactions
  ipr.anno.gene.tr = as(ipr.anno.gene, 'transactions')
  ipr.anno.cluster.tr = as(ipr.anno.cluster, 'transactions')
  ipr.anno.gene.rand.tr = as(ipr.anno.gene.rand, 'transactions')
  ipr.anno.cluster.rand.tr = as(ipr.anno.cluster.rand, 'transactions')
  
  # identify eqivalent domains and merge domain concept
  require("arules");
  require("arulesViz")
  data("Adult")
  rules <- apriori(ipr.anno.gene.tr, 
                   parameter = list(supp = 1/500, conf = 0.7,
                                    target = "maximally frequent itemsets"))
  summary(rules)
  
  cluster.rules <- apriori(ipr.anno.cluster.tr, 
                           parameter = list(supp = 1/500, conf = 0.7,
                                            target = "maximally frequent itemsets"))
  summary(cluster.rules)
  # unique coding of domain combinations in each gene
  # domain clustering rule finding
}

proMap2hints <- function(pMap, gff.file='Afu3g01340_Afu3g01340.gff', out.file = 'proMap_hints.gff', 
                         log.file = 'log.txt',
                         geneID2cdsID = function(x){paste(x, '-P', sep='')}, append = F, version = 3){
  # Yong Fuga Li, 20141216
  # create proMap hints gff file based on proMap results and gff file of the genes
  # anno = tryCatch(read.gff3(gff.file, format='gff3'), error = function(e){read.gff3(gff.file, format='gff')}, finally = NULL)
  anno = import.gff(gff.file) # 20160502
  idx.CDS = (anno$type=='CDS')
  anno = anno[idx.CDS, ]
  anno = sort.intervals(anno);   
  i = order(as.character(anno$ID)); anno = anno[i,]; # there are overlapping genes, so sort by gene names to avoid wrong orders, 20121216
  anno = anno[anno$ID %in% geneID2cdsID(names(pMap))]
  out = anno
  out$source = 'proMap2hints'
  out$phase = '.';
  out@elementMetadata$Score = 0;
  out@elementMetadata$coverage = 0;
  out@elementMetadata$mult = 0;
  
  # create introns hints as manual hints
  to.keep = vector('logical', length = length(anno)) | T
  for (g in names(pMap)){
    i = range(which(anno$ID == geneID2cdsID(g)))
    to.keep[i[2]] = F
    i = seq2(i[1], i[2]-1, 1)
    if (!pMap[[g]]$nHits){
      to.keep[i] = F
      next
    }
    if (!length(i))
      next
    if (g == 'AN9494'){
      1
    }
    width = anno@ranges@start[i+1] - (anno@ranges@start[i]+anno@ranges@width[i])
    if (any(width<=0)){
      write(paste('gene', g, 'contains introns with negative size.'), file = log.file, append = T)
      width[width<=0] = 1;
    }
    out@ranges[i] = IRanges(anno@ranges@start[i]+anno@ranges@width[i], width = width)
    i.confident.intron = (pMap[[g]]$intron$coverage * pMap[[g]]$nHits >=5 & pMap[[g]]$intron$coverage > 0.3) | pMap[[g]]$intron$match.score > 3
    out@elementMetadata$mult[i] =  round(pMap[[g]]$intron$coverage * pMap[[g]]$nHits,1)
    out@elementMetadata$coverage[i] =  round(pMap[[g]]$intron$coverage,3)
    out@elementMetadata$Score[i] =  round(pMap[[g]]$intron$match.score,1)
    to.keep[i[!i.confident.intron]] = F
    out$ID[i] = g
  }
  out = out[to.keep]
  out$type = 'intron'
  elementMetadata(out) = data.frame(out@elementMetadata[,c('source', 'type', 'score', 'phase', 'score', 'Score', 'mult', 'coverage')], grp = out@elementMetadata[,'ID'], src = 'M',pri=4)
  if (version == 3){
    out$pri = out$pri - ((out$coverage < 0.3) + (out$mult < 5) + (out$Score < 3)) # 20160818 -- assigned different priorities to evidences of different reliability
  }
  ### convert exon hints to manual hints
  
  export(out, out.file, format = 'gff3', append=append) 
}

get.CDS <- function(gene.IDs = c('AN5093', 'AN5092'), gff.file, DNA.fasta.file, geneID2cdsID=function(x){paste(x, '-P', sep='')}){
  # get the CDS sequences for a list of genes based on genome sequences and gff format file
  # v2. 20150406, retirve all seuqences when gene.IDs = NULL or empty
  require('Biostrings')
  require('rtracklayer')
  require('GenomicFeatures')
  
  
  fa=import(DNA.fasta.file, 'fasta', type='DNA')
  names(fa) <- sub('^([^ ]+) .+$','\\1', names(fa))
  
  # anno = tryCatch(read.gff3(gff.file, format='gff3'), error = function(e){read.gff3(gff.file, format='gff')}, finally = NULL)
  anno = import.gff(gff.file)
  idx.CDS = (anno$type=='CDS')
  anno = anno[idx.CDS, ]
  anno = sort.intervals(anno)
  
  DNAseq = import(DNA.fasta.file, 'fasta', type='DNA')
  CDSs = data.frame('seq'=c(), 'Exons'=c(), 'CDSspan(nt)'=c(), 'from'=c(), 'to'=c())
  
  if (is.null(gene.IDs) | !length(gene.IDs)){ # 20150406
    gene.IDs = anno$ID;
    gene.IDs = unique(gene.IDs[!is.na(gene.IDs)])
    geneID2cdsID = identity;
  }
  
  # i.strand = which(regexpr('strand',names(anno@elementMetadata@listData))>0)
  idx.parent = tolower(colnames(anno@elementMetadata)) == 'parent'
  for (g in gene.IDs){
    # i = anno$ID == geneID2cdsID(g)
    i = (anno$ID == geneID2cdsID(g)) | (as.character(anno@elementMetadata[[which(idx.parent)]]) == g) # 20150819 use the CDS's ID or the parents ID to match
    ranges = data.frame(chr = anno@seqnames[i], anno@ranges[i], ID = g, strands = anno@strand[i])
    gs = getDNA.subseq(DNA.fasta.file, locs = ranges[,1:3])
    seq = paste(gs, collapse='')
    
    # strand = anno@elementMetadata@listData[[i.strand]][i]
    strand = anno@strand[i]
    ustrand = unique(strand)
    if (length(ustrand)>1){
      stop(paste('One both strands', paste(strand, collapse = '')))
    }
    if (ustrand=='-'){
      seq = as.character(reverseComplement(DNAString(seq)))
      # ranges = ranges[rev(1:nrow(ranges)),] # 20141215
    }
    exon.sizes = ranges[,3]-ranges[,2]+1;
    if (ustrand=='-'){
      exon.sizes = rev(exon.sizes)
    }
    CDSs = rbind(CDSs, data.frame ('seq' = seq, 'Exons'=nrow(ranges), 
                                   'CDSspan(nt)' = 1 - ranges[1,2] + ranges[nrow(ranges),3],
                                   'exon.sizes' = paste(exon.sizes, collapse = ','),
                                   'from' = ranges[1,2], 'to' = ranges[nrow(ranges),3],
                                   'chr' = as.character(anno@seqnames[i])[1])) # 20141210, add exon.sizes
  }
  rownames(CDSs) = gene.IDs
  return(CDSs)
}

get.CDS.errorCorrection <- function(g = 'AN2596', DNA.fasta.file='/Users/yongli/Universe/data/NPgenome/Aspergillus/A_nidulans_FGSC_A4_current_chromosomes.fasta',
                                    gff.file="/Users/yongli/Universe/data/NPgenome/Aspergillus/A_nidulans_FGSC_A4_current_features.gff",
                                    new.exon = data.frame(chr=c(), start=c(), end=c()),
                                    error.correction = data.frame(from=c(), to=c())){
  # Yong Fuga Li, 20141014
  fa=import(DNA.fasta.file, 'fasta', type='DNA')
  names(fa) <- sub('^([^ ]+) .+$','\\1', names(fa))
  
  gff.format = sub('^.*\\.([^\\.]*$)', '\\1', gff.file)
  # anno = read.gff3(gff.file, format=gff.format)
  anno = import.gff(gff.file) # 20160502
  idx.CDS = (anno$type=='CDS')
  anno = anno[idx.CDS, ]
  anno = sort.intervals(anno)
  
  DNAseq = import(DNA.fasta.file, 'fasta', type='DNA')
  CDSs = data.frame('seq'=c(), 'Exons'=c(), 'CDSspan(nt)'=c())
  
  i = anno$ID == paste(g, '-P', sep='')
  # i.strand = which(regexpr('strand',names(anno@elementMetadata@listData))>0)
  if (sum(i)>0){
    # strand = anno@elementMetadata@listData[[i.strand]][i] # 20160502
    strand = anno@strand[i]
    ranges = data.frame(chr = anno@seqnames[i], anno@ranges[i], ID = g, strands = anno@strand[i])
    exons = rbind(ranges[,1:3], new.exon);    
  }else{
    strand = new.exon[,'strand']
    exons = rbind(new.exon);
  }
  # exons = ranges[,1:3]
  gs = getDNA.subseq(DNA.fasta.file, locs = exons[,1:3])
  seq = paste(gs, collapse='')
  for (j in seq2(1,nrow(error.correction),1)){
    seq1 = sub(error.correction[j,1], error.correction[j,2], seq)
    if (seq1 == seq)
      stop(paste('Could not replace from', error.correction[j,1], error.correction[j,2]))
    seq = seq1
  }
  
  if (unique(strand)=='-')
    seq = as.character(reverseComplement(DNAString(seq)))
  prot0 = as.character(translate(DNAString(seq),if.fuzzy.codon = 'X'))
  prot1 = as.character(translate(DNAString(seq, start=2),if.fuzzy.codon = 'X'))
  prot2 = as.character(translate(DNAString(seq, start=3),if.fuzzy.codon = 'X'))  
  
  return(list(CDS=seq, exon.seqs = gs, protein=list(frame1=prot0, frame2=prot1, frame3=prot2), exons = exons))
}

plot.genes <- function(gff.file="/Users/yongli/Universe/data/NPgenome/Aspergillus/A_nidulans_FGSC_A4_current_features.gff",
                       start.gene = 'AN2596', end.gene = 'AN2612', tag = 'S700', 
                       out.file = paste(tag, 'NPGC_Struct.pdf', sep=''),
                       cluster.anno = matrix(0,0,0),width=12, height=4, gamma = 0.4, rotation=40, extra.nt = 3000,
                       class2colors = list(oxidoreductase='red', P450='pink', monooxygenase='red', 
                                           hydrolase='orange', aldolase='orange', 
                                           unknown='grey', transporter='blue',other='grey', DUF='black', 
                                           acyltransferases='green', methyltransferase='green', transferase='green')){
  # visaulize gene structure and functions (labels) in a genome regions
  # cluster.anno: a list of gene annotation in the cluster
  #     cluster.anno = data.frame(ID, function, synthesize, class)
  # 'Gviz', 'genoPlotR', 'GenomeGraphs'
  # 'http://genometools.org'
  # YFL, 20141017
  require('Gviz')
  require('grid')
  require('GenomicFeatures')
  
  #   gff.file="/Users/yongli/Universe/data/NPgenome/Aspergillus/A_nidulans_FGSC_A4_current_features.gff"
  #   start.gene = 'AN2596'; end.gene = 'AN2612'
  gff.format = sub('^.*\\.([^\\.]*$)', '\\1', gff.file)
  
  # anno = read.gff3(gff.file, gff.format)
  anno = import.gff(gff.file) # 20160502
  anno = sort.intervals(anno)
  IDs = sub('^([^\\-]*)\\-?.*$', '\\1',anno$ID)
  gene.range.i = c(min(which(IDs==start.gene)), max(which(IDs==end.gene)))
  #   anno$feature = anno$type 
  chr = as.character(anno@seqnames[gene.range.i[1]])
  if (chr != as.character(anno@seqnames[gene.range.i[2]]))
    stop('Two genes in different chromosome')
  
  # anno.sub = cbind(as.data.frame(anno[gene.range.i[1]:gene.range.i[2]]), group=IDs[gene.range.i[1]:gene.range.i[2]])
  anno.sub = as.data.frame(anno[gene.range.i[1]:gene.range.i[2]])
  
  # anno.sub0 = anno[gene.range.i[1]:gene.range.i[2]]; anno.sub0$group=IDs[gene.range.i[1]:gene.range.i[2]]
  anno.gene = anno.sub[anno.sub$type=='gene',]; rownames(anno.gene) = anno.gene$ID
  st.nt = min(anno.gene$start)-extra.nt; en.nt = max(anno.gene$end)+extra.nt
  rownames(cluster.anno) = cluster.anno[,'ID']; cluster.anno = mat.fill.row(cluster.anno, row.names = as.character(anno.gene$ID), default = 'other')
  # tracks
  synthesize.genes = rownames(cluster.anno)[cluster.anno[,'synthesize']=='T']
  axisTrack <- GenomeAxisTrack(range=reduce(IRanges(start=anno.gene[synthesize.genes, 'start'], end = anno.gene[synthesize.genes, 'end'], names = rep('synthesize', length(synthesize.genes)))))
  options(ucscChromosomeNames=FALSE)
  gene.track.color = AnnotationTrack(start = anno.gene$start, width = anno.gene$width, chromosome = chr, 
                                     strand = as.character(anno.gene$strand),
                                     feature =  cluster.anno[as.character(anno.gene$ID),'class'], genome = "Asp. Nidulans", name = 'Genes')    
  gene.track = AnnotationTrack(start = anno.gene$start, width = anno.gene$width, chromosome = chr, 
                               strand = as.character(anno.gene$strand),
                               feature =  cluster.anno[as.character(anno.gene$ID),'function'], genome = "Asp. Nidulans", name = 'Genes')
  txDB = GenomicFeatures::makeTranscriptDbFromGFF(gff.file)
  txTr <- GeneRegionTrack(txDB, chromosome = as.character(anno.gene$seqnames[1]), group = IDs,
                          start = st.nt, end = en.nt, name='Transcripts')
  class2color[setdiff(names(class2color))]
  # plotTracks(axisTrack, from = st.nt, to = en.nt, synthesize='green')
  # feature(txTr)
  if (!is.null(out.file)){
    pdf(out.file, width = width,height = height)
  }
  
  
  plotTracks(c(axisTrack, gene.track, txTr), featureAnnotation = 'feature', just.feature = 'below', 
             fontcolor.feature='#555555', fontsize.feature=10, 
             rotation.item=rotation, transcriptAnnotation = 'group', shape='arrow', just.group = 'left')
  if (1){
    a = 0.15  
    grid.newpage()
    #     pushViewport(viewport(height=a, y=1, just="top"))
    #     grid.rect(gp=gpar(col="grey"))    
    #     plotTracks(axisTrack, add=TRUE, from = st.nt, to = en.nt)  
    #     popViewport(1)
    #     pushViewport(viewport(height=(1-a)*gamma, y=(1-a)*(1-gamma), just="bottom"))  
    #     grid.rect(gp=gpar(col="white"))      
    #     plotTracks(c(gene.track), featureAnnotation = 'feature', fontcolor.feature='#555555', fontsize=10, 
    #                rotation.item=rotation, add=TRUE, from = st.nt, to = en.nt)
    #     popViewport(1)    
    pushViewport(viewport(height=(1-a)*gamma+a, y=(1-a)*(1-gamma), just="bottom"))  
    grid.rect(gp=gpar(col="white"))
    do.call(plotTracks, c(list(trackList= c(axisTrack, gene.track.color), featureAnnotation = 'feature', fontcolor.feature='#555555', fontsize=10, 
                               rotation.item=rotation, add=TRUE, from = st.nt, to = en.nt), class2colors))
    # plotTracks(c(axisTrack, gene.track), featureAnnotation = 'feature', fontcolor.feature='#555555', fontsize=10, 
    #            rotation.item=rotation, add=TRUE, from = st.nt, to = en.nt)
    plotTracks(c(axisTrack, gene.track), featureAnnotation = 'feature', fontcolor.feature='#555555', fontsize=10, 
               rotation.item=rotation, add=TRUE, from = st.nt, to = en.nt)
    
    popViewport(1)    
    
    pushViewport(viewport(height=(1-a)*(1-gamma), y=0, just="bottom"))
    grid.rect(gp=gpar(col="white"))   
    plotTracks(txTr, add=TRUE, transcriptAnnotation = 'group', shape='box', just.group = 'left', from = st.nt, to = en.nt)
    popViewport(1)    
  }
  if (!is.null(out.file)){
    dev.off()    
  }
  #   st = 4463337 - 1000; en = 4543536 + 1000; out.image = 'cluster.png'
  #   system(paste('gt gff3 -tity', gff.file, 'tidy.gff'))
  #   system(paste('gt sketch -seqid', chr,'-start', st, '-end', en, out.image, 'tidy.gff'))
}

codon.optimizer <- function(CDS.list.file='Batch1_CDS.txt', format = c('tab', 'fasta'), CDS=NULL, 
                            N.codon.types.to.change = 6, out.file = NULL, restriction.sites = c(BsaI='GGTCTC', BsaI.rc='GAGACC',
                                                                                                AarI='CACCTGC', AarI.rc = 'GCAGGTG',
                                                                                                polyA8='AAAAAAAA', polyC8 = 'CCCCCCCC',
                                                                                                polyG5='GGGGG', polyT8 = 'TTTTTTTT'),
                            repeats = c(rep.CCAGAG='CCAGAGC'),# provide unites of the repeats
                            tag = '',
                            genetic.table=1, host.species='4932', #left.extra='CCCGGG', right.extra='CCCGGG',
                            left.extra='GATCAGCGGCCGC', right.extra='CCCGGGAACAC'){ # 20141217, use Not1 and XmaI sites
  # V1:
  #   20141023, YF Li
  # V2:
  #   20141210: add format
  #   20141211: fix a bug the one replacement creates another sites, 
  #             and allow using all alternative codons from most frequent to rarest for all codons in a site to be removed
  # V3:
  #   20141230: add repeats removal function by random sampling of codons for the same aa
  # 20151014: change the default of N.codon.types.to.change to 4, in yeast, their are 3 codons (CGG, CGA, CGC) with freq < 3/1000, 
  # and another one (CTC) with freq < 6/100. These 4 also shows the highest codon usage ratios between A. nidulans and S. cerevisiae.
  # V4: to do -- perform codon harmonization
  
  extra.nt = max(c(1,nchar(restriction.sites))) - 1 # extra nt to include on a site to be changes, so that we can check to make sure no new restriction sites are created
  require(RCurl)
  require(xlsx)
  require('Biostrings')
  
  for (i in seq2(1, length(repeats),1)){ # get two unites of the repeats
    repeats[i] = paste(repeats[i], repeats[i], sep='')
  }
  
  if (sub('^.*\\.([^\\.]+)$', '\\1', CDS.list.file) %in% c('fna', 'fasta'))
    format = 'fasta'
  print(paste('Input format', format, sep=''))
  format = match.arg(format);
  if (format=='fasta'){
    fa = read.fasta(fasta.files = CDS.list.file, type = 'DNA')    
    CDS.list.file.bck = CDS.list.file;
    CDS.list.file = paste(CDS.list.file, '.txt', sep='')
    write.table(fa[,'seq', drop=F], file = CDS.list.file, row.names = T, col.names=F, quote = F, sep='\t')
  }
  if (is.null(CDS)){
    CDS = read.table(CDS.list.file,header = F, sep='\t',as.is = T)    
  }
  if (!is.null(CDS.list.file) & is.null(out.file))
    out.file = paste('optimized', tag, '_N',N.codon.types.to.change, '_gt', genetic.table, '_h',host.species, '_', sub('\\..+', '.xls', CDS.list.file), sep='')
  codon.usage.table.url = paste('http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=', host.species, '&aa=', genetic.table, '&style=GCG', sep='');
  a <- getURL(codon.usage.table.url)
  gtable <- read.table(text=regexpr.match('PRE\\>[\\s\\S]*\\.\\.\\n \\n([\\s\\S]*)\\n\\n\\<\\/PRE',a, perl=T)[[1]], header = F, as.is=T,strip.white = T)
  colnames(gtable) = c('AA', 'Codon', 'Number', 'per1000', 'Fraction')
  rownames(gtable) = gtable[,'Codon']
  ii = which.max.by(gtable[,'Number'], gtable[,'AA'])
  gtable.max = gtable[ii,]; rownames(gtable.max) = gtable.max[,'AA']
  gtable$freqNorm = gtable$Number/gtable.max[gtable$AA, 'Number'] # normalized by the most frequent aa
  
  # get rare codons
  rare.codon = sort.by(gtable, by = gtable$freqNorm)[seq2(1, N.codon.types.to.change,1),]
  gtable$toCodon = gtable$Codon;
  gtable[as.character(rare.codon$Codon), 'toCodon'] = gtable.max[as.character(rare.codon$AA), 'Codon']
  rare.codon = sort.by(gtable, by = gtable$freqNorm)[seq2(1, N.codon.types.to.change,1),]  
  cat('rare codons to be changed')
  print(rare.codon)  
  
  # best alternative codons for restriction site optimization
  gtable$bestAlternativeCodon = gtable.max[ gtable$AA, 'Codon']
  ii = which.max.n.by(gtable[,'Number'], gtable[,'AA'], n = 2)
  gtable.max$bestAlternativeCodon = gtable$Codon[ii[gtable.max$AA]]
  gtable.max$bestAlternativeCodon[gtable.max$bestAlternativeCodon %in% rare.codon$Codon] = NA; # if the second best codon is rare than do there is no second best
  gtable.max$bestAlternativeCodon[is.na(gtable.max$bestAlternativeCodon)] = gtable.max$Codon[is.na(gtable.max$bestAlternativeCodon)]
  gtable[as.character(gtable.max$Codon), 'bestAlternativeCodon'] =  gtable.max$bestAlternativeCodon
  gtable[,c('per1000Alt', 'freqNormAlt')] = gtable[as.character(gtable$bestAlternativeCodon),c('per1000','freqNorm')]
  
  # optimize codons
  cat('\noptimizing codons & Removing restriction sites\n')
  colnames(CDS) = c('name', 'CDS')
  n.site.corrected = 0;
  n.protein.corrected = 0;
  for (i in 1:nrow(CDS)){
    cat('\n')
    cat(CDS$name[i],'\t')
    l = nchar(CDS[i,2])
    CDS[i, 'newCDS'] = paste(gtable[substring(CDS[i,2], first = seq(1, l, 3), last = seq(3, l, 3)), 'toCodon'], collapse = '')
    
    has.restriction.site = F
    notes = ''
    for (rr in seq2(1, length(restriction.sites),1)){
      r = restriction.sites[[rr]]; r.name = names(restriction.sites)[rr];
      m <- m00 <- as.matrix(matchPattern(r, DNAString(CDS[i, 'newCDS']), fixed=F)@ranges);
      r = as.character(r)
      m[,2] = m[,1]+ m[,2]-1; m0 = m;
      m[,1] = floor((m[,1]-1)/3)*3+1 # extend to cover whole codons
      m[,2] = ceiling((m[,2])/3)*3 # extend to cover whole codons
      m0 = m0 - m[,1] + 1; # match in the local coordiate
      if (length(m00>0)){
        has.restriction.site = T
        cat(r.name, r, 'site:', nrow(m00), '\t')
        notes = paste(notes, r.name, r, 'site:', nrow(m00), '    ', sep=' ')
      }
      #        if (i == 77){
      #          cat('here we are')
      #        }
      for (j in seq2(1,nrow(m),1)){
        local.seq = substring(CDS[i,'newCDS'], m[j,1], m[j,2]);
        local.seq.left = substring(CDS[i,'newCDS'], m[j,1]-extra.nt, m[j,1]-1); 
        local.seq.right = substring(CDS[i,'newCDS'], m[j,2]+1, m[j,2]+extra.nt); 
        
        ll = nchar(local.seq)
        Cs = substring(local.seq, first = seq(1, ll, 3), last = seq(3, ll, 3)) # codons        
        # to.table = gtable[Cs,]; to.table$index = 1:nrow(to.table)
        # to.table = sort.by(to.table, by = to.table$per1000Alt, decreasing=T);
        # to.table = sort.by(to.table, by = to.table$freqNormAlt, decreasing=T); # sort the codons by 
        # to.try = to.table$index[to.table$bestAlternativeCodon != to.table$Codon]; # index of the codons to change in the orders of codon prefrence
        to.table = c()
        for (tt in 1:length(Cs)){
          indx = gtable$AA == gtable[Cs[tt], 'AA']  & gtable$Codon != Cs[tt]
          if (any(indx))
            to.table = rbind(to.table, cbind(gtable[indx,c('AA', 'Codon', 'freqNorm', 'per1000')], index = tt))
        }
        to.table = sort.by(to.table, by = to.table$freqNorm, decreasing=T); # sort by the codons freq 
        succeed = F;
        for (t in 1:nrow(to.table)){ # to.try){
          CsNew = Cs;
          #CsNew[t] = to.table[CsNew[t], 'bestAlternativeCodon']
          CsNew[to.table$index[t]] = to.table[t, 'Codon']
          local.seq.new = paste(CsNew, collapse = '')
          matched.any = F
          for (r1 in restriction.sites){ # 20141211
            if (length(matchPattern(r1, DNAString(paste(local.seq.left, local.seq.new, local.seq.right, sep='')), fixed=F))){
              matched.any=T
              break
            }
          }
          if (!matched.any){
            succeed = T
            n.site.corrected = n.site.corrected + 1;
            break;
          }
        }
        if (!succeed){cat(CDS$name[i],'\t')
          warning(paste('\nFailed to remove site ', r, ' location ', m00[j,1], ' in sequence ', CDS$name[i], sep=''))
        }else{
          CDS[i, 'newCDS'] = paste(substring(CDS[i, 'newCDS'], 1, m[j,1]-1), local.seq.new, substring(CDS[i, 'newCDS'], m[j,2]+1, l), sep='')
        }        
      }
    }
    
    # handling repeats
    for (rr in seq2(1, length(repeats),1)){
      r = repeats[[rr]]; r.name = names(repeats)[rr];
      m <- m00 <- as.matrix(matchPattern(r, DNAString(CDS[i, 'newCDS']), fixed=F)@ranges);
      r = as.character(r)
      m[,2] = m[,1]+ m[,2]-1; m0 = m;
      m[,1] = floor((m[,1]-1)/3)*3+1 # extend to cover whole codons
      m[,2] = ceiling((m[,2])/3)*3 # extend to cover whole codons
      m0 = m0 - m[,1] + 1; # match in the local coordiate
      if (length(m00>0)){
        has.restriction.site = T
        cat(r.name, r, 'site:', nrow(m00), '\t')
        notes = paste(notes, r.name, r, 'site:', nrow(m00), '    ', sep=' ')
      }
      #        if (i == 77){
      #          cat('here we are')
      #        }
      for (j in seq2(1,nrow(m),1)){
        local.seq = substring(CDS[i,'newCDS'], m[j,1], m[j,2]);
        local.seq.left = substring(CDS[i,'newCDS'], m[j,1]-extra.nt, m[j,1]-1); 
        local.seq.right = substring(CDS[i,'newCDS'], m[j,2]+1, m[j,2]+extra.nt); 
        
        ll = nchar(local.seq)
        Cs = substring(local.seq, first = seq(1, ll, 3), last = seq(3, ll, 3)) # codons        
        # to.table = gtable[Cs,]; to.table$index = 1:nrow(to.table)
        # to.table = sort.by(to.table, by = to.table$per1000Alt, decreasing=T);
        # to.table = sort.by(to.table, by = to.table$freqNormAlt, decreasing=T); # sort the codons by 
        # to.try = to.table$index[to.table$bestAlternativeCodon != to.table$Codon]; # index of the codons to change in the orders of codon prefrence
        to.table = c()
        CsNew = Cs;
        for (reps in 1:10){# try 10 times to get one meet restriction site criteria  
          for (tt in 1:length(Cs)){
            indx = gtable$AA == gtable[Cs[tt], 'AA']
            if (!any(indx))
              next
            candidates = gtable[indx,'Codon'];
            to.use = sample(1:sum(indx),1, T, prob=gtable[indx, 'freqNorm']/sum(gtable[indx, 'freqNorm']))
            CsNew[tt] =  candidates[to.use] 
          }
          local.seq.new = paste(CsNew, collapse = '')
          matched.any = F
          for (r1 in restriction.sites){ # 20141211
            if (length(matchPattern(r1, DNAString(paste(local.seq.left, local.seq.new, local.seq.right, sep='')), fixed=F))){
              matched.any=T
              break
            }
          }
          if (!matched.any){
            succeed = T
            n.site.corrected = n.site.corrected + 1;
            break;
          }
        }
        if (!succeed){cat(CDS$name[i],'\t')
          warning(paste('\nFailed to remove site ', r, ' location ', m00[j,1], ' in sequence ', CDS$name[i], sep=''))
        }else{
          CDS[i, 'newCDS'] = paste(substring(CDS[i, 'newCDS'], 1, m[j,1]-1), local.seq.new, substring(CDS[i, 'newCDS'], m[j,2]+1, l), sep='')
        }        
      }
    }
    
    n.protein.corrected = n.protein.corrected + has.restriction.site
    oldSeq = strsplit(CDS$CDS[i], '')[[1]]; newSeq = strsplit(CDS$newCDS[i], '')[[1]]
    CDS[i, 'Nchanged'] = sum(oldSeq != newSeq)
    CDS[i, 'Nchanged%'] = round(CDS[i, 'Nchanged']/l*100,1)
    CDS[i, 'CG%_old'] = round(sum(oldSeq %in% c('C','G'))/l*100,1)
    CDS[i, 'CG%_new'] = round(sum(newSeq %in% c('C','G'))/l*100,1)
    CDS[i, 'CAI_old'] = exp(mean(log(gtable[substring(CDS[i,'CDS'], first = seq(1, l, 3), last = seq(3, l, 3)), 'freqNorm'])))
    CDS[i, 'CAI_new'] = exp(mean(log(gtable[substring(CDS[i,'newCDS'], first = seq(1, l, 3), last = seq(3, l, 3)), 'freqNorm'])))
    CDS[i, 'sites removed'] = notes
  }
  
  ## confirm protein seq
  cat('\nconfirm protein sequences\n')
  for (i in 1:nrow(CDS)){
    if (translate(DNAString(CDS[i, 'newCDS']),if.fuzzy.codon = 'X') != translate(DNAString(CDS[i, 'CDS']),if.fuzzy.codon = 'X'))
      stop('Protein sequence changed')
  }
  
  cat('n.protein.corrected', n.protein.corrected, '\n')
  cat('n.site.corrected', n.site.corrected, '\n')
  
  CDS$CDSwExtra = paste(left.extra, CDS$newCDS, right.extra, sep='')
  CDS$length = nchar(CDS$CDSwExtra)
  write.table(x = cbind(CDS$name, CDS$CDSwExtra), paste('optimized_',CDS.list.file,sep=''), sep='\t', row.names = F, col.names = F, append = F, quote = F)
  # write.table(x = cbind(CDS$name, CDS$newCDS), paste('new_',CDS.list.file,sep=''), sep='\t', row.names = F, col.names = F, append = F, quote = F)
  # seqlist = cbind(seq = CDS$newCDS, ID = CDS$name);
  seqlist = cbind(seq = CDS$CDSwExtra, ID = CDS$name); # 20160616
  rownames(seqlist) = CDS$name;
  write.fasta(seqlist, paste('optimized_', CDS.list.file.bck, sep=''))
  if (!is.null(out.file))
    write.xlsx(CDS,file = out.file, row.names = F)
  invisible(CDS)
}

predict.genes <- function(genome){
  
}

annotate.functions <- function(genome, genes, do=c('iprscan', 'blastp')){
  # Yong Fuga Li, 
  # 20141027
  
}

FUN.select.KU = function(x){ # 20141126
  y = matrix(F, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x));
  t = as.numeric(sub('(.*)\\%','\\1',as.character(x$Top_nonself_Hit_identity.percent))); t[is.na(t)] = 0
  if ('Existing.Anno' %in% colnames(y))      
    y[, 'Existing.Anno'] = regexpr(pattern='polyketide|alkaloid|terpenoid|terpene|nonribosomal peptide', x$Existing.Anno, ignore.case = T)>0
  if ('domains' %in% colnames(y))            
    y[, 'domains'] = (regexpr(pattern='polyketide', x$domains, ignore.case = T)>0 & as.numeric(as.character(x$length)) > 800) | # polyketide
      regexpr(pattern='alkaloid', x$domains, ignore.case = T)>0 | # alkaloid
      regexpr(pattern='terpenoid|terpene', x$domains, ignore.case = T)>0 | # terpenoid
      regexpr(pattern='nonribosomal peptide', x$domains, ignore.case = T)>0 | (regexpr(pattern='Adenylation|ACP', x$domains, ignore.case = T)>0 & regexpr(pattern='Condensation', x$domains, ignore.case = T)>0 & regexpr(pattern='Phosphopantetheine', x$domains, ignore.case = T)>0) # required domains for NRPS or PKS
  return(y)
}

is.KU <- function(anno.txt='', domain.txt=''){
  # 20150415, Yong Fuga Li
  y = regexpr(pattern='polyketide|alkaloid|terpenoid|terpene|nonribosomal peptide|secondary metabo', anno.txt, ignore.case = T)>0
  y = y | regexpr(pattern='polyketide', domain.txt, ignore.case = T)>0 | # polyketide
    regexpr(pattern='alkaloid', domain.txt, ignore.case = T)>0 | # alkaloid
    regexpr(pattern='terpenoid|terpene', domain.txt, ignore.case = T)>0 | # terpenoid
    regexpr(pattern='nonribosomal peptide', domain.txt, ignore.case = T)>0 | (regexpr(pattern='Adenylation|ACP', domain.txt, ignore.case = T)>0 & regexpr(pattern='Condensation', domain.txt, ignore.case = T)>0 & regexpr(pattern='Phosphopantetheine', domain.txt, ignore.case = T)>0) # required domains for NRPS or PKS
  
}

FUN.select.maybeKU = function(x){ # 20141126
  y = matrix(F, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x));
  t = as.numeric(sub('(.*)\\%','\\1',as.character(x$Top_nonself_Hit_identity.percent))); t[is.na(t)] = 0
  if ('Existing.Anno' %in% colnames(y))      
    y[, 'Existing.Anno'] = regexpr(pattern='polyketide|alkaloid|terpenoid|terpene|nonribosomal peptide', x$Existing.Anno, ignore.case = T)>0
  if ('domains' %in% colnames(y))            
    y[, 'domains'] = regexpr(pattern='(polyketide|acyl carrier protein)', x$domains, ignore.case = T)>0 | # polyketide
      regexpr(pattern='alkaloid', x$domains, ignore.case = T)>0 | # alkaloid
      regexpr(pattern='(terpenoid|terpene|geranyl diphosphate|farnesyl diphosphate)', x$domains, ignore.case = T)>0 | # terpenoid
      regexpr(pattern='nonribosomal peptide', x$domains, ignore.case = T)>0 | (regexpr(pattern='Adenylation|ACP', x$domains, ignore.case = T)>0 & regexpr(pattern='Condensation', x$domains, ignore.case = T)>0 & regexpr(pattern='Phosphopantetheine', x$domains, ignore.case = T)>0) # required domains for NRPS or PKS
  if ('top.5.hits' %in% colnames(y))      
    y[, 'top.5.hits'] = regexpr(pattern='polyketide|alkaloid|terpenoid|terpene|nonribosomal peptide', x$top.5.hits, ignore.case = T)>0
  return(y)
}

FUN.select.promising = function(x){
  y = matrix(F, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x));
  if ('CS' %in% colnames(y)){
    y[,'CS'] = as.numeric(as.character(x[,'CS']))>3      
  }
  return(y)
}
FUN.select.boundary = function(x){
  y = matrix(F, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x));
  if ('express' %in% colnames(y)){
    y[,'express'] = as.numeric(as.character(x$express))>9      
  }
  if ('CS' %in% colnames(y)){
    y[,'CS'] = as.numeric(as.character(x[,'CS'])) < 0.5      
  }  
  return(y)
}

FUN.select.special = function(x){
  y = matrix(F, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x));
  if ('name' %in% colnames(y))      
    y[, 'name'] = regexpr(pattern='llm|laeA', x$name, ignore.case = T)>0
  if ('Existing.Anno' %in% colnames(y))      
    y[, 'Existing.Anno'] = regexpr(pattern='llm|laeA|molyb', x$Existing.Anno, ignore.case = T)>0
  if ('domains' %in% colnames(y))      
    y[, 'domains'] = regexpr(pattern='laeA|molyb', x$domains, ignore.case = T)>0
  return(y)
}

FUN.select.interestingv1 = function(x){
  y = matrix(F, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x));
  t = as.numeric(sub('(.*)\\%','\\1',as.character(x$Top_nonself_Hit_identity.percent))); t[is.na(t)] = 0
  if ('length' %in% colnames(y))
    y[, 'length'] = as.numeric(as.character(x$length)) > 1000
  if ('Top_nonself_Hit_identity.percent' %in% colnames(y))      
    y[, 'Top_nonself_Hit_identity.percent'] = t > 50 | t < 25
  if ('Existing.Anno' %in% colnames(y))      
    y[, 'Existing.Anno'] = regexpr(pattern='polyketide|alkaloid|terpenoid|terpene|nonribosomal peptide', x$Existing.Anno, ignore.case = T)>0
  if ('domains' %in% colnames(y))            
    y[, 'domains'] = regexpr(pattern='polyketide|alkaloid|terpenoid|terpene|nonribosomal peptide', x$domains, ignore.case = T)>0
  if ('top.5.hits' %in% colnames(y))      
    y[, 'top.5.hits'] = regexpr(pattern='polyketide|alkaloid|terpenoid|terpene|nonribosomal peptide', x$top.5.hits, ignore.case = T)>0
  return(y)
}

FUN.select.interesting = function(x){
  y = matrix(F, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x));
  t = as.numeric(sub('(.*)\\%','\\1',as.character(x$Top_nonself_Hit_identity.percent))); t[is.na(t)] = 0
  if ('length' %in% colnames(y))
    y[, 'length'] = as.numeric(as.character(x$length)) > 800
  if ('Top_nonself_Hit_identity.percent' %in% colnames(y))      
    y[, 'Top_nonself_Hit_identity.percent'] = t > 75 | (t < 25 & t > 0)
  return(y)
}


FUN.select.oxidoreductase = function(x){
  y = matrix(F, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x));
  if ('Existing.Anno' %in% colnames(y)) 
    y[, 'Existing.Anno'] = regexpr(pattern='oxidoreductase|P450|oxidase|dehydrogenase|oxygenase|reductase', x$Existing.Anno, ignore.case = T)>0
  if ('domains' %in% colnames(y))  
    y[, 'domains'] = regexpr(pattern='oxidoreductase|P450|oxidase|dehydrogenase|oxygenase|reductase', x$domains, ignore.case = T)>0
  if ('top.5.hits' %in% colnames(y))  
    y[, 'top.5.hits'] = regexpr(pattern='oxidoreductase|P450|oxidase|dehydrogenase|oxygenase|reductase', x$top.5.hits, ignore.case = T)>0
  return(y)
}

FUN.select.boring = function(x){
  y = matrix(F, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x));
  if ('Existing.Anno' %in% colnames(y)) 
    y[, 'Existing.Anno'] = regexpr(pattern='secondary metab', x$Existing.Anno, ignore.case = T)>0
  return(y)
}

FUN.select.warning = function(x){
  y = matrix(F, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x));
  if ('average.intron.size' %in% colnames(y))
    y[, 'average.intron.size'] = as.numeric(as.character(x$average.intron.size)) > 100       
  if ('average.exon.size' %in% colnames(y))
    y[, 'average.exon.size'] = as.numeric(as.character(x$average.exon.size)) < 100       
  return(y)
}

FUN.select.catabolism = function(x){
  y = matrix(F, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x));
  if ('Existing.Anno' %in% colnames(y))
    y[, 'Existing.Anno'] = regexpr(pattern='transferase|synthase|synthetase|ligase', x$Existing.Anno, ignore.case = T)>0
  if ('domains' %in% colnames(y))
    y[, 'domains'] = regexpr(pattern='transferase|synthase|synthetase|ligase', x$domains, ignore.case = T)>0
  if ('top.5.hits' %in% colnames(y))
    y[, 'top.5.hits'] = regexpr(pattern='transferase|synthase|synthetase|ligase', x$top.5.hits, ignore.case = T)>0
  return(y)    
}

xlsx.extractSheet.NPGC <- function(xlsx.file='Pexpansum_MC29w20p0.005_DeepAnno.xlsx', header = T, na.strings = '|',extra.genes = 5){
  # YF Li, 20141126
  # FUN.select: a function to select the cells and return a logical matrix
  options(java.parameters = "-Xmx4g" )
  require(XLConnect)
  
  xlsx.KU.file = sub(pattern = '.([^\\.]+$)', replacement = '_KU.\\1', xlsx.file)
  xlsx.maybeUU.file = sub(pattern = '.([^\\.]+$)', replacement = '_maybeUU.\\1', xlsx.file)
  xlsx.UU.file = sub(pattern = '.([^\\.]+$)', replacement = '_UU.\\1', xlsx.file)
  
  wb <- XLConnect::loadWorkbook(xlsx.file)  
  n = length(XLConnect::getSheets(wb))
  is.KU <- is.maybeUU <- is.UU <- vector(mode = 'logical', length = n)
  jgc()  
  for (i in 1:n){
    x = XLConnect::readWorksheet(wb, sheet = i) # 20141126
    x[x==na.strings] = NA
    is.KU[i] = any(FUN.select.KU(x)[(extra.genes+1):(nrow(x)-extra.genes),], na.rm = T)
    is.maybeUU[i] = any(FUN.select.maybeKU(x)[(extra.genes+1):(nrow(x)-extra.genes),], na.rm = T)
  }
  is.maybeUU = is.maybeUU & ! is.KU  
  is.UU = ! (is.KU | is.maybeUU)    
  cat('Number of KU ', sum(is.KU))
  cat('\nNumber of maybe KU ', sum(is.maybeUU))
  cat('\nNumber of UU ', sum(is.UU))  
  
  extract.sheets <- function(xlsx.file, xlsx.out.file, to.keep, name.from = 'UU', name.to = 'KU'){
    wb = XLConnect::loadWorkbook(xlsx.file)  
    sheetnames = XLConnect::getSheets(wb)  
    XLConnect::clearSheet(wb, sheet = sheetnames[!to.keep])  
    XLConnect::removeSheet(wb, sheet = sheetnames[!to.keep])  
    sheetnames = XLConnect::getSheets(wb)    
    renameSheet(wb, sheet = sheetnames, sub(name.from, name.to, sheetnames))
    XLConnect::saveWorkbook(wb, xlsx.out.file)    
  }
  
  extract.sheets(xlsx.file, xlsx.KU.file, is.KU, name.from = 'UU', name.to = 'KU'); jgc()
  extract.sheets(xlsx.file, xlsx.maybeUU.file, is.maybeUU, name.from = 'UU', name.to = 'maybeUU'); jgc()
  extract.sheets(xlsx.file, xlsx.UU.file, is.UU, name.from = 'UU', name.to = 'UU'); jgc()
  
  xlsx.color.NPGC(xlsx.KU.file)
  xlsx.color.NPGC(xlsx.maybeUU.file)
  xlsx.color.NPGC(xlsx.UU.file)
}

gff.id.change <- function(gff.file, anno = NULL,
                          in.info = c(feature = 'CDS', id.type = 'ID'), in.ids = NULL, 
                          out.info = c(feature = 'gene', id.type = 'ID'),
                          extra.nt = 2500,
                          out.type = c('nt', 'id')){
  # 20151012, YF LI
  out.type =  match.arg(out.type)
  if (out.type == 'nt')
    out.info = c(feature = 'gene', id.type = 'ID')
  if (is.null(anno)){
    gff.format = sub('^.*\\.([^\\.]*$)', '\\1', gff.file)
    # anno = read.gff3(gff.file, format=gff.format)  
    anno = import.gff(gff.file) # 20160502
  }
  anno.in = sort.intervals(anno[anno$type==in.info['feature'], ])
  anno.out = sort.intervals(anno[anno$type==out.info['feature'], ])
  anno.in = anno.in[sort(match(in.ids,  sub('\\.\\d$', '', anno.in@elementMetadata[,in.info['id.type']])))]
  m = findOverlaps(anno.in, subject = anno.out)
  out.ids = anno.out@elementMetadata[m@subjectHits,out.info['id.type']]
  if (out.type == 'nt'){
    locs = geneRanges2ntRanges(anno, out.ids, extra.nt) 
    return(locs)
  }else{
    return(out.ids)    
  }
}


gff.match <- gff.mapping <- function(gff.file = 'cS818_augoHintstop1.gff', gff.reference = 'cS818.gff', tag = 'mapped_', gff.to.out = paste(tag, gff.file, sep=''), 
                                     geneID.re = '^([^\\.]+)(?:\\..+)?$', # extra gene IDs from element (e.g. exon) IDs
                                     match.by = 'gene', format='gff3', geneID2cdsID=identity){
  # Yong Fuga Li, 20141215
  # 20141231: add geneID2cdsID
  
  g = import(gff.file)
  g.ref = import(gff.reference)
  if (!any(match.by %in% g$type)) # 20151002
    stop(paste(g$type, 'is not a feature type for ', gff.file))
  if (!any(match.by %in% g.ref$type))
    stop(paste(g.ref$type, 'is not a feature type for ', gff.reference))
  
  g.slim = g[g$type==match.by]
  g.ref.slim = g.ref[g.ref$type==match.by]
  m = findOverlaps(g.slim, subject = g.ref.slim)
  g.ref.IDs = geneID2cdsID(g.ref.slim$ID[m@subjectHits])
  g.IDs = g.slim$ID[m@queryHits]
  g.IDs.extra = setdiff(g.slim$ID, g.IDs)
  g.IDs = c(g.IDs, g.IDs.extra)
  g.ref.IDs = c(g.ref.IDs, g.IDs.extra)
  ID.map = unlist(as.list(by(g.ref.IDs, regexpr.match(geneID.re, g.IDs), FUN = function(x){paste(unique(x), collapse = '~')})))
  locs = regexpr.match.loc(geneID.re, g$ID)
  for (i in 1: length(locs)){
    if (is.na(g$ID[i]))
      next
    g$ID[i] = paste(substr(g$ID[i], 1, locs[[i]][1,1]-1), ID.map[substr(g$ID[i], locs[[i]][1,1], locs[[i]][1,2])], substr(g$ID[i], locs[[i]][1,2]+1,nchar(g$ID[i])), sep='')
  }
  pa = unlist.multi(g$Parent);
  locs = regexpr.match.loc(geneID.re, pa)
  for (i in 1: length(locs)){
    if (is.na(pa[i]))
      next
    g$Parent[[i]] = paste(substr(pa[i], 1, locs[[i]][1,1]-1), ID.map[substr(pa[i], locs[[i]][1,1], locs[[i]][1,2])], substr(pa[i], locs[[i]][1,2]+1,nchar(pa[i])), sep='')
  }  
  export(g, gff.to.out, format = format)
}


NPGC.wclustering <- enzyme.wclustering <- function(gff.file, iprscan.tab.file = NULL, chromosome.specific=F, 
                                                   anno = NULL,
                                                   ipr.anno = NULL,
                                                   pep.fasta.file=pep.fasta.file, filter.proteins = T, min.protein.length = 150,
                                                   gene.definition = c('gene', 'transcript', 'mRNA'), proteinID = 'ID', 
                                                   annotation.by = c('OR', 'desc', 'domain'), 
                                                   tag = 'A_nidulans_FGSC_A4', window.size = 20, log.scale = F, 
                                                   simu.rep = 5, 
                                                   method = c('TXTwe', 'TXTw', 'MC29e', 'MC29'),
                                                   enzyme.weight.file = paste('/Users/yongli/Universe/write/Project_Current/t.NPbioinformatics/enzyme_weighting/', method, '.txt', sep=''),
                                                   multi.match = c('max', 'mean'),
                                                   prediction.file='Top.Clusters', min.contig.len=window.size/2,
                                                   compare.against =c('simulation','theoretical'),
                                                   p.value.cutoff = 0.005,
                                                   outformat=c('csv', 'tab')){
  # statistical analysis of weighted enzyme clustering in a genome
  # chromosome.specific: estimate chromosome specific enzyme probability estimation
  # simu.rep: simulated gene sequences
  # Yong Fuga Li, 20141220, modified from enzyme.wclustering
  
  compare.against = match.arg(compare.against)
  gene.definition = match.arg(gene.definition)  # 20141125
  outformat = match.arg(outformat)
  annotation.by = match.arg(annotation.by)  # 20141125
  method = match.arg(method);
  require('rtracklayer')
  require('genomeIntervals')
  require(lattice)
  
  if (is.null(anno)){
    gff.format = sub('^.*\\.([^\\.]*$)', '\\1', gff.file)
    # anno = read.gff3(gff.file, format=gff.format)
    anno = import.gff(gff.file) # 20160502
  }
  chrs = as.character(unique(anno@seqnames))
  
  ## keep genes only
  idx.gene = (anno$type==gene.definition) # 20141125
  anno = anno[idx.gene, ]
  anno = sort.intervals(anno)
  colnames(anno@elementMetadata) = toupper(colnames(anno@elementMetadata)) # 20141125
  if (!is.null(anno$NOTE)){ # 20141125
    desc.fname = 'NOTE'
  }else if (!is.null(anno$DESCRIPTION)){
    desc.fname = 'DESCRIPTION'  
  }else{
    warning('No description or Note field for the annotation of genes')      
    desc.fname = 'NOTE'
    anno$NOTE = ''
  }
  
  #   ### remove short chrs, 20141222
  #   to.remove = vector('logical', length(anno))
  #   for (i in 1:length(chrs)){
  #     chr = chrs[i]
  #     is.in.chr = as.vector(anno@seqnames==chr)
  #     if (sum(is.in.chr) < min.contig.len){
  #       to.remove[is.in.chr] = T
  #     }
  #   }
  
  # read ipr anno: 20141125
  if (is.null(ipr.anno)){
    ipr.anno = iprscan.flat(iprscan.tab.file, na.strings = c('-', 'NA', 'NULL'))    
  }
  ipr.anno = mat.fill.row(t(t(ipr.anno)), row.names = anno@elementMetadata[,toupper(proteinID)], default = '')[,1]
  names(ipr.anno) = anno$ID
  if (annotation.by %in% 'desc'){
    annotation.text = as.character(as.vector(anno@elementMetadata[[toupper(desc.fname)]]))    
  }else if(annotation.by %in% 'domain'){
    annotation.text = as.character(as.vector(ipr.anno));
  }else if(annotation.by %in% c('OR')){
    annotation.text = paste(as.character(as.vector(anno@elementMetadata[[toupper(desc.fname)]])), as.character(as.vector(ipr.anno)))    
  }
  
  # filter proteins by length
  cat('removing short unannotated protein\n')
  if(filter.proteins){
    len.pep = nchar(read.fasta(pep.fasta.file)[anno@elementMetadata[,toupper(proteinID)], 'seq']);    
    i.fake.protein = as.character(as.vector(ipr.anno)) == '' & (as.character(as.vector(anno@elementMetadata[[toupper(desc.fname)]])) == '' |
                                                                  as.character(as.vector(anno@elementMetadata[[toupper(desc.fname)]])) == 'Protein of unknown function')
    cat('Number of unannotated protein', sum(i.fake.protein), '\n')
    i.fake.protein = i.fake.protein & (len.pep < min.protein.length)
    anno = anno[!i.fake.protein];
    ipr.anno = ipr.anno[!i.fake.protein];
    annotation.text = annotation.text[!i.fake.protein];
    cat('removed', sum(i.fake.protein), 'proteins\n')
  }
  chrs = as.character(unique(anno@seqnames))
  
  ### read weight file
  w = read.table(file = enzyme.weight.file, header = T, sep = '\t')
  #   pat = paste('(', paste(w$enzyme, collapse = '|'), ')', sep='')
  #   # get gene weights
  #   ma = regexpr.match(pat = pat, txt = annotation.text, perl=T, ignore.case=T)
  cat('Calculating gene weights\n')
  is.enzyme = vector('numeric', length(annotation.text))
  if (multi.match=='max')
    is.enzyme = is.enzyme - Inf;
  nmatches = 0;
  for (p in 1:nrow(w)){
    if (multi.match=='max'){
      ma = c(NA, 1)[(regexpr(w$enzyme[p], text = annotation.text, perl=T)>0)+1];
      is.enzyme = rowMax(cbind(is.enzyme,  ma * w$AverageScore[p]), na.rm=T)    
    }else{
      ma = (regexpr(w$enzyme[p], text = annotation.text, perl=T)>0)
      is.enzyme = is.enzyme + ma * w$AverageScore[p]    
      nmatches = nmatches + ma
    }
  }
  if (multi.match=='mean'){
    is.enzyme = is.enzyme/(1E-10+nmatches)
  }else{
    is.enzyme[is.infinite(is.enzyme)]  = 0
  }
  
  # get local max of enzyme weights in sliding windows
  cat('Identify local max cluster scores\n')
  epsilon = 1E-10
  L.gene = list()  
  chr.ranges = matrix(0, length(chrs), ncol = 2, dimnames = list(chrs,c('start', 'end')))
  labels.succ.local.all = is.enzyme
  for (i in 1:length(chrs)){
    chr = as.character(chrs[i])
    is.in.chr = as.vector(anno@seqnames==chr)
    chr.ranges[chr, c('start', 'end')] =  range(which(is.in.chr)) 
    L.gene[[chr]] = sum(is.in.chr)# number of genes in this chromosome
    seq = is.enzyme[is.in.chr] 
    if (L.gene[[chr]] < window.size){
      labels.succ.local.all[is.in.chr] = 0
      next
    }
    labels.succ.local = label.successes.local.max(seq,window.size) # only keep the local max that are greater than 0
    labels.succ.local[labels.succ.local<0] = 0;
    labels.succ.local.all[is.in.chr] = labels.succ.local
  }
  labels.succ.local.all.wZeros = labels.succ.local.all
  labels.succ.local.all = labels.succ.local.all[labels.succ.local.all>epsilon]
  
  # simulations: get local max of enzyme weights in sliding windows
  labels.succ.local.all.simus = list();
  for (r in 1:simu.rep){
    txt = paste('iteration', r);
    cat(txt)
    labels.succ.local.all.simu = is.enzyme;
    is.enzyme.simu = is.enzyme;
    if (!chromosome.specific){
      is.enzyme.simu = is.enzyme.simu[sample.int(length(is.enzyme.simu))]
    }
    for (i in 1:length(chrs)){
      chr = chrs[i]
      idx = chr.ranges[chr,1]:chr.ranges[chr,2]
      if (L.gene[[chr]] < window.size){
        labels.succ.local.all.simu[idx] = 0;
        next
      }
      if (chromosome.specific){
        seq.simu = is.enzyme[idx]         
        seq.simu = seq.simu[sample.int(length(seq.simu))]
      } else{
        seq.simu = is.enzyme.simu[idx]         
      }
      labels.succ.local = label.successes.local.max(seq.simu, window.size) # only keep the local max that are greater than 0
      labels.succ.local[labels.succ.local<0] = 0;
      labels.succ.local.all.simu[idx] = labels.succ.local
    }
    labels.succ.local.all.simus[[r]] = labels.succ.local.all.simu[labels.succ.local.all.simu>epsilon];
    cat(paste(rep('\b',nchar(txt)), collapse = ''))
  }
  
  #   dat = rbind(data.frame(score = unlist(labels.succ.local.all.simus), data ='simulation'),
  #               data.frame(score = labels.succ.local.all, data = 'real genome'))
  # 
  #   hist.by(dat$score, by = dat$data, by.name = '')
  pdf(paste(tag, '_TrueClustersEstimates.pdf', sep=''), 5,4)
  dd = distribution.diff(sample=labels.succ.local.all, null.samples=labels.succ.local.all.simus, tag = '')
  dev.off()  
  
  ################ output top predictions
  anno.df = as.data.frame(anno)
  for (i in 1:length(anno.df)){
    if (class(anno.df[[i]])!='integer')
      anno.df[[i]] = unlist2(anno.df[[i]])
  }
  anno.df[, 'score'] = labels.succ.local.all.wZeros
  anno.df[, 'p.value'] = dd$score2pvalue(labels.succ.local.all.wZeros)
  anno.df[, 'fdr'] = dd$score2fdr(labels.succ.local.all.wZeros)
  anno.df[, '#true'] = dd$score2ntrue(labels.succ.local.all.wZeros)
  
  # mark the whole clusters
  anno.df[, 'cluster.ID'] = ''  
  l = window.size;
  succ.loc.count = 0;
  gene.ranges = c();
  for (i in which(anno.df[, 'score']>0)){
    succ.loc.count = succ.loc.count+1;
    st = max((i-l+1),chr.ranges[anno.df$seqnames[i],'start'])
    anno.df[st:i, 'score'] = rowMax(cbind(anno.df[st:i, 'score'], anno.df[i, 'score']))
    anno.df[st:i, 'p.value'] = rowMin(cbind(anno.df[st:i, 'p.value'], anno.df[i, 'p.value']))
    anno.df[st:i, 'cluster.ID'] =  paste(anno.df[st:i, 'cluster.ID'], paste('S', succ.loc.count,sep=''))
    gene.ranges = rbind(gene.ranges, c(start=anno.df$ID[st], end=anno.df$ID[i], 
                                       ID = paste('S', succ.loc.count,sep=''), p.value=anno.df[i, 'p.value']))
  }
  gene.ranges = as.data.frame(gene.ranges); gene.ranges$p.value = as.numeric(as.character(gene.ranges$p.value))
  
  # select top window and run clusters
  to.output.windows = anno.df[,'p.value'] < p.value.cutoff; 
  write.table(gene.ranges, file = paste(tag,'_geneRanges_all.tsv', sep=''), sep='\t')
  gene.ranges = gene.ranges[gene.ranges$p.value < p.value.cutoff,]
  write.table(gene.ranges, file = paste(tag,'_geneRanges_filtered.tsv', sep=''), sep='\t')
  
  # how many top clusters are included?
  s.names =  anno.df[to.output.windows, 'cluster.ID']
  s.names = strsplit(paste(s.names,collapse=' '), '\\s+',perl=T)[[1]];
  uc = unique.count(s.names)
  n.clusters.localwindows = sum(uc$counts.unique==window.size)
  
  out.names = c(intersect(c('seqnames', 'start', 'end', 'ID', 'Note', 'orf_classification', 'Gene'),colnames(anno.df)),
                colnames(anno.df)[ncol(anno.df)-5+c(5,1:4)])
  if (outformat=='csv'){
    write.table(anno.df[to.output.windows,out.names], file=paste('cluster.anno.', tag, '.p', p.value.cutoff, '.NWindowClusters',n.clusters.localwindows, '.csv',sep=''),sep=',', row.names=F)
  }else if (outformat=='tab'){
    write.table(anno.df[to.output.windows,out.names], file=paste('cluster.anno.', tag, '.p', p.value.cutoff, '.NWindowClusters',n.clusters.localwindows, '.tab',sep=''),sep='\t', row.names=F, quote = F)
  }
  
  # write clean per cluster output, 20140611
  write.NPGC <- function(anno.df, i.new.NPG = to.output.windows, window.size=window.size, 
                         file.out=paste('cluster.anno.clean', tag, '.p', p.value.cutoff, 
                                        '.NWindowClusters',n.clusters.localwindows, '.tab',sep='')){
    # 20140613
    is.SM = regexpr(pattern='secondary metab', text = as.character(as.vector(anno.df$Note)), perl=T, ignore.case=T)>0
    is.PKS = regexpr(pattern='polyketide synthase', text = as.character(as.vector(anno.df$Note)), perl=T, ignore.case=T)>0
    
    all.SID = anno.df$cluster.ID[i.new.NPG]
    all.SID = strsplit(paste(all.SID,collapse=' '), '\\s+',perl=T)[[1]];
    uc = unique.count(all.SID)
    cluster.names = names(uc$counts.unique[uc$counts.unique==window.size])      
    
    clean.table = matrix('',nrow=length(cluster.names),ncol=8,
                         dimnames=list(cluster.names, c('cluster ID', 'chr', 'coordinate', 'gene range', 'min distance to SM genes', 'closest SM gene(s)', 'p-value', 'cluster gene annotations')));
    n.correct.cluster = 0;
    for (nc in cluster.names){
      i.match = regexpr(paste(nc,'(\\s|$)',sep=''), anno.df$cluster.ID)>0      
      ## get closest SM
      chr = unique(anno.df$seqnames[i.match])
      loc.SM = t(which(is.SM & anno.df$seqnames==chr))
      loc.cluster = t(t(which(i.match)))
      dist.to.SM = repmat(loc.cluster,1,length(loc.SM)) - repmat(loc.SM, length(loc.cluster),1) 
      min.dist.to.SM = min(c(Inf, abs(dist.to.SM)))
      #if (min.dist.to.SM)
      if (!min.dist.to.SM) # 20140720
        n.correct.cluster = n.correct.cluster + 1
      closest.SM = which(abs(dist.to.SM)==min.dist.to.SM,arr.ind=T)
      if (!is.null(closest.SM) && length(closest.SM)>0){
        min.dist.to.SM = paste(dist.to.SM[closest.SM], collapse='...')
        closest.SM = loc.SM[closest.SM[,2]]        
      }
      
      # cluster coordinates
      min.l = min(c(anno.df$start[i.match], anno.df$end[i.match], Inf))
      max.l = max(c(anno.df$start[i.match], anno.df$end[i.match], Inf))
      
      # cluster gene ranges
      first.gene = anno.df$ID[min(which(i.match))]
      last.gene = anno.df$ID[max(which(i.match))]
      
      # cluster all gene annotations;
      cluster.anno = paste(anno.df$ID[i.match], anno.df$Note[i.match], sep='|', collapse='\t')
      matchedSM.anno = paste(anno.df$ID[closest.SM], anno.df$Note[closest.SM], sep='|', collapse='...')
      clean.table[nc, ] = c(nc,chr, paste(min.l, '-', max.l), 
                            paste(first.gene, '-', last.gene), min.dist.to.SM, 
                            matchedSM.anno, min(anno.df[i.match,'p.value']), cluster.anno)
      
    }  
    write(x='#Some of the predicted clusters are overlapping. They may indicate a larger cluster if the clusters significantly overlap (according to the coordiates in column 3).', file=file.out, append=F)
    write(x='#Column 5 gives the distance of the cluster to the closest known secondary metabolite genes', file=file.out, append=T)
    write(x='#Column 5, 0 means known SM genes are within the predicted cluster', file=file.out, append=T)
    write(x='#Column 6 gives the gene names and annotations of the closest SM gene(s)', file=file.out, append=T)
    write(x='#Column 5 and column 6, when there are multiple closest SM genes, they are separated by ...', file=file.out, append=T)
    write(x='#Column 8+ gives the gene names and annotations of the genes in the predicted cluster', file=file.out, append=T)
    write(x=paste('#Estimated No. true NP gene clusters:',dd$n.pos), file=file.out, append=T)
    suppressWarnings(write.table(clean.table, file=file.out,sep='\t', row.names=F, quote = F, append=T))
    # n.SM.cluster = sum((diff(which(is.SM))>1) | (diff.str(anno.df$seqnames[is.SM])))+1
    # number of known SM gene clusters cannot be determined accurately
    out = c(sum(is.SM),sum(is.PKS), sum(i.new.NPG & is.SM), 
            sum(i.new.NPG & is.PKS), n.correct.cluster);
    names(out) = c('#known SM genes', '#known PKSs', '#matched SM genes', '#matched PKS genes',
                   '#matched SM clusters')
    return(out)
  }
  
  a = write.NPGC(anno.df, i.new.NPG = to.output.windows, window.size=window.size,
                 file.out=paste('cluster.annoCompact.', tag, '.p', p.value.cutoff, '.NWindowClusters',n.clusters.localwindows, '.tab',sep=''))
  n.unknowns = sum(regexpr(pattern='Protein of unknown function', text = annotation.text, perl=T)>0) # 20140529
  n.genes = length(anno)
  cat('number of genes', n.genes, '\n')
  cat('number with unknown function', n.unknowns, '\n')
  return(gene.ranges)
}

indicible.promoters <- function(gseID, platformID = ""){
  # load series and platform data from GEO
  gset <- getGEO(gseID, GSEMatrix =TRUE)
  if (length(gset) > 1 & (!is.null(platformID) & platformID != "")) idx <- grep(platformID, attr(gset, "names")) else idx <- 1    
  gset <- gset[[idx]]
  fvarLabels(gset) <- make.names(fvarLabels(gset))
  desc = gset@phenoData@data
  geMat = exprs.gene(gset, ID.type = 'symbol', remove.unmappable = T, coding.only = F)
  colnames(geMat) = desc$title
  write.table(geMat,file = paste(gseID, '.xls', sep=''), sep='\t', col.names = NA)
  SD = rowSds(geMat)
  hist(SD)
  
  # clustering of gene and 
  r = cor.wrap(t(geMat['ADH2',,drop=F]), t(geMat))
  g.inducible = sort(r[,r>0.4],decreasing = T)
  write.table(t(g.inducible), file = paste(gseID, '_ADH2_associated.xls', sep=''), sep='\t')
  pdf(paste(gseID, 'clustering.pdf',sep=''),width = 200,height = 200)
  heatmap.quick.geMat(geMat,col.by.all = T, id.type = 'symbol')
  heatmap.quick.geMat(geMat[g.inducible,],col.by.all = T,sd.cutoff = 0)
  dev.off()
  
  # ICA modules of genes
  
  ### ICA do  
}

get.codon.usage <- function(DNA.file, gff.file, using.existing.file = T){
  # 20120507, Yong Fuga Li
  cu.file = sub('.[^\\.]*$', '.cut', x = gff.file) # codon usage table
  if (!using.existing.file | !file.exists(cu.file)){
    system(paste('getAnnoFasta.pl --seqfile=', DNA.file, ' ',  gff.file, sep=''))
    cds.file = list.files(pattern=sub('.[^\\.]*$', '.*.codingseq', x = gff.file))
    # cds.file = list.files(pattern=paste('Penex1.filtered_proteins.FilteredModels1', '.*\\.codingseq', sep=''))
    cdss = read.fasta(cds.file)
    gcounts = matrix(data = 0,64, nrow(cdss), dimnames = list(codon=names(GENETIC_CODE), gene = rownames(cdss)))
    codes = rownames(gcounts)
    for (i in 1:nrow(cdss)){
      s = toupper(cdss[i,'seq']);
      l = nchar(s);  
      ss = substring(s, first = seq(1, l, 3), last = seq(3, l, 3))
      uc = unique.count(ss)$counts.unique
      gcounts[,i] = mat.fill.row(uc, codes)
    }
    write.table(gcounts,file = cu.file, quote = F, sep = '\t', row.names = T, col.names = NA)    
  }else{
    gcounts = as.matrix(read.table(file = cu.file, header = T, row.names = 1, check.names=F))
  }
  
  return(gcounts)
}

sdf2smiles <- function (sdf) 
{
  require(ChemmineR)
  if (!class(sdf) == "SDFset") {
    stop("reference compound must be a compound of class \"SDFset\"")
  }
  if (1){ #(.haveOB()) {
    sdfstrList = as(as(sdf, "SDFstr"), "list")
    defs = paste(Map(function(x) paste(x, collapse = "\n"), 
                     sdfstrList), collapse = "\n")
    t = Reduce(rbind, strsplit(unlist(strsplit(convertFormat("SDF", 
                                                             "SMI", defs), "\n", fixed = TRUE)), "\t", fixed = TRUE))
    if (class(t) == "character") {
      smiles = t[1]
      names(smiles) = t[2]
    }
    else {
      smiles = t[, 1]
      if (ncol(t)>=2)
        names(smiles) = t[, 2]
    }
    return(smiles)
  }
  else {
    message("ChemmineOB not found, falling back to web service version. This will be much slower")
    sdf2smilesWeb(sdf)
  }
}

best.blast.hits <- function(from.file = 'GCA_000264905.1_Stehi1_protein.faa', 
                            from.gff.file = NULL, 
                            to.file='Stehi1_GeneCatalog_proteins_20101026.aa.fasta', 
                            from.IDs = c('EIM85216',  'EIM85220', '2015KU8'),
                            gene.definition = 'CDS', id.type = 'protein_id'){ # for NCBI gff files
  # 20151001: get the best blast hits for a set of sequences
  # 20151007: add from gff file, so that the internal protein IDs can be extracted accurately from the start and end protein
  fa.from = read.fasta(fasta.files = from.file, type = 'AA')
  proIDs = sub('\\.\\d$', '', rownames(fa.from)) # remove version numbers
  from.IDs = sub('\\.\\d$', '', from.IDs)
  rownames(fa.from) = proIDs
  if (is.null(from.gff.file)){
    i = sort(match(from.IDs[1:2], proIDs))
    all.IDs = proIDs[i[1]:i[2]]      
  }else{
    all.IDs = geneRanges2allGenes(from.gff.file, from.IDs[1:2], gene.definition, id.type)
  }
  all.IDs = sub('\\.\\d$', '', all.IDs)
  
  fa.from.select = paste(from.IDs[3], '.from.fa', sep='')
  write.fasta(fa.from[all.IDs,], fa.from.select)
  blastp.asn.file = paste(from.IDs[3], '.blastp.asn', sep='')
  blastp.xml.file = paste(from.IDs[3], '.blastp.xml', sep='')
  blastp.hitList = paste(from.IDs[3], '.blastp.hits', sep='')
  no.top.hits = 1
  system(paste('blastp -query', fa.from.select, '-num_threads 6 -subject', to.file,  '-outfmt 11 -out', blastp.asn.file, '-evalue 1 -max_target_seqs ', no.top.hits))
  # system(paste('blast_formatter -archive', blastp.asn.file, '-outfmt 5 -out', blastp.xml.file, '-max_target_seqs ', no.top.hits))      
  system(paste('blast_formatter -archive ', blastp.asn.file, ' -outfmt \'6 qseqid sseqid length pident mismatch gaps\'  -out', blastp.hitList, '-max_target_seqs 1'))  
  hits = read.table(blastp.hitList, header = F, comment.char = '')
  colnames(hits) = c('qseqid', 'sseqid', 'length', 'pident', 'mismatch', 'gaps')
  hits = as.data.frame(lapply(hits, FUN = function(x){unlist(as.list(by(x, paste(hits$qseqid, hits$sseqid, sep=':'), # 20151016, handle the case where there are more than 1 hits per target protein
                                                                        function(y){
                                                                          if (is.character(y) | is.factor(y))
                                                                            return(unique(as.character(y)))
                                                                          else if(is.numeric(y))
                                                                            return(sum(y))
                                                                          else
                                                                            return(y[1])})))}))
  rownames(hits)  = hits$qseqid
  hits = hits[all.IDs,] # 20151023: keep original gene orders
  hits$pident = 1-(hits$mismatch+hits$gaps)/hits$length
  hits$paln = hits$length/sapply(fa.from[all.IDs,1], nchar)[rownames(hits)] # percentage of query sequences aligned
  write.table(hits, blastp.hitList, sep='\t', row.names = F, quote = F)
  return(hits)
}

lineage.map <- function(queries = c('Aspergillus','Penicillium', 'Epichloe', 'Fusarium'), Rank = 'class', SubTree = 'Fungi'){
  # convert a list of taxonomy queries to a specific rank
  # Yong Fuga Li, 20151106
  queries.old = queries;
  queries = unique(queries)
  out.rank = ncbiTaxonomy(paste(SubTree, '[SubTree] AND ', Rank, '[Rank]', sep=''), FALSE)
  lineage = ncbiTaxonomy(queries, FALSE)
  a = unlist.dupenames(sapply(out.rank$name, function(x) grep(x,lineage$lineage)))
  out = cbind(from = lineage$name[a], to = names(a))
  rownames(out) = out[,'from']
  out = mat.fill.row(out, queries.old)
  return(out)
}
taxonomy.lineage.overlap.v1 <- function(species1, species2){
  # 20151006, YF Li
  require('genomes')
  species1 = unique(species1)
  species2 = unique(species2)  
  lineage1 = ncbiTaxonomy(species1, FALSE)
  lineage2 = ncbiTaxonomy(species2, FALSE)
  n1 = nrow(lineage1); n2 = nrow(lineage2)
  overlaps = matrix(1, nrow = n1, ncol= n2, dimnames = list(lineage1$name, lineage2$name))
  common = matrix('', nrow = n1, ncol= n2, dimnames = list(lineage1$name, lineage2$name))
  for (i in 1:n1){
    for (j in 1:n2){
      common[i,j] = substr(lineage1$lineage[i], 1, lcprefix(lineage1$lineage[i], lineage2$lineage[j]))
      overlaps[i,j] = length(strsplit(common[i,j], split = '; ')[[1]])
    }
  }
  return(list(common.lineage = common, overlaps = overlaps))
}

taxonomy.lineage.overlap <- function(species1, species2=NULL){
  # 20151006, YF Li
  # 20150324: allow species2 to be null
  require('genomes')
  Sys.setenv(email='yonli@umail.iu.edu')
  species1 = sort(unique(species1))
  #lineage1 = ncbiTaxonomy(species1, FALSE)
  # lineage1 = t(sapply(species1, function(x) ncbiTaxonomy(x, FALSE)))
  lineage1 = t(sapply(species1, function(x){y = ncbiTaxonomy(x, FALSE); y[5] = paste(y[5],'; ',x, sep=''); return(y)})) # 20160324, the lineage does not contain the last level, adding it to obtain better species similarity matrix
  
  if (all(is.null(species2))){
    species2 = species1;
    lineage2 = lineage1;
  }else{
    species2 = sort(unique(species2)) 
    #lineage2 = ncbiTaxonomy(species2, FALSE)
    lineage2 = t(sapply(species2, function(x) {y = ncbiTaxonomy(x, FALSE); y[5] = paste(y[5],'; ',x, sep=''); return(y)})) # 20160324,
  }
  
  n1 = nrow(lineage1); n2 = nrow(lineage2)
  overlaps = matrix(1, nrow = n1, ncol= n2, dimnames = list(species1, species2))
  common = matrix('', nrow = n1, ncol= n2, dimnames = list(species1, species2))
  Len1 = vector('numeric', length = nrow(lineage1))
  Len2 = vector('numeric', length = nrow(lineage2))
  for (i in 1:n1){
    for (j in 1:n2){
      if (0){
        common[i,j] = substr(lineage1[i,'lineage'], 1, lcprefix(lineage1[i,'lineage'][[1]], lineage2[j,'lineage'][[1]]))
        overlaps[i,j] = length(strsplit(common[i,j], split = '; ')[[1]])
      }else{ # 20160324
        L1 = strsplit(lineage1[i,'lineage'][[1]], split = '; ')[[1]]
        L2 = strsplit(lineage2[j,'lineage'][[1]], split = '; ')[[1]]
        Len1[i] = length(L1); Len2[j] = length(L2);
        l = min(Len1[i], Len2[j])
        overlaps[i,j] = min(which(c(L1[1:l]!=L2[1:l], T)))-1
        common[i,j] = paste(L1[seq2(1, overlaps[i,j],1)], collapse = '; ')
      }
    }
  }
  # 20160324 & 20160527-- compute normalized similarity
  # max.overlaps = diag(1/sqrt(diag(overlaps)))
  sim = diag(1/sqrt(Len1), nrow = length(Len1)) %*% overlaps %*% diag(1/sqrt(Len2), nrow = length(Len2));
  rownames(sim) <- rownames(overlaps);
  colnames(sim) <- colnames(overlaps);
  
  return(list(common.lineage = common, overlaps = overlaps, similarity = sim))
}
select.ModelSpecies.v1 <- function(query.species){
  # 20151006, YF Li
  # version 1, do so maually using NCBI common Tree functionality
  # ref: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3901538/
  # ref: http://www.ncbi.nlm.nih.gov/books/NBK179288/
  # ref: http://www.chnosz.net/manual/taxonomy.html
  require('taxize')
  require('genomes')
  q.lineage = ncbiTaxonomy(query.species, FALSE)
  # system(paste('epost -db ', from.db, ' -id ', paste(uids[((k-1)*max.query+1):(min(k*max.query,length(uids)))], collapse = ','), '| ', ' elink -target ', to.db, ' -cmd neighbor | xtract -pattern LinkSet -element Id > ', out.file, sep=''))   # -block Stat -match @type:PubMed -element  @count        
  
  query.genus =  sub('^(\\S+)\\s.*$', '\\1',query.species) # sapply(strsplit(query.species, ' '), rbind)[1,]
  # Cryptococcus neoformans gattii ==> Cryptococcus gattii
  # Cryptococcus neoformans neoformans ==> Cryptococcus neoformans
  model.species = strsplit('Homo sapiens, Drosophila melanogaster, Arabidopsis thaliana, Brugia malayi, Aedes aegypti, Tribolium castaneum, Schistosoma mansoni, Tetrahymena thermophila, Galdieria sulphuraria, Zea mays, Toxoplasma gondii, Caenorhabditis elegans, Caenorhabditis elegans , Aspergillus fumigatus, Aspergillus nidulans, Aspergillus nidulans, Aspergillus oryzae, Aspergillus terreus, Botrytis cinerea, Candida albicans, Candida guilliermondii, Candida tropicalis, Chaetomium globosum, Coccidioides immitis, Coprinus cinereus, Coprinus cinereus, Nicotiana attenuata, Cryptococcus gattii, Cryptococcus neoformans, Debaryomyces hansenii, Encephalitozoon cuniculi, Eremothecium gossypii, Fusarium graminearum, Fusarium graminearum, Histoplasma capsulatum, Histoplasma capsulatum, Kluyveromyces lactis, Laccaria bicolor, Petromyzon marinus, Leishmania tarentolae, Lodderomyces elongisporus, Magnaporthe grisea, Neurospora crassa, Neurospora crassa, Phanerochaete chrysosporium, Phanerochaete chrysosporium, Pichia stipitis, Rhizopus oryzae, Saccharomyces cerevisiae, Saccharomyces cerevisiae, Saccharomyces cerevisiae, Schizosaccharomyces pombe, Thermoanaerobacter tengcongensis, Trichinella spiralis, Ustilago maydis, Ustilago maydis, Yarrowia lipolytica, Nasonia vitripennis, Solanum lycopersicum, Chlamydomonas reinhardtii, Amphimedon queenslandica, Pneumocystis jirovecii, Triticum aestivum, Gallus gallus, Danio rerio, Escherichia coli, Staphylococcus aureus', split = ', ')[[1]]
  model.genus = sub('^(\\S+)\\s.*$', '\\1',model.species)
  cat('\nOverlapping with model species:')
  cat(intersect(tolower(query.species), tolower(model.species)));
  cat('\nOverlapping with model genus:')
  cat(intersect(toupper(query.genus), toupper(model.genus)));
  cat('\n')
  q.species = paste(c(query.species, unique(model.species)), collapse = ' OR ')
  q.genus = paste(c(query.genus, unique(model.genus)), collapse = ' OR ')
  
  URL1 = paste('http://www.ncbi.nlm.nih.gov/taxonomy/?term=', q.species, sep='')
  URL2 = paste('http://www.ncbi.nlm.nih.gov/taxonomy/?term=', q.genus, sep='')
  return(list(URL.species = URL1, URL.genus = URL2))
}

select.ModelSpecies <- function(query.species, simplify = T){
  # 20151006, YF Li
  # version 2: do so automatically based on the largest linear overlaps
  # ref: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3901538/
  # ref: http://www.ncbi.nlm.nih.gov/books/NBK179288/
  # ref: http://www.chnosz.net/manual/taxonomy.html
  require('taxize')
  require('genomes')
  query.species = unique(query.species)
  query.genus =  sub('^(\\S+)\\s.*$', '\\1',query.species) # sapply(strsplit(query.species, ' '), rbind)[1,]
  query.genus = unique(query.genus)
  # Cryptococcus neoformans gattii ==> Cryptococcus gattii
  # Cryptococcus neoformans neoformans ==> Cryptococcus neoformans
  model.species = strsplit('Homo sapiens, Drosophila melanogaster, Arabidopsis thaliana, Brugia malayi, Aedes aegypti, Tribolium castaneum, Schistosoma mansoni, Tetrahymena thermophila, Galdieria sulphuraria, Zea mays, Toxoplasma gondii, Caenorhabditis elegans, Aspergillus fumigatus, Aspergillus nidulans, Aspergillus oryzae, Aspergillus terreus, Botrytis cinerea, Candida albicans, Candida guilliermondii, Candida tropicalis, Chaetomium globosum, Coccidioides immitis, Coprinus cinereus, Coprinus cinereus, Nicotiana attenuata, Cryptococcus gattii, Cryptococcus neoformans, Debaryomyces hansenii, Encephalitozoon cuniculi, Eremothecium gossypii, Fusarium graminearum, Fusarium graminearum, Histoplasma capsulatum, Histoplasma capsulatum, Kluyveromyces lactis, Laccaria bicolor, Petromyzon marinus, Leishmania tarentolae, Lodderomyces elongisporus, Magnaporthe grisea, Neurospora crassa, Neurospora crassa, Phanerochaete chrysosporium, Phanerochaete chrysosporium, Pichia stipitis, Rhizopus oryzae, Saccharomyces cerevisiae, Saccharomyces cerevisiae, Saccharomyces cerevisiae, Schizosaccharomyces pombe, Thermoanaerobacter tengcongensis, Trichinella spiralis, Ustilago maydis, Ustilago maydis, Yarrowia lipolytica, Nasonia vitripennis, Solanum lycopersicum, Chlamydomonas reinhardtii, Amphimedon queenslandica, Pneumocystis jirovecii, Triticum aestivum, Gallus gallus, Danio rerio, Escherichia coli, Staphylococcus aureus', split = ', ')[[1]]    
  if (simplify){
    model.species = strsplit('Homo sapiens, Drosophila melanogaster, Arabidopsis thaliana, Brugia malayi, Aedes aegypti, Tribolium castaneum, Schistosoma mansoni, Tetrahymena thermophila, Galdieria sulphuraria, Zea mays, Toxoplasma gondii, Caenorhabditis elegans, Aspergillus nidulans, Botrytis cinerea, Candida albicans, Candida guilliermondii, Candida tropicalis, Chaetomium globosum, Coccidioides immitis, Coprinus cinereus, Coprinus cinereus, Nicotiana attenuata, Cryptococcus gattii, Cryptococcus neoformans, Debaryomyces hansenii, Encephalitozoon cuniculi, Eremothecium gossypii, Fusarium graminearum, Fusarium graminearum, Histoplasma capsulatum, Histoplasma capsulatum, Kluyveromyces lactis, Laccaria bicolor, Petromyzon marinus, Leishmania tarentolae, Lodderomyces elongisporus, Magnaporthe grisea, Neurospora crassa, Neurospora crassa, Phanerochaete chrysosporium, Phanerochaete chrysosporium, Pichia stipitis, Rhizopus oryzae, Saccharomyces cerevisiae, Saccharomyces cerevisiae, Saccharomyces cerevisiae, Schizosaccharomyces pombe, Thermoanaerobacter tengcongensis, Trichinella spiralis, Ustilago maydis, Ustilago maydis, Yarrowia lipolytica, Nasonia vitripennis, Solanum lycopersicum, Chlamydomonas reinhardtii, Amphimedon queenslandica, Pneumocystis jirovecii, Triticum aestivum, Gallus gallus, Danio rerio, Escherichia coli, Staphylococcus aureus', split = ', ')[[1]]      
  }else{
  }
  model.species = unique(model.species)
  model.genus = sub('^(\\S+)\\s.*$', '\\1',model.species)
  model.genus = unique(model.genus)
  cat('\nOverlapping with model species:')
  cat(intersect(tolower(query.species), tolower(model.species)));
  cat('\nOverlapping with model genus:')
  cat(intersect(toupper(query.genus), toupper(model.genus)));
  cat('\n')
  overlap.mat = taxonomy.lineage.overlap(query.species, model.species)
  matches = cbind(# query.species = rownames(overlap.mat$overlaps), 
    best.model.species = colnames(overlap.mat$overlaps)[apply(overlap.mat$overlaps, 1, FUN = which.max)],
    overlap = apply(overlap.mat$overlaps, 1, FUN = max))
  best.model.species = apply(overlap.mat$overlaps, 1, FUN = function(x)names(which(x==max(x))))
  if (is.matrix(best.model.species)){
    best.model.species = unlist.dupenames(as.list(as.data.frame(best.model.species, stringsAsFactors = F)))
  }
  all.matches = cbind(best.model.species = unlist.dupenames(best.model.species),
                      overlap = unlist.dupenames(apply(overlap.mat$overlaps, 1, FUN = function(x) x[x==max(x)])))
  return(list(matches, all.matches))
}

### HGT detection using phymltest (PhyML in ape): http://www.inside-r.org/packages/cran/ape/docs/phymltest
### ref: http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003035
### To confirm the non-metazoan origin of the sequences with hU≥30 and at least one significant metazoan hit, each transcript meeting these conditions was translated and aligned using ClustalW2 to the output (the best hits for each of the five taxa) of the previous blastx analysis. Each alignment was then trimmed to exclude regions where only one of the sequences was present, and phylogenetic trees were built in PhyML from amino-acids sequences using a JTT model [15];

gff.addgene <- function(gff.file, gene.definition = 'CDS', ID = 'proteinId', out.file){
  # 20151110, add gene features to gff files
  # Yong Fuga Li
  # anno = read.gff3(meta$gff.file[i])
  anno = import.gff(gff.file) # 20160502
  anno.c = anno[anno$type==gene.definition]
  
  a = tapply(anno.c, anno.c$proteinId, FUN = function(x){y = range(x);
  y@elementMetadata = x@elementMetadata[1,];
  return(y)}, simplify=T)
  anno.g = c()
  for (i in a){
    if (!length(anno.g))
      anno.g = i
    else
      anno.g = c(anno.g, i)
  }
  anno.g$type = 'gene'
  anno.g$phase = '.'
  anno.g$exonNumber = NA
  anno = c(anno, anno.g)
  anno = sort.intervals(anno)
  export.gff3(anno, out.file)
}

NPscanv0 <- function(gff.files = NULL, iprscan.tab.files = NULL){
  # HMM based learning and prediction of NPGC based on domain annotaiton information
  require(HMM) # ref: http://web.stanford.edu/class/stats366/hmmR2.html;
  require(hmm.discnp)
  require(Rhmm)
  # require(NHMM)
  require(G1DBN)
  # require(HiddenMarkov)
  # require('ebdbNet')
  # dynamic naive bayes; dynamic bayesian network, generalized HMM
  # http://stackoverflow.com/questions/17696547/hidden-markov-models-package-in-r
  # http://a-little-book-of-r-for-bioinformatics.readthedocs.org/en/latest/src/chapter10.html
  # Dec.2015, Jan. 2016
  setwd('/Users/yongli/Universe/write/Project_Current/9.O.NPbioinformatics/AllFungalGenomes/Aspergillus_Binary')
  gff.file = 'Aspzo1.filtered_proteins.GeneCatalog_2013_12_03_15_38.gff3' # Aspergillus zonatus v1.0
  iprscan = 'Aspzo1_GeneCatalog_proteins_20121010_IPR.tab'
  if (0){
    setwd('/Users/yongli/Universe/write/Project_Current/9.O.NPbioinformatics/AllFungalGenomes/Aspergillus_Binary')
    gff.file = 'Aspnid1.filtered_proteins.AspGD_genes.gff3' # Aspergillus nidulans from AspGD
    iprscan = 'Aspnid1_GeneCatalog_proteins_20110130_IPR.tab'
  }
  if (0){
    setwd('/Users/yongli/Universe/write/Project_Current/9.O.NPbioinformatics/AllFungalGenomes/Pleurotus/')
    gff.file = 'PleosPC15_2.filtered_proteins.FilteredModels1.gff3' # Pleurotus ostreatus PC15 v2.0
    iprscan = 'PleosPC15_2_domaininfo_FilteredModels1.tab'
    signalP = 'PleosPC15_2_signalp_FilteredModels1.tab'
    ECpathway = 'PleosPC15_2_ecpathwayinfo_FilteredModels1.tab'
  }
  
  # anno = read.gff3(gff.file)
  anno = import.gff(gff.file) # 20160502
  ipr = iprscan.flat(iprscan, out.type = 'table')
  n.hidden = 2;
  
  
}

plot.hmm <- function(hmm = training.trace[[4]]$hmm){
  # 20160203
  require(lattice)
  require(latticeExtra)
  require(gridExtra)
  require(gplots)
  color <- colorRampPalette(c('white','red'))(256) 
  x = log(hmm$transProbs); colnames(x) = NULL
  a = levelplot(t(x), col.regions=color,aspect = 'fill',
                ylab = 'From', xlab = 'To') #,scales = list(x = list(rot = 90),alternating=1))
  x = log(hmm$emissionProbs);
  colnames(x) = NULL
  b = levelplot(t(x), ylab = '', aspect = 'fill',#scales = list(x = list(rot = 90),alternating=1),
                xlab = 'Domains', col.regions=color)
  # Combination via `c.trellis`
  d = grid.arrange(a,b, nrow=1, ncol=2,widths = c(0.45,0.8), heights = 1)
  # print(d)
  #   layout_matrix = cbind(c(1,2,2,2))
  #   comb_levObj <- c(a, b, layout = c(2, 1), merge.legends = T)
  #   print(comb_levObj)
}


NPscanv1 <- function(species = 'Aspergillus nidulans from AspGD', gff.file=NULL, iprscan.tab.file = NULL, 
                     bin.file = NULL, out.file = 'A.nidu.xls',
                     domain.order = F, pseudocount = 0.5,
                     data.root = '/Users/yongli/Dropbox/NPGC/NPGCquery_data',
                     gene.definition = c('gene', 'transcript'), proteinID = 'ID'){
  # find a window of size 15 or less that meet the gene function query criteria
  # YF Li
  # 20141028, 20141111
  # 20150625-27: fixed a bug that occurs when no cluster is found, changed the interface and added species, bin.file
  #              separate the gff and iprscan file parsing and the querying code
  require(seqHMM)
  require(HMM) # ref: http://web.stanford.edu/class/stats366/hmmR2.html;
  require(HiddenMarkov)
  require(depmixS4)
  require(CRF)
  require('R2HTML')
  require('xlsx')
  root = getwd()
  setwd(data.root)
  if (!is.null(species)){
    meta = read.table('NPGCquery_meta.txt',header = T,as.is = T, sep= '\t', row.names = 1)
    bin.file = meta[species, 'bin.file']
    gff.file=NULL; iprscan.tab.file = NULL;
    proteinID = meta[species, 'proteinID']
    gene.definition = meta[species, 'gene.definition']
  }
  if (!is.null(bin.file)){
    load(bin.file)
  }else{
    bin.file = paste('NPGCquery', gff.file, '.RData', sep='')
    get.NPGC.query.bin(gff.file=gff.file, iprscan.tab.file = iprscan.tab.file, bin.file = bin.file,
                       gene.definition = gene.definition, proteinID = proteinID)
    load(bin.file)
  }
  if (0){
  domains = length(unique(ipr.tab$ipr.acc[ipr.tab$analysis %in% c('HMMPfam', 'HMMSmart', 'HMMTigr', 'HMMPIR', 'superfamily', 'BlastProDom')]))
  motifs = length(unique(ipr.tab$ipr.acc[ipr.tab$analysis %in% c('ProfileScan', 'FPrintScan', 'ScanRegExp')]))
  }
  # to.keep = ipr.tab$analysis %in% c('HMMPfam', 'HMMSmart', 'HMMTigr', 'HMMPIR', 'superfamily', 'BlastProDom')
  to.keep = ipr.tab$ipr.acc != ''
  ipr.tab = ipr.tab[to.keep,]
  emit.symbols = sort(unique(ipr.tab$ipr.acc))
  
  ##### prepare training data
  # ignor domain orders
  if (domain.order)
    stop('domain.order not implemented yet')
  seqs = c()
  accs = by(ipr.tab, ipr.tab$ID, FUN = function(x)unique(x$ipr.acc))
  nDomPerGene = sapply(accs[anno@elementMetadata[,toupper(proteinID)]], FUN = length)
  hist(nDomPerGene, xlab = '#Domains per Gene')
  for (chr in unique(anno@seqnames)){
    for (g in anno@elementMetadata[anno@seqnames == chr,toupper(proteinID)]){
      seqs = c(seqs, 'b', accs[[g]], 'e')   # add 'b' and 'e' to indicate protein start and end
    }
  }
  domCounts = unique.count(seqs)$counts.unique[emit.symbols]
  hist(log(domCounts), xlab = 'log(#Genes per Domain)')
  
  ################
  ### initialize model
  ################
  # define enzymes
  domain.anno = paste(ipr.tab$ipr.desc, ipr.tab$signature.desc, sep='~~~');
  domain.anno = by(domain.anno, INDICES = ipr.tab$ipr.acc, FUN = function(x){paste(unique(x), collapse = ';')})
  
  is.enzyme.EC6 = regexpr(pattern='(oxidoreductase|transferase|hydrolase|lyase|isomerase|ligase)', text =paste(ipr.tab$ipr.desc, ipr.tab$signature.desc), perl=T, ignore.case=T) > 0
  is.enzyme.EC6 = by(is.enzyme.EC6, ipr.tab$ipr.acc, FUN = any)
  is.enzyme.EC6 = is.enzyme.EC6[emit.symbols]
  is.enzyme.MC29e = regexpr(pattern='(oxidoreductase|hydrolase|dehydrogenase|synthase|reductase|transferase|methyltransferase|oxidase|synthetase|monooxygenase|isomerase|dehydratase|decarboxylase|deaminase|O\\-methyltransferase|transaminase|hydratase|acetyltransferase|N\\-acetyltransferase|dioxygenase|aminotransferase|O\\-acyltransferase|esterase|N\\-methyltransferase|acyltransferase|aldolase|O\\-acetyltransferase|cyclase|catalase|hydroxylase|P450|transporter|transcription factor)', text =paste(ipr.tab$ipr.desc, ipr.tab$signature.desc), perl=T, ignore.case=T) > 0 
  is.enzyme.MC29e = by(is.enzyme.MC29e, ipr.tab$ipr.acc, FUN = any)
  is.enzyme.MC29e = is.enzyme.MC29e[emit.symbols]
  is.enzyme = is.enzyme.MC29e
  
  NPG.initialProfile = domCounts * is.enzyme
  NPG.initialProfile = (NPG.initialProfile + pseudocount)/sum(NPG.initialProfile+pseudocount) 
  initialProfile = (domCounts + pseudocount)/sum(domCounts+pseudocount)
  nH = 6; nE = length(emit.symbols);
  HMM = initHMM(States = c("NPG.b","NPG.d","NPG.e","OG.b","OG.d","OG.e"), 
                Symbols = c('b','e',emit.symbols), 
                startProbs = c(.25,0,0.25,0.25,0,0.25),
                transProbs = rbind(t(c(0,0.9,0.1,0,0,0)),
                                   t(c(0,0.8,0.2,0,0,0)),
                                   t(c(0.9,0,0,0.1,0,0)),
                                   t(c(0,0,0,0,0.9,0.1)),
                                   t(c(0,0,0,0,0.8,0.2)),
                                   t(c(0.1,0,0,0.9,0,0))),
                emissionProbs = rbind(t(c(1, rep(0,nE+1))),
                                      t(c(0,0,NPG.initialProfile)),
                                      t(c(0,1,rep(0, nE))),
                                      t(c(1, rep(0, nE+1))),
                                      t(c(0,0,initialProfile)),
                                      t(c(0,1,rep(0, nE)))))
  # HMM.ori
  # HMM.dom
  # HMM.dom.ori
  
  ################
  # unsupervised learning
  ################
  step = 2; max.steps = 50; delta=1E-5
  training.trace = list()
  training.trace[['0']] = list(hmm = HMM, difference = Inf)
  s = 1
  while (1){
    cat((s-1)*step, training.trace[[paste(s-1)]]$difference)
    # training.trace[[paste(s)]] = baumWelch(training.trace[[paste(s-1)]]$hmm, observation=seqs, maxIterations=step, delta=delta, pseudoCount=0)
    training.trace[[paste(s)]] = viterbiTraining(training.trace[[paste(s-1)]]$hmm, observation=seqs, maxIterations=step, delta=delta, pseudoCount=pseudocount)
    if (all(training.trace[[paste(s)]]$difference < delta))
      break
    cat('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b')
    s = s+1
  }
  require('lattice')
  plot.hmm(training.trace[[1]]$hmm)
  plot.hmm(training.trace[[s]]$hmm)
  top.NPG.domains = sort(training.trace[[4]]$hmm$emissionProbs['NPG.d',], decreasing = T)[1:100]
  top.NPG.domains =  cbind(domain.anno[names(top.NPG.domains)], top.NPG.domains)
  top.OG.domains = sort(training.trace[[4]]$hmm$emissionProbs['OG.d',], decreasing = T)[1:100]
  top.OG.domains =  cbind(domain.anno[names(top.OG.domains)], top.OG.domains)
  top.diff.domains = sort(training.trace[[4]]$hmm$emissionProbs['NPG.d',]/training.trace[[4]]$hmm$emissionProbs['OG.d',], decreasing = T)[1:100]
  top.diff.domains =  cbind(domain.anno[names(top.diff.domains)], top.diff.domains)
  step = 2; max.steps = 50; delta=1E-9
  training.traceBW = list()
  training.traceBW[['0']] = training.trace[[length(training.trace)]]
  s = 1
  while (1){
    cat((s-1)*step, training.traceBW[[paste(s-1)]]$difference)
    training.traceBW[[paste(s)]] = baumWelch(training.traceBW[[paste(s-1)]]$hmm, observation=seqs, maxIterations=step, delta=delta, pseudoCount=0.5)
    if (all(training.traceBW[[paste(s)]]$difference < delta))
      break
    cat('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b')
    s = s+1
  }
  posterior(hmm, observation)
  viterbi(hmm, observation)
  image(vt$hmm$transProbs)
  image(HMM$transProbs)
  require(gplots)
  image(log(vt$hmm$emissionProbs),col=greenred(256))
  image(log(HMM$emissionProbs),col=greenred(256))
  
  # label KUs
  # semisupervised EM learning
  # Predictions
  
  #### output
  cat('\n#Identified clusters: ', nrow(gene.ranges))
  to.keep.extend = extend.index(core.regions, window.extend, sides='both', do.unique=T)
  to.keep.extend = to.keep.extend[to.keep.extend<=length(anno) & to.keep.extend>=1]
  anno$PARENT[1] ==c()
  is.enzyme.all[] = c('', 'Yes')[is.enzyme.all+1]
  out = cbind(chr = as.character(anno@seqnames)[], gene=anno$ID, 'protein ID' = anno@elementMetadata[,toupper(proteinID)], Existing.Anno = anno@elementMetadata[,toupper(desc.fname)],
              is.enzyme.all, domains = ipr.anno)[to.keep.extend,]
  
  rownames(geneID2clusterID) = geneID2clusterID[,1];
  out = cbind(out, clusterID = mat.fill.row(geneID2clusterID, rownames(out), '')[,2])
  
  write.xlsx(out, out.file)
  # HTML(out, 'out.html')
  return(out)
}

NPscan <- function(genome.ID = c('A.nidu_AspGD', 'Aspergillus flavus NRRL3357',
                               'Aspergillus tubingensis v1.0', 'Penicillium expansum',
                               'Trichoderma virens'),
                   nNPG= 2, # number of NPG cluster types;
                   nOG = 1, # number of other gene cluster types;
                   init.truth = c('IPR',  'MC29e', 'EC6', 'note', 'domain', 'both'),
                   init.expand = 1, # add the neighbors of inital genes as potential SM genes
                   init.decay.rate = 1-1/8, init.expand.combine = c('prob', 'max'), # max, use max when combine two probability diffused to a location
                   init.weighted = F,
                   predict.extra = F, # output the prediction results for initial HMM and final HMM removing pseudocounts
                   eval.truth = c('note', 'domain', 'both', 'IPR'),
                   do.viterbi = T, remove.pseudocounts = F,
                   EM.method = c('null', 'depmixS4', 'seqHMM', 'HMM'), # null: no EM step is done
                   domain.order = F, pseudocount = 0.1, annotation.by = 'OR',
                   data.root = '/Users/yongli/Dropbox/NPGC/NPGCquery_data',
                   dom.freq.cutoff = 2, # only the emitted domains with frequency above this are retained for training
                   gff.file=NULL, iprscan.tab.file = NULL, 
                   bin.file = NULL, out.tag =sub(' ', '_', genome.ID), 
                   out.file = paste(out.tag, '.xls', sep=''), remove.glycan=T,
                   gene.definition = c('gene', 'transcript'), proteinID = 'ID'){
  # YF Li
  # 20160202-20160204
  # v2: HMM is slow using seqHMM instead
  # v3: 20160413, add init.truth, evaluation.truth, 
  # 'IPR+/-5': mark proteins based on IPR domains and then expand to +/- 5 genes as SM genes
  # v4. 20160421: allow multiple chromosome, and multi-genomes
  
  require('R.matlab')
  require(seqHMM)
  require(HMM) # ref: http://web.stanford.edu/class/stats366/hmmR2.html;
  #   require(HiddenMarkov)
  require(depmixS4) # depmixS4 is super fast compared to HMM and seqHMM, which are similar in speed, 
  # although HMM provide viterbi training, which is order of magnitude faster
  # require(CRF)
  require('R2HTML')
  require('xlsx')
  root = getwd()
  EM.method = match.arg(EM.method)
  eval.truth = match.arg(eval.truth)
  init.truth = match.arg(init.truth)
  init.expand.combine = match.arg(init.expand.combine)
  
  ################
  ###Prepare data
  ################
  if (0){ # compared to the JGI version, AspGD gff contains gene annotation
    genome.ID = 'A.nidu_AspGD'
    setwd('/Users/yongli/Universe/write/Project_Current/9.O.NPbioinformatics/Nidulans.SlidingWindow/Annotation')
    DNA.file='A_nidulans_FGSC_A4_current_chromosomes.fasta'
    gff.file="A_nidulans_FGSC_A4_current_features.gff"
    pep.fasta.file = "A_nidulans_FGSC_A4_current_orf_trans_all.fasta"
    iprscan.tab.file = 'A_nidulans_FGSC_A4_iprscan.out.txt';
    proteinID = 'ID'; cds2gene = function(x)(sub('-T$', '', x, perl=T))
  }
  
  setwd(data.root)
  if (!is.null(genome.ID)){
    meta = read.table('NPGCquery_meta.txt',header = T,as.is = T, sep= '\t', row.names = 1)
    bin.file = meta[genome.ID, 'bin.file']
    gff.file=NULL; iprscan.tab.file = NULL;
    proteinID = meta[genome.ID, 'proteinID']
    gene.definition = meta[genome.ID, 'gene.definition']
  }
  ipr.anno.all = c(); ipr.tab.all = data.frame()
  if (!is.null(bin.file)){
    first = 1;
    for (b in bin.file){
      cat('loading', b, '\n')
      load(b)
      if (first){
        first = 0;
        anno.all = anno;
      }else{
        anno.all = c(anno.all, anno);
      }
      ipr.anno.all = c(ipr.anno.all, ipr.anno);
      ipr.tab.all = rbind(ipr.tab.all, ipr.tab[c('ID', 'analysis', "ipr.acc", "ipr.desc", "signature.acc", "signature.desc")])
    }
    
  }else{
    if (length(genome.ID)>1)
      stop('Only allow one genome if binary files is not used!\n')
    bin.file = paste(genome.ID, '.RData', sep='')
    get.NPGC.query.bin(gff.file=gff.file, iprscan.tab.file = iprscan.tab.file, bin.file = bin.file,
                       gene.definition = 'gene', proteinID = proteinID)
    load(bin.file)
  }
  
  chr = as.character(anno@seqnames)[]; gene.ID = anno$ID; prot.ID = anno@elementMetadata[,toupper(proteinID)];
  anno.txt = unlist(anno@elementMetadata@listData[[toupper(desc.fname)]])
  domain.txt = as.character(as.vector(ipr.anno))
  if (annotation.by %in% 'desc'){
    annotation.text = anno.txt
  }else if(annotation.by %in% 'domain'){
    annotation.text = domain.txt;
  }else if(annotation.by %in% c('OR')){
    annotation.text = paste(anno.txt, domain.txt)    
  }
  names(annotation.text) = names(ipr.anno)
  if (eval.truth == 'note'){ # 20160413
    t.sm.evaluation = is.KU(anno.txt)
  }else if (eval.truth == 'domain'){
    t.sm.evaluation = is.KU(domain.txt =  domain.txt)
  }else if (eval.truth == 'both'){
    t.sm.evaluation = is.KU(anno.txt, domain.txt)
  }
  
  if (init.truth == 'note'){ # 20160413
    is.ku = is.KU(anno.txt)
  }else if (init.truth == 'domain'){
    is.ku = is.KU(domain.txt =  domain.txt)
  }else if (init.truth == 'both'){
    is.ku = is.KU(anno.txt, domain.txt)
  }
  
  if (0){
    domains = length(unique(ipr.tab$ipr.acc[ipr.tab$analysis %in% c('HMMPfam', 'HMMSmart', 'HMMTigr', 'HMMPIR', 'superfamily', 'BlastProDom')]))
    motifs = length(unique(ipr.tab$ipr.acc[ipr.tab$analysis %in% c('ProfileScan', 'FPrintScan', 'ScanRegExp', 'ProSitePatterns', 'ProSiteProfiles', 'PRINTS', 'ScanRegExp')]))
    locations = c('TMHMM', 'SignalP_EUK', 'SignalP_GRAM_NEGATIVE', 'SignalP_GRAM_POSITIVE')
    to.keep = ipr.tab$analysis %in% c('HMMPfam', 'HMMSmart', 'HMMTigr', 'HMMPIR', 'superfamily', 'BlastProDom')
    ipr.tab = ipr.tab[to.keep,]
  }
  to.keep = toupper(ipr.tab$analysis) %in% toupper(c('HMMPfam', 'HMMSmart', 'HMMTigr', 'HMMPIR', 'superfamily', 'BlastProDom', 'Pfam', 'SMART', 'TIGRFAM','ProDom','PIRSF', 'Hamap', 'Gene3D'))
  to.keep = to.keep & ipr.tab$ipr.acc != ''
  ipr.tab = ipr.tab[to.keep,]

  ################
  ##### remove low frequency domains; 20160404
  ################
  accs0 = by(ipr.tab, ipr.tab$ID, FUN = function(x){y = unique(x$ipr.acc); y[!is.na(y)]})
  accs0 = do.call(list, accs0)
  domCounts0 = unique.count(unlist(accs0))$counts.unique
  domCounts0 = domCounts0[domCounts0 >= dom.freq.cutoff]
  doms.tokeep = names(domCounts0)
  accs = by(ipr.tab, ipr.tab$ID, FUN = function(x){y = unique(x$ipr.acc); y[!is.na(y) & y %in% doms.tokeep]})
  accs = do.call(list, accs)
  # accs[setdiff(names(ipr.anno), names(accs))] = NA
  # accs = accs[names(ipr.anno)]
  idx = !is.na(ipr.tab$ipr.acc)
  emit.symbols = sort(doms.tokeep)
  # emit.symbols = emit.symbols[!is.na(doms.tokeep)]
  
  ### domain based NPG labeling, 20160413
  sm.ipr <- function(ipr.acc=NULL, file = '/Users/yongli/Dropbox/NPGC/NPGCquery_data/SM.domains_manualAnno_v2.xlsx'){
    SM.doms = read.xlsx2(file,sheetIndex = 1, as.is =T)
    rownames(SM.doms) = SM.doms$id
    SM.doms$BGC.confidence = as.numeric(as.character(SM.doms$BGC.confidence))
    # sum(SM.doms$BGC.confidence>0)
    if (0){
      ipr.acc.confident = as.character(SM.doms$id[SM.doms$BGC.confidence==1])
      ipr.acc.maybe = as.character(SM.doms$id[SM.doms$BGC.confidence==0.5])
    }
    if (is.null(ipr.acc)){
      ipr.acc = SM.doms$id
    }
    x = SM.doms[ipr.acc,'BGC.confidence'];
    x[is.na(x)] = 0;
    names(x) = ipr.acc
    return(x)
  }
  
  ################
  ##### obtaining domain dependencies by association rule mining, 20160310
  ################
  if (0){
    require(arules)
    require("arules");
    require("arulesViz")
    data("Adult")
    accs.tr = as(accs, 'transactions')
    rules <- apriori(accs.tr, 
                     parameter = list(support=1/1000, conf = 1,minlen = 2, maxlen= 2,
                                      target = 'rules')) # "maximally frequent itemsets"))
    summary(rules)
    inspect(rules)
    # rules@items@itemsetInfo
    rules1 <- subset(rules, subset = lift > 2)
  }
  
  ################
  ##### prepare training data
  ################
  # ignor domain orders
  warning('Some domains has not ipr.acc')
  nDomPerGene = sapply(accs[anno@elementMetadata[,toupper(proteinID)]], FUN = length)
  hist(nDomPerGene, xlab = '#Domains per Gene')
  if (domain.order)
    stop('domain.order not implemented yet')
  seqs = c()
  seqs.nr = c() # non-redundant domain annotations
  geneseqs = c()
  for (chr in unique(anno@seqnames)){
    for (g in anno@elementMetadata[anno@seqnames == chr,toupper(proteinID)]){
      seqs = c(seqs, 'b', accs[[g]], 'e')   # add 'b' and 'e' to indicate protein start and end
      geneseqs = c(geneseqs, rep(g, length(accs[[g]])+2))
    }
  }
  emit.symbols = emit.symbols[emit.symbols%in%unique(seqs)] # some in iprscan file are not in gff file...
  domCounts = unique.count(seqs)$counts.unique[emit.symbols]
  
  pdf('Domain_prevalence.pdf',5,4)
  hist(log(domCounts), xlab = 'log(#Genes per Domain)', main ='')
  dev.off()
  
  ################
  ### initialize model
  ################
  # define enzymes
  domain.anno = paste(ipr.tab$ipr.desc, ipr.tab$signature.desc, sep='~~~');
  domain.anno = by(domain.anno, INDICES = ipr.tab$ipr.acc, FUN = function(x){paste(unique(x), collapse = ';')})
  
  if (remove.glycan){ # 20160311
    head(ipr.tab)
    is.glycan = regexpr(pattern='(glyco|galacto|fructo|gluco)', text =paste(ipr.tab$ipr.desc, ipr.tab$signature.desc), perl=T, ignore.case=T) > 0 
    is.othercontaminants = regexpr(pattern='(kinase|proteasome)', text =paste(ipr.tab$ipr.desc, ipr.tab$signature.desc), perl=T, ignore.case=T) > 0 
    # is.enzyme.KUnotGlycan = (ipr.tab$ID %in% anno@elementMetadata$ID[is.ku]) & !is.glycan & !is.othercontaminants
    # is.enzyme.KUnotGlycan = by(is.enzyme.KUnotGlycan[idx], ipr.tab$ipr.acc[idx], FUN = any)
    # is.enzyme.KUnotGlycan = is.enzyme.KUnotGlycan[emit.symbols]
    is.glycanContam = by((is.glycan | is.othercontaminants)[idx], ipr.tab$ipr.acc[idx], FUN = any)
    is.glycanContam = is.glycanContam[emit.symbols]
    is.glycan = by(is.glycan[idx], ipr.tab$ipr.acc[idx], FUN = any)
    is.glycan = is.glycan[emit.symbols]
  }else{
    is.glycan <- is.glycanContam <- zeros(length(emit.symbols))
  }
  

  is.SM = sm.ipr(emit.symbols); # domain based SM gene annotation
  is.SM = max.by(is.SM[ipr.tab$ipr.acc], ipr.tab$ID, min = 0)[ipr.tab$ID] # keep weights
  if (init.truth == 'IPR'){
    is.enzyme = sm.ipr(emit.symbols);
    if (1){ # expand from domains to genes, 20160414
      if (init.weighted){
        is.enzyme = max.by(is.enzyme[ipr.tab$ipr.acc], ipr.tab$ID, min = 0)[ipr.tab$ID] # keep weights
      }else{
        is.enzyme = ipr.tab$ID %in% ipr.tab$ID[ipr.tab$ipr.acc %in% names(which(is.enzyme>0))]
      }
    }
  }else if (init.truth == 'MC29e'){
    is.enzyme = regexpr(pattern='(oxidoreductase|hydrolase|dehydrogenase|synthase|reductase|transferase|methyltransferase|oxidase|synthetase|monooxygenase|isomerase|dehydratase|decarboxylase|deaminase|O\\-methyltransferase|transaminase|hydratase|acetyltransferase|N\\-acetyltransferase|dioxygenase|aminotransferase|O\\-acyltransferase|esterase|N\\-methyltransferase|acyltransferase|aldolase|O\\-acetyltransferase|cyclase|catalase|hydroxylase|P450|transporter|transcription factor)', text =paste(ipr.tab$ipr.desc, ipr.tab$signature.desc), perl=T, ignore.case=T) > 0 
    if (1){ # expand from domains to genes, 20160414
      is.enzyme = ipr.tab$ID %in% ipr.tab$ID[is.enzyme]
    }
  }else if (init.truth == 'EC6'){
    is.enzyme = regexpr(pattern='(oxidoreductase|transferase|hydrolase|lyase|isomerase|ligase)', text =paste(ipr.tab$ipr.desc, ipr.tab$signature.desc), perl=T, ignore.case=T) > 0
    if (1){ # expand from domains to genes, 20160414
      is.enzyme = ipr.tab$ID %in% ipr.tab$ID[is.enzyme]
    }
  }else if (init.truth %in% c('note', 'domain', 'both')){ # based on KU and SM specific keywords
    is.enzyme = ipr.tab$ID %in% anno@elementMetadata$ID[is.ku]
  }else{
  }
  if (eval.truth=='IPR'){# expand from domains to genes, 20160414
    t.sm.evaluation = sm.ipr(emit.symbols);
    t.sm.evaluation = vector.fill(max.by(t.sm.evaluation[ipr.tab$ipr.acc], ipr.tab$ID, min = 0) > 0, prot.ID)
  }
  
  if (init.expand>0){ # expand to neighbor genes
    is.enzyme = vector.fill(max.by(is.enzyme, ipr.tab$ID, min = 0), prot.ID);is.enzyme[is.na(is.enzyme)] = 0 # to gene level
    if (init.expand.combine == 'prob'){
      is.enzyme = diffuse.by(is.enzyme, init.expand, decay.rate = init.decay.rate, combine.fun = function(x,y)1-(1-x)*(1-y)) # diffuse to neighbor genes
    }else if (init.expand.combine=='max'){
      is.enzyme = diffuse.by(is.enzyme, init.expand, decay.rate = init.decay.rate, combine.fun = max2) # diffuse to neighbor genes
    }
    is.enzyme = is.enzyme[ipr.tab$ID] # to domains in the genes
    is.enzyme[is.na(is.enzyme)] = 0
  }
  #is.enzyme = is.enzyme.MC29e; tag = 'MC29e'
  #is.enzyme = is.enzyme.EC6; tag = 'EC6'
  
  if (init.weighted){
    is.enzyme = by(is.enzyme[idx], ipr.tab$ipr.acc[idx], FUN = sum) # use max to keep the weights if provided, as in 'IPR' method for init.truth
    is.enzyme = is.enzyme[emit.symbols]
    is.enzyme[is.na(is.enzyme)] = 0
    NPG.initialProfile = is.enzyme
    NPG.initialProfile = (NPG.initialProfile + pseudocount)/sum(NPG.initialProfile+pseudocount) 
    NPG.initialProfile = NPG.initialProfile[emit.symbols];
    is.enzyme = is.enzyme > 0 # use binary for remaining initialization
  }else{
    is.enzyme = by(is.enzyme[idx], ipr.tab$ipr.acc[idx], FUN = max) # use max to keep the weights if provided, as in 'IPR' method for init.truth
    is.enzyme = is.enzyme[emit.symbols]
    is.enzyme[is.na(is.enzyme)] = 0
    is.enzyme = is.enzyme > 0 # use binary for all initialization
    NPG.initialProfile = domCounts * is.enzyme
    NPG.initialProfile = (NPG.initialProfile + pseudocount)/sum(NPG.initialProfile+pseudocount) 
    NPG.initialProfile = NPG.initialProfile[emit.symbols];
  }
  
  NPG.initialProfile.noGlycan = domCounts * (is.enzyme & !is.glycanContam)
  NPG.initialProfile.noGlycan = (NPG.initialProfile.noGlycan + pseudocount)/sum(NPG.initialProfile.noGlycan+pseudocount) 
  NPG.initialProfile.noGlycan = NPG.initialProfile.noGlycan[emit.symbols];
  
  NPG.initialProfile.glycan = domCounts * is.glycan
  NPG.initialProfile.glycan = (NPG.initialProfile.glycan + pseudocount)/sum(NPG.initialProfile.glycan+pseudocount) 
  NPG.initialProfile.glycan = NPG.initialProfile.glycan[emit.symbols];
  
  initialProfile = (domCounts + pseudocount)/sum(domCounts+pseudocount)
  initialProfile = initialProfile[emit.symbols];
  
  nTypes = nNPG + nOG
  
  out.tag = paste(out.tag, 'it',toupper(init.truth),'ie',init.expand,'iw',c('F','T')[init.weighted+1], 'et', toupper(eval.truth), 
                  'iec', toupper(init.expand.combine), 'idr', signif(init.decay.rate,2),
                  'rc', c('F','T')[remove.glycan+1],'rp', c('F','T')[remove.pseudocounts+1],
                  'npg', nNPG, 'og', nOG,  'domfr',
                  dom.freq.cutoff, 'pc',pseudocount, sep='')
  nH = nTypes * 3; nE = length(emit.symbols);
  States = c()
  emissionProbs = c()
  for (i in 1:nNPG){
    States = c(States, paste('NPG', i, c('.b', '.d', '.e'), sep=''))
    if (i == 1){
      NPG.profile = (NPG.initialProfile.noGlycan + runif(nE)*mean(NPG.initialProfile.noGlycan) *0.5) * exp(rnorm(nE)*0.1)
    }else if (i==2){
      NPG.profile = (NPG.initialProfile.glycan + runif(nE)*mean(NPG.initialProfile.glycan) *0.5) * exp(rnorm(nE)*0.1)
    }else{
      NPG.profile = (NPG.initialProfile + runif(nE)*mean(NPG.initialProfile) *0.5) * exp(rnorm(nE)*0.1)
    }
    NPG.profile = NPG.profile/sum(NPG.profile)
    emissionProbs = rbind(emissionProbs, t(c(1, rep(0,nE+1))),
                          t(c(0,0,NPG.profile)),
                          t(c(0,1,rep(0, nE))))
  }
  for (i in 1:nOG){
    States = c(States, paste('OG', i, c('.b', '.d', '.e'), sep=''))
    OG.profile = (initialProfile + runif(nE)*mean(initialProfile) *0.5) * exp(rnorm(nE)*0.1)
    OG.profile = OG.profile/sum(OG.profile)
    emissionProbs = rbind(emissionProbs, 
                          t(c(1, rep(0, nE+1))),
                          t(c(0,0,OG.profile)),
                          t(c(0,1,rep(0, nE))))
  }
  startProbs = c()
  transProbs = c()
  p.i = 0.9 # intra-state transition probability
  for (i in 1:nTypes){
    startProbs = c(startProbs, c(1/nTypes, 0, 0))
    block.intra = rbind(t(c(0,0.9,0.1)),t(c(0,0.8,0.2)),t(c(p.i,0,0)))
    block.inter = rbind(t(c(0,0,0)),t(c(0,0,0)),t(c((1-p.i)/(nTypes-1),0,0)))
    transProbs1 = c()
    for (j in 1: nTypes){
      if (i == j){
        transProbs1 = cbind(transProbs1, block.intra)
      }else{
        transProbs1 = cbind(transProbs1, block.inter)
      }
    }
    transProbs = rbind(transProbs, transProbs1)
  }
  rownames(transProbs) <- colnames(transProbs) <- States
  
  if(0){#
    nH = 6; nE = length(emit.symbols);
    HMM = initHMM(States = c("NPG.b","NPG.d","NPG.e","OG.b","OG.d","OG.e"), 
                  Symbols = c('b','e',emit.symbols), 
                  startProbs = c(.25,0,0.25,0.25,0,0.25),
                  transProbs = rbind(t(c(0,0.9,0.1,0,0,0)),
                                     t(c(0,0.8,0.2,0,0,0)),
                                     t(c(0.9,0,0,0.1,0,0)),
                                     t(c(0,0,0,0,0.9,0.1)),
                                     t(c(0,0,0,0,0.8,0.2)),
                                     t(c(0.1,0,0,0.9,0,0))),
                  emissionProbs = rbind(t(c(1, rep(0,nE+1))),
                                        t(c(0,0,NPG.initialProfile)),
                                        t(c(0,1,rep(0, nE))),
                                        t(c(1, rep(0, nE+1))),
                                        t(c(0,0,initialProfile)),
                                        t(c(0,1,rep(0, nE)))))
  }
  
  ################
  ### viterbi training
  ################
  HMM0 = initHMM(States = States, 
                 Symbols = c('b','e',emit.symbols), 
                 startProbs = startProbs,
                 transProbs = transProbs,
                 emissionProbs = emissionProbs)
  if (do.viterbi){
    step = 2; max.steps = 50; delta=1E-5
    training.trace = list()
    training.trace[['0']] = list(hmm = HMM0, difference = Inf)
    s = 1
    ptm <- proc.time();
    while (1){
      cat((s-1)*step, training.trace[[paste(s-1)]]$difference)
      # training.trace[[paste(s)]] = baumWelch(training.trace[[paste(s-1)]]$hmm, observation=seqs, maxIterations=step, delta=delta, pseudoCount=0)
      training.trace[[paste(s)]] = viterbiTraining(training.trace[[paste(s-1)]]$hmm, observation=seqs, maxIterations=step, delta=delta, pseudoCount=0.1)
      if (all(training.trace[[paste(s)]]$difference < delta) | s > max.steps) # | 
          # abs(training.trace[[paste(s)]]$difference - training.trace[[paste(s-1)]]$difference)[2] < delta |
          # s > 15) # no longer improves
        break
      cat('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b')
      s = s+1
    }
    print(proc.time()-ptm)
    s = s+1
    # write.csv(anno@elementMetadata[is.ku,], 'KU.csv')
    # write.csv(ipr.anno[is.ku], 'KU_ipr.csv')
    ### prediction
    statePred = HMM::viterbi(training.trace[[s]]$hmm, seqs)
    postP = HMM::posterior(training.trace[[s]]$hmm, seqs)
    output.NPscan.hmm(HMM0, training.trace[[s]]$hmm, domain.anno, nNPG, nOG,
                      geneseqs, seqs, statePred,postP,anno,
                      anno.txt, domain.txt, t.sm.evaluation, out.tag= paste(out.tag,'_viterbi', sep=''))
    
    if (predict.extra){
      statePred0 = HMM::viterbi(HMM0, seqs)
      postP0 = HMM::posterior(HMM0, seqs)
      output.NPscan.hmm(HMM0, training.trace[[s]]$hmm, domain.anno, nNPG, nOG,
                        geneseqs, seqs, statePred0,postP0,anno,
                        anno.txt, domain.txt, t.sm.evaluation, out.tag= paste(out.tag,'_viterbi_init', sep=''))
      
      i.be = c(seq(1,nH,3),seq(3,nH,3)); # begin and end states remove pseudo counts;
      finalHMM = training.trace[[s]]$hmm
      finalHMM$emissionProbs[i.be, ] = HMM0$emissionProbs[i.be,]
      finalHMM$emissionProbs[HMM0$emissionProbs==0] = 0;
      finalHMM$emissionProbs = 1/rowSums(finalHMM$emissionProbs) * finalHMM$emissionProbs
      finalHMM$transProbs[HMM0$transProbs==0] = 0
      finalHMM$transProbs = 1/rowSums(finalHMM$transProbs) * finalHMM$transProbs
      statePred.rp = HMM::viterbi(finalHMM, seqs)
      postP.rp = HMM::posterior(finalHMM, seqs)
      output.NPscan.hmm(HMM0, finalHMM, domain.anno, nNPG, nOG,
                        geneseqs, seqs, statePred.rp, postP.rp, anno,
                        anno.txt, domain.txt, t.sm.evaluation, out.tag= paste(out.tag,'_viterbi_rp', sep=''))
    }
    
    ################
    ### remove undesired transition and emmision due to pseudocounts
    ################
    HMM.nsc =  training.trace[[s]]$hmm
    HMM.nsc$transProbs = HMM.nsc$transProbs * (training.trace[[1]]$hmm$transProbs>0)
    HMM.nsc$transProbs = 1/rowSums(HMM.nsc$transProbs) * HMM.nsc$transProbs
    i.be = c(seq(1, nH, 3), seq(3, nH, 3)) # begin and end state, remove pseudo count information
    HMM.nsc$emissionProbs[i.be,] = HMM.nsc$emissionProbs[i.be,] * (training.trace[[1]]$hmm$emissionProbs[i.be,] > 0)
    HMM.nsc$emissionProbs[i.be,] = 1/rowSums(HMM.nsc$emissionProbs[i.be,]) * HMM.nsc$emissionProbs[i.be,]
  }else{
    HMM.nsc = HMM0
  }
  
  ################
  ### EM training, based on the viterbi training output, if desired
  ################
  if (EM.method == 'seqHMM'){
    seq.formated <- seqdef(t(seqs), 1:length(seqs),
                           labels = c('b','e',emit.symbols))
    HMM = build_hmm(state_names = States, 
                    observations = seq.formated,
                    initial_probs = HMM.nsc$startProbs,
                    transition_probs = HMM.nsc$transProbs,
                    emission_probs = HMM.nsc$emissionProbs)
    # HMM.reinit = build_hmm(state_names = States, # initialize using vertabi training output
    #                     observations = seq.formated,
    #                     initial_probs = training.trace[[s]]$hmm$startProbs,
    #                     transition_probs = training.trace[[s]]$hmm$transProbs,
    #                     emission_probs = training.trace[[s]]$hmm$emissionProbs)
    # alphabet(seqs)[1:10]
    # fit.HMM <- fit_model(HMM, threads=3, control_em = list(restart = list(times = 0)))
    fit.HMM <- fit_model(HMM,control_em = list(maxeval = 100, restart = list(times = 0)),
                         global_step=T, control_global = list(maxtime=1000),
                         local_step=T)
    # plot.hmm(fit.HMM$model)
    statePred.EM = hidden_paths(fit.HMM$model)
    statePred.EM = as.character(unlist(as.list(statePred.EM)))
    postP.EM = posterior_probs(fit.HMM$model)
    postP.EM = postP.EM[,,1]
    output.NPscan.hmm(seqHMM2HMM(HMM), seqHMM2HMM(fit.HMM$model), domain.anno, nNPG, nOG,
                      geneseqs, seqs, statePred.EM,postP.EM, anno, is.SM,
                      anno.txt, domain.txt, t.sm.evaluation, out.tag=paste(out.tag,'_viterbi_seqHMM_EM', sep=''))
  }else if (EM.method == 'depmixS4'){ # speed similar to viterbi according to testing: analysis.HMM.speedComparison
    if (0){
      depmix0 <- depmix(list(obs~1), data=data.frame(obs = seqs),nstates=nH,
                        family=list(multinomial('identity')))
      depmix0@prior@parameters$coefficients = HMM.nsc$startProbs
      depmix0@init = t(HMM.nsc$startProbs)
      depmix0@trDens[] = HMM.nsc$transProbs
      for (i in 1:length(depmix0@transition)){
        depmix0@transition[[i]]@parameters$coefficients = HMM.nsc$transProbs[i,]
      }
      for (i in 1:length(depmix0@response)){
        depmix0@response[[i]][[1]]@parameters$coefficients = HMM.nsc$emissionProbs[i,]
      }
    }
    dmHMM0.viterbi <- HMM2depmix(HMM.nsc, seqs)
    set.seed(3)
    ptm <- proc.time()
    dmHMM.viterbi.EM <- fit(dmHMM0.viterbi, emc = em.control(rand=F)) # no random start, otherwise, em.depmix gives an error message Starting values not feasible; please provide them"
    proc.time()-ptm
    ptm <- proc.time()
    dmHMM.viterbi.viterbi <- fit(dmHMM0.viterbi, emc = em.control(rand=F,classification='hard')) # no random start, otherwise, em.depmix gives an error message Starting values not feasible; please provide them"
    proc.time()-ptm
    # user    system   elapsed 
    # 33425.32  87016.23 128767.56 9 (1.5 days)
    
    if (0){
      dmHMM0 <- HMM2depmix(HMM0, seqs)
      xx <- depmix2HMM(dmHMM0)
      all(xx$Symbols == HMM0$Symbols)
      all(xx$emissionProbs == HMM0$emissionProbs)
      set.seed(3)
      ptm <- proc.time()
      dmHMM.EM <- fit(dmHMM0, emc = em.control(rand=F))
      proc.time()-ptm
    }
    if (0){
      #     iteration 0 logLik: -126822.8 
      #     iteration 5 logLik: -125822.1 
      #     iteration 10 logLik: -125609 
      #     iteration 15 logLik: -125575.9 
      #     iteration 20 logLik: -125543.3 
      #     iteration 25 logLik: -125514 
      #     iteration 30 logLik: -125504.1 
      #     iteration 35 logLik: -125495 
      #     iteration 40 logLik: -125475.7 
      #     iteration 45 logLik: -125464.7 
      #     iteration 50 logLik: -125459.2 
      #     iteration 55 logLik: -125456.1 
      #     iteration 60 logLik: -125455.2 
      #     iteration 65 logLik: -125453.6 
      #     iteration 70 logLik: -125453.3 
      #     iteration 75 logLik: -125453.1 
      #     iteration 80 logLik: -125453 
      #     iteration 85 logLik: -125452.9 
      #     iteration 90 logLik: -125452.8 
      #     iteration 95 logLik: -125452.7 
      #     iteration 100 logLik: -125452.6 
      #     iteration 105 logLik: -125452.6 
      #     iteration 110 logLik: -125452.6 
      #     iteration 115 logLik: -125452 
      #     iteration 120 logLik: -125450.8 
      #     iteration 125 logLik: -125450.7 
      #     iteration 130 logLik: -125450.5 
      #     iteration 135 logLik: -125448.6 
      #     iteration 140 logLik: -125448.6 
      #     iteration 145 logLik: -125448.6 
      #     iteration 150 logLik: -125448.6 
      #     iteration 155 logLik: -125448.6 
      #     converged at iteration 156 with logLik: -125448.5      
    }

    post = depmixS4::posterior(dmHMM.viterbi.EM)
    
    statePred.depmix = post[,1]
    statePred.depmix = HMM0$States[statePred.depmix] # statePred.depmix = HMM::viterbi(depmix2HMM(dmHMM.viterbi.EM), seqs)
    postP.depmix = forwardbackward(dmHMM.viterbi.EM)$gamma # postP.depmix = HMM::posterior(depmix2HMM(dmHMM.viterbi.EM), seqs)
    colnames(postP.depmix) = names(dmHMM.viterbi.EM@prior@parameters$coefficients)
    postP.depmix = t(postP.depmix);
    # postP.depmix = post[,2:ncol(post)] this is wrong
    
    output.NPscan.hmm(HMM.nsc, depmix2HMM(dmHMM.viterbi.EM), domain.anno,
                      geneseqs, seqs, statePred.depmix,postP.depmix,anno,
                      anno.txt, domain.txt, t.sm.evaluation, out.tag=paste(out.tag,'_viterbi_depmixEM', sep=''))

    class(dmHMM.viterbi.viterbi) = 'depmix.fitted'
    post = depmixS4::posterior(dmHMM.viterbi.viterbi)
    
    statePred.depmix = post[,1]
    statePred.depmix = HMM0$States[statePred.depmix] # statePred.depmix = HMM::viterbi(depmix2HMM(dmHMM.viterbi.EM), seqs)
    postP.depmix = forwardbackward(dmHMM.viterbi.viterbi)$gamma # postP.depmix = HMM::posterior(depmix2HMM(dmHMM.viterbi.EM), seqs)
    colnames(postP.depmix) = names(dmHMM.viterbi.viterbi@prior@parameters$coefficients)
    postP.depmix = t(postP.depmix);
    # postP.depmix = post[,2:ncol(post)] this is wrong
    
    output.NPscan.hmm(HMM.nsc, depmix2HMM(dmHMM.viterbi.viterbi), domain.anno, nNPG, nOG,
                      geneseqs, seqs, statePred.depmix,postP.depmix,anno,
                      anno.txt, domain.txt, t.sm.evaluation, out.tag=paste(out.tag,'_viterbi_depmixviterbi', sep=''))
    
    
  }else if (EM.method == 'HMM'){ # HMM
    step = 2; max.steps = 50; delta=1E-9
    training.traceBW = list()
    training.traceBW[['0']] = list(hmm = HMM.nsc, difference = Inf)
    s = 1
    while (1){
      cat((s-1)*step, training.traceBW[[paste(s-1)]]$difference)
      training.traceBW[[paste(s)]] = baumWelch(training.traceBW[[paste(s-1)]]$hmm, observation=seqs, maxIterations=step, delta=delta, pseudoCount=0.5)
      if (all(training.traceBW[[paste(s)]]$difference < delta))
        break
      cat('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b')
      s = s+1
    }
    statePred.EM = HMM::viterbi(training.traceBW[[s]]$hmm, seqs)
    postP.EM = HMM::posterior(training.traceBW[[s]]$hmm, seqs)
    output.NPscan.hmm(training.traceBW[[1]]$hmm, training.traceBW[[s]]$hmm, domain.anno, nNPG, nOG,
                      geneseqs, seqs, statePred.EM, postP.EM,anno,
                      anno.txt, domain.txt, t.sm.evaluation, out.tag= paste(out.tag,'_viterbi', sep=''))
    
  }else if (EM.method == 'null'){
    # nothing is done here
  }
  # HMM.ori
  # HMM.dom
  # HMM.dom.ori
  
  if (0){
    image(vt$hmm$transProbs)
    image(HMM$transProbs)
    require(gplots)
    image(log(vt$hmm$emissionProbs),col=greenred(256))
    image(log(HMM$emissionProbs),col=greenred(256))
    
    # label KUs
    # semisupervised EM learning
    # Predictions
    
    #### output
    cat('\n#Identified clusters: ', nrow(gene.ranges))
    to.keep.extend = extend.index(core.regions, window.extend, sides='both', do.unique=T)
    to.keep.extend = to.keep.extend[to.keep.extend<=length(anno) & to.keep.extend>=1]
    anno$PARENT[1] ==c()
    is.enzyme.all[] = c('', 'Yes')[is.enzyme.all+1]
    out = cbind(chr = chr, gene=gene.ID, 'protein ID' = prot.ID, Existing.Anno = anno@elementMetadata[,toupper(desc.fname)],
                is.enzyme.all, domains = ipr.anno)[to.keep.extend,]
    
    rownames(geneID2clusterID) = geneID2clusterID[,1];
    out = cbind(out, clusterID = mat.fill.row(geneID2clusterID, rownames(out), '')[,2])
    
    write.xlsx(out, out.file)
    
    # HTML(out, 'out.html')
    return(out)
  }
}

viterbiTraining.depmix <- function(dmHMM0, seqs){
  # implement fast viterbiTraining based on depmixS4:::viterbi
  class(dmHMM0) = 'depmix.fitted'
  s = depmixS4:::viterbi(seqs, dmHMM0)
  s = depmixS4::posterior(dmHMM0)
  logLik(dmHMM0)
}

mut.gene <- function(seq, mut.model){
  # implement fast viterbiTraining based on depmixS4:::viterbi
  
}

HMM2depmix <- function(HMM, seqs=NULL){
  # 20160331, initialize a HMM model based on an HMM model in package HMM
  nH = length(HMM$States)
  if (is.null(seqs)){
    seqs = HMM$Symbols
  }
  depmix0 <- depmix(list(obs~1), data=data.frame(obs = seqs),nstates=nH,
                    respstart = setNames(t(HMM$emissionProbs)[1:length(HMM$emissionProbs)], 
                                         rep(colnames(HMM$emissionProbs), 
                                             times = nrow(HMM$emissionProbs))),
                    trstart = t(HMM$transProbs), 
                    instart = t(HMM$startProbs),
                    family=list(multinomial('identity')))
  
  # names(depmix0@prior@parameters$coefficients) = names(HMM$startProbs)
  return(depmix0)
}

depmix2HMM <- function(depmixHMM){
  # 20160331
  emissionProbs = c()
  for (i in 1:length(depmixHMM@response)){
    emissionProbs = rbind(emissionProbs, depmixHMM@response[[i]][[1]]@parameters$coefficients)
  }

  hmm.model =   initHMM(States = names(depmixHMM@prior@parameters$coefficients),
                        Symbols = colnames(emissionProbs), 
                        startProbs = depmixHMM@init,
                        transProbs = t(depmixHMM@trDens[1,,]),
                        emissionProbs = emissionProbs)
  return(hmm.model)
}

seqHMM2HMM <- function(seqhmm.model = fit.HMM$model){
  # 20160330, YFL
  hmm.model =   initHMM(States = seqhmm.model$state_names, 
                        Symbols = seqhmm.model$symbol_names, 
                        startProbs = seqhmm.model$initial_probs,
                        transProbs = seqhmm.model$transition_probs,
                        emissionProbs = seqhmm.model$emission_probs)
  return(hmm.model)
}

output.NPscan.hmm <- function(initialHMM, finalHMM, domain.anno,nNPG, nOG,
                              geneseqs, seqs, statePred, postP,anno,
                              anno.txt, domain.txt, truth = t.sm.evaluation, out.tag){
  # 20160508, v2, output cluster types, 
  require('lattice')
  pdf(paste('HMM_learning',out.tag, '.pdf', sep=''), 7,3)
  plot.hmm(initialHMM)
  plot.hmm(finalHMM)
  dev.off()
  append=F
  out.file = paste('NPScan_Emission', out.tag, '.xlsx', sep='')
  for (i in 1:nNPG){
    tag = paste('NPG', i, '.d', sep='')
    top.NPG.domains = sort(initialHMM$emissionProbs[tag,], decreasing = T)
    top.NPG.domains =  cbind(domain.anno[names(top.NPG.domains)], top.NPG.domains)
    write.xlsx2(top.NPG.domains, append=append, sheetName = paste(tag, '_init', sep=''), file = out.file)
    tag = paste('NPG', i, '.d', sep='')
    append = T
    top.NPG.domains = sort(finalHMM$emissionProbs[tag,], decreasing = T)
    top.NPG.domains =  cbind(domain.anno[names(top.NPG.domains)], top.NPG.domains)
    write.xlsx2(top.NPG.domains, append=append, sheetName = paste(tag, '_Final', sep=''), file = out.file)
  }
  for (i in 1:nOG){
    tag = paste('OG', i, '.d', sep='')
    top.NPG.domains = sort(initialHMM$emissionProbs[tag,], decreasing = T)
    top.NPG.domains =  cbind(domain.anno[names(top.NPG.domains)], top.NPG.domains)
    write.xlsx2(top.NPG.domains, append=append, sheetName = paste(tag, '_init', sep=''), file = out.file)
    tag = paste('OG', i, '.d', sep='')
    top.NPG.domains = sort(finalHMM$emissionProbs[tag,], decreasing = T)
    top.NPG.domains =  cbind(domain.anno[names(top.NPG.domains)], top.NPG.domains)
    write.xlsx2(top.NPG.domains, append=append, sheetName = paste(tag, '_Final', sep=''), file = out.file)
  }
  # top.diff.domains = sort(initialHMM$emissionProbs['NPG.d',]/initialHMM$emissionProbs['OG.d',], decreasing = T)[1:100]
  # top.diff.domains =  cbind(domain.anno[names(top.diff.domains)], top.diff.domains)
  # write.csv(top.NPG.domains, 'top.NPG1.domains_init.csv')
  # write.csv(top.OG1.domains, 'top.OG1.domains_init.csv')
  # write.csv(top.OG2.domains, 'top.OG2.domains_init.csv')
  # write.csv(top.diff.domains, 'top.diff.domains_init.csv')
  
  # top.NPG.domains = sort(finalHMM$emissionProbs[tag,], decreasing = T)   
  nH = nrow(postP)/3
  
  ##### 20160508, computing the KU domains
  Type = vector(mode = 'character', length(statePred));
  ipr2sm = IPR2SMtype();
  i.domains = !(seqs%in%c('b','e'))
  SMTypes = colnames(ipr2sm)
  i.KU.domains = !is.na(match(seqs, rownames(ipr2sm)));
  SMTypeMat = as.matrix(ipr2sm[seqs[i.KU.domains],])
  Type[i.KU.domains] = apply(SMTypeMat, MARGIN = 1, function(x)paste(SMTypes[which(x==1)], ':',1, collapse = '; ', sep=''))
  geneID = cumsum(seqs == 'b');
  Type.gene = vector(mode = 'character', max(geneID))
  a = tapply(1:sum(i.KU.domains), geneID[i.KU.domains], function(x){
    x= colSums(SMTypeMat[x,,drop=F]);
    # cat(x)
    paste(SMTypes[which(x>0)], ':',x[x>0], collapse = '; ', sep='')
  })
  Type.gene[as.numeric(names(a))]  = a;
  # assign cluster ID
  NP.Gene.pred = regexpr('^OG.*b', statePred[seqs=='b'])<=0
  cluster.start =  diff(c(F, NP.Gene.pred)) == 1;
  chr = as.character(anno@seqnames)
  chr.start = (chr != c('begin',chr[1:(length(chr)-1)]) )
  cluster.start[chr.start & NP.Gene.pred] = T
  cluster.ID =  cumsum(cluster.start);
  cluster.ID[!NP.Gene.pred] = 0
  cluster.size = cluster.ID;
  cluster.size[cluster.ID!=0] = (unique.count(cluster.ID[cluster.ID!=0])$counts.unique)[as.character(cluster.ID[cluster.ID!=0])]
  cluster.max.p = cluster.ID;
  ii = which(cluster.ID!=0)
  cluster.max.p[cluster.ID!=0] = unlist(tapply(t(postP)[seqs=='b',][cbind(ii, match((statePred[seqs=='b'])[ii],rownames(postP)))], INDEX = cluster.ID[cluster.ID!=0], max))[as.character(cluster.ID[cluster.ID!=0])]
  
  ## aggregate gene NP type to the cluster level
  Type.cluster = vector(mode = 'character', max(geneID));
  gene.SMTypeMat = do.call(rbind, tapply(1:sum(i.KU.domains), geneID[i.KU.domains], function(x){
    x= colSums(SMTypeMat[x,,drop=F])}, simplify = T));  
  a = tapply(1:nrow(gene.SMTypeMat), cluster.ID[as.numeric(rownames(gene.SMTypeMat))], function(x){
    x= colSums(gene.SMTypeMat[x,,drop=F]);
    paste(SMTypes[which(x>0)], ':',x[x>0], collapse = '; ', sep='')
  })
  a= unlist(a); a = a[setdiff(names(a), '0')]
  i.matched = as.character(cluster.ID) %in% names(a)
  Type.cluster[i.matched] = a[as.character(cluster.ID[i.matched])]
  Type.cluster[cluster.ID!=0 & !i.matched] = 'UU'
  Type.cluster[cluster.ID!=0 & !i.matched & regexpr('^NPG1.*b', statePred[seqs=='b'])<=0] = 'UU(NPG2)'
  
  pred.domains =  data.frame(Gene = geneseqs,
                        Feature = seqs,
                        Annotation = domain.anno[seqs],
                        State = statePred,
                        NP.Type = Type,
                        Posterior = t(postP))
  if (0){ # domain level predictions
    out.file = paste('NPScan_DomainPred_', out.tag, '.csv', sep='')
    write.csv(pred.domains, row.names = F,  file = out.file)
  }
  gene.out.file = paste('NPScan_GenePred2_', out.tag, '.csv', sep='')
  pred.gene = data.frame(chr = as.character(anno@seqnames), 
                         pred.domains[pred.domains[,'Feature']=='b',c(1,4, seq(6, ncol(pred.domains),3))], 
                         Gene.NP.type = Type.gene,
                         Cluster.NP.type = Type.cluster,
                         Cluster.ID = cluster.ID, 
                         Cluster.size = cluster.size,
                         Cluster.p.max = cluster.max.p,
                         known.KU = truth*1,
                         Gene.Anno = anno.txt, 
                         Domains = domain.txt)
  pdf(paste('perf_', out.tag, '.pdf', sep=''),3,3) # 20160404
  for (i in 3+(1:nH)){
    s = sub('Posterior\\.(.+)\\.b', '\\1', colnames(pred.gene)[i])
    print(hist.by(log10(pred.gene[,i]/(1-pred.gene[,i])), c('Other genes', 'True NP genes')[1+truth], xlab = paste('log odds ', s, sep=''), by.name = ''))
  }
  for (i in 3+(1:nH)){
    print(hist.by(pred.gene[,i], c('Other genes', 'True NP genes')[1+truth], xlab = paste('probability ', s, sep=''), by.name = ''))
  }
  dev.off()
  
  write.csv(pred.gene, row.names = F, file = gene.out.file)

  ### visualization of gene probability
  pdf(paste('prob_', out.tag, 'plot.pdf', sep=''), 10,6)
  dat = data.frame()
  for (i in 3+(1:(nH-1))){
    dat = rbind(dat, data.frame(gene = 1:nrow(pred.gene), score = pred.gene[[i]], 
                                type = sub('Posterior\\.(.*)\\.b', '\\1', colnames(pred.gene)[i]),
                                chr = pred.gene$chr,
                                is.KU = pred.gene$known.KU))
  }
  print(ggplot(data = dat) +  geom_line(mapping = aes(x=gene, y=score, color=type), alpha=0.3) + 
          geom_point(aes(x=gene, y=score),shape = 1,data = dat[dat$type=='NPG1' & dat$is.KU,]) + 
          facet_wrap(~chr, nrow=4,  scales="free") + labs(color = 'scores')+ 
          theme_bw() +  theme(panel.grid.major = element_blank()))
  print(ggplot(data = dat) +  geom_line(mapping = aes(x=gene, y=log10(score/(1-score)), color=type), alpha=0.3) + 
          geom_point(aes(x=gene, y=log10(score/(1-score))),shape = 1,data = dat[dat$type=='NPG1' & dat$is.KU,]) + 
          facet_wrap(~chr, nrow=4,  scales="free") + labs(color = 'scores') +
          theme_bw() +  theme(panel.grid.major = element_blank()))
  dev.off()
  # save(list = c('initialHMM', 'finalHMM', 'seqs', 'domain.anno', 'nNPG','nOG','geneseqs',
  #               'seqs', 'statePred', 'postP', 'anno','anno.txt', 'domain.txt', 'truth', 'out.tag'), file = paste(out.tag, '.RData', sep=''))
  # save(list = c('initialHMM', 'finalHMM', 'seqs', 'statePred', 'postP'), file = paste(out.tag, '.RData', sep=''))
  
}

IPR2SMtype <- function(file = '/Users/yongli/Dropbox/NPGC/NPGCquery_data/SM.domains_manualAnno_v2.xlsx'){
  # 20160508
  SM.doms = read.xlsx2(file,sheetIndex = 1, as.is =T)
  SM.doms$BGC.confidence = as.numeric(as.character(SM.doms$BGC.confidence))
  SMtypes = unique(unlist(strsplit(as.character(SM.doms$NP_Class), split = '; ')))
  SMtypes = c(SMtypes, 'UU')
  #   IPR2SM = matrix(0, nrow = sum(SM.doms$NP_Class!=''), ncol = length(SMtypes), 
  #                   dimnames = c(SM.doms$id[SM.doms$NP_Class!=''],SMtypes));
  seq = strsplit(as.character(SM.doms$NP_Class[SM.doms$NP_Class!='']),  split = '; ');
  names(seq) = SM.doms$id[SM.doms$NP_Class!='']
  IPR2SM = seq2mat(seq, alphabet = SMtypes)
  return(IPR2SM)
}

NPscan.postprocess <- function(files = dir(pattern = '^NPScan_GenePred2_*'), tag = '', remove.TF=F, remove.transporter=F,
                               meta.file = '/Users/yongli/Dropbox/NPGC/NPGCquery_data/NPGCquery_meta.txt',
                               cluster.info.file = paste(tag, 'cluster_info.xlsx', sep=''),
                               length.cutoff = 5, length.cutoff.max = 25, p.cutoff = 0.99, extra.gene = 1, Walsh.only = F, verbose = F){
  # 20160508
  # 20160517: add semi-UU (without condensation and Keto-synthase), 
  #           and highlight special protein types
  #           "radical.SAM"/'(FAD|Flavin)' & "oxygenase"/ IPR005123 - Fe(II) oxygenase, IPR014030/IPR014031 - KS domain, IPR001242 -- condensation
  # 20160526: add length.cutoff.max and Walsh.only
  # note that the cutoffs only applies to UUs and semiUUs
  cluster.info = c()
  for (f in files){
    dat = read.csv(f);
    if (0){
      idx = dat$Cluster.NP.type == 'UU' & dat$Cluster.size >= length.cutoff & dat$Cluster.size <= length.cutoff.max & dat$Cluster.p.max >=p.cutoff;
      idx = dat$Cluster.ID %in% unique(dat$Cluster.ID[idx])
      UU = dat[extend.index(which(idx), n = extra.gene),];
      idx = !(dat$Cluster.NP.type %in% c('UU', 'UU(NPG2)')) & dat$Cluster.NP.type!=''
      idx = dat$Cluster.ID %in% unique(dat$Cluster.ID[idx])
      KU = dat[extend.index(which(idx), n = extra.gene),]
    }
    idx = dat$Cluster.NP.type == 'UU(NPG2)' & dat$Cluster.size >= length.cutoff & dat$Cluster.size <= length.cutoff.max & dat$Cluster.p.max >=p.cutoff
    idx = dat$Cluster.ID %in% unique(dat$Cluster.ID[idx])
    UU.NNPG2 = dat[extend.index(which(idx), n = extra.gene), ];
    
    file = '/Users/yongli/Dropbox/NPGC/NPGCquery_data/SM.domains_manualAnno_v2.xlsx'
    SM.doms = read.xlsx2(file,sheetIndex = 1, as.is =T)
    core.enzyme.IPR = paste(SM.doms$id[SM.doms$SufficientFor!=''], collapse = '|')
    is.semiUU =  as.character(dat$Gene.NP.type)!='' & regexpr(core.enzyme.IPR,as.character(dat$Domains), perl = T)<0
    is.KU =  as.character(dat$Gene.NP.type)!='' & regexpr(core.enzyme.IPR,as.character(dat$Domains), perl = T)>0
    is.semiUU = by(is.semiUU, INDICES = dat$Cluster.ID, FUN = sum) > 0
    is.KU = by(is.KU, INDICES = dat$Cluster.ID, FUN = sum) >0
    is.semiUU = is.semiUU & ! is.KU
    is.UU = by(dat$Cluster.NP.type == 'UU', INDICES = dat$Cluster.ID, FUN = sum) > 0
    
    
    if (verbose){
      cat('\nCore enzymes motifs but not core enzyme domains:\n')
      print(dat[as.character(dat$Gene.NP.type)=='' & regexpr(core.enzyme.IPR,as.character(dat$Domains), perl = T)>0,])
      
      cat('\nCore enzymes (motifs) outside clusters:\n')
      print(dat[regexpr(core.enzyme.IPR,as.character(dat$Domains), perl = T)>0 & dat$Cluster.ID == 0,])
    }
    
    if (Walsh.only){ # only select clusters that contain understudied oxidoreductase
      is.Walsh = regexpr('(radical.SAM|IPR005123)',as.character(dat$Domains), perl = T)>0 | (regexpr('(FAD|Flavin)',as.character(dat$Domains), perl = T)>0 & regexpr('oxygenase',as.character(dat$Domains), perl = T)>0)
      is.Walsh = by(is.Walsh, INDICES = dat$Cluster.ID, FUN = sum) >0
      is.KU = is.KU & is.Walsh;
      is.semiUU = is.semiUU & is.Walsh;
      is.UU = is.UU & is.Walsh;
    }
    
    is.enzyme = vector(mode = 'logical', length = nrow(dat)) | T
    if (remove.TF)
      is.enzyme = is.enzyme & !regexpr('Transcription factor',as.character(dat$Domains), perl = T, ignore.case = T)>0
    if (remove.transporter)
      is.enzyme = is.enzyme & !regexpr('Major facilitator superfamily|transporter',as.character(dat$Domains), perl = T, ignore.case = T)>0
    cluster.size = by(is.enzyme, dat$Cluster.ID, sum);
    to.keep = cluster.size >= length.cutoff & cluster.size <= length.cutoff.max & by(dat$Cluster.p.max >= p.cutoff, INDICES = dat$Cluster.ID, FUN = sum) > 0
    UU = dat[extend.index(which(dat$Cluster.ID %in% setdiff(as.numeric(names(is.KU)[is.UU & to.keep]),0)), n = extra.gene),];
    semiUU = dat[extend.index(which(dat$Cluster.ID %in% setdiff(as.numeric(names(is.KU)[is.semiUU & to.keep]),0)), n = extra.gene),];
    KU = dat[extend.index(which(dat$Cluster.ID %in% setdiff(as.numeric(names(is.KU)[is.KU]),0)), n = extra.gene),];
    
    cat('KU cluster:', sum(is.KU & to.keep), '\nsemi-UU cluster:', sum(is.semiUU & to.keep), '\nUU:', sum(is.UU & to.keep), '\n')
    # is.KU1 =  regexpr('(PK|NRP)',as.character(dat$Gene.NP.type), perl = T) < 0 & regexpr('(IPR014030|IPR014031|IPR001242)',as.character(dat$Domains), perl = T)>0
    #     write.csv(UU, row.names = F, file = paste('UUselect_', f, sep=''))
    #     write.csv(UU.NNPG2, row.names = F, file = paste('UUselect.NPG2_', f, sep=''))
    #     write.csv(KU, row.names = F, file = paste('KUselect_', f, sep=''))
    if (nrow(UU)){
      write.xlsx2(UU, sheetName = 'UU', row.names = F, file = paste('UUselect_', tag, sub('it.*.csv', '.xlsx', f), sep=''))
      xlsx.color.NPscan(paste('UUselect_', tag, sub('it.*.csv', '.xlsx', f), sep=''))
      unlink(paste('UUselect_', tag, sub('it.*.csv', '.xlsx', f), sep=''))
    }
    if (nrow(UU.NNPG2)){
      write.xlsx2(UU.NNPG2, sheetName = 'UU.NPG2', row.names = F, file = paste('GlycoUUSelect_', tag, sub('it.*.csv', '.xlsx', f), sep=''))
      xlsx.color.NPscan(paste('GlycoUUSelect_',tag,  sub('it.*.csv', '.xlsx', f), sep=''))
      unlink(paste('GlycoUUSelect_', tag, sub('it.*.csv', '.xlsx', f), sep=''))
    }
    if (nrow(semiUU)){
      write.xlsx2(semiUU, sheetName = 'semiUU', row.names = F, file = paste('semiUUselect_', tag, sub('it.*.csv', '.xlsx', f), sep=''))
      xlsx.color.NPscan(paste('semiUUselect_', tag, sub('it.*.csv', '.xlsx', f), sep=''))
      unlink(paste('semiUUselect_', tag, sub('it.*.csv', '.xlsx', f), sep=''))
    }
    if (nrow(KU)){
      write.xlsx2(KU, sheetName = 'KU', row.names = F, file = paste('KUselect_', tag, sub('it.*.csv', '.xlsx', f), sep=''))
      xlsx.color.NPscan(paste('KUselect_', tag, sub('it.*.csv', '.xlsx', f), sep=''))
      unlink(paste('KUselect_', tag, sub('it.*.csv', '.xlsx', f), sep=''))
    }
    
    meta = read.csv(meta.file, sep = '\t', as.is = T)
    j = which(sub(' ', '_', meta$genome.ID) == sub('NPScan_GenePred2_(.+)it.+ie.+$', '\\1', f) | 
              sub(' ', '_', meta$species) == sub('NPScan_GenePred2_(.+)it.+ie.+$', '\\1', f))
    
    in.house.genome = '';
    jgi.genome = ''
    if (meta$Source[j] == 'JGI'){
      jgi.genome =  meta$genome.ID[j];
    }
    if (meta$Source[j] != 'JGI'){
      idx = c("gff.file", "iprscan.tab.file", "pep.fasta.file", "DNA.file", "gene.definition", "proteinID") # c(3,4,6,8,9,10)
      for (i in c("gff.file", "iprscan.tab.file", "pep.fasta.file", "DNA.file")){
        meta[meta$folder!='',i] = paste(meta$folder[meta$folder!=''], meta[meta$folder!='',i], sep='/')
      }
      in.house.genome = paste(idx, ' = \'', meta[j,idx], '\'', sep='', collapse = '; ')
    }
    if (nrow(UU)){
      cluster.locs = data.frame(first = which.first.by(UU$Cluster.ID), last = which.last.by(UU$Cluster.ID))
      cluster.locs = cluster.locs[setdiff(rownames(cluster.locs), '0'),]
      cluster.info = rbind(cluster.info,
                           data.frame(ClusterID = paste('UU',UU$Cluster.ID[cluster.locs$first], sep=''),  # , meta$genome.ID[j], '|'
                                      type = 'UU',
                                      'GenBank Genome' = '', 
                                      'JGI Genome' = jgi.genome,
                                      'Same Genome' = '',
                                      'In House Genome'= in.house.genome,
                                      species = sub('(\\S+ \\S+).*$','\\1', meta$species[j]),
                                      'First Protein' = UU$Gene[cluster.locs$first], 
                                      'Last Protein' = UU$Gene[cluster.locs$last]))      
    }
    if (nrow(semiUU)){
      cluster.locs = data.frame(first = which.first.by(semiUU$Cluster.ID), last = which.last.by(semiUU$Cluster.ID))
      cluster.locs = cluster.locs[setdiff(rownames(cluster.locs), '0'),]
      cluster.info = rbind(cluster.info,
                           data.frame(ClusterID = paste('semiUU', semiUU$Cluster.ID[cluster.locs$first], sep=''),  # meta$genome.ID[j], '|', 
                                      type = 'semiUU',
                                      'GenBank Genome' = '', 
                                      'JGI Genome' = jgi.genome,
                                      'Same Genome' = '',
                                      'In House Genome'= in.house.genome,
                                      species = sub('(\\S+ \\S+).*$','\\1', meta$species[j]),
                                      'First Protein' = semiUU$Gene[cluster.locs$first], 
                                      'Last Protein' = semiUU$Gene[cluster.locs$last]))      
      
    }

    ### 20160527: output cluster info
    # ClusterID	GenBank Genome	JGI Genome	Same Genome	In House Genome	species	First Protein	Last Protein
    # Ca157				iprscan.tab.file = 'CA_K87_contig_Anidulans.faa.tsv'; gff.file = 'CA_K87_contig_Anidulans.gff'; DNA.file = 'CA_K87_contig.fasta'; pep.fasta.file = 'CA_K87_contig_Anidulans.faa'; gene.definition = 'transcript'; proteinID = 'ID'	Calcarisporium arbuscula	g7062.t1	g7069.t1
    # Afu1g17740		Aspfu1			Aspergillus fumigatus	Afu1g17700	Afu1g17750
  }
  cluster.info$ClusterID = paste(cluster.info$JGI.Genome, cluster.info$GenBank.Genome, '.', cluster.info$ClusterID, sep='')
  write.xlsx2(cluster.info, row.names = F, sheetName = 'cluster.info', file = cluster.info.file)
  write.xlsx2(cluster.info[cluster.info$type == 'semiUU',], row.names = F, sheetName = 'cluster.info', file = sub('.xls', '_semiUU.xls', cluster.info.file))
  write.xlsx2(cluster.info[cluster.info$type == 'UU',], row.names = F, sheetName = 'cluster.info', file = sub('.xls', '_UU.xls', cluster.info.file))
}

xlsx.color.NPscan <- function(xlsx.file = 'nidulans.deepAnno.all.xlsx', out.file=paste('colored_', xlsx.file, sep='')){
  # Yong Fuga Li, 20141004
  xlsx.color(xlsx.file = xlsx.file, FUN.select = FUN.select.semiUU.NPScan, fill.color = 'purple', out.file = out.file, na.strings='|')
  xlsx.color(xlsx.file = out.file, FUN.select = FUN.select.Walsh.NPScan, border.color = 'blue', out.file = out.file, na.strings='|')
}

FUN.select.semiUU.NPScan = function(x){
  # semi-UU (without condensation and Keto-synthase), 
  y = matrix(F, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x));
  y[,'Domains'] <- y[,'Gene.NP.type'] <- regexpr('(PK|NRP)',as.character(x$Gene.NP.type), perl = T)>=0 & regexpr('(IPR014030|IPR014031|IPR001242)',as.character(x$Domains), perl = T)<0
  return(y)
}

FUN.select.Walsh.NPScan <- function(x){
  #   #           and highlight special protein types
  #           "radical.SAM"/'(FAD|Flavin)' & "oxygenase"/ IPR005123 - Fe(II) oxygenase, IPR014030/IPR014031 - KS domain, IPR001242 -- condensation
  y = matrix(F, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x));
  y[,'Domains'] <- y[,'Gene.NP.type'] <- regexpr('(radical.SAM|IPR005123)',as.character(x$Domains), perl = T)>0 | (regexpr('(FAD|Flavin)',as.character(x$Domains), perl = T)>0 & regexpr('oxygenase',as.character(x$Domains), perl = T)>0)
  return(y)
}

xlsx.color.mergedDeepAnno <- function(xlsx.file = 'nidulans.deepAnno.all.xlsx', out.file=paste('colored_', xlsx.file, sep='')){
  # Yong Fuga Li, 20141004
  xlsx.color(xlsx.file = xlsx.file, include.header=T, FUN.select = function(x){y = matrix(T, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x)); y},
             border.color = 'grey', 
             out.file = out.file, na.strings='|')  # change global style
  xlsx.color(xlsx.file = out.file, include.header=T, FUN.select = function(x){y = matrix(T, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x)); y},
             font=list(color = NULL, heightInPoints=12, name='Arial', isItalic=F, isBold=F, isStrikeout=F, underline=NULL), 
             out.file = out.file, na.strings='|')  # change global style
  xlsx.color(xlsx.file = out.file, header = F,include.header=F, FUN.select = function(x){y = matrix(F, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x)); y[1,]=T; y},
             font=list(color = NULL, heightInPoints=12, name='Arial', isItalic=F, isBold=T, isStrikeout=F, underline=NULL), 
             out.file = out.file, na.strings='|')  # bold headers
  xlsx.color(xlsx.file = out.file, FUN.select = FUN.select.KU.mergedDeepAnno, fill.color = 'purple', out.file = out.file, na.strings='|')
  xlsx.color(xlsx.file = out.file, FUN.select = FUN.select.Walsh.mergedDeepAnno, border.color = 'blue', out.file = out.file, na.strings='|')
}

FUN.select.KU.mergedDeepAnno = function(x){
  # semi-UU (without condensation and Keto-synthase), 
  y = matrix(F, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x));
  y[,'IPR_Domain_Annotation'] <- is.KU.ipr(x[,'IPR_Domain_Annotation'])
  y[,'Manual_Annotation'] <- is.KU(x[,'Manual_Annotation'])
  return(y)
}

FUN.select.Walsh.mergedDeepAnno <- function(x){
  #   #           and highlight special protein types
  #           "radical.SAM"/'(FAD|Flavin)' & "oxygenase"/ IPR005123 - Fe(II) oxygenase, IPR014030/IPR014031 - KS domain, IPR001242 -- condensation
  y = matrix(F, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x));
  y[,'IPR_Domain_Annotation'] <- regexpr('(radical.SAM|IPR005123)',as.character(x$IPR_Domain_Annotation), perl = T)>0 | (regexpr('(FAD|Flavin)',as.character(x$IPR_Domain_Annotation), perl = T)>0 & regexpr('oxygenase',as.character(x$IPR_Domain_Annotation), perl = T)>0)
  return(y)
}

is.KU.ipr <- function(txt){
  ipr2sm = IPR2SMtype();
  iprIDs = paste(names(which(apply(ipr2sm[,1:10], MARGIN = 1, sum)>0)), collapse ='|')
  return(regexpr(iprIDs, txt, perl=T)>0)
}

merge.deepAnno <- function(clusterIDs = c('semiUU174', 'semiUU204', 'semiUU559',
                                           'UU1', 'UU10', 'UU29', 'UU48'),
                            out.file = 'merged.xlsx',
                            root = '/Users/yongli/Dropbox/NPGC/NPGCquery_data'){
  # Yong Fuga Li, 20160604
  require(xlsx)
  require('XLConnect')
  setwd(root)
  append=F
  for (i in clusterIDs){
    f = paste('colored_', i, '.xlsx', sep='');
    dat = read.xlsx2(f, sheetIndex = 1)
    i.range = which(dat$cluster.boundary == 'Boundary')
    dat.out = dat[i.range[1]:i.range[2], c('name', 'length',	'Existing.Anno', 'domains')];
    dat.out$Existing.Anno = sub('Uncharacterized ORF; ', '', dat.out$Existing.Anno)
    dat.out$Existing.Anno = sub('^Ortholog of .*$', '', dat.out$Existing.Anno)
    dat.out$Existing.Anno = sub('^.*description:\\"([^\\"]+)\\".*$', '\\1', dat.out$Existing.Anno)
    dat.out$name = sub('transcript:','', dat.out$name)
    colnames(dat.out) = c('Gene', 'Length', 'Manual_Annotation', 'IPR_Domain_Annotation')
    write.xlsx2(dat.out, file = out.file, sheetName = i, append = append, row.names = F)
    append = T
  }
  
  wb <- loadWorkbook(out.file, create = TRUE)
  cs <- createCellStyle(wb)
  # Specify to wrap the text
  setWrapText(cs, wrap = TRUE)
  for (cID in clusterIDs){
    setColumnWidth(wb,sheet=cID,column=1,width=256*9)
    setColumnWidth(wb,sheet=cID,column=2,width=256*8)
    setColumnWidth(wb,sheet=cID,column=3,width=256*8*4)
    setColumnWidth(wb,sheet=cID,column=4,width=256*8*16)
    for (r in 1:getLastRow(wb,cID)){
      setCellStyle(wb, sheet = cID, row = r, col = 3, 
                   cellstyle = cs)
      setCellStyle(wb, sheet = cID, row = r, col = 4,
                   cellstyle = cs)
    }
    saveWorkbook(wb)
  }
  
  xlsx.color.mergedDeepAnno(out.file)
  # in excels, select all tabs together and then print "entire workbook" as pdf in landscape mode with 65% size.
}

get.NPScan.nchr <- function(NPScan.files = dir(pattern = 'NPScan_GenePred2_.*_viterbi.csv')){
  # get the number of chromosomes in a genome, 20160802
  nchr = vector('numeric', length(NPScan.files));
  names(nchr) = NPScan.files
  for (f in NPScan.files){
    dat = read.csv(f)
    nchr[f] = length(unique(dat$chr))
  }
  return(nchr)
}
