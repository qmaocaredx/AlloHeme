ow = options('warn')
options(warn = 1)


read.samplesheet <- function(manifest, is.csv = F, 
                             requiredFields.SampleSpecific = c('SampleType', 'InputDNA', 'DNAType', 'GeneticRelationship', 'Baselines', 'ExclusionFile')){
  #### read illumina sample sheet csv files
  # treat , as space to handle extra , in sample sheet csv files
  #############################
  
  dat = readLines(con = manifest, warn = F)
  
  ### remove empty lines
  dat = dat[regexpr("^[\\s,]*$", dat, perl = T) <= 0]
  
  ## remove comments
  dat = sub('#.*','',dat)
  
  ### get section information
  section.boundary = regexpr("^[\\s,]*\\[.+\\][\\s,]*$", dat, perl = T) > 0
  section.names = sub("^[\\s,]*\\[(.+)\\][\\s,]*$", '\\1',dat[section.boundary]);
  
  sections = data.frame(section.names = section.names, 
                        section.start = which(section.boundary)+1,
                        section.end = c(which(section.boundary)[-1]-1, length(dat)))
  
  # parse each sections
  out = list()
  for (i in 1:nrow(sections)){
    idx = seq2(sections$section.start[i], sections$section.end[i], 1);
    if(length(idx)){
      if (sections$section.names[i] =='Data'){
        out[[i]] = 
          read.csv(textConnection(dat[idx]), header = T, strip.white = T, blank.lines.skip = TRUE, stringsAsFactors = F)
        if (!all(requiredFields.SampleSpecific %in% names(out[[i]]))){
          stop(paste('Require the following columns in the data section: ', paste(requiredFields.SampleSpecific, collapse = ', '), sep =''))
        }  
        if ("ExclusionFile" %in% names(out[[i]])){
          out[[i]]$ExclusionFile[is.na(out[[i]]$ExclusionFile)] = '';
        }
        out[[i]] = out[[i]][,colSums(!is.na(out[[i]]))>0, drop=F]
      }else if (sections$section.names[i] == 'Reads'){
        out[[i]] = 
          read.csv(textConnection(dat[idx]), header = F, strip.white = T, blank.lines.skip = TRUE, as.is = T, stringsAsFactors = F)
        colnames(out[[i]])[1] = 'Length'
        out[[i]] = out[[i]][,colSums(!is.na(out[[i]]))>0, drop=F]
      }else{
        out[[i]] = 
          read.csv(textConnection(dat[idx]), header = F, row.names = 1, strip.white = T, blank.lines.skip = TRUE, as.is = T, stringsAsFactors = F)
        out[[i]] = out[[i]][,colSums(!is.na(out[[i]]))>0, drop=F]
        x = (out[[i]][,1])
        names(x) = rownames(out[[i]])
        out[[i]] = as.list(x)
      }
    }else{
      out[[i]] = NA
    }
  }
  names(out) = sections$section.names
  print(out)
  return(out)
}

plot.AF <- function(dat, tag = ''){
  eps = .Machine$double.eps
  n.tot = dat$Ref_Allele_Counts + dat$Alt_Allele_Counts + 2*eps
  n.minor = apply(cbind(dat$Ref_Allele_Counts, dat$Alt_Allele_Counts), MARGIN = 1, min)
  r = (n.minor+eps)/n.tot;
  xmax = min(max(r[r<0.35], na.rm = T)*1.1,0.5)
  col = c('black', 'red')[(dat$Ref_Allele_Counts>dat$Alt_Allele_Counts) + 1]
  plot(r, n.tot, main = tag, xlab = 'AF', ylab = 'total read', xlim = c(0,xmax), col = col)
}

plot.AF.ref <- function(dat, tag = ''){
  suppressPackageStartupMessages(require(scales))
  eps = .Machine$double.eps
  n.tot = dat$Ref_Allele_Counts + dat$Alt_Allele_Counts + 2*eps
  r = (dat$Ref_Allele_Counts+eps)/n.tot;

  d1 = density(r[r<0.5], bw = 0.001)
  d2 = density(1-r[r>0.5], bw = 0.001)
  ylim = c(0, max(d1$y, d2$y))
  xlim = c(0, max(d1$x, d2$x))
  plot(d1, main = "", xlab = "", ylim = ylim, xlim = xlim)
  par(new = T)
  plot(d2, col = 'red', xlab = 'AF (red ref > 0.5)', main = "", ylim = ylim, xlim = xlim)
  par(new = F)
  plot(min2(r, 1-r), n.tot, col = alpha(c('black', 'red')[(r>0.5)+1], 0.5), main = tag,
       xlab = 'AF (red: Ref > 0.5)', ylab = 'total read', xlim = c(0,0.5))
}

normalize.by.total <- function(x){
  tot = colSums(x);
  tot = tot/mean(tot)
  return(t(t(x)/tot))
}

get.geneticRelatedness.prior <- function(relationship = c('parent', 'child', 'sibling', 
                                                          'uncle/aunt', 'nephew/niece', 'unrelated'), pi = 0.5){
  # 20160928, YFLi
  relationship = match.arg(relationship, several.ok = T)
  priors = matrix(0, 6, 9, dimnames = list(relationship = c('parent', 'child', 'sibling', 
                                                             'uncle/aunt', 'nephew/niece', 'unrelated'),
                                            genotype.combination = c('00', '01', '02',
                                                                     '10', '11', '12',
                                                                     '20', '21', '22')));
  
  if (0){
    priors['parent', ] <- priors['child', ] <- c(8/64, 8/64, 0, 8/64, 16/64, 8/64, 0, 8/64, 8/64);
    priors['sibling', ] = c(9/64, 6/64, 1/64, 6/64, 20/64, 6/64, 1/64, 6/64, 9/64);
    priors['uncle/aunt', ] <- priors['nephew/niece', ] <- c(6/64, 8/64, 2/64, 8/64, 16/64, 8/64, 2/64, 8/64, 6/64);
    priors['unrelated', ] = c(4/64, 8/64, 4/64, 8/64, 16/64, 8/64, 4/64, 8/64, 4/64);
  }else{
    priors['parent', ] <- priors['child', ] <- c((1-pi)^3, pi * (1-pi)^2, 0, 
                                                 pi * (1-pi)^2, pi * (1-pi), pi^2 * (1-pi), 
                                                 0, pi^2 * (1-pi), pi^3);
    priors['sibling', ] = c((1-pi)^2*(1-pi/2)^2, pi * (1-pi)^2 * (1-pi/2), pi^2 * (1-pi)^2/4, 
                            pi * (1-pi)^2*(1-pi/2), pi * (1-pi) * (1+pi-pi^2), pi^2 * (1-pi) * (1/2+pi/2), 
                            pi^2 * (1-pi)^2/4, pi^2 * (1-pi) * (1/2 + pi/2), pi^2*(1/2+pi/2)^2);
    priors['uncle/aunt', ] <- 
      priors['nephew/niece', ] <- c((1-pi)^3*(1-pi/2), pi * (1-pi)^2 * (3/2-pi), pi^2 * (1-pi)^2/2, 
                                    pi * (1-pi)^2*(3/2-pi), pi * (1-pi) * (1/2+2*pi-2*pi^2), pi^2 * (1-pi) * (1/2+pi), 
                                    pi^2 * (1-pi)^2/2, pi^2 * (1-pi) * (1/2 + pi), pi^3*(1/2+pi/2));
    priors['unrelated', ] = c((1-pi)^4, 2*pi * (1-pi)^3, pi^2 * (1-pi)^2, 
                              2*pi * (1-pi)^3, 4*pi^2 * (1-pi)^2, 2*pi^3 * (1-pi), 
                              pi^2 * (1-pi)^2, 2*pi^3 * (1-pi), pi^4);
  }
  priors[relationship,]
}

ErrGT <- function(lambda, GT){
  return(((1- lambda) * GT + lambda * (2- GT)) / 2)
}


expectedAF <- function(GT, lambda, b){
  G = ErrGT(lambda, GT)
  print("in expected AF")
  print(GT)
  return(cbind(p = G %*% t(t(b)), q = (1-G) %*% t(t(b))))
}


runID2Folder <- function(runIDs){
  folders = c()
  for (i in runIDs){
    folders = c(folders, 
                system(paste('ls -d /Volumes/ngulbahce/AnalysisScripts/AlloHeme/Prototype/NGversion/Example_Input/Reports/', i, '*', sep = ''),intern=T))
  }
  names(folders) = runIDs
  return(folders)
}


getGT <- function(multi.runIDs = list(HCT116 = 'A000H00WVM', RKO = 'A000H00WVM'), # allow specifying multiple samples per each baseline
                  multi.ctFileREs = list(HCT116 = '(HCT116|HCT-116).*.txt|tsv', RKO = 'RKO.*.txt|tsv'),
                  multi.files = NULL,
                  default.loci = NULL){
  # user specify cell line names and number of cell lines
  if (length(multi.runIDs)!= length(multi.ctFileREs) || any(names(multi.runIDs)!= names(multi.ctFileREs))){
    stop('Please specify same list of celltypes in multi.runIDs and multi.ctFileREs!')
  }
  if (!is.null(multi.files)){
    celltypes = names(multi.files)
  }else{
    celltypes = names(multi.runIDs)
  }
  GTrep = c()
  ctIDs = c()
  for (celltype in celltypes){
    if ((!is.null(multi.files))){
      files = multi.files[[celltype]]
      runIDs = NULL;
      ctFileREs = NULL
    }else{
      runIDs = multi.runIDs[[celltype]]
      ctFileREs = multi.ctFileREs[[celltype]]
      files = NULL
    }
    
    
    if (all(is.null(runIDs) || is.na(runIDs)) && all(is.null(files)||is.na(files))){ # contributor of unknown genotype
      ctIDs = c(ctIDs, celltype)
      GT = rep(NA, length(default.loci))
      names(GT) = default.loci
      GTrep = cbind(GTrep, GT);
    }else{ # contributor of known genotype
      if (all(is.null(files))){
        if (length(runIDs)>1 & length(ctFileREs)==1){
          ctFileREs = rep(ctFileREs, length(runIDs))
        }
        runFolders = runID2Folder(runIDs)
        for (i in 1:length(runIDs)){
          ctFiles = list.files(runFolders[i], pattern = ctFileREs[i], full.names = T)
          ctIDs = c(ctIDs, rep(celltype, length(ctFiles)));
          if (!length(ctFiles)){
            stop(runIDs, 'Baseline file not found.')
          }
          for (ctFile in ctFiles){
            dat = read.csv(ctFile, sep = '\t', skip = 1)
            dat$Genotype[dat$Genotype == './.'] = NA
            GT = sapply(strsplit(as.character(dat$Genotype), '/'), function(x) 2- sum(as.numeric(x)))
            names(GT) = dat$Locus_ID
            GTrep = cbind.union(GTrep, GT, default = NA)
          }
        }
      }else{
        ctFiles = files
        ctIDs = c(ctIDs, rep(celltype, length(ctFiles)));
        if (!length(ctFiles)){
          stop(ctFiles, 'Baseline file not found.')
        }
        for (ctFile in ctFiles){
          dat = read.csv(ctFile, sep = '\t', skip = 1)
          dat$Genotype[dat$Genotype == './.'] = NA
          GT = sapply(strsplit(as.character(dat$Genotype), '/'), function(x) 2- sum(as.numeric(x)))
          names(GT) = dat$Locus_ID
          GTrep = cbind.union(GTrep, GT, default = NA)
        }
      }
    }
  }
  colnames(GTrep) = ctIDs;
  uctIDs = unique(ctIDs);
  GTs = c();
  GTpasses = c();
  for (ctID in uctIDs){
    consensus = rowMeans(GTrep[,ctIDs == ctID, drop=F], na.rm = T)
    consistency = apply(GTrep[,ctIDs == ctID, drop=F] == consensus, 1, all, na.rm=T)
    GTs = cbind.union(GTs, consensus, default = NA)
    GTpasses = cbind.union(GTpasses, consistency)
  }
  colnames(GTpasses) <- colnames(GTs) <- uctIDs
  failed = which(rowSums(!GTpasses)>0)
  if (length(failed)){
    warning('Some SNPs received conflicting genotyping in different samples')
    cat(failed)
  }
  GTs = round(GTs)
  GTs[is.nan(GTs)] = NA
  return(list(GTs = GTs, GTpasses = GTpasses))
}


PUGT.expandGT <- function(GTs, prior = NULL){
  # 1) expand partially unknown genotypes
  # 2) compute the posterior probability
  D = ncol(GTs)
  contributorNames = colnames(GTs)
  GTexp = matrix(0, nrow(GTs), 3^D)
  allGTs = names(prior)
  print("prior")
  print(prior)
  print("allGTs")
  print(allGTs)
  colnames(GTexp) = allGTs
  rownames(GTexp) = rownames(GTs)
  print("GTexp before")
  print(GTexp)
  for (g in 1:nrow(GTs)){
    print("class of GTs")
    print(class(GTs))
    G1 = GTs[g,] # expand degenerate genotypes for each loci
    print("G1")
    print(G1)
    allGT.ob = list()
    for (d in 1:D){
      if (is.na(G1[d])){
        Gs = 0:2
      }else{
        Gs = G1[d]
      }
      allGT.ob[[d]] = paste(contributorNames[d], Gs, sep = ':')
    }
    print("allGT.ob")
    print(allGT.ob)
    allGT.obs = apply(as.matrix(expand.grid(allGT.ob)), 1, FUN = paste, collapse = '|')
    print("allGT.obs")
    print(allGT.obs)
    priorNormalize = prior[allGT.obs];
    priorNormalize = priorNormalize/sum(priorNormalize);
    GTexp[g, allGT.obs] = priorNormalize
  }
  print("GTexp after")
  print(GTexp)
  return(GTexp)
}


PUGT.Prior <- function(contributorNames = contributors, relationships = 'unrelated', pi = 0.5){
  # compute the prior for 3D their corresponding prior probabilities. 
  # Some of the probabilities can be 0
  if (relationships != 'unrelated'){
    stop('relationship', relationships, 'is not implemented')
  }
  
  allGT <- allGTID <- c()
  prior1 = c(pi^2, 2*pi*(1-pi), (1-pi)^2)
  print(prior1)
  allprior = c()
  for (d in 1:length(contributorNames)){
    allGTID = cbind(allGTID, paste(contributorNames[d], 0:2, sep = ':'))
    allGT = cbind(allGT, 0:2)
    #allprior = cbind(allprior, prior1)
    allprior = get.geneticRelatedness.prior(relationship="child",pi=0.5)
  }
  print(allprior)
  colnames(allGT) = contributorNames
  allGTIDs = apply(as.matrix(expand.grid(as.data.frame(allGTID))), 1, FUN = paste, collapse = '|')
  allpriors = cbind(prior = apply(as.matrix(expand.grid(as.data.frame(allprior))), 1, FUN = prod), as.matrix(expand.grid(as.data.frame(allGT))))
  rownames(allpriors) = allGTIDs;
  print("allpriors")
  print(allpriors)

  return(allpriors)
}

PUGT.LL <- function(beta, r, lambda, model = c('EPS','S','PS','ES'), 
                    expandGT, # population allele frequency, assume to be the same for all SNPs
                    prior, ni, n1, n2){
  # log base 2
  suppressPackageStartupMessages(require(extraDistr))
  eAF = expectedAF(GT = prior[,2:ncol(prior)], lambda = lambda, b = beta) # compute the expected allele fraction for 3^D possible genotype. Return a vector of length 3^D.
  n = n1 + n2;
  print("n1,n2")
  print(n1)
  print(n2) 
  mu.n = mean(n);
  if (model == 'S'){
    BB <- function(n1, n2, a, b) log2(a)*n1 + log2(b)*n2; # binomial model
  }else if (model == 'PS'){
    BB <- function(n1, n2, a, b) (lbeta(n1 + ni*(1+r)/(1-r)*a, n2+ni*(1+r)/(1-r)*b) - lbeta(ni*(1+r)/(1-r)*a, ni*(1+r)/(1-r)*b))/log(2);
  }else if (model == 'ES'){
    BB <- function(n1, n2, a, b) (lbeta(n1 + (ni-1)*a, n2+(ni-1)*b) - lbeta((ni-1)*a, (ni-1)*b))/log(2);
  }else if (model == 'EPS'){
    print("r in EPS")
    print(r)
    BB <- function(n1, n2, a, b) (lbeta(n1 + ni*(1+r)/2*a, n2+ni*(1+r)/2*b) - lbeta(ni*(1+r)/2*a, ni*(1+r)/2*b))/log(2);
  }
  
  coeff <- (lgamma(n+1)-lgamma(n1+1)-lgamma(n2+1))/log(2);  ## add the constant coefficents to get correct likelihood values
  
  ll.mat = matrix(-Inf, nrow = length(n1), ncol = nrow(eAF))
  for (i in 1:nrow(eAF)){
    ll.mat[,i] = BB(n1, n2, eAF[i,1], eAF[i,2])
  }
  colnames(ll.mat) = rownames(eAF);
  rownames(ll.mat) = rownames(expandGT)
  
  ll.mat[expandGT==0] = -Inf;
  ll.mat[expandGT>0] = log2(expandGT[expandGT>0]) + ll.mat[expandGT>0]; ## add the genetic prior
  ll = rowLogSumExp(ll.mat);
  ll = ll + coeff
  return(sum(ll))
}


argmin <- function(x, fun){
  min.fun = Inf;
  if (!is.matrix(x)){
    x = matrix(x, nrow = length(x), ncol = 1)
  }
  for (i in 1:nrow(x)){
    xx = x[i,]
    y = fun(xx)
    if (y < min.fun){
      min.x = xx;
      min.fun = y;
    }
  }
  return(min.x)
}


mapping.bias <- function(dat){
  r = dat$Ref_Allele_Counts/(dat$Ref_Allele_Counts+dat$Alt_Allele_Counts)
  cutoffs = c(0.05, 0.1, 0.25, 0.5)
  braw <- br <- cutoffs; names(braw) <- names(br) <- cutoffs
  for (c1 in cutoffs){
    braw[paste(c1)] = sum(dat$Ref_Allele_Counts[r<c1], na.rm = T)/sum(dat$Alt_Allele_Counts[r > 1- c1], na.rm = T)
    br[paste(c1)] = sum(r[r<c1], na.rm = T)/sum(1-r[r > 1- c1], na.rm = T)
  }
  braw = c(braw, allRef_Alt = sum(dat$Ref_Allele_Counts, na.rm=T)/sum(dat$Alt_Allele_Counts, na.rm=T), 
           lambda.Alt2Ref = sum(dat$Ref_Allele_Counts[r < 0.25], na.rm=T) /sum(dat$Alt_Allele_Counts[r < 0.25], na.rm=T),
           lambda.Ref2Alt = sum(dat$Alt_Allele_Counts[r > 0.75], na.rm=T)/sum(dat$Ref_Allele_Counts[r > 0.25], na.rm=T))
  return(list(braw = braw, br = br))
}


estimate.noise <- function(dat){
  # estimate the sequencing error (point error)
  epsilon = .Machine$double.eps
  is.noise = dat$Other_Allele_Counts/(rowSums(dat[,2:4])+epsilon) < 0.01
  return(sum(dat$Other_Allele_Counts[is.noise]) /(sum(colSums(dat[is.noise,2:4]))+epsilon) / 2 + epsilon) 
  # 1/2 total error rate to 2 SNP alleles are the same as mutation rate from 1 allele to another allele
  # + epsilon to avoid zeros sequencing error rate, which leads to numerical errors
}


target.QC <- function(dat, cutoff = 50){
  # cutoff was 50, 
  # it was then changed to 100, 
  epsilon = .Machine$double.eps
  is.3allele = dat$Other_Allele_Counts/(rowSums(dat[,2:4])+epsilon) > 0.01 # the 3rd allele is more than 1%
  is.lowcoverage = rowSums(dat[,2:3]) < cutoff # less than 50 read
  return(!(is.3allele | is.lowcoverage))
}


get.ni <- function(H, h, input.DNA.ng, genome.weight){
  # convert input DNA amount and length to effective input DNA amount
  e_cfDNA = (H-h+1)/H
  ni = e_cfDNA * input.DNA.ng/genome.weight # amplifiable copies of cfDNA
  names(ni) = NULL
  return(ni)
}


BMT.pred <- function(dat, GTdat, 
                     prior.relationship = c('parent', 'child', 'sibling',  # donor's relationship to recipient
                                            'uncle/aunt', 'nephew/niece', 'unrelated'),
                     prior.pi = 0.5, 
                     exclusion.list = NULL,
                     DNA.type = c('cfDNA', 'horizon', 'gDNA'),
                     H = c(cfDNA =  165, horizon = 160, gDNA = 1E5), # average cfDNA or gDNA length
                     h = 110, # average amplicon length
                     model = c('KGT.NaiveLM', 
                               'KGT.IterLM',
                               'KGT.seq',
                               'KGT.EPS',
                               'PUGT.EPS',
                               'PUGT.S',
                               'PUGT.seq', # same as PUGT.S
                               'PUGT.PS',
                               'PUGT.ES'),
                     logsig.transform = F, 
                     c = 1, r = 1, optimize.c = T,
                     input.DNA.ng = 8, # ng
                     genome.weight = 3.59E-3, # ng/haploid
                     duplicate.adjust = F, error.modeling = c('estimate.apriori','predict', 'none')){
  # 20160919
  suppressPackageStartupMessages(require(stats4))
  prior.relationship = match.arg(prior.relationship)
  DNA.type = match.arg(DNA.type)
  target.keep = target.QC(dat);
  CoverageMu = mean((dat$Alt_Allele_Counts + dat$Ref_Allele_Counts)[target.keep],na.rm=T)
  QC.metrics = c(CoverageMu = CoverageMu, 
                 TargetPassed = sum(target.keep),
                 SeqErr = estimate.noise(dat),
                 MappingBias.bfration = mapping.bias(dat)$br,
                 MappingBias.braw = mapping.bias(dat)$braw,
                 PerTargetLL = 0,
                 FisherCV = 0,
                 ExpectedFisherCV = 0)
  if (QC.metrics['TargetPassed'] < 20){
    warning(paste('# Targets after filtering', QC.metrics['TargetPassed'], '< 20'))
  }else  if (QC.metrics['CoverageMu'] < 100){
    warning(paste('# Targets after filtering', QC.metrics['CoverageMu'], '< 100'))
  }
  
  
  dat = dat[target.keep,]
  model = match.arg(model)
  if (model == 'PS'){
    model = 'PCR-seq'
  }else if (model == 'ES'){
    model = 'cfDNA-seq'
  }else if (model == 'EPS'){
    model = 'cfDNA-PCR-seq'
  }
  error.modeling = match.arg(error.modeling)
  
  ni = get.ni(H[DNA.type], h, input.DNA.ng, genome.weight)
  cat('Input ', input.DNA.ng, 'ng; Copy Number ', round(ni,digits = 1), '\n\n')
  
  eps = .Machine$double.eps
  if (error.modeling == 'predict'){
    stop('Not implemented for KGT.NaiveLM')
  }else if (error.modeling == 'estimate.apriori'){
    lambda = estimate.noise(dat)
    print("Error estimate via other counts")
    print(lambda)
  }else{
    lambda = 1E-14
  }
  
  if (regexpr('PUGT.',model)>0){
    GTs = GTdat$GTs
    contributors = colnames(GTs)
    prior = PUGT.Prior(contributorNames = contributors, relationships = 'unrelated', pi = 0.5) ## if the genotypes is fully known, the prior will have no impact
    expandGT = PUGT.expandGT(GTs, prior[,1]) 
    n1 =  dat$Ref_Allele_Counts + eps; n2 = dat$Alt_Allele_Counts + eps
    names(n1) <- names(n2) <- dat$Locus_ID
    loci = intersect(rownames(expandGT), names(n1));
    expandGT = expandGT[loci,]
    n1 = n1[loci]
    n2 = n2[loci]
    npar = ncol(GTs)
    
    if (1){ # initialization with a global search
      if (npar > 5){
        warning('number of contributor greater than 5, it will be slow to find a good inital point for optimization')
      }
      if (npar == 2){
        x = c(-80, -40, -20, seq(-10, 10, 0.1), 20, 40, 80)
      }else{
        x = c(seq(-9, 9, 3.4));
      }
      xgrid = as.data.frame(matrix(x, nrow = length(x), ncol = npar -1));
      xgrid = as.matrix(cbind(r = 1.222358, expand.grid(xgrid)))
      min.x = argmin(xgrid, function(x)-PUGT.LL(b = softmax(x[2:npar]), r = logsig(x[1]), lambda = lambda, 
                                                model = sub('PUGT.', '',model), 
                                                expandGT, # population allele frequency, assume to be the same for all SNPs
                                                prior, ni, n1 =n1, n2 = n2))
    }
    soptim = optim(min.x, fn=function(x) -PUGT.LL(b = softmax(x[2:npar]), r = logsig(x[1]), lambda = lambda, 
                                                  model = sub('PUGT.', '',model), 
                                                  expandGT, # population allele frequency, assume to be the same for all SNPs
                                                  prior, ni, n1 =n1, n2 = n2),
                   method='BFGS', hessian=TRUE)
    beta <- softmax(soptim$par[2:npar]) # logsig(alpha$par[1])
    r = logsig(soptim$par[1])
    J <- MASS::ginv(gradient.softmax(soptim$par[2:npar])) # tranformation matrix 
    J = rbind(0,cbind(0,J)); J[1,1] = 1/(logsig(soptim$par[1])*logsig(-soptim$par[1])) # Do not use J[1,1] = 1/(r*(1-r)), so as to avoid the numeric rounding error of 1-r = 1 when r is small
    J[is.infinite((J))] = 10E10
    soptim$hessian = t(J) %*% soptim$hessian %*% J
    stderr = sqrt(MASS::ginv(soptim$hessian)[which(diag(nrow = nrow(soptim$hessian))>0)])
    stderr[is.na(stderr)] = 0
    
    logL = -soptim$value
    names(beta) <-  colnames(GTs)
    names(stderr) <- c('r', colnames(GTs))
  }else{
    stop(model, 'not implemented')
  }
  
  L_info = get.NumInfoLoci(dat, prior.relationship, GTs = GTs, beta = beta)
  Leff = beta * (1-beta) / stderr[names(beta)] ^2 / (ni * (1-exp(-CoverageMu/ni)))  # effective number of loci
  out = c(beta = beta, logL = logL, 
          CI95.low = beta + qnorm(0.025) * stderr[names(beta)], 
          CI95.high = beta + qnorm(0.975) * stderr[names(beta)],
          lambda = lambda, ni = ni, c = c, PCR.amplication.rate = r, stderr = stderr,
          CV_b = stderr[names(beta)]/beta, L_info = L_info, L_eff = Leff, CoE = Leff/L_info);
  
  QC.metrics['PerTargetLL'] = out['logL']/QC.metrics['TargetPassed'];
  return(list(pred = out, dat.filter = dat, QC = QC.metrics))
  
}


get.NumInfoLoci <- function(dat, relationship, GTs = NULL, beta = NULL){
  # compute the number of informative loci based on the known genotypes or
  # estimate the number of informative loci based on genetic relationship
  ceoff = c('sibling' = 0.21875, 'parent' = 0.25,  'child' = 0.25, 'uncle' = 0.3125, 'aunt' = 0.3125, 
            'nephew' = 0.3125, 'niece' = 0.3125, grandparent = NA, grandchild = NA, cousin =  NA,halfsibling = NA,
            'unrelated' = 0.375)
  if (! relationship %in% names(ceoff)){
    stop('unknown relationship type')
  }
  if (relationship %in% c('grandparent','grandchild', 'cousin', 'halfsibling')){
    stop(relationship, 'not implemented')
  }
  
  if (is.null(GTs)){ 
    Linfo = c(estimate = nrow(dat) * ceoff[relationship],
              exact = NA); # note that this is an estimate of the number of effective loci
  }else{ # exact number of informative SNPs
    idx.major = which.max(beta) # major contributor
    # when  some of the genotype is unknown (NAs), the # will be NA automatically
    LInfo_exact = colMeans(((GTs[,idx.major] == 0) * (GTs != 0)) + ((GTs[,idx.major] == 2) * (GTs != 2)), na.rm = T) * nrow(GTs)
    Linfo = c(estimate = nrow(dat) * ceoff[relationship],
              exact = LInfo_exact[setdiff(1:length(LInfo_exact), idx.major)])
  }
}
