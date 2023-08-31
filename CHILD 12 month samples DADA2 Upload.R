library(BiocManager)
#BiocManager::install()
library(msa)
library(dada2); packageVersion("dada2")
library(plyr)
library(reshape2)
library(ggplot2)
library(Biostrings)
library(doParallel)
library(ggalluvial)
library(ggrepel)
library(ggdendro)
library(MASS)
library(rentrez)
library(dendextend)
library(R.utils)
library(dplyr)
library(ggnewscale)
library(vegan)
library(taxize)
library(RGCCA)
library(httr)
library(Rfast)

# ==== collapseNoMismatch2 ====
collapseNoMismatch2 = function (seqtab, minOverlap = 20, orderBy = "abundance", identicalOnly = FALSE, 
                                vec = TRUE, band = -1, verbose = FALSE) {
  dupes <- duplicated(colnames(seqtab))
  if (any(dupes)) {
    st <- seqtab[, !dupes, drop = FALSE]
    for (i in which(dupes)) {
      sq <- colnames(seqtab)[[i]]
      st[, sq] <- st[, sq] + seqtab[, i]
    }
    seqtab <- st
  }
  if (identicalOnly) {
    return(seqtab)
  }
  
  unqs.srt <- sort(getUniques(seqtab), decreasing = TRUE)
  seqs <- names(unqs.srt)
  seqs.prefix <- substr(seqs, 1, minOverlap)
  # print(length(unique(seqs.prefix)))
  seqs.unused <- rep.int(TRUE,length(seqs))
  seqs.out <- rep.int(FALSE,length(seqs))
  collapsed <- matrix(0L, nrow = nrow(seqtab), ncol = ncol(seqtab))
  colnames(collapsed) <- colnames(seqtab)
  rownames(collapsed) <- rownames(seqtab)
  
  while (any(seqs.unused)) {
    # first query is always added
    queryNum = which.max(seqs.unused)
    print(queryNum)
    query = seqs[queryNum]
    
    query.prefix <- seqs.prefix[queryNum]
    seqs.out[queryNum] <-TRUE
    seqs.unused[queryNum] <- FALSE
    collapsed[, query] <- seqtab[, query]
    if(!any(seqs.unused)){
      break
    }
    collapseCandidates = rep.int(FALSE,length(seqs))
    collapseCandidates[seqs.unused] = grepl(
      query.prefix, seqs[seqs.unused], fixed = TRUE) | 
      sapply(seqs.prefix[seqs.unused], function(x) grepl(x, query, fixed = TRUE)
      )
    seqs.used = FALSE
    if(any(collapseCandidates)){
      seqs.used = sapply(seqs[collapseCandidates], function(ref){
        if (nwhamming(query, ref, vec = vec, band = band) == 0) {
          collapsed[, query] <<- collapsed[, query] + seqtab[, ref]
          return(TRUE)
        }else{
          return(FALSE)
        }
      })
      
      seqs.unused[collapseCandidates] = !seqs.used
    }
  }
  
  collapsed <- collapsed[,seqs.out,drop = FALSE]
  if (!is.null(orderBy)) {
    if (orderBy == "abundance") {
      collapsed <- collapsed[, order(colSums(collapsed),
                                     decreasing = TRUE), drop = FALSE]
    }
    else if (orderBy == "nsamples") {
      collapsed <- collapsed[, order(colSums(collapsed >
                                               0), decreasing = TRUE), drop = FALSE]
    }
  }
  collapsed <- collapsed[, order(colSums(collapsed), decreasing = TRUE),
                         drop = FALSE]
  if (verbose) 
    message("Output ", ncol(collapsed), " collapsed sequences out of ", 
            ncol(seqtab), " input sequences.")
  collapsed
}


# ==== multiplot ====
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


# ==== DADA2 Processing Pipeline function ====
process = function(path,
                   fnFs, fnRs, sample.names, errPoolName, orientFR.split = FALSE,
                   trimLeftSelect, truncLenSelect, filter.matchIDs = FALSE,
                   ErrModelMonotonicity = FALSE,
                   OMEGA_A = getDadaOpt(option = "OMEGA_A"), pool = FALSE, poolList = character(0),
                   priorsF = character(0), priorsR = character(0),
                   
                   plotToggle = FALSE,
                   printPriorHead = FALSE,
                   
                   preFilteredToggle = TRUE,
                   preCalcErrToggle = TRUE,
                   preDerepToggle = TRUE,
                   preDenoiseToggle = TRUE,
                   
                   saveFiltered = TRUE,
                   saveDereplicated = TRUE,
                   saveDenoised = TRUE){
  # ---- filter and trimming ----
  cat("Filter and trimming...\n")
  dir.create(path = paste0(path,"/filtered/"), showWarnings = TRUE)
  
  # Set destination to filtered/ subdirectory
  filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
  filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
  
  if(!preFilteredToggle){
    #filter and trim function
    #NOTE: consider relaxing maxEE if needed
    out <- filterAndTrimWinPara(fnFs, filtFs, fnRs, filtRs, truncLen = truncLenSelect, trimLeft = trimLeftSelect,
                                maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                                compress=TRUE, multithread=TRUE, matchIDs=filter.matchIDs)
    
    # creating sequence tracking table
    track <- cbind(out, NA, NA, NA, NA)
    colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
    rownames(track) <- sample.names
    
    print(track)
    
  }else{
    # creating sequence tracking table
    track <- cbind(rep(NA,length(sample.names)), rep(NA,length(sample.names)), rep(NA,length(sample.names)), rep(NA,length(sample.names)))
    colnames(track) <- c("denoisedF", "denoisedR", "merged", "nonchim")
    rownames(track) <- sample.names
  }
  
  #plotting quality score profiles
  if(plotToggle){
    plotQualityProfile(filtFs,aggregate = TRUE)
    plotQualityProfile(filtRs,aggregate = TRUE)
  }
  
  # ---- generating error models ----
  cat("\n")
  cat("Generating error models...\n")
  dir.create(path = paste0(path,"/Error Models"), showWarnings = TRUE)
  
  if(!preCalcErrToggle){
    if(orientFR.split){
      #learning error models for split samples
      orientFerrF <- learnErrors(filtFs[grep(".orientF", filtFs, fixed = TRUE)], multithread=TRUE)
      orientFerrR <- learnErrors(filtRs[grep(".orientF", filtRs, fixed = TRUE)], multithread=TRUE)
      orientRerrF <- learnErrors(filtFs[grep(".orientR", filtFs, fixed = TRUE)], multithread=TRUE)
      orientRerrR <- learnErrors(filtRs[grep(".orientR", filtRs, fixed = TRUE)], multithread=TRUE)
      
      if(ErrModelMonotonicity){
        #set error values for Q-scores<40 to equal error values for Q-scores=40
        #orientFerrF
        orientFerrFOutMono = (getErrors(orientFerrF))
        orientFerrFOutMono = apply(orientFerrFOutMono, 2, FUN = function(x){
          y = x
          print(y)
          y[y<orientFerrFOutMono[,40]] = orientFerrFOutMono[,40][y<orientFerrFOutMono[,40]]
          print(y)
          return(y)
        })
        orientFerrF$err_out = orientFerrFOutMono
        
        #orientFerrR
        orientFerrROutMono = (getErrors(orientFerrR))
        orientFerrROutMono = apply(orientFerrROutMono, 2, FUN = function(x){
          y = x
          print(y)
          y[y<orientFerrROutMono[,40]] = orientFerrROutMono[,40][y<orientFerrROutMono[,40]]
          print(y)
          return(y)
        })
        orientFerrR$err_out = orientFerrROutMono
        
        #orientRerrF
        orientRerrFOutMono = (getErrors(orientRerrF))
        orientRerrFOutMono = apply(orientRerrFOutMono, 2, FUN = function(x){
          y = x
          print(y)
          y[y<orientRerrFOutMono[,40]] = orientRerrFOutMono[,40][y<orientRerrFOutMono[,40]]
          print(y)
          return(y)
        })
        orientRerrF$err_out = orientRerrFOutMono
        
        #orientRerrR
        orientRerrROutMono = (getErrors(orientRerrR))
        orientRerrROutMono = apply(orientRerrROutMono, 2, FUN = function(x){
          y = x
          print(y)
          y[y<orientRerrROutMono[,40]] = orientRerrROutMono[,40][y<orientRerrROutMono[,40]]
          print(y)
          return(y)
        })
        orientRerrR$err_out = orientRerrROutMono
      }
      
      #saving error model files as .rds
      saveRDS(orientFerrF,file = paste0(path,"/Error Models/",errPoolName,".orientF","_errF.rds"))
      saveRDS(orientFerrR,file = paste0(path,"/Error Models/",errPoolName,".orientF","_errR.rds"))
      saveRDS(orientRerrF,file = paste0(path,"/Error Models/",errPoolName,".orientR","_errF.rds"))
      saveRDS(orientRerrR,file = paste0(path,"/Error Models/",errPoolName,".orientR","_errR.rds"))
      
    }else{
      #learning error models
      errF <- learnErrors(filtFs, multithread=TRUE)
      errR <- learnErrors(filtRs, multithread=TRUE)
      
      if(ErrModelMonotonicity){
        #set error values for Q-scores<40 to equal error values for Q-scores=40
        #ErrF
        errFOutMono = (getErrors(errF))
        errFOutMono = apply(errFOutMono, 2, FUN = function(x){
          y = x
          print(y)
          y[y<errFOutMono[,40]] = errFOutMono[,40][y<errFOutMono[,40]]
          print(y)
          return(y)
        })
        
        errF$err_out = errFOutMono
        
        #ErrR
        errROutMono = (getErrors(errR))
        errROutMono = apply(errROutMono, 2, FUN = function(x){
          y = x
          print(y)
          y[y<errROutMono[,40]] = errROutMono[,40][y<errROutMono[,40]]
          print(y)
          return(y)
        })
        
        errR$err_out = errROutMono
      }
      
      #saving error model files as .rds
      saveRDS(errF,file = paste0(path,"/Error Models/",errPoolName,"_errF.rds"))
      saveRDS(errR,file = paste0(path,"/Error Models/",errPoolName,"_errR.rds"))
    }    
  }else{
    if(orientFR.split){
      orientFerrF = readRDS(file = paste0(path,"/Error Models/",errPoolName,".orientF_errF.rds"))
      orientFerrR = readRDS(file = paste0(path,"/Error Models/",errPoolName,".orientF_errR.rds"))
      orientRerrF = readRDS(file = paste0(path,"/Error Models/",errPoolName,".orientR_errF.rds"))
      orientRerrR = readRDS(file = paste0(path,"/Error Models/",errPoolName,".orientR_errR.rds"))
    }else{
      errF = readRDS(file = paste0(path,"/Error Models/",errPoolName,"_errF.rds"))
      errR = readRDS(file = paste0(path,"/Error Models/",errPoolName,"_errR.rds"))
    }
  }
  
  #plotting error model summary graphs
  if(plotToggle){
    if(orientFR.split){
      print(plotErrors(orientFerrF, nominalQ=TRUE))
      print(plotErrors(orientFerrR, nominalQ=TRUE))
      print(plotErrors(orientRerrF, nominalQ=TRUE))
      print(plotErrors(orientRerrR, nominalQ=TRUE))
    }else{
      print(plotErrors(errF, nominalQ=TRUE))
      print(plotErrors(errR, nominalQ=TRUE))
    }
  }
  
  # ---- dereplicating filtered FASTQ files ----
  cat("\n")
  cat("Dereplicating sequences...\n")
  dir.create(path = paste0(path,"/dereplicated"), showWarnings = TRUE)
  
  if(!preDerepToggle){
    # dereplication is done multicore parallel, if possible
    # generating parallel backend
    nodes <- detectCores()
    cl <- makeCluster(nodes, type = "PSOCK")
    registerDoParallel(cl)
    
    # dereplicating fwd reads
    derepFs = aaply(filtFs,1,function(x, filePath = path){
      sampleName = laply(strsplit(x, "filtered/", fixed = TRUE), function(y) y[2])
      sampleName = laply(strsplit(sampleName, "_F_filt.fastq.gz", fixed = TRUE), function(y) y[1])
      
      derepF <- derepFastq(x, verbose=TRUE)
      saveRDS(derepF,file = paste0(filePath,"/dereplicated/",sampleName,"_derepF.rds"))
      return(paste0(filePath,"/dereplicated/",sampleName,"_derepF.rds"))
    }, filePath = path, .parallel = TRUE, .paropts = list(.packages = "dada2"))
    
    # dereplicating rev reads
    derepRs = aaply(filtRs,1,function(x, filePath = path){
      sampleName = laply(strsplit(x, "filtered/", fixed = TRUE), function(y) y[2])
      sampleName = laply(strsplit(sampleName, "_R_filt.fastq.gz", fixed = TRUE), function(y) y[1])
      
      derepR <- derepFastq(x, verbose=TRUE)
      saveRDS(derepR,file = paste0(filePath,"/dereplicated/",sampleName,"_derepR.rds"))
      return(paste0(filePath,"/dereplicated/",sampleName,"_derepR.rds"))
    }, filePath = path, .parallel = TRUE, .paropts = list(.packages = "dada2"))
    
    #stopping parallel backend
    stopCluster(cl)
    
    # Name the derep-class objects by the sample names
    names(derepFs) <- sample.names
    names(derepRs) <- sample.names
    
    # Remove filtered files if saveFiltered == FALSE and the process() generated filtered files (preFilteredToggle == FALSE)
    if(!saveFiltered & !preFilteredToggle){
      file.remove(filtFs,filtRs)
    }
  }else{
    derepFs = paste0(path,"/dereplicated/",sample.names,"_derepF.rds")
    derepRs = paste0(path,"/dereplicated/",sample.names,"_derepR.rds")
  }
  
  # ---- denoising dereplicated files ----
  cat("\n")
  cat("Denoising sequences...\n")
  dir.create(path = paste0(path,"/dada denoised"), showWarnings = TRUE)
  
  if(!preDenoiseToggle){
    dadaFs = data.frame(name = rep(NA,length(fnFs)),
                        count = rep(NA,length(fnFs)))
    
    dadaRs = data.frame(name = rep(NA,length(fnRs)),
                        count = rep(NA,length(fnRs)))
    
    # if pool != FALSE, denoise pooled samples according to poolList,
    # else set outPoolListNamesMatch (indices of non-pooled samples) to all samples
    if(pool != FALSE){
      
      # check if poolList is NULL
      if(is.null(poolList)){
        stop(("If pooling, poolList must not be NULL."))
      }
      
      # all samples in sample.names will be processed; 
      # samples in poolList will be pooled according to poolList,
      # samples not in poolList will default to non-pooled denoising
      
      # determining samples in/out poolList
      inPoolListNamesMatch = unlist(poolList, use.names = FALSE)
      
      outPoolListNamesMatch = which(!(sample.names %in% inPoolListNamesMatch))
      inPoolListNamesMatch = which(sample.names %in% inPoolListNamesMatch)
      
      a_ply(names(poolList),1, function(dadaPoolName){
        cat(paste0("denoising pool '",dadaPoolName,"'...\n"))
        
        poolSampleNames = poolList[[dadaPoolName]]
        poolSampleNamesMatch = match(poolSampleNames,sample.names)
        
        poolDerepFs = derepFs[poolSampleNamesMatch]
        poolDerepRs = derepRs[poolSampleNamesMatch]
        
        print(poolSampleNames)
        
        if(!isEmpty(priorsF)){
          priorMatchDf = adply(poolDerepFs,1, .id = NULL, function(x){
            sampleName = laply(strsplit(x, "dereplicated/", fixed = TRUE), function(y) y[2])
            sampleName = laply(strsplit(sampleName, "_derepF.rds", fixed = TRUE), function(y) y[1])
            
            derepFFile = readRDS(x)
            
            if(printPriorHead){
              # printing some dereplicated fwd reads and fwd priors to help troubleshoot trimming priors
              cat("head of dereplicated Fwd sequences:\n")
              print(head(names(derepFFile$uniques)))
              cat("head of prior Fwd sequences:\n")
              print(head(priorsF[[dadaPoolName]]))
              
              # matching fwd priors and fwd dereplicated reads for troubleshooting/sanity check
              cat("\n")
              cat(paste0(sum(names(derepFFile$uniques) %in% priorsF[[dadaPoolName]]), " Fwd dereplicated priors matches out of ", length(priorsF[[dadaPoolName]]), " Fwd priors\n"))
              cat(paste0(sum(derepFFile$uniques[names(derepFFile$uniques) %in% priorsF[[dadaPoolName]]]), " priors matches found out of ", sum(derepFFile$uniques), " Fwd sequences\n"))
            }
            
            derepInPriorsCount = sum(names(derepFFile$uniques) %in% priorsF[[dadaPoolName]])
            priorsCount = length(priorsF[[dadaPoolName]])
            derepInPriorsReadCount = sum(derepFFile$uniques[names(derepFFile$uniques) %in% priorsF[[dadaPoolName]]])
            readCount = sum(derepFFile$uniques)
            
            # same data as above in tabular format
            df = data.frame(sampleName = sampleName,
                            derepInPriorsCount = derepInPriorsCount,
                            priorsCount = priorsCount,
                            derepInPriorsReadCount = derepInPriorsReadCount,
                            readCount = readCount)
            
            return(df)
          },.progress = "text")
          print(priorMatchDf)
        }else{
          cat("no fwd priors\n")
        }
        
        # Loading Derep files into a list
        poolDerepFFileList = alply(poolDerepFs,1,function(x){
          derepFFile = readRDS(x)
        })
        
        # Setting error model to be used
        if(orientFR.split){
          is.orientR = grepl(".orientR", poolSampleNames, fixed = TRUE)
          if(!(all(is.orientR) | all(!is.orientR))){
            stop("If orientFR.split and pooling, pools should contain only samples with same orientation")
          }
          is.orientR = all(is.orientR)
          
          if(!is.orientR){
            dadaErrF = orientFerrF
          }else{
            dadaErrF = orientRerrF
          }
        }else{
          dadaErrF = errF
        }
        
        # Setting priors to be used
        if(!isEmpty(priorsF)){
          dadaPriorsF = priorsF[[dadaPoolName]]
          if(is.null(dadaPriorsF)){
            dadaPriorsF = character(0)
          }
        }else{
          dadaPriorsF = character(0)
        }
        
        #DADA denoising
        dadaFList <- dada(poolDerepFFileList, err=dadaErrF, priors = dadaPriorsF, pool=pool, OMEGA_A = OMEGA_A, multithread=TRUE)
        
        for(i in 1:length(dadaFList)){
          saveRDS(dadaFList[[i]],file = paste0(path,"/dada denoised/",poolSampleNames[i],"_dadaF.rds"))
          dadaFs[poolSampleNamesMatch[i],1] <<- paste0(path,"/dada denoised/",poolSampleNames[i],"_dadaF.rds")
          dadaFs[poolSampleNamesMatch[i],2] <<- sum(getUniques(dadaFList[[i]]))
        }
        
        track[,"denoisedF"] = as.numeric(dadaFs[,2])
        cat("\n")
        print(track[poolSampleNamesMatch,])
        
        if(!isEmpty(priorsR)){
          priorMatchDf = adply(poolDerepRs,1, .id = NULL, function(x){
            sampleName = laply(strsplit(x, "dereplicated/", fixed = TRUE), function(y) y[2])
            sampleName = laply(strsplit(sampleName, "_derepR.rds", fixed = TRUE), function(y) y[1])
            
            derepRFile = readRDS(x)
            
            if(printPriorHead){
              # printing some dereplicated rev reads and rev priors to help troubleshoot trimming priors
              cat("\n")
              cat("head of dereplicated Rev sequences:\n")
              print(head(names(derepRFile$uniques)))
              cat("head of prior Rev sequences:\n")
              print(head(priorsR[[dadaPoolName]]))
              
              # matching rev priors and rev dereplicated reads for troubleshooting/sanity check
              cat("\n")
              cat(paste0(sum(names(derepRFile$uniques) %in% priorsR[[dadaPoolName]]), " Rev dereplicated priors matches out of ", length(priorsR[[dadaPoolName]]), " Rev priors\n"))
              cat(paste0(sum(derepRFile$uniques[names(derepRFile$uniques) %in% priorsR[[dadaPoolName]]]), " priors matches found out of ", sum(derepRFile$uniques), " Rev sequences\n"))
            }
            
            derepInPriorsCount = sum(names(derepRFile$uniques) %in% priorsR[[dadaPoolName]])
            priorsCount = length(priorsR[[dadaPoolName]])
            derepInPriorsReadCount = sum(derepRFile$uniques[names(derepRFile$uniques) %in% priorsR[[dadaPoolName]]])
            readCount = sum(derepRFile$uniques)
            
            # same data as above in tabular format
            df = data.frame(sampleName = sampleName,
                            derepInPriorsCount = derepInPriorsCount,
                            priorsCount = priorsCount,
                            derepInPriorsReadCount = derepInPriorsReadCount,
                            readCount = readCount)
            
            return(df)
          },.progress = "text")
          print(priorMatchDf)
        }else{
          cat("no rev priors\n")
        }
        
        # Loading Derep files into a list
        poolDerepRFileList = alply(poolDerepRs,1,function(x){
          derepRFile = readRDS(x)
        })
        
        # Setting error model to be used
        if(orientFR.split){
          is.orientR = grepl(".orientR", poolSampleNames, fixed = TRUE)
          if(!(all(is.orientR) | all(!is.orientR))){
            stop("If orientFR.split and pooling, pools should contain only samples with same orientation")
          }
          is.orientR = all(is.orientR)
          
          if(!is.orientR){
            dadaErrR = orientFerrR
          }else{
            dadaErrR = orientRerrR
          }
        }else{
          dadaErrR = errR
        }
        
        # Setting priors to be used
        if(!isEmpty(priorsR)){
          dadapriorsR = priorsR[[dadaPoolName]]
          if(is.null(dadapriorsR)){
            dadapriorsR = character(0)
          }
        }else{
          dadapriorsR = character(0)
        }
        
        #DADA denoising
        dadaRList <- dada(poolDerepRFileList, err=dadaErrR, priors = dadapriorsR, pool=pool, OMEGA_A = OMEGA_A, multithread=TRUE)
        
        for(i in 1:length(dadaRList)){
          saveRDS(dadaRList[[i]],file = paste0(path,"/dada denoised/",poolSampleNames[i],"_dadaR.rds"))
          dadaRs[poolSampleNamesMatch[i],1] <<- paste0(path,"/dada denoised/",poolSampleNames[i],"_dadaR.rds")
          dadaRs[poolSampleNamesMatch[i],2] <<- sum(getUniques(dadaRList[[i]]))
        }
        
        track[,"denoisedR"] = as.numeric(dadaRs[,2])
        cat("\n")
        print(track[poolSampleNamesMatch,])
      })
      
      cat("\n")
      print(track[inPoolListNamesMatch,])
    }else{
      outPoolListNamesMatch = 1:length(sample.names)
    }
    
    # If any samples are non-pooled, denoise them
    if(length(outPoolListNamesMatch) > 0){
      if(!isEmpty(priorsF)){
        priorMatchDf = adply(derepFs[outPoolListNamesMatch],1, .id = NULL, function(x){
          sampleName = laply(strsplit(x, "dereplicated/", fixed = TRUE), function(y) y[2])
          sampleName = laply(strsplit(sampleName, "_derepF.rds", fixed = TRUE), function(y) y[1])
          
          derepFFile = readRDS(x)
          
          if(printPriorHead){
            # printing some dereplicated fwd reads and fwd priors to help troubleshoot trimming priors
            cat("\n")
            cat("head of dereplicated Fwd sequences:")
            print(head(names(derepFFile$uniques)))
            cat("head of prior Fwd sequences:")
            print(head(priorsF[[sampleName]]))
            
            # matching fwd priors and fwd dereplicated reads for troubleshooting/sanity check
            cat("\n")
            cat(paste0(sum(names(derepFFile$uniques) %in% priorsF[[sampleName]]), " Fwd dereplicated priors matches out of ", length(priorsF[[sampleName]]), " Fwd priors\n"))
            cat(paste0(sum(derepFFile$uniques[names(derepFFile$uniques) %in% priorsF[[sampleName]]]), " priors matches found out of ", sum(derepFFile$uniques), " Fwd sequences\n"))
          }
          
          derepInPriorsCount = sum(names(derepFFile$uniques) %in% priorsF[[sampleName]])
          priorsCount = length(priorsF[[sampleName]])
          derepInPriorsReadCount = sum(derepFFile$uniques[names(derepFFile$uniques) %in% priorsF[[sampleName]]])
          readCount = sum(derepFFile$uniques)
          
          # same data as above in tabular format
          df = data.frame(sampleName = sampleName,
                          derepInPriorsCount = derepInPriorsCount,
                          priorsCount = priorsCount,
                          derepInPriorsReadCount = derepInPriorsReadCount,
                          readCount = readCount)
          
          return(df)
        },.progress = "text")
        print(priorMatchDf)
      }else{
        cat("no fwd priors\n")
      }
      
      #DADA denoising
      dadaFs[outPoolListNamesMatch,] = aaply(derepFs[outPoolListNamesMatch],1, .drop = FALSE,function(x){
        sampleName = laply(strsplit(x, "dereplicated/", fixed = TRUE), function(y) y[2])
        sampleName = laply(strsplit(sampleName, "_derepF.rds", fixed = TRUE), function(y) y[1])
        
        derepFFile = readRDS(x)
        
        # Setting error model to be used
        if(orientFR.split){
          is.orientR = grepl(".orientR", sampleName, fixed = TRUE)
          if(!is.orientR){
            dadaErrF = orientFerrF
          }else{
            dadaErrF = orientRerrF
          }
        }else{
          dadaErrF = errF
        }
        
        # Setting priors to be used
        if(!isEmpty(priorsF)){
          dadaPriorsF = priorsF[[sampleName]]
          if(is.null(dadaPriorsF)){
            dadaPriorsF = character(0)
          }
        }else{
          dadaPriorsF = character(0)
        }
        
        cat("\n")
        dadaF <- dada(derepFFile, err=dadaErrF, priors = dadaPriorsF, pool=FALSE, OMEGA_A = OMEGA_A, multithread=TRUE)
        saveRDS(dadaF,file = paste0(path,"/dada denoised/",sampleName,"_dadaF.rds"))
        
        name = paste0(path,"/dada denoised/",sampleName,"_dadaF.rds")
        count = sum(getUniques(dadaF))
        return(c(name,count))
      },.progress = "text")
      
      track[,"denoisedF"] = as.numeric(dadaFs[,2])
      print(track[outPoolListNamesMatch,])
      
      if(!isEmpty(priorsR)){
        priorMatchDf = adply(derepRs,1, .id = NULL, function(x){
          sampleName = laply(strsplit(x, "dereplicated/", fixed = TRUE), function(y) y[2])
          sampleName = laply(strsplit(sampleName, "_derepR.rds", fixed = TRUE), function(y) y[1])
          
          derepRFile = readRDS(x)
          
          if(printPriorHead){
            # printing some dereplicated rev reads and rev priors to help troubleshoot trimming priors
            cat("\n")
            cat("head of dereplicated Rev sequences:\n")
            print(head(names(derepRFile$uniques)))
            cat("head of prior Rev sequences:\n")
            print(head(priorsR[[sampleName]]))
            
            # matching rev priors and rev dereplicated reads for troubleshooting/sanity check
            cat(" ")
            cat(paste0(sum(names(derepRFile$uniques) %in% priorsR[[sampleName]]), " Rev dereplicated priors matches out of ", length(priorsR[[sampleName]]), " Rev priors\n"))
            cat(paste0(sum(derepRFile$uniques[names(derepRFile$uniques) %in% priorsR[[sampleName]]]), " priors matches found out of ", sum(derepRFile$uniques), " Rev sequences\n"))
          }
          
          derepInPriorsCount = sum(names(derepRFile$uniques) %in% priorsR[[sampleName]])
          priorsCount = length(priorsR[[sampleName]])
          derepInPriorsReadCount = sum(derepRFile$uniques[names(derepRFile$uniques) %in% priorsR[[sampleName]]])
          readCount = sum(derepRFile$uniques)
          
          # same data as above in tabular format
          df = data.frame(sampleName = sampleName,
                          derepInPriorsCount = derepInPriorsCount,
                          priorsCount = priorsCount,
                          derepInPriorsReadCount = derepInPriorsReadCount,
                          readCount = readCount)
          
          return(df)
        },.progress = "text")
        print(priorMatchDf)
      }else{
        cat("no rev priors\n")
      }
      
      #DADA denoising
      dadaRs[outPoolListNamesMatch,] = aaply(derepRs[outPoolListNamesMatch],1, .drop = FALSE,function(x){
        sampleName = laply(strsplit(x, "dereplicated/", fixed = TRUE), function(y) y[2])
        sampleName = laply(strsplit(sampleName, "_derepR.rds", fixed = TRUE), function(y) y[1])
        
        derepRFile = readRDS(x)
        
        # Setting error model to be used
        if(orientFR.split){
          is.orientR = grepl(".orientR", sampleName, fixed = TRUE)
          if(!is.orientR){
            dadaErrR = orientFerrR
          }else{
            dadaErrR = orientRerrR
          }
        }else{
          dadaErrR = errR
        }
        
        # Setting priors to be used
        if(!isEmpty(priorsR)){
          dadaPriorsR = priorsR[[sampleName]]
          if(is.null(dadaPriorsR)){
            dadaPriorsR = character(0)
          }
        }else{
          dadaPriorsR = character(0)
        }
        
        cat("\n")
        dadaR <- dada(derepRFile, err=dadaErrR, priors = dadaPriorsR, pool=FALSE, OMEGA_A = OMEGA_A, multithread=TRUE)
        saveRDS(dadaR,file = paste0(path,"/dada denoised/",sampleName,"_dadaR.rds"))
        
        name = paste0(path,"/dada denoised/",sampleName,"_dadaR.rds")
        count = sum(getUniques(dadaR))
        return(c(name,count))
      },.progress = "text")
      
      track[,"denoisedR"] = as.numeric(dadaRs[,2])
      print(track[outPoolListNamesMatch,])
    }
    print(track)
    
  }else{
    dadaFs = data.frame(name = paste0(path,"/dada denoised/",sample.names,"_dadaF.rds"))
    dadaRs = data.frame(name = paste0(path,"/dada denoised/",sample.names,"_dadaR.rds"))
  }
  
  # ---- merging denoised F and R files ----
  cat("\n")
  cat("Merging F and R sequences...\n")
  dir.create(path = paste0(path,"/outputs"), showWarnings = TRUE)
  
  mergerArgs = data.frame(dadaFs = dadaFs[,1], derepFs, dadaRs = dadaRs[,1], derepRs, stringsAsFactors = FALSE)
  
  # merging fwd and rev sequences is done multicore parallel, if possible
  # generating parallel backend
  nodes <- detectCores()
  cl <- makeCluster(nodes, type = "PSOCK")
  registerDoParallel(cl)
  
  #merging fwd and rev sequences
  mergers = mlply(mergerArgs, function(dadaFs, derepFs, dadaRs, derepRs){
    dadaF = readRDS(dadaFs)
    derepF = readRDS(derepFs)
    dadaR = readRDS(dadaRs)
    derepR = readRDS(derepRs)
    
    merger <- mergePairs(dadaF, derepF, dadaR, derepR, maxMismatch = 0, verbose=TRUE)
    return(merger)
  }, .parallel = TRUE, .paropts = list(.packages = "dada2"))
  
  #stopping parallel backend
  stopCluster(cl)
  
  names(mergers) = sample.names
  
  # Remove dereplicated files if saveDereplicated == FALSE and the process() generated dereplicated files (preDerepToggle == FALSE)
  if(!saveDereplicated & !preDerepToggle){
    file.remove(derepFs,derepRs)
  }
  
  # Remove denoised files if saveDenoised == FALSE and the process() generated denoised files (preDenoiseToggle == FALSE)
  if(!saveDenoised & !preDenoiseToggle){
    file.remove(dadaFs[,1],dadaRs[,1])
  }
  
  #saving mergers
  for(i in 1:length(sample.names)){
    saveRDS(mergers[[i]],file = paste0(path,"/outputs/",sample.names[i],"_mergers.rds"))
  }
  
  track[,"merged"] = sapply(mergers, function(x) sum(getUniques(x)))
  print(track)
  
  # saving sequence tables
  for(i in 1:length(sample.names)){
    seqtabIndiv <- makeSequenceTable(mergers[[i]])
    saveRDS(seqtabIndiv,file = paste0(path,"/outputs/",sample.names[i],"_seqtab.rds"))
  }
  
  cat("Removing bimeras...\n")
  if(orientFR.split){
    seqtab.orientF <- makeSequenceTable(mergers[grep(".orientF", sample.names, fixed = TRUE)])
    saveRDS(seqtab.orientF,file = paste0(path,"/outputs/",errPoolName,".orientF_seqtab.rds"))
    
    seqtab.orientR <- makeSequenceTable(mergers[grep(".orientR", sample.names, fixed = TRUE)])
    saveRDS(seqtab.orientR,file = paste0(path,"/outputs/",errPoolName,".orientR_seqtab.rds"))
    
    #consensus-based chimera removal (sequences which appear to be composed of two parent sequences)
    seqtab.orientF.nochim <- removeBimeraDenovo(seqtab.orientF, method="consensus", multithread=TRUE, verbose=TRUE)
    seqtab.orientR.nochim <- removeBimeraDenovo(seqtab.orientR, method="consensus", multithread=TRUE, verbose=TRUE)
    saveRDS(seqtab.orientF.nochim,file = paste0(path,"/outputs/",errPoolName,".orientF_seqtab.nochim.rds"))
    saveRDS(seqtab.orientR.nochim,file = paste0(path,"/outputs/",errPoolName,".orientR_seqtab.nochim.rds"))
    
    rownames(seqtab.orientF.nochim) = gsub(".orientF","",rownames(seqtab.orientF.nochim))
    rownames(seqtab.orientR.nochim) = gsub(".orientR","",rownames(seqtab.orientR.nochim))
    colnames(seqtab.orientR.nochim) = as.character(reverseComplement(DNAStringSet(colnames(seqtab.orientR.nochim))))
    
    seqtab.nochim <- mergeSequenceTables(tables = list(seqtab.orientF.nochim,seqtab.orientR.nochim), repeats = "sum")
    seqtab.nochim <- removeBimeraDenovo(seqtab.nochim, method="consensus", multithread=TRUE, verbose=TRUE)
    saveRDS(seqtab.nochim,file = paste0(path,"/outputs/",errPoolName,"_seqtab.nochim.rds"))
    
    track[grep(".orientF", sample.names, fixed = TRUE),"nonchim"] = rowSums(seqtab.orientF.nochim)
    track[grep(".orientR", sample.names, fixed = TRUE),"nonchim"] = rowSums(seqtab.orientR.nochim)
  }else{
    seqtab <- makeSequenceTable(mergers)
    saveRDS(seqtab,file = paste0(path,"/outputs/",errPoolName,"_seqtab.rds"))
    
    #consensus-based chimera removal (sequences which appear to be composed of two parent sequences)
    seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
    saveRDS(seqtab.nochim,file = paste0(path,"/outputs/",errPoolName,"_seqtab.nochim.rds"))
    
    track[,"nonchim"] = rowSums(seqtab.nochim)
  }
  
  cat("\n")
  cat("Done.\n")
  return(track)
}


# ==== getQualityProfile function ====
getQualityProfile = function (fl, n = 5e+05, aggregate = FALSE) 
{
  statdf <- data.frame(Cycle = integer(0), Mean = numeric(0), 
                       Q25 = numeric(0), Q50 = numeric(0), Q75 = numeric(0), 
                       Cum = numeric(0), file = character(0))
  anndf <- data.frame(minScore = numeric(0), label = character(0), 
                      rclabel = character(0), rc = numeric(0), file = character(0))
  FIRST <- TRUE
  for (f in fl[!is.na(fl)]) {
    srqa <- ShortRead::qa(f, n = n)
    df <- srqa[["perCycle"]]$quality
    rc <- sum(srqa[["readCounts"]]$read)
    if (rc >= n) {
      rclabel <- paste("Reads >= ", n)
    }
    else {
      rclabel <- paste("Reads: ", rc)
    }
    means <- rowsum(df$Score * df$Count, df$Cycle)/rowsum(df$Count, 
                                                          df$Cycle)
    get_quant <- function(xx, yy, q) {
      xx[which(cumsum(yy)/sum(yy) >= q)][[1]]
    }
    q25s <- by(df, df$Cycle, function(foo) get_quant(foo$Score, 
                                                     foo$Count, 0.25), simplify = TRUE)
    q50s <- by(df, df$Cycle, function(foo) get_quant(foo$Score, 
                                                     foo$Count, 0.5), simplify = TRUE)
    q75s <- by(df, df$Cycle, function(foo) get_quant(foo$Score, 
                                                     foo$Count, 0.75), simplify = TRUE)
    cums <- by(df, df$Cycle, function(foo) sum(foo$Count), 
               simplify = TRUE)
    if (!all(sapply(list(names(q25s), names(q50s), names(q75s), 
                         names(cums)), identical, rownames(means)))) {
      stop("Calculated quantiles/means weren't compatible.")
    }
    if (FIRST) {
      plotdf <- cbind(df, file = basename(f))
      FIRST <- FALSE
    }
    else {
      plotdf <- rbind(plotdf, cbind(df, file = basename(f)))
    }
    statdf <- rbind(statdf, data.frame(Cycle = as.integer(rownames(means)), 
                                       Mean = means, Q25 = as.vector(q25s), Q50 = as.vector(q50s), 
                                       Q75 = as.vector(q75s), Cum = 10 * as.vector(cums)/rc, 
                                       file = basename(f)))
    anndf <- rbind(anndf, data.frame(minScore = min(df$Score), 
                                     label = basename(f), rclabel = rclabel, rc = rc, 
                                     file = basename(f)))
  }
  anndf$minScore <- min(anndf$minScore)
  if (aggregate) {
    plotdf.summary <- aggregate(Count ~ Cycle + Score, plotdf, 
                                sum)
    plotdf.summary$label <- paste(nrow(anndf), "files (aggregated)")
    means <- rowsum(plotdf.summary$Score * plotdf.summary$Count, 
                    plotdf.summary$Cycle)/rowsum(plotdf.summary$Count, 
                                                 plotdf.summary$Cycle)
    q25s <- by(plotdf.summary, plotdf.summary$Cycle, function(foo) get_quant(foo$Score, 
                                                                             foo$Count, 0.25), simplify = TRUE)
    q50s <- by(plotdf.summary, plotdf.summary$Cycle, function(foo) get_quant(foo$Score, 
                                                                             foo$Count, 0.5), simplify = TRUE)
    q75s <- by(plotdf.summary, plotdf.summary$Cycle, function(foo) get_quant(foo$Score, 
                                                                             foo$Count, 0.75), simplify = TRUE)
    cums <- by(plotdf.summary, plotdf.summary$Cycle, function(foo) sum(foo$Count), 
               simplify = TRUE)
    statdf.summary <- data.frame(Cycle = as.integer(rownames(means)), 
                                 Mean = means, Q25 = as.vector(q25s), Q50 = as.vector(q50s), 
                                 Q75 = as.vector(q75s), Cum = 10 * as.vector(cums)/rc)
    p = statdf.summary
  }
  else {
    p = statdf
  }
  p
}

# ==== plotQualityProfile2 function ====
# test = getQualityProfile(fnFs,aggregate = TRUE)
# test2 = getQualityProfile(fnRs,aggregate = TRUE)
# test3 = data.frame(CycleF = test$Cycle, MeanF = test$Mean)
# test3 = test3 %>%
#   left_join(data.frame(CycleF = 254-test2$Cycle, MeanR = test2$Mean)) %>%
#   mutate(diff = abs(MeanF - MeanR))
# test3$cSumDiff[which(test3$diff==min(test3$diff[!is.na(test3$diff)])):tail(which(!is.na(test3$diff)),n = 1)] = cumsum(test3$diff[which(test3$diff==min(test3$diff[!is.na(test3$diff)])):tail(which(!is.na(test3$diff)),n = 1)])
# test3$cSumDiff[which(test3$diff==min(test3$diff[!is.na(test3$diff)])):head(which(!is.na(test3$diff)),n = 1)] = cumsum(test3$diff[which(test3$diff==min(test3$diff[!is.na(test3$diff)])):head(which(!is.na(test3$diff)),n = 1)])
# test3$div = test3$cSumDiff/(abs(test3$CycleF-which(test3$diff==min(test3$diff[!is.na(test3$diff)])))+1)

# ggplot()+
#   geom_line(data = test, aes(x = Cycle, y = Mean), color = "red")+
#   geom_line(data = test2, aes(x = 254-Cycle, y = Mean), color = "blue")+
#   geom_vline(xintercept = trimLeftSelect[1], color = "black") + 
#   geom_vline(xintercept = truncLenSelect[1], color = "red") +
#   geom_vline(xintercept = 254-trimLeftSelect[2], color = "black") + 
#   geom_vline(xintercept = 254-truncLenSelect[2], color = "blue")+
#   geom_line(data = test3, aes(x = CycleF, y = diff), color = "black")+
#   geom_line(data = test3, aes(x = CycleF, y = div), color = "red")

plotQualityProfile2 = function(fnFs,fnRs,primerLenF,primerLenR,seqlen,trimLeftSelect,truncLenSelect){
  qpF = getQualityProfile(fnFs,aggregate = TRUE)
  qpR = getQualityProfile(fnRs,aggregate = TRUE)
  qp = data.frame(CycleF = qpF$Cycle, MeanF = qpF$Mean)
  qp = qp %>%
    left_join(data.frame(CycleF = seqlen-qpR$Cycle, MeanR = qpR$Mean)) %>%
    mutate(diff = abs(MeanF - MeanR))
  qp$cSumDiff = NA
  qp$cSumDiff[which(qp$diff==min(qp$diff[!is.na(qp$diff)])):tail(which(!is.na(qp$diff)),n = 1)] = cumsum(qp$diff[which(qp$diff==min(qp$diff[!is.na(qp$diff)])):tail(which(!is.na(qp$diff)),n = 1)])
  qp$cSumDiff[which(qp$diff==min(qp$diff[!is.na(qp$diff)])):head(which(!is.na(qp$diff)),n = 1)] = cumsum(qp$diff[which(qp$diff==min(qp$diff[!is.na(qp$diff)])):head(which(!is.na(qp$diff)),n = 1)])
  qp$div = qp$cSumDiff/(abs(qp$CycleF-which(qp$diff==min(qp$diff[!is.na(qp$diff)])))+1)
  
  # QualityProfileF = QualityProfileF + 
  #   geom_vline(xintercept = trimLeftSelect[1]) + 
  #   geom_vline(xintercept = truncLenSelect[1]) + 
  #   geom_vline(xintercept = (seqlen + primerLenF + primerLenR) - (trimLeftSelect[2] + truncLenSelect[2]), color = "red")
  # QualityProfileR = QualityProfileR + 
  #   geom_vline(xintercept = trimLeftSelect[2]) + 
  #   geom_vline(xintercept = truncLenSelect[2]) + 
  #   geom_vline(xintercept = (seqlen + primerLenF + primerLenR) - (trimLeftSelect[1] + truncLenSelect[1]), color = "red")
  
  plot = ggplot()+
    geom_line(data = qpF, aes(x = Cycle, y = Mean), color = "red")+
    geom_line(data = qpR, aes(x = 254-Cycle, y = Mean), color = "blue")+
    geom_vline(xintercept = trimLeftSelect[1], color = "black") + 
    geom_vline(xintercept = truncLenSelect[1], color = "red") +
    geom_vline(xintercept = seqlen-trimLeftSelect[2], color = "black") + 
    geom_vline(xintercept = seqlen-truncLenSelect[2], color = "blue")+
    geom_line(data = qp, aes(x = CycleF, y = diff), color = "black")+
    geom_line(data = qp, aes(x = CycleF, y = div), color = "red")
  
  print(plot)
  
  print(paste0("Fwd and Rev read overlap is ",truncLenSelect[1] - ((seqlen + primerLenF + primerLenR) - (truncLenSelect[2]))," bases"))
}


# ==== seq2fasta conversion function ====
seq2fasta = function(seqNames, sequences, fileName, savePath){
  fasta = paste0(">",seqNames,"\n",sequences,"\n")
  dir.create(path = savePath, showWarnings = TRUE)
  write(fasta,file = paste0(savePath,"/",fileName,".txt"), append = FALSE)
}
# ==== BLAST_16S search function ====
BLAST_16S = function(fastafileName,
                     RID = NULL,
                     
                     remote = FALSE, 
                     blast_path = NULL,
                     db_16S = "/Users/lab/Desktop/Christopher_Yau/BLAST_databases/16S_ribosomal_RNA/16S_ribosomal_RNA", 
                     parallel = TRUE,
                     
                     wait = FALSE,
                     remoteWaitTime = 10,
                     
                     RID_out = TRUE,
                     RID_outDirName = if(is.null(fastafileName)){
                       NULL
                     }else{
                       dirname(fastafileName)
                     },
                     outFile = NULL,
                     
                     topHit = FALSE,
                     taxonomy = FALSE){
  # ---- Retrieving BLAST result(s) by RIDs ----
  if(!is.null(RID)){
    print(paste0("Retrieving BLAST result(s) RID(s) ",paste0(RID, collapse = ", "),"..."))
    BLASTsearchTable = data.frame("RID" = RID)
    
    if(is.null(outFile)){
      if(is.null(RID_outDirName)){
        stop("RID_outDirName or outFile must be provided.")
      }
      outFile = paste0(RID_outDirName,"/",paste0(RID, collapse = "_"),"_BLAST_out.txt")
    }
  }
  # ---- Performing BLAST search ----
  else{
    # setting output file
    if(is.null(outFile)){
      outFile = paste0(fastafileName,"_BLAST_out.txt")
    }
    if(!remote){
      # ---- Local BLAST search ----
      # setting parallel parameter
      if(parallel){
        ncpu = detectCores()
      }else{
        ncpu = 1
      }
      
      print("16S BLAST will be run locally; output will be written to:")
      print(outFile)
      
      # setting parameters for local BLAST query
      if(!taxonomy){
        outfmt = "6 qseqid saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore"
      }else{
        outfmt = "6 qseqid saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids"
      }
      
      blast_args = c("-db", db_16S,
                     "-query", shQuote(fastafileName),
                     "-num_threads", ncpu,
                     "-outfmt", shQuote(outfmt)
      )
      
      if(is.null(blast_path)){
        blast_path = Sys.which("blastn")
        blast_path = blast_path[which(blast_path != "")[1]]
      }
      
      # BLAST commmand
      cmd = system2(blast_path, args = blast_args, stdout = outFile, 
                    stderr = FALSE, wait = wait)
      
      if(!wait){
        # if not waiting, run BLAST command, resulting in outFile with no header
        return(cmd)
      }
      # ---- Wait for local BLAST output ----
      else{
        # if waiting, run BLAST command, wait for outFile, overwrite outFile with outFile with header
        cmd
        blastResults = read.delim(outFile, header =  FALSE, sep = "\t",stringsAsFactors =  FALSE)
        if(!taxonomy){
          colnames(blastResults) = c("query", "subject","per_id", "alignLen", "mismatches", "gapOpens", "q_start", "q_end", "s_start", "s_end", "evalue", "bitScore")
        }else{
          colnames(blastResults) = c("query", "subject","per_id", "alignLen", "mismatches", "gapOpens", "q_start", "q_end", "s_start", "s_end", "evalue", "bitScore", "taxid")
        }
        
        write.table(blastResults,
                    file = outFile, 
                    quote = FALSE,
                    sep = "\t",
                    row.names = FALSE)
        
        # If topHits, filter for top hits
        if(topHit){
          blastResults = blastResults[order(blastResults$evalue),]
          blastResults = blastResults[!duplicated(blastResults$query),]
          blastResults = blastResults[order(blastResults$query),]
          
          # if taxonomy, query NCBI Taxonomy for taxonomic information
          if(taxonomy){
            taxidtermList = split(unique(blastResults$taxid), rep(1:ceiling(length(unique(blastResults$taxid))/500), each = 500, length.out = length(unique(blastResults$taxid))))
            
            lineageList = llply(taxidtermList, function(x){
              taxonomy = classification(x, db = 'ncbi', callopts=list(http_version = 0L))
              return(taxonomy)
            })
            
            lineage = ldply(lineageList, .id = NULL,function(x){
              df = ldply(x, .id = NULL, function(y){
                df2 = y$name[match(c("superkingdom","phylum","class","order","family","genus","species","strain"),y$rank)]
                if(is.na(df2[8])){
                  df2[8] = df2[7]
                }
                df2[which(is.na(df2))] = paste0(df2[which(is.na(df2))-1],"_unclassified")
                df2 = data.frame(t(df2))
                colnames(df2) = c("BLASTsuperkingdom","BLASTphylum","BLASTclass","BLASTorder","BLASTfamily","BLASTgenus","BLASTspecies", "BLASTmatch")
                return(df2)
              })
              return(df)
            })
            
            blastResults = data.frame(blastResults, 
                                      lineage[match(blastResults$taxid,unique(blastResults$taxid)),]
            )
          }
        }
        return(blastResults)
      }
    }else{
      # ---- Remote BLAST search ----
      # determine size of the FASTA file, split into chunks <900000 bytes in size if needed
      if(file.size(fastafileName)>1000000){
        print("FASTA file exceeds max size for web-BLAST; chunking FASTA file...")
        fastafile = readFasta(fastafileName)
        fastafileChunks = split(fastafile, rep(1:ceiling(file.size(fastafileName)/900000),
                                               each = ceiling(length(fastafile)/ceiling(file.size(fastafileName)/900000)),
                                               length.out = length(fastafile)))
      }else{
        fastafileChunks = list(readFasta(fastafileName))
      }
      
      # send BLAST searches, receive RIDs and ETAs
      BLASTsearchTable = ldply(fastafileChunks, .id = NULL, function(x){
        xChar = paste0(">",as.character(ShortRead::id(x)),"\n",as.character(sread(x)),"\n", collapse = "")
        
        BLASTsearch = POST("https://blast.ncbi.nlm.nih.gov/Blast.cgi", body = list(
          PROGRAM = "blastn",
          MEGABLAST= "on",
          QUERY = xChar,
          DATABASE = "rRNA_typestrains/16S_ribosomal_RNA",
          CMD = "Put")
        )
        
        BLASTsearchResponse = httr::content(BLASTsearch, as = "text")
        RID = regexec("RID = (.*?)\n",BLASTsearchResponse)
        RID = substr(BLASTsearchResponse,start = RID[[1]][2], stop = RID[[1]][2]+ attr(RID[[1]],"match.length")[2] -1)
        print(paste0("BLAST search ID is ",RID))
        
        ETA = regexec("RTOE = (.*?)\n",BLASTsearchResponse)
        ETA = substr(BLASTsearchResponse,start = ETA[[1]][2], stop = ETA[[1]][2]+ attr(ETA[[1]],"match.length")[2] -1)
        ETA = as.numeric(ETA)
        
        print(paste0("Estimated BLAST run time is ",ETA," seconds..."))
        print(paste0("Submitting next NCBI BLAST query in ", remoteWaitTime ," seconds..."))
        
        nextQuerypb <- txtProgressBar(min = 0, max = as.numeric(remoteWaitTime), style = 3)
        for(i in 1:remoteWaitTime){
          Sys.sleep(1)
          # update progress bar
          setTxtProgressBar(nextQuerypb, i)
        }
        close(nextQuerypb)
        
        return(data.frame("RID" = RID, "ETA" = ETA))
      })
      
      # if not waiting, return RIDs and ETAs
      if(!wait){
        return(BLASTsearchTable)
      }
      
      # ---- Wait for max ETA to elapse ----
      maxETA = max(BLASTsearchTable$ETA)
      
      print(paste0("Max estimated BLAST run time is ",maxETA," seconds..."))
      
      maxETApb <- txtProgressBar(min = 0, max = maxETA, style = 3)
      for(i in 1:maxETA){
        Sys.sleep(1)
        # update progress bar
        setTxtProgressBar(maxETApb, i)
      }
      close(maxETApb)
      print("Max estimated BLAST run time elapsed...")
    }
  }
  # ---- Check remote BLAST search status(es) ----
  BLASTStatus = "WAITING"
  
  while(any(BLASTStatus=="WAITING")){
    print(paste0("Querying NCBI BLAST for results in ", remoteWaitTime ," seconds..."))
    BLASTStatuspb <- txtProgressBar(min = 0, max = as.numeric(remoteWaitTime), style = 3)
    for(i in 1:remoteWaitTime){
      Sys.sleep(1)
      # update progress bar
      setTxtProgressBar(BLASTStatuspb, i)
    }
    close(BLASTStatuspb)
    
    remoteWaitTime = remoteWaitTime*2
    
    BLASTStatus = aaply(BLASTsearchTable$RID, .margins = 1, function(x){
      BLASTSearchInfo = POST("https://blast.ncbi.nlm.nih.gov/Blast.cgi",body = list(
        FORMAT_OBJECT="SearchInfo",
        RID=x,
        CMD="Get")
      )
      
      BLASTSearchInfoResponse = httr::content(BLASTSearchInfo, as = "text")
      BLASTStatus = regexec("Status=(.*?)\n",BLASTSearchInfoResponse)
      BLASTStatus = substr(BLASTSearchInfoResponse,start = BLASTStatus[[1]][2], stop = BLASTStatus[[1]][2]+ attr(BLASTStatus[[1]],"match.length")[2] -1)
      
      return(BLASTStatus)
    })
    print(data.frame(RID = BLASTsearchTable$RID, status = BLASTStatus))
    if(!wait){
      break
    }
  }
  if(!wait & any(BLASTStatus=="WAITING")){
    print("BLAST search still waiting...")
    return()
  }
  
  # ---- Retrieve remote BLAST search results ----
  dirName = dirname(outFile)
  print("BLAST search ready.")
  print(data.frame(RID = BLASTsearchTable$RID, status = BLASTStatus))
  
  blast_formatter_path = Sys.which("blast_formatter")
  blast_formatter_path = blast_formatter_path[which(blast_formatter_path != "")[1]]
  
  # setting parameters for blast_formatter
  if(!taxonomy){
    outfmt = "6 qseqid saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore"
  }else{
    outfmt = "6 qseqid saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids"
  }
  
  blastResults = alply(BLASTsearchTable$RID[BLASTStatus == "READY"], .margins = 1,function(rid){
    blast_formatter_args = c("-rid", rid,
                             "-outfmt", shQuote(outfmt)
    )
    blastResult = system2(blast_formatter_path, args = blast_formatter_args, stdout = TRUE, 
                          stderr = FALSE, wait = TRUE)
    blastResult = paste0(blastResult, collapse = "\n")
    blastResult = read.delim(text = blastResult, header =  FALSE, sep = "\t",stringsAsFactors =  FALSE)
    
    if(!taxonomy){
      colnames(blastResult) = c("query", "subject","per_id", "alignLen", "mismatches", "gapOpens", "q_start", "q_end", "s_start", "s_end", "evalue", "bitScore")
    }else{
      colnames(blastResult) = c("query", "subject","per_id", "alignLen", "mismatches", "gapOpens", "q_start", "q_end", "s_start", "s_end", "evalue", "bitScore", "taxid")
    }
    
    if(RID_out){
      write.table(blastResult,
                  file = paste0(dirName,"/",rid,"-Alignment-HitTable",".txt"), 
                  quote = FALSE,sep = "\t",
                  row.names = FALSE)
    }
    return(blastResult)
  })
  
  blastResults = do.call(rbind, blastResults)
  
  write.table(blastResults,
              file = outFile, 
              quote = FALSE,sep = "\t",
              row.names = FALSE)
  
  # If topHits, filter for top hits
  if(topHit){
    blastResults = blastResults[order(blastResults$evalue),]
    blastResults = blastResults[!duplicated(blastResults$query),]
    blastResults = blastResults[order(blastResults$query),]
    
    # if taxonomy, query NCBI Taxonomy for taxonomic information
    if(taxonomy){
      taxidtermList = split(unique(blastResults$taxid), rep(1:ceiling(length(unique(blastResults$taxid))/500), each = 500, length.out = length(unique(blastResults$taxid))))
      
      lineageList = llply(taxidtermList, function(x){
        taxonomy = classification(x, db = 'ncbi', callopts=list(http_version = 0L))
        return(taxonomy)
      })
      
      lineage = ldply(lineageList, .id = NULL,function(x){
        df = ldply(x, .id = NULL, function(y){
          df2 = y$name[match(c("superkingdom","phylum","class","order","family","genus","species","strain"),y$rank)]
          if(is.na(df2[8])){
            df2[8] = df2[7]
          }
          df2[which(is.na(df2))] = paste0(df2[which(is.na(df2))-1],"_unclassified")
          df2 = data.frame(t(df2))
          colnames(df2) = c("BLASTsuperkingdom","BLASTphylum","BLASTclass","BLASTorder","BLASTfamily","BLASTgenus","BLASTspecies", "BLASTmatch")
          return(df2)
        })
        return(df)
      })
      
      blastResults = data.frame(blastResults, 
                                lineage[match(blastResults$taxid,unique(blastResults$taxid)),]
      )
    }
  }
  return(blastResults)
  
}

# ==== fecal seq. data derep object generation ====
path <- "./Alessandra_16s/CHILD 12 month samples DADA2"
seqlen = 253

seqMetadata = read.delim("./Alessandra_16s/CHILD 12 month samples DADA2/input/metadata/CHILD_12month_seq_metadata.txt", stringsAsFactors = FALSE)

# fnFs <- sort(list.files(file.path(path,"input/fastq"), pattern="_pass_1.fastq.gz", full.names = TRUE))
# fnRs <- sort(list.files(file.path(path,"input/fastq"), pattern="_pass_2.fastq.gz", full.names = TRUE))

# rename files to sample names and standard FASTQ naming format
# SRR_Fs = sapply(strsplit(basename(fnFs), "_pass_"), `[`, 1)
# newFnFNames = seqMetadata$Forward_read_filename[match(SRR_Fs,seqMetadata$Run)]
# file.path(path,"input/fastq",newFnFNames)
# 
# SRR_Rs = sapply(strsplit(basename(fnRs), "_pass_"), `[`, 1)
# newFnRNames = seqMetadata$Reverse_read_filename[match(SRR_Rs,seqMetadata$Run)]
# 
# file.rename(from = fnFs,to = file.path(path,"input/fastq",newFnFNames))
# file.rename(from = fnRs,to = file.path(path,"input/fastq",newFnRNames))

fnFs <- sort(list.files(file.path(path,"input/fastq"), pattern="R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(file.path(path,"input/fastq"), pattern="R2_001.fastq.gz", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Sequencing_run in seqMetadata is pools

primerLenF = 0
primerLenR = 0
trimLeft = c(0,0)
trimLeftSelect = c(primerLenF+trimLeft[1],primerLenR+trimLeft[2])
truncLenSelect = c(157,116)

priorsF = character(0)
priorsR = character(0)

flowcellNamesSub = data.frame(fnFs,fnRs,sample.names,flowcell = seqMetadata$Sequencing_run[match(sample.names,seqMetadata$Sample.Name)])
write.table(flowcellNamesSub, "./Alessandra_16s/CHILD 12 month samples DADA2/input/metadata/CHILD_12month_seq_flowcells.txt",quote = FALSE,sep = "\t",row.names = FALSE)

uniqueFlowcellNames = unique(flowcellNamesSub$flowcell)

fnFsTest = (flowcellNamesSub %>%
              group_by(flowcell) %>%
              slice_sample(n = 5))[c("flowcell","fnFs","fnRs")]


plotQualityProfile(fnFsTest$fnFs[1:10], aggregate = TRUE)


plotQualityProfileTest = plotQualityProfile3(fnFs = fnFsTest$fnFs,
                                             fnRs = fnFsTest$fnRs,
                                             primerLenF,primerLenR,seqlen,trimLeftSelect,truncLenSelect, lenLimF = 240, lenLimR = 240, overlapLen = 20)


poolList = dlply(flowcellNamesSub,.(flowcell), function(x){
  out = x$sample.names
})

process(path = path,
        flowcellNamesSub$fnFs, flowcellNamesSub$fnRs, flowcellNamesSub$sample.names, errPoolName = "CHILD_12M", orientFR.split = FALSE,
        trimLeftSelect = trimLeftSelect, truncLenSelect = truncLenSelect, filter.matchIDs = FALSE,
        ErrModelMonotonicity = FALSE,
        OMEGA_A = getDadaOpt(option = "OMEGA_A"), pool = "pseudo", poolList = poolList,
        priorsF = character(0), priorsR = character(0),
        
        plotToggle = FALSE,
        printPriorHead = FALSE,
        
        preFilteredToggle = TRUE,
        preCalcErrToggle = TRUE,
        preDerepToggle = TRUE,
        preDenoiseToggle = TRUE,
        
        saveFiltered = TRUE,
        saveDereplicated = TRUE,
        saveDenoised = TRUE)

# ==== POST-DADA2 ====
# ==== merge Sequence Tables ====
seqtabList = list.files(paste0(path,"/outputs"), pattern="^Run.*?_seqtab.nochim\\.rds", full.names = TRUE)
#merge Sequence Tables
seqtab.nochim = mergeSequenceTables(tables = seqtabList)

saveRDS(seqtab.nochim,paste0(path,"/outputs/seqtab.nochim.RDS"))

seqtab.nochim = readRDS(paste0(path,"/outputs/seqtab.nochim.RDS"))
seqtab.nochim = seqtab.nochim[rowSums(seqtab.nochim)>0,]

seqCounts = data.frame(sample = rownames(seqtab.nochim), counts = rowSums(seqtab.nochim))

ggplot()+
  geom_histogram(data = seqCounts,bins = 100, aes(x = counts))+
  scale_x_log10()

ggplot()+
  geom_histogram(data = seqCounts,bins = 100, aes(x = counts))+
  scale_x_log10()+
  facet_wrap(~flowcell)

# ==== generate ASV table ====
ASV_Table = data.frame(ASV_Num = seq_along(colnames(seqtab.nochim)), sequences = colnames(seqtab.nochim),t(seqtab.nochim), stringsAsFactors = FALSE)
row.names(ASV_Table) = NULL
ASV_Table$ASV_Num = sprintf("%05d", ASV_Table$ASV_Num)
ASV_Table$ASV_Num = paste0("ASV_", ASV_Table$ASV_Num)
#ASV_Table = read.delim(paste0(path,"/CHILD_12M ASV_Table.txt"),sep = "\t", stringsAsFactors = FALSE)
#write.table(ASV_Table,file = paste0(path,"/CHILD_12M ASV_Table.txt"),sep = "\t", row.names = FALSE)

# ---- BLAST against NCBI for species ----
seq2fasta(ASV_Table$ASV_Num, ASV_Table$sequences, "All_ASVs_isolates", file.path(path,"ASV FASTAs"))

blastResults = BLAST_16S(fastafileName = paste0(path,"/ASV FASTAs/All_ASVs_isolates.txt"),
                         remote = FALSE,
                         wait = TRUE,
                         topHit = TRUE,
                         taxonomy = TRUE)

ASV_Table = data.frame(ASV_Table[,1:2],
                       blastResults[match(ASV_Table$ASV_Num, blastResults$query),c("BLASTsuperkingdom","BLASTphylum","BLASTclass","BLASTorder","BLASTfamily","BLASTgenus","BLASTspecies","BLASTmatch")],
                       BLASTper_id = blastResults$per_id[match(ASV_Table$ASV_Num,blastResults$query)],
                       ASV_Table[,3:ncol(ASV_Table)])

write.table(ASV_Table,file = file.path(path,"CHILD_12M ASV_Table.txt"),sep = "\t", row.names = FALSE)
ASV_Table = read.delim(file.path(path,"CHILD_12M ASV_Table.txt"))
