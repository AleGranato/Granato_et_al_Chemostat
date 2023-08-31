library(BiocManager)
library(devtools)
# devtools::install_github("benjjneb/dada2") # change the ref argument to get other versions
library(dada2); packageVersion("dada2")
library(plyr)
library(reshape2)
library(ggplot2)
library(Biostrings)
library(doParallel)
library(ggalluvial)
# library(ggrepel)
# library(ggdendro)
# library(MASS)
library(rentrez)
# library(dendextend)
library(R.utils)
library(dplyr)
library(ggnewscale)
library(vegan)
library(taxize)
library(RGCCA)
library(biomformat)
library(ShortRead)
library(Rfast)
library(httr)
library(ComplexHeatmap)
library(scales)

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
  print("Filter and trimming...")
  dir.create(path = paste0(path,"/filtered/"), showWarnings = TRUE)
  
  # Set destination to filtered/ subdirectory
  filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
  filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
  
  if(!preFilteredToggle){
    #filter and trim function
    #NOTE: consider relaxing maxEE if needed
    out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = truncLenSelect, trimLeft = trimLeftSelect,
                         maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                         compress=TRUE, multithread=TRUE, matchIDs=filter.matchIDs) # On Windows set multithread=FALSE
    
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
  print("Generating error models...")
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
  print("dereplicating sequences...")
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
      sampleName = laply(strsplit(sampleName, "_", fixed = TRUE), function(y) y[1])
      
      derepF <- derepFastq(x, verbose=TRUE)
      saveRDS(derepF,file = paste0(filePath,"/dereplicated/",sampleName,"_derepF.rds"))
      return(paste0(filePath,"/dereplicated/",sampleName,"_derepF.rds"))
    }, filePath = path, .parallel = TRUE, .paropts = list(.packages = "dada2"))
    
    # dereplicating rev reads
    derepRs = aaply(filtRs,1,function(x, filePath = path){
      sampleName = laply(strsplit(x, "filtered/", fixed = TRUE), function(y) y[2])
      sampleName = laply(strsplit(sampleName, "_", fixed = TRUE), function(y) y[1])
      
      derepR <- derepFastq(x, verbose=TRUE)
      saveRDS(derepR,file = paste0(filePath,"/dereplicated/",sampleName,"_derepR.rds"))
      name = paste0(filePath,"/dereplicated/",sampleName,"_derepR.rds")
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
  print("denoising sequences...")
  dir.create(path = paste0(path,"/dada denoised"), showWarnings = TRUE)
  
  if(!preDenoiseToggle){
    if(pool == FALSE){
      if(!isEmpty(priorsF)){
        priorMatchDf = adply(derepFs,1, .id = NULL, function(x){
          sampleName = laply(strsplit(x, "dereplicated/", fixed = TRUE), function(y) y[2])
          sampleName = laply(strsplit(sampleName, "_", fixed = TRUE), function(y) y[1])
          
          derepFFile = readRDS(x)
          
          if(printPriorHead){
            # printing some dereplicated fwd reads and fwd priors to help troubleshoot trimming priors
            print(" ")
            print("head of dereplicated Fwd sequences:")
            print(head(names(derepFFile$uniques)))
            print("head of prior Fwd sequences:")
            print(head(priorsF[[sampleName]]))
            
            # matching fwd priors and fwd dereplicated reads for troubleshooting/sanity check
            print(" ")
            print(paste0(sum(names(derepFFile$uniques) %in% priorsF[[sampleName]]), " Fwd dereplicated priors matches out of ", length(priorsF[[sampleName]]), " Fwd priors"))
            print(paste0(sum(derepFFile$uniques[names(derepFFile$uniques) %in% priorsF[[sampleName]]]), " priors matches found out of ", sum(derepFFile$uniques), " Fwd sequences"))
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
        print("no fwd priors")
      }
      
      #DADA denoising
      dadaFs = aaply(derepFs,1, .drop = FALSE,function(x){
        sampleName = laply(strsplit(x, "dereplicated/", fixed = TRUE), function(y) y[2])
        sampleName = laply(strsplit(sampleName, "_", fixed = TRUE), function(y) y[1])
        
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
        }else{
          dadaPriorsF = character(0)
        }
        
        dadaF <- dada(derepFFile, err=dadaErrF, priors = dadaPriorsF, pool=FALSE, OMEGA_A = OMEGA_A, multithread=TRUE)
        saveRDS(dadaF,file = paste0(path,"/dada denoised/",sampleName,"_dadaF.rds"))
        
        name = paste0(path,"/dada denoised/",sampleName,"_dadaF.rds")
        count = sum(getUniques(dadaF))
        return(c(name,count))
      },.progress = "text")
      
      track[,"denoisedF"] = as.numeric(dadaFs[,2])
      print(track)
      
      if(!isEmpty(priorsR)){
        priorMatchDf = adply(derepRs,1, .id = NULL, function(x){
          sampleName = laply(strsplit(x, "dereplicated/", fixed = TRUE), function(y) y[2])
          sampleName = laply(strsplit(sampleName, "_", fixed = TRUE), function(y) y[1])
          
          derepRFile = readRDS(x)
          
          if(printPriorHead){
            # printing some dereplicated rev reads and rev priors to help troubleshoot trimming priors
            print(" ")
            print("head of dereplicated Rev sequences:")
            print(head(names(derepRFile$uniques)))
            print("head of prior Rev sequences:")
            print(head(priorsR[[sampleName]]))
            
            # matching rev priors and rev dereplicated reads for troubleshooting/sanity check
            print(" ")
            print(paste0(sum(names(derepRFile$uniques) %in% priorsR[[sampleName]]), " Rev dereplicated priors matches out of ", length(priorsR[[sampleName]]), " Rev priors"))
            print(paste0(sum(derepRFile$uniques[names(derepRFile$uniques) %in% priorsR[[sampleName]]]), " priors matches found out of ", sum(derepRFile$uniques), " Rev sequences"))
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
        print("no rev priors")
      }
      
      #DADA denoising
      dadaRs = aaply(derepRs,1, .drop = FALSE,function(x){
        sampleName = laply(strsplit(x, "dereplicated/", fixed = TRUE), function(y) y[2])
        sampleName = laply(strsplit(sampleName, "_", fixed = TRUE), function(y) y[1])
        
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
        }else{
          dadaPriorsR = character(0)
        }
        
        dadaR <- dada(derepRFile, err=dadaErrR, priors = dadaPriorsR, pool=FALSE, OMEGA_A = OMEGA_A, multithread=TRUE)
        saveRDS(dadaR,file = paste0(path,"/dada denoised/",sampleName,"_dadaR.rds"))
        
        name = paste0(path,"/dada denoised/",sampleName,"_dadaR.rds")
        count = sum(getUniques(dadaR))
        return(c(name,count))
      },.progress = "text")
      
      track[,"denoisedR"] = as.numeric(dadaRs[,2])
      print(track)
      
    }else{
      dadaFs = data.frame(name = rep(NA,length(fnFs)),
                          count = rep(NA,length(fnFs)))
      
      dadaRs = data.frame(name = rep(NA,length(fnRs)),
                          count = rep(NA,length(fnRs)))
      
      a_ply(names(poolList),1, function(dadaPoolName){
        print(paste0("denoising pool '",dadaPoolName,"'..."))
        
        poolSampleNames = poolList[[dadaPoolName]]
        poolSampleNamesMatch = match(poolSampleNames,sample.names)
        
        poolDerepFs = derepFs[poolSampleNamesMatch]
        poolDerepRs = derepRs[poolSampleNamesMatch]
        
        print(poolSampleNames)
        
        if(!isEmpty(priorsF)){
          priorMatchDf = adply(poolDerepFs,1, .id = NULL, function(x){
            sampleName = laply(strsplit(x, "dereplicated/", fixed = TRUE), function(y) y[2])
            sampleName = laply(strsplit(sampleName, "_", fixed = TRUE), function(y) y[1])
            
            derepFFile = readRDS(x)
            
            if(printPriorHead){
              # printing some dereplicated fwd reads and fwd priors to help troubleshoot trimming priors
              print(" ")
              print("head of dereplicated Fwd sequences:")
              print(head(names(derepFFile$uniques)))
              print("head of prior Fwd sequences:")
              print(head(priorsF[[dadaPoolName]]))
              
              # matching fwd priors and fwd dereplicated reads for troubleshooting/sanity check
              print(" ")
              print(paste0(sum(names(derepFFile$uniques) %in% priorsF[[dadaPoolName]]), " Fwd dereplicated priors matches out of ", length(priorsF[[dadaPoolName]]), " Fwd priors"))
              print(paste0(sum(derepFFile$uniques[names(derepFFile$uniques) %in% priorsF[[dadaPoolName]]]), " priors matches found out of ", sum(derepFFile$uniques), " Fwd sequences"))
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
          print("no fwd priors")
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
        }else{
          dadaPriorsF = character(0)
        }
        
        dadaFList <- dada(poolDerepFFileList, err=dadaErrF, priors = dadaPriorsF, pool=pool, OMEGA_A = OMEGA_A, multithread=TRUE)
        
        for(i in 1:length(dadaFList)){
          saveRDS(dadaFList[[i]],file = paste0(path,"/dada denoised/",poolSampleNames[i],"_dadaF.rds"))
          dadaFs[poolSampleNamesMatch[i],1] <<- paste0(path,"/dada denoised/",poolSampleNames[i],"_dadaF.rds")
          dadaFs[poolSampleNamesMatch[i],2] <<- sum(getUniques(dadaFList[[i]]))
        }
        
        track[,"denoisedF"] = as.numeric(dadaFs[,2])
        print(track)
        
        if(!isEmpty(priorsR)){
          priorMatchDf = adply(poolDerepRs,1, .id = NULL, function(x){
            sampleName = laply(strsplit(x, "dereplicated/", fixed = TRUE), function(y) y[2])
            sampleName = laply(strsplit(sampleName, "_", fixed = TRUE), function(y) y[1])
            
            derepRFile = readRDS(x)
            
            if(printPriorHead){
              # printing some dereplicated rev reads and rev priors to help troubleshoot trimming priors
              print(" ")
              print("head of dereplicated Rev sequences:")
              print(head(names(derepRFile$uniques)))
              print("head of prior Rev sequences:")
              print(head(priorsR[[dadaPoolName]]))
              
              # matching rev priors and rev dereplicated reads for troubleshooting/sanity check
              print(" ")
              print(paste0(sum(names(derepRFile$uniques) %in% priorsR[[dadaPoolName]]), " Rev dereplicated priors matches out of ", length(priorsR[[dadaPoolName]]), " Rev priors"))
              print(paste0(sum(derepRFile$uniques[names(derepRFile$uniques) %in% priorsR[[dadaPoolName]]]), " priors matches found out of ", sum(derepRFile$uniques), " Rev sequences"))
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
          print("no rev priors")
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
        }else{
          dadapriorsR = character(0)
        }
        
        dadaRList <- dada(poolDerepRFileList, err=dadaErrR, priors = dadapriorsR, pool=pool, OMEGA_A = OMEGA_A, multithread=TRUE)
        
        for(i in 1:length(dadaRList)){
          saveRDS(dadaRList[[i]],file = paste0(path,"/dada denoised/",poolSampleNames[i],"_dadaR.rds"))
          dadaRs[poolSampleNamesMatch[i],1] <<- paste0(path,"/dada denoised/",poolSampleNames[i],"_dadaR.rds")
          dadaRs[poolSampleNamesMatch[i],2] <<- sum(getUniques(dadaRList[[i]]))
        }
        
        track[,"denoisedR"] = as.numeric(dadaRs[,2])
        print(track)
      })
    }
  }else{
    dadaFs = data.frame(name = paste0(path,"/dada denoised/",sample.names,"_dadaF.rds"))
    dadaRs = data.frame(name = paste0(path,"/dada denoised/",sample.names,"_dadaR.rds"))
  }
  
  # ---- merging denoised F and R files ----
  print("merging F and R sequences...")
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
  
  return(track)
}



# ==== getAggregateQualityScores function ====
getAggregateQualityScores = function (fl, n = 5e+05){
  plotdf = lapply(fl[!is.na(fl)], function(x){
    srqa <- ShortRead::qa(x, n = n)
    df <- srqa[["perCycle"]]$quality
  })
  plotdf = do.call(rbind, plotdf)
  
  plotdf.summary <- aggregate(Count ~ Cycle + Score, plotdf,
                              sum)
  means <- rowsum((1-10^(plotdf.summary$Score/-10)) * plotdf.summary$Count, 
                  plotdf.summary$Cycle)/rowsum(plotdf.summary$Count, 
                                               plotdf.summary$Cycle)
  
  statdf.summary <- data.frame(Cycle = as.integer(rownames(means)), 
                               Mean = means)
  return(statdf.summary)
}
# ==== plotAggregateLengths function ====
plotAggregateLengths = function (fl, n = 5e+05, table = TRUE){
  plotdf = lapply(fl[!is.na(fl)], function(x){
    sr <- ShortRead::FastqSampler(x, n = n)
    df <- yield(sr)
    df <- width(df)
  })
  plotdf = do.call(c, plotdf)
  
  if(table){
    print(table(plotdf))
  }else{
    plot = ggplot()+
      geom_histogram(bins = max(plotdf)-min(plotdf)+1,aes(x = plotdf))+
      scale_x_continuous(breaks = seq.int(min(plotdf), max(plotdf)))
    print(plot)
  }
}
# ==== plotQualityProfile3 function ====
plotQualityProfile3 = function(fnFs,fnRs,primerLenF,primerLenR,seqlen,trimLeftSelect,truncLenSelect, lenLimF, lenLimR){
  print("Getting aggregate quality scores...")
  qpF = getAggregateQualityScores(fnFs)
  qpR = getAggregateQualityScores(fnRs)
  # qpF = qpF1
  # qpR = qpR1
  
  qpF = qpF[1:lenLimF,]
  qpR = qpR[1:lenLimR,]
  
  qp = data.frame(CycleF = qpF$Cycle, MeanF = qpF$Mean) %>%
    full_join(data.frame(CycleF = seqlen-qpR$Cycle+1, CycleR = qpR$Cycle, MeanR = qpR$Mean))
  
  qp = qp[order(qp$CycleF),]
  
  FRoverlap = which(!is.na(qp$MeanF) & !is.na(qp$MeanR))
  
  # multiply log10(MeanF) and log10(MeanR) by 10^6, then truncate to integers
  # perform remaining calculations with transformed values to avoid floating point number errors
  #divide all values by 10^6 for display and output at the end
  qp = qp %>%
    mutate(MeanF_ = trunc(log10(MeanF)*10^6)) %>%
    mutate(MeanR_ = trunc(log10(MeanR)*10^6))
  
  # calculate cummulative sum of log10(mean error-free rate)
  qp$csumMeanF = NA
  qp$csumMeanF[which(!is.na(qp$MeanF))] = cumsum(qp$MeanF_[which(!is.na(qp$MeanF))])
  qp$csumMeanR = NA
  qp$csumMeanR[rev(which(!is.na(qp$MeanR)))] = cumsum(qp$MeanR_[rev(which(!is.na(qp$MeanR)))])
  
  print("Starting length n subarrays...")
  # find length n subarray with maximum subarray sum
  subarray_csumMean = aaply(1:(length(FRoverlap)), 1, function(n){
    csumMean = rep_len(NA, nrow(qp))
    csumMean[FRoverlap[(FRoverlap-FRoverlap[1]+1-n)>=0]] = qp$csumMeanF[FRoverlap[(FRoverlap-FRoverlap[1]+1-n)>=0]] + qp$csumMeanR[(FRoverlap-(n-1))[(FRoverlap-FRoverlap[1]+1-n)>=0]]
    return(csumMean)
  })
  subarray_csumMean = t(subarray_csumMean)/10^6
  
  subarray_csumMeanMax = adply(1:ncol(subarray_csumMean), 1, function(n){
    data.frame(
      n = n,
      score = max(subarray_csumMean[,n],na.rm = TRUE),
      CycleF = qp$CycleF[which.max(subarray_csumMean[,n])],
      CycleR = qp$CycleR[which.max(subarray_csumMean[,n])-n+1]
    )
  })
  
  plot = ggplot()+
    geom_line(data = qp, color = "red", aes(x = CycleF, y = 10^(csumMeanF/10^6)))+
    geom_line(data = qp, color = "blue", aes(x = CycleF, y = 10^(csumMeanR/10^6)))+
    scale_y_continuous(name = "cumulative error-free prob.")+
    geom_vline(xintercept = trimLeftSelect[1], color = "black") + 
    geom_vline(xintercept = truncLenSelect[1], color = "pink") +
    geom_vline(xintercept = seqlen-trimLeftSelect[2], color = "black") + 
    geom_vline(xintercept = seqlen-truncLenSelect[2], color = "lightblue")#+
  # geom_line(data = qp, color = "black", aes(x = CycleF, y = sumMean + mean(qp$MeanF[FRoverlap] + qp$MeanR[FRoverlap])/2 - 15))+
  # scale_y_continuous(sec.axis = sec_axis(~ . -mean(qp$MeanF[FRoverlap] + qp$MeanR[FRoverlap])/2 +15, name = "over-Mean Q"))+
  # geom_line(data = qp, aes(x = CycleF, y = csumMean), color = "green")+
  # geom_vline(xintercept = autotrimLeftSelect[1], color = "red") + 
  # geom_vline(xintercept = seqlen-autotrimLeftSelect[2], color = "blue")
  
  print(plot)
  
  # 
  # print(paste0("Over-mean Q score for supplied parameters is ",
  #              qp$csumMean[which(qp$CycleF == truncLenSelect[1])] - qp$csumMean[which(qp$CycleR == truncLenSelect[2])]
  # ))
  # print(paste0("Fwd and Rev read overlap for supplied parameters is ",
  #              truncLenSelect[1] - ((seqlen + primerLenF + primerLenR) - (truncLenSelect[2])),
  #              " bases"
  # ))
  # invisible(qp)
  invisible(subarray_csumMeanMax)
}


# ==== seq2fasta conversion function ====
seq2fasta = function(seqNames, sequences, fileName, savePath){
  fasta = paste0(">",seqNames,"\n",sequences)
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
# ==== Kalign_pim function ====
#REST API access to Kalign
Kalign_pim = function(userEmail,
                      fastafileName,
                      
                      waitTime = 10,
                      
                      outDir = dirname(fastafileName),
                      outFile = NULL){
  
  if(is.null(outFile)){
    outFile = paste0(outDir,"/",basename(fastafileName),"_Kalign_out.txt")
  }
  params = list(
    email = userEmail,
    stype = "dna",
    format = "fasta",
    sequence = readChar(fastafileName, file.info(fastafileName)$size))
  
  kalignRunResponse <- POST("https://www.ebi.ac.uk/Tools/services/rest/kalign/run", body = params)
  kalignJobID = content(kalignRunResponse)
  
  # check job status
  kalignJobStatusResponse = "Checking"
  while (kalignJobStatusResponse != "FINISHED") {
    kalignJobStatuspb <- txtProgressBar(min = 0, max = waitTime, style = 3)
    for(i in 1:waitTime){
      Sys.sleep(1)
      # update progress bar
      setTxtProgressBar(kalignJobStatuspb, i)
    }
    close(kalignJobStatuspb)
    kalignJobStatus = GET(paste0("https://www.ebi.ac.uk/Tools/services/rest/kalign/status/",kalignJobID))
    kalignJobStatusResponse = content(kalignJobStatus)
  }
  
  # GET response
  kalignJob_pimResponse = GET(paste0("https://www.ebi.ac.uk/Tools/services/rest/kalign/result/",kalignJobID,"/pim"))
  kalignJob_pim = content(kalignJob_pimResponse)
  
  # format response
  kalignJob_pim = read.delim(header = FALSE, skip = 6, sep = "", text = kalignJob_pim, stringsAsFactors = FALSE)
  kalignJob_pim = kalignJob_pim[,-c(1:2)]
  seqNames = ShortRead::id(readFasta(fastafileName))
  rownames(kalignJob_pim) = seqNames
  colnames(kalignJob_pim) = seqNames
  
  write.table(kalignJob_pim, 
              file = outFile,
              quote = FALSE,
              sep = "\t",
              row.names = TRUE,
              col.names = FALSE
  )
  return(kalignJob_pim)
}
# ==== primer trimming ====
# library(ShortRead)
# 
# # Forward:  CCTACGGGAGGCAGCAG
# # Reverse:  CTACHVGGGTWTCTAAT
primerTrim = function(path, fnFs, fnRs, primerF, primerR, linkerFlen = 20, linkerRlen = 20){
  primerFseq = DNAString(paste0(paste0(rep.int("N",linkerFlen), collapse = ""),primerF, collapse = ""))
  primerRseq = DNAString(paste0(paste0(rep.int("N",linkerRlen), collapse = ""),primerR, collapse = ""))
  
  dir.create(path = file.path(path,"primer trimmed"), showWarnings = TRUE)
  
  fnFs2 = file.path(path,"primer trimmed",basename(fnFs))
  fnRs2 = file.path(path,"primer trimmed",basename(fnRs))
  
  for(i in 1:length(fnRs2)){
    R1Stream = FastqStreamer(fnFs[i])
    R2Stream = FastqStreamer(fnRs[i])
    on.exit(close(R1Stream))
    on.exit(close(R2Stream))
    
    repeat{
      R1Chunk <- yield(R1Stream)
      R2Chunk <- yield(R2Stream)
      if (
        length(R1Chunk) == 0 |
        length(R2Chunk) == 0){
        break
      }
      
      trimF = trimLRPatterns(Lpattern = primerFseq, subject = R1Chunk, Lfixed = FALSE, ranges = FALSE)
      trimR = trimLRPatterns(Lpattern = primerRseq, subject = R2Chunk, Lfixed = FALSE, ranges = FALSE)
      
      #write to R1 file
      writeFastq(file = fnFs2[i],
                 object = trimF,
                 mode = "a"
      )
      
      #write to R2 file
      writeFastq(file = fnRs2[i],
                 object = trimR,
                 mode = "a"
      )
    }
  }
}

# ==== pseudo-log10 function ====
pseudo_log_breaks = function(n = 5, sigma = 1, base = 10){
  force(n)
  force(sigma)
  force(base)
  n_default = n
  
  function(x, n = n_default) {
    raw_rng <- suppressWarnings(range(x, na.rm = TRUE))
    if (any(!is.finite(raw_rng))) {
      return(numeric())
    }
    pos_x = raw_rng[raw_rng>0]
    neg_x = raw_rng[raw_rng<0]
    
    if(length(pos_x)>0){
      pos_breaks = log_breaks(n, base)(c(sigma,pos_x))
    }else{
      pos_breaks = numeric()
    }
    
    if(length(neg_x)>0){
      neg_breaks = -(log_breaks(n, base)(c(sigma,abs(neg_x))))
    }else{
      neg_breaks = numeric()
    }
    
    all_breaks = c(rev(neg_breaks),pos_breaks)
    
    return(c(rev(neg_breaks),pos_breaks))
  }
}

pseudo_log10_trans = function (sigma = 1, base = 10){
  force(sigma)
  force(base)
  trans_new(name = "pseudo_log10",
            transform = function(x) asinh(x/(2 * sigma))/log(base), 
            inverse = function(x) 2 * sigma * sinh(x * log(base)),
            breaks = pseudo_log_breaks(sigma = sigma,base = base)
  )
}
# when using pseudo_log10_trans transform:
# the scale will generally have reasonable resolution from (-Inf,-sigma)|(sigma,Inf)
# with ~linear behaviour between -2*sigma and 2*sigma, (~ y = 0.1913878/sigma * x) (slope = asinh(1)/(log(10)*2))
# with ~0.4 scale units between -sigma and sigma

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
# ==== unflipping R1 and R2 reads (repartitioning R1 and R2 reads) ====
unFlip = function(path,
                  fnFs, fnRs,
                  id.sep = "\\s"){
  dir.create(path = paste0(path,"/unflipped/"), showWarnings = TRUE)
  
  test = adply(1:length(fnFs), .margins = 1, .id = NULL, .fun = function(i){
    flippedOutR1 = 0
    readCount = 0
    
    R1Stream = FastqStreamer(fnFs[i])
    R2Stream = FastqStreamer(fnRs[i])
    
    on.exit(close(R1Stream))
    on.exit(close(R2Stream))
    
    repeat{
      R1Chunk <- yield(R1Stream)
      R2Chunk <- yield(R2Stream)
      
      if(length(R1Chunk) == 0 | length(R2Chunk) == 0){
        break
      }
      
      # extract sequence identifier field
      R1seqid = sapply(strsplit(as.character(id(R1Chunk)), id.sep), `[`, 2)
      R2seqid = sapply(strsplit(as.character(id(R2Chunk)), id.sep), `[`, 2)
      
      flipSelectR1 = grep("2",R1seqid, fixed = TRUE)
      flipSelectR2 = grep("1",R2seqid, fixed = TRUE)
      
      flippedOutR1 = flippedOutR1 + length(flipSelectR1)
      readCount = readCount + length(R1seqid)
      
      if(!identical(flipSelectR1,flipSelectR2)){
        stop(print("R1 & R2 reads are not matched."))
      }
      
      tempChunk = R1Chunk[flipSelectR1]
      R1Chunk[flipSelectR1] = R2Chunk[flipSelectR2]
      R2Chunk[flipSelectR1] = tempChunk
      rm(tempChunk)
      
      # write to R1 file
      writeFastq(file = paste0(path,"/unflipped/",basename(fnFs[i])),
                 object = R1Chunk,
                 mode = "a"
      )
      
      # write to R2 file
      writeFastq(file = paste0(path,"/unflipped/",basename(fnRs[i])),
                 object = R2Chunk,
                 mode = "a"
      )
    }
    return(data.frame(fnF = fnFs[i],readCount,flippedOutR1))
  })
}


# ==== mixedOrientSortTrim ====
mixedOrientSortTrim = function(seqPath, fnR1s, fnR2s, 
                               primerF, primerR, 
                               maxMismatchesF, maxMismatchesR,
                               minMatchLenF,minMatchLenR,
                               linkerFlen = 20, linkerRlen = 20){
  primerFseq = DNAString(paste0(paste0(rep.int("N",linkerFlen), collapse = ""),primerF, collapse = ""))
  primerRseq = DNAString(paste0(paste0(rep.int("N",linkerRlen), collapse = ""),primerR, collapse = ""))
  
  maxMismatchesF = c(maxMismatchesF,rep.int(maxMismatchesF[length(maxMismatchesF)],linkerFlen))
  maxMismatchesR = c(maxMismatchesR,rep.int(maxMismatchesR[length(maxMismatchesR)],linkerRlen))
  
  dir.create(path = file.path(seqPath,"mixed orientation sort and trimmed"), showWarnings = TRUE)
  
  sampleNames = sapply(strsplit(basename(fnR1s),"_"), `[`, 1)
  suffixR1Names = sapply(strsplit(basename(fnR1s[1]),sampleNames), `[`, 2)
  suffixR2Names = sapply(strsplit(basename(fnR2s[1]),sampleNames), `[`, 2)
  
  fnR1Fs = file.path(seqPath,"mixed orientation sort and trimmed",paste0(sampleNames,".orientF",suffixR1Names))
  fnR2Rs = file.path(seqPath,"mixed orientation sort and trimmed",paste0(sampleNames,".orientF",suffixR2Names))
  
  fnR1Rs = file.path(seqPath,"mixed orientation sort and trimmed",paste0(sampleNames,".orientR",suffixR1Names))
  fnR2Fs = file.path(seqPath,"mixed orientation sort and trimmed",paste0(sampleNames,".orientR",suffixR2Names))
  
  for(i in 1:length(fnR1s)){
    R1Stream = FastqStreamer(fnR1s[i])
    R2Stream = FastqStreamer(fnR2s[i])
    on.exit(close(R1Stream))
    on.exit(close(R2Stream))
    
    repeat{
      R1Chunk <- yield(R1Stream)
      R2Chunk <- yield(R2Stream)
      if (
        length(R1Chunk) == 0 |
        length(R2Chunk) == 0){
        break
      }
      
      trimR1_F = trimLRPatterns(Lpattern = primerFseq, subject = R1Chunk, max.Lmismatch = maxMismatchesF, Lfixed = FALSE, ranges = TRUE)
      trimR2_R = trimLRPatterns(Lpattern = primerRseq, subject = R2Chunk, max.Lmismatch = maxMismatchesR, Lfixed = FALSE, ranges = TRUE)
      
      trimR1_R = trimLRPatterns(Lpattern = primerRseq, subject = R1Chunk, max.Lmismatch = maxMismatchesR, Lfixed = FALSE, ranges = TRUE)
      trimR2_F = trimLRPatterns(Lpattern = primerFseq, subject = R2Chunk, max.Lmismatch = maxMismatchesF, Lfixed = FALSE, ranges = TRUE)
      
      
      keepR1_F <- (trimR1_F@start >= minMatchLenF | trimR2_R@start >= minMatchLenR) & trimR1_F@width >0 & trimR2_R@width >0
      keepR1_R <- (trimR2_F@start >= minMatchLenF | trimR1_R@start >= minMatchLenR) & trimR2_F@width >0 & trimR1_R@width >0 & !keepR1_F
      
      
      R1ChunkF = R1Chunk[keepR1_F]
      R2ChunkR = R2Chunk[keepR1_F]
      
      R1ChunkR = R1Chunk[keepR1_R]
      R2ChunkF = R2Chunk[keepR1_R]
      
      
      R1ChunkF <- narrow(R1ChunkF,start = trimR1_F@start[keepR1_F], width = trimR1_F@width[keepR1_F])
      R2ChunkR <- narrow(R2ChunkR,start = trimR2_R@start[keepR1_F], width = trimR2_R@width[keepR1_F])
      
      R1ChunkR <- narrow(R1ChunkR,start = trimR1_R@start[keepR1_R], width = trimR1_R@width[keepR1_R])
      R2ChunkF <- narrow(R2ChunkF,start = trimR2_F@start[keepR1_R], width = trimR2_F@width[keepR1_R])
      
      
      # write to .orientF files
      writeFastq(file = fnR1Fs[i],
                 object = R1ChunkF,
                 mode = "a"
      )
      
      writeFastq(file = fnR2Rs[i],
                 object = R2ChunkR,
                 mode = "a"
      )
      
      # write to .orientR files
      writeFastq(file = fnR1Rs[i],
                 object = R1ChunkR,
                 mode = "a"
      )
      
      writeFastq(file = fnR2Fs[i],
                 object = R2ChunkF,
                 mode = "a"
      )
    }
  }
}

# ==== sortFromFlipped ====
sortFromFlipped = function(seqPath,
                           fnFs, fnRs,
                           id.sep = "\\s"){
  dir.create(path = file.path(seqPath,"mixed orientation sort and trimmed"), showWarnings = TRUE)
  
  sampleNames = sapply(strsplit(basename(fnFs),"_"), `[`, 1)
  suffixR1Names = sapply(strsplit(basename(fnFs[1]),sampleNames), `[`, 2)
  suffixR2Names = sapply(strsplit(basename(fnRs[1]),sampleNames), `[`, 2)
  
  fnR1Fs = file.path(seqPath,"mixed orientation sort and trimmed",paste0(sampleNames,".orientF",suffixR1Names))
  fnR2Rs = file.path(seqPath,"mixed orientation sort and trimmed",paste0(sampleNames,".orientF",suffixR2Names))
  
  fnR1Rs = file.path(seqPath,"mixed orientation sort and trimmed",paste0(sampleNames,".orientR",suffixR1Names))
  fnR2Fs = file.path(seqPath,"mixed orientation sort and trimmed",paste0(sampleNames,".orientR",suffixR2Names))
  
  for(i in 1:length(fnFs)){
    FStream = FastqStreamer(fnFs[i])
    RStream = FastqStreamer(fnRs[i])
    
    on.exit(close(FStream))
    on.exit(close(RStream))
    
    repeat{
      FChunk <- yield(FStream)
      RChunk <- yield(RStream)
      
      if(length(FChunk) == 0 | length(RChunk) == 0){
        break
      }
      
      # extract sequence identifier field
      Fseqid = sapply(strsplit(as.character(id(FChunk)), id.sep), `[`, 2)
      Rseqid = sapply(strsplit(as.character(id(RChunk)), id.sep), `[`, 2)
      
      SelectR1F = grep("1",Fseqid, fixed = TRUE)
      SelectR2R = grep("2",Rseqid, fixed = TRUE)
      
      SelectR1R = grep("1",Rseqid, fixed = TRUE)
      SelectR2F = grep("2",Fseqid, fixed = TRUE)
      
      FChunkR1 = FChunk[SelectR1F]
      RChunkR2 = RChunk[SelectR2R]
      
      RChunkR1 = RChunk[SelectR1R]
      FChunkR2 = FChunk[SelectR2F]
      
      
      # write to .orientF files
      writeFastq(file = fnR1Fs[i],
                 object = FChunkR1,
                 mode = "a"
      )
      
      writeFastq(file = fnR2Rs[i],
                 object = RChunkR2,
                 mode = "a"
      )
      
      # write to .orientR files
      writeFastq(file = fnR1Rs[i],
                 object = RChunkR1,
                 mode = "a"
      )
      
      writeFastq(file = fnR2Fs[i],
                 object = FChunkR2,
                 mode = "a"
      )
    }
  }
}



# ==== match read IDs ====
readMatch = function(path,
                     fnR1s, fnR2s,
                     id.sep = "\\s"){
  dir.create(path = paste0(path,"/reads matched/"), showWarnings = TRUE)
  
  for(i in 1:length(fnR1s)){
    R1Stream = FastqStreamer(fnR1s[i])
    R2Stream = FastqStreamer(fnR2s[i])
    
    on.exit(close(R1Stream))
    on.exit(close(R2Stream))
    
    R1ChunkRemainder = ShortReadQ()
    R2ChunkRemainder = ShortReadQ()
    
    repeat{
      R1Chunk <- yield(R1Stream)
      R2Chunk <- yield(R2Stream)
      
      R1ChunkRemainder = ShortRead::append(R1ChunkRemainder,R1Chunk)
      R2ChunkRemainder = ShortRead::append(R2ChunkRemainder,R2Chunk)
      
      R1ChunkSize = length(R1Chunk)
      R2ChunkSize = length(R2Chunk)
      
      rm(R1Chunk)
      rm(R2Chunk)
      
      # extract sequence identifier field
      R1seqid = sapply(strsplit(as.character(id(R1ChunkRemainder)), id.sep), `[`, 1)
      R2seqid = sapply(strsplit(as.character(id(R2ChunkRemainder)), id.sep), `[`, 1)
      
      seqidMatches = match(R1seqid,R2seqid, nomatch = 0)
      
      writeR1Chunk = R1ChunkRemainder[which(seqidMatches!=0)]
      writeR2Chunk = R2ChunkRemainder[seqidMatches[seqidMatches!=0]]
      
      R1ChunkRemainder = R1ChunkRemainder[which(seqidMatches==0)]
      R2ChunkRemainder = R2ChunkRemainder[-(seqidMatches[seqidMatches!=0])]
      
      # write to R1 file
      writeFastq(file = file.path(path,"reads matched",basename(fnR1s[i])),
                 object = writeR1Chunk,
                 mode = "a"
      )
      
      # write to R2 file
      writeFastq(file = file.path(path,"reads matched",basename(fnR2s[i])),
                 object = writeR2Chunk,
                 mode = "a"
      )
      
      if(R1ChunkSize == 0 & R2ChunkSize == 0){
        break
      }
    }
  }
}

# ==== fecal seq. data derep object generation ====
path <- "./Alessandra_16s/NS0_NS1 analysis"
seqlen = 254

priorsPrimerFseq = "GTGCCAGCMGCCGCGGTAA"
priorsPrimerRseq = "GGACTACHVGGGTWTCTAA"

#  ==== Fecal & chemostat isolates cleaning for priors ==== 
NS0_F_C_Isolates = read.delim(file.path(path,"input/isolates/NS0 F_C_Isolates.txt"),header = TRUE,sep = "\t", stringsAsFactors = FALSE)
NS1_F_C_Isolates = read.delim(file.path(path,"input/isolates/NS1 F_C_Isolates.txt"),header = TRUE,sep = "\t", stringsAsFactors = FALSE)

# "length"
# "isolateAlignF_score"
# "isolateAlignR_score"
# "V4Seq"
# "matchMissing"
# "N_NumV4Seq"
# "rcV4Seq"
# "trimLengthF"
# "trimLengthR"

# ---- NS0 ----
NS0_F_C_IsolatesCleaned = NS0_F_C_Isolates
NS0_F_C_IsolatesCleaned$sangerSeq = gsub("[^ACGTN]", "", NS0_F_C_IsolatesCleaned$sangerSeq)
NS0_F_C_IsolatesCleaned$length = nchar(NS0_F_C_IsolatesCleaned$sangerSeq)

rcPrimerFseq = as.character(reverseComplement(DNAString(priorsPrimerFseq)))

mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = FALSE)

isolateAlignF = pairwiseAlignment(NS0_F_C_IsolatesCleaned$sangerSeq, rcPrimerFseq, type = "overlap", substitutionMatrix = mat)
isolateAlignR = pairwiseAlignment(NS0_F_C_IsolatesCleaned$sangerSeq, priorsPrimerRseq, type = "overlap", substitutionMatrix = mat)

NS0_F_C_IsolatesCleaned$isolateAlignF_score = score(isolateAlignF)
NS0_F_C_IsolatesCleaned$isolateAlignR_score = score(isolateAlignR)

table(score(isolateAlignF))
table(score(isolateAlignR))

isolateAlignFtest = data.frame(seq = as.character(alignedPattern(isolateAlignF)), score = score(isolateAlignF))
isolateAlignFtest = isolateAlignFtest[!duplicated(isolateAlignFtest$seq),]

isolateAlignRtest = data.frame(seq = as.character(alignedPattern(isolateAlignR)), score = score(isolateAlignR))
isolateAlignRtest = isolateAlignRtest[!duplicated(isolateAlignRtest$seq),]

V4SeqEnd = NS0_F_C_IsolatesCleaned$length
V4SeqEnd[nchar(as.character(alignedPattern(isolateAlignF)))>10] = (start(pattern(isolateAlignF))-1)[nchar(as.character(alignedPattern(isolateAlignF)))>10]

V4SeqStart = rep(1,nrow(NS0_F_C_IsolatesCleaned))
V4SeqStart[nchar(as.character(alignedPattern(isolateAlignR)))>10] = (end(pattern(isolateAlignR))+1)[nchar(as.character(alignedPattern(isolateAlignR)))>10]

any(!(V4SeqEnd>V4SeqStart))

NS0_F_C_IsolatesCleaned$V4Seq = as.character(narrow(DNAStringSet(NS0_F_C_IsolatesCleaned$sangerSeq), start = V4SeqStart, end = V4SeqEnd))

write.table(NS0_F_C_IsolatesCleaned$V4Seq, "./Alessandra_16s/NS0_NS1 analysis/input/isolates/NS0_V4seq.txt", row.names = FALSE, quote = FALSE)

V4SeqManual = read.delim("./Alessandra_16s/NS0_NS1 analysis/input/isolates/NS0_V4seqManualClean.txt", header = FALSE)
NS0_F_C_IsolatesCleaned$V4SeqManual = V4SeqManual$V1
NS0_F_C_IsolatesCleaned$N_NumV4Seq = laply(gregexpr("N", NS0_F_C_IsolatesCleaned$V4Seq),function(x) sum(x!=-1))
NS0_F_C_IsolatesCleaned$N_NumV4SeqManual = laply(gregexpr("N", NS0_F_C_IsolatesCleaned$V4SeqManual),function(x) sum(x!=-1))
NS0_F_C_IsolatesCleaned$rcV4Seq = as.character(reverseComplement(DNAStringSet(NS0_F_C_IsolatesCleaned$V4Seq)))
NS0_F_C_IsolatesCleaned$rcV4SeqManual = as.character(reverseComplement(DNAStringSet(NS0_F_C_IsolatesCleaned$V4SeqManual)))
write.table(NS0_F_C_IsolatesCleaned, file.path(path,"input/isolates/NS0 F_C_IsolatesCleaned.txt"), row.names = FALSE, quote = FALSE, sep = "\t")

# ---- NS1 ----
NS1_F_C_IsolatesCleaned = NS1_F_C_Isolates
NS1_F_C_IsolatesCleaned$sangerSeq = gsub("[^ACGTN]", "", NS1_F_C_IsolatesCleaned$sangerSeq)
NS1_F_C_IsolatesCleaned$length = nchar(NS1_F_C_IsolatesCleaned$sangerSeq)

rcPrimerFseq = as.character(reverseComplement(DNAString(priorsPrimerFseq)))

mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = FALSE)

isolateAlignF = pairwiseAlignment(NS1_F_C_IsolatesCleaned$sangerSeq, rcPrimerFseq, type = "overlap", substitutionMatrix = mat)
isolateAlignR = pairwiseAlignment(NS1_F_C_IsolatesCleaned$sangerSeq, priorsPrimerRseq, type = "overlap", substitutionMatrix = mat)

NS1_F_C_IsolatesCleaned$isolateAlignF_score = score(isolateAlignF)
NS1_F_C_IsolatesCleaned$isolateAlignR_score = score(isolateAlignR)

table(score(isolateAlignF))
table(score(isolateAlignR))

isolateAlignFtest = data.frame(seq = as.character(alignedPattern(isolateAlignF)), score = score(isolateAlignF))
isolateAlignFtest = isolateAlignFtest[!duplicated(isolateAlignFtest$seq),]

isolateAlignRtest = data.frame(seq = as.character(alignedPattern(isolateAlignR)), score = score(isolateAlignR))
isolateAlignRtest = isolateAlignRtest[!duplicated(isolateAlignRtest$seq),]

V4SeqEnd = NS1_F_C_IsolatesCleaned$length
V4SeqEnd[nchar(as.character(alignedPattern(isolateAlignF)))>10] = (start(pattern(isolateAlignF))-1)[nchar(as.character(alignedPattern(isolateAlignF)))>10]

V4SeqStart = rep(1,nrow(NS1_F_C_IsolatesCleaned))
V4SeqStart[nchar(as.character(alignedPattern(isolateAlignR)))>10] = (end(pattern(isolateAlignR))+1)[nchar(as.character(alignedPattern(isolateAlignR)))>10]

any(!(V4SeqEnd>V4SeqStart))

NS1_F_C_IsolatesCleaned$V4Seq = as.character(narrow(DNAStringSet(NS1_F_C_IsolatesCleaned$sangerSeq), start = V4SeqStart, end = V4SeqEnd))

write.table(NS1_F_C_IsolatesCleaned$V4Seq, "./Alessandra_16s/NS0_NS1 analysis/input/isolates/NS1_V4seq.txt", row.names = FALSE, quote = FALSE)

V4SeqManual = read.delim("./Alessandra_16s/NS0_NS1 analysis/input/isolates/NS1_V4seqManualClean.txt", header = FALSE)
NS1_F_C_IsolatesCleaned$V4SeqManual = V4SeqManual$V1
NS1_F_C_IsolatesCleaned$N_NumV4Seq = laply(gregexpr("N", NS1_F_C_IsolatesCleaned$V4Seq),function(x) sum(x!=-1))
NS1_F_C_IsolatesCleaned$N_NumV4SeqManual = laply(gregexpr("N", NS1_F_C_IsolatesCleaned$V4SeqManual),function(x) sum(x!=-1))
NS1_F_C_IsolatesCleaned$rcV4Seq = as.character(reverseComplement(DNAStringSet(NS1_F_C_IsolatesCleaned$V4Seq)))
NS1_F_C_IsolatesCleaned$rcV4SeqManual = as.character(reverseComplement(DNAStringSet(NS1_F_C_IsolatesCleaned$V4SeqManual)))
write.table(NS1_F_C_IsolatesCleaned, file.path(path,"input/isolates/NS1 F_C_IsolatesCleaned.txt"), row.names = FALSE, quote = FALSE, sep = "\t")


# ==== extracting priors ====
NS0_F_C_IsolatesCleaned = read.delim(file.path(path,"input/isolates/NS0 F_C_IsolatesCleaned.txt"),header = TRUE,sep = "\t", stringsAsFactors = FALSE)
NS1_F_C_IsolatesCleaned = read.delim(file.path(path,"input/isolates/NS1 F_C_IsolatesCleaned.txt"),header = TRUE,sep = "\t", stringsAsFactors = FALSE)

# NS0_priorsFecal = NS0_F_C_IsolatesCleaned$rcV4Seq[NS0_F_C_IsolatesCleaned$source == "fecal"]
# NS0_priorsChemostat = NS0_F_C_IsolatesCleaned$rcV4Seq[NS0_F_C_IsolatesCleaned$source == "chemostat"]
# NS0_priorsDC = NS0_F_C_IsolatesCleaned$rcV4Seq[NS0_F_C_IsolatesCleaned$DC_inoculum == TRUE]
# 
# NS1_priorsFecal = NS1_F_C_IsolatesCleaned$rcV4Seq[NS1_F_C_IsolatesCleaned$source == "fecal"]
# NS1_priorsChemostat = NS1_F_C_IsolatesCleaned$rcV4Seq[NS1_F_C_IsolatesCleaned$source == "chemostat"]
# NS1_priorsDC = NS1_F_C_IsolatesCleaned$rcV4Seq[NS1_F_C_IsolatesCleaned$DC_inoculum == TRUE]

NS0_priors = NS0_F_C_IsolatesCleaned$rcV4Seq
NS0_priorsDC = NS0_F_C_IsolatesCleaned$rcV4Seq[NS0_F_C_IsolatesCleaned$DC_inoculum == TRUE]

NS1_priors = NS1_F_C_IsolatesCleaned$rcV4Seq
NS1_priorsDC = NS1_F_C_IsolatesCleaned$rcV4Seq[NS1_F_C_IsolatesCleaned$DC_inoculum == TRUE]

#  ==== 2019 09 05 sequencing run samples ==== 
fnR1s <- sort(list.files(file.path(path,"input/sequencing 20190905/raw"), pattern="R1_001.fastq", full.names = TRUE))
fnR2s <- sort(list.files(file.path(path,"input/sequencing 20190905/raw"), pattern="R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnR1s), "_R"), `[`, 1)

# separate R1 and R2 reads into R1_F, R1_R & R2_F, R2_R
maxMismatchesF = c(1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
maxMismatchesF = rev(maxMismatchesF)
minMatchLenF = 12

maxMismatchesR = c(1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0)
maxMismatchesR = rev(maxMismatchesR)
minMatchLenR = 4

seqPath = "./Alessandra_16s/NS0_NS1 analysis/input/sequencing 20190905"

mixedOrientSortTrim(seqPath, fnR1s, fnR2s,
                    priorsPrimerFseq, priorsPrimerRseq,
                    maxMismatchesF, maxMismatchesR,
                    minMatchLenF,minMatchLenR,
                    linkerFlen = 20, linkerRlen = 20)

fnR1s <- sort(list.files(file.path(path,"input/sequencing 20190905/mixed orientation sort and trimmed"), pattern="R1_001.fastq", full.names = TRUE))
fnR2s <- sort(list.files(file.path(path,"input/sequencing 20190905/mixed orientation sort and trimmed"), pattern="R2_001.fastq", full.names = TRUE))

# primerLenF = 19
# primerLenR = 20
# trimLeft = c(8,0)
primerLenF = 0
primerLenR = 0
trimLeft = c(0,0)

trimLeftSelect = c(primerLenF+trimLeft[1],primerLenR+trimLeft[2])
truncLenSelect = c(195,79)

priorsR1 = list(NS0_priorsDC,NS0_priorsDC,NS0_priorsDC,NS0_priorsDC,NS0_priorsDC,NS0_priorsDC,NS0_priorsDC,NS0_priorsDC,NS0_priorsDC,NS0_priorsDC,NS0_priorsDC,NS0_priorsDC,NS0_priorsDC,
                NS0_priors,
                NS0_priors,NS0_priors,NS0_priors,NS0_priors,NS0_priors,NS0_priors,
                NS1_priorsDC,NS1_priorsDC,NS1_priorsDC,NS1_priorsDC,
                NS1_priors
                )

priorsR1 = c(priorsR1, lapply(priorsR1, function(x) as.character(reverseComplement(DNAStringSet(x)))))
names(priorsR1) = c(paste0(sample.names,".orientF"),paste0(sample.names,".orientR"))

priorsR2 = llply(priorsR1, function(x){
  y = as.character(reverseComplement(DNAStringSet(x)))
  y = substr(y,trimLeftSelect[2]+1,truncLenSelect[2])
  y = y[laply(gregexpr("N", y),function(z) sum(z!=-1) == 0)]
  y = unique(y)
})

priorsR1 = llply(priorsR1, function(x){
  y = substr(x,trimLeftSelect[1]+1,truncLenSelect[1])
  y = y[laply(gregexpr("N", y),function(z) sum(z!=-1) == 0)]
  y = unique(y)
})

# priorsF = character(0)
# priorsR = character(0)

#Legend:
# green = mean
# orange = median
# dashed orange = 25th and 75th quantiles
# first black line = start of filtered reads
# second black line = end of filtered reads
# red line = overlap with opposite reads

plotQualityProfile(fnR1s[7], aggregate = FALSE)
plotQualityProfile(fnR2s[7], aggregate = FALSE)
plotAggregateLengths(fnR1s[1:20])
plotAggregateLengths(fnR2s[1:20])
test = plotQualityProfile3(fnR1s[1:20],fnR2s[1:20],primerLenF,primerLenR,seqlen,trimLeftSelect,truncLenSelect, lenLimF = 236, lenLimR = 236, overlapLen = 20)

sample.names <- sapply(strsplit(basename(fnR1s), "_R"), `[`, 1)

# process priors only
process(path,
        fnR1s, fnR2s, sample.names, errPoolName = "20190905", orientFR.split = TRUE,
        trimLeftSelect, truncLenSelect, filter.matchIDs = FALSE,
        ErrModelMonotonicity = FALSE,
        OMEGA_A = getDadaOpt(option = "OMEGA_A"), pool = FALSE, poolList = NULL,
        priorsF = priorsR1, priorsR = priorsR2,
        
        plotToggle = FALSE,
        printPriorHead = FALSE,
        
        preFilteredToggle = TRUE,
        preCalcErrToggle = TRUE,
        preDerepToggle = TRUE,
        preDenoiseToggle = FALSE,
        
        saveFiltered = TRUE,
        saveDereplicated = TRUE,
        saveDenoised = TRUE)

#  ==== 2019 05 03 sequencing run samples ==== 
fnFs <- sort(list.files(file.path(path,"input/sequencing 20190503/raw"), pattern="R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(file.path(path,"input/sequencing 20190503/raw"), pattern="R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

seqPath = "./Alessandra_16s/NS0_NS1 analysis/input/sequencing 20190503"

sortFromFlipped(seqPath, fnFs, fnRs, id.sep = "\\s")

fnR1s <- sort(list.files(file.path(path,"input/sequencing 20190503/mixed orientation sort and trimmed"), pattern="R1_001.fastq", full.names = TRUE))
fnR2s <- sort(list.files(file.path(path,"input/sequencing 20190503/mixed orientation sort and trimmed"), pattern="R2_001.fastq", full.names = TRUE))

# match read IDs
readMatch(seqPath, fnR1s, fnR2s, id.sep = "\\s")

fnR1s <- sort(list.files(file.path(path,"input/sequencing 20190503/reads matched"), pattern="R1_001.fastq", full.names = TRUE))
fnR2s <- sort(list.files(file.path(path,"input/sequencing 20190503/reads matched"), pattern="R2_001.fastq", full.names = TRUE))

primerLenF = 0
primerLenR = 0
trimLeft = c(0,0)
trimLeftSelect = c(primerLenF+trimLeft[1],primerLenR+trimLeft[2])
truncLenSelect = c(158,116)

priorsR1 = list(NS1_priorsDC,NS1_priorsDC,NS1_priorsDC,NS1_priorsDC,NS1_priorsDC,NS1_priorsDC,NS1_priorsDC,NS1_priorsDC,NS1_priorsDC,NS1_priorsDC,NS1_priorsDC,NS1_priorsDC,NS1_priorsDC,NS1_priorsDC,NS1_priorsDC,NS1_priorsDC,NS1_priorsDC,
                NS1_priors,NS1_priors,NS1_priors,NS1_priors,NS1_priors,NS1_priors,NS1_priors,NS1_priors
)

# orientF.R2 and orientR.R1 priors need to start at base #2; 3' end is trimmed 1 base in
priorsR2 = c(lapply(priorsR1, function(x) as.character(narrow(DNAStringSet(x), start = 1, end = width(DNAStringSet(x))-1))),
             lapply(priorsR1, function(x) as.character(reverseComplement(DNAStringSet(x)))))
names(priorsR2) = c(paste0(sample.names,".orientF"),paste0(sample.names,".orientR"))
priorsR1 = c(priorsR1, lapply(priorsR1, function(x) as.character(narrow(reverseComplement(DNAStringSet(x)), start = 2, end = width(reverseComplement(DNAStringSet(x)))))))
names(priorsR1) = c(paste0(sample.names,".orientF"),paste0(sample.names,".orientR"))

priorsR2 = llply(priorsR2, function(x){
  y = as.character(reverseComplement(DNAStringSet(x)))
  y = substr(y,trimLeftSelect[2]+1,truncLenSelect[2])
  y = y[laply(gregexpr("N", y),function(z) sum(z!=-1) == 0)]
  y = unique(y)
})

priorsR1 = llply(priorsR1, function(x){
  y = substr(x,trimLeftSelect[1]+1,truncLenSelect[1])
  y = y[laply(gregexpr("N", y),function(z) sum(z!=-1) == 0)]
  y = unique(y)
})

# priorsF = character(0)
# priorsR = character(0)

#Legend:
# green = mean
# orange = median
# dashed orange = 25th and 75th quantiles
# first black line = start of filtered reads
# second black line = end of filtered reads
# red line = overlap with opposite reads

# quality scores are binned, check for need to enforce monotonicity
plotQualityProfile(fnR1s[1:20], aggregate = TRUE)
plotQualityProfile(fnR2s[1:20], aggregate = TRUE)
plotAggregateLengths(fnR1s[1:20])
plotAggregateLengths(fnR2s[1:20])
test = plotQualityProfile3(fnR1s[1:20],fnR2s[1:20],primerLenF,primerLenR,seqlen,trimLeftSelect,truncLenSelect, lenLimF = 222, lenLimR = 222)

sample.names <- sapply(strsplit(basename(fnR1s), "_R"), `[`, 1)

# process priors only
process(path,
        fnR1s, fnR2s, sample.names, errPoolName = "20190503", orientFR.split = TRUE,
        trimLeftSelect, truncLenSelect, filter.matchIDs = FALSE,
        ErrModelMonotonicity = FALSE,
        OMEGA_A = getDadaOpt(option = "OMEGA_A"), pool = FALSE, poolList = NULL,
        priorsF = priorsR1, priorsR = priorsR2,
        
        plotToggle = FALSE,
        printPriorHead = FALSE,
        
        preFilteredToggle = TRUE,
        preCalcErrToggle = TRUE,
        preDerepToggle = TRUE,
        preDenoiseToggle = FALSE,
        
        saveFiltered = TRUE,
        saveDereplicated = TRUE,
        saveDenoised = TRUE)

#  ==== manual cross-run pseudopooling ==== 
x20190905_seqtab.nochim = readRDS("./Alessandra_16s/NS0_NS1 analysis/priorsOnly_outputs/20190905_seqtab.nochim.rds")
x20190503_seqtab.nochim = readRDS("./Alessandra_16s/NS0_NS1 analysis/priorsOnly_outputs/20190503_seqtab.nochim.rds")

# 20190503 sequences may be short 1 base on the 3' end
# Pad out with "N"
colnames(x20190503_seqtab.nochim) = paste0(colnames(x20190503_seqtab.nochim),"N")

#merge Sequence Tables
seqtab = mergeSequenceTables(x20190905_seqtab.nochim, x20190503_seqtab.nochim)
seqtab.nochim = removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)
seqtab.nochim.presence = seqtab.nochim>0

seqtab.nochim.presence.NS0_F_C = seqtab.nochim.presence[grep("NS0-F",rownames(seqtab.nochim.presence), fixed = TRUE),]
seqtab.nochim.presence.NS1_F_C = seqtab.nochim.presence[grep("NS1-F",rownames(seqtab.nochim.presence), fixed = TRUE),]
seqtab.nochim.presence.NS0_DC = seqtab.nochim.presence[grep("NS0-DC",rownames(seqtab.nochim.presence), fixed = TRUE),]
seqtab.nochim.presence.NS1_DC = seqtab.nochim.presence[grep("NS1-DC",rownames(seqtab.nochim.presence), fixed = TRUE),]

# PSEUDO_PREVALENCE == 2

seqtab.nochim.presence.NS0_F_C = names(colSums(seqtab.nochim.presence.NS0_F_C)>1)
seqtab.nochim.presence.NS1_F_C = names(colSums(seqtab.nochim.presence.NS1_F_C)>1)
seqtab.nochim.presence.NS0_DC = names(colSums(seqtab.nochim.presence.NS0_DC)>1)
seqtab.nochim.presence.NS1_DC = names(colSums(seqtab.nochim.presence.NS1_DC)>1)

NS0_priorsDC = c(NS0_priorsDC,seqtab.nochim.presence.NS0_DC)
NS0_priors = c(NS0_priors,seqtab.nochim.presence.NS0_F_C)
NS1_priorsDC = c(NS1_priorsDC,seqtab.nochim.presence.NS1_DC)
NS1_priors = c(NS1_priors,seqtab.nochim.presence.NS1_F_C)

# 2019 09 05 sequencing run priors
priorsR1 = list(NS0_priorsDC,NS0_priorsDC,NS0_priorsDC,NS0_priorsDC,NS0_priorsDC,NS0_priorsDC,NS0_priorsDC,NS0_priorsDC,NS0_priorsDC,NS0_priorsDC,NS0_priorsDC,NS0_priorsDC,NS0_priorsDC,
                NS0_priors,
                NS0_priors,NS0_priors,NS0_priors,NS0_priors,NS0_priors,NS0_priors,
                NS1_priorsDC,NS1_priorsDC,NS1_priorsDC,NS1_priorsDC,
                NS1_priors
)

priorsR1 = c(priorsR1, lapply(priorsR1, function(x) as.character(reverseComplement(DNAStringSet(x)))))
names(priorsR1) = c(paste0(sample.names,".orientF"),paste0(sample.names,".orientR"))

priorsR2 = llply(priorsR1, function(x){
  y = as.character(reverseComplement(DNAStringSet(x)))
  y = substr(y,trimLeftSelect[2]+1,truncLenSelect[2])
  y = y[laply(gregexpr("N", y),function(z) sum(z!=-1) == 0)]
  y = unique(y)
})

priorsR1 = llply(priorsR1, function(x){
  y = substr(x,trimLeftSelect[1]+1,truncLenSelect[1])
  y = y[laply(gregexpr("N", y),function(z) sum(z!=-1) == 0)]
  y = unique(y)
})

sample.names <- sapply(strsplit(basename(fnR1s), "_R"), `[`, 1)

process(path,
        fnR1s, fnR2s, sample.names, errPoolName = "20190905", orientFR.split = TRUE,
        trimLeftSelect, truncLenSelect, filter.matchIDs = FALSE,
        ErrModelMonotonicity = FALSE,
        OMEGA_A = getDadaOpt(option = "OMEGA_A"), pool = FALSE, poolList = NULL,
        priorsF = priorsR1, priorsR = priorsR2,
        
        plotToggle = FALSE,
        printPriorHead = FALSE,
        
        preFilteredToggle = TRUE,
        preCalcErrToggle = TRUE,
        preDerepToggle = TRUE,
        preDenoiseToggle = FALSE,
        
        saveFiltered = TRUE,
        saveDereplicated = TRUE,
        saveDenoised = TRUE)



# 2019 05 03 sequencing run priors
priorsR1 = list(NS1_priorsDC,NS1_priorsDC,NS1_priorsDC,NS1_priorsDC,NS1_priorsDC,NS1_priorsDC,NS1_priorsDC,NS1_priorsDC,NS1_priorsDC,NS1_priorsDC,NS1_priorsDC,NS1_priorsDC,NS1_priorsDC,NS1_priorsDC,NS1_priorsDC,NS1_priorsDC,NS1_priorsDC,
                NS1_priors,NS1_priors,NS1_priors,NS1_priors,NS1_priors,NS1_priors,NS1_priors,NS1_priors
)

# orientF.R2 and orientR.R1 priors need to start at base #2; 3' end is trimmed 1 base in
priorsR2 = c(lapply(priorsR1, function(x) as.character(narrow(DNAStringSet(x), start = 1, end = width(DNAStringSet(x))-1))),
             lapply(priorsR1, function(x) as.character(reverseComplement(DNAStringSet(x)))))
names(priorsR2) = c(paste0(sample.names,".orientF"),paste0(sample.names,".orientR"))
priorsR1 = c(priorsR1, lapply(priorsR1, function(x) as.character(narrow(reverseComplement(DNAStringSet(x)), start = 2, end = width(reverseComplement(DNAStringSet(x)))))))
names(priorsR1) = c(paste0(sample.names,".orientF"),paste0(sample.names,".orientR"))

priorsR2 = llply(priorsR2, function(x){
  y = as.character(reverseComplement(DNAStringSet(x)))
  y = substr(y,trimLeftSelect[2]+1,truncLenSelect[2])
  y = y[laply(gregexpr("N", y),function(z) sum(z!=-1) == 0)]
  y = unique(y)
})

priorsR1 = llply(priorsR1, function(x){
  y = substr(x,trimLeftSelect[1]+1,truncLenSelect[1])
  y = y[laply(gregexpr("N", y),function(z) sum(z!=-1) == 0)]
  y = unique(y)
})

sample.names <- sapply(strsplit(basename(fnR1s), "_R"), `[`, 1)

process(path,
        fnR1s, fnR2s, sample.names, errPoolName = "20190503", orientFR.split = TRUE,
        trimLeftSelect, truncLenSelect, filter.matchIDs = FALSE,
        ErrModelMonotonicity = FALSE,
        OMEGA_A = getDadaOpt(option = "OMEGA_A"), pool = FALSE, poolList = NULL,
        priorsF = priorsR1, priorsR = priorsR2,
        
        plotToggle = FALSE,
        printPriorHead = FALSE,
        
        preFilteredToggle = TRUE,
        preCalcErrToggle = TRUE,
        preDerepToggle = TRUE,
        preDenoiseToggle = FALSE,
        
        saveFiltered = TRUE,
        saveDereplicated = TRUE,
        saveDenoised = TRUE)

#  ==== manual cross-run pseudopooling based on dada files (not run) ==== 
dadaFileList = list.files("./Alessandra_16s/NS0_NS1 analysis/dada denoised",full.names = TRUE)

dadaFileList[grep("NS0-F.*orientF_dadaF",dadaFileList, fixed = FALSE)]
dadaFileList[grep("NS0-F.*orientR_dadaF",dadaFileList, fixed = FALSE)]
dadaFileList[grep("NS0-F.*orientF_dadaR",dadaFileList, fixed = FALSE)]
dadaFileList[grep("NS0-F.*orientR_dadaR",dadaFileList, fixed = FALSE)]

dadaFileList[grep("NS1-F.*orientF_dadaF",dadaFileList, fixed = FALSE)]
dadaFileList[grep("NS1-F.*orientR_dadaF",dadaFileList, fixed = FALSE)]
dadaFileList[grep("NS1-F.*orientF_dadaR",dadaFileList, fixed = FALSE)]
dadaFileList[grep("NS1-F.*orientR_dadaR",dadaFileList, fixed = FALSE)]

dadaFileList[grep("NS0-DC.*orientF_dadaF",dadaFileList, fixed = FALSE)]
dadaFileList[grep("NS0-DC.*orientR_dadaF",dadaFileList, fixed = FALSE)]
dadaFileList[grep("NS0-DC.*orientF_dadaR",dadaFileList, fixed = FALSE)]
dadaFileList[grep("NS0-DC.*orientR_dadaR",dadaFileList, fixed = FALSE)]

dadaFileList[grep("NS1-DC.*orientF_dadaF",dadaFileList, fixed = FALSE)]
dadaFileList[grep("NS1-DC.*orientR_dadaF",dadaFileList, fixed = FALSE)]
dadaFileList[grep("NS1-DC.*orientF_dadaR",dadaFileList, fixed = FALSE)]
dadaFileList[grep("NS1-DC.*orientR_dadaR",dadaFileList, fixed = FALSE)]

dadaFileListSearch = c("NS0-F.*orientF_dadaF",
                        "NS0-F.*orientR_dadaF",
                        "NS0-F.*orientF_dadaR",
                        "NS0-F.*orientR_dadaR",
                        
                        "NS1-F.*orientF_dadaF",
                        "NS1-F.*orientR_dadaF",
                        "NS1-F.*orientF_dadaR",
                        "NS1-F.*orientR_dadaR",
                        
                        "NS0-DC.*orientF_dadaF",
                        "NS0-DC.*orientR_dadaF",
                        "NS0-DC.*orientF_dadaR",
                        "NS0-DC.*orientR_dadaR",
                        
                        "NS1-DC.*orientF_dadaF",
                        "NS1-DC.*orientR_dadaF",
                        "NS1-DC.*orientF_dadaR",
                        "NS1-DC.*orientR_dadaR")

test = alply(dadaFileListSearch,1, function(searchTerm){
  dadaFiles = alply(dadaFileList[grep(searchTerm,dadaFileList, fixed = FALSE)], 1, function(fileName){
    dadaFile = readRDS(fileName)
  })
  names(dadaFiles) = basename(dadaFileList[grep(searchTerm,dadaFileList, fixed = FALSE)])
  seqtab = makeSequenceTable(dadaFiles)
})

names(test) = c("NS0_F_C.orientF.R1",
                "NS0_F_C.orientR.R1",
                "NS0_F_C.orientF.R2",
                "NS0_F_C.orientR.R2",
                
                "NS1_F_C.orientF.R1",
                "NS1_F_C.orientR.R1",
                "NS1_F_C.orientF.R2",
                "NS1_F_C.orientR.R2",
                
                "NS0_DC.orientF.R1",
                "NS0_DC.orientR.R1",
                "NS0_DC.orientF.R2",
                "NS0_DC.orientR.R2",
                
                "NS1_DC.orientF.R1",
                "NS1_DC.orientR.R1",
                "NS1_DC.orientF.R2",
                "NS1_DC.orientR.R2")
  
llply(test, function(x){
  nchar(colnames(x))
})

# ==== POST-DADA2 ====
x20190905_seqtab.nochim = readRDS("./Alessandra_16s/NS0_NS1 analysis/outputs/20190905_seqtab.nochim.rds")
x20190503_seqtab.nochim = readRDS("./Alessandra_16s/NS0_NS1 analysis/outputs/20190503_seqtab.nochim.rds")

# x20190905_seqtab.nochim = readRDS("./Alessandra_16s/NS0_NS1 analysis/priorsOnly_outputs/20190905_seqtab.nochim.rds")
# x20190503_seqtab.nochim = readRDS("./Alessandra_16s/NS0_NS1 analysis/priorsOnly_outputs/20190503_seqtab.nochim.rds")

#rename seqtab.nochim row names
rownames(x20190905_seqtab.nochim) = gsub("-",".",rownames(x20190905_seqtab.nochim), fixed = TRUE)
rownames(x20190503_seqtab.nochim) = gsub("-",".",rownames(x20190503_seqtab.nochim), fixed = TRUE)

#merge Sequence Tables
seqtab = mergeSequenceTables(x20190905_seqtab.nochim, x20190503_seqtab.nochim)
seqtab.nochim = removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)

#isolates processing
isolateSeqData = read.delim(file.path(path,"input/isolates/isolateSeqData20221219.txt"), sep = "\t",stringsAsFactors = FALSE)

isolates_seqtab = matrix(data = as.numeric(
  c(NS0_F_C_IsolatesCleaned$source == "fecal",rep_len(0,nrow(NS1_F_C_IsolatesCleaned)),
    NS0_F_C_IsolatesCleaned$source == "chemostat",rep_len(0,nrow(NS1_F_C_IsolatesCleaned)),
    NS0_F_C_IsolatesCleaned$DC_inoculum,rep_len(0,nrow(NS1_F_C_IsolatesCleaned)),
    rep_len(0,nrow(NS0_F_C_IsolatesCleaned)),NS1_F_C_IsolatesCleaned$source == "fecal",
    rep_len(0,nrow(NS0_F_C_IsolatesCleaned)),NS1_F_C_IsolatesCleaned$source == "chemostat",
    rep_len(0,nrow(NS0_F_C_IsolatesCleaned)),NS1_F_C_IsolatesCleaned$DC_inoculum
  )
),
nrow = 6, ncol = length(NS0_F_C_IsolatesCleaned$rcV4Seq)+length(NS1_F_C_IsolatesCleaned$rcV4Seq), byrow = TRUE)


# isolates_seqtab = matrix(data = as.numeric(
#   c(NS0_F_C_IsolatesCleaned$source == "fecal",rep_len(0,nrow(NS1_F_C_IsolatesCleaned)),
#     NS0_F_C_IsolatesCleaned$source == "chemostat",rep_len(0,nrow(NS1_F_C_IsolatesCleaned)),
#     NS0_F_C_IsolatesCleaned$DC_inoculum,rep_len(0,nrow(NS1_F_C_IsolatesCleaned)),
#     rep_len(0,nrow(NS0_F_C_IsolatesCleaned)),NS1_F_C_IsolatesCleaned$source == "fecal",
#     rep_len(0,nrow(NS0_F_C_IsolatesCleaned)),NS1_F_C_IsolatesCleaned$source == "chemostat",
#     rep_len(0,nrow(NS0_F_C_IsolatesCleaned)),NS1_F_C_IsolatesCleaned$DC_inoculum
#   )
# ),
# nrow = 6, ncol = length(NS0_F_C_IsolatesCleaned$rcV4Seq)+length(NS1_F_C_IsolatesCleaned$rcV4Seq), byrow = TRUE)

colnames(isolates_seqtab) = c(NS0_F_C_IsolatesCleaned$rcV4Seq,NS1_F_C_IsolatesCleaned$rcV4Seq)
rownames(isolates_seqtab) = c("NS0.isoF","NS0.isoC","NS0.isoDC","NS1.isoF","NS1.isoC","NS1.isoDC")
isolates_seqtab = collapseNoMismatch2(isolates_seqtab,verbose = TRUE)

seqtab.nochim = mergeSequenceTables(seqtab.nochim,isolates_seqtab)
seqtab.nochim_collapse = collapseNoMismatch2(seqtab.nochim,verbose = TRUE)

#reorder samples
seqtab.nochim_collapse = seqtab.nochim_collapse[c("NS0.F",
                                                  "NS0.FCd01",
                                                  "NS0.FCd05",
                                                  "NS0.FCd09",
                                                  "NS0.FCd13",
                                                  "NS0.FCd17",
                                                  "NS0.FCd21",
                                                  "NS0.DCd00",
                                                  "NS0.DCV1d01",
                                                  "NS0.DCV1d05",
                                                  "NS0.DCV1d09",
                                                  "NS0.DCV1d13",
                                                  "NS0.DCV1d17",
                                                  "NS0.DCV1d21",
                                                  "NS0.DCV2d01",
                                                  "NS0.DCV2d05",
                                                  "NS0.DCV2d09",
                                                  "NS0.DCV2d13",
                                                  "NS0.DCV2d17",
                                                  "NS0.DCV2d21",
                                                  "NS1.F",
                                                  "NS1.FCd01",
                                                  "NS1.FCd05",
                                                  "NS1.FCd09",
                                                  "NS1.FCd13",
                                                  "NS1.FCd17",
                                                  "NS1.FCd21",
                                                  "NS1.FCd25",
                                                  "NS1.FCd27",
                                                  "NS1.DCd00",
                                                  "NS1.DCV1d01",
                                                  "NS1.DCV1d05",
                                                  "NS1.DCV1d09",
                                                  "NS1.DCV1d09R",
                                                  "NS1.DCV1d13",
                                                  "NS1.DCV1d17",
                                                  "NS1.DCV1d21",
                                                  "NS1.DCV2d01",
                                                  "NS1.DCV2d05",
                                                  "NS1.DCV2d09",
                                                  "NS1.DCV2d09R",
                                                  "NS1.DCV2d13",
                                                  "NS1.DCV2d17",
                                                  "NS1.DCV2d21",
                                                  "NS1.DCV3d01",
                                                  "NS1.DCV3d05",
                                                  "NS1.DCV3d09R",
                                                  "NS1.DCV3d13",
                                                  "NS1.DCV3d17",
                                                  "NS1.DCV3d21",
                                                  "NS0.isoF",
                                                  "NS0.isoC",
                                                  "NS0.isoDC",
                                                  "NS1.isoF",
                                                  "NS1.isoC",
                                                  "NS1.isoDC"),]

sampleNum = 50
sample_isoNum = 56

ASV_Table = data.frame(ASV_Num = seq_along(colnames(seqtab.nochim_collapse)), sequences = colnames(seqtab.nochim_collapse),t(seqtab.nochim_collapse), stringsAsFactors = FALSE)
row.names(ASV_Table) = NULL
# taxa = assignTaxonomy(seqtab.nochim, "./silva_nr_v132_train_set.fa.gz", multithread=TRUE)
# taxa = addSpecies(taxa, "./silva_species_assignment_v132.fa.gz")
# ASV_Table = data.frame(ASV_Table,taxa, stringsAsFactors = FALSE)
ASV_Table$ASV_Num = sprintf("%03d", ASV_Table$ASV_Num)
ASV_Table$ASV_Num = paste0("ASV_", ASV_Table$ASV_Num)
#ASV_Table = read.delim(file.path(path,"NS0_NS1 ASV_Table.txt"),sep = "\t", stringsAsFactors = FALSE)
#write.table(ASV_Table,file = file.path(path,"NS0_NS1 ASV_Table.txt"),sep = "\t", row.names = FALSE)
#write.table(ASV_Table,file = file.path(path,"NS0_NS1 ASV_TableNon_pooled.txt"),sep = "\t", row.names = FALSE)
# ==== Export seqs.fna (sequences), table.biom (counts), and metadata.tsv (metadata) for picrust2 ====
install.packages("castor")
library(castor)
#seqs.fna
seq2fasta(ASV_Table$ASV_Num, ASV_Table$sequences, "seqs", file.path(path,"picrust2_inputs"))
file.rename(file.path(path,"picrust2_inputs/seqs.txt"),file.path(path,"picrust2_inputs/seqs.fna"))
#table.biom
table.biom = ASV_Table[,12:(sampleNum+1)]
row.names(table.biom) = ASV_Table[,1]
colnames(table.biom) = colnames(ASV_Table)[12:(sampleNum+1)]
table.biom = make_biom(table.biom)
write_biom(table.biom, file.path(path,"picrust2_inputs/table.biom"))
#metadata.tsv
metadata.tsv = data.frame(sample_16S = colnames(ASV_Table)[12:(sampleNum+1)])
metadata.tsv = left_join(metadata.tsv,sample_Table)

write.table(metadata.tsv, file = file.path(path,"picrust2_inputs/metadata.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)

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

write.table(ASV_Table,file = file.path(path,"NS0_NS1 ASV_Table.txt"),sep = "\t", row.names = FALSE)

# ---- k-Align for pairwise identity matrix----
# idMat = Kalign_pim(userEmail = "christopher.yau@mail.utoronto.ca",
#                    fastafileName = paste0(path,"/ASV FASTAs/All_ASVs_isolates.txt"))
"./Alessandra_16s/NS0_NS1 analysis/ASV FASTAs/All_ASVs_isolates.txt_Kalign_out.txt"
idMat = read.delim(paste0(path,"/ASV FASTAs/All_ASVs_isolates.txt_Kalign_out.txt"),header = FALSE, stringsAsFactors = FALSE)

idMat = (100-idMat)/100
idMatDist = as.dist(idMat)

# ASVcluster = hclust(idMatDist, method = "single")
ASVcluster = hclust(idMatDist, method = "complete")
# ASVcluster = hclust(idMatDist, method = "average")
# ASVcluster = hclust(idMatDist, method = "ward.D")
# ASVcluster = hclust(idMatDist, method = "ward.D2")
#plot(ASVcluster,hang = -1)
ASVDendro = as.dendrogram(ASVcluster)
plot(ASVDendro)

#convert to within-sample relative abundance
ASV_TableRelAb = ASV_Table
ASV_TableRelAb[,12:(sample_isoNum+12-1)] = colwise(function(x) x/sum(x))(ASV_TableRelAb[,12:(sample_isoNum+12-1)])

# ---- Adding collapse metadata to ASV table and cleanup ---- 
# write.table(ASV_Table,file = paste0(path,"/NS0_NS1 ASV_Table.txt"),sep = "\t", row.names = FALSE)
ASV_TableSave = ASV_Table
ASV_Table = ASV_TableSave

#convert to within-sample relative abundance
ASV_TableRelAb = ASV_Table
ASV_TableRelAb[,12:(sample_isoNum+12-1)] = colwise(function(x) x/sum(x))(ASV_TableRelAb[,12:(sample_isoNum+12-1)])

# ---- Collapse to taxa (closest BLAST-match species)----
ASV_Table_TaxaCol = ddply(ASV_TableRelAb, .(BLASTspecies), function(x){
  countSum = colSums(x[,12:(sample_isoNum+12-1)])
  df = c(ASV_Num = x$ASV_Num[1], countSum, x[1,(sample_isoNum+2):length(x)])
  df = data.frame(df, stringsAsFactors = FALSE)
  return(df)
})

# ---- Collapse by sequence identity----
ASV_Table_idCol = ddply(ASV_Table, .(idCluster), function(x){
  if(nrow(x)>1){
    countSum = colSums(x[,12:(sample_isoNum+12-1)])
    clusterIdMat = idMat[as.character(x$ASV_Num),as.character(x$ASV_Num)]
    minDistSums = which.min(colSums(clusterIdMat))
    df = c(ASV_Num = x$ASV_Num[minDistSums], countSum, x[minDistSums,(sample_isoNum+2):length(x)])
    df = data.frame(df, clustMemNum = nrow(x), stringsAsFactors = FALSE)
    return(df)
  }else{
    df = data.frame(x, clustMemNum = 1, stringsAsFactors = FALSE)
    return(df)
  }
})
ASV_TableRelAb_idCol = ddply(ASV_TableRelAb, .(idCluster), function(x){
  if(nrow(x)>1){
    countSum = colSums(x[,12:(sample_isoNum+12-1)])
    clusterIdMat = idMat[as.character(x$ASV_Num),as.character(x$ASV_Num)]
    minDistSums = which.min(colSums(clusterIdMat))
    df = c(ASV_Num = x$ASV_Num[minDistSums], countSum, x[minDistSums,(sample_isoNum+2):length(x)])
    df = data.frame(df, clustMemNum = nrow(x), stringsAsFactors = FALSE)
    return(df)
  }else{
    df = data.frame(x, clustMemNum = 1, stringsAsFactors = FALSE)
    return(df)
  }
})

# ---- Sample Table setup for plots ----
sample_Table = data.frame(sample_16S = colnames(ASV_Table)[12:(sampleNum+1)], stringsAsFactors = FALSE)
sample_Table = data.frame(sample_Table,
                          source = substr(sample_Table$sample_16S,1,3), 
                          type = substr(sample_Table$sample_16S,5,6), 
                          vessel = sapply(strsplit(sample_Table$sample_16S ,"V|(d.*)", fixed = FALSE), function(x) x[2]), 
                          day = as.numeric(sapply(strsplit(sample_Table$sample_16S ,"d|R", fixed = FALSE), function(x) x[2])))
sample_Table$day[is.na(sample_Table$day)] = 0
# ---- ASV Table setup for plots ----
#select collapse form
# graphData = ASV_TableSave
graphData = ASV_TableRelAb
# graphData = ASV_Table_TaxaCol
# graphData = ASV_Table_idCol
graphData = ASV_TableRelAb_idCol
