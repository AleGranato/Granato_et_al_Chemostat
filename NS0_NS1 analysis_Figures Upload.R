library(BiocManager)
# install.packages("BiocManager")
# BiocManager::install(version = "3.16")
# BiocManager::install("msa")

library(ComplexHeatmap)
library(scales)
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
library(clValid)
library(gtable)
library(cowplot)
library(httr)
library(ShortRead)
library(ggdendro)
library(dendextend)
library(DECIPHER)

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
          print("Filtering for top hits...")
          blastResults = blastResults[order(blastResults$evalue),]
          blastResults = blastResults[!duplicated(blastResults$query),]
          blastResults = blastResults[order(blastResults$query),]
          
          # if taxonomy, query NCBI Taxonomy for taxonomic information
          if(taxonomy){
            print("Querying NCBI Taxonomy for taxonomic information...")
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
  kalignJobID = httr::content(kalignRunResponse)
  
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
    kalignJobStatusResponse = httr::content(kalignJobStatus)
  }
  
  # GET response
  kalignJob_pimResponse = GET(paste0("https://www.ebi.ac.uk/Tools/services/rest/kalign/result/",kalignJobID,"/pim"))
  kalignJob_pim = httr::content(kalignJob_pimResponse)
  
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

# ==== DECIPHER_clusterASVs function ====
library("DECIPHER")
DECIPHER_clusterASVs = function(seqs,
                                identityCutoff = 1,
                                testBounds = FALSE, upperBound = 1, lowerBound = 1, testIncrement = 0.1,
                                ncores = NULL){
  
  # detect cores; uses all cores if ncores = NULL
  if(is.null(ncores)){
    ncores <- detectCores()
  }
  
  dna <- Biostrings::DNAStringSet(seqs)
  
  # align sequences and calculate distance matrix
  cat("Aligning sequences...\n")
  aln <- DECIPHER::AlignSeqs(dna, processors = ncores)
  cat("Calculating distance matrix...\n")
  dmat <- DECIPHER::DistanceMatrix(aln, processors = ncores)
  
  if(!testBounds){
    cat(paste0("Clustering at ",identityCutoff * 100,"% identity...\n"))
    clusters <- DECIPHER::TreeLine(
      myDistMatrix = dmat,
      method = "complete",
      cutoff = 1-identityCutoff, # use cutoff = 0.03 for a 97% OTU
      type = "clusters",
      processors = ncores)
    
  }else{
    cat(paste0("Clustering between",lowerBound*100,"-",upperBound*100,"% identity, incrementing by ",testIncrement*100,"%...\n"))
    clusters <- DECIPHER::TreeLine(
      myDistMatrix=dmat,
      method = "complete",
      cutoff = 1-seq(lowerBound, upperBound, by = testIncrement),
      type = "clusters",
      processors = ncores)
  }
  return(clusters)
}

# ==== Loading ASV table ====
setwd("~/Desktop/Christopher_Yau/R Data Analysis/Dada2/Alessandra_16s/NS0_NS1 analysis")
path <- "."

sampleNum = 50
sample_isoNum = 56

# load pseudo-pooled ASV table
ASV_Table = read.delim(paste0(path,"/NS0_NS1 ASV_Table.txt"),sep = "\t", stringsAsFactors = FALSE)

# test = ASV_Table[,12:61]
# sum(rowSums(test)>0)
# sort(colSums(test))

#convert to within-sample relative abundance
ASV_TableRelAb = ASV_Table
ASV_TableRelAb[,12:(sample_isoNum+12-1)] = colwise(function(x) x/sum(x))(ASV_TableRelAb[,12:(sample_isoNum+12-1)])
# ASV_TableRelAbSave = ASV_TableRelAb
# write.table(ASV_TableRelAb,file = paste0(path,"/NS0_NS1 ASV_TableRelAb.txt"),sep = "\t", row.names = FALSE)

#convert to within-sample relative abundance
ASV_TableRelAb = ASV_Table
ASV_TableRelAb[,12:(sample_isoNum+12-1)] = colwise(function(x) x/sum(x))(ASV_TableRelAb[,12:(sample_isoNum+12-1)])


# load non-pooled ASV table
ASV_TableNon_pooled = read.delim(paste0(path,"/NS0_NS1 ASV_TableNon_pooled.txt"),sep = "\t", stringsAsFactors = FALSE)
#convert to within-sample relative abundance
ASV_TableNon_pooledRelAb = ASV_TableNon_pooled
ASV_TableNon_pooledRelAb[,3:(sample_isoNum+3-1)] = colwise(function(x) x/sum(x))(ASV_TableNon_pooledRelAb[,3:(sample_isoNum+3-1)])


# load full-pooled ASV table
ASV_TableFull_pooled = read.delim(paste0(path,"/NS0_NS1 ASV_TableFull_pooled.txt"),sep = "\t", stringsAsFactors = FALSE)
#convert to within-sample relative abundance
ASV_TableFull_pooledRelAb = ASV_TableFull_pooled
ASV_TableFull_pooledRelAb[,3:(sample_isoNum+3-1)] = colwise(function(x) x/sum(x))(ASV_TableFull_pooledRelAb[,3:(sample_isoNum+3-1)])

# load pseudo-pooled 7 participant ASV table
ASV_Table_7par = read.delim(paste0(path,"/NS0_NS1 ASV_Table_reseq.txt"),sep = "\t", stringsAsFactors = FALSE)
#convert to within-sample relative abundance
ASV_Table_7parRelAb = ASV_Table_7par
ASV_Table_7parRelAb[,12:ncol(ASV_Table_7parRelAb)] = colwise(function(x) x/sum(x))(ASV_Table_7parRelAb[,12:ncol(ASV_Table_7parRelAb)])


# load pseudo-pooled resequenced ASV table for MS co-analysis
ASV_Table_reseq = read.delim(paste0(path,"/DIABIMMUNE Chemostat MS_reseq ASV_Table.txt"),sep = "\t", stringsAsFactors = FALSE)
# keep NS0 and NS1 samples
ASV_Table_reseq = cbind(ASV_Table_reseq[,1:11],ASV_Table_reseq[,grep("NS0|NS1",colnames(ASV_Table_reseq), fixed = FALSE)])
#convert to within-sample relative abundance
ASV_Table_reseqRelAb = ASV_Table_reseq
ASV_Table_reseqRelAb[,12:ncol(ASV_Table_reseqRelAb)] = colwise(function(x) x/sum(x))(ASV_Table_reseqRelAb[,12:ncol(ASV_Table_reseqRelAb)])

# ==== Loading isolates table ====
path <- "./"

# load NS0 isolates table
NS0_F_C_IsolatesCleaned = read.delim(file.path(path,"input/isolates/NS0 F_C_IsolatesCleaned.txt"),sep = "\t", stringsAsFactors = FALSE)

# load NS1 isolates table
NS1_F_C_IsolatesCleaned = read.delim(file.path(path,"input/isolates/NS1 F_C_IsolatesCleaned.txt"),sep = "\t", stringsAsFactors = FALSE)

# ==== load NMR metabolites data ====
NMRmetabolite_Table = read.delim(file = "./input/NS0_NS1_NMRmetabolites_media_names.txt", stringsAsFactors = FALSE, check.names = FALSE)
colnames(NMRmetabolite_Table)[6] = "day_NMR"
NMRmetabolite_Table_Media = cbind.data.frame(NMRmetabolite_Table[,1:6],scale(NMRmetabolite_Table[,7:ncol(NMRmetabolite_Table)], center = NMRmetabolite_Table[1,7:ncol(NMRmetabolite_Table)], scale = FALSE))
# NMRmetabolite_Table_Media = NMRmetabolite_Table
NMRmetaboliteMelt = melt(NMRmetabolite_Table,measure.vars = 7:ncol(NMRmetabolite_Table), variable.name = "metabolite", value.name = "conc", stringsAsFactors = FALSE)
NMRmetaboliteMelt$vessel[is.na(NMRmetaboliteMelt$vessel)] = 1

NMRmetabolite_MediaMelt = melt(NMRmetabolite_Table_Media,measure.vars = 7:ncol(NMRmetabolite_Table_Media), variable.name = "metabolite", value.name = "conc", stringsAsFactors = FALSE)
NMRmetabolite_MediaMelt$vessel[is.na(NMRmetabolite_MediaMelt$vessel)] = 1

# ==== load MS metabolites data ====
# FC MS data
MSmetabolite_TableFC_clean = read.delim(file = "./input/NS0_NS1_MSmetabolitesFC_clean_names.txt", stringsAsFactors = FALSE, check.names = FALSE)
MSmetabolite_TableFC_clean_Media = cbind.data.frame(MSmetabolite_TableFC_clean[,1:7],scale(MSmetabolite_TableFC_clean[,8:ncol(MSmetabolite_TableFC_clean)], center = colMeans(MSmetabolite_TableFC_clean[1:4,8:ncol(MSmetabolite_TableFC_clean)]), scale = FALSE))

# MSmetabolite_TableFC_clean_Zscore = MSmetabolite_TableFC_clean[!is.na(MSmetabolite_TableFC_clean$type),]
# MSmetabolite_TableFC_clean_Zscore = cbind.data.frame(MSmetabolite_TableFC_clean_Zscore[,1:7],scale(MSmetabolite_TableFC_clean_Zscore[,8:ncol(MSmetabolite_TableFC_clean_Zscore)], center = TRUE, scale = TRUE))
# 
# MSmetabolite_TableFC_clean_MediaFold = cbind.data.frame(MSmetabolite_TableFC_clean[,1:7],sweep(MSmetabolite_TableFC_clean[,8:ncol(MSmetabolite_TableFC_clean)],2,colMeans(MSmetabolite_TableFC_clean[1:4,8:ncol(MSmetabolite_TableFC_clean)]),'/'))

# NMRmetabolite_Table_Media = NMRmetabolite_Table
MSmetabolite_TableFC_cleanMelt = MSmetabolite_TableFC_clean
MSmetabolite_TableFC_cleanMelt = MSmetabolite_TableFC_cleanMelt[!is.na(MSmetabolite_TableFC_cleanMelt$source),]
MSmetabolite_TableFC_cleanMelt = MSmetabolite_TableFC_cleanMelt[,c(1:7,(8:ncol(MSmetabolite_TableFC_cleanMelt))[colSums(MSmetabolite_TableFC_cleanMelt[,8:ncol(MSmetabolite_TableFC_cleanMelt)])!=0])]
MSmetabolite_TableFC_cleanMelt = melt(MSmetabolite_TableFC_cleanMelt,measure.vars = 8:ncol(MSmetabolite_TableFC_cleanMelt), variable.name = "metabolite", value.name = "conc", stringsAsFactors = FALSE)
MSmetabolite_TableFC_cleanMelt$vessel[is.na(MSmetabolite_TableFC_cleanMelt$vessel)] = 1

MSmetabolite_TableFC_clean_MediaMelt = melt(MSmetabolite_TableFC_clean_Media,measure.vars = 8:ncol(NMRmetabolite_Table_Media), variable.name = "metabolite", value.name = "conc", stringsAsFactors = FALSE)
MSmetabolite_TableFC_clean_MediaMelt$vessel[is.na(MSmetabolite_TableFC_clean_MediaMelt$vessel)] = 1

# MSmetabolite_TableFC_clean_ZscoreMelt = melt(MSmetabolite_TableFC_clean_Zscore,measure.vars = 8:ncol(MSmetabolite_TableFC_clean_Zscore), variable.name = "metabolite", value.name = "conc", stringsAsFactors = FALSE)
# MSmetabolite_TableFC_clean_ZscoreMelt$vessel[is.na(MSmetabolite_TableFC_clean_ZscoreMelt$vessel)] = 1

# DC MS data
MSmetabolite_TableDC_clean = read.delim(file = "./input/NS0_NS1_MSmetabolitesDC_clean_names.txt", stringsAsFactors = FALSE, check.names = FALSE)
MSmetabolite_TableDC_clean_Media = cbind.data.frame(MSmetabolite_TableDC_clean[,1:7],scale(MSmetabolite_TableDC_clean[,8:ncol(MSmetabolite_TableDC_clean)], center = colMeans(MSmetabolite_TableDC_clean[1:4,8:ncol(MSmetabolite_TableDC_clean)]), scale = FALSE))

# MSmetabolite_TableDC_clean_Zscore = MSmetabolite_TableDC_clean[!is.na(MSmetabolite_TableDC_clean$type),]
# MSmetabolite_TableDC_clean_Zscore = cbind.data.frame(MSmetabolite_TableDC_clean_Zscore[,1:7],scale(MSmetabolite_TableDC_clean_Zscore[,8:ncol(MSmetabolite_TableDC_clean_Zscore)], center = TRUE, scale = TRUE))
# 
# MSmetabolite_TableDC_clean_MediaFold = cbind.data.frame(MSmetabolite_TableDC_clean[,1:7],sweep(MSmetabolite_TableDC_clean[,8:ncol(MSmetabolite_TableDC_clean)],2,colMeans(MSmetabolite_TableDC_clean[1:4,8:ncol(MSmetabolite_TableDC_clean)]),'/'))

# NMRmetabolite_Table_Media = NMRmetabolite_Table
MSmetabolite_TableDC_cleanMelt = melt(MSmetabolite_TableDC_clean,measure.vars = 8:ncol(MSmetabolite_TableDC_clean), variable.name = "metabolite", value.name = "conc", stringsAsFactors = FALSE)
MSmetabolite_TableDC_cleanMelt$vessel[is.na(MSmetabolite_TableDC_cleanMelt$vessel)] = 1

MSmetabolite_TableDC_clean_MediaMelt = melt(MSmetabolite_TableDC_clean_Media,measure.vars = 8:ncol(NMRmetabolite_Table_Media), variable.name = "metabolite", value.name = "conc", stringsAsFactors = FALSE)
MSmetabolite_TableDC_clean_MediaMelt$vessel[is.na(MSmetabolite_TableDC_clean_MediaMelt$vessel)] = 1

# MSmetabolite_TableDC_clean_ZscoreMelt = melt(MSmetabolite_TableDC_clean_Zscore,measure.vars = 8:ncol(MSmetabolite_TableDC_clean_Zscore), variable.name = "metabolite", value.name = "conc", stringsAsFactors = FALSE)
# MSmetabolite_TableDC_clean_ZscoreMelt$vessel[is.na(MSmetabolite_TableDC_clean_ZscoreMelt$vessel)] = 1


# ---- Sample Table setup for plots ----
sample_Table = data.frame(sample_16S = c(colnames(ASV_Table)[12:(sampleNum+12-1)],colnames(ASV_Table_7par)[12:ncol(ASV_Table_7par)]), stringsAsFactors = FALSE)
sample_Table = sample_Table[!duplicated(sample_Table$sample_16S), ,drop = FALSE]

sample_Table = full_join(sample_Table,NMRmetabolite_Table[,c(1:6)])

# sample_Table$source[is.na(sample_Table$source)] = substr(sample_Table$sample_16S,1,3)[is.na(sample_Table$source)]
sample_Table$source[is.na(sample_Table$source)] = sapply(strsplit(sample_Table$sample_16S,"\\."), function(x) x[1])[is.na(sample_Table$source)]

# sample_Table$type[is.na(sample_Table$type)] = substr(sample_Table$sample_16S,5,6)[is.na(sample_Table$type)]
sample_Table$type[is.na(sample_Table$type)] = substr(sapply(strsplit(sample_Table$sample_16S,"\\."), function(x) x[2]),1,2)[is.na(sample_Table$type)]

sample_Table$vessel[is.na(sample_Table$vessel)] = sapply(strsplit(sample_Table$sample_16S ,"V|(d.*)", fixed = FALSE), function(x) x[2])[is.na(sample_Table$vessel)]
sample_Table$vessel[is.na(sample_Table$vessel)&!is.na(sample_Table$sample_16S)] = 0
sample_Table$vessel[grep("DCd00",sample_Table$sample_16S)] = "inoculum"

sample_Table$day_16S = as.numeric(sapply(strsplit(sample_Table$sample_16S ,"d|R", fixed = FALSE), function(x) x[2]))
sample_Table$day_16S[is.na(sample_Table$day_16S)&!is.na(sample_Table$sample_16S)] = 0
sample_Table = sample_Table[,c(1:5,7,6)]



# ==== Figure uniformity setup ====
# Set uniform phylum colors
# Chose a RcolorBrewer palette
library(RColorBrewer)

# Accent
# Dark2
# Paired
# Pastel1
# Pastel2
# Set1
# Set2
# Set3

phylumColors = brewer.pal(9, "Set1")
names(phylumColors) = c(
  "Actinomycetota",
  "Bacillota",
  "Bacteroidota",
  "Campylobacterota",
  "Fusobacteriota",
  "Lentisphaerota",
  "Pseudomonadota",
  "Thermodesulfobacteriota",
  "Verrucomicrobiota"
)

show_col(phylumColors, ncol = 3)

# hue_pal()(7) ggplot default
# phylumColors = c(
#   "Actinobacteria" = "#F8766D",
#   "Bacteroidetes" = "#C49A00",
#   "Firmicutes" = "#53B400",
#   "Fusobacteria" = "#00C094",
#   "Lentisphaerae" = "#00B6EB",
#   "Proteobacteria" = "#A58AFF",
#   "Verrucomicrobia" = "#FB61D7"
# )

# ---- absencePresenceLevels setup ----
absencePresenceLevels2 = c(
  # "FALSE.FALSE",
  "TRUE.FALSE",
  "FALSE.TRUE",
  "TRUE.TRUE"
)

absencePresence2LegendData = data.frame(levelName = absencePresenceLevels2,
                                        y = c(2,1,0),
                                        x1 = c(1,0,1),
                                        x2 = c(0,1,1))

absencePresenceLevels2_0 = c(
  "FALSE.FALSE",
  "TRUE.FALSE",
  "FALSE.TRUE",
  "TRUE.TRUE"
)

absencePresence2_0LegendData = data.frame(levelName = absencePresenceLevels2_0,
                                        y = c(3,2,1,0),
                                        x1 = c(0,1,0,1),
                                        x2 = c(0,0,1,1))

absencePresenceLevels3 = c(
  # "FALSE.FALSE.FALSE",
  "TRUE.FALSE.FALSE",
  "FALSE.TRUE.FALSE",
  "FALSE.FALSE.TRUE",
  "TRUE.TRUE.FALSE",
  "TRUE.FALSE.TRUE",
  "FALSE.TRUE.TRUE",
  "TRUE.TRUE.TRUE"
)

absencePresence3LegendData = data.frame(levelName = absencePresenceLevels3,
                                        y = c(6,5,4,3,2,1,0),
                                        x1 = c(1,0,0,1,1,0,1),
                                        x2 = c(0,1,0,1,0,1,1),
                                        x3 = c(0,0,1,0,1,1,1))

# ---- Themes and sample names ----
titleSize = 26
titleSize.axis = 30
legendSize = 28
textSize = 28 #25
textSize.x = 28
textSize.y = 30
textSize.point = 16
textSize.point.legend = 25
textSize.point.select = 34
titleSize.axis.select = 36
# textSize.point2 = 14
legendTextSize.fecal2 = 16
legendTextSize.chemostat2 = 13
legendTextSize.chemostat3 = 20
legendTextSize.day2 = 30
legendTextSize.day3 = 16

alluvial.theme_x0 = theme_bw()+
  theme(axis.title.x = element_blank(),
        # axis.title.x = element_text(size=titleSize, face = "bold", family="Helvetica"),
        axis.title.y = element_text(size=titleSize.axis, face = "bold", family="Helvetica"),
        axis.text.x = element_text(size=textSize.x, face = "bold", family="Helvetica"),
        axis.text.y = element_text(size=textSize.y, face = "bold", family="Helvetica"),
        title = element_text(size=textSize, face = "bold", family="Helvetica"),
        strip.text.x = element_text(size=titleSize, face = "bold", family="Helvetica"),
        strip.text.y = element_text(size=titleSize, face = "bold", family="Helvetica"),
        legend.title = element_text(size=legendSize, face = "bold", family="Helvetica"),
        legend.text = element_text(size=textSize, face = "bold", family="Helvetica"),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())

alluvial.theme_x1 = theme_bw()+
  theme(# axis.title.x = element_blank(),
        axis.title.x = element_text(size=titleSize.axis, face = "bold", family="Helvetica"),
        axis.title.y = element_text(size=titleSize.axis, face = "bold", family="Helvetica"),
        axis.text.x = element_text(size=textSize.x, face = "bold", family="Helvetica"),
        axis.text.y = element_text(size=textSize.y, face = "bold", family="Helvetica"),
        title = element_text(size=textSize, face = "bold", family="Helvetica"),
        strip.text.x = element_text(size=titleSize, face = "bold", family="Helvetica"),
        strip.text.y = element_text(size=titleSize, face = "bold", family="Helvetica"),
        legend.title = element_text(size=legendSize, face = "bold", family="Helvetica"),
        legend.text = element_text(size=textSize, face = "bold", family="Helvetica"),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())

alluvial.theme_x0_legend0 = alluvial.theme_x0+
  theme(legend.position = "none")

alluvial.theme_x1_legend0 = alluvial.theme_x1+
  theme(legend.position = "none")

legend.fecal2.theme = theme_void()+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=legendTextSize.fecal2, face = "bold", family="Helvetica", vjust = 0),
        plot.title = element_text(size=legendSize, face = "bold", family="Helvetica", hjust = 0.66)
  )

legend.chemostat2.theme = theme_void()+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=legendTextSize.chemostat2, face = "bold", family="Helvetica", vjust = 0),
        plot.title = element_text(size=legendSize, face = "bold", family="Helvetica", hjust = 0.66)
  )

legend.chemostat3.theme = theme_void()+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=legendTextSize.chemostat3, face = "bold", family="Helvetica", vjust = 0),
        plot.title = element_text(size=legendSize, face = "bold", family="Helvetica", hjust = 0.5)
  )

legend.day2.theme = theme_void()+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=legendTextSize.day2, face = "bold", family="Helvetica", vjust = 0),
        plot.title = element_text(size=legendSize, face = "bold", family="Helvetica", hjust = 0.5)
  )

legend.day3.theme = theme_void()+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=legendTextSize.day3, face = "bold", family="Helvetica", vjust = 0),
        plot.title = element_text(size=legendSize, face = "bold", family="Helvetica", hjust = 0.66)
  )

PCoA.theme = theme_bw()+
  theme(axis.title.y = element_text(size=titleSize.axis, face = "bold", family="Helvetica"),
        axis.title.x = element_text(size=titleSize.axis, face = "bold", family="Helvetica"),
        axis.text.x = element_text(size=textSize.x, face = "bold", family="Helvetica"),
        axis.text.y = element_text(size=textSize.y, face = "bold", family="Helvetica"),
        title = element_text(size=textSize, face = "bold", family="Helvetica"),
        strip.text.x = element_text(size=titleSize, face = "bold", family="Helvetica"),
        strip.text.y = element_text(size=titleSize, face = "bold", family="Helvetica"),
        legend.title = element_text(size=legendSize, face = "bold", family="Helvetica"),
        legend.text = element_text(size=textSize, face = "bold", family="Helvetica"),
        plot.title = element_text(hjust = 0.5))

line.theme = theme_bw()+
  theme(axis.title.y = element_text(size=titleSize.axis, face = "bold", family="Helvetica"),
        axis.title.x = element_text(size=titleSize.axis, face = "bold", family="Helvetica"),
        axis.text.x = element_text(size=textSize.x, face = "bold", family="Helvetica"),
        axis.text.y = element_text(size=textSize.y, face = "bold", family="Helvetica"),
        title = element_text(size=textSize, face = "bold", family="Helvetica"),
        strip.text.x = element_text(size=titleSize, face = "bold", family="Helvetica"),
        strip.text.y = element_text(size=titleSize, face = "bold", family="Helvetica"),
        legend.title = element_text(size=legendSize, face = "bold", family="Helvetica"),
        legend.text = element_text(size=textSize, face = "bold", family="Helvetica"),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())

point.theme = theme_bw()+
  theme(axis.title.y = element_text(size=titleSize.axis, face = "bold", family="Helvetica"),
        axis.title.x = element_text(size=titleSize.axis, face = "bold", family="Helvetica"),
        axis.text.x = element_text(size=textSize.point, face = "bold", family="Helvetica"),
        axis.text.y = element_text(size=textSize.point, face = "bold", family="Helvetica"),
        title = element_text(size=textSize, face = "bold", family="Helvetica"),
        strip.text.x = element_text(size=titleSize, face = "bold", family="Helvetica"),
        strip.text.y = element_text(size=titleSize, face = "bold", family="Helvetica"),
        legend.title = element_text(size=legendSize, face = "bold", family="Helvetica"),
        legend.text = element_text(size=textSize.point.legend, face = "bold", family="Helvetica"),
        plot.title = element_text(hjust = 0.5))

point.theme_legendtl = point.theme +
  theme(legend.title = element_blank(),
        legend.position = c(0.05, 0.95),
        legend.justification = c(0, 1),
        legend.background = element_blank(),
        legend.key = element_blank())

point.theme_legendtr = point.theme +
  theme(legend.title = element_blank(),
        legend.position = c(0.95, 0.95),
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        legend.key = element_blank())


point.select.theme = theme_bw()+
  theme(axis.title.y = element_text(size=titleSize.axis.select, face = "bold", family="Helvetica"),
        axis.title.x = element_text(size=titleSize.axis.select, face = "bold", family="Helvetica"),
        axis.text.x = element_text(size=textSize.point.select, face = "bold", family="Helvetica"),
        axis.text.y = element_text(size=textSize.point.select, face = "bold", family="Helvetica"),
        title = element_text(size=textSize, face = "bold", family="Helvetica"),
        strip.text.x = element_text(size=titleSize, face = "bold", family="Helvetica"),
        strip.text.y = element_text(size=titleSize, face = "bold", family="Helvetica"),
        legend.title = element_text(size=legendSize, face = "bold", family="Helvetica"),
        legend.text = element_text(size=textSize.point.select, face = "bold", family="Helvetica"),
        plot.title = element_text(hjust = 0.5))

point.select.theme_legendtl = point.select.theme +
  theme(legend.title = element_blank(),
        legend.position = c(0.05, 0.95),
        legend.justification = c(0, 1),
        legend.background = element_blank(),
        legend.key = element_blank())

point.select.theme_legendtr = point.select.theme +
  theme(legend.title = element_blank(),
        legend.position = c(0.95, 0.95),
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        legend.key = element_blank())
  
histogram.theme = theme_bw()+
  theme(# axis.title.x = element_blank(),
    axis.title.x = element_text(size=titleSize.axis, face = "bold", family="Helvetica"),
    axis.title.y = element_text(size=titleSize.axis, face = "bold", family="Helvetica"),
    axis.text.x = element_text(size=textSize.x, face = "bold", family="Helvetica"),
    axis.text.y = element_text(size=textSize.y, face = "bold", family="Helvetica"),
    title = element_text(size=textSize, face = "bold", family="Helvetica"),
    strip.text.x = element_text(size=titleSize, face = "bold", family="Helvetica"),
    strip.text.y = element_text(size=titleSize, face = "bold", family="Helvetica"),
    legend.title = element_text(size=legendSize, face = "bold", family="Helvetica"),
    legend.text = element_text(size=textSize, face = "bold", family="Helvetica"),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_blank(), 
    panel.grid.minor.x = element_blank(),
    # legend.title = element_blank(),
    legend.position = c(0.95, 0.95),
    legend.justification = c(1, 1),
    legend.background = element_blank(),
    legend.key = element_blank())

# point.theme2 = theme_bw()+
#   theme(axis.title.y = element_text(size=titleSize, face = "bold", family="Helvetica"),
#         axis.title.x = element_text(size=titleSize, face = "bold", family="Helvetica"),
#         axis.text.x = element_text(size=textSize.point2, face = "bold", family="Helvetica"),
#         axis.text.y = element_text(size=textSize.point2, face = "bold", family="Helvetica"),
#         title = element_text(size=textSize, face = "bold", family="Helvetica"),
#         strip.text.x = element_text(size=titleSize, face = "bold", family="Helvetica"),
#         strip.text.y = element_text(size=titleSize, face = "bold", family="Helvetica"),
#         legend.title = element_text(size=legendSize, face = "bold", family="Helvetica"),
#         legend.text = element_text(size=textSize.point2, face = "bold", family="Helvetica"),
#         plot.title = element_text(hjust = 0.5))

sampleNameLevelsOrig = colnames(ASV_TableRelAb)[12:ncol(ASV_TableRelAb)]
sampleNameLevels = colnames(ASV_TableRelAb)[12:ncol(ASV_TableRelAb)]
sampleNameLevels = gsub(pattern = ".F",
            replacement = "\nFecal",
            x = sampleNameLevels,
            fixed = TRUE)
sampleNameLevels = gsub(pattern = "d",
            replacement = "\nday ",
            x = sampleNameLevels,
            fixed = FALSE)
sampleNameLevels = gsub(pattern = ".DC",
            replacement = "\nDefined Community",
            x = sampleNameLevels,
            fixed = TRUE)
sampleNameLevels = gsub(pattern = "V",
            replacement = "\nvessel ",
            x = sampleNameLevels,
            fixed = TRUE)

sampleNameLevels.3lines = sampleNameLevels
sampleNameLevels.3lines = gsub(pattern = "lC",
                        replacement = "l\nChemostat",
                        x = sampleNameLevels.3lines,
                        fixed = TRUE)

sampleNameLevels.3lines = gsub(pattern = "\nday 17",
                        replacement = "",
                        x = sampleNameLevels.3lines,
                        fixed = TRUE)

sampleNameLevels.3lines = gsub(pattern = "NS0\nDefined Community\nday 00",
                               replacement = "NS0\nIsolate\nLibrary",
                               x = sampleNameLevels.3lines,
                               fixed = TRUE)

sampleNameLevels.3lines = gsub(pattern = "NS1\nDefined Community\nday 00",
                        replacement = "NS1\nIsolate\nLibrary",
                        x = sampleNameLevels.3lines,
                        fixed = TRUE)

names(sampleNameLevels.3lines) = sampleNameLevelsOrig


sampleNameLevels = gsub(pattern = "lC",
            replacement = "l Chemostat",
            x = sampleNameLevels,
            fixed = TRUE)

sampleNameLevels = gsub(pattern = "\nday 17",
                        replacement = "",
                        x = sampleNameLevels,
                        fixed = TRUE)

sampleNameLevels = gsub(pattern = "NS0\nDefined Community\nday 00",
                        replacement = "NS0\nIsolate Library",
                        x = sampleNameLevels,
                        fixed = TRUE)

sampleNameLevels = gsub(pattern = "NS1\nDefined Community\nday 00",
                        replacement = "NS1\nIsolate Library",
                        x = sampleNameLevels,
                        fixed = TRUE)

names(sampleNameLevels) = sampleNameLevelsOrig

# #### Export setup ####
# newDir = gsub("-","",Sys.Date(), fixed = TRUE)
# dir.create(paste0(path,"Figure Outputs/",newDir,"/"))
newDir = "20240517"

# #### Tables ####
# ---- Tables for reference ----
# ---- sample alpha diversity table ----
# select non-pooled samples
ASV_TableSelectNonPooled = ASV_TableNon_pooled[,c("NS0.F","NS1.F","NS0.FCd17","NS1.FCd17")]
rownames(ASV_TableSelectNonPooled) = ASV_TableNon_pooled$sequences
colnames(ASV_TableSelectNonPooled) = paste0("nonpooled_",colnames(ASV_TableSelectNonPooled))
ASV_TableSelectNonPooled = t(ASV_TableSelectNonPooled)

ASV_TableSeqTab = ASV_Table[12:ncol(ASV_Table)]
rownames(ASV_TableSeqTab) = ASV_Table$sequences
ASV_TableSeqTab = t(ASV_TableSeqTab)

ASV_TableAlphaDiversityMerged = mergeSequenceTables(ASV_TableSelectNonPooled, ASV_TableSeqTab)
ASV_TableAlphaDiversityMerged_collapse = collapseNoMismatch2(ASV_TableAlphaDiversityMerged,verbose = TRUE)

alphaDiversityTable = data.frame(sample = rownames(ASV_TableAlphaDiversityMerged_collapse),
                                 ASV_count = rowSums(ASV_TableAlphaDiversityMerged_collapse >0),
                                 chao1 = t(estimateR(x = ASV_TableAlphaDiversityMerged_collapse))[,"S.chao1"],
                                 shannon = diversity(x = ASV_TableAlphaDiversityMerged_collapse,
                                                     index = "shannon"))

write.table(alphaDiversityTable, file = paste0(path,"Figure Outputs/",newDir,"/alphaDiversityTable.txt"), 
            quote = FALSE,
            sep = "\t",
            row.names =  FALSE, 
            col.names = TRUE)

# ---- Tables for manuscript ----
# ---- Supplementary Table 1: ASVs present in NS0.F, NS1.F, or NS0.FCd17, but absent in NS1.FCd17 ----
graphData = ASV_TableRelAb
graphData = graphData[,c(1:11,which(colnames(graphData) %in% c("NS0.F", "NS1.F" , "NS0.FCd17", "NS1.FCd17")))]
graphData = graphData[rowSums(graphData[,12:ncol(graphData)])>0,]
graphData = graphData[graphData$NS1.FCd17 == 0,]

graphData = graphData[,c(1,4,7,9,11,12:(ncol(graphData)-1))]

write.table(graphData, file = paste0(path,"Figure Outputs/",newDir,"/Supp_Table1.txt"), 
            quote = FALSE,
            sep = "\t",
            row.names =  FALSE, 
            col.names = TRUE)

# ---- Supplementary Table 3,4 cleanup: Sanger sequences of NS0 and NS1 isolates ----
# load NS0 isolates table
NS0_F_C_Isolates = readxl::read_xlsx(path = file.path(path,"input/Supplementary_table3.xlsx"))

# load NS1 isolates table
NS1_F_C_Isolates = readxl::read_xlsx(path = file.path(path,"input/Supplementary_table4.xlsx"))

# this sequence was reverse complemented
NS0_F_C_Isolates$`16S sequence (V3-V6 region)`[NS0_F_C_Isolates$`Bacterial strain ID` =="T1D_NS_0 10 FAA AN"] = as.character(reverseComplement(DNAString(NS0_F_C_Isolates$`16S sequence (V3-V6 region)`[NS0_F_C_Isolates$`Bacterial strain ID` =="T1D_NS_0 10 FAA AN"])))

# remove any "-" characters from sequencing
NS0_F_C_Isolates$`16S sequence (V3-V6 region)` = gsub("-","",NS0_F_C_Isolates$`16S sequence (V3-V6 region)`,fixed = TRUE)
NS1_F_C_Isolates$`16S sequence (V3-V6 region)` = gsub("-","",NS1_F_C_Isolates$`16S sequence (V3-V6 region)`,fixed = TRUE)

write.table(NS0_F_C_Isolates, file = paste0(path,"Figure Outputs/",newDir,"/Supp_Table3.txt"), 
            quote = FALSE,
            sep = "\t",
            row.names =  FALSE, 
            col.names = TRUE)

write.table(NS1_F_C_Isolates, file = paste0(path,"Figure Outputs/",newDir,"/Supp_Table4.txt"), 
            quote = FALSE,
            sep = "\t",
            row.names =  FALSE, 
            col.names = TRUE)

# ---- ASV vs Metabolite correlation tables ----
# ---- Supplementary Table 8: ASV vs NMR metabolite correlations, NS1DCV1-3d01-d17, without d1 ----
X_16S = data.frame(t(ASV_TableRelAb[,12:(sampleNum+12-1)])) %>% tibble::rownames_to_column(var = "sample_16S")
colnames(X_16S)[2:ncol(X_16S)] = ASV_TableRelAb$ASV_Num

X_16S = X_16S[grep("NS1.DC", X_16S$sample_16S, fixed = TRUE),]

# remove d01
X_16S = X_16S[!grepl("d01", X_16S$sample_16S, fixed = TRUE),]
# remove >d17
X_16S = X_16S[substr(X_16S$sample_16S, nchar(X_16S$sample_16S)-1,nchar(X_16S$sample_16S))<=17,]

X_metabolites = semi_join(NMRmetabolite_Table,X_16S, by = c(sample_16S = "sample_16S"))
X_16S = semi_join(X_16S,NMRmetabolite_Table, by = c(sample_16S = "sample_16S"))
X_16S = X_16S[,c(1,which(colSums(X_16S[2:ncol(X_16S)])!=0)+1)]
X_16S = X_16S[,colSums(X_16S != 0)>1]

X_metabolites = X_metabolites[match(X_16S$sample_16S,X_metabolites$sample_16S),]

cortest = adply(7:26, 1, .id = NULL, function(x){
  adply(2:ncol(X_16S),1, .id = NULL, function(y){
    test = cor.test(X_metabolites[,x],X_16S[,y])
    # test = cor.test(X_metabolites[,x],X_16S[,y], method = "spearman")
    df = data.frame("metabolite" = colnames(X_metabolites)[x],
                    "ASV" = colnames(X_16S)[y],
                    "cor" = test$estimate,
                    "pval" = test$p.value)
  })
})

ggplot()+
  geom_histogram(data = cortest, bins = 100, aes(x = pval))

cortest$BLASTspecies = ASV_TableRelAb$BLASTspecies[match(cortest$ASV,ASV_TableRelAb$ASV_Num)]
cortest$BLASTper_id = ASV_TableRelAb$BLASTper_id[match(cortest$ASV,ASV_TableRelAb$ASV_Num)]
cortest = cortest[,c(1,2,5,6,3,4)]
cortest$FDR = p.adjust(cortest$pval, method = "fdr")
cortest = cortest[order(cortest$pval),]
cortest$FDRrank = rank(cortest$FDR, ties.method = "max")
cortest$FDNum = cortest$FDR * cortest$FDRrank

write.table(cortest, file = paste0(path,"/Figure Outputs/",newDir,"/Supp_Table8.txt"), 
            quote = FALSE,
            sep = "\t",
            row.names =  FALSE, 
            col.names = TRUE)

# ---- Supplementary Table 9: ASV vs MS metabolite correlations, NS1DCV1-3d01-d17, without d1 ----
X_16S = data.frame(t(ASV_TableRelAb[,12:(sampleNum+12-1)])) %>% tibble::rownames_to_column(var = "sample_16S")
colnames(X_16S)[2:ncol(X_16S)] = ASV_TableRelAb$ASV_Num

X_16S = X_16S[grep("NS1.DC", X_16S$sample_16S, fixed = TRUE),]

# remove d00 and d01
X_16S = X_16S[!grepl("d00", X_16S$sample_16S, fixed = TRUE),]
X_16S = X_16S[!grepl("d01", X_16S$sample_16S, fixed = TRUE),]
# remove >d17
X_16S = X_16S[substr(X_16S$sample_16S, nchar(X_16S$sample_16S)-1,nchar(X_16S$sample_16S))<=17,]

X_metabolites = semi_join(MSmetabolite_TableDC_clean,X_16S, by = c(sample_16S = "sample_16S"))
X_16S = semi_join(X_16S,MSmetabolite_TableDC_clean, by = c(sample_16S = "sample_16S"))
X_16S = X_16S[,c(1,which(colSums(X_16S[2:ncol(X_16S)])!=0)+1)]
X_16S = X_16S[,colSums(X_16S != 0)>1]

X_metabolites = X_metabolites[match(X_16S$sample_16S,X_metabolites$sample_16S),]
X_metabolites[,8:ncol(X_metabolites)] = scale(X_metabolites[,8:ncol(X_metabolites)])

X_metabolites = cbind(X_metabolites[,1:7],X_metabolites[,8:ncol(X_metabolites)][,apply(X_metabolites[,8:ncol(X_metabolites)],2,function(x) !any(is.nan(x)))])

cortest = adply(8:ncol(X_metabolites), 1, .id = NULL, function(x){
  adply(2:ncol(X_16S),1, .id = NULL, function(y){
    test = cor.test(X_metabolites[,x],X_16S[,y])
    # test = cor.test(X_metabolites[,x],X_16S[,y], method = "spearman")
    df = data.frame("metabolite" = colnames(X_metabolites)[x],
                    "ASV" = colnames(X_16S)[y],
                    "cor" = test$estimate,
                    "pval" = test$p.value)
  })
})

ggplot()+
  geom_histogram(data = cortest, bins = 100, aes(x = pval))

cortest$BLASTspecies = ASV_TableRelAb$BLASTspecies[match(cortest$ASV,ASV_TableRelAb$ASV_Num)]
cortest$BLASTper_id = ASV_TableRelAb$BLASTper_id[match(cortest$ASV,ASV_TableRelAb$ASV_Num)]
cortest = cortest[,c(1,2,5,6,3,4)]
cortest$FDR = p.adjust(cortest$pval, method = "fdr")
cortest = cortest[order(cortest$pval),]
cortest$FDRrank = rank(cortest$FDR, ties.method = "max")
cortest$FDNum = cortest$FDR * cortest$FDRrank

write.table(cortest, file = paste0(path,"/Figure Outputs/",newDir,"/Supp_Table9.txt"), 
            quote = FALSE,
            sep = "\t",
            row.names =  FALSE, 
            col.names = TRUE)

# #### Figures ####
# ==== Full Phylum legend ====
graphData = ASV_TableRelAb[1:9,1:11]
graphData$BLASTphylum = c(
  "Actinomycetota",
  "Bacillota",
  "Bacteroidota",
  "Campylobacterota",
  "Fusobacteriota",
  "Lentisphaerota",
  "Pseudomonadota",
  "Thermodesulfobacteriota",
  "Verrucomicrobiota"
)

graphData$data = 1

graphDataMelt = melt(graphData, measure.vars = 12, variable.name = "Sample", value.name = "Relative_Abundance", stringsAsFactors = FALSE)

phylumLegendSource = ggplot() +
  geom_stratum(data = graphDataMelt, linetype = "solid",  aes(x = Sample, stratum = ASV_Num, y = Relative_Abundance*100, fill = BLASTphylum))+
  scale_x_discrete(name = "Chemostat Culture (day)")+
  scale_y_continuous(name = "Relative Abundance (%)")+
  scale_fill_manual(name = "Phylum", values = phylumColors)+
  alluvial.theme_x1+
  theme_void()

phylum.legend.gTable = get_legend(phylumLegendSource)
# phylum.legend.gTable = phylum.legend.gTable$grobs[[1]]


# ==== Supplementary Figures ====
# ---- Supp_1A: scatter plot, selected non-pooled vs pseudo-pooled samples ----
# keep NS0.F,NS1.F,NS0.FCd17,NS1.FCd17
ASV_TableSelectPooled = ASV_Table[,c("NS0.F","NS1.F","NS0.FCd17","NS1.FCd17")]
rownames(ASV_TableSelectPooled) = ASV_Table$sequences
colnames(ASV_TableSelectPooled) = paste0("pooled_",colnames(ASV_TableSelectPooled))
ASV_TableSelectPooled = t(ASV_TableSelectPooled)

ASV_TableSelectNonPooled = ASV_TableNon_pooled[,c("NS0.F","NS1.F","NS0.FCd17","NS1.FCd17")]
rownames(ASV_TableSelectNonPooled) = ASV_TableNon_pooled$sequences
colnames(ASV_TableSelectNonPooled) = paste0("nonpooled_",colnames(ASV_TableSelectNonPooled))
ASV_TableSelectNonPooled = t(ASV_TableSelectNonPooled)

ASV_TableSelectMerged = mergeSequenceTables(ASV_TableSelectNonPooled, ASV_TableSelectPooled)
ASV_TableSelectMerged_collapse = collapseNoMismatch2(ASV_TableSelectMerged,verbose = TRUE)

ASV_TableMerged = data.frame(ASV_Num = seq_along(colnames(ASV_TableSelectMerged_collapse)), sequences = colnames(ASV_TableSelectMerged_collapse),t(ASV_TableSelectMerged_collapse), stringsAsFactors = FALSE)
row.names(ASV_TableMerged) = NULL
ASV_TableMerged$ASV_Num = sprintf("%03d", ASV_TableMerged$ASV_Num)
ASV_TableMerged$ASV_Num = paste0("ASV_", ASV_TableMerged$ASV_Num)

seq2fasta(ASV_TableMerged$ASV_Num, ASV_TableMerged$sequences, "pseudo_pool_compare", file.path(path,"ASV FASTAs"))

blastResults = BLAST_16S(fastafileName = paste0(path,"/ASV FASTAs/pseudo_pool_compare.txt"),
                         remote = FALSE,
                         wait = TRUE,
                         topHit = TRUE,
                         taxonomy = TRUE)

ASV_TableMerged = data.frame(ASV_TableMerged[,1:2],
                       blastResults[match(ASV_TableMerged$ASV_Num, blastResults$query),c("BLASTsuperkingdom","BLASTphylum","BLASTclass","BLASTorder","BLASTfamily","BLASTgenus","BLASTspecies","BLASTmatch")],
                       BLASTper_id = blastResults$per_id[match(ASV_TableMerged$ASV_Num,blastResults$query)],
                       ASV_TableMerged[,3:ncol(ASV_TableMerged)])

write.table(ASV_TableMerged,file = file.path(path,"pseudo_pool_compare ASV_Table.txt"),sep = "\t", row.names = FALSE)
ASV_TableMerged = read.delim(file.path(path,"pseudo_pool_compare ASV_Table.txt"))

#convert to within-sample relative abundance
ASV_TableMergedRelAb = ASV_TableMerged
ASV_TableMergedRelAb[,12:ncol(ASV_TableMergedRelAb)] = colwise(function(x) x/sum(x))(ASV_TableMergedRelAb[,12:ncol(ASV_TableMergedRelAb)])

ASV_TableMergedRelAbMelt = melt(ASV_TableMergedRelAb,id.vars = 1:11)
ASV_TableMergedRelAbMelt$pooled = sapply(strsplit(as.character(ASV_TableMergedRelAbMelt$variable),"_"),function(x) x[1])
ASV_TableMergedRelAbMelt$sample = sapply(strsplit(as.character(ASV_TableMergedRelAbMelt$variable),"_"),function(x) x[2])
ASV_TableMergedRelAbMelt$sample = factor(ASV_TableMergedRelAbMelt$sample)
ASV_TableMergedRelAbMelt$variable = NULL

ASV_TableMergedRelAbMeltCast = dcast(ASV_TableMergedRelAbMelt, formula = ...~pooled, value.var = "value")
levels(ASV_TableMergedRelAbMeltCast$sample) = sampleNameLevels[match(levels(ASV_TableMergedRelAbMeltCast$sample),colnames(ASV_TableRelAb)[12:ncol(ASV_TableRelAb)])]
textData = ddply(ASV_TableMergedRelAbMeltCast, .(sample), function(x){
  nnp = sum(x$nonpooled > 0)
  npp = sum(x$pooled > 0)
  df = data.frame(nnp = nnp, npp = npp)
})

Supp_1A.plot.point = ggplot()+
  geom_point(data = ASV_TableMergedRelAbMeltCast, aes(x = nonpooled*100, y = pooled*100))+
  geom_abline(color = "grey90",slope = 1, intercept = 0)+
  scale_x_log10(name = "non-pooled\nRelative Abundance (%)", labels = label_scientific(digits = 1))+
  scale_y_log10(name = "pseudo-pooled\nRelative Abundance (%)", labels = label_scientific(digits = 1))+
  geom_text(data = textData, size = 7, parse = TRUE, aes(x = 10, y = 0.001, hjust = "inward", vjust = "inward", label = paste("n[np] ==",nnp)))+
  geom_text(data = textData, size = 7, parse = TRUE, aes(x = 0.001, y = 10, hjust = "inward", vjust = "inward",label = paste("n[pp] ==",npp)))+
  facet_wrap(~sample)+
  coord_fixed()+
  point.theme
# ---- Supp_1B: scatter plot, selected full-pooled vs pseudo-pooled samples ----
# keep NS0.F,NS1.F,NS0.FCd17,NS1.FCd17
ASV_TableSelectPooled = ASV_Table[,c("NS0.F","NS1.F","NS0.FCd17","NS1.FCd17")]
rownames(ASV_TableSelectPooled) = ASV_Table$sequences
colnames(ASV_TableSelectPooled) = paste0("pseudopooled_",colnames(ASV_TableSelectPooled))
ASV_TableSelectPooled = t(ASV_TableSelectPooled)

ASV_TableSelectFullPooled = ASV_TableFull_pooled[,c("NS0.F","NS1.F","NS0.FCd17","NS1.FCd17")]
rownames(ASV_TableSelectFullPooled) = ASV_TableFull_pooled$sequences
colnames(ASV_TableSelectFullPooled) = paste0("fullpooled_",colnames(ASV_TableSelectFullPooled))
ASV_TableSelectFullPooled = t(ASV_TableSelectFullPooled)

ASV_TableSelectMerged = mergeSequenceTables(ASV_TableSelectFullPooled, ASV_TableSelectPooled)
ASV_TableSelectMerged_collapse = collapseNoMismatch2(ASV_TableSelectMerged,verbose = TRUE)

ASV_TableMerged = data.frame(ASV_Num = seq_along(colnames(ASV_TableSelectMerged_collapse)), sequences = colnames(ASV_TableSelectMerged_collapse),t(ASV_TableSelectMerged_collapse), stringsAsFactors = FALSE)
row.names(ASV_TableMerged) = NULL
ASV_TableMerged$ASV_Num = sprintf("%03d", ASV_TableMerged$ASV_Num)
ASV_TableMerged$ASV_Num = paste0("ASV_", ASV_TableMerged$ASV_Num)

seq2fasta(ASV_TableMerged$ASV_Num, ASV_TableMerged$sequences, "full_pool_compare", file.path(path,"ASV FASTAs"))

blastResults = BLAST_16S(fastafileName = paste0(path,"/ASV FASTAs/full_pool_compare.txt"),
                         remote = FALSE,
                         wait = TRUE,
                         topHit = TRUE,
                         taxonomy = TRUE)

ASV_TableMerged = data.frame(ASV_TableMerged[,1:2],
                             blastResults[match(ASV_TableMerged$ASV_Num, blastResults$query),c("BLASTsuperkingdom","BLASTphylum","BLASTclass","BLASTorder","BLASTfamily","BLASTgenus","BLASTspecies","BLASTmatch")],
                             BLASTper_id = blastResults$per_id[match(ASV_TableMerged$ASV_Num,blastResults$query)],
                             ASV_TableMerged[,3:ncol(ASV_TableMerged)])

write.table(ASV_TableMerged,file = file.path(path,"full_pool_compare ASV_Table.txt"),sep = "\t", row.names = FALSE)
ASV_TableMerged = read.delim(file.path(path,"full_pool_compare ASV_Table.txt"))

#convert to within-sample relative abundance
ASV_TableMergedRelAb = ASV_TableMerged
ASV_TableMergedRelAb[,12:ncol(ASV_TableMergedRelAb)] = colwise(function(x) x/sum(x))(ASV_TableMergedRelAb[,12:ncol(ASV_TableMergedRelAb)])

ASV_TableMergedRelAbMelt = melt(ASV_TableMergedRelAb,id.vars = 1:11)
ASV_TableMergedRelAbMelt$pooled = sapply(strsplit(as.character(ASV_TableMergedRelAbMelt$variable),"_"),function(x) x[1])
ASV_TableMergedRelAbMelt$sample = sapply(strsplit(as.character(ASV_TableMergedRelAbMelt$variable),"_"),function(x) x[2])
ASV_TableMergedRelAbMelt$sample = factor(ASV_TableMergedRelAbMelt$sample)
ASV_TableMergedRelAbMelt$variable = NULL

ASV_TableMergedRelAbMeltCast = dcast(ASV_TableMergedRelAbMelt, formula = ...~pooled, value.var = "value")
levels(ASV_TableMergedRelAbMeltCast$sample) = sampleNameLevels[match(levels(ASV_TableMergedRelAbMeltCast$sample),colnames(ASV_TableRelAb)[12:ncol(ASV_TableRelAb)])]
textData = ddply(ASV_TableMergedRelAbMeltCast, .(sample), function(x){
  nfp = sum(x$fullpooled > 0)
  npp = sum(x$pseudopooled > 0)
  df = data.frame(nfp = nfp, npp = npp)
})

Supp_1B.plot.point = ggplot()+
  geom_point(data = ASV_TableMergedRelAbMeltCast, aes(x = fullpooled*100, y = pseudopooled*100))+
  geom_abline(color = "grey90",slope = 1, intercept = 0)+
  scale_x_log10(name = "full-pooled\nRelative Abundance (%)", labels = label_scientific(digits = 1))+
  scale_y_log10(name = "pseudo-pooled\nRelative Abundance (%)", labels = label_scientific(digits = 1))+
  geom_text(data = textData, size = 7, parse = TRUE, aes(x = 10, y = 0.001, hjust = "inward", vjust = "inward", label = paste("n[fp] ==",nfp)))+
  geom_text(data = textData, size = 7, parse = TRUE, aes(x = 0.001, y = 10, hjust = "inward", vjust = "inward",label = paste("n[pp] ==",npp)))+
  facet_wrap(~sample)+
  coord_fixed()+
  point.theme

# ---- Supp_2A-J: Alluvial plot, all samples F, FCd01-d17, colour by F, FCd01-d17 presence ----
# graphData = ASV_TableRelAb
graphData = ASV_Table_7parRelAb
graphData = graphData[,c(1:11,grep(".F",colnames(graphData), fixed = TRUE))]
# subset graphData by sample
graphDataList = alply(c("NS0","NS1","S2","S3","NS4","S5","NS6"), 1, function(sampleName){
  # graphDataList = alply(c("NS0","NS1"), 1, function(sampleName){
  sampleName.F = paste0(sampleName,".F")
  
  output = graphData[,c(1:11,grep(sampleName,colnames(graphData), fixed = TRUE))]
  
  # order by relative abundance in F sample and NS0.FC samples
  output = output[order(output[,sampleName.F], decreasing = TRUE),]
  output = output[order(rowSums(output[,grep("FC",colnames(output), fixed = TRUE)])),]
  
  #order by phylum
  # output = output[order(output$BLASTphylum),]
  
  output$presenceGroupColor = interaction(output[,sampleName.F]!=0,
                                          rowSums(output[,grep("FC",colnames(output), fixed = TRUE)])!=0)
  
  #order by sample presence
  output = output[order(output$presenceGroupColor),]
  
  # generate percentage data for absence/presence data
  legendData1 = data.frame("day.0" = output[,sampleName.F],
                           "days.1_17" = rowMeans(output[,grep("FC",colnames(output), fixed = TRUE)]))
  legendData = aggregate(legendData1, by = list("presenceGroupColor" = interaction(legendData1$day.0!=0,
                                                                                   legendData1$days.1_17!=0)),
                         FUN = function(x) sum(x))
  
  legendData = legendData[,1:2]
  
  legendData$label = paste0(sprintf("%.1f",round(legendData$day.0*100,1)),"%")
  legendData$label[legendData$day.0==0] = ""
  
  legendDataMelt = melt(legendData, measure.vars = 2, variable.name = "Sample", value.name = "Relative_Abundance", stringsAsFactors = FALSE)
  legendDataMelt = legendDataMelt[legendDataMelt$Relative_Abundance!=0,]
  legendDataMelt$day_16S = 0
  
  # generate legend data for absence/presence tile legend
  legendDataTile = aggregate(legendData1, by = list("levelName" = output$presenceGroupColor), FUN = function(x) sum(x)*100)
  
  legendDataTile = legendDataTile[-1,]
  legendDataTile = legendDataTile[match(absencePresence2LegendData$levelName,legendDataTile$levelName),]
  
  legendDataTile$y = c(2,1,0)
  
  legendDataTileMelt = melt(legendDataTile, id.vars = c(1,ncol(legendDataTile)), variable.name = "x")
  legendDataTileMelt$filled = as.numeric(legendDataTileMelt$value>0)
  
  legendDataTileMelt$x = gsub(".","\n",legendDataTileMelt$x, fixed = TRUE)
  legendDataTileMelt$x = gsub("_","-",legendDataTileMelt$x, fixed = TRUE)
  
  legendTile = ggplot()+
    geom_tile(data = legendDataTileMelt[legendDataTileMelt$filled == 1,],
              color = "black",
              size = 0.5,
              aes(x = x, y = y, fill = levelName))+
    # geom_text(data = legendDataMelt[legendDataMelt$filled == 1,],
    #           size = 5.5,
    #           aes(x = x, y = y, label = value))+
    geom_tile(data = legendDataTileMelt[legendDataTileMelt$filled == 0,],
              color = "black",
              fill = "white",
              size = 0.5,
              aes(x = x, y = y))+
    # geom_tile(data = legendDataTileMelt,
    #           color = "black",
    #           height = 0.70,
    #           width = 0.70,
    #           size = 0.25,
    #           aes(x = "", y = y, fill = levelName))+
    scale_x_discrete(position = "top")+   
    scale_y_continuous(expand = expansion(add = 0.1))+
    ggtitle("Identified on:")+
    coord_fixed()+
    legend.day2.theme
  
  
  #renumber ASVs based on reordering
  output$ASV_Num = seq_along(output$ASV_Num)
  output$ASV_Num = sprintf("%04d", output$ASV_Num)
  
  outputMelt = melt(output, measure.vars = 12:(ncol(output)-1), variable.name = "Sample", value.name = "Relative_Abundance", stringsAsFactors = FALSE)
  outputMelt = outputMelt[outputMelt$Relative_Abundance!=0,]
  
  outputMelt$day_16S = sample_Table$day_16S[match(outputMelt$Sample,sample_Table$sample_16S)]
  
  # return(outputMelt)
  
  #plots
  #colored by sample presence
  Supp_2.plot.alluvial.ap = ggplot() + 
    geom_stratum(data = outputMelt, linetype = "solid",  aes(x = factor(day_16S), stratum = ASV_Num, y = Relative_Abundance*100, fill = presenceGroupColor))+
    geom_errorbar(data = legendDataMelt, stat = "stratum", position = position_nudge(x = -0.3), width = 0.1, size = 1, aes(x = factor(day_16S), stratum = presenceGroupColor, y = Relative_Abundance*100))+
    geom_text(data = legendDataMelt, stat = "stratum", nudge_x = -0.375, angle = 90, size = 6, aes(x = factor(day_16S), stratum = presenceGroupColor, y = Relative_Abundance*100, label = label, vjust = 0, hjust = 0.5))+
    geom_flow(data = outputMelt, color = "black",aes(x = factor(day_16S),stratum = ASV_Num, alluvium = ASV_Num, y = Relative_Abundance*100, fill = presenceGroupColor))+
    scale_x_discrete(name = "Chemostat Culture (day)", expand = expansion(add = c(0.6,0.4)))+
    scale_y_continuous(name = "Relative Abundance (%)")+
    # scale_fill_manual(values = absencePresenceColors3)+
    alluvial.theme_x1_legend0
  
  #colored by phylum
  Supp_2.plot.alluvial.ph = ggplot() +
    geom_stratum(data = outputMelt, linetype = "solid",  aes(x = factor(day_16S), stratum = ASV_Num, y = Relative_Abundance*100, fill = BLASTphylum))+
    geom_flow(data = outputMelt, color = "black",aes(x = factor(day_16S),stratum = ASV_Num, alluvium = ASV_Num, y = Relative_Abundance*100, fill = BLASTphylum))+
    scale_x_discrete(name = "Chemostat Culture (day)", expand = expansion(add = c(0.6,0.4)))+
    scale_y_continuous(name = "Relative Abundance (%)")+
    scale_fill_manual(name = "Phylum", values = phylumColors)+
    alluvial.theme_x1_legend0
  
  return(list(Supp_2.plot.alluvial.ap,Supp_2.plot.alluvial.ph,legendTile))
})
names(graphDataList) = c("NS0","NS1","S2","S3","NS4","S5","NS6")

Supp_2A.plot.alluvial.ap = graphDataList[[3]][[1]]
Supp_2B.plot.alluvial.ap = graphDataList[[4]][[1]]
Supp_2C.plot.alluvial.ap = graphDataList[[5]][[1]]
Supp_2D.plot.alluvial.ap = graphDataList[[6]][[1]]
Supp_2E.plot.alluvial.ap = graphDataList[[7]][[1]]

Supp_2.legend.ap = graphDataList[[1]][[3]]

Supp_2F.plot.alluvial.ph = graphDataList[[3]][[2]]
Supp_2G.plot.alluvial.ph = graphDataList[[4]][[2]]
Supp_2H.plot.alluvial.ph = graphDataList[[5]][[2]]
Supp_2I.plot.alluvial.ph = graphDataList[[6]][[2]]
Supp_2J.plot.alluvial.ph = graphDataList[[7]][[2]]

Supp_2.ph.legend.gTable = phylum.legend.gTable

install.packages("patchwork")
library("patchwork")

test = 
(
  (
    (
      Supp_2A.plot.alluvial.ap|
      Supp_2B.plot.alluvial.ap|
      Supp_2C.plot.alluvial.ap|
      Supp_2D.plot.alluvial.ap|
      Supp_2E.plot.alluvial.ap|
      free(Supp_2.legend.ap)
      )+plot_layout(widths = c(1,1,1,1,1,0.4), axis_titles = "collect"))/
  (
    (
      Supp_2F.plot.alluvial.ph| 
      Supp_2G.plot.alluvial.ph| 
      Supp_2H.plot.alluvial.ph| 
      Supp_2I.plot.alluvial.ph| 
      Supp_2J.plot.alluvial.ph| 
      Supp_2.ph.legend.gTable
      )+plot_layout(widths = c(1,1,1,1,1,0.4), axis_titles = "collect"))/ 
Supp_2K.plot.PCoA
)+plot_layout(heights = c(1,1,3))

ggsave(paste0(path,"Figure Outputs/",newDir,"/test.pdf"),
       plot = test,
       device=cairo_pdf, # add this to embed the Calibri font
       width = 45, height = 30,
       units = "in", # other options c("in", "cm", "mm"),
       dpi = 300,
       limitsize = FALSE)


test = plot_grid(Supp_2A.plot.alluvial.ap,
               Supp_2B.plot.alluvial.ap,
               Supp_2C.plot.alluvial.ap,
               Supp_2D.plot.alluvial.ap,
               Supp_2E.plot.alluvial.ap,
               Supp_2.legend.ap,
               Supp_2F.plot.alluvial.ph,
               Supp_2G.plot.alluvial.ph,
               Supp_2H.plot.alluvial.ph,
               Supp_2I.plot.alluvial.ph,
               Supp_2J.plot.alluvial.ph,
               Supp_2.ph.legend.gTable
               ,
          nrow = 2,
          ncol = 6,
          rel_widths = c(1,1,1,1,1,0.4),
          rel_heights	= c(1,1)
          )

test = plot_grid(test,
                 Supp_2K.plot.PCoA,
                 nrow = 2,
                 ncol = 1,
                 rel_widths = c(1,1),
                 rel_heights	= c(1,1.5)
)


plot_grid(Supp_2A.plot.alluvial.ap,
               Supp_2B.plot.alluvial.ap,
               Supp_2C.plot.alluvial.ap,
               Supp_2D.plot.alluvial.ap,
               Supp_2E.plot.alluvial.ap,
nrow = 1,
ncol = 5,
rel_widths = c(1,1,1,1,1),
rel_heights	= c(1)
)

# ---- Supp_2K: PCoA plot, all samples F, FCd13, colour by source, shape by type, path by source ----
graphData = ASV_Table_7parRelAb
graphData = graphData[,c(grep("F$",colnames(graphData),fixed = FALSE),grep("FCd13",colnames(graphData),fixed = FALSE))]
BC_dis = vegdist(t(graphData))

#PCoA
PCoA = cmdscale(BC_dis, k = 2, eig = TRUE)

Plot_PCoA = sample_Table %>%
  right_join(tibble::rownames_to_column(data.frame(PCoA$points), var = "sample_16S") %>%
               rename(PCo1 = X1,
                      PCo2 = X2))
Plot_PCoA$type = c("F" = "Fecal", "FC" = "Chemostat")[Plot_PCoA$type]

PCoA.percentvariance = PCoA$eig/sum(PCoA$eig) *100

Supp_2K.plot.PCoA = ggplot()+
  geom_point(data = Plot_PCoA, size = 3, aes(x = PCo1, y =PCo2, color = source, shape = type))+
  geom_path(data = Plot_PCoA,aes(x = PCo1, y =PCo2, group = source))+
  scale_x_continuous(name = paste0("PCo1 (",signif(PCoA.percentvariance[1],3),"%)"))+
  scale_y_continuous(name = paste0("PCo2 (",signif(PCoA.percentvariance[2],3),"%)"))+
  scale_shape_discrete(name = "sample type")+
  coord_fixed()+
  PCoA.theme
# ---- Supp_3: Dendrogram and Boxplot, NS0 and NS1 isolates, phylogenetic relationship and relative abundance in CHILD study samples ---- 
# # load DIABIMMUNE combined ASV table
# DIABIMMUNE_ASV_Table = read.delim(paste0(path,"/DIABIMMUNE_combined ASV_Table.txt"),sep = "\t", stringsAsFactors = FALSE)
# DIABIMMUNE_metadata = read.delim(paste0(path,"/DIABIMMUNE_combined metadata.txt"),sep = "\t", stringsAsFactors = FALSE)

# load CHILD 12 month ASV table
CHILD_12M_ASV_Table = read.delim(paste0(path,"/CHILD_12M ASV_Table.txt"),sep = "\t", stringsAsFactors = FALSE)

CHILD_ASV_seqtab = as.matrix(t(CHILD_12M_ASV_Table[,12:ncol(CHILD_12M_ASV_Table)]))
CHILDSampleCount = nrow(CHILD_ASV_seqtab)
colnames(CHILD_ASV_seqtab) = CHILD_12M_ASV_Table$sequences

graphData = ASV_Table
# select samples
graphData = graphData[,c(1:11,which(colnames(graphData) %in% c("NS0.DCd00","NS1.DCd00")))]
# remove ASVs not in DC
graphData = graphData[rowSums(graphData[,12:ncol(graphData)])>0,]
# make isolate ASVs most abundant to ensure they are kept post-collapse
graphData = data.frame(graphData,"Isolates_presence" = ((rowSums(graphData[,12:ncol(graphData)])>0)+0)*10^7)
graphData_seqtab = as.matrix(t(graphData[,12:ncol(graphData)]))
colnames(graphData_seqtab) = graphData$sequences
rownames(graphData_seqtab) = c("NS0.DCd00","NS1.DCd00","Isolates_presence")

# merge CHILD study and isolates seqtabs
merged_seqtab = mergeSequenceTables(CHILD_ASV_seqtab, graphData_seqtab)
merged_seqtab_collapsed = collapseNoMismatch2(merged_seqtab,verbose = TRUE, minOverlap = 100)

# reformat as ASV table
mergedASV_Table = data.frame(ASV_Num = seq_along(colnames(merged_seqtab_collapsed)), sequences = colnames(merged_seqtab_collapsed),t(merged_seqtab_collapsed), stringsAsFactors = FALSE)
row.names(mergedASV_Table) = NULL

mergedASV_Table$ASV_Num = sprintf("%04d", mergedASV_Table$ASV_Num)
mergedASV_Table$ASV_Num = paste0("ASV_", mergedASV_Table$ASV_Num)

seq2fasta(mergedASV_Table$ASV_Num, mergedASV_Table$sequences, "CHILD_prevalence", file.path(path,"ASV FASTAs"))

blastResults = BLAST_16S(fastafileName = paste0(path,"/ASV FASTAs/CHILD_prevalence.txt"),
                         remote = FALSE,
                         wait = TRUE,
                         topHit = TRUE,
                         taxonomy = TRUE)

mergedASV_Table = data.frame(mergedASV_Table[,1:2],
                             blastResults[match(mergedASV_Table$ASV_Num, blastResults$query),c("BLASTsuperkingdom","BLASTphylum","BLASTclass","BLASTorder","BLASTfamily","BLASTgenus","BLASTspecies","BLASTmatch")],
                             BLASTper_id = blastResults$per_id[match(mergedASV_Table$ASV_Num,blastResults$query)],
                             mergedASV_Table[,3:ncol(mergedASV_Table)])

write.table(mergedASV_Table,file = file.path(path,"CHILD_prevalence ASV_Table.txt"),sep = "\t", row.names = FALSE)
CHILD_mergedASV_Table = read.delim(file.path(path,"CHILD_prevalence ASV_Table.txt"))

# prevalence in the CHILD study samples = # of samples where ASV is >0 / # of samples
CHILD_mergedASV_Table$CHILD_prevalence = 
  rowSums(CHILD_mergedASV_Table[,grep(pattern = "X\\d",x = colnames(CHILD_mergedASV_Table), fixed = FALSE)]>0) / 
  sum(grepl(pattern = "X\\d",x = colnames(CHILD_mergedASV_Table), fixed = FALSE))

# CHILD_mergedASV_Table = CHILD_mergedASV_Table %>% dplyr::relocate(Isolates_presence, .before = X7104822814)

# subset out ASVs in isolates
CHILD_mergedASV_Table_isoSubset = CHILD_mergedASV_Table[CHILD_mergedASV_Table$Isolates_presence>0,]
CHILD_mergedASV_Table_nonisoSubset = CHILD_mergedASV_Table[CHILD_mergedASV_Table$Isolates_presence==0,]

# load NS0 isolates table
NS0_F_C_Isolates = readxl::read_xlsx(path = file.path(path,"input/Supplementary_table3.xlsx"))

# load NS1 isolates table
NS1_F_C_Isolates = readxl::read_xlsx(path = file.path(path,"input/Supplementary_table4.xlsx"))

CheckTableNS0 = data.frame(seqName = NS0_F_C_Isolates$`Bacterial strain ID`,
                           source = "NS0",
                           isoNum = paste0("NS0",1:nrow(NS0_F_C_Isolates)),
                           seq = as.character(reverseComplement(DNAStringSet(NS0_F_C_Isolates$`16S sequence (V3-V6 region)`))),
                           NS0count = 1,
                           NS1count = 0,
                           count = 1)

CheckTableNS0$seq[CheckTableNS0$seqName =="T1D_NS_0 10 FAA AN"] = as.character(reverseComplement(DNAString(CheckTableNS0$seq[CheckTableNS0$seqName =="T1D_NS_0 10 FAA AN"])))

CheckTableNS1 = data.frame(seqName = NS1_F_C_Isolates$`Bacterial strain ID`,
                           source = "NS1",
                           isoNum = paste0("NS1",1:nrow(NS1_F_C_Isolates)),
                           seq = as.character(reverseComplement(DNAStringSet(NS1_F_C_Isolates$`16S sequence (V3-V6 region)`))),
                           NS0count = 0,
                           NS1count = 1,
                           count = 1)

CheckTableSanger = rbind(CheckTableNS0,CheckTableNS1)
CheckTableSanger$seq = gsub("-","",CheckTableSanger$seq,fixed = TRUE)

CheckTableIsoSubset = data.frame(seqName = CHILD_mergedASV_Table_isoSubset$ASV_Num,
                                 source = "iso",
                                 isoNum = CHILD_mergedASV_Table_isoSubset$ASV_Num,
                                 seq = CHILD_mergedASV_Table_isoSubset$sequences,
                                 NS0count = CHILD_mergedASV_Table_isoSubset$NS0.DCd00,
                                 NS1count = CHILD_mergedASV_Table_isoSubset$NS1.DCd00,
                                 count = CHILD_mergedASV_Table_isoSubset$NS0.DCd00+CHILD_mergedASV_Table_isoSubset$NS1.DCd00)

CheckTable = rbind(CheckTableSanger,CheckTableIsoSubset)

seq2fasta(CheckTable$isoNum, CheckTable$seq, "CheckTable", file.path(path,"ASV FASTAs"))

blastResults = BLAST_16S(fastafileName = paste0(path,"/ASV FASTAs/CheckTable.txt"),
                         remote = FALSE,
                         wait = TRUE,
                         topHit = TRUE,
                         taxonomy = TRUE)

CheckTable = data.frame(CheckTable[,1:4],
                        blastResults[match(CheckTable$isoNum, blastResults$query),c("BLASTsuperkingdom","BLASTphylum","BLASTclass","BLASTorder","BLASTfamily","BLASTgenus","BLASTspecies","BLASTmatch")],
                        BLASTper_id = blastResults$per_id[match(CheckTable$isoNum,blastResults$query)],
                        CheckTable[,5:ncol(CheckTable)])

write.table(CheckTable,file = file.path(path,"CheckTable ASV_Table.txt"),sep = "\t", row.names = FALSE)
CheckTable = read.delim(file.path(path,"CheckTable ASV_Table.txt"))

CheckTableClust = DECIPHER_clusterASVs(CheckTable$seq, testBounds = TRUE,
                                       lowerBound = 0.80, upperBound = 1, testIncrement = 0.001)

CheckTableClust.sanger_not_iso_Count = apply(CheckTableClust,2, function(col){
  sanger = col[CheckTable$source %in% c("NS0","NS1")]
  iso = col[CheckTable$source %in% c("iso")]
  
  sum(!(sanger %in% iso))
})

CheckTableClust.sanger_not_iso_Species = apply(CheckTableClust,2, function(col){
  sanger = CheckTable$source %in% c("NS0","NS1")
  iso = CheckTable$source %in% c("iso")
  
  sangerCol = col[sanger]
  isoCol = col[iso]
  out = rep(NA, nrow(CheckTable))
  out[col %in% (sangerCol[!(sangerCol %in% isoCol)])] = 
    CheckTable$BLASTspecies[col %in% (sangerCol[!(sangerCol %in% isoCol)])]
  out
})

CheckTableClust.iso_not_sanger_Count = apply(CheckTableClust,2, function(col){
  sanger = CheckTable$source %in% c("NS0","NS1")
  iso = CheckTable$source %in% c("iso")
  
  sangerCol = col[sanger]
  isoCol = col[iso]
  
  sum(CheckTable$count[col %in% (isoCol[!(isoCol %in% sangerCol)])])
  
})

# filter isoSubset ASVs: remove ASV as spurious if they don't cluster with an isolate Sanger sequence at 0.941 identity
# filterClustCutoffCol = CheckTableClust$cluster_0_082
filterClustCutoffCol = CheckTableClust$cluster_0_059

filterClustNames = CheckTable$seqName[filterClustCutoffCol %in% 
                                        (filterClustCutoffCol[CheckTable$source %in% c("iso")][!(filterClustCutoffCol[CheckTable$source %in% c("iso")] %in% 
                                                                                                   filterClustCutoffCol[CheckTable$source %in% c("NS0","NS1")])])]

CHILD_mergedASV_Table_isoSubsetFiltered = CHILD_mergedASV_Table_isoSubset[!(CHILD_mergedASV_Table_isoSubset$ASV_Num %in%filterClustNames),]

# recombine CHILD_mergedASV_Table
CHILD_mergedASV_Table = rbind(CHILD_mergedASV_Table_isoSubsetFiltered,CHILD_mergedASV_Table_nonisoSubset)

# collapse to species level by BLASTspecies names,
# keeping metadata from ASV present in NS0 and NS1 isolates and most prevalent in CHILD study samples as "representative"
CHILD_mergedASV_Table_SpeciesCol = CHILD_mergedASV_Table[order(CHILD_mergedASV_Table$CHILD_prevalence, decreasing = TRUE),]
CHILD_mergedASV_Table_SpeciesCol = CHILD_mergedASV_Table_SpeciesCol[order(CHILD_mergedASV_Table_SpeciesCol$Isolates_presence, decreasing = TRUE),]

CHILD_mergedASV_Table_SpeciesCol = data.frame(CHILD_mergedASV_Table_SpeciesCol[which(!duplicated(CHILD_mergedASV_Table_SpeciesCol$BLASTspecies)),1:11],
                                              rowsum(CHILD_mergedASV_Table_SpeciesCol[12:ncol(CHILD_mergedASV_Table_SpeciesCol)],
                                                     group = CHILD_mergedASV_Table_SpeciesCol$BLASTspecies,
                                                     reorder = FALSE)
)

# recalculate prevalence in the CHILD study samples for species collapsed ASVs
# prevalence in the CHILD study samples = # of samples where ASV is >0 / # of samples
CHILD_mergedASV_Table_SpeciesCol$CHILD_prevalence = 
  rowSums(CHILD_mergedASV_Table_SpeciesCol[,grep(pattern = "X\\d",x = colnames(CHILD_mergedASV_Table_SpeciesCol), fixed = FALSE)]>0) / 
  sum(grepl(pattern = "X\\d",x = colnames(CHILD_mergedASV_Table_SpeciesCol), fixed = FALSE))

# remove species not in isolates
CHILD_mergedASV_Table_SpeciesCol = CHILD_mergedASV_Table_SpeciesCol[CHILD_mergedASV_Table_SpeciesCol$Isolates_presence>0,]

# check # of ASVs in CHILD_mergedASV_Table_SpeciesCol which share taxonomic labels with isolate Sanger sequences
sum((CHILD_mergedASV_Table_SpeciesCol$BLASTspecies %in% CheckTable$BLASTspecies[grep("NS_0|G37020",CheckTable$seqName)]))
sum((CHILD_mergedASV_Table_SpeciesCol$BLASTgenus %in% CheckTable$BLASTgenus[grep("NS_0|G37020",CheckTable$seqName)]))
sum((CHILD_mergedASV_Table_SpeciesCol$BLASTfamily %in% CheckTable$BLASTfamily[grep("NS_0|G37020",CheckTable$seqName)]))
sum((CHILD_mergedASV_Table_SpeciesCol$BLASTorder %in% CheckTable$BLASTorder[grep("NS_0|G37020",CheckTable$seqName)]))
nrow(CHILD_mergedASV_Table_SpeciesCol)

# CHILD_mergedASV_Table_SpeciesCol = CHILD_mergedASV_Table_SpeciesCol %>% dplyr::relocate(Isolates_presence, .before = X7104822814)

# Kalign 16S sequences to construct tree
seq2fasta(CHILD_mergedASV_Table_SpeciesCol$ASV_Num, CHILD_mergedASV_Table_SpeciesCol$sequences, "CHILD_prevalence_filtered", file.path(path,"ASV FASTAs"))

# align to construct tree (other options?)
# idMat = Kalign_pim(userEmail = "christopher.yau@mail.utoronto.ca",
#                    fastafileName = paste0(path,"/ASV FASTAs/CHILD_prevalence_filtered.txt"))

idMat = read.delim(paste0(path,"/ASV FASTAs/CHILD_prevalence_filtered.txt_Kalign_out.txt"),header = FALSE, row.names = 1, stringsAsFactors = FALSE)

idMat = (100-idMat)/100
idMatDist = as.dist(idMat)

isolatescluster = hclust(idMatDist, method = "complete")

# reorder ASV table and add treeOrder column
CHILD_mergedASV_Table_SpeciesCol = CHILD_mergedASV_Table_SpeciesCol[isolatescluster$order,]
CHILD_mergedASV_Table_SpeciesCol$treeOrder = paste0("treeOrder_", sprintf("%03d", seq_along(CHILD_mergedASV_Table_SpeciesCol$ASV_Num)))
CHILD_mergedASV_Table_SpeciesCol = CHILD_mergedASV_Table_SpeciesCol %>% 
  relocate(treeOrder)

# convert CHILD study sample count into relative abundance
CHILD_mergedASV_Table_SpeciesCol[,grep("X\\d" ,colnames(CHILD_mergedASV_Table_SpeciesCol), fixed = FALSE)] = sweep(CHILD_mergedASV_Table_SpeciesCol[,grep("X\\d" ,colnames(CHILD_mergedASV_Table_SpeciesCol), fixed = FALSE)],
                                                                                                                   2,
                                                                                                                   colSums(CHILD_mergedASV_Table_SpeciesCol[,grep("X\\d" ,colnames(CHILD_mergedASV_Table_SpeciesCol), fixed = FALSE)]),
                                                                                                                   FUN = '/')

# make ggdend for dendrogram
merged_ggdend = as.dendrogram(isolatescluster) %>%
  set_labels(labels = CHILD_mergedASV_Table_SpeciesCol$BLASTspecies)
merged_ggdend = as.ggdend(merged_ggdend)

# melt table for boxplot data
CHILD_mergedASV_Table_SpeciesColBoxplotMelt = reshape2::melt(CHILD_mergedASV_Table_SpeciesCol,
                                                             measure.vars = grep("X\\d" ,colnames(CHILD_mergedASV_Table_SpeciesCol), fixed = FALSE),
                                                             variable.name = "sample",
                                                             value.name = "relative_abundance")

# for Kalign-based dendrogram plot
Supp_3.dendrogram = ggplot()+
  # plotting branches
  geom_segment(data = merged_ggdend$segments,
               linewidth = 0.1,
               # y is moved into negatives to match boxplot y-axis scale (for spacing)
               aes(x = x, y = y-1, xend = xend, yend = yend-1))+ 
  # plotting leaves
  geom_point(data = merged_ggdend$nodes[!is.na(merged_ggdend$nodes$leaf),],
             shape = 22, size = 3.5,
             # y is moved into negatives to match boxplot y-axis scale (for spacing)
             aes(x = x, y = y-1, color = CHILD_mergedASV_Table_SpeciesCol$BLASTphylum, fill = CHILD_mergedASV_Table_SpeciesCol$BLASTphylum))+
  # plotting species labels
  geom_text(data = merged_ggdend$labels,
            # nudge_y to avoid overlap with leaves
            hjust = 1, vjust = 0.5, size = 3.5, nudge_y = -0.02, angle = 90,
            # y is moved into negatives to match boxplot y-axis scale (for spacing)
            aes(x = x, y = y-1, label = merged_ggdend$labels$label))+
  ylab("log10(rel. abundance)")+ #placeholder text to match boxplot y-axis scale (for spacing)
  # 0.6 units is default expansion for x-axis for discrete scales in ggplot2, matches boxplot x-axis (for spacing)
  scale_x_continuous(expand = expansion(add = c(0.6, 0.6)))+
  # 0.7 unit expansion below y-axis prevents species labels from being cut off in combined plot (manually set)
  # "-1" break is shown to match boxplot y-axis scale (for spacing)
  scale_y_continuous(expand = expansion(add = c(0.7,0)), breaks = c(-1))+
  scale_color_manual(name = "Phylum", values = phylumColors)+
  scale_fill_manual(name = "Phylum", values = phylumColors,
                    guide = guide_legend(override.aes = list(size = 10)))+
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        # y-axis elements are coloured white to act as spacers to match boxplot
        axis.title.y = element_text(color = "white"),
        axis.text.y = element_text(color = "white"),
        axis.ticks.y = element_line(color = "white"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        # removing bottom margin to place dendrogram and boxplot closer
        plot.margin = margin(5.5,5.5,0,5.5, unit = "pt"))+
  theme(legend.title = element_text(size=legendSize, face = "bold", family="Helvetica"),
        legend.text = element_text(size=textSize, face = "bold", family="Helvetica"),
        legend.key = element_blank())

Supp_3.legend.gTable = get_legend(Supp_3.dendrogram)
Supp_3.dendrogram = Supp_3.dendrogram + 
  theme(legend.position = "none")

Supp_3.dendrogram.CoordFlip =
  ggplot()+
  # plotting branches
  geom_segment(data = merged_ggdend$segments,
               linewidth = 0.1,
               # y is moved into negatives to match boxplot y-axis scale (for spacing)
               aes(y = x, x = y-1, yend = xend, xend = yend-1))+ 
  # plotting leaves
  geom_point(data = merged_ggdend$nodes[!is.na(merged_ggdend$nodes$leaf),],
             shape = 22, size = 3.5,
             # y is moved into negatives to match boxplot y-axis scale (for spacing)
             # 0.05 is subtracted from y to keep leaves from overlapping dendrogram lines (manually set)
             aes(y = x, x = y-1, color = CHILD_mergedASV_Table_SpeciesCol$BLASTphylum, fill = CHILD_mergedASV_Table_SpeciesCol$BLASTphylum))+
  # plotting species labels
  geom_text(data = merged_ggdend$labels,
            # nudge_y to avoid overlap with leaves
            hjust = 0, vjust = 0.5, size = 3.5, nudge_x = 0.02, angle = 0,
            # y is moved into negatives to match boxplot y-axis scale (for spacing)
            aes(y = x, x = y-1, label = merged_ggdend$labels$label))+
  xlab("log10(rel. abundance)")+ #placeholder text to match boxplot y-axis scale (for spacing)
  # 0.6 units is default expansion for x-axis for discrete scales in ggplot2, matches boxplot x-axis (for spacing)
  scale_y_continuous(expand = expansion(add = c(0.6, 0.6)))+
  # 0.7 unit expansion below y-axis prevents species labels from being cut off in combined plot (manually set)
  # "-1" break is shown to match boxplot y-axis scale (for spacing)
  scale_x_reverse(expand = expansion(add = c(0, 0.7)), breaks = c(-1))+
  scale_color_manual(name = "Phylum", values = phylumColors)+
  scale_fill_manual(name = "Phylum", values = phylumColors,
                    guide = guide_legend(override.aes = list(size = 10)))+
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        # y-axis elements are coloured white to act as spacers to match boxplot
        axis.title.x = element_text(color = "white"),
        axis.text.x = element_text(color = "white"),
        axis.ticks.x = element_line(color = "white"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        # removing right margin to place dendrogram and boxplot closer
        plot.margin = margin(5.5,0,5.5,5.5, unit = "pt"))+
  theme(legend.title = element_text(size=legendSize, face = "bold", family="Helvetica"),
        legend.text = element_text(size=textSize, face = "bold", family="Helvetica"),
        legend.key = element_blank())

Supp_3.CoordFlip.legend.gTable = get_legend(Supp_3.dendrogram.CoordFlip)
Supp_3.dendrogram.CoordFlip = Supp_3.dendrogram.CoordFlip + 
  theme(legend.position = "none")


Supp_3.boxplot = ggplot()+
  geom_boxplot(data = CHILD_mergedASV_Table_SpeciesColBoxplotMelt, size = 0.1, outlier.size = 0.1,
               aes(x = treeOrder, y = log10(relative_abundance)))+
  ylab("log10(rel. abundance)")+
  theme(plot.background = element_blank(),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        # removing top margin to place dendrogram and boxplot closer
        plot.margin = margin(0,5.5,5.5,5.5, unit = "pt"),
        legend.position = "none")

Supp_3.boxplot.CoordFlip =
  ggplot()+
  geom_boxplot(data = CHILD_mergedASV_Table_SpeciesColBoxplotMelt, size = 0.1, outlier.size = 0.1,
               aes(y = treeOrder, x = log10(relative_abundance)))+
  scale_x_reverse()+
  xlab("log10(rel. abundance)")+
  theme(plot.background = element_blank(),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        # removing left margin to place dendrogram and boxplot closer
        plot.margin = margin(5.5,5.5,5.5,0, unit = "pt"),
        legend.position = "none")

Supp_3.plot_grid.Norm = plot_grid(Supp_3.dendrogram,
                                  Supp_3.boxplot,
                                  ncol = 1, rel_heights = c(5,3))

Supp_3.plot_grid.CoordFlip = plot_grid(Supp_3.dendrogram.CoordFlip,
                             Supp_3.boxplot.CoordFlip, 
                             ncol = 2, rel_widths = c(5,3))

# ---- Supp_4: Alluvial plot, NS1 DCd00, DCV1-3 d01-d17, colour by DCd00, DCd01-d17 presence and phylum ----
graphData = ASV_TableRelAb
colnames(graphData)[colnames(graphData)=="NS1.DCV3d09R"] = "NS1.DCV3d09"
colnames(graphData)[grep("d21", colnames(graphData), fixed = TRUE)] = "X"
graphData = graphData[,c(1:11,grep("NS1\\.DCd00|NS1\\.DCV.d..$",colnames(graphData), fixed = FALSE))]

# subset graphData by sample
graphDataList = alply(c("V1","V2","V3"), 1, function(sampleName){
  sampleName.DCd00 = "NS1.DCd00"
  
  output = graphData[,c(1:11,
                        grep(sampleName.DCd00,colnames(graphData), fixed = TRUE),
                        grep(sampleName,colnames(graphData), fixed = TRUE))]
  
  # order by relative abundance in F sample and NS0.FC samples
  output = output[order(output[,sampleName.DCd00], decreasing = TRUE),]
  output = output[order(rowSums(output[,grep(sampleName,colnames(output), fixed = TRUE)])),]
  
  #order by phylum
  output = output[order(output$BLASTphylum),]
  
  output$presenceGroupColor = interaction(output[,sampleName.DCd00]!=0,
                                          rowSums(output[,grep(sampleName,colnames(output), fixed = TRUE)])!=0)
  
  #order by sample presence
  output = output[order(output$presenceGroupColor),]
  
  # generate percentage data for absence/presence data
  legendData1 = data.frame("day.0" = output[,sampleName.DCd00],
                           "days.1_17" = rowMeans(output[,grep(sampleName,colnames(output), fixed = TRUE)]))
  legendData = aggregate(legendData1, by = list("presenceGroupColor" = interaction(legendData1$day.0!=0,
                                                                                   legendData1$days.1_17!=0)),
                         FUN = function(x) sum(x))
  
  legendData = legendData[,1:2]
  
  legendData$label = paste0(sprintf("%.1f",round(legendData$day.0*100,1)),"%")
  legendData$label[legendData$day.0==0] = ""
  
  legendDataMelt = melt(legendData, measure.vars = 2, variable.name = "Sample", value.name = "Relative_Abundance", stringsAsFactors = FALSE)
  legendDataMelt = legendDataMelt[legendDataMelt$Relative_Abundance!=0,]
  legendDataMelt$day_16S = 0
  
  # generate legend data for absence/presence tile legend
  legendDataTile = aggregate(legendData1, by = list("levelName" = output$presenceGroupColor), FUN = function(x) sum(x)*100)
  
  legendDataTile = legendDataTile[-1,]
  legendDataTile = legendDataTile[match(absencePresence2LegendData$levelName,legendDataTile$levelName),]
  
  legendDataTile$y = c(2,1,0)
  
  legendDataTileMelt = melt(legendDataTile, id.vars = c(1,ncol(legendDataTile)), variable.name = "x")
  legendDataTileMelt$filled = as.numeric(legendDataTileMelt$value>0)
  
  legendDataTileMelt$x = gsub(".","\n",legendDataTileMelt$x, fixed = TRUE)
  legendDataTileMelt$x = gsub("_","-",legendDataTileMelt$x, fixed = TRUE)
  
  legendTile = ggplot()+
    geom_tile(data = legendDataTileMelt[legendDataTileMelt$filled == 1,],
              color = "black",
              size = 0.5,
              aes(x = x, y = y, fill = levelName))+
    # geom_text(data = legendDataMelt[legendDataMelt$filled == 1,],
    #           size = 5.5,
    #           aes(x = x, y = y, label = value))+
    geom_tile(data = legendDataTileMelt[legendDataTileMelt$filled == 0,],
              color = "black",
              fill = "white",
              size = 0.5,
              aes(x = x, y = y))+
    # geom_tile(data = legendDataTileMelt,
    #           color = "black",
    #           height = 0.70,
    #           width = 0.70,
    #           size = 0.25,
    #           aes(x = "", y = y, fill = levelName))+
    scale_x_discrete(position = "top")+   
    scale_y_continuous(expand = expansion(add = 0.1))+
    ggtitle("Identified on:")+
    coord_fixed()+
    legend.day2.theme
  
  #renumber ASVs based on reordering
  output$ASV_Num = seq_along(output$ASV_Num)
  output$ASV_Num = sprintf("%04d", output$ASV_Num)
  
  outputMelt = melt(output, measure.vars = 12:(ncol(output)-1), variable.name = "Sample", value.name = "Relative_Abundance", stringsAsFactors = FALSE)
  outputMelt = outputMelt[outputMelt$Relative_Abundance!=0,]
  
  outputMelt$day_16S = sample_Table$day_16S[match(outputMelt$Sample,sample_Table$sample_16S)]
  
  # return(outputMelt)
  
  #plots
  #colored by sample presence
  Supp_4.plot.alluvial.ap = ggplot() + 
    geom_stratum(data = outputMelt, linetype = "solid",  aes(x = factor(day_16S), stratum = ASV_Num, y = Relative_Abundance*100, fill = presenceGroupColor))+
    geom_errorbar(data = legendDataMelt, stat = "stratum", position = position_nudge(x = -0.3), width = 0.1, size = 1, aes(x = factor(day_16S), stratum = presenceGroupColor, y = Relative_Abundance*100))+
    geom_text(data = legendDataMelt, stat = "stratum", nudge_x = -0.375, angle = 90, size = 6, aes(x = factor(day_16S), stratum = presenceGroupColor, y = Relative_Abundance*100, label = label, vjust = 0, hjust = 0.5))+
    geom_flow(data = outputMelt, color = "black",aes(x = factor(day_16S),stratum = ASV_Num, alluvium = ASV_Num, y = Relative_Abundance*100, fill = presenceGroupColor))+
    scale_x_discrete(name = "Chemostat Culture (day)", expand = expansion(add = c(0.6,0.4)))+
    scale_y_continuous(name = "Relative Abundance (%)")+
    # scale_fill_manual(values = absencePresenceColors3)+
    alluvial.theme_x1_legend0
  
  #colored by phylum
  Supp_4.plot.alluvial.ph = ggplot() +
    geom_stratum(data = outputMelt, linetype = "solid",  aes(x = factor(day_16S), stratum = ASV_Num, y = Relative_Abundance*100, fill = BLASTphylum))+
    geom_flow(data = outputMelt, color = "black",aes(x = factor(day_16S),stratum = ASV_Num, alluvium = ASV_Num, y = Relative_Abundance*100, fill = BLASTphylum))+
    scale_x_discrete(name = "Chemostat Culture (day)", expand = expansion(add = c(0.6,0.4)))+
    scale_y_continuous(name = "Relative Abundance (%)")+
    scale_fill_manual(name = "Phylum", values = phylumColors)+
    alluvial.theme_x1_legend0
  
  return(list(Supp_4.plot.alluvial.ap,Supp_4.plot.alluvial.ph,legendTile))
})
names(graphDataList) = c("V1","V2","V3")

Supp_4A.plot.alluvial.ap = graphDataList[[1]][[1]]
Supp_4B.plot.alluvial.ph = graphDataList[[1]][[2]]

Supp_4C.plot.alluvial.ap = graphDataList[[2]][[1]]
Supp_4D.plot.alluvial.ph = graphDataList[[2]][[2]]

Supp_4E.plot.alluvial.ap = graphDataList[[3]][[1]]
Supp_4F.plot.alluvial.ph = graphDataList[[3]][[2]]

Supp_4.legend.ap = graphDataList[[1]][[3]]
Supp_4.ph.legend.gTable = phylum.legend.gTable

# ---- Supp_5A: heatmap of ASV vs NMR metabolite correlations, NS1DCV1-3d05-d17, without d1 ----
X_16S = data.frame(t(ASV_TableRelAb[,12:(sampleNum+12-1)])) %>% tibble::rownames_to_column(var = "sample_16S")
colnames(X_16S)[2:ncol(X_16S)] = ASV_TableRelAb$ASV_Num

X_16S = X_16S[grep("NS1.DC", X_16S$sample_16S, fixed = TRUE),]

# remove d01
X_16S = X_16S[!grepl("d01", X_16S$sample_16S, fixed = TRUE),]
# remove >d17
X_16S = X_16S[substr(X_16S$sample_16S, nchar(X_16S$sample_16S)-1,nchar(X_16S$sample_16S))<=17,]

X_metabolites = semi_join(NMRmetabolite_Table,X_16S, by = c(sample_16S = "sample_16S"))
X_16S = semi_join(X_16S,NMRmetabolite_Table, by = c(sample_16S = "sample_16S"))
X_16S = X_16S[,c(1,which(colSums(X_16S[2:ncol(X_16S)])!=0)+1)]
X_16S = X_16S[,colSums(X_16S != 0)>1]

X_metabolites = X_metabolites[match(X_16S$sample_16S,X_metabolites$sample_16S),]

cortest = cor(X_metabolites[,7:26],X_16S[,2:ncol(X_16S)])
colnames(cortest) = ASV_TableRelAb$BLASTmatch[match(colnames(cortest),ASV_TableRelAb$ASV_Num)]

Supp_5A.heatplot = Heatmap(cortest, 
                         name = "ASVs vs. Metabolites", #title of legend
                         row_title = "Metabolites", column_title = "ASVs",
                         show_column_names = FALSE,
                         show_row_names = TRUE,
                         row_names_gp = gpar(fontsize = 15), # Text size for row names
                         column_names_gp = gpar(fontsize = 3), # Text size for column names
                         clustering_method_rows = "complete", clustering_method_columns = "complete",
                         row_km = 1, column_km = 6,
                         show_parent_dend_line = FALSE,
                         column_gap = unit(0.005, "npc"),
                         col = circlize::colorRamp2(breaks = c(-1,0,1),
                                                    colors = c("#4575B4", "#FFFFBF", "#D73027")),
                         heatmap_legend_param = list(
                           legend_width = unit(12, "cm"),
                           direction = "horizontal"
                         )
)

Supp_5A.heatplot = ComplexHeatmap::draw(object = Supp_5A.heatplot, heatmap_legend_side = "top")

# ---- Supp_5B,C: scatter plot of selected ASV vs NMR metabolite correlations, NS1 ----
X_16S = data.frame(t(ASV_TableRelAb[,12:(sampleNum+12-1)])) %>% tibble::rownames_to_column(var = "sample_16S")
colnames(X_16S)[2:ncol(X_16S)] = ASV_TableRelAb$ASV_Num

X_16S = X_16S[grep("NS1.DC", X_16S$sample_16S, fixed = TRUE),]

# remove d01
X_16S = X_16S[!grepl("d01", X_16S$sample_16S, fixed = TRUE),]
# remove >d17
X_16S = X_16S[substr(X_16S$sample_16S, nchar(X_16S$sample_16S)-1,nchar(X_16S$sample_16S))<=17,]

X_metabolites = semi_join(NMRmetabolite_Table,X_16S, by = c(sample_16S = "sample_16S"))
X_16S = semi_join(X_16S,NMRmetabolite_Table, by = c(sample_16S = "sample_16S"))
X_16S = X_16S[,c(1,which(colSums(X_16S[2:ncol(X_16S)])!=0)+1)]

X_metabolites = X_metabolites[match(X_16S$sample_16S,X_metabolites$sample_16S),]

cortest = cor(X_metabolites[,7:26],X_16S[,2:ncol(X_16S)])
cortestMelt = melt(cortest)
cortestMelt$abundance = apply(X_16S[,2:ncol(X_16S)], 2, function(x) sum(x))[cortestMelt$Var2]
cortestMelt$conc = apply(X_metabolites[,7:ncol(X_metabolites)], 2, function(x) sum(x))[cortestMelt$Var1]
cortestMelt$abundanceVar = apply(X_16S[,2:ncol(X_16S)], 2, function(x) sd(x))[cortestMelt$Var2]
cortestMelt$concVar = apply(X_metabolites[,7:ncol(X_metabolites)], 2, function(x) sd(x))[cortestMelt$Var1]
cortestMelt$metric = cortestMelt$value * cortestMelt$abundanceVar * cortestMelt$concVar
# cortestMelt$metric = cortestMelt$value * cortestMelt$abundanceVar

graphData1 = data.frame(correlation = "positive", ASV = "ASV_008", X_16S = X_16S$ASV_008*100, metabolite = "Acetate", X_metabolites = X_metabolites$Acetate)
graphData2 = data.frame(correlation = "negative", ASV = "ASV_002", X_16S = X_16S$ASV_002*100, metabolite = "Acetate", X_metabolites = X_metabolites$Acetate)
graphData = rbind(graphData1,graphData2)

slope1 = summary(lm(X_metabolites~X_16S, data = graphData1))$coefficients[2,1]
rsq1 = summary(lm(X_metabolites~X_16S, data = graphData1))$r.squared

graphData1$label = paste0("slope ==",signif(slope1,3),"*';'~R^2 ==", signif(rsq1,3))

slope2 = summary(lm(X_metabolites~X_16S, data = graphData2))$coefficients[2,1]
rsq2 = summary(lm(X_metabolites~X_16S, data = graphData2))$r.squared

graphData2$label = paste0("slope ==",signif(slope2,3),"*';'~R^2 ==", signif(rsq2,3))


# Supp_5B.plot = ggplot()+
#   geom_point(data = graphData, aes(x = X_16S, y = X_metabolites, color = ASV))+
#   facet_wrap(~correlation, scales = "free")

Supp_5B.plot.point = ggplot()+
  geom_point(data = graphData1, size = 5, aes(x = X_16S, y = X_metabolites))+
  xlab(bquote(atop(italic(.(ASV_TableRelAb$BLASTspecies[ASV_TableRelAb$ASV_Num==graphData1$ASV[1]])),"% relative abundance")))+
  ylab(paste0("[",graphData1$metabolite[1],"]"," (mM)"))+
  new_scale_color()+
  geom_smooth(data = graphData1, method = "lm", se = FALSE, aes(x = X_16S, y = X_metabolites, color = label))+
  scale_colour_discrete(name = "",
                        labels = label_parse())+
  point.select.theme_legendtl

Supp_5C.plot.point = ggplot()+
  geom_point(data = graphData2, size = 5, aes(x = X_16S, y = X_metabolites))+   
  xlab(bquote(atop(italic(.(ASV_TableRelAb$BLASTspecies[ASV_TableRelAb$ASV_Num==graphData2$ASV[1]])),"\n% relative abundance")))+   
  ylab(paste0("[",graphData2$metabolite[1],"]"," (mM)"))+
  new_scale_color()+
  geom_smooth(data = graphData2, method = "lm", se = FALSE, aes(x = X_16S, y = X_metabolites, color = label))+
  scale_colour_discrete(name = "",
                        labels = label_parse())+
  point.select.theme_legendtr

# ==== Figure 2 ====
# ---- A,B,C,D: Alluvial plot, NS0 and NS1 samples F, FCd01-d17, colour by F, FCd01-d17 presence and phylum----
graphData = ASV_TableRelAb
# graphData = graphData[,c(1:11,grep(".F",colnames(graphData), fixed = TRUE))]
graphData = graphData[,c(1:11,which(colnames(graphData) %in% c("NS0.F","NS0.FCd01", "NS0.FCd05", "NS0.FCd09", "NS0.FCd13", "NS0.FCd17",
                                                               "NS1.F","NS1.FCd01", "NS1.FCd05", "NS1.FCd09", "NS1.FCd13", "NS1.FCd17")))]

# subset graphData by sample
graphDataList = alply(c("NS0","NS1"), 1, function(sampleName){
  sampleName.F = paste0(sampleName,".F")
  
  output = graphData[,c(1:11,grep(sampleName,colnames(graphData), fixed = TRUE))]
  
  # order by relative abundance in F sample and NS0.FC samples
  output = output[order(output[,sampleName.F], decreasing = TRUE),]
  output = output[order(rowSums(output[,grep("FC",colnames(output), fixed = TRUE)])),]
  
  #order by phylum
  output = output[order(output$BLASTphylum),]
  
  output$presenceGroupColor = interaction(output[,sampleName.F]!=0,
                                          rowSums(output[,grep("FC",colnames(output), fixed = TRUE)])!=0)
  
  #order by sample presence
  output = output[order(output$presenceGroupColor),]
  
  # generate percentage data for absence/presence data
  legendData1 = data.frame("day.0" = output[,sampleName.F],
                           "days.1_17" = rowMeans(output[,grep("FC",colnames(output), fixed = TRUE)]))
  legendData = aggregate(legendData1, by = list("presenceGroupColor" = interaction(legendData1$day.0!=0,
                                                                                   legendData1$days.1_17!=0)),
                         FUN = function(x) sum(x))
  
  legendData = legendData[,1:2]
  
  legendData$label = paste0(sprintf("%.1f",round(legendData$day.0*100,1)),"%")
  legendData$label[legendData$day.0==0] = ""
  
  legendDataMelt = melt(legendData, measure.vars = 2, variable.name = "Sample", value.name = "Relative_Abundance", stringsAsFactors = FALSE)
  legendDataMelt = legendDataMelt[legendDataMelt$Relative_Abundance!=0,]
  legendDataMelt$day_16S = 0
  
  # generate legend data for absence/presence tile legend
  legendDataTile = aggregate(legendData1, by = list("levelName" = output$presenceGroupColor), FUN = function(x) sum(x)*100)
  
  legendDataTile = legendDataTile[-1,]
  legendDataTile = legendDataTile[match(absencePresence2LegendData$levelName,legendDataTile$levelName),]
  
  legendDataTile$y = c(2,1,0)
  
  legendDataTileMelt = melt(legendDataTile, id.vars = c(1,ncol(legendDataTile)), variable.name = "x")
  legendDataTileMelt$filled = as.numeric(legendDataTileMelt$value>0)
 
  legendDataTileMelt$x = gsub(".","\n",legendDataTileMelt$x, fixed = TRUE)
  legendDataTileMelt$x = gsub("_","-",legendDataTileMelt$x, fixed = TRUE)
  
  legendTile = ggplot()+
    geom_tile(data = legendDataTileMelt[legendDataTileMelt$filled == 1,],
              color = "black",
              size = 0.5,
              aes(x = x, y = y, fill = levelName))+
    # geom_text(data = legendDataMelt[legendDataMelt$filled == 1,],
    #           size = 5.5,
    #           aes(x = x, y = y, label = value))+
    geom_tile(data = legendDataTileMelt[legendDataTileMelt$filled == 0,],
              color = "black",
              fill = "white",
              size = 0.5,
              aes(x = x, y = y))+
    # geom_tile(data = legendDataTileMelt,
    #           color = "black",
    #           height = 0.70,
    #           width = 0.70,
    #           size = 0.25,
    #           aes(x = "", y = y, fill = levelName))+
    scale_x_discrete(position = "top")+   
    scale_y_continuous(expand = expansion(add = 0.1))+
    ggtitle("Identified on:")+
    coord_fixed()+
    legend.day2.theme
  
  #renumber ASVs based on reordering
  output$ASV_Num = seq_along(output$ASV_Num)
  output$ASV_Num = sprintf("%04d", output$ASV_Num)
  
  outputMelt = melt(output, measure.vars = 12:(ncol(output)-1), variable.name = "Sample", value.name = "Relative_Abundance", stringsAsFactors = FALSE)
  outputMelt = outputMelt[outputMelt$Relative_Abundance!=0,]
  
  outputMelt$day_16S = sample_Table$day_16S[match(outputMelt$Sample,sample_Table$sample_16S)]
  
  # return(outputMelt)
  
  #plots
  #colored by sample presence
  Fig2.plot.alluvial.ap = ggplot() + 
    geom_stratum(data = outputMelt, linetype = "solid",  aes(x = factor(day_16S), stratum = ASV_Num, y = Relative_Abundance*100, fill = presenceGroupColor))+
    geom_errorbar(data = legendDataMelt, stat = "stratum", position = position_nudge(x = -0.3), width = 0.1, size = 1, aes(x = factor(day_16S), stratum = presenceGroupColor, y = Relative_Abundance*100))+
    geom_text(data = legendDataMelt, stat = "stratum", nudge_x = -0.375, angle = 90, size = 6, aes(x = factor(day_16S), stratum = presenceGroupColor, y = Relative_Abundance*100, label = label, vjust = 0, hjust = 0.5))+
    geom_flow(data = outputMelt, color = "black",aes(x = factor(day_16S),stratum = ASV_Num, alluvium = ASV_Num, y = Relative_Abundance*100, fill = presenceGroupColor))+
    scale_x_discrete(name = "Chemostat Culture (day)", expand = expansion(add = c(0.6,0.4)))+
    scale_y_continuous(name = "Relative Abundance (%)")+
    # scale_fill_manual(values = absencePresenceColors3)+
    alluvial.theme_x1_legend0
  
  #colored by phylum
  Fig2.plot.alluvial.ph = ggplot() +
    geom_stratum(data = outputMelt, linetype = "solid",  aes(x = factor(day_16S), stratum = ASV_Num, y = Relative_Abundance*100, fill = BLASTphylum))+
    geom_flow(data = outputMelt, color = "black",aes(x = factor(day_16S),stratum = ASV_Num, alluvium = ASV_Num, y = Relative_Abundance*100, fill = BLASTphylum))+
    scale_x_discrete(name = "Chemostat Culture (day)", expand = expansion(add = c(0.6,0.4)))+
    scale_y_continuous(name = "Relative Abundance (%)")+
    scale_fill_manual(name = "Phylum", values = phylumColors)+
    alluvial.theme_x1_legend0
  
  return(list(Fig2.plot.alluvial.ap,Fig2.plot.alluvial.ph,legendTile))
})
names(graphDataList) = c("NS0","NS1")

Fig2A.plot.alluvial.ap = graphDataList[[1]][[1]]
Fig2B.plot.alluvial.ph = graphDataList[[1]][[2]]

Fig2C.plot.alluvial.ap = graphDataList[[2]][[1]]
Fig2D.plot.alluvial.ph = graphDataList[[2]][[2]]

Fig2.legend.ap = graphDataList[[1]][[3]]
Fig2.ph.legend.gTable = phylum.legend.gTable

# ==== Figure 3 ====
# ---- A,B,C,D: Alluvial plot, NS0 and NS1 samples F, FCd17, colour by F, FCd17 presence and phylum ----
graphData = ASV_TableRelAb
# graphData = graphData[,c(1:11,grep(".F",colnames(graphData), fixed = TRUE))]
graphData = graphData[,c(1:11,which(colnames(graphData) %in% c("NS0.F","NS0.FCd17",
                                                               "NS1.F","NS1.FCd17")))]

# subset graphData by sample
graphDataList = alply(c("NS0","NS1"), 1, function(sampleName){
  sampleName.F = paste0(sampleName,".F")
  sampleName.FCd17 = paste0(sampleName,".FCd17")
  
  output = graphData[,c(1:11,grep(sampleName,colnames(graphData), fixed = TRUE))]
  
  # order by relative abundance in F sample and FCd17 samples
  output = output[order(output[,sampleName.F], decreasing = TRUE),]
  output = output[order(output[,sampleName.FCd17], decreasing = TRUE),]
  
  #order by phylum
  output = output[order(output$BLASTphylum),]
  
  output$presenceGroupColor = interaction(output[,sampleName.F]!=0,
                                          output[,sampleName.FCd17]!=0)
  
  #order by sample presence
  output = output[order(output$presenceGroupColor),]
  
  # generate percentage data for absence/presence data
  legendData1 = data.frame(output[,c(sampleName.F,sampleName.FCd17)])
  
  # return(legendData1)
  
  legendData = aggregate(legendData1, by = list("presenceGroupColor" = interaction(legendData1[,sampleName.F]!=0,
                                                                                   legendData1[,sampleName.FCd17]!=0)),
                         FUN = function(x) sum(x))
  
  legendData = legendData[,1:2]
  
  legendData$label = paste0(sprintf("%.1f",round(legendData[,sampleName.F]*100,1)),"%")
  legendData$label[legendData[,sampleName.F]==0] = ""
  
  legendDataMelt = melt(legendData, measure.vars = 2, variable.name = "Sample", value.name = "Relative_Abundance", stringsAsFactors = FALSE)
  legendDataMelt = legendDataMelt[legendDataMelt$Relative_Abundance!=0,]
  legendDataMelt$day_16S = 0
  legendDataMelt$Sample = sampleNameLevels.3lines[as.character(legendDataMelt$Sample)]
  # return(legendDataMelt)
  
  # generate legend data for absence/presence tile legend
  legendDataTile = aggregate(legendData1, by = list("levelName" = output$presenceGroupColor), FUN = function(x) sum(x)*100)
  
  legendDataTile = legendDataTile[-1,]
  legendDataTile = legendDataTile[match(absencePresence2LegendData$levelName,legendDataTile$levelName),]
  
  legendDataTile$y = c(2,1,0)
  
  # return(legendDataTile)
  
  legendDataTileMelt = melt(legendDataTile, id.vars = c(1,ncol(legendDataTile)), variable.name = "x")
  legendDataTileMelt$filled = as.numeric(legendDataTileMelt$value>0)
  
  levels(legendDataTileMelt$x) = sampleNameLevels.3lines[levels(legendDataTileMelt$x)]
  
  # legendDataTileMelt$x = gsub(".","\n",legendDataTileMelt$x, fixed = TRUE)
  # legendDataTileMelt$x = gsub("_","-",legendDataTileMelt$x, fixed = TRUE)
  # return(legendDataTileMelt)
  
  legendTile = ggplot()+
    geom_tile(data = legendDataTileMelt[legendDataTileMelt$filled == 1,],
              color = "black",
              size = 0.5,
              aes(x = x, y = y, fill = levelName))+
    # geom_text(data = legendDataMelt[legendDataMelt$filled == 1,],
    #           size = 5.5,
    #           aes(x = x, y = y, label = value))+
    geom_tile(data = legendDataTileMelt[legendDataTileMelt$filled == 0,],
              color = "black",
              fill = "white",
              size = 0.5,
              aes(x = x, y = y))+
    # geom_tile(data = legendDataTileMelt,
    #           color = "black",
    #           height = 0.70,
    #           width = 0.70,
    #           size = 0.25,
    #           aes(x = "", y = y, fill = levelName))+
    scale_x_discrete(position = "top")+   
    scale_y_continuous(expand = expansion(add = 0.1))+
    ggtitle("Identified in:")+
    coord_fixed()+
    legend.day2.theme
  
  #renumber ASVs based on reordering
  output$ASV_Num = seq_along(output$ASV_Num)
  output$ASV_Num = sprintf("%04d", output$ASV_Num)
  
  outputMelt = melt(output, measure.vars = 12:(ncol(output)-1), variable.name = "Sample", value.name = "Relative_Abundance", stringsAsFactors = FALSE)
  outputMelt = outputMelt[outputMelt$Relative_Abundance!=0,]
  
  outputMelt$day_16S = sample_Table$day_16S[match(outputMelt$Sample,sample_Table$sample_16S)]
  
  levels(outputMelt$Sample) = sampleNameLevels.3lines[levels(outputMelt$Sample)]
  
  # return(outputMelt)
  
  #plots
  #colored by sample presence
  Fig3.plot.alluvial.ap = ggplot() + 
    geom_stratum(data = outputMelt, linetype = "solid",  aes(x = Sample, stratum = ASV_Num, y = Relative_Abundance*100, fill = presenceGroupColor))+
    geom_errorbar(data = legendDataMelt, stat = "stratum", position = position_nudge(x = -0.3), width = 0.1, size = 1, aes(x = Sample, stratum = presenceGroupColor, y = Relative_Abundance*100))+
    geom_text(data = legendDataMelt, stat = "stratum", nudge_x = -0.375, angle = 90, size = 6, aes(x = Sample, stratum = presenceGroupColor, y = Relative_Abundance*100, label = label, vjust = 0, hjust = 0.5))+
    geom_flow(data = outputMelt, color = "black",aes(x = Sample,stratum = ASV_Num, alluvium = ASV_Num, y = Relative_Abundance*100, fill = presenceGroupColor))+
    scale_x_discrete(name = "Chemostat Culture (day)", expand = expansion(add = c(0.6,0.4)))+
    scale_y_continuous(name = "Relative Abundance (%)")+
    # scale_fill_manual(values = absencePresenceColors3)+
    alluvial.theme_x1_legend0
  
  #colored by phylum
  Fig3.plot.alluvial.ph = ggplot() +
    geom_stratum(data = outputMelt, linetype = "solid",  aes(x = Sample, stratum = ASV_Num, y = Relative_Abundance*100, fill = BLASTphylum))+
    geom_flow(data = outputMelt, color = "black",aes(x = Sample,stratum = ASV_Num, alluvium = ASV_Num, y = Relative_Abundance*100, fill = BLASTphylum))+
    scale_x_discrete(name = "Chemostat Culture (day)")+
    scale_y_continuous(name = "Relative Abundance (%)")+
    scale_fill_manual(name = "Phylum", values = phylumColors)+
    alluvial.theme_x1_legend0
  
  return(list(Fig3.plot.alluvial.ap,Fig3.plot.alluvial.ph,legendTile))
})
names(graphDataList) = c("NS0","NS1")

Fig3A.plot.alluvial.ap = graphDataList[[1]][[1]]
Fig3B.plot.alluvial.ph = graphDataList[[1]][[2]]

Fig3C.plot.alluvial.ap = graphDataList[[2]][[1]]
Fig3D.plot.alluvial.ph = graphDataList[[2]][[2]]

Fig3.legend.ap = graphDataList[[1]][[3]]
Fig3.legend.ap = Fig3.legend.ap + scale_x_discrete(position = "top", labels=c("Fecal", "Fecal\nChemostat"))

Fig3.ph.legend.gTable = phylum.legend.gTable

# ==== Figure 4 ====
# ---- A,B,C,D: Alluvial plot, NS0 and NS1 samples F, FCd17, DCd00, colour by F, FCd17, DCd00 presence, and phylum ----
graphData = ASV_TableRelAb
graphData = graphData[,c(1:11,which(colnames(graphData) %in% c("NS0.F","NS0.DCd00","NS0.FCd17",
                                                               "NS1.F","NS1.DCd00","NS1.FCd17")))]

# subset graphData by sample
graphDataList = alply(c("NS0","NS1"), 1, function(sampleName){
  sampleName.F = paste0(sampleName,".F")
  sampleName.DCd00 = paste0(sampleName,".DCd00")
  sampleName.FCd17 = paste0(sampleName,".FCd17")
  
  output = graphData[,c(1:11,grep(sampleName,colnames(graphData), fixed = TRUE))]
  # return(output)
  #order by phylum
  # output = output[order(output[,sampleName.FCd17]),]
  output = output[order(output$BLASTphylum),]
  
  output[,sampleName.DCd00][output[,sampleName.DCd00]!=0] = 1/sum(output[,sampleName.DCd00]!=0)
  
  output$presenceGroupColor = interaction(output[,sampleName.F]!=0,
                                          output[,sampleName.DCd00]!=0
  )
  
  output$presenceGroupColor2 = interaction(output[,sampleName.F]!=0,
                                          output[,sampleName.DCd00]!=0,
                                          output[,sampleName.FCd17]!=0
  )
  
  #order by sample presence
  output = output[order(output$presenceGroupColor2),]
  output = output[order(output$presenceGroupColor),]
  
  # generate percentage data for absence/presence data
  legendData1 = data.frame(output[,c(sampleName.F,sampleName.DCd00)])
  legendData = aggregate(legendData1, by = list("presenceGroupColor" = interaction(output[,sampleName.F]!=0,
                                                                                   output[,sampleName.DCd00]!=0)),
                         FUN = function(x) sum(x))
  # return(legendData)
  legendData = legendData[,1:2]
  
  # return(legendData)
  
  legendData$label = paste0(sprintf("%.1f",round(legendData[,sampleName.F]*100,1)),"%")
  legendData$label[legendData[,sampleName.F]==0] = ""
  
  legendDataMelt = melt(legendData, measure.vars = 2, variable.name = "Sample", value.name = "Relative_Abundance", stringsAsFactors = FALSE)
  legendDataMelt = legendDataMelt[legendDataMelt$Relative_Abundance!=0,]
  legendDataMelt$day_16S = 0
  
  
  legendDataMelt$Sample = as.character(legendDataMelt$Sample)
  legendDataMelt$Sample = sampleNameLevels.3lines[legendDataMelt$Sample]
  # return(legendDataMelt)
  
  
  # generate legend data for absence/presence tile legend
  legendData2 = data.frame(output[,c(sampleName.F,sampleName.DCd00,sampleName.FCd17)])
  legendDataTile = aggregate(legendData2, by = list("levelName" = output$presenceGroupColor2), FUN = function(x) sum(x)*100)
  # return(legendDataTile)
  legendDataTile = legendDataTile[-1,]
  legendDataTile = legendDataTile[match(absencePresence3LegendData$levelName,legendDataTile$levelName),]
  
  legendDataTile$y = c(6,5,4,3,2,1,0)
  
  legendDataTileMelt = melt(legendDataTile, id.vars = c(1,ncol(legendDataTile)), variable.name = "x")
  legendDataTileMelt$x = as.character(legendDataTileMelt$x)
  legendDataTileMelt$x = sampleNameLevels.3lines[legendDataTileMelt$x]
  legendDataTileMelt$x = factor(legendDataTileMelt$x)
  legendDataTileMelt$x = factor(as.character(legendDataTileMelt$x),levels(legendDataTileMelt$x)[c(1,3,2)])
  
  legendDataTileMelt$filled = as.numeric(legendDataTileMelt$value>0)
  
  # legendDataTileMelt$x = gsub(".","\n",legendDataTileMelt$x, fixed = TRUE)
  # legendDataTileMelt$x = gsub("_","-",legendDataTileMelt$x, fixed = TRUE)
  # return(legendDataTileMelt)
  legendTile = ggplot()+
    geom_tile(data = legendDataTileMelt[legendDataTileMelt$filled == 1,],
              color = "black",
              size = 0.5,
              aes(x = x, y = y, fill = levelName))+
    # geom_text(data = legendDataMelt[legendDataMelt$filled == 1,],
    #           size = 5.5,
    #           aes(x = x, y = y, label = value))+
    geom_tile(data = legendDataTileMelt[legendDataTileMelt$filled == 0,],
              color = "black",
              fill = "white",
              size = 0.5,
              aes(x = x, y = y))+
    # geom_tile(data = legendDataTileMelt,
    #           color = "black",
    #           height = 0.70,
    #           width = 0.70,
    #           size = 0.25,
    #           aes(x = "", y = y, fill = levelName))+
    scale_x_discrete(position = "top")+   
    scale_y_continuous(expand = expansion(add = 0.1))+
    ggtitle("Identified in:")+
    coord_fixed()+
    legend.chemostat3.theme
  
  #renumber ASVs based on reordering
  output$ASV_Num = seq_along(output$ASV_Num)
  output$ASV_Num = sprintf("%04d", output$ASV_Num)
  
  outputMelt = melt(output, measure.vars = 12:(ncol(output)-2), variable.name = "Sample", value.name = "Relative_Abundance", stringsAsFactors = FALSE)
  outputMelt = outputMelt[outputMelt$Relative_Abundance!=0,]
  
  outputMelt$day_16S = sample_Table$day_16S[match(outputMelt$Sample,sample_Table$sample_16S)]
  
  levels(outputMelt$Sample) = sampleNameLevels.3lines[levels(outputMelt$Sample)]
  outputMelt$Sample = factor(outputMelt$Sample, levels = levels(outputMelt$Sample)[c(1,3,2)])
  # return(outputMelt)
  
  #plots
  #colored by sample presence
  Fig4.plot.alluvial.ap = ggplot() + 
    geom_stratum(data = outputMelt, linetype = "solid",  aes(x = Sample, stratum = ASV_Num, y = Relative_Abundance*100, fill = presenceGroupColor2))+
    geom_errorbar(data = legendDataMelt, stat = "stratum", position = position_nudge(x = -0.3), width = 0.1, size = 1, aes(x = Sample, stratum = presenceGroupColor, y = Relative_Abundance*100))+
    geom_text(data = legendDataMelt, stat = "stratum", nudge_x = -0.375, angle = 90, size = 6, aes(x = Sample, stratum = presenceGroupColor, y = Relative_Abundance*100, label = label, vjust = 0, hjust = 0.5))+
    geom_flow(data = outputMelt, color = "black",aes(x = Sample,stratum = ASV_Num, alluvium = ASV_Num, y = Relative_Abundance*100, fill = presenceGroupColor2))+
    scale_x_discrete(name = "Chemostat Culture (day)", expand = expansion(add = c(0.6,0.6)))+
    scale_y_continuous(name = "Relative Abundance (%)")+
    # scale_fill_manual(values = absencePresenceColors3)+
    alluvial.theme_x1_legend0
  
  #colored by phylum
  Fig4.plot.alluvial.ph = ggplot() +
    geom_stratum(data = outputMelt, linetype = "solid",  aes(x = Sample, stratum = ASV_Num, y = Relative_Abundance*100, fill = BLASTphylum))+
    geom_flow(data = outputMelt, color = "black",aes(x = Sample,stratum = ASV_Num, alluvium = ASV_Num, y = Relative_Abundance*100, fill = BLASTphylum))+
    scale_x_discrete(name = "Chemostat Culture (day)", expand = expansion(add = c(0.6,0.6)))+
    scale_y_continuous(name = "Relative Abundance (%)")+
    scale_fill_manual(name = "Phylum", values = phylumColors)+
    alluvial.theme_x1_legend0
  
  return(list(Fig4.plot.alluvial.ap,Fig4.plot.alluvial.ph,legendTile))
})

names(graphDataList) = c("NS0","NS1")

Fig4A.plot.alluvial.ap = graphDataList[[1]][[1]]
Fig4B.plot.alluvial.ph = graphDataList[[1]][[2]]

Fig4C.plot.alluvial.ap = graphDataList[[2]][[1]]
Fig4D.plot.alluvial.ph = graphDataList[[2]][[2]]

Fig4.legend.ap = graphDataList[[1]][[3]]
Fig4.legend.ap = Fig4.legend.ap + scale_x_discrete(position = "top", labels=c("Fecal", "Isolate\nLibrary", "Fecal\nChemostat"))

Fig4.ph.legend.gTable = phylum.legend.gTable

# ---- E: Histogram, percentage of species with binned prevalences in CHILD study samples, faceted by in CHILD vs. in CHILD and isolates ----
# load CHILD 12 month ASV table
CHILD_12M_ASV_Table = read.delim(paste0(path,"/CHILD_12M ASV_Table.txt"),sep = "\t", stringsAsFactors = FALSE)

CHILD_ASV_seqtab = as.matrix(t(CHILD_12M_ASV_Table[,12:ncol(CHILD_12M_ASV_Table)]))
CHILDSampleCount = nrow(CHILD_ASV_seqtab)
colnames(CHILD_ASV_seqtab) = CHILD_12M_ASV_Table$sequences

graphData = ASV_Table
# select samples
graphData = graphData[,c(1:11,which(colnames(graphData) %in% c("NS0.DCd00","NS1.DCd00")))]
# remove ASVs not in DC
graphData = graphData[rowSums(graphData[,12:ncol(graphData)])>0,]
graphData = data.frame(graphData[,c(1:11)],"Isolates_presence" = (rowSums(graphData[,12:ncol(graphData)])>0)+0)
graphData_seqtab = as.matrix(t(graphData[,12:ncol(graphData)]))
colnames(graphData_seqtab) = graphData$sequences
rownames(graphData_seqtab) = "Isolates_presence"

# merge CHILD study and isolates seqtabs
merged_seqtab = mergeSequenceTables(CHILD_ASV_seqtab, graphData_seqtab)
merged_seqtab_collapsed = collapseNoMismatch2(merged_seqtab,verbose = TRUE)

# reformat as ASV table
mergedASV_Table = data.frame(ASV_Num = seq_along(colnames(merged_seqtab_collapsed)), sequences = colnames(merged_seqtab_collapsed),t(merged_seqtab_collapsed), stringsAsFactors = FALSE)
row.names(mergedASV_Table) = NULL

mergedASV_Table$ASV_Num = sprintf("%04d", mergedASV_Table$ASV_Num)
mergedASV_Table$ASV_Num = paste0("ASV_", mergedASV_Table$ASV_Num)

seq2fasta(mergedASV_Table$ASV_Num, mergedASV_Table$sequences, "CHILD_prevalence", file.path(path,"ASV FASTAs"))

blastResults = BLAST_16S(fastafileName = paste0(path,"/ASV FASTAs/CHILD_prevalence.txt"),
                         remote = FALSE,
                         wait = TRUE,
                         topHit = TRUE,
                         taxonomy = TRUE)

mergedASV_Table = data.frame(mergedASV_Table[,1:2],
                             blastResults[match(mergedASV_Table$ASV_Num, blastResults$query),c("BLASTsuperkingdom","BLASTphylum","BLASTclass","BLASTorder","BLASTfamily","BLASTgenus","BLASTspecies","BLASTmatch")],
                             BLASTper_id = blastResults$per_id[match(mergedASV_Table$ASV_Num,blastResults$query)],
                             mergedASV_Table[,3:ncol(mergedASV_Table)])

write.table(mergedASV_Table,file = file.path(path,"CHILD_prevalence ASV_Table.txt"),sep = "\t", row.names = FALSE)
CHILD_mergedASV_Table = read.delim(file.path(path,"CHILD_prevalence ASV_Table.txt"))

# subset out ASVs in isolates
CHILD_mergedASV_Table_isoSubset = CHILD_mergedASV_Table[CHILD_mergedASV_Table$Isolates_presence>0,]
CHILD_mergedASV_Table_nonisoSubset = CHILD_mergedASV_Table[CHILD_mergedASV_Table$Isolates_presence==0,]

# load NS0 isolates table
NS0_F_C_Isolates = readxl::read_xlsx(path = file.path(path,"input/Supplementary_table3.xlsx"))

# load NS1 isolates table
NS1_F_C_Isolates = readxl::read_xlsx(path = file.path(path,"input/Supplementary_table4.xlsx"))

CheckTableNS0 = data.frame(seqName = NS0_F_C_Isolates$`Bacterial strain ID`,
                           source = "NS0",
                           isoNum = paste0("NS0",1:nrow(NS0_F_C_Isolates)),
                           seq = as.character(reverseComplement(DNAStringSet(NS0_F_C_Isolates$`16S sequence (V3-V6 region)`))),
                           NS0count = 1,
                           NS1count = 0,
                           count = 1)

CheckTableNS0$seq[CheckTableNS0$seqName =="T1D_NS_0 10 FAA AN"] = as.character(reverseComplement(DNAString(CheckTableNS0$seq[CheckTableNS0$seqName =="T1D_NS_0 10 FAA AN"])))

CheckTableNS1 = data.frame(seqName = NS1_F_C_Isolates$`Bacterial strain ID`,
                           source = "NS1",
                           isoNum = paste0("NS1",1:nrow(NS1_F_C_Isolates)),
                           seq = as.character(reverseComplement(DNAStringSet(NS1_F_C_Isolates$`16S sequence (V3-V6 region)`))),
                           NS0count = 0,
                           NS1count = 1,
                           count = 1)

CheckTableSanger = rbind(CheckTableNS0,CheckTableNS1)
CheckTableSanger$seq = gsub("-","",CheckTableSanger$seq,fixed = TRUE)

CheckTableIsoSubset = data.frame(seqName = CHILD_mergedASV_Table_isoSubset$ASV_Num,
                                 source = "iso",
                                 isoNum = CHILD_mergedASV_Table_isoSubset$ASV_Num,
                                 seq = CHILD_mergedASV_Table_isoSubset$sequences,
                                 NS0count = CHILD_mergedASV_Table_isoSubset$NS0.DCd00,
                                 NS1count = CHILD_mergedASV_Table_isoSubset$NS1.DCd00,
                                 count = CHILD_mergedASV_Table_isoSubset$NS0.DCd00+CHILD_mergedASV_Table_isoSubset$NS1.DCd00)

CheckTable = rbind(CheckTableSanger,CheckTableIsoSubset)

seq2fasta(CheckTable$isoNum, CheckTable$seq, "CheckTable", file.path(path,"ASV FASTAs"))

blastResults = BLAST_16S(fastafileName = paste0(path,"/ASV FASTAs/CheckTable.txt"),
                         remote = FALSE,
                         wait = TRUE,
                         topHit = TRUE,
                         taxonomy = TRUE)

CheckTable = data.frame(CheckTable[,1:4],
                        blastResults[match(CheckTable$isoNum, blastResults$query),c("BLASTsuperkingdom","BLASTphylum","BLASTclass","BLASTorder","BLASTfamily","BLASTgenus","BLASTspecies","BLASTmatch")],
                        BLASTper_id = blastResults$per_id[match(CheckTable$isoNum,blastResults$query)],
                        CheckTable[,5:ncol(CheckTable)])


CheckTableClust = DECIPHER_clusterASVs(CheckTable$seq, testBounds = TRUE,
                                       lowerBound = 0.80, upperBound = 1, testIncrement = 0.001)

CheckTableClust.sanger_not_iso_Count = apply(CheckTableClust,2, function(col){
  sanger = col[CheckTable$source %in% c("NS0","NS1")]
  iso = col[CheckTable$source %in% c("iso")]
  
  sum(!(sanger %in% iso))
})

CheckTableClust.sanger_not_iso_Species = apply(CheckTableClust,2, function(col){
  sanger = CheckTable$source %in% c("NS0","NS1")
  iso = CheckTable$source %in% c("iso")
  
  sangerCol = col[sanger]
  isoCol = col[iso]
  out = rep(NA, nrow(CheckTable))
  out[col %in% (sangerCol[!(sangerCol %in% isoCol)])] = 
    CheckTable$BLASTspecies[col %in% (sangerCol[!(sangerCol %in% isoCol)])]
  out
})

CheckTableClust.iso_not_sanger_Count = apply(CheckTableClust,2, function(col){
  sanger = CheckTable$source %in% c("NS0","NS1")
  iso = CheckTable$source %in% c("iso")
  
  sangerCol = col[sanger]
  isoCol = col[iso]
  
  sum(CheckTable$count[col %in% (isoCol[!(isoCol %in% sangerCol)])])
  
})

# filter isoSubset ASVs: remove ASV as spurious if they don't cluster with an isolate Sanger sequence at 0.941 identity
# filterClustCutoffCol = CheckTableClust$cluster_0_082
filterClustCutoffCol = CheckTableClust$cluster_0_059

filterClustNames = CheckTable$seqName[filterClustCutoffCol %in% 
                                        (filterClustCutoffCol[CheckTable$source %in% c("iso")][!(filterClustCutoffCol[CheckTable$source %in% c("iso")] %in% 
                                                                                                   filterClustCutoffCol[CheckTable$source %in% c("NS0","NS1")])])]

CHILD_mergedASV_Table_isoSubsetFiltered = CHILD_mergedASV_Table_isoSubset[!(CHILD_mergedASV_Table_isoSubset$ASV_Num %in%filterClustNames),]

# recombine CHILD_mergedASV_Table
CHILD_mergedASV_Table = rbind(CHILD_mergedASV_Table_isoSubsetFiltered,CHILD_mergedASV_Table_nonisoSubset)

# collapse to species level by BLASTspecies names,
# keeping metadata from ASV present in NS0 and NS1 isolates and most prevalent in CHILD study samples as "representative"
CHILD_mergedASV_Table_SpeciesCol = CHILD_mergedASV_Table[order(CHILD_mergedASV_Table$CHILD_prevalence, decreasing = TRUE),]
CHILD_mergedASV_Table_SpeciesCol = CHILD_mergedASV_Table_SpeciesCol[order(CHILD_mergedASV_Table_SpeciesCol$Isolates_presence, decreasing = TRUE),]

CHILD_mergedASV_Table_SpeciesCol = data.frame(CHILD_mergedASV_Table_SpeciesCol[which(!duplicated(CHILD_mergedASV_Table_SpeciesCol$BLASTspecies)),1:11],
                                              rowsum(CHILD_mergedASV_Table_SpeciesCol[12:ncol(CHILD_mergedASV_Table_SpeciesCol)],
                                                     group = CHILD_mergedASV_Table_SpeciesCol$BLASTspecies,
                                                     reorder = FALSE)
)

# collapse to species level by BLASTspecies names, keeping metadata from ASV with the highest counts as "representative"
# seqCountOrder = order(rowSums(CHILD_mergedASV_Table[12:ncol(CHILD_mergedASV_Table)]), decreasing = TRUE)

# CHILD_mergedASV_Table_SpeciesCol = CHILD_mergedASV_Table[seqCountOrder,]
# CHILD_mergedASV_Table_SpeciesCol = data.frame(CHILD_mergedASV_Table_SpeciesCol[which(!duplicated(CHILD_mergedASV_Table_SpeciesCol$BLASTspecies)),1:11],
#                                               rowsum(CHILD_mergedASV_Table_SpeciesCol[12:ncol(CHILD_mergedASV_Table_SpeciesCol)],
#                                                      group = CHILD_mergedASV_Table_SpeciesCol$BLASTspecies,
#                                                      reorder = FALSE)
# )

# recalculate prevalence in the CHILD study samples for species collapsed ASVs
# prevalence in the CHILD study samples = # of samples where ASV is >0 / # of samples
CHILD_mergedASV_Table_SpeciesCol$CHILD_prevalence = 
  rowSums(CHILD_mergedASV_Table_SpeciesCol[,grep(pattern = "X\\d",x = colnames(CHILD_mergedASV_Table_SpeciesCol), fixed = FALSE)]>0) / 
  sum(grepl(pattern = "X\\d",x = colnames(CHILD_mergedASV_Table_SpeciesCol), fixed = FALSE))

# creating graph data with 2 groups: all species vs. species in isolates
graphData1 = CHILD_mergedASV_Table_SpeciesCol[CHILD_mergedASV_Table_SpeciesCol$Isolates_presence>0,]
graphData1$group = "In CHILD study and isolates library"
graphData2 = CHILD_mergedASV_Table_SpeciesCol
graphData2$group = "In CHILD study"
graphData = rbind(graphData1,graphData2)

Fig4E.histogram = ggplot()+
  geom_histogram(data = graphData, bins = 100, position = "dodge", aes(x = CHILD_prevalence,y = after_stat(density)*after_stat(width)*100, fill = group))+
  geom_density(data = graphData,bounds = c(0, 1), aes(x = CHILD_prevalence, color = group))+
  scale_color_discrete(guide = "none")+
  scale_fill_discrete(name = "Species presence")+
  scale_x_continuous(name = "Prevalence in CHILD study samples")+
  scale_y_continuous(name = "Species in group (%)")+
  histogram.theme

Fig4E.histogramPseudoLog.sigma1=
ggplot()+
  geom_histogram(data = graphData, bins = 100, position = "dodge", aes(x = CHILD_prevalence,y = after_stat(density)*after_stat(width)*100, fill = group))+
  geom_density(data = graphData,bounds = c(0, 1), trim = TRUE, adjust = 5, aes(x = CHILD_prevalence, color = group))+
  scale_color_discrete(guide = "none")+
  scale_fill_discrete(name = "Species presence")+
  scale_x_continuous(name = "Prevalence in CHILD study samples")+
  scale_y_continuous(name = "Species in group (%)",
                     trans=scales::pseudo_log_trans(base = 10, sigma = 1),
                     breaks = c(40,20,10,5,2.5,1,0.5,0))+
  histogram.theme

Fig4E.histogramPseudoLog.sigma0.2=
  ggplot()+
  geom_histogram(data = graphData, bins = 100, position = "dodge", aes(x = CHILD_prevalence,y = after_stat(density)*after_stat(width)*100, fill = group))+
  geom_density(data = graphData,bounds = c(0, 1), trim = TRUE, adjust = 5, aes(x = CHILD_prevalence, color = group))+
  scale_color_discrete(guide = "none")+
  scale_fill_discrete(name = "Species presence")+
  scale_x_continuous(name = "Prevalence in CHILD study samples")+
  scale_y_continuous(name = "Species in group (%)",
                     trans=scales::pseudo_log_trans(base = 10, sigma = 0.2),
                     breaks = c(40,20,10,5,2.5,1,0.5, 0.2,0))+
  histogram.theme

Fig4E.histogramPseudoLog.sigma0.1=
  ggplot()+
  geom_histogram(data = graphData, bins = 100, position = "dodge", aes(x = CHILD_prevalence,y = after_stat(density)*after_stat(width)*100, fill = group))+
  geom_density(data = graphData,bounds = c(0, 1), trim = TRUE, adjust = 5, aes(x = CHILD_prevalence, color = group))+
  scale_color_discrete(guide = "none")+
  scale_fill_discrete(name = "Species presence")+
  scale_x_continuous(name = "Prevalence in CHILD study samples")+
  scale_y_continuous(name = "Species in group (%)",
                     trans=scales::pseudo_log_trans(base = 10, sigma = 0.1),
                     breaks = c(40,20,10,5,2.5,1,0.5, 0.2,0.1,0))+
  histogram.theme

# Fig4E.histogramAxisBreak=
#   ggplot()+
#   geom_histogram(data = graphData, bins = 100, position = "dodge", aes(x = CHILD_prevalence,y = after_stat(density)*after_stat(width)*100, fill = group))+
#   geom_density(data = graphData,bounds = c(0, 1), trim = TRUE, adjust = 3, aes(x = CHILD_prevalence, color = group))+
#   scale_color_discrete(guide = "none")+
#   scale_fill_discrete(name = "Species presence")+
#   scale_x_continuous()+
#   scale_y_continuous(name = "Species in group (%)")+
#   xlab("Prevalence in CHILD study samples")+
#   scale_y_break(breaks = c(4,4),  
#                 scales = 0.25, 
#                 # ticklabels = c(46,47),
#                 expand = FALSE)+
#   histogram.theme+
#   theme(legend.title = element_blank(),
#       legend.position = c(1, 1),
#       legend.background = element_blank(),
#       legend.key = element_blank())
# 
# test = get_legend(Fig4E.histogramAxisBreak)
# 
# Fig4E.histogramAxisBreak = Fig4E.histogramAxisBreak+
#   theme(legend.position = "none")


# ==== Figure 5 ====
# ---- A: scatter plot, NS1 DCV1 vs DCV2 vs DCV3, all days ----
graphData = ASV_TableRelAb[,c(1:11,grep("NS1\\.DCV.d..$",colnames(ASV_TableRelAb), fixed = FALSE))]

graphDataMelt = melt(graphData,measure.vars = 12:ncol(graphData), variable.name = "sample")

graphDataMelt = left_join(graphDataMelt,sample_Table[,c(1,5,6)], by = c("sample" = "sample_16S")) 
graphDataMelt$vessel = paste0("Vessel ",graphDataMelt$vessel)

graphDataMeltCast = dcast(graphDataMelt[,c(1,13:15)], formula = ... ~ vessel)

graphData = alply(3:(ncol(graphDataMeltCast)-1), 1, function(x){
  df = melt(graphDataMeltCast,id.vars = 1:x, measure.vars = (x+1):ncol(graphDataMeltCast))
  df = data.frame(df[,1:2],
                  var1 = colnames(df)[x],
                  value1 = df[,x],
                  var2 = df$variable,
                  value2 = df$value)
})
graphData = do.call(rbind,graphData)

graphData1 = graphData[graphData$value1!=0 & graphData$value2!=0,]
graphData1 = graphData1[!is.na(graphData1$ASV_Num),]
graphData1$group = interaction(graphData1$var1,graphData1$var2)

slopes = ddply(graphData1, .(group), function(x){
  model = lm(log10(value2)~log10(value1), data = x)
  slope = summary(model)$coefficients[2,1]
})

rsq = ddply(graphData1, .(group), function(x){
  model = lm(log10(value2)~log10(value1), data = x)
  rsq = summary(model)$r.squared
})

levels(graphData1$group) =  c(paste0("slope ==",signif(slopes$V1[1],3),"*'0;'~R^2 ==", signif(rsq$V1[1],3)),
                              paste0(""),
                              paste0("slope ==",signif(slopes$V1[2],3),"*';'~R^2 ==", signif(rsq$V1[2],3)),
                              paste0("slope ==",signif(slopes$V1[3],3),"*';'~R^2 ==", signif(rsq$V1[3],3))
)

Fig5A.base = ggplot()+
  geom_point(data = graphData, shape = 21, aes(x = value1*100, y = value2*100, fill = factor(day_16S)))+
  geom_smooth(data = graphData1, method = "lm", se = FALSE, aes(x = value1*100, y = value2*100, group = group, color = group))+
  scale_x_log10(name = "Relative Abundance (%)", labels = label_scientific(digits = 1))+
  scale_y_log10(name = "Relative Abundance (%)", labels = label_scientific(digits = 1))+
  scale_color_discrete(name = "Linear Regression",
                       labels = label_parse(), guide = guide_legend(order = 1))+
  scale_fill_discrete(name = "Sample day", guide = guide_legend(order = 2, ncol=2))+
  facet_grid(var2~var1, switch = "both")+
  coord_fixed()+
  point.theme+
  theme(legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_text(size = 26),
        legend.text = element_text(size = 23))

# Turn ggplots into grobs
Fig5A_gTable = ggplotGrob(Fig5A.base)

Fig5A_gTable$grobs[[4]] = nullGrob()
Fig5A_gTable = gtable_add_grob(Fig5A_gTable, grobs = Fig5A_gTable$grobs[[22]], t = 7, l = 8)
Fig5A_gTable$grobs[[22]] = nullGrob()

Fig5A_gTable = Fig5A_gTable[,-12]

Fig5A.plot.gTable.point = Fig5A_gTable

grid.newpage()
grid.draw(Fig5A.plot.gTable.point)





# ---- B: PCoA plot, NS1 DCd00, DCV1-3d01-d17, colour by day, shape by source ----
#Bray-Curtis Dissimilarity
graphData = ASV_TableRelAb

BC_dis = vegdist(t(graphData[,c("NS1.DCd00",
                                "NS1.DCV1d01",
                                "NS1.DCV1d05",
                                "NS1.DCV1d09R",
                                "NS1.DCV1d13",
                                "NS1.DCV1d17",
                                "NS1.DCV2d01",
                                "NS1.DCV2d05",
                                "NS1.DCV2d09R",
                                "NS1.DCV2d13",
                                "NS1.DCV2d17",
                                "NS1.DCV3d01",
                                "NS1.DCV3d05",
                                "NS1.DCV3d09R",
                                "NS1.DCV3d13",
                                "NS1.DCV3d17")]))

#PCoA
PCoA = cmdscale(BC_dis, k = 2, eig = TRUE)
Plot_PCoA = sample_Table %>%
  right_join(tibble::rownames_to_column(data.frame(PCoA$points), var = "sample_16S") %>%
               rename(PCo1 = X1,
                      PCo2 = X2)) 
PCoA.percentvariance = PCoA$eig/sum(PCoA$eig) *100

Fig5B.plot.PCoA = ggplot()+
  geom_point(data = Plot_PCoA, size = 3,aes(x =PCo1, y =PCo2, shape = factor(vessel), color = factor(day_16S)))+
  scale_x_continuous(name = paste0("PCo1 (", signif(PCoA.percentvariance[1],3),"%)"))+
  scale_y_continuous(name = paste0("PCo2 (", signif(PCoA.percentvariance[2],3),"%)"))+
  scale_shape_manual(name = "source", values = c(15:18))+
  scale_color_discrete(name = "day")+
  coord_fixed()+
  PCoA.theme


# ---- C: scatter plot, NS1 FCd17, mean(DCV1-3d17) ----
graphData = ASV_TableRelAb
# graphData = graphData[,c(1:11,grep("NS1.FCd21",colnames(graphData), fixed = TRUE), grep("NS1\\.DCV.d21",colnames(graphData), fixed = FALSE))]
graphData$NS1.DCd17 = rowMeans(graphData[,c("NS1.DCV1d17","NS1.DCV2d17","NS1.DCV3d17")])
# graphData = graphData[graphData$NS1.F!=0 |graphData$NS1.FCd17!=0 | graphData$NS1.DCd17!=0,]
graphData = graphData[graphData$NS1.DCd00!=0,]
# graphData = graphData[rowSums(graphData[,c("NS1.F", "NS1.FCd01", "NS1.FCd05", "NS1.FCd09", "NS1.FCd13", "NS1.FCd17",
#                                            "NS1.DCd00")])!=0,]

graphData = graphData[,c(1:11,which(colnames(graphData) %in% c("NS1.F","NS1.FCd17","NS1.DCd00","NS1.DCd17")))]

graphData$presenceGroupColor = interaction(graphData$NS1.F!=0,
                                           graphData$NS1.FCd17!=0
)

# graphData$presenceGroupColor = interaction(graphData$NS1.F!=0,
#                                            rowSums(graphData[,c("NS1.FCd01", "NS1.FCd05", "NS1.FCd09", "NS1.FCd13")])!=0,
#                                            graphData$NS1.FCd17!=0
# )

graphData1 = graphData[graphData$NS1.FCd17!=0 & graphData$NS1.DCd17!=0,]

slopes = summary(lm(log10(NS1.DCd17)~log10(NS1.FCd17), data = graphData1))$coefficients[2,1]

rsq = summary(lm(log10(NS1.DCd17)~log10(NS1.FCd17), data = graphData1))$r.squared

graphData1$label = paste0("slope ==",signif(slopes,3),"*';'~R^2 ==", signif(rsq,3))

# generate legend data
# generate count data for absence/presence data
legendData1 = graphData[,c("NS1.F", "NS1.FCd17","NS1.DCd00")]

legendData = aggregate(legendData1, by = list("levelName" = graphData$presenceGroupColor), FUN = function(x) sum(x>0))
legendData = legendData[match(absencePresence2_0LegendData$levelName,legendData$levelName),]
legendData$y = c(3,2,1,0)

legendDataMelt = melt(legendData, id.vars = c(1,ncol(legendData)), variable.name = "x")
legendDataMelt$filled = as.numeric(legendDataMelt$value>0)

legendDataMelt$value = paste0("n = ",legendDataMelt$value)
legendDataMelt$value[legendDataMelt$x != "NS1.DCd00"] = ""

levels(legendDataMelt$x) = sampleNameLevels.3lines[levels(legendDataMelt$x)]

Fig5C.legend = ggplot()+
  geom_tile(data = legendDataMelt[legendDataMelt$filled == 1,],
            color = "black",
            size = 0.5,
            aes(x = x, y = y, fill = levelName))+
  geom_text(data = legendDataMelt[legendDataMelt$filled == 1,],
            size = 13,
            aes(x = x, y = y, label = value))+
  geom_tile(data = legendDataMelt[legendDataMelt$filled == 0,],
            color = "black",
            fill = "white",
            size = 0.5,
            aes(x = x, y = y))+
  # geom_tile(data = legendDataMelt,
  #           color = "black",
  #           height = 0.70,
  #           width = 0.70,
  #           size = 0.25,
  #           aes(x = "", y = y, fill = levelName))+
  scale_x_discrete(position = "top", limits = c(levels(legendDataMelt$x)))+
  scale_y_continuous(expand = expansion(add = 0.1))+
  ggtitle("Identified in:")+
  coord_fixed()+
  legend.chemostat3.theme

# base plot
Fig5C.base = ggplot()+
  geom_point(data = graphData, aes(x = NS1.FCd17*100, y = NS1.DCd17*100, color = presenceGroupColor))+
  scale_colour_discrete(guide = "none")+
  new_scale_color()+
  geom_smooth(data = graphData1, method = "lm", se = FALSE, aes(x = NS1.FCd17*100, y = NS1.DCd17*100, color = label))+
  scale_colour_discrete(name = "", type = "black",
                        labels = label_parse())+
  geom_line(color = "grey90", aes(x = c(seq(0,0.001,0.0001),0.01,0.1,1,10), y = c(seq(0,0.001,0.0001),0.01,0.1,1,10)))+
  scale_x_continuous(name = "fecal chemostat\nRelative Abundance (%)", labels = label_scientific(digits = 1),
                     trans=scales::pseudo_log_trans(base = 10, sigma = 0.001), breaks = c(0,0.01, 0.1, 1, 10), expand = expansion(mult = (13.5/13*1.05)-1))+
  scale_y_continuous(name = "defined community chemostat\nRelative Abundance (%)", labels = label_scientific(digits = 1),
                     trans=scales::pseudo_log_trans(base = 10, sigma = 0.001), breaks = c(0, 0.001, 0.01, 0.1, 1, 10), expand = expansion(mult = (28.5/28*1.05)-1))+
  coord_fixed()+
  point.theme_legendtl

Fig5C_gTable = ggplotGrob(Fig5C.base)
# add row for top histogram
Fig5C_gTable = gtable_add_rows(Fig5C_gTable, unit(0.25, "null"), pos = 6)

# add row for right histogram axis
Fig5C_gTable = gtable_add_rows(Fig5C_gTable, unit(0, "null"), pos = 8)

# add rows for xlab-b and space for right histogram axis title
Fig5C_gTable = gtable_add_rows(Fig5C_gTable, unit(0.5, "grobheight",Fig5C_gTable$grobs[[12]]), pos = 10)
Fig5C_gTable = gtable_add_rows(Fig5C_gTable, unit(0.5, "grobheight",Fig5C_gTable$grobs[[12]]), pos = 10)

# add column for right histogram
Fig5C_gTable = gtable_add_cols(Fig5C_gTable, unit(0.25, "null"), pos = 5)
# add column for top histogram axis
Fig5C_gTable = gtable_add_cols(Fig5C_gTable, unit(0, "null"), pos = 4)

# Fig5C_gTable$heights
# Fig5C_gTable$widths

# histograms on t and r
Fig5C.histt = ggplot()+
  geom_histogram(data = graphData, position = "dodge", bins = 13, center = 0, aes(x = NS1.FCd17*100, fill = presenceGroupColor))+
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10, sigma = 0.001), breaks = c(0, 0.01, 0.1, 1, 10))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10, sigma = 1), breaks = c(0, 10, 25, 50, 100))+
  scale_fill_discrete(guide = "none")+
  theme_bw()+
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.x = element_blank())

Fig5C.histt_gTable = ggplotGrob(Fig5C.histt)

# add top histogram
Fig5C_gTable = gtable_add_grob(Fig5C_gTable, Fig5C.histt_gTable$grobs[[6]], t = 7, l = 6)
# add top histogram axis
Fig5C_gTable = gtable_add_grob(Fig5C_gTable, Fig5C.histt_gTable$grobs[[3]], t = 7, l = 5)
# add top histogram axis title
Fig5C_gTable = gtable_add_grob(Fig5C_gTable, Fig5C.histt_gTable$grobs[[13]], t = 7, l = 4)

# plot(Fig5C_gTable)

Fig5C.histr = ggplot()+
  geom_histogram(data = graphData, position = "dodge", bins = 28, center = 0, aes(x = NS1.DCd17*100, fill = presenceGroupColor))+
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10, sigma = 0.001), breaks = c(0, 0.001, 0.01, 0.1, 1, 10))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10, sigma = 1), breaks = c(0, 10, 25, 50, 100))+
  scale_fill_discrete(guide = "none")+
  coord_flip()+
  theme_bw()+
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank())


Fig5C.histr_gTable = ggplotGrob(Fig5C.histr)

# add right histogram
Fig5C_gTable = gtable_add_grob(Fig5C_gTable, Fig5C.histr_gTable$grobs[[6]], t = 8, l = 7)
# add right histogram axis
Fig5C_gTable = gtable_add_grob(Fig5C_gTable, Fig5C.histr_gTable$grobs[[7]], t = 9, l = 7)

# move xlab-b
Fig5C_gTable = gtable_add_grob(Fig5C_gTable, Fig5C_gTable$grobs[[12]], t = 11, b = 12, l = 6)
Fig5C_gTable$grobs[[12]] = nullGrob()

# add right histogram axis title
Fig5C_gTable = gtable_add_grob(Fig5C_gTable, Fig5C.histr_gTable$grobs[[12]], t = 10, b = 11, l = 7)

Fig5C.plot.gTable.point = Fig5C_gTable

grid.newpage()
grid.draw(Fig5C.plot.gTable.point)


# ==== Figure 6 ====
# ---- A: NMR metabolite line plot, both FCd03-d16, facet by source, colour by metabolite ----
graphData = NMRmetaboliteMelt

graphData = graphData[graphData$type == "FC" & !is.na(graphData$type),]
# only using data up to d17, after day 1
graphData = graphData[graphData$day_NMR<=17,]
graphData = graphData[graphData$day_NMR>1,]

# capture average order of metabolites on last day in display
metaboliteOrder = ddply(graphData[graphData$day_NMR == 16,], .(metabolite), function(x){
  mean(x$conc)
})
colnames(metaboliteOrder)[2] = "mean"
metaboliteOrder = metaboliteOrder[order(metaboliteOrder$mean,decreasing = TRUE),]
# reorder levels to match mean metabolite concentrations
graphData$metabolite = factor(graphData$metabolite, levels = metaboliteOrder$metabolite)

# reorder colors to improve contrasts
metaboliteReorder = as.vector(cbind(seq.int(1,21,3),seq.int(2,21,3),seq.int(3,21,3)))
metaboliteReorder = metaboliteReorder[-21]

Fig6A.plot.line = ggplot()+
  geom_line(data = graphData, size = 1.5, aes(x = day_NMR, y = conc, group = interaction(source,type,vessel,metabolite), color = metabolite))+
  scale_x_continuous(name = "Day")+
  scale_y_log10(name = "[Metabolite] (mM)")+
  scale_color_manual(name = "Metabolite", values = (hue_pal()(20))[metaboliteReorder])+
  facet_wrap(~source)+
  line.theme


# ---- B: MS metabolite heatmap, NS0 FCd02-d16, facet by source, colour by metabolite ----
X_metabolites = MSmetabolite_TableFC_clean
# X_metabolites = MSmetabolite_TableFC_clean_Zscore

X_metabolites = X_metabolites[X_metabolites$source == "NS0" & !is.na(X_metabolites$type),]
X_metabolites = X_metabolites[X_metabolites$day_MS<=17,]
X_metabolites = X_metabolites[X_metabolites$day_MS>1,]

graphData = X_metabolites[,8:ncol(X_metabolites)]
graphData = do.call(rbind,apply(graphData,1,function(x) x/graphData[1,]))
graphData = log10(graphData)
graphData = t(graphData)
colnames(graphData) = X_metabolites$sample_MSMetabolites
graphData = graphData[apply(graphData,1,function(x) !any(is.nan(x))),]
# graphData = graphData[apply(graphData,1,function(x) sd(x)!=0),]
graphData = graphData[apply(graphData,1,function(x) !any(is.infinite(x))),]
colnames(graphData) = X_metabolites$day_MS

graphDataHclust = hclust(dist(1-cor(t(graphData))), method = "complete")

colBreak = max(c(abs(min(graphData)),abs(max(graphData))))

Fig6B.heatplot = Heatmap(graphData,
                         name = "log10(fold change) from day 2", #title of legend
                         row_title = "Metabolites", column_title = "Timepoints (day)",
                         show_column_names = TRUE,
                         show_row_names = TRUE,
                         row_names_gp = gpar(fontsize = 10), # Text size for row names
                         column_title_side = c("bottom"),
                         column_names_gp = gpar(fontsize = 15), # Text size for column names
                         column_names_rot = 0, # column names rotation
                         clustering_method_rows = "complete",
                         cluster_columns = FALSE,
                         row_km = 8,
                         show_parent_dend_line = FALSE,
                         column_gap = unit(0.005, "npc"),
                         col = circlize::colorRamp2(breaks = c(-colBreak,0,colBreak),
                                                    colors = c("#4575B4", "#FFFFBF", "#D73027")),
                         heatmap_legend_param = list(
                           legend_width = unit(12, "cm"),
                           direction = "horizontal"
                         )
)

Fig6B.heatplot = ComplexHeatmap::draw(object = Fig6B.heatplot, heatmap_legend_side = "top")

# ---- C: MS metabolite heatmap, NS1 FCd02-d16, facet by source, colour by metabolite ----
X_metabolites = MSmetabolite_TableFC_clean
# X_metabolites = MSmetabolite_TableFC_clean_Zscore

X_metabolites = X_metabolites[X_metabolites$source == "NS1" & !is.na(X_metabolites$type),]
X_metabolites = X_metabolites[X_metabolites$day_MS<=17,]
X_metabolites = X_metabolites[X_metabolites$day_MS>1,]

graphData = X_metabolites[,8:ncol(X_metabolites)]
graphData = do.call(rbind,apply(graphData,1,function(x) x/graphData[1,]))
graphData = log10(graphData)
graphData = t(graphData)
colnames(graphData) = X_metabolites$sample_MSMetabolites
graphData = graphData[apply(graphData,1,function(x) !any(is.nan(x))),]
# graphData = graphData[apply(graphData,1,function(x) sd(x)!=0),]
graphData = graphData[apply(graphData,1,function(x) !any(is.infinite(x))),]
colnames(graphData) = X_metabolites$day_MS

graphDataHclust = hclust(dist(1-cor(t(graphData))), method = "complete")

colBreak = max(c(abs(min(graphData)),abs(max(graphData))))

Fig6C.heatplot = Heatmap(graphData,
                         name = "log10(fold change) from day 2", #title of legend
                         row_title = "Metabolites", column_title = "Timepoints (day)",
                         show_column_names = TRUE,
                         show_row_names = TRUE,
                         row_names_gp = gpar(fontsize = 10), # Text size for row names
                         column_title_side = c("bottom"),                          
                         column_names_gp = gpar(fontsize = 15), # Text size for column names
                         column_names_rot = 0, # column names rotation
                         clustering_method_rows = "complete",
                         cluster_columns = FALSE,
                         row_km = 8,
                         show_parent_dend_line = FALSE,
                         column_gap = unit(0.005, "npc"),
                         col = circlize::colorRamp2(breaks = c(-colBreak,0,colBreak),
                                                    colors = c("#4575B4", "#FFFFBF", "#D73027")),
                         heatmap_legend_param = list(
                           legend_width = unit(12, "cm"),
                           direction = "horizontal"
                         )
)

Fig6C.heatplot = ComplexHeatmap::draw(object = Fig6C.heatplot, heatmap_legend_side = "top")


# ---- D: NMR metabolite scatter plot, NS1 DCV1 vs DCV2 vs DCV3 all days ----
X_metabolites = NMRmetabolite_Table

X_metabolites = X_metabolites[X_metabolites$source == "NS1" & !is.na(X_metabolites$type) & X_metabolites$type == "DC",]

X_metabolitesMelt = reshape2::melt(X_metabolites,measure.vars = 7:ncol(X_metabolites), variable.name = "metabolite")
X_metabolitesMeltCast = dcast(X_metabolitesMelt[,3:8], formula = ... ~ vessel)

colnames(X_metabolitesMeltCast)[5:7] = c("Vessel 1",
                                         "Vessel 2",
                                         "Vessel 3")

graphData = alply(5:(ncol(X_metabolitesMeltCast)-1), 1, function(x){
  df = melt(X_metabolitesMeltCast,id.vars = 1:x, measure.vars = (x+1):ncol(X_metabolitesMeltCast))
  df = data.frame(df[,1:4],
                  var1 = colnames(df)[x],
                  value1 = df[,x],
                  var2 = df$variable,
                  value2 = df$value)
})
graphData = do.call(rbind,graphData)

graphData1 = graphData[graphData$value1!=0 & graphData$value2!=0,]
graphData1$group = interaction(graphData1$var1,graphData1$var2)

slopes = ddply(graphData1, .(group), function(x){
  model = lm(log10(value2)~log10(value1), data = x)
  slope = summary(model)$coefficients[2,1]
})

rsq = ddply(graphData1, .(group), function(x){
  model = lm(log10(value2)~log10(value1), data = x)
  rsq = summary(model)$r.squared
})

levels(graphData1$group) =  c(paste0("slope ==",signif(slopes$V1[1],3),"*'0;'~R^2 ==", signif(rsq$V1[1],3)),
                              paste0(""),
                              paste0("slope ==",signif(slopes$V1[2],3),"*';'~R^2 ==", signif(rsq$V1[2],3),"*'0'"),
                              paste0("slope ==",signif(slopes$V1[3],3),"*';'~R^2 ==", signif(rsq$V1[3],3))
)

Fig6D.base = ggplot()+
  geom_point(data = graphData, shape = 21, aes(x = value1, y = value2, fill = factor(day_NMR)))+
  geom_smooth(data = graphData1, method = "lm", se = FALSE, aes(x = value1, y = value2, group = group, color = group))+
  scale_x_log10(name = "[metabolite] (mM)", labels = label_scientific(digits = 1))+
  scale_y_log10(name = "[metabolite] (mM)", labels = label_scientific(digits = 1))+
  scale_color_discrete(name = "Linear Regression",
                       labels = label_parse())+
  scale_fill_discrete(name = "Sample day", guide = guide_legend(ncol=2))+
  facet_grid(var2~var1, switch = "both")+
  coord_fixed()+
  point.theme+
  theme(legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_text(size = 26),
        legend.text = element_text(size = 23)
  )

# Turn ggplots into grobs
Fig6D_gTable = ggplotGrob(Fig6D.base)

is.gtable(Fig6D_gTable$grobs[[22]])

Fig6D_gTable$grobs[[4]] = nullGrob()

# Fig6a_gTable$grobs[[4]] = Fig6a_gTable$grobs[[22]][2:6,2:4]

Fig6D_gTable = gtable_add_grob(Fig6D_gTable, grobs = Fig6D_gTable$grobs[[22]], 
                               t = 7, l = 8,
                               b = 7, r = 8,
                               name = "new_guide_box")

Fig6D_gTable$grobs[[22]] = nullGrob()

Fig6D_gTable = Fig6D_gTable[,-12]

Fig6D.plot.gTable.point = Fig6D_gTable

grid.newpage()
grid.draw(Fig6D.plot.gTable.point)

# ---- E: MS metabolite scatter plot, NS1 DCV1 vs DCV2 vs DCV3 all days ----
X_metabolites = MSmetabolite_TableDC_clean
# X_metabolites = MSmetabolite_TableDC_clean_Zscore

X_metabolites = X_metabolites[X_metabolites$source == "NS1" & !is.na(X_metabolites$type),]
X_metabolites = X_metabolites[X_metabolites$day_MS<=17,]
# X_metabolites = X_metabolites[X_metabolites$day_MS>1,]

X_metabolitesMelt = reshape2::melt(X_metabolites,measure.vars = 8:ncol(X_metabolites), variable.name = "metabolite")
X_metabolitesMeltCast = dcast(X_metabolitesMelt[,4:9], formula = ... ~ vessel)

colnames(X_metabolitesMeltCast)[5:7] = c("Vessel 1",
                                         "Vessel 2",
                                         "Vessel 3")

graphData = alply(5:(ncol(X_metabolitesMeltCast)-1), 1, function(x){
  df = melt(X_metabolitesMeltCast,id.vars = 1:x, measure.vars = (x+1):ncol(X_metabolitesMeltCast))
  df = data.frame(df[,1:4],
                  var1 = colnames(df)[x],
                  value1 = df[,x],
                  var2 = df$variable,
                  value2 = df$value)
})
graphData = do.call(rbind,graphData)

graphData1 = graphData[graphData$value1!=0 & graphData$value2!=0,]
graphData1$group = interaction(graphData1$var1,graphData1$var2)

slopes = ddply(graphData1, .(group), function(x){
  model = lm(log10(value2)~log10(value1), data = x)
  slope = summary(model)$coefficients[2,1]
})

rsq = ddply(graphData1, .(group), function(x){
  model = lm(log10(value2)~log10(value1), data = x)
  rsq = summary(model)$r.squared
})

levels(graphData1$group) =  c(paste0("slope ==",signif(slopes$V1[1],3),"*';'~R^2 ==", signif(rsq$V1[1],3)),
                              paste0(""),
                              paste0("slope ==",signif(slopes$V1[2],3),"*';'~R^2 ==", signif(rsq$V1[2],3)),
                              paste0("slope ==",signif(slopes$V1[3],3),"*'.004;'~R^2 ==", signif(rsq$V1[3],3),"*'0'")
)

Fig6E.base = ggplot()+
  geom_point(data = graphData, shape = 21, aes(x = value1, y = value2, fill = factor(day_MS)))+
  geom_smooth(data = graphData1, method = "lm", se = FALSE, aes(x = value1, y = value2, group = group, color = group))+
  scale_x_log10(name = "Arbitrary Units", labels = label_scientific(digits = 1))+
  scale_y_log10(name = "Arbitrary Units", labels = label_scientific(digits = 1))+
  scale_color_discrete(name = "Linear Regression",
                       labels = label_parse())+
  scale_fill_discrete(name = "Sample day", guide = guide_legend(ncol=2))+
  facet_grid(var2~var1, switch = "both")+
  coord_fixed()+
  point.theme+
  theme(legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_text(size = 26),
        legend.text = element_text(size = 23)
        )

# Turn ggplots into grobs
Fig6E_gTable = ggplotGrob(Fig6E.base)

is.gtable(Fig6E_gTable$grobs[[22]])

Fig6E_gTable$grobs[[4]] = nullGrob()

# Fig6_gTable$grobs[[4]] = Fig6_gTable$grobs[[22]][2:6,2:4]

Fig6E_gTable = gtable_add_grob(Fig6E_gTable, grobs = Fig6E_gTable$grobs[[22]], 
                              t = 7, l = 8,
                              b = 7, r = 8,
                              name = "new_guide_box")

Fig6E_gTable$grobs[[22]] = nullGrob()

Fig6E_gTable = Fig6E_gTable[,-12]

Fig6E.plot.gTable.point = Fig6E_gTable

grid.newpage()
grid.draw(Fig6E.plot.gTable.point)

#==== save plots ====
#save all .plot.alluvial.ap PDF files to ./Figure Outputs
a_ply(ls()[grepl(".plot.alluvial.ap",ls(),fixed = TRUE)],1,function(x){
  ggsave(paste0(path,"Figure Outputs/",newDir,"/",x,".pdf"),
         plot = eval(parse(text = x)),
         device=cairo_pdf, # add this to embed the Calibri font
         width = 7.126, height = 10,
         units = "in", # other options c("in", "cm", "mm"),
         dpi = 300)
})

#save all .legend PDF files to ./Figure Outputs
a_ply(ls()[grepl(".legend",ls(),fixed = TRUE)],1,function(x){
  ggsave(paste0(path,"Figure Outputs/",newDir,"/",x,".pdf"),
         plot = eval(parse(text = x)),
         device=cairo_pdf, # add this to embed the Calibri font
         width = 2.874, height = 10,
         units = "in", # other options c("in", "cm", "mm"),
         dpi = 300)
})

#save all .plot.alluvial.ph PDF files to ./Figure Outputs
a_ply(ls()[grepl(".plot.alluvial.ph",ls(),fixed = TRUE)],1,function(x){
  ggsave(paste0(path,"Figure Outputs/",newDir,"/",x,".pdf"),
         plot = eval(parse(text = x)),
         device=cairo_pdf, # add this to embed the Calibri font
         width = 7.126, height = 10,
         units = "in", # other options c("in", "cm", "mm"),
         dpi = 300)
})

#save all .plot.PCoA PDF files to ./Figure Outputs
a_ply(ls()[grepl(".plot.PCoA",ls(),fixed = TRUE)],1,function(x){
  ggsave(paste0(path,"Figure Outputs/",newDir,"/",x,".pdf"),
         plot = eval(parse(text = x)),
         device=cairo_pdf, # add this to embed the Calibri font
         width = 10, height = 10,
         units = "in", # other options c("in", "cm", "mm"),
         dpi = 300)
})

#save all .plot.line PDF files to ./Figure Outputs
a_ply(ls()[grepl(".plot.line",ls(),fixed = TRUE)],1,function(x){
  ggsave(paste0(path,"Figure Outputs/",newDir,"/",x,".pdf"),
         plot = eval(parse(text = x)),
         device=cairo_pdf, # add this to embed the Calibri font
         width = 10, height = 10,
         units = "in", # other options c("in", "cm", "mm"),
         dpi = 300)
})

#save all .heatplot PDF files to ./Figure Outputs
a_ply(ls()[grepl(".heatplot",ls(),fixed = TRUE)],1,function(x){
  pdf(paste0(path,"Figure Outputs/",newDir,"/",x,".pdf"), width = 10, height = 10)
  ComplexHeatmap::draw(object = get(x))
  dev.off()
})

#save all .plot.point PDF files to ./Figure Outputs
a_ply(ls()[grepl(".plot.point",ls(),fixed = TRUE)],1,function(x){
  ggsave(paste0(path,"Figure Outputs/",newDir,"/",x,".pdf"),
         plot = eval(parse(text = x)),
         device=cairo_pdf, # add this to embed the Calibri font
         width = 10, height = 10,
         units = "in", # other options c("in", "cm", "mm"),
         dpi = 300)
})

#save all .plot.gTable.point PDF files to ./Figure Outputs
a_ply(ls()[grepl(".plot.gTable.point",ls(),fixed = TRUE)],1,function(x){
  grid.newpage()
  cairo_pdf(filename = paste0(path,"Figure Outputs/",newDir,"/",x,".pdf"),
            width = 10, height = 10,
            bg = "transparent")
  grid.draw(eval(parse(text = x)))
  dev.off()
})

#save all .legend.gTable PDF files to ./Figure Outputs
a_ply(ls()[grepl(".legend.gTable",ls(),fixed = TRUE)],1,function(x){
  grid.newpage()
  cairo_pdf(filename = paste0(path,"Figure Outputs/",newDir,"/",x,".pdf"),
            width = 10, height = 10,
            bg = "transparent")
  grid.draw(eval(parse(text = x)))
  dev.off()
})

#save all .plot_grid.Norm PDF files to ./Figure Outputs
a_ply(ls()[grepl(".plot_grid.Norm",ls(),fixed = TRUE)],1,function(x){
  grid.newpage()
  cairo_pdf(filename = paste0(path,"Figure Outputs/",newDir,"/",x,".pdf"),
            width = 20, height = 8,
            bg = "transparent")
  grid.draw(eval(parse(text = x)))
  dev.off()
})

a_ply(ls()[grepl(".plot_grid.CoordFlip",ls(),fixed = TRUE)],1,function(x){
  grid.newpage()
  cairo_pdf(filename = paste0(path,"Figure Outputs/",newDir,"/",x,".pdf"),
            width = 8, height = 20,
            bg = "transparent")
  grid.draw(eval(parse(text = x)))
  dev.off()
})

#save all .histogram PDF files to ./Figure Outputs
a_ply(ls()[grepl(".histogram",ls(),fixed = TRUE)],1,function(x){
  grid.newpage()
  cairo_pdf(filename = paste0(path,"Figure Outputs/",newDir,"/",x,".pdf"),
            width = 10, height = 10,
            bg = "transparent")
  grid.draw(eval(parse(text = x)))
  dev.off()
})

# write.table(CHILD_mergedASV_Table_SpeciesCol[,1:12],file = paste0(path,"Figure Outputs/",newDir,"/","NS0_NS1_isolates speciesCol.txt"), sep = "\t")


#save all graphDataList files to ./Figure Outputs
# a_ply(names(graphDataList),1,function(x){
#   grid.newpage()
#   cairo_pdf(filename = paste0(path,"Figure Outputs/",newDir,"/",x,".pdf"),
#             width = 10, height = 10,
#             bg = "transparent")
#   grid.draw(graphDataList[[x]][[1]])
#   dev.off()
# })
# 

