#!/bin/R
###################################################################################################
## Set log
###################################################################################################
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

###################################################################################################
## Load required packages
###################################################################################################
libs <- paste0(Sys.getenv("CONDA_PREFIX"),"/lib/R/library")
.libPaths(libs)
#.libPaths()
print("Loading required R libraries")
if(!require(genomeIntervals)){
  BiocManager::install("GenomeInfoDbData",update=F,ask=F)
  require(genomeIntervals)
}
library(ggplot2)
library(gtools)
library(gtable)
library(data.table)
library(reshape)
library(grid)
library(grDevices)
library(stringr)
library(xtable)
#binny libraries:
library(caTools)
library(fpc)
library(FNN)
library(RColorBrewer)
library(scales)
library(diptest)
library(mixtools)
library(gclus)

###################################################################################################
## Input file locations from snakemake
###################################################################################################

cluster_dir         <- snakemake@input[["outdir"]]
binny_dir           <- snakemake@params[["binnydir"]]
MG.depth_file       <- snakemake@input[["mgdepth"]]
coords_file         <- snakemake@input[["vizbin"]]
annot_file          <- snakemake@input[["gff"]]
function_script     <- snakemake@params[["plot_functions"]]
pk 		    <- snakemake@config[["binning"]][["binny"]][["pk"]]
nn 		    <- snakemake@config[["binning"]][["binny"]][["nn"]]
all_outputs <- snakemake@output

###################################################################################################
## Initialize functions for various calculations and normalizations
###################################################################################################
print("START: Reading functions")
source(function_script)
print("DONE: Reading functions")

###################################################################################################
## Read in the necessary input files
###################################################################################################
print("START: Reading data")

## MG depth (bedtools coverageBed output)
print("Read in MG depth file")
MG.depth <- read.table(MG.depth_file, colClasses=c("factor", "numeric"))
colnames(MG.depth) <- c("contig", "MG_depth")

## gff annotation file
print("Read in gff3 annotation file")
annot <- readGff3(annot_file, isRightOpen = FALSE, quiet=TRUE)

# extract only interesting information, i.e start and end of genes
print("Processing gff3 annotation file")
dbs <- setdiff(unique(unlist(lapply(parseGffAttributes(annot),names))),c("ID","inference","locus_tag","product","note","Name","gene"))
annot.1 <- as.data.frame(
  cbind(
    as.character(annot@annotation$seq_name),
    annot@.Data,
    getGffAttribute(annot, dbs)
  )
)
colnames(annot.1)[1:3] <- c("contig", "start", "end")

# Remove rows with no start and stop positions
if (nrow(annot.1[is.na(annot.1$start),]) != 0) {
  annot.1 <- annot.1[-is.na(annot.1$start),]
} else if (nrow(annot.1[is.na(annot.1$end),]) != 0) {
  annot.1 <- annot.1[-is.na(annot.1$end),]
} else {
  print("Genes are OK!")
}

# Replace non-annotated genes with NA
annot.1[,2:3] <- sapply(annot.1[2:3], as.character)
annot.1[,2:3] <- sapply(annot.1[2:3], as.numeric)

# save.image(paste0(binny_dir,"/binny_WS.Rdata"))

# create annotation table with essential genes
if(!"essential" %in% colnames(annot.1)){
  # if no essential genes are found, binny is ended
  for(o in setdiff(unlist(all_outputs),paste0(binny_dir,"/binny_WS.Rdata"))){
    system(paste0("touch ",o))
  }
  print("no essential genes. Ending Binny.")
  system("touch binny.error")
}else{
  print("Creating annotation table")
  annot.E.prelim <- annot.1[,c("contig","start","end","essential")]
  annot.E.prelim <- annot.E.prelim[!is.na(annot.E.prelim$essential),]
  annot.E.prelim$essential <- as.character(annot.E.prelim$essential)
  annot.E.prelim$contig <- as.character(annot.E.prelim$contig)
  # save.image(paste0(binny_dir,"/binny_WS.Rdata"))
  
# read vizbin coordinates
  
  print("Reading in vizbin coordinates")
  coords <- read.table(coords_file, colClasses=c("factor", "numeric", "numeric"),
                       sep="\t", col.names=c("contig", "x", "y"))
  
  if(nrow(coords)==0 | length(intersect(annot.E.prelim$contig,coords$contig))==0 ){
    # if no coordinates exist, or non of the essential genes are on contigs with coordinates, binny ends.
    for(o in setdiff(unlist(all_outputs),paste0(binny_dir,"/binny_WS.Rdata"))){
      system(paste0("touch ",o))
    }
    system("touch binny.error")
    print("no vizbin coordinates for essential genes. Ending Binny")
  }else{
    print("DONE: Reading data")
    # save.image(paste0(binny_dir,"/binny_WS.Rdata"))
    
    ###################################################################################################
    ## Merge coordinates, essential gene annotations, and metagenomic depth
    ###################################################################################################
    print("START: Incorporating VizBin, depth, and annotation data")
    
    # depth and coordinates per contig
    vb_dat <- merge(MG.depth, coords, by=c("contig"), all=T, incomparables=NA)
    vb_dat <- vb_dat[!is.na(vb_dat$x),]
    vb_dat$MG_depth[is.na(vb_dat$MG_depth)] <- 0
    # save.image(paste0(binny_dir,"/binny_WS.Rdata"))
    
    # per-contig data and essential genes
    contigInfo <- vb_dat
    contigInfo$contig <- as.character(contigInfo$contig)
    annot.E <- annot.E.prelim[annot.E.prelim$contig %in% vb_dat$contig,]
    essAnno <- aggregate(annot.E$essential,list(annot.E$contig),function(x) paste(x,sep=";",collapse=";"))
    colnames(essAnno) <- c("contig","essentialGene")
    # save.image(paste0(binny_dir,"/binny_WS.Rdata"))
    
    contigInfo <- merge(contigInfo,essAnno,by.x=1,by.y=1,all.x=T)
    contigInfo$essentialGene[is.na(contigInfo$essentialGene)] <- "notEssential"
    # some genes have more than one model, unify these here.
    hmm2 <- data.frame("gene"=rep(c("glyS","proS","pheT","rpoC"),each=2),
                       "HMM"=paste("TIGR",c("00388","00389","00408","00409","00471","00472","02386","02387"),sep=""),
                        stringsAsFactors=F)
    for(hmm in hmm2$HMM) contigInfo$essentialGene <- gsub(hmm,hmm2$gene[hmm2$HMM==hmm],contigInfo$essentialGene)    

    print("DONE: Incorporating VizBin and annotation data")
    # save.image(paste0(binny_dir,"/binny_WS.Rdata"))
    
    ###################################################################################################
    ## binny - 1st iteration
    ###################################################################################################
    
    print("START: Running binny on VizBin data")
    
    #### find which columns contain coordinates
    coco <- c(which(colnames(contigInfo)=="x"),which(colnames(contigInfo)=="y"))
    ### automated estimate of reachability distance (based on number of neighbouring points)
    skn <- sort(knn.dist(contigInfo[,coco],pk)[,pk]) #sort distance of neighbouring points
    sdkn <- runsd(skn,10) #calculate running standard deviation between 10 neighbouring points
    est <- sort(skn)[min(which(sdkn>quantile(sdkn,0.975)&skn>median(skn)))] #find the first jump in distances at the higher end
    #### document estimated reachability distances
    write.table(t(c("scan","reachabilityDist")),
                paste0(binny_dir,paste("/reachabilityDistanceEstimates",pk,nn,"tsv",sep=".")),
                row.names=F,col.names=F,quote=F,sep="\t")
    write.table(t(c("first",est)),
                paste0(binny_dir,paste("/reachabilityDistanceEstimates",pk,nn,"tsv",sep=".")),
                row.names=F,col.names=F,quote=F,sep="\t",append=T)
    ### run DBscan with reachibility distance
    cdb <- dbscan(contigInfo[,coco],est,pk)

    # save.image(paste0(binny_dir,"/binny_WS.Rdata"))
    ### count contigs and essential genes in clusters and plot coordinates with cluster membership
    pdf(paste0(binny_dir,paste("/scatterPlot1",pk,nn,"pdf",sep=".")))
    plot(contigInfo[,coco],pch=16,cex=0.25,ann=F,axes=F)
    j <- 1
    cdbTab <- data.frame("cluster"=names(table(cdb$cluster)),"contigs"=0,"numberEss"=0,"uniqueEss"=0,stringsAsFactors=F)
    for(i in names(table(cdb$cluster))) {
      points(contigInfo[cdb$cluster==i,coco],pch=16,cex=0.4,col=c("grey",colors(distinct=T))[j])
      cdbTab$contigs[cdbTab$cluster==i] <- table(cdb$cluster)[names(table(cdb$cluster))==i]
      cdbTab$numberEss[cdbTab$cluster==i] <- length(unlist(sapply(contigInfo$essentialGene[contigInfo$essentialGene!="notEssential"&
                                                                                             cdb$cluster==i],
                                                                  function(x) unlist(strsplit(x,split=";")))))
      cdbTab$uniqueEss[cdbTab$cluster==i] <- length(unique(unlist(sapply(contigInfo$essentialGene[contigInfo$essentialGene!="notEssential"&
                                                                                                    cdb$cluster==i],
                                                                         function(x) unlist(strsplit(x,split=";"))))))
      j<-j+1
    }
    box()
    dev.off()
    
    ### initialize tables for cluster membership
    write.table(cdbTab,paste0(binny_dir,paste("/clusterFirstScan",pk,nn,"tsv", sep=".")),
                sep="\t",row.names=F,quote=F)
    write.table(t(c("clusterName","cutoff")),paste0(binny_dir,paste("/bimodalClusterCutoffs",pk,nn,"tsv",sep=".")),
                sep="\t",row.names=F,col.names=F,quote=F)
    
    ### assign bin names to contigs 
    #### names are "N" for contigs recognized as noise by DBscan
    #### all other bins keep their DBscan number and get a letter at the front
    #### "E" for bins without essential genes
    #### "C" for bins with less than or equal to 20% duplicated essential genes
    #### "D" for bins with more than 20% duplicated essential genes
    clusterRes <- data.frame("contig"=contigInfo$contig,
                             "cluster"="x",
                             stringsAsFactors=F)
    #### N:
    clusterRes$cluster[cdb$cluster==0] <- "N" #DBscan calls noise clusters "0"
    #### E:
    emptyClus <- cdbTab$cluster[cdbTab$uniqueEss==0&cdbTab$cluster!=0]
    clusterRes$cluster[cdb$cluster %in% emptyClus] <- paste("E",cdb$cluster[cdb$cluster %in% emptyClus],sep="")
    #### C:
    closeClus <- cdbTab$cluster[cdbTab$numberEss/cdbTab$uniqueEss<=1.2&cdbTab$cluster!=0]
    clusterRes$cluster[cdb$cluster %in% closeClus] <- paste("C",cdb$cluster[cdb$cluster %in% closeClus],sep="")
    #### D:
    duClus <- cdbTab$cluster[cdbTab$numberEss/cdbTab$uniqueEss>1.2&cdbTab$cluster!=0]
    clusterRes$cluster[cdb$cluster %in% duClus] <- paste("D",cdb$cluster[cdb$cluster %in% duClus],sep="")
    
    ### cut clusters by metagenomic coverage depth, document cut-offs
    #### initialise output file
    write.table(t(c("cluster","cutoff")),
                paste0(binny_dir,paste("/bimodalClusterCutoffs",pk,nn,"tsv",sep=".")),
                sep="\t",col.names=F,row.names=F,quote=F)
    ### put each cluster with detected duplications into depth-cutting:
    for(ds in unique(clusterRes$cluster[grep("D",clusterRes$cluster)])){
      clusterRes <- muClus(ds,clusterRes,paste0(binny_dir,paste("/bimodalClusterCutoffs",pk,nn,"tsv",sep=".")),1)
    }
    
    ###################################################################################################
    ## binny - iterations, for bins with more than 20% duplicated essential genes
    ###################################################################################################
    
    ### 2nd iteration:
    #### - reachability estimates are based on nth nearest neighbour (independent of number of neighbouring points)
    #### - neighbouring points are increased by 2
    pk2 <- pk + 2
    pdf(paste0(binny_dir,paste("/scatterPlots2",pk,nn,"pdf",sep=".")))
    clusterRes <- binny_iteration(pk2, cRes=clusterRes, cInfo=contigInfo, coco=coco, pk=pk, nn=nn, binny_dir=binny_dir, cluster_dir=cluster_dir)
    dev.off()
    # save.image(paste0(binny_dir,"/binny_WS.Rdata"))
    
    ###third iteration of clustering on clusters with more than 20% duplicated essential genes
    #### - reachability estimates are based on nth nearest neighbour (independent of number of neighbouring points)
    #### - neighbouring points are increased by 2 again
    pk4 <- pk + 4
    pdf(paste0(binny_dir,paste("/scatterPlots3",pk,nn,"pdf",sep=".")))
    clusterRes <- binny_iteration(pk4, cRes=clusterRes, cInfo=contigInfo, coco=coco, pk=pk, nn=nn, binny_dir=binny_dir, cluster_dir=cluster_dir)
    dev.off()
    # save.image(paste0(binny_dir,"/binny_WS.Rdata"))
    
    #fourth iteration of clustering on clusters with more than 20% duplicated essential genes
    ### - reachability estimates are based on nth nearest neighbor (independent of number of neighbouring points)
    ### - neighboring points are increased by 2 again
    pk6 <- pk + 6
    pdf(paste0(binny_dir,paste("/scatterPlots4",pk,nn,"pdf",sep=".")))
    clusterRes <- binny_iteration(pk4, cRes=clusterRes, cInfo=contigInfo, coco=coco, pk=pk, nn=nn, binny_dir=binny_dir, cluster_dir=cluster_dir)
    dev.off()
    # save.image(paste0(binny_dir,"/binny_WS.Rdata"))
    
    ###################################################################################################
    ## Assignment of final names to clusters
    ###################################################################################################
    ### based on completeness, the C is replaced by:
    #### T=Tier
    #### "C90P92": more than 96/107 essential genes, less than 115 essential genes in total (<8% duplicated genes, >90% complete)
    #### "C66P80": more than 70/107 essential genes (>66% complete), less than 20% in duplicate
    #### "C50P80": more than 53/107 essential genes (>50% complete), less than 20% in duplicate
    #### "C33P80": more than 35/107 essential genes (>33% complete), less than 20% in duplicate
    #### "C15P80": more than 15/107 essential genes (>14% complete), less than 20% in duplicate
    ### which leaves: 
    #### "C": at least 1/107 essential genes (>= 1% complete), less than 20% in duplicate
    #### "E": no essential genes
    #### "B": at least 1/107 essential genes (>= 1% complete), at least 20% in duplicate
    #### "N": noise
    for(clus in grep("C",unique(clusterRes$cluster),value=T)){ #iterate over all bins excep: overcomplete (B), 0void of essential genes (E), or noise (N)
      uni <- length(unique(unlist(sapply(contigInfo$essentialGene[contigInfo$essentialGene!="notEssential"&
                                                             clusterRes$cluster==clus],
                                  function(x) unlist(strsplit(x,split=";")))))) # number of different essential genes
      ess <- length(unlist(sapply(contigInfo$essentialGene[contigInfo$essentialGene!="notEssential"&
                                                                    clusterRes$cluster==clus],
                                         function(x) unlist(strsplit(x,split=";"))))) # number of all essential genes
      
      # Name bins after completion and purity estimates
      completeness <- round((uni/115)*100, digits = 0)
      purity <- round((uni/ess)*100, digits = 0)                      
      clusterRes$cluster[clusterRes$cluster == clus] <- gsub("C",
                                                             paste('C', completeness, '_P', purity, '_', sep = ''),
                                                             clus)
      
      # if(uni>96&ess<115){
      #   clusterRes$cluster[clusterRes$cluster == clus] <- gsub("C","C90P92",clus)
      # } else if(uni>70){
      #   clusterRes$cluster[clusterRes$cluster == clus] <- gsub("C","C66P80",clus)
      # } else if(uni>53){
      #   clusterRes$cluster[clusterRes$cluster == clus] <- gsub("C","C50P80",clus)
      # } else if(uni>35){
      #   clusterRes$cluster[clusterRes$cluster == clus] <- gsub("C","C33P80",clus)
      # } else if(uni>15){
      #   clusterRes$cluster[clusterRes$cluster == clus] <- gsub("C","C15P80",clus)
      # }
    }
    
    ###################################################################################################
    ## Save results
    ###################################################################################################
    ### change bin names that went into long refinement loops
    clusterRes$cluster <- gsub("(2.){8,}","2222.",gsub("(1.){8,}","1111.",
                                                       clusterRes$cluster))
    ### save bins to contig table for all contigs
    write.table(clusterRes,
                paste0(binny_dir,paste("/contigs2clusters",pk,nn,"tsv",sep=".")),
                sep="\t",row.names=F,quote=F)
    # write.table(clusterRes[grep("^[PGOL]",clusterRes$cluster),],
    #            "scaffold2bin.tsv",
    #            sep="\t",row.names=F,quote=F)
    
    ### save used functions and values + result in data object
    save(list=c("contigInfo","coco","pk","nn","essAnno","find_cutoff","muClus","clusterRes"),
         file=paste0(binny_dir,paste("/clusteringWS",pk,nn,"Rdata",sep=".")))
    ### save just bins to contig table as R object
    saveRDS(clusterRes,paste0(binny_dir,paste("/contigs2clusters",pk,nn,"RDS",sep=".")))
    # save.image(paste0(binny_dir,"/binny_WS.Rdata"))
    print("Saved tables. Plotting final map now.")
    
    ### plot final map
    essPal <- colorRampPalette(brewer.pal(11,"Spectral"))(109)[109:1]
    png(paste0(binny_dir,paste("/finalClusterMap",pk,nn,"png",sep=".")),width=17,height=17,pointsize=7,res=300,units = "cm")
    layout(matrix(1:2,nrow=1,ncol=2),c(5/6,1/6))
    par(mar=rep(1,4),pty="s")
    plot(contigInfo[,coco],pch=16,cex=0.25,ann=F,axes=F,bty="o")
    cI <- data.frame("essNum"=vector(mode="numeric",length=length(unique(clusterRes$cluster))),
                     "medX"=vector(mode="numeric",length=length(unique(clusterRes$cluster))),
                     stringsAsFactors=F)
    i <- 1
    for(clus in unique(clusterRes$cluster)){
      uni <- length(unique(unlist(lapply(contigInfo$essentialGene[contigInfo$essentialGene!="notEssential"&
                                                                    clusterRes$cluster==clus],
                                         function(x) unlist(strsplit(x,split=";"))))))
      cI$essNum[i] <- uni
      rownames(cI)[i] <- clus
      cI$medX[i] <- median(contigInfo[clusterRes$cluster==clus,coco[1]])
      cI$medY[i] <- median(contigInfo[clusterRes$cluster==clus,coco[2]])
      cI$high[i] <- max(contigInfo[clusterRes$cluster==clus,coco[2]])
      rm(uni)
      i <- i +1
    }
    cI2 <- cI[rownames(cI)!="N",][order(cI$essNum[rownames(cI)!="N"]),]
    cI2 <- rbind(cI[rownames(cI)=="N",],cI2)
    for(i in 1:nrow(cI2)){
      clus <- rownames(cI2)[i]
      if(clus=="N") pc <- "grey" else if(grepl("E",clus)) pc <- "black" else if(grepl("B",clus)) pc <- brewer.pal(8,"Set3")[8] else{
        pc <- essPal[cI2$essNum[i]]
      }
      points(contigInfo[clusterRes$cluster==clus,coco],cex=0.3,pch=16,col=pc)
      if(cI2$essNum[i]>15&clus!="N") text(cI2$medX[i],cI2$high[i]+0.5,labels=clus,cex=1.0,font=2,col=pc,pos=3)
      rm(pc)
    }
    par(mar=c(3,3,3,1),pty="m")
    plotcolors(cbind(c(colorRampPalette(brewer.pal(11,"Spectral"))(107),"black"),c(colorRampPalette(brewer.pal(11,"Spectral"))(107),"black")))
    # save.image(paste0(binny_dir,"/binny_WS.Rdata"))
    hi <- hist(cI$essNum,breaks=0:107,plot=F)
    a <- apply(cbind(hi$counts,hi$breaks[-length(hi$breaks)]+1),1,function(x) lines(c(0.5,(x[1]*1.9/max(hi$counts[-1]))+0.5),c(rep(x[2],2))))
    abline(h=c(101.5,71.5,51.5,31.5,1.5),col=c(rep("grey20",4),"grey80"),lty=2)
    axis(1,seq(from=0.5,to=2.4,length.out=3),labels=c(0,max(hi$counts[-1])/2,max(hi$counts[-1])))
    axis(2,seq(from=1,to=11000/100,length.out=6),labels=c(0,20,40,60,80,100),las=1)
    mtext("% completeness",2,2)
    mtext("clusters",1,2)
    dev.off()
    print("binny done")
    
    ### remove intermediate results
    rm(list=c("contigInfo","coco","pk","nn","clusterRes",
              "essPal","cI","cI2","hi","a","clus","ess",
              "emptyClus","closeClus","cdbTab","skn","sdkn","est",
              "pk2","pk4","pk6"))
    
    
    ####################################################################
    ## Save the R workspace
    ####################################################################
    
    print("START: Saving R image: binny_WS.Rdata")
    save.image(paste0(binny_dir,"/binny_WS.Rdata"))
    print("DONE: Saving R image: binny_WS.Rdata")
  }}
