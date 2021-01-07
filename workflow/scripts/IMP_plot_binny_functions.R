#!/bin/R

###################################################################################################
## This file contains the functions required for binny
###################################################################################################

### cut-off between two modes
find_cutoff <- function(data,k=2,proba=0.5) {
  model <- try(normalmixEM(x=data, k=k,maxrestarts=50,verb=F),silent=T)
  if(class(model)!="try-error"){
    i <- which.min(model$mu)
    f <- function(x) {
      proba - (model$lambda[i]*dnorm(x, model$mu[i], model$sigma[i]) /(model$lambda[1]*dnorm(x, model$mu[1], model$sigma[1]) + model$lambda[2]*dnorm(x, model$mu[2], model$sigma[2])))
    }
    lower <- min(model$mu)
    upper <- max(model$mu)
    if(f(lower)*f(upper)>0){
      return(median(data))
    }else{
      return(uniroot(f=f, lower=lower, upper=upper)$root)  # Careful with division by zero if changing lower and upper
    }
  }else{
    return(median(data))
  }
}

### recursive function to split clusters according to depth
muClus <- function(clusterName,cRes,resFile,recDepth){
  print(recDepth)
  covcol <- which(colnames(contigInfo)=="MG_depth")
  dInfo <- contigInfo[cRes$cluster==clusterName,]
  if(dip.test(log10(1+dInfo[,covcol][dInfo$essentialGene!="notEssential"]))$p.value<0.05){
    cuto <- find_cutoff(log10(1+dInfo[,covcol][dInfo$essentialGene!="notEssential"]))
    write.table(t(c(clusterName,cuto)),resFile,append=T,row.names=F,col.names=F,quote=F,sep="\t")
    subs1 <- dInfo[log10(1+dInfo[,covcol])<cuto,]
    num1 <- length(unlist(sapply(subs1$essentialGene[subs1$essentialGene!="notEssential"],function(x) unlist(strsplit(x,split=";")))))
    uni1 <- length(unique(unlist(sapply(subs1$essentialGene[subs1$essentialGene!="notEssential"],function(x) unlist(strsplit(x,split=";"))))))
    if(uni1==0) {
      cRes$cluster[cRes$contig %in% subs1$contig] <- paste(gsub("D","E",cRes$cluster[cRes$contig %in% subs1$contig]),1,sep=".")
    } else if(num1/uni1<=1.2){
      cRes$cluster[cRes$contig %in% subs1$contig] <- paste(gsub("D","C",cRes$cluster[cRes$contig %in% subs1$contig]),1,sep=".")
    } else if(recDepth < 900){
      cRes$cluster[cRes$contig %in% subs1$contig] <- paste(cRes$cluster[cRes$contig %in% subs1$contig],1,sep=".")
      cRes <- muClus(unique(cRes$cluster[cRes$contig %in% subs1$contig]),cRes,resFile,recDepth+1)
    }else{
      cRes$cluster[cRes$contig %in% subs1$contig] <- paste(gsub("D","B",cRes$cluster[cRes$contig %in% subs1$contig]),1,sep=".")
    }
    subs2 <- dInfo[log10(1+dInfo[,covcol])>=cuto,]
    num2 <- length(unlist(sapply(subs2$essentialGene[subs2$essentialGene!="notEssential"],function(x) unlist(strsplit(x,split=";")))))
    uni2 <- length(unique(unlist(sapply(subs2$essentialGene[subs2$essentialGene!="notEssential"],function(x) unlist(strsplit(x,split=";"))))))
    if(uni2==0) {
      cRes$cluster[cRes$contig %in% subs2$contig] <- paste(gsub("D","E",cRes$cluster[cRes$contig %in% subs2$contig]),2,sep=".")
    } else if(num2/uni2<=1.2){
      cRes$cluster[cRes$contig %in% subs2$contig] <- paste(gsub("D","C",cRes$cluster[cRes$contig %in% subs2$contig]),2,sep=".")
    } else if(recDepth < 900){
      cRes$cluster[cRes$contig %in% subs2$contig] <- paste(cRes$cluster[cRes$contig %in% subs2$contig],2,sep=".")
      cRes <- muClus(unique(cRes$cluster[cRes$contig %in% subs2$contig]),cRes,resFile,recDepth+1)
    } else{
      cRes$cluster[cRes$contig %in% subs2$contig] <- paste(gsub("D","B",cRes$cluster[cRes$contig %in% subs2$contig]),2,sep=".")
    }
  } else {
    cRes$cluster[cRes$contig %in% dInfo$contig] <- gsub("^.","B",cRes$cluster[cRes$contig %in% dInfo$contig])
  }
  return(cRes)
}

### one binny iteration:
binny_iteration <- function(pk_i, cRes=clusterRes, cInfo=contigInfo, coco=coco, pk=pk, nn=nn, binny_dir=binny_dir, cluster_dir=cluster_dir){
  for(bb in unique(cRes$cluster[grep("B",cRes$cluster)])){
    bbInfo <- cInfo[cRes$cluster==bb,]
    #### estimate EST for current cluster:
    if(nrow(bbInfo)>nn){ #check if there are enough contigs
      sknBB <- sort(knn.dist(bbInfo[,coco],nn)[,nn])
      if(length(sknBB>10)){ #check if there are enough contigs
        sdknBB <- runsd(sknBB,10)
        plot(sdknBB,sknBB)
        plot(sdknBB)
        abline(h=quantile(sdknBB,0.975))
        estBB <- sort(sknBB)[min(which(sdknBB>quantile(sdknBB,0.975)&sknBB>median(sknBB)))]
        if(is.na(estBB)){
          estBB <- sort(sknBB)[max(sdknBB)]
          print(paste(bb,": EST could not be estimated."))
        }
        plot(sknBB)
        abline(h=estBB)
        
        #### record EST:
        write.table(t(c(bb,estBB)),
                    paste0(binny_dir,paste("/reachabilityDistanceEstimates",pk,nn,"tsv",sep=".")),
                    row.names=F,col.names=F,quote=F,sep="\t",append=T)
        #### run DBscan on current cluster:
        BBcdb <- dbscan(bbInfo[,coco],estBB,pk_i)
        
        #### plot current cluster with new subclusters:
        plot(bbInfo[,coco],pch=16,cex=0.25,ann=F)
        j <- 1
        BBcdbTab <- data.frame("cluster"=names(table(BBcdb$cluster)),"contigs"=0,"numberEss"=0,"uniqueEss"=0,stringsAsFactors=F)
        for(i in names(table(BBcdb$cluster))) { #iterate over subclusters
          points(bbInfo[BBcdb$cluster==i,coco],pch=16,cex=0.4,col=c("grey",colors(distinct=T))[j])
          BBcdbTab$contigs[BBcdbTab$cluster==i] <- table(BBcdb$cluster)[names(table(BBcdb$cluster))==i]
          BBcdbTab$numberEss[BBcdbTab$cluster==i] <- length(unlist(sapply(bbInfo$essentialGene[bbInfo$essentialGene!="notEssential"&
                                                                                                 BBcdb$cluster==i],
                                                                          function(x) unlist(strsplit(x,split=";")))))
          BBcdbTab$uniqueEss[BBcdbTab$cluster==i] <- length(unique(unlist(sapply(bbInfo$essentialGene[bbInfo$essentialGene!="notEssential"&
                                                                                                        BBcdb$cluster==i],
                                                                                 function(x) unlist(strsplit(x,split=";"))))))
          j<-j+1
        }
        #### document clusters with simplified names (if they went into loops)
        write.table(BBcdbTab,
                    paste0(cluster_dir,paste("/blob1Cluster",pk,nn,gsub("(2.){8,}","2222.",gsub("(1.){8,}","1111.",
                                                                                                bb)),"tsv",sep=".")),
                    sep="\t",row.names=F,quote=F)
        
        #### assign bins to contigs
        ##### use some temporary dataframes: for noise
        BBclusterRes <- data.frame("contig"=bbInfo$contig,"cluster"="x",stringsAsFactors=F)
        BBclusterRes$cluster[BBcdb$cluster==0] <- "N" # noise as N
        ##### write into table that should be returned in the end
        cRes$cluster[cRes$contig %in% BBclusterRes$contig[BBclusterRes$cluster=="N"]] <- "N"
        
        ##### use some temporary dataframes: for clusters without essential genes
        BBemptyClus <- BBcdbTab$cluster[BBcdbTab$uniqueEss==0&BBcdbTab$cluster!=0]
        BBclusterRes$cluster[BBcdb$cluster %in% BBemptyClus] <- paste("E",BBcdb$cluster[BBcdb$cluster %in% BBemptyClus],sep="")
        BBemptyCont <- BBclusterRes$contig[BBcdb$cluster %in% BBemptyClus]
        ##### write into table that should be returned in the end
        cRes$cluster[cRes$contig %in% BBemptyCont] <- paste(gsub("B","E",bb),
                                                            gsub("E","",unlist(sapply(cRes$contig[cRes$contig %in% BBemptyCont],
                                                                                      function(x)BBclusterRes$cluster[BBclusterRes$contig==x]))),
                                                            sep=".")
        
        ##### use some temporary dataframes: for clusters with essential genes and less than 20% duplication
        BBcloseClus <- BBcdbTab$cluster[BBcdbTab$numberEss/BBcdbTab$uniqueEss<=1.2&BBcdbTab$cluster!=0]
        BBclusterRes$cluster[BBcdb$cluster %in% BBcloseClus] <- paste("C",BBcdb$cluster[BBcdb$cluster %in% BBcloseClus],sep="")
        BBcloseCont <- BBclusterRes$contig[BBcdb$cluster %in% BBcloseClus]
        if(length(BBcloseCont>0)) print(paste("new cluster in",bb,"\n"))
        ##### write into table that should be returned in the end
        cRes$cluster[cRes$contig %in% BBcloseCont] <- paste(gsub("B","C",bb), #replace old B with C, as they're not anymore overcomplete
                                                            gsub("C","",unlist(sapply(cRes$contig[cRes$contig %in% BBcloseCont],
                                                                                      function(x)BBclusterRes$cluster[BBclusterRes$contig==x]))),
                                                            sep=".")
        
        ##### use some temporary dataframes: for clusters with more than 20% duplication of essential genes
        BBduClus <- BBcdbTab$cluster[BBcdbTab$numberEss/BBcdbTab$uniqueEss>1.2&BBcdbTab$cluster!=0]
        BBclusterRes$cluster[BBcdb$cluster %in% BBduClus] <- paste("D",BBcdb$cluster[BBcdb$cluster %in% BBduClus],sep="")
        BBduCont <- BBclusterRes$contig[BBcdb$cluster %in% BBduClus]
        uniD <- unique(paste(gsub("B","D",bb),
                             gsub("D","",unlist(sapply(clusterRes$contig[clusterRes$contig %in% BBduCont],
                                                       function(x)BBclusterRes$cluster[BBclusterRes$contig==x]))),
                             sep="."))
        ##### write into table that should be returned in the end
        cRes$cluster[cRes$contig %in% BBduCont] <- paste(gsub("B","D",bb), #replace old B for new D to go into splitting by depth
                                                         gsub("D","",unlist(sapply(cRes$contig[cRes$contig %in% BBduCont],
                                                                                   function(x)BBclusterRes$cluster[BBclusterRes$contig==x]))),
                                                         sep=".")
        for(bds in uniD){
          cRes <- muClus(bds,cRes,paste0(binny_dir,paste("/bimodalClusterCutoffs",pk,nn,"tsv",sep=".")),1) #overwriting cRes
        }
      }else{
        print(paste("Cluster",bb,"too small for knn-approach.\n"))
      }
    }else{
      print(paste("Cluster",bb,"too small for knn-approach.\n"))
    }
  }
  cRes #return result for all clusters
}
