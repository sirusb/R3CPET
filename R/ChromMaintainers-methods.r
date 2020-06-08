.valid.ChromMaintainers <- function(x){
  
  if(class(x@maintainers) != "HLDAResult")
    return("the maintainers slot should an HLDAResult object");
  
  if(!is.matrix(x@topEdges))
    return("topEdges should be a matrix")
  
  if(!is.matrix(x@topNodes))
    return("topNodes should be a matrix")
  
  if(ncol(x@topEdges) != ncol(x@topNodes))
      return("topEdges and topNodes should have the same number of networks")
      
  if(!is.list(x@networks)){
    return("networks should be a list")   
  }    
  else{
    ## if there are some elements check if they are of class igraph
    if(length(x@networks) >0){
      alligraph <- all(sapply(x@networks, function(elem) class(elem) == "igraph"))
      
      if(!alligraph)
        return("the networks slot should be a list of igraph objects")
    }
  }
                                 
  return(TRUE)
}
                              
.ConvertToHDA<-function(Nets,tfspace){
  
  #tenPercent<- floor(length(Nets)/10);
  Documents<-list();
  
  for(i in 1:length(Nets)){
    termCount<-length(unique(Nets[[i]]));
    
    if(termCount >0){
      counts<-table(Nets[[i]]);
      pos<-match(names(counts), tfspace);
      ord<-order(pos);
      netName<-names(Nets)[i];
      if(is.null(netName)){ netName<-i; }
      counts<-counts[ord];
      names(counts)<-pos[ord];
      res<-matrix(0,nrow=2,ncol=length(counts))
      res[1,]<-as.numeric(names(counts))-1; #Should be 0 indexed 
      res[2,]<-counts
      Documents[[netName]]<-res;
      
    }
        
  }
  
  return(Documents);
}


.plot.clustOrderHeatmap<-function(cluster,data, path, W=2048, H=1024){
  elementsOrder<-order(cluster);
  data<-data[elementsOrder,];
  cluster<-cluster[elementsOrder];
  
  lbls<-paste("cluster",as.numeric(sort(unique(cluster))),sep="");
  annot<-data.frame(Cluster=factor(cluster,labels = lbls));  
  rownames(annot)<-rownames(data);
  
  Var1<- sample(colors()[100:700],length(lbls))
  names(Var1)<-lbls;
  ann_colors<-list(Cluster=Var1)
  
  if(path != ""){
    filename <- file.path(path,"ClusHeatmap.png");
    png(filename,height=H,width=W)  
  }
  
  p <- pheatmap(t(data), cluster_rows=FALSE,cluster_cols=FALSE,border_color=NA, show_colnames=FALSE,
                color= colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")) (200),
                annotation= annot, annotation_colors=ann_colors)
  p
  if(path != ""){
    dev.off()    
  }
  return(p)
}


.plot.TopicsSimilarity <- function(topics){

  simMat<-matrix(0,nrow=ncol(topics),ncol=ncol(topics));
  
  for(i in 1:(ncol(topics)-1)){
    for(j in (i+1):ncol(topics)){
      simMat[i,j]<-length(intersect(topics[,i],topics[,j]))/ length(unique(c(topics[,i],topics[,j])));
      simMat[j,i]<-simMat[i,j];
    }
  }
  colnames(simMat)<-paste("Network",1:ncol(topics),sep="")
  rownames(simMat)<-paste("Network",1:ncol(topics),sep="")
  diag(simMat) <- 1
  p <- pheatmap(simMat)  
  invisible( list(plot=p, simMat = simMat) )
}

## Get the the list of gene promoters near to the given regions
.GetClusterInfo<-function(object){

  if(!class(object) %in% "GRanges"){
    stop("object should be of class GRanges");
  }
  
  requireNamespace("TxDb.Hsapiens.UCSC.hg19.knownGene")
  requireNamespace("org.Hs.eg.db")
  
  hg19.known<-TxDb.Hsapiens.UCSC.hg19.knownGene;  
  genesPromoters<-subsetByOverlaps(promoters(hg19.known,2500,2500),object)
  
  ucscIDs<-elementMetadata(genesPromoters)$tx_name	
  ucsc2Entrez<-toTable(org.Hs.egUCSCKG);
  pos<-match(ucscIDs,ucsc2Entrez$ucsc_id);
  pos<-pos[!is.na(pos)];
  EntrezIDs<-ucsc2Entrez$gene_id[pos];
  #hgnc <- EntrezToHGNC(EntrezIDs)
  res<- select(org.Hs.eg.db, EntrezIDs, c("GENENAME", "SYMBOL"))		
  res<-res[!duplicated(res),]
  return(res);  
  #hgnc <- hgnc[!(duplicated),]
  #return(hgnc)
}

## Creates a directory that contains a the list of genes for each cluster
.get.ClusterInvolvedGenes<-function(hdaRes,data,path="ClustersGenes"){
    
  message(paste("creating directory",path))
  dir.create(path, showWarnings = FALSE)
  
  clus<-sort(unique(getClusters(hdaRes)));
  message("processing clusters ....")
  for(i in clus){
    clusRegions <-  getRegionsIncluster(hdaRes,data, cluster= i)
    clusInfo<-.GetClusterInfo(clusRegions);
    fname<-file.path(path, paste( c("cluster",i,"_genes.txt"),collapse="") )
    write.table(clusInfo,file=fname,row.names=FALSE,quote=FALSE, sep="\t")
  }
}

.get.NetworksGenes<- function(hdaRes, data, path){
    
  requireNamespace("ChIPpeakAnno")
  data("TSS.human.GRCh37", package="ChIPpeakAnno", envir=environment())
  TSS.human.GRCh37 <- get("TSS.human.GRCh37", envir= environment())

  message(paste("creating directory",path))
  dir.create(path, showWarnings = FALSE)
  
  nets <-1:ncol(topEdges(hdaRes))
  message("processing networks ....")
  for(net in nets){    
    
    NetworkRegions <- getRegionsInNetwork(hdaRes,data,net)
    if(length(NetworkRegions)>0){
      #networkInfo<-.GetClusterInfo(NetworkRegions);
      tmp.rd <- GRanges( gsub("chr","", as.character(seqnames(NetworkRegions))), 
                        IRanges(start(NetworkRegions), end(NetworkRegions)))
        
      ##ensemble <- useMart("ensembl")
      ## hsp <- useDataset(mart = ensemble, dataset = "hsapiens_gene_ensembl")
      
      tmp.anno <- annotatePeakInBatch(tmp.rd,
                                      featureType = "TSS",
                                      AnnotationData = TSS.human.GRCh37,
                                      ## mart= hsp,
                                      select = "all", 
                                      PeakLocForDistance = "middle")
      
      res <- with(tmp.anno, { subset(as.data.frame(tmp.anno),
             abs(distancetoFeature) <= 2500 | insideFeature =="includeFeature" )})
      
      if(nrow(res)>0){
        conv <- EnsemblToHGNC(res$feature)
        pos <- match(res$feature, conv$ensembl_gene_id)
        
        res$name <- conv$hgnc_symbol[pos]
        res$space <- paste("chr",as.character(res$space))
        res <- res[,c("feature","name")]      
        fname<-file.path(path, paste( c("Network",net,"_genes.txt"),collapse="") )
        write.table(res,file=fname,row.names=FALSE,quote=FALSE, sep="\t")
      }
    }
  }
}

.buildNetFromEdges<-function(edgesList){

  g<-graph.empty()
  for(e in edgesList){
    vertices <-unlist(strsplit(e,split="_"))
    NotIn <- which(!vertices %in% V(g)$name)
    #Add nodes
    if(length(NotIn)>0){
      for(n in NotIn){
        g <- g+ vertex(vertices[n])
      }
    }  	
    #Add edges
    g[vertices[1],vertices[2],directed=FALSE]<-1		
  }
  g<-as.undirected(g)
  return(g)
}

## TO Draw basier like edges
## inspired from http://is-r.tumblr.com/post/38459242505/beautiful-network-diagrams-with-ggplot2

.edgeMaker <- function(whichRow, len = 100, curved = TRUE,adjacencyList,layoutCoordinates){
  
  fromC <- as.matrix( layoutCoordinates[adjacencyList[whichRow, 1],1:2 ] ) # Origin
  toC <- as.matrix( layoutCoordinates[adjacencyList[whichRow, 2],1:2 ] )  # Terminus
  
  # Add curve:
  graphCenter <- colMeans(layoutCoordinates[,1:2])  # Center of the overall graph
  bezierMid <- unlist(c(fromC[1], toC[2]))  # A midpoint, for bended edges
  distance1 <- sum((graphCenter - bezierMid)^2)
  if(distance1 < sum((graphCenter - unlist(c(toC[1], fromC[2])))^2)){
    bezierMid <- c(toC[1], fromC[2])
  }  # To select the best Bezier midpoint
  bezierMid <- (fromC + toC + bezierMid) / 3  # Moderate the Bezier midpoint
  if(curved == FALSE){bezierMid <- (fromC + toC) / 2}  # Remove the curve
  
  edge <- data.frame(bezier(c(fromC[1], bezierMid[1], toC[1]),  # Generate
                            c(fromC[2], bezierMid[2], toC[2]),  # X & y
                            evaluation = len))  # Bezier path coordinates
  edge$Sequence <- 1:len  # For size and colour weighting in plot
  edge$Group <- paste(adjacencyList[whichRow, 1:2], collapse = ">")
  return(edge)
}

## inspired from http://is-r.tumblr.com/post/38459242505/beautiful-network-diagrams-with-ggplot2
.plotNetwork<-function(g,layout.fct=layout.kamada.kawai, title=""){
  
  plotcord <- data.frame(layout.fct(g) )
  petnet<-NULL;
  
  if(nrow(get.edgelist(g))>0){
    edglist <- melt(as.matrix(get.adjacency(g)))
    edglist <- edglist[edglist$value > 0, ]
    edglist[,1]<-factor(edglist[,1],levels=V(g)$name)
    edglist[,2]<-factor(edglist[,2],levels=V(g)$name)
    edges <- data.frame(plotcord[edglist[,1],], plotcord[edglist[,2],])
    colnames(edges) <-  c("X1","Y1","X2","Y2")
    edges$midX  <- (edges$X1 + edges$X2) / 2
    edges$midY  <- (edges$Y1 + edges$Y2) / 2    
    
    allEdges <- lapply(1:nrow(edglist), .edgeMaker, len = 500, curved = TRUE, 
                       adjacencyList= edglist, layoutCoordinates= plotcord)
    allEdges <- do.call(rbind, allEdges)
        
       
    pnet <- with(allEdges,{ ggplot(allEdges) + geom_path(aes(x = x, y = y, group = Group,size = -Sequence, colour= Sequence))})
    
  }
  else{
    pnet <- ggplot() 		
  }
  plotcord$type <- as.factor(V(g)$type)
  plotcord$name <-V(g)$name;
  
  pnet <- pnet + geom_point(aes_string(x='X1', y='X2'), data=plotcord,size = 10, pch=21, color="#e34a33", fill="#fdbb84") 
  pnet <- pnet + scale_colour_gradient(low = gray(0), high = gray(9/10), guide = "none")
  pnet <- pnet + geom_text(data = plotcord, aes_string(x='X1',y='X2', label = 'name'),size=2,family="Courier", fontface="bold")  
  pnet <- pnet + scale_size(range = c(1/10, 1), guide = "none") 
  pnet <- pnet + theme(panel.background = element_blank()) # + theme(legend.position="none")
  pnet <- pnet + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) 
  pnet <- pnet + theme( legend.background = element_rect(colour = NA))
  pnet <- pnet + theme(panel.background = element_rect(fill = "white", colour = "black")) 
  pnet <- pnet + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
  pnet <- pnet + ggtitle(title)
    
  return(pnet);
}

.plotAllNet <- function(networks,layoutfct=layout.kamada.kawai, file="AllGraphs.pdf"){
  
  if(length(networks) == 0 || !all(sapply(networks, is.igraph)))
    stop("the networks member should a list of igraph object")
      
  message("plotting networks")
  plots <- list()
  if(file != ""){
    message(paste("plots will be also available on file", file))
    pdf(file, width=14, height=7)
  }
  for(i in 1:length(networks)){
    subplot<-.plotNetwork(networks[[i]],layoutfct, paste("Network",i));  	
    #print(subplot, vp=vplayout( ceiling(i/netPerRow), ifelse(i %% netPerRow==0,netPerRow,i %% netPerRow)) );  	
    if(file != "") plot(subplot)    
    plots[[i]] <- subplot    
  }
  
  if(file != ""){
    dev.off()
  }
  return(plots)
}

.annotateNodesExpression<-function(graphs,RPKMS){
  
  if(length(graphs) == 0)
    stop("Please generate the igraph objects first. Check the GenerateNetworks method.")
  for(i in 1:length(graphs)){
    g<-graphs[[i]]
    pos<- match(V(g)$name,RPKMS[,1])
    V(g)[!is.na(pos)]$RPKM<-RPKMS[pos[!is.na(pos)],2]
    graphs[[i]]<-g
  }
  
  return(graphs)
}

###################################################################################
##
##        ChromatinMaintainers-methods
##
####################################################################################


setMethod("clusterInteractions", signature = c(object="ChromMaintainers"),
          function(object, method=c("clues","sota"), nbClus=20 ){
           
			cat("clusterInteractions : checking\n")
            if(is.null(object@maintainers@docPerTopic) || 0 %in% dim(object@maintainers@docPerTopic))
              stop("The docPerTopic matrix should not be empty")
            
			cat("clusterInteractions : reading args\n")			
            method <- match.arg(method)            
                        
            if(method == "clues"){
			  cat("using clues\n")
              requireNamespace("clues")
			 }
            else{
              if(nbClus <= 0 || is.null(nbClus))
                stop("nbClus should be a positive number")              
			   cat("using clValid\n")
              requireNamespace("clValid")
            }
            
            clusRes <- switch(method,
                              clues = clues(object@maintainers@docPerTopic) ,
                              sota = sota(object@maintainers@docPerTopic,maxCycles=nbClus-1) )
            object@clusRes <- clusRes;
            message(paste("DNA interactions have been clustered into",length(unique(getClusters(object))), "cluster"))
            return(object)
          })

setMethod("plot3CPETRes", signature = c(object="ChromMaintainers"),
          function(object, path="", W=14, H=7 , type=c("heatmap","clusters","curve","avgCurve","netSim", "networks"),
                   byEdge=TRUE, layoutfct=layout.kamada.kawai, ...){
            
            type <- match.arg(type)
            ## we can do getClusters(object) but sometime we get 
            ## some weired erros.
            clusters<- getClusters(object) #@clusRes$mem
            
            p<- NULL;
            if(type== "heatmap"){
              if(is.null(slot(object@maintainers,"docPerTopic")))
                stop('No infered networks were found, please check the Method InferNetworks')
             p <- .plot.clustOrderHeatmap(clusters, slot(object@maintainers,"docPerTopic"), path, W, H)
            }
            else{
              if(type %in% c("clusters","curve","avgCurve") && is.null(object@clusRes))
                stop("No clustering results found, please check method cluster")
              
              par(mar = rep(2, 4))
              
              if(type == "clusters"){                
                if("sota" %in% class(object@clusRes))
                  stop("clusters plot are only supported for clues objects")
                else
                 p<- plotClusters(object@clusRes$y,object@clusRes$mem)
              }
              else{
                if(type == "curve"){
                  if("sota" %in% class(object@clusRes))
                    p<- plot(object@clusRes)
                  else
                    p <- plotCurves(object@clusRes$y,object@clusRes$mem)
                }
                else
                  if(type == "avgCurve"){
                    if("sota" %in% class(object@clusRes)){
                       message("curves and avgCurves are plotted for clues objects")
                       p <- plot(object@clusRes)
                    }
                    else
                      p <- plotAvgCurves(object@clusRes$y,object@clusRes$mem)
                  }
                else{
                  if(type == "netSim"){
                    if(byEdge){
                      p <- .plot.TopicsSimilarity(slot(object,"topEdges"))
                    }
                    else
                      p <- .plot.TopicsSimilarity(slot(object,"topNodes"))
                  }
                  else{
                    if(type== "networks"){
                      library(ggplot2)
                     if(path == "") path = "AllGraphs.pdf"                     
                     p<-  .plotAllNet(networks(object), layoutfct,path)
                    }
                  }
                    
                }
              }              
            }
                                
            invisible(p)
          })

setMethod("getClusters",signature= c(object= "ChromMaintainers"),
          function(object){
            if(is.null(object@clusRes))
              return(NA)
            
            clusters <- c();
            if("clues" %in% class(object@clusRes)){              
              clusters <- object@clusRes$mem
            }
            else{              
              clusters <- object@ClusRes$clust
            }            
            return(clusters)
            
          })

setMethod("getRegionsIncluster", signature=c(hdaRes="ChromMaintainers", data="ChiapetExperimentData",
                                             cluster="numeric"),
          function(hdaRes,data, cluster=1, ...){
                      
            if(is.null(hdaRes@clusRes))
              stop("You need to do the clustering first, check the cluster method")
            
            clusters <- hdaRes@clusRes$mem
            clusElements<-which(clusters == cluster);          
            if(length(clusElements) <= 0){
              warning("The provided cluster does not exist")              
              return(new("GRanges"))
            }
            
            petNames <-rownames(slot(hdaRes@maintainers,"docPerTopic"))
            
            regionRoot <- paste("PET#",petNames[clusElements],sep="")
            
            pos <- which( gsub("\\.\\d|PET#","",pet(data)$PET_ID) %in% petNames[clusElements])
            return(pet(data)[pos])                  
})


setMethod("getRegionsInNetwork", signature=c(hdaRes="ChromMaintainers", data="ChiapetExperimentData",
                                             net="numeric"),
          function(hdaRes,data, net=1,thr=0.5, ...){
            
            if(is.null(hdaRes@clusRes))
              stop("You need to do the clustering first, check the cluster method")
            
            
            if(ncol(slot(hdaRes@maintainers,"docPerTopic")) < net){
              warning("The provided net does not exist")              
              return(new("GRanges"))
            }
            
			      maxes <- apply(slot(hdaRes@maintainers,"docPerTopic"),1,function(x) which(x==max(x))[1])
            topInter <- which(maxes == net)
            #topInter <- which(slot(hdaRes@maintainers,"docPerTopic")[,net] >= thr)
            
            petNames <-rownames(slot(hdaRes@maintainers,"docPerTopic"))[topInter]
            regionRoot <- paste("PET#",petNames,sep="")      
            pos <- which( gsub("\\.\\d|PET#","",pet(data)$PET_ID) %in% petNames)
            return(pet(data)[pos])                  
          })

setMethod("outputGenesPerClusterToDir", signature=c(hdaRes="ChromMaintainers", data="ChiapetExperimentData"),
          function(hdaRes,data,path="ClustersGenes", ...){
            .get.ClusterInvolvedGenes(hdaRes,data,path)
})

setMethod("outputGenesPerNetworkToDir", signature=c(hdaRes="ChromMaintainers", data="ChiapetExperimentData"),
          function(hdaRes, data, path="NetworksGenes", ...){
           .get.NetworksGenes(hdaRes, data, path) 
          })

## TODO: Don't use a lot of packages 
setMethod("visualizeCircos", signature = c(object= "ChromMaintainers", data= "ChiapetExperimentData",  cluster="numeric"),
          function(object, data, cluster = 1, chrLenghts = NULL){
            
			requireNamespace("biovizBase")
			requireNamespace("ggbio")
            interactions <- getRegionsIncluster(object, data, cluster = cluster)
            
            if(length(interactions) ==0)
              return(NA)
            if(is.null(chrLenghts)){
              data("hg19Ideogram", package = "biovizBase", envir= environment())            
              hg19Ideogram <- get("hg19Ideogram", envir = environment())
              hg19Ideogram <- hg19Ideogram[ as.character(seqnames(hg19Ideogram)) %in% seqlevels(interactions) ]
              hg19Ideo <- keepSeqlevels(hg19Ideogram, seqlevels(interactions))
              seqlengths(interactions) <- seqlengths(hg19Ideo)
            }
            else{
              if(length(names(chrLenghts)) == 0 || !is.numeric(chrLenghts))
                stop("chrLenghts should be a names numeric vector")
              
              if(! all(seqlevels(interactions) %in% names(chrLenghts)) )
                 stop("some chromosomes are missing from chrLenghts")
                 
              pos <- match(seqlevels(interactions), names(chrLenghts))
              
              seqlengths(interactions) <- chrLenghts[pos]
            }
            
            ## get left-side interactions
            leftID <- grep("PET#\\d+\\.1",interactions$PET_ID)
            RightID <- grep("PET#\\d+\\.2",interactions$PET_ID)
            
            circos <- interactions[leftID]
            values(circos)$to.gr <- interactions[RightID]
            
            
            p <- ggplot() + layout_circle(hg19Ideo, geom = "ideo", fill = "#9ecae1", color="#636363", radius = 30,trackWidth = 4)
            p <- p + layout_circle(hg19Ideo, geom = "text", aes(label = seqnames), 
                                   vjust = 0,radius = 32, trackWidth = 7)
            
            p <- p + layout_circle(circos, geom = "link", linked.to = "to.gr",radius = 29, trackWidth = 1, color="#f03b20")            
            p <- p + ggtitle(paste("Interactions in cluster", cluster)) + 
                     theme(plot.title = element_text(lineheight=.8, face="bold"))
            plot(p)
            invisible(list(circos = circos,plot = p))
          })

setMethod("topEdges", signature = c(object = "ChromMaintainers"),
          function(object){
            return(object@topEdges)
          })

setMethod("topNodes", signature = c(object = "ChromMaintainers"),
          function(object){
            return(object@topNodes)
          })

setMethod("networks", signature= c(object= "ChromMaintainers"), 
          function(object){
            return(slot(object,"networks"))
          })

setMethod("updateResults", signature=c(object="ChromMaintainers",nets="NetworkCollection", thr="numeric"),
          function(object,nets,thr=0.5){
            if(!is.null(nets)){
                object@topEdges<- .print.topwords(object@maintainers@wordsPerTopic,as.matrix(TF(nets)),thr)
                object@topNodes <- .getNodesList(object@topEdges) 
                ## if the networks were previously generated then update them
                if(length(object@networks) > 0){
                    object<- GenerateNetworks(object)
                }
                return(object)
            }
            else{
                warning("a NetworkCollection object should be specified")
            }
          }
         )            

setMethod("GenerateNetworks", signature = c(object = "ChromMaintainers"),
           function(object,...) {
             
             ## if one of the dimensions is zero we consider it as non valid
             if(0 %in% dim(topEdges(object)) || is.null(topEdges(object)) )
               stop("No topEdge reults found")
             
             motifs <- colnames(wordsPerTopic(object@maintainers))             
             subgraphs<-list()
             topics <- topEdges(object)
             for(i in 1:ncol(topics)){
               edgesList<-unique(topics[,i])
               edgesList<- edgesList[edgesList!=""]
               g <- .buildNetFromEdges(edgesList)
               V(g)$type="co-factor";
               tfs<-which(V(g)$name %in% motifs);
               if(length(tfs)>0){
                 V(g)$type[tfs]<-"TF";
               }
               subgraphs[[i]]<-g;
             } 
             names(subgraphs) <- paste("Network", 1:length(subgraphs),sep="")
             object@networks <- subgraphs;
             return(object)
           })


setMethod("annotateExpression", signature= c(object = "ChromMaintainers", RPKMS = "data.frame"),
          function(object, RPKMS){
            if(ncol(RPKMS) <2)
              stop("a data.frame with at least 2 columns should be provided")
            
            if( all(is.na(as.numeric(RPKMS[,2]) )) )
                stop("The second column should have numeric values")
            
            object@networks <-  .annotateNodesExpression(networks(object),RPKMS)
            return(object)
          })


setMethod("show", signature=c(object="ChromMaintainers"),
          function(object){
            cat("class:", class(object),"\n")
            cat("HLDA Results:\n")
            cat("------------\n")
            print(object@maintainers)            
          })
                
setValidity("ChromMaintainers",.valid.ChromMaintainers)
                
## An S3 user freindly method
ChromMaintainers<- function( maintainers,topEdges,topNodes, clusRes = NULL, networks = list()){    
  return( new("ChromMaintainers", maintainers = maintainers,topEdges = topEdges,
              topNodes = topNodes, clusRes = clusRes, networks = networks) )
}
