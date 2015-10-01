## This method uses the David webserive to annotate the genes controlled by each cluster
## between each request the method sleeps for 5 secs
.DavidGOAnalysis.folder <- function(folder,pval=0.05,GOlimit=10){
  GOFolder <- file.path(folder,"GO")
  dir.create(GOFolder ,showWarnings=FALSE)
  files<-list.files(folder,pattern=".txt",full.names=TRUE);
  dav<-list()
  for(f in files){
    message(paste("GO enrichment for file", f))
    genesList<-read.csv(f,sep="\t");		
    DavRes<-.DoDAVIDGOAnnotation.list(as.character(unique(genesList[,1])),pval,GOlimit);
    pngName<-gsub(".txt","_GO.png",f);
    
    dav[[basename(f)]]<-DavRes;
    if(nrow(DavRes)){      
      with(DavRes,{
	      ggplot(DavRes,aes(y=-log(pvalue),x=GO, fill=GO))+ 
		geom_bar(position="identity") + coord_flip();
	})
      ggsave(pngName,width=50.8,height=28.575,units="cm")
      fname = file.path(GOFolder, basename(gsub(".txt","_GO.csv",f)) )
      write.table(DavRes,sep="\t",row.names=FALSE,quote=FALSE,col.names=TRUE,file=fname)
    }
    Sys.sleep(5)		
  }
  return(dav)
}


.DoDAVIDGOAnnotation.list <-function(genesList,pval=0.05,GOlimit=5,IDtype="ENTREZ_GENE_ID"){
  
  
  #Groups<-c();
  Annot<-c();
  pvalue<-c();
  counts<-c();  
  genes<-c();
  genesList<-genesList[genesList!=""];
  
  #Get GO enrichment from DAVID
  dav<-.DAVIDQuery(ids=genesList,type=IDtype,annot="GOTERM_BP_ALL",tool="chartReport",URLlengthLimit = 600400);
  
  if(length(dav$downloadFileName)>0){
    dav$DAVIDQueryResult$Benjamini<-as.double(dav$DAVIDQueryResult$Benjamini);
    sigGO<-with(dav,{ subset(dav$DAVIDQueryResult[-1,], Benjamini <= pval)});
    
    #if we have some significant annotation
    if(nrow(sigGO)>0){
      ord<-order(sigGO$Benjamini)				
      sigGO<-sigGO[ord,]                
      if(nrow(sigGO)>GOlimit ) {sigGO<-sigGO[1:GOlimit,] }                               				
      Annot<- paste(sigGO$Term,"(",as.integer(sigGO$Count),")",sep="")
      sigGO$Count<-as.integer(sigGO$Count)
      counts<-sigGO$Count/length(genesList)
      #gr<-paste(clus,"(",length(genesList),")",sep="")
      #Groups<-c(Groups,rep(gr,length(sigGO$Term)));
      pvalue<-c(pvalue,sigGO$Benjamini);
      genes<- c(genes, sigGO$Genes)
    }           
  }    
  
  #result<-data.frame(Cluster=factor(Groups),GO=Annot,pvalue=as.double(pvalue),Count=as.double(counts));
  result<-data.frame(GO=Annot,pvalue=as.double(pvalue),Count=as.double(counts), Genes=genes);
  return(result);
}


.DoDAVIDGOAnnotation <-function(topics,pval=0.05,GOlimit=5){
  
  clusters<-1:ncol(topics)
  
  Groups<-c();
  Annot<-c();
  pvalue<-c();
  counts<-c();
  for(clus in clusters){
    message(paste("GO enrichment for network", clus))
    #Get the genes in that cluster
    inCluster<-as.character(topics[,clus]);
    inCluster<-inCluster[inCluster!=""];
    
    #Get GO enrichment from DAVID
    dav<- .DAVIDQuery(ids=inCluster,type="OFFICIAL_GENE_SYMBOL",annot="GOTERM_BP_ALL",tool="chartReport",URLlengthLimit = 102400);
    
    if(length(dav$downloadFileName)>0){
      dav$DAVIDQueryResult$Benjamini<-as.double(dav$DAVIDQueryResult$Benjamini);
      sigGO<-with(dav, {subset(dav$DAVIDQueryResult[-1,], Benjamini <= pval)})
      
      #if we have some significant annotation
      if(nrow(sigGO)>0){
        if(nrow(sigGO)>GOlimit ) {sigGO<-sigGO[1:GOlimit,] }
        #anot<-as.character( paste( sigGO$Term, "(" , sigGO$Count, ")", sep="") )
        #Annot<-c(Annot,anot)          
        ord<-order(sigGO$Benjamini)
        sigGO<-sigGO[ord,]
        Annot<-c(Annot,sigGO$Term)
        sigGO$Count<-as.integer(sigGO$Count)
        counts<-c(counts,sigGO$Count/length(inCluster))
        gr<-paste(clus,"(",length(inCluster),")",sep="")
        Groups<-c(Groups,rep(gr,length(sigGO$Term)));
        pvalue<-c(pvalue,sigGO$Benjamini);
      }           
    }
    Sys.sleep(5)
  }
  
  result<-data.frame(Cluster=factor(Groups),GO=Annot,pvalue=as.double(pvalue),Count=as.double(counts));
  return(result);
}

.DavidGOAnalysis <- function(topics,pval=0.05,GOlimit=5, path=""){
    
  DavRes<-.DoDAVIDGOAnnotation(topics,pval);  
  
  if(nrow(DavRes)){
    
    p <- ggplot(DavRes, aes_string(x='Cluster', y='GO',colour='pvalue')) + 
      geom_point(aes_string(size='Count'))+ 
      scale_colour_gradient(low="red", high="blue")
    
    
    if(path != ""){
      pngName<-file.path(path,"NetworksEnrich.pdf");
      ggsave(pngName,width=50.8,height=28.575,units="cm")
    }
    plot(p)
  }                            
  return(list(GO=DavRes, plot= p))
}



#######################################################################################
##
##       S4 Methods
##
#######################################################################################


setMethod("GOEnrich.folder",signature = c(folder="character"),
          function(folder, fdr=0.05,GOlimit=20){
            res <- .DavidGOAnalysis.folder(folder,pval=fdr,GOlimit=GOlimit)
            invisible(res)
          })

setMethod("GOEnrich.networks",signature= c(object = "ChromMaintainers"),
          function(object, fdr=0.05, GOlimit= 5,path=""){          
            
            res <- .DavidGOAnalysis(slot(object,"topNodes"),fdr,GOlimit, path)            
            invisible(res)
          })
