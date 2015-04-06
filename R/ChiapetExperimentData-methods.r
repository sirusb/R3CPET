## checks and loads the pet paramter into the object x
## if pet is a GRanges it will loaded directly
## otherwise it will loaded from a file
.check.pet <- function(pet, x, IsBed, petHasHeader, dist){
  ## if it is a GRanges check it has the correct fields
  if(class(pet) == "GRanges"){
    if(! "PET_ID" %in% names(mcols(pet))){
      stop("The pet object should contain a metadata column named PET_ID") 
    }
    else{
      ## if the PET_ID column exists check that it is conform to the naming convention
      
      nbLeft  <- length( grep('PET#\\d+\\.1', pet$PET_ID ))
      nbRight <- length( grep('PET#\\d+\\.2', pet$PET_ID ))
      
      if(nbLeft != nbRight)
        stop("The the number of Left regions is diffrent from the right regions \n make sure that the regions are named correctly")
      
      ## if valid put it in the pet slot
      x@pet <- pet
    }
  }
  else{
    if(is.character(pet) && pet!="") x <- loadPETs(x,petFile = pet,IsBed= IsBed, 
                                        header = petHasHeader, dist = dist)        
  }
  
  return(x)
}

## checks and loads the tfbs paramter into the object x
## if tfbs is a GRanges it will loaded directly
## otherwise it will loaded from a file
.check.tfbs <- function(tfbs, x, tfbsHasHeader){
  ## if it is a GRanges check it has the correct fields
  if(class(tfbs) == "GRanges"){
    if(! "TF" %in% names(mcols(tfbs))){
      stop("The tfbs object should contain a metadata column named TF") 
    }    
    
    ## if valid put it in the tfbs slot
    x@tfbs = tfbs
  }
  else{
    if(is.character(tfbs) && tfbs!="") x <- loadTFBS(x, tfbsFile = tfbs, 
                                         header= tfbsHasHeader)
  }  
  return(x)
}


.ReadPET<-function(file, header=TRUE){
  
  #Read file  
  PET<-read.table(file,sep="\t",comment.char="+",header=header);
  
  #Check if the fields exist
  if(ncol(PET) < 4){
    stop("Some columns are missing, please verify the file format");
  }  
  
  #Chek if the first line is the chromosome position
  chr<-grep("chr",PET[,1])
  if(length(chr)!= nrow(PET)){
    dif<-setdiff(1:nrow(PET),chr);
    stop(paste(length(dif)," elements are not chromosome names :",head(dif)))
  }
  
  if(!is.numeric(PET[,2])){
    stop("The second column should be numeric");
  }
  
  if(!is.numeric(PET[,3])){
    stop("The Third column should be numeric");
  }
  
  nbLeft <-length( grep('PET#\\d+\\.1', as.character(PET[,4])) )
  nbRight <- length( grep('PET#\\d+\\.2', as.character(PET[,4])) )
  
  if(nbLeft != nbRight)
    stop("The the number of Left regions is diffrent from the right regions \n make sure that the regions are named correctly")
  
  PET.ranges <- GRanges(seqnames=Rle(PET[,1]),ranges = IRanges(start=PET[,2],end=PET[,3]),
                        PET_ID=as.character(PET[,4]))
  return(PET.ranges);
}




## Reads a BED file that contains the peak positions in the first three columns 
## and the name of the TF in the fourt, 
## here we suppose that all the TF files are combined into one file
.ReadTFBS<-function(file, header=FALSE){
  
  TFBS<-read.table(file,sep="\t", header= header);
  #Check if the file has 4 columns
  if(ncol(TFBS)!= 4){
    stop("Some columns are missing from the file")
  }
  
  #Check the type of each column
  chr<-grep("chr",TFBS[,1]);
  if(length(chr)!= nrow(TFBS)){
    dif<-setdiff(1:nrow(TFBS),chr);
    stop(paste(length(dif)," elements are not chromosome names :",head(dif)))
  }
  
  #Check if the second and third colum are numeric
  if(!is.numeric(TFBS[,2])){
    stop("The second column is not numeric");
  }
  
  if(!is.numeric(TFBS[,3])){
    stop("The third column is not numeric");
  }
  
  if(!is.factor(TFBS[,4])){
    warning("The fourth column should be a character that contains the TF name");
  }
  
  TFbed.ranges<-GRanges(seqnames=Rle(TFBS[,1]),ranges = IRanges(start=TFBS[,2],end=TFBS[,3]),TF=TFBS[,4])
  
  return(TFbed.ranges);
}

.getNuclearPPI<-function(ppiType=c("HPRD","Biogrid"), filter=FALSE,term ="GO:0005634", 
                         annot=NULL, RPKM= NULL, threshold=1){
  
  ppiType <- match.arg(ppiType)
  
  PPI<-c();
  #Load a PPI network
  if(ppiType == "HPRD"){
    data("HPRD", envir= environment());    
    PPI<-get("PPI.HPRD", envir=environment());
    message(paste("loading HPRD network, with", vcount(PPI),"nodes"))
  }
  else{
    data("Biogrid", envir= environment())
    PPI<-get("PPI.Biogrid",envir= environment());
    message(paste("loading Biogrid network, with", vcount(PPI),"nodes"))
  }
  
  
  if(filter){
    PPI <- .filterPPI(PPI, term, annot, RPKM, threshold)
  }
  
  return(PPI)
  
}

.filterPPI <-function(PPI,term="GO:0005634", annot=NULL, RPKM= NULL, threshold=1){
  
  if(term == "" & ( nrow(RPKM) ==0 || is.null(RPKM))) 
    stop("No filtering criteria was provided")
  
  
  ## if no annotation was provided we use the one provided by the package
  if(is.null(annot)){
    data("geneLocations",envir= environment());
    annot <- get("geneLocations.nucleus",envir= environment());
  }
  else{   
    if(class(annot) != "data.frame")
      stop("the annot parameter should be a data.frame object")
    
    if(! "cellular_component_term" %in% colnames(annot))
      stop("please make sure that annot contains a column named cellular_component_term")
  }
  
  ## Get the list of genes that have this GO term
  pos<-grep(term,annot$cellular_component_term)
  
  if(length(pos) == 0){
    warning(paste("no filtering done, could not find '", term,"' in the annotation table",sep=""))
    return(PPI)
  }
  
  inNucleus<-which(V(PPI)$name %in%  as.character(annot$geneSymbol[pos]))
  PPI_nucleus<-induced.subgraph(PPI, V(PPI)[inNucleus])
  ## Remove isolated genes
  isolated <- which(degree(PPI_nucleus) < 1)
  if(length(isolated) > 0){
    #warning(paste(length(isolated),"isolated genes with term ",term,"have been filtered \n"))
    #warning(paste("the isolated genes are :", V(PPI_nucleus)$name[isolated],sep="", collapse="\t"))
    PPI_nucleus<- delete.vertices(PPI_nucleus, isolated)
  }
  
  if(!is.null(RPKM)){
    PPI_nucleus <- .filterByExpression(PPI_nucleus,RPKM, threshold)
  }
  
  return(PPI_nucleus)
}


## Filters a given PPI given the RPKM expression of its genes
.filterByExpression<-function(PPI,RPKM, threshold=1){
  
  if(class(RPKM) != "data.frame"){
    stop("The RPKM paramter should be a data.frame object")
  }
  
  if(ncol(RPKM) <2)
    stop("the RPKM dataset should at least contain two columns")
  
  if(ncol(RPKM) == 2){
    colnames(RPKM) <- c("GeneSymbol", "RPKM")
  }
  else{
    if(!all( c("GeneSymbol","RPKM") %in% colnames(RPKM) ) )
      stop("Provide a two columns table or set your dataset to contain a column named GeneSymbol and RPKM")
  }
  
  pos<- match(V(PPI)$name,RPKM[,"GeneSymbol"])  
  if(length(pos) == 0){
    warning("None of the provided GeneSymbol corresponds to the PPI names, no filtering done") 
  }
  else{
    V(PPI)[!is.na(pos)]$RPKM<-RPKM[pos[!is.na(pos)],"RPKM"]  	
    PPI<-induced.subgraph(PPI,V(PPI)[which(V(PPI)$RPKM >= threshold)])
  }
  
  return(PPI)
}

## Calculates all the shorest paths connecting the provided nodes in the given PPI network
.calculate.all.shortest.paths<-function(PPI,nodes){  
  paths<-list();
  
  nodes<-nodes[which(nodes %in% V(PPI)$name)]
  #Calculate the distance between each node and the others and save in a list
  for(i in 1:length(nodes)){
    sourceName<-nodes[i];
    allPaths<- get.all.shortest.paths(PPI, V(PPI)[sourceName], V(PPI)[nodes]);
    #Get the destinations names
    desNames<-V(PPI)$name[sapply(allPaths$res,function(x) tail(x,1))];
    #Melt the paths with the same source together.
    names(allPaths$res)<-desNames;
    #Get the list or genes with multi-path
    idx<-which(duplicated(desNames)==FALSE);  	
    
    p<-sapply(1:length(idx),function(i) {
      res<-list();
      if(i==length(idx)) {
        if(idx[i]<length(allPaths$res)){
          res<-unique(unlist(allPaths$res[ idx[i]:length(allPaths$res)]))
          #work around: Just to get the source element at the end of the list
          res<-c(res,tail(allPaths$res[[idx[i]]],1)) 
        }
        else{
          res<-allPaths$res[[ idx[i] ]]
        }
      } 
      else { 
        if( idx[i+1]-idx[i] >1) {
          res<-unique(unlist(allPaths$res[ idx[i]:(idx[i+1]-1) ]))
          res<-c(res,tail(allPaths$res[[idx[i]]],1)) 
        } 
        else{ 
          res<-allPaths$res[[ idx[i] ]]
        }
      }
      
      return(res);
    }					
    )
    names(p)<-V(PPI)$name[sapply(p,function(x) tail(x,1))];
    paths[[sourceName]]<- p;
    #names(paths[[sourceName]])<-nodes;
  }
  
  return(paths)
}

## create a data.table that indexes the TF associated with each DNA region
.create.datatables_chip<-function(PETS, TFBS,minOverlap=50){
  
  
  #Create the PETs table
  print("Creating PET table")
  
  petNames <- as.character(PETS$PET_ID)
  
  petDT<-data.table(petID=1:length(petNames), 
                    petName= petNames, 
                    rootName=substr(petNames,1,nchar(petNames)-1)
  );
  
  #Because the table is large we put it as a reference using setkey
  print("Sorting Table");
  with(petDT, {setkey(petDT,rootName)}); 
  
  print("Get TF associated with each PETs")
  PET.TFs<-findOverlaps(PETS,TFBS,minoverlap=minOverlap)  
  
  print("Creating Motifs table")
  tfsName<-toupper(sort(unique(TFBS[subjectHits(PET.TFs),]@elementMetadata$TF)));	
  motifsDT<-data.table(motifID=1:length(tfsName),motifName= tfsName);
  print("Sorting Table");
  with(motifsDT, {setkey(motifsDT,motifID)});
  
  print("Creating hasMotif table")	
  #Estimate the size to prealocate memory so we gain in execution time
  pos<-match(toupper(TFBS[subjectHits(PET.TFs),]@elementMetadata$TF),tfsName)
  
  hasMotifDT<-data.table(petID=queryHits(PET.TFs), motifID=pos);	
  
  #put this table also as a reference to gain in access speed
  print("Sorting Table");
  with(hasMotifDT, {setkey(hasMotifDT,petID)});
  
  res<-list(PET=petDT, motifs= motifsDT, hasMotif=hasMotifDT);
  
  return(res);
  
}

## Builds network and returning a list of edges
.build.PETNetworks.byedge<-function(PET,Motifs,hasMotifs,PPI){
  
  NetsList<-list();
  
  pets<-unique(with(PET,{PET[,rootName]}));
  possibleMotifs<-unique(with(Motifs, {Motifs[,motifName]}));
  
  all.paths<-.calculate.all.shortest.paths(PPI,possibleMotifs);
  
  f<-function(p,paths,bgNetwork){  	
     rootName=petID=motifID=motifName=J=NULL
    #Get the two pairs that interact
    p.pair<-PET[rootName == p,petID];
    
    tf1<- Motifs[ J( hasMotifs[ J(p.pair[1])] [, motifID]  )] [ , motifName];
    tf2<- Motifs[ J( hasMotifs[ J(p.pair[2])] [, motifID]  )] [ , motifName];
    
    tf1<-unique(tf1);
    tf2<-unique(tf2);
    
    if(is.na(tf1) || is.na(tf2)){
      return(NULL);
    }
    
    #Get the list of the involved proteins
    net<-.get.network.byedge(paths,tf1,tf2, bgNetwork);
    return(net)
  }
  
  #res<-lapply(pets,f,paths=all.paths, bgNetwork=PPI)
  ## parallelize the computation so the users will not wait for a long time.
  ## works on all OS
  nbCores <- detectCores(logical=TRUE)
  message(paste(nbCores, "cores detected on your computer"))
  cl <- makeCluster(getOption("cl.cores", nbCores))
  message("")
  res <- parLapply(cl,pets, f,paths=all.paths, bgNetwork=PPI)
  stopCluster(cl)
  
  ##res<-mclapply(pets,f,paths=all.paths, bgNetwork=PPI, mc.cores= 1L, mc.preschedule=TRUE)
  names(res)<-pets;
  return(res);
}


## Calculates the global frequency of a TF interaction in the collection of networks
.get.TFsStat<-function(PetNet){
  TFs<-c()  
  
  for(i in 1:length(PetNet)){
    freq<-table(PetNet[[i]]);
    tfs<- names(freq);
    #if the TFs already exist increment their appearance frequency
    exist<-tfs[which(tfs %in% names(TFs))];
    TFs[exist]<-TFs[exist]+1;
    
    #For the elements who are not included yet, include them
    notExist<-setdiff(tfs,exist);
    newNames<-c(names(TFs),notExist);
    TFs<-c(TFs,rep(1,length(notExist)));
    names(TFs)<- newNames;		
  }
  
  return(TFs);
}

.get.network.byedge<-function(paths,source,dist,PPI){
  
  source<-source[which(source %in% names(paths))];
  
  
  if(length(source)==0){		
    return(NULL);
  }
  
  #lst<-c(source,dist);
  lst<-c();
  dist<-dist[ which(dist %in% V(PPI)$name) ]
  
  for(s in source){
    
    net<-unique(as.integer(unlist(paths[[s]][dist])));
    net<-net[!is.null(net)];
    
    if(length(net)>0){
      subNet<-simplify(induced.subgraph(PPI,V(PPI)[net]));
      res<-apply(get.edges(graph=subNet,E(subNet)),1,function(x){paste(sort(V(subNet)[x]$name),collapse="_")})
      lst<-c(lst , res);	
    }				
  }	
  
  return(lst);
}

.filter.singnifTF<-function(petNets,selectedTF){
    
    for(i in 1:length(petNets)){  
      if(length(petNets[[i]])>0){
        toUse<- which(petNets[[i]] %in% selectedTF);
        #petNets[[i]] <- setdiff(petNets[[i]],toRemove);		
        petNets[[i]]<- petNets[[i]][toUse];
      }
    }	
    return(petNets);
}

## Remove the edges that apear too much of seldomly
.RemoveOutliers<-function(petNets, minFreq=0.25,maxFreq=0.75){
  
  TF_freq.nucleus <- .get.TFsStat(PetNet= petNets)
  
  if(minFreq> 1 || maxFreq>1){
    stop("The pprovided probabilities should be between 0 and 1")
  }
  
  if(minFreq >= maxFreq){
    stop("The minimum frequency should be smaller than the maximum frequency")
  }
  
  prob<-c(0,minFreq,maxFreq,1)
  quant<- quantile(TF_freq.nucleus, probs= prob )
  TF_filtered.nucleus <- TF_freq.nucleus[which( TF_freq.nucleus >= quant[2] ) ]
  TF_filtered.nucleus <- TF_filtered.nucleus[which( TF_filtered.nucleus <= quant[3] ) ]  
  petNets.filtered <-.filter.singnifTF(petNets, names( TF_filtered.nucleus ) )
  
  return(petNets.filtered);
}

.get.TFSpace<-function(PetNet){
  
  TFs<-c();
  setSizes<-c();
  for(i in 1:length(PetNet)){
    TFs<-unique(c( TFs, PetNet[[i]]));
    setSizes<-c(setSizes,length(PetNet[[i]]));
  }
  
  res<-list(TF=TFs,size=setSizes);
  return(res);
}

.valid.buildNetworks <-function(object, minFreq, maxFreq){
  
  if(!is.numeric(minFreq) || !is.numeric(maxFreq))
    stop("minFreq and maxFreq should be numeric")
  
  ## check if all the slots are full
  if(length(object@pet) == 0)
    stop("no interaction data was loaded, check the method loadPET")
  
  if(length(object@pet) == 0)
    stop("no TF binding site locations were loaded, check the method loadTFBS")
  
  if(length(object@.dt) == 0)              
    stop("Data was not indexed, check the method createIndexes")
  else{
    if(! all(names(object@.dt) %in% c("PET","motifs","hasMotif")) )
      stop("The data.tables PET, motifs and hasMotif should be defined, check creatIndexes")
  }
}


###################################################################
##
##     Methods definition
##
######################################################################

setMethod("loadPETs",signature=c(object = "ChiapetExperimentData", petFile = "character"),
          function(object, petFile, IsBed = TRUE, header = TRUE,dist = 1000){            
            if(petFile==""){
              stop("no file was specified")
            }
            
            if(!file.exists(petFile)){
              stop(paste("could not find file ", petFile))
            }
            
            
            tryCatch({              
              ## if it is already a PET file we read it
              if(missing(IsBed) || (!missing(IsBed) & IsBed)){
                if(missing(header)) header= TRUE                
                object@pet <- .ReadPET(petFile,header)
              }
              else{
                if(!missing(dist) & dist <= 0){
                  stop("dist should be > 0")
                }                
                PET<- CreateCenteredBED(petFile, header, dist)              
                object@pet <- GRanges(seqnames=Rle(PET[,1]),
                                      ranges = IRanges(start=as.numeric(PET[,2]),end= as.numeric(PET[,3])),
                                      PET_ID=as.character(PET[,4]))
              }
              message(paste(length(object@pet)/2, "interacting DNA regions loaded"))              
              return(object)
              
            },error = function(err){
              stop(paste( "error when loading interactions \n",conditionMessage(err),sep=""))
            }
            )
            
          })

setMethod("pet", signature = c(object = "ChiapetExperimentData"), 
          function(object){
            return( object@pet)
          })
              
setReplaceMethod("pet", "ChiapetExperimentData",
          function(object, value){
            if(!"PET_ID" %in% names(mcols(object@pet)))
              stop("There is no metadata colum named PET_ID")
            object@pet <- value
            object
            })

setMethod("loadTFBS", signature=c(object="ChiapetExperimentData", tfbsFile="character"),
          def = function(object, tfbsFile,header = FALSE, ...){                        
            
            if(tfbsFile==""){
              stop("no file was specified")
            }
            
            if(!file.exists(tfbsFile)){
              stop(paste("could not find file ", tfbsFile))
            }
            
            tryCatch({
              object@tfbs <- .ReadTFBS(tfbsFile, header)
              message(paste("a total of ", length(object@tfbs), "binding sites for",
                            length(unique(object@tfbs$TF)),"TF  were loaded" ))
              return(object)
              
            },error = function(err){
              stop(paste( "error when loading TF binding sites \n",conditionMessage(err),sep=""))
            })
          })


setMethod("tfbs", signature = c(object="ChiapetExperimentData"),
          function(object){
            return(object@tfbs)
          })
              
setReplaceMethod("tfbs", "ChiapetExperimentData",
                 function(object, value){
                   if(!"TF" %in% names(mcols(object@tfbs)))
                     stop("There is no metadata colum named TF")
                   object@tfbs <- value
                   object
                 })
              
## loads the PPI, by default no filtering is done, how ever the user can give some custome
## filtering paramters
setMethod("loadPPI", signature=c(object="ChiapetExperimentData"),
          function(object,type=c("HPRD","Biogid"),customPPI= NULL, 
                   filter = FALSE, term ="GO:0005634", annot=NULL, RPKM= NULL, threshold=1 ){                        
            
            if(!is.null(customPPI)){
              if(class(customPPI)=="igraph"){
                object@ppi <- customPPI
              }
              else{
                if(class(customPPI) != "character")
                  stop("the customPPI parameter can be a path or an igraph object")
                
                PPI <- read.graph(customPPI, format="ncol")
                if(filter){
                  PPI <- .filterPPI(PPI, term, annot, RPKM, threshold)
                }
                object@ppi <- PPI
              }
            }
            else{
              object@ppi <- .getNuclearPPI(type, filter, term, annot, RPKM, threshold)
            }
            return(object)
          })

setMethod("ppi", signature= c(object = "ChiapetExperimentData"),
          function(object){
            return(object@ppi)
          })


setMethod("createIndexes", signature(object="ChiapetExperimentData"),
          def= function(object, minOverlap=50){            
            if(length(object@pet) == 0)
              stop("no interaction data was loaded, check the method loadPET")
            
            if(length(object@pet) == 0)
              stop("no TF binding site locations were loaded, check the method loadTFBS")
            PETS <- object@pet
            TFBS <- object@tfbs
            object@.dt <- .create.datatables_chip(PETS,TFBS, minOverlap)
            return(object)
          })


setMethod("buildNetworks", signature= c(object="ChiapetExperimentData"),
          function(object, minFreq = 0.25, maxFreq = 0.75){
            
            .valid.buildNetworks(object, minFreq, maxFreq)
            
            print(paste("building",nrow(object@.dt$PET)/2, "networks"))
            petNets<- .build.PETNetworks.byedge(PET=object@.dt$PET,Motifs=object@.dt$motifs,
                                                hasMotifs=object@.dt$hasMotif,object@ppi)
            
            print("Filtering networks")
            petNets.filtered<- .RemoveOutliers(petNets, minFreq, maxFreq)
            print("get the list of the involved interactions")
            tfSpace<- .get.TFSpace(PetNet=petNets.filtered)
            empty<-which(tfSpace$size==0)
            NonemptyNet<-which(tfSpace$size>0)
            petNets.filtered<-petNets.filtered[NonemptyNet]  
            
            netCollection <- new('NetworkCollection', 
                                 networks = petNets.filtered,
                                 sizes = tfSpace$size[NonemptyNet],
                                 TFCollection = tfSpace$TF)
            
          })

setMethod("visualizeInteractions", signature = c(object= "ChiapetExperimentData", range= "GRanges"),
          function(object, range){
            
			requireNamespace("biovizBase")
            interactions <- subsetByOverlaps(pet(object), range)          
            if(length(interactions) == 0)
              return("NA");
            ## to avoid the check notes
            data("hg19Ideogram", package = "biovizBase", envir=environment())            
            hg19Ideogram <- get("hg19Ideogram", envir=environment());
            hg19Ideo <- keepSeqlevels(hg19Ideogram, seqlevels(pet(object)))
            seqlengths(object@pet) <- seqlengths(hg19Ideo)
                
            
            roots <- unique(gsub("\\.\\d+","",interactions$PET_ID))
            
            ## get left and right side interactions
            leftID  <- match(paste(roots,".1",sep=""), pet(object)$PET_ID)
            rightID <- match(paste(roots,".2",sep=""), pet(object)$PET_ID)
            
                        
            circos <- pet(object)[leftID]
            values(circos)$to.gr <- pet(object)[rightID]
                        
            p <- ggplot() + layout_circle(hg19Ideo, geom = "ideo", fill = "#9ecae1", color="#636363", radius = 30,trackWidth = 4)
            p <- p + layout_circle(hg19Ideo, geom = "text", aes(label = seqnames), 
                                   vjust = 0,radius = 32, trackWidth = 7)
            
            p <- p + layout_circle(circos, geom = "link", linked.to = "to.gr",radius = 29, trackWidth = 1, color="#f03b20")            
            p <- p + ggtitle("Interactions in the selected range") + 
              theme(plot.title = element_text(lineheight=.8, face="bold"))
            plot(p)
            invisible(list(circos = circos,plot = p))
})


setMethod("plotTrack", signature = c(object= "ChiapetExperimentData", range= "GRanges"),
          function(object, range){
			
			requireNamespace("biovizBase")
            interactions <- subsetByOverlaps(pet(object), range)
            
            if(length(interactions) == 0)
              return("NA");
            
            data("hg19Ideogram", package = "biovizBase",envir=environment())    
            hg19Ideogram <- get("hg19Ideogram", envir=environment())
            data("hg19IdeogramCyto", package = "biovizBase", envir=environment())
            hg19IdeogramCyto <- get("hg19IdeogramCyto", envir=environment())
                              
            roots <- unique(gsub("\\.\\d+","",interactions$PET_ID))
            
            ## get left and right side interactions
            leftID  <- match(paste(roots,".1",sep=""), pet(object)$PET_ID)
            rightID <- match(paste(roots,".2",sep=""), pet(object)$PET_ID)
            
            intra <- which( as.character(seqnames(pet(object)[leftID])) == as.character( seqnames(pet(object)[rightID]) ) )
            
            leftPos <- leftID[intra]
            rightPos <- rightID[intra]
            
            interactions <- GRanges(as.character(seqnames(pet(object)[leftPos])), 
                                    IRanges( start(pet(object)[leftPos]), start(pet(object)[rightPos]) )                                    
                                   )
            
            hg19Ideo <- keepSeqlevels(hg19Ideogram, seqlevels(interactions))
            seqlengths(interactions) <- seqlengths(hg19Ideo)
            
            
            p1 <- autoplot(interactions, geom="arch", color="#67a9cf",size=1) +theme_bw()
            p2<-plotIdeogram(hg19IdeogramCyto, as.character(seqnames(interactions)[1]),
                             zoom.region = c(min(start(interactions)),
                                             max(end(interactions))))
            
            tks <- tracks(p2, p1, heights=c(0.2,1)) + xlim(start(range),end(range))
            tks
            invisible(tks)
          })



setMethod("show",'ChiapetExperimentData', 
          function(object){
            cat("class:", class(object),"\n")
            if(!is.null(object@pet))
              cat(length(object@pet)/2, " interacting regions","\n")
            else
              cat("no interaction data have been loaded yet\n")
            
            if(!is.null(object@tfbs))
              cat(length(unique(object@tfbs@elementMetadata$TF)), "TF used","\n")
            else
              cat("no TF binding data have been loaded yet\n")
            
            if(!is.null(object@ppi))
              cat(paste("Background PPI: \n nodes:",vcount(object@ppi),"\t edges:", ecount(object@ppi), "\n"))
            else
              cat("no PPI has been loaded yet\n")
            
            if(!is.null(object@.dt) && length(object@.dt) > 0)
              cat("indexes tables have been created \n")
            else
              cat("No indexes created yet \n")
          })
                
                
                
## A user-frindly S3 method for the construsctor
ChiapetExperimentData <- function(pet='', tfbs='', ppi=NULL, 
                                  ## loadPETs options
                                  IsBed=TRUE, petHasHeader=FALSE, dist=1000,
                                  ## loadTFBS options
                                  tfbsHasHeader=FALSE,
                                  ## loadPPI options
                                  ppiType=c("HPRD","Biogid"),
                                  filter=FALSE, term="GO:0005634", annot=NULL,
                                  RPKM= NULL, threshold=1
                                  ){
   
  obj <- new("ChiapetExperimentData")
  obj <- .check.pet(pet, obj, IsBed, petHasHeader, dist)
  obj <- .check.tfbs(tfbs,obj, tfbsHasHeader)
  obj <- loadPPI(obj,type =  ppiType, customPPI = ppi, filter = filter,
                 term = term, annot = annot, RPKM = RPKM,
                 threshold = threshold)  
  
  return(obj)
}
