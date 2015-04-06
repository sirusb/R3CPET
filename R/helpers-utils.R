#This function used in preprocessing, it will split each interaction into two lines
#one contains information about the head and the other about the tail,
#each PET will be named

.checkChIAPETFormat<- function(file, header){
  
  if(!file.exists(file)){
    stop(paste("couldn't find file", file))
  }
  
  if(!is.logical(header)){
    stop("the varaible header should be logical")
  }
  
  interactionFile <- read.table(file,sep="\t",header=TRUE);
  
  if(nrow(interactionFile) < 6){
    stop("At least 6 columns file should be provided \n Check help for more info")
  }
  
  if(length(grep("chr",interactionFile[,1])) == 0 || 
       (!is.factor(interactionFile[,1]) & 
          !is.character(interactionFile[,1]))){
    stop("The first column should specify the left chromatin name (exp: chr1)")
  }
  
  if(length(grep("chr",interactionFile[,4])) == 0 || 
       (!is.factor(interactionFile[,4]) & 
          !is.character(interactionFile[,4]))){
    stop("The first column should specify the right chromatin name (exp: chr1)")
  }
  
  if(!is.numeric(interactionFile[,2]) || !is.numeric(interactionFile[,3])){
    stop("The start end stop coodrinates of the left region should be numeric")
  }
  
  if(!is.numeric(interactionFile[,5]) || !is.numeric(interactionFile[,6])){
    stop("The start end stop coodrinates of the right region should be numeric")
  }
  
  ## if every thing is ok, return the interaction table after renaming it
  colnames(interactionFile)<- c("chromleft","startleft","endleft",
                                "chromright","startright","endright")
  return(interactionFile)
}

.createBEDFile.Centered<-function(file, header=TRUE,dist=1000){
    
  if(!is.numeric(dist)){
    stop("dist should be numeric")
  }
  
  ## Read the file    
  chiapetFile<- .checkChIAPETFormat(file, header)
  
  ## First remove all interactions with chrM
  
  pos<-which(as.character(chiapetFile$chromleft) == "chrM");
  pos<-c(pos,which(as.character(chiapetFile$chromright) == "chrM"));
  
  pos<-unique(pos);
  
  if(length(pos)>0){
    chiapetFile<- chiapetFile[setdiff(1:nrow(chiapetFile),pos),];
  }
  
  #Create the letft pets BED file
  petName<-paste("PET#",1:nrow(chiapetFile),".1",sep="");
  nbElem<-2*nrow(chiapetFile)
  allpets<-data.frame(chr=rep("",nbElem),start=c(0,nbElem),end=c(0,nbElem),
                      ID=rep("",nbElem))    
  
  ## The data.frame factor problem :(
  allpets$chr = as.character(allpets$chr)
  allpets$ID  = as.character(allpets$ID)
                      
  PetCenter<- floor((chiapetFile$startleft+ chiapetFile$endleft)/2);
  
  firstPetIndex<-seq(1,nrow(allpets),2) #Put it in the odd rows
  secondPetIndex<-seq(2,nrow(allpets),2) #Put in the even rows
  
  allpets[firstPetIndex,1]<-as.character(chiapetFile$chromleft);
  allpets[firstPetIndex,2]<-PetCenter-dist;
  allpets[firstPetIndex,2][as.numeric(allpets[firstPetIndex,2]) < 0] <-0;
  allpets[firstPetIndex,3]<-PetCenter+dist;
  allpets[firstPetIndex,4]<-petName;
  
  #Add the right side pets
  
  petName<-paste("PET#",1:nrow(chiapetFile),".2",sep="");
  
  PetCenter<- floor((chiapetFile$startright + chiapetFile$endright)/2);
  
  allpets[secondPetIndex,1]<-as.character(chiapetFile$chromright);
  allpets[secondPetIndex,2]<-PetCenter-dist;
  allpets[secondPetIndex,2][ as.numeric( allpets[secondPetIndex,2] ) < 0] <-0;
  allpets[secondPetIndex,3]<-PetCenter + dist;
  allpets[secondPetIndex,4]<-petName;
   
  return(allpets)
}


setMethod("CreateCenteredBED","character",
          def=.createBEDFile.Centered
)
