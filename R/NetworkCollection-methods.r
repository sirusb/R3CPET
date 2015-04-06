## TODO: check if the names of the edges are all of the form A_B
## need to be done in a fast way (maybe use Rcpp)
.valid.NetworkCollection.networks <- function(x){
    
  if(!is.null(x@networks) && !all(sapply(x@networks, is.character)) )
    return(paste("netowrks should be a list of edges objects"))
  else
    return(TRUE)
}

.valid.NetworkCollection.sizes <- function(x){
    
  if(!is.null(x@sizes)){
    if(length(x@sizes) != length(x@networks)){
      return(paste("The sizes and networks properties should have the same length"))
    }
    else{    
      sameSize= TRUE
      ##check if the sizes are correct
      for(i in 1:length(x@sizes)){
        if(x@sizes[i] != length( x@networks[[i]] ) ){
          sameSize=FALSE
          break;
        }
      }
      
      if(!sameSize)
        return(paste("The size of network",i, "was is not correct"))
      else
        return(TRUE)
    }
  }
  else
    return(TRUE)
}

## we don't check on this on the creation time, as it can be time consuming
.valid.NetworkCollection.TFCollection <- function(x){
  return(TRUE)
}

.valid.NetworkCollection <- function(x){  
  msg =.valid.NetworkCollection.networks(x)
   if(!is.logical(msg) || !msg)
     return(msg)
  else{
    msg = .valid.NetworkCollection.sizes(x)
    if(!is.logical(msg) || !msg)
      return(msg)
    else
      return(.valid.NetworkCollection.TFCollection(x))      
  }    
}

## get the set of enriched edges per network
.print.topwords <- function(words, vocab, thr=0.5) 
{
  
  num.topics <- nrow(words) 
  topics <- NULL 
  for (k in seq(num.topics))
  {
    prob <- words[k,]
    total <- sum(prob)
    prob <- prob/total
    s <- sort.int(x=prob, decreasing=TRUE, index.return=TRUE)
    psum<-0;
    pos<-1;
    repeat{  					
      psum<- psum +  s$x[pos];
      if(psum >= thr) {break;}			
      pos<- pos+1;			
    }
    top.idx <- s$ix[1:pos]
    topic.prob <- prob[top.idx]
    topic.words <- as.character( vocab[top.idx]);
    if(!is.null(topics) && nrow(topics) < length(topic.words)){
      diff<- length(topic.words) - nrow(topics);
      toAdd<-matrix("",ncol=ncol(topics),nrow=diff);
      topics<- rbind(topics,toAdd);
    }
    else{
      if(!is.null(topics) && nrow(topics) > length(topic.words)){
        
        diff<-nrow(topics) - length(topic.words);
        topic.words<-c(topic.words,rep("",diff))
      }			
    }
    
    
    topics <- cbind(topics, topic.words) 
    #head <- paste(head, sprintf("%50d", k), sep="")
  }
  colnames(topics)<-paste("Topic_",1:ncol(topics),sep="");   
  return(topics);
}


## covert the list of edges into a list of nodes
.getNodesList<-function(topedges){
  
  topwords<- NULL;
  
  for(i in 1:ncol(topedges)){
    lst<-unique(unlist(strsplit(topedges[,i],split="_")))
    
    #if there are some elements, check the size
    if(! is.null(topwords)){
      #if less, add some empty elements
      if(length(lst) <= nrow(topwords)){
        dif<-nrow(topwords)-length(lst);
        lst<-c(lst,rep("",dif));
      }
      else{
        dif<-length(lst)-nrow(topwords);
        toAdd<-matrix("",ncol=ncol(topwords),nrow=dif);
        topwords<- rbind(topwords,toAdd);
      }
    }  		
    topwords<-cbind(topwords,lst);
  }
  
  colnames(topwords)<-colnames(topedges);
  return(topwords);
}


######################################################################################
##
##      S4 Methods
##
#####################################################################################

setMethod('InferNetworks',signature= c(object="NetworkCollection"),
          function(object,thr =0.5,max_iter = 500L, max_time = 3600L, ...){
            
            petNets <-networks(object)
            HDA_DATA<- .ConvertToHDA(petNets ,tfspace=TF(object))
            print("Estimating the number of topics");
            HLDA_res<- RunHLDA(HDA_DATA, max_iter = max_iter, max_time = max_time);
            
            #Normalize the frequencies
            HLDA_res$topicPerDoc <-t(apply(HLDA_res$topicPerDoc ,1,function(x){x/sum(x)}))
            colnames(HLDA_res$topicPerDoc)<- paste("Topic", 1:ncol(HLDA_res$topicPerDoc),sep="")
            rownames(HLDA_res$topicPerDoc)<-gsub("[\\.|PET#]","",names(petNets))                                                    
            Res<- HLDAResult(HLDA_res$topicPerDoc, HLDA_res$wordsPerTopic, HLDA_res$Betas)
            
            HLDA_res$topics_edge<- .print.topwords(HLDA_res$wordsPerTopic,as.matrix(TF(object)),thr)
            HLDA_res$topics_TF <- .getNodesList(HLDA_res$topics_edge)          
            
            return(ChromMaintainers(Res, HLDA_res$topics_edge, 
                                    HLDA_res$topics_TF, structure(list(),class= "cluesOrSota")
            )
            )
          })

setMethod("show","NetworkCollection",
          function(object){
            cat("class", class(object),"\n")
            cat(length(object@networks), "network loaded \n")
            cat(length(object@TFCollection), "diffrent edges has been used \n")
          })


setMethod("networks","NetworkCollection",
          function(object){
            return(object@networks)
          })

setMethod("sizes","NetworkCollection",
          function(object){
            return(object@sizes)
          })

setMethod("TF","NetworkCollection",
          function(object){
            return(object@TFCollection)
          })

setValidity("NetworkCollection", .valid.NetworkCollection)

NetworkCollection <- function(networks, sizes,  TFCollection){
  return(new('NetworkCollection', 
      networks = networks,
      sizes = sizes,
      TFCollection = TFCollection)
  )
}



