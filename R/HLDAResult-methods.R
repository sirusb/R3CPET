.valid.HLDAResult <-function(x){
  
  if(!is.matrix(x@docPerTopic))  
    return("docPerTopic should be a matrix objet")
  
  if(!is.matrix(x@wordsPerTopic))
    return("wordsPerTopic should be a matrix object")
  
  if(ncol(x@docPerTopic) != nrow(x@wordsPerTopic))
    return("The number of columns in docPerTopic should be equal to the number of rows in wordsPerTopic")
  
  if(!is.numeric(x@betas) || length(x@betas) != ncol(x@docPerTopic))
    return("The beta vector should have the same length as the number of topics")
  
  return(TRUE)
}



###############################################################################
##
##              S4 Methods
##
###############################################################################

setValidity("HLDAResult",.valid.HLDAResult)

## a user friendly S3 function
HLDAResult<-function(docPerTopic,wordsPerTopic ,betas){  
  
  return(new("HLDAResult",docPerTopic = docPerTopic, wordsPerTopic = wordsPerTopic ,betas = betas))
}

setMethod("wordsPerTopic","HLDAResult",
           function(object) {
             return(slot(object,"wordsPerTopic"))
           })

setMethod("docPerTopic", "HLDAResult",
          function(object){
            return( slot(object,"docPerTopic"))
          })

setMethod("betas", "HLDAResult",
          function(object){
            return( slot(object, "betas"))
          })

setMethod("show", "HLDAResult", 
          function(object){
            if(!is.null(slot(object,"docPerTopic"))){
              cat(nrow(slot(object,"docPerTopic")), "element have been classified into",ncol(slot(object,"docPerTopic")), "topics \n")
              cat("Number of different words", ncol(slot(object,"wordsPerTopic")), "\n")
            }
            else
              cat("No Results yet\n")
})
