.preparareData<- function(petFile, tfbsFile,petIsBed){

  print(paste("Reading PET file:",petFile))
  PETs<-.ReadPET(petFile);
  K562_BED<-.ReadTFBS(tfbsFile);
  DTs<-.create.datatables_chip(PETs,K562_BED)  
}


setMethod("PrepareData", signature(petFile="character",tfbsFile="character", petIsBed="logical"), 
          def= function(petFile,tfbsFile, petIsBed = TRUE)
          {
            
            if(petFile == ""){
              stop("no interaction file was specified in the variable petFile")
            }
            if(!file.exists(petFile)){
              stop(paste("could not open file", petFile))
            }
            
            if(tfbsFile==""){
              stop("no TF binding sites were specified in the variable tfbsFile")
            }
            
            if(!file.exists(tfbsFile)){
              stop(paste("could not open file", tfbsFile))
            }
            
            tryCatch({
              .preparareData(petFile, tfbsFile, petIsBed)
            },error= function(err){
              stop(conditionMessage(err))
            })
          })