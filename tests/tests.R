test_create_ChiapetExperimentData <- function(){
    x <- ChiapetExperimentData()
    checkTrue(class(x) == "ChiapetExperimentData", 
                "No problem creating ChiapetExperimentData ")
}

test_interactions_file <- function(){
    petFile <- file.path(system.file("example",package="R3CPET"),
                         "HepG2_interactions.txt")
    
    chechTrue(file.exists(petFile))
}

test_TFBS_file <- function(){
    tfFile <- file.path(system.file("example",package="R3CPET"),
                        "HepG2_TF.txt.gz")
    chechTrue(file.exists(tfFile))
}

test_loadPETS <- function(){
    x <- ChiapetExperimentData()
    checkEquals(class(x),"ChiapetExperimentData")
    petFile <- file.path(system.file("example",package="R3CPET"),
                         "HepG2_interactions.txt")
    
    test_interactions_file()
    
    x <- loadPETs(x,petFile=petFile, IsBed=FALSE)
    
    checkTrue(length(pet(x)) >0, "PETs can be loadded")
}

test_loadPETS <- function(){
    x <- ChiapetExperimentData()
    checkEquals(class(x),"ChiapetExperimentData")
    tfFile <- file.path(system.file("example",package="R3CPET"),
                        "HepG2_TF.txt.gz")
    
    test_TFBS_file()
    
    x <- loadTFBS(x,tfbsFile= tfFile)
    
    checkTrue(length(tfbs(x)) >0, "TFBS can be loadded")
}


test_createIndex <- function(){
    x <- ChiapetExperimentData()
    
    
    tfFile <- file.path(system.file("example",package="R3CPET"),
                        "HepG2_TF.txt.gz")
    x <- loadTFBS(x,tfbsFile= tfFile)
    
    petFile <- file.path(system.file("example",package="R3CPET"),
                         "HepG2_interactions.txt")
    x <- loadPETs(x,petFile=petFile, IsBed=FALSE)
    
    x<- createIndexes(x)
    
    checkEquals(length(x@.dt),3)
    checkIdentical(names(x@.dt), c("PET","motifs", "hasMotif"))
    
    for(i in 1:3) checkTrue("data.table" %in% class(x@.dt[[i]]) )
}