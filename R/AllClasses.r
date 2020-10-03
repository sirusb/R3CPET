setOldClass("igraph")

setClass("ChiapetExperimentData", 
         representation(
           pet="GRanges",
           tfbs="GRanges",
           ppi="igraph",
           .dt="list"),
         prototype = prototype(pet=GRanges(), tfbs = GRanges(),
                          ppi= NULL, .dt =list())
         )

setClass("NetworkCollection",         
         representation(
            networks = "list",
            sizes = "numeric",
            TFCollection = "character"
           ),
         prototype= prototype(networks=NULL, sizes = numeric(), TFCollection = character())
        )



setClass("HLDAResult",                  
         representation(
           docPerTopic = "matrix",
           wordsPerTopic = "matrix",
           betas = "numeric"           
           ),
         prototype= prototype(docPerTopic=NULL, wordsPerTopic = NULL, betas= numeric())
         )

setOldClass("sota")
setClassUnion("cluesOrSota", c("sota","NULL"))

setClass("ChromMaintainers",         
        representation(
          maintainers = "HLDAResult",
           topEdges = "matrix",
           topNodes = "matrix",
           networks = "list",
           clusRes = "cluesOrSota"          
           )
        )
