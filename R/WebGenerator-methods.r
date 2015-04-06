setMethod("createServer", signature = c(x = "ChiapetExperimentData", nets="NetworkCollection",
                                        hlda="ChromMaintainers"),
          function(x,nets,hlda){
            requireNamespace("shiny")
            serverDir <- system.file("www", package="R3CPET")
            tmpDir <- file.path(tempdir(),"3CPET_Shiny")
            
            ## Get the path the server files
            rScripts <- list.files(serverDir,"*\\.R", full.names=TRUE)
            jsScripts <- list.files(serverDir,"*\\.js", full.names=TRUE)
            ## create the temporary folder to hold the server
            dir.create(tmpDir,showWarnings= FALSE)
            
            ## Copy the scripts
            file.copy(rScripts, tmpDir)
            file.copy(jsScripts, tmpDir)
            
            ## Create the data directory
            dataFolder <- file.path(tmpDir,"data")
            dir.create(dataFolder, showWarnings= FALSE)
                        
            ## Save the objects there
            save(x,file= file.path(dataFolder, "x.RData"))
            save(hlda,file= file.path(dataFolder, "hlda.RData"))
            save(nets,file= file.path(dataFolder,"nets.RData"))
                        
            
            ## Run the application
            runApp(tmpDir)
          })
