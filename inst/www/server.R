library(shiny)
library(R3CPET)
library(ggplot2)
library(rCharts)

load("data/x.RData")
load("data/hlda.RData")
load("data/nets.RData")

constructNetwork <- function(func){
    reactive({
        res <- func()        
        i <- res$idx
        hlda <- res$hlda
        if(i > length(networks(hlda)))
            stop("The provided network index is not correct")
            
        nodes <- V(networks(hlda)[[i]])$name
        
        edgeList <-get.edgelist(hlda@networks[[i]])
        
        links <- t(apply(edgeList,1, function(x) {c(match(x[1], nodes)-1, match(x[2], nodes)-1, 1) }))        
        colnames(links) <- c("source","target", "weight")
        data <- list(names=nodes, links= links)
        return(data)
    })
}

df <- as.data.frame(pet(x))


shinyServer(function(input, output, session) {

    print("Hellow")
    left <- df[grep("PET#\\d+\\.1",df$PET_ID),]
    right <- df[grep("PET#\\d+\\.2",df$PET_ID),]    
    observe({      
        ## #################################################################
        ##          Raw data display
        ## #################################################################
        
        if(input$DataType == 'RawTab'){
        if(input$Raw == "PET"){
            
             if(input$petStat == 'petInfo'){
              observe({
                isolate(
                 output$distPlot <- renderChart({   
                 y <-as.data.frame(table(df$seqnames))                                                                  				 
                 p <- rPlot(x = "Var1", y = "Freq", data = y, type = "bar", color="Var1")     
                     
                p$guides(x = list(title = "Chromosome"))
                p$guides(y = list(title = "Number of interactions"))
                        
                p$addParams(dom = 'distPlot',
                 title = "Number of interactions per chromosome")                 
                 return(p)
                 }) 
                ) ## isolate
                }) ## observe
                 
                 observe({
                isolate( output$petInfoTab <- renderTable({
                    inter <- data.frame(left= left$seqnames, right = right$seqnames)

                    y<- data.frame("Total Interactions" = c(nrow(left)),
                                   "Intra-chromosome interaction" = nrow(subset(inter, left == right)),
                                   "Inter-chromosome interactions" = nrow(subset(inter, left != right))
                                  )
                    y
                    })
                  ) ## isolate
                }) ## observe
             }
            else if(input$petStat == 'petHeatmap'){
                    observe({
                    isolate( output$heatPlot <- renderChart({

                        inter <- data.frame(table(data.frame(left= left$seqnames, right = right$seqnames)))
                        p <- rPlot(x = 'left', y = 'right', color = 'Freq', data = inter, type = 'tile')
                        #p$guides(x = list(title = "Left"), ticks = unique(inter$left))
                        #p$guides(y = list(title = "Right"), ticks = unique(inter$right))                        
                        p$addParams(dom = 'heatPlot',title = "Interaction Frequency heatmap")                 
                        return(p)                      

                    })
                    ) ## isolate
                    }) ## observe
                }                
                
                else if(input$petStat == 'span'){

                    observe({
                        isolate({ span <- abs(left$start - right$start);                
                                output$slider <- renderUI({sliderInput("petSpanSlider", "span", 
                                                                       min= min(span), max= max(span), 
                                                                       value= floor( 0.5 * (max(span)+ min(span))), ,step=500)  
                                                          })                
                                })

                       if(! is.null(input$petSpanSlider)){                
                        output$petSpanPlot <- renderChart({ 
                            dense = density(span[span <=  input$petSpanSlider ])
                            dense = data.frame(dense$x, dense$y)
                            p <- hPlot(x = "dense.x", y = "dense.y", data = dense, type = "line")
                            p$params$xAxis[[1]]$title$text <- "Span"
                            p$params$yAxis[[1]]$title$text <- "Density"

                            p$params$title$text <- "Density plot of the interactions span";
                            p$params$dom <- 'petSpanPlot';    
                            return(p)
                            #plot(density(span[span <=  input$petSpanSlider ]), xlab="span",main="Interaction Span distribution") 
                        })
                       }
                  })
                }
              else if(input$petStat == 'circos'){
                  observe({
                      input$btnDispCircos;
                      plot(input$btnDispCircos)
                      isolate({
                      if(input$circosChrPos != ""){
                          pos <- unlist(strsplit(input$circosChrPos, split=":"))          
                          if(length(pos) ==2 && grep("chr",pos[1]) ==1){
                              chr <- pos[1];
                              pos <- unlist(strsplit(pos[2],split="-"))                            
                              gr <- GRanges(chr, IRanges(as.numeric(pos[1]),as.numeric(pos[2])))                              
                              print("plotting Circos")
                              output$circosPlot<- renderPlot({  visualizeInteractions(x,gr) })              
                              # TDDO: try to fix this
                              #print("plotting Track")                              
                              #output$trackPlot <- renderPlot({ trks <-plotTrack(x,gr)+ xlim(start(gr),end(gr));                                                                                                                              
                              #                                 print(trks@xlim)
                              #                                 print(trks)})
                          }
                      }
                  })
                  })
             }
        } ## input$Raw == "PET"
        else{
            if(input$Raw == "TFBS"){                
                output$tfbsStatTable <- renderDataTable({
                    dt <- x@.dt$hasMotif[,length(unique(petID)),by=motifID ]
                    setnames(dt, "V1", "nb_Interacting_Regions");                                        
                    motifNames <- x@.dt$motifs[motifID %in% dt[,motifID], motifID, motifName]                                        
                    dt <- merge(motifNames,dt, by="motifID")
                    freq <- as.data.frame(table(tfbs(x)$TF)) 
                    names(freq) <- c("motifName", "Total_Binding_Site")
                    dt <- merge(dt,freq, by="motifName")
                    dt
                })
                
                
                output$tfbsPerRegionPlot <- renderChart({
                    dt <- x@.dt$hasMotif[,length(unique(motifID)),by=petID ]
                    dt <- as.data.frame(dt)                    
                    p <- rPlot(x = "bin(V1,1)", y = "count(petID)", data = dt, type = "bar")
                    #hist(dt[,V1],freq=T,xlab="Number of TFs per regions",
                    #     main="Histogram of the number of TFs per interacting region")
                    p$addParams(dom = 'tfbsPerRegionPlot')
                    return(p)
                })
            }
            else{
                if(input$Raw == "PPI"){
                    
                    output$ppiTableOutput <- renderTable({
                        ppiInfo <- data.frame("Nodes" = vcount(ppi(x)) , "Edges" = ecount(ppi(x)),
                                              "diameter" = diameter(ppi(x)), "density" = graph.density(ppi(x))                                              
                                             ) 
                        ppiInfo
                    })
                    
                    output$ppiDegreeDist <- renderChart({
                        y <- degree.distribution(ppi(x))
                        dd <- data.frame(degree= 0:(length(y)-1), frequency = y)
                        p <- rPlot(frequency ~ degree, data=dd, type="point")
                        p$guides(x = list(max = length(y) + 5, title = ""))
                        p$addParams(dom = 'ppiDegreeDist')
                        
                        return(p)
                    })
                    
                }
            }
        }
      }
      ## #################################################################
      ##          Results display
      ## #################################################################
      else if(input$DataType == "ResultsTab" ){          
        if(input$Results == "Nets"){                        
		    observe({
                print('Updating the list of clusters')
                clusList <- list();
				for(i in sort(unique(getClusters(hlda))))
					clusList[[sprintf("Cluster %d",i)]] <- i				
                    
				updateSelectInput(session,"ClusListInput",choices = clusList, selected=1);	
                
				print("Updating the list of networks")
				netList <- list();
				for(i in 1:ncol(topNodes(hlda)))
					netList[[sprintf("Network %d",i)]] <- i				
                    
				updateSelectInput(session,"Maintainers",choices = netList, selected=1);							                                            						                
			})			
			if(input$netStat == "allNetSizeDist"){
				observe({								
					#isolate({
						
					output$MainterNetsDist <- renderChart({ 
                        print("Rendring chart")
                        hlda <-  updateResults(object=hlda, nets, input$cutThresholdInput)
						 netSizes <- c()
						 for(i in 1:length(networks(hlda))){
							#netSizes <- c(netSizes, rep(i,vcount(networks(hlda)[[i]])))
							netSizes <- c(netSizes, vcount(networks(hlda)[[i]]))
						 }
						 
						 sdf <-data.frame(network = factor(paste("Net",1:length(networks(hlda)), sep="")), Freq= netSizes) 						 
						 p <- rPlot(x = "network", y = "Freq", data =sdf, type = "bar", color="network")     
						
						p$guides(x = list(title = "Chromatin maintainer network"))
						p$guides(y = list(title = "Number of nodes"))
								
						p$addParams(dom = 'MainterNetsDist',
									title = "Chromatin maintainers networks sizes")                 
						 return(p)
					#}) 
					})
                }) ## observe				
			}
			else{
				if(input$netStat == "hldaNets"){                       
					if(input$rbtnNets == "cmnList"){
						output$hldaNetTableOutput <- renderDataTable({     
                            print("rendring table")
                            hlda <-  updateResults(object=hlda, nets, input$cutThresholdInput)
							return(topNodes(hlda))
						})
					}
					else                
					if(input$rbtnNets == "visNets"){
					
					output$hldaNetOutput <- constructNetwork(function(){
                        print("rendring netwrok")
                        hlda <-  updateResults(object=hlda, nets, input$cutThresholdInput)
                        input$cutThresholdInput                                                
						input$rbtnNets						
						i <- as.numeric(input$Maintainers)
						return(list(idx=i, hlda=hlda));
						})
					}
				}
			}
        }
        else{
            if(input$Results == "clusters"){				
                if(input$clusterInfo == "clusHeatmap"){
					output$clusHeatmapOutput<- renderPlot({  plot3CPETRes(hlda,type="heatmap") }) 
				}
				else{					
					if(input$clusterInfo == "clusCurves"){
						output$clusCurvesOutput<- renderPlot({ plot3CPETRes(hlda,type="avgCurve") }) 
					}
					else{
						if(input$clusterInfo == "genePerCluster"){
                            output$genesPerClusStatTable <- renderDataTable({
							     
						  })
						}
                        else{
                            if(input$clusterInfo == "clusCircos"){
                                output$clusCircosOutput <- renderPlot({                                    
                                    visualizeCircos(hlda,x, cluster=as.numeric(input$ClusListInput))                                    
                                })
                            }
                        }
					}
				}				
            }
        }
      }
                
    })  
})