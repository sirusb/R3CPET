library(shiny)
library(rCharts)

## https://gist.github.com/xiaodaigh/7012930
actionButton <- function(inputId, label, btn.style = "primary" , css.class = "") {
  if ( btn.style %in% c("primary","info","success","warning","danger","inverse","link")) {
    btn.css.class <- paste("btn",btn.style,sep="-")
  } else btn.css.class = ""
  
  tags$button(id=inputId, type="button", class=paste("btn action-button",btn.css.class,css.class,collapse=" "), label)
}
    
reactiveNetwork <- function (outputId) 
{
  HTML(paste("<div id=\"", outputId, "\" class=\"shiny-network-output\"><svg /></div>", sep=""))
}
    
# Define UI for application that plots random distributions 
shinyUI(pageWithSidebar(

  # Application title
  headerPanel("Visualize 3CPET results"),

  # Sidebar with a slider input for number of observations
  sidebarPanel( 
      tabsetPanel(
            tabPanel("Raw Data",
                wellPanel(
                    h6("Select row data to display."),
                    selectInput("Raw", "",
                               c("Interactions" = "PET",
                                "TFBS"= "TFBS",
                                "PPI"= "PPI"))
                    ,                    
                    ## Panel to manipulate interactions                
                    conditionalPanel(
                        condition = "input.Raw == 'PET'",                                                
                        h6("Statistics:"),
                        selectInput("petStat", "",
                                   c("General info" = "petInfo",
                                     "Heatmap" = "petHeatmap",
                                     "Interaction Span" = "span",
                                     "circos map" = "circos")
                                    ),                        
                        
                        ## Panel to display circos map
                        conditionalPanel(
                            condition = "input.petStat == 'circos'",
                            textInput("circosChrPos","Position:",
                                     "chr1:1000000-15000000"),
                            actionButton("btnDispCircos","generate circos")
                           )                        
                        ),
                        
                     ## ToDO: Panel to manipulate TFBS
                     conditionalPanel(
                         condition = "input.Raw == 'TFBS'",                         
                              h6("Statistics:"),
                              selectInput("tfbsStat", "",
                                   c("General info" = "tfbsInfo",
                                     "AdditionalStat" = "tfbsHist"
                                     )
                                )                         
                     ) ##  conditional TFBS
                    #,                    
                    # conditionalPanel(
                    #     condition = "input.Raw == 'PPI'",                        
                    #     h6("Select an option"),
                    #     selectInput("ppiStat", "",
                    #                  c("General info" = "ppiInfo",
                    #                    "Topological info" = "ppiTopo")                               
                    #        )                         
                     #) ##  conditional PPI
                    ),
                    br(),br(),br(),br(),
                    value="RawTab"
                ), ## Raw Data
                     
                tabPanel("Results",
                    h6("Select results to display"),
                    selectInput("Results", "",
                               c("Networks" = "Nets",
                                "Clusters"= "clusters")),
                    conditionalPanel(
                         condition = "input.Results == 'Nets'",                         
                              h6("Statistics:"),
                              selectInput("netStat", "",
                                   c("All networks size distribution" = "allNetSizeDist",
                                     "Chromatin Maintainers" = "hldaNets",
                                     "GO analysis" ="netsGO"
                                     )
                                ),
                              sliderInput("cutThresholdInput",label="Cut-off threshold",
                                          min=0, max=1, value=0.5),
                        
                        ## displayed when the Chromatin maintainer networks are selected
                        conditionalPanel(
                            condition = "input.netStat == 'hldaNets'",
                            radioButtons("rbtnNets", "Info to display :",
                             c("Elements list" = "cmnList",
                               "Visualize networks" = "visNets")
                               ),
                            
							  conditionalPanel(
								condition="input.netStat == 'hldaNets' && input.rbtnNets == 'visNets'",
								h6("Select network to display."),
								selectInput("Maintainers", "",
								   c())
							  )                            
                        )
                     ),
					 conditionalPanel(
						condition= "input.Results == 'clusters'",
						h6("Cluster information:"),
						selectInput("clusterInfo","",
							c("Heatmap" = "clusHeatmap",
							  "Enrichment profile" = "clusCurves",
							  "Genes per cluster" = "genePerCluster",
							  "Cluster circos" = "clusCircos"						  
							),                        
						),
                        conditionalPanel(
                            condition= "input.clusterInfo == 'clusCircos'",
                            h6("Select cluster to analyse:"),
                            selectInput("ClusListInput", "",
                                       c())
                            )
                     ),                    
                   br(),br(),br(),br(),
                   value="ResultsTab"
                ),  ## Results
          id="DataType"
        ) ## tabsetPanel
),## sidebarPanel

  # Show a plot of the generated distribution
  mainPanel(
  includeHTML("d3.min.js"),
  includeHTML("graph.js"),
      
  wellPanel(
    conditionalPanel(
        condition = "input.DataType == 'RawTab' && input.petStat == 'petInfo' && input.Raw == 'PET'",
        h3("General info"),
        tableOutput("petInfoTab"),
        #h3("Number of interacting DNA regions per chromosome"),
        #plotOutput("distPlot", width="90%")        
        showOutput("distPlot", "polycharts")        
    ),
      
    conditionalPanel(
        condition = "input.DataType == 'RawTab' && input.petStat == 'petHeatmap' && input.Raw == 'PET'",        
        h3("Interaction heatmap"),
        showOutput("heatPlot", "polycharts")        
    ),
      
    conditionalPanel(
        condition = "input.DataType == 'RawTab' && input.petStat == 'span' && input.Raw == 'PET'",
        #h3("Intra-Interactions span distribution"),
        #plotOutput("petSpanPlot"),               
        showOutput("petSpanPlot", "highcharts"),  
        uiOutput("slider")
    ),
    
    conditionalPanel(
        condition = "input.DataType == 'RawTab' && input.petStat == 'circos' && input.Raw == 'PET'",
        h3("Circos map"),
        plotOutput("circosPlot",height="500px"),
        #TODO: try to fix this
        #h3("Interaction in the region"),
        #plotOutput("trackPlot")        
    ),
      
    conditionalPanel(
        condition = "input.DataType == 'RawTab' && input.Raw == 'TFBS' && input.tfbsStat == 'tfbsInfo' ",
        h3("Some statistics about the TF binding site"),        
        dataTableOutput("tfbsStatTable")               
    ),
      
    conditionalPanel(
        condition = "input.DataType == 'RawTab' && input.Raw == 'TFBS' && input.tfbsStat == 'tfbsHist' ",
        h3("Additional statistics"),
        showOutput("tfbsPerRegionPlot", "polycharts")        
    ),
      
    conditionalPanel(
        condition = "input.DataType == 'RawTab' && input.Raw == 'PPI'",
        h3("Used PPI general info"),
        tableOutput("ppiTableOutput"),
        h3("Degree distribution"),
        showOutput("ppiDegreeDist", "polycharts")
    ),
     
	 conditionalPanel(
        condition = "input.DataType == 'ResultsTab' && input.Results == 'Nets' && input.netStat == 'allNetSizeDist'",
        h3("Number of nodes per network"),
		showOutput("MainterNetsDist", "polycharts")        
    ),
	 
    conditionalPanel(
        condition = "input.DataType == 'ResultsTab' && input.Results == 'Nets' && input.netStat == 'hldaNets' && input.rbtnNets == 'cmnList' ",
        h3("Member of the chromatin maintainer networks"),		
        dataTableOutput("hldaNetTableOutput")
    ),
      
    conditionalPanel(
        condition = "input.DataType == 'ResultsTab' && input.Results == 'Nets' && input.netStat == 'hldaNets' && input.rbtnNets == 'visNets' ",
        h3("chromatin maintainer networks"),
        
        reactiveNetwork("hldaNetOutput")
    ),
	
    conditionalPanel(
        condition = "input.DataType == 'ResultsTab' && input.Results == 'clusters' && input.clusterInfo == 'clusHeatmap' ",
        h3("DNA enrichment profile"),
        plotOutput("clusHeatmapOutput")
    ),
    
	conditionalPanel(
        condition = "input.DataType == 'ResultsTab' && input.Results == 'clusters' && input.clusterInfo == 'clusCurves' ",
        h3("DNA enrichment profile"),
        plotOutput("clusCurvesOutput")
    ),
	
	conditionalPanel(
        condition = "input.DataType == 'ResultsTab' && input.Results == 'clusters' && input.clusterInfo == 'genePerCluster' ",
        h3("Genes per cluster"),
        dataTableOutput("genesPerClusStatTable")
    ),
    conditionalPanel(
        condition = "input.DataType == 'ResultsTab' && input.Results == 'clusters' && input.clusterInfo == 'clusCircos' ",
        h3("Circos map of the cluster interactions"),
        plotOutput("clusCircosOutput", height="60%")
    )
  , height="100%")
  )
))