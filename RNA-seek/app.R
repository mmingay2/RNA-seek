library(shiny)
library(shinydashboard)
library(datasets)
library(ggplot2)
library(shinythemes)
library(tidyr)
library(reshape2)
library(data.table)
library(dplyr)
library(DT)

##Set theme

commonTheme = list(labs(color="Density",fill="Density"),
                   theme_bw(),
                   theme(legend.position=c(0,1),
                         legend.justification=c(0,1)))

##Load and process Data 

rpks <- as.data.frame(fread("./rpks.txt"))
allg <- readRDS(file = "./msigdesc.rds")
msigdb <- readRDS("./msigdb2.rds")
datab <- readRDS(file="./vitclists.RDS")

##Create User Interface object which designs the webpage and allows for user input through input$"inputID"

ui <- navbarPage(title="Vitamin C induced epigenetic remodelling", theme = shinythemes::shinytheme("sandstone"),
                 tabPanel(title = "RNA-seek",
                          h4(textOutput("apptitle1")),
                          h4(textOutput("apptitle2")),
                          h4(tags$a(href="https://www.nature.com/articles/leu2017171", "Click here to read my paper and learn more")),
                          br(),
                          fluidPage(
                            sidebarLayout(
                            sidebarPanel(
                              selectInput("fillm", "Select an msigDB term to highlight it's associated genes:", c(allg, names(datab)), multiple = F, selected = allg[4748])
                            ),
                            # Show the caption and plot of the requested variable against mpg
                            mainPanel(
                              plotOutput("rnaPlot", 
                                         width = "125%",  height = "125%",
                                         hover = "plot_hover",
                                         click = "plot_click"),
                              br(),
                              h4(textOutput("tabletitle")),
                              br(),
                              dataTableOutput("infob")
                          )
                        )
                      )
               )
)


#Create a server object that generates output$ plots or tables that are reactitve '({})' to changes in the UI (input$*) to update the output (output$*)

server <- shinyServer(function(input, output) {
  
  plata <- reactive({ 
    x <- log(1+rpks[,paste0("unt_1")])
    y <- log(1+rpks[,paste0("vitc_1")])
    z <- rpks$gene %in% c(msigdb[[input$fillm]], toupper(datab[[input$fillm]]$V1))
    #z <- rpks$gene %in% toupper(data[[input$fillm]]$V1)
    gene <- rpks$gene
    ensid <- rpks$ensid
    df <- data.frame(gene, ensid, x,y,z)
    colnames(df) <- c("gene", "ensid", "Untreated", "VitaminC", "InGroup")
    df
  })
  
  
  output$rnaPlot <- renderPlot({
    plata() %>%
      ggplot(aes(x=Untreated, y=VitaminC, colour=InGroup, size=InGroup, alpha=InGroup))+
      # ggplot(data=rpks, aes(x=log(unt_2), y=log(vitc_2),colour=get(input$fillx))+
      geom_point(position="jitter")+
      #   geom_point(data=subset(rpks, input$fillx != "n.s"), aes(x=log(get(paste0("unt_",input$exp_num))), y=log(get(paste0("vitc_",input$exp_num))), colour=get(input$fillx)), position="jitter", size=1.2)+
      theme_bw()+
      xlab("Untreated")+
      ylab("Vitamin C")+
      scale_size_manual(values=c(0.5,1))+
      scale_alpha_manual(values=c(0.25, 1))+
      guides(alpha=FALSE, size=FALSE)+
      scale_colour_manual(values=c("gray65", "red"), guide = guide_legend(title = "ingroup?"))+
      geom_abline(slope=1, color="gray10")+
      ggtitle("Click Points on Plot to view Gene ID and RPKM")
    #coord_cartesian(xlim=c(-2,4), ylim=c(-2,4))
  }, width = 700, height = 400)
  
  output$infob <- renderDataTable({
    plata() %>% nearPoints(input$plot_click, maxpoints = 25, threshold = 10) %>% mutate(fc=log2(((10^VitaminC)-1)/((10^Untreated)-1))) %>% dplyr::rename("log(1+RPKM) VitC treated" = "VitaminC", 
                                                                      "log(1+RPKM) untreated" = "Untreated", 
                                                                      "Gene Name" = "gene", 
                                                                      "Ensembl ID" = "ensid", 
                                                                      "RPKM log2 Fold Change (VitC/Unt)" = "fc",
                                                                      "In MSigDB list?" = "InGroup")
  })
  
  output$tabletitle <- renderText("Cursor click the scatterplot above to display genes nearby the click position")
  output$apptitle1 <- renderText("Use this application to explore pairwise gene expression differences between 2 RNA-seq libraries")
  output$apptitle2 <- renderText("These libraries come from Untreated and Vitamin C treated (1mM) IDH1 R132H expression mouse bone marrow cells")
})

shinyApp(ui = ui, server = server)

