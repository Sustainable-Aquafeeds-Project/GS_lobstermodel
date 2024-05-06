#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
library(shiny)
library(plotly)
library(DT)
library(tidyverse)
library(pbapply)
library(ggrepel)

#source code with function to run (this is a script that only contains functions not analyses of model)
source("GS_model_functions.R")
#get input data
#wirtz<-read.csv("wirtzinput.csv")
#wirtz$Protein_intake<-NA
#wirtz$Lipid_intake<-NA
#obswirtz<-read_csv(("wirztobs.csv")) 

#input data of all experiments are in their individual calibration scripts



#newinput<-data.frame(G=0.0057,fV=0.5,betaV=input$scalebetaV*0.655,betaH=input$scalebetaH*0.53)
# newinput<-data.frame(G=0.0057,fV=0.5,betaV=0.655,betaH=0.53)
# 
# newinput$Protein_intake<-NA
# newinput$Lipid_intake<-NA
# newinput[1,c(5:6)]<-GSreversemodel(input=newinput ,tau=input$tau,eta=input$eta,kstarN=input$kstarN,thetaZ=input$thetaZ,thetaV=input$thetaV)

shiny_gs <- function(expinput=,obs=obslandman1) {

  ui=fluidPage(
    # Application title
    titlePanel("Geometric Stoichiometry Model: Calibration"),
    fluidRow(
      column(4, wellPanel(
        sliderInput("tau", "Biomass turnover, d-1:", min = 0.01, max = 0.1, value = 0.02,
                    step = 0.01),
        sliderInput("xi", "Other Basal Metabolism, d-1:", min = 1e-5, max = 1e-3, value =1e-4,
                    step = 1e-5),
        sliderInput("eta", "Specific Dynamic Action:", min = 0.01, max = 0.1, value = 0.03,
                     step = 0.01),
        sliderInput("kstarN", "Net N Synthesis Efficiency:", min = 0.26, max = 1.0, value = 0.9,
                    step = 0.1),  
        sliderInput("thetaZ", "C:N of Consumer:", min = 2, max = 6, value = 5.2,
                                              step = 0.1),
        sliderInput("thetaV", "C:N of Protein:", min = 2, max = 6, value = 3.7,
                    step = 0.1),
        sliderInput("phi", "Penalty for Protein:", min = 0, max = 1, value = 0,
                    step = 0.1)
          )),
      column(6,
             plotOutput("plot1", width = 600, height = 600)
      )),

    
  )



  server = function(input,output,expinput=landman1,obs=obslandman1) {
    
    output$plot1 <- renderPlot({
      # set up params using values given, need check 
      out<-data.matrix(unlist(pbapply(X=expinput,1,GSreversemodel,tau=input$tau,xi=input$xi,eta=input$eta,kstarN=input$kstarN,thetaZ=input$thetaZ,thetaV=input$thetaV, phi=input$phi,simplify=T)))
      expinput[,c(5:6)]<-cbind(out[which(rownames(out)=="Protein_intake")],out[which(rownames(out)=="Lipid_intake")])
      coefs<-lm(data=obs,Lipid_intake ~Protein_intake-1)
      # add new point on the plot - at fV=0.04, G = 0.006, and rescale betaV and betaH
      # need to work out how to enter in a data table interactively
      newpoints<-GSreversemodel(input = c(0.0057,0.5,input$scalebetaV*0.655,input$scalebetaH*0.53),tau=input$tau,eta=input$eta,kstarN=input$kstarN,thetaZ=input$thetaZ,thetaV=input$thetaV)
      
       ggplot(expinput,aes(x=Protein_intake,y= Lipid_intake,group=G)) + 
        scale_color_continuous(type="viridis") +
        geom_path(aes(color=G)) +
        geom_point(size=3.5,data=obs,aes(color=G))+ 
        geom_abline(slope=unlist(coefs[1]), intercept=0,linetype=2) + 
        lims(x=c(0,0.11),y=c(0,0.08))+ 
        labs(size=4,x=expression("Protein intake mol C mol C"^"-1"*"day"^"-1"),y=expression("Lipid intake mol C mol C"^"-1"*"day"^"-1"),color="growth rate")+
        theme(axis.title = element_text(size=14))+
        theme(legend.position = c(.98, .5),
              legend.justification = c("right", "top"),
              legend.box.just = "right",
              legend.margin = margin(6, 6, 6, 6)) +
        geom_text_repel(size=4,data=obs,aes(label = Feed)) 
      
      
      })
     

  }
  shinyApp(ui, server)
}

shiny_gs()
