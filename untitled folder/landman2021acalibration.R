#### Calibration work-flow for Lobster GS model application

# Step 0. Create reverse mode input files for each experimental dataset
# Step 1. Run app.R to explore parameters to tune model to observations in reverse mode
# Step 2. Determine fV by fitting regression of Lipid vs Protein intakes across feeds for the experimental observation and run as inputs in forward mode
# Step 3. Using fV and model inputs estimate the best parameter set for each feed (lowest error between pred and obs) for the selected parameter(s) to be tuned by using optimization or grid-search# Step 3. Using fV and model inputs estimate the best parameter set for each feed (lowest error between pred and obs) for the selected parameter(s) to be tuned by using optimization or grid-search
# Step 4. Repeat Steps 0-3  for each experiment and make plot of pred vs obs across all model and data. 
# Step 5. Using best parameter sets for all experiments, run model in forward mode to predict N waste vs % Protein in Feed (and could predict other outputs such as growth, respiration etc., to explore further)

## Step 0 - get data

#source code with function to run (this is a script that only contains functions not analyses of model)
source("GS_model_functions.R")
#get input data from Landm0an et al 2021a
landman1<-read.csv("landman2021input.csv")
landman1$Protein_intake<-NA
landman1$Lipid_intake<-NA
obslandman1<-read.csv("landman2021aobs.csv")
input<-read.csv("landman2021aobs.csv") 

## Step 1 - run app
#source("app.R")
#shiny_gs()
## interactively, this indicates kstarN as candidate parameter

## Step 2 - get fV values

#use linear regression or exact observed intake values for each feed

coefs<-lm(data=input,Lipid_intake ~Protein_intake-1)
predLI<-predict(coefs)
coefs
predLI

# or for each data point:


out<-pbapply(X=as.matrix(input[,c(2,3,5,6)]),1,GSforwardmodel,kstarN=0.9,phi=0)
fV<-c(out[[1]]$fV,out[[2]]$fV,out[[3]]$fV,out[[4]]$fV,out[[5]]$fV,out[[6]]$fV)

# # add to obsinput

input$fV<-fV

## Step 3 - get best parameter set for each feed type
# for more than one variable to estimate need to use optim

D1est<-optim(c(0.9,0.5),getError,input=input,feed="D1",paramname="both",method="L-BFGS-B",upper=c(1.0,1.0),lower=c(0.00001,0.00001))
D2est<-optim(c(0.9,0.5),getError,input=input,feed="D2",paramname="both",method="L-BFGS-B",upper=c(1.0,1.0),lower=c(0.00001,0.00001))
D3est<-optim(c(0.9,0.5),getError,input=input,feed="D3",paramname="both",method="L-BFGS-B",upper=c(1.0,1.0),lower=c(0.00001,0.00001))
D4est<-optim(c(0.9,0.5),getError,input=input,feed="D4",paramname="both",method ="L-BFGS-B",upper=c(1.0,1.0),lower=c(0.00001,0.00001))
D5est<-optim(c(0.9,0.5),getError,input=input,feed="D5",paramname="both",method ="L-BFGS-B",upper=c(1.0,1.0),lower=c(0.00001,0.00001))
D6est<-optim(c(0.9,0.5),getError,input=input,feed="D6",paramname="both",method ="L-BFGS-B",upper=c(1.0,1.0),lower=c(0.00001,0.00001))
pars<-rbind(D1est$par,D2est$par,D3est$par,D4est$par,D5est$par,D6est$par)


### check to see if these make sense 

feed=c("D1","D2","D3","D4","D5","D6")
input$pred_Protein_intake<-NA
input$pred_Lipid_intake<-NA
for (i in 1:6){
  finput<-data.frame(G=input[input$Feed==feed[i],"G"],fV=input[input$Feed==feed[i],"fV"],betaV=input[input$Feed==feed[i],"betaV"],betaH=input[input$Feed==feed[i],"betaH"],input[input$Feed==feed[i],"Lipid_intake"],input[input$Feed==feed[i],"Protein_intake"])
  # used optimised
  #input[i,c(8:9)]<-GSreversemodel(finput,kstarN = kstarNEst[i])
  input[i,c(8:9)]<-GSreversemodel(finput,kstarN = pars[i,1],phi = pars[i,2])
}

#add pars kstarN and phi to dataframe
input$kstarN<-pars[,1]
input$phi<-pars[,2]

library(readr)
write_csv(input,"~/OneDrive - University of Tasmania/Sowdamini_R_PhDWork/Salmon_Geometric_Stoichiometry/landman1predictions.csv")

linput1<-read.csv("landman1modelpredictions.csv")

  landman2021aplot<-ggplot(input,aes(x=Protein_intake,y= Lipid_intake,group=G)) + 
      scale_color_continuous(type="viridis") +
      # geom_path(aes(color=G)) +
      geom_point(size=3.5,aes(color=G))+ 
      geom_point(size=5,shape="X",data=input,aes(x=pred_Protein_intake,y=pred_Lipid_intake,color=G))+ 
      geom_abline(slope=unlist(coefs[1]), intercept=0,linetype=2)+ 
      geom_path(data=linput1,aes(x=Protein_intake,y=Lipid_intake,color=G))+ 
      lims(x=c(0,0.11),y=c(0,0.06))+ 
      labs(size=4,x=expression("Protein intake mol C mol C"^"-1"*"day"^"-1"),y=expression("Lipid intake mol C mol C"^"-1"*"day"^"-1"),color="growth rate")+
      theme(axis.title = element_text(size=14))+
      theme(legend.position = c(.98, .98),
            legend.justification = c("right", "top"),
            legend.box.just = "right",
            legend.margin = margin(6, 6, 6, 6)) +
      geom_text_repel(size=3,data=input,aes(label = Feed)) + annotate("text", x=0.025, y=0.007, label= "fV=0.22", angle= 0)+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
  landman2021aplot
#renaming pred_protein_intake and pred_Lipid_intake as IV and IH for saving df to run in forward mode
names(input)[8]<-"Protein_intake"
names(input)[9]<-"Lipid_intake"
#save df as RDS just in case you need to go back and ref and you've cleared the environment
library(readr)
#saveRDS(input,"~/OneDrive - University of Tasmania/Sowdamini_R_PhDWork/Salmon_Geometric_Stoichiometry/wirtzpredictions4fwdmode.RDS")
# below automatically saves in your  project directory
saveRDS(input,"newlandman1input.RDS")

#wasteplot 
landman1fwdmode<-readRDS("newlandman1input.RDS")
rownames<-landman1fwdmode$Feed
landman1fwdmode<-as.matrix(landman1fwdmode[,c(5,6,8,9,10,11)])
rownames(landman1fwdmode)<-rownames
out4<-pbapply(X=landman1fwdmode,1,GSforwardmodel)
out5<-c(out4[[1]]$Pellet_N,out4[[2]]$Pellet_N,out4[[3]]$Pellet_N,out4[[4]]$Pellet_N,out4[[5]]$Pellet_N,out4[[6]]$Pellet_N)
landman1fwdmode$Nwaste<-out5
landman1fwdmode$CP_feed<-c(60.5,59.77,60.57,61.34,61.02, 61.82)
landman1fwdmode<-data.frame(landman1fwdmode)
l1waste<-ggplot(landman1fwdmode,aes(x=CP_feed,y=Nwaste,color=rownames))+geom_point(color="black")+ geom_text_repel(size=4,color="black",data=landman1fwdmode,aes(label=rownames))+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(legend.position="none")+
  lims(x=c(52,67),y=c(0,0.010))+ scale_color_manual(values = wes_palette("Royal1"))
  
