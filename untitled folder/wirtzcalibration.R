 #### Calibration work-flow for Lobster GS model application
##Wirtz et al. 2022 data
# Step 0. Create reverse mode input files for each experimental dataset
# Step 1. Run app.R to explore parameters to tune model to observations in reverse mode
# Step 2. Determine fV by fitting regression of Lipid vs Protein intakes across feeds for the experimental observation and run as inputs in forward mode
# Step 3. Using fV and model inputs estimate the best parameter set for each feed (lowest error between pred and obs) for the selected parameter(s) to be tuned by using optimization or grid-search# Step 3. Using fV and model inputs estimate the best parameter set for each feed (lowest error between pred and obs) for the selected parameter(s) to be tuned by using optimization or grid-search
# Step 4. Repeat Steps 0-3  for each experiment and make plot of pred vs obs across all model and data. 
# Step 5. Using best parameter sets for all experiments, run model in forward mode to predict N waste vs % Protein in Feed (and could predict other outputs such as growth, respiration etc., to explore further)

#### clear workspace
rm(list=ls())
install.packages("tidyverse")
library(conflicted) 
library(shiny)
library(plotly)
library(DT)
library(tidyverse)
library(pbapply)
library(ggrepel)

##### Experiment 1 - Wirtz data

## Step 0 - get data

#source code with function to run (this is a script that only contains functions not analyses of model)
source("GS_model_functions.R")
#get input data
wirtz<-read.csv("wirtzinput.csv")
wirtz$Protein_intake<-NA
wirtz$Lipid_intake<-NA
input<-read_csv(("wirztobs.csv"))
obswirtz<-read_csv(("wirztobs.csv"))

## Step 1 - run app
#source("app.R")
#shiny_gs()
## interactively, this indicates kstarN as candidate parameter

## Step 2 - get fV values

#use linear regression or exact observed intake values for each feed

coefs<-lm(data=input,Lipid_intake ~Protein_intake-1)
predLI<-predict(coefs)

# or for each data point:

out<-pbapply(X=as.matrix(input[,c(2,3,5,6)]),1,GSforwardmodel,kstarN=0.9,phi=0)
fV<-c(out[[1]]$fV,out[[2]]$fV,out[[3]]$fV,out[[4]]$fV)

# # add to obsinput

input$fV<-fV

## Step 3 - get best parameter set for each feed type

## added function getError to  GS_model_functions

# for more than one variable to estimate need to use optim

FMest<-optim(c(0.9,0.5),getError,input=input,feed="FM",paramname="both",method="L-BFGS-B",upper=c(1.0,1.0),lower=c(0.00001,0.00001))
KMest<-optim(c(0.9,0.5),getError,input=input,feed="KM",paramname="both",method="L-BFGS-B",upper=c(1.0,1.0),lower=c(0.00001,0.00001))
SBMest<-optim(c(0.9,0.5),getError,input=input,feed="SBM",paramname="both",method="L-BFGS-B",upper=c(1.0,1.0),lower=c(0.00001,0.00001))
SWMest<-optim(c(0.9,0.5),getError,input=input,feed="SWM",paramname="both",method ="L-BFGS-B",upper=c(1.0,1.0),lower=c(0.00001,0.00001))
pars<-rbind(FMest$par, KMest$par,SBMest$par,SWMest$par)

### check to see if these make sense 
feed=c("FM","KM","SBM","SWM")
input$pred_Protein_intake<-NA
input$pred_Lipid_intake<-NA
for (i in 1:4) {
  finput<-data.frame(G=input[input$Feed==feed[i],"G"],fV=input[input$Feed==feed[i],"fV"],betaV=input[input$Feed==feed[i],"betaV"],betaH=input[input$Feed==feed[i],"betaH"], input[input$Feed==feed[i],"Lipid_intake"],input[input$Feed==feed[i],"Protein_intake"])
  # used optimised
  #input[i,c(8:9)]<-GSreversemodel(finput,kstarN = kstarNEst[i])
  input[i,c(8:9)]<-GSreversemodel(finput,kstarN = pars[i,1],phi = pars[i,2])
  # or used guessed values from shiny
  # input[i,c(8:9)]<-GSreversemodel(finput,kstarN = c(1,0.95,1,0.6))
}
#adding kstarN and phi to input dataframe
input$kstarN<-pars[,1]
input$phi<-pars[,2]

#importing model predictions for fv=0-1 and swm with 0 penalty datasets
input1<-read.csv("wirtzmodelpredictions1204.csv")
input2<-read.csv("swm0penaltyoutput.csv")

#plotting input arrays, comparing calibrated vs uncalibrated vs expt
library(viridisLite)
library(viridis)
library(ggplot2)
wirtzplot <- ggplot(input,aes(x = Protein_intake, y = Lipid_intake, group = G)) +
  geom_point(size = 3.5, aes(color = G))+
  geom_point(size = 5, shape = "X", data = input, aes(x = pred_Protein_intake, y = pred_Lipid_intake, color = G)) + 
  geom_abline(slope = unlist(coefs[1]), intercept = 0, linetype = 2) +
  geom_path(data = input1, aes(x = Protein_intake, y = Lipid_intake, color = G)) +  scale_color_continuous(type="viridis")+
  geom_line(linetype = "dashed", color = "#FDE725FF", data = input2, aes(x = Protein_intake, y = Lipid_intake, group = G)) + 
  lims(x = c(0, 0.11), y = c(0, 0.06)) + 
  labs(size = 4, x = expression("Protein intake mol C mol C"^"-1"*"day"^"-1"), y = expression("Lipid intake mol C mol C"^"-1"*"day"^"-1"), color = "growth rate") +
  theme(axis.title = element_text(size = 14)) +
  theme(legend.position = c(.98, .98),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)) +
  geom_text_repel(size = 3, data = input, aes(label = Feed)) + 
  annotate("text", x = 0.008, y = 0.008, label = "fV=0.095", angle = 0)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

# not too bad! the only feed with phi > 0 is SWM
# add lines for all values of fV but using these 

# select only the calibrated modelcolumns need for forward mode now

#renaming pred_protein_intake and pred_Lipid_intake as IV and IH for saving df to run in forward mode
names(input)[8]<-"Protein_intake"
names(input)[9]<-"Lipid_intake"
#save df as RDS just in case you need to go back and ref and you've cleared the environment
library(readr)
#saveRDS(input,"~/OneDrive - University of Tasmania/Sowdamini_R_PhDWork/Salmon_Geometric_Stoichiometry/wirtzpredictions4fwdmode.RDS")
# below automatically saves in your  project directory
saveRDS(input,"newwirtzinput.RDS")


#wirtzfwdmode<-input
#now we attempt to do wirtz waste plot in forward mode
wirtzfwdmode<-readRDS("newwirtzinput.RDS")
rownames<-wirtzfwdmode$Feed
#select  which columns are needed to run in forwrad mode
wirtzfwdmode<-as.matrix(wirtzfwdmode[,c(5,6,8,9,10,11)])
rownames(wirtzfwdmode)<-rownames
# names(wirtzfwdmode)[1]<-"Protein_intake"
# names(wirtzfwdmode)[2]<-"Lipid_intake"
# library(pbapply)
# library(ggrepel)
# library(ggplot2)
# library(viridis)
#here, remember to modify and check the function,check values of kstarN and phi in fwd mode. 

# below the columsn being selected  needed to be fixed:
out2<-pbapply(X=wirtzfwdmode,1,GSforwardmodel)
# test run holding kstarN and phi constant
test2<-pbapply(X=wirtzfwdmode[,c(1,2,3,4)],1,GSforwardmodel,kstarN=0.9,phi=0)

#out2<-pbapply(X=wirtzfwdmode,1,GSforwardmodel)
out3<-c(out2[[1]]$Pellet_N,out2[[2]]$Pellet_N,out2[[3]]$Pellet_N,out2[[4]]$Pellet_N)

# compare with test values
test3<-c(test2[[1]]$Pellet_N,test2[[2]]$Pellet_N,test2[[3]]$Pellet_N,test2[[4]]$Pellet_N)
# waste looks the same

outg<-c(out2[[1]]$Growth_N,out2[[2]]$Growth_N,out2[[3]]$Growth_N,out2[[4]]$Growth_N)
testg<-c(test2[[1]]$Growth_N,test2[[2]]$Growth_N,test2[[3]]$Growth_N,test2[[4]]$Growth_N)
# growth is very different

#now convert  to  df
wirtzfwdmode<-data.frame(wirtzfwdmode)
wirtzfwdmode$Nwaste<-out3
wirtzfwdmode$CP_feed<-c(63.75,65.1,60.05,60.55)
wirtzfwdmode$Feed<-rownames

#Plotting 
# could not read in the file below using read.csv
library(RColorBrewer)
install.packages("wesanderson")
library(wesanderson)
library(ggplot2)
library(viridis)

wasteexpt<-read_csv("wastecalculations.csv")
wwaste<-ggplot(wirtzfwdmode,aes(x=CP_feed,y=Nwaste,color=Feed))+geom_point(size=8)+
  geom_point(shape="x",size=8,data=wasteexpt,aes(x=CP_feed,y=Nwaste,color=Feed))+ geom_text_repel(size=7.5,data=wasteexpt,aes(label = Feed))+ 
  labs(x = expression("Percentage of crude protein in feeds"), y = expression("N waste output mol N mol C"^"-1"*"day"^"-1")) + 
  theme(axis.title = element_text(size = 20),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(legend.position="none")+
  lims(x=c(52,67),y=c(0,0.010))
#+scale_color_manual(values = wes_palette("GrandBudapest1"))

#axis.title.x=element_blank(),axis.title.y=element_blank(),
##I see that the plot has remained unchanged
#trying to see if the df is the issue here
dim(wirtzfwdmode)

#check if the model picks up the right values for kstarN and phi 
pbapply(wirtzfwdmode[, c("kstarN", "phi")], 1, mean)
# Define a function to print the values of kstarN and phi
printValues <- function(x) {
  cat("kstarN =", x["kstarN"], "\n")
  cat("phi =", x["phi"], "\n\n")
}

#looks ok

# Call pbapply with the printValues function
pbapply(wirtzfwdmode, 1, printValues)
#the df doesnt seem to be issue here then


#calculating diff b/w pred and observed waste
