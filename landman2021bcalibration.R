#### Calibration work-flow for Lobster GS model application

# Step 0. Create reverse mode input files for each experimental dataset
# Step 1. Run app.R to explore parameters to tune model to observations in reverse mode
# Step 2. Determine fV by fitting regression of Lipid vs Protein intakes across feeds for the experimental observation and run as inputs in forward mode
# Step 3. Using fV and model inputs estimate the best parameter set for each feed (lowest error between pred and obs) for the selected parameter(s) to be tuned by using optimization or grid-search# Step 3. Using fV and model inputs estimate the best parameter set for each feed (lowest error between pred and obs) for the selected parameter(s) to be tuned by using optimization or grid-search
# Step 4. Repeat Steps 0-3  for each experiment and make plot of pred vs obs across all model and data. 
# Step 5. Using best parameter sets for all experiments, run model in forward mode to predict N waste vs % Protein in Feed (and could predict other outputs such as growth, respiration etc., to explore further)
library(tidyverse)
library(pbapply)
library(ggrepel)
## Step 0 - get data

#source code with function to run (this is a script that only contains functions not analyses of model)
source("GS_model_functions.R")
#get input data from Landman et al 2021a
landman2<-read.csv("landman2021binputs.csv")
landman2$Protein_intake<-NA
landman2$Lipid_intake<-NA
obslandman2<-read.csv("landman2021bobs.csv")
input<-read.csv("landman2021bobs.csv") 

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
fV<-c(out[[1]]$fV,out[[2]]$fV,out[[3]]$fV,out[[4]]$fV,out[[5]]$fV,out[[6]]$fV,out[[7]]$fV)

# # add to obsinput

input$fV<-fV

## Step 3 - get best parameter set for each feed type
# for more than one variable to estimate need to use optim

BM0est<-optim(c(0.9,0.5,0.9,0.7),getError,input=input,feed="BM0%",paramname="withbetas",method="L-BFGS-B",upper=c(1.0,1.0,1.0,1.0),lower=c(0.00001,1e-10,1e-10,1e-10))
BM16est<-optim(c(0.9,0.5,0.9,0.7),getError,input=input,feed="BM1.6%",paramname="withbetas",method="L-BFGS-B",upper=c(1.0,1.0,1.0,1.0),lower=c(0.00001,1e-10,1e-10,1e-10))
BM31est<-optim(c(0.9,0.5,0.9,0.7),getError,input=input,feed="BM3.1%",paramname="withbetas",method="L-BFGS-B",upper=c(1.0,1.0,1.0,1.0),lower=c(0.00001,1e-10,1e-10,1e-10))
BM63est<-optim(c(0.9,0.5,0.9,0.7),getError,input=input,feed="BM6.3%",paramname="withbetas",method ="L-BFGS-B",upper=c(1.0,1.0,1.0,1.0),lower=c(0.00001,1e-10,1e-10,1e-10))
BM125est<-optim(c(0.9,0.5,0.9,0.7),getError,input=input,feed="BM12.5%",paramname="withbetas",method ="L-BFGS-B",upper=c(1.0,1.0,1.0,1.0),lower=c(0.00001,1e-10,1e-10,1e-10))
BM250est<-optim(c(0.9,0.5,0.9,0.7),getError,input=input,feed="BM25%",paramname="withbetas",method ="L-BFGS-B",upper=c(1.0,1.0,1.0,1.0),lower=c(0.00001,1e-10,1e-10,1e-10))
BMHSest<-optim(c(0.9,0.5,0.9,0.7),getError,input=input,feed="BMHS",paramname="withbetas",method ="L-BFGS-B",upper=c(1.0,1.0,1.0,1.0),lower=c(0.00001,1e-10,1e-10,1e-10))

pars<-rbind(BM0est$par,BM16est$par,BM31est$par,BM63est$par,BM125est$par,BM250est$par,BMHSest$par)
colnames(pars)<-c("kstarN","phi","betaV","betaH")
rownames(pars)<-input$Feed

#replace input with calibrated values
input[,"betaV"]<-pars[,"betaV"]
input[,"betaH"]<-pars[,"betaH"]
input$phi<-input$kstarN<-NA
input[,"kstarN"]<-pars[,"kstarN"]
input[,"phi"]<-pars[,"phi"]
 
## replace matrix 10 times for each value of fV 0:1
M<-as.matrix(input[,-1])
input_array<-matrix(rep(t(M),10),ncol=dim(M)[2],byrow=TRUE)
colnames(input_array)<-names(input[,-1])
input_array<-data.frame(Feed=input[,1],input_array)
# replace fVs with 0.1 to 1 instead of calibrated values
input_array[,"fV"]<-rep(seq(0.1,1,0.1),each=7)

### check to see if these make sense 

# create full array of inputs with varying fV form 0:1

# write over observed Protein and Lipid intake values with predicted ones
for (i in 1:dim(input_array)[1]) input_array[i,c("pred_Protein_intake","pred_Lipid_intake")]<-unlist(GSreversemodel(input_array[i,],kstarN=input_array[i,"kstarN"],phi=input_array[i,"phi"],betaV=input_array[i,"betaV"],betaH=input_array[i,"betaH"]))



#### get predicted calibrated values for plotting
#feed=c("BM0","BM16","BM31","BM63","BM125","BM250","BMHS")
input$pred_Protein_intake<-NA
input$pred_Lipid_intake<-NA
for (i in 1:7){
  #finput<-data.frame(G=input[input$Feed==feed,"G"],fV=input[input$Feed==feed,"fV"],betaV=input[input$Feed==feed,"betaV"],betaH=input[input$Feed==feed,"betaH"], Protein_intake=input[input$Feed==feed,"Protein_intake"],Lipid_intake=input[input$Feed==feed,"Lipid_intake"])
  #finput<-data.frame(G=input[input$Feed==feed[i],"G"],fV=input[input$Feed==feed[i],"fV"],betaV=input[input$Feed==feed[i],"betaV"],betaH=input[input$Feed==feed[i],"betaH"],input[input$Feed==feed[i],"Lipid_intake"],input[input$Feed==feed[i],"Protein_intake"])
  # used optimised
  #input[i,c(8:9)]<-GSreversemodel(finput,kstarN = kstarNEst[i])
  input[i,c(10:11)]<-GSreversemodel(as.matrix(input[i,-1]),kstarN = pars[i,"kstarN"],phi = pars[i,"phi"],betaV =  pars[i,"betaV"],betaH =  pars[i,"betaH"])
}


library(ggplot2)
library(ggrepel)
library(readr)
write_csv(input,"~/OneDrive - University of Tasmania/Sowdamini_R_PhDWork/Salmon_Geometric_Stoichiometry/landman2predictions.csv")

linput2<-read.csv("landman2modelpredictions.csv")
input<-as.data.frame(input)
landman2021bplot<-ggplot(input,aes(x=Protein_intake,y= Lipid_intake,group=G)) + 
  scale_color_continuous(type="viridis") +
  geom_point(size=3.5)+ 
  geom_point(data=input,size=5,shape="X",aes(x=pred_Protein_intake,y=pred_Lipid_intake,color=G))+ 
  geom_abline(slope=unlist(coefs[1]), intercept=0,linetype=2) +
  geom_line(data=linput2,aes(x=Protein_intake,y=Lipid_intake,color=G))+ 
  lims(x=c(0,0.11),y=c(0,0.06))+
  labs(size=4,x=expression("Protein intake mol C mol C"^"-1"*"day"^"-1"),y=expression("Lipid intake mol C mol C"^"-1"*"day"^"-1"),color="growth rate")+
  theme(axis.title = element_text(size=14))+
  theme(legend.position = c(0.98,0.98),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)) +
  geom_text_repel(size=3,data=input,aes(label = Feed))+annotate("text", x=0.01, y=0.006, label= "fV=0.18", angle= 0)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

#renaming pred_protein_intake and pred_Lipid_intake as IV and IH for saving df to run in forward mode
names(input)[10]<-"Protein_intake"
names(input)[11]<-"Lipid_intake"
#save df as RDS just in case you need to go back and ref and you've cleared the environment
library(readr)
#saveRDS(input,"~/OneDrive - University of Tasmania/Sowdamini_R_PhDWork/Salmon_Geometric_Stoichiometry/wirtzpredictions4fwdmode.RDS")
# below automatically saves in your  project directory
saveRDS(input,"newlandman2input.RDS")
#wasteplot 
landman2fwdmode<-readRDS("newlandman2input.RDS")
rownames<-landman2fwdmode$Feed
landman2fwdmode<-as.matrix(landman2fwdmode[,c(5,6,8,9,10,11)])
rownames(landman2fwdmode)<-rownames
out6<-pbapply(X=landman2fwdmode,1,GSforwardmodel)
out7<-c(out6[[1]]$Pellet_N,out6[[2]]$Pellet_N,out6[[3]]$Pellet_N,out6[[4]]$Pellet_N,out6[[5]]$Pellet_N,out6[[6]]$Pellet_N,out6[[7]]$Pellet_N)
landman2fwdmode$Nwaste<-out7
landman2fwdmode$CP_feed<-c(62.31,62.68,61.19,60.82,61.44,61.49,54.31)
landman2fwdmode<-as.data.frame(landman2fwdmode)
l2waste<-ggplot(landman2fwdmode,aes(x=CP_feed,y=Nwaste,color=rownames))+geom_point()+ geom_text_repel(size=4,color="black",data=landman2fwdmode,aes(label=rownames))+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(legend.position="none")+
  lims(x=c(52,67),y=c(0,0.010))+ scale_color_manual(values = wes_palette("Royal1"))


