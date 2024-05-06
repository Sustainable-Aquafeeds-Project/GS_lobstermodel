library(tidyverse)
#source code with function to run (this is a script that only contains functions not analyses of model)
GSreversemodel<- function(input,tau=0.04,xi=0.00001,eta=0.0259,kstarN=NULL,thetaZ=5.2,thetaV=3.7,phi=NULL,betaV=NULL,betaH=NULL){
  # Default Parameters
  # tau <- 0.04      # biomass turnover, d-1. CHANGING VALUES TO FIT MODEL TO DATA
  # xi <- 0.00001    # other basal metabolism, d-1 model default
  # eta <- 0.0259       # SDA from Wang et al 2021 spiny lobster data for CD60 DIET
  # betaV <- 0.69     # absorption efficiency: protein from Anderson et al 2005, changed it to input because this varies in data
  # kstarN <- 0.9     # max. N synthesis efficiency model default
  # thetaZ <- 5.2 #5.5     # C:N of consumer (Q: would this be wet or dry ratio?) #data obtained from Ikeda et al. 2011
  # 
  # thetaV <- 3.7     # C:N of protein (Q: would this be wet or dry ratio?)
  # phi <- 0     # penalty set to 0 to start off with, can change later
  
  # Specify G (growth rate, d-1) and fV
  # These settings set to match the forward model for comparison
  G <- input[,"G"] 
  fV <- input[,"fV"] # 0 <= fV <= 1

  if (is.null(betaV)) betaV<-input[,"betaV"]     
  if (is.null(betaH)) betaH <- input[,"betaH"] 
  if (is.null(kstarN)) kstarN<-input[,"kstarN"]   
  if (is.null(phi)) phi <- input[,"phi"]
  
   # absorption efficiency: lipid from Landman et al., 2021
  
  # Step 1 Calculate demand for absorbed C
  DN <- (tau+G)/kstarN         # Eq. 1
  DV <- thetaV/thetaZ*DN        # Eq. 2
  DX <- (thetaZ-thetaV)/thetaZ*(tau+G)     # Eq. 3
  
  # Step 3 Calculation of IV and IH (before step 2)
  if (fV < 1) {
    a <- fV/(1-phi*fV)
    numerator <- (DV+a*(DX+xi))/(betaV-eta*a)+(DX+xi)/eta
    denominator <- betaH*(1-phi*fV)/(eta*(1-fV))-eta*a/(betaV-eta*a)-1
    IH <- numerator/denominator            # Eq. 8
    IV <- betaH*IH*(1-phi*fV)/(eta*(1-fV))-(DX+xi)/eta-IH   # Eq. 9
  } else {
    IH <- 0.
    IV <- (DV+(DX+xi)/(1-phi))/(betaV-eta/(1-phi))
  }
  
  # Step 2 Energetic costs, including penalty
  DA <- DX + xi + eta*(IV+IH)
  if (phi > 0) {
    Omega <- DA*(1/(1-phi*fV)-1)
  } else {
    Omega <- 0.
  }
  
  # Step 4 C and N budgets
  R <- tau+(1-kstarN)*DV+xi+eta*(IV+IH)+Omega     # Eq. 11
  E <- tau/thetaZ+(1-kstarN)*DV/thetaV+fV/thetaV*(DX+xi+eta*(IV+IH)+Omega)
  # Eq. 12
  WC <- (1-betaV)*IV+(1-betaH)*IH            # Eq. 13
  WN <- (1-betaV)*IV/thetaV                  # Eq. 14
  
  IC <- IV+IH       # C intake, mol C mol C-1 d-1
  IN <- IV/thetaV   # N intake, mol N mol C-1 d-1
  Ctot <- G+R+WC    # C budget check
  GN <- G/thetaZ    # growth, mol N mol C-1 d-1
  Ntot <- GN+E+WN   # N budget check
  
  # print("         ")
  # print(c("Intake V (IV)",IV))
  # print(c("Intake H (IH)",IH))
  # print("         ")
  # print("C budget:")
  # print(c("Intake C (IV+IH)",IC))
  # print(c("Growth (G)",G))
  # print(c("Respiration (R)",R))
  # print(c("Pellet C (WC)",WC))
  # print(c("C sinks (G+R+WC)",Ctot))
  # print("         ")
  # print("N budget:")
  # print(c("Intake N (IV+IH)",IN))
  # print(c("Growth (GN)",GN))
  # print(c("Excretion (E)",E))
  # print(c("Pellet N (WN)",WN))
  # print(c("N sinks (G+R+WC)",Ntot))
  # 
  #return(data.frame(Protein_intake=IV,Lipid_intake=IH, protein_absorption_efficiency=betaV,lipid_absorption_efficiency=betaH, protein_absorption=betaV*IV, lipid_absorption=betaH*IH, Respiration=R, Excretion=E, Pellet_C=WC, Pellet_N=WN, Growth_G = G, Growth_N=GN )) #changed IH= carb intake to lipid intake
  return(list(Protein_intake=as.numeric(IV),Lipid_intake=as.numeric(IH))) 
}



#Growth forward mode

GSforwardmodel <- function(input,tau=0.04,xi=0.00001,eta=0.0259,thetaZ=5.2,thetaV=3.7,kstarN=NULL,phi=NULL,betaV=NULL,betaH=NULL){
#(input=input){
  # Parameters
 #tau <- 0.04      # biomass turnover, d-1. CHANGING VALUES TO FIT MODEL TO DATA
 #xi <- 0.00001    # other basal metabolism, d-1 model default
 #eta <- 0.0259       # SDA from Wang et al 2021 spiny lobster data for CD60 DIET
 #betaV <- 0.69     # absorption efficiency: protein from Anderson et al 2005, changed it to input because this varies in data
 #kstarN <- 0.9     # max. N synthesis efficiency model default
 #thetaZ <- 5.2 #5.5     # C:N of consumer (Q: would this be wet or dry ratio?) #data obtained from Ikeda et al. 2011
 #thetaV <- 3.7     # C:N of protein (Q: would this be wet or dry ratio?)
 #phi <- 0     # penalty set to 0 to start off with, can change later
  
  #Input intake
  IV <- input["Protein_intake"]
  IH <- input["Lipid_intake"]

  
  if (is.null(betaV)) betaV<-input["betaV"]     
  if (is.null(betaH)) betaH <- input["betaH"] 
  if (is.null(kstarN)) kstarN<-input["kstarN"]   
  if (is.null(phi)) phi <- input["phi"]

  #IV <- input[,"IV"]        # intake protein, mol C mol C-1 d-1
  #IH <- input[,"IH"]
  #if (is.null(betaV)) betaV<-input[3]     
  #if (is.null(betaH)) betaH <- input[4] 
  #if (is.null(kstarN)) kstarN<-input[5]   
  #if (is.null(phi)) phi <- input[6]
  
  #attempting to check if the function is being passed on correctly
  print(kstarN)
  print(phi)
  # intake carbohydrate, mol C mol C-1 d-1
  # Step 1 Utilisation of carbohydrate
  IVm <- thetaV*tau/(betaV*kstarN*thetaZ)  # Eq. S1
  DXm <- (thetaZ-thetaV)/thetaZ*tau        # Eq. S2
  IstarHm <- (eta*IVm+DXm+xi)/(betaH-eta)  # Eq. S3
  
  InVG <- thetaV/(betaV*kstarN*thetaZ)     # Eq. S4
  DnXG <- (thetaZ-thetaV)/thetaZ           # Eq. S5
  InstarHG <- (eta*InVG+DnXG)/(betaH-eta)  # Eq. S6
  
  thetastarHVG <- InstarHG/InVG            # Eq. S7
  thetastarHV <- (IstarHm+(IV-IVm)*thetastarHVG)/IV   # Eq. S8
  IstarH <- IV*thetastarHV                 # Eq. S9
  IHU <- min(IH,IstarH)                    # Eq. S10
  CX <- (betaH-eta)*(IH-IHU)               # Eq. S11
  
  # Step 2 Calculation of fV and G
  if (IH < IstarH){
    a <- betaV*IV
    b <- thetaV/(kstarN*(thetaZ-thetaV))
    cc <- xi + eta*(IV+IHU)
    d <- betaH*IHU
    fV <- (-a+b*d-b*cc)/(-a+phi*b*d-d-b*cc)  # Eq. S12
  } else {
    fV <- 0
  }
  
   if (fV < 1) {
    DCX <- (betaH*IHU*(1-phi*fV))/(1-fV)-(xi+eta*(IV+IHU))     # Eq. S13
    G <- thetaZ*DCX/(thetaZ-thetaV)-tau      # Eq. S14a
  } else {
    numerator <- betaV*IV-(xi+eta*(IV+IH))/(1-phi)
    denominator <- thetaV/(kstarN*thetaZ)+(thetaZ-thetaV)/((1-phi)*thetaZ)
    G <- numerator/denominator-tau           # Eq. S14b
  }
  
  
  # Calculate Omega
  DX <- (thetaZ-thetaV)/thetaZ*(tau+G)       # Eq. 3
  DA <- DX + xi + eta*(IV+IH)                # Eq. 4
  if (phi > 0) {
    Omega <- DA*(1/(1-phi*fV)-1)             # Eq. 5
  } else {
    Omega <- 0.
  }
  
  
  # Step 3 C and N budgets
  DV <- thetaV/thetaZ*(tau+G)/kstarN         # Eq. 2
  R <- tau+(1-kstarN)*DV+xi+eta*(IV+IH)+Omega+CX     # Eq. S15
  E <- tau/thetaZ+(1-kstarN)*DV/thetaV+fV/thetaV*(DX+xi+eta*(IV+IH)+Omega)
  # Eq. S16
  WC <- (1-betaV)*IV+(1-betaH)*IH            # Eq. S17
  WN <- (1-betaV)*IV/thetaV                  # Eq. S18
  
  IC <- IV+IH       # C intake, mol C mol C-1 d-1
  IN <- IV/thetaV   # N intake, mol N mol C-1 d-1
  Ctot <- G+R+WC    # C budget check
  GN <- G/thetaZ    # growth, mol N mol C-1 d-1
  Ntot <- GN+E+WN   # N budget check
  
  #Outputs
  
  # print("         ")
  # print(c("Intake V (IV; specified)",IV))
  # print(c("Intake H (IH; specified)",IH))
  # print(c("fV; calculated)",fV))
  # print("         ")
  # print("C budget:")
  # print(c("Intake C (IV+IH)",IC))
  # print(c("Growth (G)",G))
  # print(c("Respiration (R)",R))
  # print(c("Pellet C (WC)",WC))
  # print(c("C sinks (G+R+WC)",Ctot))
  # print("         ")
  # print("N budget:")
  # print(c("Intake N (IV+IH)",IN))
  # print(c("Growth (GN)",GN))
  # print(c("Excretion (E)",E))
  # print(c("Pellet N (WN)",WC))
  # print(c("N sinks (G+R+WC)",Ntot))
  
return(data.frame(growth=G,Respiration=R,Excretion=E,fV=fV,protein_absorption_efficiency=betaV,lipid_absorption_efficiency=betaH, protein_absorption=betaV*IV, lipid_absorption=betaH*IH, Pellet_C=WC, Pellet_N=WN, Growth_G = G, Growth_N=GN))
  
} #end of function


##  Error function to compare model Iv and Ih predictions with observations
getError<-function(val=c(0.3,0), input=input,feed="BM0%",paramname="both"){
  
finput<-data.frame(G=input[input$Feed==feed,"G"],fV=input[input$Feed==feed,"fV"],betaV=input[input$Feed==feed,"betaV"],betaH=input[input$Feed==feed,"betaH"], Protein_intake=input[input$Feed==feed,"Protein_intake"],Lipid_intake=input[input$Feed==feed,"Lipid_intake"])
  
  if (paramname=="kstarN") out<-GSreversemodel(finput,kstarN = val)
  if (paramname=="phi") out<-GSreversemodel(finput,phi = val)
  if (paramname=="both") out<-GSreversemodel(finput,kstarN=val[1],phi = val[2])
  if (paramname=="withbetas") out<-GSreversemodel(finput,kstarN=val[1],phi = val[2], betaV=val[3],betaH=val[4])
  
   ## calculate and output error, convert observed from tonnes to grams (per m2 per yr)
  squared_error <- c((out$Protein_intake- finput$Protein_intake)^2,(out$Lipid_intake- finput$Lipid_intake)^2)
  
  rmse<-sqrt(sum(squared_error,na.rm=T)/sum(!is.na(squared_error)))
  
  return(rmse)
  
}

