############################################################################
#
# Covid simulation code
#
# This is a discrete-time, age stratified, age-of-infection, compartmental model for COVID
# Data is customized to Israel and updated to the end of December 2020
#
# Main details of the model are described in 
#
# Katriel, "Stochastic discrete-time age-of-infection epidemic models." International Journal of Biomathematics 6, no. 01 (2013): 1250066.
#
# Yaari R, Katriel G, Stone L, Mendelson E, Mandelboim M, Huppert A. 2016 Model-based reconstruction of an epidemic using multiple datasets: understanding influenza A/H1N1 pandemic dynamics in Israel. 
# J. R. Soc. Interface 13: 20160099. http://dx.doi.org/10.1098/rsif.2016.0099
#
# Developed by Nir Gavish (Technion), Amit Huppert (Gertner Institute), Haggai Katriel (Ord Braude), and Rami Yaari (Gertner Institute)
# 
############################################################################

# This function prepares the inputs for the simulations and runs it
main <- function() {
  ####################################################### 
  # Parameters, initial conditions and other definitions
  ####################################################### 
  
  # Gather parameters & extract out common parameters to avoid using parms$
  parms=gatherParameters();
  m<-parms$m # Number of age groups
  N<-parms$N # Size of age groups
  ag_dist<-parms$ag_dist # Size of age groups
  d<-parms$d # effective gtd support
  Tmax<-parms$Tmax # Simulation time
  
  # Initial condition
  Reff0 <- 1.25 # Effective reproduction number in the beginning of the simulation
  S0 <- 0.9*parms$N # Assuming 10% of the population has already recovered
  importedCases <- matrix(0,nrow=Tmax,ncol=m) # Initialized here with no imported cases
  
  # Initialize infected and detected cases in the last 21 days using actual data of detected cases in the last 21 days 
  i0 <- read.csv("i0matrix_init.csv", header = TRUE,colClasses='numeric'); 
  i0 <- round(sapply(1:m,function(i) i0[,i]/parms$det_rates[i])) # Need to fix to include detection rates!
  i0 <- i0*Reff0 #compensating for infection-detection lag when using detected to intialize infected 
  i0 <- win_smooth_mat(i0,7)
  d0 <- nrow(i0) 
  
  # Initialize severe cases using actual data
  sev0 <- as.matrix(read.csv("sev0matrix_init.csv", header = TRUE,colClasses='numeric'))
  
  ## New vaccinations in the last 28 days 
  v_ag_dist <- c(0,0,0,1753,6320,6320,11107,11107,15871,15871,19061,19061,73023,73023,57465,102413)/parms$N
  vacHist <- c(rep(0,21),7.7,23.7,42.6,60.5,74.1,38.7,32) %*% t(v_ag_dist) 
  vacHist <- vacHist*1e3;
  
  ## Vaccinated susceptibles
  susV20 <- rep(0,m)
  
  ################################################## 
  # Contact matrices
  ################################################## 
  # We use POLYMOD-based contact matrices adapted to Israel, see https://doi.org/10.1371/journal.pcbi.1005697 for more details.
  # Additional treatment of ours accounts to changes during the spread of the COVID pandemic.
  C <- getContactMatrices(parms$ag_dist)
  beta <- C$home+0.33*C$work+0.20*C$school+C$other
  beta[15,] <- beta[15,]*0.75 # Decreased contacts of 70-74 ag yields far better agreement with actual data
  
  ngm <- t(t(beta)*parms$ag_dist) # Computing the next generation matrix (ngm) from the contact matrix
  ngm <- ngm*parms$susceptibilityFactor # Accounting for decreased susceptibility & infectivity in children
  ngm <- t(t(ngm)*parms$infectivityFactor)
  # Normalized by initial condition so that Reff0*ngm corresponds to effective reproduction number Reff0
  ngm_polymod <-ngm/max(Re(eigen(t(ngm)*S0/N)$values)) 
  
  ################################################## 
  # Vaccination assumptions and strategies
  ################################################## 
  Tv1 <- 7 # Duration of phase 0 
  Tv2 <- 28-Tv1 # Duration of phase 1
  eff1 <- 0.9 # Vaccine efficacy to preventing disease during phase 1 (more than Tv1 days and before it reaches full effficacy)
  eff2 <- 0.95 # Vaccine efficacy to preventing disease during phase 2 (after in reaches full efficacy Tv1+Tv2 days after first dose)
  Seff1 <- 0.36 # Vaccine efficacy to preventing infection during phase 1
  Seff2 <- 0.55 # Vaccine efficacy to preventing infection during phase 2
  vaccineEfficacyAssumption = c(Tv1=Tv1,Tv2=Tv2,eff1=eff1,eff2=eff2,Seff1=Seff1,Seff2=Seff2)
  
  # Vaccination schedule during simulation time
  vaccineRate <- 150000
  vaccinationSchedule <- matrix(0,m,Tmax+d0);
  
  ag_dist_vaccine<-0*parms$N;
  ag_dist_vaccine[13:16]<-parms$N[13:16];
  ag_dist_vaccine<-ag_dist_vaccine/sum(ag_dist_vaccine)
  vaccine_distribution_matrix<-matrix(rep(ag_dist_vaccine,7),nrow=m);
  
  vaccinationSchedule[,(d0+1):(d0+7)] <- vaccineRate*vaccine_distribution_matrix # 40000 every day
  
  ##################
  # Run simulation
  ##################
  simulationResults <- run_ag_aoi_ex_model(parms,ngm_polymod,1.25,S0,i0,vacHist,susV20,sev0,importedCases,vaccinationSchedule,vaccineEfficacyAssumption) 
  
  # An example of a simple plot of the detected cases over time
  plot(rowSums(simulationResults$detected))
}

#########################################################
#
# This is the function running the simulation itself.  
# 
#########################################################
run_ag_aoi_ex_model <- function(parms,ngm,R0,S0,i0,vacHist,sus_V20,sev0,importedCases,vaccinationSchedule,vaccineEfficacyAssumption,ssn=0) {
  
  # Extracting common parameters to avoid using parm$ in code
  Tmax<-parms$Tmax # Simulation time
  m<-parms$m # Number of age groups
  N<-parms$N # Size of age groups 
  d<-parms$d # effective gtd (generation time distribution) support
  
  # Extracting parameters regarding assumptions of vaccine efficacy to avoid using vaccineEfficacyAssumption$ in code
  Tv1<-vaccineEfficacyAssumption[['Tv1']]
  Tv2<-vaccineEfficacyAssumption[['Tv2']]
  eff1<-vaccineEfficacyAssumption[['eff1']]
  eff2<-vaccineEfficacyAssumption[['eff2']]
  Seff1<-vaccineEfficacyAssumption[['Seff1']]
  Seff2<-vaccineEfficacyAssumption[['Seff2']]
  
  # Verification of parameters
  d0 <- nrow(i0)
  d1 <- length(parms$Pdet)
  d2 <- length(parms$Phosp)
  d3 <- length(parms$Pmodp)
  d4 <- length(parms$Psevp)
  d5 <- length(parms$Presp)
  d6 <- length(parms$Pdead)
  dv <- Tv1+Tv2
  
  stopifnot(d<=d0,d<=d1,d<=d2,d<=d3,d<=d4,d<=d5,d<=d6)
  stopifnot(length(parms$det_rates)==m,length(parms$hosp_rates)==m,length(parms$modp_rates)==m,
            length(parms$sevp_rates)==m,length(parms$resp_rates)==m,length(parms$death_rates)==m)
  stopifnot(is.matrix(i0),ncol(i0)==m)
  stopifnot(length(S0)==m)
  stopifnot(nrow(importedCases)==Tmax,ncol(importedCases)==m)
  stopifnot(is.matrix(ngm),nrow(ngm)==m,ncol(ngm)==m)

  # Compartments
  sus <- matrix(0,Tmax+d0,m)             # Susceptible (S)
  infected <- matrix(0,Tmax+d0,m)        # newly infected (i)
  sus_V2 <- matrix(0,Tmax+d0,m)          # Vaccinated at least Tv1+Tv2 days after vaccination (vaccine reached it's full efficacy)
  v0 <- matrix(0,Tmax+dv,m)              # New vaccinated (v0) 'susceptibles' 
  v0_infected <- matrix(0,Tmax+d0,m)     # Newly infected less than Tv1 days after vaccinated
  v1_infected <- matrix(0,Tmax+d0,m)     # Newly infected less than Tv1+Tv2 days after vaccinated
  v2_infected <- matrix(0,Tmax+d0,m)     # Newly infected at least Tv1+Tv2 days after vaccinated
  hosp <- matrix(0,Tmax+d0,m)            # Hospitalized
  modp <- matrix(0,Tmax+d0,m)            # Moderate condition
  sevp <- matrix(0,Tmax+d0,m)            # Severe condition
  resp <- matrix(0,Tmax+d0,m)            # Need respirator
  dead <- matrix(0,Tmax+d0,m)            # Deceased
  
  detected <- matrix(0,Tmax+d0,m)        # Detected     
  scaled_detected <- matrix(0,Tmax+d0,m) # auxilary variable for estimating morbidity after vaccination

  # Auxilary variables for number of susceptibles at different stages after vaccine administration
  Sv0 <- rep(0,m)
  Sv1 <- rep(0,m)
  
  # Initialization
  infected[1:d0,] <- i0
  v0[1:dv,] <- vacHist
  detected[1:d0,] <- round(sapply(1:m,function(i) i0[,i]*parms$det_rates[i]))
  scaled_detected[1:d0,] <- detected[1:d0,] # Assume in initialization that there are no infected people who have been vaccinated
  sevp[(d0-8):d0,] <- sev0
  sus[d0,] <- S0*(1-colSums(v0)/N) # This initialization does not yet take into account people vaccinated more than Tv1+Tv2 days ago
  v0_infected[1:d0,] <- 0*i0
  v1_infected[1:d0,] <- 0*i0
  v2_infected[1:d0,] <- 0*i0
  ngmt <- R0*ngm
  
  # Pre-allocate additional variables
  Reff <- rep(0,Tmax)
  
  # Main loop marching in time
  for(t in (d0+1):(Tmax+d0)) {
    t1 <- max(1,(t-d1+1)):t
    t2 <- max(1,(t-d2+1)):t
    t3 <- max(1,(t-d3+1)):t
    t4 <- max(1,(t-d4+1)):t
    t5 <- max(1,(t-d5+1)):t
    t6 <- max(1,(t-d6+1)):t
    
    # For each age group j
    for(j in 1:m) {
      # Compute number of susceptible vaccinated individuals in phase 0,1
      Sv0[j]=sum(v0[(t+dv-d0-1):(t+dv-d0-Tv1),j]) # sum of all those that were vaccinated in the last Tv1 days
      Sv1[j]=sum(v0[(t+dv-d0-1-Tv1):(t+dv-d0-Tv2-Tv1),j]) # sum of all those that were vaccinated Tv1->(Tv1+Tv2) days ago
      
      # Force of infection   
      # Here we assume that the everyone have the same infectivity (symptomatic,a-symptomatic,vaccinated at all stages)
      Finf=sum(sapply(1:m, function(k) ngmt[j,k]*sum(parms$gtd[d:1]*(infected[(t-d):(t-1),k]+v0_infected[(t-d):(t-1),k]+v1_infected[(t-d):(t-1),k]+v2_infected[(t-d):(t-1),k]))))/N[j]
      
      # Infections 
      infected[t,j] <- sus[t-1,j]*Finf+importedCases[t-d0,j]
      v0_infected[t,j] <- Sv0[j]*Finf
      v1_infected[t,j] <- (1-Seff1)*Sv1[j]*Finf
      v2_infected[t,j] <- (1-Seff2)*sus_V2[t-1,j]*Finf 
      
      # Disease dynamics
      detected[t,j] <- sum(parms$det_rates[j]*parms$Pdet[length(t1):1]*infected[t1,j])
      # scaled detected is an auxilary variable to account for the efficacy of the vaccine in preventing disease (but perhaps not infection)
      scaled_detected[t,j] <- sum(parms$det_rates[j]*parms$Pdet[length(t1):1]*(infected[t1,j]+v0_infected[t1,j]+(1-eff1)*v1_infected[t1,j]+(1-eff2)*v2_infected[t1,j])) 
      hosp[t,j] <- sum(parms$hosp_rates[j]*parms$Phosp[length(t2):1]*scaled_detected[t2,j])
      modp[t,j] <- sum(parms$modp_rates[j]*parms$Pmodp[length(t3):1]*scaled_detected[t3,j])
      sevp[t,j] <- sum(parms$sevp_rates[j]*parms$Psevp[length(t4):1]*scaled_detected[t4,j])
      resp[t,j] <- sum(parms$resp_rates[j]*parms$Presp[length(t5):1]*scaled_detected[t5,j])
      dead[t,j] <- sum(parms$death_rates[j]*parms$Pdead[length(t6):1]*scaled_detected[t6,j])
      
      # Book-keeping
      sus[t,j] <- sus[t-1,j]-infected[t,j]
      
      # Remove infected vaccinated from count of historically vaccinated assuming number of infected is proportional to the number of vaccinated every day
      if (Sv0[j]>0) {v0[(t+dv-d0-1):(t+dv-d0-Tv1),j] <- v0[(t+dv-d0-1):(t+dv-d0-Tv1),j]*(1-v0_infected[t,j]/Sv0[j])}
      if (Sv1[j]>0) {v0[(t+dv-d0-1-Tv1):(t+dv-d0-Tv2-Tv1),j] <- v0[(t+dv-d0-1-Tv1):(t+dv-d0-Tv2-Tv1),j]*(1-v1_infected[t,j]/Sv1[j])}
      
      # New susceptible vaccinated 
      sus_V2[t,j] <- sus_V2[t-1,j]+v0[t+dv-d0-Tv1-Tv2,j]-v2_infected[t,j]
      
      # Administrate vaccinations
      v0[t+dv-d0,j] <- vaccinationSchedule[j,t] # vac counts newly vaccinated
      sus[t,j] <- sus[t,j]-v0[t+dv-d0,j] # Here we assume that all vaccines go to susceptibles
    }
    # Compute the effective reproduction number
    sus_aux <- (sus[t-1,]+Sv0+(1-Seff1)*Sv1+(1-Seff2)*sus_V2[t-1,])/N
    Reff[t-d0] <- max(Re(eigen(t(ngmt)*sus_aux)$values))
  }
  
  # Wrap up the simulation data 
  inf <- as.matrix(infected[(d0+1):(Tmax+d0),])
  det <- as.matrix(detected[(d0+1):(Tmax+d0),])
  hosp <- as.matrix(hosp[(d0+1):(Tmax+d0),])
  modp <- as.matrix(modp[(d0+1):(Tmax+d0),])
  sevp <- as.matrix(sevp[(d0+1):(Tmax+d0),])
  resp <- as.matrix(resp[(d0+1):(Tmax+d0),])
  dead <- as.matrix(dead[(d0+1):(Tmax+d0),])
  sus <- as.matrix(sus[(d0+1):(Tmax+d0),])
  vac <- as.matrix(sus_V2[(d0+1):(Tmax+d0),])
  
  # Compute number of overall patients at severe/critical conditions hospitalized at a given time under the assumption that a patient resides on average 9 days in such a situtation
  # This fits well the statistics in the last few months, but will likely decrease toward 7-8 days as the median age of severe cases become younger
  sevpAux <- rowSums(sevp)
  sevpAccumulated <- sapply(1:Tmax, function(k) sum(sevpAux[(d0+k-8):(d0+k)]))
  
  return (list(infected=inf,detected=det,hosp=hosp,modp=modp,sevp=sevp,resp=resp,dead=dead,sus=sus,
               Reff=Reff,vac=vac,sevpAccumulated=sevpAccumulated))
}

#########################################################
#
# Gathering parameters
#
#########################################################
gatherParameters <- function() {
  
  ######################
  # Gather parameters
  ######################
  
  # Generation time distribution 
  gtd <- load_gtd() # generation time distribution
  d <- length(gtd)
  
  ## Demography
  # 16 age groups (compatible with POLYMOD matrices) 
  ag_dist <- c(0.0994,0.0948,0.0846,0.0771,0.0719,0.0665,0.0660,0.0652,
               0.0630,0.0566,0.0478,0.0437,0.0413,0.0411,0.0318,0.0492)
  ag_labels <- c('0-4','5-9','10-14','15-19','20-24','25-29','30-34','35-39',
                 '40-44','45-49','50-54','55-59','60-64','65-69','70-74','75+')
  
  N <- 9e6*ag_dist
  m <- length(N)
  
  # Simulation time
  Tmax <- 100
  
  # Relative probability of different events (detection,hospitalization,..) from day of infection to recovery
  Pdet <- c(0.042,0.061,0.083,0.105,0.120,0.121,0.108,0.089,0.070,0.054,
            0.041,0.031,0.023,0.017,0.012,0.008,0.006,0.004,0.003,0.001,0.001)
  
  Phosp <- c(0.452,0.097,0.055,0.051,0.050,0.050,0.050,0.047,0.039,0.032,0.025,0.018,0.015,0.012,0.007,
             0,0,0,0,0,0,0)
  
  Pmodp <- c(0.171,0.141,0.106,0.068,0.057,0.053,0.051,0.050,0.044,0.043,0.039,0.035,0.030,0.024,0.020,
             0.015,0.013,0.010,0.009,0.007,0.007,0.007)
  
  Psevp <- c(0.187,0.156,0.073,0.072,0.063,0.063,0.063,0.061,0.053,0.049,0.039,0.028,0.018,0.015,0.012,
             0.010,0.009,0.007,0.007,0.006,0.005,0.004)
  
  Presp <- c(0.086,0.069,0.057,0.057,0.057,0.069,0.069,0.071,0.066,0.063,0.054,0.041,0.036,0.030,0.029,
             0.023,0.022,0.017,0.014,0.010,0.009,0.008,0.008,0.007,0.007,0.007,0.005,0.005,0.004)
  
  Pdead <- c(0.031,0.033,0.036,0.040,0.043,0.046,0.047,0.049,0.050,0.049,0.048,0.046,0.043,0.040,0.037,
             0.034,0.030,0.027,0.025,0.023,0.020,0.019,0.017,0.016,0.014,0.013,0.012,0.011,0.010,0.010,
             0.009,0.009,0.009,0.008,0.008,0.007,0.007,0.006,0.005,0.005,0.004,0.004)
  
  # Rate of occurence of different events (detection,hospitalization,..) per age group
  det_rates <- c(rep(0.5,m-1),0.7) # Higher detection rate in age group 75+ in agreement with data
  
  hosp_rates <- c(0.021,0.004,0.004,0.007,0.018,0.026,0.032,0.036,0.042,
                  0.051,0.068,0.092,0.128,0.180,0.248,0.420)
  
  modp_rates <- c(0,0,0,0.001,0.003,0.005,0.009,0.014,0.019,0.029,0.042,
                  0.063,0.090,0.134,0.188,0.348)
  
  sevp_rates <- c(0,0,0,0.001,0.001,0.003,0.005,0.008,0.010,0.018,0.028,
                  0.042,0.063,0.100,0.147,0.287)
  
  resp_rates <- c(0,0,0,0,0,0,0.001,0.001,0.001,0.003,0.005,0.007,0.014,0.027,0.035,0.057)
  
  death_rates <- c(0,0,0,0,0,0,0,0,0,0.001,0.003,0.006,0.012,0.024,0.049,0.164)
  
  # Kids (age-dependent susceptibility and infectivity)
  susceptibilityFactor <- c(0.7,0.7,0.7,1,1,1,1,1,1,1,1,1,1,1,1,1) # To account for reduced susceptiblity in children 
  infectivityFactor <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)  # To account for reduced infectivity in children
  
  return(list(gtd=gtd,d=d,ag_dist=ag_dist,ag_labels=ag_labels,N=N,m=m,Tmax=Tmax,
              Pdet=Pdet,Phosp=Phosp,Pmodp=Pmodp,Psevp=Psevp,Presp=Presp,Pdead=Pdead,
              det_rates=det_rates,hosp_rates=hosp_rates,modp_rates=modp_rates,sevp_rates=sevp_rates,
              resp_rates=resp_rates,death_rates=death_rates,susceptibilityFactor=susceptibilityFactor,
              infectivityFactor=infectivityFactor))
}

load_gtd <- function(mean_P=4.5, sd_P=3.5) {
  x <- 1:21
  shape_P <- (mean_P/sd_P)^2
  scale_P <- mean_P/shape_P 
  P <- pgamma(x,shape=shape_P,scale=scale_P)-pgamma(x-1,shape=shape_P,scale=scale_P)
  P <- P/sum(P)
  d <- length(P)
  mean_gtd <- sum((1:d-0.5)*P)
  sd_gtd <- sqrt(sum((1:d-0.5)^2*P)-mean_gtd^2)
  return(P)
}

getContactMatrices <- function(ag_dist) {
  # Read matrix data from file
  data_dir <- getwd();
  home <- read.csv(paste(data_dir,"/home.csv",sep=""), header = FALSE);
  other <- read.csv(paste(data_dir,"/other.csv",sep=""), header = FALSE);
  school <- read.csv(paste(data_dir,"/school.csv",sep=""), header = FALSE);
  work <- read.csv(paste(data_dir,"/work.csv",sep=""), header = FALSE);
  
  # Age distribution corresponding to the contact matrices loaded
  NiCij <- c(0.102383712,0.0948894219,0.08539818,0.0784777,0.071752378,0.069193072,
             0.067747007448,0.064545005,0.062525105,0.05431925,0.047020072,0.044736207866136,
             0.042452342969943,0.039181481183937,0.026924357017433,0.0484546)
  
  # Normalize according to updated age distributions
  work <- t(t(work)/NiCij*ag_dist)
  home <- t(t(home)/NiCij*ag_dist)
  other <- t(t(other)/NiCij*ag_dist)
  school <- t(t(school)/NiCij*ag_dist)
  
  return(list("home"=home,"work"=work,"other"=other,"school"=school));
}

win_smooth <- function(dat, window, smooth_right=T) {
  
  i1 <- (window-1)/2+1
  l <- length(dat)
  dats <- dat
  fact <- sum(dats[(l-window+1):l])/sum(dats[(l-window):(l-1)])
  dats[(l+1):(l+i1-1)] <- dats[l]*(fact^(1:(i1-1)))
  dats[i1:l] <- zoo::rollmean(dats,window)
  dats <- dats[1:l]
  if(!smooth_right) {
    i2 <- l-i1
    dats[i2:l] <- dat[i2:l]
  }
  return (dats)
}

win_smooth_mat <- function(dat, window, smooth_right=T) {
  
  i1 <- (window-1)/2+1
  l <- nrow(dat)
  dats <- dat
  for(j in 1:ncol(dat)) {
    datj <- dats[,j]
    fact <- sum(datj[(l-window+1):l])/sum(datj[(l-window):(l-1)])
    datj[(l+1):(l+i1-1)] <- datj[l]*(fact^(1:(i1-1)))
    dats[i1:l,j] <- zoo::rollmean(datj,window)
    if(!smooth_right) {
      i2 <- l-i1
      dats[i2:l,j] <- dat[i2:l,j]
    }
  }
  
  return (dats)
}
library(ggtext)
library(ggpubr)
library(reshape2)
main()
