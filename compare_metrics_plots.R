# This script compares the solution obtain from the max param models of the 
# old and new (addition of gut lumen) structures

# This script creates plots for the maximum and GA problems



########################
#======================
#  ***  FUNCTIONS  ***
#======================
#######################

#=====================
#1. Create parameters  
#=====================
create.params <- function(comp_names, w){
  
  # List with names of all possible compartments
  all_comps <- list("RoB"="RoB","Heart"="Heart", "Kidneys"="Kidneys", "Brain"="Brain", "Spleen"="Spleen",
                    "Lungs"="Lungs", "Liver"="Liver", "Uterus"="Uterus", "Bone"="Bone", "Adipose"="Adipose", "Skin"="Skin", "Muscles"="Muscles",
                    "GIT"="GIT") # List with names of all possible compartments
  
  ### Density of tissues/organs
  d_tissue <- 1 #g/ml
  d_skeleton <- 1.92 #g/ml
  d_adipose <- 0.940 #g/ml
  
  Q_total <- (1.54*w^0.75)*60 # Total Cardiac Output (ml/h)
  
  Total_Blood <- 0.06*w+0.77 # Total blood volume (ml)
  
  fr_ad <- 0.0199*w + 1.644 # w in g,  Brown et al.1997 p.420. This equation gives the  adipose % of body weight 
  
  #read data from excel
  fractions <- openxlsx::read.xlsx("Rat physiological parameters.xlsx", sheet = 1, colNames = T, rowNames = T)
  fractions <- as.matrix(sapply(fractions, as.numeric))
  rownames(fractions) <- all_comps
  
  #Tissue weight fraction 
  Tissue_fractions <- fractions[,1]/100 # % of BW. Na values refers to the volume of the rest organs(RoB)
  Tissue_fractions[10] <- fr_ad/100
  #Regional blood flow fraction
  Regional_flow_fractions <- fractions[,2]/100 # % of total cardiac output
  #Capillary volume fractions (fractions of tissue volume)
  Capillary_fractions <- fractions[,3] # of tissue volume
  #Macrophage content as fraction tissue volume for each tissue/organ
  Macrophage_fractions <- fractions[,4] 
  
  W_tis <- rep(0,length(comp_names))
  V_tis <- rep(0,length(comp_names))
  V_cap <- rep(0,length(comp_names))
  W_macro <- rep(0,length(comp_names))  #one more for blood compartment
  Q <- rep(0,length(comp_names))
  
  
  for (i in 1:length(comp_names)) {
    control <- comp_names[i]
    
    Tissue_fractions[i] <- ifelse(is.na(control), NA, Tissue_fractions[i])
    Regional_flow_fractions[i] <- ifelse(is.na(control), NA, Regional_flow_fractions[i])
    Capillary_fractions[i] <- ifelse(is.na(control), NA, Capillary_fractions[i])
    Macrophage_fractions[i] <- ifelse(is.na(control), NA, Macrophage_fractions[i])
    
    ### Calculation of tissue weights  
    W_tis[i] <- w*Tissue_fractions[i]
    
    
    ###Calculation of tissue volumes
    
    if (i==9){
      V_tis[i] <- W_tis[i]/d_skeleton
    } else if(i==10){
      V_tis[i] <- W_tis[i]/d_adipose
    } else{
      V_tis[i] <- W_tis[i]/d_tissue 
    }
    
    ###Calculation of capillary volumes
    V_cap[i] <- V_tis[i]*Capillary_fractions[i]
    
    ###Volume of macrophage contents
    W_macro[i] <- W_tis[i]*Macrophage_fractions[i]
    
    ###Calculation of regional blood flows
    Q[i] <- Q_total*Regional_flow_fractions[i]
  }
  
  #Vm_ven <- 0.01*Vven #macrophage content in veins
  #Vm_art <- 0.01*Vart #0.02*Vart #macrophage content in arteries
  
  ### Calculations for "Soft tissue" compartment
  W_tis[1] <- w - sum(W_tis[2:length(W_tis)], na.rm = TRUE)
  V_tis[1] <- W_tis[1]/d_adipose     
  Q[1] <- Q_total - sum(Q[2:length(Q)],na.rm = TRUE) + Q[6]
  V_cap[1] <- V_tis[1]*Capillary_fractions[1] #Total_Blood - Vven - Vart - sum(V_cap[2:length(V_cap)], na.rm = TRUE)
  W_macro[1] <- W_tis[1]*Macrophage_fractions[1]
  
  parameters <- matrix(c(W_tis[],V_tis[],V_cap[],Q[],W_macro[]), ncol = 5)
  colnames(parameters) <- c("W_tis", "V_tis", "V_cap", "Q", "W_macro")
  rownames(parameters) <- all_comps
  
  Vven=0.64*Total_Blood
  Vart=0.15*Total_Blood
  Wm_ven=0.01*Vven
  Wm_art=0.01*Vart
  
  return(c(
    "Q_total"=Q_total, "V_blood"=Total_Blood, "V_ven"=Vven, "V_art"=Vart,
    
    "w_rob"=parameters[1,1], "w_ht"=parameters[2,1], "w_ki"=parameters[3,1], "w_spl"=parameters[5,1], "w_lu"=parameters[6,1], "w_li"=parameters[7,1], "w_bone"=parameters[9,1], "w_git"=parameters[13,1],
    
    "V_tis_rob"=parameters[1,2], "V_tis_ht"=parameters[2,2], "V_tis_ki"=parameters[3,2], "V_tis_spl"=parameters[5,2], "V_tis_lu"=parameters[6,2], "V_tis_li"=parameters[7,2], "V_tis_bone"=parameters[9,2], "V_tis_git"=parameters[13,2], 
    
    "V_cap_rob"=parameters[1,3], "V_cap_ht"=parameters[2,3], "V_cap_ki"=parameters[3,3], "V_cap_spl"=parameters[5,3], "V_cap_lu"=parameters[6,3], "V_cap_li"=parameters[7,3], "V_cap_bone"=parameters[9,3], "V_cap_git"=parameters[13,3],
    
    "Q_rob"=parameters[1,4], "Q_ht"=parameters[2,4], "Q_ki"=parameters[3,4], "Q_spl"=parameters[5,4], "Q_lu"=parameters[6,4], "Q_li"=parameters[7,4], "Q_bone"=parameters[9,4], "Q_git"=parameters[13,4]
    
  ))
}

# Physiological parameters units
# V_blood, V_ven, V_art (ml): Volume of total blood, venous blood and arterial blood
# w_i (g):                    mass of tissue or organ "i"
# V_tis_i (ml):                volume of tissue or organ "i"
# V_cap_i (ml):                volume of capillary blood in tissue "i"
# Q_i, Q_total (ml/h):        regional blood flow of tissue or organ "i"




#=======================
#2. Create system matrix  
#=======================
#--------------------------------------------------------------------------------------------------
# "create_ODE_matrix()" creates the matrix with the coefficients of the state variables of the 
# desired ODE system. It takes as input the values of parameters and returns the matrix.
#--------------------------------------------------------------------------------------------------

create_ODE_matrix <- function(phys_pars, fit_pars, position){
  with( as.list(phys_pars),{
    #============================
    #Indexing of state variables
    #============================
    # x1 <- M_ven
    # x2 <- M_art
    
    # Capillaries       | Tissue  
    # x3 <- M_cap_ht    | x4 <- M_ht
    # x5 <- M_cap_lu    | x6 <- M_lu
    # x7 <- M_cap_li    | x8 <- M_li
    # x9 <- M_cap_spl   | x10 <- M_spl
    # x11 <- M_cap_ki   | x12 <- M_ki
    # x13 <- M_cap_git  | x14 <- M_git
    # x15 <- M_cap_bone | x16 <- M_bone
    # x17 <- M_cap_rob  | x18 <- M_rob 
    # x19 <- M_feces    | x20 <- M_urine
    # x21 <- M_lumen
    
    # Matrix "A" contains the coefficients of the ODEs system of the PBPK. The ODEs system contains 20 variables, 
    # so the dimensions of matrix A are 20x20. Each row of the matrix represents the differential equation of each 
    # state variable x_i and each column represents the value of the coefficient of each state variable x_j in ODE 
    # of each x_i. The indexing of state variables is analytically presented in the table "Indexing of state variables".
    P_ht <- fit_pars[position[1]]
    P_lu <- fit_pars[position[2]]
    P_li <- fit_pars[position[3]]
    P_spl <- fit_pars[position[4]]
    P_ki <- fit_pars[position[5]]
    P_git <- fit_pars[position[6]]
    P_bone <- fit_pars[position[7]]
    P_rob <- fit_pars[position[8]]
    
    x_ht <- fit_pars[position[9]]
    x_lu <- fit_pars[position[10]]
    x_li <- fit_pars[position[11]]
    x_spl <- fit_pars[position[12]]
    x_ki <- fit_pars[position[13]]
    x_git <- fit_pars[position[14]]
    x_bone <- fit_pars[position[15]]
    x_rob <- fit_pars[position[16]]
    
    CLE_f <- fit_pars[length(fit_pars)-2]
    CLE_u <- fit_pars[length(fit_pars)-1]
    CLE_h <- fit_pars[length(fit_pars)]
    
    A <- matrix(c(rep(0,21^2)), 
                nrow = 21)
    rownames(A) <- c("M_ven", "M_art",
                     "M_cap_ht" ,"M_ht",
                     "M_cap_lu" ,"M_lu",
                     "M_cap_li" ,"M_li",
                     "M_cap_spl" ,"M_spl",
                     "M_cap_ki" ,"M_ki",
                     "M_cap_git" ,"M_git",
                     "M_cap_bone" ,"M_bone",
                     "M_cap_rob" ,"M_rob", 
                     "M_feces"  ,"M_urine",
                     "M_lumen")
    colnames(A) <- rownames(A)
    
    #Venous
    A[1,1]<- -Q_total/V_ven; A[1,3]<-Q_ht/V_cap_ht; A[1,7]<-(Q_spl+Q_li+Q_git)/V_cap_li; A[1,11]<-Q_ki/V_cap_ki;
    A[1,15]<-Q_bone/V_cap_bone; A[1,17]<-Q_rob/V_cap_rob
    
    #Arterial
    A[2,2]<- -Q_total/V_art; A[2,5]<-Q_total/V_cap_lu
    
    #Heart - Capillaries
    A[3,2]<- Q_ht/V_art; A[3,3] <- -Q_ht/V_cap_ht -x_ht*Q_ht/V_cap_ht; A[3,4]<- x_ht*Q_ht/(w_ht*P_ht) 
    #Heart - Tissue
    A[4,3]<- x_ht*Q_ht/V_cap_ht; A[4,4] <- - x_ht*Q_ht/(w_ht*P_ht)
    
    #Lungs- Capillaries
    A[5,1] <- Q_total/V_ven; A[5,5] <- -(Q_total/V_cap_lu + x_lu*Q_total/V_cap_lu); A[5,6] <- x_lu*Q_total/(w_lu*P_lu)
    #Lungs - Tissue
    A[6,5] <- x_lu*Q_total/V_cap_lu; A[6,6] <- -x_lu*Q_total/(w_lu*P_lu)
    
    #Liver - capillaries
    A[7,2] <- Q_li/V_art; A[7,7]<- -Q_li/V_cap_li - Q_spl/V_cap_li - Q_git/V_cap_li  - x_li*Q_li/V_cap_li; 
    A[7,8] <- x_li*Q_li/(w_li*P_li); A[7,9] <- Q_spl/V_cap_spl; A[7,13] <- Q_git/V_cap_git
    #Liver - Tissue
    A[8,7]<-x_li*Q_li/V_cap_li; A[8,8]<- - x_li*Q_li/(w_li*P_li) - CLE_h
    
    #Spleen - Capillaries
    A[9,2] <- Q_spl/V_art; A[9,9]<- -Q_spl/V_cap_spl - x_spl*Q_spl/V_cap_spl; A[9,10] <- x_spl*Q_spl/(w_spl*P_spl)
    #Spleen - Tissue
    A[10,9] <- x_spl*Q_spl/V_cap_spl; A[10,10]<- -x_spl*Q_spl/(w_spl*P_spl)
    
    # Kidneys - Capillaries
    A[11,2] <- Q_ki/V_art; A[11,11] <- -Q_ki/V_cap_ki -x_ki*Q_ki/V_cap_ki - CLE_u; A[11,12] <- x_ki*Q_ki/(w_ki*P_ki)
    #Kidneys -Tissue
    A[12,11] <- x_ki*Q_ki/V_cap_ki; A[12,12] <- - x_ki*Q_ki/(w_ki*P_ki)
    
    #Git - Capillaries
    A[13,2] <- Q_git/V_art; A[13,13] <- - Q_git/V_cap_git - x_git*Q_git/V_cap_git; A[13,14] <- x_git*Q_git/(w_git*P_git)
    #Git - Tissue
    A[14,13] <- x_git*Q_git/V_cap_git; A[14,14] <- - x_git*Q_git/(w_git*P_git) 
    #Git - Lumen
    A[21,8] <- CLE_h; A[21,21] <- - CLE_f
    
    #Bone - Capillaries
    A[15,2] <- Q_bone/V_art; A[15,15]<- -Q_bone/V_cap_bone -x_bone*Q_bone/V_cap_bone; A[15,16] <- x_bone*Q_bone/(w_bone*P_bone)
    #Bone - Tissue
    A[16,15] <- x_bone*Q_bone/V_cap_bone; A[16,16] <- - x_bone*Q_bone/(w_bone*P_bone)
    
    #RoB - Capillaries
    A[17,2] <- Q_rob/V_art; A[17,17] <- - Q_rob/V_cap_rob - x_rob*Q_rob/V_cap_rob; A[17,18] <- x_rob*Q_rob/(w_rob*P_rob)
    #RoB - Tissue
    A[18,17] <- x_rob*Q_rob/V_cap_rob; A[18,18] <- - x_rob*Q_rob/(w_rob*P_rob)
    
    #Feces 
    A[19,21] <- CLE_f
    
    #Urine
    A[20,11] <- CLE_u
    
    return(A)
  })
}

#====================
#3. Matrix exponent 
#====================
#--------------------------------------------------------------------------------------------------
# "Solve_exp_matrix()" is a function that solves the ODE system using the matrix "x" (which 
# contains the coefficients of the system), "time" which is the desired time points to 
# be calculated and "y_init" is the initial values of the state variables.
#--------------------------------------------------------------------------------------------------

solve_exp_matrix <- function(x, time, y_init, phys_pars){
  with( as.list(phys_pars),{
    if(!is.matrix(x)){
      stop("x must be a NxN matrix")
    }
    
    if(!is.numeric(y_init)){
      stop("y_init must be a numeric vector")
    }
    
    if(dim(x)[1] != dim(x)[2]){
      stop("Matrix x must be NxN")
    }
    
    if(dim(x)[1] != length(y_init)){
      stop("Dimension of y_init must be equal to dimension of matrix x")
    }
    
    
    y_t  <- matrix(data=NA, nrow = nrow(x), ncol = length(time))
    colnames(y_t) <- as.character(time)
    
    y_t[,1] <- y_init
    for (t in 2:length(time)) {
      solution_t <- expm::expm(x*time[t])%*%y_init
      y_t[,t] <- solution_t
    }
    rownames(y_t) <- rownames(x)
    
    y_t <- data.frame(t(y_t))
    
    # Transform TiO2 masses to concentrations
    concentrations <- cbind(time,
                            (y_t$M_ven + y_t$M_art + y_t$M_cap_ht + y_t$M_cap_lu +
                               y_t$M_cap_li+y_t$M_cap_spl+
                               y_t$M_cap_ki+ y_t$M_cap_git +
                               y_t$M_cap_bone + y_t$M_cap_rob)/V_blood,
                            
                            y_t$M_ht/w_ht,
                            y_t$M_lu/w_lu,
                            y_t$M_li/w_li,
                            y_t$M_spl/w_spl,
                            y_t$M_ki/w_ki,
                            (y_t$M_git+y_t$M_lumen)/w_git,
                            y_t$M_bone/w_bone,
                            y_t$M_feces,
                            y_t$M_urine)
    colnames(concentrations) <- c("Time","C_blood", "C_ht", "C_lu", "C_li",
                                  "C_spl", "C_ki", "C_git", "C_bone", "Feces", "Urine")
    
    #return(list(y_t, concentrations))
    return(data.frame(concentrations))
  })
}


fitness.metric <- function(observed, predicted, comp.names =NULL){
  # Check if the user provided the correct input format
  if (!is.list(observed) || !is.list(predicted)){
    stop(" The observations and predictions must be lists")
  }
  # Check if the user provided equal length lists
  if (length(observed) != length(predicted)){
    stop(" The observations and predictions must have the same compartments")
  }
  Ncomp <- length(observed) # Number of compartments
  I <- rep(NA, Ncomp) # Compartment discrepancy index
  N_obs <- rep(NA, Ncomp) #Number of observations per compartment
  #loop over the compartments
  for (i in 1:Ncomp){
    Et <- 0 #relative error with observations
    St <- 0  #relative error with simulations
    N <- length(observed[[i]]) # number of observations for compartment i
    # Check if observations and predictions have equal length
    if(N != length(predicted[[i]])){
      stop(paste0("Compartment ",i," had different length in the observations and predictions"))
    }
    N_obs[i] <- N # populate the N_obs vector
    for (j in 1:N){
      # sum of relative squared errors (error = observed - predicted)
      Et <- Et + ( abs(observed[[i]][j] - predicted[[i]][j])  / observed[[i]][j] )  ^2
      St <- St + ( abs(observed[[i]][j] - predicted[[i]][j])  / predicted[[i]][j] )  ^2
    }
    
    # root mean of the square of observed values
    RMEt <- sqrt(Et/N)
    # root mean of the square of simulated values
    RMSt <- sqrt( St/N)
    
    I[i] <- (RMEt + RMSt)/2   
  }
  # Total number of observations
  Ntot <- sum(N_obs)
  # Initialise the consolidated discrepancy index
  Ic <-0
  for (i in 1:Ncomp){
    # Give weight to compartments with more observations (more information)
    Ic <- Ic +  I[i]* N_obs[i]/Ntot
  }
  # Name the list of compartment discrepancy indices
  if ( !is.null(comp.names)){
    names(I) <- comp.names
  }else if (!is.null(names(observed))){
    names(I) <- names(observed)
  } else if (!is.null(names(predicted)) && is.null(comp.names) ){
    names(I) <- names(predicted)
  }
  return(Ic)
  #return(list(Total_index = Ic, Compartment_index= I))
}

# Function to estimate the percentage of percent of model-predicted concentrations
# falling within twofold of the corresponding observed concentrations
two.fold <- function(predictions, observations, times=NULL){
 y_obs <- unlist(observations)
 y_pred <- unlist(predictions)
 # Total number of observations
 N<- length(y_obs)
 # Counter for counting how many observations lie within two fold from the data
 counter <- 0
 for ( i in 1:N){
   if ((y_pred[i]<=2*y_obs[i]) & (y_pred[i]>=0.5*y_obs[i])){
     counter <- counter + 1
   }
 }
 twofold_percentage <- (counter/N)*100
 return(twofold_percentage)
}

#  absolute average fold error
AAFE <- function(predictions, observations, times=NULL){
  y_obs <- unlist(observations)
  y_pred <- unlist(predictions)
  # Total number of observations
  N<- length(y_obs)
  log_ratio <- rep(NA, N) 
  for ( i in 1:N){
    log_ratio[i] <- abs(log((y_pred[i]/y_obs[i]), base = 10))
  }
  aafe <- 10^(sum(log_ratio)/N) 
  return(aafe)
}

#  R-squared between predictions and observations
r.squared <- function(predictions, observations, times=NULL){
  y_pred <- unlist(predictions)
  y_obs <- unlist(observations)

  lm.model <- lm(y_obs~y_pred)
  r_squared <- summary(lm.model)$r.squared 

  return(r_squared)
}

#  Root-mean-square deviation
rmsd <- function(predictions, observations, times=NULL){
  y_obs <- unlist(observations)
  y_pred <- unlist(predictions)
  # Total number of observations
  N<- length(y_obs)
  summation <- 0
  for ( i in 1:N){
    summation <- summation + (y_obs[i]-y_pred[i])^2
  }
  rmsd <- sqrt(summation/N)
  
  return(rmsd)
}


pbpk.index <- function(observed, predicted, comp.names =NULL){
  # Check if the user provided the correct input format
  if (!is.list(observed) || !is.list(predicted)){
    stop(" The observations and predictions must be lists")
  }
  # Check if the user provided equal length lists
  if (length(observed) != length(predicted)){
    stop(" The observations and predictions must have the same compartments")
  }
  Ncomp <- length(observed) # Number of compartments
  I <- rep(NA, Ncomp) # Compartment discrepancy index
  N_obs <- rep(NA, Ncomp) #Number of observations per compartment
  #loop over the compartments
  for (i in 1:Ncomp){
    et <- 0 # errors
    Et <-0  # experimental
    N <- length(observed[[i]]) # number of observations for compartment i
    # Check if observations and predictions have equal length
    if(N != length(predicted[[i]])){
      stop(paste0("Compartment ",i," had different length in the observations and predictions"))
    }
    N_obs[i] <- N # populate tne N_obs vector
    for (j in 1:N){
      # sum of absolute squared errors (error = observed - predicted)
      et <- et + (abs(observed[[i]][j] - predicted[[i]][j]))^2
      # Sum of squared observed values
      Et <- Et + (observed[[i]][j])^2
    }
    # root mean square of the absolute error
    RMet2 <-sqrt(et/N)
    # root mean of the square of observed values
    RMEt2 <- sqrt(Et/N)
    
    I[i] <- RMet2/RMEt2   
  }
  # Total number of observations
  Ntot <- sum(N_obs)
  # Initialise the consolidated discrepancy index
  Ic <-0
  for (i in 1:Ncomp){
    Ic <- Ic +  I[i]* N_obs[i]/Ntot
  }
  # Name the list of compartment discrepancy indices
  if ( !is.null(comp.names)){
    names(I) <- comp.names
  }else if (!is.null(names(observed))){
    names(I) <- names(observed)
  } else if (!is.null(names(predicted)) && is.null(comp.names) ){
    names(I) <- names(predicted)
  }
  return(Ic)
  #return(list(Total_index = Ic, Compartment_index= I))
}

#======================
#5. Objective function  
#======================

obj.func <- function(params, ...){
  
  dots <- list(...)
  with(as.list(dots),{
    
      # Create the matrix of the system  
      A <- create_ODE_matrix(phys_pars = phys_pars, fit_pars =exp(params),  position = position )
      # Solve the ODE system using the exponential matrix method  
      solution <-  solve_exp_matrix(x = A, time = sample_time, y_init = y_init,phys_pars = phys_pars )
      
    concentrations <- solution[solution$Time %in% time_points, 2:(dim(solution)[2]-2)]
    excr_solution <-  data.frame(solution$Time, solution$Feces, solution$Urine)
    excr_solution <- excr_solution[solution$Time %in% excretion_time_points, c(2:3)]
    
    observed <- list()
    predicted <- list()
    
    for (i in 1:(length(concentrations))) {
      observed[[i]] <- df[,i]
      predicted[[i]] <- concentrations[,i]
    }
    observed[[i+1]] <- excretion[,1] #feces
    observed[[i+2]] <- excretion[,2] #urine
    predicted[[i+1]] <- excr_solution[,1] #feces
    predicted[[i+2]] <- excr_solution[,2] #urine
    
    if(w_version == "r.squared"){
      discrepancy <- - r.squared(observed, predicted)
    }else if(w_version == "PBPK_index"){
      discrepancy <- pbpk.index(observed, predicted)
    }else if(w_version == "proposed_metric"){
      discrepancy <- fitness.metric(observed, predicted)
    }else if(w_version == "two.fold"){
      discrepancy <- - two.fold(observed, predicted)
    }else if(w_version == "AAFE"){
      discrepancy <- AAFE(observed, predicted)
    }else if(w_version == "rmsd"){
      discrepancy <- rmsd(observed, predicted)
    }
    return(discrepancy)
  })
}


#===============
# Load data  
#===============
setwd("C:/Users/ptsir/Documents/GitHub/PBPK_Genetic_Algorithm")

dose_kg <- 10 # mg/kg rat body
mass <- 250 # g  
dose <- dose_kg*mass/1000 # mg TiO2

# Load raw data from paper Xie et al.2011
df <- openxlsx::read.xlsx("TiO2_iv_rat.xlsx", sheet = 1, colNames = T, rowNames = T) # TiO2 NPs %ID/g of tissue  (Table 1)
excretion <- openxlsx::read.xlsx("Cummulative_Excretion.xlsx", sheet = 2, colNames = T, rowNames = F) # accumulated excretory rate, expressed as %ID
excretion_time <- round(excretion[,1])*24 # hours
excretion <- excretion[,c(2:3)]

# Transform to (mg of NPs)/(g of tissue)
df <- (df/100)*dose
df$Intestine <- df$Intestine +df$Stomach
colnames(df)[which(names(df)=="Intestine")] <- "Git"
df <- subset(df, select = -c(Stomach, Brain))

df[5,1] <- 1e-05

excretion <- (excretion/100)*dose

### Important!!! each compartment has a specific index vectors Tissue_fractions, Regional_flow_fractions, Capillary_fractions and cannot be changed
# The index of each compartment:
#Rest of Body (rob) --> 1
#Heart (ht) --> 2
#Kidneys (ki) --> 3
#Brain (br) --> 4
#Spleen (spl) --> 5
#Lungs (lu) --> 6
#Liver (li) --> 7
#Uterus (ut) --> 8
#Bone (bone) --> 9
#Adipose (ad) --> 10
#Skin (skin) --> 11
#Muscles (mu) --> 12
#Gastrointestinal track (GIT) --> 13


#### If any of these compartments don not exist in pbpk, just give it the value NA in compartments vector, example: "Heart" = NA and it will remove it 
#### from the equilibriums and the corresponding V_tis, V_cap, Q will be equal to NA.


compartments <- list( "RoB"="RoB","Heart"="Heart", "Kidneys"="Kidneys", "Brain"= NA, "Spleen"="Spleen",
                      "Lungs"="Lungs", "Liver"="Liver", "Uterus"=NA, "Bone"="Bone", "Adipose"=NA, "Skin"=NA, "Muscles"=NA, "GIT"="GIT") #used as input in function, compartments that are used in pbpk


# Nelder-Mead from dfoptim package
y_init <- c(dose, rep(0,20))
time_points <- c(1,3,7, 15, 30)*24 # hours
excretion_time_points <- excretion_time
sample_time <- seq(0, 30*24, 1)
# Initialise vector of physiological parameters
phys_pars <- create.params(compartments,mass)

#---------------------------
N_p <- 8
N_x <- 8
grouping <- c(1:8,1:8)
# Define size of P and X groups
P_groups <- length(unique(grouping[1:N_p]))  # sample size
X_groups <- length(unique(grouping[(N_p+1):(N_p+N_x)]))  # sample size
# set.seed(0)
# Initilise parameter values
fitted <- rep(NA,P_groups+X_groups+3)
# Initialise naming vectors
pnames <- rep(NA, P_groups)
xnames <- rep(NA, X_groups)

#Define names for P and X groups
for (i in 1:P_groups){
  pnames[i] <- paste0("P", as.character(unique(grouping[1:N_p])[i]))
}
for (j in 1:X_groups){
  xnames[j] <- paste0("X", as.character(unique(grouping[(N_p+1):(N_p+N_x)])[j]))
}
# Define the total parameter vector names
names(fitted) <- c(pnames, xnames,"CLE_f", "CLE_u")
# Variable for keeping which value in the fitted params vector corresponds to each coefficient
position = rep(NA, length(grouping))
for (i in 1:(length(position))){
  if(i<=8){
    position[i] <- which(names(fitted) == paste0("P", as.character(grouping[i])))
  }else{
    position[i] <- which(names(fitted) == paste0("X", as.character(grouping[i])))
  }
}
# Some initialisations fail to obtain solution, so resample until you do
fitted[] <- c(log(exp(runif(P_groups, 3,6))),log(exp(runif(X_groups+3, -3,1))))


MAX <- 800
w_version <- "r.squared"  
# Run the Nelder Mead algorithmm to estimate the parameter values
nm_optimizer_max_r<- dfoptim::nmk(par = fitted, fn = obj.func,
                                     control = list(maxfeval=MAX, trace=T), y_init = y_init,
                                     time_points = time_points,
                                     excretion_time_points =  excretion_time_points,
                                     sample_time = sample_time,
                                     phys_pars = phys_pars, 
                                     position = position )
max_params_r<- exp(nm_optimizer_max_r$par)

w_version <- "PBPK_index"  
# Run the Nelder Mead algorithmm to estimate the parameter values
nm_optimizer_max_pbpk <- dfoptim::nmk(par = fitted, fn = obj.func,
                                      control = list(maxfeval=MAX, trace=T), y_init = y_init,
                                      time_points = time_points,
                                      excretion_time_points =  excretion_time_points,
                                      sample_time = sample_time,
                                      phys_pars = phys_pars, 
                                      position = position )
max_params_pbpk <- exp(nm_optimizer_max_pbpk$par)

w_version <- "proposed_metric"  
# Run the Nelder Mead algorithmm to estimate the parameter values
nm_optimizer_max_new<- dfoptim::nmk(par = fitted, fn = obj.func,
                                       control = list(maxfeval=MAX, trace=T), y_init = y_init,
                                       time_points = time_points,
                                       excretion_time_points =  excretion_time_points,
                                       sample_time = sample_time,
                                       phys_pars = phys_pars, 
                                       position = position )
max_params_new <- exp(nm_optimizer_max_new$par)

w_version <- "two.fold"  
# Run the Nelder Mead algorithmm to estimate the parameter values
nm_optimizer_max_two<- dfoptim::nmk(par = fitted, fn = obj.func,
                                    control = list(maxfeval=MAX, trace=T), y_init = y_init,
                                    time_points = time_points,
                                    excretion_time_points =  excretion_time_points,
                                    sample_time = sample_time,
                                    phys_pars = phys_pars, 
                                    position = position )
max_params_two <- exp(nm_optimizer_max_two$par)

w_version <- "AAFE"  
# Run the Nelder Mead algorithmm to estimate the parameter values
nm_optimizer_max_aafe<- dfoptim::nmk(par = fitted, fn = obj.func,
                                    control = list(maxfeval=MAX, trace=T), y_init = y_init,
                                    time_points = time_points,
                                    excretion_time_points =  excretion_time_points,
                                    sample_time = sample_time,
                                    phys_pars = phys_pars, 
                                    position = position )
max_params_aafe <- exp(nm_optimizer_max_aafe$par)

w_version <- "rmsd"  
# Run the Nelder Mead algorithmm to estimate the parameter values
nm_optimizer_max_rmsd<- dfoptim::nmk(par = fitted, fn = obj.func,
                                     control = list(maxfeval=MAX, trace=T), y_init = y_init,
                                     time_points = time_points,
                                     excretion_time_points =  excretion_time_points,
                                     sample_time = sample_time,
                                     phys_pars = phys_pars, 
                                     position = position )
max_params_rmsd <- exp(nm_optimizer_max_rmsd$par)

# Create the matrix of the system  
A_r <- create_ODE_matrix(phys_pars = phys_pars, fit_pars = max_params_r,  position = position )
A_pbpk <- create_ODE_matrix(phys_pars = phys_pars, fit_pars = max_params_pbpk,  position = position )
A_new <- create_ODE_matrix(phys_pars = phys_pars, fit_pars = max_params_new,  position = position )
A_two <- create_ODE_matrix(phys_pars = phys_pars, fit_pars = max_params_two,  position = position )
A_aafe <- create_ODE_matrix(phys_pars = phys_pars, fit_pars = max_params_aafe,  position = position )
A_rmsd<- create_ODE_matrix(phys_pars = phys_pars, fit_pars = max_params_rmsd,  position = position )


# Solve the ODE system using the exponential matrix method  
solution_r <-  as.data.frame(solve_exp_matrix(x = A_r, time = sample_time, 
                                                    y_init = y_init,phys_pars = phys_pars ))
names(solution_r) <- c("Time","Blood", "Heart", "Lungs", "Liver", "Spleen",
                         "Kidneys","Git", "Bone",  "Feces", "Urine")
solution_pbpk <-  as.data.frame(solve_exp_matrix(x = A_pbpk, time = sample_time, 
                                                        y_init = y_init,phys_pars = phys_pars ))
names(solution_pbpk) <- c("Time","Blood", "Heart", "Lungs", "Liver", "Spleen",
                           "Kidneys","Git", "Bone",  "Feces", "Urine")
solution_new <-  as.data.frame(solve_exp_matrix(x = A_new, time = sample_time, 
                                                          y_init = y_init,phys_pars = phys_pars ))
names(solution_new) <- c("Time","Blood", "Heart", "Lungs", "Liver", "Spleen",
                            "Kidneys","Git", "Bone",  "Feces", "Urine")
solution_two <-  as.data.frame(solve_exp_matrix(x = A_two, time = sample_time, 
                                                y_init = y_init,phys_pars = phys_pars ))
names(solution_two) <- c("Time","Blood", "Heart", "Lungs", "Liver", "Spleen",
                         "Kidneys","Git", "Bone",  "Feces", "Urine")
solution_aafe <-  as.data.frame(solve_exp_matrix(x = A_aafe, time = sample_time, 
                                                y_init = y_init,phys_pars = phys_pars ))
names(solution_aafe) <- c("Time","Blood", "Heart", "Lungs", "Liver", "Spleen",
                         "Kidneys","Git", "Bone",  "Feces", "Urine")
solution_rmsd <-  as.data.frame(solve_exp_matrix(x = A_rmsd, time = sample_time, 
                                                 y_init = y_init,phys_pars = phys_pars ))
names(solution_rmsd) <- c("Time","Blood", "Heart", "Lungs", "Liver", "Spleen",
                          "Kidneys","Git", "Bone",  "Feces", "Urine")


# Create a single data frame to hold the observation data 
observations <- data.frame( Time =c(24,  72, 168, 360, 720), excretion, df)

library(ggplot2)

# Defining the linetype and colour of each curve
ltp <- c("R-squared" = "twodash","RMSD" ="dotted", "PBKOF" = "solid", "AAFE" = "longdash","PBPK index" = "dashed")
cls <-  c("R-squared" = "#56B4E9", "RMSD" = "#E69F00", "PBKOF" ="#000000", "AAFE" = "#009E73", "PBPK index" ="#CC79A7",
          "Observations" = "#D55E00")


create.plots <- function(compartment){  
  excreta <- compartment %in% c("Feces", "Urine")
  ggplot(data = solution_r)+
    geom_line( aes_string(x= "Time", y= rlang::expr(!!compartment), 
                          color = '"R-squared"',linetype = '"R-squared"'),  size=1.5,alpha = 0.7) +
    geom_line(data=solution_pbpk, aes_string(x= "Time", y= rlang::expr(!!compartment),
                                             color = '"PBPK index"',linetype ='"PBPK index"'), size=1.5,alpha = 0.9) +
    geom_line(data=solution_new, aes_string(x= "Time", y= rlang::expr(!!compartment),
                                            color =  '"PBKOF"',linetype =  '"PBKOF"'), size=1.5,alpha = 0.7) +
    geom_line(data=solution_aafe, aes_string(x= "Time", y= rlang::expr(!!compartment), 
                                             color = '"AAFE"',linetype ='"AAFE"'), size=1.5,alpha = 0.7) +
    geom_line(data=solution_rmsd, aes_string(x= "Time", y= rlang::expr(!!compartment), 
                                             color = '"RMSD"',linetype = '"RMSD"'), size=1.5,alpha = 0.7) +
    geom_point(data=observations, aes_string(x="Time", y= rlang::expr(!!compartment), 
                                             color='"Observations"'), size=4)+
    labs(title = rlang::expr(!!compartment), 
         y = ifelse(excreta,"TiO2 (mg)","TiO2 (mg/g tissue)" ),
         x = "Time (hours)")+
    theme(plot.title = element_text(hjust = 0.5))+
  scale_y_continuous(trans='log10')+
    scale_color_manual("",values=cls)+
    scale_linetype_manual("Metrics", values=ltp)
  
}
plots <- lapply(names(observations)[2:length(observations)],create.plots)
p1 <-  plots[[1]]
p2 <-  plots[[2]]
p3 <-  plots[[3]]
p4 <-  plots[[4]]
p5 <-  plots[[5]]
p6 <-  plots[[6]]
p7 <-  plots[[7]]
p8 <-  plots[[8]]
p9 <-  plots[[9]]
p10 <-  plots[[10]]
#gridExtra::grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8, p9,p10, nrow = 4)
#gridExtra::grid.arrange(p5,p6,p7,p8,nrow = 2)
#gridExtra::grid.arrange(p9,p10,nrow = 2)

ggpubr::ggarrange(p1, p2, p3, p4,p5,p6,p7,p8, p9,p10, ncol=3, nrow=4, 
                  common.legend = TRUE, legend="right")

