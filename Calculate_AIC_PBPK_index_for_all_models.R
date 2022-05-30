# This script calculates the AIC and PBPK index of the maximum, minimum and GA problem


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
      
      # Matrix "A" contains the coefficients of the ODEs system of the PBPK. The ODEs system contains 20 variables, 
      # so the dimensions of matrix A are 20x20. Each row of the matrix represents the differential equation of each 
      # state variable x_i and each column represents the value of the coefficient of each state variable x_j in ODE 
      # of each x_i. The indexing of state variables is analytically presented in the table "Indexing of state variables".
      
      # Numbering of parameters
      #1:ht, 2:lu, 3:li, 4:spl, 5:ki, 6:git, 7:bone, 8:rob
      
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
      
      CLE_f <- fit_pars[length(fit_pars)-1]
      CLE_u <- fit_pars[length(fit_pars)]
      
      
      A <- matrix(c(rep(0,20^2)), 
                  nrow = 20)
      rownames(A) <- c("M_ven", "M_art",
                       "M_cap_ht" ,"M_ht",
                       "M_cap_lu" ,"M_lu",
                       "M_cap_li" ,"M_li",
                       "M_cap_spl" ,"M_spl",
                       "M_cap_ki" ,"M_ki",
                       "M_cap_git" ,"M_git",
                       "M_cap_bone" ,"M_bone",
                       "M_cap_rob" ,"M_rob", 
                       "M_feces"  ,"M_urine")
      colnames(A) <- rownames(A)
      
      #Venous
      A[1,1]<- -Q_total/V_ven; A[1,3]<-Q_ht/V_cap_ht; A[1,7]<-(Q_spl+Q_li)/V_cap_li; A[1,11]<-Q_ki/V_cap_ki; A[1,13]<-Q_git/V_cap_git;
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
      A[7,2] <- Q_li/V_art; A[7,7]<- -Q_li/V_cap_li - Q_spl/V_cap_li - x_li*Q_li/V_cap_li; 
      A[7,8] <- x_li*Q_li/(w_li*P_li); A[7,9] <- Q_spl/V_cap_spl
      #Liver - Tissue
      A[8,7]<-x_li*Q_li/V_cap_li; A[8,8]<- - x_li*Q_li/(w_li*P_li)
      
      #Spleen - Capillaries
      A[9,2] <- Q_spl/V_art; A[9,9]<- -Q_spl/V_cap_spl - x_spl*Q_spl/V_cap_spl; A[9,10] <- x_spl*Q_spl/(w_spl*P_spl)
      #Spleen - Tissue
      A[10,9] <- x_spl*Q_spl/V_cap_spl; A[10,10]<- -x_spl*Q_spl/(w_spl*P_spl)
      
      # Kidneys - Capillaries
      A[11,2] <- Q_ki/V_art; A[11,11] <- -Q_ki/V_cap_ki -x_ki*Q_ki/V_cap_ki; A[11,12] <- x_ki*Q_ki/(w_ki*P_ki)
      #Kidneys -Tissue
      A[12,11] <- x_ki*Q_ki/V_cap_ki; A[12,12] <- - x_ki*Q_ki/(w_ki*P_ki) -CLE_u
      
      #Git - Capillaries
      A[13,2] <- Q_git/V_art; A[13,13] <- - Q_git/V_cap_git - x_git*Q_git/V_cap_git; A[13,14] <- x_git*Q_git/(w_git*P_git)
      #Git - Tissue
      A[14,13] <- x_git*Q_git/V_cap_git; A[14,14] <- - x_git*Q_git/(w_git*P_git) - CLE_f
      
      #Bone - Capillaries
      A[15,2] <- Q_bone/V_art; A[15,15]<- -Q_bone/V_cap_bone -x_bone*Q_bone/V_cap_bone; A[15,16] <- x_bone*Q_bone/(w_bone*P_bone)
      #Bone - Tissue
      A[16,15] <- x_bone*Q_bone/V_cap_bone; A[16,16] <- - x_bone*Q_bone/(w_bone*P_bone)
      
      #RoB - Capillaries
      A[17,2] <- Q_rob/V_art; A[17,17] <- - Q_rob/V_cap_rob - x_rob*Q_rob/V_cap_rob; A[17,18] <- x_rob*Q_rob/(w_rob*P_rob)
      #RoB - Tissue
      A[18,17] <- x_rob*Q_rob/V_cap_rob; A[18,18] <- - x_rob*Q_rob/(w_rob*P_rob)
      
      #Feces 
      A[19,14] <- CLE_f
      
      #Urine
      A[20,12] <- CLE_u
      
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
                              (y_t$M_ven + y_t$M_art)/V_blood,
                              y_t$M_ht/w_ht,
                              y_t$M_lu/w_lu,
                              y_t$M_li/w_li,
                              y_t$M_spl/w_spl,
                              y_t$M_ki/w_ki,
                              y_t$M_git/w_git,
                              y_t$M_bone/w_bone,
                              y_t$M_feces,
                              y_t$M_urine)
      colnames(concentrations) <- c("Time","C_blood", "C_ht", "C_lu", "C_li",
                                    "C_spl", "C_ki", "C_git", "C_bone", "Feces", "Urine")
      
      return(data.frame(concentrations))
    })
  }
  
  #===============
  #4. PBPK INDEX  
  #===============
  
  ############# Calculate PBPK indices #############
  # pbpk.index a function the returns the compartment and consolidated (Total) discrepancy index
  # of a PBPK model, given some experimental data. It follows the paper of Krishnan et al.1995.
  # observed: list of vectors containing the experimental data
  # predictions: list of vectors containing the predicted data
  # names of the compartments
  
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
      St <- 0  # simulated
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
        St <- St +  (predicted[[i]][j])^2
      }
      # root mean square of the absolute error
      RMet2 <-sqrt(et/N)
      # root mean of the square of observed values
      RMEt2 <- sqrt(Et/N)
      # root mean of the square of simulated values
      RMSt2 <- sqrt(St/N)
      
      I[i] <- ( (RMet2/RMEt2) + (RMet2/RMSt2))/2   
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
  #===============
  #4. Improved fitness metric  
  #===============
  
  ############# Calculate PBPK indices #############
  # pbpk.index a function the returns the compartment and consolidated (Total) discrepancy index
  # of a PBPK model, given some experimental data. It follows the paper of Krishnan et al.1995.
  # observed: list of vectors containing the experimental data
  # predictions: list of vectors containing the predicted data
  # names of the compartments
  
  
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
      MEt <- Et/N
      # root mean of the square of simulated values
      MSt <- St/N
      
      I[i] <- (MEt + MSt)/2   
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
      
      discrepancy <- fitness.metric(observed, predicted)
      
      return(discrepancy)
    })
  }
  
  #=====================================
  #6. Calculate Residual Sum of Squares  
  #=====================================
  
  RSS <- function(predictions, observations, times=NULL){
    
    if(!is.list(observations)){
      stop("Observations must be given as a list")
    }
    
    for (i in 1:length(observations)) {
      if(!is.matrix(observations[[i]])){
        stop("Each element of observations list must bea 2-column matrix")
      }
      if(ncol(observations[[i]]) != 2){
        stop("Each element of observations must be a 2-column matrix")
      }
    }
    
    # Checking if all the compartments have been measured for the same time points
    for (i in 1:(length(observations)-1)) {
      if(all(observations[[i]][,1] == observations[[i+1]][,1])){
        different_times <- FALSE
      }else{
        # If at least one compartmenthas different time points, break
        different_times <- TRUE
        break
      }
    }
    
    # Set the observation times vector if all observations have the same time points
    if(different_times == FALSE){
      observations_time <- observations[[1]][,1]
    }
    
    # If user provided time points but observations have different time points, then times can't be used  
    if (!is.null(times) & (different_times == TRUE)){
      warning("parameter 'times' will not be used because different time vectors have
            been detected in the observations provided")
    }
    
    # If user provided time points that are not part of the observation time points, then times can't be used  
    if (!is.null(times) & (different_times == FALSE)){
      if (sum(times%in%observations_time)<length(times)){
      warning("parameter 'times' will not be used because it contains time points 
             that are not part of the observations ")
        times <- observations_time
      }
    }
    
    # If user did not provide time points, then keep the time points of the observations
    if (is.null(times) & (different_times == FALSE)){
     times <- observations_time
    }
    
    predicted <- list()
    if(different_times == TRUE){ 
      # if the data time points for each compartment are different, ignore times parameter and keep all the values from the data
      for(i in colnames(predictions)[2:dim(predictions)[2]]){
        predicted[[i]] <- predictions[which(predictions$Time %in% observations[[i]][,1]) ,i]
      }
    }else{ # if the data time points for each compartment are the same, use times parameter and keep data only for these moments
      for(i in colnames(predictions)[2:dim(predictions)[2]]){
        predicted[[i]] <- predictions[which(predictions$Time %in% times),i]
        observations[[i]] <- observations[[i]][which(observations[[i]][,1] %in% times ),]
      }
    }
    
    observed <- list()
    for (i in 1:length(observations)) {
      # If user provided only one point the observations have collapsed into a vector
      if(length(times)==1){
        observed[[i]] <- observations[[i]][2] #drop the column of time for each compartment and keep ony the data
      }else{
        observed[[i]] <- observations[[i]][,2] #drop the column of time for each compartment and keep ony the data
      }
    }
    
    res <- list() 
    for (i in 1:length(observed)) { # loop for each compartment
      res[[i]] <- observed[[i]] - predicted[[i]] # calculate the residuals of each compartment and store them to lists
    }
    names(res) <- names(predictions[2:length(names(predictions))])
    
    return(sum((unlist(res))^2)) # Unlist all residuals and sum their squared values
  }
  #=================
  #7. Calculate AIC  
  #=================
  #=====================
  # Akaike information criteria corrected for small sample size
  # n = Number of total observations 
  # k = Number of model parameters
  # Observations is a list of lists, each sublist consist of a 2 col-matrix,
  # with the first column being time and the second the actual observations. 
  # Predictions is a dataframe with time as its first column and the rest
  # of the columns being the tissues of interest. The names of the predictions and
  # observations should be the same.
  AICc <- function(k, predictions, observations, n = NULL, times=NULL){
    
    # calculate n in case it is not given
    if(is.null(n) & is.null(times)){
      n <- 0 
      for (i in 1:length(observations)) {
        n <- n + dim(observations[[i]])[1]
      }
    }else if(is.null(n) & !is.null(times)){
      n <- length(observations) * lengths(times)
    }
    print(paste0("RSS is: ",RSS(predictions,observations,times)))
    AICc <- n*log(RSS(predictions,observations,times)/n) + 2*k + (2*k*(k+1))/(n-k-1)
    return(AICc)
  }
 
  ####################################
  #===================================
  #  *** MAX PROBLEM ***  
  #===================================
  #####################################
  
  # Nelder-Mead from dfoptim package
  y_init <- c(dose, rep(0,19))
  time_points <- c(1,3,7, 15, 30)*24 # hours
  excretion_time_points <- excretion_time
  sample_time <- seq(0, 30*24, 1)
  # Initialise vector of physiological parameters
  phys_pars <- create.params(compartments,mass)
  
  #---------------------------
  # Define fitting parameters 
  #---------------------------
  N_p <-8 #   Number of partition coefficients
  N_x <- 8#   Number of permeability coefficients
  # Convert the binary encoding to integer
  grouping <- c(1:8,1:8)
  # Define size of P and X groups
  P_groups <- length(unique(grouping[1:N_p]))  # sample size
  X_groups <- length(unique(grouping[(N_p+1):(N_p+N_x)]))  # sample size
 # set.seed(0)
  # Initilise parameter values
  fitted <- log(exp(runif(P_groups+X_groups+2, -2,2)))
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
  # Run the Nelder Mead algorithmm to estimate the parameter values
  nm_optimizer_max <- dfoptim::nmk(par = fitted, fn = obj.func,
                                   control = list(maxfeval=2000, trace=T), y_init = y_init,
                                   time_points = time_points,
                                   excretion_time_points =  excretion_time_points,
                                   sample_time = sample_time,
                                   phys_pars = phys_pars, 
                                   position = position )

  # Extract the converged parameter values in the log space
  params <- nm_optimizer_max$par
  # Create the matrix of the system  
  A <- create_ODE_matrix(phys_pars = phys_pars, fit_pars =exp(params),  position = position )
  # Solve the ODE system using the exponential matrix method  
  solution <-  solve_exp_matrix(x = A, time = sample_time, y_init = y_init,phys_pars = phys_pars )
  
  observed <- list()
  for (i in 1:(length(df))) {
    observed[[i]] <- cbind(time_points, df[,i])
  }
  observed[[i+1]] <-  cbind(excretion_time_points,excretion[,1]) #feces
  observed[[i+2]] <-  cbind(excretion_time_points,excretion[,2]) #urine
  names(observed) <- c(names(df), names(excretion))
  
  predicted <- solution
  names(predicted) <- c("Times",names(df), "Feces", "Urine")
  
  #Obtain AIC for predictions vs observations
  AIC_result <- AICc(k =length(params), predictions = predicted, observations = observed)
  # GA solves a maximisation problem, and best model gives minimum AIC, so take opposite of AIC
  fit_value_max <- AIC_result
  print(paste0(" AIC value is ", fit_value_max))
  
  print(paste0("PBPK index is ", nm_optimizer_max$value))
  
  
  
  ####################################
  #===================================
  #  *** MIN PROBLEM ***  
  #===================================
  #####################################
  #---------------------------
  # Define fitting parameters 
  #---------------------------
  # Convert the binary encoding to integer
  grouping <- rep(1,16)
  # Define size of P and X groups
  P_groups <- length(unique(grouping[1:N_p]))  # sample size
  X_groups <- length(unique(grouping[(N_p+1):(N_p+N_x)]))  # sample size
  # set.seed(0)
  # Initilise parameter values
  fitted <- log(exp(runif(P_groups+X_groups+2, -2,2)))
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
  # Run the Nelder Mead algorithmm to estimate the parameter values
  nm_optimizer_min <- dfoptim::nmk(par = fitted, fn = obj.func,
                                   control = list(maxfeval=100, trace=T), y_init = y_init,
                                   time_points = time_points,
                                   excretion_time_points =  excretion_time_points,
                                   sample_time = sample_time,
                                   phys_pars = phys_pars, 
                                   position = position )
  
  # Extract the converged parameter values in the log space
  params <- nm_optimizer_min$par
  # Create the matrix of the system  
  A <- create_ODE_matrix(phys_pars = phys_pars, fit_pars =exp(params),  position = position )
  # Solve the ODE system using the exponential matrix method  
  solution <-  solve_exp_matrix(x = A, time = sample_time, y_init = y_init,phys_pars = phys_pars )
  
  observed <- list()
  for (i in 1:(length(df))) {
    observed[[i]] <- cbind(time_points, df[,i])
  }
  observed[[i+1]] <-  cbind(excretion_time_points,excretion[,1]) #feces
  observed[[i+2]] <-  cbind(excretion_time_points,excretion[,2]) #urine
  names(observed) <- c(names(df), names(excretion))
  
  predicted <- solution
  names(predicted) <- c("Times",names(df), "Feces", "Urine")
  
  #Obtain AIC for predictions vs observations
  AIC_result <- AICc(k =length(params), predictions = predicted, observations = observed)
  # GA solves a maximisation problem, and best model gives minimum AIC, so take opposite of AIC
  fit_value_min <- AIC_result
  print(paste0(" AIC value is ", fit_value_min))
  
  print(paste0("Fitness metric is ", nm_optimizer_min$value))
  
  
  
  ####################################
  #===================================
  #  *** BEST BINARY PROBLEM ***  
  #===================================
  #####################################
  #---------------------------
  # Define fitting parameters 
  #---------------------------
  #==================
  #5.Binary mapping 
  #==================
  # Function for mapping the binary number to integer
  # Since with 4 digits numbers from 0-15 can be mapped and here we have 8 
  # different compartments, every two integers correspond to one compartment
  bin2int <- function(bin_seq){
    int <- GA::binary2decimal(bin_seq)
    if(int == 0 || int == 1){
      return(1)
    }else if(int == 2 || int == 3){
      return(2)
    }else if(int == 4 || int == 5){
      return(3)
    }else if(int == 6 || int == 7){
      return(4)
    }else if(int == 8 || int == 9){
      return(5)
    }else if(int == 10 || int == 11){
      return(6)
    }else if(int == 12 || int == 13){
      return(7)
    }else if(int == 14 || int == 15){
      return(8)
    }
  }
  
  
  #=============================
  #6. Convert binary to grouping  
  #=============================
  # Function for converting binary into integer (from )
  decode_ga_bin <- function(binary_num)
  { 
    # Convert binary encoding to gray encoding to avoid the Hamming cliff problem
    gray_num <- GA::gray2binary(binary_num) 
    gray_num <- binary_num 
    
    #Four digit binary encodes up to 15, if we are past 13, assign the value 13
    
    # Partition coefficient grouping
    P1 <-bin2int(gray_num[1:4])
    P2 <-bin2int(gray_num[5:8])
    P3 <-bin2int(gray_num[9:12])
    P4 <-bin2int(gray_num[13:16])
    P5 <-bin2int(gray_num[17:20])
    P6 <-bin2int(gray_num[21:24])
    P7 <-bin2int(gray_num[25:28])
    P8 <-bin2int(gray_num[29:32])
    
    
    # Permeability coefficient grouping
    X1 <-bin2int(gray_num[33:36])
    X2 <-bin2int(gray_num[37:40])
    X3 <-bin2int(gray_num[41:44])
    X4 <-bin2int(gray_num[45:48])
    X5 <-bin2int(gray_num[49:52])
    X6 <-bin2int(gray_num[53:56])
    X7 <-bin2int(gray_num[57:60])
    X8 <-bin2int(gray_num[61:64])
    
    out <- structure(c(P1,P2,P3,P4,P5,P6,P7,P8, X1,X2,X3,X4,X5,
                       X6,X7,X8), names = c("P1","P2","P3","P4",
                                            "P5","P6", "P7", "P8", "X1",
                                            "X2", "X3", "X4", "X5", "X6", "X7", "X8"))
    return(out)
  }
  
  # Nelder-Mead from dfoptim package
  y_init <- c(dose, rep(0,19))
  time_points <- c(1,3,7, 15, 30)*24 # hours
  excretion_time_points <- excretion_time
  sample_time <- seq(0, 30*24, 1)
  # Initialise vector of physiological parameters
  phys_pars <- create.params(compartments,mass)
  
  #---------------------------
  # Define fitting parameters 
  #---------------------------
  N_p <-8 #   Number of partition coefficients
  N_x <- 8#   Number of permeability coefficients
  #---------------------------
  # Define fitting parameters 
  #---------------------------
  GA_results_bin <- GA_results
  # Convert the binary encoding to integer
  grouping <- decode_ga_bin(GA_results_bin@solution[1,])
  # Define size of P and X groups
  P_groups <- length(unique(grouping[1:N_p]))  # sample size
  X_groups <- length(unique(grouping[(N_p+1):(N_p+N_x)]))  # sample size
  # set.seed(0)
  # Initilise parameter values
  fitted <- log(exp(runif(P_groups+X_groups+2, -2,2)))
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
  # Run the Nelder Mead algorithmm to estimate the parameter values
  nm_optimizer_bin <- dfoptim::nmk(par = fitted, fn = obj.func,
                                   control = list(maxfeval=2000, trace=T), y_init = y_init,
                                   time_points = time_points,
                                   excretion_time_points =  excretion_time_points,
                                   sample_time = sample_time,
                                   phys_pars = phys_pars, 
                                   position = position )
  
  # Extract the converged parameter values in the log space
  params <- nm_optimizer_bin$par
  # Create the matrix of the system  
  A <- create_ODE_matrix(phys_pars = phys_pars, fit_pars =exp(params),  position = position )
  # Solve the ODE system using the exponential matrix method  
  solution <-  solve_exp_matrix(x = A, time = sample_time, y_init = y_init,phys_pars = phys_pars )
  
  observed <- list()
  for (i in 1:(length(df))) {
    observed[[i]] <- cbind(time_points, df[,i])
  }
  observed[[i+1]] <-  cbind(excretion_time_points,excretion[,1]) #feces
  observed[[i+2]] <-  cbind(excretion_time_points,excretion[,2]) #urine
  names(observed) <- c(names(df), names(excretion))
  
  predicted <- solution
  names(predicted) <- c("Times",names(df), "Feces", "Urine")
  
  #Obtain AIC for predictions vs observations
  AIC_result <- AICc(k =length(params), predictions = predicted, observations = observed)
  # GA solves a maximisation problem, and best model gives minimum AIC, so take opposite of AIC
  print(paste0(" AIC value is ", AIC_result))
  
  print(paste0("PBPK index is ", nm_optimizer_bin$value))
  
  
  ####################################
  #===================================
  #  *** BEST REAL PROBLEM ***  
  #===================================
  #####################################
  #---------------------------
  # Define fitting parameters 
  #---------------------------
  # Function for converting binary into integer (from )
  decode_ga_real <- function(real_num)
  { 
    # Partition coefficient grouping
    P1 <- floor(real_num[1])
    P2 <- floor(real_num[2])
    P3 <- floor(real_num[3])
    P4 <- floor(real_num[4])
    P5 <- floor(real_num[5])
    P6 <- floor(real_num[6])
    P7 <- floor(real_num[7])
    P8 <- floor(real_num[8])
    
    
    # Permeability coefficient grouping
    X1 <- floor(real_num[9])
    X2 <- floor(real_num[10])
    X3 <- floor(real_num[11])
    X4 <- floor(real_num[12])
    X5 <- floor(real_num[13])
    X6 <- floor(real_num[14])
    X7 <- floor(real_num[15])
    X8 <- floor(real_num[16])
    
    out <- structure(c(P1,P2,P3,P4,P5,P6,P7,P8, X1,X2,X3,X4,X5,
                       X6,X7,X8), names = c("P1","P2","P3","P4",
                                            "P5","P6", "P7", "P8", "X1",
                                            "X2", "X3", "X4", "X5", "X6", "X7", "X8"))
    return(out)
  }

  GA_results_real <- GA_results
  # Convert the binary encoding to integer
  grouping <- decode_ga_real(GA_results_real@solution[1,])
  # Define size of P and X groups
  P_groups <- length(unique(grouping[1:N_p]))  # sample size
  X_groups <- length(unique(grouping[(N_p+1):(N_p+N_x)]))  # sample size
  # set.seed(0)
  # Initilise parameter values
  fitted <- log(exp(runif(P_groups+X_groups+2, -2,2)))
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
  # Run the Nelder Mead algorithmm to estimate the parameter values
  nm_optimizer_real <- dfoptim::nmk(par = fitted, fn = obj.func,
                                   control = list(maxfeval=2000, trace=T), y_init = y_init,
                                   time_points = time_points,
                                   excretion_time_points =  excretion_time_points,
                                   sample_time = sample_time,
                                   phys_pars = phys_pars, 
                                   position = position )
  
  # Extract the converged parameter values in the log space
  params <- nm_optimizer_bin$par
  # Create the matrix of the system  
  A <- create_ODE_matrix(phys_pars = phys_pars, fit_pars =exp(params),  position = position )
  # Solve the ODE system using the exponential matrix method  
  solution <-  solve_exp_matrix(x = A, time = sample_time, y_init = y_init,phys_pars = phys_pars )
  
  observed <- list()
  for (i in 1:(length(df))) {
    observed[[i]] <- cbind(time_points, df[,i])
  }
  observed[[i+1]] <-  cbind(excretion_time_points,excretion[,1]) #feces
  observed[[i+2]] <-  cbind(excretion_time_points,excretion[,2]) #urine
  names(observed) <- c(names(df), names(excretion))
  
  predicted <- solution
  names(predicted) <- c("Times",names(df), "Feces", "Urine")
  
  #Obtain AIC for predictions vs observations
  AIC_result <- AICc(k =length(params), predictions = predicted, observations = observed)
  # GA solves a maximisation problem, and best model gives minimum AIC, so take opposite of AIC
  fit_value_real <- -AIC_result
  print(paste0("Minus AIC value is ", fit_value_real))
  
  print(paste0("PBPK index is ", nm_optimizer_real$value))
  
  