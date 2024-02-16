# Set the directory from where all relevant data all loaded
setwd("C:/Users/ptsir/Documents/GitHub/PBPK_Genetic_Algorithm/Kreyling/NLOPTR")

# Load results from genetic algorithm grouping
load("FPG_nloptr.RData")
GA_results_FPG <- GA_results@solution[1,]  
load("PNG_nloptr.RData")
GA_results_PNG <- GA_results@solution[1,]  
load("SPPCG_nloptr.RData")
GA_results_SPPCG <- GA_results@solution[1,]  
# Clear all results except from the groupings
rm(list=ls()[! ls() %in% c("GA_results_FPG", "GA_results_PNG","GA_results_SPPCG")])


 #####################################
  ### Function to create Parameters ###
  #####################################
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
    
    W_tis <- rep(0,length(comp_names))
    V_tis <- rep(0,length(comp_names))
    V_cap <- rep(0,length(comp_names))
    Q <- rep(0,length(comp_names))
    
    # The following values were calculated by dividing the %ID/ g tissue with the %ID w/o free 48 from Table 2 of Kreyling et al. (2017)
    # Thus, they represent the average mass, in grams, of the respective tissues in each time group.
    liver_expw <- mean(c(8.57, 8.92, 9.30, 8.61, 9.20))
    spleen_expw <- mean(c(0.93, 0.75, 0.97, 0.68, 0.71))
    kidneys_expw <- mean(c(2.27, 2.36, 2.44, 2.11, 2.26))
    lungs_expw <- mean(c(1.87, 1.60, 1.80, 1.48, 1.31))
    heart_expw <- mean(c(0.89, 1.00, 1.00, 1.00, 0.88))
    blood_expw <- mean(c(16.52, 17.45, 15.33, 18.50, 18.00))
    carcass_expw <- mean(c(206.00, 203.33, 184.00, 202.00, 203.75))
    skeleton_expw <- mean(c(26.15, 27.50, 25.56, 25.79, 25.26))
    soft_tissues <- mean(c(228.57, 253.85, 214.29, 225.93, 231.04))
    
    
    ### Calculation of tissue weights  
    W_tis[2] <- heart_expw
    W_tis[3] <- kidneys_expw
    W_tis[5] <- spleen_expw
    W_tis[6] <- lungs_expw
    W_tis[7] <- liver_expw
    W_tis[9] <- skeleton_expw
    W_tis[13] <- Tissue_fractions[13]*w
    
    for (i in 1:length(comp_names)) {
      control <- comp_names[i]
      
      Regional_flow_fractions[i] <- ifelse(is.na(control), NA, Regional_flow_fractions[i])
      Capillary_fractions[i] <- ifelse(is.na(control), NA, Capillary_fractions[i])
      
      
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
      
      
      ###Calculation of regional blood flows
      Q[i] <- Q_total*Regional_flow_fractions[i]
    }
    
    
    ### Calculations for "Soft tissue" compartment
    W_tis[1] <- w - sum(W_tis[2:length(W_tis)], na.rm = TRUE)-Total_Blood
    V_tis[1] <- W_tis[1]/d_adipose     
    Q[1] <- Q_total - sum(Q[2:length(Q)],na.rm = TRUE) + Q[6]
    V_cap[1] <- V_tis[1]*Capillary_fractions[1] #Total_Blood - Vven - Vart - sum(V_cap[2:length(V_cap)], na.rm = TRUE)
    
    parameters <- matrix(c(W_tis[],V_tis[],V_cap[],Q[]), ncol = 4)
    colnames(parameters) <- c("W_tis", "V_tis", "V_cap", "Q")
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
  #===============================================
  #2. Function to create initial values for ODEs 
  #===============================================
  
  create.inits <- function(parameters, dose){
    with( as.list(parameters),{
      M_ht<-0; M_lu<-0; M_li<-0; M_spl<-0; M_ki<-0; M_git<-0; M_bone<-0; M_rob<-0;
      
      M_cap_ht<-0; M_cap_lu<-0; M_cap_li<-0; M_cap_spl<-0; M_cap_ki<-0; M_cap_git<-0; M_cap_bone<-0; M_cap_rob<-0;
      
      M_lumen <- 0;
      M_ven <- dose; M_art<-0
      M_feces<-0; M_urine<-0 
      
      return(c("M_ht" = M_ht, "M_lu" = M_lu, 
               "M_li" = M_li, "M_spl" = M_spl, 
               "M_ki" = M_ki, "M_git" = M_git, 
               "M_bone" = M_bone,"M_rob"=M_rob,
               
               "M_cap_ht" = M_cap_ht, "M_cap_lu" = M_cap_lu, 
               "M_cap_li" = M_cap_li, "M_cap_spl" = M_cap_spl, 
               "M_cap_ki" = M_cap_ki, "M_cap_git" = M_cap_git, 
               "M_cap_bone" = M_cap_bone,"M_cap_rob"=M_cap_rob,
               
               "M_lumen" = M_lumen,
               "M_ven" = M_ven, "M_art" = M_art, "M_feces" = M_feces, "M_urine" = M_urine))
      
    })
  }
  
  #==============
  #3. ODEs System
  #==============
  ode.func <- function(time, inits, params){
    position <- params[37:(37+15)]
    fit_pars <- params[(37+16):length(params)]
    
    with(as.list(c(inits, params)),{
      
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
      CLE_h <- fit_pars[length(fit_pars)]
      CLE_u <- 0
      
      
      # Concentrations (mg of NPs)/(g of wet tissue)
      C_ht <- M_ht/w_ht
      C_cap_ht <- M_cap_ht/V_cap_ht
      C_lu <- M_lu/w_lu
      C_cap_lu <- M_cap_lu/V_cap_lu
      C_li <- M_li/w_li
      C_cap_li <- M_cap_li/V_cap_li
      C_spl <- M_spl/w_spl
      C_cap_spl <- M_cap_spl/V_cap_spl
      C_ki <- M_ki/w_ki
      C_cap_ki <- M_cap_ki/V_cap_ki
      C_git <- M_git/w_git
      C_cap_git <- M_cap_git/V_cap_git
      C_bone <- M_bone/w_bone
      C_cap_bone <- M_cap_bone/V_cap_bone
      C_rob <- M_rob/w_rob
      C_cap_rob <- M_cap_rob/V_cap_rob
      
      C_ven <- M_ven/V_ven
      C_art <- M_art/V_art
      
      # Heart
      dM_cap_ht <- Q_ht*(C_art - C_cap_ht) - x_ht*Q_ht*(C_cap_ht - C_ht/P_ht)
      dM_ht <- x_ht*Q_ht*(C_cap_ht - C_ht/P_ht) 
      
      # Lungs
      dM_cap_lu <- Q_total*(C_ven - C_cap_lu) - x_lu*Q_total*(C_cap_lu - C_lu/P_lu)
      dM_lu <-  x_lu*Q_total*(C_cap_lu - C_lu/P_lu)
      
      # Liver 
      # Liver 
      dM_cap_li <- Q_li*(C_art - C_cap_li) + Q_spl*(C_cap_spl - C_cap_li) + Q_git*(C_cap_git - C_cap_li) -
        x_li*(Q_li+Q_git+Q_spl)*(C_cap_li - C_li/P_li)
      dM_li <- x_li*(Q_li+Q_git+Q_spl)*(C_cap_li - C_li/P_li) - CLE_h*M_li
      
      # Spleen
      dM_cap_spl <- Q_spl*(C_art - C_cap_spl) - x_spl*Q_spl*(C_cap_spl - C_spl/P_spl)
      dM_spl <- x_spl*Q_spl*(C_cap_spl - C_spl/P_spl) 
      
      # Kidneys
      dM_cap_ki <- Q_ki*(C_art - C_cap_ki) - x_ki*Q_ki*(C_cap_ki - C_ki/P_ki)- CLE_u*M_cap_ki
      dM_ki <- x_ki*Q_ki*(C_cap_ki - C_ki/P_ki) 
      
      # GIT - Gastrointestinal Track
      dM_cap_git <- Q_git*(C_art - C_cap_git) - x_git*Q_git*(C_cap_git - C_git/P_git)
      dM_git <- x_git*Q_git*(C_cap_git - C_git/P_git) 
      dM_lumen <- CLE_h*M_li - CLE_f *M_lumen 
      
      # Bone
      dM_cap_bone <- Q_bone*(C_art - C_cap_bone) - x_bone*Q_bone*(C_cap_bone - C_bone/P_bone)
      dM_bone <- x_bone*Q_bone*(C_cap_bone - C_bone/P_bone) 
      
      
      # RoB - Rest of Body
      dM_cap_rob <- Q_rob*(C_art - C_cap_rob) - x_rob*Q_rob*(C_cap_rob - C_rob/P_rob)
      dM_rob <- x_rob*Q_rob*(C_cap_rob - C_rob/P_rob) 
      
      # Urine
      dM_urine <- CLE_u*M_cap_ki
      
      # Feces
      dM_feces <- CLE_f*M_lumen
      
      # Venous Blood
      dM_ven <- Q_ht*C_cap_ht + (Q_li + Q_spl+Q_git)*C_cap_li + Q_ki*C_cap_ki +
        Q_bone*C_cap_bone + Q_rob*C_cap_rob - Q_total*C_ven
      
      # Arterial Blood
      dM_art <-  Q_total*C_cap_lu - Q_total*C_art
      
      Blood_total <- M_ven + M_art + M_cap_ht + M_cap_lu +M_cap_li+M_cap_spl+
        M_cap_ki+ M_cap_git+M_cap_bone+M_cap_rob
      Blood <- Blood_total/(V_blood)
      
      C_soft <- (M_git+M_lumen+M_rob)/(w_git + w_rob)
      
      list(c("dM_ht" = dM_ht, "dM_lu" = dM_lu, 
             "dM_li" = dM_li, "dM_spl" = dM_spl, 
             "dM_ki" = dM_ki, "dM_git" = dM_git, 
             "dM_bone" = dM_bone,"dM_rob"=dM_rob,
             
             "dM_cap_ht" = dM_cap_ht, "dM_cap_lu" = dM_cap_lu, 
             "dM_cap_li" = dM_cap_li, "dM_cap_spl" = dM_cap_spl, 
             "dM_cap_ki" = dM_cap_ki, "dM_cap_git" = dM_cap_git, 
             "dM_cap_bone" = dM_cap_bone,"dM_cap_rob"=dM_cap_rob,
             
             "dM_lumen" = dM_lumen,
             "dM_ven" = dM_ven, "dM_art" = dM_art, "dM_feces" = dM_feces, "dM_urine" = dM_urine),
           
           "Blood"=Blood,
           "C_ht"=C_ht, "C_lu"=C_lu, "C_li"=C_li, "C_spl"=C_spl,
           "C_ki"=C_ki,  "C_bone"=C_bone, "C_soft"=C_soft,
           "Feces"=M_feces)
    })
  }
  #======================
  #3. Objective function  
  #======================
  obj.func <- function(par,phys_pars, position,  sample_time, inits, time_points, excretion_time_points){
    
    params <- c(phys_pars, position, exp(par))
    solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                        y = inits, parms = params, 
                                        method="lsodes",rtol = 1e-3, atol = 1e-3))
    
    concentrations <- data.frame(solution$time, solution$C_li, solution$C_spl, solution$C_ki,
                                 solution$C_lu, solution$C_ht, solution$Blood,
                                 solution$C_bone,  solution$C_soft)
    
    concentrations <- concentrations[solution$time %in%time_points, 2:dim(concentrations)[2]]
    excr_solution <-  data.frame(solution$time, solution$Feces)
    excr_solution <- excr_solution[solution$time %in% excretion_time_points,2]
    observed <- list()
    predicted <- list()
    
    for (i in 1:length(concentrations)) {
      observed[[i]] <- df[,i]
      predicted[[i]] <- concentrations[,i]
    }
    observed[[i+1]] <- excretion #feces
    predicted[[i+1]] <- excr_solution #feces
    
    discrepancy <- PBKOF(observed, predicted)
    
    return(discrepancy)
  }
  
  #===============
  #4. Improved fitness metric  
  #===============
  
  ############# Calculate PBPK indices #############
  # PBKOF function the returns the compartment and consolidated (Total) discrepancy index
  # of a PBPK model, given some experimental data. It follows the paper of Krishnan et al.1995.
  # observed: list of vectors containing the experimental data
  # predictions: list of vectors containing the predicted data
  # names of the compartments
  
  
  PBKOF <- function(observed, predicted, comp.names =NULL){
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
  
# Create functions for decoding the results of the Genetic Algorithms grouping 
# Function for decoding FPG and PNG schemes
  decode_ga <- function(real_num){ 
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
  
  # Function for decoding SPPCG scheme
  decode_ga_sppcg <- function(real_num)
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
    
    
    out <- structure(c(P1,P2,P3,P4,P5,P6,P7,P8), names = c("P1","P2","P3","P4",
                                                           "P5","P6", "P7", "P8"))
    return(out)
  }
  
  #=============================
  #8. Create position  
  #=============================  
  # Function for creating the position from which to draw each param from the fitted params vector
  create.position <- function(grouping){
    #---------------------------
    # Define fitting parameters 
    #---------------------------
    N_p <-8 #   Number of partition coefficients
    N_x <- 8#   Number of permeability coefficients  
    # Define size of P and X groups
    P_groups <- length(unique(grouping[1:N_p]))  # sample size
    X_groups <- length(unique(grouping[(N_p+1):(N_p+N_x)]))  # sample size
    # set.seed(0)
    # Initilise parameter values
    fitted <- rep(NA, P_groups+X_groups+2)
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
    names(fitted) <- c(pnames, xnames,"CLE_f",  "CLE_h")
    # Variable for keeping which value in the fitted params vector corresponds to each coefficient
    position = rep(NA, length(grouping))
    for (i in 1:(length(position))){
      if(i<=8){
        position[i] <- which(names(fitted) == paste0("P", as.character(grouping[i])))
      }else{
        position[i] <- which(names(fitted) == paste0("X", as.character(grouping[i])))
      }
    }
    fitted[] <-  c(log(rep(10,P_groups)),log(rep(0.01,X_groups)), log(0.16), log(4.5e-05))
    
    return(list("position"=position,"fitted"=fitted, 'P_groups' = P_groups, X_groups = X_groups))
  }
  
  
  create.position_sppcg <- function(grouping){
    #---------------------------
    # Define fitting parameters 
    #---------------------------
    N_group <- 8 #   Number of groups fro partition and permeability coefficients
    # Define size of P and X groups
    P_groups <- length(unique(grouping))  # sample size
    X_groups <- P_groups
    # set.seed(0)
    # Initilise parameter values
    fitted <- rep(NA, P_groups+X_groups+2)
    # Initialise naming vectors
    pnames <- rep(NA, P_groups)
    xnames <- rep(NA, X_groups)
    
    #Define names for P and X groups
    for (i in 1:P_groups){
      pnames[i] <- paste0("P", as.character(unique(grouping[1:N_group])[i]))
    }
    for (j in 1:X_groups){
      xnames[j] <- paste0("X", as.character(unique(grouping[1:N_group])[j]))
    }
    # Define the total parameter vector names
    names(fitted) <- c(pnames, xnames,"CLE_f",  "CLE_h")
    # Variable for keeping which value in the fitted params vector corresponds to each coefficient
    position = rep(NA, 2*length(grouping))
    for (i in 1:(length(position)/2)){
      position[i] <- which(names(fitted) == paste0("P", as.character(grouping[i])))
      position[i+8] <- position[i] + P_groups
    }
    fitted[] <-  c(log(rep(10,P_groups)),log(rep(0.01,X_groups)), log(0.16), log(4.5e-05))
    
    return(list("position"=position,"fitted"=fitted, 'P_groups' = P_groups, X_groups = X_groups))
  }
  
  #===============
  # Load data  
  #===============
  
  dose <- 18.15 # ug # Since results are in %ID we used a random dose that is within the dose range given to the rats
  mass <- 263 #g, female Wistar Kyoto rats
  
  # Load raw data from paper Kreyling et al.2017, which are given in %ID/g tissue
  df_percent <- openxlsx::read.xlsx("Kreyling-IV-data.xlsx", sheet = 6, colNames = T, rowNames = T) # TiO2 NPs %ID/g of tissue  (Table 1)
  excretion_percent <- openxlsx::read.xlsx("Kreyling-IV-data.xlsx", sheet = 3, colNames = T, rowNames = F) # accumulated excretory rate, expressed as %ID
  # Drop the first time points because the graph is supposed to be cumulative dose but the cumulative feces in day 1 are less that the first time points 
  excretion_time <- round(excretion_percent[3:5,1])*24 # hours
  # Convert doses from %ID to masses
  df_all <- (df_percent/100) * dose # Concentrations in (ug of NPs)/(g of tissue)
  # Drop unused compartments
  df <- df_all[, !(names(df_all) %in% c("Uterus", "Brain", "Carcass"))]
  excretion <- (excretion_percent[3:5,2]/100) * dose
  
  
  
  compartments <- list( "RoB"="RoB","Heart"="Heart", "Kidneys"="Kidneys", "Brain"= NA, "Spleen"="Spleen",
                        "Lungs"="Lungs", "Liver"="Liver", "Uterus"= NA, "Bone"="Bone", "Adipose"=NA, "Skin"=NA, "Muscles"=NA, "GIT"="GIT") #used as input in function, compartments that are used in pbpk
  
  
  # Nelder-Mead from dfoptim package
  time_points <- c(1,4,24, 7*24, 28*24) # hours
  excretion_time_points <- excretion_time
  sample_time <- seq(0, 28*24, 1)
  # Initialise vector of physiological parameters
  phys_pars <- create.params(compartments,mass)
  inits <- create.inits(phys_pars, dose)
  
  # Create the parameter grouping for all five estimation schemes
  grouping_MANG <- c(1:8, 1:8)
  grouping_MING <- c(rep(1,8), rep(1,8))
  grouping_FPG <- decode_ga(GA_results_FPG)  
  grouping_PNG <- decode_ga(GA_results_PNG)  
  grouping_SPPCG <- decode_ga_sppcg(GA_results_SPPCG)  
  
  
  # Create the position vector to match the ODE parameters with the fitted parameter values
  
  grouping_position_MANG <- create.position(grouping_MANG)
  position_MANG <- grouping_position_MANG$position
  P_groups_MANG <- grouping_position_MANG$P_groups
  X_groups_MANG <- grouping_position_MANG$X_groups
  fitted_MANG <- grouping_position_MANG$fitted
  
  grouping_position_MING <- create.position(grouping_MING)
  position_MING <- grouping_position_MING$position
  P_groups_MING <- grouping_position_MING$P_groups
  X_groups_MING <- grouping_position_MING$X_groups
  fitted_MING <- grouping_position_MING$fitted
  
  
  grouping_position_FPG <- create.position(grouping_FPG)
  position_FPG <- grouping_position_FPG$position
  P_groups_FPG <- grouping_position_FPG$P_groups
  X_groups_FPG <- grouping_position_FPG$X_groups
  fitted_FPG <- grouping_position_FPG$fitted
  
  
  grouping_position_PNG <- create.position(grouping_PNG)
  position_PNG <- grouping_position_PNG$position
  P_groups_PNG <- grouping_position_PNG$P_groups
  X_groups_PNG <- grouping_position_PNG$X_groups
  fitted_PNG <- grouping_position_PNG$fitted
  
  
  grouping_position_SPPCG <- create.position_sppcg(grouping_SPPCG)
  position_SPPCG <- grouping_position_SPPCG$position
  P_groups_SPPCG <- grouping_position_SPPCG$P_groups
  X_groups_SPPCG <- grouping_position_SPPCG$X_groups
  fitted_SPPCG <- grouping_position_SPPCG$fitted
  
  # Options for the NLOPT algorithm
  opts <- list( "algorithm" = "NLOPT_LN_NEWUOA",
                "xtol_rel" = 1e-06,
                "ftol_rel" = 0.0,
                "ftol_abs" = 0.0,
                "xtol_abs" = 0.0 ,
                "maxeval" = 1000) 

  #=======
  # MANG #
  #=======
  try(
    optimizer_MANG <- nloptr::nloptr( x0= fitted_MANG,
                                      eval_f = obj.func,
                                      lb	= rep(-15, length(fitted_MANG)),
                                      ub = rep(9, length(fitted_MANG)),
                                      opts = opts,
                                      phys_pars = phys_pars, position = position_MANG, sample_time = sample_time,
                                      inits = inits,
                                      time_points = time_points,
                                      excretion_time_points =  excretion_time_points),
    silent = TRUE
  )
  if(exists("optimizer_MANG")){
    params_MANG<- exp(optimizer_MANG$solution)
    value_MANG <- optimizer_MANG$objective
  }else{
    params_MANG<- NULL
    value_MANG <- NULL
  }
  
  #=======
  # MING #
  #=======
  try(
    optimizer_MING <- nloptr::nloptr( x0= fitted_MING,
                                      eval_f = obj.func,
                                      lb	= rep(-15, length(fitted_MING)),
                                      ub = rep(9, length(fitted_MING)),
                                      opts = opts,
                                      phys_pars = phys_pars, position = position_MING, sample_time = sample_time,
                                      inits = inits,
                                      time_points = time_points,
                                      excretion_time_points =  excretion_time_points),
    silent = TRUE
  )
  if(exists("optimizer_MING")){
    params_MING<- exp(optimizer_MING$solution)
    value_MING <- optimizer_MING$objective
  }else{
    params_MING<- NULL
    value_MING <- NULL
  }
  
  #=======
  # FPG #
  #=======
  try(
    optimizer_FPG <- nloptr::nloptr( x0= fitted_FPG,
                                     eval_f = obj.func,
                                     lb	= rep(-15, length(fitted_FPG)),
                                     ub = rep(9, length(fitted_FPG)),
                                     opts = opts,
                                     phys_pars = phys_pars, position = position_FPG, sample_time = sample_time,
                                     inits = inits,
                                     time_points = time_points,
                                     excretion_time_points =  excretion_time_points),
    silent = TRUE
  )
  if(exists("optimizer_FPG")){
    params_FPG<- exp(optimizer_FPG$solution)
    value_FPG <- optimizer_FPG$objective
  }else{
    params_FPG<- NULL
    value_FPG <- NULL
  }
  
  
  #=======
  # PNG #
  #=======
  try(
    optimizer_PNG <- nloptr::nloptr( x0= fitted_PNG,
                                     eval_f = obj.func,
                                     lb	= rep(-15, length(fitted_PNG)),
                                     ub = rep(9, length(fitted_PNG)),
                                     opts = opts,
                                     phys_pars = phys_pars, position = position_PNG, sample_time = sample_time,
                                     inits = inits,
                                     time_points = time_points,
                                     excretion_time_points =  excretion_time_points),
    silent = TRUE
  )
  if(exists("optimizer_PNG")){
    params_PNG<- exp(optimizer_PNG$solution)
    value_PNG <- optimizer_PNG$objective
  }else{
    params_PNG<- NULL
    value_PNG <- NULL
  }
  
  
  #========
  # SPPCG #
  #========
  try(
    optimizer_SPPCG <- nloptr::nloptr( x0= fitted_SPPCG,
                                       eval_f = obj.func,
                                       lb	= rep(-15, length(fitted_SPPCG)),
                                       ub = rep(9, length(fitted_SPPCG)),
                                       opts = opts,
                                       phys_pars = phys_pars, position = position_SPPCG, sample_time = sample_time,
                                       inits = inits,
                                       time_points = time_points,
                                       excretion_time_points =  excretion_time_points),
    silent = TRUE
  )
  if(exists("optimizer_SPPCG")){
    params_SPPCG<- exp(optimizer_SPPCG$solution)
    value_SPPCG <- optimizer_SPPCG$objective
  }else{
    params_SPPCG<- NULL
    value_SPPCG <- NULL
  }


# Number of parameter values per estimation scheme
L_MANG <- length(params_MANG)
L_MING <- length(params_MING)
L_FPG <- length(params_FPG)
L_PNG <- length(params_PNG)
L_SPPCG <- length(params_SPPCG)

#Create one matrix to hold the value of each parameter for each model per initialization
N_par <- 18 
N_model <- 5 #MANG, MING, FPG, PNG,SPPCG
pmatrix <- matrix(rep(NA, N_par*N_model),nrow = N_model)
colnames(pmatrix) <- c("P_ht", "P_lu", "P_li", "P_spl","P_ki", "P_git",
                       "P_bone", "P_rob", "x_ht", "x_lu","x_li", "x_spl",
                       "x_ki", "x_git", "x_bone","x_rob", "CLE_f",  "CLE_h" )
rownames(pmatrix) <-  c("MANG",  "MING", "FPG", "PNG","SPPCG")

# Populate the matrix with the values
pmatrix[1,1:16] <-  unlist(lapply(1:16, function(x) {params_MANG[position_MANG[x]]}))
pmatrix[1,17:18] <- params_MANG[c(length(params_MANG)-1, length(params_MANG))]
pmatrix[2,1:16] <-  unlist(lapply(1:16, function(x) {params_MING[position_MING[x]]}))
pmatrix[2,17:18] <- params_MING[c(length(params_MING)-1, length(params_MING))]
pmatrix[3,1:16] <-  unlist(lapply(1:16, function(x) {params_FPG[position_FPG[x]]}))
pmatrix[3,17:18] <- params_FPG[c(length(params_FPG)-1, length(params_FPG))]
pmatrix[4,1:16] <-  unlist(lapply(1:16, function(x) {params_PNG[position_PNG[x]]}))
pmatrix[4,17:18] <- params_PNG[c(length(params_PNG)-1, length(params_PNG))]
pmatrix[5,1:16] <-  unlist(lapply(1:16, function(x) {params_SPPCG[position_SPPCG[x]]}))
pmatrix[5,17:18] <- params_SPPCG[c(length(params_SPPCG)-1, length(params_SPPCG))]

View(t(signif(pmatrix,3)))

# Create solutions for each model
parameters_MANG <- c(phys_pars, position_MANG, params_MANG)
sol_MANG <-  data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                           y = inits, parms = parameters_MANG, 
                           method="bdf",rtol = 1e-5, atol = 1e-5))

solution_MANG <- data.frame(sol_MANG$time, sol_MANG$C_li, sol_MANG$C_spl, sol_MANG$C_ki,
                            sol_MANG$C_lu, sol_MANG$C_ht, sol_MANG$Blood,
                            sol_MANG$C_bone,  sol_MANG$C_soft, sol_MANG$Feces)
names(solution_MANG) <- c("Time", "Liver",  "Spleen",
                         "Kidneys", "Lungs", "Heart","Blood", "Bone",  "Rob","Feces")



###
parameters_MING <- c(phys_pars, position_MING, params_MING)
sol_MING <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                           y = inits, parms = parameters_MING, 
                           method="bdf",rtol = 1e-5, atol = 1e-5))

solution_MING <- data.frame(sol_MING$time, sol_MING$C_li, sol_MING$C_spl, sol_MING$C_ki,
                            sol_MING$C_lu, sol_MING$C_ht, sol_MING$Blood,
                            sol_MING$C_bone,  sol_MING$C_soft, sol_MING$Feces)
names(solution_MING) <- c("Time", "Liver",  "Spleen",
                          "Kidneys", "Lungs", "Heart","Blood", "Bone",  "Rob","Feces")


###
parameters_FPG <- c(phys_pars, position_FPG, params_FPG)
sol_FPG <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                           y = inits, parms = parameters_FPG, 
                           method="bdf",rtol = 1e-5, atol = 1e-5))

solution_FPG <- data.frame(sol_FPG$time, sol_FPG$C_li, sol_FPG$C_spl, sol_FPG$C_ki,
                            sol_FPG$C_lu, sol_FPG$C_ht, sol_FPG$Blood,
                            sol_FPG$C_bone,  sol_FPG$C_soft, sol_FPG$Feces)
names(solution_FPG) <- c("Time", "Liver",  "Spleen",
                          "Kidneys", "Lungs", "Heart","Blood", "Bone",  "Rob","Feces")

###
parameters_PNG <- c(phys_pars, position_PNG, params_PNG)
sol_PNG <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                           y = inits, parms = parameters_PNG, 
                           method="bdf",rtol = 1e-5, atol = 1e-5))

solution_PNG <- data.frame(sol_PNG$time, sol_PNG$C_li, sol_PNG$C_spl, sol_PNG$C_ki,
                            sol_PNG$C_lu, sol_PNG$C_ht, sol_PNG$Blood,
                            sol_PNG$C_bone,  sol_PNG$C_soft, sol_PNG$Feces )
names(solution_PNG) <- c("Time", "Liver",  "Spleen",
                          "Kidneys", "Lungs", "Heart","Blood", "Bone",  "Rob","Feces")

##
parameters_SPPCG <- c(phys_pars, position_SPPCG, params_SPPCG)
sol_SPPCG <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                           y = inits, parms = parameters_SPPCG, 
                           method="bdf",rtol = 1e-5, atol = 1e-5))

solution_SPPCG <- data.frame(sol_SPPCG$time, sol_SPPCG$C_li, sol_SPPCG$C_spl, sol_SPPCG$C_ki,
                            sol_SPPCG$C_lu, sol_SPPCG$C_ht, sol_SPPCG$Blood,
                            sol_SPPCG$C_bone,  sol_SPPCG$C_soft, sol_SPPCG$Feces)
names(solution_SPPCG) <- c("Time", "Liver",  "Spleen",
                      "Kidneys", "Lungs", "Heart","Blood", "Bone",  "Rob","Feces")

# Create a single data frame to hold the observation data 
observations <- data.frame( Time =c(1,4,24, 7*24, 28*24), df,  c(NA,NA,excretion))
names(observations) <- c("Time", "Liver",  "Spleen",
                         "Kidneys", "Lungs", "Heart","Blood", "Bone",  "Rob","Feces")

###########################################
#==========================================
#               Model Metrics             #
#==========================================
###########################################
#
#Function to estimate the percentage of percent of model-predicted concentrations
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

# Function for estimating all metrics
metric.print <- function(predictions, observations ){

  concentrations <- predictions[predictions$Time %in% observations$Time, ]
  observed <- list()
  predicted <- list()
  
  for (i in 2:(length(concentrations))) {
    if(names(concentrations)[i] == "Feces"){
      observed[[i-1]] <- observations[3:5,i]
      predicted[[i-1]] <- concentrations[3:5,i]
    }else{
      observed[[i-1]] <- observations[,i]
      predicted[[i-1]] <- concentrations[,i]
    }
  }
  
  print(paste0("pbpk index of model is: ", pbpk.index(observed, predicted)))
  print(paste0("R-squared of model is: ", r.squared(observed, predicted)))
  print(paste0("AAFE  of model is: ", AAFE(observed, predicted)))
  print(paste0("RMSD of model is ", rmsd(observed, predicted)))
  print(paste0("2-fold % of model is: ", two.fold(observed, predicted)))
  print(paste0("PBKOF of model is: ", PBKOF(observed, predicted)))
  
}

metric.print(solution_MANG, observations)
metric.print(solution_MING, observations)
metric.print(solution_FPG, observations)
metric.print(solution_PNG, observations)
metric.print(solution_SPPCG, observations)


#####################
#    PLOTS  ########
######################

library(ggplot2)

# Defining the linetype and colour of each curve
ltp <- c("MANG" = "twodash", "FPG" = "solid", "PNG" = "dotted","SPPCG" = "dashed"
         ,"MING" = "dotdash" )
cls <-  c("MANG" = "#56B4E9",  "FPG" ="#000000", "PNG" = "#009E73", "SPPCG" ="#CC79A7",
          "Observations" = "#D55E00"
          , "MING" = "#E69F00")


create.plots <- function(compartment){  
  excreta <- compartment %in% c("Feces", "Urine")
  ggplot(data = solution_MANG)+
             geom_line( aes_string(x= "Time", y= rlang::expr(!!compartment), 
                                   color = '"MANG"',linetype = '"MANG"'),  size=1.5,alpha = 0.7) +
             geom_line(data=solution_FPG, aes_string(x= "Time", y= rlang::expr(!!compartment),
                                                      color = '"FPG"',linetype ='"FPG"'), size=1.5,alpha = 0.7) +
             geom_line(data=solution_MING, aes_string(x= "Time", y= rlang::expr(!!compartment),
                                                      color = '"MING"',linetype ='"MING"'), size=1.5,alpha = 0.7) +
             geom_line(data=solution_PNG, aes_string(x= "Time", y= rlang::expr(!!compartment),
                                                      color =  '"PNG"',linetype =  '"PNG"'), size=1.5,alpha = 0.7) +
             geom_line(data=solution_SPPCG, aes_string(x= "Time", y= rlang::expr(!!compartment), 
                                                      color = '"SPPCG"',linetype ='"SPPCG"'), size=1.5,alpha = 0.7) +
             geom_point(data=observations, aes_string(x="Time", y= rlang::expr(!!compartment), 
                                                      color='"Observations"'), size=4)+
             labs(title = rlang::expr(!!compartment), 
                  y = ifelse(excreta,expression("TiO2 (" * mu * "g)"),expression("TiO2 (" * mu* "g/g tissue)" )),
                  x = "Time (hours)")+
             theme(plot.title = element_text(hjust = 0.5))+
             {if(compartment %in% c("Blood", "Kidneys", "Bone", "Rob", "Lungs", "Heart" ))scale_y_continuous(trans='log10')}+
             scale_color_manual("", values=cls)+
             scale_linetype_manual("Models", values=ltp) +
             theme(legend.justification = "top" ,
               legend.key.size = unit(1.5, 'cm'),  
                   legend.title = element_text(size=14),
                   axis.title=element_text(size=14),
                   legend.text = element_text(size=14))+
    theme(plot.margin=grid::unit(c(0.25,0.25,0.25,0.25), "cm"))
           
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
#gridExtra::grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8, p9,p10, nrow = 4)
#gridExtra::grid.arrange(p5,p6,p7,p8,nrow = 2)
#gridExtra::grid.arrange(p9,p10,nrow = 2)


final_plot<- ggpubr::ggarrange(p1, p2, p3, p4,p5,p6,p7,p8, p9, ncol=3, nrow=4, 
                               common.legend = TRUE, legend="right")+
  theme(plot.margin=grid::unit(c(0.25,0.25,0.25,0.25), "cm"))

ggsave("model_benchmark_low.png",plot = final_plot,
       device='png', dpi=320,
       width = 15,
       height = 12,
       units = "in")
dev.off()

#save.image(file = "Benchmarking.RData")
