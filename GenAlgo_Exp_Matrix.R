library(expm)
setwd("C:/Users/vassi/Documents/GitHub/PBPK_Genetic_Algorithm")

dose_kg <- 10 # mg/kg rat body
mass <- 250 # g  
dose <- dose_kg*mass/1000 # mg TiO2

# Load raw data from paper Xie et al.2011
df <- openxlsx::read.xlsx("TiO2_iv_rat.xlsx", sheet = 1, colNames = T, rowNames = T) # TiO2 NPs %ID/g of tissue  (Table 1)
excretion <- openxlsx::read.xlsx("Cummulative_Excretion.xlsx", sheet = 1, colNames = T, rowNames = F) # accumulated excretory rate, expressed as %ID
excretion <- excretion[,c(2:3)]

# Transform to (mg of NPs)/(g of tissue)
df <- (df/100)*dose
excretion <- (excretion/100)*dose

####################################################################
###################             DATA             ###################
####################################################################


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


####################
### User's INPUT ###
####################
#### If any of these compartments don not exist in pbpk, just give it the value NA in compartments vector, example: "Heart" = NA and it will remove it 
#### from the equilibriums and the corresponding V_tis, V_cap, Q will be equal to NA.


compartments <- list( "RoB"="RoB","Heart"="Heart", "Kidneys"="Kidneys", "Brain"= NA, "Spleen"="Spleen",
                      "Lungs"="Lungs", "Liver"="Liver", "Uterus"=NA, "Bone"="Bone", "Adipose"=NA, "Skin"=NA, "Muscles"=NA, "GIT"="GIT") #used as input in function, compartments that are used in pbpk


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

params<-create.params(compartments,mass)

x_gen <- 0.001 # random value - unitless
P_gen <- 0.001 # random value - unitless
CLE_f <- 0.001
CLE_u <- 0.001

x_ht <- x_gen
x_lu <- x_gen
x_li <- x_gen
x_spl <- x_gen
x_ki <- x_gen
x_git <- x_gen 
x_bone <- x_gen
x_rob <- x_gen

P_ht <- P_gen
P_lu <- P_gen
P_li <- P_gen
P_spl <- P_gen
P_ki <- P_gen
P_git <- P_gen 
P_bone <- P_gen
P_rob <- P_gen


#--------------------------------------------------------------------------------------------------
# "create_ODE_matrix()" creates the matrix with the coefficients of the state variables of the 
# desired ODE system. It takes as input the values of parameters and returns the matrix.
#--------------------------------------------------------------------------------------------------
create_ODE_matrix <- function(parameters){
  with( as.list(parameters),{
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

A <- create_ODE_matrix(params)

y_init <- c(dose, rep(0,19))

sample_time <- c(0,1,3,7, 15, 30)*24 # hours


#--------------------------------------------------------------------------------------------------
# "Solve_exp_matrix()" is a function that solves the ODE system using the matrix "x" (which 
# contains the coefficients of the system), "time" which is the desired time points to 
# be calculated and "y_init" is the initial values of the state variables.
#--------------------------------------------------------------------------------------------------

Solve_exp_matrix <- function(x, time, y_init){
  
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
    solution_t <- expm(x*time[t])%*%y_init
    y_t[,t] <- solution_t
  }
  rownames(y_t) <- rownames(A)
  return(t(y_t))
}

start_time <- Sys.time()
solution <-  Solve_exp_matrix(x = A, time = sample_time, y_init = y_init)
end_time <- Sys.time()

Matrix_solution_duration <- end_time - start_time
Matrix_solution_duration