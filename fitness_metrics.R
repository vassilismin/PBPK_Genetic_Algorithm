# predictions <- solution
# names(predictions) <- c("Time","Blood","Heart", "Lungs", "Liver", "Spleen", "Kidneys", "Git", "Bone", "Feces", "Urine")
# #observations <- df
# 
# observations <- list()
# time_points <- c(24,72,168,360,720)
# 
# for (i in 1:dim(df)[2]) {
#   observations[[i]] <- matrix(c(time_points, df[,i]), nrow = length(time_points), byrow = FALSE)  
#   names(observations)[i] <- names(df)[i]
#   colnames(observations[[i]]) <- c("Time", names(df)[i])
# }
# observations[[i+1]] <- matrix(c(excretion_time_points, excretion[,1]), nrow = length(excretion_time_points), byrow = FALSE) # feces
# names(observations)[i+1] <- "Feces"
# colnames(observations[[i+1]]) <- c("Time", "Feces")
# observations[[i+2]] <- matrix(c(excretion_time_points, excretion[,2]), nrow = length(excretion_time_points), byrow = FALSE) # urine
# names(observations)[i+2] <- "Urine"
# colnames(observations[[i+2]]) <- c("Time", "Urine")


#=====================
#   RSS Function
#=====================

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
      different_times <- TRUE
      break
    }
  }
  if (!is.null(times) && (different_times == TRUE)){
    warning("parameter 'times' will not be used because different time vectors have
            been detected in the observations provided")
  }
  
  predicted <- list()
  if(different_times){ # if the data time points for each compartment are different, ignore times parameter and keep all the values from the data
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
    observed[[i]] <- observations[[i]][,2] #drop the column of time for each compartment and keep ony the data
  }
  
  res <- list() 
  for (i in 1:length(observed)) { # loop for each compartment
    res[[i]] <- observed[[i]] - predicted[[i]] # calculate the residuals of each compartment and store them to lists
  }
  names(res) <- names(observed)
   
  return(sum((unlist(res))^2)) # Unlist all residuals and sum their squared values
}


#RSS(predictions, observations)

#=====================
#   AICc Function
#=====================
# Akaike information criteria corrected for small sample size
# n = Number of total observations 
# k = Number of model parameters

AICc <- function(n = NULL,k, predictions, observations, times=NULL){
  
  if(is.null(n) & !is.null(times)){
    stop("The number of total observations n must be given")
  }
  
  # calculate n in case it is not given
  if(is.null(n) & is.null(times)){
    n <- 0 
    for (i in 1:length(observations)) {
      n <- n + dim(observations[[i]])[1]
    }
  }
  
  return(-2*log(RSS(predictions,observations,times)/n) + 2*k + (2*k*(k+1))/(n-k-1))
}

