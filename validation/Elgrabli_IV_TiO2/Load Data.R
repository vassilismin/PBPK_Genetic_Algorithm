# Set the directory of the excel file with the Elgrabli data
setwd("C:/Users/vassi/Documents/LAB")

# Load the mean values 
# The units of the values are micro-g Ti per organ
mean_values <- openxlsx::read.xlsx("Elgrabli_IV_TiO2.xlsx", sheet = "Mean", colNames = T, rowNames = T)

# Take the Urine data
Urine_mean <- mean_values["Urine"]
# Remove the urine data from the main dataframe
mean_values <- mean_values[,- 9]

# Take the control values of eah compartment.
control_values <- mean_values[1,]
# Remove the control values row from the main dataframe
mean_values <- mean_values[-1,]

# Set the rownames of the main dataframe
# (They are the time points of observations given in hours)
rownames(mean_values) <- c("0.1667", "1", "24", "168", "672", "1344") 

# Subtract the control values from each compartment and time point
for (i in 1:dim(mean_values)[1]) {
  mean_values[i,] <- mean_values[i,] - control_values
}
# Transform all negative values to zeros
mean_values[mean_values < 0] <- 0

# Replace the Not-Detected values with the half of the minimum value of the dataset
mean_values[mean_values == "ND"] <- min(mean_values, na.rm = T)/2

# Load the sd values
sd_values <- openxlsx::read.xlsx("Elgrabli_IV_TiO2.xlsx", sheet = "SD", colNames = T, rowNames = T)

# Take the Urine data
Urine_sd <- sd_values["Urine"]
# Remove the urine data from the main dataframe
sd_values <- sd_values[,- 9]

# Remove the control values row from the main dataframe
sd_values <- sd_values[-1,]

# Set the rownames of the sd dataframe
# (They are the time points of observations given in hours)
rownames(sd_values) <- c("0.1667", "1", "24", "168", "672", "1344") 


# Create an extra df for the urine data
# This values are the absolute Ti measured in urine at 
# each time point (not cumulative amount).
Urine <- data.frame(cbind(Urine_mean, Urine_sd))
# Subtract the control value
urine_control <- Urine[1,1]
for(i in 1:dim(Urine)[1]){
  Urine[i,1] <- Urine[i,1] - urine_control 
}
# Drop the control row
Urine <- Urine[-c(1,6,7),]

# Calculate the cumulative amount of Ti in urine 
Urine[,1] <- cumsum(Urine[,1])
# Do we need to do the same on sd values ??????

# Final step: Transform the Ti to TiO2 mass.
frac_Ti <- 59.9/100 #fraction g of Ti in g of TiO2 (given in g)
mean_values <- mean_values/frac_Ti
sd_values <- sd_values/frac_Ti
Urine <- Urine/frac_Ti 
rownames(Urine) <- c("0.1667", "1", "24", "168") 
colnames(Urine) <- c("Mean", "SD")


###########################
Body_mass <- 200 # g
Dose_per_rat <- 1.7 # mg

data_list <- list("Mean" = mean_values,
                  "SD" = sd_values,
                  "Urine" = Urine,
                  "Body_mass" = Body_mass,
                  "Dose_per_rat" = Dose_per_rat)

data_list