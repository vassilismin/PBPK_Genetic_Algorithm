chromosome <- rbinom(n=64, size=1, prob=0.5) 
print(binary_test)
print(decode_ga(binary_test))

par <- decode_ga(binary_test)
P_groups <- length(unique(par[1:13]))  # sample size
X_groups <- length(unique(par[14:length(par)]))  # sample size

c <- par[2]  # acceptance number

-----------------------------------------------------------------------------------------------
  
  myf <- function(...){
    dots <- list(...)
    
    with(as.list(dots),{
      print(a)
    })
  }

k <- c(0,1,0,1)	
print(GA::binary2decimal(GA::gray2binary(k)))
print(GA::binary2decimal(k))


data_list <- list("1" = runif(16,1,8.99),
                  "2" = runif(16,1,8.99),
                  "3" = runif(16,1,8.99),
                  "4" = runif(16,1,8.99),
                  "5" = runif(16,1,8.99),
                  "6" = runif(16,1,8.99),
                  "7" = runif(16,1,8.99),
                   "8" = runif(16,1,8.99))

start_time = Sys.time()
clus <- parallel::makeCluster(parallel::detectCores())
parallel::clusterExport(clus,"ga_fitness")
re <-parallel::parLapply(clus,data_list,ga_fitness)
end_time = Sys.time()
print(end_time-start_time)
# Close cluster
parallel::stopCluster(clus)


