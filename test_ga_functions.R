binary_test <- x <- rbinom(n=104, size=1, prob=0.5) 
print(binary_test)
print(decode_ga(binary_test))

par <- decode_ga(binary_test)
P_groups <- length(unique(par[1:13]))  # sample size
X_groups <- length(unique(par[14:length(par)]))  # sample size

c <- par[2]  # acceptance number

-----------------------------------------------------------------------------------------------
  
  myf <- function(parameters){
    with(as.list(parameters),{
      with(as.list(pars),{
        
      print(grouping)
    })
    })
    }