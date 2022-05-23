binary_test <- x <- rbinom(n=104, size=1, prob=0.5) 
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

start_time = Sys.time()
count = 0
re <- ga_fitness(runif(16,1,8.99))
end_time = Sys.time()
print(end_time-start_time)
print(re)