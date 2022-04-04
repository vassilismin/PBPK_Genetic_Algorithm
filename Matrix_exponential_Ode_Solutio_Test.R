library(expm)
library(deSolve)

Solve_exp_matrix <- function(x, time, y_init){
  
  
  y_t  <- matrix(data=NA, nrow = nrow(x), ncol = length(time))
  colnames(y_t) <- as.character(time)
  
  y_t[,1] <- y_init
  for (t in 2:length(time)) {
    solution_t <- expm(x*time[t])%*%y_init
    y_t[,t] <- solution_t
  }
  
  return(y_t)
}

# Test the function over a simple 2x2 system 
# dx/dt = -0.5x + 0.5y
# dy/dt = 0.5x - 0.5y

A <- matrix(c(-0.5, 0.5, 0.5 ,-0.5), nrow = 2, byrow = TRUE)
inits <- c(0.1, 0.2)
time <- seq(0,10)

exp_solution <- Solve_exp_matrix(A, time, inits)
t(exp_solution)

### Numerical Solution

params <- c(A)
inits_ode <- c("x"=0.1, "y"=0.2)

ode.func <-  function(time, inits, params){
  with(as.list(c(inits, params)),{
    dx <- params[1]*x + params[2]*y
    
    dy <- params[3]*x + params[4]*y
    
    list(c(dx,dy))
  })
}
ode(times=time, func=ode.func, y=inits_ode, parms=params, 
    method="bdf", rtol=1e-5, atol=1e-5)
