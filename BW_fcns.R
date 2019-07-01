suppressMessages(library(tidyverse))

### Baum Welch Algorithm
## Inputs: training data from cdcfluview, lower and uppper bounds for variance priors, 
##         error cutoff (tolerance) value for convergence, number of random starts (i.e. iterations
##         to test convergence)
## Output: List with sublists of MLEs for params rho, trans_matrix, sigma_0, sigma_1, and initial states
baum_welch = function(train.dat, var_lbound, var_ubound, err_cutoff, n_rand_start,
                      rinit_states = TRUE){
  
  # initialize variables
  N = unique(train.dat$season)
  S = unique(train.dat$week)
  revS = rev(S)
  
  a = var_lbound
  b = var_ubound
  error = err_cutoff
  
  out_params = list()
  
  for (n in 1:n_rand_start){
    set.seed(n*var_ubound)
    
    rho = runif(1, -1,1)
    a_0.0 = rbeta(1, 0.5,0.5)
    a_1.1 = rbeta(1, 0.5, 0.5)
    
    init.state = c(0,1)
    if (rinit_states){
      pi0 = runif(1)
      init.state = c(pi0,  1-pi0)
    }
   
    P = matrix(c(a_0.0, 1- a_0.0, 1- a_1.1, a_1.1), 2, 2, byrow = TRUE)
    
    theta_l = runif(1, a, b)
    theta_m1 = runif(1, theta_l, b)
    theta_m2 = runif(1, theta_m1, b)
    theta_h = runif(1, theta_m2, b)
    
    sigma_0 = runif(1, theta_l, theta_m1)
    sigma_1 = runif(1, theta_m2, theta_h)
    
    sigma_0 = rep(sigma_0, length(N))
    sigma_1 = rep(sigma_1, length(N))
    
    cnt = 0
    logic_vec = c()
    
    while(sum(logic_vec) !=5){
      logic_vec = c()
      
      ######################
      ###### E-step #######
      #####################
      
      # Normalized Forward-Backward Algorithm
      fwd_mat = back_mat = array(NA, c(2,max(S),length(N))) # Create fwd and back matrices
      
      # Create gamma and xi arrays
      gamma  = array(NA, c(2,max(S),length(N)))
      xi  = array(NA, c(4,max(S)-1,length(N)))
      
      # Create normalization weight arrays
      eta = array(NA, c(1,max(S),length(N)))
      
      for (i in 1:length(N)){
        diff_ts = train.dat %>% filter(season == N[i]) %>% 
          select(diff_weighted) %>% unlist()
        
        # Initialize Forward-Backward Matrices
        fwd_mat[1,1, i] = init.state[1] * dnorm(diff_ts[1], 0, sd = sigma_0[i])
        fwd_mat[2,1, i] = init.state[2] * dnorm(diff_ts[1], 0, sd = sigma_1[i])
        back_mat[1,52, i] = back_mat[2,52, i] = 1
        
        eta[,1,i] = 1/sum(fwd_mat[, 1, i]) # Normlization weight at t = 1
        fwd_mat[, 1, i] = fwd_mat[, 1, i] * eta[1,1, i] # Normalized values
        
        ### Forward steps
        for (t in S[-1]){ # For each week
          fwd_mat[1, t, i] = (fwd_mat[1, t-1, i]*P[1,1] + fwd_mat[2, t-1, i]*P[2,1]) * dnorm(diff_ts[t], 0,  
                                                                                             sd = sigma_0[i]) 
          fwd_mat[2, t, i] = (fwd_mat[1, t-1, i]*P[1,2] + fwd_mat[2, t-1, i]*P[2,2]) * dnorm(diff_ts[t], 
                                                                                             rho*diff_ts[t-1], 
                                                                                             sd = sigma_1[i])
          eta[, t, i] = 1/sum(fwd_mat[, t, i])
          fwd_mat[, t, i] = fwd_mat[, t, i] * eta[, t, i] # Normalized alphas
        }
        
        
        ## Backward Steps
        for (t.1 in revS[-1]){
          back_mat[1, t.1, i] = (back_mat[1, t.1+1, i]*P[1,1]*dnorm(diff_ts[t.1+1], 0,  
                                                                    sd = sigma_0[i])) + 
            (back_mat[2, t.1+1, i]*P[1,2]*dnorm(diff_ts[t.1+1], rho*diff_ts[t.1], 
                                                sd = sigma_1[i]))
          
          back_mat[2, t.1, i] = (back_mat[1, t.1+1, i]*P[2,1]*dnorm(diff_ts[t.1+1], 0,  
                                                                    sd = sigma_0[i])) + 
            (back_mat[2, t.1+1, i]*P[2,2] * dnorm(diff_ts[t.1+1], rho*diff_ts[t.1], 
                                                  sd = sigma_1[i]))
          
          back_mat[, t.1, i] = back_mat[, t.1, i] *eta[, t.1 + 1, i] # Normalized betas
        }
        
        ## Gamma weights
        for(t in S){
          gamma[1, t, i] = fwd_mat[1, t, i]*back_mat[1, t, i]
          gamma[2, t, i] = fwd_mat[2, t, i]*back_mat[2, t, i]
          gamma[, t, i] = gamma[, t, i]/ sum(gamma[, t, i])
        }
        
        ## Xi weights
        for(t in 1:(length(S) -1)){
          # Xi starting state i = 0 to j = 0,1
          xi[1, t, i] = (gamma[1,t,i]*P[1,1]*back_mat[1, t+1, i]*eta[, t+1, i]*dnorm(diff_ts[t+1], 0,
                                                                                     sd = sigma_0[i]))/back_mat[1,t,i]
          xi[2, t, i] = (gamma[1,t,i]*P[1,2]*back_mat[2, t+1, i]*eta[, t+1, i]*dnorm(diff_ts[t+1], rho*diff_ts[t],
                                                                                     sd = sigma_1[i]))/back_mat[1,t,i]
          
          # Xi starting state i = 1 to j = 0,1
          xi[3, t, i] = (gamma[2,t,i]*P[2,1]*back_mat[1, t+1,i]*eta[, t+1, i]*dnorm(diff_ts[t+1], 0,
                                                                                    sd = sigma_0[i]))/back_mat[2,t,i]
          xi[4, t, i] = (gamma[2,t,i]*P[2,2]*back_mat[2, t+1,i]*eta[, t+1, i]*dnorm(diff_ts[t+1], rho*diff_ts[t],
                                                                                    sd = sigma_1[i]))/back_mat[2,t,i]
        }
        
      } #### End of E-Step
    
      #####################
      ###### M-step #######
      #####################
      
      ## pi_k updates
      new.init.state = c(sum(gamma[1,1,])/sum(gamma[,1,]), sum(gamma[2,1,])/sum(gamma[,1,]))
      logic_vec = c(logic_vec, all(abs(init.state - new.init.state) < error))
      
      ## Transition prob updates
      P.new = matrix(NA,2,2)
      P.new[1,1] = sum(xi[1,1:51,])/sum(xi[1:2,1:51,])
      P.new[1,2] = sum(xi[2,1:51,])/sum(xi[1:2,1:51,])
      P.new[2,1] = sum(xi[3,1:51,])/sum(xi[3:4,1:51,])
      P.new[2,2] = sum(xi[4,1:51,])/sum(xi[3:4,1:51,])
      logic_vec = c( logic_vec, all(abs(c(P - P.new)) < error))
      
      
      
      ## Rho updates
      rho.num = rho.denom = rho.year.num = rho.year.denom = c()
      for (i in 1:length(N)){
        diff_ts = train.dat %>% filter(season == N[i]) %>% 
          select(diff_weighted) %>% unlist()
        for (t in S[-1]){
          rho.num = c(rho.num, gamma[2,t,i]*diff_ts[t]*diff_ts[t-1])
          rho.denom = c(rho.denom, diff_ts[t-1]^2)
        }
        
        rho.year.num[i] = sum(rho.num)
        rho.year.denom[i] = sum(rho.denom)
      }
      
      rho.new = sum(rho.year.num)/sum(rho.year.denom)
      logic_vec = c(logic_vec, abs(c(rho - rho.new))< error)
      
      
      ## Sigma^2.0 updates
      sigma_0.new = rep(NA, length(N))
      for (i in 1:length(N)){
        diff_ts = train.dat %>% filter(season == N[i]) %>% 
          select(diff_weighted) %>% unlist()
        
        sigma_0.new[i] = sum(gamma[1,,i]*(diff_ts^2))/sum(gamma[1,,i])
      }
      
      logic_vec = c(logic_vec, all(abs(c(sigma_0 - sqrt(sigma_0.new)))< error))
      
      
      ### Sigma^2.1 updates
      sigma_1.new = rep(NA, length(N))
      for (i in 1:length(N)){
        diff_ts = train.dat %>% filter(season == N[i]) %>% 
          select(diff_weighted) %>% unlist()
        sigma_1.num = c()
        sigma_1.num = c(sum(gamma[2,,i]*(diff_ts[1]^2))) # t =1 for k = 1 no AR1, instead normal centerd at 0
        
        for (t in S[-1]){
          sigma_1.num = c(sigma_1.num, gamma[2,,i]*(diff_ts[t] - rho.new*diff_ts[t-1])^2)
        }
        sigma_1.new[i] = sum(sigma_1.num)/sum(gamma[2,,i])
      }
      logic_vec = c(logic_vec, all(abs(c(sigma_1 - sqrt(sigma_1.new)))< error))
      
    
      ## Update params based on MLEs
      init.state = new.init.state
      P = P.new
      rho = rho.new
      sigma_0 = sqrt(sigma_0.new)
      sigma_1 = sqrt(sigma_1.new)
      
      cnt = cnt +1
    }# End of M-Step
  
    
    out_params[[n]] = list(iter = cnt, st_states = init.state, 
                         rho = rho, trans_mat = c(P),
                         std_dev0 = sigma_0, std_dev1 = sigma_1,
                         gamma = gamma)
  }
  
  return(out_params)
  
}
