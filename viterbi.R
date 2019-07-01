suppressMessages(library(tidyverse))


viterbi_decoding = function(test.dat, start_states, trans_mat, std_dev0, std_dev1, rho){
  
  N = unique(test.dat$season)
  S = unique(test.dat$week)
  
  fwd_mat = array(NA, c(2,max(S),length(N))) # Create fwd  matrices
  viterb = array(NA, c(2,max(S),length(N)))
  viterb[,1,] = c(1,1)
  decode_path = list()
  
  new.init.state = start_states %>% as.vector()
  P = matrix(trans_mat, 2, 2, byrow = TRUE)
  sigma_0 = std_dev0
  sigma_1 = std_dev1
  
  # For each year
  for (i in 1:length(N)){
    
    diff_ts = test.dat %>% filter(season == N[i]) %>% 
      select(diff_weighted) %>% unlist()
    
    # Initialize Forward-Backward Matrices
    fwd_mat[1,1, i] = log(new.init.state[1]) + dnorm(diff_ts[1], 0, sd = sigma_0[i], log = T)
    fwd_mat[2,1, i] = log(new.init.state[2]) + dnorm(diff_ts[1], 0, sd = sigma_1[i], log = T)
    
    ### Forward steps
    # For each week
    for (t in S[-1]){
      
      tmp_t0 = c(fwd_mat[1, t-1, i] + log(P[1,1]), fwd_mat[2, t-1, i] + log(P[2,1]))
      fwd_mat[1, t, i] = max(tmp_t0) + dnorm(diff_ts[t], 0,sd = sigma_0[i], log = T)
      viterb[1,t,i] = which(tmp_t0 == max(tmp_t0))
      
      tmp_t1 = c(fwd_mat[1, t-1, i] + log(P[1,2]),fwd_mat[2, t-1, i] + log(P[2,2]))
      
      fwd_mat[2, t, i] = max(tmp_t1) +  dnorm(diff_ts[t], rho*diff_ts[t-1], sd = sigma_1[i], log = T)
      viterb[2,t,i] = which(tmp_t1 == max(tmp_t1))
      
    }
    
    path_id = which(fwd_mat[,max(S),i] == max(fwd_mat[,max(S),i]) )
    decode_path[[i]] = data_frame(states = viterb[path_id,,i] - 1,
                                  p_st0 = fwd_mat[1,,i],
                                  p_st1 = fwd_mat[2,,i])
    
  }
  return(decode_path)
  
}