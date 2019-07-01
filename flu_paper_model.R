library(rjags)


paper_model = function(dat, var_lbound, var_ubound, burn_in, mc_iter, nthin, n_chains){
  N = unique(dat$season)
  S = unique(dat$week) 
  diff.rates = comp =  lambda= matrix(NA, max(S), max(N))
  
  for (i in N){
    diff.rates[,i] = dat %>% filter(season == i) %>% select(diff_weighted) %>%  unlist()
  }
  
  
  DSmodel <- "model{
  for (j in 1:nyear){
    dif.rates[1, j] ~ dnorm(0,tau[1, j])
    tau[1, j] <- pow(lambda[comp[1, j], j],-2)
  }
  
  for (j in 1:nyear)
  {
    for (i in 2:nweek[j]) {
      dif.rates[i, j] ~ dnorm(mu[i, j],tau[i, j])
      tau[i, j] <- pow(lambda[comp[i, j], j],-2)
      
      mu[i, j] <- ro*dif.rates[i-1, j]*equals(comp[i, j],2)
    }
  }
  
  ro ~ dunif(-1,1)
  
  for (j in 1:nyear) {
    comp[1, j] ~ dcat(P0[])
    lambda[1, j] ~ dunif(linf,lmed1)
    lambda[2, j] ~ dunif(lmed2,lsup)
  }
  
  linf ~ dunif(a,b)
  lmed1 ~ dunif(linf,b)
  lmed2 ~ dunif(lmed1,b)
  lsup ~ dunif(lmed2,b)
  
  for (j in 1:nyear)
  {
    for (i in 2:nweek[j]){
    comp[i, j] ~ dcat(P.mat[comp[i-1, j], ])
    }
  }
  P0[1]<-0.5
  P0[2]<-0.5
  P.mat[1,1] ~ dbeta(0.5,0.5)
  P.mat[2,2] ~ dbeta(0.5,0.5)
  P.mat[1,2]<- 1-P.mat[1,1]
  P.mat[2,1]<- 1-P.mat[2,2]
}"

  write_file(DSmodel, "/Users/ecoronado/Documents/Spring2019/STA531/final_project/DSModel.bugs")
  
  jags <- jags.model("/Users/ecoronado/Documents/Spring2019/STA531/final_project/DSModel.bugs",
                     data = list('nyear' = max(N),
                                 'nweek' = rep(max(S),max(N)),
                                 "dif.rates" = diff.rates,
                                 a = var_lbound,
                                 b = var_ubound),
                     n.chains = n_chains, quiet = TRUE)
  
  
  update(jags, burn_in)
  
  true_parms = jags.samples(jags,
                            c('ro', 'lambda', "P.mat", "linf", "lmed1", "lmed2", "lsup",
                              "comp"),
                            mc_iter, thin = nthin)
  
  
}

