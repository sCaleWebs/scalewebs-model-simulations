## functions for data simulation

simulate_degrees <- function(n_basal, S, degree_sd_logit = .3){
  mean_p <- (n_basal - 1)/(S - 1)
  
  consumer_probs_logit <- plogis(qlogis(mean_p) + rnorm(S - n_basal, mean = 0, sd = degree_sd_logit))
  
  degrees <- rbinom(S - n_basal, S - 1, prob = consumer_probs_logit) + 1
  
  return(degrees)
}