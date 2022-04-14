## functions for data simulation

simulate_degrees <- function(n_basal, S, degree_sd_logit = .3){
  mean_p <- (n_basal - 1)/(S - 1)
  
  consumer_probs_logit <- plogis(qlogis(mean_p) + rnorm(S - n_basal, mean = 0, sd = degree_sd_logit))
  
  degrees <- rbinom(S - n_basal, S - 1, prob = consumer_probs_logit) + 1
  
  return(degrees)
}

# function to take a bunch of degrees and return a matrix
simulate_network <- function(S, n_basal, degrees){
  consumers_adj <- degrees |> 
    map(~ rep(0:1, times = c(S-.x, .x))) |> 
    map(sample, replace = FALSE) |> 
    bind_cols()
  
  final_sp <- S-n_basal
  
  uneaten_resources <- rowSums(consumers_adj[,-final_sp]) == 0
  consumers_adj[,final_sp] <- as.numeric(uneaten_resources)
  
  species_names <- c(paste0("sp", 1:(S - n_basal)),paste0("Res", 1:n_basal))
  
  adj_mat <- consumers_adj |> 
    cbind(matrix(0, nrow = S, ncol = n_basal)) |> 
    set_names(species_names)

  return(adj_mat)
}
  

simulate_site_env <- function(nsites = 45, 
                              xmin = -3, xmax = 3, 
                              bar_phi = -1, 
                              beta_phi = .3){
  
  x <- runif(nsites, min = xmin, max = xmax)
  
  phi <- bar_phi + beta_phi * x
  
  return(data.frame(x, phi))
  
}

sim_ab <- function(corr_ab = .5, sd_a = .4, sd_b = 1, avg_a = 2, avg_b = -3){
  corr_mat <- matrix(c(1, corr_ab, corr_ab, 1), nrow = 2)
  sds <- c(sd_a, sd_b)
  
  sig <- diag(sds) %*% corr_mat %*% diag(sds)
  
  consumer_a_b <- MASS::mvrnorm(S - n_basal, mu = c(avg_a, avg_b), Sigma = sig)
  
  return(consumer_a_b)
}

# put togeterh in fake data
create_fake_data <- function(nsites, n_basal = 3, S, site_df, consumer_a_b){
  expand_grid(site = paste0("site", 1:nsites),
              sp = paste0("sp", 1:(S - n_basal))) |> 
    mutate(
      env_x = site_df$x[site],
      logit_prob = consumer_a_b[sp, 1] * (site_df$phi[site] - consumer_a_b[sp, 2])
    )
}


# plotting ----------------------------------------------------------------

plot_degree_dist <- function(d, S){
  tibble(degrees = d) |> 
    count(degrees) |> 
    ggplot(aes(x = degrees, y = n)) + 
    geom_bar(stat = "identity") + 
    coord_cartesian(xlim = c(0, S)) + 
    labs(x = "Degrees of nodes",
         y = "Frequency")
}


plot_network_mat <- function(adj_mat){
  adj_mat |> 
    mutate(prey_name = species_names) |> 
    pivot_longer(-prey_name, values_to = "interact", names_to = "consumer_name") |> 
    mutate(consumer_name = forcats::fct_reorder(consumer_name, interact, .fun = sum),
           prey_name = factor(prey_name, levels = levels(consumer_name))) |> 
    ggplot(aes(x = prey_name, y = consumer_name, fill = as.factor(interact))) + 
    geom_tile() + 
    scale_fill_manual(values = c("white", "black")) + 
    guides(fill = "none") + 
    coord_fixed()
}