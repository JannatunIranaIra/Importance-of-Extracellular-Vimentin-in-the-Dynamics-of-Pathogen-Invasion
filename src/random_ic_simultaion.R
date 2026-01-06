#dynamics of the model for different gamma, lambda1 value (series of phase potraits try)
library(deSolve)
library(ggplot2)
library(dplyr)
library(egg)


# Generate 10,000 random initial conditions
set.seed(42)
n_sims <- 100
init_conditions <- data.frame(
  P = 1530,#runif(n_sims, 1, 100),
  V = 0.10,#runif(n_sims, 1, 100),
  Pv = 10000,#runif(n_sims, 0.1, 100),
  I = 998.1644#runif(n_sims, 0, 1000)  # kc = 1000, ensure kc - I > 0
)

# Store all results
all_results <- list()

for (i in 1:n_sims) {
  #parameter values
  parameters <- list(lambda1 = 2.5e0, #loop
                     b = 0.061,
                     delta1 = 0.01,
                     delta2 = 0.05,
                     pi = 2, 
                     k = 0.1,
                     kc = 10^3,
                     rho1 = 0.5, 
                     gamma1 = 1e-05,  #loop
                     beta1 = 1.5 *10^(-4),
                     beta2 = 10^(-4),
                     rho2 = 1.8,
                     sigma1 = 10^(-3), 
                     rhoM = 0.0024)
  
  # initial conditions
  ic <- as.numeric(init_conditions[i, ])
  state <- c(P = ic[1],
             V = ic[2],
             Pv = ic[3],
             I = ic[4])
  
  # ODE system
  Model<-function(t, state, parameters) {
    with(as.list(c(state, parameters)),{
      # rate of change
      dP <- - parameters$lambda1*P*V + parameters$delta2 *parameters$pi*I- parameters$rho1*P
      dV <- - parameters$lambda1*P*V + parameters$gamma1*(parameters$kc-I)*V+ parameters$beta2*P+ parameters$beta1*P*parameters$kc - parameters$rho2*V
      dPv <- parameters$lambda1*P*V - parameters$sigma1*Pv*(parameters$kc-I)- parameters$rhoM*Pv
      dI <- parameters$sigma1*Pv*(parameters$kc-I) - parameters$delta2 *I
      
      # return the rate of change
      list(c(dP, dV, dPv, dI))
    }) # end with(as.list ...
  }
  
  # Time span
  times <- seq(0, 1500, by = 1)
  
  result <- tryCatch({
    ode(y = state, times = times, func = Model, parms = parameters)
    }, error = function(e) NULL)
  
  if (!is.null(result)) {
    df <- as.data.frame(result)
    df$sim_id <- i
    all_results[[i]] <- df
  }


}
Equi=matrix(data=0,nrow = 100,ncol = 4 )
for (i in 1:length(all_results)) {
  #
  temp = as.matrix(all_results[[i]][dim(all_results[[i]])[1],c("P","V","Pv","I")])
  #print(length(temp))
  Equi[i,]=temp
}
Equi<-bind_rows(all_results)

#plot(Equi[])
# Combine into one data frame
trajectories <- bind_rows(all_results)

# Plot a subset (e.g. first 50) of simulations
sampled_ids <- sample(unique(trajectories$sim_id), 50)
p1 <- ggplot(filter(trajectories, sim_id %in% sampled_ids), aes(x = time)) +
  geom_line(aes(y = P, group = sim_id), alpha = 0.5, color = "steelblue") +
  labs(title = "Sampled Pathogen Trajectories",
       y = "Pathogen", x = "Time") +
  theme_minimal()

p2 <- ggplot(filter(trajectories, sim_id %in% sampled_ids), aes(x = time)) +
  geom_line(aes(y = V, group = sim_id), alpha = 0.5, color = "darkgreen") +
  labs(title = "Sampled Vimentin Trajectories",
       y = "Vimentin", x = "Time") +
  theme_minimal()

p3 <- ggplot(filter(trajectories, sim_id %in% sampled_ids), aes(x = time)) +
  geom_line(aes(y = Pv, group = sim_id), alpha = 0.5, color = "purple") +
  labs(title = "Sampled P-V Complex Trajectories",
       y = "P-V Complex", x = "Time") +
  theme_minimal()

p4 <- ggplot(filter(trajectories, sim_id %in% sampled_ids), aes(x = time)) +
  geom_line(aes(y = I, group = sim_id), alpha = 0.5, color = "red") +
  labs(title = "Sampled Infected Cell Trajectories (50 out of 10,000)",
       y = "Infected Cell", x = "Time") +
  theme_minimal()

# Print all plots
ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)

# phase potraits
 ggplot(filter(trajectories, sim_id %in% sampled_ids), aes(x = P, y = I)) +
   geom_path(arrow = arrow(type = "closed", length = unit(0.1, "inches"))) +
   geom_point(alpha = 0.3, color = "darkred") +
  labs(title = "Pathogen vs Infected phase potrait",
       y = "Infected", x = "Pathogen") +
  theme_minimal()

 ggplot(filter(trajectories, sim_id %in% sampled_ids), aes(x = P, y = V)) +
   geom_path(arrow = arrow(type = "closed", length = unit(0.1, "inches"))) +
   geom_point(alpha = 0.3, color = "green") +
   labs(title = "Pathogen vs Vimentin phase potrait",
        y = "Vimentin", x = "Pathogen") +
   theme_minimal()
 ggplot(filter(trajectories, sim_id %in% sampled_ids), aes(x = V, y = I)) +
   geom_path(arrow = arrow(type = "closed", length = unit(0.1, "inches"))) +
   geom_point(alpha = 0.3, color = "darkblue") +
   labs(title = "Vimentin vs Infected phase potrait",
        y = "Infected", x = "Vimentin") +
   theme_minimal()




