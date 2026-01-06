# Necessary libraries
library(polynom)
library(ggplot2)
library(dplyr)
library(deSolve)
library(egg)

# Phase portrait investigation
# Function for solving nonzero, non-negative and real I^* 
find_Equilibruim <- function(A_coeff,B_coeff,C_coeff,D_coeff, tolerance = 1e-5){
  
  # Compute the discriminant
  discriminant <- 18 * A_coeff * B_coeff * C_coeff * D_coeff - 4 * B_coeff^3 * D_coeff + B_coeff^2 * C_coeff^2 - 4 * A_coeff * C_coeff^3 - 27 * A_coeff^2 * D_coeff^2
  
  poly <- polynomial(c(D_coeff,C_coeff,B_coeff,A_coeff)) # replaced polyroot function with polynomial+solve command
  solution = solve(poly)
  
  # avoiding very small imaginary parts
  for (ii1 in 1:length(solution)) {
    if (abs(Im(solution[ii1])) < tolerance) {
      solution[ii1] <- Re(solution[ii1])
    }
  }
  if (discriminant > 0) {
    # print("All roots are real and distinct.")
    solution=solution[which(Re(solution)>0)]
    return(solution)
  } else if (discriminant == 0) {
    # print("At least two roots are real and repeated.")
    solution <- solution[!duplicated(Re(solution))]
    solution  = solution[which(Re(solution)>0)]
    return(solution)
  } else {
    # print("One root is real, and two are complex conjugates.")
    solution <- as.complex(solution)
    solution <- solution[!(solution %in% solution[duplicated(Re(solution))])]
    solution <- solution[which(Im(solution)==0)]
    solution=solution[which(Re(solution)>0)]
    return(solution)
  }
}

# output results variables
para_value <- c()
equi_ife <- c()
equi_ie <- c()
nature <- c()
count <- 0
eq_points <- data.frame(Para=numeric(), P_star=numeric(), V_star=numeric(), Pv_star=numeric(), equi_I=numeric())

#5e6, 5.1e6, 0.44, 0.43, .00079, .00081
critical_point <- c(0.70,   #stable ife near switch point
                    1.15,     #unstable ife near switch point
                    5.0,
                    6.4,
                    6.6,
                    7.0,
                    11.01)     #unstable ife

# Loop to simulate equilibria and analyze
for (j in 1:length(critical_point)){
  #parameter values
  parameters <- list(lambda1 =2.5e0, 
                     b = 0.061,
                     delta1 = 0.01,
                     delta2 = 10^3,
                     pi = critical_point[j], 
                     k = 47,
                     kc = 10^3,
                     rho1 = 0.8, 
                     gamma1 = 2e-02, # gamma = 2e-03 te sob stable r gamma = 2e-01 a unstable stable milano
                     beta1 = 1.5 *10^(-4),
                     beta2 = 10^(-4),
                     rho2 = 2.2,
                     sigma1 = 10^(-3), 
                     rhoM = 0.0024)
  # # Some expressions used for Equilibrium values
  # L <- (parameters$beta1 * parameters$kc + parameters$beta2) / parameters$rho1
  # TT <- parameters$k * parameters$kc + parameters$delta2 - (parameters$b - parameters$delta1)
  # Tb <- parameters$k * parameters$kc + parameters$delta2 + parameters$delta1
  #
  # Coefficients
  A_coeff <- with(parameters, {
    sigma1 * ((1 - pi) * lambda1 * delta2 * (1 + ((1 - pi) / rho1) * (beta2 + beta1 * kc)) + rho1 * sigma1 * gamma1)
  })
  
  B_coeff <- with(parameters, {
    -lambda1 * delta2 * (
      (1 - pi) * (2 * kc * sigma1 + rhoM) +
        rhoM * sigma1 +
        ((1 - pi) / rho1) * (beta2 + beta1 * kc) * ((2 * kc * sigma1 * (1 - pi) + rhoM) + rhoM * sigma1)
    ) -
      rho1 * sigma1 * (
        (gamma1 * kc - rho2) * sigma1 + gamma1 * (2 * sigma1 * kc + rhoM)
      )
  })
  
  C_coeff <- with(parameters, {
    lambda1 * delta2 * (kc * (1 - pi) + rhoM) * (
      (sigma1 * kc + rhoM) +
        (1 / rho1) * (beta2 + beta1 * kc) * (sigma1 * kc * (1 - pi) + rhoM)
    ) +
      rho1 * sigma1 * (
        (gamma1 * kc - rho2) * (2 * sigma1 * kc + rhoM) +
          kc * gamma1 * (sigma1 + rhoM)
      )
  })
  
  D_coeff <- with(parameters, {
    -rho1 * sigma1 * kc * (sigma1 * kc + rhoM) * (gamma1 * kc - rho2)
  })
  # Find equilibrium points
  I_star <- find_Equilibruim(A_coeff, B_coeff, C_coeff, D_coeff, tolerance = 1e-5)
  I_star <- c(0,I_star)
  equi_I <- c()
  # Check equilibrium values before plotting
  for (j1 in 1:length(I_star)) {
    equi_I <- Re(I_star[j1])
    # Compute equilibrium points
    # Pv_star <- (-(parameters$b - parameters$delta1 - parameters$delta2 - parameters$k * parameters$kc) * equi_I) /
    # (parameters$sigma1 * (parameters$kc - equi_I))
    Pv_star <- with(parameters, { (delta2 * equi_I) / (sigma1 * (kc - equi_I)) })
    
    # P_star <- equi_I / (parameters$rho1 * parameters$sigma1 * (parameters$kc - equi_I)) *
    #   ((parameters$sigma1 * (parameters$kc - equi_I) + parameters$rhoM) * parameters$b +
    #      (parameters$delta1 + parameters$delta2 + parameters$k * parameters$kc) *
    #      (parameters$sigma1 * (parameters$kc - equi_I) * (parameters$pi - 1) - parameters$rhoM))
    P_star <- with(parameters, {  (-delta2 * equi_I / (rho1 * sigma1 * (kc - equi_I))) *
        ((1 - pi) * sigma1 * (kc - equi_I) + rhoM) })
    
    # V_star <- ((parameters$beta2 + parameters$beta1 * parameters$kc) * equi_I) /
    #   (parameters$rho1 * parameters$sigma1 * (parameters$kc - equi_I) *
    #      (parameters$lambda1 * P_star + parameters$rho2 - parameters$gamma1 * (parameters$kc - equi_I))) *
    #   ((parameters$sigma1 * (parameters$kc - equi_I) + parameters$rhoM) * parameters$b +
    #      (parameters$delta1 + parameters$delta2 + parameters$k * parameters$kc) *
    #      (parameters$sigma1 * (parameters$kc - equi_I) * (parameters$pi - 1) - parameters$rhoM))
    numerator <- with(parameters, { -delta2 * equi_I * (1 + (rhoM / (sigma1 * (kc - equi_I)))) +
        (beta2 + beta1 * kc) * P_star })
    denominator <- with(parameters, { (gamma1 * kc - rho2 - gamma1 * equi_I) })
    V_star <- numerator / denominator
    
    if (any(is.nan(P_star) | is.nan(V_star) | is.nan(Pv_star) | is.nan(equi_I))) {
      next
    } else {
      if (equi_I>= 0 & Pv_star >= 0 & P_star>= 0 & V_star>=0) {
        count <- count +1
        eq_points <- rbind(eq_points, c(critical_point[j], P_star, V_star, Pv_star, equi_I))
        para_value[count] <- critical_point[j]
        if (equi_I == 0){
          equi_ife[count] <- equi_I
          equi_ie[count] <- NA
          # Calculate jacobian at IFE
          agg_jacobian = matrix(c(-parameters$rho1, 0, 0, (parameters$delta2)*parameters$pi,
                                  parameters$beta2+ parameters$beta1*parameters$kc,  parameters$gamma1 * parameters$kc-parameters$rho2, 0, 0,
                                  0, 0, -parameters$sigma1*(parameters$kc)-parameters$rhoM, 0,
                                  0,0, parameters$sigma1*(parameters$kc), -parameters$delta2),
                                nrow = 4, ncol = 4, byrow = TRUE )
          
        } else {
          equi_ie[count] <- equi_I
          equi_ife[count] <- NA
          #Calculate jacobian at IE
          agg_jacobian = matrix(c(-parameters$lambda1*V_star-parameters$rho1, -parameters$lambda1*P_star, 0, (parameters$delta2)*parameters$pi,
                                  -parameters$lambda1*V_star+ parameters$beta2+ parameters$beta1*parameters$kc, -parameters$lambda1*P_star+ parameters$gamma1*(parameters$kc-equi_ie[count])-parameters$rho2, 0, -parameters$gamma1*V_star,
                                  parameters$lambda1*V_star, parameters$lambda1*P_star, -parameters$sigma1*(parameters$kc-equi_ie[count])-parameters$rhoM, parameters$sigma1*Pv_star,
                                  0,0, parameters$sigma1*(parameters$kc-equi_ie[count]), -parameters$sigma1*Pv_star-parameters$delta2),
                                nrow = 4, ncol = 4, byrow = TRUE )
        }
        ev <- eigen(agg_jacobian)
        evalues <- ev$values
        if (all(Re(evalues) <0)){
          nature[count] <- "stable"
        } else {
          nature[count] <- "unstable"
        }
      }
    }
  }
}


for (j in 1:length(critical_point)){
  #parameter values
  parameters <- list(lambda1 =2.5e0, 
                     b = 0.061,
                     delta1 = 0.01,
                     delta2 = 10^3,
                     pi = critical_point[j], 
                     k = 47,
                     kc = 10^3,
                     rho1 = 0.8, 
                     gamma1 = 2e-02, # gamma = 2e-03 te sob stable r gamma = 2e-01 a unstable stable milano
                     beta1 = 1.5 *10^(-4),
                     beta2 = 10^(-4),
                     rho2 = 2.2,
                     sigma1 = 10^(-3), 
                     rhoM = 0.0024)
  
  
  #Solving ODEs
  state <- c(P = 1,
             V = 1e-1,
             Pv = 0,
             I = 0)
  
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
  
  times <- seq(0, 150, by = 0.1)
  out <- ode(y = state, times = times, func = Model, parms = parameters) #, atol = 1e-6, rtol = 1e-6, hmin = 1e-6, method = "lsoda"
  # Convert ODE output to data frame
  out_df <- as.data.frame(out)
  
  match_indices <- which(eq_points[,1] == critical_point[j])
  
  pdata_points <- c()
  vdata_points <- c()
  pvdata_points <- c()
  idata_points <- c()
  if (length(match_indices) > 0) {
    pdata_points <- eq_points[match_indices, 2]
    vdata_points <- eq_points[match_indices, 3]
    pvdata_points <- eq_points[match_indices, 4]
    idata_points <- eq_points[match_indices, 5]
    
    # Create a data frame for ggplot2
    eq_data <- data.frame(P = pdata_points, V = vdata_points, Pv = pvdata_points, I = idata_points)
  } else {
    eq_data <- data.frame(P = numeric(0), V = numeric(0), Pv = numeric(0), I = numeric(0))
  }
  
  
  # 2D phase portraits with ggplot2
  p1 <- ggplot(out_df, aes(x = P, y = V)) +
    geom_path(arrow = arrow(type = "closed", length = unit(0.1, "inches"))) +
    geom_point(data = eq_data, aes(x = P, y = V), color = "red", size = 3, shape = 16) +
    # xlim(0,max(1.1 * eq_data$P)) +
    # ylim(0,max(1.1 * eq_data$V)) +
    labs(x = "Pathogen (P)", y = "Vimentin (V)") +
    theme_minimal()
  
  p2 <- ggplot(out_df, aes(x = P, y = Pv)) +
    geom_path(arrow = arrow(type = "closed", length = unit(0.1, "inches"))) +
    geom_point(data = eq_data, aes(x = P, y = Pv), color = "red", size = 3, shape = 16) +
    # xlim(0,max(1.1 * eq_data$P)) +
    # ylim(0,max(1.1 * eq_data$Pv)) +
    labs(x = "Pathogen (P)", y = "P-V complex") +
    theme_minimal()
  
  p3 <- ggplot(out_df, aes(x = P, y = I)) +
    geom_path(arrow = arrow(type = "closed", length = unit(0.1, "inches"))) +
    geom_point(data = eq_data, aes(x = P, y = I), color = "red", size = 3, shape = 16) +
    # xlim(0,max(1.1 * eq_data$P)) +
    # ylim(0,max(1.1 * eq_data$I)) +
    labs(x = "Pathogen (P)", y = "Infected Cells (I)") +
    theme_minimal()
  
  p4 <- ggplot(out_df, aes(x = V, y = Pv)) +
    geom_path(arrow = arrow(type = "closed", length = unit(0.1, "inches"))) +
    geom_point(data = eq_data, aes(x = V, y = Pv), color = "red", size = 3, shape = 16) +
    # xlim(0,max(1.1 * eq_data$V)) +
    # ylim(0,max(1.1 * eq_data$Pv)) +
    labs(x = "Vimentin (V)", y = "P-V complex") +
    theme_minimal()
  
  p5 <- ggplot(out_df, aes(x = V, y = I)) +
    geom_path(arrow = arrow(type = "closed", length = unit(0.1, "inches"))) +
    geom_point(data = eq_data, aes(x = V, y = I), color = "red", size = 3, shape = 16) +
    # xlim(0,max(1.1 * eq_data$V)) +
    # ylim(0,max(1.1 * eq_data$I)) +
    labs(x = "Vimentin (V)", y = "Infected Cells (I)") +
    theme_minimal()
  
  p6 <- ggplot(out_df, aes(x = Pv, y = I)) +
    geom_path(arrow = arrow(type = "closed", length = unit(0.1, "inches"))) +
    geom_point(data = eq_data, aes(x = Pv, y = I), color = "red", size = 3, shape = 16) +
    # xlim(0,max(1.1 * eq_data$Pv)) +
    # ylim(0,max(1.1 * eq_data$I)) +
    labs(x = "P-V complex", y = "Infected Cells (I)") +
    theme_minimal()
  
  # Print all plots
  ggarrange(p1, p2, p3, p4, p5, p6, ncol = 3, nrow = 2)
}

# #3D plots
library(plotly)
# 
# critical_point <- c(5e-4,     # stable ife
#                     0.0017,   #stable ife near switch point
#                     0.0019,     #unstable ife near switch point
#                     500,
#                     45000,
#                     46000,
#                     100000)     #unstable ife
# 
for (j in 1:length(critical_point)){
  #parameter values
  parameters <- list(lambda1 =2.5e0, 
                     b = 0.061,
                     delta1 = 0.01,
                     delta2 = 10^3,
                     pi = critical_point[j], 
                     k = 47,
                     kc = 10^3,
                     rho1 = 0.8, 
                     gamma1 = 2e-02, # gamma = 2e-03 te sob stable r gamma = 2e-01 a unstable stable milano
                     beta1 = 1.5 *10^(-4),
                     beta2 = 10^(-4),
                     rho2 = 2.2,
                     sigma1 = 10^(-3), 
                     rhoM = 0.0024)


  #Solving ODEs
  state <- c(P = 1,
             V = 1e-1,
             Pv = 0,
             I = 0)

  Model<-function(t, state, parameters) {
    with(as.list(c(state, parameters)),{
      # rate of change
      dP <- - parameters$lambda1*P*V + (parameters$delta1 + parameters$delta2 + parameters$k* parameters$kc)*parameters$pi*I- parameters$rho1*P
      dV <- - parameters$lambda1*P*V + parameters$gamma1*(parameters$kc-I)*V+ parameters$beta2*P+ parameters$beta1*P*parameters$kc - parameters$rho2*V
      dPv <- parameters$lambda1*P*V - parameters$sigma1*Pv*(parameters$kc-I)- parameters$rhoM*Pv
      dI <- parameters$sigma1*Pv*(parameters$kc-I) +(- parameters$delta2+ parameters$b- parameters$delta1- parameters$k* parameters$kc)*I

      # return the rate of change
      list(c(dP, dV, dPv, dI))
    }) # end with(as.list ...
  }

  times <- seq(0, 100, by = 0.001)
  out <- ode(y = state, times = times, func = Model, parms = parameters, atol = 1e-6, rtol = 1e-6, hmin = 1e-6, method = "lsoda")
  # Convert ODE output to data frame
  out_df <- as.data.frame(out)

  match_indices <- which(eq_points[,1] == critical_point[j])

  pdata_points <- c()
  vdata_points <- c()
  pvdata_points <- c()
  idata_points <- c()
  if (length(match_indices) > 0) {
    pdata_points <- eq_points[match_indices, 2]
    vdata_points <- eq_points[match_indices, 3]
    pvdata_points <- eq_points[match_indices, 4]
    idata_points <- eq_points[match_indices, 5]

    # Create a data frame for ggplot2
    eq_data <- data.frame(P = pdata_points, V = vdata_points, Pv = pvdata_points, I = idata_points)
  } else {
    eq_data <- data.frame(P = numeric(0), V = numeric(0), Pv = numeric(0), I = numeric(0))
  }


  f1 <- plot_ly(out_df, x = ~P, y = ~V, z = ~Pv, type = "scatter3d", mode = "lines",
                line = list(width = 3, color = "blue")) %>%
    add_trace(x = eq_data$P, y = eq_data$V, z = eq_data$Pv,
              type = 'scatter3d', mode = 'markers',
              marker = list(size = 5, color = 'red', opacity = 0.8)) %>%
    layout(scene = list(
      xaxis = list(title = "Pathogen (P)"),
      yaxis = list(title = "Vimentin (V)"),
      zaxis = list(title = "P-V complex (Pv)")))

  f2 <- plot_ly(out_df, x = ~P, y = ~V, z = ~I, type = "scatter3d", mode = "lines",
                line = list(width = 3, color = "blue")) %>%
    add_trace(x = eq_data$P, y = eq_data$V, z = eq_data$I,
              type = 'scatter3d', mode = 'markers',
              marker = list(size = 5, color = 'red', opacity = 0.8)) %>%
    layout(scene = list(
      xaxis = list(title = "Pathogen (P)"),
      yaxis = list(title = "Vimentin (V)"),
      zaxis = list(title = "Infected Cells (I)")))

  f3 <- plot_ly(out_df, x = ~V, y = ~Pv, z = ~I, type = "scatter3d", mode = "lines",
                line = list(width = 3, color = "blue")) %>%
    add_trace(x = eq_data$V, y = eq_data$Pv, z = eq_data$I,
              type = 'scatter3d', mode = 'markers',
              marker = list(size = 5, color = 'red', opacity = 0.8)) %>%
    layout(scene = list(
      xaxis = list(title = "Vimentin (V)"),
      yaxis = list(title = "P-V complex (Pv)"),
      zaxis = list(title = "Infected Cells (I)")))

  f1
  f2
  f3
}



