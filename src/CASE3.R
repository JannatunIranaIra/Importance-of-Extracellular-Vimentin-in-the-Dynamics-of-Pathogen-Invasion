# Necessary libraries
library(polynom)
library(ggplot2)
library(dplyr)
library(deSolve)
library(egg)

#generating parameter set
bodyframe <-seq(from=0,to=1,length.out= 1001)
# m8 = 1e-7*bodyframe
# m7 = 1e-6*bodyframe
# m6 = 1e-5*bodyframe
# m5 = 1e-4*bodyframe
m4 = 1e-3*bodyframe
m3 = 1e-2*bodyframe
m2 = 1e-1*bodyframe
m1 = bodyframe
p0 = 10*bodyframe
p1 = 10^2*bodyframe
p2 = 10^3*bodyframe
p3 = 10^4*bodyframe
# p4 = 10^5*bodyframe
# p5 = 10^6*bodyframe
# p6 = 10^7*bodyframe
all1 = unique(c(m4,m3,m2,m1,p0,p1,p2,p3))
all1=all1[2:length(all1)] #gamma1



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
equi_P <- c()
equi_V <- c()
equi_Pv <- c()
nature <- c()
count_comp <- c()
resilience_vinf <- c()
reactivity_v0 <- c()
count_re0 <- c()
eig1 <- c()
eig2 <- c()
eig3 <- c()
eig4 <- c()
count <- 0

# Loop to simulate equilibria and analyze
for (j in 1:length(all1)){
  #parameter values
  parameters <- list(lambda1 =2.5e0, 
                     b = 0.061,
                     delta1 = 0.01,
                     delta2 = 10^4,
                     pi = 10, 
                     k = 47,
                     kc = 10^3,
                     rho1 = 0.5, 
                     gamma1 = all1[j], 
                     beta1 = 1.5 *10^(-4),
                     beta2 = 10^(-4),
                     rho2 = 0.8,
                     sigma1 = 10^(-3), 
                     rhoM = 0.0024)
  
  
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
  I_star <- find_Equilibruim(A_coeff,B_coeff,C_coeff,D_coeff, tolerance = 1e-5)
  I_star <- c(0,I_star)
  equi_I <- c()
  # Check equilibrium values before plotting
  for (j1 in 1:length(I_star)) {
    
    equi_I <- Re(I_star[j1])
    
    # Compute equilibrium points
    Pv_star <- with(parameters, { (delta2 * equi_I) / (sigma1 * (kc - equi_I)) })
    
    P_star <- with(parameters, {  (-delta2 * equi_I / (rho1 * sigma1 * (kc - equi_I))) *
        ((1 - pi) * sigma1 * (kc - equi_I) + rhoM) })
    
    numerator <- with(parameters, { -delta2 * equi_I * (1 + (rhoM / (sigma1 * (kc - equi_I)))) +
        (beta2 + beta1 * kc) * P_star })
    denominator <- with(parameters, { (gamma1 * kc - rho2 - gamma1 * equi_I) })
    V_star <- numerator / denominator
    
    # Debugging - Check values
    # cat("Equilibrium", j1, ":", "P_star =", P_star, "V_star =", V_star, "Pv_star =", Pv_star, "I_star =", equi_I, "\n")
    if (any(is.nan(P_star) | is.nan(V_star) | is.nan(Pv_star) | is.nan(equi_I))) {
      next
    } else {
      if (equi_I>= 0 & Pv_star >= 0 & P_star>= 0 & V_star>=0) {
        count <- count +1
        para_value[count] <- all1[j]
        equi_Pv[count] <- Pv_star
        equi_V[count] <- V_star
        equi_P[count] <- P_star
        
        if (equi_I == 0){
          equi_ife[count] <- equi_I
          equi_ie[count] <- NA
        } else {
          equi_ie[count] <- equi_I
          equi_ife[count] <- NA
        }
        #Calculate jacobian at IE
        agg_jacobian = matrix(c(-parameters$lambda1*V_star-parameters$rho1, -parameters$lambda1*P_star, 0, (parameters$delta2)*parameters$pi,
                                -parameters$lambda1*V_star+ parameters$beta2+ parameters$beta1*parameters$kc, -parameters$lambda1*P_star+ parameters$gamma1*(parameters$kc-equi_I)-parameters$rho2, 0, -parameters$gamma1*V_star,
                                parameters$lambda1*V_star, parameters$lambda1*P_star, -parameters$sigma1*(parameters$kc-equi_I)-parameters$rhoM, parameters$sigma1*Pv_star,
                                0,0, parameters$sigma1*(parameters$kc-equi_I), -parameters$sigma1*Pv_star-parameters$delta2),
                              nrow = 4, ncol = 4, byrow = TRUE )
        evalues <- eigen(agg_jacobian)$values
        eig1[count] <- evalues[1]
        eig2[count] <- evalues[2]
        eig3[count] <- evalues[3]
        eig4[count] <- evalues[4]
        
        # Number of complex eigenvalues
        temp_comp <- 0
        for (i3 in 1:length(evalues)){
          if (Im(evalues[i3]) != 0){
            temp_comp <- temp_comp+1
          }
        }
        count_comp[count] <- temp_comp
        
        if (all(Re(evalues) <0)){
          nature[count] <- "stable"
          eig_dominant <- max(abs(evalues))
          # Resilience of the system
          resilience_vinf[count] <- - Re(eig_dominant)
          
          # Hermitian for the Reactivity
          Hermi = 0.5 * (agg_jacobian + t(agg_jacobian) )
          ev_H <- eigen(Hermi)
          evalues_H <- ev_H$values
          dom_H <- - max(evalues_H)
          # Reactivity of the system
          temp_v0 <- dom_H
          for (i4 in 1:length(evalues_H)){
            if (dom_H == -evalues_H[i4] ){
              temp_v0 <- evalues_H[i4]
            }
          }
          reactivity_v0[count] <- temp_v0
          
          count_re0[count] <- NA
          
        } else {
          nature[count] <- "unstable"
          resilience_vinf[count] <- NA
          reactivity_v0[count] <- NA
          
          # counting the number of eigenvalues where Re(eigenvalue)=0
          unstable_evalues <- evalues
          temp_r0 <- 0
          for (i5 in 1:length(unstable_evalues)) {
            if(abs(Re(unstable_evalues[i5])) < 10^(-6) ){
              temp_r0 <- temp_r0 +1
              critical_point <- para_value[count]
              critical_equi <- equi_I
              print(critical_point)
            }
          }
          count_re0[count]<- temp_r0
        }
        
        
      }
    }
  }
}

# Data frame
df1 <- data.frame(
  parameter = para_value,
  IFE = equi_ife,
  IE = equi_ie,
  P = equi_P,
  V = equi_V,
  Pv = equi_Pv,
  nature = nature,
  complex_eigenvalues = count_comp,
  resilience = resilience_vinf,
  reactivity = reactivity_v0,
  real_eigen_zero = count_re0,
  eigenvalue1  = eig1,
  eigenvalue2 = eig2,
  eigenvalue3 = eig3,
  eigenvalue4 = eig4
)

df <- data.frame(
  parameter = para_value,
  Equ = c(equi_ife, equi_ie),
  nature = nature,
  complex_eigenvalues = count_comp,
  resilience = resilience_vinf,
  reactivity = reactivity_v0,
  real_eigen_zero = count_re0,
  eigenvalue1  = eig1,
  eigenvalue2 = eig2,
  eigenvalue3 = eig3,
  eigenvalue4 = eig4
)

# two separate bifurcation diagram
# Diagram 1
p1=ggplot(data = df1, aes(x = (parameter),
                          y = (IFE),
                          color=nature,
                          # size = factor(nature),
                          shape = factor(complex_eigenvalues))) +
  geom_point() +
  # geom_point(aes(x= critical_point,y= critical_equi),colour="green", size = 3.5) +          # Equilibrium points
  # geom_line(aes(group = 1),       # Line connecting points
  # linetype = "dashed", color = "gray")  +
  scale_color_manual(values = c(stable = "blue",
                                unstable = "red"))+
  labs(title = "Bifurcation Diagram",
       x = "Parameter",
       y = "Equilibrium")+
  # scale_y_log10()  +
  scale_x_log10()  +
  theme_minimal()
print(p1)

# Diagram 2
p2=ggplot(data = df1, aes(x = (parameter),
                          y = (IE),
                          color=nature,
                          size = factor(complex_eigenvalues),
                          shape = factor(nature))) +
  geom_point() +
  # geom_point(aes(x= critical_point,y= critical_equi),colour="green", size = 3.5) +          # Equilibrium points
  # geom_line(aes(group = 1),       # Line connecting points
  # linetype = "dashed", color = "gray")  +
  scale_color_manual(values = c(stable = "blue",
                                unstable = "red"))+
  labs(title = "Bifurcation Diagram",
       x = "Parameter",
       y = "Equilibrium")+
  # scale_y_log10()  +
  scale_x_log10()  +
  theme_minimal()
print(p2)

# # All together bifurcation diagram
p3=ggplot(data = df, aes(x = parameter, 
                         y = Equ, 
                         color = nature, 
                         size = factor(complex_eigenvalues), 
                         shape = factor(nature))) +
  geom_point(position = position_jitter(width = 0.05, height = 0.05), alpha = 0.4) +
  scale_color_manual(values = c(stable = "blue", unstable = "red")) +
  # scale_size_manual(values = c("0" = 2, "1" = 4)) +
  labs(title = "Bifurcation Diagram",
       x = "Parameter",
       y = "Equilibrium") +
  scale_x_log10() +
  theme_minimal()
print(p3)




# Phase portrait investigation
# output results variables
para_value <- c()
equi_ife <- c()
equi_ie <- c()
nature <- c()
count <- 0
eq_points <- data.frame(Para=numeric(), P_star=numeric(), V_star=numeric(), Pv_star=numeric(), equi_I=numeric())

#5e6, 5.1e6, 0.44, 0.43, .00079, .00081
critical_point <- c(1.0e-06,     # stable ife
                    0.00079,   #stable ife near switch point
                    0.00081,     #unstable ife near switch point
                    0.0081,
                    0.4,
                    0.44,
                    50000,
                    5.0e6,
                    5.1e6)     #unstable ife

# Loop to simulate equilibria and analyze
for (j in 1:length(critical_point)){
  #parameter values
  parameters <- list(lambda1 =2.5e0, 
                     b = 0.061,
                     delta1 = 0.01,
                     delta2 = 10^4,
                     pi = 10, 
                     k = 47,
                     kc = 10^3,
                     rho1 = 0.5, 
                     gamma1 = critical_point[j], 
                     beta1 = 1.5 *10^(-6),
                     beta2 = 10^(-6),
                     rho2 = 1.8,
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
                     delta2 = 10^4,
                     pi = 10,
                     k = 47,
                     kc = 10^3,
                     rho1 = 0.5,
                     gamma1 = critical_point[j],
                     beta1 = 1.5 *10^(-4),
                     beta2 = 10^(-4),
                     rho2 = 0.8,
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

  times <- seq(0, 1500, by = 0.1)
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

#3D plots
library(plotly)

critical_point <- c(5e-4,     # stable ife
                    0.0017,   #stable ife near switch point
                    0.0019,     #unstable ife near switch point
                    500,
                    45000,
                    46000,
                    100000)     #unstable ife

for (j in 1:length(critical_point)){
  #parameter values
  parameters <- list(lambda1 =2.5e0,
                     b = 0.061,
                     delta1 = 0.01,
                     delta2 = 10^4,
                     pi = 1,
                     k = 47,
                     kc = 10^3,
                     rho1 = 0.5,
                     gamma1 = critical_point[j],
                     beta1 = 1.5 *10^(-7),
                     beta2 = 10^(-7),
                     rho2 = 1.8,
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










