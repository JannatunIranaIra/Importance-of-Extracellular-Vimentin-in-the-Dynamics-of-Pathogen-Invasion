# Computing equilibrium (both zero + non-zero) 
# and analyzing the equilibria by computing its nature, resilience, reactivity
# counting complex eigenvalues and where the real part of the eigenvalues are zero
# two separate bifurcation diagram
# Phase planes 2-D + 3-D
# eigenvalue vs parameter plot to check Hopf bifurcation


# Necessary libraries
library(polynom)
library(ggplot2)
library(dplyr)

#generating parameter set
bodyframe <-seq(from=0,to=1,length.out= 1001)
m8 = 1e-7*bodyframe
m7 = 1e-6*bodyframe
m6 = 1e-5*bodyframe
m5 = 1e-4*bodyframe
m4 = 1e-3*bodyframe
m3 = 1e-2*bodyframe
m2 = 1e-1*bodyframe
m1 = bodyframe
p0 = 10*bodyframe
p1 = 10^2*bodyframe
p2 = 10^3*bodyframe
p3 = 10^4*bodyframe
p4 = 10^5*bodyframe
p5 = 10^6*bodyframe
p6 = 10^7*bodyframe
all1 = unique(c(m4,m3,m2,m1,p0,p1,p2,p3))
all1=all1[2:length(all1)] #gamma1


# # Function for solving nonegative real I^* (doesn't work)
# find_Equilibruim <- function(A_coeff,B_coeff,C_coeff,D_coeff){
#   solution = polyroot(c(A_coeff,B_coeff,C_coeff,D_coeff))
#   
#   for (i2 in 1:length(solution)){
#     if (Im(solution[i2]) == 0 && Re(solution[i2])>= 0) {
#       return(solution[i2])
#     }
#     
#   }
# }

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
    
    # Debugging - Check values
    # cat("Equilibrium", j1, ":", "P_star =", P_star, "V_star =", V_star, "Pv_star =", Pv_star, "I_star =", equi_I, "\n")
    if (any(is.nan(P_star) | is.nan(V_star) | is.nan(Pv_star) | is.nan(equi_I))) {
      # stop("NaN detected in equilibrium calculations!")
      next
    } else {
      
      if (equi_I>= 0 & Pv_star >= 0 & P_star>= 0 & V_star>=0) {
        count <- count +1
        para_value[count] <- all1[j]
        
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
            if(Re(unstable_evalues[i5]) == 0 ){
              temp_r0 <- temp_r0 +1
              critical_point <- para_value[count]
              critical_equi <- equi_I
              print(unstable_evalues[i5])
              print(para_value[count])
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
                          size = factor(nature),
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
print(p2)

# All together bifurcation diagram
p3=ggplot(data = df, aes(x = parameter, 
                         y = Equ, 
                         color = nature, 
                         size = factor(nature), 
                         shape = factor(complex_eigenvalues))) +
  geom_point(position = position_jitter(width = 0.05, height = 0.05), alpha = 0.4) +
  scale_color_manual(values = c(stable = "blue", unstable = "red")) +
  # scale_size_manual(values = c("0" = 2, "1" = 4)) +
  ylim(0,0.051) +
  labs(title = "Bifurcation Diagram",
       x = "Parameter",
       y = "Equilibrium") +
  scale_x_log10() +
  theme_minimal()
print(p3)

# Phase planes 2-D + 3-D investigation
# library(deSolve)
# # critical_point <- c(critical_point, 8.0e-08, 9.0e-08, 1e+02, 19000)
# for (j in 1:length(critical_point)){
#   #parameter values
#   parameters <- list(lambda1 =2.5e0, 
#                      b = 0.061,
#                      delta1 = 0.01,
#                      delta2 = 10^4,
#                      pi = 20, 
#                      k = 47,
#                      kc = 10^3,
#                      rho1 = 0.5, 
#                      gamma1 = 8.0e-08, 
#                      beta1 = 1.5 *10^(-6),
#                      beta2 = 10^(-6),
#                      rho2 = 2.2,
#                      sigma1 = 10^(-5), 
#                      rhoM = 0.0024)
#   
#   #Solving ODEs
#   state <- c(P = 10,
#              V = 10^-4,
#              Pv = 0,
#              I = 0)
#   
#   Model<-function(t, state, parameters) {
#     with(as.list(c(state, parameters)),{
#       # rate of change
#       dP <- - parameters$lambda1*P*V + (parameters$delta1 + parameters$delta2 + parameters$k* parameters$kc)*parameters$pi*I- parameters$rho1*P
#       dV <- - parameters$lambda1*P*V + parameters$gamma1*(parameters$kc-I)*V+ parameters$beta2*P+ parameters$beta1*P*parameters$kc - parameters$rho2*V
#       dPv <- parameters$lambda1*P*V - parameters$sigma1*Pv*(parameters$kc-I)- parameters$rhoM*Pv
#       dI <- parameters$sigma1*Pv*(parameters$kc-I) +(- parameters$delta2+ parameters$b- parameters$delta1- parameters$k* parameters$kc)*I
#       
#       # return the rate of change
#       list(c(dP, dV, dPv, dI))
#     }) # end with(as.list ...
#   }
#   
#   times <- seq(0, 100, by = 0.1)
#   out <- ode(y = state, times = times, func = Model, parms = parameters, atol = 1e-14, rtol = 1e-14)
#   
#   # Some expressions used for Equilibrium values
#   L <- (parameters$beta1 * parameters$kc + parameters$beta2) / parameters$rho1
#   TT <- parameters$k * parameters$kc + parameters$delta2 - (parameters$b - parameters$delta1)
#   Tb <- parameters$k * parameters$kc + parameters$delta2 + parameters$delta1
#   
#   # Coefficients
#   A_coeff <- ( 
#     ( (2*L+1)*parameters$pi*Tb + (parameters$gamma1* parameters$rho1 )/parameters$lambda1  )*TT - parameters$pi^2*L*Tb^2 - (L+1)*TT^2 )* parameters$lambda1*parameters$sigma1^2
#   
#   B_coeff <- (
#     ( -(2*L+1)*parameters$pi*Tb + (parameters$gamma1* parameters$rho1 )/parameters$lambda1  )* (2*parameters$sigma1*parameters$kc +parameters$rhoM )*parameters$lambda1 *TT*parameters$sigma1
#     -(parameters$gamma1*parameters$kc - parameters$rho2)*parameters$rho1*TT*parameters$sigma1^2 
#     +2*parameters$pi^2*L*parameters$kc*Tb^2*parameters$sigma1^2*parameters$lambda1 + 2*(L+1)*(parameters$sigma1* parameters$kc +parameters$rhoM)*TT^2*parameters$lambda1 *parameters$sigma1)
#   
#   C_coeff <- ( ( (2*L +1)*parameters$pi*Tb +(parameters$gamma1* parameters$rho1) /(parameters$lambda1) )* (parameters$sigma1*parameters$kc +parameters$rhoM)*parameters$sigma1 *parameters$kc* TT *parameters$lambda1 
#                + (parameters$gamma1 *parameters$kc - parameters$rho2) * (2* parameters$sigma1 *parameters$kc +parameters$rhoM) * parameters$sigma1* parameters$rho1 *TT - parameters$pi^2 *parameters$kc^2 *parameters$sigma1^2 *parameters$lambda1 *L*Tb^2 
#                - (L+1) *(parameters$sigma1 *parameters$kc +parameters$rhoM)^2 *parameters$lambda1 *TT^2)
#   
#   D_coeff <- -(
#     (parameters$gamma1 * parameters$kc - parameters$rho2) * parameters$rho1*TT * (parameters$sigma1 * parameters$kc + parameters$rhoM )* parameters$sigma1* parameters$kc
#   )
#   
#   # Find equilibrium points
#   I_star <- find_Equilibruim(A_coeff, B_coeff, C_coeff, D_coeff, tolerance = 1e-5)
#   equi_I <- c()
#   for (j1 in 1:length(I_star)) {
#     equi_I <- I_star[j1]
#   # Calculate other variables at I*
#   Pv_star = (-(parameters$b-parameters$delta1- parameters$delta2 - parameters$k*parameters$kc)*equi_I[j1])/(parameters$sigma1*(parameters$kc-equi_I[j1]))
#   P_star = equi_I[j1]/(parameters$rho1*parameters$sigma1*(parameters$kc-equi_I[j1]))*((parameters$sigma1*(parameters$kc-equi_I[j1]) + parameters$rhoM)*parameters$b +(parameters$delta1+ parameters$delta2 + parameters$k*parameters$kc)*(parameters$sigma1*(parameters$kc-equi_I[j1])*(parameters$pi-1)-parameters$rhoM))
#   V_star = ((parameters$beta2+parameters$beta1*parameters$kc)*equi_I[j1])/(parameters$rho1*parameters$sigma1*(parameters$kc-equi_I[j1])*(parameters$lambda1*P_star + parameters$rho2 - parameters$gamma1 *(parameters$kc-equi_I[j1])))*((parameters$sigma1*(parameters$kc-equi_I[j1]) + parameters$rhoM)*parameters$b +(parameters$delta1+parameters$delta2 + parameters$k*parameters$kc)*(parameters$sigma1*(parameters$kc-equi_I[j1])*(parameters$pi-1)  - parameters$rhoM))
#   }
#     
#     par(mfrow=c(2,3))
#     plot(out[,2],out[,3],type = 'l', xlab = "Pathogen", ylab = "Vimentin" )
#     for (j1 in 1:length(I_star)) {
#       # Adding the equilibrium points (x=critical_point, y= critical_equi)
#     points(critical_point[j],I_star[j1], col='red', pch=16)
#     }
#     plot(out[,2],out[,4],type = 'l', xlab = "Pathogen", ylab = "P-V complex" )
#     for (j1 in 1:length(I_star)) {
#       # Adding the equilibrium points (x=critical_point, y= critical_equi)
#       points(critical_point[j],I_star[j1], col='red', pch=16)
#     }
#     plot(out[,2],out[,5],type = 'l', xlab = "Pathogen", ylab = "Infected cell" )
#     for (j1 in 1:length(I_star)) {
#       # Adding the equilibrium points (x=critical_point, y= critical_equi)
#       points(critical_point[j],I_star[j1], col='red', pch=16)
#     }
#     plot(out[,3],out[,4],type = 'l', xlab = "Vimentin", ylab = "P-V complex" )
#     # points(critical_point[j],critical_equi, col='red', pch=16)
#     plot(out[,3],out[,5],type = 'l', xlab = "Vimentin", ylab = "Infected cell" )
#     # points(critical_point[j],critical_equi, col='red', pch=16)
#     plot(out[,4],out[,5],type = 'l', xlab = "P-V complex", ylab = "Infected cell" )
#     # points(critical_point[j],critical_equi, col='red', pch=16)
# }
# 
# 
#   
#   
# # }


# eigenvalue vs parameter plot to check Hopf bifurcation
plot(Re(df1$eigenvalue1),Im(df1$eigenvalue1))
plot(Re(df1$eigenvalue2),Im(df1$eigenvalue2))
plot(Re(df1$eigenvalue3),Im(df1$eigenvalue3))
plot(Re(df1$eigenvalue4),Im(df1$eigenvalue4))
# persp(df1$parameter, Re(df1$eigenvalue1),Im(df1$eigenvalue1))
# install.packages("rgl")
library(rgl)
plot3d(x = Re(df1$eigenvalue1), y = Im(df1$eigenvalue1), z = df1$parameter,
       size = 5,
       xlab = "Real Part", ylab = "Imaginary Part", zlab = "Parameter")

# Open a new 3D window
open3d()

# Plot 3D Points inside a Box
plot3d(x = Re(df1$eigenvalue1), y = Im(df1$eigenvalue1), z = df1$parameter,
       size = 1.5, type = "s",
       xlab = "Real Part", ylab = "Imaginary Part", zlab = "Parameter")

# Add Grid Box for Better Visualization
grid3d(c("x", "y", "z"))

library(plotly)
plot_ly(x = Re(df1$eigenvalue1), y = Im(df1$eigenvalue1), z = df1$parameter,
        # type = "scatter3d", mode = "markers",
        marker = list(size = 5, color = df1$parameter, colorscale = "Viridis")) %>%
  layout(scene = list(
    xaxis = list(title = "Real Part"),
    yaxis = list(title = "Imaginary Part"),
    zaxis = list(title = "Parameter"),
    title = "3D Eigenvalue Plot"
  ))

plot_ly(x = Re(df1$eigenvalue2), y = Im(df1$eigenvalue2), z = df1$parameter,
        # type = "scatter3d", mode = "markers",
        marker = list(size = 5, color = df1$parameter, colorscale = "Viridis")) %>%
  layout(scene = list(
    xaxis = list(title = "Real Part"),
    yaxis = list(title = "Imaginary Part"),
    zaxis = list(title = "Parameter"),
    title = "3D Eigenvalue Plot"
  ))

plot_ly(x = Re(df1$eigenvalue3), y = Im(df1$eigenvalue3), z = df1$parameter,
        # type = "scatter3d", mode = "markers",
        marker = list(size = 5, color = df1$parameter, colorscale = "Viridis")) %>%
  layout(scene = list(
    xaxis = list(title = "Real Part"),
    yaxis = list(title = "Imaginary Part"),
    zaxis = list(title = "Parameter"),
    title = "3D Eigenvalue Plot"
  ))

plot_ly(x = Re(df1$eigenvalue4), y = Im(df1$eigenvalue4), z = df1$parameter,
        # type = "scatter3d", mode = "markers",
        marker = list(size = 5, color = df1$parameter, colorscale = "Viridis")) %>%
  layout(scene = list(
    xaxis = list(title = "Real Part"),
    yaxis = list(title = "Imaginary Part"),
    zaxis = list(title = "Parameter"),
    title = "3D Eigenvalue Plot"
  ))

