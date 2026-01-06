# Bifurcation for pi
# Necessary libraries
library(polynom)
library(ggplot2)
library(dplyr)

#generating parameter set
all1= seq(0,15, length.out = 1001) #pi



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
                     delta2 = 10^3,
                     pi = all1[j], 
                     k = 47,
                     kc = 10^3,
                     rho1 = 0.8, 
                     gamma1 = 8e-02, # gamma = 2e-03 te sob stable r gamma = 2e-01 a unstable stable milano
                     beta1 = 1.5 *10^(-4),
                     beta2 = 10^(-4),
                     rho2 = 2.2,
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


# Data frame to print
df2 <- data.frame(
  parameter = para_value,
  IFE = equi_ife,
  IE = equi_ie,
  P = equi_P,
  V = equi_V,
  Pv = equi_Pv,
  nature = nature
)

# Install if not already installed
install.packages("xtable")

# Load package
library(xtable)

# Convert to LaTeX
print(xtable(df2[455:475,]), type = "latex")



# dynamics of the model for different gamma, lambda1 value (series of phase potraits try)
library(polynom)
library(deSolve)
library(ggplot2)
library(dplyr)
library(egg)

# Pi values
all1= seq(0,15, length.out = 101) #pi

count <- 0



# Define the directory
plotout_dir <- "C:/Users/zanna/Downloads/Final Project - MATH 8996/pi_phase_plots"
# dataout_dir <- "C:/Users/zanna/Downloads/Final Project - MATH 8996" 
# 
# Create the directory if it doesn't exist
if (!dir.exists(plotout_dir)) {
dir.create(plotout_dir)
}



# for (i in 1:length(all1)) {
for (j in 1:length(all1)){
  count <- count+1
  
  #parameter values
  parameters <- list(lambda1 =2.5e0, 
                     b = 0.061,
                     delta1 = 0.01,
                     delta2 = 10^3,
                     pi = all1[j], 
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
  state <- c(P = 10,
             V = 10^-4,
             Pv = 0,
             I = 0)  
  
  #The main model
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
  
  
  times <- seq(0, 100, by = 0.01)
  out <- ode(y = state, times = times, func = Model, parms = parameters, atol = 1e-6, rtol = 1e-6, hmin = 1e-6, method = "lsoda")
  
  
  I_star = out[length(out[,1]),5]
  # if (I_star < 1e-10){I_star =0}
  # Pv_star <- with(parameters, { (delta2 * I_star) / (sigma1 * (kc - I_star)) })
  # P_star <- with(parameters, {  (-delta2 * I_star / (rho1 * sigma1 * (kc - I_star))) *
  #     ((1 - pi) * sigma1 * (kc - I_star) + rhoM) })
  # numerator <- with(parameters, { -delta2 * I_star * (1 + (rhoM / (sigma1 * (kc - I_star)))) +
  #     (beta2 + beta1 * kc) * P_star })
  # denominator <- with(parameters, { (gamma1 * kc - rho2 - gamma1 * I_star) })
  # V_star <- numerator / denominator
  
  
  # Set the file name
  filename <- file.path(plotout_dir, paste0("plot_", count, ".png"))
  # data_output_file <- file.path(dataout_dir, "output.csv")
  
  
  # Open a graphics device
  png(filename)
  
  #plotting 
  par(mfrow = c(2, 3), mar = c(4.5, 4.5, 3, 2), oma = c(0, 0, 3, 0))
  
  plot(out[,1],out[,2],type = 'l', xlab = "time", ylab = "Pathogen")
  plot(out[,1],out[,3],type = 'l', xlab = "time", ylab = "Vimentin" )
  plot(out[,1],out[,4],type = 'l', xlab = "time", ylab = "P-V complex" )
  plot(out[,1],out[,5],type = 'l', xlab = "time", ylab = "Infected cell" )
  plot(out[,3],out[,5],type = 'l', xlab = "Vimentin", ylab = "Infected cell" )
  # lines(V_star,I_star,type = 'p',col="red")
  plot(out[,2],out[,5],type = 'l', xlab = "Pathogen", ylab = "Infected cell" )
  # lines(P_star,Pv_star,type = 'p',col="red")
  #Title of plots
  #mtext(paste("lambda=", all1[i], ", gamma =", all2[j]), side = 3, line = -2, outer = TRUE, cex = 1.1 )
  mtext(sprintf("pi= %f", parameters$pi, all1[j]), side = 3, line = -2, outer = TRUE, cex = 1.1 )
  # mtext(sprintf("eigenvalues= %s, %s, %s, %s",  values[1],values[2],values[3],values[4]), side = 1, line = -2, outer = TRUE, cex = 1.1 )
  # complex1 <- paste0(Re(values[1]), " + ", Im(values[1]), "i")
  # complex2 <- paste0(Re(values[2]), " + ", Im(values[2]), "i")
  # complex3 <- paste0(Re(values[3]), " + ", Im(values[3]), "i")
  # complex4 <- paste0(Re(values[4]), " + ", Im(values[4]), "i")
  # mtext(sprintf("eigenvalues= %s, %s, %s, %s", complex1,complex2,complex3,complex4), side = 1, line = -2, outer = TRUE, cex = 1.1)
  
  # Close the graphics device
  dev.off()
  
  # #Adding lambda and gamma value in the data list
  # gamma_out <- rep(all2[j],nrow(out))
  # lambda_out <- rep(all1[i], nrow(out))
  # modified_out <- cbind(out, gamma_out, lambda_out)
  # global_out <<- rbind(global_out,modified_out)
  # write.csv(global_out, data_output_file)
  
}
# }


cat("Plots saved in:", plotout_dir)
# cat("CSV file saved in:", data_output_file) 
