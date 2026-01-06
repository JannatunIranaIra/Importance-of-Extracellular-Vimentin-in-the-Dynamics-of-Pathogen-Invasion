# Phase planes 2-D + 3-D investigation
library(polynom)
library(deSolve)
library(ggplot2)
library(dplyr)

#generating parameter set
bodyframe <-seq(from=0,to=1,length.out= 101)
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
all1 = unique(c(m6,m5,m4,m3,m2,m1,p0,p1,p2))
all1=all1[2:length(all1)] #gamma1

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


critical_point <- c(0.0001,0.02)
# 8.0e-08, 9.0e-08, 1e+02, 19000)
for (j in 1:length(critical_point)){
  #parameter values
  parameters <- list(lambda1 =2.5e0, 
                     b = 0.061,
                     delta1 = 0.01,
                     delta2 = 10^4,
                     pi = 5, 
                     k = 47,
                     kc = 10^3,
                     rho1 = 0.5, 
                     gamma1 = all1[j], 
                     beta1 = 1.5 *10^(-6),
                     beta2 = 10^(-6),
                     rho2 = 1.8,
                     sigma1 = 10^(-4), 
                     rhoM = 0.0024)
  
  #Solving ODEs
   state <- c(P = 10,
              V = 10^-4,
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

   times <- seq(0, 100, by = 0.01)
   out <- ode(y = state, times = times, func = Model, parms = parameters, atol = 1e-8, rtol = 1e-12)

   # par(mfrow=c(2,3))
   # plot(out[,2],out[,3],type = 'l', xlab = "Pathogen", ylab = "Vimentin" )
   # plot(out[,2],out[,4],type = 'l', xlab = "Pathogen", ylab = "P-V complex" )
   # plot(out[,2],out[,5],type = 'l', xlab = "Pathogen", ylab = "Infected cell" )
   # plot(out[,3],out[,4],type = 'l', xlab = "Vimentin", ylab = "P-V complex" )
   # plot(out[,3],out[,5],type = 'l', xlab = "Vimentin", ylab = "Infected cell" )
   # plot(out[,4],out[,5],type = 'l', xlab = "P-V complex", ylab = "Infected cell" )

  # Some expressions used for Equilibrium values
    L <- (parameters$beta1 * parameters$kc + parameters$beta2) / parameters$rho1
    TT <- parameters$k * parameters$kc + parameters$delta2 - (parameters$b - parameters$delta1)
    Tb <- parameters$k * parameters$kc + parameters$delta2 + parameters$delta1

    # Coefficients
    A_coeff <- (
      ( (2*L+1)*parameters$pi*Tb + (parameters$gamma1* parameters$rho1 )/parameters$lambda1  )*TT - parameters$pi^2*L*Tb^2 - (L+1)*TT^2 )* parameters$lambda1*parameters$sigma1^2

    B_coeff <- (
      ( -(2*L+1)*parameters$pi*Tb + (parameters$gamma1* parameters$rho1 )/parameters$lambda1  )* (2*parameters$sigma1*parameters$kc +parameters$rhoM )*parameters$lambda1 *TT*parameters$sigma1
      -(parameters$gamma1*parameters$kc - parameters$rho2)*parameters$rho1*TT*parameters$sigma1^2
      +2*parameters$pi^2*L*parameters$kc*Tb^2*parameters$sigma1^2*parameters$lambda1 + 2*(L+1)*(parameters$sigma1* parameters$kc +parameters$rhoM)*TT^2*parameters$lambda1 *parameters$sigma1)

    C_coeff <- ( ( (2*L +1)*parameters$pi*Tb +(parameters$gamma1* parameters$rho1) /(parameters$lambda1) )* (parameters$sigma1*parameters$kc +parameters$rhoM)*parameters$sigma1 *parameters$kc* TT *parameters$lambda1
                 + (parameters$gamma1 *parameters$kc - parameters$rho2) * (2* parameters$sigma1 *parameters$kc +parameters$rhoM) * parameters$sigma1* parameters$rho1 *TT - parameters$pi^2 *parameters$kc^2 *parameters$sigma1^2 *parameters$lambda1 *L*Tb^2
                 - (L+1) *(parameters$sigma1 *parameters$kc +parameters$rhoM)^2 *parameters$lambda1 *TT^2)

    D_coeff <- -(
      (parameters$gamma1 * parameters$kc - parameters$rho2) * parameters$rho1*TT * (parameters$sigma1 * parameters$kc + parameters$rhoM )* parameters$sigma1* parameters$kc
    )
  # Find equilibrium points
    I_star <- find_Equilibruim(A_coeff, B_coeff, C_coeff, D_coeff, tolerance = 1e-5)
    I_star <- c(0,I_star)
    equi_I <- c()

    
  # Check equilibrium values before plotting
  for (j1 in 1:length(I_star)) {
    equi_I <- Re(I_star[j1])

    # Compute equilibrium points
    Pv_star <- (-(parameters$b - parameters$delta1 - parameters$delta2 - parameters$k * parameters$kc) * equi_I) /
      (parameters$sigma1 * (parameters$kc - equi_I))
    P_star <- equi_I / (parameters$rho1 * parameters$sigma1 * (parameters$kc - equi_I)) *
      ((parameters$sigma1 * (parameters$kc - equi_I) + parameters$rhoM) * parameters$b +
         (parameters$delta1 + parameters$delta2 + parameters$k * parameters$kc) *
         (parameters$sigma1 * (parameters$kc - equi_I) * (parameters$pi - 1) - parameters$rhoM))
    V_star <- ((parameters$beta2 + parameters$beta1 * parameters$kc) * equi_I) /
      (parameters$rho1 * parameters$sigma1 * (parameters$kc - equi_I) *
         (parameters$lambda1 * P_star + parameters$rho2 - parameters$gamma1 * (parameters$kc - equi_I))) *
      ((parameters$sigma1 * (parameters$kc - equi_I) + parameters$rhoM) * parameters$b +
         (parameters$delta1 + parameters$delta2 + parameters$k * parameters$kc) *
         (parameters$sigma1 * (parameters$kc - equi_I) * (parameters$pi - 1) - parameters$rhoM))

    # Debugging - Check values
    # cat("Equilibrium", j1, ":", "P_star =", P_star, "V_star =", V_star, "Pv_star =", Pv_star, "I_star =", equi_I, "\n")
    if (any(is.nan(P_star) | is.nan(V_star) | is.nan(Pv_star) | is.nan(equi_I))) {
      # stop("NaN detected in equilibrium calculations!")
      next
    } else {
      

    if (equi_I>= 0 & Pv_star >= 0 & P_star>= 0 & V_star>=0) {
      count <- count +1
      # para_value[count] <- all1[j]
      Pv_star [count] <- (-(parameters$b - parameters$delta1 - parameters$delta2 - parameters$k * parameters$kc) * equi_I) /
        (parameters$sigma1 * (parameters$kc - equi_I))
 


      }
    }
    par(mfrow = c(2, 3))  # Keep multi-panel layout
    plot(out[, 2], out[, 3], type = 'l', xlab = "Pathogen", ylab = "Vimentin")
    plot(out[, 2], out[, 4], type = 'l', xlab = "Pathogen", ylab = "P-V complex")
    plot(out[, 2], out[, 5], type = 'l', xlab = "Pathogen", ylab = "Infected cell")
    plot(out[, 3], out[, 4], type = 'l', xlab = "Vimentin", ylab = "P-V complex")
    plot(out[, 3], out[, 5], type = 'l', xlab = "Vimentin", ylab = "Infected cell")
    plot(out[, 4], out[, 5], type = 'l', xlab = "P-V complex", ylab = "Infected cell")
  }
}

# Data frame
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



  # Redraw the phase plots with equilibrium points
  par(mfrow = c(2, 3))  # Keep multi-panel layout
  plot(log10(out[, 2]), log10(out[, 3]), type = 'l', xlim = c(log10(min(out[,2],P_star)), log10(max(out[,2],P_star))), ylim = c(log10(min(out[,3],V_star)), log10(max(out[,3],V_star))) ,xlab = "Pathogen", ylab = "Vimentin")
  for (j1 in 1:length(I_star)) {
    points(log10(P_star[j1]), log10(V_star[j1]), col = "red", pch = 19, cex = 1.5)
  }

  plot(log10(out[, 2]), log10(out[, 4]), type = 'l', xlim = c(log10(min(out[,2],P_star)), log10(max(out[,2],P_star))), ylim = c(log10(min(out[,4],Pv_star)), log10(max(out[,4],Pv_star))) ,xlab = "Pathogen", ylab = "P-V complex")
  for (j1 in 1:length(I_star)) {
    points(log10(P_star[j1]), log10(Pv_star[j1]), col = "red", pch = 19, cex = 1.5)
  }

  plot(log10(out[, 2]), log10(out[, 5]), type = 'l', xlim = c(log10(min(out[,2],P_star)), log10(max(out[,2],P_star))), ylim = c(log10(min(out[,5],equi_I)), log10(max(out[,5],equi_I))) ,xlab = "Pathogen", ylab = "Infected")
  for (j1 in 1:length(I_star)) {
    points(log10(P_star[j1]), log10(equi_I), col = "red", pch = 19, cex = 1.5)
  }

  plot(log10(out[, 3]), log10(out[, 4]), type = 'l', xlim = c(log10(min(out[,3],V_star)), log10(max(out[,3],V_star))), ylim = c(log10(min(out[,4],Pv_star)), log10(max(out[,4],Pv_star))) ,xlab = "Vimentin", ylab = "P-V complex")
  for (j1 in 1:length(I_star)) {
    points(log10(V_star[j1]), log10(Pv_star[j1]), col = "red", pch = 19, cex = 1.5)
  }

  plot(log10(out[, 3]), log10(out[, 5]), type = 'l', xlim = c(log10(min(out[,3],V_star)), log10(max(out[,3],V_star))), ylim = c(log10(min(out[,5],equi_I)), log10(max(out[,5],equi_I))) ,xlab = "Vimentin", ylab = "Infected")
  for (j1 in 1:length(I_star)) {
    points(log10(V_star[j1]), log10(equi_I), col = "red", pch = 19, cex = 1.5)
  }

  plot(log10(out[, 4]), log10(out[, 5]), type = 'l', xlim = c(log10(min(out[,4],Pv_star)), log10(max(out[,4],Pv_star))), ylim = c(log10(min(out[,5],equi_I)), log10(max(out[,5],equi_I))) ,xlab = "P-V complex", ylab = "Infected")
  for (j1 in 1:length(I_star)) {
    points(log10(Pv_star[j1]), log10(equi-I), col = "red", pch = 19, cex = 1.5)
  }





  