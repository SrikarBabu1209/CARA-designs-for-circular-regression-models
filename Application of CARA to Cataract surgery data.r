rm(list = ls())

library(readxl)
library(dplyr)
library(circular)
library(numDeriv)

clean_sheet <- function(path, sheetname, trt_label){
  
  raw <- read_excel(path,sheet = sheetname,col_names = FALSE)
  
  #Header rows
  h1 <- as.character(unlist(raw[3, ]))
  h2 <- as.character(unlist(raw[4, ]))
  
  #Create names
  cn <- ifelse(is.na(h1) | h1 == "",h2,paste(h1, h2, sep = "_"))
  
  #Remove header rows
  dat <- raw[-c(1,2,3,4), ]
  colnames(dat) <- cn
  rownames(dat) <- NULL
  
  #Add treatment
  dat$treatment <- trt_label
  return(dat)
}
file_path <- file.choose()
SICS  <- clean_sheet(file_path, "SICS (A)", 1)

SNARE <- clean_sheet(file_path, "SNARE (B)", 2)

CONV  <- clean_sheet(file_path, "CONV (C)", 3)

TORS  <- clean_sheet(file_path, "TORS (D)", 4)

all_dat <- bind_rows(SICS, SNARE, CONV, TORS)
dim(all_dat)


cataract <- all_dat %>%
  
  transmute(
    
    patient_id = PID,
    
    treatment = factor(treatment,levels = c(1,2,3,4),labels = c("SICS","SNARE","CONV","TORS")),
    
    age = as.numeric(Age),
    
    sex = Sex,
    
    surgeon = Surgeon,
    
    eye = Eye,
    
    #Pre-operative astigmatism
    
    preop_k = as.numeric(`K...14`),
    
    preop_axis = as.numeric(`Axis...15`),
    
    #Surgery induced astigmatism (1 Month)
    
    sia_mag = as.numeric(`SIA_K...36`),
    
    sia_axis = as.numeric(`AXIS...37`))

names(cataract)
head(cataract)
summary(cataract)
table(cataract$treatment)
any(is.na(cataract))


#Cleaning data

cataract_clean <- cataract %>%
  
  filter(!is.na(age),!is.na(preop_k),!is.na(preop_axis),!is.na(sia_mag),!is.na(sia_axis),!is.na(treatment))

#Converting Treatment to Factor

cataract_clean$treatment <- as.factor(cataract_clean$treatment)

#Summary Table

summary_table <- data.frame(
  
  Variable = c("Age","Preop K","Preop Axis","SIA Magnitude","SIA Axis"),
  Mean = c(mean(cataract_clean$age),mean(cataract_clean$preop_k),mean(cataract_clean$preop_axis),mean(cataract_clean$sia_mag),mean(cataract_clean$sia_axis)),
  SD = c(sd(cataract_clean$age),sd(cataract_clean$preop_k),sd(cataract_clean$preop_axis),sd(cataract_clean$sia_mag),sd(cataract_clean$sia_axis))
  )

summary_table[, -1] <- round(summary_table[, -1],3)

cat("\n========================================\n")
cat("Summary Table\n")
cat("========================================\n")

print(summary_table)


#Double Angle Transformation

cataract2 <- cataract_clean %>%
  mutate(
    preop_axis2 =(2 * preop_axis) %% 360,
    sia_axis2 =(2 * sia_axis) %% 360)


#Circular Variables

preop_circ <- circular(cataract2$preop_axis2,units = "degrees",modulo = "2pi")
sia_circ <- circular(cataract2$sia_axis2,units = "degrees",modulo = "2pi")


#Converting to Radians
theta_x <- as.numeric(conversion.circular(preop_circ,units = "radians"))
theta_y <- as.numeric(conversion.circular(sia_circ,units = "radians"))

#Complex Representation

X <- exp(1i * theta_x)
Y <- exp(1i * theta_y)

#Standardizing Covariates

cataract2$age_std <- as.numeric(scale(cataract2$age))
cataract2$k_std <- as.numeric(scale(cataract2$preop_k))

invlogit <- function(x){1 / (1 + exp(-x))}


#Wrapped Cauchy Negative Loglikelihood

wc_negloglik <- function(par, X, Y){
  r_raw   <- par[1]
  phi     <- par[2]
  rho_raw <- par[3]
  r <- 0.75 * invlogit(r_raw)
  rho <- 0.80 * invlogit(rho_raw)
  beta <- r * exp(1i * phi)
  eta <- (X + beta) / (1 + Conj(beta) * X)
  mu <- Arg(eta)
  y <- Arg(Y)
  
  eps <- 1e-8

  ll <- sum(log(1 - rho^2 + eps)-log(1 +rho^2 -2 * rho * cos(y - mu) +eps))
  ridge <- 0.20 * r^2
  return(  -(ll - ridge))
}

#Fitting CIRCULAR-CIRCULAR MODEL

fit_wc <- optim(par = c(0,0,0),fn = wc_negloglik,X = X,Y = Y,method = "BFGS",control = list(maxit = 3000,reltol = 1e-8))

r_hat <- 0.75 * invlogit(fit_wc$par[1])

phi_hat <- fit_wc$par[2]

rho_hat <- 0.80 * invlogit(fit_wc$par[3])

beta_hat <- r_hat * exp(1i * phi_hat)


#Circular-Circular Results

cc_results <- data.frame(Parameter = c("Modulus","Argument","Rho"),Estimate = c(Mod(beta_hat),Arg(beta_hat),rho_hat))

cc_results$Estimate <- round(cc_results$Estimate,4)

cat("\n========================================\n")
cat("Circular-Circular Model\n")
cat("========================================\n")

print(cc_results)

eta_hat <- (X + beta_hat) / (1 + Conj(beta_hat) * X)
mu_hat <- Arg(eta_hat)


residuals_wc <- atan2(sin(theta_y - mu_hat),cos(theta_y - mu_hat))

res_num <- as.numeric(residuals_wc)

residual_table <- data.frame(Statistic = c("Mean","SD","Circular Concentration"),Value = c(mean(res_num),sd(res_num),mean(cos(res_num))))

residual_table$Value <- round(residual_table$Value,4)

cat("\n========================================\n")
cat("Residual Diagnostics\n")
cat("========================================\n")

print(residual_table)

#Linear-Circular Negative Loglikelihood

lc_negloglik <- function(par,x1,x2,y){
  mu0     <- par[1]
  b1      <- par[2]
  b2      <- par[3]
  rho_raw <- par[4]
  rho <- 0.75 * invlogit(rho_raw)
  eta <- b1 * x1 +b2 * x2
  mu <- mu0 + 2 * atan(eta)
  eps <- 1e-8
  ll <- sum(log(1 - rho^2 + eps)-log(1 +rho^2 -2 * rho * cos(y - mu) +eps))
  ridge <- 0.20 * (b1^2 + b2^2)
 return(-(ll - ridge))
}

#Fitting Global Model

fit_lc <- optim(par = c(0,0,0,0),
  fn = lc_negloglik,
  x1 = cataract2$age_std,
  x2 = cataract2$k_std,
  y = theta_y,
  method = "BFGS",
  control = list(maxit = 3000,reltol = 1e-8))

mu0_hat <- fit_lc$par[1]
b1_hat <- fit_lc$par[2]
b2_hat <- fit_lc$par[3]
rho_lc <- 0.75 * invlogit(fit_lc$par[4])

lc_results <- data.frame(Parameter = c("mu0","b1","b2","rho"),Estimate = c(mu0_hat,b1_hat,b2_hat,rho_lc))

lc_results$Estimate <- round(lc_results$Estimate,4)

cat("\n========================================\n")
cat("Global Linear-Circular Model\n")
cat("========================================\n")

print(lc_results)


treatment_results <- data.frame()
xi_vec <- numeric()

for(trt in levels(cataract2$treatment)){
  dat_trt <- cataract2[cataract2$treatment == trt,]
  if(nrow(dat_trt) < 10){next}
  
  y_trt <- as.numeric(conversion.circular(circular(dat_trt$sia_axis2,units = "degrees",modulo = "2pi"),units = "radians"))
  
  fit_trt <- optim(par = c(0,0,0,0),
    fn = lc_negloglik,
    x1 = dat_trt$age_std,
    x2 = dat_trt$k_std,
    y = y_trt,
    method = "BFGS",
    control = list(maxit = 2000))
  
  mu0_t <- fit_trt$par[1]
  b1_t <- fit_trt$par[2]
  b2_t <- fit_trt$par[3]
  rho_t <- 0.75 * invlogit(fit_trt$par[4])

    temp_table <- data.frame(Treatment = trt,mu0 = mu0_t,b1 = b1_t,b2 = b2_t,rho = rho_t)
  
  treatment_results <- rbind(treatment_results,temp_table)
  
  H <- hessian(func = lc_negloglik,x = fit_trt$par,x1 = dat_trt$age_std,x2 = dat_trt$k_std,y = y_trt)
  
  H <- H + diag(0.05,nrow(H))
  eigvals <- Re(eigen(H)$values)
  
  if(any(eigvals <= 0)){
    xi_val <- 1} else {
    vcov_mat <- solve(H)
    xi_val <- mean(diag(vcov_mat))}
  
  #Stabilization
  xi_val <- median(c(0.25, xi_val, 5))
  xi_vec[trt] <- xi_val
}

treatment_results[, -1] <- round(treatment_results[, -1],4)

cat("\n========================================\n")
cat("Treatmentwise Linear-Circular Fits\n")
cat("========================================\n")

print(treatment_results)

xi_df <- data.frame(Treatment = names(xi_vec),xi_k = as.numeric(xi_vec))
xi_df$xi_k <- round(xi_df$xi_k,4)

cat("\n========================================\n")
cat("Asymptotic Precision\n")
cat("========================================\n")

print(xi_df)

#Failure Probabilities
delta <- pi / 2
failure_tab <- cataract2 %>%
  mutate(y_rad = theta_y,
    failure =abs(y_rad) > delta) %>%
  group_by(treatment) %>%
  summarise(qk = mean(failure))

failure_tab <- as.data.frame(failure_tab)


colnames(failure_tab)[1] <- "Treatment"

failure_tab$qk <- pmax(failure_tab$qk,0.05)

failure_tab$qk <- round(failure_tab$qk,4)

cat("\n========================================\n")
cat("Failure Probabilities\n")
cat("========================================\n")

print(failure_tab)

alloc_tab <- merge(failure_tab,xi_df,by = "Treatment")

#Stabilized CARA Allocation

alloc_tab$weight <- sqrt(alloc_tab$xi_k / alloc_tab$qk)

#Smooth Shrinkage

alloc_tab$weight <- sqrt(alloc_tab$weight)

alloc_tab$pi_star <- alloc_tab$weight / sum(alloc_tab$weight)

alloc_tab$qk <- round(alloc_tab$qk,4)

alloc_tab$xi_k <- round(alloc_tab$xi_k,4)

alloc_tab$weight <- round(alloc_tab$weight,4)

alloc_tab$pi_star <- round(alloc_tab$pi_star,4)

cat("\n========================================\n")
cat("Stabilized CARA Allocation\n")
cat("========================================\n")

print(alloc_tab[order(-alloc_tab$pi_star),])

final_summary <- data.frame(Quantity = c("CC Modulus","CC Rho","LC Rho","Residual Concentration"),
  Value = c(Mod(beta_hat),rho_hat,rho_lc,mean(cos(res_num))))

final_summary$Value <- round(final_summary$Value,4)

cat("\n========================================\n")
cat("Final Summary\n")
cat("========================================\n")

print(final_summary)

cat("\n========================================\n")
cat("ANALYSIS COMPLETE\n")
cat("========================================\n")
  

