install.packages(c("readxl", "dplyr", "circular", "CircStats"))
library(readxl)
library(dplyr)
library(circular)
library(CircStats)


#Importing and Cleaning Raw Excel Sheets

clean_sheet <- function(path, sheetname, trt_label){
  raw <- read_excel(path, sheet = sheetname, col_names = FALSE)
  
  #Extracting header rows
  h1 <- as.character(unlist(raw[3, ]))
  h2 <- as.character(unlist(raw[4, ]))
  
  cn <- ifelse(is.na(h1) | h1 == "", h2, paste(h1, h2, sep = "_"))
  
  #Removing raw header rows
  dat <- raw[-c(1, 2, 3, 4), ]
  colnames(dat) <- cn
  rownames(dat) <- NULL
  
  #Adding treatment label
  dat$treatment <- trt_label
  return(dat)
}


cat("Select Excel data file in the pop-up window...\n")
file_path <- file.choose()

SICS  <- clean_sheet(file_path, "SICS (A)", 1)
SNARE <- clean_sheet(file_path, "SNARE (B)", 2)
CONV  <- clean_sheet(file_path, "CONV (C)", 3)
TORS  <- clean_sheet(file_path, "TORS (D)", 4)

#Binding all groups
all_dat <- bind_rows(SICS, SNARE, CONV, TORS)


#Transmute and Clean Primary Variables
cataract <- all_dat %>%
  transmute(
    patient_id = PID,
    treatment  = factor(treatment, levels = c(1, 2, 3, 4), labels = c("SICS", "SNARE", "CONV", "TORS")),
    age        = as.numeric(Age),
    sex        = Sex,
    surgeon    = Surgeon,
    eye        = Eye,
    #Pre-operative astigmatism
    preop_k    = as.numeric(`K...14`),
    preop_axis = as.numeric(`Axis...15`),
    #Surgery induced astigmatism (1 Month)
    sia_mag    = as.numeric(`SIA_K...36`),
    sia_axis   = as.numeric(`AXIS...37`)
  )

#Removing cases with missing values in key circular orientations
cataract_clean <- cataract %>%
  filter(
    !is.na(preop_axis),
    !is.na(sia_axis),
    !is.na(sia_mag)
  ) %>%
  mutate(
    #Converting axial angles [0, 180) to doubled angles [0, 360)
    preop_axis2 = (2 * preop_axis) %% 360,
    sia_axis2   = (2 * sia_axis) %% 360,
    
    #Defining circular class objects
    preop_circ  = circular(preop_axis2, units = "degrees", modulo = "2pi"),
    sia_circ    = circular(sia_axis2, units = "degrees", modulo = "2pi"),
    
    x_rad       = as.numeric(preop_circ) * pi / 180,
    y_rad       = as.numeric(sia_circ) * pi / 180,
    
    #Standardizing linear covariate (Age) to align with simulation specifications
    age_std     = (age - mean(age)) / sd(age)
  )


#Exploratory Circular Data Analysis (EDA)
eda_summary <- cataract_clean %>%
  group_by(treatment) %>%
  summarise(
    n             = n(),
    mean_preop    = mean(preop_circ),
    mean_sia      = mean(sia_circ),
    rho_preop     = rho.circular(preop_circ),
    rho_sia       = rho.circular(sia_circ),
    .groups       = "drop"
  )
print(eda_summary)


#Circular-Circular Möbius Regression Fitting

#Using polar coordinates to enforce constraint r = |beta| < 1
mobius_polar_loglik <- function(params, x, y) {
  r    <- params[1]
  phi  <- params[2]
  rho  <- 1 / (1 + exp(-params[3])) #Logit link function to force rho in [0, 1)
  beta <- r * exp(1i * phi)
  X_complex <- exp(1i * x)
  eta <- (X_complex + beta) / (1 + Conj(beta) * X_complex)
  mu  <- Arg(eta) %% (2 * pi)
  
  log_dens <- log(1 - rho^2) - log(2 * pi) - log(1 + rho^2 - 2 * rho * cos(y - mu))
  return(-sum(log_dens))
}

fit_mobius_group <- function(data_sub) {
  #Constrained optimization using L-BFGS-B (r bounded in [0, 0.95] to prevent boundary failure)
  lower_bounds <- c(0, 0, -5)
  upper_bounds <- c(0.95, 2 * pi, 5)
  init_val     <- c(0.2, pi, 0)
  
  opt <- optim(
    init_val, 
    mobius_polar_loglik, 
    x = data_sub$x_rad, 
    y = data_sub$y_rad, 
    method = "L-BFGS-B",
    lower = lower_bounds,
    upper = upper_bounds,
    hessian = TRUE
  )
  
  r_est    <- opt$par[1]
  phi_est  <- opt$par[2]
  beta_est <- r_est * exp(1i * phi_est)
  rho_est  <- 1 / (1 + exp(-opt$par[3]))
  
  #Calculating per-observation variance parameter (xi) using the localized Hessian
  n_k <- nrow(data_sub)
  hessian_sub <- opt$hessian[1:2, 1:2]
  
  if (all(eigen(hessian_sub)$values > 0)) {
    info_beta <- solve(hessian_sub / n_k)
    xi_est    <- mean(diag(info_beta))
  } else {
    xi_est    <- 1 - rho_est  #Falling back to circular variance if boundary optimization is unstable
  }
  
  return(list(beta = beta_est, r = r_est, phi = phi_est, rho = rho_est, xi = xi_est))
}

treatments     <- levels(cataract_clean$treatment)
mobius_results <- list()

for (trt in treatments) {
  sub_data <- filter(cataract_clean, treatment == trt)
  fit <- fit_mobius_group(sub_data)
  mobius_results[[trt]] <- fit
  cat(sprintf("Treatment: %-5s | Beta: %7.4f + %7.4fi (Mod: %.4f) | Rho: %.4f | Xi (Precision): %.4f\n", 
              trt, Re(fit$beta), Im(fit$beta), Mod(fit$beta), fit$rho, fit$xi))
}


#Linear-Circular Regression Model Fitting

#Joint log-likelihood
lin_circ_loglik_real <- function(params, x, y, trt_vec) {
  #params: beta_1, beta_2, beta_3, beta_4, mu0, logit_rho
  beta <- params[1:4]
  mu0  <- params[5]
  rho  <- 1 / (1 + exp(-params[6]))
  
  ll <- 0
  for (i in 1:length(y)) {
    mu_i <- (mu0 + 2 * atan(beta[trt_vec[i]] * x[i])) %% (2 * pi)
    ll <- ll + log(1 - rho^2) - log(1 + rho^2 - 2 * rho * cos(y[i] - mu_i))
  }
  return(-ll)
}

x_val   <- cataract_clean$age_std
y_val   <- cataract_clean$y_rad
trt_val <- as.numeric(cataract_clean$treatment) #SICS=1, SNARE=2, CONV=3, TORS=4

#Setting stable starting parameters for joint optimization
init_val <- c(rep(0, 4), mean(y_val), 0)

op <- optim(
  init_val, 
  lin_circ_loglik_real, 
  x = x_val, 
  y = y_val, 
  trt_vec = trt_val,
  method = "L-BFGS-B",
  lower = c(rep(-2, 4), 0, -5),
  upper = c(rep(2, 4), 2 * pi, 5),
  hessian = TRUE
)

#Extracting joint parameters
beta_est <- op$par[1:4]
mu0_est  <- op$par[5]
rho_est  <- 1 / (1 + exp(-op$par[6]))

for (k in 1:4) {
  cat(sprintf("Treatment: %-5s | Intercept (mu0): %.4f rad | Slope (Beta): %7.4f | Global Rho: %.4f\n", 
              treatments[k], mu0_est, beta_est[k], rho_est))
}


#CARA Adaptive Target Allocation Estimation

q_fun_exact <- function(mu, rho) {
  delta <- pi / 4 #clinical window threshold (45 degrees)
  f <- function(theta) {
    (1 - rho^2) / (2 * pi * (1 + rho^2 - 2 * rho * cos(theta - mu)))
  }
  integrate(f, lower = delta, upper = 2 * pi - delta)$value
}

#Extracting parameter variances (Xi) from the joint parameter Hessian
H <- (op$hessian + t(op$hessian)) / 2
V <- tryCatch(
  solve(H),
  error = function(e) diag(1, 6)
)
Vdiag <- pmax(diag(V)[1:4], 1e-6) #Variance factor for treatment effects (xi_k)

#Computing analytical failure probability q_k at standard baseline (x = 0)
q_est <- numeric(4)
for (k in 1:4) {
  mu_predicted <- (mu0_est + 2 * atan(beta_est[k] * 0)) %% (2 * pi)
  q_est[k]     <- q_fun_exact(mu_predicted, rho_est)
}
q_est <- pmax(q_est, 1e-6)

#CARA allocation calculation: target allocation is proportional to sqrt(xi / q)
p_raw   <- sqrt(Vdiag / q_est)
pi_star <- p_raw / sum(p_raw)

#Organizing findings into table format
results_synchronized <- data.frame(
  Treatment            = treatments,
  Beta_Estimate        = round(beta_est, 4),
  Failure_Prob_qk      = round(q_est, 4),
  Precision_Factor_xik = round(Vdiag, 4),
  Target_Allocation_pi = round(pi_star, 4),
  CR_Allocation        = 0.2500
)

print(results_synchronized)

