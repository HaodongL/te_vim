# Check CATE estimates
repo_path = "~/Repo/te_vim/"
source(paste0(repo_path, "R/simu/simu_dgd.R")) #Used for the current examples
source(paste0(repo_path, "R/est_function/sl3_config.R"))
source(paste0(repo_path, "R/est_function/fit_para.R"))
source(paste0(repo_path, "R/est_function/vim.R"))

# Q, g
set.seed(1994)
N <- 5e2 #size of generated data
df <- generate_data_simple(N)

ws = c('X2')
cv = F
dr = T
lfm_linear = FALSE
max.it = 600

Q_bounds = c(0.001, 0.999)
g_bounds = c(0.025, 0.975)
tau_bounds = c(-1+1e-3, 1-1e-3)
tau_s_bounds = c(-1+1e-3, 1-1e-3)
gamma_s_bounds = c(1e-6, 1-1e-6)

y_l <- min(df$Y)
y_u <- max(df$Y)
df$Y <- scale01(df$Y, y_l, y_u)


# Q_bounds = c(-1e4, 1e4)
# g_bounds = c(0.025, 0.975)
# tau_bounds = c(-1e4, 1e4)
# tau_s_bounds = c(-1e4, 1e4)
# gamma_s_bounds = c(-1e4, 1e4)
# y_l <- 0
# y_u <- 1

# fit Q, g
# fit tau, tau_s, gamma_s
all_covar = setdiff(names(df), 'Y')
all_w = setdiff(all_covar, 'A')
all_wsc = setdiff(all_w, ws)
y <- df[["Y"]]
a <- df[["A"]]

res_Qg <- fitSL(df, sl_Q, sl_g, Q_bounds, g_bounds) # big object rm later
Q_fit <- res_Qg$Q_fit # big object rm later
g_fit <- res_Qg$g_fit # big object rm later

Q0_task <- make_sl3_Task(data = df %>% mutate(A=0), covariates = all_covar)
Q1_task <- make_sl3_Task(data = df %>% mutate(A=1), covariates = all_covar)

Qbar0W <- Q_fit$predict(Q0_task)
Qbar1W <- Q_fit$predict(Q1_task)
# bound Q if need
if (!is.null(Q_bounds)){
  Qbar0W <- bound(Qbar0W, Q_bounds)
  Qbar1W <- bound(Qbar1W, Q_bounds)
}
QbarAW <- ifelse(df$A == 1, Qbar1W, Qbar0W)
gn <- g_fit$predict()
# bound g
gn <- bound(gn, g_bounds)


# use true g
# gn <- plogis(0.1*df$X1*df$X2-0.4*df$X1)



po = (y - QbarAW)*(2*a - 1)/gn + Qbar1W - Qbar0W
tau_fit <- fit_x(df = df, sl_x = sl_x, po = po, outcome = 'po', para = 'tau')


# T cate
f_cate_t <- function(w1, w2){
  df <- data.frame("A" = 0, "X1" = w1, "X2" = w2)
  Q0_task <- make_sl3_Task(data = df, covariates = all_covar)
  Q1_task <- make_sl3_Task(data = df %>% mutate(A=1), covariates = all_covar)
  
  Qbar0W <- Q_fit$predict(Q0_task)
  Qbar1W <- Q_fit$predict(Q1_task)
  
  tau <- Qbar1W - Qbar0W
  tau <- bound(tau, tau_bounds)
  tau_t <- tau*(y_u - y_l)
  return(tau_t)
}


# DR cate
f_cate_dr <- function(w1, w2){
  df <- data.frame("A" = 0, "X1" = w1, "X2" = w2)
  Q0_task <- make_sl3_Task(data = df, covariates = all_covar)
  Q1_task <- make_sl3_Task(data = df %>% mutate(A=1), covariates = all_covar)
  
  Qbar0W <- Q_fit$predict(Q0_task)
  Qbar1W <- Q_fit$predict(Q1_task)
  
  tau_task <- make_sl3_Task(data = df, covariates = all_w)
  tau <- tau_fit$predict(tau_task)
  tau <- bound(tau, tau_bounds)
  tau_dr <- tau*(y_u - y_l)
  return(tau_dr)
}


####

####

f_cate <- function (x, y) {
  CATE <- x^2*(x+7/5) + (5*y/3)^2
  return (CATE)
}


# png(filename= paste0("~/Repo/te_vim/tnp/cate_", N,".png"),
#     width = 2048,
#     height = 1024,
#     res = 180,
#     pointsize = 10)
par(mfrow = c(1, 3))


w1 <- seq(-1, 1, length= 30)
w2 <- w1
tau <- outer(w1, w2, f_cate)
# z[is.na(z)] <- 1
op <- par(bg = "white")
persp(w1, w2, tau, theta = 30, phi = 25, expand = 0.5, col = "#669bbc",
      main="True CATE function")


# cate t
w1 <- seq(-1, 1, length= 30)
w2 <- w1
tau_n <- outer(w1, w2, f_cate_t)
mse = round(mean((tau_n - tau)^2), 3)
# z[is.na(z)] <- 1
op <- par(bg = "white")
persp(w1, w2, tau_n, theta = 30, phi = 25, expand = 0.5, col = "#669bbc",
      main= paste0("CATE estimates with T-learner ", 
                   "(MSE = ", mse, ")"))


# cate dr
w1 <- seq(-1, 1, length= 30)
w2 <- w1
tau_n <- outer(w1, w2, f_cate_dr)
mse = round(mean((tau_n - tau)^2), 3)
# z[is.na(z)] <- 1
op <- par(bg = "white")
persp(w1, w2, tau_n, theta = 30, phi = 25, expand = 0.5, col = "#669bbc",
      main= paste0("CATE estimates with DR-learner ", 
                   "(MSE = ", mse, ")"))

# dev.off()



