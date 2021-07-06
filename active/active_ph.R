
library(bipartite)
library(parallel)

##setwd("~/Desktop/Planiliza haematocheilus (Lh)_arreglado/general ")

####----------------------  ACTIVE PARASITE COMMUNITY ---------------------------
###---------------------------  SEA OF AZOV   ----------------------------------

azov_active<- read.table("azov_active.txt", header = T, row.names = 1)

#------------- REAL VALUES OF THE NETWORK (AZOV) --------------------------------------

t <- proc.time()
network_active_azov <- networklevel(azov_active, index = c("connectance", "weighted NODF"))
modularity_active_azov <- computeModules(azov_active)@likelihood
proc.time()-t

#save the results
index_active_azov <- data.frame((t(network_active_azov)), modularity_active_azov)
write.csv2(index_active_azov,file = "realvalues_azov.csv", row.names = FALSE)

#--------------------- BOOTSTRAP (AZOV) ------------------------------------------------

t <- proc.time()

bootSxSP <- function (azov_active) {
  a <- sample(seq_len(nrow(azov_active)), nrow(azov_active), replace = TRUE)
  azov_active <- azov_active[a, ]
  return(azov_active)
}
A_R <- replicate(2, bootSxSP(azov_active), simplify = FALSE)

cores <- detectCores()
cl <- makeForkCluster(cores)

netlevelazov <- parSapply(cl, A_R, function(x) networklevel (x, index = c("connectance", "weighted NODF")))
modularity_azov<- parSapply(cl, A_R, function(x)computeModules(x)@likelihood)

stopCluster(cl)
proc.time()-t

# save the results 
boot_active_azov <- data.frame(t(netlevelazov), modularity_azov)
write.csv2(boot_active_azov, file = "bootnet_active_azov.csv", row.names = FALSE)

#---------------------  NULL MODEL (AZOV)  --------------------------------------------

t <- proc.time()

## with swap.web algorithm (for WNDOF and modularity)
nullazov_active_swap <- nullmodel(azov_active, N=1000, method = 2) #method 2/"swap.web"

## with vaznull algorithm (for WNODF only)
nullazov_active_vaz <- nullmodel(azov_active, N=1000, method = 3) #method 3/"vaznull"

cores <- detectCores()
cl <- makeForkCluster(cores)

null_netlevelazov_s <- parSapply(cl, nullazov_active_swap, function(x) networklevel (x, index = c("weighted NODF")))
null_modularazov_s <-parSapply(cl, nullazov_active_swap, function(x) computeModules(x)@likelihood)

null_netlevelazov_v <- parSapply(cl, nullazov_active_vaz, function(x) networklevel (x, index = c("weighted NODF")))

stopCluster(cl)
proc.time()-t

# save the results 
null_swap_active_azov <- data.frame((t(data.frame(null_netlevelazov_s))), null_modularazov_s)
write.csv2(null_swap_active_azov, file = "null_swap_active_azov.csv", row.names = FALSE)

null_vaz_active_azov <- t(null_netlevelazov_v)
write.csv2(null_vaz_active_azov, file = "null_vaz_active_azov.csv", row.names = FALSE)

# ------------------------ MEAN AND STANDARD DEVIATION OF REAL VALUES (CONNECTANCE INDEX) (AZOV) -------------------------

# Necessary values: 
## Network real value
index_active_azov[,1]
## Bootstrap values 
boot_active_azov[,1]

a_a_value_c <- index_active_azov[,1]
a_a_mean_c <- mean(boot_active_azov[,1])
a_a_sd_c <- sd(boot_active_azov[,1])
a_a_ci_c <- quantile(boot_active_azov[,1], c(0.025,0.975))

### ------------------------ MEAN AND STANDARD DEVIATION OF STANDARISED VALUES (WNODF AND MODULARITY INDICES) (AZOV) -------------------------

## STANDARDISE THE REAL VALUES OF THE NETWORK, NULLMODEL VALUES AND DETERMINE THE PERCENTILES 

# WNODF values
index_active_azov [,2] #real value
boot_active_azov [,2] #bootstrap values 

stan_WNODF  <- function (x) {
  mean <- mean(null_swap_active_azov[,1])
  sd <- sd(null_swap_active_azov[,1])
  y <- (x)
  result <- (y-mean)/sd
}

a_a_value_wnodf <- sapply(index_active_azov[,2], stan_WNODF) #real standardised value 
boot_stanWNODF_azov <- as.data.frame(sapply(boot_active_azov[,3], stan_WNODF)) #bootstrap standardised values 

hist(boot_stanWNODF_azov)

a_a_mean_wnodf <- mean(boot_stanWNODF_azov)
a_a_sd_wnodf <- sd(boot_stanWNODF_azov)
a_a_ci_wnodf <- quantile(boot_stanWNODF_azov, c(0.025,0.975))

#write.csv2(boot_stan_azov, file = "boot_stan_azov_g.csv")


# Q VALUES 
index_active_azov [,3] #real value
boot_active_azov [,3] #bootstrap values 

stan_Q  <- function (x) {
  mean <- mean(null_swap_active_azov[,2])
  sd <- sd(null_swap_active_azov[,2])
  y <- (x)
  result <- (y-mean)/sd
}

a_a_value_q <- sapply(index_active_azov[,3], stan_Q) #real standardised value 
boot_stanQ_azov <- as.data.frame(sapply(boot_active_azov[,7], stan_Q)) #bootstrap standardised values 

hist(boot_stanQ_azov)

a_a_mean_q <- mean(boot_stanQ_azov)
a_a_sd_q <- sd(boot_stanQ_azov)
a_a_ci_q <- quantile(boot_stanQ_azov, c(0.025,0.975))

#write.csv2(boot_stan_azov, file = "boot_stan_azov_g.csv")



###--------------------------------  SEA OF JAPAN   -------------------------------------------

japan_active<- read.table("japan_active.txt", header = T, row.names = 1)

#----------------- REAL VALUES OF THE NETWORK (JAPAN) --------------------------------------

t <- proc.time()
network_active_japan <- networklevel(japan_active, index = c("connectance", "weighted NODF"))
modularity_active_japan <- computeModules(japan_active)@likelihood
proc.time()-t

#save the results
index_active_japan <- data.frame((t(data.frame(network_active_japan))), modularity_active_japan)
write.csv2(index_active_japan,file = "realvalues_japan.csv")

#--------------------- BOOTSTRAP (JAPAN) ------------------------------------------------

t <- proc.time()

bootSxSP <- function (japan_active) {
  a <- sample(seq_len(nrow(japan_active)), nrow(japan_active), replace = TRUE)
  japan_active <- japan_active[a, ]
  return(japan_active)
}
A_R <- replicate(1000, bootSxSP(japan_active), simplify = FALSE)

cores <- detectCores()
cl <- makeForkCluster(cores)

netleveljapan <- parSapply(cl, A_R, function(x) networklevel (x, index = c("connectance", "weighted NODF")))
modularity_japan<- parSapply(cl, A_R, function(x)computeModules(x)@likelihood)

stopCluster()
proc.time()-t

# save the results 
boot_active_japan <- data.frame((t(data.frame(netleveljapan))), modularity_japan)
write.csv2(boot_active_japan, file = "bootnet_active_japan.csv", row.names = FALSE)

#---------------------  NULL MODEL (JAPAN)  --------------------------------------------

t <- proc.time()

## with swap.web algorithm (for WNDOF and modularity)
nulljapan_active_swap <- nullmodel(japan_active, N=1000, method = 2) #method 2/"swap.web"

## with vaznull algorithm (for WNODF only)
nulljapan_active_vaz <- nullmodel(japan_active, N=1000, method = 3) #method 3/"vaznull"

cores <- detectCores()
cl <- makeForkCluster(cores)

null_netleveljapan_s <- parSapply(cl, nulljapan_active_swap, function(x) networklevel (x, index = c("weighted NODF")))
null_modularjapan_s <-parSapply(cl, nulljapan_active_swap, function(x) computeModules(x)@likelihood)

null_netleveljapan_v <- parSapply(cl, nulljapan_active_vaz, function(x) networklevel (x, index = c("weighted NODF")))

stopCluster(cl)
proc.time()-t

# save the results 
null_swap_active_japan <- data.frame((t(data.frame(null_netleveljapan_s))), null_modularjapan_s)
write.csv2(null_swap_active_japan, file = "null_swap_active_japan.csv", row.names = FALSE)

null_vaz_active_japan <- data.frame((t(data.frame(null_netleveljapan_v))), null_modularjapan_s)
write.csv2(null_vaz_active_japan, file = "null_vaz_active_japan.csv", row.names = FALSE)

# ------------------------ MEAN AND STANDARD DEVIATION OF REAL VALUES (CONNECTANCE INDEX) (JAPAN) -------------------------

# Necessary values: 
## Network real value
index_active_japan[,1]
## Bootstrap values 
boot_active_japan[,1]

j_a_value_c <- index_active_japan[,1]
j_a_mean_c <- mean(boot_active_japan[,1])
j_a_sd_c <- sd(boot_active_japan[,1])
j_a_ci_c <- quantile(boot_active_japan[,1], c(0.025,0.975))

### ------------------------ MEAN AND STANDARD DEVIATION OF STANDARISED VALUES (WNODF AND MODULARITY INDICES) (JAPAN) -------------------------

## STANDARDISE THE REAL VALUES OF THE NETWORK, NULLMODEL VALUES AND DETERMINE THE PERCENTILES 

# WNODF values
index_active_japan [,2] #real value
boot_active_japan [,2] #bootstrap values 

stan_WNODF  <- function (x) {
  mean <- mean(null_swap_active_japan[,1])
  sd <- sd(null_swap_active_japan[,1])
  y <- (x)
  result <- (y-mean)/sd
}

j_a_value_wnodf <- sapply(index_active_japan[,2], stan_WNODF) #real standardised value 
boot_stanWNODF_japan <- as.data.frame(sapply(boot_active_japan[,3], stan_WNODF)) #bootstrap standardised values 

hist(boot_stanWNODF_japan)

j_a_mean_wnodf <- mean(boot_stanWNODF_japan)
j_a_mean_sd <- sd(boot_stanWNODF_japan)
j_a_ci_wndof <- quantile(boot_stanWNODF_japan, c(0.025,0.975))

#write.csv2(boot_stan_japan, file = "boot_stan_japan_g.csv")


# Q VALUES 
index_active_japan [,3] #real value
boot_active_japan [,3] #bootstrap values 

stan_Q  <- function (x) {
  mean <- mean(null_swap_active_japan[,2])
  sd <- sd(null_swap_active_japan[,2])
  y <- (x)
  result <- (y-mean)/sd
}

j_a_value_q <- sapply(index_active_japan[,3], stan_Q) #real standardised value 
boot_stanQ_japan <- as.data.frame(sapply(boot_active_japan[,7], stan_Q)) #bootstrap standardised values 

hist(boot_stanQ_japan)

j_a_mean_q <- mean(boot_stanQ_japan)
j_a_sd_q <- sd(boot_stanQ_japan)
j_a_ci_q <- quantile(boot_stanQ_japan, c(0.025,0.975))

#write.csv2(boot_stan_japan, file = "boot_stan_japan_g.csv")


#---------------------- MANN - WHITNEY TEST   -----------------------

# CONNECTANCE
wilcox.test(index_active_azov[,1], index_active_japan[,1], conf.int = TRUE)

# WNODF
wilcox.test(null_swap_active_azov[,1], null_swap_active_japan[,1], conf.int = TRUE)


# MODULARIDAD
wilcox.test(null_swap_active_azov[,2], null_swap_active_japan[,2], conf.int = TRUE)

wilcox.test((quantile(boot_stanq_azov, c(0.025,0.975))),(quantile(boot_stanq_japan, c(0.025,0.975))))



# SEA OF AZOV
## C VALUES 
a_active_value_c <- index_active_azov[,1] #real value
a_active_mean_c <- mean(boot_active_azov[,1])
a_active_sd_c <- sd(boot_active_azov[,1])
a_active_ci_c <- quantile(boot_active_azov[,1], c(0.025,0.975))

## WNODF VALUES 
a_active_value_wnodf <- sapply(index_active_azov[,2], stan_WNODF) #real standardized value 
a_active_mean_wnodf <- mean(boot_stanWNODF_azov)
a_active_sd_wnodf <- sd(boot_stanWNODF_azov)
a_active_ci_wnodf <- quantile(boot_stanWNODF_azov, c(0.025,0.975))

## Q VALUES 
a_active_value_q <- sapply(index_active_azov[,3], stan_Q) #real standardized value
a_active_mean_q <- mean(boot_stanQ_azov)
a_active_sd_q <- sd(boot_stanQ_azov)
a_active_ci_q <- quantile(boot_stanQ_azov, c(0.025,0.975))


# SEA OF JAPAN 
## C VALUES 
j_active_value_c <- index_active_japan[,1] #real value
j_active_mean_c <- mean(boot_active_japan[,1])
j_active_sd_c <- sd(boot_active_japan[,1])
j_active_ci_c <- quantile(boot_active_japan[,1], c(0.025,0.975))

## WNODF VALUES 
j_active_value_wnodf <- sapply(index_active_japan[,2], stan_WNODF) #real standardized value 
j_active_mean_wnodf <- mean(boot_stanWNODF_japan)
j_active_sd_wndof <- sd(boot_stanWNODF_japan)
j_active_ci_wndof <- quantile(boot_stanWNODF_japan, c(0.025,0.975))

## Q VALUES 
j_active_value_q <- sapply(index_active_japan[,3], stan_Q) #real standardized value 
j_active_mean_q <- mean(boot_stanQ_japan)
j_active_sd_q <- sd(boot_stanQ_japan)
j_active_ci_q <- quantile(boot_stanQ_japan, c(0.025,0.975))


## Create a data frame with all the index results 

a_active_connectance <- c(a_active_value_c, a_active_mean_c, a_active_sd_c, a_active_ci_c)
j_active_connectance <- c(j_active_value_c, j_active_mean_c, j_active_sd_c, j_active_ci_c)

a_active_wnodf <-c(a_active_value_wnodf, a_active_mean_wnodf, a_active_sd_wnodf, a_active_ci_wnodf)
j_active_wnodf <-c(j_active_value_wnodf, j_active_mean_wnodf, j_active_sd_wnodf, j_active_ci_wnodf)

a_active_q <-c(a_active_value_q, a_active_mean_q, a_active_sd_q, a_active_ci_q)
j_active_q <- c(j_active_value_q, j_active_mean_q, j_active_sd_q, j_active_ci_q)

active_values <- rbind(a_active_connectance, j_active_connectance, a_active_wnodf, j_active_wnodf, a_active_q, j_active_q)
colnames(active_values) <- c("value", "mean", "sd", "ci")
active_values

