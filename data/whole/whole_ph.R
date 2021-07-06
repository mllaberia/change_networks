
library(bipartite)
library(parallel)
setwd("~/Desktop/Planiliza haematocheilus (Lh)/whole")


####----------------------  WHOLE PARASITE COMMUNITY ---------------------------
###---------------------------  SEA OF AZOV   ----------------------------------

azov_whole<- read.table("azov_whole.txt", header = T, row.names = 1)

#------------- REAL VALUES OF THE NETWORK (AZOV) --------------------------------------

t <- proc.time()
network_whole_azov <- networklevel(azov_whole, index = c("connectance", "weighted NODF"))
modularity_whole_azov <- computeModules(azov_whole)@likelihood
proc.time()-t

#save the results
index_whole_azov <- data.frame((t(network_whole_azov)), modularity_whole_azov)
write.csv2(index_whole_azov,file = "realvalues_azov.csv", row.names = FALSE)

#--------------------- BOOTSTRAP (AZOV) ------------------------------------------------

t <- proc.time()

bootSxSP <- function (azov_whole) {
  a <- sample(seq_len(nrow(azov_whole)), nrow(azov_whole), replace = TRUE)
  azov_whole <- azov_whole[a, ]
  return(azov_whole)
}
A_R <- replicate(1000, bootSxSP(azov_whole), simplify = FALSE)

cores <- detectCores()
cl <- makeForkCluster(20)

netlevelazov <- parSapply(cl, A_R, function(x) networklevel (x, index = c("connectance", "weighted NODF")))
modularity_azov<- parSapply(cl, A_R, function(x)computeModules(x)@likelihood)

stopCluster(cl)
proc.time()-t

# save the results 
boot_whole_azov <- data.frame(t(netlevelazov), modularity_azov)
write.csv2(boot_whole_azov, file = "bootnet_whole_azov.csv", row.names = FALSE)

#---------------------  NULL MODEL (AZOV)  --------------------------------------------

t <- proc.time()

## with swap.web algorithm (for WNDOF and modularity)
nullazov_whole_swap <- nullmodel(azov_whole, N=1000, method = 2) #method 2/"swap.web"

## with vaznull algorithm (for WNODF only)
nullazov_whole_vaz <- nullmodel(azov_whole, N=1000, method = 3) #method 3/"vaznull"

cores <- detectCores()
cl <- makeForkCluster(20)

null_netlevelazov_s <- parSapply(cl, nullazov_whole_swap, function(x) networklevel (x, index = c("weighted NODF")))
null_modularazov_s <-parSapply(cl, nullazov_whole_swap, function(x) computeModules(x)@likelihood)

null_netlevelazov_v <- parSapply(cl, nullazov_whole_vaz, function(x) networklevel (x, index = c("weighted NODF")))

stopCluster(cl)
proc.time()-t

# save the results 
null_swap_whole_azov <- data.frame(data.frame(null_netlevelazov_s), null_modularazov_s)
write.csv2(null_swap_whole_azov, file = "null_swap_whole_azov.csv", row.names = FALSE)

null_vaz_whole_azov <- data.frame(null_netlevelazov_v)
write.csv2(null_vaz_whole_azov, file = "null_vaz_whole_azov.csv", row.names = FALSE)

# ------------------------ MEAN AND STANDARD DEVIATION OF REAL VALUES (CONNECTANCE INDEX) (AZOV) -------------------------

# Necessary values: 
## Network real value
index_whole_azov[,1]
## Bootstrap values 
boot_whole_azov[,1]

a_w_value_c <- index_whole_azov[,1]
a_w_mean_c <- mean(boot_whole_azov[,1])
a_w_sd_c <- sd(boot_whole_azov[,1])
a_w_ci_c <- quantile(boot_whole_azov[,1], c(0.025,0.975))

### ------------------------ MEAN AND STANDARD DEVIATION OF STANDARISED VALUES (WNODF AND MODULARITY INDICES) (AZOV) -------------------------

## STANDARDISE THE REAL VALUES OF THE NETWORK, NULLMODEL VALUES AND DETERMINE THE PERCENTILES 

# WNODF values
index_whole_azov [,2] #real value
boot_whole_azov [,2] #bootstrap values 

stan_WNODF  <- function (x) {
  mean <- mean(null_swap_whole_azov[,1])
  sd <- sd(null_swap_whole_azov[,1])
  y <- (x)
  result <- (y-mean)/sd
}

a_w_value_wnodf <- sapply(index_whole_azov[,2], stan_WNODF) #real standardised value 
boot_stanWNODF_azov <- as.data.frame(sapply(boot_whole_azov[,3], stan_WNODF)) #bootstrap standardised values 

hist(boot_stanWNODF_azov)

a_w_mean_wnodf <- mean(boot_stanWNODF_azov)
a_w_sd_wnodf <- sd(boot_stanWNODF_azov)
a_w_ci_wnodf <- quantile(boot_stanWNODF_azov, c(0.025,0.975))

#write.csv2(boot_stan_azov, file = "boot_stan_azov_g.csv")


# Q VALUES 
index_whole_azov [,3] #real value
boot_whole_azov [,3] #bootstrap values 

stan_Q  <- function (x) {
  mean <- mean(null_swap_whole_azov[,2])
  sd <- sd(null_swap_whole_azov[,2])
  y <- (x)
  result <- (y-mean)/sd
}

a_w_value_q <- sapply(index_whole_azov[,3], stan_Q) #real standardised value 
boot_stanQ_azov <- as.data.frame(sapply(boot_whole_azov[,7], stan_Q)) #bootstrap standardised values 

hist(boot_stanQ_azov)

a_w_mean_q <- mean(boot_stanQ_azov)
a_w_sd_q <- sd(boot_stanQ_azov)
a_w_ci_q <- quantile(boot_stanQ_azov, c(0.025,0.975))

#write.csv2(boot_stan_azov, file = "boot_stan_azov_g.csv")



###--------------------------------  SEA OF JAPAN   -------------------------------------------

japan_whole<- read.table("japan_whole.txt", header = T, row.names = 1)

#----------------- REAL VALUES OF THE NETWORK (JAPAN) --------------------------------------

t <- proc.time()
network_whole_japan <- networklevel(japan_whole, index = c("connectance", "weighted NODF"))
modularity_whole_japan <- computeModules(japan_whole)@likelihood
proc.time()-t

#save the results
index_whole_japan <- data.frame((t(data.frame(network_whole_japan))), modularity_whole_japan)
write.csv2(index_whole_japan,file = "realvalues_japan.csv")

#--------------------- BOOTSTRAP (JAPAN) ------------------------------------------------

t <- proc.time()

bootSxSP <- function (japan_whole) {
  a <- sample(seq_len(nrow(japan_whole)), nrow(japan_whole), replace = TRUE)
  japan_whole <- japan_whole[a, ]
  return(japan_whole)
}
A_R <- replicate(1000, bootSxSP(japan_whole), simplify = FALSE)

cores <- detectCores()
cl <- makeForkCluster(20)

netleveljapan <- parSapply(cl, A_R, function(x) networklevel (x, index = c("connectance", "weighted NODF")))
modularity_japan<- parSapply(cl, A_R, function(x)computeModules(x)@likelihood)

stopCluster()
proc.time()-t

# save the results 
boot_whole_japan <- data.frame((t(data.frame(netleveljapan))), modularity_japan)
write.csv2(boot_whole_japan, file = "bootnet_whole_japan.csv", row.names = FALSE)

#---------------------  NULL MODEL (JAPAN)  --------------------------------------------

t <- proc.time()

## with swap.web algorithm (for WNDOF and modularity)
nulljapan_whole_swap <- nullmodel(japan_whole, N=1000, method = 2) #method 2/"swap.web"

## with vaznull algorithm (for WNODF only)
nulljapan_whole_vaz <- nullmodel(japan_whole, N=1000, method = 3) #method 3/"vaznull"

cores <- detectCores()
cl <- makeForkCluster(20)

null_netleveljapan_s <- parSapply(cl, nulljapan_whole_swap, function(x) networklevel (x, index = c("weighted NODF")))
null_modularjapan_s <-parSapply(cl, nulljapan_whole_swap, function(x) computeModules(x)@likelihood)

null_netleveljapan_v <- parSapply(cl, nulljapan_whole_vaz, function(x) networklevel (x, index = c("weighted NODF")))

stopCluster(cl)
proc.time()-t

# save the results 
null_swap_whole_japan <- data.frame(data.frame(null_netleveljapan_s), null_modularjapan_s)
write.csv2(null_swap_whole_japan, file = "null_swap_whole_japan.csv", row.names = FALSE)

null_vaz_whole_japan <- data.frame(null_netleveljapan_v)
write.csv2(null_vaz_whole_japan, file = "null_vaz_whole_japan.csv", row.names = FALSE)

# ------------------------ MEAN AND STANDARD DEVIATION OF REAL VALUES (CONNECTANCE INDEX) (JAPAN) -------------------------

# Necessary values: 
## Network real value
index_whole_japan[,1]
## Bootstrap values 
boot_whole_japan[,1]

j_w_value_c <- index_whole_japan[,1]
j_w_mean_c <- mean(boot_whole_japan[,1])
j_w_sd_c <- sd(boot_whole_japan[,1])
j_w_ci_c <- quantile(boot_whole_japan[,1], c(0.025,0.975))

### ------------------------ MEAN AND STANDARD DEVIATION OF STANDARISED VALUES (WNODF AND MODULARITY INDICES) (JAPAN) -------------------------

## STANDARDISE THE REAL VALUES OF THE NETWORK, NULLMODEL VALUES AND DETERMINE THE PERCENTILES 

# WNODF values
index_whole_japan [,2] #real value
boot_whole_japan [,2] #bootstrap values 

stan_WNODF  <- function (x) {
  mean <- mean(null_swap_whole_japan[,1])
  sd <- sd(null_swap_whole_japan[,1])
  y <- (x)
  result <- (y-mean)/sd
}

j_w_value_wnodf <- sapply(index_whole_japan[,2], stan_WNODF) #real standardised value 
boot_stanWNODF_japan <- as.data.frame(sapply(boot_whole_japan[,3], stan_WNODF)) #bootstrap standardised values 

hist(boot_stanWNODF_japan)

j_w_mean_wnodf <- mean(boot_stanWNODF_japan)
j_w_mean_sd <- sd(boot_stanWNODF_japan)
j_w_ci_wndof <- quantile(boot_stanWNODF_japan, c(0.025,0.975))

#write.csv2(boot_stan_japan, file = "boot_stan_japan_g.csv")


# Q VALUES 
index_whole_japan [,3] #real value
boot_whole_japan [,3] #bootstrap values 

stan_Q  <- function (x) {
  mean <- mean(null_swap_whole_japan[,2])
  sd <- sd(null_swap_whole_japan[,2])
  y <- (x)
  result <- (y-mean)/sd
}

j_w_value_q <- sapply(index_whole_japan[,3], stan_Q) #real standardised value 
boot_stanQ_japan <- as.data.frame(sapply(boot_whole_japan[,7], stan_Q)) #bootstrap standardised values 

hist(boot_stanQ_japan)

j_w_mean_q <- mean(boot_stanQ_japan)
j_w_sd_q <- sd(boot_stanQ_japan)
j_w_ci_q <- quantile(boot_stanQ_japan, c(0.025,0.975))

#write.csv2(boot_stan_japan, file = "boot_stan_japan_g.csv")


#---------------------- MANN - WHITNEY TEST   -----------------------

# CONNECTANCE
wilcox.test(index_whole_azov[,1], index_whole_japan[,1], conf.int = TRUE)

# WNODF
wilcox.test(null_swap_whole_azov[,1], null_swap_whole_japan[,1], conf.int = TRUE)


# MODULARIDAD
wilcox.test(null_swap_whole_azov[,2], null_swap_whole_japan[,2], conf.int = TRUE)

wilcox.test((quantile(boot_stanq_azov, c(0.025,0.975))),(quantile(boot_stanq_japan, c(0.025,0.975))))



# SEA OF AZOV
## C VALUES
a_whole_value_c <- index_whole_azov[,1] #real value
a_whole_mean_c <- mean(boot_whole_azov[,1])
a_whole_sd_c <- sd(boot_whole_azov[,1])
a_whole_ci_c <- quantile(boot_whole_azov[,1], c(0.025,0.975))

## WNODF VALUES 
a_whole_value_wnodf <- sapply(index_whole_azov[,2], stan_WNODF) #real standardized value 
a_whole_mean_wnodf <- mean(boot_stanWNODF_azov)
a_whole_sd_wnodf <- sd(boot_stanWNODF_azov)
a_whole_ci_wnodf <- quantile(boot_stanWNODF_azov, c(0.025,0.975))

## Q VALUES 
a_whole_value_q <- sapply(index_whole_azov[,3], stan_Q) #real standardized value
a_whole_mean_q <- mean(boot_stanQ_azov)
a_whole_sd_q <- sd(boot_stanQ_azov)
a_whole_ci_q <- quantile(boot_stanQ_azov, c(0.025,0.975))


# SEA OF JAPAN 
## C VALUES 
j_whole_value_c <- index_whole_japan[,1] #real value
j_whole_mean_c <- mean(boot_whole_japan[,1])
j_whole_sd_c <- sd(boot_whole_japan[,1])
j_whole_ci_c <- quantile(boot_whole_japan[,1], c(0.025,0.975))

## WNODF VALUES 
j_whole_value_wnodf <- sapply(index_whole_japan[,2], stan_WNODF) #real standardized value 
j_whole_mean_wnodf <- mean(boot_stanWNODF_japan)
j_whole_sd_wndof <- sd(boot_stanWNODF_japan)
j_whole_ci_wndof <- quantile(boot_stanWNODF_japan, c(0.025,0.975))

## Q VALUES 
j_whole_value_q <- sapply(index_whole_japan[,3], stan_Q) #real standardized value 
j_whole_mean_q <- mean(boot_stanQ_japan)
j_whole_sd_q <- sd(boot_stanQ_japan)
j_whole_ci_q <- quantile(boot_stanQ_japan, c(0.025,0.975))


## Create a data frame with all the index results 

a_whole_connectance <- c(a_whole_value_c, a_whole_mean_c, a_whole_sd_c, a_whole_ci_c)
j_whole_connectance <- c(j_whole_value_c, j_whole_mean_c, j_whole_sd_c, j_whole_ci_c)

a_whole_wnodf <-c(a_whole_value_wnodf, a_whole_mean_wnodf, a_whole_sd_wnodf, a_whole_ci_wnodf)
j_whole_wnodf <-c(j_whole_value_wnodf, j_whole_mean_wnodf, j_whole_sd_wnodf, j_whole_ci_wnodf)

a_whole_q <-c(a_whole_value_q, a_whole_mean_q, a_whole_sd_q, a_whole_ci_q)
j_whole_q <- c(j_whole_value_q, j_whole_mean_q, j_whole_sd_q, j_whole_ci_q)

whole_values <- rbind(a_whole_connectance, j_whole_connectance, a_whole_wnodf, j_whole_wnodf, a_whole_q, j_whole_q)
colnames(whole_values) <- c("value", "mean", "sd", "ci")
whole_values

