###### library
library(Matrix)
library(igraph)
library(genlasso)
library(dplyr)
library(meta)
library(netmeta)
library(ggplot2)
library(scales)
library(MASS)


###### helper functions 

# Formula to calculate AICc by given rss
AICc_rss = function(rss, K, n){
  2 * K + n * log(rss/n) + 2 * K *(K + 1)/(n - K - 1) 
}

# Formula to calculate BICc by given rss
BICc_rss = function(rss, K, n){
  log(n) * K + n * log(rss/n) + 2 * K *(K + 1)/(n - K - 1) 
}


# Return the plot of fitted model by adding an extra 0 line
fusedLasso_gr = function(fit){
  fit$beta <- rbind(rep(0,  dim(fit$beta)[2]), fit$beta)
  fit$bls <- c(0, fit$bls)
  # fit$lambda <- log(fit$lambda)
  # fit$lambda <- (fit$lambda)^0.5
  cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73","#CC79A7")
  plot(fit, numbers = TRUE, ylab = "Coordinates of d", col = cbp1)
  # plot(fit, numbers = TRUE, xlab = expression(paste("log(", lambda,") ")))
  # plot(fit, numbers = TRUE, xlab = expression(sqrt(lambda)))
}


# Return the plot of fitted model by adding an extra 0 line with log scale
fusedLasso_gr_log = function(fit){
  fit$beta <- rbind(rep(0,  dim(fit$beta)[2]), fit$beta)
  fit$bls <- c(0, fit$bls)
  fit$lambda <- log(fit$lambda)
  plot(fit, numbers = TRUE, xlab = expression(paste("log(", lambda,") ")), 
       ylab = "Coordinates of d")
}

# Return the plot of fitted model by adding an extra 0 line with log scale
fusedLasso_gr_cubic = function(fit){
  fit$beta <- rbind(rep(0,  dim(fit$beta)[2]), fit$beta)
  fit$bls <- c(0, fit$bls)
  fit$lambda <- (fit$lambda)^(1/3)
  plot(fit, numbers = TRUE, xlab = expression(lambda^{frac(1, 3)}), 
       ylab = "Coordinates of d")
}


# Return the row that have the smallest AICc
#   rss: residual sum of square; 
#        length(rss) = total number of rows
#   K: number of model parameters; same length as `rss`

best_row = function(rss, Ks, n){
  best_value = NA
  row = 0
  # total number of rows
  n = length(rss)
  for (i in 1:n){
    K = Ks[i]
    BICc = BICc_rss(rss[i], K, n) 
    # first round
    if (is.na(best_value)){
      best_value = BICc
      row = i
    } else if (BICc < best_value){
      best_value = BICc
      row = i
    }
  }
  row
}


# Return a matrix X with 
#   nrow = length(treatment1), ncol = length(total_levels)
#   according to the contrast-based info
X_matrix = function(treatment1, treatment2, ref, n) {
  
  # factorize the treatment
  treat1 = factor(treatment1)
  treat2 = factor(treatment2)
  
  # calculate the total categories in treatment 1 & 2
  total_levels = sort(unique(c(levels(treat1), levels(treat2))))
  
  # create a 0-matrix by given m x n
  X = matrix(0, nrow = length(treatment1), ncol = length(total_levels))
  
  # find the clo index of the reference
  ref_index = which(total_levels == ref)
  
  # t1: col position of treatment 1 (need to be -1)
  t1 = factor(treatment1, levels = total_levels, labels = seq(1, n, 1))
  # t2: col position of treatment 2 (need to be 1)
  t2 = factor(treatment2, levels = total_levels, labels = seq(1, n, 1))
  # insert 1, -1 into the 0-matrix
  for (i in 1:nrow(X)){
    X[i, c(t1[i], t2[i])] = c(1, -1)
  }
  
  # rename the columns
  colnames(X) = total_levels
  
  # remove the ref row in X
  X = X[, -ref_index]
  
  return(X)
}


# Return the study that is multi-arm
multi_study = function(dataset){
  dataset$studlab = as.factor(dataset$studlab)
  # gives name of all labs
  labs = levels(dataset$studlab)
  output = data.frame()
  for (i in labs){
    study_arms = filter(dataset, studlab == i)
    # find the multi-arm study, add to the output data frame
    if (nrow(study_arms) > 2){output = rbind(output, study_arms)}
  }
  output
}


# return the name of pairs (numbers) that is pooled by the given
#    lambda (index) and fused lasso (fit) model
pairwise_pool_names <- function(index, fit){
  if (index > ncol(fit$beta)) {return(c())}
  zero_i = c()
  values = gene_dd(fit$beta[,index])
  names = c(12, 13, 14, 15, 23, 24, 25, 34, 35, 45)
  for (i in 1:length(values)) {
    if (values[i] == 0){
      zero_i = append(zero_i, names[i])
    }
  }
  zero_i
}


# return a vector that includes all pairwise d_n given d_1, 
#   n = number of treatments in the study.
#   e.g. given c(d12,d13,d14), will return c(d12,d13,d14,d23,d24,d34)
gene_dd = function(d1k){
  n = length(d1k) + 1
  dd = c()
  for (i in 1:(n - 2)){
    di = c()
    for (j in (i + 1):(n - 1)){
      di = append(di, d1k[j] - d1k[i])
    }
    dd = append(dd, di)
  }
  # adding all the d_1k 
  dd = append(dd, d1k, 0)
  dd
}

prob_sum<- function(sim_sum){
  n = 1000
  
  ### P(||I||_1 > 1)
  I_geq1_yes = sim_sum[which(sim_sum$path_indicator == "yes"),]
  I_geq1_no = sim_sum[which(sim_sum$path_indicator == "no"),]
  p_I_yes = nrow(I_geq1_yes)/1000
  p_I_no = nrow(I_geq1_no)/1000
  
  ### P(FP_iota = 1 | ||I||_1 > 1)
  # P(yes|yes)
  FP_iota_1_ygy = I_geq1_yes[which(I_geq1_yes$best_indicator == "yes"),]
  p_fp_ygy = nrow(FP_iota_1_ygy)/nrow(I_geq1_yes)
  # P(no|yes)
  FP_iota_1_ngy = I_geq1_yes[which(I_geq1_yes$best_indicator == "no"),]
  p_fp_ngy = nrow(FP_iota_1_ngy)/nrow(I_geq1_yes)
  # P(yes|no)
  FP_iota_1_ygn = I_geq1_no[which(I_geq1_no$best_indicator == "yes"),]
  p_fp_ygn = nrow(FP_iota_1_ygn)/nrow(I_geq1_no)
  # P(no|no)
  FP_iota_1_ngn = I_geq1_no[which(I_geq1_no$best_indicator == "no"),]
  p_fp_ngn = nrow(FP_iota_1_ngn)/nrow(I_geq1_no)
  
  p_ygy = p_fp_ygy * p_I_yes
  p_ngy = p_fp_ngy * p_I_yes
  p_ygn = p_fp_ygn * p_I_no
  p_ngn = p_fp_ngn * p_I_no
  
  ### update data
  yes_yes = c(p_I_yes, p_fp_ygy, p_ygy)
  yes_no = c(p_I_yes, p_fp_ngy, p_ngy)
  no_yes = c(p_I_no, p_fp_ygn, p_ygn)
  no_no = c(p_I_no, p_fp_ngn, p_ngn)
  
  
  ### pairwise prob
  if (p_ygy == 0|is.na(p_ygy)){yes_yes = append(yes_yes, rep(0,10))} 
  if(p_ngy == 0|is.na(p_ngy)){yes_no = append(yes_no, rep(0,10))} 
  if(p_ygn == 0|is.na(p_ygn)){no_yes = append(no_yes, rep(0,10))} 
  if(p_ngn == 0|is.na(p_ngn)){no_no = append(no_no, rep(0,10))} 
  for (i in 12:21){
    if(length(yes_yes) < 13){
      yes_yes = append(yes_yes, sum(FP_iota_1_ygy[,i])/nrow(FP_iota_1_ygy))
    }
    if(length(yes_no) < 13){
      yes_no = append(yes_no, sum(FP_iota_1_ngy[,i])/nrow(FP_iota_1_ngy))
    }
    if(length(no_yes) < 13){
      no_yes = append(no_yes, sum(FP_iota_1_ygn[,i])/nrow(FP_iota_1_ygn))
    }
    if(length(no_no) < 13){
      no_no = append(no_no, sum(FP_iota_1_ngn[,i])/nrow(FP_iota_1_ngn))
    }
  }
  
  result = data.frame(rbind(yes_yes, yes_no, no_yes, no_no))
  colnames(result) = c("P(||I||_1 > 1)", "P(FP_iota = 1 | ||I||_1 > 1)", "p", 
                       "1-2", "1-3", "1-4", "1-5", "2-3", "2-4", "2-5", "3-4", "3-5", "4-5")
  result
}






