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
library(gemtc)


###### helper functions 
contains <- function(x,y) {
  x %in% y
}

# Formula to calculate AICc by given rss
AICc_rss = function(rss, K, n){
  2 * K + rss + 2 * K *(K + 1)/(n - K - 1) 
}

# Formula to calculate BICc by given rss
BICc_rss = function(rss, K, n){
  log(n) * K + rss + 2 * K *(K + 1)/(n - K - 1) 
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

# Return the plot of fitted model by adding an extra 0 line with cubic scale
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

best_row_aic = function(rss, Ks, n){
  best_value = NA
  row = 0
  # total number of rows
  n = length(rss)
  for (i in 1:n){
    K = Ks[i]
    AICc = AICc_rss(rss[i], K, n) 
    # first round
    if (is.na(best_value)){
      best_value = AICc
      row = i
    } else if (AICc < best_value){
      best_value = AICc
      row = i
    }
  }
  row
}



# Return the row that have the smallest BICc
#   rss: residual sum of square; 
#        length(rss) = total number of rows
#   K: number of model parameters; same length as `rss`

best_row_bic = function(rss, Ks, n){
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
# Treatments 1,2 and ref names given in string form
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


# Return a matrix X with (number of studies) x (number of treatments)
#   Treatments 1,2 and ref names given in int or double form
X_matrix_num <- function(t1, t2, ref, n){
  X = matrix(0, nrow = length(t1), ncol = n)
  for (i in 1:nrow(X)){
    X[i, c(t1[i], t2[i])] = c(1, -1)
  }
  colnames(X) = seq(1, n, 1)
  X = X[, -ref]
  X
}


# Given total number of distinct treatments, return the corresponding 
#   vector to draw the fully connected graph.
#   i.e given n = 5, it will return c(1,2,1,3,1,4,2,3,2,4,3,4)

n_treats_graph_generator <- function(n){
  if(n <= 2) {return(1)}
  
  n = n - 1
  e = c()
  for (i in 1:n){
    for (j in (i + 1):n){
      e = append(e, c(i, j))
    }
  }
  e = head(e, -4)
  e
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

##### helper for simulation study 
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


# A helper function that can transfer the info in cells 
#   of the netleague from string into numbers in form of 
#   c(estimate, left ci bound, right ci bound)
league_str_num <- function(a){
  # find the index of all whitespaces
  w_space_i = unlist(gregexpr(' ', a))
  # index of character where end the estimator 
  est_end = w_space_i[1] - 1
  # index of character where start and end the left ci bound
  ci_l_start = w_space_i[1] + 2
  if(length(w_space_i) == 4){ci_l_end = w_space_i[3] - 2}else{ci_l_end = w_space_i[2] - 2}
  # index of character where start and end the right ci bound
  ci_r_start = w_space_i[length(w_space_i)] + 1
  ci_r_end = nchar(a) - 1
  
  c(as.numeric(substr(a,1, est_end)),
    as.numeric(substr(a,ci_l_start, ci_l_end)), 
    as.numeric(substr(a,ci_r_start, ci_r_end)))
}

simulti_p_aicc <- function(X, bb, fac){
  ### e
  # generate 1000 epsilon
  size = 1000                                      
  meanvector = rep(0, nrow(X))                                 
  # matrix with all var on diagonal
  var_cov = diag((pksn_multi$seTE*fac)^2)
  # find the indexes that stands for multi-arm 
  multi_index = which(pksn_multi$studlab == unique(as.vector(p_study$studlab)))
  # update the cov terms in the var_cov matrix 
  var_cov[multi_index[1],multi_index[2]] = cov_Guttman1997 * (fac^2)
  var_cov[multi_index[2],multi_index[1]] = cov_Guttman1997 * (fac^2)
  e <- mvrnorm(n = size, mu = meanvector, Sigma = var_cov)
  
  ### y
  y_ = matrix(rep(as.numeric(X %*% bb), 1000), ncol = nrow(X), byrow=TRUE) + e
  
  sim_sum = data.frame()
  for (i in 1:1000) {
    y = y_[i,]
    n = length(y)
    K = ncol(X)
    # Choleski Decomposition of var_cov matrix
    w = chol(solve(var_cov))
    gr <- graph(c(1,2,1,3,1,4,2,3,2,4,3,4), directed=FALSE)
    fit = fusedlasso(y=w%*%y, X=w%*%X, graph = gr, gamma=1)
    # fusedLasso_gr(fit) (graph code)
    
    ### df, lambda, rss, AICc
    rss = summary(fit)[,3]
    Ks = summary(fit)[,1]
    beta = fit$beta[, best_row_aic(rss, Ks, n)]
    # Calculate the rss for the full model 
    mod = glm(w%*%y ~ w%*%X - 1)
    rss_full = deviance(mod)
    # update variables 
    rss = c(rss, rss_full)
    Ks = c(Ks, K)
    # adding the full model to the summary
    sum_table = rbind(summary(fit), c(K, 0, rss_full))
    # Add corresponding AICc values into the table 
    sum_table_a = data.frame(cbind(sum_table, AICc = AICc_rss(rss, Ks, n)))
    rownames(sum_table_a) = 1:nrow(sum_table_a)
    
    
    # based on each case, give the correct df
    #   (0,0,0,0), all same, correct df = 0
    #   (0.5,1,3,5), all diff, correct df = 4
    #   (1,1,-2,-2), correct df = 2
    #   (1,1,1,-2), correct df = 1
    if (all(bb == rep(0, 4))){correct_df = 0}
    if (all(bb == c(0.5,1,3,5))){correct_df = 4}
    if (all(bb == c(1,1,-2,-2))){correct_df = 2}
    if (all(bb == c(1,1,1,-2))){correct_df = 2}
    
    #### delta_AIC 
    #   (the difference between the correct model and the selected one)
    # if df = 0 gives the smallest AICc, then only the smallest lambda is added
    select_index = max(which(sum_table_a$AICc == min(sum_table_a$AICc)))
    select = sum_table_a[select_index,]
    # under the correct df, find the one that has the smallest AICc 
    df_index = which(sum_table_a$df == correct_df)
    df_rows = sum_table_a[df_index,]
    correct = df_rows[max(which(df_rows$AICc == min(df_rows$AICc))),]
    delta_AICc = abs(select$AICc - correct$AICc)
    
    #### Path check
    # based on each case, give the correct pooling pairs
    if (all(bb == c(0.5,1,3,5))){cpair = c()}
    if (all(bb == c(1,1,-2,-2))){cpair = c(23, 45)}
    if (all(bb == c(1,1,1,-2))){cpair = c(23, 24, 34)}
    # under correct df, see if there exist a model that has correct path
    # path indicator (path_indicator) is set to be no as default
    path_indicator = "no"
    if (correct_df == 0 | correct_df == 4){path_indicator = "yes"}
    # (1,1,-2,-2)
    if (all(bb == c(1,1,-2,-2))){
      j = 1
      while (path_indicator == "no" & j <= length(df_index)) {
        pool_pairs_model = pairwise_pool_names(df_index[j], fit)
        #print(i)
        #print(pool_pairs_model)
        if (length(pool_pairs_model) == 2){
          if (all(pool_pairs_model == cpair)){path_indicator = "yes"}
        }
        j = j + 1
      }
    }
    # (1,1,1,-2)
    if (all(bb == c(1,1,1,-2))){
      j = 1
      while (path_indicator == "no" & j <= length(df_index)) {
        pool_pairs_model = pairwise_pool_names(df_index[j], fit)
        #print(i)
        #print(pool_pairs_model)
        if (length(pool_pairs_model) == 3){
          if (all(pool_pairs_model == cpair)){path_indicator = "yes"}
        }
        j = j + 1
      }
    }
    
    
    ### false pooling within the best model 
    pair_values = rep(0, 10)
    names = c(12, 13, 14, 15, 23, 24, 25, 34, 35, 45)
    pool_best = pairwise_pool_names(select_index, fit)
    # indicator of FP, default is no, which means no FP
    best_indicator = "no"
    if (correct_df == 0){
      best_indicator = "no"
    }else if (correct_df == 4){
      if (!is.null(pool_best)){best_indicator = "yes"}
    }else{
      ss = append(cpair, pool_best)
      # if further pooling are there compared to the correct pairs
      if (length(unique(ss)) > length(cpair)) {best_indicator = "yes"}
    }
    # update 1-2 1-3 1-4 1-5 2-3 2-4 2-5 3-4 3-5 4-5
    for (k in pool_best){
      pair_values[which(names == k)] = 1
    }
    pair_values = matrix(pair_values, nrow = 1)
    
    
    
    ### (estimated dd - true dd)^2 for best model and full model 
    # true values 
    true_dd = gene_dd(bb)
    # full model 
    esti_full = gene_dd(as.vector(coef(mod)))
    full_diff2 = matrix((true_dd - esti_full)^2, nrow = 1)
    # best model 
    if (select_index == nrow(sum_table_a)){ #full model is selected
      best_diff2 = full_diff2
    }else{
      esti_best = gene_dd(fit$beta[,select_index])
      best_diff2 = matrix((true_dd - esti_best)^2, nrow = 1)
    }
    
    ########## check 1-5
    est_full_1_5 = as.vector(coef(mod))[4]
    if (select_index == nrow(sum_table_a)){ #full model is selected
      est_best_1_5 = est_full_1_5
    }else{
      est_best_1_5 = fit$beta[,select_index][4]
    }
    
    real_1_5 = bb[4]
    
    
    
    # add the row that has the smallest AICc values into the summary table
    sim_sum = rbind(sim_sum, 
                    cbind(select, correct, delta_AICc, path_indicator, 
                          best_indicator, pair_values, full_diff2, best_diff2,
                          est_full_1_5, est_best_1_5, real_1_5))
    if(i%%200==0){print(i)}
  }
  rownames(sim_sum) = 1:nrow(sim_sum)
  colnames(sim_sum)[c(1, 5)] = c("select_df", "correct_df")
  colnames(sim_sum)[c(12:21)] = c("1-2_poll","1-3", "1-4", "1-5", "2-3", "2-4",
                                  "2-5", "3-4", "3-5", "4-5_poll")
  colnames(sim_sum)[c(22:31)] = c("1-2_full","1-3", "1-4", "1-5", "2-3", "2-4",
                                  "2-5", "3-4", "3-5", "4-5_full")
  colnames(sim_sum)[c(32:41)] = c("1-2_best","1-3", "1-4", "1-5", "2-3", "2-4",
                                  "2-5", "3-4", "3-5", "4-5_best")
  
  
  sim_sum
}




simulti_p_bicc <- function(X, bb, fac){
  ### e
  # generate 1000 epsilon
  size = 1000                                      
  meanvector = rep(0, nrow(X))                                 
  # matrix with all var on diagonal
  var_cov = diag((pksn_multi$seTE*fac)^2)
  # find the indexes that stands for multi-arm 
  multi_index = which(pksn_multi$studlab == unique(as.vector(p_study$studlab)))
  # update the cov terms in the var_cov matrix 
  var_cov[multi_index[1],multi_index[2]] = cov_Guttman1997 * (fac^2)
  var_cov[multi_index[2],multi_index[1]] = cov_Guttman1997 * (fac^2)
  e <- mvrnorm(n = size, mu = meanvector, Sigma = var_cov)
  
  ### y
  y_ = matrix(rep(as.numeric(X %*% bb), 1000), ncol = nrow(X), byrow=TRUE) + e
  
  sim_sum = data.frame()
  for (i in 1:1000) {
    y = y_[i,]
    n = length(y)
    K = ncol(X)
    # Choleski Decomposition of var_cov matrix
    w = chol(solve(var_cov))
    gr <- graph(c(1,2,1,3,1,4,2,3,2,4,3,4), directed=FALSE)
    fit = fusedlasso(y=w%*%y, X=w%*%X, graph = gr, gamma=1)
    # fusedLasso_gr(fit) (graph code)
    
    ### df, lambda, rss, BICc
    rss = summary(fit)[,3]
    Ks = summary(fit)[,1]
    beta = fit$beta[, best_row_bic(rss, Ks, n)]
    # Calculate the rss for the full model 
    mod = glm(w%*%y ~ w%*%X - 1)
    rss_full = deviance(mod)
    # update variables 
    rss = c(rss, rss_full)
    Ks = c(Ks, K)
    # adding the full model to the summary
    sum_table = rbind(summary(fit), c(K, 0, rss_full))
    # Add corresponding BICc values into the table 
    sum_table_a = data.frame(cbind(sum_table, BICc = BICc_rss(rss, Ks, n)))
    rownames(sum_table_a) = 1:nrow(sum_table_a)
    
    
    # based on each case, give the correct df
    #   (0,0,0,0), all same, correct df = 0
    #   (0.5,1,3,5), all diff, correct df = 4
    #   (1,1,-2,-2), correct df = 2
    #   (1,1,1,-2), correct df = 1
    if (all(bb == rep(0, 4))){correct_df = 0}
    if (all(bb == c(0.5,1,3,5))){correct_df = 4}
    if (all(bb == c(1,1,-2,-2))){correct_df = 2}
    if (all(bb == c(1,1,1,-2))){correct_df = 2}
    
    #### delta_BIC 
    #   (the difference between the correct model and the selected one)
    # if df = 0 gives the smallest BICc, then only the smallest lambda is added
    select_index = max(which(sum_table_a$BICc == min(sum_table_a$BICc)))
    select = sum_table_a[select_index,]
    # under the correct df, find the one that has the smallest BICc 
    df_index = which(sum_table_a$df == correct_df)
    df_rows = sum_table_a[df_index,]
    correct = df_rows[max(which(df_rows$BICc == min(df_rows$BICc))),]
    delta_BICc = abs(select$BICc - correct$BICc)
    
    #### Path check
    # based on each case, give the correct pooling pairs
    if (all(bb == c(0.5,1,3,5))){cpair = c()}
    if (all(bb == c(1,1,-2,-2))){cpair = c(23, 45)}
    if (all(bb == c(1,1,1,-2))){cpair = c(23, 24, 34)}
    # under correct df, see if there exist a model that has correct path
    # path indicator (path_indicator) is set to be no as default
    path_indicator = "no"
    if (correct_df == 0 | correct_df == 4){path_indicator = "yes"}
    # (1,1,-2,-2)
    if (all(bb == c(1,1,-2,-2))){
      j = 1
      while (path_indicator == "no" & j <= length(df_index)) {
        pool_pairs_model = pairwise_pool_names(df_index[j], fit)
        #print(i)
        #print(pool_pairs_model)
        if (length(pool_pairs_model) == 2){
          if (all(pool_pairs_model == cpair)){path_indicator = "yes"}
        }
        j = j + 1
      }
    }
    # (1,1,1,-2)
    if (all(bb == c(1,1,1,-2))){
      j = 1
      while (path_indicator == "no" & j <= length(df_index)) {
        pool_pairs_model = pairwise_pool_names(df_index[j], fit)
        #print(i)
        #print(pool_pairs_model)
        if (length(pool_pairs_model) == 3){
          if (all(pool_pairs_model == cpair)){path_indicator = "yes"}
        }
        j = j + 1
      }
    }
    
    
    ### false pooling within the best model 
    pair_values = rep(0, 10)
    names = c(12, 13, 14, 15, 23, 24, 25, 34, 35, 45)
    pool_best = pairwise_pool_names(select_index, fit)
    # indicator of FP, default is no, which means no FP
    best_indicator = "no"
    if (correct_df == 0){
      best_indicator = "no"
    }else if (correct_df == 4){
      if (!is.null(pool_best)){best_indicator = "yes"}
    }else{
      ss = append(cpair, pool_best)
      # if further pooling are there compared to the correct pairs
      if (length(unique(ss)) > length(cpair)) {best_indicator = "yes"}
    }
    # update 1-2 1-3 1-4 1-5 2-3 2-4 2-5 3-4 3-5 4-5
    for (k in pool_best){
      pair_values[which(names == k)] = 1
    }
    pair_values = matrix(pair_values, nrow = 1)
    
    
    
    ### (estimated dd - true dd)^2 for best model and full model 
    # true values 
    true_dd = gene_dd(bb)
    # full model 
    esti_full = gene_dd(as.vector(coef(mod)))
    full_diff2 = matrix((true_dd - esti_full)^2, nrow = 1)
    # best model 
    if (select_index == nrow(sum_table_a)){ #full model is selected
      best_diff2 = full_diff2
    }else{
      esti_best = gene_dd(fit$beta[,select_index])
      best_diff2 = matrix((true_dd - esti_best)^2, nrow = 1)
    }
    
    ########## check 1-5
    est_full_1_5 = as.vector(coef(mod))[4]
    if (select_index == nrow(sum_table_a)){ #full model is selected
      est_best_1_5 = est_full_1_5
    }else{
      est_best_1_5 = fit$beta[,select_index][4]
    }
    
    real_1_5 = bb[4]
    
    
    
    # add the row that has the smallest BICc values into the summary table
    sim_sum = rbind(sim_sum, 
                    cbind(select, correct, delta_BICc, path_indicator, 
                          best_indicator, pair_values, full_diff2, best_diff2,
                          est_full_1_5, est_best_1_5, real_1_5))
    if(i%%200==0){print(i)}
  }
  rownames(sim_sum) = 1:nrow(sim_sum)
  colnames(sim_sum)[c(1, 5)] = c("select_df", "correct_df")
  colnames(sim_sum)[c(12:21)] = c("1-2_poll","1-3", "1-4", "1-5", "2-3", "2-4",
                                  "2-5", "3-4", "3-5", "4-5_poll")
  colnames(sim_sum)[c(22:31)] = c("1-2_full","1-3", "1-4", "1-5", "2-3", "2-4",
                                  "2-5", "3-4", "3-5", "4-5_full")
  colnames(sim_sum)[c(32:41)] = c("1-2_best","1-3", "1-4", "1-5", "2-3", "2-4",
                                  "2-5", "3-4", "3-5", "4-5_best")
  
  
  sim_sum
}

