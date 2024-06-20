
# Given the estimated beta, return group of treatments that have
#   the same effect
group_dd <- function(est_beta){
  # find the estimate treatment diff between all pairs of treatment
  est_dd = gene_dd(est_beta)
  # find the corresponding treatment numbers 
  pair_treats = n_treats_graph_generator(length(est_beta) + 2)
  t1 = 1
  t2 = 2
  p_list = c()
  # p_list in form of c(1,2,1,3,2,3) if no diff between 1,2,3
  for(i in 1:length(est_dd)){
    if(est_dd[i] == 0){
      p_list = append(p_list, pair_treats[c(t1, t2)]) 
    }
    # each time goes up two numbers since how the pair_treats is defined
    t1 = t1 + 2
    t2 = t2 + 2
  }
  
  # group all treatments that have same effect
  group_index = 1
  group = list()
  # Add the first pair of pooled treatments
  group[[1]] = p_list[1:2]
  
  # if more than 1 pair
  if(length(p_list) > 2){
    t1 = 3
    t2 = 4
    # if the smaller value in the pair is in the group, 
    # then add the lager one as well
    for (i in 2:(length(p_list)/2)){
      for (index in 1:group_index) {
        if (p_list[t1] %in% group[[index]]){
          group[[index]] = unique(append(group[[index]], p_list[t2]))
          # only create a new subgroup if the number is in none of the sublist
        }else if(!any(sapply(group, contains, x = p_list[t1]))){
          group_index = group_index + 1
          group[[group_index]] = p_list[t1:t2]
        }
      }
      # each time goes up two numbers since how the p_list is defined
      t1 = t1 + 2
      t2 = t2 + 2
      }
  }
  #print(group)
  group
}




# Given dataf is a generate data, with TE, seTE, studlab, t1, t2
#   n is the total number of distinct treatments.
# Return a list with two elements. The first one is a list indicate group 
#   of pooled treatments or list() if no treatments are pooled. The second 
#   element is a NMA model after pool or 0 (if all treatments are pooled)

GLF_aicc_2arm <- function(dataf, N.treat){
  y = dataf$TE
  sigma = dataf$seTE
  X = X_matrix_num(dataf$t1, dataf$t2, ref = 1, n = N.treat) 
  w = diag(1/sigma)
  n = length(y)
  K = ncol(X)
  e = n_treats_graph_generator(N.treat)
  gr = graph(e, directed=FALSE)
  fit = fusedlasso(y=w%*%y, X=w%*%X, graph = gr, gamma=1)
  #fusedLasso_gr(fit)
  #print(dataf)
  
  ### df, lambda, rss, AICc
  rss = summary(fit)[,3]
  Ks = summary(fit)[,1]
  
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
  
  min_row = max(which(sum_table_a$AICc == min(sum_table_a$AICc)))
  
  
  
  
  
  # if full model is selected, then no treatments are pooled
  if(sum_table_a$lambda[min_row] == 0){
    list(list(), netmeta(TE, seTE, t1, t2, studlab,
                    data = dataf, sm = "MD", reference = 1))
    
    
    
    
    
    # rename the pooled treatments, and create new NMA model
  }else{
    best_beta = fit$beta[, min_row]
    t_group_list = group_dd(best_beta)
    
    update = dataf
    # rename treatments that are pooled as the smallest number in the group
    for (k in 1:length(t_group_list)){
      # change t1
      c1 = sapply(dataf$t1, contains, y = t_group_list[[k]])
      if(any(c1)){
        update$t1[which(c1)] = min(t_group_list[[k]])
      }
      # change t2
      c2 = sapply(dataf$t2, contains, y = t_group_list[[k]])
      if(any(c2)){
        update$t2[which(c2)] = min(t_group_list[[k]])
      }
    }
    
    # drop rows that compare the same treatments (t1=t2)
    update = update[-which(update$t1 == update$t2),]
    
    ### create new NMA model based on the updated info
    # if all treatments are pooled
    if(nrow(update) == 0){
      return(list(t_group_list, 0))
    }else{
      list(t_group_list, netmeta(TE, seTE, t1, t2, studlab,
                                 data = update, sm = "MD", reference = 1))
    }
    
  }
}




# Given smry_info, like summary(net)$common$TE, list of pooled 
#   treatments and number of treatments, return all d_1i's
t_vs_ref = function(smry_info, N.treat, group){
  gene = rep(0, N.treat)
  # treatments that are given in the nma model
  nma_t_index = as.numeric(rownames(data.frame(smry_info[1, ])))[-1]
  gene[nma_t_index] = smry_info[1, ][-1]
  
  for (subgroup in group) {
    gene[subgroup] = gene[min(subgroup)]
  }
  
  # remove d11
  # gene = gene[-1]
  gene
}



rankogram_check <- function(net){
  tryCatch(rankogram(net),
           error = function(e){
             "error"
           })
}


### Generate the full lower bound table from the pooled one
gen_full_lower <-function(N.treat, smry, group_i){
  
  full_lower = data.frame(matrix(rep(0, N.treat*N.treat), nrow = N.treat))
  rownames(full_lower) = colnames(full_lower) = 1:N.treat
  row_names = as.numeric(rownames(data.frame(smry$random$lower)))
  
  row_group_num = 0
  for (row in 1:N.treat) {
    for (subgroup in group_i) {
      if(row %in% subgroup){row_group_num = min(subgroup)}
    }
    # no pooled group
    if(row_group_num == 0){row_group_num = row}
    col_group_num = 0
    for (col in 1:N.treat) {
      for (subgroup in group_i) {
        if(col %in% subgroup){col_group_num = min(subgroup)}
      }
      if(col_group_num == 0){col_group_num = col}
      
      c1=which(row_names == row_group_num)
      c2=which(row_names == col_group_num)
      
      full_lower[row,col] = smry$random$lower[c1,c2]
      col_group_num = 0
    }
    row_group_num = 0
  }
  full_lower
}




