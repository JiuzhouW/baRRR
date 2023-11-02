########### functions ###########

####################### helper functions #######################

# generate V or V-shaped zero matrix
V_generator = function(nc=10,nr=100,nonzero=1){
  if(nonzero==1){
    V = matrix(rnorm(nc*nr),ncol = nc)
  }else if(nonzero==0){
    V = matrix(rep(0,nc*nr),ncol = nc)
  }
  return(V)
}

# combine Y matrix for each module 
Y_combiner = function(Yi_list,mod){
  Y = NULL
  for (i in 1:length(mod)) {
    if(mod[i]==1){
      Y = cbind(Y,Yi_list[[i]])
    }else if(mod[i]==0){
      msize = dim(Yi_list[[i]])
      Y = cbind(Y,matrix(rep(0,msize[1]*msize[2]),ncol = msize[2]))
    }
  }
  return(Y)
}


# combine Y matrix without zero-like matrix
Y_combiner_nonzero = function(Yi_list,mod){
  Y = NULL
  for (i in 1:length(mod)) {
    if(mod[i]==1){
      Y = cbind(Y,Yi_list[[i]])
    }
  }
  return(Y)
}

############### data generation ################# 

data_generate = function(
  seed = 1,
  n_sample = c(200,100),
  n_feature = c(100,50),
  # number of features for Y, here assume Y is only multi-cohort
  # not multi-view
  p_y = 30,
  #R_b = 10, estimation parameter
  R_b_true = 1,
  #R_s = 5, estimation parameter
  R_s_true = 5,
  correlation = 0,
  orth_gen = FALSE,
  Binvolved = TRUE,
  Sinvolved = TRUE,
  # length of cohorts and views should match 
  modules_B_cohort = list(c(1,1),c(1,1),c(1,1),c(1,0),c(0,1),c(1,0),c(1,0),c(0,1),c(0,1)),
  modules_B_view = list(c(1,1),c(1,0),c(0,1),c(1,1),c(1,1),c(1,0),c(0,1),c(1,0),c(0,1)),
  modules_S_cohort = modules_B_cohort,
  modules_S_view = modules_B_view,
  # control the level of signal in each module
  sd_B = c(1,1,1,1,1,1,1,1,1),
  # correlates to generation of S
  # by default, make the structure of S same as that of B
  sd_S = sd_B
){
  # initialize 
  n_cohort = length(n_sample)
  n_sample_sum = sum(n_sample)
  n_view = length(n_feature)
  n_feature_sum = sum(n_feature)
  n_mod_B = length(modules_B_cohort)
  n_mod_S = length(modules_S_cohort)
  
  # not for data generation but 
  # prepare for optimization in S
  col_index = list()
  row_index = list()
  temp=0
  for(i in 1:n_cohort){
    col_index[[i]] = (temp+1):(temp+n_sample[i])
    temp = temp + n_sample[i]
  }
  temp=0
  for(i in 1:n_view){
    row_index[[i]] = (temp+1):(temp+n_feature[i])
    temp = temp + n_feature[i]
  }
  
  if(Sinvolved){
    # indicate sample index in each S module
    col_S_mod_index = list()
    for(i in 1:n_mod_S){
      temp = c()
      indicator = modules_S_cohort[[i]]
      for (j in 1:n_cohort) {
        if(indicator[j]==1){
          temp = c(temp,col_index[[j]])
        }
      }
      col_S_mod_index[[i]] = temp
    }
    
    # indicate feature index in each S module
    row_S_mod_index = list()
    for(i in 1:n_mod_S){
      temp = c()
      indicator = modules_S_view[[i]]
      for (j in 1:n_view) {
        if(indicator[j]==1){
          temp = c(temp,row_index[[j]])
        }
      }
      row_S_mod_index[[i]] = temp
    }
  }
  
  
  if(Binvolved){
    # indicate feature index in each B module
    col_B_mod_index = list()
    for(i in 1:n_mod_B){
      temp = c()
      indicator = modules_B_cohort[[i]]
      for (j in 1:n_cohort) {
        if(indicator[j]==1){
          temp = c(temp,col_index[[j]])
        }
      }
      col_B_mod_index[[i]] = temp
    }
    
    row_B_mod_index = list()
    for(i in 1:n_mod_B){
      temp = c()
      indicator = modules_B_view[[i]]
      for (j in 1:n_view) {
        if(indicator[j]==1){
          temp = c(temp,row_index[[j]])
        }
      }
      row_B_mod_index[[i]] = temp
    }
    
  }
  
  
  
  # Generate Y only if we include B, 
  # the data generating process should be different from optimization
  # covariate effect only has one view
  # maybe extend this to multiple views
  
  library(MASS)
  set.seed(seed)
  
  # individual Y_i, regressors in each group
  Y_org_list = list()
  sd_Y_org = NULL
  
  if(correlation==0){
    for (i in 1:n_cohort) {
      Y_org_list[[i]] = matrix(rnorm(n_sample[i]*p_y),ncol = n_sample[i])
      sd_Y_org = c(sd_Y_org,sd(Y_org_list[[i]]))
    }
  }else if(correlation>0){
    sigma = matrix(rep(correlation,p_y^2),ncol = p_y)+diag(1-correlation,p_y)
    for (i in 1:n_cohort) {
      Y_org_list[[i]] = t(mvrnorm(n_sample[i],rep(0,p_y),Sigma = sigma))
      sd_Y_org = c(sd_Y_org,sd(Y_org_list[[i]]))
    }
  }
  
  
  
  # individually orthogonalize each Yi
  if(orth_gen){
    for(i in 1:n_cohort){
      Y_org_list[[i]] = t(svd(Y_org_list[[i]])$v)
      # important step
      # rescale orthogonalized Y to have the same sd(Y)
      Y_org_list[[i]] = Y_org_list[[i]]/sd(Y_org_list[[i]])*sd_Y_org[i]
    }
  }
  
  #combined modulized Y to generate X_tot
  # another way to handle this is by generate X separately and combine together
  Y_gen_list = list()
  
  for (i in 1:n_mod_B) {
    Y_gen_list[[i]] = Y_combiner(Y_org_list,modules_B_cohort[[i]])
  }
  
  
  #initialize UB and VB, therefore B
  
  UB_list = list()
  VB_list = list()
  B_list = list()
  
  
  # scaling will be different if including S
  if(Sinvolved){
    for(i in 1:n_mod_B){
      
      UB_list[[i]] = matrix(rnorm(n_feature_sum*R_b_true),ncol = R_b_true)
      UB_list[[i]][-row_B_mod_index[[i]],] = 0
      VB_list[[i]] = matrix(rnorm(p_y*R_b_true),ncol = R_b_true)
      # calculate the current sd without blocks of zeros
      temp = UB_list[[i]]%*%t(VB_list[[i]])%*%Y_gen_list[[i]]
      temp_nonzero = temp[row_B_mod_index[[i]], col_B_mod_index[[i]]]
      B_list[[i]] = sd_B[i]*UB_list[[i]]%*%t(VB_list[[i]])/sd(temp_nonzero)
    }
  }else if(!Sinvolved){
    for(i in 1:n_mod_B){
      
      UB_list[[i]] = matrix(rnorm(n_feature_sum*R_b_true),ncol = R_b_true)
      UB_list[[i]][-row_B_mod_index[[i]],] = 0
      VB_list[[i]] = matrix(rnorm(p_y*R_b_true),ncol = R_b_true)
      # calculate the current sd without blocks of zeros
      temp = UB_list[[i]]%*%t(VB_list[[i]])
      temp_nonzero = temp[row_B_mod_index[[i]], ]
      B_list[[i]] = sd_B[i]*UB_list[[i]]%*%t(VB_list[[i]])/sd(temp_nonzero)
    }
  }
  
  # check the correlation after standardization
  # if orth_gen == 0, then it corresponds with the original since multiply a number to the Yi does not affect corr
  # orthonormalize Y will automatically make the correlation to be 0
  #cor(t(Y_org_list[[1]]))
  #cor(t(Y_list[[1]][,1:n_sample[1]]))
  
  # generate S if involved
  
  if(Sinvolved){
    #initialize U and V, therefore S
    
    U_list = list()
    V_list = list()
    S_list = list()
    
    for(i in 1:n_mod_S){
      
      #U_s may contain parts of zero matrices
      temp = NULL
      for (j in 1:length(modules_S_view[[i]])) {
        temp = rbind(temp,V_generator(nc=R_s_true,nr=n_feature[j],nonzero=modules_S_view[[i]][j]))
      }
      U_list[[i]] = temp
      
      #V_s may contain parts of zero matrices
      temp = NULL
      for (j in 1:length(modules_S_cohort[[i]])) {
        #print(modules_S_cohort[[i]][j])
        temp = rbind(temp,V_generator(nc=R_s_true,nr=n_sample[j],nonzero=modules_S_cohort[[i]][j]))
      }
      V_list[[i]] = temp
      
      temp = U_list[[i]]%*%t(V_list[[i]])
      temp_nonzero = temp[row_S_mod_index[[i]], col_S_mod_index[[i]]]
      S_list[[i]] = sd_S[i]*U_list[[i]]%*%t(V_list[[i]])/sd(temp_nonzero)
    }
  }else if(!Sinvolved){
    row_S_mod_index = "not involved"
    col_S_mod_index = "not involved"
    S_list = list()
    
    for(i in 1:n_mod_S){
      S_list[[i]] = matrix(0, ncol = n_sample_sum, nrow = n_feature_sum)
    }
  }
  
  
  X_tot = 0
  S_tot = 0
  BY_tot = 0
  for (i in 1:n_mod_B) {
    BY_tot = BY_tot + B_list[[i]]%*%Y_gen_list[[i]] 
  }
  for(i in 1:n_mod_S){
    S_tot = S_tot + S_list[[i]]
  }
  E_tot = matrix(rnorm(n_feature_sum*n_sample_sum),nrow = n_feature_sum)
  X_tot = S_tot + BY_tot + E_tot
  
  return(list(X_tot = X_tot,
              BY_tot = BY_tot,
              S_tot = S_tot,
              Y_org_list = Y_org_list,
              B_list = B_list,
              S_list = S_list,
              row_B_mod_index = row_B_mod_index,
              row_S_mod_index = row_S_mod_index,
              col_S_mod_index = col_S_mod_index,
              col_index = col_index,
              row_index = row_index
  ))
}


# Generate best penalties for auxiliary structures S 
# based on the random matrix theory mentioned in the manuscripts

lambda_S_gen = function(n_feature,modules_S_view,n_sample,modules_S_cohort){
  all_lambda = c()
  for (i in 1:length(modules_S_view)) {
    temp = sqrt(sum(n_sample[modules_S_cohort[[i]]==1])) + sqrt(sum(n_feature[modules_S_view[[i]]==1]))
    all_lambda = c(all_lambda,temp)
  }
  return(all_lambda)
}


# Generate best penalties for covariate effects B 
# based on the random matrix theory mentioned in the manuscripts

lambda_B_gen = function(n_feature, modules_B_view, p_y){
  all_lambda = c()
  for (i in 1:length(modules_B_view)) {
    temp = sqrt(sum(n_feature[modules_B_view[[i]]==1])) + sqrt(p_y)
    all_lambda = c(all_lambda,temp)
  }
  return(all_lambda)
}



######################parameter initialization#################

starters_generate = function(
  # set the maximum rank for estimators
  R_b = 5,
  R_s = 10,
  # set the initial values as 0
  Binvolved = TRUE,
  Sinvolved = TRUE,
  seed = 1,
  n_feature_sum = 150,
  n_sample_sum = 300,
  p_y = 10,
  n_mod_B = 1,
  n_mod_S = 1
){
  set.seed(seed)
  
  UB_s_list = list()
  VB_s_list = list()
  B_s_list = list()
  U_s_list = list()
  V_s_list = list()
  S_s_list = list()
  
  if(Binvolved){
    for(i in 1:n_mod_B){
      UB_s_list[[i]] = matrix(rnorm(n_feature_sum*R_b),ncol = R_b)
      # turn blocks of uninvolved views into zero
      # notice row_S_mod_index is made from S
      # UB_s_list[[i]][-row_B_mod_index[[i]],] = 0
      VB_s_list[[i]] = matrix(rnorm(p_y*R_b),ncol = R_b)
      B_s_list[[i]] = UB_s_list[[i]] %*% t(VB_s_list[[i]])
    }
  }else if(!Binvolved){
    for(i in 1:n_mod_B){
      #U_s is universal for one module
      UB_s_list[[i]] = matrix(0,nrow = n_feature_sum, ncol = R_b)
      VB_s_list[[i]] = matrix(0,nrow = p_y, ncol = R_b)
      B_s_list[[i]] = matrix(0, nrow = n_feature_sum, ncol = p_y)
    }
  }
  
  
  if(Sinvolved){
    for(i in 1:n_mod_S){
      U_s_list[[i]] = matrix(rnorm(n_feature_sum*R_s),ncol = R_s)
      V_s_list[[i]] = matrix(rnorm(n_sample_sum*R_s),ncol = R_s)
      S_s_list[[i]] = U_s_list[[i]]%*%t(V_s_list[[i]])
    }
  }else if(!Sinvolved){
    for(i in 1:n_mod_S){
      U_s_list[[i]] = matrix(rep(0,n_feature_sum*R_s),ncol = R_s)
      V_s_list[[i]] = matrix(rep(0,n_sample_sum*R_s),ncol = R_s)
      S_s_list[[i]] = U_s_list[[i]]%*%t(V_s_list[[i]])
    }
  }
  
  return(list(
    UB_s_list = UB_s_list,
    VB_s_list = VB_s_list,
    B_s_list = B_s_list,
    U_s_list = U_s_list,
    V_s_list = V_s_list,
    S_s_list = S_s_list
  ))
}


################## assign missingness ###########

missing_assigner = function(
  col_index = list(1:5,6:10,11:15),
  row_index = list(1:12,13:20),
  missing_prop = 0.2
){
  row_cut_off = 0
  for (i in 1:length(row_index)) {
    row_cut_off[i+1] = tail(row_index[[i]],1)
  }
  n_tot_row = tail(row_cut_off,1)
  
  col_cut_off = 0
  for (i in 1:length(col_index)) {
    col_cut_off[i+1] = tail(col_index[[i]],1)
  }
  n_tot_col = tail(col_cut_off,1)
  
  
  
  
  ########### partial column missing ##########
  # each list is selected index for the current cohort
  sel_col_list = list()
  # record the remaining choose of column index
  candidate_col_index_list = list()
  for (i in 1:length(col_index)) {
    candidate_col_index_list[[i]] = (col_cut_off[i]+1):col_cut_off[i+1]
  }
  
  
  for (i in 1:length(row_index)) {
    index_per_j = NULL
    for (j in 1:length(col_index)) {
      temp = candidate_col_index_list[[j]]
      if(length(temp)==1){
        # to avoid "sample" takes place from 1:x
        index_per_ij = temp
      }else{
        index_per_ij = sample(temp,
                              size = round(length(col_index[[j]])*missing_prop),
                              replace = FALSE)
      }
      index_per_j = c(index_per_j,index_per_ij)
      candidate_col_index_list[[j]] = temp[! temp %in% index_per_j]
    }
    sel_col_list[[i]] = index_per_j
  }
  miss_col_index = NULL
  for (i in 1:length(row_index)) {
    for (current_column in sel_col_list[[i]]) {
      temp = (current_column-1)*n_tot_row + row_cut_off[i]
      miss_col_index = append(miss_col_index,
                              (temp+1):(temp + length(row_index[[i]])))
    }
  }
  
  
  
  
  ##### total column missing #######
  sel_col_tot = NULL
  # the candidate_col_index is used 
  for (j in 1:length(col_index)) {
    if(length(candidate_col_index_list[[j]]) == 1){
      temp = candidate_col_index_list[[j]]
    }else{
      temp = sample(candidate_col_index_list[[j]],
                    size = round(length(col_index[[j]])*missing_prop),
                    replace = FALSE)
    }
    sel_col_tot = c(sel_col_tot,temp)
  }
  
  miss_col_tot_index = NULL
  for (current_column in sel_col_tot) {
    temp = (current_column-1)*n_tot_row 
    miss_col_tot_index = append(miss_col_tot_index,
                                (temp + 1):(temp + n_tot_row))
  }  
  
  
  
  
  ######### missing partial rows #########
  # each list is selected index for the current view
  sel_row_list = list()
  
  # record the remaining choose of column index
  candidate_row_index_list = list()
  for (i in 1:length(row_index)) {
    candidate_row_index_list[[i]] = (row_cut_off[i]+1):row_cut_off[i+1]
  }
  
  
  for (j in 1:length(col_index)) {
    index_per_i = NULL
    for (i in 1:length(row_index)) {
      temp = candidate_row_index_list[[i]]
      if(length(temp)==1){
        # to avoid "sample" takes place from 1:x
        index_per_ij = temp
      }else{
        index_per_ij = sample(temp,
                              size = round(length(row_index[[i]])*missing_prop),
                              replace = FALSE)
      }
      index_per_i = c(index_per_i,index_per_ij)
      candidate_row_index_list[[i]] = temp[! temp %in% index_per_i]
    }
    sel_row_list[[j]] = index_per_i
  }
  
  
  miss_row_index = NULL
  for (i in 1:length(col_index)) {
    current_rows = sel_row_list[[i]]
    for (j in col_index[[i]]) {
      miss_row_index = append(miss_row_index,
                              current_rows + (j-1)*n_tot_row)
    }
  }
  
  
  
  
  
  ########## missing entries ###########
  n_entries = n_tot_col*n_tot_row
  temp = 1:n_entries
  candidate_entry_index = temp[! temp %in% c(miss_col_index,
                                             miss_col_tot_index,
                                             miss_row_index)]
  miss_entry_index = sample(candidate_entry_index, 
                            size = round(missing_prop*n_entries),
                            replace = FALSE)
  
  
  
  # find common elements from multiple vectors
  intersect(miss_col_index,miss_col_tot_index)
  # check: should be integer(0)
  miss_part_col_row_index = intersect(miss_col_index, miss_row_index)
  miss_full_col_row_index = intersect(miss_col_tot_index, miss_row_index)
  miss_part_col_index = miss_col_index[! miss_col_index %in% miss_part_col_row_index]
  miss_full_col_index = miss_col_tot_index[! miss_col_tot_index %in% miss_full_col_row_index]
  miss_only_row_index = miss_row_index[! miss_row_index %in% c(miss_full_col_row_index,miss_part_col_row_index)]
  miss_tot_index = unique(c(miss_row_index,miss_col_index,miss_col_tot_index,miss_entry_index))
  
  
  return(list(
    miss_col_index = miss_col_index,
    miss_col_tot_index = miss_col_tot_index,
    miss_row_index = miss_row_index,
    miss_entry_index = miss_entry_index,
    # only full col, no row
    miss_full_col_index = miss_full_col_index,
    # both full col and row
    miss_full_col_row_index = miss_full_col_row_index,
    miss_part_col_index = miss_part_col_index,
    miss_part_col_row_index = miss_part_col_row_index,
    miss_only_row_index = miss_only_row_index,
    miss_tot_index = miss_tot_index
  ))
}


################### optimization ####################

# concatenate all matrices in a module, including zeros

mod_refine = function(X,mod_view,mod_cohort,n_sample,n_feature){
  ###########################################################################################################################
  ## @parameters: 
  ##    n_sample: numeric vector, number of samples for each cohort
  ##    X: matrix, full concatenated outcome matrix without zeros 
  ##    mod: numeric vector, current module information, e.g. c(1,0,1)
  ## ESTIMATION:
  ##    No estimation in the data generation function.
  ## @returns:
  ##    temp: matrix, concatenated outcome matrix for current module with zeros
  ###########################################################################################################################  
  
  temp = X
  
  sample_index_upper = 0
  sample_index_lower = 1
  for(i in 1:length(mod_cohort)){
    sample_index_upper = sample_index_upper + n_sample[i]
    
    feature_index_upper = 0
    feature_index_lower = 1
    
    for(j in 1:length(mod_view)){
      feature_index_upper = feature_index_upper + n_feature[j]
      
      #print(mod_cohort[i]*mod_view[j])
      
      if(mod_cohort[i]*mod_view[j]==0){
        temp[feature_index_lower:feature_index_upper,
             sample_index_lower:sample_index_upper] = 0
      }
      
      feature_index_lower = feature_index_lower + n_feature[j]
    }
    sample_index_lower = sample_index_lower + n_sample[i]
  }
  return(temp)
}


# replace all negative numbers with zeros in a vector

thres = function(vector){
  new_vec = NULL
  for (i in vector) {
    new_vec = c(new_vec,max(0,i))
  }
  return(new_vec)
}

# calculate current loss
# computationally slow

loss_calculator = function(X_tot,n_mod_B,n_mod_S,lambdaBs,lambdaSs,B_s_list,S_s_list,Y_list){
  X_res = X_tot 
  for(i in 1:n_mod_B){
    X_res = X_res - B_s_list[[i]]%*%Y_list[[i]] 
  }
  for(i in 1:n_mod_S){
    X_res = X_res - S_s_list[[i]]
  }
  
  loss = 0.5*sum(X_res^2) 
  for(i in 1:n_mod_B){
    loss = loss + lambdaBs[i]*sum(abs(svd(B_s_list[[i]])$d)) 
  }
  for(i in 1:n_mod_S){
    loss = loss + lambdaSs[i]*sum(abs(svd(S_s_list[[i]])$d))
  }
  return(loss)
}


impute_BS = function(
  X_tot_incom,
  Y_org_list,
  B_s_list,
  S_s_list,
  modules_B_cohort,
  modules_B_view,
  modules_S_cohort,
  modules_S_view,
  n_sample,
  n_feature,
  lambdaBs,
  lambdaSs,
  row_B_mod_index,
  row_S_mod_index,
  col_S_mod_index,
  Binvolved = TRUE,
  Sinvolved = TRUE,
  orth_sol = TRUE,
  loss_comp = FALSE,
  max_iter = 100,
  bound = 1e-5
){
  # get the matrix index of missing positions
  missing_pos = is.na(X_tot_incom)
  
  
  # initialize 
  
  n_sample_sum = sum(n_sample)
  n_feature_sum = sum(n_feature)
  n_mod_B = length(modules_B_cohort)
  n_mod_S = length(modules_S_cohort)
  p_y = dim(B_s_list[[1]])[2]
  # p_x denotes the total number of features in multiple views
  p_x = dim(B_s_list[[1]])[1]
  
  
  Y_list = list()
  
  for (i in 1:n_mod_B) {
    Y_list[[i]] = Y_combiner(Y_org_list,modules_B_cohort[[i]])
  }
  
  Y_SVD_list = list()
  
  for(i in 1:length(Y_list)){
    Y_SVD_list[[i]] = svd(Y_list[[i]])
  }
  
  # standardize Y and B if do not orthogonalize data
  # aim to make each row of Y has norm one
  if(!orth_sol){
    trans4B = list()
    for(i in 1:n_mod_B){
      reg_const = sqrt(rowSums(Y_list[[i]]^2))
      #B_list[[i]] = B_list[[i]]%*%diag(reg_const,ncol = length(reg_const))
      Y_list[[i]] = diag(1/reg_const,ncol = length(reg_const))%*%Y_list[[i]]
      # The estimated B should multiply the same scale
      trans4B[[i]] = diag(1/reg_const,ncol = length(reg_const))
    }
  }else if(orth_sol){
    trans4B = list()
    for(i in 1:length(Y_list)){
      Y_list[[i]] = t(Y_SVD_list[[i]]$v)
      trans4B[[i]] = diag(Y_SVD_list[[i]]$d^(-1))%*%t(Y_SVD_list[[i]]$u)
    }
  }
  
  # Use x_tot for update
  X_tot = X_tot_incom
  
  # initialize some values for missing
  # fill the NA with 0
  X_tot[missing_pos] = 0
  # # fill the NA with row means
  # for(i in 1:nrow(X_tot)){
  #   X_tot[i,is.na(X_tot[i,])] = mean(X_tot[i,], na.rm = TRUE)
  # }
  # current values that fill NA
  current_filling = X_tot[missing_pos]
  
  
  time = 0
  loss_list = NULL
  diff_list = NULL
  S_s_list_prev = S_s_list
  B_s_list_prev = B_s_list
  
  
  
  ##############soft-thresholding############
  
  for (iter in 1:max_iter) {
    time = time + 1
    criter = 0
    
    #construct residual
    X_res = X_tot
    
    if(Binvolved){
      for(i in 1:n_mod_B){
        X_res = X_res - B_s_list[[i]]%*%Y_list[[i]] 
      }
    }
    if(Sinvolved){
      for (i in 1:n_mod_S) {
        X_res = X_res - S_s_list[[i]]
      }
    }
    
    
    # update B
    if(Binvolved){
      for(i in 1:n_mod_B){
        X_res = X_res + B_s_list[[i]]%*%Y_list[[i]]
        #03/01/23:X_res does not need mod_refine since B does not have zero blocks
        #09/01/23:this is wrong, we might have zero blocks in B as well
        #refine X_res
        X_res_re = mod_refine(X_res, mod_cohort = modules_B_cohort[[i]],
                              mod_view = modules_B_view[[i]],
                              n_sample = n_sample, n_feature = n_feature)
        XresY_svd = svd(X_res_re%*%t(Y_list[[i]]))
        B_s_list[[i]] = XresY_svd$u%*%diag(thres(XresY_svd$d-lambdaBs[i]))%*%
          t(XresY_svd$v)
        
        criter = criter + sum((B_s_list[[i]]-B_s_list_prev[[i]])^2)
        
        X_res = X_res - B_s_list[[i]]%*%Y_list[[i]]
      }
    }
    
    
    
    # update S
    if(Sinvolved){
      for (i in 1:n_mod_S) {
        
        X_res = X_res + S_s_list[[i]]
        
        # only SVD for parameters of interests
        
        X_res2compose = X_res[row_S_mod_index[[i]],col_S_mod_index[[i]]]
        S_new_svd_2 = svd(X_res2compose)
        S_new2 = S_new_svd_2$u%*%diag(thres(S_new_svd_2$d - lambdaSs[i]))%*%t(S_new_svd_2$v)
        
        S_new2back = matrix(0, nrow = n_feature_sum, ncol = n_sample_sum)
        
        S_new2back[row_S_mod_index[[i]],col_S_mod_index[[i]]] = S_new2
        S_s_list[[i]] = S_new2back
        
        X_res = X_res - S_s_list[[i]]
        
        criter = criter + sum((S_s_list[[i]]-S_s_list_prev[[i]])^2)
      }
    }
    
    if(loss_comp){
      loss_list = c(loss_list,loss_calculator(X_tot,n_mod_B,
                                              n_mod_S,lambdaBs,lambdaSs,B_s_list,S_s_list,Y_list))
    }
    S_s_list_prev = S_s_list
    B_s_list_prev = B_s_list
    
    # generate new X_tot
    X_tot_new = 0
    if(Binvolved){
      for (i in 1:n_mod_B) {
        X_tot_new = X_tot_new + B_s_list[[i]]%*%Y_list[[i]] 
      }
    }
    if(Sinvolved){
      for(i in 1:n_mod_S){
        X_tot_new = X_tot_new + S_s_list[[i]]
      }
    }
    
    # replace the NA with updates
    X_tot[missing_pos] = X_tot_new[missing_pos]
    # current difference with previous NA filling values
    diff_filling = sum((current_filling-X_tot_new[missing_pos])^2)
    diff_list = c(diff_list, diff_filling)
    # update the current filling values 
    current_filling = X_tot_new[missing_pos]
    
    if(max(criter,diff_filling)<bound){
      break
    }
    
  }
  
  
  # transform back to access the true scale of B
  for(i in 1:length(Y_list)){
    B_s_list[[i]] = B_s_list[[i]]%*%trans4B[[i]]
  }
  
  
  return(list(
    X_tot_com = X_tot,
    B_s_list = B_s_list,
    B_s_list_prev = B_s_list_prev,
    S_s_list = S_s_list,
    S_s_list_prev = S_s_list_prev,
    Y_list = Y_list,
    time = time,
    criterion = criter,
    loss_list = loss_list
  ))
}


################### impute with UV decomposition ##########



impute_UV = function(
  X_tot_incom,
  Y_org_list,
  UB_s_list,
  VB_s_list,
  B_s_list,
  U_s_list,
  V_s_list,
  S_s_list,
  modules_B_cohort,
  modules_B_view,
  modules_S_cohort,
  modules_S_view,
  n_sample,
  n_feature,
  lambdaBs,
  lambdaSs,
  #row_B_mod_index,
  #row_S_mod_index,
  #col_S_mod_index,
  Binvolved = TRUE,
  Sinvolved = TRUE,
  orth_sol = TRUE,
  loss_comp = FALSE,
  max_iter = 100,
  bound = 1e-5
){
  # get the matrix index of missing positions
  missing_pos = is.na(X_tot_incom)
  
  # initialize 
  
  n_sample_sum = sum(n_sample)
  n_feature_sum = sum(n_feature)
  n_mod_B = length(modules_B_cohort)
  n_mod_S = length(modules_S_cohort)
  R_b = dim(UB_s_list[[1]])[2]
  p_y = dim(VB_s_list[[1]])[1]
  # p_x denotes the total number of features in multiple views
  p_x = dim(UB_s_list[[1]])[1]
  R_s = dim(U_s_list[[1]])[2]
  
  Y_list = list()
  
  for (i in 1:n_mod_B) {
    Y_list[[i]] = Y_combiner(Y_org_list,modules_B_cohort[[i]])
  }
  
  Y_SVD_list = list()
  
  for(i in 1:length(Y_list)){
    Y_SVD_list[[i]] = svd(Y_list[[i]])
  }
  
  # standardize Y and B if do not orthogonalize data
  # aim to make each row of Y has norm one
  if(!orth_sol){
    trans4B = list()
    for(i in 1:n_mod_B){
      reg_const = sqrt(rowSums(Y_list[[i]]^2))
      #B_list[[i]] = B_list[[i]]%*%diag(reg_const,ncol = length(reg_const))
      Y_list[[i]] = diag(1/reg_const,ncol = length(reg_const))%*%Y_list[[i]]
      # The estimated B should multiply the same scale
      trans4B[[i]] = diag(1/reg_const,ncol = length(reg_const))
    }
  }else if(orth_sol){
    trans4B = list()
    for(i in 1:length(Y_list)){
      Y_list[[i]] = t(Y_SVD_list[[i]]$v)
      trans4B[[i]] = diag(Y_SVD_list[[i]]$d^(-1))%*%t(Y_SVD_list[[i]]$u)
    }
  }
  
  X_tot = X_tot_incom
  
  # initialize some values for missing
  # fill the NA with row means
  X_tot[missing_pos] = 0
  # # fill the NA with row means
  # for(i in 1:nrow(X_tot)){
  #   X_tot[i,is.na(X_tot[i,])] <- mean(X_tot[i,], na.rm = TRUE)
  # }
  # current values that fill NA
  current_filling = X_tot[missing_pos]
  
  time = 0
  loss_list = NULL
  diff_list = 0
  S_s_list_prev = S_s_list
  B_s_list_prev = B_s_list
  
  
  for (iter in 1:max_iter) {
    
    time = time + 1
    criter = 0
    
    #construct residual
    X_res = X_tot
    if(Binvolved){
      for(i in 1:n_mod_B){
        X_res = X_res - B_s_list[[i]]%*%Y_list[[i]] 
      }
    }
    if(Sinvolved){
      for (i in 1:n_mod_S) {
        X_res = X_res - S_s_list[[i]]
      }
    }
    
    
    # update B
    if(Binvolved){
      for(i in 1:n_mod_B){
        X_res = X_res + B_s_list[[i]]%*%Y_list[[i]]
        
        X_res_re = mod_refine(X_res, mod_cohort = modules_B_cohort[[i]],
                              mod_view = modules_B_view[[i]],
                              n_sample = n_sample, n_feature = n_feature)
        
        UB_s_list[[i]] = X_res_re%*%t(Y_list[[i]])%*%VB_s_list[[i]]%*%
          solve(t(VB_s_list[[i]])%*%
                  crossprod(t(Y_list[[i]]))%*%
                  VB_s_list[[i]]+lambdaBs[i]*diag(R_b))
        
        
        inv = solve(kronecker(crossprod(UB_s_list[[i]]),crossprod(t(Y_list[[i]]))) + lambdaBs[i]*diag(R_b*p_y))
        V_s_vec = inv%*%as.vector(Y_list[[i]]%*%t(X_res_re)%*%UB_s_list[[i]])
        VB_s_list[[i]] = matrix(V_s_vec, ncol = R_b, byrow = FALSE)
        
        
        B_s_list[[i]] = UB_s_list[[i]]%*%t(VB_s_list[[i]])
        
        criter = criter + sum((B_s_list[[i]]-B_s_list_prev[[i]])^2)
        
        X_res = X_res - B_s_list[[i]]%*%Y_list[[i]]
      }
    }
    
    
    # update S
    if(Sinvolved){
      for (i in 1:n_mod_S) {
        
        X_res = X_res + S_s_list[[i]]
        X_res_re = mod_refine(X_res, mod_cohort = modules_S_cohort[[i]],
                              mod_view = modules_S_view[[i]],
                              n_sample = n_sample, n_feature = n_feature)
        
        # refine X_s to satisfy the module structure, ie. some parts to be 0
        U_s_list[[i]] = X_res_re%*%V_s_list[[i]]%*%solve(t(V_s_list[[i]])%*%V_s_list[[i]] + lambdaSs[i]*diag(R_s))
        
        V_s_list[[i]] = t(X_res_re)%*%U_s_list[[i]]%*%solve(t(U_s_list[[i]])%*%U_s_list[[i]] + lambdaSs[i]*diag(R_s))
        S_s_list[[i]] = U_s_list[[i]]%*%t(V_s_list[[i]])
        
        X_res = X_res - S_s_list[[i]]
        
        criter = criter + sum((S_s_list[[i]]-S_s_list_prev[[i]])^2)
      }
    }
    
    if(loss_comp){
      loss_list = c(loss_list,loss_calculator(X_tot,n_mod_B,
                                              n_mod_S,lambdaBs,lambdaSs,B_s_list,S_s_list,Y_list))
    }
    S_s_list_prev = S_s_list
    B_s_list_prev = B_s_list
    
    # generate new X_tot
    X_tot_new = 0
    if(Binvolved){
      for (i in 1:n_mod_B) {
        X_tot_new = X_tot_new + B_s_list[[i]]%*%Y_list[[i]] 
      }
    }
    if(Sinvolved){
      for(i in 1:n_mod_S){
        X_tot_new = X_tot_new + S_s_list[[i]]
      }
    }
    
    # replace the NA with updates
    X_tot[missing_pos] = X_tot_new[missing_pos]
    # current difference with previous NA filling values
    diff_filling = sum((current_filling-X_tot_new[missing_pos])^2)
    diff_list = c(diff_list, diff_filling)
    # update the current filling values 
    current_filling = X_tot_new[missing_pos]
    
    if(max(criter,diff_filling)<bound){
      break
    }
    
  }
  
  # transform back to access the true scale of B
  for(i in 1:length(Y_list)){
    B_s_list[[i]] = B_s_list[[i]]%*%trans4B[[i]]
  }
  
  return(list(
    X_tot_com = X_tot,
    B_s_list = B_s_list,
    B_s_list_prev = B_s_list_prev,
    S_s_list = S_s_list,
    S_s_list_prev = S_s_list_prev,
    UB_s_list = UB_s_list,
    VB_s_list = VB_s_list,
    U_s_list = U_s_list,
    V_s_list = V_s_list,
    Y_list = Y_list,
    time = time,
    criterion = criter,
    loss_list = loss_list
  ))
}

################### Opt with BS soft-thresholding ##########



opt_BS = function(
  X_tot,
  Y_org_list,
  B_s_list,
  S_s_list,
  modules_B_cohort,
  modules_B_view,
  modules_S_cohort,
  modules_S_view,
  n_sample,
  n_feature,
  lambdaBs,
  lambdaSs,
  row_B_mod_index,
  row_S_mod_index,
  col_S_mod_index,
  Binvolved = TRUE,
  Sinvolved = TRUE,
  orth_sol = TRUE,
  loss_comp = FALSE,
  max_iter = 100,
  bound = 1e-5
){
  # initialize 
  
  n_sample_sum = sum(n_sample)
  n_feature_sum = sum(n_feature)
  n_mod_B = length(modules_B_cohort)
  n_mod_S = length(modules_S_cohort)
  p_y = dim(B_s_list[[1]])[2]
  # p_x denotes the total number of features in multiple views
  p_x = dim(B_s_list[[1]])[1]
  
  
  Y_list = list()
  
  for (i in 1:n_mod_B) {
    Y_list[[i]] = Y_combiner(Y_org_list,modules_B_cohort[[i]])
  }
  
  Y_SVD_list = list()
  
  for(i in 1:length(Y_list)){
    Y_SVD_list[[i]] = svd(Y_list[[i]])
  }
  
  # standardize Y and B if do not orthogonalize data
  # aim to make each row of Y has norm one
  if(!orth_sol){
    trans4B = list()
    for(i in 1:n_mod_B){
      reg_const = sqrt(rowSums(Y_list[[i]]^2))
      #B_list[[i]] = B_list[[i]]%*%diag(reg_const,ncol = length(reg_const))
      Y_list[[i]] = diag(1/reg_const,ncol = length(reg_const))%*%Y_list[[i]]
      # The estimated B should multiply the same scale
      trans4B[[i]] = diag(1/reg_const,ncol = length(reg_const))
    }
  }else if(orth_sol){
    trans4B = list()
    for(i in 1:length(Y_list)){
      Y_list[[i]] = t(Y_SVD_list[[i]]$v)
      trans4B[[i]] = diag(Y_SVD_list[[i]]$d^(-1))%*%t(Y_SVD_list[[i]]$u)
    }
  }
  
  
  
  time = 0
  loss_list = NULL
  S_s_list_prev = S_s_list
  B_s_list_prev = B_s_list
  
  
  
  ##############soft-thresholding############
  
  for (iter in 1:max_iter) {
    time = time + 1
    criter = 0
    
    #construct residual
    X_res = X_tot
    
    if(Binvolved){
      for(i in 1:n_mod_B){
        X_res = X_res - B_s_list[[i]]%*%Y_list[[i]] 
      }
    }
    if(Sinvolved){
      for (i in 1:n_mod_S) {
        X_res = X_res - S_s_list[[i]]
      }
    }

    
    # update B
    if(Binvolved){
      for(i in 1:n_mod_B){
        X_res = X_res + B_s_list[[i]]%*%Y_list[[i]]
        #03/01/23:X_res does not need mod_refine since B does not have zero blocks
        #09/01/23:this is wrong, we might have zero blocks in B as well
        #refine X_res
        X_res_re = mod_refine(X_res, mod_cohort = modules_B_cohort[[i]],
                              mod_view = modules_B_view[[i]],
                              n_sample = n_sample, n_feature = n_feature)
        XresY_svd = svd(X_res_re%*%t(Y_list[[i]]))
        B_s_list[[i]] = XresY_svd$u%*%diag(thres(XresY_svd$d-lambdaBs[i]))%*%
          t(XresY_svd$v)
        
        criter = criter + sum((B_s_list[[i]]-B_s_list_prev[[i]])^2)
        
        X_res = X_res - B_s_list[[i]]%*%Y_list[[i]]
      }
    }
    
    
    
    # update S
    if(Sinvolved){
      for (i in 1:n_mod_S) {
        
        X_res = X_res + S_s_list[[i]]
        
        # only SVD for parameters of interests
        
        X_res2compose = X_res[row_S_mod_index[[i]],col_S_mod_index[[i]]]
        S_new_svd_2 = svd(X_res2compose)
        S_new2 = S_new_svd_2$u%*%diag(thres(S_new_svd_2$d - lambdaSs[i]))%*%t(S_new_svd_2$v)
        
        S_new2back = matrix(0, nrow = n_feature_sum, ncol = n_sample_sum)
        
        S_new2back[row_S_mod_index[[i]],col_S_mod_index[[i]]] = S_new2
        S_s_list[[i]] = S_new2back
        
        X_res = X_res - S_s_list[[i]]
        
        criter = criter + sum((S_s_list[[i]]-S_s_list_prev[[i]])^2)
      }
    }
    
    if(loss_comp){
      loss_list = c(loss_list,loss_calculator(X_tot,n_mod_B,
                                              n_mod_S,lambdaBs,lambdaSs,B_s_list,S_s_list,Y_list))
    }
    S_s_list_prev = S_s_list
    B_s_list_prev = B_s_list
    
    if(criter<bound){
      break
    }
    
  }
  
  
  # transform back to access the true scale of B
  for(i in 1:length(Y_list)){
    B_s_list[[i]] = B_s_list[[i]]%*%trans4B[[i]]
  }
  
  
  return(list(
    B_s_list = B_s_list,
    B_s_list_prev = B_s_list_prev,
    S_s_list = S_s_list,
    S_s_list_prev = S_s_list_prev,
    Y_list = Y_list,
    time = time,
    criterion = criter,
    loss_list = loss_list
  ))
}


################### Opt with UV decomposition ##########



opt_UV = function(
  X_tot,
  Y_org_list,
  UB_s_list,
  VB_s_list,
  B_s_list,
  U_s_list,
  V_s_list,
  S_s_list,
  modules_B_cohort,
  modules_B_view,
  modules_S_cohort,
  modules_S_view,
  n_sample,
  n_feature,
  lambdaBs,
  lambdaSs,
  #row_B_mod_index,
  #row_S_mod_index,
  #col_S_mod_index,
  Binvolved = TRUE,
  Sinvolved = TRUE,
  orth_sol = TRUE,
  loss_comp = FALSE,
  max_iter = 100,
  bound = 1e-5
){
  # initialize 
  
  n_sample_sum = sum(n_sample)
  n_feature_sum = sum(n_feature)
  n_mod_B = length(modules_B_cohort)
  n_mod_S = length(modules_S_cohort)
  R_b = dim(UB_s_list[[1]])[2]
  p_y = dim(VB_s_list[[1]])[1]
  # p_x denotes the total number of features in multiple views
  p_x = dim(UB_s_list[[1]])[1]
  R_s = dim(U_s_list[[1]])[2]
  
  Y_list = list()
  
  for (i in 1:n_mod_B) {
    Y_list[[i]] = Y_combiner(Y_org_list,modules_B_cohort[[i]])
  }
  
  Y_SVD_list = list()
  
  for(i in 1:length(Y_list)){
    Y_SVD_list[[i]] = svd(Y_list[[i]])
  }
  
  # standardize Y and B if do not orthogonalize data
  # aim to make each row of Y has norm one
  if(!orth_sol){
    trans4B = list()
    for(i in 1:n_mod_B){
      reg_const = sqrt(rowSums(Y_list[[i]]^2))
      #B_list[[i]] = B_list[[i]]%*%diag(reg_const,ncol = length(reg_const))
      Y_list[[i]] = diag(1/reg_const,ncol = length(reg_const))%*%Y_list[[i]]
      # The estimated B should multiply the same scale
      trans4B[[i]] = diag(1/reg_const,ncol = length(reg_const))
    }
  }else if(orth_sol){
    trans4B = list()
    for(i in 1:length(Y_list)){
      Y_list[[i]] = t(Y_SVD_list[[i]]$v)
      trans4B[[i]] = diag(Y_SVD_list[[i]]$d^(-1))%*%t(Y_SVD_list[[i]]$u)
    }
  }
  
  
  
  time = 0
  loss_list = NULL
  S_s_list_prev = S_s_list
  B_s_list_prev = B_s_list
  
  
  for (iter in 1:max_iter) {
    
    time = time + 1
    criter = 0
    
    #construct residual
    X_res = X_tot
    
    if(Binvolved){
      for(i in 1:n_mod_B){
        X_res = X_res - B_s_list[[i]]%*%Y_list[[i]] 
      }
    }
    if(Sinvolved){
      for (i in 1:n_mod_S) {
        X_res = X_res - S_s_list[[i]]
      }
    }

    
    # update B
    if(Binvolved){
      for(i in 1:n_mod_B){
        X_res = X_res + B_s_list[[i]]%*%Y_list[[i]]
        
        X_res_re = mod_refine(X_res, mod_cohort = modules_B_cohort[[i]],
                              mod_view = modules_B_view[[i]],
                              n_sample = n_sample, n_feature = n_feature)
        
        UB_s_list[[i]] = X_res_re%*%t(Y_list[[i]])%*%VB_s_list[[i]]%*%
          solve(t(VB_s_list[[i]])%*%
                  crossprod(t(Y_list[[i]]))%*%
                  VB_s_list[[i]]+lambdaBs[i]*diag(R_b))
        
        
        inv = solve(kronecker(crossprod(UB_s_list[[i]]),crossprod(t(Y_list[[i]]))) + lambdaBs[i]*diag(R_b*p_y))
        V_s_vec = inv%*%as.vector(Y_list[[i]]%*%t(X_res_re)%*%UB_s_list[[i]])
        VB_s_list[[i]] = matrix(V_s_vec, ncol = R_b, byrow = FALSE)
        
        
        B_s_list[[i]] = UB_s_list[[i]]%*%t(VB_s_list[[i]])
        
        criter = criter + sum((B_s_list[[i]]-B_s_list_prev[[i]])^2)
        
        X_res = X_res - B_s_list[[i]]%*%Y_list[[i]]
      }
    }
    
    
    # update S
    if(Sinvolved){
      for (i in 1:n_mod_S) {
        
        X_res = X_res + S_s_list[[i]]
        X_res_re = mod_refine(X_res, mod_cohort = modules_S_cohort[[i]],
                              mod_view = modules_S_view[[i]],
                              n_sample = n_sample, n_feature = n_feature)
        
        # refine X_s to satisfy the module structure, ie. some parts to be 0
        U_s_list[[i]] = X_res_re%*%V_s_list[[i]]%*%solve(t(V_s_list[[i]])%*%V_s_list[[i]] + lambdaSs[i]*diag(R_s))
        
        V_s_list[[i]] = t(X_res_re)%*%U_s_list[[i]]%*%solve(t(U_s_list[[i]])%*%U_s_list[[i]] + lambdaSs[i]*diag(R_s))
        S_s_list[[i]] = U_s_list[[i]]%*%t(V_s_list[[i]])
        
        X_res = X_res - S_s_list[[i]]
        
        criter = criter + sum((S_s_list[[i]]-S_s_list_prev[[i]])^2)
      }
    }
    
    if(loss_comp){
      loss_list = c(loss_list,loss_calculator(X_tot,n_mod_B,
                                              n_mod_S,lambdaBs,lambdaSs,B_s_list,S_s_list,Y_list))
    }
    S_s_list_prev = S_s_list
    B_s_list_prev = B_s_list
    
    if(criter<bound){
      break
    }
    
  }
  
  # transform back to access the true scale of B
  for(i in 1:length(Y_list)){
    B_s_list[[i]] = B_s_list[[i]]%*%trans4B[[i]]
  }
  
  return(list(
    B_s_list = B_s_list,
    B_s_list_prev = B_s_list_prev,
    S_s_list = S_s_list,
    S_s_list_prev = S_s_list_prev,
    UB_s_list = UB_s_list,
    VB_s_list = VB_s_list,
    U_s_list = U_s_list,
    V_s_list = V_s_list,
    Y_list = Y_list,
    time = time,
    criterion = criter,
    loss_list = loss_list
  ))
}

