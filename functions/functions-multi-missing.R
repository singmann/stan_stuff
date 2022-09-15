
prep_data <- function(df) {
  df_use <- df[rowSums(is.na(df)) != ncol(df),]
  
  if (any(colMeans(is.na(df_use)) == 1)) {
    stop("one column has only NAs.", call. = FALSE)
  }
  
  colnum <- seq_len(ncol(df))
  
  subsets <- unlist(lapply(colnum, combn, 
                           x = colnum, simplify = FALSE), 
                    recursive = FALSE)
  subsets <- rev(subsets)
  df_out <- vector("list", length = length(subsets))
  #names(df_out) <- paste0("d_", vapply(subsets, paste, collapse = "_", FUN.VALUE = ""))
  names(df_out) <- paste0("d_", seq_len(length(subsets)))
  
  df_loop <- df_use
  
  for (i in seq_len(length(subsets))) {
    tmp <- select_pick_omit(df_loop, subsets[[i]])
    df_out[[i]] <- tmp$use
    df_loop <- tmp$df_loop
  }
  dims <- lapply(df_out, dim)
  dims <- do.call("rbind", dims)
  attr(dims, "dimnames") <- NULL
  
  which_var <- subsets
  names(which_var) <- paste0("whichvar_", seq_len(length(subsets)))
  
  c(
    nvars = ncol(df_use), nsets = length(subsets), 
    ddim = list(dims), which_var, df_out
  )
}


make_data_declaration <- function(dat_prep, file = "") {
  cat("data {\n", file = file)
  cat("  int nvars;\n", sep = "", file = file)
  cat("  int nsets;\n", sep = "", file = file)
  cat("  int ddim[nsets, 2]", ";\n", sep = "", file = file)
  
  ## variables for mapping columns
  for (i in grep("^whichvar_", names(dat_prep))) {
    # cat("  vector[", length(dat_prep[[i]]), "] ", names(dat_prep)[i], ";\n",
    #     sep = "", file = file)
    if (length(dat_prep[[i]]) > 1) {
      cat("  int ", names(dat_prep)[i], "[", length(dat_prep[[i]]), "];\n",
          sep = "", file = file) 
    } else {
      cat("  int ", names(dat_prep)[i], ";\n",
          sep = "", file = file)
    }
  }
  
  c <- 1
  for (i in grep("^d_", names(dat_prep))) {
    cat("  vector[ddim[", c, ", 2]] ",
        names(dat_prep)[i], "[ddim[", c, ", 1]];\n", sep = "", file = file)
    # cat("  array[ddim[", c, ", 2]] vector[ddim[", c, ", 1]] ", 
    #     names(dat_prep)[i], ";\n", sep = "", file = file)
    c <- c+1
  }
  cat("}\n", file = file)
}

make_parameters_block <- function(dat_prep, file = "") {
  cat("parameters {\n", file = file)
  cat("  corr_matrix[nvars] Omega;\n", sep = "", file = file)
  cat("  vector<lower=0>[nvars] sigma;\n", sep = "", file = file)
  cat("  vector[nvars] mu;\n", sep = "", file = file)
  cat("}\n", file = file)
}

make_model_block <- function(dat_prep, file = "", 
                             lkj_eta = 1, 
                             t_nu = 3, t_scale = 3.5) {
  cat("model {\n", file = file)
  cat("  Omega ~ lkj_corr(", lkj_eta,");\n", sep = "", file = file)
  cat("  sigma ~ student_t(", t_nu , ", 0, ", t_scale, ");\n\n", 
      sep = "", file = file)
  
  cat("  d_1 ~ multi_normal(mu, quad_form_diag(Omega, sigma));\n", 
      sep = "", file = file)
  #browser()
  subsets <- dat_prep[paste0("whichvar_", seq_len(dat_prep$nsets))][-1]
  length_subsets <- vapply(subsets, length, 1)
  
  vectors <- subsets[length_subsets == 1]
  cor_subsets <- subsets[length_subsets > 1]
  
  for (i in seq_along(cor_subsets)) {
    tmp_num <- strsplit(names(cor_subsets)[[i]], "_")[[1]][2]
    cat("  d_", tmp_num, " ~ multi_normal(", 
        "mu[", names(cor_subsets)[[i]], "], ", 
        "quad_form_diag(Omega_", tmp_num, 
        ", sigma[", names(cor_subsets)[[i]], "]));\n", 
        sep = "", file = file)
    
  }
  
  for (i in seq_along(vectors)) {
    tmp_num <- strsplit(names(vectors)[[i]], "_")[[1]][2]
    if (dat_prep$ddim[as.numeric(tmp_num), 1] == 0L) next
    cat("  d_", tmp_num, "[1] ~ normal(", 
        "mu[", names(vectors)[[i]], "], ", 
        "sqrt(sigma[", names(vectors)[[i]], "]));\n", 
        sep = "", file = file)
    # y[ns1] ~ normal(mu[1], sqrt(Sigma[1, 1]));
    
  }
  
  cat("}\n", file = file)
}

make_trans_parameters_block <- function(dat_prep, file = "") {
  cat("transformed parameters {\n", file = file)
  subsets <- dat_prep[paste0("whichvar_", seq_len(dat_prep$nsets))]
  
  temp_cor <- subsets[-1]
  temp_cor <- temp_cor[vapply(temp_cor, length, 0) != 1]
  
  for (i in seq_along(temp_cor)) {
    cat("  corr_matrix[", length(temp_cor[[i]]), "] Omega_", i+1,";\n", 
        sep = "", file = file)  
  }
  
  for (i in seq_along(temp_cor)) {
    tmp_subsets <- unlist(lapply(seq_len(length(temp_cor[[i]])), combn, 
                                 x = temp_cor[[i]], simplify = FALSE), 
                          recursive = FALSE)
    tmp_pairs <- tmp_subsets[vapply(tmp_subsets, length, 1) == 2]
    for (j in seq_len(length(temp_cor[[i]]))) {
      cat("  Omega_", i+1, "[", j, ",", j, "] = 1.0", ";\n", 
          sep = "", file = file)
    }
    for (j in seq_along(tmp_pairs)) {
      tar_x <- which(temp_cor[[i]] == tmp_pairs[[j]][[1]])
      tar_y <- which(temp_cor[[i]] == tmp_pairs[[j]][[2]])
      cat("  Omega_", i+1, "[", tar_x, ",", tar_y, "] = ",
          "Omega[", tmp_pairs[[j]][[1]], ",", tmp_pairs[[j]][[2]], "]", 
          ";\n", 
          sep = "", file = file)
      cat("  Omega_", i+1, "[", tar_y, ",", tar_x, "] = ",
          "Omega[", tmp_pairs[[j]][[2]], ",", tmp_pairs[[j]][[1]], "]", 
          ";\n", 
          sep = "", file = file)
    }
    
  }
  cat("}\n", file = file)
}


select_pick_omit <- function(data, select) {
  dat_out <- na.omit(data[select])
  dat_remain <- data[attr(dat_out, "na.action"),]
  attr(dat_out, "na.action") <- NULL
  list(use = dat_out, df_loop = dat_remain)
}

