####Functions to predict recovery using models#####
####with and without blood based biomarkers#####

######################################################
#
# Fit models with and without chosen biomarker
#
######################################################

fit_mymodels_bio <- function(df, mymodels){
  
  models_imp <- vector(mode = "list", length = length(mymodels))
  names(models_imp) <- mymodels
  
  if("bio_only" %in% mymodels){
    models_imp[["bio_only"]] <-
      glm(Recovery ~ 
            bio, 
          data = df,
          family = "binomial")
  }
  
  if("plain_upfrontEM" %in% mymodels){
    
    models_imp[["plain_upfrontEM"]] <- 
      glm(Recovery ~ 
            Sex +
            Education*Age +
            HxMH +
            InjViolenceVictimAlcohol +
            GCSScoreBaselineDerived +
            PTAover1h, 
          data = df,
          family = "binomial")
    
    models_imp[["bio_upfrontEM"]] <- 
      glm(Recovery ~ 
            Sex +
            Education*Age +
            HxMH +
            InjViolenceVictimAlcohol +
            GCSScoreBaselineDerived +
            PTAover1h +
            bio, 
          data = df,
          family = "binomial")
  }
  
  
  if("plain_headsmart" %in% mymodels){
    
    models_imp[["plain_headsmart"]] <- 
      glm(Recovery ~ 
            Age +
            HxMH +
            acute_RPQHeadaches +
            acute_RPQPoorConcentration +
            acute_RPQLightSensitivity, 
          data = df,
          family = "binomial")
    
    models_imp[["bio_headsmart"]] <- 
      glm(Recovery ~ 
            Age +
            HxMH +
            acute_RPQHeadaches +
            acute_RPQPoorConcentration +
            acute_RPQLightSensitivity +
            bio, 
          data = df,
          family = "binomial")
  }
  
  if("plain_centerEM" %in% mymodels){
    
    models_imp[["plain_centerEM"]] <- 
      glm(Recovery ~ 
            Age +
            ISS +
            GCSScoreBaselineDerived +
            Sex +
            HxMH +
            ASA +
            Cause +
            acute_RPQTotalScore, 
          data = df,
          family = "binomial")
    
    models_imp[["bio_centerEM"]] <- 
      glm(Recovery ~ 
            Age +
            ISS +
            GCSScoreBaselineDerived +
            Sex +
            HxMH +
            ASA +
            Cause +
            acute_RPQTotalScore +
            bio, 
          data = df,
          family = "binomial")
  }
  
  return(models_imp)
  
}

####################################################
#
# save myres as a pretty table in word
#
####################################################

save_result_as_pretty_table_bio <- function(myres, mytitle) {
  #paste point estimates with corresponding 95% confidence interval into the same cell
  temp <- paste0(round_fun(myres[["pooled_estimate"]]), 
                 "\n(", 
                 round_fun(myres[["pooled_ll"]]), 
                 "-", 
                 round_fun(myres[["pooled_ul"]]), 
                 ")"
  )
  
  dim(temp) <- c(nrow(myres[["pooled_estimate"]]), ncol(myres[["pooled_estimate"]]))
  colnames(temp) <- colnames(myres[["pooled_estimate"]])
  rownames(temp) <- rownames(myres[["pooled_estimate"]])
  
  #format model name display
  temp <- 
    temp %>%
    as.data.frame() 
  
  temp$model <-
    row.names(temp)
  
  temp <- 
    temp %>%
    select(auc,delta_auc, nag, delta_nag, brier, sens, spec, ppv, npv, hl, lrt, model)
  
  #adjust p-value (note only one p per two models, so NAs present
  temp[c(1:4), "lrt"] <- p.adjust(temp[c(1:4), "lrt"], method = "fdr")
  
  
  #reorder rows
  myorder <- c("bio_only",
               "plain_upfrontEM",
               "bio_upfrontEM",
               "plain_headsmart",
               "bio_headsmart",
               "plain_centerEM",
               "bio_centerEM")
  temp$model <-
    factor(temp$model,
           levels = myorder)
  temp <- 
    temp %>%
    arrange(model)
  
  #remove CI for H-L p-value as non-sensical
  temp <-
    temp %>%
    separate(col = "hl", 
             into = "hl",
             sep = "\n") %>%
    separate(col = "lrt", 
             into = "lrt",
             sep = "\n")
  
  #for p-values (HL and lrt) display as "<0.001" if appropriate
  temp$hl  <- ifelse(temp$hl < 0.001, "<0.001", temp$hl)
  temp$lrt <- ifelse(temp$lrt < 0.001, "<0.001", temp$lrt)
  temp$hl  <- ifelse(temp$hl > 0.999, ">0.999", temp$hl)
  
  temp$lrt <-
    ifelse(temp$lrt == "NA", "", temp$lrt)
  temp$delta_auc <-
    ifelse(temp$delta_auc == "NA\n(NA-NA)", "", temp$delta_auc)
  temp$delta_nag <-
    ifelse(temp$delta_nag == "NA\n(NA-NA)", "", temp$delta_nag)
  
  row.names(temp) <- temp$model
  temp <- temp %>% select(-model)
  temp <- t(temp) %>% as.data.frame()
  
  temp$Metric <- row.names(temp)
  temp <- 
    temp %>%
    select(Metric, everything())
  
  #select only those metrics I want to display
  temp <- 
    temp %>%
    filter(Metric %in% c("auc", "nag", "sens", "spec", "ppv", "npv", "lrt"))
  
  #change order and label of metrics displayed
  temp$Metric <-
    factor(temp$Metric,
           levels = c("auc", "delta_auc", "nag", "delta_nag", "brier", "sens", "spec", "ppv", "npv", "hl", "lrt"),
           labels = c("Area under\nthe curve", 
                      "Incremental AUC (logit scale)",
                      "Variation\nexplained (%)", 
                      "Incremental variation\nexplained (%)",
                      "Brier score", 
                      "Sensitivity", 
                      "Specificity", 
                      "Positive predictive value",
                      "Negative predictive value",
                      "Hosmer-\nLemeshow", 
                      "Likelihood ratio\ntest p-value"))
  
  gose5 <-
    temp %>% 
    select(Metric, all_of(myorder))
  
  
  #create flextable object for GOSE < 5 results
  ft5 <-
    gose5 %>%
    regulartable() %>%
    set_header_labels("plain_upfrontEM" = "w/o bio",
                      "plain_headsmart" = "w/o bio",
                      "plain_centerEM" = "w/o bio",
                      "bio_upfrontEM" = "with bio",
                      "bio_headsmart" = "with bio",
                      "bio_centerEM" = "with bio",
                      "Metric" = " ",
                      "bio_only" = " ") %>%
    add_header_row(values = c("Metric","bio only", "UPFRONT-ED", "HeadSMART", "CENTER-ED"),
                   colwidths = c(1, 1, 2, 2, 2))%>%
    vline(j = c(1, 2, 4, 6, 8))%>%
    vline_left() %>%
    vline_right() %>%
    align(align = "right", part = "all")%>%
    align(align = "center", part = "header", i = c(1:2)) %>%
    align(align = "left", j = c(1), part = "all") %>% 
    padding(padding = 0, part = "all", padding.right = 3, padding.left = 3, padding.top = 3, padding.bottom = 3) %>%
    fontsize(size = 10, part = "all") %>%
    fix_border_issues(part = "all") 
  
  save_as_docx(ft5, path = paste0("Figures_bio/performance_", mytitle, ".docx"))
  
  
}

#############################################################
# display model pooled coefficients for
# all models with and without DTI
# fitted to the original (imputed but not boottsrapped) data
# and save as pretty table
#############################################################
save_model_coefficients_bio <- function(orig_coef_list, mybio){
  
  temp <- lapply(orig_coef_list, df_fun)
  
  for (i in 1:length(temp)){
    temp[[i]]$Model <- names(temp)[[i]]
  }
  
  temp <- bind_rows(temp)
  #temp <-
  #  temp %>%
  #  rename("pooled_estimate" = Estimate, 
  #         "pooled_se" = SE_within)
  temp$LL <- temp$pooled_estimate - I(1.96*temp$pooled_se)
  temp$UL <- temp$pooled_estimate + I(1.96*temp$pooled_se)
  temp$Sig <- ifelse(temp$LL <0 & temp$UL >0, "ns", "signif")
  
  temp$pooled_estimate <- format(as.numeric(round(temp$pooled_estimate, 3)), nsmall = 3)
  temp$pooled_se <- format(as.numeric(round(temp$pooled_se, 3)), nsmall = 3)
  temp$UL <- format(as.numeric(round(temp$UL, 3)), nsmall = 3)
  temp$LL <- format(as.numeric(round(temp$LL, 3)), nsmall = 3)
  
  temp <- 
    temp %>%
    select(Model, Variable, pooled_estimate, LL, UL, Sig)
  temp$Model<-
    factor(temp$Model,
           levels = c("bio_only",
                      "plain_headsmart",
                      "plain_upfrontEM",
                      "plain_centerEM",
                      "bio_headsmart",
                      "bio_upfrontEM",
                      "bio_centerEM"),
           labels = c("Bio only",
                      "HeadSMART",
                      "UPFRONT-ED",
                      "CENTER-ED",
                      "HeadSMART with bio",
                      "UPFRONT-ED with bio",
                      "CENTER-ED with bio"))
  
  
  #turn table into flextable object
  ft_model_coefs <- 
    temp %>%
    regulartable() %>%
    set_header_labels("pooled_estimate" = "Pooled beta") %>%
    vline_left() %>%
    vline_right() %>%
    vline(j = 1) %>%
    merge_v(j = 1) %>%
    fix_border_issues(part = "all") %>%
    bold(part = "header") %>%
    align(j = c(3,4,5,6), align = "right") %>%
    padding(padding = 1, padding.top = 3) %>%
    fontsize(size = 10, part = "all") %>%
    valign(valign = "bottom", j = c(2:6))
  
  
  #add horizontal line after merged cell
  row_loc <- rle(cumsum(ft_model_coefs$body$spans$columns[,1] ))$values
  ft_model_coefs <- 
    ft_model_coefs %>% 
    border(border.bottom = fp_border_default(),
           i=row_loc, 
           j = 1:6, 
           part="body") 
  ft_model_coefs <- 
    ft_model_coefs %>% 
    border(border.bottom = fp_border_default(),
           i = ft_model_coefs$body$spans$columns[,1] > 1, 
           j = 1, 
           part="body") %>% 
    border(border.bottom = fp_border_default(), 
           border.top = fp_border_default(),
           part = "header")
  
  ft_model_coefs
  
  
  #save to word document
  save_as_docx(ft_model_coefs, path = paste0("Figures_bio/modelcoefs_", mybio, ".docx"))
  
}


#################################################
#
# Mega function for a chosen biomarker
# reference models are fitted with and without
# the biomarker and optimism-corrected
# performance measures are calculated
#
# all intermediate objects are saved as RDS
# a final results table (but no graphs)
# are saved
# graphs could later be created from RDS objects
##################################################

biomarker_megafunction <- function(orig_pat_list, imp_boot_data_list, mybio){
  ####################################
  # Assign the biomarker of interest
  # to the variable bio
  # in each bootstrap sample
  ####################################
  for (i in 1:m) {
    
    orig_pat_list[[i]]$bio <- orig_pat_list[[i]][, mybio]
    for (b in 1:B){
      imp_boot_data_list[[i]][[b]]$bio <- imp_boot_data_list[[i]][[b]][, mybio]
    }
  }
  
  saveRDS(imp_boot_data_list, paste0("RDS_bio/imp_boot_data_list_", mybio, ".rds"))
  saveRDS(orig_pat_list, paste0("RDS_bio/orig_pat_list_", mybio, ".rds"))
  
  ############################################################################
  # For original data
  # fit prognostic models with and without DTI
  # derive unadjusted performance metrics 
  #(optimism-adjustment will be done later)
  ############################################################################
  
  #prepare model names
  mymodels <<- 
    c( "bio_only",            
       "plain_upfrontEM",    
       "plain_headsmart",
       "plain_centerEM",
       "bio_upfrontEM",      
       "bio_headsmart",
       "bio_centerEM")      
  
  #fit models to each of the 10 imputed datasets
  orig_models_by_imp <- vector(mode = "list", length = m)
  
  for (i in 1:m) {
    orig_models_by_imp[[i]] <- fit_mymodels_bio(orig_pat_list[[i]], mymodels)
  }
  
  
  #extract coefficients for all models across the imputed datasets
  orig_coef_list              <- vector(mode = "list", length = length(mymodels))
  names(orig_coef_list)       <- mymodels
  orig_models_by_model        <- vector(mode = "list", length = length(mymodels))
  names(orig_models_by_model) <- mymodels
  
  for (i in 1:length(mymodels)) {
    #reorganise nested list by model type
    orig_models_by_model[[i]] <- lapply(orig_models_by_imp, `[[`, i)
    #pool coefficients using Rubin's rules to report in paper
    orig_coef_list[[i]]       <- pool_my_coefs(orig_models_by_model[[i]]) 
    
  }
  
  # test model performance with and without bio
  # in each imputed dataset of the original data
  orig_perform_list <- vector(mode = "list", length = m)
  
  for (i in 1:m){
    orig_perform_list[[i]] <- test_mymodels(orig_pat_list[[i]], orig_models_by_imp[[i]])
  }
  
  #I am not pooling performance metrics now, only after optimism-correction
  
  saveRDS(orig_coef_list, paste0("RDS_bio/orig_coef_list_", mybio, ".rds"))
  
  
  #########################################################################
  # For each bootstrap sample within each imputed dataset
  #
  # a) fit models with and without DTI to bootstrap sample
  # b) test models from b) on bootstrap sample
  # c) test models from b) on original data
  # d) calculate optimism for this bootstrap sample
  #########################################################################
  
  #For each bootstrap sample within each imputed dataset I want to store
  # 1. mean performance of boot_models on boot_data (for SE calculation across bootstrap samples)
  # 2. optimism of performance (comparing performance of boot_models on boot_data vs on original_data)
  bootonboot_perform_list <- rec.list(c(m, B))
  bootonboot_prob_list    <- rec.list(c(m, B))
  optimism_list           <- rec.list(c(m, B))
  
  for (i in 1:m){
    
    for (b in 1:B){
      
      myboot   <- imp_boot_data_list[[i]][[b]]
      
      #a) fit models with and without DTI to bootstrap sample
      boot_models <- fit_mymodels_bio(imp_boot_data_list[[i]][[b]], mymodels)
      
      #b) test models from b) on bootstrap sample
      bootonboot_perform_list[[i]][[b]]  <- test_mymodels(imp_boot_data_list[[i]][[b]], boot_models)[["res"]]
      
      #c) test models from b) on original data
      bootonorig_perform <- test_mymodels(orig_pat_list[[i]], boot_models)[["res"]]
      
      #d) calculate optimism for this bootstrap sample
      optimism_list[[i]][[b]] <- bootonboot_perform_list[[i]][[b]] - bootonorig_perform
      
    }
  }
  
  
  saveRDS(optimism_list, paste0("RDS_bio/optimism_list_", mybio,".rds"))
  saveRDS(bootonboot_perform_list, paste0("RDS_bio/bootonboot_perform_list_", mybio,".rds"))
  
  ###################################################################
  # Calculate optimism corrected model performance
  # for each imputed dataset
  ###################################################################
  
  adj_perform_list <- vector(mode = "list", length = m)
  boot_se_list     <- vector(mode = "list", length = m)
  
  #for each imputed dataset
  for (i in 1:m) {
    
    #calculate mean optimism across bootstrap samples
    arr_opt                 <- array( unlist(optimism_list[[i]]) , 
                                      c(nrow(optimism_list[[i]][[1]]),
                                        ncol(optimism_list[[i]][[1]]),
                                        B) 
    )
    mean_optimism           <- apply(arr_opt, 1:2, mean)
    colnames(mean_optimism) <- colnames(optimism_list[[i]][[1]])
    rownames(mean_optimism) <- mymodels
    
    adj_perform_list[[i]]   <- orig_perform_list[[i]][["res"]] - mean_optimism
    
    #remember NOT to optimism correct R2 again, as the rms package has already done this
    adj_perform_list[[i]][, "nag"]       <- orig_perform_list[[i]][["res"]][, "nag"]
    adj_perform_list[[i]][, "delta_nag"] <- orig_perform_list[[i]][["res"]][, "delta_nag"]
    
    #calculate standard errors
    arr_se                 <- array( unlist(bootonboot_perform_list[[i]]) , 
                                     c(nrow(bootonboot_perform_list[[i]][[1]]),
                                       ncol(bootonboot_perform_list[[i]][[1]]),
                                       B) 
    )
    boot_se_list[[i]]      <- apply(arr_se, 1:2, sd)
    colnames(boot_se_list[[i]]) <- colnames(bootonboot_perform_list[[i]][[1]])
    rownames(boot_se_list[[i]]) <- mymodels
    
  }
  
  #Save model objects
  saveRDS(adj_perform_list, paste0("RDS_bio/adj_perform_list_", mybio, ".rds"))
  saveRDS(boot_se_list, paste0("RDS_bio/boot_se_list", mybio, ".rds"))
  
  
  #################################################
  # Pool optimism-corrected performance estimates
  # and their standard errors
  # across the m imputed datasets
  #################################################
  pooled_results <- pool_myresults(adj_perform_list, boot_se_list, m, mymodels)
  
  
  save_result_as_pretty_table_bio(pooled_results, paste0("pooled_results_", mybio))
  
  
  #################################################
  # extract probabilities of incomplete recovery
  # for mxB datasets
  #################################################
  
  #extract the probabilities of incomplete recovery as calculated by each model
  #for each imputed dataset (in the original unbootstrapped data)
  prob_list <- lapply(orig_perform_list, `[[`, 2)
  prob_list <- lapply(prob_list, t)
  
  #add in a column with the corresponding observed outcome
  for (i in 1:length(prob_list)){
    prob_list[[i]] <- cbind((orig_pat_list[[i]]$Recovery), prob_list[[i]])
    colnames(prob_list[[i]]) <- c("Recovery", mymodels)
    prob_list[[i]] <- as.data.frame(prob_list[[i]])
  }
  
  #save model object
  saveRDS(prob_list, paste0("RDS_bio/prob_list", mybio, ".rds"))
  
  
  #############################################################
  # display model pooled coefficients for
  # all models with and without DTI
  # fitted to the original (imputed but not bootstrapped) data
  #############################################################
  
  save_model_coefficients_bio(orig_coef_list, mybio)
  
}




