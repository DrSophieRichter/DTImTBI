#Functions for script DTI_mTBI#

###########################################
#rendering functions for table1
###########################################

#set up display of continuous variables
my.render.cont <- function(x) 
{
  #display 2 sig figs:
  with(stats.apply.rounding(stats.default(x), digits=2), 
       #in first column display variable name and "Median (IQR)":
       c("", "Median (Min-Max)" = 
           #in second column calculate and display:
           sprintf("%s (%s %s %s)", 
                   MEDIAN,  
                   round(min(x, na.rm = T),1) , 
                   "-", 
                   round(max(x, na.rm = T),1))))   #calculations requested
}


#set up display of categorical variables
my.render.cat <- function(x) 
{
  #in first column display variable name:
  c("", 
    #in second column calculate and display:
    sapply(stats.default(x), function(y) with(y, sprintf("%d (%0.0f %%)", 
                                                         FREQ, 
                                                         round(PCT, 0)))))
}


#######################################
#Now adjust values for age
########################################

adjust_for_age <- function(mydata) {
  
  #I know dti in wm changes with age in a quadratic manner:
  #I know age adjustment needs to be ROI specific: 
  
  #Convert DTI columns to numeric
  temp <- 
    mydata %>%
    select(contains("ROI"))
  mycols <- colnames(temp)
  mydata[mycols] <- sapply(mydata[mycols],as.numeric)
  
  for (i in 1:length(mycols)){
    
    mycol <- mycols[i]
    
    #model dti ~ poly(age,2)
    fit <- lm(get(mycol) ~ poly(Age_at_this_scan, 2), data = mydata %>% filter(Category != "TBI"))
    
    #get predictions for each patient using their actual age (e.g. 0.6 for 20 yo)
    pred_real_age <- predict(fit, newdata = mydata)
    
    #get predictions for each patient using 43 years of age (e.g. 0.3 for 50 yo)
    mydata50                  <- mydata
    mydata50$Age_at_this_scan <- 43
    pred_50                   <- predict(fit, newdata = mydata50)
    
    #get the ratio standard/actual (e.g. 0.5)
    delta_pred <- pred_50 - pred_real_age
    
    #adjusted dti = observed dti * ratio (e.g. 0.5 * 0.5 = 0.25)
    newname <- paste0("adj_", mycol)
    mydata[, newname] <- mydata[,mycol] + delta_pred
    
  }
  
  return(mydata)
  
}
#note I checked model assumption on individual ROIs before rolling out as a function


#######################################
#Sensitivity analysis
#do NOT adjust values for age
########################################

do_not_adjust_for_age <- function(mydata) {
  
  #I know dti in wm changes with age in a quadratic manner:
  #I know age adjustment needs to be ROI specific: 
  
  #Convert DTI columns to numeric
  temp <- 
    mydata %>%
    select(contains("ROI"))
  mycols <- colnames(temp)
  mydata[mycols] <- sapply(mydata[mycols],as.numeric)
  
  for (i in 1:length(mycols)){
    
    mycol <- mycols[i]
    
    #just change name without adjusting for age
    newname <- paste0("adj_", mycol)
    mydata[, newname] <- mydata[,mycol]
    
  }
  
  return(mydata)
  
}

###################################
#Harmonize imaging data
###################################

harmonize_data <- function(mydata){
  
  #batch(Scanner identifier)
  mydata$Scanner <- droplevels(mydata$Scanner)
  batch <- as.factor(mydata$Scanner)
  
  #dat (imaging features)
  fa_dat <- 
    mydata %>%
    select(contains("FA_ROI")) %>%
    t()
  
  md_dat <- 
    mydata %>%
    select(contains("MD_ROI")) %>%
    t()
  
  
  #mod (biological variables to preserve)
  mydata$Age_squared <- mydata$Age_at_this_scan^2
    mod <- 
      model.matrix(~ Age_at_this_scan + 
                     Age_squared +
                     Sex +
                     GCSScoreBaselineDerived +
                     MR_abnormality +
                     Time_to_scan, 
                   mydata)
  
  
  #harmonize fa data
  fa_res <- neuroCombat(dat = fa_dat, 
                        batch = batch, 
                        mod = mod, 
                        verbose = FALSE, 
                        mean.only = FALSE)
  
  fa_harmo <- 
    cbind(as.character(mydata$Scan_ID), t(fa_res[[1]])) %>% 
    as.data.frame()
  colnames(fa_harmo)[1] <- c("Scan_ID")
  
  #harmonize md data
  md_res <- neuroCombat(dat = md_dat, 
                        batch = batch, 
                        mod = mod, 
                        verbose = FALSE, 
                        mean.only = FALSE)
  
  md_harmo <- 
    cbind(as.character(mydata$Scan_ID), t(md_res[[1]])) %>% 
    as.data.frame()
  colnames(md_harmo)[1] <- c("Scan_ID")
  
  #replace unharmonized with harmonized data
  mydata <-
    mydata %>%
    select(-contains("ROI"))
  mydata <-
    merge(mydata, fa_harmo, by = "Scan_ID")
  mydata <-
    merge(mydata, md_harmo, by = "Scan_ID")
  
  return(mydata)
  
}

#####################################################################
#scale and center patient dti data
#to yield standardized (comparable) lasso coefficents later
#else MD coefficients are going to be orders of magnitude higher
#####################################################################

scale_and_center <- function(pat_list){
  
  #identify dti variable names
  temp <-
    pat_list[[1]] %>%
    select(contains("adj"))
  cols <-  colnames(temp)
  
  #standardize outcome and patient data
  for (i in 1:length(pat_list)) {
    
    #select dti variables
    mypat        <- 
      pat_list[[i]][, c(cols)]
    
    #scale and center data
    pre_proc_val    <- preProcess(mypat[,cols], method = c("center", "scale"))
    mypat[,cols]    <- predict(pre_proc_val, mypat[,cols])
    colnames(mypat) <- paste0("stan_", colnames(mypat))
    
    #add scaled and centred data to pat_list
    pat_list[[i]]   <- bind_cols(pat_list[[i]], mypat)
  }
  
  return(pat_list)
  
}

################################################
#
# choose optimal panel for ROIs of FA and MD
# using lasso regression
# using multiple imputed datasets
#
################################################

select_rois_multimp <- function(pat_list, myoutcome) {
  
  m <- length(pat_list)
  
  #initialize vector to store lambdas from m imputed datasets
  my_lambdas <- vector()
  
  #variables to consider for selection
  temp <-
    pat_list[[1]] %>%
    select(contains("stan_adj"))
  cols = colnames(temp)
  
  #1.	Find one lambda per imputed dataset (using 10 x 10fold cross-validation). 
  #This will yield 10 lambdas
  for (pat_num in 1:length(pat_list)) {
    
    #select one imputed dataset
    mypat        <- 
      pat_list[[pat_num]][, c(myoutcome, cols)]
    
    #make model matrix
    x <-  as.matrix(mypat[, -1])
    y <-  mypat[,1]
    
    #determine lambda using 10x repeated 10fold cross-validation
    #note glmnet standardizes variables by default
    MSEs <- NULL
    for (i in 1:10){
      cv   <- cv.glmnet(x = x, y = y, alpha=1, nfolds=10)  
      MSEs <- cbind(MSEs, cv$cvm)
    }
    rownames(MSEs)        <- cv$lambda
    lambda_min            <- as.numeric(names(which.min(rowMeans(MSEs))))
    my_lambdas[[pat_num]] <- lambda_min
  }
  
  
  #2.	I take the mean lambda as my optimal_lambda
  optimal_lambda <- mean(my_lambdas) 
  
  #3.	I stack all 10 imputed datasets
  mystack        <- 
    bind_rows(pat_list) %>% 
    select(all_of(myoutcome), contains("stan_adj"))
  
  #make model matrix
  x <-  as.matrix(mystack[, -1])
  y <-  mystack[,1]
  
  #4.	I fit one lasso model to the stacked data, using 
  #    a.	Weights = 0.1
  #    b.	Lambda = optimal_lambda
  myweights   <- rep((1/m), nrow(mystack))
  lasso_model <- glmnet(x, y, 
                        alpha = 1, #for lasso
                        lambda = optimal_lambda, #mean of 10 imputed datasets
                        weights = myweights #weighted to not inflate sample size
  )
  
  #Extract coefficients
  mycoefs <- 
    coef(lasso_model) %>%
    as.matrix() %>%
    as.data.frame()
  
  colnames(mycoefs) <- myoutcome
  
  res        <- list(lasso_model, mycoefs)
  names(res) <- c("model", "coefs")
  
  return(res)
}

################################################
#
# choose optimal panel for ROIs of FA and MD
# using lasso regression
# using s single dataset
#
################################################

select_rois <- function(data, myoutcome) {
  
  #variables to consider for selection
  temp <-
    data %>%
    select(contains("stan_adj"))
  cols = colnames(temp)
  
  #Find lambda (using 10 x 10fold cross-validation)
  
  #select relevant columns from my data
  mypat        <- 
    data[, c(myoutcome, cols)]
  
  #make model matrix
  x <-  as.matrix(mypat[, -1])
  y <-  mypat[,1]
  
  #determine lambda using 10x repeated 10fold cross-validation
  #note glmnet standardizes variables by default
  MSEs <- NULL
  for (i in 1:10){
    cv   <- cv.glmnet(x = x, y = y, alpha=1, nfolds=10)  
    MSEs <- cbind(MSEs, cv$cvm)
  }
  rownames(MSEs)        <- cv$lambda
  lambda_min            <- as.numeric(names(which.min(rowMeans(MSEs))))
  
  
  #fit lasso model to the data, using lambda_min
  lasso_model <- glmnet(x, y, 
                        alpha = 1, #for lasso
                        lambda = lambda_min
  )
  
  #Extract coefficients
  mycoefs <- 
    coef(lasso_model) %>%
    as.matrix() %>%
    as.data.frame()
  
  colnames(mycoefs) <- myoutcome
  
  res        <- list(lasso_model, mycoefs)
  names(res) <- c("model", "coefs")
  
  return(res)
}



######################################################
#
# Fit models with and without DTI_score
#
######################################################

fit_mymodels <- function(df, mymodels){
  
  models_imp <- vector(mode = "list", length = length(mymodels))
  names(models_imp) <- mymodels
  
  if("DTI_only" %in% mymodels){
    models_imp[["DTI_only"]] <-
      glm(Recovery ~ 
            DTI_score, 
          data = df,
          family = "binomial")
  }

  if("plain_upfrontEM" %in% mymodels){
    
    models_imp[["plain_upfrontEM"]] <- 
      glm(Recovery ~ 
            Sex +
            Education*Age_at_this_scan +
            HxMH +
            InjViolenceVictimAlcohol +
            GCSScoreBaselineDerived +
            PTAover1h, 
          data = df,
          family = "binomial")
    
    models_imp[["dti_upfrontEM"]] <- 
      glm(Recovery ~ 
            Sex +
            Education*Age_at_this_scan +
            HxMH +
            InjViolenceVictimAlcohol +
            GCSScoreBaselineDerived +
            PTAover1h +
            DTI_score, 
          data = df,
          family = "binomial")
  }
 
  if("plain_upfrontPLUS" %in% mymodels){
    
    models_imp[["plain_upfrontPLUS"]] <- 
      glm(Recovery ~ 
            Education +
            HxMH +
            InjViolenceVictimAlcohol +
            GCSScoreBaselineDerived +
            PTAover1h +
            wk2_GAD7TotalScore +
            wk2_PHQ9TotlScre +
            wk2_RPQTotalScore, 
          data = df,
          family = "binomial")
    
    models_imp[["dti_upfrontPLUS"]] <- 
      glm(Recovery ~ 
            Education +
            HxMH +
            InjViolenceVictimAlcohol +
            GCSScoreBaselineDerived +
            PTAover1h +
            wk2_GAD7TotalScore +
            wk2_PHQ9TotlScre +
            wk2_RPQTotalScore +
            DTI_score, 
          data = df,
          family = "binomial")
  }
  
  if("plain_headsmart" %in% mymodels){
    
    models_imp[["plain_headsmart"]] <- 
      glm(Recovery ~ 
            Age_at_this_scan +
            HxMH +
            acute_RPQHeadaches +
            acute_RPQPoorConcentration +
            acute_RPQLightSensitivity, 
          data = df,
          family = "binomial")
    
    models_imp[["dti_headsmart"]] <- 
      glm(Recovery ~ 
            Age_at_this_scan +
            HxMH +
            acute_RPQHeadaches +
            acute_RPQPoorConcentration +
            acute_RPQLightSensitivity +
            DTI_score, 
          data = df,
          family = "binomial")
  }
  
  if("plain_centerEM" %in% mymodels){
    
    models_imp[["plain_centerEM"]] <- 
      glm(Recovery ~ 
            Age_at_this_scan +
            ISS +
            GCSScoreBaselineDerived +
            Sex +
            HxMH +
            ASA +
            Cause +
            acute_RPQTotalScore,
          data = df,
          family = "binomial")
    
    models_imp[["dti_centerEM"]] <- 
      glm(Recovery ~ 
            Age_at_this_scan +
            ISS +
            GCSScoreBaselineDerived +
            Sex +
            HxMH +
            ASA +
            Cause +
            acute_RPQTotalScore +
            DTI_score, 
          data = df,
          family = "binomial")
  }
  
  if("plain_centerPLUS" %in% mymodels){
    
    models_imp[["plain_centerPLUS"]] <- 
      glm(Recovery ~ 
            Age_at_this_scan +
            ISS +
            GCSScoreBaselineDerived +
            HxMH +
            Cause +
            wk2_RPQTotalScore +
            wk2_PCL5TotalScore,
          data = df,
          family = "binomial")
    
    models_imp[["dti_centerPLUS"]] <- 
      glm(Recovery ~ 
            Age_at_this_scan +
            ISS +
            GCSScoreBaselineDerived +
            HxMH +
            Cause +
            wk2_RPQTotalScore +
            wk2_PCL5TotalScore +
            DTI_score, 
          data = df,
          family = "binomial")
  }
  
    return(models_imp)
  
}

#############################################################
#
# Return NA if matrix singular for rms::val.prob
#
#############################################################

get_nag <- function(j, prob_imp, df){
  tryCatch(
    {nag  <- rms::val.prob(p = prob_imp[[j]], y = df$Recovery, pl = FALSE)
    nag <- nag["R2"][[1]]*100
    return(nag)
    
    },
    error = function(e){
      NA
    }
  )
}

##############################################################
#
# Test the chosen models on the chosen dataset
#
##############################################################


test_mymodels <- function(df, models_imp){
  
  #Make predictions on imp data
  
  prob_imp             <- vector(mode = "list", length = length(mymodels))
  names(prob_imp)      <- mymodels
  predclass_imp        <- vector(mode = "list", length = length(mymodels))
  names(predclass_imp) <- mymodels
  
  for (j in 1:length(mymodels)){
    #probability of incomplete recovery
    prob_imp[[j]]  <- predict(object = models_imp[[j]],
                              newdata = df,
                              type = "response")
    #predicted recovery (0 or 1)
    predclass_imp[[j]]  <- ifelse(prob_imp[[j]] <= 0.5, 0, 1)
  }
  
  #Now get all the performance measures I want to report
  # so that they can be corrected for optimism
  # AUC
  # Nagelkerke R2
  # sensitivity & specificty
  # PPV and NPV
  # p-value LRT
  
  mycols <- c("auc", "hl", "brier" ,"nag","sens", "spec","ppv", "npv", "lrt", "delta_auc", "delta_nag")
  
  res <- matrix(, ncol = length(mycols), nrow = length(mymodels))
  colnames(res) <- mycols
  rownames(res) <- mymodels
  
  for (j in 1:length(mymodels)){
    
    #AUC
    res[j, "auc"] <- ci.auc(roc(df$Recovery ~ prob_imp[[j]], quiet = TRUE))[2]
    
    #Hosmer-Lemeshow test
    res[j, "hl"]  <- NA
    
    #Brier score
    res[j, "brier"]     <- DescTools::BrierScore(resp = df$Recovery,
                                                 pred = prob_imp[[j]])
    
    #Nagelkerke R2
    res[j, "nag"] <- get_nag(j, prob_imp, df)
    
    #Sensitivity
    res[j, "sens"] <- sensitivity(factor(predclass_imp[[j]]), 
                                  factor(df$Recovery),
                                  positive = "1")
    #Specificity
    res[j, "spec"] <- specificity(factor(predclass_imp[[j]]), 
                                  factor(df$Recovery), 
                                  negative = "0")
    
    #PPV
    res[j, "ppv"] <- posPredValue(factor(predclass_imp[[j]]), 
                                  factor(df$Recovery),
                                  positive = "1")
    #NPV
    res[j, "npv"] <- negPredValue(factor(predclass_imp[[j]]), 
                                  factor(df$Recovery), 
                                  negative = "0")
    
    #LRT (comparing nested models)
    a <- (length(mymodels) + 1)/2
    b <- (length(mymodels) - 1)/2
    c <- 0.0000000001 # constant in case AUC = 1, logit of which would be Inf
    if (j <= a){ 
      res[j, "lrt"] <- NA
      res[j, "delta_auc"] <- NA
      res[j, "delta_nag"] <- NA
    } else {
      temp          <- lmtest::lrtest(models_imp[[j-b]], models_imp[[j]]) 
      res[j, "lrt"] <- temp[2, "Pr(>Chisq)"] 
      res[j, "delta_auc"] <- LaplacesDemon::logit(I(res[j, "auc"] - c)) - LaplacesDemon::logit(I(res[j-b, "auc"] - c))
      res[j, "delta_nag"] <- res[j, "nag"] - res[j-b, "nag"]
    }
    
  }
  
  prob            <- bind_rows(prob_imp)
  res_list        <- list(res, prob)
  names(res_list) <- c("res", "prob")
  return(res_list)
  
}



##################################################################################
#
# Use Rubin's rules to pool results
# across m imputed datasets
#
##################################################################################

pool_myresults <- function(estimate_list, se_list, m, mymodels){
  #Pool results across imputed datasets
  #using Rubin's rules
  mycols <- colnames(estimate_list[[1]])
  mymodels <- rownames(estimate_list[[1]])
  
  #####For point estimates that is the same as taking the mean
  
  arr_estimate <- array(unlist(estimate_list) , c(length(mymodels),length(mycols),m) )
  pooled_estimate <- apply(arr_estimate, 1:2, mean, na.rm = TRUE)
  colnames(pooled_estimate) <- mycols
  rownames(pooled_estimate) <- mymodels
  
  
  #####For SE Rubin's rules are more complicated
  
  #first calculate within imputation variance
  squ_fun <- function(x){
    x^2
  }
  temp <- lapply(se_list, squ_fun)
  
  arr_within <- array(unlist(temp) , c(length(mymodels),length(mycols),m) )
  within_variance <- apply(arr_within, 1:2, mean, na.rm = TRUE)
  colnames(within_variance) <- mycols
  rownames(within_variance) <- mymodels
  
  
  #then calculate between imputation variance
  sub_fun <- function(x){
    x - pooled_estimate
  }
  
  delta_estimate <- lapply(estimate_list, sub_fun)
  delta_estimate_squ <- lapply(delta_estimate, squ_fun)
  arr_between <- array(unlist(delta_estimate_squ) , c(length(mymodels),length(mycols),m) )
  
  spec_fun <- function(x){
    sum(x, na.rm = TRUE)/(length(x)-1)
  }
  between_variance <- apply(arr_between, 1:2, spec_fun)
  colnames(between_variance) <- mycols
  rownames(between_variance) <- mymodels
  
  
  #Then calculate total variance and pooled standard error
  total_variance <- within_variance + between_variance + I(between_variance/m)
  pooled_se      <- sqrt(total_variance)
  
  
  #calculate the limits of the 95% confidence interval
  pooled_ll <- pooled_estimate - I(1.96*pooled_se)
  pooled_ul <- pooled_estimate + I(1.96*pooled_se)
  
  #set sensible limits for confidence intervals of
  #auc (0-1)
  pooled_ll[,"auc"] <- ifelse(pooled_ll[,"auc"] < 0, 0, pooled_ll[,"auc"])
  pooled_ul[,"auc"] <- ifelse(pooled_ul[,"auc"] > 1, 1, pooled_ul[,"auc"])
  #nag (-100 - +100)
  pooled_ll[,"nag"] <- ifelse(pooled_ll[,"nag"] < -100, -100, pooled_ll[,"nag"])
  pooled_ul[,"nag"] <- ifelse(pooled_ul[,"nag"] > 100, 100, pooled_ul[,"nag"])
  
  #collect results in a list
  pooled_results <- list(pooled_estimate, pooled_ll, pooled_ul, pooled_se)
  names(pooled_results) <- c("pooled_estimate", "pooled_ll", "pooled_ul", "pooled_se")
  
  return(pooled_results)
}


#####################################################
#
# Round figures in results table
#
#####################################################

round_fun <- function(df){
  #round p-values to 3 digits
  #auc, slope, int 2
  df[,c("hl", "lrt")] <- 
    format(round(df[,c("hl", "lrt")],3), 
           nsmall = 3)
  df[, c("auc", "brier", "delta_auc", "sens", "spec", "ppv", "npv")] <- 
    format(round(as.numeric(df[, c("auc", "brier", "delta_auc","sens", "spec", "ppv", "npv")]), 2),
           nsmall = 2)
  df[, c("nag", "delta_nag")] <-
    format(round(as.numeric(df[, c("nag", "delta_nag")]), 0),
           nsmall = 0)
  df <- 
    df %>%
    trimws()
  return(df)
}

####################################################
#
# save myres as a pretty table in word
#
####################################################

save_result_as_pretty_table <- function(myres, mytitle) {
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
  temp[c(1:6), "lrt"] <- p.adjust(temp[c(1:6), "lrt"], method = "fdr")
  
  
  #reorder rows
  myorder <- c("DTI_only",
               "plain_upfrontEM",
               "dti_upfrontEM",
               "plain_upfrontPLUS",
               "dti_upfrontPLUS",
               "plain_headsmart",
               "dti_headsmart",
               "plain_centerEM",
               "dti_centerEM",
               "plain_centerPLUS",
               "dti_centerPLUS")
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
    set_header_labels("plain_upfrontEM" = "w/o DTI",
                      "plain_upfrontPLUS" = "w/o DTI",
                      "plain_headsmart" = "w/o DTI",
                      "plain_centerEM" = "w/o DTI",
                      "plain_centerPLUS" = "w/o DTI",
                      "dti_upfrontEM" = "with DTI",
                      "dti_upfrontPLUS" = "with DTI",
                      "dti_headsmart" = "with DTI",
                      "dti_centerEM" = "with DTI",
                      "dti_centerPLUS" = "with DTI",
                      "Metric" = " ",
                      "DTI_only" = " ") %>%
    add_header_row(values = c("Metric",
                              "DTI only", 
                              "UPFRONT-ED", 
                              "UPFRONT-PLUS", 
                              "HeadSMART", 
                              "CENTER-ED",
                              "CENTER-PLUS"),
                   colwidths = c(1, 1, 2, 2, 2, 2, 2))%>%
    vline(j = c(1, 2, 4, 6, 8, 10, 12))%>%
    vline_left() %>%
    vline_right() %>%
    align(align = "right", part = "all")%>%
    align(align = "center", part = "header", i = c(1:2)) %>%
    align(align = "left", j = c(1), part = "all") %>% 
    padding(padding = 0, part = "all", padding.right = 3, padding.left = 3, padding.top = 3, padding.bottom = 3) %>%
    fontsize(size = 10, part = "all") %>%
    fix_border_issues(part = "all") 
  
  save_as_docx(ft5, path = paste0("Figures/performance_", mytitle, ".docx"))
  
  
}

###############################################################
# 
# Use prob_list to make a panel
# of calibration plots
#
###############################################################

make_calibration_plots <- function(prob_list, mytitle, midline_dot){
  
  #calculate the mean probability of recovery across the 10 imputed datasets
  #so that the calibration plot is based on the same nuber of patients
  # as the original sample (so that standard errors are not shrunk)
  arr_prob <- 
    array( unlist(prob_list) , 
           c(nrow(prob_list[[1]]), ncol(prob_list[[1]]), m) )
  
  myprobs  <- 
    apply(arr_prob, 1:2, mean) %>% 
    as.data.frame()
  colnames(myprobs) <- colnames(prob_list[[1]])
  
  #make individual plots and store them in a list
  calplot_list <- lapply(c(2:ncol(myprobs)), FUN = function(x) {
    
    predvar    <- names(myprobs)[x]
    myoutcome  <- "Recovery"
    
    plot_cal(mydata = myprobs,
             predvar = predvar,
             myoutcome = myoutcome,
             mycolor = "red",
             midline_dot = midline_dot)
  })
  
  names(calplot_list) <- colnames(myprobs)[-1]
  
  
  #arrange plots in a panel and save image
  calplot_panel <- 
    grid.arrange(top  = textGrob("Model without DTI                                Model with DTI", 
                                 gp = gpar(fontsize = 15)),
                 left = textGrob("CENTER-PLUS                           CENTER-ED                          HeadSMART                            UPFRONT-PLUS                        UPFONT-ED", 
                                 gp = gpar(fontsize = 15), 
                                 rot = 90),
                 calplot_list[["plain_upfrontEM"]],
                 calplot_list[["dti_upfrontEM"]], 
                 calplot_list[["plain_upfrontPLUS"]], 
                 calplot_list[["dti_upfrontPLUS"]], 
                 calplot_list[["plain_headsmart"]], 
                 calplot_list[["dti_headsmart"]],
                 calplot_list[["plain_centerEM"]], 
                 calplot_list[["dti_centerEM"]],
                 calplot_list[["plain_centerPLUS"]], 
                 calplot_list[["dti_centerPLUS"]],
      ncol=2, nrow=5, 
      widths = c(4,4),
      heights = c(4,4,4,4,4))
  
  ggsave(filename = paste0("Figures/calplot_", mytitle, ".pdf"), 
         plot = calplot_panel,
         dpi = 300,
         width = 18,
         height = 37,
         units = "cm",
         limitsize = FALSE)
  
}


###########################################################
#
# Make calibration plot
#
###########################################################

plot_cal <- function(mydata, predvar, myoutcome, mycolor, midline_dot){
  temp <- mutate(mydata, bin = ntile(get(predvar), 10)) %>% 
    # Bin prediction into 10ths
    group_by(bin) %>%
    mutate(n = n(), # Get ests and CIs
           bin_pred = mean(get(predvar)), 
           bin_prob = mean(as.numeric(get(myoutcome))), 
           se = sqrt((bin_prob * (1 - bin_prob)) / n), 
           ul = bin_prob + 1.96 * se, 
           ll = bin_prob - 1.96 * se) %>%
    ungroup()
  
  mod <- lm(bin_prob ~ bin_pred, data = temp)
  intercept <- coef(mod)[["(Intercept)"]]
  intercept <- format(round(intercept,3), nsmall = 3)
  slope <- coef(mod)[["bin_pred"]]
  slope <- format(round(slope,3), nsmall = 3)
  
  
  if (midline_dot == "yes") {
    
    slope <- gsub(".", "·", slope, fixed=TRUE)
    intercept <- gsub(".", "·", intercept, fixed=TRUE)
    
    g1 <- ggplot(data = temp, aes(x = bin_pred, y = bin_prob, ymin = ll, ymax = ul)) +
      geom_pointrange(size = 0.5, color = "black") +
      scale_y_continuous(breaks = seq(-0.2, 1.2, by = 0.2),
                         labels = scales::label_number(decimal.mark = "·")) +
      coord_cartesian(ylim = c(0, 1),) +
      scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2),
                         labels = scales::label_number(decimal.mark = "·")) +
      geom_abline() + # 45 degree line indicating perfect calibration
      geom_smooth(method = "lm", se = FALSE, linetype = "dashed", 
                  color = "black", formula = y~-1 + x) + 
      # straight line fit through estimates
      geom_smooth(data = temp,
                  aes(x = get(predvar), y = as.numeric(get(myoutcome))), 
                  color = mycolor, se = FALSE, method = "loess") + 
      # loess fit through estimates
      xlab("") +
      ylab("Observed Probability") +
      theme_minimal() +
      annotate("text", 
               x= 0.03, y= 0.85, hjust = 0,
               label= paste("slope =", slope)) + 
      annotate("text", 
               x = 0.03, y=0.95, hjust = 0, 
               label = paste("intercept =", intercept))
    
    # The distribution plot  
    g2 <- 
      ggplot(data = temp, aes(x = get(predvar), y = after_stat(density))) +
      geom_histogram(fill = "black", bins = 30) +
      scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2),
                         labels = scales::label_number(decimal.mark = "·")) +
      scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, by = 2)) +
      xlab("Predicted Probability") +
      ylab("%") +
      theme_minimal() +
      theme(panel.grid.minor = element_blank())
    
  } else {
    g1 <- ggplot(data = temp, aes(x = bin_pred, y = bin_prob, ymin = ll, ymax = ul)) +
      geom_pointrange(size = 0.5, color = "black") +
      scale_y_continuous(breaks = seq(-0.2, 1.2, by = 0.2)) +
      coord_cartesian(ylim = c(0, 1),) +
      scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
      geom_abline() + # 45 degree line indicating perfect calibration
      geom_smooth(method = "lm", se = FALSE, linetype = "dashed", 
                  color = "black", formula = y~-1 + x) + 
      # straight line fit through estimates
      geom_smooth(data = temp,
                  aes(x = get(predvar), y = as.numeric(get(myoutcome))), 
                  color = mycolor, se = FALSE, method = "loess") + 
      # loess fit through estimates
      xlab("") +
      ylab("Observed Probability") +
      theme_minimal() +
      annotate("text", 
               x= 0.03, y= 0.85, hjust = 0,
               label= paste("slope =", slope)) + 
      annotate("text", 
               x = 0.03, y=0.95, hjust = 0, 
               label = paste("intercept =", intercept))
    
    # The distribution plot  
    g2 <- 
      ggplot(data = temp, aes(x = get(predvar), y = after_stat(density))) +
      geom_histogram(fill = "black", bins = 30) +
      scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
      scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, by = 2)) +
      xlab("Predicted Probability") +
      ylab("%") +
      theme_minimal() +
      theme(panel.grid.minor = element_blank())
    
  }
  
  # Combine them    
  myplot <- cowplot::plot_grid(g1, NULL, g2, align = "hv", nrow = 3, rel_heights = c(3/4, -0.06, 1/4))
  
  return(myplot)
}

####################################################################
#
# Create panels of ROC curves and save as images
#
####################################################################

make_roc_plots <- function(prob_list, myres, mytitle, midline_dot) {
  
  myprobs <- bind_rows(prob_list)
  
  #Make individual ROC curves and store in a list
  roc_list <- vector(mode = "list", length = length(mymodels))
  names(roc_list) <- mymodels
  
  for (i in 1:length(mymodels)){
    roc_list[[i]] <- roc(unlist(myprobs[, "Recovery"]), 
                         unlist(myprobs[, mymodels[i]]))
  }
  
  aucs <- formatC(myres[["pooled_estimate"]][, "auc"], digits = 2, format = "f")
  lls  <- formatC(myres[["pooled_ll"]][, "auc"], digits = 2, format = "f")
  uls  <- formatC(myres[["pooled_ul"]][, "auc"], digits = 2, format = "f")
  aucs <- paste0(aucs, " (", lls, "-", uls, ")")
  
  
  rocplot_list <- vector(mode = "list", length = (length(mymodels)/2 + 1))
  names(rocplot_list) <- mymodels[str_detect(mymodels, "plain")==FALSE]
  
  for (i in 2:length(rocplot_list)){
    
    with    <- i+5
    without <- i
    rocplot_list[[i]] <- 
      make_one_roc_plot(
        roc_with = roc_list[[with]],
        roc_without = roc_list[[without]],
        color_with = "orange",
        color_without = "darkcyan",
        auc_with = aucs[with],
        auc_without = aucs[without],
        midline_dot = midline_dot)
    
  }
  
  rocplot_list[[1]] <-
    make_one_roc_plot_DTI_only(
      roc_with = roc_list[[1]],
      color_with = "orange",
      auc_with = aucs[1],
      midline_dot = midline_dot)
  
  #arrange roc plots in a panel and save image
  rocplot_panel <- 
    grid.arrange(
      arrangeGrob(rocplot_list[[1]], top = "DTI only"),
      arrangeGrob(rocplot_list[[2]], top = "UPFRONT-ED"), 
      arrangeGrob(rocplot_list[[3]], top = "UPFRONT-PLUS"),
      arrangeGrob(rocplot_list[[4]], top = "HeadSMART"),
      arrangeGrob(rocplot_list[[5]], top = "CENTER-ED"),
      arrangeGrob(rocplot_list[[6]], top = "CENTER-PLUS"),
      ncol=3, nrow=2)
  
  ggsave(filename = paste0("Figures/rocplot_", mytitle,".pdf"), 
         plot = rocplot_panel,
         dpi = 300,
         width = 20,
         height = 13,
         units = "cm")
  
}


###########################################################
#
# Make ROC plot 
# two helper function to draw individual ROC plots
#
###########################################################

make_one_roc_plot <- function(roc_with, roc_without, color_with, color_without, auc_with, auc_without, midline_dot){
  roclist <- list("roc_with" =  roc_with, "roc_without" =  roc_without)
  
  myplot <- 
    ggroc(roclist, legacy.axes= T, size=0.6, aes = "colour") +
    theme_classic(base_size = 6) +
    theme(legend.justification=c(1,0), 
          legend.position=c(0.98,0.08), 
          legend.text=element_text(size=6),
          legend.title=element_text(size=6),
          legend.key.size = unit(0.6,"line"),
          plot.title = element_text(size = 9, hjust = 0.5),
          axis.title.x = element_text(size = 6),
          axis.title.y = element_text(size = 6),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6)) +
    geom_abline(intercept = 0, 
                slope = 1, 
                color = "darkgrey", 
                linetype = "dashed") 
  
  if (midline_dot == "yes") {
    
    auc_with    <- gsub(".", "·", auc_with, fixed=TRUE)
    auc_without <- gsub(".", "·", auc_without, fixed=TRUE)
    
    myplot <- 
      myplot +
      scale_y_continuous(labels = scales::label_number(decimal.mark = "·")) +
      scale_x_continuous(labels = scales::label_number(decimal.mark = "·")) +
      scale_colour_manual(name   = "AUC (95% CI)",
                          values = c("roc_with" =  color_with, 
                                     "roc_without" = color_without),
                          labels = c(paste0("with DTI: ", auc_with),
                                     paste0("without DTI: ", auc_without))) 
  } else {
    myplot <- 
      myplot +
      scale_colour_manual(name   = "AUC (95% CI)",
                          values = c("roc_with" =  color_with, 
                                     "roc_without" = color_without),
                          labels = c(paste0("with DTI: ", auc_with),
                                     paste0("without DTI: ", auc_without)))
  }
  
  return(myplot)
}


###

make_one_roc_plot_DTI_only <- function(roc_with, color_with, color_without, auc_with, auc_without, midline_dot){
  roclist <- list("roc_with" =  roc_with)
  
  myplot <- 
    ggroc(roclist, legacy.axes= T, size=0.6, aes = "colour") +
    theme_classic(base_size = 6) +
    theme(legend.justification=c(1,0), 
          legend.position=c(0.98,0.08), 
          legend.text=element_text(size=6),
          legend.title=element_text(size=6),
          legend.key.size = unit(0.6,"line"),
          plot.title = element_text(size = 9, hjust = 0.5),
          axis.title.x = element_text(size = 6),
          axis.title.y = element_text(size = 6),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6)) +
    geom_abline(intercept = 0, 
                slope = 1, 
                color = "darkgrey", 
                linetype = "dashed")
  
  if (midline_dot == "yes") {
    
    auc_with    <- gsub(".", "·", auc_with, fixed=TRUE)
    
    myplot <- 
      myplot +
      scale_y_continuous(labels = scales::label_number(decimal.mark = "·")) +
      scale_x_continuous(labels = scales::label_number(decimal.mark = "·")) +
      scale_colour_manual(name   = "AUC (95% CI)",
                          values = c("roc_with" =  color_with),
                          labels = c(paste0("with DTI: ", auc_with))
      ) 
  } else {
    myplot <- 
      myplot +
      scale_colour_manual(name   = "AUC (95% CI)",
                          values = c("roc_with" =  color_with),
                          labels = c(paste0("with DTI: ", auc_with))
      )
  }
  
  
  return(myplot)
}


###########################################################
#
# Format pooled coefficients
#
###########################################################

format_pooled_coefs <- function(pooled_coefs, model){
  
  temp <- format(round(pooled_coefs, 3), nsmall = 3)
  temp <- temp %>% as.data.frame()
  
  temp$Variable <- rownames(temp)
  temp$Model    <- model
  temp <- 
    temp %>%
    select(Model, Variable, everything())
  
  colnames(temp) <- gsub("Pr(>|z|)", "p-value", colnames(temp), fixed = TRUE)
  temp[, "p-value"] <- ifelse(as.numeric(temp[, "p-value"]) < 0.001, "<0.001", temp[, "p-value"])
  temp[, "p-value"] <- as.character(temp[, "p-value"])
  
  return(temp)
}

###########################################################
#
#Pool coefficients across imputed datasets
#
###########################################################

pool_my_coefs <- function(m_list){
  
  m <- length(m_list)
  
  #for each imputed dataset
  estimate_list   <- vector(mode = "list", length = m)
  se_within_list  <- vector(mode = "list", length = m)
  
  for (i in 1:m) {
    df <- summary(m_list[[i]])$coefficients %>% as.matrix()
    
    #extract the mean estimate
    estimate_list[[i]]           <- t(t(df[,1]))
    colnames(estimate_list[[i]]) <- "Estimate"
    
    #extract SE error of mean estimate (=within imputation variance)
    se_within_list[[i]] <- t(t(df[,2]))
    colnames(se_within_list[[i]]) <- "SE_within"
  }
  
  #combine with Rubin's rules
  #pool mean estimates
  arr_estimate <- array(unlist(estimate_list) , c(nrow(estimate_list[[1]]),ncol(estimate_list[[1]]),m) )
  pooled_estimate <- apply(arr_estimate, 1:2, mean)
  colnames(pooled_estimate) <- colnames(estimate_list[[1]])
  rownames(pooled_estimate) <- rownames(estimate_list[[1]])
  
  #calculate mean within imputation variance
  squ_fun <- function(x){
    x^2
  }
  
  temp <- lapply(se_within_list, squ_fun)
  arr_within <- array(unlist(temp) , c(nrow(se_within_list[[1]]),ncol(se_within_list[[1]]),m) )
  within_variance <- apply(arr_within, 1:2, mean)
  colnames(within_variance) <- colnames(se_within_list[[1]])
  rownames(within_variance) <- rownames(se_within_list[[1]])
  within_variance
  
  #calculate between imputation variance
  sub_fun <- function(x){
    x - pooled_estimate
  }
  
  delta_estimate     <- lapply(estimate_list, sub_fun)
  delta_estimate_squ <- lapply(delta_estimate, squ_fun)
  arr_between        <- array(unlist(delta_estimate_squ) , c(nrow(estimate_list[[1]]),ncol(estimate_list[[1]]),m) )
  
  spec_fun <- function(x){
    sum(x)/(m-1)
  }
  between_variance <- apply(arr_between, 1:2, spec_fun)
  colnames(between_variance) <- colnames(estimate_list[[1]])
  rownames(between_variance) <- rownames(estimate_list[[1]])
  between_variance
  
  #Then calculate total variance and pooled standard error
  total_variance <- within_variance + between_variance + I(between_variance/m)
  pooled_se      <- sqrt(total_variance)
  
  mycoef           <- cbind(pooled_estimate, pooled_se)
  colnames(mycoef) <- c("pooled_estimate", "pooled_se")
  
  return(mycoef)
  
}


##############################################################
#
# extract coefficients from a model object
#
##############################################################

coef_fun <- function(x){
  coefs <- summary(x)$coefficients
  return(coefs)
}


 

############################################
# 
# count NAs per columns and display
# those columns which exceed a percentage 
# threshold of missingness
#
############################################

count_na <- function(x, t) {
  na_count <- sapply(x, function(y) sum(length(which(is.na(y))))) %>% as.data.frame()
  colnames(na_count) <- "count"
  na_count$perc <- na_count$count/nrow(x)*100
  na_count <-
    na_count %>%
    filter(perc >= t)
  return(na_count)
}

##################################################
#
# Add real ROI names to lasso coefficients
#
#################################################

get_roi_names <- function(lasso, rois) {
  lasso$Metric <- ifelse(str_detect(lasso$X, "FA"), "FA", 
                         ifelse(str_detect(lasso$X, "MD"), "MD", NA))
  lasso$ROI_num <- gsub("stan_adj_FA_", "", lasso$X)
  lasso$ROI_num <- gsub("stan_adj_MD_", "", lasso$ROI_num)
  
  lasso <- merge(lasso, rois, by = "ROI_num", all.x = TRUE, all.y = FALSE )
  lasso$Coefficient <- 
    format(round(lasso$Recovery, 3), 3) %>% 
    as.numeric()
  lasso$Atlas <- ifelse(lasso$ROI_num == "(Intercept)", "(Model intercept)", as.character(lasso$Atlas))
  lasso <-
    lasso %>%
    select(Atlas, 
           Region = Name,
           Side,
           Metric,
           Coefficient)
  
  lasso$Atlas <- factor(lasso$Atlas,
                        levels = c("(Model intercept)",
                                   "JHU",
                                   "Extension WM",
                                   "Extension Thalamus",
                                   "Extension Brainstem"))
  lasso <-
    lasso %>%
    arrange(Atlas, Region, Side, Metric)
  
  return(lasso)
}


#######################
# initiate nested list
#######################

rec.list <- function(len){
  if(length(len) == 1){
    vector("list", len)
  } else {
    lapply(1:len[1], function(...) rec.list(len[-1]))
  }
}



#############################################
# make lasso summary table
# i.e. which tracts where selected how often
#############################################

make_lasso_summary_table <- function(imp_boot_lasso_list, lasso_result_orig_list) {
  #collate coefficients from all lasso models 
  #from all bootstrap samples 
  #within all m imputed datasets
  lasso_coef_list <- vector(mode = "list", length = m)
  
  for (i in 1:m){
    lasso_coef_list[[i]] <- lapply(imp_boot_lasso_list[[i]], `[[`, 2)
  }
  
  lasso_coef_table <- bind_cols(lasso_coef_list)
  lasso_coef_table$Count <- rowSums(lasso_coef_table!=0)
  lasso_coef_table$Tract <- rownames(lasso_coef_table)
  lasso_coef_table <-
    lasso_coef_table %>%
    select(Tract, Count) %>%
    filter(Tract != "(Intercept)")
  
  #identify tracts selected in original data (in at least 1 of m imputed datasets)
  lasso_tract_list <- vector(mode = "list", length = m)
  for (i in 1:m){
    lasso_tract_list <- lapply(lasso_result_orig_list, `[[`, 2)
  }
  lasso_tract_table <- bind_cols(lasso_tract_list)
  
  #select only tracts that were selected in at least half m imputed datasets
  lasso_tract_table$Count <- rowSums(lasso_tract_table != 0)
  t <- m/2
  lasso_tract_table <-
    lasso_tract_table %>%
    filter(Count >= t)
  
  lasso_coef_table <-
    lasso_coef_table %>%
    filter(Tract %in% rownames(lasso_tract_table))
  
  #replace tract abbreviations with real names
  lasso_coef_table$Tract <- gsub("stan_adj_", "", lasso_coef_table$Tract)
  lasso_coef_table$Tract <- gsub("_R", " R", lasso_coef_table$Tract)
  lasso_coef_table <- separate(lasso_coef_table, Tract, into = c("Metric", "ROI"), sep = " ")
  
  rois <- read.csv("/Users/sophierichter/Documents/Medicine/PhD/Trajectory_paper/Raw_data/JHU_updated/Atlas_Oct_2021_JHU_Extra_FINAL_Sophie.csv")
  rois$ROI_num <- paste0("ROI_", rois$ROI_num)
  rois <- 
    rois %>%
    select(ROI_num, Side, Name)
  
  lasso_coef_table <- merge(lasso_coef_table, rois, by.x = "ROI", by.y = "ROI_num", all.x = TRUE, all.y = FALSE)
  
  #display the frequency with which each of the tracts 
  #selected in the original data
  #was also selected in mxB bootstrap samples
  lasso_coef_table <- 
    lasso_coef_table %>%
    select(Name, Side, Metric, Count)
  lasso_coef_table$Percent <-
    (lasso_coef_table$Count/(m*B))*100 
  lasso_coef_table$Percent <- round(lasso_coef_table$Percent, 0)
  lasso_coef_table$Frequency <-
    paste0(lasso_coef_table$Count, " (", lasso_coef_table$Percent, "%)")
  lasso_coef_table <-
    lasso_coef_table %>%
    arrange(Name) %>%
    select(Tract = Name,
           Side,
           Metric,
           Frequency)
  
  #make a pretty flextable
  ft_lasso <- 
    lasso_coef_table %>%
    regulartable()
  
  save_as_docx(ft_lasso, path = paste0("Figures/lasso_table.docx"))
  
}


######################################
# Helper function to make
# biomarker performance table
# for one timepoint
######################################

format_tab <- function(tab, mysource, myvar){
  n <- colSums(!is.na(mysource))[[myvar]]
  colnames(tab) <- c("est", "ll", "ul")
  
  
  fun <- function(x){
    format(round(x,2), nsmall = 2)
  }
  tab <- apply(tab, 2, fun)
  tab <- as.data.frame(tab)
  tab$res <- paste0(tab$est, " (", tab$ll, "-", tab$ul, ")")
  tab$res <- gsub(" (  NA-  NA)", "", tab$res, fixed = TRUE)
  tab$res <- gsub(" (  NA-   NA)", "", tab$res, fixed = TRUE)
  tab$perc <- round((as.numeric(tab$res)/n *100),0)
  
  #calculated avoided scans
  x <- sum(mysource[, myvar] < as.numeric(tab[["cutoff", "est"]]), na.rm = TRUE)
  y <- round(x/n*100, 0)
  avoided <- paste0(x, " (", y, "%)")
  
  tab[["FN", "res"]] <- paste0(round(as.numeric(tab[["FN", "res"]]), 0), " (", tab[["FN", "perc"]], "%)")
  tab[["FP", "res"]] <- paste0(round(as.numeric(tab[["FP", "res"]]), 0), " (", tab[["FP", "perc"]], "%)")
  tab <- tab %>% select(res)
  tab <- t(tab) %>% as.data.frame()
  tab$Avoided <- avoided
  tab$N <- n
  tab <-
    tab %>%
    select(N, 
           cutoff,
           Se,
           Sp,
           Avoided,
           FP,
           FN)
  
  return(tab)
}


##########################################
# Helper function
# to turn list elements into dataframes
##########################################

df_fun <- function(x){
  t <- as.data.frame(x)
  t$Variable <- rownames(t)
  return(t)
}


#############################################################
# display model pooled coefficients for
# all models with and without DTI
# fitted to the original (imputed but not boottsrapped) data
# and save as pretty table
#############################################################
save_model_coefficients <- function(orig_coef_list){
  
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
           levels = c("DTI_only",
                      "plain_headsmart",
                      "plain_upfrontEM",
                      "plain_upfrontPLUS",
                      "plain_centerEM",
                      "plain_centerPLUS",
                      "dti_headsmart",
                      "dti_upfrontEM",
                      "dti_upfrontPLUS",
                      "dti_centerEM",
                      "dti_centerPLUS"),
           labels = c("DTI only",
                      "HeadSMART",
                      "UPFRONT-ED",
                      "UPFRONT-PLUS",
                      "CENTER-ED",
                      "CENTER-PLUS",
                      "HeadSMART with DTI",
                      "UPFRONT-ED with DTI",
                      "UPFRONT-PLUS with DTI",
                      "CENTER-ED with DTI",
                      "CENTER-PLUS with DTI"))
  
  
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
  save_as_docx(ft_model_coefs, path = paste0("Figures/modelcoefs.docx"))
  
}




