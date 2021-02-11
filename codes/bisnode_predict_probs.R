#########################################################################################
# Prepared for Gabor's Data Analysis
#
# Data Analysis for Business, Economics, and Policy
# by Gabor Bekes and  Gabor Kezdi
# Cambridge University Press 2021
#
# gabors-data-analysis.com 
#
# License: Free to share, modify and use for educational purposes. 
# 	Not to be used for commercial purposes.


# This code is co-authored with Zsuzsa Holler and Jeno Pal
# Chapter 17
# CH17A
# using the bisnode-firmd dataset
# version 0.9 2020-09-10
#########################################################################################



# ------------------------------------------------------------------------------------------------------
#### SET UP
# CLEAR MEMORY
rm(list=ls())

# Import libraries
library(haven)
library(glmnet)
library(purrr)
library(margins)
library(skimr)
library(kableExtra)
library(Hmisc)
library(cowplot)
library(gmodels) 
library(lspline)
library(sandwich)
library(modelsummary)

library(rattle)
library(caret)
library(pROC)
library(ranger)
library(rpart)
library(partykit)
library(rpart.plot)
library(viridis)
library(xgboost)
library(tidyverse)




# Loading and preparing data ----------------------------------------------

# Use R format so it keeps factor definitions

data <- readRDS("C:/Users/T450s/Desktop/programming/git/business_growth_prediction/data/clean/bisnode_firms_clean.rds")
data <- data.frame(data)

#get helper functions
source('https://raw.githubusercontent.com/zsomborh/business_growth_prediction/main/codes/helper.R')

#summary
datasummary_skim(data, type='numeric', histogram = TRUE)
# datasummary_skim(data, type="categorical")


# Define variable sets ----------------------------------------------
# (making sure we use ind2_cat, which is a factor)

rawvars <-  c("curr_assets", "curr_liab", "extra_exp", "extra_inc", "extra_profit_loss", "fixed_assets",
              "inc_bef_tax", "intang_assets", "inventories", "liq_assets", "material_exp", "personnel_exp",
              "profit_loss_year", "sales", "share_eq", "subscribed_cap")
qualityvars <- c("balsheet_flag", "balsheet_length", "balsheet_notfullyear")
engvar <- c("total_assets_bs", "fixed_assets_bs", "liq_assets_bs", "curr_assets_bs",
            "share_eq_bs", "subscribed_cap_bs", "intang_assets_bs", "extra_exp_pl",
            "extra_inc_pl", "extra_profit_loss_pl", "inc_bef_tax_pl", "inventories_pl",
            "material_exp_pl", "profit_loss_year_pl", "personnel_exp_pl")
engvar2 <- c("extra_profit_loss_pl_quad", "inc_bef_tax_pl_quad",
             "profit_loss_year_pl_quad", "share_eq_bs_quad", 'labor_avg_mod_sq')
engvar3 <- c(grep("*flag_low$", names(data), value = TRUE),
             grep("*flag_high$", names(data), value = TRUE),
             grep("*flag_error$", names(data), value = TRUE),
             grep("*flag_zero$", names(data), value = TRUE))
d1 <-  c("d1_sales_mil_log_mod", "d1_sales_mil_log_mod_sq",
         "flag_low_d1_sales_mil_log", "flag_high_d1_sales_mil_log")
hr <- c("female", "ceo_age", "flag_high_ceo_age", "flag_low_ceo_age",
        "flag_miss_ceo_age", "ceo_count", "labor_avg_mod",
        "flag_miss_labor_avg", "foreign_management")
firm <- c("age", "age2", "new", "ind2_cat", "m_region_loc", "urban_m")

# interactions for logit, LASSO
interactions1 <- c("ind2_cat*age", "ind2_cat*age2",
                   "ind2_cat*d1_sales_mil_log_mod", "ind2_cat*sales_mil_log",
                   "ind2_cat*ceo_age", "ind2_cat*foreign_management",
                   "ind2_cat*female",   "ind2_cat*urban_m", "ind2_cat*labor_avg_mod")
interactions2 <- c("sales_mil_log*age", "sales_mil_log*female",
                   "sales_mil_log*profit_loss_year_pl", "sales_mil_log*foreign_management")


X1 <- c("sales_mil_log", "sales_mil_log_sq", "d1_sales_mil_log_mod", "profit_loss_year_pl", "ind2_cat")
X2 <- c("sales_mil_log", "sales_mil_log_sq", "d1_sales_mil_log_mod", "profit_loss_year_pl", "fixed_assets_bs","share_eq_bs","curr_liab_bs",   "curr_liab_bs_flag_high", "curr_liab_bs_flag_error",  "age","foreign_management" , "ind2_cat")
X3 <- c("sales_mil_log", "sales_mil_log_sq", firm, engvar,                   d1)
X4 <- c("sales_mil_log", "sales_mil_log_sq", firm, engvar, engvar2, engvar3, d1, hr, qualityvars)
X5 <- c("sales_mil_log", "sales_mil_log_sq", firm, engvar, engvar2, engvar3, d1, hr, qualityvars, interactions1, interactions2)




# for LASSO
logitvars <- c("sales_mil_log", "sales_mil_log_sq", engvar, engvar2, engvar3, d1, hr, firm, qualityvars, interactions1, interactions2)

# for RF (no interactions, no modified features)
rfvars  <-  c("sales_mil", "d1_sales_mil_log", rawvars, hr, firm, qualityvars)


# Check simplest model X1
ols_modelx1 <- lm(formula(paste0("fast_g ~", paste0(X1, collapse = " + "))),
                  data = data)
summary(ols_modelx1)

glm_modelx1 <- glm(formula(paste0("fast_g ~", paste0(X1, collapse = " + "))),
                   data = data, family = "binomial")
summary(glm_modelx1)


# Check model X2
glm_modelx2 <- glm(formula(paste0("fast_g ~", paste0(X2, collapse = " + "))),
                   data = data, family = "binomial")
summary(glm_modelx2)

#calculate average marginal effects (dy/dx) for logit
mx2 <- margins(glm_modelx2)

sum_table <- summary(glm_modelx2) %>%
    coef() %>%
    as.data.frame() %>%
    select(Estimate) %>%
    mutate(factor = row.names(.)) %>%
    merge(summary(mx2)[,c("factor","AME")])

kable(x = sum_table, format = "latex", digits = 3,
      col.names = c("Variable", "Coefficient", "dx/dy"),
      caption = "Average Marginal Effects (dy/dx) for Logit Model") %>%
    cat(.,file= paste0(output,"AME_logit_X2.tex"))

sum_table

summary(glm_modelx2) %>%
    coef() %>%
    as.data.frame() %>%  select(Estimate)
# baseline model is X4 (all vars, but no interactions) -------------------------------------------------------

ols_model <- lm(formula(paste0("fast_g ~", paste0(X4, collapse = " + "))),
                data = data)
summary(ols_model)

glm_model <- glm(formula(paste0("fast_g ~", paste0(X3, collapse = " + "))),
                 data = data, family = "binomial")
summary(glm_model)

#calculate average marginal effects (dy/dx) for logit
# vce="none" makes it run much faster, here we do not need variances

m <- margins(glm_model, vce = "none")

sum_table2 <- summary(glm_model) %>%
    coef() %>%
    as.data.frame() %>%
    select(Estimate, `Std. Error`) %>%
    mutate(factor = row.names(.)) %>%
    merge(summary(m)[,c("factor","AME")])

kable(x = sum_table2, format = "latex", digits = 3,
      col.names = c("Variable", "Coefficient", "SE", "dx/dy"),
      caption = "Average Marginal Effects (dy/dx) for Logit Model") 

saveRDS(sum_table2, 'logit_3_summary.rds')
# separate datasets -------------------------------------------------------

set.seed(8)

train_indices <- as.integer(createDataPartition(data$fast_g, p = 0.8, list = FALSE))
data_train <- data[train_indices, ]
data_holdout <- data[-train_indices, ]

dim(data_train)
dim(data_holdout)

Hmisc::describe(data$fast_g_f)
Hmisc::describe(data_train$fast_g_f)
Hmisc::describe(data_holdout
                $fast_g_f)

glm_model <- glm(formula(paste0("fast_g ~", paste0(X3, collapse = " + "))),
                 data = data_train, family = "binomial")
summary(glm_model)
#######################################################x
# PART I PREDICT PROBABILITIES
# Predict logit models ----------------------------------------------
#######################################################x

# 5 fold cross-validation

train_control <- trainControl(
    method = "cv",
    number = 5,
    classProbs = TRUE,
    summaryFunction = twoClassSummaryExtended,
    savePredictions = TRUE
)

# Train Logit Models ----------------------------------------------

logit_model_vars <- list("X1" = X1, "X2" = X2, "X3" = X3, "X4" = X4, "X5" = X5)

CV_RMSE_folds <- list()
logit_models <- list()


for (model_name in names(logit_model_vars)) {
    
    features <- logit_model_vars[[model_name]]
    
    set.seed(7)
    glm_model <- train(
        formula(paste0("fast_g_f ~", paste0(features, collapse = " + "))),
        method = "glm",
        data = data_train,
        family = binomial,
        trControl = train_control
    )
    
    logit_models[[model_name]] <- glm_model
    # Calculate RMSE on test for each fold
    CV_RMSE_folds[[model_name]] <- glm_model$resample[,c("Resample", "RMSE")]
    
}

CV_RMSE_folds
# Logit lasso -----------------------------------------------------------

lambda <- 10^seq(-1, -4, length = 10)
grid <- expand.grid("alpha" = 1, lambda = lambda)

set.seed(9)
system.time({
    logit_lasso_model <- train(
        formula(paste0("fast_g_f ~", paste0(logitvars, collapse = " + "))),
        data = data_train,
        method = "glmnet",
        preProcess = c("center", "scale"),
        family = "binomial",
        trControl = train_control,
        tuneGrid = grid,
        na.action=na.exclude
    )
})

tuned_logit_lasso_model <- logit_lasso_model$finalModel
best_lambda <- logit_lasso_model$bestTune$lambda
logit_models[["LASSO"]] <- logit_lasso_model
lasso_coeffs <- as.matrix(coef(tuned_logit_lasso_model, best_lambda))

CV_RMSE_folds[["LASSO"]] <- logit_lasso_model$resample[,c("Resample", "RMSE")]


#############################################x
# PART I
# No loss fn
########################################

# Draw ROC Curve and calculate AUC for each folds --------------------------------
CV_AUC_folds <- list()

for (model_name in names(logit_models)) {
    
    auc <- list()
    model <- logit_models[[model_name]]
    for (fold in c("Fold1", "Fold2", "Fold3", "Fold4", "Fold5")) {
        cv_fold <-
            model$pred %>%
            filter(Resample == fold)
        
        roc_obj <- roc(cv_fold$obs, cv_fold$fast)
        auc[[fold]] <- as.numeric(roc_obj$auc)
    }
    
    CV_AUC_folds[[model_name]] <- data.frame("Resample" = names(auc),
                                             "AUC" = unlist(auc))
}


# For each model: average RMSE and average AUC for models ----------------------------------

CV_RMSE <- list()
CV_AUC <- list()

for (model_name in names(logit_models)) {
    CV_RMSE[[model_name]] <- mean(CV_RMSE_folds[[model_name]]$RMSE)
    CV_AUC[[model_name]] <- mean(CV_AUC_folds[[model_name]]$AUC)
}

CV_AUC
# We have 6 models, (5 logit and the logit lasso). For each we have a 5-CV RMSE and AUC.
# We pick our preferred model based on that. -----------------------------------------------

nvars <- lapply(logit_models, FUN = function(x) length(x$coefnames))
nvars[["LASSO"]] <- sum(lasso_coeffs != 0)

logit_summary1 <- data.frame("Number of predictors" = unlist(nvars),
                             "CV RMSE" = unlist(CV_RMSE),
                             "CV AUC" = unlist(CV_AUC))


saveRDS(logit_summary1, 'logit_summary.rds')
    
# Take best model and estimate RMSE on holdout  -------------------------------------------

best_logit_no_loss <- logit_models[["X3"]]

logit_predicted_probabilities_holdout <- predict(best_logit_no_loss, newdata = data_holdout, type = "prob")
data_holdout[,"best_logit_no_loss_pred"] <- logit_predicted_probabilities_holdout[,"fast"]
RMSE(data_holdout[, "best_logit_no_loss_pred", drop=TRUE], data_holdout$fast_g)

# discrete ROC (with thresholds in steps) on holdout -------------------------------------------------
thresholds <- seq(0.05, 0.75, by = 0.05)

cm <- list()
true_positive_rates <- c()
false_positive_rates <- c()
for (thr in thresholds) {
    holdout_prediction <- ifelse(data_holdout[,"best_logit_no_loss_pred"] < thr, "not_fast", "fast") %>%
        factor(levels = c("not_fast", "fast"))
    cm_thr <- confusionMatrix(holdout_prediction,data_holdout$fast_g_f)$table
    cm[[as.character(thr)]] <- cm_thr
    true_positive_rates <- c(true_positive_rates, cm_thr["fast", "fast"] /
                                 (cm_thr["fast", "fast"] + cm_thr["not_fast", "fast"]))
    false_positive_rates <- c(false_positive_rates, cm_thr["fast", "not_fast"] /
                                  (cm_thr["fast", "not_fast"] + cm_thr["not_fast", "not_fast"]))
}



tpr_fpr_for_thresholds <- tibble(
    "threshold" = thresholds,
    "true_positive_rate" = true_positive_rates,
    "false_positive_rate" = false_positive_rates
)
discrete_roc_plot <- ggplot(
    data = tpr_fpr_for_thresholds,
    aes(x = false_positive_rate, y = true_positive_rate, color = threshold)) +
    labs(x = "False positive rate (1 - Specificity)", y = "True positive rate (Sensitivity)") +
    geom_point(size=2, alpha=0.8) +
    scale_color_viridis(option = "D", direction = -1) +
    scale_x_continuous(expand = c(0.01,0.01), limit=c(0,1), breaks = seq(0,1,0.1)) +
    scale_y_continuous(expand = c(0.01,0.01), limit=c(0,1), breaks = seq(0,1,0.1)) +
    theme_minimal() +
    theme(legend.position ="right") +
    theme(legend.title = element_text(size = 4), 
          legend.text = element_text(size = 4),
          legend.key.size = unit(.4, "cm")) 
discrete_roc_plot


# continuous ROC on holdout with best model (Logit 3) -------------------------------------------

roc_obj_holdout <- roc(data_holdout$fast_g, data_holdout$best_logit_no_loss_pred)

roc_holdout_logit <- createRocPlot(roc_obj_holdout, "best_logit_no_loss_roc_plot_holdout")

saveRDS(roc_holdout_logit, 'roc_holdout_logit.rds')

# Confusion table with different tresholds ----------------------------------------------------------

# default: the threshold 0.5 is used to convert probabilities to binary classes
logit_class_prediction <- predict(best_logit_no_loss, newdata = data_holdout)
summary(logit_class_prediction)

# confusion matrix: summarize different type of errors and successfully predicted cases
# positive = "yes": explicitly specify the positive case
cm_object1 <- confusionMatrix(logit_class_prediction, data_holdout$fast_g_f, positive = "fast")
cm1 <- cm_object1$table
cm1

# we can apply different thresholds

# 0.5 same as before
holdout_prediction <-
    ifelse(data_holdout$best_logit_no_loss_pred < 0.5, "not_fast", "fast") %>%
    factor(levels = c("not_fast", "fast"))
cm_object1b <- confusionMatrix(holdout_prediction,data_holdout$fast_g_f)
cm1b <- cm_object1b$table
cm1b

# a sensible choice: mean of predicted probabilities

mean_predicted_default_prob <- mean(data_holdout$best_logit_no_loss_pred)
mean_predicted_default_prob
holdout_prediction <-
    ifelse(data_holdout$best_logit_no_loss_pred < mean_predicted_default_prob, "not_fast", "fast") %>%
    factor(levels = c("not_fast", "fast"))
cm_object2 <- confusionMatrix(holdout_prediction,data_holdout$fast_g_f)
cm2 <- cm_object2$table
cm2



# Calibration curve -----------------------------------------------------------
# how well do estimated vs actual event probabilities relate to each other?


logit_cal_plot <- create_calibration_plot(data_holdout,
                        prob_var = "best_logit_no_loss_pred", 
                        actual_var = "fast_g",
                        n_bins = 10)

saveRDS(logit_cal_plot , 'logit_cal_plot.rds')
#############################################x
# PART II.
# We have a loss function
########################################

# Introduce loss function
# relative cost of of a false negative classification (as compared with a false positive classification)


# Initial invesment is 10m - investment horizon is 2 years
# Cost from False positive - missing out on the investment - this will be equal to the fast-growth threshold specified in cleaner

sp <- read_csv('C:/Users/T450s/Desktop/programming/git/business_growth_prediction/data/raw/sp-500-historical-annual-returns.csv', skip = 15)
sp$year <- substr(sp$date,0,4)
expected_earning <- (mean(sp[sp$year >= 2009,]$value)/ 100)*2
two_year_exp_earn <-(expected_earning * (1+expected_earning)) 

# Cost of false negative - average sales growth for 
worst_80_perc <- quantile(data[data$fast_g == 0,]$g, .20, na.rm=T)
worst_80_perc
FP=10 * two_year_exp_earn
FN=-10 * worst_80_perc[[1]] 

cost = FN/FP
# the prevalence, or the proportion of cases in the population (n.cases/(n.controls+n.cases))
prevelance = sum(data_train$fast_g)/length(data_train$fast_g)

# Draw ROC Curve and find optimal threshold with loss function --------------------------

best_tresholds <- list()
expected_loss <- list()
logit_cv_rocs <- list()
logit_cv_threshold <- list()
logit_cv_expected_loss <- list()

for (model_name in names(logit_models)) {
    
    model <- logit_models[[model_name]]
    colname <- paste0(model_name,"_prediction")
    
    best_tresholds_cv <- list()
    expected_loss_cv <- list()
    
    for (fold in c("Fold1", "Fold2", "Fold3", "Fold4", "Fold5")) {
        cv_fold <-
            model$pred %>%
            filter(Resample == fold)
        
        roc_obj <- roc(cv_fold$obs, cv_fold$fast)
        best_treshold <- coords(roc_obj, "best", ret="all", transpose = FALSE,
                                best.method="youden", best.weights=c(cost, prevelance))
        best_tresholds_cv[[fold]] <- best_treshold$threshold
        expected_loss_cv[[fold]] <- (best_treshold$fp*FP + best_treshold$fn*FN)/length(cv_fold$fast)
    }
    
    # average
    best_tresholds[[model_name]] <- mean(unlist(best_tresholds_cv))
    expected_loss[[model_name]] <- mean(unlist(expected_loss_cv))
    
    # for fold #5
    logit_cv_rocs[[model_name]] <- roc_obj
    logit_cv_threshold[[model_name]] <- best_treshold
    logit_cv_expected_loss[[model_name]] <- expected_loss_cv[[fold]]
    
}

logit_summary2 <- data.frame("Avg of optimal thresholds" = unlist(best_tresholds),
                             "Threshold for Fold5" = sapply(logit_cv_threshold, function(x) {x$threshold}),
                             "Avg expected loss" = unlist(expected_loss),
                             "Expected loss for Fold5" = unlist(logit_cv_expected_loss))

# Create plots based on Fold5 in CV ----------------------------------------------

for (model_name in names(logit_cv_rocs)) {
    
    r <- logit_cv_rocs[[model_name]]
    best_coords <- logit_cv_threshold[[model_name]]
    createLossPlot(r, best_coords,
                   paste0(model_name, "_loss_plot"))
    createRocPlotWithOptimal(r, best_coords,
                             paste0(model_name, "_roc_plot"))
}

# Pick best model based on average expected loss ----------------------------------

best_logit_with_loss <- logit_models[["X4"]]
best_logit_optimal_treshold <- best_tresholds[["X4"]]

logit_predicted_probabilities_holdout <- predict(best_logit_with_loss, newdata = data_holdout, type = "prob")
data_holdout[,"best_logit_with_loss_pred"] <- logit_predicted_probabilities_holdout[,"fast"]

# ROC curve on holdout
roc_obj_holdout <- roc(data_holdout$fast_g, data_holdout[, "best_logit_with_loss_pred", drop=TRUE])

# Get expected loss on holdout
holdout_treshold <- coords(roc_obj_holdout, x = best_logit_optimal_treshold, input= "threshold",
                           ret="all", transpose = FALSE)
expected_loss_holdout <- (holdout_treshold$fp*FP + holdout_treshold$fn*FN)/length(data_holdout$fast_g)
expected_loss_holdout

# Confusion table on holdout with optimal threshold
holdout_prediction <-
    ifelse(data_holdout$best_logit_with_loss_pred < best_logit_optimal_treshold, "not_fast", "fast") %>%
    factor(levels = c("not_fast", "fast"))
cm_object3 <- confusionMatrix(holdout_prediction,data_holdout$fast_g_f)
cm3 <- cm_object3$table


#################################################
# PREDICTION WITH RANDOM FOREST
#################################################

# -----------------------------------------------
# RANDOM FOREST GRAPH EXAMPLE
# -----------------------------------------------

data_for_graph <- data_train
levels(data_for_graph$fast_g_f) <- list("not_fast" = "not_fast", "fast" = "fast")


set.seed(8)
rf_for_graph <-
    rpart(
        formula = fast_g_f ~ sales_mil + profit_loss_year+ foreign_management,
        data = data_for_graph,
        control = rpart.control(cp = 0.00010, minbucket = 87)
    )

rpart.plot(rf_for_graph, tweak=1, digits=2, extra=107, under = TRUE)


#################################################
# Probability forest
# Split by gini, ratio of 1's in each tree, average over trees
#################################################

# 5 fold cross-validation

train_control <- trainControl(
    method = "cv",
    n = 5,
    classProbs = TRUE, # same as probability = TRUE in ranger
    summaryFunction = twoClassSummaryExtended,
    savePredictions = TRUE,
    verboseIter = TRUE
)


tune_grid <- expand.grid(
    .mtry = c(5, 6, 7),
    .splitrule = "gini",
    .min.node.size = c(10, 15)
)

# getModelInfo("ranger")
set.seed(7)
rf_model_p <- train(
    formula(paste0("fast_g_f ~ ", paste0(rfvars , collapse = " + "))),
    method = "ranger",
    data = data_train,
    tuneGrid = tune_grid,
    trControl = train_control,
    num.threads = 2
)


best_mtry <- rf_model_p$bestTune$mtry
best_min_node_size <- rf_model_p$bestTune$min.node.size


# Get average (ie over the folds) RMSE and AUC ------------------------------------
CV_RMSE_folds[["rf_p"]] <- rf_model_p$resample[,c("Resample", "RMSE")]

auc <- list()
for (fold in c("Fold1", "Fold2", "Fold3", "Fold4", "Fold5")) {
    cv_fold <-
        rf_model_p$pred %>%
        filter(Resample == fold)
    
    roc_obj <- roc(cv_fold$obs, cv_fold$fast)
    auc[[fold]] <- as.numeric(roc_obj$auc)
}
CV_AUC_folds[["rf_p"]] <- data.frame("Resample" = names(auc),
                                     "AUC" = unlist(auc))

CV_RMSE[["rf_p"]] <- mean(CV_RMSE_folds[["rf_p"]]$RMSE)
CV_AUC[["rf_p"]] <- mean(CV_AUC_folds[["rf_p"]]$AUC)

CV_RMSE


# CREATE AND SAVE ROC PLOTS 

rf_pred_prob_holdout <- predict(rf_model, newdata = data_holdout, type = "prob")
data_holdout[,"rf_pred1"] <- logit_predicted_probabilities_holdout[,"fast"]

roc_obj_holdout <- roc(data_holdout$fast_g, data_holdout$rf_pred1)
roc_holdout_rf <- createRocPlot(roc_obj_holdout, "rf_roc_plot_holdout")
roc_holdout_rf
saveRDS(roc_holdout_rf, 'roc_holdout_rf.rds')

rf_cal_plot <- create_calibration_plot(data_holdout,
                                          prob_var = "rf_pred1", 
                                          actual_var = "fast_g",
                                          n_bins = 10)
rf_cal_plot

saveRDS(rf_cal_plot , 'rf_cal_plot.rds')

# Naive Approach - threshold is mean of predictions 

mean_predicted_default_prob <- mean(data_holdout$rf_pred1)

holdout_prediction <-
    ifelse(data_holdout$rf_pred1 < mean_predicted_default_prob, "not_fast", "fast") %>%
    factor(levels = c("not_fast", "fast"))
cm_object2 <- confusionMatrix(holdout_prediction,data_holdout$fast_g_f)

saveRDS(cm_object2,'naive_confusion.rds')

# Now use loss function and search for best thresholds and expected loss over folds -----
best_tresholds_cv <- list()
expected_loss_cv <- list()

for (fold in c("Fold1", "Fold2", "Fold3", "Fold4", "Fold5")) {
    cv_fold <-
        rf_model_p$pred %>%
        filter(mtry == best_mtry,
               min.node.size == best_min_node_size,
               Resample == fold)
    
    roc_obj <- roc(cv_fold$obs, cv_fold$fast)
    best_treshold <- coords(roc_obj, "best", ret="all", transpose = FALSE,
                            best.method="youden", best.weights=c(cost, prevelance))
    best_tresholds_cv[[fold]] <- best_treshold$threshold
    expected_loss_cv[[fold]] <- (best_treshold$fp*FP + best_treshold$fn*FN)/length(cv_fold$fast)
}

# average
best_tresholds[["rf_p"]] <- mean(unlist(best_tresholds_cv))
expected_loss[["rf_p"]] <- mean(unlist(expected_loss_cv))


rf_summary <- data.frame("CV RMSE" = CV_RMSE[["rf_p"]],
                         "CV AUC" = CV_AUC[["rf_p"]],
                         "Avg of optimal thresholds" = best_tresholds[["rf_p"]],
                         "Threshold for Fold5" = best_treshold$threshold,
                         "Avg expected loss" = expected_loss[["rf_p"]],
                         "Expected loss for Fold5" = expected_loss_cv[[fold]])



# Create plots - this is for Fold5

best_loss <- createLossPlot(roc_obj, best_treshold, "rf_p_loss_plot")
createRocPlotWithOptimal(roc_obj, best_treshold, "rf_p_roc_plot")
saveRDS(best_loss, 'best_loss.rds')
# Create confusion matrix 


holdout_prediction_best <-
    ifelse(data_holdout$rf_pred1 < best_treshold$threshold, "not_fast", "fast") %>%
    factor(levels = c("not_fast", "fast"))

cm_object3<- confusionMatrix(holdout_prediction_best,data_holdout$fast_g_f)

saveRDS(cm_object3,'smart_confusion.rds')

# Take model to holdout and estimate RMSE, AUC and expected loss ------------------------------------


rf_predicted_probabilities_holdout <- predict(rf_model_p, newdata = data_holdout, type = "prob")
data_holdout$rf_p_prediction <- rf_predicted_probabilities_holdout[,"fast"]
RMSE(data_holdout$rf_p_prediction, data_holdout$fast_g)

# Estimate losses 


loss_logit <- as.data.frame(as.matrix(cm_object1))[1,'fast'] * FN +
    as.data.frame(as.matrix(cm_object1))[2,'not_fast'] * FP


loss_naive <- as.data.frame(as.matrix(cm_object2))[1,'fast'] * FN +
    as.data.frame(as.matrix(cm_object2))[2,'not_fast'] * FP


loss_smart <- as.data.frame(as.matrix(cm_object3))[1,'fast'] * FN +
    as.data.frame(as.matrix(cm_object3))[2,'not_fast'] * FP

loss_naive

loss_smart/loss_logit

# ROC curve on holdout
roc_obj_holdout <- roc(data_holdout$fast_g, data_holdout[, "rf_p_prediction", drop=TRUE])
roc_obj_holdout
# AUC
as.numeric(roc_obj_holdout$auc)

# Get expected loss on holdout with optimal threshold
holdout_treshold <- coords(roc_obj_holdout, x = best_tresholds[["rf_p"]] , input= "threshold",
                           ret="all", transpose = FALSE)
expected_loss_holdout <- (holdout_treshold$fp*FP + holdout_treshold$fn*FN)/length(data_holdout$fast_g)
expected_loss_holdout

#################################################
# Classification forest
# Split by Gini, majority vote in each tree, majority vote over trees
#################################################
# Show expected loss with classification RF and default majority voting to compare

train_control <- trainControl(
    method = "cv",
    n = 5
)
train_control$verboseIter <- TRUE

set.seed(7)
rf_model_f <- train(
    formula(paste0("fast_g_f ~ ", paste0(rfvars , collapse = " + "))),
    method = "ranger",
    data = data_train,
    tuneGrid = tune_grid,
    trControl = train_control
)

data_train$rf_f_prediction_class <-  predict(rf_model_f,type = "raw")
data_holdout$rf_f_prediction_class <- predict(rf_model_f, newdata = data_holdout, type = "raw")

#We use predicted classes to calculate expected loss based on our loss fn
fp <- sum(data_holdout$rf_f_prediction_class == "fast" & data_holdout$fast_g_f == "not_fast")
fn <- sum(data_holdout$rf_f_prediction_class == "not_fast" & data_holdout$fast_g_f == "fast")
(fp*FP + fn*FN)/length(data_holdout$fast_g)


# Summary results ---------------------------------------------------

nvars[["rf_p"]] <- length(rfvars)

summary_results <- data.frame("Number of predictors" = unlist(nvars),
                              "CV RMSE" = unlist(CV_RMSE),
                              "CV AUC" = unlist(CV_AUC),
                              "CV threshold" = unlist(best_tresholds),
                              "CV expected Loss" = unlist(expected_loss))

model_names <- c("Logit X1", "Logit X4",
                 "Logit LASSO","RF probability")
summary_results <- summary_results %>%
    filter(rownames(.) %in% c("X1", "X4", "LASSO", "rf_p"))
rownames(summary_results) <- model_names

kable(x = summary_results, format = "latex", booktabs=TRUE,  digits = 3, row.names = TRUE,
      linesep = "", col.names = c("Number of predictors", "CV RMSE", "CV AUC",
                                  "CV threshold", "CV expected Loss")) %>%
    cat(.,file= paste0(output, "summary_results.tex"))

summary_results


# 5 fold cross-validation

train_control <- trainControl(
    method = "cv",
    number = 5,
    classProbs = TRUE, # same as probability = TRUE in ranger
    summaryFunction = twoClassSummaryExtended,
    savePredictions = TRUE,
    verboseIter = TRUE, 
    allowParallel = TRUE
)


gbm_grid <-  expand.grid(nrounds=c(350),
                         max_depth = c(2, 4, 6),
                         eta = c(0.01,0.05, 0.1),
                         gamma = c(0.01,1),
                         colsample_bytree = c(0.75),
                         subsample = c(0.50),
                         min_child_weight = c(0))

gbm_grid


set.seed(7)
xgb_model <- train(
    formula(paste0("fast_g_f ~ ", paste0(rfvars , collapse = " + "))),
    method = "xgbTree",
    #metric = 'NormalizedGini',
    data = data_train,
    tuneGrid = gbm_grid,
    trControl = train_control,
    num.threads = 3
)


xgb_model <- readRDS('C:/Users/T450s/Desktop/programming/git/business_growth_prediction/models/xgb_model.rds')

CV_RMSE_folds[["xgb"]] <- xgb_model$resample[,c("Resample", "RMSE")]

auc <- list()
for (fold in c("Fold1", "Fold2", "Fold3", "Fold4", "Fold5")) {
    cv_fold <-
        xgb_model$pred %>%
        filter(Resample == fold)
    
    roc_obj <- roc(cv_fold$obs, cv_fold$fast)
    auc[[fold]] <- as.numeric(roc_obj$auc)
}
CV_AUC_folds[["xgb"]] <- data.frame("Resample" = names(auc),
                                     "AUC" = unlist(auc))

CV_RMSE[["xgb"]] <- mean(CV_RMSE_folds[["xgb"]]$RMSE)
CV_AUC[["xgb"]] <- mean(CV_AUC_folds[["xgb"]]$AUC)

nvars[["xgb"]] <- length(rfvars)

all_results <- data.frame("Number of predictors" = unlist(nvars),
           "CV RMSE" = unlist(CV_RMSE),
           "CV AUC" = unlist(CV_AUC))

all_results <- all_results[rownames(all_results) %in% c('X1', 'X3', 'rf_p', 'xgb'),]

saveRDS(logit_models, 'logit_models.rds')
saveRDS(logit_lasso_model, 'logit_lasso.rds')
saveRDS(rf_model_p, 'rf_model_p.rds')
saveRDS(rf_model_f, 'rf_model_f.rds')
saveRDS(xgb_model, 'xgb_model.rds')
saveRDS(all_results, 'all_results.rds')
