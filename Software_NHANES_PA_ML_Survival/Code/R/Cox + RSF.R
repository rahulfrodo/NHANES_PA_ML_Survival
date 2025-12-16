###########################
## Libraries
###########################

rm(list = ls())    
Sys.setenv(LANGUAGE = "en")
gc() # free unused memory

library(dplyr)
library(survival)
library(refund)
library(caret)
library(MASS)
library(grpreg)
library(mgcv)
library(randomForestSRC)
library(SurvMetrics)
library(ggplot2)
library(scales)
library(patchwork)
library(survminer)
library(ggpubr)

############################
## The DataSet
############################

load("NHANES_11_14_survPA_TD_Final.RData") 

##
library(haven)
library(foreign)
bmidfG<- read.xport("BMX_G.xpt")
bmidfH<- read.xport("BMX_H.xpt")
bmidfcomb<-rbind(bmidfG,bmidfH)
bmidfbase<-bmidfcomb[bmidfcomb$SEQN%in%baseline$SEQN,]
baseline$height<-bmidfbase$BMXHT

plot(baseline$height,baseline$RIDAGEYR,xlab="Height",ylab="Age" ,main="NHANES 2011-14 (older than 50 years)")
abline(h=75,col="2",lty=2)
abline(v=180,col="2",lty=2)

cor(baseline$height,baseline$RIDAGEYR,use="pairwise.complete.obs")
#-0.1126176 significant?
x<-baseline$height
y<-baseline$RIDAGEYR
cor(x, y, method = c("pearson"),use="pairwise.complete.obs")
cor.test(x, y, method=c("pearson"),use="pairwise.complete.obs")
#significant association# -0.1126176 



##
time_mort<-baseline$permth_exm/12 #age in month to year-unit
baseline$time_mort<-time_mort
baseline$Age<-baseline$RIDAGEYR
baseline$Race<-baseline$RIDRETH3

df <- baseline %>%
  dplyr::select(Age, Race, BMI, sex, Mobility, mortstat, diabetes.y, poverty_level, 
                Asthma, Arthritis, heart_failure, coronary_heart_disease, angina, 
                stroke, thyroid, bronchitis, cancer, time_mort, actmat) %>%
  filter(!is.na(mortstat) &
           !is.na(time_mort) &
           !is.na(BMI))

# Compute row means of activity matrix
row_means <- rowMeans(df$actmat, na.rm = TRUE)

# Attach it to the original dataframe
df$act_mean <- row_means

# Group by sex and death, then count occurrences
df_summary <- df %>%
  group_by(sex, mortstat) %>%
  summarise(count = n(), .groups = "drop")

print(df_summary)

mean(df$act_mean, na.rm = TRUE)
tapply(df$act_mean, df$mortstat, mean)

##############################
## KM analysis
##############################

# Create age groups
km_df <- df %>%
  mutate(age_group = cut(Age, 
                         breaks = c(50, 60, 70, Inf), 
                         labels = c("50-60", "60-70", "70-above"),
                         right = FALSE))


# KM plot by age group
fit_age <- survfit(Surv(time = km_df$time_mort, event = km_df$mortstat) ~ age_group, data = km_df)
plot_age <- ggsurvplot(fit_age, data = km_df, pval = TRUE, conf.int = FALSE,
                       risk.table = TRUE, title = "KM Curve by Age Group")

# KM plot by sex
fit_sex <- survfit(Surv(time = km_df$time_mort, event = km_df$mortstat) ~ sex, data = km_df)
plot_sex <- ggsurvplot(fit_sex, data = km_df, pval = TRUE, conf.int = FALSE,
                       risk.table = TRUE, title = "KM Curve by Sex")

# KM plot by race
fit_race <- survfit(Surv(time = km_df$time_mort, event = km_df$mortstat) ~ Race, data = km_df)
plot_race <- ggsurvplot(fit_race, data = km_df, pval = TRUE, conf.int = FALSE,
                        risk.table = TRUE, title = "KM Curve by Race")

# KM plot by mobility
fit_mobility <- survfit(Surv(time = km_df$time_mort, event = km_df$mortstat) ~ Mobility, data = km_df)
plot_mobility <- ggsurvplot(fit_mobility, data = km_df, pval = TRUE, conf.int = FALSE,
                            risk.table = TRUE, title = "KM Curve by Mobility")

# Combine the plots: 2 columns × 2 rows
panel_plot <- ggarrange(
  plot_age$plot, plot_sex$plot,
  plot_race$plot, plot_mobility$plot,
  ncol = 2, nrow = 2,
  labels = c("A", "B", "C", "D")
)

# To display the plot
print(panel_plot)


##############################
## FPCA & Data Preprocessing
##############################

fpca_result <- fpca.sc(Y = df$actmat, pve = 0.99)  #can use fast alternative#
scoremat<-fpca_result$scores

dfrfs <- subset(df, select = -c(actmat))

# Automatically add columns for all principal components
for (i in 1:ncol(scoremat)) {
  dfrfs[[paste0("PC", i)]] <- scoremat[, i]
}
# Columns to convert to numeric
cols_to_convert <- c("Race", "sex", "Mobility")

# Convert specified columns to numeric
dfrfs[cols_to_convert] <- lapply(dfrfs[cols_to_convert], as.factor)

#################### Cox_Regression #############################

# Set seed for reproducibility
set.seed(123)

# Create train-test split (e.g., 80% train, 20% test)
sample_indices <- sample(1:nrow(dfrfs), size = 0.8 * nrow(dfrfs))
train_data <- dfrfs[sample_indices, ]
test_data  <- dfrfs[-sample_indices, ]

cox_model <- coxph(Surv(time_mort, mortstat) ~ Age + Race + BMI + sex + 
                     Mobility + diabetes.y + poverty_level + 
                     Asthma + Arthritis + heart_failure + 
                     coronary_heart_disease + angina + stroke + thyroid + 
                     bronchitis + cancer  + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + 
                     PC7 + PC8 + PC9, data = dfrfs) #do in whole data

summary(cox_model)
#POSITIVE NEGATIVE SIGN INTERPRET##
##JUST COMBINE BETA AND EIGENFUNCC FOR SOME PLOT###
betascore<-coef(cox_model)[21:29]
betafunc<-fpca_result$efunctions%*%betascore
S<-c(0:1439)/60
plot(S,betafunc,type="l",ylab="Beta(Time of Day)",xlab="Time of Day")
abline(h=0,col="grey",lty=2)

##
# Compute risk scores
risk_cox <- rowSums(predict(cox_model, test_data, type = "terms"))

# Compute C-index
cal_c(marker = risk_cox, Stime = test_data$time_mort, 
      status = test_data$mortstat)


###
cox_model2 <- coxph(Surv(time_mort, mortstat) ~ Age + Race + BMI + sex + 
                     Mobility + diabetes.y + poverty_level + 
                     Asthma + Arthritis + heart_failure + 
                     coronary_heart_disease + angina + stroke + thyroid + 
                     bronchitis + cancer  + act_mean, data = dfrfs) #do in whole data

summary(cox_model2)


#

################## RSF ########################

library(randomForestSRC)
library(SurvMetrics)

rsf_tuned <- rfsrc(Surv(time_mort,mortstat) ~ Age + Race + BMI + sex + 
                     Mobility + diabetes.y + poverty_level + 
                     Asthma + Arthritis + heart_failure + 
                     coronary_heart_disease + angina + stroke + thyroid + 
                     bronchitis + cancer  + PC1 + PC2 + PC3 + PC4 + PC5 + 
                     PC6 + PC7 + PC8 + PC9,
                   data = dfrfs, ntree = 500, mtry = 3, 
                   nodesize = 15, nsplit = 10, importance = TRUE)

# Step 4: View Model Summary
print(rsf_tuned)

# Plot the OOB survival curve
plot(rsf_tuned)

# Load required library
library(ggplot2)

# Extract variable importance
var_imp <- rsf_tuned$importance

# Convert to data frame
var_imp_df <- data.frame(
  Variable = names(var_imp),
  Importance = var_imp
)

# Sort by importance
var_imp_df <- var_imp_df[order(var_imp_df$Importance, decreasing = TRUE), ]

# Plot using ggplot2
ggplot(var_imp_df, aes(x = reorder(Variable, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Variable Importance (Random Survival Forest)",
       x = "Variable",
       y = "Importance") +
  theme_minimal()


#################################
## 100 Split 
#################################

# Set path to your folder containing the CSV files
data_dir <- "D:/Final_Year_Project/splits_final"

cal_c <- function(marker, Stime, status){
  utimes <- sort(unique(Stime[status==1]))
  num <- denom <- 0
  for(ut in seq_along(utimes)){
    ti    <- utimes[ut]
    inx_i <- which(Stime == ti & status==1)
    inx_j <- which(Stime > ti)
    n_case_t    <- length(inx_i)
    n_control_t <- length(inx_j)
    for(i in seq_along(inx_i)){
      num   <- num + sum( (marker[inx_j] > marker[inx_i[i]] ) ) + 0.5*sum(marker[inx_j] == marker[inx_i[i]])
    }
    denom <- denom + n_case_t*n_control_t
  }
  1-num/denom
}
##############################
## Cox PH Model
##############################


## Mean Act based
calculate_C_index_cox <- function(n_splits = 100, data_path = "./splits", 
                                  train_prefix = "train_split_", 
                                  test_prefix = "test_split_",
                                  factor_cols = c("Race", "sex", "Mobility", 
                                                  "diabetes.y", "Asthma", 
                                                  "Arthritis", 
                                                  "heart_failure", 
                                                  "coronary_heart_disease", 
                                                  "angina", "stroke", 
                                                  "thyroid", "bronchitis", 
                                                  "cancer")) {
  c_index_values <- numeric(n_splits)
  
  for (i in 1:n_splits) {
    train_file <- file.path(data_path, paste0(train_prefix, i, ".csv"))
    test_file  <- file.path(data_path, paste0(test_prefix, i, ".csv"))
    
    df_train <- read.csv(train_file)
    df_test  <- read.csv(test_file)
    
    # Convert specified columns to factors in both train and test sets
    for (col in factor_cols) {
      df_train[[col]] <- as.factor(df_train[[col]])
      df_test[[col]]  <- as.factor(df_test[[col]])
    }
    
    # Fit Cox model
    cox_model <- coxph(Surv(time_mort, mortstat) ~ Age + Race + BMI + sex + 
                         Mobility + diabetes.y + poverty_level + 
                         Asthma + Arthritis + heart_failure + 
                         coronary_heart_disease + angina + stroke + thyroid + 
                         bronchitis + cancer + act_mean, data = df_train)
    
    # Compute risk scores
    risk_cox <- rowSums(predict(cox_model, df_test, type = "terms"))
    
    # Compute C-index
    c_index_values[i] <- cal_c(marker = risk_cox, Stime = df_test$time_mort, 
                               status = df_test$mortstat)
  }
  
  return(c_index_values)
}

c_index_cox_MeanAct <- calculate_C_index_cox(data_path = data_dir)
mean_c_index_cox_MeanAct <- mean(c_index_cox_MeanAct)

## PC based

calculate_C_index_cox <- function(n_splits = 100, data_path = "./splits", 
                                  train_prefix = "train_split_", 
                                  test_prefix = "test_split_",
                                  factor_cols = c("Race", "sex", "Mobility", 
                                                  "diabetes.y", "Asthma", 
                                                  "Arthritis", 
                                                  "heart_failure", 
                                                  "coronary_heart_disease", 
                                                  "angina", "stroke", 
                                                  "thyroid", "bronchitis", 
                                                  "cancer")) {
  c_index_values <- numeric(n_splits)
  
  for (i in 1:n_splits) {
    train_file <- file.path(data_path, paste0(train_prefix, i, ".csv"))
    test_file  <- file.path(data_path, paste0(test_prefix, i, ".csv"))
    
    df_train <- read.csv(train_file)
    df_test  <- read.csv(test_file)
    
    # Convert specified columns to factors in both train and test sets
    for (col in factor_cols) {
      df_train[[col]] <- as.factor(df_train[[col]])
      df_test[[col]]  <- as.factor(df_test[[col]])
    }
    
    # Fit Cox model
    cox_model <- coxph(Surv(time_mort, mortstat) ~ Age + Race + BMI + sex + 
                         Mobility + diabetes.y + poverty_level + 
                         Asthma + Arthritis + heart_failure + 
                         coronary_heart_disease + angina + stroke + thyroid + 
                         bronchitis + cancer + PC1 + PC2 + PC3 + PC4 + PC5 + 
                         PC6 + PC7 + PC8 + PC9, data = df_train)
    
    # Compute risk scores
    risk_cox <- rowSums(predict(cox_model, df_test, type = "terms"))
    
    # Compute C-index
    c_index_values[i] <- cal_c(marker = risk_cox, Stime = df_test$time_mort, 
                               status = df_test$mortstat)
  }
  
  return(c_index_values)
}

c_index_cox_FPC <- calculate_C_index_cox(data_path = data_dir)
mean_c_index_cox_FPC <- mean(c_index_cox_FPC)

################################
## RSF
################################

## Mean Act based

calculate_C_index_rsf <- function(n_splits = 100, data_path = "./splits", 
                                  train_prefix = "train_split_", 
                                  test_prefix = "test_split_",
                                  factor_cols = c("Race", "sex", "Mobility", 
                                                  "diabetes.y", "Asthma", 
                                                  "Arthritis", 
                                                  "heart_failure", 
                                                  "coronary_heart_disease", 
                                                  "angina", "stroke", 
                                                  "thyroid", "bronchitis", 
                                                  "cancer")) {
  c_index_values <- numeric(n_splits)
  
  for (i in 1:n_splits) {
    # File paths
    train_file <- file.path(data_path, paste0(train_prefix, i, ".csv"))
    test_file  <- file.path(data_path, paste0(test_prefix, i, ".csv"))
    
    # Read train and test data
    df_train <- read.csv(train_file)
    df_test  <- read.csv(test_file)
    
    # Convert specified columns to factors
    for (col in factor_cols) {
      df_train[[col]] <- as.factor(df_train[[col]])
      df_test[[col]]  <- as.factor(df_test[[col]])
    }
    
    # Fit RSF model
    rsf_tuned <- rfsrc(Surv(time_mort, mortstat) ~ Age + Race + BMI + sex + 
                         Mobility + diabetes.y + poverty_level + 
                         Asthma + Arthritis + heart_failure + 
                         coronary_heart_disease + angina + stroke + thyroid + 
                         bronchitis + cancer + act_mean,
                       data = df_train, ntree = 500, mtry = 3, 
                       nodesize = 15, nsplit = 10, importance = TRUE)
    
    # Determine median death time index for C-index computation
    dis_time <- sort(unique(df_train$time_mort[df_train$mortstat == 1]))
    med_index <- median(1:length(dis_time))
    
    # Predict survival
    pred <- predict(rsf_tuned, newdata = df_test, times = dis_time, se = FALSE)
    surv_probs <- pred$survival
    
    # Compute C-index
    surv_obj <- Surv(df_test$time_mort, df_test$mortstat)
    c_index_values[i] <- Cindex(surv_obj, predicted = surv_probs[, med_index])
  }
  
  return(c_index_values)
}

c_index_rsf_MeanAct <- calculate_C_index_rsf(data_path = data_dir)
mean_c_index_rsf_MeanAct <- mean(c_index_rsf_MeanAct)

## PC based

calculate_C_index_rsf <- function(n_splits = 100, data_path = "./splits", 
                                  train_prefix = "train_split_", 
                                  test_prefix = "test_split_",
                                  factor_cols = c("Race", "sex", "Mobility", 
                                                  "diabetes.y", "Asthma", 
                                                  "Arthritis", 
                                                  "heart_failure", 
                                                  "coronary_heart_disease", 
                                                  "angina", "stroke", 
                                                  "thyroid", "bronchitis", 
                                                  "cancer")) {
  c_index_values <- numeric(n_splits)
  
  for (i in 1:n_splits) {
    # File paths
    train_file <- file.path(data_path, paste0(train_prefix, i, ".csv"))
    test_file  <- file.path(data_path, paste0(test_prefix, i, ".csv"))
    
    # Read train and test data
    df_train <- read.csv(train_file)
    df_test  <- read.csv(test_file)
    
    # Convert specified columns to factors
    for (col in factor_cols) {
      df_train[[col]] <- as.factor(df_train[[col]])
      df_test[[col]]  <- as.factor(df_test[[col]])
    }
    
    # Fit RSF model
    rsf_tuned <- rfsrc(Surv(time_mort, mortstat) ~ Age + Race + BMI + sex + 
                         Mobility + diabetes.y + poverty_level + 
                         Asthma + Arthritis + heart_failure + 
                         coronary_heart_disease + angina + stroke + thyroid + 
                         bronchitis + cancer + PC1 + PC2 + PC3 + PC4 + PC5 + 
                         PC6 + PC7 + PC8 + PC9,
                       data = df_train, ntree = 500, mtry = 3, 
                       nodesize = 15, nsplit = 10, importance = TRUE)
    
    # Determine median death time index for C-index computation
    dis_time <- sort(unique(df_train$time_mort[df_train$mortstat == 1]))
    med_index <- median(1:length(dis_time))
    
    # Predict survival
    pred <- predict(rsf_tuned, newdata = df_test, times = dis_time, se = FALSE)
    surv_probs <- pred$survival
    
    # Compute C-index
    surv_obj <- Surv(df_test$time_mort, df_test$mortstat)
    c_index_values[i] <- Cindex(surv_obj, predicted = surv_probs[, med_index])
  }
  
  return(c_index_values)
}

c_index_rsf_FPC <- calculate_C_index_rsf(data_path = data_dir)
mean_c_index_rsf_FPC <- mean(c_index_rsf_FPC)


##########################
## Result Dataframe
##########################

results_df <- data.frame(
  c_index_cox_MeanAct = c_index_cox_MeanAct,
  c_index_cox_FPC = c_index_cox_FPC,
  c_index_rsf_MeanAct = c_index_rsf_MeanAct,
  c_index_rsf_FPC = c_index_rsf_FPC
)
