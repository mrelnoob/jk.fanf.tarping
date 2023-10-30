# ******************************************************************************************** #
# ******************************************************************************************** #
# ********************************** MODELLING eradication *********************************** #
# ******************************************************************************************** #
# ******************************************************************************************** #


# ================================================= #
##### Data preparation for modelling 'eradication' #####
# ================================================= #

### Sourcing previous scripts:
source(file = here::here("R/01_03_exploration_functions.R")) # I need to source my previous scripts (R files)
# to be able to import my cleaned dataset (with my 'model_datasets()' function) and perform my analyses.

.pardefault <- par() # To save the default graphical parameters (in case I want to restore them).

readr::read_csv(here::here("data", "erad.csv"), col_names = TRUE, col_types =
                  readr::cols(
                    manager_id = readr::col_factor(),
                    xp_id = readr::col_factor(),
                    eff_eradication = readr::col_factor(c("0", "1")),
                    eradication = readr::col_factor(c("0", "1")),
                    goals = readr::col_factor(),
                    geomem = readr::col_factor(c("0", "1")),
                    maxveg = readr::col_factor(c("0", "1")),
                    uprootexcav = readr::col_factor(c("0", "1")),
                    fully_tarped = readr::col_factor(c("0", "1")),
                    stripsoverlap_ok = readr::col_factor(c("0", "1")),
                    tarpfix_multimethod = readr::col_factor(c("0", "1")),
                    tarpfix_pierced = readr::col_factor(c("0", "1")),
                    plantation = readr::col_factor(c("0", "1")),
                    repairs = readr::col_factor(c("0", "1")),
                    add_control = readr::col_factor(c("0", "1")),
                    pb_fixation = readr::col_factor(c("0", "1")),
                    pb_durability = readr::col_factor(c("0", "1")))) -> eff
summary(eff)



##### Final pre-modelling assumption checks ##### (run only when required)
# --------------------------------------------- #

# ### Testing the relevance of the random effect structure:
# m0.glm <- stats::glm(eradication ~ fully_tarped, data = eff, family = binomial)
# m0.glmer <- lme4::glmer(eradication ~ fully_tarped + (1|manager_id), data = eff, family = binomial)
# aic.glm <- AIC(logLik(m0.glm))
# aic.glmer <- AIC(logLik(m0.glmer))
#
# # Likelihood Ratio Test:
# null.id <- -2 * logLik(m0.glm) + 2 * logLik(m0.glmer)
# pchisq(as.numeric(null.id), df=1, lower.tail=F)
# rm(m0.glm, m0.glmer, aic.glm, aic.glmer, null.id)
# # The Likelihood Ratio Tests are NOT significant so the use of the random effect structure may not be
# # necessary! However, further tests on the model residuals may indicate otherwise.


# IMPORTANT NOTE:
# Since we'll use a regularization modelling method (i.e. penalized regression) to avoid overfitting
# and deal with our "complete separation" problem, we will not use model selection and multimodel inference
# for this response variable.
# Consequently, we will assess further logistic regression assumptions using the most parsimonious full model
# we could think of, i.e. including all available predictors we thought should be important to predict
# our outcome (i.e. the eradication or near-eradication of knotweeds after tarping), based on our knowledge
# of the biology and ecology of the plant.





# =========================================================================== #
##### Cross-validated Ridge logistic regression with bootstrap validation #####
# =========================================================================== #

##### Comparing Ridge model with and without interactions #####
# ----------------------------------------------------------- #

### Data preparation:

mydata <- eff[,c("eradication", "add_control", "distance", "fully_tarped", "geomem", "obstacles",
                 "plantation", "pb_fixation", "pb_durability", "repairs", "slope", "stand_surface",
                 "stripsoverlap_ok", "sedicover_height", "tarping_duration", "uprootexcav")]
# Create data matrices (as accepted by glmnet):
x <- stats::model.matrix(eradication~., mydata)[,-1] # Matrix of potential predictors (WITHOUT interactions)
f <- as.formula(eradication~add_control+distance+fully_tarped+geomem+obstacles+plantation+pb_fixation
                +pb_durability+repairs+slope+stand_surface+stripsoverlap_ok+sedicover_height+tarping_duration
                +uprootexcav+pb_fixation*repairs+pb_fixation*add_control+add_control*fully_tarped
                +distance*fully_tarped+distance*uprootexcav+fully_tarped*geomem+stripsoverlap_ok*pb_durability
                +uprootexcav*pb_fixation+pb_fixation*pb_durability+uprootexcav*plantation+repairs*plantation+0)
x_int <- stats::model.matrix(f, mydata)[,-1] # Matrix of potential predictors (WITH interactions)
y <- as.numfactor(x = mydata$eradication) %>% as.matrix()
# NOTE: to create all possible interactions, the formula should be: f <- eradication~.^2 (for second-order
# interactions, or ^3 for third order)!



### Ridge regression:
# Find the optimal value of lambda that minimizes the cross-validation error:
set.seed(653)
cv.ridge <- glmnet::cv.glmnet(x = x, y = y, alpha = 0, family = "binomial",
                              type.measure = "deviance", nfolds = 10)
cv.ridge_int <- glmnet::cv.glmnet(x = x_int, y = y, alpha = 0, family = "binomial",
                              type.measure = "deviance", nfolds = 10)
plot(cv.ridge)
plot(cv.ridge_int)
# glmnet::coef.glmnet(object = cv.ridge, s = cv.ridge$lambda.min) # To have look at the coefficients
# glmnet::coef.glmnet(object = cv.ridge_int, s = cv.ridge$lambda.min) # To have look at the coefficients

# Compute the full Ridge model:
ridge.model <- glmnet::glmnet(x = x, y = y, alpha = 0, family = "binomial")
ridge.model_int <- glmnet::glmnet(x = x_int, y = y, alpha = 0, family = "binomial")



### Evaluate the performance of the full Ridge model (WITHOUT interactions):
pred.class <- predict(object = ridge.model, newx = x,
                         s = min(ridge.model$lambda), type = "class") %>% as.factor()
pred.prob  <- predict(ridge.model, newx = x, s = min(ridge.model$lambda), type = "response")

MLmetrics::LogLoss(y_pred = pred.prob, y_true = y) # 0.447
MLmetrics::AUC(y_pred = pred.prob, y_true = y) # 0.865
MLmetrics::R2_Score(y_pred = pred.prob, y_true = y) # 0.352
stats::deviance(ridge.model) # 75.95

conf.mat <- MLmetrics::ConfusionMatrix(y_pred = pred.class, y_true = y) # Or table(y, pred.class)
error_rate <- (conf.mat[2]+conf.mat[3])/length(y)
error_rate # Mean error rate = 0.188
tpr <- (conf.mat[4])/(conf.mat[2]+conf.mat[4])
tpr # True positive rate = 0.7

# # To plot the ROC curve:
# pred <- ROCR::prediction(predictions = pred.prob, labels = mydata$eradication)
# auc <- ROCR::performance(prediction.obj = pred, measure = "auc")@y.values[[1]][1]
# perf <- ROCR::performance(prediction.obj = pred, measure = "tpr","fpr")
# plot(perf, col="navyblue", cex.main=1,
#      main= paste("Logistic Regression ROC Curve: AUC =", round(auc,3)))
# abline(a=0, b = 1, col='darkorange1')



### Evaluate the performance of the full Ridge model (WITH interactions):
pred.class_int <- predict(object = ridge.model_int, newx = x_int,
                         s = min(ridge.model_int$lambda), type = "class") %>% as.factor()
pred.prob_int  <- predict(ridge.model_int, newx = x_int, s = min(ridge.model_int$lambda), type = "response")

MLmetrics::LogLoss(y_pred = pred.prob_int, y_true = y) # 0.399
MLmetrics::AUC(y_pred = pred.prob_int, y_true = y) # 0.894
MLmetrics::R2_Score(y_pred = pred.prob_int, y_true = y) # 0.43
stats::deviance(ridge.model_int) # 67.878

conf.mat_int <- MLmetrics::ConfusionMatrix(y_pred = pred.class_int, y_true = y) # Or table(y, pred.class_int)
error_rate_int <- (conf.mat_int[2]+conf.mat_int[3])/length(y)
error_rate_int # Mean error rate = 0.176
tpr_int <- (conf.mat_int[4])/(conf.mat_int[2]+conf.mat_int[4])
tpr_int # True positive rate = 0.73

# # To plot the ROC curve:
# pred <- ROCR::prediction(predictions = pred.prob_int, labels = mydata$eradication)
# auc <- ROCR::performance(prediction.obj = pred, measure = "auc")@y.values[[1]][1]
# perf <- ROCR::performance(prediction.obj = pred, measure = "tpr","fpr")
# plot(perf, col="navyblue", cex.main=1,
#      main= paste("Logistic Regression ROC Curve: AUC =", round(auc,3)))
# abline(a=0, b = 1, col='darkorange1')



### Enhanced (optimism) bootstrap comparison:
# Create a function suitable for boot that will return the optimism estimates for statistics testing
# models against the full original sample:
compare_opt <- function(orig_data, i){
  # Create the resampled data
  train_data <- orig_data[i, ]

  x.train <- stats::model.matrix(eradication~., train_data)[,-1] # Matrix of potential predictors (WITHOUT interactions)
  f <- as.formula(eradication~add_control+distance+fully_tarped+geomem+obstacles+plantation+pb_fixation
                  +pb_durability+repairs+slope+stand_surface+stripsoverlap_ok+sedicover_height+tarping_duration
                  +uprootexcav+pb_fixation*repairs+pb_fixation*add_control+add_control*fully_tarped
                  +distance*fully_tarped+distance*uprootexcav+fully_tarped*geomem+stripsoverlap_ok*pb_durability
                  +uprootexcav*pb_fixation+pb_fixation*pb_durability+uprootexcav*plantation+repairs*plantation+0)
  x.train_int <- stats::model.matrix(f, train_data)[,-1] # Matrix of potential predictors (WITH interactions)
  y.train <- as.numfactor(x = train_data$eradication) %>% as.matrix() # Same for Y

  # Run the entire modelling process:
  model_full <- glmnet::glmnet(x = x.train, y = y.train, alpha = 0, family = "binomial")
  model_full_int <- glmnet::glmnet(x = x.train_int, y = y.train, alpha = 0, family = "binomial")

  # Predict the values on the trained, resampled data:
  train_pred.class <- predict(object = model_full, newx = x.train,
                              s = min(model_full$lambda), type = "class") %>% as.factor()
  train_pred.prob  <- predict(object = model_full, newx = x.train,
                              s = min(model_full$lambda), type = "response")
  train_pred.class_int <- predict(object = model_full_int, newx = x.train_int,
                              s = min(model_full_int$lambda), type = "class") %>% as.factor()
  train_pred.prob_int  <- predict(object = model_full_int, newx = x.train_int,
                              s = min(model_full_int$lambda), type = "response")

  # Predict the values on the original, unresampled data:
  full_pred.class <- predict(object = model_full, newx = x,
                             s = min(model_full$lambda), type = "class") %>% as.factor()
  full_pred.prob  <- predict(object = model_full, newx = x,
                             s = min(model_full$lambda), type = "response")
  full_pred.class_int <- predict(object = model_full_int, newx = x_int,
                             s = min(model_full_int$lambda), type = "class") %>% as.factor()
  full_pred.prob_int  <- predict(object = model_full_int, newx = x_int,
                             s = min(model_full_int$lambda), type = "response")

  train_conf.mat <- MLmetrics::ConfusionMatrix(y_pred = train_pred.class, y_true = y.train)
  train_error_rate <- (train_conf.mat[2]+train_conf.mat[3])/length(y.train)
  train_tpr <- (train_conf.mat[4])/(train_conf.mat[2]+train_conf.mat[4])
  full_conf.mat <- MLmetrics::ConfusionMatrix(y_pred = full_pred.class, y_true = y)
  full_error_rate <- (full_conf.mat[2]+full_conf.mat[3])/length(y)
  full_tpr <- (full_conf.mat[4])/(full_conf.mat[2]+full_conf.mat[4])

  train_conf.mat_int <- MLmetrics::ConfusionMatrix(y_pred = train_pred.class_int, y_true = y.train)
  train_error_rate_int <- (train_conf.mat_int[2]+train_conf.mat_int[3])/length(y.train)
  train_tpr_int <- (train_conf.mat_int[4])/(train_conf.mat_int[2]+train_conf.mat_int[4])
  full_conf.mat_int <- MLmetrics::ConfusionMatrix(y_pred = full_pred.class_int, y_true = y)
  full_error_rate_int <- (full_conf.mat_int[2]+full_conf.mat_int[3])/length(y)
  full_tpr_int <- (full_conf.mat_int[4])/(full_conf.mat_int[2]+full_conf.mat_int[4])


  # Return a vector of summary optimism results
  results <- c(
    boot_LogLoss = MLmetrics::LogLoss(y_pred = train_pred.prob, y_true = y.train) -
      MLmetrics::LogLoss(y_pred = full_pred.prob, y_true = y),
    boot_AUC = MLmetrics::AUC(y_pred = train_pred.prob, y_true = y.train) -
      MLmetrics::AUC(y_pred = full_pred.prob, y_true = y),
    boot_R2 = MLmetrics::R2_Score(y_pred = train_pred.prob, y_true = y.train) -
      MLmetrics::R2_Score(y_pred = full_pred.prob, y_true = y),
    boot_m.error = train_error_rate - full_error_rate,
    boot_tpr = train_tpr - full_tpr,
    boot_LogLoss_int = MLmetrics::LogLoss(y_pred = train_pred.prob_int, y_true = y.train) -
      MLmetrics::LogLoss(y_pred = full_pred.prob_int, y_true = y),
    boot_AUC_int = MLmetrics::AUC(y_pred = train_pred.prob_int, y_true = y.train) -
      MLmetrics::AUC(y_pred = full_pred.prob_int, y_true = y),
    boot_R2_int = MLmetrics::R2_Score(y_pred = train_pred.prob_int, y_true = y.train) -
      MLmetrics::R2_Score(y_pred = full_pred.prob_int, y_true = y),
    boot_m.error_int = train_error_rate_int - full_error_rate_int,
    boot_tpr_int = train_tpr_int - full_tpr_int
  )
  return(results)
}



### Perform bootstrapping and return optimism-corrected results:
res_opt <- boot::boot(data = mydata, statistic = compare_opt, R = 1000) # For large datasets, use a parallel
# cluster computing (it's quite straightforward).

# Calculate the results for the original model:
original <- c(
  MLmetrics::LogLoss(y_pred = pred.prob, y_true = y),
  MLmetrics::AUC(y_pred = pred.prob, y_true = y),
  MLmetrics::R2_Score(y_pred = pred.prob, y_true = y),
  error_rate,
  tpr,
  MLmetrics::LogLoss(y_pred = pred.prob_int, y_true = y),
  MLmetrics::AUC(y_pred = pred.prob_int, y_true = y),
  MLmetrics::R2_Score(y_pred = pred.prob_int, y_true = y),
  error_rate_int,
  tpr_int
)
res_original <- as.data.frame(x = original) %>%
  dplyr::rename("Initial value" = original)
res_original$metrics <- c("LogLoss", "AUC",
                          "R2", "Mean Error Rate", "True Positive Rate",
                          "LogLoss (w/ inter.)", "AUC (w/ inter.)",
                          "R2 (w/ inter.)", "Mean Error Rate (w/ inter.)",
                          "True Positive Rate (w/ inter.)")

# Compute the mean bootstrapped results and the enhanced results:
optimism <- apply(na.omit(res_opt$t), 2, mean)
corrected_results <- original - optimism %>% as.data.frame()
res_corrected <- dplyr::rename(.data = corrected_results, "Optimism-corrected value" = .)
res_corrected$metrics <- c("LogLoss", "AUC",
                           "R2", "Mean Error Rate", "True Positive Rate",
                           "LogLoss (w/ inter.)", "AUC (w/ inter.)",
                           "R2 (w/ inter.)", "Mean Error Rate (w/ inter.)",
                           "True Positive Rate (w/ inter.)")
# Note that, for an unknown reason, I didn't manage to make the random number generator work as I want it to, so
# the results vary every time you launch the code again (the seed is not set appropriately). It's not very important
# as they always are in the same range of results but I should fix it one day.



### Export result tables
readr::write_csv2(x = res_original, file = here::here("output", "tables", "Orig.predic.metrics_eradication.csv"))
readr::write_csv2(x = res_corrected, file = here::here("output", "tables", "Correct.predic.metrics_eradication.csv"))





##### Predicting new (near-)eradication events #####
# ------------------------------------------------ #

### Generate new data
new.dat <- mydata %>%
  dplyr::mutate_if(is.numeric, mean) %>%
  dplyr::mutate(add_control = 1) %>%
  dplyr::mutate(fully_tarped = 1) %>%
  dplyr::mutate(geomem = 1) %>%
  dplyr::mutate(plantation = 0) %>%
  dplyr::mutate(pb_fixation = 0) %>%
  dplyr::mutate(pb_durability = 0) %>%
  dplyr::mutate(repairs = 1) %>%
  dplyr::mutate(stripsoverlap_ok = 1) %>%
  dplyr::mutate(uprootexcav = 1)
new.dat <- new.dat[1:50,]
summary(new.dat) # Generated a new dataset simulating knotweed stands with average values for numeric
# variables and that was fully covered by a geomembrane after having been uprooted and with at least 40cm
# of fabric strips overlap, the presence of additional control, repairs, no plantation, and no particular
# problems.
new.dat_nf <- new.dat %>% dplyr::mutate(fully_tarped = 0)

# For a range of distances and 3 stand surfaces, when FULLY TARPED
distance <- seq(0.1, 5, length.out = 50)
new.dist <- new.dat
new.dist$distance <- distance
new.dist5 <- dplyr::mutate(.data = new.dist,
                           stand_surface = 5)
new.dist50 <- dplyr::mutate(.data = new.dist,
                           stand_surface = 50)
new.dist500 <- dplyr::mutate(.data = new.dist,
                           stand_surface = 500)
new.dist5 <- stats::model.matrix(f, new.dist5) # Matrix of potential predictors to be used by predict().
new.dist50 <- stats::model.matrix(f, new.dist50)
new.dist500 <- stats::model.matrix(f, new.dist500)

# For a range of distances and 3 stand surfaces, when NOT FULLY TARPED:
distance <- seq(0.1, 5, length.out = 50)
new.dist <- new.dat_nf
new.dist$distance <- distance
new.dist_nf5 <- dplyr::mutate(.data = new.dist,
                           stand_surface = 5)
new.dist_nf50 <- dplyr::mutate(.data = new.dist,
                            stand_surface = 50)
new.dist_nf500 <- dplyr::mutate(.data = new.dist,
                             stand_surface = 500)
new.dist_nf5 <- stats::model.matrix(f, new.dist_nf5) # Matrix of potential predictors to be used by predict().
new.dist_nf50 <- stats::model.matrix(f, new.dist_nf50)
new.dist_nf500 <- stats::model.matrix(f, new.dist_nf500)



### Predict probabilities of (near-)eradication for the new data:
# Predictions for a range of distance values (according to some stand surfaces)
pred.prob_dist5  <- predict(ridge.model_int, newx = new.dist5, s = min(ridge.model_int$lambda), type = "response")
pred.prob_dist50  <- predict(ridge.model_int, newx = new.dist50, s = min(ridge.model_int$lambda), type = "response")
pred.prob_dist500  <- predict(ridge.model_int, newx = new.dist500, s = min(ridge.model_int$lambda), type = "response")

pred.prob_dist_nf5  <- predict(ridge.model_int, newx = new.dist_nf5, s = min(ridge.model_int$lambda), type = "response")
pred.prob_dist_nf50  <- predict(ridge.model_int, newx = new.dist_nf50, s = min(ridge.model_int$lambda), type = "response")
pred.prob_dist_nf500  <- predict(ridge.model_int, newx = new.dist_nf500, s = min(ridge.model_int$lambda), type = "response")



### Plot predictions and export the plot:
# In JPEG format
grDevices::jpeg(filename = here::here("output", "plots", "Figure_1.jpeg"))

graphics::par(cex.lab=1.3, font.lab=2, bty = "n", fg = "gray35",
    col.axis = "gray35", col.lab = "gray20", cex = 0.9, tcl = -0.3,
    mgp = c(2, 0.6, 0.1), oma = c(1, 0, 1, 0), lab = c(5, 10, 7))
graphics::plot(pred.prob_dist5~new.dist5[,2], ylim = c(0,1), type = "n",
               xlab = "Distance (m)", ylab = "(Near-) eradication probability")
graphics::lines(pred.prob_dist5~new.dist5[,2], data = new.dist5, col = "gold", lwd = 3, lty = 1,
                panel.first = {grid(col="lavender",nx = 5,ny = 9, lty = 6)})
graphics::lines(pred.prob_dist50~new.dist50[,2], data = new.dist50, col = "sandybrown", lwd = 3, lty = 1)
graphics::lines(pred.prob_dist500~new.dist500[,2], data = new.dist500, col = "chocolate4", lwd = 3, lty = 1)
graphics::lines(pred.prob_dist_nf5~new.dist_nf5[,2], data = new.dist_nf5, col = "gold", lwd = 3, lty = 4)
graphics::lines(pred.prob_dist_nf50~new.dist_nf50[,2], data = new.dist_nf50, col = "sandybrown", lwd = 3, lty = 4)
graphics::lines(pred.prob_dist_nf500~new.dist_nf500[,2], data = new.dist_nf500, col = "chocolate4", lwd = 3, lty = 4)
graphics::legend(x = 3, y = 0.6,
                 legend = c("Fully tarped", "Not fully tarped", "Stand surface 5 m2", "Stand surface 50 m2",
                            "Stand surface 500 m2"),
                 col = c("black", "black", "gold", "sandybrown", "chocolate4"), lty = c(1,4,1,1,1),
                 lwd = 3, bty = "n", cex = 1.1)
dev.off()

# In PDF format:
grDevices::pdf(file = here::here("output", "plots", "Figure_1.pdf"))

graphics::par(cex.lab=1.3, font.lab=2, bty = "n", fg = "gray35",
              col.axis = "gray35", col.lab = "gray20", cex = 0.9, tcl = -0.3,
              mgp = c(2, 0.6, 0.1), oma = c(1, 0, 1, 0), lab = c(5, 10, 7))
graphics::plot(pred.prob_dist5~new.dist5[,2], ylim = c(0,1), type = "n",
               xlab = "Distance (m)", ylab = "(Near-) eradication probability")
graphics::lines(pred.prob_dist5~new.dist5[,2], data = new.dist5, col = "gold", lwd = 3, lty = 1,
                panel.first = {grid(col="lavender",nx = 5,ny = 9, lty = 6)})
graphics::lines(pred.prob_dist50~new.dist50[,2], data = new.dist50, col = "sandybrown", lwd = 3, lty = 1)
graphics::lines(pred.prob_dist500~new.dist500[,2], data = new.dist500, col = "chocolate4", lwd = 3, lty = 1)
graphics::lines(pred.prob_dist_nf5~new.dist_nf5[,2], data = new.dist_nf5, col = "gold", lwd = 3, lty = 4)
graphics::lines(pred.prob_dist_nf50~new.dist_nf50[,2], data = new.dist_nf50, col = "sandybrown", lwd = 3, lty = 4)
graphics::lines(pred.prob_dist_nf500~new.dist_nf500[,2], data = new.dist_nf500, col = "chocolate4", lwd = 3, lty = 4)
graphics::legend(x = 3, y = 0.6,
                 legend = c("Fully tarped", "Not fully tarped", "Stand surface 5 m2", "Stand surface 50 m2",
                            "Stand surface 500 m2"),
                 col = c("black", "black", "gold", "sandybrown", "chocolate4"), lty = c(1,4,1,1,1),
                 lwd = 3, bty = "n", cex = 1.1)
dev.off()
