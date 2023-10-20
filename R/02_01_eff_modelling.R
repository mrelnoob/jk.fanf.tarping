### Sourcing previous scripts:
source(file = here::here("R/01_03_exploration_functions.R")) # I need to source my previous scripts (R files)
# to be able to import my cleaned dataset (with my 'model_datasets()' function) and perform my analyses.



# ------------------------------------------------------ #
##### Data preparation for modelling 'effectiveness' #####
# ------------------------------------------------------ #

##### * 1.1. Data import and preparation ---------------------------------------
# ---------------------------------------------------------------------------- #

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
                    pb_durability = readr::col_factor(c("0", "1")))) %>%
  dplyr::mutate(effectiveness = effectiveness/10) %>% # For 'effectiveness' to be look like a percentage that could
  # be modelled using a Beta Regression Model.
  dplyr::mutate(effectiveness = ifelse(effectiveness == 1, 0.999, effectiveness)) %>%
  dplyr::mutate(distance_cent = scale(log2(distance+1), center = TRUE, scale = FALSE)) %>%
  dplyr::mutate(st_surface_cent = scale(log2(stand_surface), center = TRUE, scale = FALSE)) %>%
  dplyr::mutate(followups_cent = scale(x = followups, center = TRUE, scale = FALSE)) -> eff
summary(eff)



### Custom functions for modelling:
# To help creating QQ-plots for beta regression models:

qq.line = function(x) {
  # following four lines from base R's qqline()
  y <- quantile(x[!is.na(x)], c(0.25, 0.75))
  x <- qnorm(c(0.25, 0.75))
  slope <- diff(y)/diff(x)
  int <- y[1L] - slope * x[1L]
  return(c(int = int, slope = slope))
}





# ---------------------------------------- #
##### Building of the candidate models #####
# ---------------------------------------- #

Cand.mod <- list()
R.ajust <- data.frame(Model=integer(0), R2=numeric(0)) # Creates an empty data.frame with 2 variables



# ### Testing the relevance of the random effect structure:
# m0.glm <- glmmTMB::glmmTMB(formula = effectiveness~1, data = eff,
#                            family = glmmTMB::beta_family(link = "logit"), REML = FALSE)
# m0.glmer <- glmmTMB::glmmTMB(formula = effectiveness~(1|manager_id), data = eff,
#                              family = glmmTMB::beta_family(link = "logit"), REML = FALSE)
# aic.glm <- AIC(logLik(m0.glm))
# aic.glmer <- AIC(logLik(m0.glmer))
#
# # Likelihood Ratio Test:
# null.id <- -2 * logLik(m0.glm) + 2 * logLik(m0.glmer)
# pchisq(as.numeric(null.id), df=1, lower.tail=F) # The Likelihood Ratio Test is NOT significant so
# # the use of the random effect structure is not necessary!
# rm(m0.glm, m0.glmer, aic.glm, aic.glmer, null.id)



##### Model 1 (null model) #####
# ------------------------------

Cand.mod[[1]] <- glmmTMB::glmmTMB(formula = effectiveness~1, data = eff,
                                  family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

# ### Model evaluation (Flutterbys method)
# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[1]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[1]])))
# QQline <- qq.line(resid(Cand.mod[[1]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[1]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals:
# dat.sim <- simulate(Cand.mod[[1]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[1]]))
#
# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[1]])
# performance::check_collinearity(Cand.mod[[1]])
# performance::r2_nakagawa(Cand.mod[[1]])
#
# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[1]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[1]])
#
# ### Computing a quasiR^2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[1]], ~1))[1] - logLik(Cand.mod[[1]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=1, R2=R2[1]))



##### Model 2 #####
# -----------------

Cand.mod[[2]] <- glmmTMB::glmmTMB(formula = effectiveness~log2(distance+1) + followups, data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation (Fluttersbys method):

# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[2]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[2]])))
# QQline <- qq.line(resid(Cand.mod[[2]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[2]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_29):
# dat.sim <- simulate(Cand.mod[[2]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[2]]))

# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[2]])
# performance::check_collinearity(Cand.mod[[2]])
# performance::r2_nakagawa(Cand.mod[[2]])
#
# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[2]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[2]])
# family(Cand.mod[[2]])$linkinv(glmmTMB::fixef(Cand.mod[[2]])$cond) # To get interpretable coefficients.

### Computing a Pseudo-R2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[2]], ~1))[1] - logLik(Cand.mod[[2]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=2, R2=R2[1]))




##### Model 3 #####
# -----------------

Cand.mod[[3]] <- glmmTMB::glmmTMB(formula = effectiveness~log2(distance+1) + fully_tarped, data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

# ### Model evaluation (Fluttersbys method):
# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[3]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[3]])))
# QQline <- qq.line(resid(Cand.mod[[3]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[3]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_3):
# dat.sim <- simulate(Cand.mod[[3]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[3]]))
#
# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[3]])
# performance::check_collinearity(Cand.mod[[3]])
# performance::r2_nakagawa(Cand.mod[[3]])
#
# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[3]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[3]])
# family(Cand.mod[[3]])$linkinv(glmmTMB::fixef(Cand.mod[[3]])$cond) # To get interpretable coefficients.

### Computing a Pseudo-R2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[3]], ~1))[1] - logLik(Cand.mod[[3]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=3, R2=R2[1]))



##### Model 4 #####
# -----------------

Cand.mod[[4]] <- glmmTMB::glmmTMB(formula = effectiveness~log2(distance+1) + obstacles, data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation (Fluttersbys method):
# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[4]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[4]])))
# QQline <- qq.line(resid(Cand.mod[[4]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[4]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_29):
# dat.sim <- simulate(Cand.mod[[4]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[4]]))

# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[4]])
# performance::check_collinearity(Cand.mod[[4]])
# performance::r2_nakagawa(Cand.mod[[4]])
#
# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[4]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[4]])
# family(Cand.mod[[4]])$linkinv(glmmTMB::fixef(Cand.mod[[4]])$cond) # To get interpretable coefficients.

### Computing a Pseudo-R2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[4]], ~1))[1] - logLik(Cand.mod[[4]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=4, R2=R2[1]))



##### Model 5 #####
# -----------------

Cand.mod[[5]] <- glmmTMB::glmmTMB(formula = effectiveness~log2(distance+1) + log2(stand_surface), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation (Fluttersbys method):
# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[5]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[5]])))
# QQline <- qq.line(resid(Cand.mod[[5]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[5]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_29):
# dat.sim <- simulate(Cand.mod[[5]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[5]]))

# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[5]])
# performance::check_collinearity(Cand.mod[[5]])
# performance::r2_nakagawa(Cand.mod[[5]])
#
# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[5]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[5]])
# family(Cand.mod[[5]])$linkinv(glmmTMB::fixef(Cand.mod[[5]])$cond) # To get interpretable coefficients.

### Computing a Pseudo-R2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[5]], ~1))[1] - logLik(Cand.mod[[5]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=5, R2=R2[1]))



##### Model 6 #####
# -----------------

Cand.mod[[6]] <- glmmTMB::glmmTMB(formula = effectiveness~log2(distance+1) + log2(sedicover_height+1), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation (Fluttersbys method):
# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[6]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[6]])))
# QQline <- qq.line(resid(Cand.mod[[6]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[6]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_29):
# dat.sim <- simulate(Cand.mod[[6]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[6]]))

# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[6]])
# performance::check_collinearity(Cand.mod[[6]])
# performance::r2_nakagawa(Cand.mod[[6]])
#
# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[6]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[6]])
# family(Cand.mod[[6]])$linkinv(glmmTMB::fixef(Cand.mod[[6]])$cond) # To get interpretable coefficients.

### Computing a Pseudo-R2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[6]], ~1))[1] - logLik(Cand.mod[[6]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=6, R2=R2[1]))



##### Model 7 #####
# -----------------

Cand.mod[[7]] <- glmmTMB::glmmTMB(formula = effectiveness~log2(distance+1) + uprootexcav, data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation (Fluttersbys method):
# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[7]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[7]])))
# QQline <- qq.line(resid(Cand.mod[[7]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[7]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_29):
# dat.sim <- simulate(Cand.mod[[7]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[7]]))

# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[7]])
# performance::check_collinearity(Cand.mod[[7]])
# performance::r2_nakagawa(Cand.mod[[7]])
#
# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[7]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[7]])
# family(Cand.mod[[7]])$linkinv(glmmTMB::fixef(Cand.mod[[7]])$cond) # To get interpretable coefficients.

### Computing a Pseudo-R2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[7]], ~1))[1] - logLik(Cand.mod[[7]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=7, R2=R2[1]))



##### Model 8 #####
# -----------------

Cand.mod[[8]] <- glmmTMB::glmmTMB(formula = effectiveness~followups + fully_tarped, data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation (Fluttersbys method):
# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[8]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[8]])))
# QQline <- qq.line(resid(Cand.mod[[8]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[8]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_29):
# dat.sim <- simulate(Cand.mod[[8]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[8]]))

# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[8]])
# performance::check_collinearity(Cand.mod[[8]])
# performance::r2_nakagawa(Cand.mod[[8]])
#
# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[8]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[8]])
# family(Cand.mod[[8]])$linkinv(glmmTMB::fixef(Cand.mod[[8]])$cond) # To get interpretable coefficients.

### Computing a Pseudo-R2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[8]], ~1))[1] - logLik(Cand.mod[[8]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=8, R2=R2[1]))



##### Model 9 #####
# -----------------

Cand.mod[[9]] <- glmmTMB::glmmTMB(formula = effectiveness~followups + pb_fixation, data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation (Fluttersbys method):
# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[9]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[9]])))
# QQline <- qq.line(resid(Cand.mod[[9]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[9]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_29):
# dat.sim <- simulate(Cand.mod[[9]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[9]]))

# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[9]])
# performance::check_collinearity(Cand.mod[[9]])
# performance::r2_nakagawa(Cand.mod[[9]])
#
# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[9]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[9]])
# family(Cand.mod[[9]])$linkinv(glmmTMB::fixef(Cand.mod[[9]])$cond) # To get interpretable coefficients.

### Computing a Pseudo-R2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[9]], ~1))[1] - logLik(Cand.mod[[9]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=9, R2=R2[1]))



##### Model 10 #####
# -----------------

Cand.mod[[10]] <- glmmTMB::glmmTMB(formula = effectiveness~followups + obstacles, data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation (Fluttersbys method):
# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[10]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[10]])))
# QQline <- qq.line(resid(Cand.mod[[10]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[10]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_29):
# dat.sim <- simulate(Cand.mod[[10]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[10]]))

# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[10]])
# performance::check_collinearity(Cand.mod[[10]])
# performance::r2_nakagawa(Cand.mod[[10]])
#
# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[10]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[10]])
# family(Cand.mod[[10]])$linkinv(glmmTMB::fixef(Cand.mod[[10]])$cond) # To get interpretable coefficients.

### Computing a Pseudo-R2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[10]], ~1))[1] - logLik(Cand.mod[[10]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=10, R2=R2[1]))



##### Model 11 #####
# -----------------

Cand.mod[[11]] <- glmmTMB::glmmTMB(formula = effectiveness~followups + log2(stand_surface), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation (Fluttersbys method):
# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[11]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[11]])))
# QQline <- qq.line(resid(Cand.mod[[11]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[11]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_29):
# dat.sim <- simulate(Cand.mod[[11]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[11]]))

# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[11]])
# performance::check_collinearity(Cand.mod[[11]])
# performance::r2_nakagawa(Cand.mod[[11]])
#
# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[11]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[11]])
# family(Cand.mod[[11]])$linkinv(glmmTMB::fixef(Cand.mod[[11]])$cond) # To get interpretable coefficients.

### Computing a Pseudo-R2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[11]], ~1))[1] - logLik(Cand.mod[[11]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=11, R2=R2[1]))



##### Model 12 #####
# -----------------

Cand.mod[[12]] <- glmmTMB::glmmTMB(formula = effectiveness~followups + log2(tarping_duration), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation (Fluttersbys method):
# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[12]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[12]])))
# QQline <- qq.line(resid(Cand.mod[[12]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[12]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_29):
# dat.sim <- simulate(Cand.mod[[12]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[12]]))

# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[12]])
# performance::check_collinearity(Cand.mod[[12]])
# performance::r2_nakagawa(Cand.mod[[12]])
#
# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[12]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[12]])
# family(Cand.mod[[12]])$linkinv(glmmTMB::fixef(Cand.mod[[12]])$cond) # To get interpretable coefficients.

### Computing a Pseudo-R2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[12]], ~1))[1] - logLik(Cand.mod[[12]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=12, R2=R2[1]))



##### Model 13 #####
# -----------------

Cand.mod[[13]] <- glmmTMB::glmmTMB(formula = effectiveness~followups_cent * pb_fixation, data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

# ### Model evaluation (Fluttersbys method):
# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[13]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[13]])))
# QQline <- qq.line(resid(Cand.mod[[13]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[13]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_17):
# dat.sim <- simulate(Cand.mod[[13]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[13]]))
#
# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[13]])
# performance::check_collinearity(Cand.mod[[13]])
# performance::r2_nakagawa(Cand.mod[[13]])
#
# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[13]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[13]])
# family(Cand.mod[[13]])$linkinv(glmmTMB::fixef(Cand.mod[[13]])$cond) # To get interpretable coefficients.

### Computing a Pseudo-R2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[13]], ~1))[1] - logLik(Cand.mod[[13]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=13, R2=R2[1]))



##### Model 14 #####
# -----------------

Cand.mod[[14]] <- glmmTMB::glmmTMB(formula = effectiveness~fully_tarped + plantation, data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation (Fluttersbys method):
# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[14]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[14]])))
# QQline <- qq.line(resid(Cand.mod[[14]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[14]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_17):
# dat.sim <- simulate(Cand.mod[[14]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[14]]))

# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[14]])
# performance::check_collinearity(Cand.mod[[14]])
# performance::r2_nakagawa(Cand.mod[[14]])
#
# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[14]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[14]])
# family(Cand.mod[[14]])$linkinv(glmmTMB::fixef(Cand.mod[[14]])$cond) # To get interpretable coefficients.

### Computing a Pseudo-R2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[14]], ~1))[1] - logLik(Cand.mod[[14]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=14, R2=R2[1]))



##### Model 15 #####
# -----------------

Cand.mod[[15]] <- glmmTMB::glmmTMB(formula = effectiveness~fully_tarped + log2(sedicover_height+1), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation (Fluttersbys method):
# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[15]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[15]])))
# QQline <- qq.line(resid(Cand.mod[[15]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[15]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_17):
# dat.sim <- simulate(Cand.mod[[15]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[15]]))
#
# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[15]])
# performance::check_collinearity(Cand.mod[[15]])
# performance::r2_nakagawa(Cand.mod[[15]])
#
# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[15]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[15]])
# family(Cand.mod[[15]])$linkinv(glmmTMB::fixef(Cand.mod[[15]])$cond) # To get interpretable coefficients.

### Computing a Pseudo-R2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[15]], ~1))[1] - logLik(Cand.mod[[15]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=15, R2=R2[1]))



##### Model 16 #####
# -----------------

Cand.mod[[16]] <- glmmTMB::glmmTMB(formula = effectiveness~fully_tarped + log2(stand_surface), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation (Fluttersbys method):
# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[16]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[16]])))
# QQline <- qq.line(resid(Cand.mod[[16]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[16]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_17):
# dat.sim <- simulate(Cand.mod[[16]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[16]]))

# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[16]])
# performance::check_collinearity(Cand.mod[[16]])
# performance::r2_nakagawa(Cand.mod[[16]])
#
# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[16]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[16]])
# family(Cand.mod[[16]])$linkinv(glmmTMB::fixef(Cand.mod[[16]])$cond) # To get interpretable coefficients.

### Computing a Pseudo-R2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[16]], ~1))[1] - logLik(Cand.mod[[16]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=16, R2=R2[1]))



##### Model 17 #####
# -----------------

Cand.mod[[17]] <- glmmTMB::glmmTMB(formula = effectiveness~fully_tarped + obstacles, data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation (Fluttersbys method):
# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[17]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[17]])))
# QQline <- qq.line(resid(Cand.mod[[17]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[17]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_17):
# dat.sim <- simulate(Cand.mod[[17]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[17]]))
#
# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[17]])
# performance::check_collinearity(Cand.mod[[17]])
# performance::r2_nakagawa(Cand.mod[[17]])

# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[17]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[17]])
# family(Cand.mod[[17]])$linkinv(glmmTMB::fixef(Cand.mod[[17]])$cond) # To get interpretable coefficients.

### Computing a Pseudo-R2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[17]], ~1))[1] - logLik(Cand.mod[[17]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=17, R2=R2[1]))



##### Model 18 #####
# -----------------

Cand.mod[[18]] <- glmmTMB::glmmTMB(formula = effectiveness~fully_tarped + log2(tarping_duration), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation (Fluttersbys method):
# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[18]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[18]])))
# QQline <- qq.line(resid(Cand.mod[[18]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[18]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_17):
# dat.sim <- simulate(Cand.mod[[18]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[18]]))
#
# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[18]])
# performance::check_collinearity(Cand.mod[[18]])
# performance::r2_nakagawa(Cand.mod[[18]])

# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[18]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[18]])
# family(Cand.mod[[18]])$linkinv(glmmTMB::fixef(Cand.mod[[18]])$cond) # To get interpretable coefficients.

### Computing a Pseudo-R2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[18]], ~1))[1] - logLik(Cand.mod[[18]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=18, R2=R2[1]))



##### Model 19 #####
# -----------------

Cand.mod[[19]] <- glmmTMB::glmmTMB(formula = effectiveness~log2(stand_surface) + geomem, data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation (Fluttersbys method):
# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[19]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[19]])))
# QQline <- qq.line(resid(Cand.mod[[19]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[19]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_17):
# dat.sim <- simulate(Cand.mod[[19]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[19]]))

# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[19]])
# performance::check_collinearity(Cand.mod[[19]])
# performance::r2_nakagawa(Cand.mod[[19]])
#
# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[19]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[19]])
# family(Cand.mod[[19]])$linkinv(glmmTMB::fixef(Cand.mod[[19]])$cond) # To get interpretable coefficients.

### Computing a Pseudo-R2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[19]], ~1))[1] - logLik(Cand.mod[[19]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=19, R2=R2[1]))



##### Model 20 #####
# -----------------

Cand.mod[[20]] <- glmmTMB::glmmTMB(formula = effectiveness~log2(stand_surface) + obstacles, data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation (Fluttersbys method):
# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[20]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[20]])))
# QQline <- qq.line(resid(Cand.mod[[20]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[20]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_17):
# dat.sim <- simulate(Cand.mod[[20]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[20]]))

# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[20]])
# performance::check_collinearity(Cand.mod[[20]])
# performance::r2_nakagawa(Cand.mod[[20]])
#
# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[20]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[20]])
# family(Cand.mod[[20]])$linkinv(glmmTMB::fixef(Cand.mod[[20]])$cond) # To get interpretable coefficients.

### Computing a Pseudo-R2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[20]], ~1))[1] - logLik(Cand.mod[[20]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=20, R2=R2[1]))



##### Model 21 #####
# -----------------

Cand.mod[[21]] <- glmmTMB::glmmTMB(formula = effectiveness~log2(stand_surface) + plantation, data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation (Fluttersbys method):
# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[21]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[21]])))
# QQline <- qq.line(resid(Cand.mod[[21]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[21]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_17):
# dat.sim <- simulate(Cand.mod[[21]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[21]]))

# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[21]])
# performance::check_collinearity(Cand.mod[[21]])
# performance::r2_nakagawa(Cand.mod[[21]])
#
# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[21]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[21]])
# family(Cand.mod[[21]])$linkinv(glmmTMB::fixef(Cand.mod[[21]])$cond) # To get interpretable coefficients.

### Computing a Pseudo-R2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[21]], ~1))[1] - logLik(Cand.mod[[21]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=21, R2=R2[1]))



##### Model 22 #####
# -----------------

Cand.mod[[22]] <- glmmTMB::glmmTMB(formula = effectiveness~log2(stand_surface) + slope, data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation (Fluttersbys method):
# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[22]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[22]])))
# QQline <- qq.line(resid(Cand.mod[[22]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[22]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_17):
# dat.sim <- simulate(Cand.mod[[22]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[22]]))

# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[22]])
# performance::check_collinearity(Cand.mod[[22]])
# performance::r2_nakagawa(Cand.mod[[22]])
#
# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[22]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[22]])
# family(Cand.mod[[22]])$linkinv(glmmTMB::fixef(Cand.mod[[22]])$cond) # To get interpretable coefficients.

### Computing a Pseudo-R2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[22]], ~1))[1] - logLik(Cand.mod[[22]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=22, R2=R2[1]))



##### Model 23 #####
# -----------------

Cand.mod[[23]] <- glmmTMB::glmmTMB(formula = effectiveness~log2(stand_surface) + log2(tarping_duration), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation (Fluttersbys method):
# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[23]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[23]])))
# QQline <- qq.line(resid(Cand.mod[[23]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[23]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_17):
# dat.sim <- simulate(Cand.mod[[23]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[23]]))

# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[23]])
# performance::check_collinearity(Cand.mod[[23]])
# performance::r2_nakagawa(Cand.mod[[23]])
#
# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[23]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[23]])
# family(Cand.mod[[23]])$linkinv(glmmTMB::fixef(Cand.mod[[23]])$cond) # To get interpretable coefficients.

### Computing a Pseudo-R2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[23]], ~1))[1] - logLik(Cand.mod[[23]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=23, R2=R2[1]))



##### Model 24 #####
# -----------------

Cand.mod[[24]] <- glmmTMB::glmmTMB(formula = effectiveness~slope + plantation, data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation (Fluttersbys method):
# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[24]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[24]])))
# QQline <- qq.line(resid(Cand.mod[[24]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[24]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_17):
# dat.sim <- simulate(Cand.mod[[24]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[24]]))
#
# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[24]])
# performance::check_collinearity(Cand.mod[[24]])
# performance::r2_nakagawa(Cand.mod[[24]])

# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[24]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[24]])
# family(Cand.mod[[24]])$linkinv(glmmTMB::fixef(Cand.mod[[24]])$cond) # To get interpretable coefficients.

### Computing a Pseudo-R2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[24]], ~1))[1] - logLik(Cand.mod[[24]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=24, R2=R2[1]))



##### Model 25 #####
# -----------------

Cand.mod[[25]] <- glmmTMB::glmmTMB(formula = effectiveness~uprootexcav + geomem, data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation (Fluttersbys method):
# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[25]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[25]])))
# QQline <- qq.line(resid(Cand.mod[[25]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[25]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_17):
# dat.sim <- simulate(Cand.mod[[25]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[25]]))

# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[25]])
# performance::check_collinearity(Cand.mod[[25]])
# performance::r2_nakagawa(Cand.mod[[25]])
#
# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[25]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[25]])
# family(Cand.mod[[25]])$linkinv(glmmTMB::fixef(Cand.mod[[25]])$cond) # To get interpretable coefficients.

### Computing a Pseudo-R2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[25]], ~1))[1] - logLik(Cand.mod[[25]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=25, R2=R2[1]))



##### Model 26 #####
# -----------------

Cand.mod[[26]] <- glmmTMB::glmmTMB(formula = effectiveness~uprootexcav + plantation, data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation (Fluttersbys method):
# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[26]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[26]])))
# QQline <- qq.line(resid(Cand.mod[[26]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[26]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_17):
# dat.sim <- simulate(Cand.mod[[26]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[26]]))

# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[26]])
# performance::check_collinearity(Cand.mod[[26]])
# performance::r2_nakagawa(Cand.mod[[26]])
#
# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[26]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[26]])
# family(Cand.mod[[26]])$linkinv(glmmTMB::fixef(Cand.mod[[26]])$cond) # To get interpretable coefficients.

### Computing a Pseudo-R2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[26]], ~1))[1] - logLik(Cand.mod[[26]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=26, R2=R2[1]))



##### Model 27 #####
# -----------------

Cand.mod[[27]] <- glmmTMB::glmmTMB(formula = effectiveness~obstacles + slope, data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation (Fluttersbys method):
# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[27]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[27]])))
# QQline <- qq.line(resid(Cand.mod[[27]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[27]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_17):
# dat.sim <- simulate(Cand.mod[[27]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[27]]))

# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[27]])
# performance::check_collinearity(Cand.mod[[27]])
# performance::r2_nakagawa(Cand.mod[[27]])
#
# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[27]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[27]])
# family(Cand.mod[[27]])$linkinv(glmmTMB::fixef(Cand.mod[[27]])$cond) # To get interpretable coefficients.

### Computing a Pseudo-R2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[27]], ~1))[1] - logLik(Cand.mod[[27]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=27, R2=R2[1]))



##### Model 28 #####
# -----------------

Cand.mod[[28]] <- glmmTMB::glmmTMB(formula = effectiveness~fully_tarped * st_surface_cent, data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

# ### Model evaluation (Fluttersbys method):
# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[28]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[28]])))
# QQline <- qq.line(resid(Cand.mod[[28]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[28]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_17):
# dat.sim <- simulate(Cand.mod[[28]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[28]]))
#
# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[28]])
# performance::check_collinearity(Cand.mod[[28]])
# performance::r2_nakagawa(Cand.mod[[28]])
#
# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[28]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[28]])
# family(Cand.mod[[28]])$linkinv(glmmTMB::fixef(Cand.mod[[28]])$cond) # To get interpretable coefficients.

### Computing a Pseudo-R2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[28]], ~1))[1] - logLik(Cand.mod[[28]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=28, R2=R2[1]))



##### Model 29 #####
# -----------------

Cand.mod[[29]] <- glmmTMB::glmmTMB(formula = effectiveness~log2(distance+1) + fully_tarped + followups, data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation (Fluttersbys method):
# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[29]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[29]])))
# QQline <- qq.line(resid(Cand.mod[[29]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[29]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_17):
# dat.sim <- simulate(Cand.mod[[29]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[29]]))

# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[29]])
# performance::check_collinearity(Cand.mod[[29]])
# performance::r2_nakagawa(Cand.mod[[29]])
#
# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[29]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[29]])
# family(Cand.mod[[29]])$linkinv(glmmTMB::fixef(Cand.mod[[29]])$cond) # To get interpretable coefficients.

### Computing a Pseudo-R2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[29]], ~1))[1] - logLik(Cand.mod[[29]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=29, R2=R2[1]))



##### Model 30 #####
# -----------------

Cand.mod[[30]] <- glmmTMB::glmmTMB(formula = effectiveness~log2(distance+1) + fully_tarped + geomem, data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation (Fluttersbys method):
# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[30]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[30]])))
# QQline <- qq.line(resid(Cand.mod[[30]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[30]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_17):
# dat.sim <- simulate(Cand.mod[[30]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[30]]))

# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[30]])
# performance::check_collinearity(Cand.mod[[30]])
# performance::r2_nakagawa(Cand.mod[[30]])
#
# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[30]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[30]])
# family(Cand.mod[[30]])$linkinv(glmmTMB::fixef(Cand.mod[[30]])$cond) # To get interpretable coefficients.

### Computing a Pseudo-R2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[30]], ~1))[1] - logLik(Cand.mod[[30]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=30, R2=R2[1]))



##### Model 31 #####
# -----------------

Cand.mod[[31]] <- glmmTMB::glmmTMB(formula = effectiveness~log2(distance+1) + fully_tarped + obstacles, data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation (Fluttersbys method):
# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[31]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[31]])))
# QQline <- qq.line(resid(Cand.mod[[31]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[31]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_17):
# dat.sim <- simulate(Cand.mod[[31]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[31]]))
#
# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[31]])
# performance::check_collinearity(Cand.mod[[31]])
# performance::r2_nakagawa(Cand.mod[[31]])
#
# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[31]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")

# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[31]])
# family(Cand.mod[[31]])$linkinv(glmmTMB::fixef(Cand.mod[[31]])$cond) # To get interpretable coefficients.

### Computing a Pseudo-R2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[31]], ~1))[1] - logLik(Cand.mod[[31]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=31, R2=R2[1]))



##### Model 32 #####
# -----------------

Cand.mod[[32]] <- glmmTMB::glmmTMB(formula = effectiveness~log2(distance+1) + fully_tarped + plantation, data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation (Fluttersbys method):
# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[32]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[32]])))
# QQline <- qq.line(resid(Cand.mod[[32]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[32]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_17):
# dat.sim <- simulate(Cand.mod[[32]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[32]]))

# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[32]])
# performance::check_collinearity(Cand.mod[[32]])
# performance::r2_nakagawa(Cand.mod[[32]])
#
# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[32]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[32]])
# family(Cand.mod[[32]])$linkinv(glmmTMB::fixef(Cand.mod[[32]])$cond) # To get interpretable coefficients.

### Computing a Pseudo-R2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[32]], ~1))[1] - logLik(Cand.mod[[32]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=32, R2=R2[1]))



##### Model 33 #####
# -----------------

Cand.mod[[33]] <- glmmTMB::glmmTMB(formula = effectiveness~log2(distance+1) + fully_tarped + slope, data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation (Fluttersbys method):
# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[33]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[33]])))
# QQline <- qq.line(resid(Cand.mod[[33]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[33]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_17):
# dat.sim <- simulate(Cand.mod[[33]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[33]]))

# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[33]])
# performance::check_collinearity(Cand.mod[[33]])
# performance::r2_nakagawa(Cand.mod[[33]])
#
# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[33]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[33]])
# family(Cand.mod[[33]])$linkinv(glmmTMB::fixef(Cand.mod[[33]])$cond) # To get interpretable coefficients.

### Computing a Pseudo-R2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[33]], ~1))[1] - logLik(Cand.mod[[33]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=33, R2=R2[1]))



##### Model 34 #####
# -----------------

Cand.mod[[34]] <- glmmTMB::glmmTMB(formula = effectiveness~log2(distance+1) + fully_tarped + log2(sedicover_height+1), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation (Fluttersbys method):
# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[34]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[34]])))
# QQline <- qq.line(resid(Cand.mod[[34]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[34]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_17):
# dat.sim <- simulate(Cand.mod[[34]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[34]]))

# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[34]])
# performance::check_collinearity(Cand.mod[[34]])
# performance::r2_nakagawa(Cand.mod[[34]])
#
# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[34]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[34]])
# family(Cand.mod[[34]])$linkinv(glmmTMB::fixef(Cand.mod[[34]])$cond) # To get interpretable coefficients.

### Computing a Pseudo-R2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[34]], ~1))[1] - logLik(Cand.mod[[34]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=34, R2=R2[1]))



##### Model 35 #####
# ------------------

Cand.mod[[35]] <- glmmTMB::glmmTMB(formula = effectiveness~log2(distance+1) + fully_tarped + log2(stand_surface), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation (Fluttersbys method):
# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[35]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[35]])))
# QQline <- qq.line(resid(Cand.mod[[35]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[35]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_17):
# dat.sim <- simulate(Cand.mod[[35]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[35]]))

# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[35]])
# performance::check_collinearity(Cand.mod[[35]])
# performance::r2_nakagawa(Cand.mod[[35]])
#
# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[35]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[35]])
# family(Cand.mod[[35]])$linkinv(glmmTMB::fixef(Cand.mod[[35]])$cond) # To get interpretable coefficients.

### Computing a Pseudo-R2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[35]], ~1))[1] - logLik(Cand.mod[[35]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=35, R2=R2[1]))



##### Model 36 #####
# -----------------

Cand.mod[[36]] <- glmmTMB::glmmTMB(formula = effectiveness~log2(distance+1) + fully_tarped + log2(tarping_duration), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation (Fluttersbys method):
# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[36]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[36]])))
# QQline <- qq.line(resid(Cand.mod[[36]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[36]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_17):
# dat.sim <- simulate(Cand.mod[[36]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[36]]))

# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[36]])
# performance::check_collinearity(Cand.mod[[36]])
# performance::r2_nakagawa(Cand.mod[[36]])
#
# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[36]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[36]])
# family(Cand.mod[[36]])$linkinv(glmmTMB::fixef(Cand.mod[[36]])$cond) # To get interpretable coefficients.

### Computing a Pseudo-R2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[36]], ~1))[1] - logLik(Cand.mod[[36]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=36, R2=R2[1]))



##### Model 37 #####
# -----------------

Cand.mod[[37]] <- glmmTMB::glmmTMB(formula = effectiveness~log2(distance+1) + log2(stand_surface) + obstacles, data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation (Fluttersbys method):
# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[37]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[37]])))
# QQline <- qq.line(resid(Cand.mod[[37]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[37]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_17):
# dat.sim <- simulate(Cand.mod[[37]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[37]]))

# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[37]])
# performance::check_collinearity(Cand.mod[[37]])
# performance::r2_nakagawa(Cand.mod[[37]])
#
# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[37]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[37]])
# family(Cand.mod[[37]])$linkinv(glmmTMB::fixef(Cand.mod[[37]])$cond) # To get interpretable coefficients.

### Computing a Pseudo-R2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[37]], ~1))[1] - logLik(Cand.mod[[37]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=37, R2=R2[1]))



##### Model 38 #####
# -----------------

Cand.mod[[38]] <- glmmTMB::glmmTMB(formula = effectiveness~log2(distance+1) + log2(stand_surface) + pb_fixation, data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation (Fluttersbys method):
# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[38]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[38]])))
# QQline <- qq.line(resid(Cand.mod[[38]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[38]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_17):
# dat.sim <- simulate(Cand.mod[[38]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[38]]))

# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[38]])
# performance::check_collinearity(Cand.mod[[38]])
# performance::r2_nakagawa(Cand.mod[[38]])
#
# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[38]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[38]])
# family(Cand.mod[[38]])$linkinv(glmmTMB::fixef(Cand.mod[[38]])$cond) # To get interpretable coefficients.

### Computing a Pseudo-R2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[38]], ~1))[1] - logLik(Cand.mod[[38]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=38, R2=R2[1]))



##### Model 39 #####
# -----------------

Cand.mod[[39]] <- glmmTMB::glmmTMB(formula = effectiveness~fully_tarped + pb_fixation + followups, data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation (Fluttersbys method):
# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[39]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[39]])))
# QQline <- qq.line(resid(Cand.mod[[39]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[39]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_17):
# dat.sim <- simulate(Cand.mod[[39]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[39]]))

# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[39]])
# performance::check_collinearity(Cand.mod[[39]])
# performance::r2_nakagawa(Cand.mod[[39]])
#
# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[39]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[39]])
# family(Cand.mod[[39]])$linkinv(glmmTMB::fixef(Cand.mod[[39]])$cond) # To get interpretable coefficients.)

### Computing a Pseudo-R2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[39]], ~1))[1] - logLik(Cand.mod[[39]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=39, R2=R2[1]))



##### Model 40 #####
# -----------------

Cand.mod[[40]] <- glmmTMB::glmmTMB(formula = effectiveness~fully_tarped + log2(stand_surface) + geomem, data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation (Fluttersbys method):
# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[40]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[40]])))
# QQline <- qq.line(resid(Cand.mod[[40]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[40]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_17):
# dat.sim <- simulate(Cand.mod[[40]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[40]]))

# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[40]])
# performance::check_collinearity(Cand.mod[[40]])
# performance::r2_nakagawa(Cand.mod[[40]])
#
# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[40]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[40]])
# family(Cand.mod[[40]])$linkinv(glmmTMB::fixef(Cand.mod[[40]])$cond) # To get interpretable coefficients.

### Computing a Pseudo-R2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[40]], ~1))[1] - logLik(Cand.mod[[40]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=40, R2=R2[1]))



##### Model 41 #####
# -----------------

Cand.mod[[41]] <- glmmTMB::glmmTMB(formula = effectiveness~fully_tarped + log2(stand_surface) + pb_fixation, data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation (Fluttersbys method):
# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[41]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[41]])))
# QQline <- qq.line(resid(Cand.mod[[41]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[41]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_17):
# dat.sim <- simulate(Cand.mod[[41]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[41]]))

# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[41]])
# performance::check_collinearity(Cand.mod[[41]])
# performance::r2_nakagawa(Cand.mod[[41]])
#
# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[41]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[41]])
# family(Cand.mod[[41]])$linkinv(glmmTMB::fixef(Cand.mod[[41]])$cond) # To get interpretable coefficients.

### Computing a Pseudo-R2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[41]], ~1))[1] - logLik(Cand.mod[[41]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=41, R2=R2[1]))



##### Model 42 #####
# -----------------

Cand.mod[[42]] <- glmmTMB::glmmTMB(formula = effectiveness~fully_tarped + log2(stand_surface) + log2(sedicover_height+1), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation (Fluttersbys method):
# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[42]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[42]])))
# QQline <- qq.line(resid(Cand.mod[[42]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[42]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_17):
# dat.sim <- simulate(Cand.mod[[42]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[42]]))

# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[42]])
# performance::check_collinearity(Cand.mod[[42]])
# performance::r2_nakagawa(Cand.mod[[42]])
#
# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[42]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[42]])
# family(Cand.mod[[42]])$linkinv(glmmTMB::fixef(Cand.mod[[42]])$cond) # To get interpretable coefficients.

### Computing a Pseudo-R2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[42]], ~1))[1] - logLik(Cand.mod[[42]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=42, R2=R2[1]))



##### Model 43 #####
# -----------------

Cand.mod[[43]] <- glmmTMB::glmmTMB(formula = effectiveness~fully_tarped + log2(stand_surface) + slope, data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation (Fluttersbys method):
# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[43]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[43]])))
# QQline <- qq.line(resid(Cand.mod[[43]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[43]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_17):
# dat.sim <- simulate(Cand.mod[[43]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[43]]))

# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[43]])
# performance::check_collinearity(Cand.mod[[43]])
# performance::r2_nakagawa(Cand.mod[[43]])
#
# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[43]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[43]])
# family(Cand.mod[[43]])$linkinv(glmmTMB::fixef(Cand.mod[[43]])$cond) # To get interpretable coefficients.

### Computing a Pseudo-R2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[43]], ~1))[1] - logLik(Cand.mod[[43]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=43, R2=R2[1]))



##### Model 44 #####
# -----------------

Cand.mod[[44]] <- glmmTMB::glmmTMB(formula = effectiveness~fully_tarped + log2(stand_surface) + uprootexcav, data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation (Fluttersbys method):
# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[44]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[44]])))
# QQline <- qq.line(resid(Cand.mod[[44]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[44]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_17):
# dat.sim <- simulate(Cand.mod[[44]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[44]]))

# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[44]])
# performance::check_collinearity(Cand.mod[[44]])
# performance::r2_nakagawa(Cand.mod[[44]])
#
# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[44]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[44]])
# family(Cand.mod[[44]])$linkinv(glmmTMB::fixef(Cand.mod[[44]])$cond) # To get interpretable coefficients.

### Computing a Pseudo-R2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[44]], ~1))[1] - logLik(Cand.mod[[44]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=44, R2=R2[1]))



##### Model 45 #####
# -----------------

Cand.mod[[45]] <- glmmTMB::glmmTMB(formula = effectiveness~fully_tarped + log2(stand_surface) + log2(tarping_duration), data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation (Fluttersbys method):
# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[45]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[45]])))
# QQline <- qq.line(resid(Cand.mod[[45]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[45]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_17):
# dat.sim <- simulate(Cand.mod[[45]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[45]]))

# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[45]])
# performance::check_collinearity(Cand.mod[[45]])
# performance::r2_nakagawa(Cand.mod[[45]])
#
# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[45]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[45]])
# family(Cand.mod[[45]])$linkinv(glmmTMB::fixef(Cand.mod[[45]])$cond) # To get interpretable coefficients.

### Computing a Pseudo-R2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[45]], ~1))[1] - logLik(Cand.mod[[45]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=45, R2=R2[1]))



##### Model 46 #####
# -----------------

Cand.mod[[46]] <- glmmTMB::glmmTMB(formula = effectiveness~fully_tarped + followups + obstacles, data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation (Fluttersbys method):
# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[46]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[46]])))
# QQline <- qq.line(resid(Cand.mod[[46]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[46]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_17):
# dat.sim <- simulate(Cand.mod[[46]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[46]]))

# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[46]])
# performance::check_collinearity(Cand.mod[[46]])
# performance::r2_nakagawa(Cand.mod[[46]])
#
# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[46]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[46]])
# family(Cand.mod[[46]])$linkinv(glmmTMB::fixef(Cand.mod[[46]])$cond) # To get interpretable coefficients.

### Computing a Pseudo-R2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[46]], ~1))[1] - logLik(Cand.mod[[46]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=46, R2=R2[1]))



##### Model 47 #####
# -----------------

Cand.mod[[47]] <- glmmTMB::glmmTMB(formula = effectiveness~fully_tarped + obstacles + pb_fixation, data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation (Fluttersbys method):
# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[47]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[47]])))
# QQline <- qq.line(resid(Cand.mod[[47]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[47]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_17):
# dat.sim <- simulate(Cand.mod[[47]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[47]]))

# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[47]])
# performance::check_collinearity(Cand.mod[[47]])
# performance::r2_nakagawa(Cand.mod[[47]])
#
# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[47]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[47]])
# family(Cand.mod[[47]])$linkinv(glmmTMB::fixef(Cand.mod[[47]])$cond) # To get interpretable coefficients.

### Computing a Pseudo-R2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[47]], ~1))[1] - logLik(Cand.mod[[47]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=47, R2=R2[1]))



##### Model 48 #####
# -----------------

Cand.mod[[48]] <- glmmTMB::glmmTMB(formula = effectiveness~fully_tarped + obstacles + slope, data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation (Fluttersbys method):
# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[48]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[48]])))
# QQline <- qq.line(resid(Cand.mod[[48]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[48]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_17):
# dat.sim <- simulate(Cand.mod[[48]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[48]]))

# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[48]])
# performance::check_collinearity(Cand.mod[[48]])
# performance::r2_nakagawa(Cand.mod[[48]])
#
# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[48]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[48]])
# family(Cand.mod[[48]])$linkinv(glmmTMB::fixef(Cand.mod[[48]])$cond) # To get interpretable coefficients.

### Computing a Pseudo-R2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[48]], ~1))[1] - logLik(Cand.mod[[48]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=48, R2=R2[1]))



##### Model 49 #####
# -----------------

Cand.mod[[49]] <- glmmTMB::glmmTMB(formula = effectiveness~fully_tarped + obstacles + uprootexcav, data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

### Model evaluation (Fluttersbys method):
# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[49]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[49]])))
# QQline <- qq.line(resid(Cand.mod[[49]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[49]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_17):
# dat.sim <- simulate(Cand.mod[[49]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[49]]))
#
# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[49]])
# performance::check_collinearity(Cand.mod[[49]])
# performance::r2_nakagawa(Cand.mod[[49]])
#
# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[49]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[49]])
# family(Cand.mod[[49]])$linkinv(glmmTMB::fixef(Cand.mod[[49]])$cond) # To get interpretable coefficients.

### Computing a Pseudo-R2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[49]], ~1))[1] - logLik(Cand.mod[[49]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=49, R2=R2[1]))



##### Model 50 #####
# -----------------

Cand.mod[[50]] <- glmmTMB::glmmTMB(formula = effectiveness~distance_cent * st_surface_cent, data = eff,
                                   family = glmmTMB::beta_family(link = "logit"), REML = FALSE)

# ### Model evaluation (Fluttersbys method):
# ggplot2::ggplot(data = NULL) + ggplot2::geom_point(ggplot2::aes(y = residuals(Cand.mod[[50]],
#                                                                               type = "pearson"), x = fitted(Cand.mod[[50]])))
# QQline <- qq.line(resid(Cand.mod[[50]], type = "pearson"))
# ggplot2::ggplot(data = NULL, ggplot2::aes(sample = resid(Cand.mod[[50]], type = "pearson"))) +
#   ggplot2::stat_qq() + ggplot2::geom_abline(intercept = QQline[1], slope = QQline[2])
# # To simulate residuals (see https://www.flutterbys.com.au/stats/tut/tut10.5a.html#h2_17):
# dat.sim <- simulate(Cand.mod[[50]], n = 250)
# par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
# resid <- NULL
# for (i in 1:nrow(dat.sim)) {
#   e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
#   resid[i] <- e(eff$effectiveness[i] + runif(250, -0.5, 0.5))
# }
# par(.pardefault) # To restore defaults graphical parameters
# plot(resid ~ fitted(Cand.mod[[50]]))
#
# ### Model evaluation (with the 'performance' package)
# performance::check_autocorrelation(Cand.mod[[50]])
# performance::check_collinearity(Cand.mod[[50]])
# performance::r2_nakagawa(Cand.mod[[50]])
#
# # To plot the observed vs. fitted values:
# plot(x = fitted(Cand.mod[[50]]), y = eff$effectiveness, xlab = "Fitted values", ylab = "Observed values")
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[50]])
# family(Cand.mod[[50]])$linkinv(glmmTMB::fixef(Cand.mod[[50]])$cond) # To get interpretable coefficients.

### Computing a Pseudo-R2:
R2 <- 1 - exp((2/nrow(eff)) * (logLik(update(Cand.mod[[50]], ~1))[1] - logLik(Cand.mod[[50]])[1])) # Methods from flutterbys.com
R.ajust <- rbind(R.ajust, data.frame(Model=50, R2=R2[1]))





# ------------------------------------- #
##### Model selection and averaging #####
# ------------------------------------- #

Model <- (1:50)

Candidate <- c("null",
               "log2(distance + 1) + followups", "log2(distance + 1) + fully_tarped", "log2(distance + 1) + obstacles",
               "log2(distance + 1) + log2(stand_surface)", "log2(distance + 1) + log2(sedicover_height + 1)", "log2(distance + 1) + uprootexcav",
               "followups + fully_tarped", "followups + pb_fixation", "followups + obstacles",
               "followups + log2(stand_surface)", "followups + tarping_duration", "followups * pb_fixation",
               "fully_tarped + plantation",
               "fully_tarped + log2(sedicover_height + 1)", # Ok n°15
               "fully_tarped + log2(stand_surface)",
               "fully_tarped + obstacles",
               "fully_tarped + tarping_duration",
               "log2(stand_surface) + geomem",
               "log2(stand_surface) + obstacles", # Ok n°20
               "log2(stand_surface) + plantation",
               "log2(stand_surface) + slope",
               "log2(stand_surface) + tarping_duration",
               "slope + plantation",
               "uprootexcav + geomem", # Ok n°25
               "uprootexcav + plantation",
               "obstacles + slope",
               "followups * plantation",
               "log2(distance + 1) + fully_tarped + followups",
               "log2(distance + 1) + fully_tarped + geomem", # Ok n°30
               "log2(distance + 1) + fully_tarped + obstacles",
               "log2(distance + 1) + fully_tarped + plantation",
               "log2(distance + 1) + fully_tarped + slope",
               "log2(distance + 1) + fully_tarped + log2(sedicover_height + 1)",
               "log2(distance + 1) + fully_tarped + log2(stand_surface)", # Ok n°35
               "log2(distance + 1) + fully_tarped + tarping_duration",
               "log2(distance + 1) + log2(stand_surface) + obstacles",
               "log2(distance + 1) + log2(stand_surface) + pb_fixation",
               "fully_tarped + pb_fixation + followups",
               "fully_tarped + log2(stand_surface) + geomem",  # Ok n°40
               "fully_tarped + log2(stand_surface) + pb_fixation",
               "fully_tarped + log2(stand_surface) + log2(sedicover_height + 1)",
               "fully_tarped + log2(stand_surface) + slope",
               "fully_tarped + log2(stand_surface) + uprootexcav",
               "fully_tarped + log2(stand_surface) + tarping_duration", # Ok n°45
               "fully_tarped + followups + obstacles",
               "fully_tarped + obstacles + pb_fixation",
               "fully_tarped + obstacles + slope",
               "fully_tarped + obstacles + uprootexcav",
               "log2(distance + 1) * log2(stand_surface)") # Ok n°50

Cand.model <- data.frame(Model, Candidate)



##### Model selection #####
# -------------------------

### Rank models based on AICc:
AICc <- MuMIn::model.sel(object = Cand.mod)

### Improved formating:
AICc.model <- as.data.frame(AICc)
AICc.model$csweigth <- cumsum(AICc.model$weight) # Add a column containing the cumulative sum of AICc weights
AICc.model$Model <- row.names(AICc.model) # Add a column containing the n° of each candidate model

AICc.model <- merge(AICc.model, Cand.model, by="Model") # CAUTION: this merger reorder the rows of the
# the table based on Model, and it removes the NULL model from the table (as it doesn't exist in Cand.model).
# Consequently, the new ordering begins at 10 (10, 11, 12, ..., 20, 21, ...)!

AICc.model <- merge(AICc.model, R.ajust, by="Model") # This merger adds the computed R2 (here: pseudo-R2)

AICc.model <- AICc.model[order(AICc.model$delta),] # To reorder rows according to delta AICc
AICc.model$Response <- "effectiveness"
AICc.model$Rank <- 1:nrow(AICc.model)
AICc.model <- AICc.model[,c("Rank", "Model", "Response", "Candidate", "df", "AICc", "delta", "weight", "R2")]
AICc.model$AICc <- format(round(AICc.model$AICc, digits=1))
AICc.model$delta <- format(round(AICc.model$delta, digits=3))
AICc.model$weight <- format(round(AICc.model$weight, digits=3))
AICc.model$R2 <- format(round(AICc.model$R2, digits=3))
colnames(AICc.model)[colnames(AICc.model) == 'Candidate'] <- 'Candidate model'
colnames(AICc.model)[colnames(AICc.model) == 'df'] <- 'k'
colnames(AICc.model)[colnames(AICc.model) == 'delta'] <- 'delta AICc'
colnames(AICc.model)[colnames(AICc.model) == 'weight'] <- 'W'
colnames(AICc.model)[colnames(AICc.model) == 'R2'] <- 'pseudo-R2'

### Table export:
readr::write_csv2(x = AICc.model, file = here::here("output", "tables", "Models_effectiveness.csv"))



##### Multimodel Inference (averaging) #####
# ------------------------------------------

Parameters <- c("Intercept", "log2(distance + 1)", "followups", "fully_tarped", "geomem", "obstacles",
                "plantation", "pb_fixation", "log2(sedicover_height + 1)", "slope", "log2(stand_surface)",
                "tarping_duration", "uprootexcav")
Var <- c("cond((Int))", "cond(log2(distance + 1))", "cond(followups)", "cond(fully_tarped1)", "cond(geomem1)",
         "cond(obstacles)", "cond(plantation1)", "cond(pb_fixation1)", "cond(log2(sedicover_height + 1))",
         "cond(slope)", "cond(log2(stand_surface))", "cond(log2(tarping_duration))", "cond(uprootexcav1)")
Var.Imp <- c("cond((Int))", "cond(log2(distance + 1))", "cond(followups)", "cond(fully_tarped)", "cond(geomem)",
             "cond(obstacles)", "cond(plantation)", "cond(pb_fixation)", "cond(log2(sedicover_height + 1))",
             "cond(slope)", "cond(log2(stand_surface))", "cond(log2(tarping_duration))", "cond(uprootexcav)")

Para.model <- data.frame(Parameters, Var, Var.Imp)

### Removing the models with interactions from the top-models:
AICc[-c(32,35,48),] -> AICc # Removing the 3 lines containing the 3 interaction models. I need to
# remove them to avoid biaising the coefficients as they, by definition, have coefficient that do not
# mean the same thing as for the other models.

### Select the top models:
#top.models <- MuMIn::get.models(AICc, cumsum(weight) <= 0.95) # To take those with a cumulative sum of
# AICc weights <= 0.95
top.models <- MuMIn::get.models(AICc, cumsum(weight) <= 1) # To take them all (we chose this option because
# it was equally interesting to highlight the explanatory variables whose "importance" was supported
# by the data and those that were not)!
# We could also select models according to their delta AICc (see Burnham & Anderson, 2002)!

### Actual model parameters averaging:
Parameter <- MuMIn::model.avg(top.models, revised.var=T, adjusted=T, fit=T)
Parameter.model <- as.data.frame(cbind(MuMIn::coefTable(Parameter), stats::confint(Parameter))) # Reports
# the conditional averaged parameters with their 95% confidence interval

### Improved formating:
Parameter.model$Var <- row.names(Parameter.model)
Parameter.model <- merge(Parameter.model, Para.model, by="Var", all=TRUE)
Imp <- as.data.frame(format(round(MuMIn::sw(Parameter), digits=2))) # Gives each predictor a
# relative "importance" value based on the sum of the model weights of the models in which the variable
# is included as a predictor
colnames(Imp) <- "Imp."
Imp$Var.Imp <- row.names(Imp)
Parameter.model <- merge(Parameter.model, Imp, by="Var.Imp", all=TRUE)

Parameter.model$'Estimate (±SE)' <- paste0(format(round(Parameter.model$Estimate, digits=3), trim=T),
                                           " (±", format(round(Parameter.model$'Std. Error', digits=3),
                                                         trim=T), ")")
Parameter.model$'(95% CI)' <- paste0("(", format(round(Parameter.model$'2.5 %', digits=3), trim=T),
                                     "; ", format(round(Parameter.model$'97.5 %', digits=3), trim=T), ")")
Parameter.model$Response <- "effectiveness"
Parameter.model$'N model' <- length(top.models)
Parameter.model <- Parameter.model[,c("Response", "Parameters", "Imp.", "Estimate (±SE)", "(95% CI)",
                                      'N model')]
Parameter.model <- Parameter.model[!is.na(Parameter.model$Imp.), ]

### Table export:
readr::write_csv2(x = Parameter.model, file = here::here("output", "tables", "Parameters_effectiveness.csv"))


