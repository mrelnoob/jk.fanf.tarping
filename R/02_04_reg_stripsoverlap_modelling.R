### Sourcing previous scripts:
source(file = here::here("R/01_03_exploration_functions.R")) # I need to source my previous scripts (R files)
# to be able to import my cleaned dataset (with my 'model_datasets()' function) and perform my analyses.



# ---------------------------------------------------------- #
##### Data preparation for modelling 'reg_stripsoverlap' #####
# ---------------------------------------------------------- #

# List of used packages (for publication or package building): here, readr, MuMIn, lme4, ggplot2,
# (broom.mixed), stats

.pardefault <- par() # To save the default graphical parameters (in case I want to restore them).

readr::read_csv(here::here("data", "roverlaps.csv"), col_names = TRUE, col_types =
                  readr::cols(
                    manager_id = readr::col_factor(),
                    xp_id = readr::col_factor(),
                    reg_stripsoverlap = readr::col_factor(c("0", "1")),
                    geomem = readr::col_factor(c("0", "1")),
                    levelling = readr::col_factor(c("0", "1")),
                    fully_tarped = readr::col_factor(c("0", "1")),
                    tarpfix_multimethod = readr::col_factor(c("0", "1")),
                    tarpfix_pierced = readr::col_factor(c("0", "1")),
                    stripfix_pierced = readr::col_factor(c("0", "1")),
                    stripfix_taped = readr::col_factor(c("0", "1")),
                    plantation = readr::col_factor(c("0", "1")),
                    repairs = readr::col_factor(c("0", "1")),
                    add_control = readr::col_factor(c("0", "1")),
                    pb_fixation = readr::col_factor(c("0", "1")),
                    pb_durability = readr::col_factor(c("0", "1")),
                    reg_elsewhere = readr::col_factor(c("0", "1")))) %>%
  dplyr::mutate(latitude = jitter(x = latitude, factor = 0.1)) %>%
  dplyr::mutate(longitude = jitter(x = longitude, factor = 0.1)) -> roverlaps # Added a very small amount of
# noise to coordinates to avoid groups with exactly similar coordinates (related to low Lat/Long resolution)
# which prevent the proper use of the DHARMa package autocorrelation test!
summary(roverlaps)



##### Final pre-modelling assumption checks ##### (run only when required)
# --------------------------------------------- #

# ### Testing the relevance of the random effect structure:
# m0.glm <- stats::glm(reg_stripsoverlap ~ sqrt(strips_overlap), data = roverlaps, family = binomial)
# m0.glmer <- lme4::glmer(reg_stripsoverlap ~ sqrt(strips_overlap) + (1|manager_id), data = roverlaps, family = binomial)
# m0.glmer2 <- lme4::glmer(reg_stripsoverlap ~ sqrt(strips_overlap) + (1|manager_id) + (1|xp_id), data = roverlaps, family = binomial)
# aic.glm <- AIC(logLik(m0.glm))
# aic.glmer <- AIC(logLik(m0.glmer))
# aic.glmer2 <- AIC(logLik(m0.glmer2))
#
# # Likelihood Ratio Test:
# null.id <- -2 * logLik(m0.glm) + 2 * logLik(m0.glmer)
# pchisq(as.numeric(null.id), df=1, lower.tail=F)
# null.id <- -2 * logLik(m0.glm) + 2 * logLik(m0.glmer2)
# pchisq(as.numeric(null.id), df=1, lower.tail=F)
# rm(m0.glm, m0.glmer, m0.glmer2, aic.glm, aic.glmer, aic.glmer2, null.id)
# The Likelihood Ratio Tests are NOT significant so the use of the random effect structure may not be
# necessary! However, further tests on the model residuals may indicate otherwise.


# ### (Re-)Assessing the linearity assumption:
# rov <- roverlaps[,c("reg_stripsoverlap", "fully_tarped", "levelling", "plantation", "obstacles",
#                     "sedicover_height", "stand_surface", "strips_overlap",
#                     "stripfix_taped", "tarpfix_multimethod")]
# model <- glm(reg_stripsoverlap ~., data = rov, family = binomial)
# # Predict the probability (p) of regrowths at strip overlaps:
# probabilities <- predict(model, type = "response")
# # Transforming predictors
# mydata <- rov %>%
#   dplyr::select_if(is.numeric) %>%
#   dplyr::mutate("stand_surface (log)" = log(stand_surface)) %>%
#   dplyr::mutate("sedicover_height (log+1)" = log(sedicover_height+1)) %>%
#   dplyr::mutate("strips_overlap (sqrt)" = sqrt(strips_overlap)) # These transformations were made to
# # linearize the relationships although they may not be necessary for sedicover_height and strips_overlap.
# predictors <- colnames(mydata)
# # Bind the logit and tidying the data for plot (ggplot2, so long format)
# mydata <- mydata %>%
#   dplyr::mutate(logit = log(probabilities/(1-probabilities))) %>%
#   tidyr::gather(key = "predictors", value = "predictor.value", -logit)
# # Create scatterplot
# ggplot2::ggplot(mydata, ggplot2::aes(y = logit, x = predictor.value))+
#   ggplot2::geom_point(size = 0.5, alpha = 0.5) +
#   ggplot2::geom_smooth(method = "loess") +
#   ggplot2::theme_bw() +
#   ggplot2::facet_wrap(~predictors, scales = "free_x")
#
#
# ### Assessing multicollinearity:
# car::vif(mod = model) # There is no signs of multicollinearity as all GVIF value are under 1.65.
# rm(rov, model, probabilities, predictors, mydata)






# ---------------------------------------- #
##### Building of the candidate models #####
# ---------------------------------------- #

Cand.mod <- list()
R.ajust <- data.frame(Model=integer(0), R2=numeric(0)) # Creates an empty data.frame with 2 variables



##### Model 1 (null model) #####
# ------------------------------

Cand.mod[[1]] <- glmmTMB::glmmTMB(formula = reg_stripsoverlap~1, data = roverlaps,
                                  family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[1]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = roverlaps$longitude, y = roverlaps$latitude, plot = TRUE) # Ok (if
# # random effects are included, otherwise Moran's I test is significant)!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[1]])
# performance::check_autocorrelation(Cand.mod[[1]])
# performance::check_collinearity(Cand.mod[[1]])
# performance::check_singularity(Cand.mod[[1]])
#
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[1]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[1]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[1]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[1]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[1]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[1]]) # Tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=1, R2=R2[[1]]))



##### Model 2 #####
# -----------------

Cand.mod[[2]] <- glmmTMB::glmmTMB(formula = reg_stripsoverlap~obstacles, data = roverlaps,
                                  family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[2]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$obstacles) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = roverlaps$longitude, y = roverlaps$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[2]])
# performance::check_autocorrelation(Cand.mod[[2]])
# performance::check_collinearity(Cand.mod[[2]])
# performance::check_singularity(Cand.mod[[2]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[2]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[2]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[2]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[2]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[2]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[2]]) # Tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=2, R2=R2[[1]]))



##### Model 3 #####
# -----------------

Cand.mod[[3]] <- glmmTMB::glmmTMB(formula = reg_stripsoverlap~fully_tarped, data = roverlaps,
                                  family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[3]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$fully_tarped) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = roverlaps$longitude, y = roverlaps$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[3]])
# performance::check_autocorrelation(Cand.mod[[3]])
# performance::check_collinearity(Cand.mod[[3]])
# performance::check_singularity(Cand.mod[[3]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[3]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[3]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[3]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[3]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[3]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[3]]) # Tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=3, R2=R2[[1]]))



##### Model 4 #####
# -----------------

Cand.mod[[4]] <- glmmTMB::glmmTMB(formula = reg_stripsoverlap~levelling, data = roverlaps,
                                  family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[4]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$levelling) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = roverlaps$longitude, y = roverlaps$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[4]])
# performance::check_autocorrelation(Cand.mod[[4]])
# performance::check_collinearity(Cand.mod[[4]])
# performance::check_singularity(Cand.mod[[4]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[4]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[4]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[4]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[4]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[4]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[4]]) # Tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=4, R2=R2[[1]]))



##### Model 5 #####
# -----------------

Cand.mod[[5]] <- glmmTMB::glmmTMB(formula = reg_stripsoverlap~plantation, data = roverlaps,
                                  family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[5]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$plantation) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = roverlaps$longitude, y = roverlaps$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[5]])
# performance::check_autocorrelation(Cand.mod[[5]])
# performance::check_collinearity(Cand.mod[[5]])
# performance::check_singularity(Cand.mod[[5]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[5]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[5]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[5]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[5]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[5]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[5]]) # Tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=5, R2=R2[[1]]))



##### Model 6 #####
# -----------------

Cand.mod[[6]] <- glmmTMB::glmmTMB(formula = reg_stripsoverlap~log2(sedicover_height + 1), data = roverlaps,
                                  family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[6]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$sedicover_height) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = roverlaps$longitude, y = roverlaps$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[6]])
# performance::check_autocorrelation(Cand.mod[[6]])
# performance::check_collinearity(Cand.mod[[6]])
# performance::check_singularity(Cand.mod[[6]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[6]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[6]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[6]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[6]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[6]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[6]]) # Tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=6, R2=R2[[1]]))



##### Model 7 #####
# -----------------

Cand.mod[[7]] <- glmmTMB::glmmTMB(formula = reg_stripsoverlap~log2(stand_surface), data = roverlaps,
                                  family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[7]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$stand_surface) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = roverlaps$longitude, y = roverlaps$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[7]])
# performance::check_autocorrelation(Cand.mod[[7]])
# performance::check_collinearity(Cand.mod[[7]])
# performance::check_singularity(Cand.mod[[7]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[7]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[7]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[7]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[7]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[7]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[7]]) # Tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=7, R2=R2[[1]]))



##### Model 8 #####
# -----------------

Cand.mod[[8]] <- glmmTMB::glmmTMB(formula = reg_stripsoverlap~strips_overlap, data = roverlaps,
                                  family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[8]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$strips_overlap) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = roverlaps$longitude, y = roverlaps$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[8]])
# performance::check_autocorrelation(Cand.mod[[8]])
# performance::check_collinearity(Cand.mod[[8]])
# performance::check_singularity(Cand.mod[[8]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[8]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[8]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[8]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[8]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[8]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[8]]) # Tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=8, R2=R2[[1]]))



##### Model 9 #####
# -----------------

Cand.mod[[9]] <- glmmTMB::glmmTMB(formula = reg_stripsoverlap~slope, data = roverlaps,
                                  family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[9]], n = 1000, plot = FALSE)
# plot(simu.resid) # Nope!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$slope) # Nope! However, considering the only marginal
# # improvement made by the inclusion of a random effect structure compared to the aggravation of overfitting
# # it induces, we decided to keep the model as is.
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = roverlaps$longitude, y = roverlaps$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[9]])
# performance::check_autocorrelation(Cand.mod[[9]])
# performance::check_collinearity(Cand.mod[[9]])
# performance::check_singularity(Cand.mod[[9]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[9]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[9]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[9]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[9]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[9]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[9]]) # Tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=9, R2=R2[[1]]))



##### Model 10 #####
# -----------------

Cand.mod[[10]] <- glmmTMB::glmmTMB(formula = reg_stripsoverlap~stripfix_taped, data = roverlaps,
                                  family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[10]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$stripfix_taped) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = roverlaps$longitude, y = roverlaps$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[10]])
# performance::check_autocorrelation(Cand.mod[[10]])
# performance::check_collinearity(Cand.mod[[10]])
# performance::check_singularity(Cand.mod[[10]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[10]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[10]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[10]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[10]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[10]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[10]]) # Tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=10, R2=R2[[1]]))



##### Model 4 #####
# -----------------

Cand.mod[[11]] <- glmmTMB::glmmTMB(formula = reg_stripsoverlap~tarpfix_multimethod, data = roverlaps,
                                  family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[11]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$tarpfix_multimethod) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = roverlaps$longitude, y = roverlaps$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[11]])
# performance::check_autocorrelation(Cand.mod[[11]])
# performance::check_collinearity(Cand.mod[[11]])
# performance::check_singularity(Cand.mod[[11]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[11]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[11]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[11]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[11]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[11]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[11]]) # Tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=11, R2=R2[[1]]))



##### Model 12 #####
# -----------------

Cand.mod[[12]] <- glmmTMB::glmmTMB(formula = reg_stripsoverlap~strips_overlap + obstacles, data = roverlaps,
                                  family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[12]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish! Note that the loess smoother are perhaps too sensible.
# DHARMa::plotResiduals(simu.resid, form = roverlaps$strips_overlap) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$obstacles) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = roverlaps$longitude, y = roverlaps$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[12]])
# performance::check_autocorrelation(Cand.mod[[12]])
# performance::check_collinearity(Cand.mod[[12]])
# performance::check_singularity(Cand.mod[[12]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[12]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[12]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[12]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[12]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[12]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[12]]) # Tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=12, R2=R2[[1]]))



##### Model 13 #####
# -----------------

Cand.mod[[13]] <- glmmTMB::glmmTMB(formula = reg_stripsoverlap~strips_overlap + fully_tarped, data = roverlaps,
                                  family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[13]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$strips_overlap) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$fully_tarped) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = roverlaps$longitude, y = roverlaps$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[13]])
# performance::check_autocorrelation(Cand.mod[[13]])
# performance::check_collinearity(Cand.mod[[13]])
# performance::check_singularity(Cand.mod[[13]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[13]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[13]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[13]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[13]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[13]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[13]]) # Tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=13, R2=R2[[1]]))



##### Model 14 #####
# -----------------

Cand.mod[[14]] <- glmmTMB::glmmTMB(formula = reg_stripsoverlap~strips_overlap + levelling, data = roverlaps,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[14]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$strips_overlap) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$levelling) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = roverlaps$longitude, y = roverlaps$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[14]])
# performance::check_autocorrelation(Cand.mod[[14]])
# performance::check_collinearity(Cand.mod[[14]])
# performance::check_singularity(Cand.mod[[14]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[14]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[14]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[14]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[14]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[14]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[14]]) # Tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=14, R2=R2[[1]]))



##### Model 15 #####
# -----------------

Cand.mod[[15]] <- glmmTMB::glmmTMB(formula = reg_stripsoverlap~strips_overlap + plantation, data = roverlaps,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[15]], n = 1000, plot = FALSE)
# plot(simu.resid) # Nope! Same as before.
# DHARMa::plotResiduals(simu.resid, form = roverlaps$strips_overlap) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$plantation) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = roverlaps$longitude, y = roverlaps$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[15]])
# performance::check_autocorrelation(Cand.mod[[15]])
# performance::check_collinearity(Cand.mod[[15]])
# performance::check_singularity(Cand.mod[[15]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[15]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[15]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[15]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[15]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[15]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[15]]) # Tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=15, R2=R2[[1]]))



##### Model 16 #####
# -----------------

Cand.mod[[16]] <- glmmTMB::glmmTMB(formula = reg_stripsoverlap~strips_overlap + log2(sedicover_height + 1), data = roverlaps,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[16]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$strips_overlap) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$sedicover_height) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = roverlaps$longitude, y = roverlaps$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[16]])
# performance::check_autocorrelation(Cand.mod[[16]])
# performance::check_collinearity(Cand.mod[[16]])
# performance::check_singularity(Cand.mod[[16]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[16]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[16]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[16]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[16]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[16]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[16]]) # Tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=16, R2=R2[[1]]))



##### Model 17 #####
# -----------------

Cand.mod[[17]] <- glmmTMB::glmmTMB(formula = reg_stripsoverlap~strips_overlap + log2(stand_surface), data = roverlaps,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[17]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$strips_overlap) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$stand_surface) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = roverlaps$longitude, y = roverlaps$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[17]])
# performance::check_autocorrelation(Cand.mod[[17]])
# performance::check_collinearity(Cand.mod[[17]])
# performance::check_singularity(Cand.mod[[17]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[17]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[17]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[17]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[17]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[17]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[17]]) # Tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=17, R2=R2[[1]]))



##### Model 18 #####
# -----------------

Cand.mod[[18]] <- glmmTMB::glmmTMB(formula = reg_stripsoverlap~strips_overlap + slope, data = roverlaps,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[18]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$strips_overlap) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$slope) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = roverlaps$longitude, y = roverlaps$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[18]])
# performance::check_autocorrelation(Cand.mod[[18]])
# performance::check_collinearity(Cand.mod[[18]])
# performance::check_singularity(Cand.mod[[18]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[18]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[18]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[18]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[18]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[18]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[18]]) # Tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=18, R2=R2[[1]]))



##### Model 19 #####
# -----------------

Cand.mod[[19]] <- glmmTMB::glmmTMB(formula = reg_stripsoverlap~strips_overlap + stripfix_taped, data = roverlaps,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[19]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$strips_overlap) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$stripfix_taped) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = roverlaps$longitude, y = roverlaps$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[19]])
# performance::check_autocorrelation(Cand.mod[[19]])
# performance::check_collinearity(Cand.mod[[19]])
# performance::check_singularity(Cand.mod[[19]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[19]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[19]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[19]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[19]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[19]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[19]]) # Tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=19, R2=R2[[1]]))



##### Model 20 #####
# -----------------

Cand.mod[[20]] <- glmmTMB::glmmTMB(formula = reg_stripsoverlap~strips_overlap + tarpfix_multimethod, data = roverlaps,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[20]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$strips_overlap) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$tarpfix_multimethod) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = roverlaps$longitude, y = roverlaps$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[20]])
# performance::check_autocorrelation(Cand.mod[[20]])
# performance::check_collinearity(Cand.mod[[20]])
# performance::check_singularity(Cand.mod[[20]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[20]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[20]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[20]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[20]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[20]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[20]]) # Tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=20, R2=R2[[1]]))



##### Model 21 #####
# -----------------

Cand.mod[[21]] <- glmmTMB::glmmTMB(formula = reg_stripsoverlap~obstacles + levelling, data = roverlaps,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[21]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$obstacles) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$levelling) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = roverlaps$longitude, y = roverlaps$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[21]])
# performance::check_autocorrelation(Cand.mod[[21]])
# performance::check_collinearity(Cand.mod[[21]])
# performance::check_singularity(Cand.mod[[21]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[21]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[21]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[21]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[21]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[21]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[21]]) # Tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=21, R2=R2[[1]]))



##### Model 22 #####
# -----------------

Cand.mod[[22]] <- glmmTMB::glmmTMB(formula = reg_stripsoverlap~obstacles + log2(stand_surface), data = roverlaps,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[22]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$obstacles) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$stand_surface) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = roverlaps$longitude, y = roverlaps$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[22]])
# performance::check_autocorrelation(Cand.mod[[22]])
# performance::check_collinearity(Cand.mod[[22]])
# performance::check_singularity(Cand.mod[[22]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[22]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[22]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[22]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[22]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[22]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[22]]) # Tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=22, R2=R2[[1]]))



##### Model 23 #####
# -----------------

Cand.mod[[23]] <- glmmTMB::glmmTMB(formula = reg_stripsoverlap~obstacles + slope, data = roverlaps,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[23]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$obstacles) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$slope) # Nope! But tolerable.
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = roverlaps$longitude, y = roverlaps$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[23]])
# performance::check_autocorrelation(Cand.mod[[23]])
# performance::check_collinearity(Cand.mod[[23]])
# performance::check_singularity(Cand.mod[[23]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[23]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[23]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[23]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[23]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[23]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[23]]) # Tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=23, R2=R2[[1]]))



##### Model 24 #####
# -----------------

Cand.mod[[24]] <- glmmTMB::glmmTMB(formula = reg_stripsoverlap~obstacles + stripfix_taped, data = roverlaps,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[24]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$obstacles) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$stripfix_taped) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = roverlaps$longitude, y = roverlaps$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[24]])
# performance::check_autocorrelation(Cand.mod[[24]])
# performance::check_collinearity(Cand.mod[[24]])
# performance::check_singularity(Cand.mod[[24]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[24]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[24]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[24]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[24]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[24]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[24]]) # Tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=24, R2=R2[[1]]))



##### Model 25 #####
# -----------------

Cand.mod[[25]] <- glmmTMB::glmmTMB(formula = reg_stripsoverlap~obstacles + tarpfix_multimethod, data = roverlaps,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[25]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$obstacles) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$tarpfix_multimethod) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = roverlaps$longitude, y = roverlaps$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[25]])
# performance::check_autocorrelation(Cand.mod[[25]])
# performance::check_collinearity(Cand.mod[[25]])
# performance::check_singularity(Cand.mod[[25]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[25]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[25]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[25]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[25]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[25]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[25]]) # Tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=25, R2=R2[[1]]))



##### Model 26 #####
# -----------------

Cand.mod[[26]] <- glmmTMB::glmmTMB(formula = reg_stripsoverlap~stripfix_taped + levelling, data = roverlaps,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[26]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$stripfix_taped) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$levelling) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = roverlaps$longitude, y = roverlaps$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[26]])
# performance::check_autocorrelation(Cand.mod[[26]])
# performance::check_collinearity(Cand.mod[[26]])
# performance::check_singularity(Cand.mod[[26]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[26]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[26]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[26]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[26]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[26]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[26]]) # Tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=26, R2=R2[[1]]))



##### Model 27 #####
# -----------------

Cand.mod[[27]] <- glmmTMB::glmmTMB(formula = reg_stripsoverlap~stripfix_taped + log2(sedicover_height + 1), data = roverlaps,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[27]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$stripfix_taped) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$sedicover_height) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = roverlaps$longitude, y = roverlaps$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[27]])
# performance::check_autocorrelation(Cand.mod[[27]])
# performance::check_collinearity(Cand.mod[[27]])
# performance::check_singularity(Cand.mod[[27]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[27]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[27]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[27]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[27]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[27]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[27]]) # Tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=27, R2=R2[[1]]))



##### Model 28 #####
# -----------------

Cand.mod[[28]] <- glmmTMB::glmmTMB(formula = reg_stripsoverlap~stripfix_taped + tarpfix_multimethod, data = roverlaps,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[28]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$stripfix_taped) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$tarpfix_multimethod) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = roverlaps$longitude, y = roverlaps$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[28]])
# performance::check_autocorrelation(Cand.mod[[28]])
# performance::check_collinearity(Cand.mod[[28]])
# performance::check_singularity(Cand.mod[[28]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[28]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[28]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[28]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[28]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[28]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[28]]) # Tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=28, R2=R2[[1]]))



##### Model 29 #####
# -----------------

Cand.mod[[29]] <- glmmTMB::glmmTMB(formula = reg_stripsoverlap~stripfix_taped + slope, data = roverlaps,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[29]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$stripfix_taped) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$slope) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = roverlaps$longitude, y = roverlaps$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[29]])
# performance::check_autocorrelation(Cand.mod[[29]])
# performance::check_collinearity(Cand.mod[[29]])
# performance::check_singularity(Cand.mod[[29]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[29]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[29]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[29]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[29]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[29]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[29]]) # Tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=29, R2=R2[[1]]))



##### Model 30 #####
# -----------------

Cand.mod[[30]] <- glmmTMB::glmmTMB(formula = reg_stripsoverlap~tarpfix_multimethod + fully_tarped, data = roverlaps,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[30]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$tarpfix_multimethod) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$fully_tarped) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = roverlaps$longitude, y = roverlaps$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[30]])
# performance::check_autocorrelation(Cand.mod[[30]])
# performance::check_collinearity(Cand.mod[[30]])
# performance::check_singularity(Cand.mod[[30]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[30]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[30]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[30]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[30]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[30]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[30]]) # Tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=30, R2=R2[[1]]))



##### Model 31 #####
# -----------------

Cand.mod[[31]] <- glmmTMB::glmmTMB(formula = reg_stripsoverlap~plantation + fully_tarped, data = roverlaps,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[31]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$plantation) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$fully_tarped) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = roverlaps$longitude, y = roverlaps$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[31]])
# performance::check_autocorrelation(Cand.mod[[31]])
# performance::check_collinearity(Cand.mod[[31]])
# performance::check_singularity(Cand.mod[[31]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[31]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[31]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[31]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[31]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[31]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[31]]) # Tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=31, R2=R2[[1]]))



##### Model 32 #####
# -----------------

Cand.mod[[32]] <- glmmTMB::glmmTMB(formula = reg_stripsoverlap~log2(sedicover_height + 1) + fully_tarped, data = roverlaps,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[32]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$sedicover_height) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$fully_tarped) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = roverlaps$longitude, y = roverlaps$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[32]])
# performance::check_autocorrelation(Cand.mod[[32]])
# performance::check_collinearity(Cand.mod[[32]])
# performance::check_singularity(Cand.mod[[32]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[32]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[32]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[32]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[32]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[32]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[32]]) # Tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=32, R2=R2[[1]]))



##### Model 33 #####
# -----------------

Cand.mod[[33]] <- glmmTMB::glmmTMB(formula = reg_stripsoverlap~log2(stand_surface) + fully_tarped, data = roverlaps,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[33]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$stand_surface) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$fully_tarped) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = roverlaps$longitude, y = roverlaps$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[33]])
# performance::check_autocorrelation(Cand.mod[[33]])
# performance::check_collinearity(Cand.mod[[33]])
# performance::check_singularity(Cand.mod[[33]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[33]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[33]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[33]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[33]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[33]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[33]]) # Tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=33, R2=R2[[1]]))



##### Model 34 #####
# -----------------

Cand.mod[[34]] <- glmmTMB::glmmTMB(formula = reg_stripsoverlap~log2(sedicover_height + 1) + levelling, data = roverlaps,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[34]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$sedicover_height) # Ok!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$levelling) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = roverlaps$longitude, y = roverlaps$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[34]])
# performance::check_autocorrelation(Cand.mod[[34]])
# performance::check_collinearity(Cand.mod[[34]])
# performance::check_singularity(Cand.mod[[34]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[34]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[34]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[34]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[34]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[34]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[34]]) # Tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=34, R2=R2[[1]]))



##### Model 35 #####
# -----------------

Cand.mod[[35]] <- glmmTMB::glmmTMB(formula = reg_stripsoverlap~plantation + log2(stand_surface), data = roverlaps,
                                   family = stats::binomial(link = "logit"), REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[35]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$plantation) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = roverlaps$stand_surface) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = roverlaps$longitude, y = roverlaps$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[35]])
# performance::check_autocorrelation(Cand.mod[[35]])
# performance::check_collinearity(Cand.mod[[35]])
# performance::check_singularity(Cand.mod[[35]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[35]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[35]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[35]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[35]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[35]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[35]]) # Tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=35, R2=R2[[1]]))





# ------------------------------------- #
##### Model selection and averaging #####
# ------------------------------------- #

Model <- (1:35)

Candidate <- c("null",
               "obstacles", "fully_tarped", "levelling", "plantation", "log2(sedicover_height + 1)",
               "log2(stand_surface)", "strips_overlap", "slope", "stripfix_taped", "tarpfix_multimethod",
               "strips_overlap + obstacles", "strips_overlap + fully_tarped", "strips_overlap + levelling",
               "strips_overlap + plantation", "strips_overlap + log2(sedicover_height + 1)",
               "strips_overlap + log2(stand_surface)", "strips_overlap + slope", "strips_overlap + stripfix_taped",
               "strips_overlap + tarpfix_multimethod", "obstacles + levelling", "obstacles + log2(stand_surface)",
               "obstacles + slope", "obstacles + stripfix_taped", "obstacles + tarpfix_multimethod",
               "stripfix_taped + levelling", "stripfix_taped + log2(sedicover_height + 1)",
               "stripfix_taped + tarpfix_multimethod", "stripfix_taped + slope",
               "tarpfix_multimethod + fully_tarped", "fully_tarped + plantation",
               "fully_tarped + log2(sedicover_height + 1)", "fully_tarped + log2(stand_surface)",
               "levelling + log2(sedicover_height + 1)", "plantation + log2(stand_surface)")

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
AICc.model$Response <- "reg_stripsoverlap"
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
readr::write_csv2(x = AICc.model, file = here::here("output", "tables", "Models_reg_stripsoverlap.csv"))


# ### Testing if the delta AICc between the "best" model and the null model is significantly large:
# m0.glm <- glmmTMB::glmmTMB(formula = reg_stripsoverlap~1, data = roverlaps,
#                            family = stats::binomial(link = "logit"), REML = FALSE)
# m1.glm <- glmmTMB::glmmTMB(formula = reg_stripsoverlap~log2(sedicover_height + 1) + levelling, data = roverlaps,
#                            family = stats::binomial(link = "logit"), REML = FALSE)
# aic.glm <- AIC(logLik(m0.glm))
# aic.glm1 <- AIC(logLik(m1.glm))
#
# # Likelihood Ratio Test:
# null.id <- -2 * logLik(m0.glm) + 2 * logLik(m1.glm)
# pchisq(as.numeric(null.id), df=1, lower.tail=F)
#
# rm(m0.glm, m1.glm, aic.glm, aic.glm1, null.id)
# # The Likelihood Ratio Tests is significant so we have support to consider that the delta AIC between the
# # top-ranked model and the "null" model is large enough to consider that the former better explains our
# # response variable than the latter!



##### Multimodel Inference (averaging) #####
# ------------------------------------------

Parameters <- c("Intercept", "fully_tarped", "obstacles", "levelling", "plantation",
                "tarpfix_multimethod", "slope", "log2(sedicover_height + 1)", "log2(stand_surface)",
                "strips_overlap", "stripfix_taped")
Var <- c("cond((Int))", "cond(fully_tarped1)", "cond(obstacles)", "cond(levelling1)",
         "cond(plantation1)", "cond(tarpfix_multimethod1)", "cond(slope)", "cond(log2(sedicover_height + 1))",
         "cond(log2(stand_surface))", "cond(strips_overlap)", "cond(stripfix_taped1)")
Var.Imp <- c("cond((Int))", "cond(fully_tarped)", "cond(obstacles)", "cond(levelling)",
             "cond(plantation)", "cond(tarpfix_multimethod)", "cond(slope)", "cond(log2(sedicover_height + 1))",
             "cond(log2(stand_surface))", "cond(strips_overlap)", "cond(stripfix_taped)")

Para.model <- data.frame(Parameters, Var, Var.Imp)

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
Parameter.model$Response <- "reg_stripsoverlap"
Parameter.model$'N model' <- length(top.models)
Parameter.model <- Parameter.model[,c("Response", "Parameters", "Imp.", "Estimate (±SE)", "(95% CI)",
                                      'N model')]
Parameter.model <- Parameter.model[!is.na(Parameter.model$Imp.), ]

### Table export:
readr::write_csv2(x = Parameter.model, file = here::here("output", "tables", "Parameters_reg_stripsoverlap.csv"))
