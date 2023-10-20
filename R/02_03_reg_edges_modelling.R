### Sourcing previous scripts:
source(file = here::here("R/01_03_exploration_functions.R")) # I need to source my previous scripts (R files)
# to be able to import my cleaned dataset (with my 'model_datasets()' function) and perform my analyses.



# --------------------------------------------------- #
##### Data preparation for modelling 'reg_edges' #####
# --------------------------------------------------- #

.pardefault <- par() # To save the default graphical parameters (in case I want to restore them).

readr::read_csv(here::here("data", "redges.csv"), col_names = TRUE, col_types =
                  readr::cols(
                    manager_id = readr::col_factor(),
                    xp_id = readr::col_factor(),
                    reg_edges = readr::col_factor(),
                    geomem = readr::col_factor(c("0", "1")),
                    maxveg = readr::col_factor(c("0", "1")),
                    uprootexcav = readr::col_factor(c("0", "1")),
                    stripsoverlap_ok = readr::col_factor(c("0", "1")),
                    tarpfix_multimethod = readr::col_factor(c("0", "1")),
                    tarpfix_pierced = readr::col_factor(c("0", "1")),
                    plantation = readr::col_factor(c("0", "1")),
                    repairs = readr::col_factor(c("0", "1")),
                    add_control = readr::col_factor(c("0", "1")),
                    pb_fixation = readr::col_factor(c("0", "1")),
                    pb_durability = readr::col_factor(c("0", "1")),
                    reg_elsewhere = readr::col_factor(c("0", "1")))) %>%
  dplyr::mutate(latitude = jitter(x = latitude, factor = 0.1)) %>%
  dplyr::mutate(longitude = jitter(x = longitude, factor = 0.1)) %>%
  dplyr::mutate(distance_cent = scale(log2(distance+1), center = TRUE, scale = FALSE),
                trench_cent = scale(log2(trench_depth+1), center = TRUE, scale = FALSE),
                st_surface_cent = scale(log2(stand_surface), center = TRUE, scale = FALSE)) -> redges # Added a very small amount of
# noise to coordinates to avoid groups with exactly similar coordinates (related to low Lat/Long resolution)
# which prevent the proper use of the DHARMa package autocorrelation test!
summary(redges)





# ---------------------------------------- #
##### Building of the candidate models #####
# ---------------------------------------- #

Cand.mod <- list()
R.ajust <- data.frame(Model=integer(0), R2=numeric(0)) # Creates an empty data.frame with 2 variables



# ### Testing the relevance of the random effect structure:
# m0.glm <- glmmTMB::glmmTMB(formula = reg_edges~1, data = redges,
#                            family = stats::binomial(link = "logit"), REML = FALSE)
# m0.glmer <- glmmTMB::glmmTMB(formula = reg_edges~(1|manager_id), data = redges,
#                              family = stats::binomial(link = "logit"), REML = FALSE)
# aic.glm <- AIC(logLik(m0.glm))
# aic.glmer <- AIC(logLik(m0.glmer))
#
# # Likelihood Ratio Test:
# null.id <- -2 * logLik(m0.glm) + 2 * logLik(m0.glmer)
# pchisq(as.numeric(null.id), df=1, lower.tail=F) # The Likelihood Ratio Test is NOT significant suggesting
# # that the use of the random effect structure is not necessary! HOWEVER, model diagnostics for subsequent
# # models have shown that failing to include "manager_id" as a random effect leads to model misspecification.
# # Consequently, and as initially planned, we included a random structure within our candiate models.
# rm(m0.glm, m0.glmer, aic.glm, aic.glmer, null.id)



##### Model 1 (null model) #####
# ------------------------------

Cand.mod[[1]] <- glmmTMB::glmmTMB(formula = reg_edges~1 + (1|manager_id), data = redges,
                                  family = stats::binomial(link = "logit"), REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[1]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[1]])
# performance::check_autocorrelation(Cand.mod[[1]])
# performance::check_collinearity(Cand.mod[[1]])
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
R2 <- performance::r2_tjur(Cand.mod[[1]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=1, R2=R2[[1]]))



##### Model 2 #####
# -----------------

Cand.mod[[2]] <- glmmTMB::glmmTMB(formula = reg_edges~log2(distance+1) + add_control + (1|manager_id), data = redges,
                                  family = stats::binomial(link = "logit"), REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[2]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$add_control) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$distance) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[2]])
# performance::check_autocorrelation(Cand.mod[[2]])
# performance::check_collinearity(Cand.mod[[2]])
# performance::check_singularity(Cand.mod[[2]]) # Singularity! Dealt with decreased tolerance for convergence
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
R2 <- performance::r2_tjur(Cand.mod[[2]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=2, R2=R2[[1]]))



##### Model 3 #####
# -----------------

Cand.mod[[3]] <- glmmTMB::glmmTMB(formula = reg_edges~log2(distance+1) + obstacles + (1|manager_id), data = redges,
                                   family = stats::binomial(link = "logit"), REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[3]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$obstacles) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$distance) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
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
R2 <- performance::r2_tjur(Cand.mod[[3]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=3, R2=R2[[1]]))



##### Model 4 #####
# -----------------

Cand.mod[[4]] <- glmmTMB::glmmTMB(formula = reg_edges~log2(distance+1) + repairs + (1|manager_id), data = redges,
                                   family = stats::binomial(link = "logit"), REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[4]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$repairs) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$distance) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[4]])
# performance::check_autocorrelation(Cand.mod[[4]])
# performance::check_collinearity(Cand.mod[[4]])
# performance::check_singularity(Cand.mod[[4]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[4]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[4]])) # Not ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[4]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[4]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[4]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[4]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=4, R2=R2[[1]]))



##### Model 5 #####
# -----------------

Cand.mod[[5]] <- glmmTMB::glmmTMB(formula = reg_edges~log2(distance+1) + pb_fixation + (1|manager_id), data = redges,
                                   family = stats::binomial(link = "logit"), REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[5]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$pb_fixation) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$distance) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[5]])
# performance::check_autocorrelation(Cand.mod[[5]])
# performance::check_collinearity(Cand.mod[[5]])
# performance::check_singularity(Cand.mod[[5]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[5]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[5]])) # Close
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[5]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[5]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[5]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[5]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=5, R2=R2[[1]]))



##### Model 6 #####
# -----------------

Cand.mod[[6]] <- glmmTMB::glmmTMB(formula = reg_edges~log2(distance+1) + log2(stand_surface) + (1|manager_id), data = redges,
                                   family = stats::binomial(link = "logit"), REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[6]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$stand_surface) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$distance) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[6]])
# performance::check_autocorrelation(Cand.mod[[6]])
# performance::check_collinearity(Cand.mod[[6]])
# performance::check_singularity(Cand.mod[[6]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[6]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[6]])) # Close
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[6]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[6]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[6]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[6]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=6, R2=R2[[1]]))




##### Model 7 #####
# -----------------

Cand.mod[[7]] <- glmmTMB::glmmTMB(formula = reg_edges~log2(distance+1) + geomem + (1|manager_id), data = redges,
                                   family = stats::binomial(link = "logit"), REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[7]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$geomem) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$distance) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[7]])
# performance::check_autocorrelation(Cand.mod[[7]])
# performance::check_collinearity(Cand.mod[[7]])
# performance::check_singularity(Cand.mod[[7]]) # Singularity! Dealt with decreased tolerance for convergence
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[7]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[7]])) # Close!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[7]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[7]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[7]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[7]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=7, R2=R2[[1]]))



##### Model 8 #####
# -----------------

Cand.mod[[8]] <- glmmTMB::glmmTMB(formula = reg_edges~log2(distance+1) + tarping_duration + (1|manager_id), data = redges,
                                   family = stats::binomial(link = "logit"), REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[8]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$tarping_duration) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$distance) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[8]])
# performance::check_autocorrelation(Cand.mod[[8]])
# performance::check_collinearity(Cand.mod[[8]])
# performance::check_singularity(Cand.mod[[8]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[8]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[8]])) # Close
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[8]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[8]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[8]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[8]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=8, R2=R2[[1]]))



##### Model 9 #####
# -----------------

Cand.mod[[9]] <- glmmTMB::glmmTMB(formula = reg_edges~log2(distance+1) + log2(trench_depth+1) + (1|manager_id), data = redges,
                                   family = stats::binomial(link = "logit"), REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[9]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$trench_depth) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$distance) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[9]])
# performance::check_autocorrelation(Cand.mod[[9]])
# performance::check_collinearity(Cand.mod[[9]])
# performance::check_singularity(Cand.mod[[9]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[9]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[9]])) # Close
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[9]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[9]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[9]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[9]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=9, R2=R2[[1]]))



##### Model 10 #####
# -----------------

Cand.mod[[10]] <- glmmTMB::glmmTMB(formula = reg_edges~distance_cent * trench_cent + (1|manager_id), data = redges,
                                   family = stats::binomial(link = "logit"), REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[10]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$trench_depth) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$distance) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[10]])
# performance::check_autocorrelation(Cand.mod[[10]])
# performance::check_collinearity(Cand.mod[[10]]) # Produces moderate correlation even though predictors are
# # centered (and/or scaled). This is unusual but not very problematic since standard errors remain acceptable
# performance::check_singularity(Cand.mod[[10]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[10]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[10]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[10]])
# Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[10]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[10]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[10]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=10, R2=R2[[1]]))



##### Model 11 #####
# -----------------

Cand.mod[[11]] <- glmmTMB::glmmTMB(formula = reg_edges~log2(trench_depth+1) + add_control + (1|manager_id), data = redges,
                                   family = stats::binomial(link = "logit"), REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[11]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$trench_depth) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$add_control) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[11]])
# performance::check_autocorrelation(Cand.mod[[11]])
# performance::check_collinearity(Cand.mod[[11]])
# performance::check_singularity(Cand.mod[[11]]) # Singularity! Dealt with decreased tolerance for convergence
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
R2 <- performance::r2_tjur(Cand.mod[[11]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=11, R2=R2[[1]]))



##### Model 12 #####
# -----------------

Cand.mod[[12]] <- glmmTMB::glmmTMB(formula = reg_edges~log2(trench_depth+1) + obstacles + (1|manager_id), data = redges,
                                   family = stats::binomial(link = "logit"), REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[12]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$trench_depth) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$obstacles) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
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
R2 <- performance::r2_tjur(Cand.mod[[12]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=12, R2=R2[[1]]))



##### Model 13 #####
# -----------------

Cand.mod[[13]] <- glmmTMB::glmmTMB(formula = reg_edges~log2(trench_depth+1) + log2(stand_surface) + (1|manager_id), data = redges,
                                   family = stats::binomial(link = "logit"), REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[13]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$trench_depth) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$stand_surface) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[13]])
# performance::check_autocorrelation(Cand.mod[[13]])
# performance::check_collinearity(Cand.mod[[13]])
# performance::check_singularity(Cand.mod[[13]]) # Singularity! Dealt with decreased tolerance for convergence
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
R2 <- performance::r2_tjur(Cand.mod[[13]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=13, R2=R2[[1]]))



##### Model 14 #####
# -----------------

Cand.mod[[14]] <- glmmTMB::glmmTMB(formula = reg_edges~log2(stand_surface) + add_control + (1|manager_id), data = redges,
                                   family = stats::binomial(link = "logit"), REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[14]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$stand_surface) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$add_control) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[14]])
# performance::check_autocorrelation(Cand.mod[[14]])
# performance::check_collinearity(Cand.mod[[14]])
# performance::check_singularity(Cand.mod[[14]]) # Singularity! Dealt with decreased tolerance for convergence
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
R2 <- performance::r2_tjur(Cand.mod[[14]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=14, R2=R2[[1]]))



##### Model 15 #####
# -----------------

Cand.mod[[15]] <- glmmTMB::glmmTMB(formula = reg_edges~log2(stand_surface) + obstacles + (1|manager_id), data = redges,
                                   family = stats::binomial(link = "logit"), REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[15]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$stand_surface) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$obstacles) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
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
R2 <- performance::r2_tjur(Cand.mod[[15]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=15, R2=R2[[1]]))



##### Model 16 #####
# -----------------

Cand.mod[[16]] <- glmmTMB::glmmTMB(formula = reg_edges~st_surface_cent * distance_cent + (1|manager_id), data = redges,
                                   family = stats::binomial(link = "logit"), REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[16]], n = 1000, plot = FALSE)
# plot(simu.resid) # Not ok! I fitted various alternative models with various packages and not red flags
# # appeared. Therefore, I'm not sure where the problem comes from... As all test are okay, I will leave
# # the model as is (it won't be among the influancial models anyway so it should not bias inferences too
# # badly).
# DHARMa::plotResiduals(simu.resid, form = redges$stand_surface) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$distance) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
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
# Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[16]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[16]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[16]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=16, R2=R2[[1]]))



##### Model 17 #####
# -----------------

Cand.mod[[17]] <- glmmTMB::glmmTMB(formula = reg_edges~log2(stand_surface) + tarping_duration + (1|manager_id), data = redges,
                                   family = stats::binomial(link = "logit"), REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[17]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$stand_surface) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$tarping_duration) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
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
R2 <- performance::r2_tjur(Cand.mod[[17]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=17, R2=R2[[1]]))



##### Model 18 #####
# -----------------

Cand.mod[[18]] <- glmmTMB::glmmTMB(formula = reg_edges~log2(trench_depth+1) + repairs + (1|manager_id), data = redges,
                                   family = stats::binomial(link = "logit"), REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[18]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$trench_depth) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$repairs) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[18]])
# performance::check_autocorrelation(Cand.mod[[18]])
# performance::check_collinearity(Cand.mod[[18]])
# performance::check_singularity(Cand.mod[[18]]) # Singularity! Dealt with decreased tolerance for convergence
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
R2 <- performance::r2_tjur(Cand.mod[[18]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=18, R2=R2[[1]]))



##### Model 19 #####
# -----------------

Cand.mod[[19]] <- glmmTMB::glmmTMB(formula = reg_edges~log2(trench_depth+1) + slope + (1|manager_id), data = redges,
                                   family = stats::binomial(link = "logit"), REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[19]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$trench_depth) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$slope) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[19]])
# performance::check_autocorrelation(Cand.mod[[19]])
# performance::check_collinearity(Cand.mod[[19]])
# performance::check_singularity(Cand.mod[[19]]) # Singularity! Dealt with decreased tolerance for convergence
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
R2 <- performance::r2_tjur(Cand.mod[[19]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=19, R2=R2[[1]]))



##### Model 20 #####
# -----------------

Cand.mod[[20]] <- glmmTMB::glmmTMB(formula = reg_edges~add_control + obstacles + (1|manager_id), data = redges,
                                   family = stats::binomial(link = "logit"), REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[20]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$add_control) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$obstacles) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
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
R2 <- performance::r2_tjur(Cand.mod[[20]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=20, R2=R2[[1]]))



##### Model 21 #####
# -----------------

Cand.mod[[21]] <- glmmTMB::glmmTMB(formula = reg_edges~add_control + pb_fixation + (1|manager_id), data = redges,
                                   family = stats::binomial(link = "logit"), REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[21]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$add_control) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$pb_fixation) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[21]])
# performance::check_autocorrelation(Cand.mod[[21]])
# performance::check_collinearity(Cand.mod[[21]])
# performance::check_singularity(Cand.mod[[21]]) # Singularity! Dealt by decreased tolerance for convergence
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
R2 <- performance::r2_tjur(Cand.mod[[21]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=21, R2=R2[[1]]))



##### Model 22 #####
# -----------------

Cand.mod[[22]] <- glmmTMB::glmmTMB(formula = reg_edges~add_control + slope + (1|manager_id), data = redges,
                                   family = stats::binomial(link = "logit"), REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[22]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$add_control) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$slope) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[22]])
# performance::check_autocorrelation(Cand.mod[[22]])
# performance::check_collinearity(Cand.mod[[22]])
# performance::check_singularity(Cand.mod[[22]]) # Singularity! Dealt with decreased tolerance for convergence
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
R2 <- performance::r2_tjur(Cand.mod[[22]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=22, R2=R2[[1]]))



##### Model 23 #####
# -----------------

Cand.mod[[23]] <- glmmTMB::glmmTMB(formula = reg_edges~add_control + tarping_duration + (1|manager_id), data = redges,
                                   family = stats::binomial(link = "logit"), REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[23]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$add_control) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$tarping_duration) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[23]])
# performance::check_autocorrelation(Cand.mod[[23]])
# performance::check_collinearity(Cand.mod[[23]])
# performance::check_singularity(Cand.mod[[23]]) # Singularity! Dealt with decreased tolerance for convergence
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
R2 <- performance::r2_tjur(Cand.mod[[23]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=23, R2=R2[[1]]))



##### Model 24 #####
# ------------------

Cand.mod[[24]] <- glmmTMB::glmmTMB(formula = reg_edges~obstacles + pb_fixation + (1|manager_id), data = redges,
                                   family = stats::binomial(link = "logit"), REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[24]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$obstacles) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$pb_fixation) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
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
R2 <- performance::r2_tjur(Cand.mod[[24]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=24, R2=R2[[1]]))



##### Model 25 #####
# -----------------

Cand.mod[[25]] <- glmmTMB::glmmTMB(formula = reg_edges~obstacles + geomem + (1|manager_id), data = redges,
                                   family = stats::binomial(link = "logit"), REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[25]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$obstacles) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$geomem) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[25]])
# performance::check_autocorrelation(Cand.mod[[25]])
# performance::check_collinearity(Cand.mod[[25]])
# performance::check_singularity(Cand.mod[[25]]) # Singularity! Dealt with decreased tolerance for convergence
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
R2 <- performance::r2_tjur(Cand.mod[[25]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=25, R2=R2[[1]]))



##### Model 26 #####
# -----------------

Cand.mod[[26]] <- glmmTMB::glmmTMB(formula = reg_edges~obstacles + repairs + (1|manager_id), data = redges,
                                   family = stats::binomial(link = "logit"), REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[26]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$obstacles) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$repairs) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
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
R2 <- performance::r2_tjur(Cand.mod[[26]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=26, R2=R2[[1]]))



##### Model 27 #####
# -----------------

Cand.mod[[27]] <- glmmTMB::glmmTMB(formula = reg_edges~slope + geomem + (1|manager_id), data = redges,
                                   family = stats::binomial(link = "logit"), REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[27]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$slope) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$geomem) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[27]])
# performance::check_autocorrelation(Cand.mod[[27]])
# performance::check_collinearity(Cand.mod[[27]])
# performance::check_singularity(Cand.mod[[27]]) # Singularity! Dealt with decreased tolerance for convergence
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
R2 <- performance::r2_tjur(Cand.mod[[27]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=27, R2=R2[[1]]))



##### Model 28 #####
# -----------------

Cand.mod[[28]] <- glmmTMB::glmmTMB(formula = reg_edges~slope + pb_fixation + (1|manager_id), data = redges,
                                   family = stats::binomial(link = "logit"), REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[28]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$slope) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$pb_fixation) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
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
R2 <- performance::r2_tjur(Cand.mod[[28]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=28, R2=R2[[1]]))



##### Model 29 #####
# -----------------

Cand.mod[[29]] <- glmmTMB::glmmTMB(formula = reg_edges~slope + tarping_duration + (1|manager_id), data = redges,
                                   family = stats::binomial(link = "logit"), REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[29]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$slope) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$tarping_duration) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
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
R2 <- performance::r2_tjur(Cand.mod[[29]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=29, R2=R2[[1]]))



##### Model 30 #####
# -----------------

Cand.mod[[30]] <- glmmTMB::glmmTMB(formula = reg_edges~log2(distance+1) + log2(trench_depth+1) + add_control
                                   + (1|manager_id), data = redges, family = stats::binomial(link = "logit"),
                                   REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[30]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$distance) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$trench_depth) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$add_control) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[30]])
# performance::check_autocorrelation(Cand.mod[[30]])
# performance::check_collinearity(Cand.mod[[30]])
# performance::check_singularity(Cand.mod[[30]]) # Singularity! Dealt with decreased tolerance for convergence
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
R2 <- performance::r2_tjur(Cand.mod[[30]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=30, R2=R2[[1]]))



##### Model 31 #####
# -----------------

Cand.mod[[31]] <- glmmTMB::glmmTMB(formula = reg_edges~log2(distance+1) + log2(trench_depth+1) + obstacles
                                   + (1|manager_id), data = redges, family = stats::binomial(link = "logit"),
                                   REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[31]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$distance) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$trench_depth) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$obstacles) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
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
R2 <- performance::r2_tjur(Cand.mod[[31]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=31, R2=R2[[1]]))



##### Model 32 #####
# -----------------

Cand.mod[[32]] <- glmmTMB::glmmTMB(formula = reg_edges~log2(distance+1) + log2(trench_depth+1) + geomem
                                   + (1|manager_id), data = redges, family = stats::binomial(link = "logit"),
                                   REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[32]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$distance) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$trench_depth) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$geomem) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[32]])
# performance::check_autocorrelation(Cand.mod[[32]])
# performance::check_collinearity(Cand.mod[[32]])
# performance::check_singularity(Cand.mod[[32]]) # Singularity! Dealt with decreased tolerance for convergence
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
R2 <- performance::r2_tjur(Cand.mod[[32]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=32, R2=R2[[1]]))



##### Model 33 #####
# -----------------

Cand.mod[[33]] <- glmmTMB::glmmTMB(formula = reg_edges~log2(distance+1) + log2(trench_depth+1) + log2(stand_surface)
                                   + (1|manager_id), data = redges, family = stats::binomial(link = "logit"),
                                   REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[33]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$distance) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$trench_depth) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$stand_surface) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
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
R2 <- performance::r2_tjur(Cand.mod[[33]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=33, R2=R2[[1]]))



##### Model 34 #####
# -----------------

Cand.mod[[34]] <- glmmTMB::glmmTMB(formula = reg_edges~log2(distance+1) + log2(trench_depth+1) + slope
                                   + (1|manager_id), data = redges, family = stats::binomial(link = "logit"),
                                   REML = FALSE)

# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[34]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$distance) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$trench_depth) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$slope) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
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
R2 <- performance::r2_tjur(Cand.mod[[34]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=34, R2=R2[[1]]))



##### Model 35 #####
# -----------------

Cand.mod[[35]] <- glmmTMB::glmmTMB(formula = reg_edges~log2(distance+1) + log2(stand_surface) + add_control
                                   + (1|manager_id), data = redges, family = stats::binomial(link = "logit"),
                                   REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[35]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$distance) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$stand_surface) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$add_control) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[35]])
# performance::check_autocorrelation(Cand.mod[[35]])
# performance::check_collinearity(Cand.mod[[35]])
# performance::check_singularity(Cand.mod[[35]]) # Singularity! Dealt with decreased tolerance for convergence
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[35]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[35]])) # Nope!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[35]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[35]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[35]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[35]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=35, R2=R2[[1]]))



##### Model 36 #####
# -----------------

Cand.mod[[36]] <- glmmTMB::glmmTMB(formula = reg_edges~log2(distance+1) + log2(stand_surface) + obstacles
                                   + (1|manager_id), data = redges, family = stats::binomial(link = "logit"),
                                   REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[36]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$distance) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$stand_surface) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$obstacles) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[36]])
# performance::check_autocorrelation(Cand.mod[[36]])
# performance::check_collinearity(Cand.mod[[36]])
# performance::check_singularity(Cand.mod[[36]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[36]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[36]])) # Close!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[36]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[36]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[36]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[36]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=36, R2=R2[[1]]))



##### Model 37 #####
# -----------------

Cand.mod[[37]] <- glmmTMB::glmmTMB(formula = reg_edges~log2(distance+1) + log2(stand_surface) + pb_fixation
                                   + (1|manager_id), data = redges, family = stats::binomial(link = "logit"),
                                   REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[37]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$distance) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$stand_surface) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$pb_fixation) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[37]])
# performance::check_autocorrelation(Cand.mod[[37]])
# performance::check_collinearity(Cand.mod[[37]])
# performance::check_singularity(Cand.mod[[37]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[37]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[37]])) # Close!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[37]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[37]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[37]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[37]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=37, R2=R2[[1]]))



##### Model 38 #####
# -----------------

Cand.mod[[38]] <- glmmTMB::glmmTMB(formula = reg_edges~log2(distance+1) + geomem + tarping_duration
                                   + (1|manager_id), data = redges, family = stats::binomial(link = "logit"),
                                   REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[38]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$distance) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$geomem) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$tarping_duration) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[38]])
# performance::check_autocorrelation(Cand.mod[[38]])
# performance::check_collinearity(Cand.mod[[38]])
# performance::check_singularity(Cand.mod[[38]]) # Singularity! Dealt with decreased tolerance for convergence
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[38]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[38]])) # Ok!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[38]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[38]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[38]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[38]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=38, R2=R2[[1]]))



##### Model 39 #####
# ------------------

Cand.mod[[39]] <- glmmTMB::glmmTMB(formula = reg_edges~log2(stand_surface) + log2(trench_depth+1) + obstacles
                                   + (1|manager_id), data = redges, family = stats::binomial(link = "logit"),
                                   REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[39]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$trench_depth) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$stand_surface) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$obstacles) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[39]])
# performance::check_autocorrelation(Cand.mod[[39]])
# performance::check_collinearity(Cand.mod[[39]])
# performance::check_singularity(Cand.mod[[39]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[39]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[39]])) # Close!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[39]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[39]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[39]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[39]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=39, R2=R2[[1]]))



##### Model 40 #####
# -----------------

Cand.mod[[40]] <- glmmTMB::glmmTMB(formula = reg_edges~log2(stand_surface) + log2(trench_depth+1) + slope
                                   + (1|manager_id), data = redges, family = stats::binomial(link = "logit"),
                                   REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[40]], n = 1000, plot = FALSE)
# plot(simu.resid) # Nope! Same as for model #24!
# DHARMa::plotResiduals(simu.resid, form = redges$trench_depth) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$stand_surface) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$slope) # Nope!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[40]])
# performance::check_autocorrelation(Cand.mod[[40]])
# performance::check_collinearity(Cand.mod[[40]])
# performance::check_singularity(Cand.mod[[40]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[40]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[40]])) # Close!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[40]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[40]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[40]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[40]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=40, R2=R2[[1]]))



##### Model 41 #####
# -----------------

Cand.mod[[41]] <- glmmTMB::glmmTMB(formula = reg_edges~log2(distance+1) + tarping_duration + add_control
                                   + (1|manager_id), data = redges, family = stats::binomial(link = "logit"),
                                   REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[41]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$distance) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$tarping_duration) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$add_control) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[41]])
# performance::check_autocorrelation(Cand.mod[[41]])
# performance::check_collinearity(Cand.mod[[41]])
# performance::check_singularity(Cand.mod[[41]]) # Singularity!
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[41]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[41]])) # Close!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[41]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[41]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[41]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[41]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=41, R2=R2[[1]]))



##### Model 42 #####
# -----------------

Cand.mod[[42]] <- glmmTMB::glmmTMB(formula = reg_edges~obstacles + repairs + add_control
                                   + (1|manager_id), data = redges, family = stats::binomial(link = "logit"),
                                   REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[42]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$obstacles) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$repairs) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$add_control) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[42]])
# performance::check_autocorrelation(Cand.mod[[42]])
# performance::check_collinearity(Cand.mod[[42]])
# performance::check_singularity(Cand.mod[[42]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[42]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[42]])) # Close!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[42]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[42]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[42]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[42]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=42, R2=R2[[1]]))



##### Model 43 #####
# -----------------

Cand.mod[[43]] <- glmmTMB::glmmTMB(formula = reg_edges~obstacles + add_control + geomem
                                   + (1|manager_id), data = redges, family = stats::binomial(link = "logit"),
                                   REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[43]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$obstacles) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$add_control) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$geomem) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[43]])
# performance::check_autocorrelation(Cand.mod[[43]])
# performance::check_collinearity(Cand.mod[[43]])
# performance::check_singularity(Cand.mod[[43]]) # Singularity!
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[43]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[43]])) # Close!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[43]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[43]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[43]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[43]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=43, R2=R2[[1]]))



##### Model 44 #####
# -----------------

Cand.mod[[44]] <- glmmTMB::glmmTMB(formula = reg_edges~obstacles + add_control + tarping_duration
                                   + (1|manager_id), data = redges, family = stats::binomial(link = "logit"),
                                   REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[44]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$obstacles) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$add_control) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$tarping_duration) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[44]])
# performance::check_autocorrelation(Cand.mod[[44]])
# performance::check_collinearity(Cand.mod[[44]])
# performance::check_singularity(Cand.mod[[44]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[44]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[44]])) # Close!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[44]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[44]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[44]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[44]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=44, R2=R2[[1]]))



##### Model 45 #####
# -----------------

Cand.mod[[45]] <- glmmTMB::glmmTMB(formula = reg_edges~obstacles + add_control + log2(trench_depth+1)
                                   + (1|manager_id), data = redges, family = stats::binomial(link = "logit"),
                                   REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[45]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$obstacles) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$add_control) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$trench_depth) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[45]])
# performance::check_autocorrelation(Cand.mod[[45]])
# performance::check_collinearity(Cand.mod[[45]])
# performance::check_singularity(Cand.mod[[45]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[45]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[45]])) # Close!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[45]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[45]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[45]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[45]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=45, R2=R2[[1]]))



##### Model 46 #####
# -----------------

Cand.mod[[46]] <- glmmTMB::glmmTMB(formula = reg_edges~obstacles + add_control + log2(stand_surface)
                                   + (1|manager_id), data = redges, family = stats::binomial(link = "logit"),
                                   REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[46]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$obstacles) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$add_control) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$stand_surface) # Ok!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[46]])
# performance::check_autocorrelation(Cand.mod[[46]])
# performance::check_collinearity(Cand.mod[[46]])
# performance::check_singularity(Cand.mod[[46]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[46]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[46]])) # Close!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[46]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[46]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[46]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[46]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=46, R2=R2[[1]]))



##### Model 47 #####
# -----------------

Cand.mod[[47]] <- glmmTMB::glmmTMB(formula = reg_edges~obstacles + add_control + slope
                                   + (1|manager_id), data = redges, family = stats::binomial(link = "logit"),
                                   REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[47]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$obstacles) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$add_control) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$slope) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[47]])
# performance::check_autocorrelation(Cand.mod[[47]])
# performance::check_collinearity(Cand.mod[[47]])
# performance::check_singularity(Cand.mod[[47]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[47]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[47]])) # Close!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[47]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[47]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[47]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[47]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=47, R2=R2[[1]]))



##### Model 48 #####
# -----------------

Cand.mod[[48]] <- glmmTMB::glmmTMB(formula = reg_edges~log2(trench_depth+1) + log2(distance+1) + repairs
                                   + (1|manager_id), data = redges, family = stats::binomial(link = "logit"),
                                   REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[48]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$trench_depth) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$distance) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$repairs) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[48]])
# performance::check_autocorrelation(Cand.mod[[48]])
# performance::check_collinearity(Cand.mod[[48]])
# performance::check_singularity(Cand.mod[[48]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[48]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[48]])) # Close!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[48]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[48]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[48]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[48]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=48, R2=R2[[1]]))



##### Model 49 #####
# -----------------

Cand.mod[[49]] <- glmmTMB::glmmTMB(formula = reg_edges~log2(trench_depth+1) + obstacles + repairs
                                   + (1|manager_id), data = redges, family = stats::binomial(link = "logit"),
                                   REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[49]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$trench_depth) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$obstacles) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$repairs) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[49]])
# performance::check_autocorrelation(Cand.mod[[49]])
# performance::check_collinearity(Cand.mod[[49]])
# performance::check_singularity(Cand.mod[[49]])
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[49]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[49]])) # Close!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[49]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[49]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[49]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[49]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=49, R2=R2[[1]]))



##### Model 50 #####
# -----------------

Cand.mod[[50]] <- glmmTMB::glmmTMB(formula = reg_edges~log2(trench_depth+1) + obstacles + geomem
                                   + (1|manager_id), data = redges, family = stats::binomial(link = "logit"),
                                   REML = FALSE)


# ### Model diagnostics:
# # Simulation-based scaled residuals (DHARMa method):
# simu.resid <- DHARMa::simulateResiduals(fittedModel = Cand.mod[[50]], n = 1000, plot = FALSE)
# plot(simu.resid) # Ok-ish!
# DHARMa::plotResiduals(simu.resid, form = redges$trench_depth) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$obstacles) # Ok!
# DHARMa::plotResiduals(simu.resid, form = redges$geomem) # Ok-ish!
# # Testing for overdispersion:
# DHARMa::testDispersion(simu.resid) # Ok!
# # Testing for spatial autocorrelation:
# DHARMa::testSpatialAutocorrelation(simulationOutput = simu.resid,
#                                    x = redges$longitude, y = redges$latitude, plot = TRUE) # Ok!
# # Using the 'performance' package:
# performance::check_distribution(Cand.mod[[50]])
# performance::check_autocorrelation(Cand.mod[[50]])
# performance::check_collinearity(Cand.mod[[50]])
# performance::check_singularity(Cand.mod[[50]]) # Singularity!
#
# ### Assessing goodness-of-fit:
# # Test of Pearson's Chi2 residuals:
# dat.resid <- sum(stats::resid(Cand.mod[[50]], type = "pearson")^2)
# 1 - stats::pchisq(dat.resid, stats::df.residual(Cand.mod[[50]])) # Close!
#
# ### Exploring the model parameters and test hypotheses:
# summary(Cand.mod[[50]])
# # Computing a pseudo-R2:
1 - (as.numeric(-2 * stats::logLik(Cand.mod[[50]]))/as.numeric(-2 * stats::logLik(
  update(Cand.mod[[50]], ~1)))) # McFadden's pseudo-R2
R2 <- performance::r2_tjur(Cand.mod[[50]]) # tjur's (pseudo-R2 for GLMMs)
R.ajust <- rbind(R.ajust, data.frame(Model=50, R2=R2[[1]]))





# ------------------------------------- #
##### Model selection and averaging #####
# ------------------------------------- #

Model <- (1:50)

Candidate <- c("null",
               "log2(distance + 1) + add_control", "log2(distance + 1) + obstacles",
               "log2(distance + 1) + repairs", "log2(distance + 1) + pb_fixation", "log2(distance + 1) + log2(stand_surface)",
               "log2(distance + 1) + geomem", "log2(distance + 1) + tarping_duration", "log2(distance + 1) + log2(trench_depth + 1)",
               "log2(distance + 1) * log2(trench_depth + 1)", "log2(trench_depth + 1) + add_control", "log2(trench_depth + 1) + obstacles",
               "log2(trench_depth + 1) + log2(stand_surface)", "log2(stand_surface) + add_control", "log2(stand_surface) + obstacles",
               "log2(stand_surface) * log2(distance + 1)", "log2(stand_surface) + tarping_duration", "log2(trench_depth + 1) + repairs",
               "log2(trench_depth + 1) + slope", "add_control + obstacles", "add_control + pb_fixation",
               "add_control + slope", "add_control + tarping_duration", "obstacles + pb_fixation",
               "obstacles + geomem", "obstacles + repairs", "slope + geomem", "slope + pb_fixation",
               "slope + tarping_duration", "log2(distance + 1) + log2(trench_depth + 1) + add_control",
               "log2(distance + 1) + log2(trench_depth + 1) + obstacles", "log2(distance + 1) + log2(trench_depth + 1) + geomem",
               "log2(distance + 1) + log2(trench_depth + 1) + log2(stand_surface)", "log2(distance + 1) + log2(trench_depth + 1) + slope",
               "log2(distance + 1) + log2(stand_surface) + add_control", "log2(distance + 1) + log2(stand_surface) + obstacles",
               "log2(distance + 1) + log2(stand_surface) + pb_fixation", "log2(distance + 1) + geomem + tarping_duration",
               "log2(stand_surface) + log2(trench_depth + 1) + obstacles", "log2(stand_surface) + log2(trench_depth + 1) + slope",
               "tarping_duration + log2(distance + 1) + add_control", "obstacles + repairs + add_control",
               "obstacles + add_control + geomem", "obstacles + add_control + tarping_duration",
               "obstacles + add_control + log2(trench_depth + 1)", "obstacles + add_control + log2(stand_surface)",
               "obstacles + add_control + slope", "log2(trench_depth + 1) + log2(distance + 1) + repairs",
               "log2(trench_depth + 1) + obstacles + repairs", "log2(trench_depth + 1) + obstacles + geomem")

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
AICc.model$Response <- "reg_edges"
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
readr::write_csv2(x = AICc.model, file = here::here("output", "tables", "Models_reg_edges.csv"))



##### Multimodel Inference (averaging) #####
# ------------------------------------------

Parameters <- c("Intercept", "add_control", "log2(distance + 1)", "geomem", "obstacles",
                "pb_fixation", "slope", "log2(stand_surface)", "tarping_duration",
                "log2(trench_depth + 1)", "repairs")
Var <- c("cond((Int))", "cond(add_control1)", "cond(log2(distance + 1))", "cond(geomem1)",
         "cond(obstacles)", "cond(pb_fixation1)", "cond(slope)", "cond(log2(stand_surface))",
         "cond(tarping_duration)", "cond(log2(trench_depth + 1))","cond(repairs1)")
Var.Imp <- c("cond((Int))", "cond(add_control)", "cond(log2(distance + 1))", "cond(geomem)",
             "cond(obstacles)", "cond(pb_fixation)", "cond(slope)", "cond(log2(stand_surface))",
             "cond(tarping_duration)", "cond(log2(trench_depth + 1))","cond(repairs)")

Para.model <- data.frame(Parameters, Var, Var.Imp)

### Select the top models:
# Removing the models with interactions before averaging:
AICc.2 <- AICc[-c(1,12),]
#top.models <- MuMIn::get.models(AICc.2, cumsum(weight) <= 0.95) # To take those with a cumulative sum of
# AICc weights <= 0.95
top.models <- MuMIn::get.models(AICc.2, cumsum(weight) <= 1) # To take them all (we chose this option because
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
Parameter.model$Response <- "reg_edges"
Parameter.model$'N model' <- length(top.models)
Parameter.model <- Parameter.model[,c("Response", "Parameters", "Imp.", "Estimate (±SE)", "(95% CI)",
                                      'N model')]
Parameter.model <- Parameter.model[!is.na(Parameter.model$Imp.), ]

### Table export:
readr::write_csv2(x = Parameter.model, file = here::here("output", "tables", "Parameters_reg_edges.csv"))
