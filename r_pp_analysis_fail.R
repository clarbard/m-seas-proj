##final attempt at r_pp 


library(mizer)
library(mizerSeasonal)
library(mizerExperimental)
library(dplyr)
library(ggplot2)

##creating single species model with defaults, but with seasonality turned on
#no_w = *20 : params_20
####create a single species model, start with presets, then turn on seasonality
non_seas <- newSingleSpeciesParams(species_name = "defaults")
pars <- non_seas
##now encode seasonality, using sampele values for mu and kappa
species_params(pars)$rdd_vonMises_r0 <- getRDD(pars)
species_params(pars)$rdd_vonMises_kappa <- 1
species_params(pars)$rdd_vonMises_mu <- 0.5
species_params(pars)$vonMises_r0 <- 100
species_params(pars)$vonMises_kappa <- species_params(pars)$rdd_vonMises_kappa
species_params(pars)$vonMises_mu <- species_params(pars)$rdd_vonMises_mu
pars <- setSeasonalReproduction(pars, release_func = "seasonalVonMisesRelease",
                                RDD = "seasonalVonMisesRDD")
s_sim <- project(pars, t_max = 50, dt = 0.001)
s_ps <- setInitialValues(pars, s_sim)
s_sim_1 <- project(s_ps, t_max = 1, dt = 0.001, t_save = 0.001)
rdi <- getTimeseries(s_sim_1, func = getRDI)
rdd <- getTimeseries(s_sim_1, func = getRDD)
rdi_df <- melt(rdi)
rdi_df$type <- "RDI"
rdd_df <- melt(rdd)
rdd_df$type <- "observed"
df <- rbind(rdi_df, rdd_df)

ggplot(df) + 
  geom_line(aes(x = time, y = value, colour = type)) +
  facet_wrap(vars(sp)) +
  scale_y_log10()
other_params(s_ps)$r_max <- rdi * rdd / (rdi - rdd)
params_default <- setRateFunction(s_ps, "RDD", "seasonalBevertonHoltRDD")



#now do a new one messing with the resource

non_seas <- newSingleSpeciesParams(species_name = "mod1")

pars <- setResource(non_seas, resource_dynamics = "resource_logistic")



##now encode seasonality, using sampele values for mu and kappa
species_params(pars)$rdd_vonMises_r0 <- getRDD(pars)
species_params(pars)$rdd_vonMises_kappa <- 1
species_params(pars)$rdd_vonMises_mu <- 0.5
species_params(pars)$vonMises_r0 <- 100
species_params(pars)$vonMises_kappa <- species_params(pars)$rdd_vonMises_kappa
species_params(pars)$vonMises_mu <- species_params(pars)$rdd_vonMises_mu
pars <- setSeasonalReproduction(pars, release_func = "seasonalVonMisesRelease",
                                RDD = "seasonalVonMisesRDD")
s_sim <- project(pars, t_max = 50, dt = 0.001)
s_ps <- setInitialValues(pars, s_sim)
s_sim_1 <- project(s_ps, t_max = 1, dt = 0.001, t_save = 0.001)
rdi <- getTimeseries(s_sim_1, func = getRDI)
rdd <- getTimeseries(s_sim_1, func = getRDD)
rdi_df <- melt(rdi)
rdi_df$type <- "RDI"
rdd_df <- melt(rdd)
rdd_df$type <- "observed"
df <- rbind(rdi_df, rdd_df)

ggplot(df) + 
  geom_line(aes(x = time, y = value, colour = type)) +
  facet_wrap(vars(sp)) +
  scale_y_log10()
other_params(s_ps)$r_max <- rdi * rdd / (rdi - rdd)
params_mod1 <- setRateFunction(s_ps, "RDD", "seasonalBevertonHoltRDD")

sim_mod1 <- project(params_mod1, t_max = 1, dt = 0.01)
plotSpectra(sim_mod1, power = 2) # setting power = 2 means can see not constant resource
plotResourceLevel(sim_mod1) #so logistic is happening

#changing resource_rate
params_mod1_2 <- setResource(params_mod1, resource_rate = 100)
sim_mod1_2 <- project(params_mod1_2, t_max = 1, dt = 0.01)
plotSpectra(sim_mod1_2, power = 2)
plotResourceLevel(sim_mod1_2)

#changing further
params_mod1_3 <- setResource(params_mod1, n = 15)
sim_mod1_3 <- project(params_mod1_3, t_max = 1, dt = 0.01)
plotSpectra(sim_mod1_3, power = 2)
plotResourceLevel(sim_mod1_3)









#changing further
params_mod1_4 <- params_mod1
params_mod1_5 <- params_mod1

resource_rate(params_mod1_4) <- 10000
resource_rate(params_mod1_5) <- 0.5



plotResourceLevel(params_mod1)
plotResourceLevel(params_mod1_4)
plotResourceLevel(params_mod1_5)



plotSpectra2(params_mod1, params_mod1_5, power =2) #looks same


plotSpectra2(params_mod1, params_mod1_5, power =2, wlim = c(30,100)) #still looks same




plotSpectra2(params_mod1, params_mod1_5, power =2, wlim = c(0,100), 
             time_range = c(0:1))






params_mod2_1 <- params_mod1

vec <- abs((1:282))

resource_rate(params_mod2_1) <- vec
plotResourceLevel(params_mod2_1)


plotSpectra(params_mod2_1, power = 2)


sim_mod1_4 <- project(params_mod1_4, t_max = 1, dt = 0.01)
plotSpectra(sim_mod1_4, power = 2)
plotResourceLevel(sim_mod1_4)
##doesn't even look different??? why?????





##try semichemostatic

non_seas <- newSingleSpeciesParams(species_name = "defaults")
pars <- non_seas
##now encode seasonality, using sampele values for mu and kappa
species_params(pars)$rdd_vonMises_r0 <- getRDD(pars)
species_params(pars)$rdd_vonMises_kappa <- 1
species_params(pars)$rdd_vonMises_mu <- 0.5
species_params(pars)$vonMises_r0 <- 100
species_params(pars)$vonMises_kappa <- species_params(pars)$rdd_vonMises_kappa
species_params(pars)$vonMises_mu <- species_params(pars)$rdd_vonMises_mu
pars <- setSeasonalReproduction(pars, release_func = "seasonalVonMisesRelease",
                                RDD = "seasonalVonMisesRDD")
s_sim <- project(pars, t_max = 50, dt = 0.001)
s_ps <- setInitialValues(pars, s_sim)
s_sim_1 <- project(s_ps, t_max = 1, dt = 0.001, t_save = 0.001)
rdi <- getTimeseries(s_sim_1, func = getRDI)
rdd <- getTimeseries(s_sim_1, func = getRDD)
rdi_df <- melt(rdi)
rdi_df$type <- "RDI"
rdd_df <- melt(rdd)
rdd_df$type <- "observed"
df <- rbind(rdi_df, rdd_df)

ggplot(df) + 
  geom_line(aes(x = time, y = value, colour = type)) +
  facet_wrap(vars(sp)) +
  scale_y_log10()
other_params(s_ps)$r_max <- rdi * rdd / (rdi - rdd)
params_default <- setRateFunction(s_ps, "RDD", "seasonalBevertonHoltRDD")



#now do a new one messing with the resource

non_seas <- newSingleSpeciesParams(species_name = "mod2")

pars <- setResource(non_seas, resource_dynamics = "resource_semichemostat")



##now encode seasonality, using sampele values for mu and kappa
species_params(pars)$rdd_vonMises_r0 <- getRDD(pars)
species_params(pars)$rdd_vonMises_kappa <- 1
species_params(pars)$rdd_vonMises_mu <- 0.5
species_params(pars)$vonMises_r0 <- 100
species_params(pars)$vonMises_kappa <- species_params(pars)$rdd_vonMises_kappa
species_params(pars)$vonMises_mu <- species_params(pars)$rdd_vonMises_mu
pars <- setSeasonalReproduction(pars, release_func = "seasonalVonMisesRelease",
                                RDD = "seasonalVonMisesRDD")
s_sim <- project(pars, t_max = 50, dt = 0.001)
s_ps <- setInitialValues(pars, s_sim)
s_sim_1 <- project(s_ps, t_max = 1, dt = 0.001, t_save = 0.001)
rdi <- getTimeseries(s_sim_1, func = getRDI)
rdd <- getTimeseries(s_sim_1, func = getRDD)
rdi_df <- melt(rdi)
rdi_df$type <- "RDI"
rdd_df <- melt(rdd)
rdd_df$type <- "observed"
df <- rbind(rdi_df, rdd_df)

ggplot(df) + 
  geom_line(aes(x = time, y = value, colour = type)) +
  facet_wrap(vars(sp)) +
  scale_y_log10()
other_params(s_ps)$r_max <- rdi * rdd / (rdi - rdd)
params_mod2 <- setRateFunction(s_ps, "RDD", "seasonalBevertonHoltRDD")


plotResourceLevel(params_mod2) #still looks logistic
params_mod2 <- setResource(params_mod2, resource_dynamics = "resource_semichemostat")


resource_dynamics(params_mod2) <- "resource_semichemostat"
plotResourceLevel(params_mod2)

##still don't see anything


#these are all just straight lines
plotSpectra(params_mod2, wlim = c(0,0.001))
plotSpectra(params_mod1_5, wlim = c(0,0.001))
plotSpectra(params_mod1, wlim = c(0,0.001))

