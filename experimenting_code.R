##new script 27th jan 2025

##This is for investigating changing no_w and dt to see if that changes where seasonal variations stop

library(mizer)
library(mizerSeasonal)
library(mizerExperimental)
library(dplyr)
library(ggplot2)


##fixing my errors



###want to try: looking at no_w coef *20, *0.1, *100
##and try looking at dt = 0.01, 0.1, 0.001





###########################################################

##creating single species model with defaults, but with seasonality turned on

#no_w = *20 : params_20
####create a single species model, start with presets, then turn on seasonality
non_seas <- newSingleSpeciesParams(species_name = "*20", no_w = log10(100/0.001) * 20 + 1)
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
params_20 <- setRateFunction(s_ps, "RDD", "seasonalBevertonHoltRDD")





##no_w *0.1: params_01
p_start <- newSingleSpeciesParams(species_name = "*0.1", no_w = log10(100/0.001) * 0.1 + 1)
pars <- p_start

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
params_01 <- setRateFunction(s_ps, "RDD", "seasonalBevertonHoltRDD")



##no_w *100: params_100
p_start <- newSingleSpeciesParams(species_name = "*100", no_w = log10(100/0.001) * 100 + 1)
pars <- p_start

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
params_100 <- setRateFunction(s_ps, "RDD", "seasonalBevertonHoltRDD")












##so we have params_01, params_20, params_100


###########################################################

##simulations (changing dt)


##simulation_1: *20 and dt = 0.01
simulation_1 <-  project(params_20, t_max = 1, dt = 0.01, t_save = 0.01)



##simulation_2: *20 and dt = 0.001
simulation_2 <-  project(params_20, t_max = 1, dt = 0.001, t_save = 0.01)

##simulation_3: *20 and dt = 0.1
simulation_3 <-  project(params_20, t_max = 1, dt = 0.1, t_save = 0.01)




##simulation_4 : *0.1 and dt = 0.01
simulation_4 <-  project(params_01, t_max = 1, dt = 0.01, t_save = 0.01)

##simulation_5: *0.1 and dt = 0.001

simulation_5 <-  project(params_01, t_max = 1, dt = 0.001, t_save = 0.01)

#simulation_6: *0.1 and dt = 0.1
simulation_6 <-  project(params_01, t_max = 1, dt = 0.1, t_save = 0.01)




##simulation_7: *100 and dt = 0.01
simulation_7 <-  project(params_100, t_max = 1, dt = 0.01, t_save = 0.01)


##simulation_8: *100 and dt = 0.001

simulation_8 <-  project(params_100, t_max = 1, dt = 0.001, t_save = 0.01)

#simulation_9: *100 and dt = 0.1
simulation_9 <-  project(params_100, t_max = 1, dt = 0.1, t_save = 0.01)




##PLOT


##make plots, and zoom into the interesting size range


plot_1<- animateSpectra(simulation_1, power = 2, wlim = c(10,60))
plot_2<- animateSpectra(simulation_2, power = 2, wlim = c(10,60))
plot_3<- animateSpectra(simulation_3, power = 2, wlim = c(10,60))
plot_4<- animateSpectra(simulation_4, power = 2, wlim = c(10,60))
plot_5<- animateSpectra(simulation_5, power = 2, wlim = c(10,60))
plot_6<- animateSpectra(simulation_6, power = 2, wlim = c(10,60))
plot_7<- animateSpectra(simulation_7, power = 2, wlim = c(10,60))
plot_8<- animateSpectra(simulation_8, power = 2, wlim = c(10,60))
plot_9<- animateSpectra(simulation_9, power = 2, wlim = c(10,60))





###now plot them with titles


plot_1 #+ ggtitle("no_w x20, dt 0.01")

plot_2 #+ ggtitle("no_w x20, dt 0.001")

plot_3 #+ ggtitle("no_w x20, dt 0.1")

plot_4 #+ ggtitle("no_w x0.1, dt 0.01")

plot_5 #+ ggtitle("no_w x0.1, dt 0.001")

plot_6 #+ ggtitle("no_w x0.1, dt 0.1")

plot_7 #+ ggtitle("no_w x100, dt 0.01")

plot_8 #+ ggtitle("no_w x100, dt 0.001")

plot_9 #+ ggtitle("no_w x100, dt 0.1")






##looking at the spectras #replace number with desired simulation
plotSpectra(simulation_1)
plotSpectra(simulation_3)






#no_w *20 and *0.1
plotSpectra2(simulation_1, simulation_4)

#no_w *20 and *100
plotSpectra2(simulation_1, simulation_7)

##no_w *0.1 and *100
plotSpectra2(simulation_4, simulation_7)


##dt = 0.01 and dt = 0.1 for no_w = *20
plotSpectra2(simulation_1, simulation_3)




###save the animmated spectras as html on desktop
#library(htmlwidgets)
# Save the plot as an HTML file
#htmlwidgets::saveWidget(plot_1, file = "~/Desktop/plot_1.html")



##plot relative difference between the no_w changes, relative to their average
#plotSpectraRelative(simulation_1, simulation_3)






plotSpectra2(simulation_1, simulation_9, power = 2, wlim = c(10,100))





















###making another simulation, trying no_w: *1,000 with dt = 0.001



##no_w *1000: params_1000
p_start <- newSingleSpeciesParams(species_name = "*1000", no_w = log10(100/0.001) * 1000 + 1)
pars <- p_start

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
params_1000 <- setRateFunction(s_ps, "RDD", "seasonalBevertonHoltRDD")



simulation_10 <-  project(params_1000, t_max = 1, dt = 0.01, t_save = 0.01)
simulation_11 <-  project(params_1000, t_max = 1, dt = 0.001, t_save = 0.01)
simulation_12 <-  project(params_1000, t_max = 1, dt = 0.1, t_save = 0.01)



plotSpectra(simulation_10)
plotSpectra(simulation_11)
plotSpectra(simulation_12)






plotSpectra2(simulation_1, simulation_11, power = 2)
plotSpectra2(simulation_1, simulation_11, power = 2, wlim = c(10,100))
plotSpectra2(simulation_1, simulation_10, power = 2, wlim = c(10,100))




plotSpectra2(simulation_1, simulation_2, power = 2, wlim = c(10,100))
plotSpectra2(simulation_4, simulation_5, power = 2, wlim = c(10,100))
plotSpectra2(simulation_7, simulation_8, power = 2, wlim = c(10,100))
plotSpectra2(simulation_10, simulation_11, power = 2, wlim = c(10,100))
             
             
             
             
             




             ###ONE MORE
##no_w *10000: params_10000
p_start <- newSingleSpeciesParams(species_name = "*10000", no_w = log10(100/0.001) * 10000 + 1)
pars <- p_start

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
params_10000 <- setRateFunction(s_ps, "RDD", "seasonalBevertonHoltRDD")









simulation_13 <-  project(params_10000, t_max = 1, dt = 0.01, t_save = 0.01)
simulation_14 <-  project(params_10000, t_max = 1, dt = 0.001, t_save = 0.01)


plotSpectra(simulation_13)
plotSpectra(simulation_14)


plotSpectra2(simulation_13, simulation_14, power = 2, wlim = c(0,60)) +
  scale_x_continuous()

plotSpectra2(simulation_11, simulation_14, power = 2, wlim = c(10,100))



plotSpectra2(simulation_1, simulation_14, power = 2, wlim = c(10,100))




















#####try taking average over a year, see if matches that of nonseasonal and others with less weight bins
nonseasonal <- project(non_seas, t_max = 1, dt = 0.01, t_save = 0.01)

plotSpectra(nonseasonal, time_range = 0:1)




##try over longer timespan

nonseasonal_long <- project(non_seas, t_max = 10, dt = 0.01, t_save = 0.01)
seasonal_default_long <- project(params_20, t_max = 10, dt = 0.01, t_save = 0.01)

#over last time step
plotSpectra2(nonseasonal_long, seasonal_default_long)
#averaged over last year
plotSpectra2(nonseasonal_long, seasonal_default_long, time_range = 9:10)
#averaged over full 10 years
plotSpectra2(nonseasonal_long, seasonal_default_long, time_range = 0:10)

#over last year and over 10 years looks same
#they are vertically shifted 

#try projecting longer
seasonal_default_longer <- project(params_20, t_max = 100, dt = 0.01, t_save = 0.01)
nonseasonal_longer <- project(non_seas, t_max = 100, dt = 0.01, t_save = 0.01)

#aveaged over last time step
plotSpectra2(nonseasonal_longer, seasonal_default_longer)
#averaged over the last year
plotSpectra2(nonseasonal_longer, seasonal_default_longer, time_range = 99:100)
#averaged over full 100 years
plotSpectra2(nonseasonal_longer, seasonal_default_longer, time_range = 0:100)



#try one last projecting even longer
seasonal_default_longest <- project(params_20, t_max = 1000, dt = 0.01, t_save = 0.01)
nonseasonal_longest <- project(non_seas, t_max = 1000, dt = 0.01, t_save = 0.01)


plotSpectra2(nonseasonal_longest, seasonal_default_longest)

plotSpectra2(nonseasonal_longest, seasonal_default_longest, time_range = 999:1000)

plotSpectra2(nonseasonal_longest, seasonal_default_longest, time_range = 0:1000)























###THE FOLLOWING IS FOR THE NUMERICAL APPROACH USING MATLAB, GETTING THE RIGHT PARAM VALUES


##FIND THE MORTALITY RATE COEFFICIENT AND EXPONENT  FOR PARAMS_20 THE DEFAULT
mort_rate_frame <- melt(getMort(params_20))
plot_mort <- ggplot(mort_rate_frame) +
  geom_line(aes(x = w_prey, y = value)) +
  scale_x_log10(name = "Weight [g]") +
  scale_y_log10(name = "Mortality rate [1/year]")
plot_mort

mm <- lm(log(mort_rate_frame$value) ~ log(mort_rate_frame$w_prey))
mm
m0 <- exp(coef(mm)[[1]])
m0



###GET THE GROWTH RATE COEFFICIENT AND EXPONENT FOR PARAMS_20 THE DEFAULT
growth_rate <- getEGrowth(params_20)
growth_rate_frame <- melt(growth_rate)
plot_growtht <- ggplot(growth_rate_frame) +
  geom_line(aes(x = w, y = value)) +
  scale_x_log10(name = "Weight [g]") +
  scale_y_log10(name = "Growth rate [g/year]")
plot_growtht
g_small_fish <- filter(growth_rate_frame, w <= 10)
lm(log(g_small_fish$value) ~ log(g_small_fish$w))








##FIND THE MORTALITY RATE COEFFICIENT AND EXPONENT  FOR PARAMS_20 THE DEFAULT
mort_rate_frame <- melt(getMort(params_1000))
plot_mort <- ggplot(mort_rate_frame) +
  geom_line(aes(x = w_prey, y = value)) +
  scale_x_log10(name = "Weight [g]") +
  scale_y_log10(name = "Mortality rate [1/year]")
plot_mort

mm <- lm(log(mort_rate_frame$value) ~ log(mort_rate_frame$w_prey))
mm
m0 <- exp(coef(mm)[[1]])
m0



###GET THE GROWTH RATE COEFFICIENT AND EXPONENT FOR PARAMS_20 THE DEFAULT
growth_rate <- getEGrowth(params_1000)
growth_rate_frame <- melt(growth_rate)
plot_growtht <- ggplot(growth_rate_frame) +
  geom_line(aes(x = w, y = value)) +
  scale_x_log10(name = "Weight [g]") +
  scale_y_log10(name = "Growth rate [g/year]")
plot_growtht
g_small_fish <- filter(growth_rate_frame, w <= 10)
lm(log(g_small_fish$value) ~ log(g_small_fish$w))











###now create the plots for looking at different averaging times for the discrepancy of start and end

plotSpectra2(nonseasonal_longest, seasonal_default_longest, power = 2, time_range= 0:1)
plotSpectra2(nonseasonal_longest, seasonal_default_longest, power = 2, time_range= 0:10)
plotSpectra2(nonseasonal_longest, seasonal_default_longest, power = 2, time_range= 0:100)
plotSpectra2(nonseasonal_longest, seasonal_default_longest, power = 2, time_range= 0:1000)
plotSpectra2(nonseasonal_longest, seasonal_default_longest, power = 2, time_range= 999:1000)













########
plotSpectra2(nonseasonal_longest, seasonal_default_longest, power = 2, time_range = 0:1)



initialN(nonseasonal_longest)
initialN(seasonal_default_longest)



#doesn't work
ns_n <- get_initial_n(non_seas)
s_n <- get_initial_n(params_20)







###Come back - find where the initial values come from in every sim; 





