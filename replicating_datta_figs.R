##Replicating the plots from Datta and Blanchard using MizerSeasonal


library(mizer)
library(mizerSeasonal)
library(mizerExperimental)
library(dplyr)
library(ggplot2)



##we start with the parameters, the same as used by Datta and Blanchard
#
##take datta_params, encode seasonality, produce the right plots
pars <- datta_params
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
s_sim_1 <- project(s_ps, t_max = 1, dt = 0.001, t_save = 0.01)
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
datta_params_seasonal <- setRateFunction(s_ps, "RDD", "seasonalBevertonHoltRDD")


#Project the models
#
##500 year projected seasonal model
datta_500yr <- project(datta_params_seasonal, t_max = 500, dt = 0.01, t_save = 0.01)
##500 year projected nonseasonal model
datta_nonseasonal_500yr <- project(datta_params, t_max = 500, dt = 0.1, t_save = 0.1)









#Plotting

##now produce the plots
plotSpectra(datta_nonseasonal_500yr, power =1, resource = FALSE,
            time_range = 499:500,total = TRUE, wlim = c(1e-03,1e05), 
            ylim= c(1e-05, 1e15))

plotSpectra(datta_nonseasonal_500yr, power =0, resource = FALSE,
              time_range = 499:500, total = TRUE, wlim = c(1e-03,1e05), 
              ylim= c(1e-05, 1e15))

plotSpectra(datta_500yr, power =0, resource = FALSE,
            time_range = 499:500, total = TRUE, wlim = c(1e-03,1e05), 
            ylim= c(1e-05, 1e15))






##try to isolate the two community spectras
#
plotSpectra2(datta_500yr, datta_nonseasonal_500yr, title1= "seasonal", title2 = "nonseasonal",
             power = 1, resource = TRUE, 
             time_range = 499:500, total = TRUE, species = c("") )




  
  
  
  
  
##plotting the individual spectras for different fish species, their seasonal
#with their nonseasonal

#sprat
plotSpectra2(datta_500yr, datta_nonseasonal_500yr, 
             name1 = "Seasonal", name2 = "Nonseasonal", power = 0, 
             species = c("Sprat") , resource = FALSE,
             time_span = 499:500)

#saithe
plotSpectra2(datta_500yr, datta_nonseasonal_500yr, 
             name1 = "Seasonal", name2 = "Nonseasonal", power = 0, 
             species = c("Saithe") , resource = FALSE,
             time_span = 499:500)

#cod
plotSpectra2(datta_500yr, datta_nonseasonal_500yr, 
             name1 = "Seasonal", name2 = "Nonseasonal", power = 0, 
             species = c("Cod") , resource = FALSE,
             time_span = 499:500)

