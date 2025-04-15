#figures for paper

library(scales)
non_seas <- newSingleSpeciesParams(species_name = "Species", no_w = log10(100/0.001) * 20 + 1)
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


##1
simulation_1 <-  project(params_20, t_max = 1, dt = 0.01, t_save = 0.01)

nonseas_sim <- project(non_seas, t_max = 1, dt = 0.01, t_save = 0.01)

plotSpectra2(nonseas_sim, simulation_1, name1 = 'Nonseasonal', name2 = 'Seasonal',
             wl = c(1e-03, 1e+02),power = 2,
 baseSize = 20) +  # Increase base font size
  theme(
    axis.text = element_text(size = 20),    # Increase axis number size
    axis.title = element_text(size = 22),   # Increase axis title size
    legend.text = element_text(size = 18),  # Increase legend text
    legend.title = element_text(size = 20), 
    panel.grid.major = element_line(size = 1.2, linetype = "solid"), # Thicker major grid
    panel.grid.minor = element_line(size = 0.8, linetype = "dashed") # Thicker minor grid
  ) 






##2
plotSpectra2(ps, ps_datta)
##

library(ChemoSpec)

plotSpectra2(ps, ps_datta, name1 = "mizerSeasonal", name2 = "Datta", baseSize = 20) +  # Increase base font size
  theme(
    axis.text = element_text(size = 20),    # Increase axis number size
    axis.title = element_text(size = 22),   # Increase axis title size
    legend.text = element_text(size = 18),  # Increase legend text
    legend.title = element_text(size = 20), 
    panel.grid.major = element_line(size = 1.2, linetype = "solid"), # Thicker major grid
    panel.grid.minor = element_line(size = 0.8, linetype = "dashed") # Thicker minor grid
  )









##3
plotSpectra2(simulation_1, simulation_13, name1= '20', name2 = '10000',
             power = 2, wlim = c(10,100),
             baseSize = 20) +  # Increase base font size
  theme(
    axis.text = element_text(size = 20),    # Increase axis number size
    axis.title = element_text(size = 22),   # Increase axis title size
    legend.text = element_text(size = 18),  # Increase legend text
    legend.title = element_text(size = 20), 
    
    panel.grid.major = element_line(size = 1.2, linetype = "solid"), # Thicker major grid
    panel.grid.minor = element_line(size = 0.8, linetype = "dashed") # Thicker minor grid
  )





##4



plotSpectra2(simulation_1, simulation_2, power = 2, wlim = c(0.25,100),ylim = c(1e-06,1e-03))
plotSpectra2(simulation_4, simulation_5, power = 2, wlim = c(0.25,100),ylim = c(1e-06,1e-03))
plotSpectra2(simulation_7, simulation_8, power = 2, wlim = c(0.25,100),ylim = c(1e-06,1e-03))
plotSpectra2(simulation_10, simulation_11, power = 2, wlim = c(0.25,100), ylim = c(1e-06,1e-03))


# Load required libraries
library(ChemoSpec)
library(ggplot2)
library(patchwork)

# Create individual plots
p2 <- plotSpectra2(simulation_1, simulation_2, name1 = "dt = 0.01", name2 = "dt = 0.001", power = 2, wlim = c(0.25,100),ylim = c(1e-06,1e-03)) + 
  ggtitle(expression(Omega == 20))  + theme(plot.title = element_text(size = 18),
    legend.position = "none", axis.text = element_text(size = 18),    # Increase axis number size
                             axis.title = element_text(size = 18),   # Increase axis title size
                             legend.text = element_text(size = 18),  # Increase legend text
                             legend.title = element_text(size = 20), 
                             panel.grid.major = element_line(size = 1.2, linetype = "solid"), # Thicker major grid
                             panel.grid.minor = element_line(size = 0.8, linetype = "dashed") # Thicker minor grid
  )
p1 <- plotSpectra2(simulation_4, simulation_5,power = 2, wlim = c(0.25,100),ylim = c(1e-06,1e-03)) + 
  ggtitle(expression(Omega == 0.1)) + theme(plot.title = element_text(size = 18),legend.position = "none", axis.text = element_text(size = 18),    # Increase axis number size
                             axis.title = element_text(size = 18),   # Increase axis title size
                             legend.text = element_text(size = 18),  # Increase legend text
                             legend.title = element_text(size = 20), 
                             panel.grid.major = element_line(size = 1.2, linetype = "solid"), # Thicker major grid
                             panel.grid.minor = element_line(size = 0.8, linetype = "dashed") # Thicker minor grid
  )
p3 <- plotSpectra2(simulation_7, simulation_8, power = 2, wlim = c(0.25,100),ylim = c(1e-06,1e-03)) + 
  ggtitle(expression(Omega == 100)) + theme(plot.title = element_text(size = 18),legend.position = "none", axis.text = element_text(size = 18),    # Increase axis number size
                             axis.title = element_text(size = 18),   # Increase axis title size
                             legend.text = element_text(size = 18),  # Increase legend text
                             legend.title = element_text(size = 20), 
                             
                             panel.grid.major = element_line(size = 1.2, linetype = "solid"), # Thicker major grid
                             panel.grid.minor = element_line(size = 0.8, linetype = "dashed") # Thicker minor grid
  )
p4 <- plotSpectra2(simulation_10, simulation_11, power = 2, wlim = c(0.25,100),ylim = c(1e-06,1e-03)) + 
  ggtitle(expression(Omega == 1000)) + theme(plot.title = element_text(size = 18),legend.position = "none", axis.text = element_text(size = 18),    # Increase axis number size
                              axis.title = element_text(size = 18),   # Increase axis title size
                              legend.text = element_text(size = 18),  # Increase legend text
                              legend.title = element_text(size = 20), 
                              
                              panel.grid.major = element_line(size = 1.2, linetype = "solid"), # Thicker major grid
                              panel.grid.minor = element_line(size = 0.8, linetype = "dashed") # Thicker minor grid
  )

forthelegend <- plotSpectra2(simulation_1, simulation_2, name1 = "dt = 0.01", name2 = "dt = 0.001",power = 2, wlim = c(0.25,100),ylim = c(1e-06,1e-03)) + 
  ggtitle(expression(Omega == 20))  + theme(plot.title = element_text(size = 18),legend.position = "bottom", axis.text = element_text(size = 18),    # Increase axis number size
                             axis.title = element_text(size = 18),   # Increase axis title size
                             legend.text = element_text(size = 22),  # Increase legend text
                             legend.title = element_text(size = 22), 
                             
                             panel.grid.major = element_line(size = 1.2, linetype = "solid"), # Thicker major grid
                             panel.grid.minor = element_line(size = 0.8, linetype = "dashed") # Thicker minor grid
  )

get_legend <- function(plot) {
  g <- ggplotGrob(plot)
  legend <- g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]
  return(legend)
}

legend <- get_legend(forthelegend) 

# Arrange plots in a 2x2 grid
plot_grid <- (p1 + p2) / (p3 + p4) / (legend)

# Display the final figure
print(plot_grid)










###5
library(ChemoSpec)
library(ggplot2)
library(patchwork)

library(mizer)
library(mizerExperimental)
library(mizerSeasonal)
library(ggplot2)

simulation_2_wspeciesname <-  project(params_20, t_max = 1, dt = 0.001, t_save = 0.001, label = "Species")

#(dt = 0.001)
#simulation_2 (20)
#simulation_5 (0.1)
#simulation_8 (100)
#simulation_11 (1000)
#simulation_14 (10000)



plotSpectra(simulation_2_wspeciesname)










# Load required libraries
library(ChemoSpec)
library(ggplot2)
library(patchwork)

# Plot with paired simulations
p1 <- plotSpectra2(simulation_5, simulation_5, power = 2, wlim = c(20, 100)) +
  ggtitle(expression(Omega == 0.1))+
  theme(
    legend.position = "none",
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18),
    panel.grid.major = element_line(size = 1),
    panel.grid.minor = element_line(size = 0.5)
  )

p2 <- plotSpectra2(simulation_2_wspeciesname, simulation_2_wspeciesname, power = 2, wlim = c(20, 100)) +
  ggtitle(expression(Omega == 20)) + 
  guides(linetype = "none", colour = guide_legend(title = NULL)) +
  theme(legend.position = c(0,0),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 18),
        panel.grid.major = element_line(size = 1),
        panel.grid.minor = element_line(size = 0.5))

p3 <- plotSpectra2(simulation_8, simulation_8, power = 2, wlim = c(20, 100)) +
  ggtitle(expression(Omega == 100)) +
  theme(legend.position = "none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 18),
        panel.grid.major = element_line(size = 1),
        panel.grid.minor = element_line(size = 0.5))

p4 <- plotSpectra2(simulation_10, simulation_10, power = 2, wlim = c(20, 100)) +
  ggtitle(expression(Omega == 1000))+
  theme(legend.position = "none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 18),
        panel.grid.major = element_line(size = 1),
        panel.grid.minor = element_line(size = 0.5))

# This one includes the legend!
p5 <- plotSpectra2(simulation_14, simulation_14, power = 2, wlim = c(20, 100)) +
  ggtitle(expression(Omega == 10000)) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    panel.grid.major = element_line(size = 1),
    panel.grid.minor = element_line(size = 0.5)
  )

# Combine all plots with p2 at the bottom to show the legend
final_plot <- (p1 + p3) / (p2 + p4) / p5 +
  plot_layout(guides = "collect") 


# Display
print(final_plot)
















##########. 
library(ChemoSpec)
library(ggplot2)
library(patchwork)

library(mizer)
library(mizerExperimental)
library(mizerSeasonal)
library(ggplot2)

s1 <- project(params_20, t_max = 1, dt = 0.1, t_save = 0.1)
s2 <- project(params_20, t_max = 1, dt = 0.01, t_save = 0.01)
s3 <- project(params_20, t_max = 1, dt = 0.001, t_save = 0.001)
s0 <- project(params_20, t_max = 1, dt = 1, t_save = 0.1)
nonseas_sim <- project(non_seas, t_max = 1, dt = 0.01, t_save = 0.01)
nonseas_sim2 <- project(non_seas, t_max = 1, dt = 0.1, t_save = 0.1)
nonseas_sim3 <- project(non_seas, t_max = 1, dt = 0.01, t_save = 0.01)
nonseas_sim4 <- project(non_seas, t_max = 1, dt = 0.001, t_save = 0.001)

ps0 <- plotSpectra2(s0, nonseas_sim, power = 2, wlim = c(1e-3, 100)) +
  ggtitle("dt = 1")+
  theme(
    legend.position = "none",
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18),
    panel.grid.major = element_line(size = 1),
    panel.grid.minor = element_line(size = 0.5)
  )

ps1 <- plotSpectra2(s1, nonseas_sim2, power = 2, wlim = c(1e-3, 100)) +
  ggtitle("dt = 0.1")+
  theme(
    legend.position = "none",
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18),
    panel.grid.major = element_line(size = 1),
    panel.grid.minor = element_line(size = 0.5)
  )
ps2 <- plotSpectra2(s2, nonseas_sim3, name1 = "Seasonal", name2 = "Nonseasonal", power = 2, wlim = c(1e-3, 100)) +
  ggtitle("dt = 0.01")+
  theme(
    legend.position = "bottom",
    legend.justification = c(0, 0),  # aligns legend to bottom-left
    legend.box.just = "left",
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    panel.grid.major = element_line(size = 1),
    panel.grid.minor = element_line(size = 0.5)
  )

ps3 <- plotSpectra2(s3, nonseas_sim4, name1 = "Seasonal", name2 = "Nonseasonal", power = 2, wlim = c(1e-3, 100)) +
  ggtitle("dt = 0.001")+

  theme(
    legend.position = "none",
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18),
    panel.grid.major = element_line(size = 1),
    panel.grid.minor = element_line(size = 0.5)
  )


final_plot <- (ps0 + ps1 + ps2 +ps3)
  plot_layout(guides = "collect") 


# Display
print(final_plot)




