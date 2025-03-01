##combining the seasonality from mizerSeasonal Delius 2024 with the growth curve
#calculations and plotting code from Datta et al. 2019
library(mizer)
library(mizerSeasonal)
library(mizerExperimental)
library(dplyr)
library(ggplot2)





#start with datta_params from mizerSeasonal, already taken to match baseModel$param

#######################################

##take datta_params, encode seasonality
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
other_params(s_ps)$r_max <- rdi * rdd / (rdi - rdd)
datta_params_seasonal <- setRateFunction(s_ps, "RDD", "seasonalBevertonHoltRDD")




################################

#work with two params objects, same parameters as baseModel (Datta et al.), one nonseasonal and the 
#other with seasonality encoded by mizerSeasonal

look_nonseas <- datta_params #nonseasonal params object
look_seas <- datta_params_seasonal #seasonal params object


################################
#declare some parameters to  use, like Datta et al.

nspp <- length(look_seas@species_params$species) #number of species
sp <- look_seas@species_params$species #species
ngrid <- length(look_seas@w) #width of weight bin grid (how many bins)
n_pp_grid <- length(look_seas@w_full) #how many weight bins in grid including resource
tmax <- 500 #they project 500 years
dt = 1/52 # time step (years)
isave <- 1 #would change, but this is what they have it equal to #saving slot

#set up grid like Datta
w <- look_seas@w
dw <- look_seas@dw
dw[ngrid] <- dw[ngrid-1] # Set final dw as same as one before
wFull <- look_seas@w_full
dwFull <- look_seas@dw_full
ngridfull<-length(wFull)
dwFull[ngridfull] <- dwFull[ngridfull-1]
nsave   <- floor(tmax/(dt*isave)) #num of time slots to save


##For reference
#param <- look_seas
#nspp <- length(look_seas@species_params$species)
#sp <- look_seas@species_params$species
#ngrid <- length(look_seas@w)
#dt <- look_seas$dt
#w <- look_seas@w
#dw <- look_seas@dw
#wFull <- look_seas@w_full
#dwFull <- look_seas@dw_full
#alpha <- sp$alpha
#eRepro <- sp$erepro
#theta <- param$interaction


######################################


#now we project the models forward in time, 500 years, gather data matrices


####Seasonal
model_sim <- project(look_seas, t_max = tmax, t_save = dt*isave) #project in time


#initialise storage arrays
gg_mat <- array(0, dim = c(nsave, nspp, ngrid))  # dimensions: time, species, weight bins
pm_mat <- array(0, dim = c(nsave, nspp, ngrid))
total_mort <- array(0, dim = c(nsave, nspp, ngrid))
nm_mat <- array(0, dim = c(nsave,n_pp_grid))


#
itimemax  <- tmax / dt 
isav <- 0  # counter for saving

for (itime in 1:itimemax) {
  gg <- getEGrowth(look_seas, n = N(model_sim)[itime,,], n_pp = NResource(model_sim)[itime,])#growth
  
  pm <- getPredMort(look_seas, n = N(model_sim)[itime,,], n_pp = NResource(model_sim)[itime,]) #predation mortality
  
  tm<- getMort(look_seas, n = N(model_sim)[itime,,], n_pp = NResource(model_sim)[itime,]) #total mortality
  
  nm <- getResourceMort(look_seas, n = N(model_sim)[itime,,], n_pp = NResource(model_sim)[itime,]) #resource mortality
  
  
  if ((itime %% isave) == 0) {
    isav <- isav + 1
    
    # Assign gg to the corresponding slice of gg_mat; likewise for rest
    gg_mat[isav, , ] <- gg  # Time slice (isav), species, weight bins
    
    pm_mat[isav, , ]<- pm
  
    total_mort[isav, , ] <- tm
    
    nm_mat[isav,  ] <- nm
  }
}
##debug checkpoint
#print(dim(gg_mat)) 





#####Nonseasonal
model_sim_nonseas <- project(look_nonseas, t_max = tmax, t_save = dt*isave) #project

#initialise
gg_mat_ns <- array(0, dim = c(nsave, nspp, ngrid))  # dimensions: time, species, weight bins
pm_mat_ns <- array(0, dim = c(nsave, nspp, ngrid))
total_mort_ns <- array(0, dim = c(nsave, nspp, ngrid))
nm_mat_ns <- array(0, dim = c(nsave, n_pp_grid))


#
itimemax  <- tmax / dt 
isav <- 0  # Counter for saving

for (itime in 1:itimemax) {
  gg <- getEGrowth(look_nonseas, n = N(model_sim_nonseas)[itime,,], n_pp = NResource(model_sim_nonseas)[itime,]) #growth
  
  pm <- getPredMort(look_nonseas, n = N(model_sim_nonseas)[itime,,], n_pp = NResource(model_sim_nonseas)[itime,]) #predation mortality
  
  tm<- getMort(look_nonseas, n = N(model_sim_nonseas)[itime,,], n_pp = NResource(model_sim_nonseas)[itime,]) #total mortality
  
  nm <- getResourceMort(look_nonseas, n = N(model_sim_nonseas)[itime,,], n_pp = NResource(model_sim_nonseas)[itime,]) #resource mortality
  
  
  if ((itime %% isave) == 0) {
    isav <- isav + 1
    
    # assign gg to the corresponding slice of gg_mat etc.
    gg_mat_ns[isav, , ] <- gg  # Time slice (isav), species, weight bins
    pm_mat_ns[isav, , ]<- pm
    total_mort_ns[isav, , ] <- tm
    nm_mat_ns[isav, ] <- nm
  }
}

############################
#now we want to do the growth curve calculations as in Datta et al.


#definitions
growth_st <- gg_mat_ns #gg_mat_ns from the nonseasonal model
growth <- gg_mat ##gg_mat from seasonal model
laststep = length(model_sim@n[,1,1])
pred_st = pm_mat_ns #PREDATION MORTALITY ON FISH SPECIES nonseasonal
pred = pm_mat #predation mortality on fish species seasonal
natmort_st = nm_mat_ns #natural mortality nonseasonal
natmort = nm_mat #natural mortality seasonal





#Now we take the following segment directly from Datta et al. code:
# LOOK AT GROWTH AND SURVIVAL
T = 15/dt; # number of time steps you want to follow cohort for
t_start = laststep - T; # setting the start time for tracking the cohort
z_st = array(0, c(12, T+1)); # row vector for following cohort weight
dc_st = array(0, c(12, T+1)); # vector for cohort survival
z = array(0, c(12, T+1)); # row vector for following cohort weight
dc = array(0, c(12, T+1)); # vector for cohort survival

# NEWBORNS OVER LIFETIME
z_st[,1] = w[1]; # log weight initially (newborn)
dc_st[,1] = model_sim_nonseas@n[t_start,,1]; # initial population in spectrum
z[,1] = w[1]; # same for seasonal system
dc[,1] = model_sim@n[t_start,,1];

for (q in seq(1,12)){ 
  for (t in seq(1,T)){ # within time period you're interested in
    zy = max(which(z[q,t] - w >= 0)) # weight bin of cohort from last time step (this will probably need updating for your code)
    z[q,t+1] = z[q,t]+dt*growth[t_start+t-1,q,zy] # using growth rate in that bin to update to z(t-t_start+1)
    dc[q,t+1] = dc[q,t]*exp(-dt*(pred[t_start+t-1,q,zy]+natmort[t_start+t-1,30+zy])) # updating amount surviving using death rate
    zy = max(which(z_st[q,t] - w >= 0)) # weight bin of cohort from last time step (this will probably need updating for your code)
    z_st[q,t+1] = z_st[q,t]+dt*growth_st[t_start+t-1,q,zy] # using growth rate in that bin to update to z(t-t_start+1)
    dc_st[q,t+1] = dc_st[q,t]*exp(-dt*(pred_st[t_start+t-1,q,zy]+natmort_st[t_start+t-1,30+zy])) # updating amount surviving using death rate
  }
}


####From here we can make the plot
s_order = c(1:3, 5, 4, 8, 7, 6, 9:12) 
specs<-look_seas@species_params$species #species
 


#Now we take the following plotting code directly from Datta et al. code:
#PLOT
par(mfrow=c(4,3),mai=c(0,0,0,0),omi=c(1,1,1,1))
for (ispec in s_order) {
  plot(dt*seq(0,T),z_st[ispec,],type="l",lty=2,lwd=2,xaxt="n",yaxt="n",ylim=c(range(z[ispec,],z_st[ispec,])),col="grey")
  lines(dt*seq(0,T),z[ispec,],type="l",lty=1,lwd=2)
  mtext(specs[ispec],3,line=-1,cex=0.8)
  
  # add axes to outer plots    
  if (is.element(ispec,c(1,5,7,10))) {
    axis(side=2, at = NULL, labels = TRUE, tick = TRUE,cex.axis=0.8,hadj=TRUE)
  }	else {
    axis(side=4, at = NULL, labels = TRUE, tick = TRUE,cex.axis=0.8,hadj=TRUE)   
  }
  if (is.element(ispec,c(10,11,12)))	axis(side=1, at = NULL, labels = TRUE, tick = TRUE,cex.axis=0.8,hadj=TRUE)
}
# add axes labels on outside
mtext("Time (years)",1,outer=T,line=3)
mtext(expression(paste("Mass (g)",sep="")),2,outer=T,line=3)



###so we have outputted the growth curves for offspring of different species hatching into spectra in
#both nonseasonal (dashed grey lines) and seasonal (solid black lines) systems, over 15 years
#but this time using the mizerSeasonal encoded gonadic mass based seasonality 
##aka here is capital breeding rather than Datta et al. income breeding


