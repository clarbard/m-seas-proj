#replicating a few more datta figures

##copying pasting from my growthcurvesreplicating script
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
model_sim <- project(look_seas, t_max = tmax, t_save = dt*isave) #project in time --seasonal
model_sim_nonseas <- project(look_nonseas, t_max = tmax, t_save = dt*isave) #project -- nonseasonal





#############
####Seasonal

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



#definitions
growth_st <- gg_mat_ns #gg_mat_ns from the nonseasonal model
growth <- gg_mat ##gg_mat from seasonal model
laststep = length(model_sim@n[,1,1])
pred_st <- pm_mat_ns #PREDATION MORTALITY ON FISH SPECIES nonseasonal
pred <- pm_mat #predation mortality on fish species seasonal
natmort_st <- nm_mat_ns #natural mortality nonseasonal
natmort <- nm_mat #natural mortality seasonal



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




#######################

#replicate figure 3


par(mfrow=c(4,3),mai=c(0,0,0,0),omi=c(1,1,1,1))
specs<-look_seas@species_params$species
s_order = c(1:3, 5, 4, 8, 7, 6, 9:12)
months = laststep - (1/dt) + round(seq(1,12)/(12*dt)) # indices in a year for monthly snapshots

for (ispec in s_order) {
  plot(w,N(model_sim)[laststep,ispec,],log="xy",typ="l",lwd=2,xaxt="n",yaxt="n",ylim=c(1e-2,1e15),col="grey")
  for (i in 1:11) lines(w,N(model_sim)[months[i],ispec,],typ="l",lty=1,lwd=2,col="grey")
  lines(w,N(model_sim)[laststep,ispec,],log="xy",typ="l",lty=2,lwd=2,xaxt="n",yaxt="n",ylim=c(1e-2,1e15))
  # for (i in 1:12) points(baseModel2$w,baseModel$N[laststep-i,baseModel$param$species$species==specs[ispec],],typ="l",col="black")
  mtext(specs[ispec],3,line=-1,cex=0.8)
  
  # add axes to outer plots    
  if (is.element(ispec,c(1,4,7,10)))	axis(side=2, at = NULL, labels = TRUE, tick = TRUE,cex.axis=0.8,hadj=TRUE)
  if (is.element(ispec,c(10,11,12)))	axis(side=1, at = NULL, labels = TRUE, tick = TRUE,cex.axis=0.8,hadj=TRUE)
}

# add axes labels on outside
mtext("Body mass (g)",1,outer=T,line=3)
mtext(expression(paste("Number density ( ", m^-3*g^-1, ")",sep="")),2,outer=T,line=3)



#################################

##figure 4

par(mfrow=c(4,3),mai=c(0,0,0,0),omi=c(1,1,1,1))
cl <- rainbow(12)
specs<-look_seas@species_params$species
s_order = c(1:3, 5, 4, 8, 7, 6, 9:12)



for (ispec in s_order) {
  bla = which(model_sim@n[laststep,ispec,] < 1) # set limit
  pick = bla[1]
  plot(w[1:pick],model_sim_nonseas@n[laststep,ispec,1:pick]/model_sim_nonseas@n[laststep,ispec,1:pick],log="xy",typ="l",lty=1,lwd=2,xaxt="n",yaxt="n",ylim=c(5e-1,5))
  for (i in 1:12) points(w[1:pick],model_sim@n[months[i],ispec,1:pick]/model_sim@n[laststep,ispec,1:pick],typ="l",col=cl[i])
  lines(w[1:pick],model_sim_nonseas@n[laststep,ispec,1:pick]/model_sim_nonseas@n[laststep,ispec,1:pick],typ="l",lty=1,lwd=2) # re-do horizontal line
  mtext(specs[ispec],3,line=-1,cex=0.8)
  
  # add axes to outer plots    
  if (is.element(ispec,c(1,5,7,10)))	axis(side=2, at = NULL, labels = TRUE, tick = TRUE,cex.axis=0.8,hadj=TRUE)
  if (is.element(ispec,c(10,11,12)))	axis(side=1, at = NULL, labels = TRUE, tick = TRUE,cex.axis=0.8,hadj=TRUE)
  if (ispec == 1) legend("topright", legend = c("J","F","M","A","M","J","J","A","S","O","N","D"), col=cl, lty=1, ncol=2, bty="n") # optional legend
}
# add axes labels on outside
mtext("Body mass (g)",1,outer=T,line=3)
mtext(expression(paste("Seasonal density / non-seasonal density ",sep="")),2,outer=T,line=3)



###################################

#survival --doesn't look good
par(mfrow=c(4,3))
for (ispec in 1:12) {
  plot(dt*seq(0,T),dc_st[ispec,],type="l",lty=2,lwd=2,xaxt="n",yaxt="n",ylim=c(1e-4,max(dc[ispec,],dc_st[ispec,])))
  lines(dt*seq(0,T),dc[ispec,],type="l",lwd=2,col="grey")
  mtext(specs[ispec],3,line=-1,cex=0.8)
  
  # add axes to outer plots    
  if (is.element(ispec,c(1,5,9)))	axis(side=2, at = NULL, labels = TRUE, tick = TRUE,cex.axis=0.8,hadj=TRUE)
  if (is.element(ispec,c(length(look_seas@species_params$species)-3,length(look_seas@species_params$species)-2,length(look_seas@species_params$species)-1,length(look_seas@species_params$species))))	axis(side=1, at = NULL, labels = TRUE, tick = TRUE,cex.axis=0.8,hadj=TRUE)
}

# add axes labels on outside
mtext("Time (years)",1,outer=T,line=3)
mtext(expression(paste("Biomass remaining",sep="")),2,outer=T,line=3)












###figure 4, spawning in a year try
par(mfrow=c(1,1))

pick = 7 # choose species (7 = sole)

t_pick = 11 # year to start 

##already have z, z_st defined previously

#taken directly from datta
r_peak_s = c(3.6047, 2.9994, 1.944, 0.40493, 1.141, 1.4257, 4.973, 0.795, 4.0495, 3.6567, 5.4732, 1.951) # reproduction peakinesses
r_peakt = c(0.54765, 0.064324, 0.3574, 0.85759, 0.38245, 0.37164, 0.48564, 0.50213, 0.22454, 0.41806, 0.31236, 0.33333) # time of peak
#can i use these?(above)

t_v=seq(t_pick/dt,t_pick/dt+1/dt) # 13 time steps

plot(seq(0,1,dt), z[pick,t_v],type="l",lty=1,lwd=4,ylab="Weight (g)",xlab="Time (years)",
     ylim=c(min(z_st[pick,t_v],z[pick,t_v]),max(z_st[pick,t_v],z[pick,t_v])),
     cex.lab=1.2, cex.axis=1.2, cex.main=1.2, cex.sub=1.2)
lines(seq(0,1,dt),z_st[pick,t_v],type="l",lty=2,lwd=2,col="grey")
lines(seq(0,1,1/1000),min(z_st[pick,t_v],z[pick,t_v])+(max(z_st[pick,t_v],z[pick,t_v])-min(z_st[pick,t_v],z[pick,t_v+12]))*
        exp(r_peak_s[pick]*cos(2*pi*(seq(0,1,1/1000)-r_peakt[pick])))/exp(r_peak_s[pick]),type="l",lty=3,lwd=2)









####Plot feeding level



#first calculate feeding level
f_ns <-  array(0, dim = c(nsave, nspp, ngrid))
f_s <-  array(0, dim = c(nsave, nspp, ngrid))
#
itimemax  <- tmax / dt 
isav <- 0  # Counter for saving
for (itime in 1:itimemax) {
  fl <- getFeedingLevel(look_nonseas, n = N(model_sim_nonseas)[itime,,], n_pp = NResource(model_sim_nonseas)[itime,]) #feeding level nonseasonal
  fl_s <- getFeedingLevel(look_seas, n = N(model_sim)[itime,,], n_pp = NResource(model_sim)[itime,]) #feeding level seasonal
  
  if ((itime %% isave) == 0) {
    isav <- isav + 1
  
    f_ns[isav, , ] <- fl  # Time slice (isav), species, weight bins
    f_s[isav, , ] <- fl_s

  }
}




##plot the seasonal feeding levels
f <- f_s

#if(is.na(meantsteps))
  y <- f_s[dim(f_s)[1],,]
#else # or mean of last meantsteps
#  y <- apply(f[(dim(f)[1]-meantsteps+1):dim(f)[1],,],c(2,3),mean)


col <- rep(rainbow(ceiling(dim(f)[2]/2),start=0,end=5/6),2)


lty <- rep(c(1,2), each = ceiling(dim(f)[2]/2))
# Get y lims
ylim <- c(0,1)
# Plot empty
plot(x=look_seas@w, y=y[1,], log="x",ylim=ylim, type="n", xlab = "mass (g)", ylab = "feeding level")
# Plot the lines
# Only plot w <= winf
for (i in 1:nspp)
  points(x=look_seas@w[look_seas@w <= look_seas@species_params$w_inf[i]], y=y[i,look_seas@w <= look_seas@species_params$w_inf[i]],
         col=col[i], type="l", lty=lty[i])


#legend(x = "bottomright" , legend = as.character(look_seas@species_params$species),
#         col= col, lty=lty, cex=0.7, ncol=2)






#plot natural mortality levels -- seasonal
##M2 is pred
#for seasonl is pred, nonseasonal pred_st
y <- pred[dim(pred)[1],,]
f <- f_s

col <- rep(rainbow(ceiling(dim(f)[2]/2),start=0,end=5/6),2)

lty <- rep(c(1,2), each = ceiling(dim(f)[2]/2))
# Get y lims
ylim <- c(0,max(y))

# Plot empty
plot(x=look_seas@w, y=y[1,], log="x",ylim=ylim, type="n", xlab = "mass (g)", ylab = "natural mortality")
# Plot the lines
for (i in 1:nspp)
  points(x=look_seas@w[look_seas@w<= look_seas@species_params$w_inf[i]], y=y[i,look_seas@w<= look_seas@species_params$w_inf[i]],
         col=col[i], lty=lty[i], type="l")

#legend(x = "bottomright" , legend = as.character(look_seas@species_params$species),
#       col= col, lty=lty, cex=0.7, ncol=2)







###plot abundance -- seasonal

N <- model_sim@n[dim(model_sim@n)[1],,]


# w is in grams
#spBiomass <- sweep(N,2,model$dw,"*") / 1e6
# Want to calculate biomass at each point (not in each size class)
# N is abundance at point, not abundance in that size class
spN <- N / 1e6
PPN <- model_sim@n_pp[dim(model_sim@n_pp)[1],] / 1e6
refSpec <- model_sim@params@resource_params$kappa*look_seas@w_full^(-model_sim@params@resource_params$lambda) / 1e6



col <- rep(rainbow(ceiling(dim(f)[2]/2),start=0,end=5/6),2)


lty <- rep(c(1,2), each = ceiling(dim(f)[2]/2))
# Set xlim and ylim values
xlim <- c(1e-2, max(look_seas@w_full)) # don't plot the full resource spectrum
ylim <- c(1e-6,max(PPN,spN,refSpec))

# Set it up so x is full spectrum
plot(x=look_seas@w_full, y=PPN, log="xy",type="n", ylim=ylim, xlim=xlim, xlab = "mass (g)", ylab = "Abundance")

# Reference spectrum:
points(x=look_seas@w_full, y= refSpec, type="l", lty=3, col = 1)
# Resource spectrum:
points(x=look_seas@w_full, y=PPN, type="l", lty=3, col = 3)

# Species spectrum
for (i in 1:dim(spN)[1])
  points(x=look_seas@w, y = spN[i,], col=col[i], type="l", lty=lty[i])



#legend(x = "bottomright" , legend = as.character(look_seas@species_params$species),
#       col= col, lty=lty, cex=0.7, ncol=2)










###total biomass --seasonal [loooks weird, the pink dashed]

trange <- 1:dim(model_sim@n)[1]

# N * w then summed across w to get total biomass through time
#spBiomass <- apply(sweep(model$N,c(3,2),model$w / 1e6,"*"),c(1,2),sum)
spBiomass <- rowSums(sweep(model_sim@n,3,look_seas@w*look_seas@dw,"*"),dims=2)

col <- rep(rainbow(ceiling(dim(f)[2]/2),start=0,end=5/6),2)
lty <- rep(c(1,2), each = ceiling(dim(f)[2]/2))

ylim <- c(1,max(spBiomass))
plot(x=trange,y=trange,type="n",ylim=ylim,log="y",xlab = "timestep", ylab = "total biomass (t)")  

for (i in 1:dim(spBiomass)[2])
  points(x=trange, y = spBiomass[trange,i], col=col[i], lty=lty[i], type="l")


#legend(x = "bottomright" , legend = as.character(look_seas@species_params$species),
 #      col= col, lty=lty, cex=0.7, ncol=2)























########Replicating figure 2
par(mfrow=c(1,1),mai=c(0,0,0,0),omi=c(1,1,1,1))

cspectrum<-apply(model_sim_nonseas@n,c(1,3),sum)
cspectrum2<-apply(model_sim@n,c(1,3),sum)


plot(w,cspectrum2[laststep,],log="xy",type="l",lwd=4,ylab="Number density",xlab="Body mass (g)",xlim=c(1e-4,1e5),ylim=c(1e-5,1e18),col="grey",
     cex.lab=1.2, cex.axis=1.2, cex.main=1.2, cex.sub=1.2)
points(w,cspectrum[laststep,],type="l",lty=2,lwd=4)
points(look_seas@w_full,model_sim@n_pp[laststep,],type="l",lty=3,lwd=2)


plot(w,cspectrum[laststep,],log="xy",typ="l",lwd=3,ylab="Number density",xlab="Body mass (g)",ylim=c(1e-5,1e15),
     cex.lab=1.2, cex.axis=1.2, cex.main=1.2, cex.sub=1.2)
for(i in 1:12) lines(w,model_sim_nonseas@n[laststep,i,],log="xy",lty=2,xlim=c(1e-5,1e5),ylim=c(1e-3,1e25),col="grey")





plot(w,cspectrum2[laststep,],log="xy",typ="l",lwd=3,ylab="Number density",xlab="Body mass (g)",ylim=c(1e-5,1e15),
     cex.lab=1.2, cex.axis=1.2, cex.main=1.2, cex.sub=1.2)
for(i in 1:12) {
  lines(w,model_sim@n[laststep,i,],log="xy",lty=2,xlim=c(1e-5,1e5),ylim=c(1e-3,1e25),col="grey")
  bla = which(model_sim@n[laststep,i,] < 1e-6)
  points(w[bla[1]-1],model_sim@n[laststep,i,bla[1]-1],pch=15, col="black")
}
