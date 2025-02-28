

testGrowthTry <- function(object, species = NULL, t_max = 1, dt = 0.01, t_save = 0.01) {
  if (!is(object, "MizerSim")) {
    stop("not a sim object")
  }
  
  
  #input simulation, take its parameters, then set initial conditions as those at final timestep
  params <- object@params
  params <- setInitialValues(params, object)
  
  #which species included (sets to all)
  species <- valid_species_arg(params, species)
  
  
  
  #project forward until t_max, with timesteps dt
  #note: dt needs to not be less than the dt used in initial 500-yr simulation projection
  #NECESSARY: t_save, otherwise will just store annually
  sim <- project(params, t_max = t_max, dt = dt,  t_save = t_save)
  
  
  
  #extract the vector of timesteps in sim
  time_steps <- as.numeric(dimnames(sim@n)[[1]])  #array of timesteps, t_min to t_max by dt(t_save)
  
  
  #initialise empty array to store the results
  ws <- array(dim = c(length(species), length(time_steps)),
              dimnames = list(Species = species, Time = time_steps))
  
  
  #extract abundance arrays from sim
  N_vals <- sim@n  # (time × species × size) #fish species
  N_pp_vals <- sim@n_pp # (time × size) #resource
  
  
  
  
  
  #function to get growth at each timestep
  get_growth_at_t <- function(t) {
    
    
    
    # Call getEGrowth with correctly shaped input
    #have to take N_t and N_pp_t at each of the t
    
    #note: need to declare t=something when doing the checks
    n_t <- (N_vals[t,,])
    n_pp_t <- (N_pp_vals[t,])
    g_t <- getEGrowth(params, n = n_t, n_pp = n_pp_t)
    
    # Interpolate growth rate across size classes 
    g_interp <- stats::approxfun(c(params@w, params@species_params$w_max[[i]]),
                       c(g_t[i, ], 0))
    
    return(g_interp)
    
    #}
    
    
  }
  
  
  #loop for all species
  for (j in seq_along(species)) {
    i <- which(params@species_params$species == species[j]) 
    
    
    #define ODE for weight growth
    myodefun <- function(t, state, parameters) {
      
      g_fn <- get_growth_at_t(t)
      
      
      return(list(g_fn(state)))
      
      
    }
    
    # Solve for weight growth over time
    ws[j, ] <- deSolve::ode(y = params@w[params@w_min_idx[i]], 
                            times = seq_along(time_steps), 
                            func = myodefun
                            )[, 2]
  }
  
  
  
  
  return(ws)
  
  
}


test_growth <- testGrowthTry(datta_500yr, t_max = 1, dt = 0.01, t_save = 0.01)


# Convert to long format for plotting
test_growth_long <- melt(test_growth, varnames = c("Species", "Time"), value.name = "Weight")

ggplot(test_growth_long, aes(x = Time, y = Weight, color = Species)) +
  geom_line() +  # Line plot
  scale_y_log10() +  # Log scale for y-axis
  labs(
    x = "Age (years)", 
    y = "Weight (g)", 
    color = "Species", 
    title = "mizerSeasonal Seasonal Growth Curves"  # Title
  ) +  # Labels and title
  theme_minimal() +  # Clean theme
  theme(
    legend.position = "none",  # Hide the legend
    plot.title = element_text(hjust = 0.5, size = 16)  # Center the title and adjust font size
  ) +
  facet_wrap(~ Species, ncol = 3)  # Facet by Species (12 species, arranged in 3 columns)


