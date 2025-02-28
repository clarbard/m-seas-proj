




# Function to compute growth rate directly for each timestep and species
get_growth_rate <- function(t, sp, state) {
  # Find the closest time step
  time_idx <- which.min(abs(time_steps - t))  # Closest timestep
  
  # Extract abundance at the current time
  n_t <- N_vals[time_idx,,]  # Abundance for each species and size class at current time
  n_pp_t <- N_pp_vals[time_idx,]  # Primary producer biomass at current time
  
  # Compute the growth rate directly (dw/dt) using the model dynamics
  g_t <- getEGrowth(params, n = n_t, n_pp = n_pp_t)[sp, ]
  #names(g_t) <- NULL
  
  # Create interpolation function for growth rate
 # g_interp <- stats::approxfun(c(params@w, params@species_params$w_max[i]),
  #                             c(g_t, 0))
  
  # Evaluate the growth rate at the current weight (state)
 # g_t_2 <- g_interp(state)  # Evaluate the interpolated function at the current state
  

  state_indx <- which(params@w == state)
  
  g_t_2 <- g_t[state_indx]
  names(g_t_2) <- NULL
  

  
  
  return(g_t_2)
}

# Main function to simulate growth
testGrowthTry <- function(object, species = NULL, t_max = 1, dt = 0.01, t_save = 0.01) {
  if (!is(object, "MizerSim")) {
    stop("not a sim object")
  }
  
  # Input simulation, extract parameters, set initial conditions
  params <- object@params
  params <- setInitialValues(params, object)
  
  species <- valid_species_arg(params, species)
  
  # Project forward with small timesteps to capture seasonality
  sim <- project(params, t_max = t_max, dt = dt, t_save = t_save)
  
  time_steps <- as.numeric(dimnames(sim@n)[[1]])
  
  # Initialize storage for species' weight over time
  ws <- array(dim = c(length(species), length(time_steps)),
              dimnames = list(Species = species, Time = time_steps))
  
  # Extract population size and primary producer biomass from simulation results
  N_vals <- sim@n  # (time × species × size)
  N_pp_vals <- sim@n_pp # (time × size)
  
  # Loop through species to calculate growth rate and simulate weight growth
  for (j in seq_along(species)) {
    i <- which(params@species_params$species == species[j])  # Get the index of the species
    
    # Define ODE function for weight growth
    myodefun <- function(t, state, parameters) {
      # Get the growth rate directly at time t
      g_t_2 <- get_growth_rate(t, i, state)
      
      # Compute dw/dt directly as the growth rate, ensuring it is a scalar
      return(list(g_t_2 * state))  # Growth rate multiplied by the current state (weight)
    }
    
    # Set initial weight for the species (scalar)
    initial_weight <- params@w[params@w_min_idx[i]]
    
    # Solve the ODE for weight growth over time with dynamically computed growth
    result <- deSolve::ode(y = initial_weight, 
                           times = time_steps, 
                           func = myodefun, 
                           parms = NULL)
    
    # Store the results (weight over time)
    ws[j, ] <- result[, 2]  # The second column contains the solution for weight
  }
  
  return(ws)
}

# Run test with seasonal variation included
test_growth <- testGrowthTry(datta_500yr, t_max = 1, dt = 0.01, t_save = 0.01)

# Convert to long format for plotting
test_growth_long <- melt(test_growth, varnames = c("Species", "Time"), value.name = "Weight")

# Plot the results
ggplot(test_growth_long, aes(x = Time, y = Weight, color = Species)) +
  geom_line() +
  scale_y_log10() +
  labs(
    x = "Age (years)", 
    y = "Weight (g)", 
    color = "Species", 
    title = "mizerSeasonal Seasonal Growth Curves"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none", 
    plot.title = element_text(hjust = 0.5, size = 16)
  ) +
  facet_wrap(~ Species, ncol = 3)

