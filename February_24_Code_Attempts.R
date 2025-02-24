###February 24th Coding Attempts


library(mizer)
library(mizerSeasonal)
library(mizerExperimental)
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggplot2)


#must have datta_500yr and datta_nonseasonal_500yr loaded [see replicating_datta_figures.R]






##starting with teh attemptt hat worked best yesterday


getGrowthCurvesRedo <- function(object, 
                                species = NULL,
                                max_age = 20,
                                # percentage = FALSE
) {
  
  ##extract params and set initial values if sim, else validate params
  #
  if (is(object, "MizerSim")) {
    params <- object@params
    params <- setInitialValues(params, object)
  } else if (is(object, "MizerParams")) {
    params <- validParams(object)
  } else {
    stop("The first argument to `getGrowthCurvesRedo()` must be a ",
         "MizerParams or a MizerSim object.")
  }
  
  
  
  ##reorder species to match listing in params
  #
  species <- valid_species_arg(params, species)
  # reorder list of species to coincide with order in params
  idx <- which(params@species_params$species %in% species)
  species <- params@species_params$species[idx]
  
  
  #1000 evenly spaced age points
  age <- seq(0, max_age, length.out = 1000)
  
  
  #initialises empty matrix for wieght/age
  ws <- array(dim = c(length(species), length(age)),
              dimnames = list(Species = species, Age = age))
  
  #get the growth rates [HERE THE EDITS?]
  g <- getEGrowth(params)
  
  
  
  for (j in seq_along(species)) {
    i <- idx[j]
    
    #interpolates growth rate as function of weight (piecewise linear approx)
    g_fn <- stats::approxfun(c(params@w, params@species_params$w_max[[i]]),
                             c(g[i, ], 0))
    
    #rate of change of weight = interpolated growth rate above
    myodefun <- function(t, state, parameters) {
      return(list(g_fn(state)))
    }
    
    #numerically solves for weight over time
    ws[j, ] <- deSolve::ode(y = params@w[params@w_min_idx[i]], 
                            times = age, func = myodefun)[, 2]
    
    
    
    
    # if (percentage) {
    #  ws[j, ] <- ws[j, ] / params@species_params$w_max[i] * 100
    #}
  }
  return(ws)
}



#get the growth curves data
test1_s <- getGrowthCurvesRedo(datta_500yr, max_age = 15)
test1_ns <- getGrowthCurvesRedo(datta_nonseasonal_500yr, max_age = 15)


# Reshape the growth curves into long format using melt
# This will create a data frame with Age, Species, and Value columns
test1_s_long <- melt(test1_s, varnames = c("Species", "Age"), value.name = "value")
test1_ns_long <- melt(test1_ns, varnames = c("Species", "Age"), value.name = "value")


























###trying to refine it


getGrowthCurvesRedo <- function(object, 
                                species = NULL,
                                max_age = 20 ) {
  
  ##extract params and set initial values if sim, else validate params
  #
  if (is(object, "MizerSim")) {
    params <- object@params
    params <- setInitialValues(params, object)
  }
  #else if (is(object, "MizerParams")) {
  #   params <- validParams(object)
  # } else {
  #   stop("The first argument to `getGrowthCurvesRedo()` must be a ",
  #        "MizerParams or a MizerSim object.")
  # }
  
  
  
  ##reorder species to match listing in params
  #
  species <- valid_species_arg(params, species)
  # reorder list of species to coincide with order in params
  idx <- which(params@species_params$species %in% species)
  species <- params@species_params$species[idx]
  
  
  #1000 evenly spaced age points
  age <- seq(0, max_age, length.out = 1000*max_age)  # More time points
  
  
  
  #initialises empty matrix for wieght/age
  ws <- array(dim = c(length(species), length(age)),
              dimnames = list(Species = species, Age = age))
  
  #get the growth rates [HERE THE EDITS?]
  g <- getEGrowth(params)
  
  
  for (j in seq_along(species)) {
    i <- idx[j]
    
    #interpolates growth rate as function of weight (piecewise linear approx)
    g_fn <- stats::approxfun(seq(min(params@w), params@species_params$w_max[[i]], length.out = 1000*max_age),
                             approx(c(params@w, params@species_params$w_max[[i]]), c(g[i, ], 0), 
                                    xout = seq(min(params@w), params@species_params$w_max[[i]], length.out = 1000*max_age))$y)
    ###adding mroe weight points
    
    #rate of change of weight = interpolated growth rate above
    myodefun <- function(t, state, parameters) {
      return(list(g_fn(state)))
    }
    
    #numerically solves for weight over time
    ws[j, ] <- deSolve::ode(y = params@w[params@w_min_idx[i]], 
                            times = seq(0, max_age, length.out = 1000*max_age),  # More time points
                            func = myodefun, 
                            method = "lsoda")[, 2]  # Use "lsoda" for adaptive step sizes
    
    
    
    
    
  }
  return(ws)
}



#get the growth curves data
test1_s <- getGrowthCurvesRedo(datta_500yr, max_age = 15)
test1_ns <- getGrowthCurvesRedo(datta_nonseasonal_500yr, max_age = 15)


# Reshape the growth curves into long format using melt
# This will create a data frame with Age, Species, and Value columns
test1_s_long <- melt(test1_s, varnames = c("Species", "Age"), value.name = "value")
test1_ns_long <- melt(test1_ns, varnames = c("Species", "Age"), value.name = "value")





























######plot
##SEASONAL
ggplot(test1_s_long, aes(x = Age, y = value, color = Species)) +
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




##NONSEASONAL
ggplot(test1_ns_long, aes(x = Age, y = value, color = Species)) +
  geom_line() +  # Line plot
  scale_y_log10() +  # Log scale for y-axis
  labs(
    x = "Age (years)", 
    y = "Weight (g)", 
    color = "Species", 
    title = "mizerSeasonal NONSeasonal Growth Curves"  # Title
  ) +  # Labels and title
  theme_minimal() +  # Clean theme
  theme(
    legend.position = "none",  # Hide the legend
    plot.title = element_text(hjust = 0.5, size = 16)  # Center the title and adjust font size
  ) +
  facet_wrap(~ Species, ncol = 3)  # Facet by Species (12 species, arranged in 3 columns)









#

library(mizer)
library(deSolve)


tester <- newSingleSpeciesParams()
plotGrowthCurves(tester)


testGrowthTry <- function(object, species = NULL, max_age = 20){
  
  
  if (is(object, "MizerParams")) {
    params <- validParams(object)
  } else {
    stop("not a params object")
  }
  
  
}










library(mizer)
library(deSolve)
library(ggplot2)
library(reshape2)
library(mizer)
library(deSolve)
library(ggplot2)
library(reshape2)


testGrowthTry <- function(object, species = NULL, t_max = 50, dt = 0.01) {
  if (!is(object, "MizerParams")) {
    stop("not a params object")
  }
  
  params <- validParams(object)
  species <- valid_species_arg(params, species)
  
  # Run the simulation with fine time resolution
  sim <- project(params, t_max = t_max, dt = dt)
  
  # Extract actual time steps available in N(sim)
  available_time_steps <- as.numeric(dimnames(N(sim))[[1]])  # Extract actual times
  time_steps <- seq(0, t_max, by = dt)  # Requested time steps
  
  # Array to store results
  ws <- array(dim = c(length(species), length(time_steps)),
              dimnames = list(Species = species, Time = time_steps))
  
  # Extract abundance arrays
  N_vals <- N(sim)  # (time × species × size)
  N_pp_vals <- NResource(sim)  # (time × size)
  
  # Check dimensions
  print(dim(N_vals))  # Debugging
  
  for (j in seq_along(species)) {
    i <- which(params@species_params$species == species[j])
    
    # Define function to get growth at each time step
    get_growth_at_t <- function(t) {
      # Ensure `t_idx` stays within available bounds
      t_idx <- which.min(abs(available_time_steps - t))
      t_idx <- max(1, min(dim(N_vals)[1], t_idx))  # Keep within limits
      
      # Extract abundance at time t_idx, ensuring correct shape
      n_t <- array(N_vals[t_idx, , ], dim = dim(N_vals)[2:3])  # Force array
      n_pp_t <- N_pp_vals[t_idx, ]  # Ensure (size vector)
      
      # Call getEGrowth with correctly shaped input
      g_t <- getEGrowth(params, n = n_t, n_pp = n_pp_t)
      
      # Interpolate growth rate across size classes
      g_interp <- approx(c(params@w, params@species_params$w_max[[i]]),
                         c(g_t[i, ], 0),
                         xout = params@w)$y
      
      return(stats::approxfun(params@w, g_interp))
    }
    
    # Define ODE for weight growth
    myodefun <- function(t, state, parameters) {
      g_fn <- get_growth_at_t(t)
      return(list(g_fn(state)))
    }
    
    # Solve for weight growth over time
    ws[j, ] <- deSolve::ode(y = params@w[params@w_min_idx[i]], 
                            times = time_steps, 
                            func = myodefun, 
                            method = "lsoda")[, 2]
  }
  
  return(ws)
}



# Test the function
#use params_1000
test_growth <- testGrowthTry(params_1000, t_max = 50, dt = 0.01)

# Convert to long format for plotting
test_growth_long <- melt(test_growth, varnames = c("Species", "Time"), value.name = "Weight")

# Plot growth curves
ggplot(test_growth_long, aes(x = Time, y = Weight, color = Species)) +
  geom_line(size = 1) +
  labs(title = "High-Resolution Growth Curves from Projection", x = "Time (years)", y = "Weight (g)") +
  theme_minimal()





test_growth <- testGrowthTry(params_10000, t_max = 50, dt = 0.01)

# Convert to long format for plotting
test_growth_long <- melt(test_growth, varnames = c("Species", "Time"), value.name = "Weight")

# Plot growth curves
ggplot(test_growth_long, aes(x = Time, y = Weight, color = Species)) +
  geom_line(size = 1) +
  labs(title = "High-Resolution Growth Curves from Projection", x = "Time (years)", y = "Weight (g)") +
  theme_minimal()


ggplot(test_growth_long, aes(x = Time, y = Weight, color = Species)) +
  geom_line(size = 0.8, alpha = 0.5) +  # Light smooth curve
  geom_point(size = 1.5, alpha = 0.8) +  # Raw data points
  scale_y_log10() +  # Set y-axis to logarithmic scale
  scale_x_log10()
labs(title = "Growth Curves with Raw Data Points (Log Scale)", 
     x = "Time (years)", 
     y = "Weight (g)") +
  theme_minimal()

























library(mizer)
library(deSolve)
library(ggplot2)
library(reshape2)
library(mizer)
library(deSolve)
library(ggplot2)
library(reshape2)


testGrowthTry <- function(object, species = NULL, t_max = 15, dt = 0.01) {
  params <- object@params
  params <- setInitialValues(params, object)
  species <- valid_species_arg(params, species)
  
  sim <- object
  
  # Extract actual time steps available in N(sim)
  available_time_steps <- as.numeric(dimnames(N(sim))[[1]])  # Extract actual times
  time_steps <- seq(0, t_max, by = dt)  # Requested time steps
  
  # Array to store results
  ws <- array(dim = c(length(species), length(time_steps)),
              dimnames = list(Species = species, Time = time_steps))
  
  # Extract abundance arrays
  N_vals <- N(sim)  # (time × species × size)
  N_pp_vals <- NResource(sim)  # (time × size)
  
  # Check dimensions
  print(dim(N_vals))  # Debugging
  
  for (j in seq_along(species)) {
    i <- which(params@species_params$species == species[j])
    
    # Define function to get growth at each time step
    get_growth_at_t <- function(t) {
      # Ensure `t_idx` stays within available bounds
      t_idx <- which.min(abs(available_time_steps - t))
      t_idx <- max(1, min(dim(N_vals)[1], t_idx))  # Keep within limits
      
      # Extract abundance at time t_idx, ensuring correct shape
      n_t <- array(N_vals[t_idx, , ], dim = dim(N_vals)[2:3])  # Force array
      n_pp_t <- N_pp_vals[t_idx, ]  # Ensure (size vector)
      
      # Call getEGrowth with correctly shaped input
      g_t <- getEGrowth(params, n = n_t, n_pp = n_pp_t)
      
      # Interpolate growth rate across size classes
      g_interp <- approx(c(params@w, params@species_params$w_max[[i]]),
                         c(g_t[i, ], 0),
                         xout = params@w)$y
      
      return(stats::approxfun(params@w, g_interp))
    }
    
    # Define ODE for weight growth
    myodefun <- function(t, state, parameters) {
      g_fn <- get_growth_at_t(t)
      return(list(g_fn(state)))
    }
    
    # Solve for weight growth over time
    ws[j, ] <- deSolve::ode(y = params@w[params@w_min_idx[i]], 
                            times = time_steps, 
                            func = myodefun, 
                            method = "lsoda")[, 2]
  }
  
  return(ws)
}



# Test the function
#use params_1000
test_growth <- testGrowthTry(simulation_1, t_max = 15, dt = 0.01)

# Convert to long format for plotting
test_growth_long <- melt(test_growth, varnames = c("Species", "Time"), value.name = "Weight")

# Plot growth curves
ggplot(test_growth_long, aes(x = Time, y = Weight, color = Species)) +
  geom_line(size = 1) +
  labs(title = "High-Resolution Growth Curves from Projection", x = "Time (years)", y = "Weight (g)") +
  theme_minimal()




ggplot(test_growth_long, aes(x = Time, y = Weight, color = Species)) +
  geom_point(size = 2, alpha = 0.6) +  # Plot raw data points
  scale_y_log10() +  # Set y-axis to logarithmic scale
  labs(title = "High-Resolution Growth Data Points (Log Scale)", 
       x = "Time (years)", 
       y = "Weight (g)") +
  theme_minimal()


ggplot(test_growth_long, aes(x = Time, y = Weight, color = Species)) +
  geom_line(size = 0.8, alpha = 0.5) +  # Light smooth curve
  geom_point(size = 1.5, alpha = 0.8) +  # Raw data points
  scale_y_log10() +  # Set y-axis to logarithmic scale
  labs(title = "Growth Curves with Raw Data Points (Log Scale)", 
       x = "Time (years)", 
       y = "Weight (g)") +
  theme_minimal()






# Test the function

test_growth <- testGrowthTry(simulation_14, t_max = 15, dt = 0.01)

# Convert to long format for plotting
test_growth_long <- melt(test_growth, varnames = c("Species", "Time"), value.name = "Weight")


ggplot(test_growth_long, aes(x = Time, y = Weight, color = Species)) +
  geom_line(size = 0.8, alpha = 0.5) +  # Light smooth curve
  geom_point(size = 1.5, alpha = 0.8) +  # Raw data points
  scale_y_log10() +  # Set y-axis to logarithmic scale
  labs(title = "Growth Curves with Raw Data Points (Log Scale)", 
       x = "Time (years)", 
       y = "Weight (g)") +
  theme_minimal()


library(mizer)
library(deSolve)
library(ggplot2)
library(reshape2)
library(mizer)
library(deSolve)
library(ggplot2)
library(reshape2)





####MESSING WITH GROWTH, TRYING TO SHOW SEASONALITY IN GROWTH CURVES
###CURRENTLY UNSUCCESSFUL
testGrowthTry <- function(object, species = NULL, t_max = 50, dt = 0.01) {
  if (!is(object, "MizerParams")) {
    stop("not a params object")
  }
  
  params <- validParams(object)
  species <- valid_species_arg(params, species)
  
  # Run the simulation with fine time resolution
  sim <- project(params, t_max = t_max, dt = dt)
  
  # Extract actual time steps available in N(sim)
  available_time_steps <- as.numeric(dimnames(N(sim))[[1]])  # Extract actual times
  time_steps <- seq(0, t_max, by = dt)  # Requested time steps
  
  # Array to store results
  ws <- array(dim = c(length(species), length(time_steps)),
              dimnames = list(Species = species, Time = time_steps))
  
  # Extract abundance arrays
  N_vals <- N(sim)  # (time × species × size)
  N_pp_vals <- NResource(sim)  # (time × size)
  
  # Check dimensions
  print(dim(N_vals))  # Debugging
  
  for (j in seq_along(species)) {
    i <- which(params@species_params$species == species[j])
    
    # Define function to get growth at each time step
    get_growth_at_t <- function(t) {
      # Ensure `t_idx` stays within available bounds
      t_idx <- which.min(abs(available_time_steps - t))
      t_idx <- max(1, min(dim(N_vals)[1], t_idx))  # Keep within limits
      
      # Extract abundance at time t_idx, ensuring correct shape
      n_t <- array(N_vals[t_idx, , ], dim = dim(N_vals)[2:3])  # Force array
      n_pp_t <- N_pp_vals[t_idx, ]  # Ensure (size vector)
      
      # Call getEGrowth with correctly shaped input
      g_t <- getEGrowth(params, n = n_t, n_pp = n_pp_t)
      
      # Interpolate growth rate across size classes
      g_interp <- approx(c(params@w, params@species_params$w_max[[i]]),
                         c(g_t[i, ], 0),
                         xout = params@w)$y
      
      return(stats::approxfun(params@w, g_interp))
    }
    
    # Define ODE for weight growth
    myodefun <- function(t, state, parameters) {
      g_fn <- get_growth_at_t(t)
      return(list(g_fn(state)))
    }
    
    # Solve for weight growth over time
    ws[j, ] <- deSolve::ode(y = params@w[params@w_min_idx[i]], 
                            times = time_steps, 
                            func = myodefun, 
                            method = "lsoda")[, 2]
  }
  
  return(ws)
}






test_growth <- testGrowthTry(datta_params_seasonal, t_max = 50, dt = 0.01)



# Reshape the growth curves into long format using melt
test_with_seasonality_long <- melt(test_growth_long, varnames = c("Species", "Age"), value.name = "value")

# Check the structure of the reshaped data
str(test_with_seasonality_long)

# Plotting the Seasonal Growth Curves
ggplot(test_with_seasonality_long, aes(x = Age, y = value, color = Species)) +
  geom_line() +  # Line plot
  scale_y_log10() +  # Log scale for y-axis
  scale_x_log10() +
  labs(
    x = "Age (years)", 
    y = "Weight (g)", 
    color = "Species", 
    title = "Growth Curves with Seasonal Variation"  # Title
  ) +  # Labels and title
  theme_minimal() +  # Clean theme
  theme(
    legend.position = "none",  # Hide the legend
    plot.title = element_text(hjust = 0.5, size = 16)  # Center the title and adjust font size
  ) +
  facet_wrap(~ Species, ncol = 3)  # Facet by Species (12 species, arranged in 3 columns)



















#try from getting growth curves
test1_s <- getGrowthCurves(datta_500yr)
test1_ns <- getGrowthCurves(datta_nonseasonal_500yr)

# Assuming test1_s and test1_ns are the growth curves obtained from getGrowthCurves()
library(reshape2)
library(ggplot2)

# Reshape the growth curves into long format using melt
# This will create a data frame with Age, Species, and Value columns
test1_s_long <- melt(test1_s, varnames = c("Species", "Age"), value.name = "value")
test1_ns_long <- melt(test1_ns, varnames = c("Species", "Age"), value.name = "value")

##SEASONAL
ggplot(test1_s_long, aes(x = Age, y = value, color = Species)) +
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

##NONSEASONAL
ggplot(test1_ns_long, aes(x = Age, y = value, color = Species)) +
  geom_line() +  # Line plot
  scale_y_log10() +  # Log scale for y-axis
  labs(
    x = "Age (years)", 
    y = "Weight (g)", 
    color = "Species", 
    title = "mizerSeasonal NONSeasonal Growth Curves"  # Title
  ) +  # Labels and title
  theme_minimal() +  # Clean theme
  theme(
    legend.position = "none",  # Hide the legend
    plot.title = element_text(hjust = 0.5, size = 16)  # Center the title and adjust font size
  ) +
  facet_wrap(~ Species, ncol = 3)  # Facet by Species (12 species, arranged in 3 columns)













library(reshape2)
library(ggplot2)
library(dplyr)

# Define the custom order for species
species_order <- c("sprat", "sandeel", "N.pout", 
                   "dab", "herring", "gurnard", 
                   "sole", "whiting", "plaice", 
                   "haddock", "cod", "saithe")

# Reshape and add Scenario labels
test1_s_long <- melt(test1_s, varnames = c("Species", "Age"), value.name = "value") %>%
  mutate(Scenario = "Seasonal")

test1_ns_long <- melt(test1_ns, varnames = c("Species", "Age"), value.name = "value") %>%
  mutate(Scenario = "Non-Seasonal")

# Combine both datasets
combined_data <- bind_rows(test1_s_long, test1_ns_long)

# Convert Species to a factor with the desired order
combined_data$Species <- factor(combined_data$Species, levels = species_order)





library(reshape2)
library(ggplot2)
library(dplyr)

# Define the custom order for species
species_order <- c("Sprat", "Sandeel", "N.pout", 
                   "Dab", "Herring", "Gurnard", 
                   "Sole", "Whiting", "Plaice", 
                   "Haddock", "Cod", "Saithe")

# Reshape and add Scenario labels
test1_s_long <- melt(test1_s, varnames = c("Species", "Age"), value.name = "value") %>%
  mutate(Scenario = "Seasonal")

test1_ns_long <- melt(test1_ns, varnames = c("Species", "Age"), value.name = "value") %>%
  mutate(Scenario = "Non-Seasonal")

# Combine both datasets
combined_data <- bind_rows(test1_s_long, test1_ns_long)

# Convert Species to a factor with the desired order
combined_data$Species <- factor(combined_data$Species, levels = species_order)

# Plot with fixed color, correct linetype mapping, and transparency adjustment
ggplot(combined_data, aes(x = Age, y = value, group = interaction(Species, Scenario), linetype = Scenario)) +
  geom_line(aes(alpha = Scenario), color = "black", size = 1) +  # Set color globally and adjust alpha
  #scale_y_log10() +  

  scale_linetype_manual(values = c("Seasonal" = "solid", "Non-Seasonal" = "dashed")) +  # Correct linetype mapping
  scale_alpha_manual(values = c("Seasonal" = 1, "Non-Seasonal" = 0.5)) +  # Adjust transparency
  labs(
    x = "Time (years)", 
    y = "Weight (g)", 
    linetype = "Scenario",
    alpha = "Scenario",
    title = "mizerSeasonal Growth Curves (Seasonal vs Non-Seasonal)"
  ) +  
  theme_minimal() +  
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    legend.position = "top"
  ) +
  facet_wrap(~ Species, ncol = 3)  # Arranges species across columns









####PRODUCING THE DATTA GROWTH CURVES REPLICAS



library(reshape2)
library(ggplot2)
library(dplyr)

# Define the custom order for species
species_order <- c("Sprat", "Sandeel", "N.pout", 
                   "Dab", "Herring", "Gurnard", 
                   "Sole", "Whiting", "Plaice", 
                   "Haddock", "Cod", "Saithe")

# Reshape and add Scenario labels
test1_s_long <- melt(test1_s, varnames = c("Species", "Age"), value.name = "value") %>%
  mutate(Scenario = "Seasonal")

test1_ns_long <- melt(test1_ns, varnames = c("Species", "Age"), value.name = "value") %>%
  mutate(Scenario = "Non-Seasonal")

# Combine both datasets
combined_data <- bind_rows(test1_s_long, test1_ns_long)

# Convert Species to a factor with the desired order
combined_data$Species <- factor(combined_data$Species, levels = species_order)

# Calculate min and max weight per species
species_limits <- combined_data %>%
  group_by(Species) %>%
  summarise(y_min = min(value, na.rm = TRUE), y_max = max(value, na.rm = TRUE))

# Plot with fixed color, correct linetype mapping, and transparency adjustment
ggplot(combined_data, aes(x = Age, y = value, group = interaction(Species, Scenario), linetype = Scenario)) +
  geom_line(aes(alpha = Scenario), color = "black", size = 1) +  # Set color globally and adjust alpha
  scale_linetype_manual(values = c("Seasonal" = "solid", "Non-Seasonal" = "dashed")) +  # Correct linetype mapping
  scale_alpha_manual(values = c("Seasonal" = 1, "Non-Seasonal" = 0.5)) +  # Adjust transparency
  labs(
    x = "Time (years)", 
    y = "Weight (g)", 
    linetype = "Scenario",
    alpha = "Scenario",
    title = "mizerSeasonal Growth Curves (Seasonal vs Non-Seasonal)"
  ) +  
  theme_minimal() +  
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    legend.position = "top"
  ) +
  facet_wrap(~ Species, ncol = 3, scales = "free_y")  # Allows unique y-axis limits for each species

