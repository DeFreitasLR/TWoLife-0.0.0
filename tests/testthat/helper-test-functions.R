# helper-test-functions.R
# Test helper functions for TWoLife package
# Updated to match CRAN-compliant examples

#' Create standard binary test landscape
#' 
#' Creates consistent binary landscape using seed 100
#' 
#' @return A landscape_params list
create_binary_test_landscape <- function() {
  set.seed(100)
  create_fractal_landscape(
    cells_per_row = 10,
    fractality = 0.5,
    habitat_proportion = 0.4,
    return_as_landscape_params = TRUE
  )
}

#' Create standard continuous test landscape
#' 
#' Creates consistent continuous landscape using seed 200
#' 
#' @return A landscape_params list
create_continuous_test_landscape <- function() {
  set.seed(200)
  create_fractal_landscape(
    cells_per_row = 10,
    fractality = 0.5,
    min_value = 0,
    max_value = 1,
    return_as_landscape_params = TRUE
  )
}

#' Run simple test simulation
#' 
#' Wrapper for consistent test simulations
#' Uses seed 300 for simple simulations
#' 
#' @param landscape Landscape params (default: continuous test landscape)
#' @param population Integer. Initial population size (default: 50)
#' @param events Integer. Maximum events (default: 100)
#' @param seed Integer. Random seed (default: 300)
#' @param ... Additional parameters passed to twolife_simulation
#' 
#' @return A twolife_result object
run_simple_test_simulation <- function(landscape = NULL,
                                       population = 50,
                                       events = 100,
                                       seed = 300,
                                       ...) {
  
  # Use default continuous landscape if not provided
  if (is.null(landscape)) {
    landscape <- create_continuous_test_landscape()
  }
  
  # Run simulation
  set.seed(seed)
  twolife_simulation(
    landscape_params = landscape,
    individual_params = list(
      initial_population_size = population,
      base_birth_rate = 0.5,
      base_mortality_rate = 0.2
    ),
    simulation_params = list(
      max_events = events
    ),
    ...
  )
}

#' Run genetic test simulation
#' 
#' Wrapper for simulations with genetic variation
#' Uses seed 400 for genetic simulations
#' 
#' @param landscape Landscape params (default: continuous test landscape)
#' @param population Integer. Initial population size (default: 50)
#' @param events Integer. Maximum events (default: 100)
#' @param seed Integer. Random seed (default: 400)
#' @param ... Additional parameters passed to twolife_simulation
#' 
#' @return A twolife_result object
run_genetic_test_simulation <- function(landscape = NULL,
                                        population = 50,
                                        events = 100,
                                        seed = 400,
                                        ...) {
  
  # Use default continuous landscape if not provided
  if (is.null(landscape)) {
    landscape <- create_continuous_test_landscape()
  }
  
  # Run simulation with genetic variation
  set.seed(seed)
  twolife_simulation(
    landscape_params = landscape,
    individual_params = list(
      initial_population_size = population,
      base_birth_rate = 0.5,
      base_mortality_rate = 0.2,
      matrix_mortality_multiplier = 3.0
    ),
    genetic_params = list(
      genotype_means = runif(population, min = 0, max = 1),
      genotype_sds = 0.5,
      sampling_points = 10,
      habitat_selection_temperatures = 0.1
    ),
    simulation_params = list(
      max_events = events
    ),
    ...
  )
}
