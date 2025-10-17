# helper-test-functions.R
# Test helper functions for TWoLife package

#' Create a simple test landscape
#' 
#' Helper function to create basic landscapes for testing
#' 
#' @param size Integer. Size of square landscape (default: 10)
#' @param habitat_prop Numeric. Proportion of habitat (default: 0.5)
#' 
#' @return A binary matrix
create_test_landscape <- function(size = 10, habitat_prop = 0.5) {
  landscape <- matrix(0, nrow = size, ncol = size)
  n_habitat <- round(size * size * habitat_prop)
  habitat_indices <- sample(size * size, n_habitat)
  landscape[habitat_indices] <- 1
  return(landscape)
}

#' Run a simple test simulation
#' 
#' Wrapper function to run simulations with sensible defaults for testing.
#' This avoids the argument matching issues with the underlying C++ function.
#' 
#' @param steps Integer. Maximum events to simulate
#' @param n Integer. Initial population size
#' @param landscape Matrix. Habitat landscape (if NULL, creates default)
#' @param ... Additional parameters passed to twolife_simulation
#' 
#' @return A twolife_result object
run_simple_test_simulation <- function(steps = 100, 
                                       n = 10, 
                                       landscape = NULL,
                                       ...) {
  
  # Create default landscape if not provided
  if (is.null(landscape)) {
    landscape <- create_test_landscape(size = 10, habitat_prop = 0.5)
  }
  
  # Prepare parameters
  landscape_params <- list(
    habitat = landscape,
    cell_size = 1.0
  )
  
  individual_params <- list(
    initial_population_size = as.integer(n),
    neighbor_radius = 2.0,
    base_birth_rate = 0.35,
    base_mortality_rate = 0.25
  )
  
  simulation_params <- list(
    max_events = as.integer(steps)
  )
  
  # Merge with any additional parameters
  dots <- list(...)
  if ("landscape_params" %in% names(dots)) {
    landscape_params <- modifyList(landscape_params, dots$landscape_params)
    dots$landscape_params <- NULL
  }
  if ("individual_params" %in% names(dots)) {
    individual_params <- modifyList(individual_params, dots$individual_params)
    dots$individual_params <- NULL
  }
  if ("simulation_params" %in% names(dots)) {
    simulation_params <- modifyList(simulation_params, dots$simulation_params)
    dots$simulation_params <- NULL
  }
  
  # Run simulation
  result <- do.call(twolife_simulation, c(
    list(
      landscape_params = landscape_params,
      individual_params = individual_params,
      simulation_params = simulation_params
    ),
    dots
  ))
  
  return(result)
}