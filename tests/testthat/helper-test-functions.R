# helper-test-functions.R
# Test helper functions for TWoLife package
# Mirrors the two canonical simulation cases from the package examples.
#
# Canonical cases:
#   Binary   — 5x5 binary landscape (seed 100), N=1000, no genetics, master_seed=45
#   Genetic  — 5x5 continuous landscape (seed 50, 0-100), N=1000, full S5 genetics, master_seed=49


#' Create standard binary test landscape
#'
#' 5x5 binary landscape, seed 100, habitat_proportion=0.4.
#' Returns a fully explicit landscape_params list (raw matrix + all keys).
#'
#' @return A landscape_params list
create_binary_test_landscape <- function() {
  set.seed(100)
  binary_matrix <- create_fractal_landscape(
    cells_per_row      = 5,
    fractality         = 0.5,
    habitat_proportion = 0.4
  )
  list(
    habitat                     = binary_matrix,
    cell_size                   = 1.0,
    boundary_condition          = 1,
    density_type                = 1,
    matrix_mortality_multiplier = 2.0,
    matrix_dispersal_multiplier = 0.5
  )
}


#' Create standard continuous test landscape
#'
#' 5x5 continuous landscape, seed 50, range [0, 100].
#' Returns a fully explicit landscape_params list (raw matrix + all keys).
#'
#' @return A landscape_params list
create_continuous_test_landscape <- function() {
  set.seed(50)
  continuous_matrix <- create_fractal_landscape(
    cells_per_row = 5,
    fractality    = 0.5,
    min_value     = 0.0,
    max_value     = 100.0
  )
  list(
    habitat                     = continuous_matrix,
    cell_size                   = 1.0,
    boundary_condition          = 1,
    density_type                = 1,
    matrix_mortality_multiplier = 2.0,
    matrix_dispersal_multiplier = 0.5
  )
}


#' Run simple test simulation
#'
#' Binary landscape, N=1000, no genetics. Mirror of the canonical `result` case.
#'
#' @param landscape   landscape_params list (default: canonical binary landscape)
#' @param population  Integer. Initial population size (default: 1000)
#' @param events      Integer. Maximum events (default: 100)
#' @param master_seed Integer. Simulation seed (default: 45)
#' @param ...         Additional parameters passed to twolife_simulation
#'
#' @return A twolife_result object
run_simple_test_simulation <- function(landscape   = NULL,
                                       population  = 1000,
                                       events      = 100,
                                       master_seed = 45,
                                       ...) {
  if (is.null(landscape)) {
    landscape <- create_binary_test_landscape()
  }

  twolife_simulation(
    landscape_params = landscape,
    individual_params = list(
      initial_population_size = population,
      neighbor_radius         = 2.0,
      vision_angle            = pi,
      step_length             = 5.0,
      base_dispersal_rate     = 0.4,
      base_birth_rate         = 0.6,
      base_mortality_rate     = 0.2,
      birth_density_slope     = 0.02,
      mortality_density_slope = 0.02
    ),
    simulation_params = list(max_events = events),
    master_seed = master_seed,
    ...
  )
}


#' Run genetic test simulation
#'
#' Continuous landscape (0-100), N=1000, full S5 genetics.
#' Mirror of the canonical `result_genetic` case.
#'
#' @param landscape   landscape_params list (default: canonical continuous landscape)
#' @param population  Integer. Initial population size (default: 1000)
#' @param events      Integer. Maximum events (default: 100)
#' @param master_seed Integer. Simulation seed (default: 49)
#' @param ...         Additional parameters passed to twolife_simulation
#'
#' @return A twolife_result object
run_genetic_test_simulation <- function(landscape   = NULL,
                                        population  = 1000,
                                        events      = 100,
                                        master_seed = 49,
                                        ...) {
  if (is.null(landscape)) {
    landscape <- create_continuous_test_landscape()
  }

  initial_genotypes <- seq(0, 100, length.out = population)

  twolife_simulation(
    landscape_params = landscape,
    individual_params = list(
      initial_population_size = population,
      neighbor_radius         = 2.0,
      vision_angle            = pi,
      step_length             = 5.0,
      base_dispersal_rate     = 0.4,
      base_birth_rate         = 0.6,
      base_mortality_rate     = 0.2,
      birth_density_slope     = 0.02,
      mortality_density_slope = 0.02
    ),
    genetic_params = list(
      genotype_means                 = initial_genotypes,
      genotype_sds                   = rep(0.4,   population),
      mutation_rates                 = rep(0.001, population),
      plasticities                   = rep(0.001, population),
      sampling_points                = rep(10,    population),
      habitat_selection_temperatures = rep(0.1,   population)
    ),
    simulation_params = list(max_events = events),
    master_seed = master_seed,
    ...
  )
}
