# COMPLETE FIXED TEST FILE
# Replace your entire tests/testthat/test-basic-functionality.R with this

test_that("twolife_simulation runs and returns expected structure", {
  # Create simple landscape
  landscape <- create_fractal_landscape(
    cells_per_row = 5,
    fractality = 0.5,
    habitat_proportion = 0.6,
    return_as_landscape_params = TRUE
  )
  
  # Run basic simulation
  result <- twolife_simulation(
    landscape_params = landscape,
    individual_params = list(
      initial_population_size = 15,
      base_birth_rate = 0.4,
      base_mortality_rate = 0.15
    ),
    simulation_params = list(max_events = 150),
    master_seed = 100
  )
  
  # Test structure
  expect_true("twolife_result" %in% class(result))
  expect_true(is.list(result))
  expect_true("summary" %in% names(result))
  expect_true("survivors" %in% names(result))
  expect_true("events" %in% names(result))
  expect_true("parameters" %in% names(result))
  
  # Test summary contents
  expect_true("final_population_size" %in% names(result$summary))
  expect_true("total_events" %in% names(result$summary))
  expect_true("status" %in% names(result$summary))
  expect_true("duration" %in% names(result$summary))
  
  # Test survivors structure
  expect_true(is.list(result$survivors))
  expect_true("x" %in% names(result$survivors))
  expect_true("y" %in% names(result$survivors))
  expect_true("id" %in% names(result$survivors))
})

test_that("Population size can be computed", {
  # Create simple landscape
  landscape <- create_fractal_landscape(
    cells_per_row = 5,
    fractality = 0.5,
    habitat_proportion = 0.6,
    return_as_landscape_params = TRUE
  )
  
  # Run simulation
  result <- twolife_simulation(
    landscape_params = landscape,
    individual_params = list(
      initial_population_size = 15,
      base_birth_rate = 0.4,
      base_mortality_rate = 0.15
    ),
    simulation_params = list(max_events = 150),
    master_seed = 101
  )
  
  # ✅ FIXED: Changed compute_population_size to population_size
  pop_trajectory <- population_size(result)
  
  # Test trajectory structure
  expect_true(is.data.frame(pop_trajectory))
  expect_true("time" %in% names(pop_trajectory))
  expect_true("population_size" %in% names(pop_trajectory))
  expect_true(nrow(pop_trajectory) > 0)
  expect_true(all(pop_trajectory$population_size >= 0))
})

test_that("Simulation handles genetic parameters", {
  landscape <- create_fractal_landscape(
    cells_per_row = 5,
    fractality = 0.5,
    habitat_proportion = 0.6,
    return_as_landscape_params = TRUE
  )
  
  result <- twolife_simulation(
    landscape_params = landscape,
    individual_params = list(initial_population_size = 10),
    genetic_params = list(
      genotype_means = runif(10, 0.3, 0.7),
      genotype_sds = 0.15,
      mutation_rates = 0.02
    ),
    simulation_params = list(max_events = 100),
    master_seed = 102
  )
  
  expect_true("genotype" %in% names(result$survivors))
  expect_true(length(result$survivors$genotype) == result$summary$final_population_size)
})

test_that("Different history detail levels work", {
  landscape <- create_fractal_landscape(
    cells_per_row = 5,
    fractality = 0.5,
    habitat_proportion = 0.6,
    return_as_landscape_params = TRUE
  )
  
  # Minimal history
  result_minimal <- twolife_simulation(
    landscape_params = landscape,
    individual_params = list(initial_population_size = 10),
    simulation_params = list(max_events = 100),
    history_detail = "minimal",
    master_seed = 103
  )
  
  # Full history
  result_full <- twolife_simulation(
    landscape_params = landscape,
    individual_params = list(initial_population_size = 10),
    simulation_params = list(max_events = 100),
    history_detail = "full",
    master_seed = 104
  )
  
  expect_true("twolife_result" %in% class(result_minimal))
  expect_true("twolife_result" %in% class(result_full))
})
