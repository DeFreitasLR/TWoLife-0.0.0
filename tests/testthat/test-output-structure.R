# test-output-structure.R
# Tests for simulation output structure using consistent parameters

test_that("output has expected fields", {
  # Run simple test simulation
  result <- run_simple_test_simulation(
    population = 50,
    events = 100,
    seed = 300
  )
  
  # Check main structure
  expect_named(result, c("summary", "survivors", "spatial", "events", "parameters"))
  
  # Check summary structure
  expect_named(result$summary, c("final_population_size", "total_events", "duration", "status"))
  
  # Check survivors structure (if any survived)
  if (result$summary$final_population_size > 0) {
    expect_true(is.data.frame(result$survivors))
    expect_true(all(c("id", "x", "y", "genotype", "phenotype", "width") %in% names(result$survivors)))
  }
  
  # Check events structure
  expect_true("times" %in% names(result$events))
  expect_true("types" %in% names(result$events))
  expect_true("individual_ids" %in% names(result$events))
  
  # Check parameters structure
  expect_true(all(c("landscape", "individual", "genetic", "simulation") %in% names(result$parameters)))
})

test_that("survivors data frame has correct structure", {
  # Run genetic simulation for more interesting survivors
  result <- run_genetic_test_simulation(
    population = 50,
    events = 100,
    seed = 400
  )
  
  # If population survived
  if (result$summary$final_population_size > 0) {
    # Check it's a data frame
    expect_s3_class(result$survivors, "data.frame")
    
    # Check number of rows matches population
    expect_equal(nrow(result$survivors), result$summary$final_population_size)
    
    # Check column types
    expect_type(result$survivors$id, "integer")
    expect_type(result$survivors$x, "double")
    expect_type(result$survivors$y, "double")
    expect_type(result$survivors$genotype, "double")
    expect_type(result$survivors$phenotype, "double")
    expect_type(result$survivors$width, "double")
    
    # Check reasonable ranges
    landscape <- create_continuous_test_landscape()
    grid_size <- nrow(landscape$habitat)
    
    expect_true(all(result$survivors$x >= 0 & result$survivors$x < grid_size))
    expect_true(all(result$survivors$y >= 0 & result$survivors$y < grid_size))
    expect_true(all(result$survivors$genotype >= 0 & result$survivors$genotype <= 1))
    expect_true(all(result$survivors$phenotype >= 0 & result$survivors$phenotype <= 1))
    expect_true(all(result$survivors$width > 0))
  }
})

test_that("events vectors have consistent lengths", {
  # Run simple simulation
  result <- run_simple_test_simulation(
    population = 50,
    events = 100,
    seed = 300
  )
  
  n_events <- length(result$events$times)
  
  # All event vectors should have same length
  expect_equal(length(result$events$types), n_events)
  expect_equal(length(result$events$individual_ids), n_events)
  
  # For standard history detail
  if ("x_coordinates" %in% names(result$events)) {
    expect_equal(length(result$events$x_coordinates), n_events)
    expect_equal(length(result$events$y_coordinates), n_events)
    expect_equal(length(result$events$genotypes), n_events)
  }
  
  # Check that event times are monotonically increasing
  if (n_events > 1) {
    expect_true(all(diff(result$events$times) >= 0))
  }
})

test_that("spatial history has correct structure", {
  landscape <- create_continuous_test_landscape()
  
  # Run with full history to ensure spatial tracking
  set.seed(400)
  result <- twolife_simulation(
    landscape_params = landscape,
    individual_params = list(
      initial_population_size = 30,
      base_birth_rate = 0.5,
      base_mortality_rate = 0.2
    ),
    simulation_params = list(max_events = 50),
    history_detail = "full"
  )
  
  # Check spatial structure exists
  expect_true("x_history" %in% names(result$spatial))
  expect_true("y_history" %in% names(result$spatial))
  expect_true("times" %in% names(result$spatial))
  
  # If there's spatial history, check consistency
  if (length(result$spatial$times) > 0) {
    expect_equal(length(result$spatial$x_history), length(result$spatial$y_history))
    expect_equal(length(result$spatial$x_history), length(result$spatial$times))
  }
})

test_that("parameters are stored correctly", {
  landscape <- create_continuous_test_landscape()
  
  # Run with specific parameters
  set.seed(300)
  result <- twolife_simulation(
    landscape_params = landscape,
    individual_params = list(
      initial_population_size = 50,
      base_birth_rate = 0.5,
      base_mortality_rate = 0.2,
      matrix_mortality_multiplier = 2.5
    ),
    genetic_params = list(
      genotype_means = rep(0.5, 50),
      genotype_sds = 0.4,
      mutation_rates = 0.01,
      plasticities = 0.05
    ),
    simulation_params = list(
      max_events = 100,
      neutral_mode = FALSE
    )
  )
  
  # Check individual parameters
  expect_equal(result$parameters$individual$initial_population_size, 50)
  expect_equal(result$parameters$individual$base_birth_rate, 0.5)
  expect_equal(result$parameters$individual$base_mortality_rate, 0.2)
  expect_equal(result$parameters$individual$matrix_mortality_multiplier, 2.5)
  
  # Check genetic parameters
  expect_equal(length(result$parameters$genetic$genotype_means), 50)
  expect_equal(result$parameters$genetic$genotype_sds, 0.4)
  expect_equal(result$parameters$genetic$mutation_rates, 0.01)
  expect_equal(result$parameters$genetic$plasticities, 0.05)
  
  # Check simulation parameters
  expect_equal(result$parameters$simulation$max_events, 100)
  expect_equal(result$parameters$simulation$neutral_mode, FALSE)
  
  # Check landscape parameters
  expect_true(is.matrix(result$parameters$landscape$habitat))
  expect_equal(nrow(result$parameters$landscape$habitat), 10)
  expect_equal(ncol(result$parameters$landscape$habitat), 10)
})

test_that("summary statistics are reasonable", {
  # Run multiple simulations
  result <- run_simple_test_simulation(
    population = 50,
    events = 100,
    seed = 300
  )
  
  # Check summary values are non-negative
  expect_true(result$summary$final_population_size >= 0)
  expect_true(result$summary$total_events >= 0)
  expect_true(result$summary$duration >= 0)
  
  # Check final population is integer
  expect_equal(result$summary$final_population_size, 
               as.integer(result$summary$final_population_size))
  
  # Check status is valid
  expect_true(result$summary$status %in% c(
    "completed",
    "max_events_reached",
    "population_extinct",
    "max_time_reached"
  ))
  
  # If population survived, check survivors match summary
  if (result$summary$final_population_size > 0) {
    expect_equal(nrow(result$survivors), result$summary$final_population_size)
  }
})
