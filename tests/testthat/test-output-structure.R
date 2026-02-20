# test-output-structure.R
# Tests for simulation output structure

test_that("output has expected fields", {
  result <- run_simple_test_simulation(steps = 100, n = 10)
  
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

test_that("survivors data frame has correct structure when population survives", {
  # Run with favorable parameters
  landscape <- matrix(1, nrow = 10, ncol = 10)  # All habitat
  
  result <- twolife_simulation(
    landscape_params = list(habitat = landscape),
    individual_params = list(
      initial_population_size = 20,
      base_birth_rate = 0.4,
      base_mortality_rate = 0.2
    ),
    simulation_params = list(max_events = 200)
  )
  
  # If population survived
  if (result$summary$final_population_size > 0) {
    expect_s3_class(result$survivors, "data.frame")
    expect_equal(nrow(result$survivors), result$summary$final_population_size)
    
    # Check column types
    expect_type(result$survivors$id, "integer")
    expect_type(result$survivors$x, "double")
    expect_type(result$survivors$y, "double")
    expect_type(result$survivors$genotype, "double")
    expect_type(result$survivors$phenotype, "double")
    expect_type(result$survivors$width, "double")
  }
})

test_that("events vectors have consistent lengths", {
  result <- run_simple_test_simulation(steps = 100, n = 10)
  
  n_events <- length(result$events$times)
  
  expect_equal(length(result$events$types), n_events)
  expect_equal(length(result$events$individual_ids), n_events)
  
  # For standard history detail
  if ("x_coordinates" %in% names(result$events)) {
    expect_equal(length(result$events$x_coordinates), n_events)
    expect_equal(length(result$events$y_coordinates), n_events)
    expect_equal(length(result$events$genotypes), n_events)
  }
})