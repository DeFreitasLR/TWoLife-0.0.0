# test-parameter-validation.R
# Tests for parameter validation and coercion

test_that("parameters are coerced to integer", {
  # Non-integer values should be accepted and coerced internally by C++
  landscape <- create_test_landscape(size = 10)
  
  result <- twolife_simulation(
    landscape_params = list(habitat = landscape),
    individual_params = list(
      initial_population_size = 5.9,  # Will be coerced to integer internally
      base_birth_rate = 0.35,
      base_mortality_rate = 0.25
    ),
    simulation_params = list(max_events = 100.1)  # Will be coerced internally
  )
  
  expect_s3_class(result, "twolife_result")
  
  # Check that parameters are stored and numeric
  expect_true(is.numeric(result$parameters$individual$initial_population_size))
  expect_equal(result$parameters$individual$initial_population_size, 5.9)
  
  # The key test: simulation should run successfully with non-integer inputs
  # (coercion happens internally in C++ code)
  expect_true("summary" %in% names(result))
  expect_true(is.numeric(result$summary$final_population_size))
})

test_that("invalid landscape parameters are caught", {
  # Missing habitat
  expect_error(
    twolife_simulation(
      landscape_params = list(),
      individual_params = list(initial_population_size = 10)
    ),
    "habitat matrix is required"
  )
  
  # Non-matrix habitat
  expect_error(
    twolife_simulation(
      landscape_params = list(habitat = c(1, 2, 3, 4)),
      individual_params = list(initial_population_size = 10)
    ),
    "must be a matrix"
  )
  
  # Empty matrix
  expect_error(
    twolife_simulation(
      landscape_params = list(habitat = matrix(nrow = 0, ncol = 0)),
      individual_params = list(initial_population_size = 10)
    ),
    "positive dimensions"
  )
})

test_that("demographic parameters are validated", {
  landscape <- create_test_landscape(size = 10)
  
  # Birth rate must be > mortality rate
  expect_error(
    twolife_simulation(
      landscape_params = list(habitat = landscape),
      individual_params = list(
        initial_population_size = 10,
        base_birth_rate = 0.2,
        base_mortality_rate = 0.3
      )
    ),
    "base_birth_rate must be > base_mortality_rate"
  )
})

test_that("history_detail parameter is validated", {
  landscape <- create_test_landscape(size = 10)
  
  # Valid values
  expect_error(
    twolife_simulation(
      landscape_params = list(habitat = landscape),
      individual_params = list(initial_population_size = 5),
      simulation_params = list(max_events = 50),
      history_detail = "minimal"
    ),
    NA
  )
  
  expect_error(
    twolife_simulation(
      landscape_params = list(habitat = landscape),
      individual_params = list(initial_population_size = 5),
      simulation_params = list(max_events = 50),
      history_detail = "standard"
    ),
    NA
  )
  
  expect_error(
    twolife_simulation(
      landscape_params = list(habitat = landscape),
      individual_params = list(initial_population_size = 5),
      simulation_params = list(max_events = 50),
      history_detail = "full"
    ),
    NA
  )
  
  # Invalid value
  expect_error(
    twolife_simulation(
      landscape_params = list(habitat = landscape),
      individual_params = list(initial_population_size = 5),
      simulation_params = list(max_events = 50),
      history_detail = "invalid"
    ),
    "history_detail must be one of"
  )
})

test_that("master_seed parameter is validated", {
  landscape <- create_test_landscape(size = 10)
  
  # Valid seed
  expect_error(
    twolife_simulation(
      landscape_params = list(habitat = landscape),
      individual_params = list(initial_population_size = 5),
      simulation_params = list(max_events = 50),
      master_seed = 12345
    ),
    NA
  )
  
  # Invalid seeds
  expect_error(
    twolife_simulation(
      landscape_params = list(habitat = landscape),
      individual_params = list(initial_population_size = 5),
      simulation_params = list(max_events = 50),
      master_seed = "not_a_number"
    ),
    "must be a single integer"
  )
  
  expect_error(
    twolife_simulation(
      landscape_params = list(habitat = landscape),
      individual_params = list(initial_population_size = 5),
      simulation_params = list(max_events = 50),
      master_seed = c(1, 2, 3)
    ),
    "must be a single integer"
  )
})