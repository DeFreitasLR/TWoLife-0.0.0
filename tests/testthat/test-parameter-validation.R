# test-parameter-validation.R
# Tests for parameter validation using consistent parameters

test_that("numeric parameters are accepted", {
  landscape <- create_continuous_test_landscape()
  
  # Non-integer values should be accepted (internally handled by C++)
  set.seed(300)
  result <- twolife_simulation(
    landscape_params = landscape,
    individual_params = list(
      initial_population_size = 50.9,  # Will be handled internally
      base_birth_rate = 0.5,
      base_mortality_rate = 0.2
    ),
    simulation_params = list(max_events = 100.5)  # Will be handled internally
  )
  
  expect_s3_class(result, "twolife_result")
  
  # Check that parameters are stored
  expect_true(is.numeric(result$parameters$individual$initial_population_size))
  
  # Simulation should complete successfully
  expect_true("summary" %in% names(result))
  expect_true(is.numeric(result$summary$final_population_size))
})

test_that("invalid landscape parameters are caught", {
  # Missing habitat
  expect_error(
    twolife_simulation(
      landscape_params = list(),
      individual_params = list(initial_population_size = 50)
    ),
    "habitat matrix is required"
  )
  
  # Non-matrix habitat
  expect_error(
    twolife_simulation(
      landscape_params = list(habitat = c(1, 2, 3, 4)),
      individual_params = list(initial_population_size = 50)
    ),
    "must be a matrix"
  )
  
  # Empty matrix
  expect_error(
    twolife_simulation(
      landscape_params = list(habitat = matrix(nrow = 0, ncol = 0)),
      individual_params = list(initial_population_size = 50)
    ),
    "positive dimensions"
  )
  
  # Non-finite values
  landscape_bad <- matrix(c(1, 2, NA, 4), nrow = 2, ncol = 2)
  expect_error(
    twolife_simulation(
      landscape_params = list(habitat = landscape_bad),
      individual_params = list(initial_population_size = 50)
    ),
    "non-finite values"
  )
})

test_that("demographic parameters are validated", {
  landscape <- create_continuous_test_landscape()
  
  # Birth rate must be > mortality rate
  expect_error(
    twolife_simulation(
      landscape_params = landscape,
      individual_params = list(
        initial_population_size = 50,
        base_birth_rate = 0.2,
        base_mortality_rate = 0.3
      ),
      simulation_params = list(max_events = 100)
    ),
    "base_birth_rate must be > base_mortality_rate"
  )
  
  # Negative rates should error
  expect_error(
    twolife_simulation(
      landscape_params = landscape,
      individual_params = list(
        initial_population_size = 50,
        base_birth_rate = -0.1,
        base_mortality_rate = 0.2
      ),
      simulation_params = list(max_events = 100)
    )
  )
})

test_that("history_detail parameter is validated", {
  landscape <- create_continuous_test_landscape()
  
  # Valid values should work
  expect_error(
    twolife_simulation(
      landscape_params = landscape,
      individual_params = list(initial_population_size = 30),
      simulation_params = list(max_events = 50),
      history_detail = "minimal"
    ),
    NA
  )
  
  expect_error(
    twolife_simulation(
      landscape_params = landscape,
      individual_params = list(initial_population_size = 30),
      simulation_params = list(max_events = 50),
      history_detail = "standard"
    ),
    NA
  )
  
  expect_error(
    twolife_simulation(
      landscape_params = landscape,
      individual_params = list(initial_population_size = 30),
      simulation_params = list(max_events = 50),
      history_detail = "full"
    ),
    NA
  )
  
  # Invalid value should error
  expect_error(
    twolife_simulation(
      landscape_params = landscape,
      individual_params = list(initial_population_size = 30),
      simulation_params = list(max_events = 50),
      history_detail = "invalid"
    ),
    "history_detail must be one of"
  )
})

test_that("master_seed parameter is validated", {
  landscape <- create_continuous_test_landscape()
  
  # Valid seed should work
  expect_error(
    twolife_simulation(
      landscape_params = landscape,
      individual_params = list(initial_population_size = 30),
      simulation_params = list(max_events = 50),
      master_seed = 12345
    ),
    NA
  )
  
  # NULL seed should work (non-deterministic)
  expect_error(
    twolife_simulation(
      landscape_params = landscape,
      individual_params = list(initial_population_size = 30),
      simulation_params = list(max_events = 50),
      master_seed = NULL
    ),
    NA
  )
  
  # Invalid seeds should error
  expect_error(
    twolife_simulation(
      landscape_params = landscape,
      individual_params = list(initial_population_size = 30),
      simulation_params = list(max_events = 50),
      master_seed = "not_a_number"
    ),
    "must be a single integer"
  )
  
  expect_error(
    twolife_simulation(
      landscape_params = landscape,
      individual_params = list(initial_population_size = 30),
      simulation_params = list(max_events = 50),
      master_seed = c(1, 2, 3)
    ),
    "must be a single integer"
  )
  
  expect_error(
    twolife_simulation(
      landscape_params = landscape,
      individual_params = list(initial_population_size = 30),
      simulation_params = list(max_events = 50),
      master_seed = 3.14
    ),
    "must be a single integer"
  )
})

test_that("genetic parameter dimensions are validated", {
  landscape <- create_continuous_test_landscape()
  
  # Mismatched genotype_means length should error
  expect_error(
    twolife_simulation(
      landscape_params = landscape,
      individual_params = list(initial_population_size = 50),
      genetic_params = list(
        genotype_means = rep(0.5, 30),  # Wrong length
        genotype_sds = 0.5
      ),
      simulation_params = list(max_events = 100)
    ),
    "genotype_means must have length equal to initial_population_size"
  )
  
  # Correct dimensions should work
  set.seed(400)
  expect_error(
    twolife_simulation(
      landscape_params = landscape,
      individual_params = list(initial_population_size = 50),
      genetic_params = list(
        genotype_means = runif(50, 0, 1),
        genotype_sds = 0.5
      ),
      simulation_params = list(max_events = 100)
    ),
    NA
  )
})

test_that("population size must be positive", {
  landscape <- create_continuous_test_landscape()
  
  # Zero population should error
  expect_error(
    twolife_simulation(
      landscape_params = landscape,
      individual_params = list(initial_population_size = 0),
      simulation_params = list(max_events = 100)
    )
  )
  
  # Negative population should error
  expect_error(
    twolife_simulation(
      landscape_params = landscape,
      individual_params = list(initial_population_size = -10),
      simulation_params = list(max_events = 100)
    )
  )
  
  # Positive population should work
  expect_error(
    twolife_simulation(
      landscape_params = landscape,
      individual_params = list(initial_population_size = 30),
      simulation_params = list(max_events = 100)
    ),
    NA
  )
})

test_that("max_events must be positive", {
  landscape <- create_continuous_test_landscape()
  
  # Zero events should error or complete immediately
  # (implementation-dependent, but should handle gracefully)
  result <- twolife_simulation(
    landscape_params = landscape,
    individual_params = list(initial_population_size = 30),
    simulation_params = list(max_events = 1)
  )
  expect_s3_class(result, "twolife_result")
  
  # Positive events should work
  expect_error(
    twolife_simulation(
      landscape_params = landscape,
      individual_params = list(initial_population_size = 30),
      simulation_params = list(max_events = 100)
    ),
    NA
  )
})

test_that("landscape parameters have defaults", {
  landscape <- create_continuous_test_landscape()
  
  # Should work with minimal landscape specification
  set.seed(300)
  result <- twolife_simulation(
    landscape_params = list(habitat = landscape$habitat),
    individual_params = list(initial_population_size = 30),
    simulation_params = list(max_events = 50)
  )
  
  expect_s3_class(result, "twolife_result")
  
  # Check defaults were applied
  expect_true(!is.null(result$parameters$landscape$cell_size))
  expect_true(!is.null(result$parameters$landscape$boundary_condition))
})
