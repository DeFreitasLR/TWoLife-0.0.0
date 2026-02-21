# test-parameter-validation.R
# Tests for parameter validation — canonical landscape structure throughout

test_that("numeric parameters are accepted", {
  landscape <- create_continuous_test_landscape()

  # Non-integer values should be accepted (internally handled by C++)
  result <- twolife_simulation(
    landscape_params = landscape,
    individual_params = list(
      initial_population_size = 1000.9,  # Will be handled internally
      neighbor_radius         = 2.0,
      vision_angle            = pi,
      step_length             = 5.0,
      base_dispersal_rate     = 0.4,
      base_birth_rate         = 0.6,
      base_mortality_rate     = 0.2,
      birth_density_slope     = 0.02,
      mortality_density_slope = 0.02
    ),
    simulation_params = list(max_events = 100.5),  # Will be handled internally
    master_seed = 45
  )

  expect_s3_class(result, "twolife_result")
  expect_true(is.numeric(result$parameters$individual$initial_population_size))
  expect_true("summary" %in% names(result))
  expect_true(is.numeric(result$summary$final_population_size))
})

test_that("invalid landscape parameters are caught", {
  # Missing habitat
  expect_error(
    twolife_simulation(
      landscape_params  = list(),
      individual_params = list(initial_population_size = 1000)
    ),
    "habitat matrix is required"
  )

  # Non-matrix habitat
  expect_error(
    twolife_simulation(
      landscape_params  = list(habitat = c(1, 2, 3, 4)),
      individual_params = list(initial_population_size = 1000)
    ),
    "must be a matrix"
  )

  # Empty matrix
  expect_error(
    twolife_simulation(
      landscape_params  = list(habitat = matrix(nrow = 0, ncol = 0)),
      individual_params = list(initial_population_size = 1000)
    ),
    "positive dimensions"
  )

  # Non-finite values
  landscape_bad <- matrix(c(1, 2, NA, 4), nrow = 2, ncol = 2)
  expect_error(
    twolife_simulation(
      landscape_params  = list(habitat = landscape_bad),
      individual_params = list(initial_population_size = 1000)
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
        initial_population_size = 1000,
        neighbor_radius         = 2.0,
        vision_angle            = pi,
        step_length             = 5.0,
        base_dispersal_rate     = 0.4,
        base_birth_rate         = 0.2,
        base_mortality_rate     = 0.3,
        birth_density_slope     = 0.02,
        mortality_density_slope = 0.02
      ),
      simulation_params = list(max_events = 100),
      master_seed = 45
    ),
    "base_birth_rate must be > base_mortality_rate"
  )

  # Negative birth rate should error
  expect_error(
    twolife_simulation(
      landscape_params = landscape,
      individual_params = list(
        initial_population_size = 1000,
        neighbor_radius         = 2.0,
        vision_angle            = pi,
        step_length             = 5.0,
        base_dispersal_rate     = 0.4,
        base_birth_rate         = -0.1,
        base_mortality_rate     = 0.2,
        birth_density_slope     = 0.02,
        mortality_density_slope = 0.02
      ),
      simulation_params = list(max_events = 100),
      master_seed = 45
    )
  )
})

test_that("history_detail parameter is validated", {
  landscape <- create_continuous_test_landscape()

  base_individual <- list(
    initial_population_size = 1000,
    neighbor_radius         = 2.0,
    vision_angle            = pi,
    step_length             = 5.0,
    base_dispersal_rate     = 0.4,
    base_birth_rate         = 0.6,
    base_mortality_rate     = 0.2,
    birth_density_slope     = 0.02,
    mortality_density_slope = 0.02
  )

  # Valid values should work
  expect_error(
    twolife_simulation(
      landscape_params  = landscape,
      individual_params = base_individual,
      simulation_params = list(max_events = 100),
      master_seed    = 45,
      history_detail = "minimal"
    ),
    NA
  )

  expect_error(
    twolife_simulation(
      landscape_params  = landscape,
      individual_params = base_individual,
      simulation_params = list(max_events = 100),
      master_seed    = 45,
      history_detail = "standard"
    ),
    NA
  )

  expect_error(
    twolife_simulation(
      landscape_params  = landscape,
      individual_params = base_individual,
      simulation_params = list(max_events = 100),
      master_seed    = 45,
      history_detail = "full"
    ),
    NA
  )

  # Invalid value should error
  expect_error(
    twolife_simulation(
      landscape_params  = landscape,
      individual_params = base_individual,
      simulation_params = list(max_events = 100),
      master_seed    = 45,
      history_detail = "invalid"
    ),
    "history_detail must be one of"
  )
})

test_that("master_seed parameter is validated", {
  landscape <- create_continuous_test_landscape()

  base_individual <- list(
    initial_population_size = 1000,
    neighbor_radius         = 2.0,
    vision_angle            = pi,
    step_length             = 5.0,
    base_dispersal_rate     = 0.4,
    base_birth_rate         = 0.6,
    base_mortality_rate     = 0.2,
    birth_density_slope     = 0.02,
    mortality_density_slope = 0.02
  )

  # Valid seed should work
  expect_error(
    twolife_simulation(
      landscape_params  = landscape,
      individual_params = base_individual,
      simulation_params = list(max_events = 100),
      master_seed = 12345
    ),
    NA
  )

  # NULL seed should work (non-deterministic)
  expect_error(
    twolife_simulation(
      landscape_params  = landscape,
      individual_params = base_individual,
      simulation_params = list(max_events = 100),
      master_seed = NULL
    ),
    NA
  )

  # Invalid seeds should error
  expect_error(
    twolife_simulation(
      landscape_params  = landscape,
      individual_params = base_individual,
      simulation_params = list(max_events = 100),
      master_seed = "not_a_number"
    ),
    "must be a single integer"
  )

  expect_error(
    twolife_simulation(
      landscape_params  = landscape,
      individual_params = base_individual,
      simulation_params = list(max_events = 100),
      master_seed = c(1, 2, 3)
    ),
    "must be a single integer"
  )

  expect_error(
    twolife_simulation(
      landscape_params  = landscape,
      individual_params = base_individual,
      simulation_params = list(max_events = 100),
      master_seed = 3.14
    ),
    "must be a single integer"
  )
})

test_that("genetic parameter dimensions are validated", {
  landscape <- create_continuous_test_landscape()

  base_individual <- list(
    initial_population_size = 1000,
    neighbor_radius         = 2.0,
    vision_angle            = pi,
    step_length             = 5.0,
    base_dispersal_rate     = 0.4,
    base_birth_rate         = 0.6,
    base_mortality_rate     = 0.2,
    birth_density_slope     = 0.02,
    mortality_density_slope = 0.02
  )

  # Mismatched genotype_means length should error (not 1 or 1000)
  expect_error(
    twolife_simulation(
      landscape_params  = landscape,
      individual_params = base_individual,
      genetic_params = list(
        genotype_means = rep(50, 30),  # Wrong length
        genotype_sds   = 0.4
      ),
      simulation_params = list(max_events = 100),
      master_seed = 49
    ),
    "must have length 1 or"
  )

  # Correct dimensions should work
  expect_error(
    twolife_simulation(
      landscape_params  = landscape,
      individual_params = base_individual,
      genetic_params = list(
        genotype_means = seq(0, 100, length.out = 1000),
        genotype_sds   = rep(0.4, 1000)
      ),
      simulation_params = list(max_events = 100),
      master_seed = 49
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
      individual_params = list(
        initial_population_size = 0,
        neighbor_radius         = 2.0,
        vision_angle            = pi,
        step_length             = 5.0,
        base_dispersal_rate     = 0.4,
        base_birth_rate         = 0.6,
        base_mortality_rate     = 0.2,
        birth_density_slope     = 0.02,
        mortality_density_slope = 0.02
      ),
      simulation_params = list(max_events = 100),
      master_seed = 45
    )
  )

  # Negative population should error
  expect_error(
    twolife_simulation(
      landscape_params = landscape,
      individual_params = list(
        initial_population_size = -10,
        neighbor_radius         = 2.0,
        vision_angle            = pi,
        step_length             = 5.0,
        base_dispersal_rate     = 0.4,
        base_birth_rate         = 0.6,
        base_mortality_rate     = 0.2,
        birth_density_slope     = 0.02,
        mortality_density_slope = 0.02
      ),
      simulation_params = list(max_events = 100),
      master_seed = 45
    )
  )

  # Positive population should work
  expect_error(
    twolife_simulation(
      landscape_params = landscape,
      individual_params = list(
        initial_population_size = 1000,
        neighbor_radius         = 2.0,
        vision_angle            = pi,
        step_length             = 5.0,
        base_dispersal_rate     = 0.4,
        base_birth_rate         = 0.6,
        base_mortality_rate     = 0.2,
        birth_density_slope     = 0.02,
        mortality_density_slope = 0.02
      ),
      simulation_params = list(max_events = 100),
      master_seed = 45
    ),
    NA
  )
})

test_that("max_events must be positive", {
  landscape <- create_continuous_test_landscape()

  base_individual <- list(
    initial_population_size = 1000,
    neighbor_radius         = 2.0,
    vision_angle            = pi,
    step_length             = 5.0,
    base_dispersal_rate     = 0.4,
    base_birth_rate         = 0.6,
    base_mortality_rate     = 0.2,
    birth_density_slope     = 0.02,
    mortality_density_slope = 0.02
  )

  # Single event should complete
  result <- twolife_simulation(
    landscape_params  = landscape,
    individual_params = base_individual,
    simulation_params = list(max_events = 1),
    master_seed = 45
  )
  expect_s3_class(result, "twolife_result")

  # Positive events should work
  expect_error(
    twolife_simulation(
      landscape_params  = landscape,
      individual_params = base_individual,
      simulation_params = list(max_events = 100),
      master_seed = 45
    ),
    NA
  )
})

test_that("landscape parameters have defaults", {
  landscape <- create_continuous_test_landscape()

  # Should work with only habitat provided (all other keys get defaults)
  result <- twolife_simulation(
    landscape_params = list(habitat = landscape$habitat),
    individual_params = list(
      initial_population_size = 1000,
      neighbor_radius         = 2.0,
      vision_angle            = pi,
      step_length             = 5.0,
      base_dispersal_rate     = 0.4,
      base_birth_rate         = 0.6,
      base_mortality_rate     = 0.2,
      birth_density_slope     = 0.02,
      mortality_density_slope = 0.02
    ),
    simulation_params = list(max_events = 100),
    master_seed = 45
  )

  expect_s3_class(result, "twolife_result")

  # Check defaults were applied
  expect_true(!is.null(result$parameters$landscape$cell_size))
  expect_true(!is.null(result$parameters$landscape$boundary_condition))
})
